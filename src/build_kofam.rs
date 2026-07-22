//! Build the KEGG KOfam reference DB — a pure-Rust provisioner mirroring the
//! layout of the shell-built `db/kofam/`. Unlike PSC there is NO external binary:
//! every step (HTTP download, gzip, tar extraction, and the prokaryotic-profile
//! concatenation) runs in-process here (no curl/gunzip/tar/cat).
//!
//! Produces, under `kofam_dir`:
//!   * `ko_list`            — per-KO adaptive bit-score thresholds + definitions
//!                            (gunzipped from `ko_list.gz`)
//!   * `profiles/`          — the extracted per-KO HMMER3 profiles + `.hal` lists
//!                            (from `profiles.tar.gz`)
//!   * `kofam_prok.hmm`     — the concatenated PROKARYOTIC KOfam profiles, built by
//!                            appending every `profiles/<KO>.hmm` named in
//!                            `profiles/prokaryote.hal` (in list order)
//!
//! This is the KofamScan search method's model DB: a single HMM library that is
//! the concatenation of exactly the prokaryotic subset KEGG defines in
//! `prokaryote.hal`. It is a FULL-tier naming DB (like PSC) — provisioned only
//! under `bactars setup-db --full` / `--only kofam`, never in the default set.
//! Heavy: ~2.4 GB of downloads and ~5-6 GB extracted; run it once.
//!
//! Source (KEGG, verified live 2026-07): `https://www.genome.jp/ftp/db/kofam/`
//! — `ko_list.gz` and `profiles.tar.gz`. The prokaryotic subset is defined by
//! `profiles/prokaryote.hal` (a plain list of `<KO>.hmm` filenames) shipped inside
//! `profiles.tar.gz`; `kofam_prok.hmm` is the raw concatenation of those files.

use crate::util_io;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// KEGG KOfam per-KO thresholds + definitions (gzip, ~0.9 MB).
const KO_LIST_URL: &str = "https://www.genome.jp/ftp/db/kofam/ko_list.gz";
/// KEGG KOfam per-KO HMMER3 profiles + `.hal` subset lists (tar.gz, ~1.5 GB).
const PROFILES_URL: &str = "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz";
/// The subset list (inside `profiles/`) naming the prokaryotic KO profiles.
const PROK_HAL: &str = "prokaryote.hal";

/// Inputs / outputs for a KOfam build.
pub struct KofamBuildConfig {
    /// Output directory (e.g. `<bundle>/kofam`).
    pub kofam_dir: PathBuf,
    /// Rebuild even if `kofam_prok.hmm` already exists.
    pub force: bool,
}

/// Build the KOfam DB end-to-end. Idempotent: skips when `kofam_prok.hmm` exists
/// unless `force`. Individual intermediate steps (ko_list, profiles extraction)
/// are also skipped when their output is already present, so an interrupted build
/// resumes cheaply. Pure-Rust: no external binary is required.
pub fn build_kofam(cfg: &KofamBuildConfig) -> Result<(), String> {
    let dir = &cfg.kofam_dir;
    std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;

    let prok_hmm = dir.join("kofam_prok.hmm");
    if prok_hmm.exists() && !cfg.force {
        eprintln!(
            "[kofam] kofam_prok.hmm already present at {} — skipping build (use --force)",
            prok_hmm.display()
        );
        return Ok(());
    }

    // --- Step 0: ko_list (download ko_list.gz + gunzip) ---
    let ko_list = dir.join("ko_list");
    if !ko_list.exists() || cfg.force {
        let gz = dir.join("ko_list.gz");
        eprintln!("[kofam][0] downloading ko_list.gz (~0.9 MB)");
        util_io::http_download(KO_LIST_URL, &gz)?;
        eprintln!("[kofam][0] gunzip ko_list.gz -> ko_list");
        util_io::gunzip_file(&gz, &ko_list)?;
        if !ko_list.exists() {
            return Err(format!("gunzip did not produce {}", ko_list.display()));
        }
        // Keep ko_list.gz alongside (matches the shell-built db/kofam layout).
    }

    // --- Step 1: profiles/ (download profiles.tar.gz + extract) ---
    let profiles = dir.join("profiles");
    let hal = profiles.join(PROK_HAL);
    if !hal.exists() || cfg.force {
        let tgz = dir.join("profiles.tar.gz");
        if !tgz.exists() || cfg.force {
            eprintln!("[kofam][1] downloading profiles.tar.gz (~1.5 GB)");
            util_io::http_download(PROFILES_URL, &tgz)?;
        }
        // The archive nests everything under a top-level `profiles/`, so extracting
        // into `kofam_dir` produces `kofam_dir/profiles/<KO>.hmm` + the `.hal` lists.
        eprintln!("[kofam][1] extracting profiles.tar.gz -> {}", profiles.display());
        util_io::extract_tar_gz(&tgz, dir)?;
        if !hal.exists() {
            return Err(format!(
                "profiles.tar.gz extracted but {} missing (is this the KOfam profiles archive?)",
                hal.display()
            ));
        }
        // Keep profiles.tar.gz alongside (matches the shell-built db/kofam layout).
    }

    // --- Step 2: build kofam_prok.hmm (concatenate the prok subset in hal order) ---
    eprintln!("[kofam][2] concatenating prokaryotic profiles -> {}", prok_hmm.display());
    let names = read_hal(&hal)?;
    if names.is_empty() {
        return Err(format!("{} listed no profiles", hal.display()));
    }
    let n = concat_profiles(&profiles, &names, &prok_hmm)?;
    eprintln!("[kofam][done] kofam_prok.hmm built with {n} prokaryotic KO profiles at {}", prok_hmm.display());
    Ok(())
}

/// Parse a KEGG `.hal` subset list into the ordered list of profile filenames it
/// names (one `<KO>.hmm` per line). Blank lines and surrounding whitespace are
/// ignored; every kept entry is a bare filename (no path component).
fn parse_hal(text: &str) -> Vec<String> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        // Guard against any path components in a `.hal` line — the concatenator
        // joins these onto the profiles dir, so only a bare filename is safe.
        .map(|l| l.rsplit(['/', '\\']).next().unwrap_or(l).to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

/// Read + parse a `.hal` file from disk into its ordered profile-filename list.
fn read_hal(hal: &Path) -> Result<Vec<String>, String> {
    let text = std::fs::read_to_string(hal).map_err(|e| format!("read {}: {e}", hal.display()))?;
    Ok(parse_hal(&text))
}

/// Concatenate `profiles/<name>` for every `name` in `names` (in order) into
/// `out_path`, streaming each profile so the multi-GB result never loads into
/// memory. Each KOfam profile is a complete `HMMER3/f … //` record, so a raw
/// byte concatenation reproduces the shell-built `kofam_prok.hmm`. Returns the
/// number of profiles concatenated.
fn concat_profiles(profiles_dir: &Path, names: &[String], out_path: &Path) -> Result<u64, String> {
    let out = File::create(out_path).map_err(|e| format!("create {}: {e}", out_path.display()))?;
    let mut w = BufWriter::with_capacity(1 << 20, out);
    let mut buf = vec![0u8; 1 << 16];
    let mut n = 0u64;
    for name in names {
        let src = profiles_dir.join(name);
        let f = File::open(&src).map_err(|e| format!("open profile {}: {e}", src.display()))?;
        let mut r = BufReader::with_capacity(1 << 16, f);
        loop {
            let got = std::io::Read::read(&mut r, &mut buf).map_err(|e| format!("read {}: {e}", src.display()))?;
            if got == 0 {
                break;
            }
            w.write_all(&buf[..got]).map_err(|e| format!("write {}: {e}", out_path.display()))?;
        }
        n += 1;
    }
    w.flush().map_err(|e| format!("flush {}: {e}", out_path.display()))?;
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufRead; // for `.lines()` in the real-data parity test

    #[test]
    fn parse_hal_orders_and_trims() {
        let text = "K00001.hmm\nK00002.hmm\n\n  K00010.hmm  \nK00003.hmm\n";
        let got = parse_hal(text);
        assert_eq!(got, vec!["K00001.hmm", "K00002.hmm", "K00010.hmm", "K00003.hmm"]);
    }

    #[test]
    fn parse_hal_strips_path_components() {
        // Defensive: a `.hal` entry must resolve to a bare filename under profiles/.
        let got = parse_hal("profiles/K00001.hmm\n../evil.hmm\nK00002.hmm\n");
        assert_eq!(got, vec!["K00001.hmm", "evil.hmm", "K00002.hmm"]);
    }

    #[test]
    fn concat_profiles_appends_in_order() {
        let dir = std::env::temp_dir().join(format!("bactars_kofam_{}", std::process::id()));
        let prof = dir.join("profiles");
        std::fs::create_dir_all(&prof).unwrap();
        // Two minimal "profiles" — the concatenator is format-agnostic (raw bytes).
        std::fs::write(prof.join("K1.hmm"), b"HMMER3/f AAA\n//\n").unwrap();
        std::fs::write(prof.join("K2.hmm"), b"HMMER3/f BBB\n//\n").unwrap();
        let out = dir.join("out.hmm");
        // hal order is K2 then K1 — the output must follow the hal order, not
        // lexicographic / directory order.
        let n = concat_profiles(&prof, &["K2.hmm".into(), "K1.hmm".into()], &out).unwrap();
        assert_eq!(n, 2);
        let got = std::fs::read(&out).unwrap();
        assert_eq!(got, b"HMMER3/f BBB\n//\nHMMER3/f AAA\n//\n");
        // Model count = number of HMMER3/ record headers = number of hal entries.
        let models = String::from_utf8(got).unwrap().matches("HMMER3/").count();
        assert_eq!(models, 2);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn concat_profiles_missing_profile_errors() {
        let dir = std::env::temp_dir().join(format!("bactars_kofam_miss_{}", std::process::id()));
        let prof = dir.join("profiles");
        std::fs::create_dir_all(&prof).unwrap();
        let out = dir.join("out.hmm");
        let e = concat_profiles(&prof, &["nope.hmm".into()], &out).unwrap_err();
        assert!(e.contains("open profile"), "unexpected error: {e}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Real-data parity: on a machine with the shell-built `db/kofam`, the number
    /// of profiles our concatenator would emit (lines in `prokaryote.hal`) must
    /// EXACTLY equal the number of `HMMER3/` records in the built `kofam_prok.hmm`.
    /// This proves the prokaryote.hal → kofam_prok.hmm mechanism matches the real
    /// build. Auto-skips where `db/kofam` is absent (CI / fresh checkout).
    #[test]
    fn prok_hal_matches_kofam_prok_hmm() {
        let kofam = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../db/kofam");
        let hal = kofam.join("profiles/prokaryote.hal");
        let prok_hmm = kofam.join("kofam_prok.hmm");
        if !hal.exists() || !prok_hmm.exists() {
            eprintln!("prok_hal_matches: db/kofam absent — skipping");
            return;
        }
        let names = read_hal(&hal).unwrap();
        let n_hal = names.len();
        // Every hal entry is a `<KO>.hmm` filename.
        assert!(names.iter().all(|n| n.ends_with(".hmm")), "hal entry not a .hmm filename");
        // Count HMMER3/ record headers in the built library (streamed).
        let f = File::open(&prok_hmm).unwrap();
        let mut models = 0usize;
        for line in BufReader::new(f).lines().map_while(Result::ok) {
            if line.starts_with("HMMER3/") {
                models += 1;
            }
        }
        assert_eq!(
            n_hal, models,
            "prokaryote.hal lists {n_hal} profiles but kofam_prok.hmm has {models} models"
        );
    }
}
