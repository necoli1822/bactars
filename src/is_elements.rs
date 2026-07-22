//! IS-element detection via rust-ise (Rust ISEScan, in-process).
//!
//! Calls `rust_ise::run` (byte-identical to the rust-ise binary) and emits
//! `Feature{IsElement}`. rust-ise internally shells out to `mmseqs2` (must be on
//! PATH) and needs an IS profile DB (`mmdb_union/…`, passed as `is_db`).

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use rust_ise::{run, IseConfig};
use std::collections::HashMap;
use std::io::{IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;

/// rust-ise's fp_control module consults `<is_db>/fpc/refset.dbtype`; strict
/// FP-control (host-recombinase false-positive filtering) is a silent **no-op**
/// unless that file exists. This is the single canonical presence check, reused
/// by `setup-db` to verify a freshly-extracted IS DB carries `fpc/`.
pub fn fpc_present(is_db: &Path) -> bool {
    is_db.join("fpc").join("refset.dbtype").is_file()
}

/// Pure, testable decision: should `detect_is` interactively prompt to download
/// the complete (fpc-carrying) IS DB?
///
/// Prompt only when the IS DB is unusable for strict mode — either its directory
/// is missing entirely, or `fpc/` is absent — AND we are safe to ask: strict is
/// on (the user did not set `BACTARS_ISSTRICT=0`) and stdin is a real terminal.
/// Non-interactive/pipeline runs (non-TTY) and strict-disabled runs never prompt,
/// so nothing can hang or block a scripted invocation.
fn should_prompt_download(db_missing: bool, fpc_missing: bool, strict: bool, is_tty: bool) -> bool {
    strict && is_tty && (db_missing || fpc_missing)
}

/// Guards the one-time-per-process fpc warning so batch runs over many genomes
/// do not repeat the same multi-line notice.
static FPC_WARNED: AtomicBool = AtomicBool::new(false);

/// Process-wide memo of the interactive IS-DB resolution: the outcome of the
/// (at-most-once) download prompt. `Some(dir)` = a freshly-fetched fpc-carrying DB
/// to search with; `None` = keep the configured DB. Memoized so a multi-genome
/// batch prompts at most once — later genomes reuse this outcome and never re-ask.
static RESOLVED_ISDB: OnceLock<Option<PathBuf>> = OnceLock::new();

/// Emit the "strict FP-control is a no-op" warning to stderr, once per process.
fn warn_fpc_missing(is_db: &str, db_missing: bool) {
    if FPC_WARNED.swap(true, Ordering::Relaxed) {
        return;
    }
    if db_missing {
        eprintln!(
            "[is][warn] IS DB directory does not exist: {is_db}\n\
             [is][warn] IS detection will have no profiles to search, and strict precision mode\n\
             [is][warn] (host-recombinase FP filtering) cannot run — it needs `fpc/refset`."
        );
    } else {
        eprintln!(
            "[is][warn] IS FP-control (strict precision mode) is a NO-OP: the IS DB at\n\
             [is][warn]   {is_db}\n\
             [is][warn] lacks an `fpc/refset` subdirectory, so rust-ise will NOT filter\n\
             [is][warn] host-recombinase false positives (e.g. fimB/fimE, other site-specific\n\
             [is][warn] recombinases) — IS-element precision will be lower than calibrated."
        );
    }
    eprintln!(
        "[is][warn] Get the complete fpc-carrying IS DB (rust-ise-isdb-fpc.tar.gz):\n\
         [is][warn]   export BACTARS_ISDB_URL=https://<host>/rust-ise-isdb-fpc.tar.gz\n\
         [is][warn]   bactars setup-db --only isdb --out DB      # extracts DB/isdb/{{mmdb_union,fpc}}\n\
         [is][warn]   bactars <genome> --is-db DB/isdb ...\n\
         [is][warn] Or set BACTARS_ISSTRICT=0 to disable strict mode and silence this warning."
    );
}

/// Cache directory for an interactively-downloaded IS DB. Prefers
/// `$XDG_CACHE_HOME/bactars/isdb`, then `$HOME/.cache/bactars/isdb`, falling back
/// to a per-user temp path so the download always has somewhere to land.
fn isdb_cache_dir() -> PathBuf {
    if let Ok(x) = std::env::var("XDG_CACHE_HOME") {
        if !x.is_empty() {
            return PathBuf::from(x).join("bactars").join("isdb");
        }
    }
    if let Ok(h) = std::env::var("HOME") {
        if !h.is_empty() {
            return PathBuf::from(h).join(".cache").join("bactars").join("isdb");
        }
    }
    std::env::temp_dir().join("bactars_isdb")
}

/// Prompt on the terminal and return true only for an explicit yes. Any read
/// error or EOF is treated as "no" so the run proceeds rather than blocks.
fn prompt_yes_no(question: &str) -> bool {
    eprint!("{question}");
    if std::io::stderr().flush().is_err() {
        return false;
    }
    let mut line = String::new();
    if std::io::stdin().read_line(&mut line).is_err() {
        return false;
    }
    let ans = line.trim().to_ascii_lowercase();
    ans == "y" || ans == "yes"
}

/// Attempt to download + extract the fpc-carrying IS DB into the cache dir,
/// returning the directory to hand to `--is-db`. Resolves the source URL from
/// `BACTARS_ISDB_URL` (or the compiled-in default); errors clearly if none is set.
fn try_download_isdb() -> Result<PathBuf, String> {
    let override_url = std::env::var("BACTARS_ISDB_URL").ok();
    let url = crate::setup_db::effective_isdb_url(override_url.as_deref()).ok_or_else(|| {
        "no IS-DB source URL configured (set BACTARS_ISDB_URL to the \
         rust-ise-isdb-fpc.tar.gz location)"
            .to_string()
    })?;
    let dir = isdb_cache_dir();
    eprintln!("[is] downloading complete IS DB from {url}");
    let is_dir = crate::setup_db::download_isdb(&url, &dir, false)?;
    if !fpc_present(&is_dir) {
        eprintln!(
            "[is][warn] downloaded IS DB at {} still lacks `fpc/refset.dbtype` \
             (not the fpc-carrying build); strict mode will remain a no-op.",
            is_dir.display()
        );
    }
    Ok(is_dir)
}

/// Detect IS elements in a genome FASTA, as `Feature{IsElement}`. `is_db` is the
/// isscan DB directory; `mmseqs` must be on PATH.
pub fn detect_is(genome_path: &str, is_db: &str, threads: usize) -> Result<Vec<Feature>, String> {
    // Intermediate files (proteome.faa, mmseqs tmp) go under a work dir we own.
    // Use a per-process dir and wipe any stale contents first — mmseqs refuses to
    // overwrite existing DBs, so a reused dir makes `mmseqs search` fail.
    let work_dir = std::env::temp_dir().join(format!("bactars_is_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&work_dir);
    // rust-ise passes `threads` straight to `mmseqs --threads`, which rejects 0;
    // resolve the "0 = all cores" convention to a concrete count here.
    let threads = if threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        threads
    };
    // Precision mode: rust-ise drops `fpFlag == "host"` calls — host-lineage FPs (e.g.
    // fimB/fimE/pinE and other site-specific recombinases) that share the DDE/tyrosine
    // fold with transposases but have NO IS-positive nt homology and no TIR/TSD. Enabled
    // by default: cross-genome validation over 5 diverse taxa (E. coli O157, Klebsiella,
    // Pseudomonas, Acinetobacter, Bacillus) showed strict drops ONLY genuine host FPs
    // (9 total, all recombinases/integrases per RefSeq) and lost ZERO real IS elements,
    // lifting IS precision (e.g. E. coli 0.85->0.90, Bacillus 0.76->0.84). Overridable
    // to off with BACTARS_ISSTRICT=0 (requires the IS DB to carry an `fpc/refset` subdir).
    let strict = std::env::var("BACTARS_ISSTRICT")
        .map(|v| !(v == "0" || v.eq_ignore_ascii_case("false")))
        .unwrap_or(true);

    // Tell the user when strict precision mode cannot actually do anything because
    // the IS DB lacks `fpc/refset` (fp_control silently no-ops without it), and —
    // on an interactive terminal only — offer to fetch the complete fpc-carrying
    // DB. Non-TTY / BACTARS_ISSTRICT=0 runs never prompt, so pipelines never hang.
    // The (possibly freshly-downloaded) DB path we ultimately search with:
    let mut is_db: PathBuf = PathBuf::from(is_db);
    if strict {
        let db_missing = !is_db.exists();
        let fpc_missing = !fpc_present(&is_db);
        if fpc_missing {
            // Resolve the download decision ONCE per process and memoize it, so an
            // interactive multi-genome batch prompts at most once — subsequent
            // genomes reuse the outcome (a downloaded dir, or `None` = keep config)
            // instead of re-prompting per genome.
            let is_db_str = is_db.to_string_lossy().into_owned();
            let is_tty = std::io::stdin().is_terminal();
            let resolved = RESOLVED_ISDB.get_or_init(move || {
                warn_fpc_missing(&is_db_str, db_missing);
                if should_prompt_download(db_missing, fpc_missing, strict, is_tty)
                    && prompt_yes_no("Download the complete IS DB now? [y/N] ")
                {
                    match try_download_isdb() {
                        Ok(new_dir) => {
                            eprintln!(
                                "[is] using freshly-fetched IS DB at {}",
                                new_dir.display()
                            );
                            Some(new_dir)
                        }
                        Err(e) => {
                            eprintln!(
                                "[is][warn] IS DB download failed ({e}); proceeding with \
                                 the existing configuration."
                            );
                            None
                        }
                    }
                } else {
                    None
                }
            });
            if let Some(new_dir) = resolved {
                is_db = new_dir.clone();
            }
        }
    }

    let cfg = IseConfig {
        threads,
        gpu: false,
        strict,
        db_dir: is_db,
    };
    let calls = run(Path::new(genome_path), &work_dir, &cfg)?;

    let mut per_contig: HashMap<String, usize> = HashMap::new();
    let mut feats = Vec::new();
    for c in calls {
        let n = per_contig.entry(c.seqid.clone()).or_insert(0);
        *n += 1;
        let (start, end) = if c.is_begin <= c.is_end {
            (c.is_begin, c.is_end)
        } else {
            (c.is_end, c.is_begin)
        };
        feats.push(Feature {
            kind: FeatureKind::IsElement,
            contig: c.seqid.clone(),
            id: format!("{}_is{}", c.seqid, n),
            start,
            end,
            strand: if c.strand == '-' { -1 } else { 1 },
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: format!("rust-ise:{}", c.fp_flag),
                accession: c.family.clone(),
                name: c.family.clone(),
                score: c.pident as f32,
                // rust-ise carries a real ISfinder homology-search E-value on each
                // call; keep it (guarding the rare unset/NaN case).
                evalue: (!c.evalue.is_nan()).then_some(c.evalue),
                ref_len: None,
            }],
            func: Functional::default(),
        });
    }
    Ok(feats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Make a unique scratch dir under the system temp dir for a test.
    fn scratch(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!(
            "bactars_is_test_{}_{}_{:?}",
            tag,
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or(0)
        ));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn fpc_present_detects_refset_dbtype() {
        // A DB dir WITHOUT fpc/ -> absent.
        let no_fpc = scratch("nofpc");
        std::fs::create_dir_all(no_fpc.join("mmdb_union")).unwrap();
        assert!(!fpc_present(&no_fpc), "no fpc/ dir should read as absent");

        // fpc/ present but no refset.dbtype file -> still absent.
        let empty_fpc = scratch("emptyfpc");
        std::fs::create_dir_all(empty_fpc.join("fpc")).unwrap();
        assert!(
            !fpc_present(&empty_fpc),
            "fpc/ without refset.dbtype should read as absent"
        );

        // A DB dir WITH fpc/refset.dbtype -> present.
        let with_fpc = scratch("withfpc");
        std::fs::create_dir_all(with_fpc.join("fpc")).unwrap();
        std::fs::write(with_fpc.join("fpc").join("refset.dbtype"), b"x").unwrap();
        assert!(fpc_present(&with_fpc), "fpc/refset.dbtype present should read as present");

        // A missing directory entirely -> absent (no panic).
        let missing = scratch("missing");
        std::fs::remove_dir_all(&missing).unwrap();
        assert!(!fpc_present(&missing));

        // Cleanup.
        for d in [&no_fpc, &empty_fpc, &with_fpc] {
            let _ = std::fs::remove_dir_all(d);
        }
    }

    #[test]
    fn should_prompt_download_tty_and_strict_gating() {
        // The happy path: strict on, on a TTY, and the DB is unusable (fpc missing
        // or dir missing) -> prompt.
        assert!(should_prompt_download(false, true, true, true), "fpc missing on TTY should prompt");
        assert!(should_prompt_download(true, true, true, true), "db missing on TTY should prompt");

        // Never prompt when stdin is not a terminal (pipeline / CI) — no hang.
        assert!(!should_prompt_download(true, true, true, false), "non-TTY must never prompt");
        assert!(!should_prompt_download(false, true, true, false), "non-TTY must never prompt");

        // Never prompt when strict is disabled (BACTARS_ISSTRICT=0), even on a TTY.
        assert!(!should_prompt_download(true, true, false, true), "strict off must never prompt");
        assert!(!should_prompt_download(false, true, false, true), "strict off must never prompt");

        // Nothing to fix (fpc present, dir present) -> no prompt even on a TTY.
        assert!(!should_prompt_download(false, false, true, true), "healthy DB must not prompt");
    }
}
