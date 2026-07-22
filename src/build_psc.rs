//! Build the PSC (Protein Similarity Cluster) reference DB from UniRef90 — a
//! pure-Rust port of `db/psc/build_psc.sh`. The ONLY external binary is mmseqs
//! (already bactars' one allowed non-Rust dependency); every other step —
//! HTTP download, gzip/tar, the NCBI-taxonomy descent, and the 32 GB streaming
//! taxonomy filter — is done in-process here (no curl/pigz/gawk/python).
//!
//! Produces, under `psc_dir`:
//!   * `taxdump/{nodes.dmp,merged.dmp}` — NCBI taxonomy (downloaded)
//!   * `accepted_taxids.txt`            — TaxIDs under Bacteria/Archaea/Viruses
//!   * `uniref90_prok.fasta`            — taxonomy-filtered UniRef90 reps
//!   * `names.tsv`                      — `id<TAB>product` for CDS naming
//!   * `psc_db`, `psc_db.*`             — mmseqs search DB + index
//!
//! This is the FULL-tier naming DB `bactars --full` uses (PSC naming + the
//! `--full` reference-length pseudogene signal). It is heavy (32 GB download,
//! ~194 GB index, tens of minutes to hours) — run it once, non-interactively.

use crate::util_io;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

/// NCBI taxonomy dump (nodes.dmp + merged.dmp), ~75 MB.
const TAXDUMP_URL: &str = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";
/// UniRef90 cluster representatives (all domains), ~32 GB gzip.
const UNIREF90_URL: &str =
    "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz";
/// Superkingdom roots kept: Bacteria, Archaea, Viruses.
const ROOTS: [u64; 3] = [2, 2157, 10239];

/// Inputs / outputs for a PSC build.
pub struct PscBuildConfig {
    /// Output directory (e.g. `<bundle>/psc`).
    pub psc_dir: PathBuf,
    /// UniRef90 gzip FASTA (`<bundle>/uniref/uniref90.fasta.gz`); downloaded if absent.
    pub uniref_gz: PathBuf,
    /// mmseqs `createindex` threads.
    pub threads: usize,
    /// Rebuild even if `psc_db` already exists.
    pub force: bool,
}

/// Build the PSC DB end-to-end. Idempotent: skips when `psc_db` exists unless
/// `force`. Individual intermediate steps are also skipped when their output is
/// already present, so an interrupted build resumes cheaply.
pub fn build_psc(cfg: &PscBuildConfig) -> Result<(), String> {
    if !crate::mmseqs::available() {
        return Err("mmseqs not found on PATH (needed to build the PSC DB)".into());
    }
    let psc = &cfg.psc_dir;
    std::fs::create_dir_all(psc).map_err(|e| format!("create {}: {e}", psc.display()))?;

    let psc_db = psc.join("psc_db");
    if psc_db.with_extension("dbtype").exists() && !cfg.force {
        eprintln!("[psc] psc_db already present at {} — skipping build (use --force)", psc_db.display());
        return Ok(());
    }

    // --- Step 0: NCBI taxdump ---
    let taxdir = psc.join("taxdump");
    let nodes = taxdir.join("nodes.dmp");
    let merged = taxdir.join("merged.dmp");
    if !nodes.exists() || cfg.force {
        std::fs::create_dir_all(&taxdir).map_err(|e| format!("create {}: {e}", taxdir.display()))?;
        let tgz = taxdir.join("taxdump.tar.gz");
        eprintln!("[psc][0] downloading NCBI taxdump (~75 MB)");
        util_io::http_download(TAXDUMP_URL, &tgz)?;
        util_io::extract_tar_gz(&tgz, &taxdir)?;
        let _ = std::fs::remove_file(&tgz);
        if !nodes.exists() {
            return Err(format!("taxdump extracted but {} missing", nodes.display()));
        }
    }

    // --- Step 1: accepted TaxID set (tree descent + merged) ---
    eprintln!("[psc][1] computing accepted TaxIDs (Bacteria/Archaea/Viruses descendants)");
    let accepted = accepted_taxids(&nodes, &merged)?;
    eprintln!("[psc][1] accepted TaxIDs: {}", accepted.len());
    // Persist for parity with the shell build / debugging.
    write_taxids(&psc.join("accepted_taxids.txt"), &accepted)?;

    // --- Step 2: UniRef90 (download if absent) + streaming taxonomy filter ---
    if !cfg.uniref_gz.exists() {
        if let Some(parent) = cfg.uniref_gz.parent() {
            std::fs::create_dir_all(parent).map_err(|e| format!("create {}: {e}", parent.display()))?;
        }
        eprintln!("[psc][2] downloading UniRef90 (~32 GB) -> {}", cfg.uniref_gz.display());
        util_io::http_download(UNIREF90_URL, &cfg.uniref_gz)?;
    }
    let prok_fa = psc.join("uniref90_prok.fasta");
    let names = psc.join("names.tsv");
    eprintln!("[psc][2] filtering UniRef90 by TaxID -> {} (+ names.tsv)", prok_fa.display());
    let (seen, kept) = filter_uniref(&cfg.uniref_gz, &accepted, &prok_fa, &names)?;
    eprintln!("[psc][2] records seen {seen}, kept {kept}");
    if kept == 0 {
        return Err("taxonomy filter kept 0 records — check TaxID parsing / accepted set".into());
    }

    // --- Step 3: mmseqs createdb + createindex ---
    let bin = crate::mmseqs::bin();
    let db = psc_db.to_string_lossy().into_owned();
    let fa = prok_fa.to_string_lossy().into_owned();
    let tmp = psc.join("tmp");
    std::fs::create_dir_all(&tmp).map_err(|e| format!("create {}: {e}", tmp.display()))?;
    let tmp_s = tmp.to_string_lossy().into_owned();
    let nthreads = if cfg.threads == 0 { "16".to_string() } else { cfg.threads.to_string() };
    eprintln!("[psc][3] mmseqs createdb");
    crate::mmseqs::run_with(&bin, &["createdb", &fa, &db, "--compressed", "1", "-v", "1"])?;
    eprintln!("[psc][3] mmseqs createindex (threads={nthreads}) — heavy, builds ~194 GB index");
    crate::mmseqs::run_with(
        &bin,
        &["createindex", &db, &tmp_s, "--threads", &nthreads, "--compressed", "1", "-v", "1"],
    )?;
    let _ = std::fs::remove_dir_all(&tmp); // scratch

    eprintln!("[psc][done] PSC DB built at {}", psc_db.display());
    Ok(())
}

/// Parse `nodes.dmp` (`taxid \t| \t parent \t| ...`) into a parent→children map,
/// descend from [`ROOTS`], then fold in `merged.dmp` obsolete ids resolving into
/// the accepted set. Returns the accepted TaxID set.
fn accepted_taxids(nodes: &Path, merged: &Path) -> Result<HashSet<u64>, String> {
    let f = File::open(nodes).map_err(|e| format!("open {}: {e}", nodes.display()))?;
    let mut children: HashMap<u64, Vec<u64>> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line.map_err(|e| format!("read nodes.dmp: {e}"))?;
        // fields separated by "\t|\t"; [0]=taxid [1]=parent
        let mut it = line.split("\t|\t");
        let taxid: u64 = match it.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let parent: u64 = match it.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if taxid != parent {
            children.entry(parent).or_default().push(taxid);
        }
    }
    let mut accepted: HashSet<u64> = HashSet::new();
    for &root in &ROOTS {
        let mut stack = vec![root];
        while let Some(t) = stack.pop() {
            if !accepted.insert(t) {
                continue;
            }
            if let Some(kids) = children.get(&t) {
                stack.extend(kids.iter().copied());
            }
        }
    }
    // merged.dmp: old_taxid \t| \t new_taxid \t|
    if let Ok(f) = File::open(merged) {
        for line in BufReader::new(f).lines() {
            let line = match line {
                Ok(l) => l,
                Err(_) => break,
            };
            let mut it = line.split("\t|\t");
            let old: Option<u64> = it.next().and_then(|s| s.trim().parse().ok());
            let new: Option<u64> = it
                .next()
                .and_then(|s| s.split('\t').next())
                .and_then(|s| s.trim().parse().ok());
            if let (Some(old), Some(new)) = (old, new) {
                if accepted.contains(&new) {
                    accepted.insert(old);
                }
            }
        }
    }
    Ok(accepted)
}

fn write_taxids(path: &Path, set: &HashSet<u64>) -> Result<(), String> {
    let mut v: Vec<u64> = set.iter().copied().collect();
    v.sort_unstable();
    let mut w = BufWriter::new(File::create(path).map_err(|e| format!("create {}: {e}", path.display()))?);
    for t in v {
        writeln!(w, "{t}").map_err(|e| format!("write taxids: {e}"))?;
    }
    w.flush().map_err(|e| e.to_string())
}

/// Stream the gzip UniRef90 FASTA, keep only records whose `TaxID=` is accepted,
/// and emit the filtered FASTA + an `id<TAB>product` names table in ONE pass —
/// never loading the 32 GB file into memory. Returns `(records_seen, records_kept)`.
fn filter_uniref(
    uniref_gz: &Path,
    accepted: &HashSet<u64>,
    out_fasta: &Path,
    out_names: &Path,
) -> Result<(u64, u64), String> {
    let gz = File::open(uniref_gz).map_err(|e| format!("open {}: {e}", uniref_gz.display()))?;
    let mut rdr = BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(gz));
    let mut fa = BufWriter::with_capacity(
        1 << 20,
        File::create(out_fasta).map_err(|e| format!("create {}: {e}", out_fasta.display()))?,
    );
    let mut nm = BufWriter::with_capacity(
        1 << 20,
        File::create(out_names).map_err(|e| format!("create {}: {e}", out_names.display()))?,
    );
    let mut buf: Vec<u8> = Vec::with_capacity(256);
    let mut keep = false;
    let (mut seen, mut kept) = (0u64, 0u64);
    loop {
        buf.clear();
        let n = rdr.read_until(b'\n', &mut buf).map_err(|e| format!("read uniref: {e}"))?;
        if n == 0 {
            break;
        }
        if buf.first() == Some(&b'>') {
            seen += 1;
            keep = header_accepted(&buf, accepted);
            if keep {
                kept += 1;
                write_name_row(&buf, &mut nm)?;
                fa.write_all(&buf).map_err(|e| format!("write fasta: {e}"))?;
            }
        } else if keep {
            fa.write_all(&buf).map_err(|e| format!("write fasta: {e}"))?;
        }
    }
    fa.flush().map_err(|e| e.to_string())?;
    nm.flush().map_err(|e| e.to_string())?;
    Ok((seen, kept))
}

/// Is a `>` header line's `TaxID=NNN` in the accepted set?
fn header_accepted(header: &[u8], accepted: &HashSet<u64>) -> bool {
    match parse_taxid(header) {
        Some(tid) => accepted.contains(&tid),
        None => false,
    }
}

/// Extract the integer after the first `TaxID=` in a header line.
fn parse_taxid(header: &[u8]) -> Option<u64> {
    const TAG: &[u8] = b"TaxID=";
    let mut i = 0;
    while i + TAG.len() <= header.len() {
        if &header[i..i + TAG.len()] == TAG {
            let mut j = i + TAG.len();
            let start = j;
            while j < header.len() && header[j].is_ascii_digit() {
                j += 1;
            }
            if j > start {
                return std::str::from_utf8(&header[start..j]).ok()?.parse().ok();
            }
            return None;
        }
        i += 1;
    }
    None
}

/// Write `id<TAB>product` for a kept header. `id` = first token without `>`;
/// `product` = text between the first space and ` n=` (UniRef metadata), or
/// `"uncharacterized protein"` when empty. Mirrors filter_uniref.awk.
fn write_name_row<W: Write>(header: &[u8], out: &mut W) -> Result<(), String> {
    // strip leading '>' and trailing newline
    let line = {
        let s = if header.first() == Some(&b'>') { &header[1..] } else { header };
        let end = s.iter().rposition(|&b| b != b'\n' && b != b'\r').map(|p| p + 1).unwrap_or(0);
        &s[..end]
    };
    let sp = line.iter().position(|&b| b == b' ');
    let id = match sp {
        Some(p) => &line[..p],
        None => line,
    };
    let mut prod: &[u8] = match sp {
        Some(p) => &line[p + 1..],
        None => b"",
    };
    // cut at the UniRef metadata trailer " n=<count> ..." (awk: / n=[0-9]+ Tax=/).
    // Require a digit after " n=" so a product name containing " n=word" is not
    // truncated; UniRef's cluster-size field is always " n=<digits>".
    if let Some(p) = find_n_meta(prod) {
        prod = &prod[..p];
    }
    out.write_all(id).map_err(|e| e.to_string())?;
    out.write_all(b"\t").map_err(|e| e.to_string())?;
    if prod.is_empty() {
        out.write_all(b"uncharacterized protein").map_err(|e| e.to_string())?;
    } else {
        out.write_all(prod).map_err(|e| e.to_string())?;
    }
    out.write_all(b"\n").map_err(|e| e.to_string())?;
    Ok(())
}

/// Index of the UniRef metadata trailer `" n=<digit>"` (the cluster-size field),
/// or `None`. Requires a digit right after `" n="` so a product name that itself
/// contains `" n=word"` is not truncated.
fn find_n_meta(prod: &[u8]) -> Option<usize> {
    let mut i = 0;
    while i + 3 < prod.len() {
        if &prod[i..i + 3] == b" n=" && prod[i + 3].is_ascii_digit() {
            return Some(i);
        }
        i += 1;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn taxid_and_name_parse() {
        let h = b">UniRef90_P12345 DNA gyrase subunit A n=42 Tax=Escherichia coli TaxID=562 RepID=GYRA_ECOLI\n";
        assert_eq!(parse_taxid(h), Some(562));
        let mut out = Vec::new();
        write_name_row(h, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "UniRef90_P12345\tDNA gyrase subunit A\n"
        );
    }

    #[test]
    fn missing_taxid_and_empty_product() {
        assert_eq!(parse_taxid(b">X something with no taxid\n"), None);
        // Product name containing " n=word" must NOT be truncated (only " n=<digit>").
        let mut out = Vec::new();
        write_name_row(b">UniRef90_A2 protein n=terminal region n=3 Tax=B TaxID=1\n", &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "UniRef90_A2\tprotein n=terminal region\n");
        // Empty product (metadata immediately after id) -> uncharacterized.
        let mut out2 = Vec::new();
        write_name_row(b">UniRef90_Q1  n=1 Tax=x TaxID=2\n", &mut out2).unwrap();
        assert_eq!(String::from_utf8(out2).unwrap(), "UniRef90_Q1\tuncharacterized protein\n");
    }

    /// Real-data parity: on a machine that already has the shell-built PSC dir,
    /// the Rust taxid descent must reproduce the python-built `accepted_taxids.txt`
    /// EXACTLY (same set). Auto-skips where the taxdump is absent (CI / fresh checkout).
    #[test]
    fn taxid_parity_vs_python_build() {
        let psc = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../db/psc");
        let nodes = psc.join("taxdump/nodes.dmp");
        let merged = psc.join("taxdump/merged.dmp");
        let expected = psc.join("accepted_taxids.txt");
        if !nodes.exists() || !expected.exists() {
            eprintln!("taxid_parity: taxdump/accepted_taxids.txt absent — skipping");
            return;
        }
        let got = accepted_taxids(&nodes, &merged).unwrap();
        let want: std::collections::HashSet<u64> = std::fs::read_to_string(&expected)
            .unwrap()
            .lines()
            .filter_map(|l| l.trim().parse().ok())
            .collect();
        assert_eq!(got.len(), want.len(), "rust {} vs python {}", got.len(), want.len());
        assert_eq!(got, want, "accepted taxid SETS differ from the python build");
    }

    #[test]
    fn descent_and_merged() {
        let dir = std::env::temp_dir().join(format!("bactars_psctax_{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let nodes = dir.join("nodes.dmp");
        let merged = dir.join("merged.dmp");
        // tree: 1<-2(Bacteria)<-10<-20 ; 1<-9999(other)
        std::fs::write(
            &nodes,
            "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n10\t|\t2\t|\tphylum\t|\n20\t|\t10\t|\tgenus\t|\n9999\t|\t1\t|\tsuperkingdom\t|\n",
        )
        .unwrap();
        std::fs::write(&merged, "30\t|\t20\t|\n").unwrap(); // 30 -> 20 (accepted)
        let acc = accepted_taxids(&nodes, &merged).unwrap();
        assert!(acc.contains(&2) && acc.contains(&10) && acc.contains(&20));
        assert!(acc.contains(&30)); // merged into accepted
        assert!(!acc.contains(&9999)); // sibling superkingdom excluded
        let _ = std::fs::remove_dir_all(&dir);
    }
}
