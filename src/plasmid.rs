//! PlasmidFinder replicon typing: identify which contigs carry a plasmid
//! replication-initiation (rep/Inc) locus, and report the replicon (Inc) type.
//!
//! This is CONTIG-LEVEL metadata (which contigs look like plasmids + their Inc
//! type), NOT a per-gene [`Feature`]. We run an mmseqs nucleotide search
//! (`--search-type 3`) of every contig against the concatenated PlasmidFinder
//! replicon-allele database and keep hits at PlasmidFinder's standard thresholds
//! (>=80% identity, >=60% coverage of the replicon allele). Returns one
//! [`PlasmidHit`] per (contig, replicon-type) best match.
//!
//! Pipeline (mirrors psc.rs / vfdb.rs): untar the DB and concatenate its `.fsa`
//! replicon files ONCE, caching the extracted dir + a concatenated FASTA next to
//! the tarball so repeat runs skip the rebuild -> write the contigs to a temp
//! FASTA -> `mmseqs createdb` (query + target, nucleotide) -> `search
//! --search-type 3` -> `convertalis` -> best replicon per (contig, type). Any
//! mmseqs / IO failure logs to stderr and returns an empty list (never crashes).

use crate::fasta::Contig;
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Minimum percent identity for a replicon call (PlasmidFinder default).
const MIN_PIDENT: f64 = 80.0;
/// Minimum fraction (0..1) of the replicon allele covered (PlasmidFinder default).
const MIN_TCOV: f64 = 0.60;

/// A detected plasmid replicon on one contig.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasmidHit {
    /// Contig the replicon locus sits on.
    pub contig: String,
    /// Replicon (Inc) type, e.g. `IncFII`, `IncHI1B(R27)`, `Rep3`.
    pub replicon_type: String,
    /// Percent identity of the best alignment (0..100).
    pub identity: f64,
    /// 1-based start of the alignment on the contig.
    pub start: i64,
    /// 1-based end of the alignment on the contig.
    pub end: i64,
}

/// Extract the replicon (Inc) type from a PlasmidFinder allele name. The naming
/// convention is `<repliconType>_<alleleNo>_[<name>]_<accession>`, so the type is
/// everything up to the first `_` (which itself may contain `/`, `(...)`, `-`).
/// Examples: `IncFII_1__CP...` -> `IncFII`; `IncHI1B(R27)_1_R27_AF250878` ->
/// `IncHI1B(R27)`; `IncB/O/K/Z_1__CU928147` -> `IncB/O/K/Z`.
pub fn replicon_type(allele: &str) -> String {
    allele.split('_').next().unwrap_or(allele).to_string()
}

/// One parsed convertalis row.
#[derive(Clone, Debug)]
struct M8Hit {
    query: String,
    target: String,
    pident: f64,
    tcov: f64,
    qstart: i64,
    qend: i64,
    bits: f32,
}

/// Parse `query,target,pident,qcov,tcov,qstart,qend,bits` into the best hit per
/// (contig, replicon-type) that clears the PlasmidFinder thresholds. Best =
/// highest bit score.
fn parse_m8_best(path: &Path) -> HashMap<(String, String), PlasmidHit> {
    let mut best: HashMap<(String, String), (f32, PlasmidHit)> = HashMap::new();
    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return HashMap::new(),
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 8 {
            continue;
        }
        let pident_raw: f64 = p[2].parse().unwrap_or(0.0);
        // mmseqs may report pident as a fraction (0..1); normalise to a percent so
        // the MIN_PIDENT (%) compare is consistent across all stages (F8).
        let pident = if pident_raw <= 1.0 { pident_raw * 100.0 } else { pident_raw };
        let hit = M8Hit {
            query: p[0].to_string(),
            target: p[1].to_string(),
            pident,
            tcov: p[4].parse().unwrap_or(0.0),
            qstart: p[5].parse().unwrap_or(0),
            qend: p[6].parse().unwrap_or(0),
            bits: p[7].parse().unwrap_or(0.0),
        };
        if hit.pident < MIN_PIDENT || hit.tcov < MIN_TCOV {
            continue;
        }
        let rtype = replicon_type(&hit.target);
        let (start, end) = if hit.qstart <= hit.qend {
            (hit.qstart, hit.qend)
        } else {
            (hit.qend, hit.qstart)
        };
        let key = (hit.query.clone(), rtype.clone());
        let ph = PlasmidHit {
            contig: hit.query.clone(),
            replicon_type: rtype,
            identity: hit.pident,
            start,
            end,
        };
        match best.get(&key) {
            Some((b, _)) if *b >= hit.bits => {}
            _ => {
                best.insert(key, (hit.bits, ph));
            }
        }
    }
    best.into_iter().map(|(k, (_, ph))| (k, ph)).collect()
}

/// Ensure the tarball is extracted and its `.fsa` replicon files concatenated to
/// `<dir>/plasmid_refs.fna` (cached). Returns the concatenated FASTA path.
fn ensure_refs(plasmid_dir: &Path) -> Result<std::path::PathBuf, String> {
    let refs = plasmid_dir.join("plasmid_refs.fna");
    if refs.exists() {
        return Ok(refs);
    }
    let tar_gz = plasmid_dir.join("plasmidfinder_db.tar.gz");
    if !tar_gz.exists() {
        return Err(format!("{} not found", tar_gz.display()));
    }
    let extract_dir = plasmid_dir.join("extracted");
    let _ = fs::remove_dir_all(&extract_dir);
    fs::create_dir_all(&extract_dir).map_err(|e| format!("mkdir extracted: {e}"))?;
    // Extract in-process (pure Rust), guarding against path traversal.
    crate::util_io::extract_tar_gz(&tar_gz, &extract_dir)?;
    // Concatenate every .fsa found under the extracted tree.
    let mut concat = String::new();
    collect_fsa(&extract_dir, &mut concat);
    if concat.is_empty() {
        return Err("no .fsa replicon files found in the PlasmidFinder tarball".to_string());
    }
    fs::write(&refs, concat).map_err(|e| format!("write plasmid_refs.fna: {e}"))?;
    Ok(refs)
}

/// Recursively append the contents of every `*.fsa` file under `dir` to `out`.
fn collect_fsa(dir: &Path, out: &mut String) {
    let Ok(entries) = fs::read_dir(dir) else { return };
    for e in entries.flatten() {
        let p = e.path();
        if p.is_dir() {
            collect_fsa(&p, out);
        } else if p.extension().and_then(|x| x.to_str()) == Some("fsa") {
            if let Ok(text) = fs::read_to_string(&p) {
                out.push_str(&text);
                if !text.ends_with('\n') {
                    out.push('\n');
                }
            }
        }
    }
}

/// Detect plasmid replicons across the contigs. `plasmid_dir` holds
/// `plasmidfinder_db.tar.gz`. Degrades gracefully: on any mmseqs / IO error it
/// logs to stderr and returns an empty list so the pipeline continues.
pub fn detect(contigs: &[Contig], plasmid_dir: &str, threads: usize) -> Result<Vec<PlasmidHit>, String> {
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    let plasmid_dir = Path::new(plasmid_dir);
    if !crate::mmseqs::available() {
        eprintln!("[plasmid] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return Ok(Vec::new());
    }
    let bin = crate::mmseqs::bin();

    let refs = match ensure_refs(plasmid_dir) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("[plasmid] cannot prepare PlasmidFinder DB under {}: {e} — skipping", plasmid_dir.display());
            return Ok(Vec::new());
        }
    };

    // Build the contig query FASTA.
    let mut query = String::new();
    for c in contigs {
        query.push('>');
        query.push_str(&c.name);
        query.push('\n');
        query.push_str(&String::from_utf8_lossy(&c.seq));
        query.push('\n');
    }

    let work = std::env::temp_dir().join(format!("bactars_plasmid_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[plasmid] cannot create temp dir {}: {e}", work.display());
        return Ok(Vec::new());
    }
    let result = run_search(&bin, &work, &refs, &query, threads);
    let _ = fs::remove_dir_all(&work);

    let best = match result {
        Ok(b) => b,
        Err(e) => {
            eprintln!("[plasmid] mmseqs pipeline failed: {e} — continuing without replicon typing");
            return Ok(Vec::new());
        }
    };

    let mut hits: Vec<PlasmidHit> = best.into_values().collect();
    // Deterministic order: contig, then replicon type.
    hits.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.replicon_type.cmp(&b.replicon_type))
    });
    if !hits.is_empty() {
        eprintln!("[plasmid] detected {} plasmid replicon(s)", hits.len());
    }
    Ok(hits)
}

/// createdb(query, refs nucleotide) -> search --search-type 3 -> convertalis ->
/// best replicon per (contig, type).
fn run_search(
    bin: &str,
    work: &Path,
    refs: &Path,
    query: &str,
    threads: usize,
) -> Result<HashMap<(String, String), PlasmidHit>, String> {
    let query_fna = work.join("query.fna");
    let query_db = work.join("queryDB");
    let target_db = work.join("targetDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");
    let nt = threads.clamp(1, 8).to_string();

    fs::write(&query_fna, query).map_err(|e| format!("write query.fna: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_fna, q_db, t_db, res_db, tmp_s, refs_s, m8_s) = (
        s(&query_fna),
        s(&query_db),
        s(&target_db),
        s(&result_db),
        s(&tmp),
        s(refs),
        s(&m8),
    );

    crate::mmseqs::run_with(bin, &["createdb", &q_fna, &q_db, "--dbtype", "2", "-v", "1"])?;
    crate::mmseqs::run_with(bin, &["createdb", &refs_s, &t_db, "--dbtype", "2", "-v", "1"])?;
    crate::mmseqs::run_with(
        bin,
        &[
            "search", &q_db, &t_db, &res_db, &tmp_s, "--threads", &nt, "--search-type", "3",
            "-e", "1e-10", "--max-seqs", "300", "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        bin,
        &[
            "convertalis", &q_db, &t_db, &res_db, &m8_s, "--format-output",
            "query,target,pident,qcov,tcov,qstart,qend,bits", "--threads", &nt, "-v", "1",
        ],
    )?;
    Ok(parse_m8_best(&m8))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn replicon_type_strips_allele_suffix() {
        assert_eq!(replicon_type("IncFII_1__CP011611"), "IncFII");
        assert_eq!(replicon_type("IncHI1B(R27)_1_R27_AF250878"), "IncHI1B(R27)");
        assert_eq!(replicon_type("IncB/O/K/Z_1__CU928147"), "IncB/O/K/Z");
        assert_eq!(replicon_type("pKPC-CAV1321_1__CP011611"), "pKPC-CAV1321");
        assert_eq!(replicon_type("Rep3_1_repA"), "Rep3");
    }

    #[test]
    fn m8_best_hit_applies_thresholds_and_collapses_by_type() {
        let dir = std::env::temp_dir().join(format!("bactars_plasmid_m8_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        let m8 = dir.join("result.m8");
        let mut fh = fs::File::create(&m8).unwrap();
        // query,target,pident,qcov,tcov,qstart,qend,bits
        // Two alleles of the same IncFII replicon on contig1; higher-bit one wins.
        writeln!(fh, "contig1\tIncFII_1__A\t95.0\t0.10\t0.90\t500\t1200\t400").unwrap();
        writeln!(fh, "contig1\tIncFII_2__B\t92.0\t0.10\t0.85\t600\t1300\t350").unwrap();
        // Below identity threshold -> rejected.
        writeln!(fh, "contig1\tIncX_1__C\t70.0\t0.10\t0.90\t10\t700\t300").unwrap();
        // Below coverage threshold -> rejected.
        writeln!(fh, "contig2\tRep3_1_r\t99.0\t0.10\t0.30\t5\t400\t250").unwrap();
        // A distinct replicon on contig2 that passes.
        writeln!(fh, "contig2\tColRNAI_1__D\t88.0\t0.10\t0.75\t900\t100\t200").unwrap();
        drop(fh);

        let best = parse_m8_best(&m8);
        let incfii = best.get(&("contig1".to_string(), "IncFII".to_string())).unwrap();
        assert_eq!(incfii.identity, 95.0);
        assert_eq!(incfii.start, 500);
        assert_eq!(incfii.end, 1200);
        assert!(!best.contains_key(&("contig1".to_string(), "IncX".to_string())));
        assert!(!best.contains_key(&("contig2".to_string(), "Rep3".to_string())));
        // Reverse-orientation coords normalised (qstart 900 > qend 100).
        let col = best.get(&("contig2".to_string(), "ColRNAI".to_string())).unwrap();
        assert_eq!(col.start, 100);
        assert_eq!(col.end, 900);
        let _ = fs::remove_dir_all(&dir);
    }
}
