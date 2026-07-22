//! VFDB virulence-factor annotation: flag CDS that are near-full-length homologs
//! of a curated VFDB core (setA) virulence-factor protein.
//!
//! For every CDS protein we run an mmseqs search against the VFDB **core** protein
//! set (`VFDB_setA_pro.fas.gz`, experimentally-verified VFs only). setA is the
//! higher-precision default: on non-pathogenic E. coli K-12 MG1655 it roughly
//! halves virulence calls (~94 → ~45) versus the full setB, dropping the
//! broadly-conserved flagellar/chemotaxis/housekeeping-secretion noise while
//! losing nothing setA itself would keep. A qualifying best hit (strong E-value AND high
//! query coverage — a virulence factor call should span most of the CDS, not a
//! short motif) ADDS a `/note` ("virulence factor: <description>") plus a
//! structured `/inference` (`VFDB:<VFG id>`) provenance record onto the CDS. It
//! never overwrites the resolved product — this is additive evidence only.
//!
//! Pipeline (mirrors psc.rs / amr_variant.rs): gunzip the DB to a plain FASTA and
//! `mmseqs createdb` it ONCE, caching both next to the `.gz` so repeat runs skip
//! the rebuild -> write the candidate CDS proteins to a temp FASTA ->
//! `mmseqs createdb` -> `search` -> `convertalis` -> best hit per query -> apply
//! the note/inference. Any mmseqs / IO failure logs to stderr and returns Ok(())
//! so the surrounding pipeline never crashes.

use crate::feature::{Feature, FeatureKind, Inference};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

/// mmseqs search sensitivity (VFDB homology naming is fine at high-ish -s).
const SEARCH_SENSITIVITY: &str = "5.7";
/// Max target sequences returned per query.
const MAX_SEQS: &str = "25";
/// Reporting E-value cutoff passed to mmseqs.
const EVALUE_CUTOFF: &str = "1e-10";
/// Accept a virulence-factor call only when the E-value is at least this good.
const MAX_EVALUE: f64 = 1e-10;
/// ...AND the alignment covers at least this fraction of the query CDS. VFDB
/// matches should be near-full-length, not a short shared domain.
const MIN_QCOV: f64 = 0.80;
/// ...AND the alignment covers at least this fraction of the TARGET VFDB protein.
/// Without a target-coverage floor a CDS that merely shares a domain with a large
/// VF protein is mislabelled; requiring most of the reference to align keeps only
/// genuine full-length homologs (the same false-positive control the sORF PSC tier
/// uses).
const MIN_TCOV: f64 = 0.70;
/// ...AND at least this percent identity. A virulence-factor assignment should be a
/// close homolog, not a distant match to a broadly-conserved family. NOTE: even the
/// curated setA core still contains genuine close homologs present in any
/// enterobacterium (iron uptake, type-1 fimbriae, curli), so a con-specific genome
/// legitimately carries several hits — VFDB reports virulence-ASSOCIATED homology
/// (presence != pathogenicity). setA is the higher-precision default; setB (full)
/// adds the flagellar/chemotaxis/predicted long tail if maximal recall is wanted.
const MIN_PIDENT: f64 = 90.0;

/// The VFDB core (setA) file layout under `vfdb_dir`: the gzipped source, the
/// gunzipped plain FASTA, and the mmseqs index. The index name (`vfdb_setA_db`)
/// is deliberately DISTINCT from the legacy setB index (`vfdb_db`) so a stale
/// setB `.dbtype` in the same directory is never silently reused.
fn vfdb_paths(vfdb_dir: &Path) -> (std::path::PathBuf, std::path::PathBuf, std::path::PathBuf) {
    (
        vfdb_dir.join("VFDB_setA_pro.fas.gz"),
        vfdb_dir.join("VFDB_setA_pro.fas"),
        vfdb_dir.join("vfdb_setA_db"),
    )
}

/// Emit a single warning the first time an m8 e-value column fails to parse.
static WARNED_EVALUE: AtomicBool = AtomicBool::new(false);
fn parse_evalue(raw: &str) -> f64 {
    match raw.parse() {
        Ok(v) => v,
        Err(_) => {
            if !WARNED_EVALUE.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "[vfdb] warning: malformed e-value '{raw}' in mmseqs m8 output — treating as +inf (further such rows silenced)"
                );
            }
            f64::INFINITY
        }
    }
}

/// One parsed VFDB header: the VFG id and a human-readable description, plus the
/// mmseqs target id (the first whitespace token of the header — what convertalis
/// reports in the `target` column) used to key hits back to this entry.
#[derive(Clone, Debug, PartialEq)]
pub struct VfdbEntry {
    /// mmseqs target id = first whitespace token, e.g. `VFG037170(gb|WP_001081754)`.
    pub target_id: String,
    /// The bare VFG accession, e.g. `VFG037170`.
    pub vfg_id: String,
    /// Gene symbol in the second parenthesis group, e.g. `plc1` (if present).
    pub gene: Option<String>,
    /// Free-text description, e.g. `phospholipase C`.
    pub description: String,
}

/// Parse a VFDB setA/setB protein header (with or without the leading `>`) into a
/// [`VfdbEntry`]. Header grammar:
///   `VFGxxxxxx(gb|ACC) (gene) description [Category (VF####) - ...] [species]`
/// Returns None for anything that does not begin with a `VFG` accession.
pub fn parse_vfdb_header(header: &str) -> Option<VfdbEntry> {
    let h = header.trim_start_matches('>').trim();
    if h.is_empty() {
        return None;
    }
    // The mmseqs target id is the first whitespace-delimited token.
    let target_id = h.split_whitespace().next().unwrap_or("").to_string();
    // VFG accession is the leading run up to the first '(' of the target id.
    let id_end = target_id.find('(').unwrap_or(target_id.len());
    let vfg_id = target_id[..id_end].trim().to_string();
    if !vfg_id.starts_with("VFG") {
        return None;
    }
    // Skip past the first parenthesis group in the FULL header (the `(gb|ACC)`).
    let first_open = h.find('(')?;
    let rest = &h[first_open..];
    let first_close = rest.find(')')?;
    let after_first = rest[first_close + 1..].trim_start();
    // Optional gene symbol is the next parenthesis group `(gene)`.
    let (gene, desc_part) = if let Some(stripped) = after_first.strip_prefix('(') {
        match stripped.find(')') {
            Some(c) => (
                Some(stripped[..c].trim().to_string()),
                stripped[c + 1..].trim_start(),
            ),
            None => (None, after_first),
        }
    } else {
        (None, after_first)
    };
    // Description runs up to the first bracketed metadata group ` [`.
    let description = match desc_part.find(" [") {
        Some(i) => desc_part[..i].trim(),
        None => desc_part.trim(),
    }
    .to_string();
    Some(VfdbEntry {
        target_id,
        vfg_id,
        gene: gene.filter(|g| !g.is_empty()),
        description,
    })
}

/// One parsed convertalis row that cleared the virulence-call gate.
#[derive(Clone, Debug)]
struct M8Hit {
    target: String,
    evalue: f64,
    bits: f32,
}

/// Parse a `query,target,pident,evalue,bits,qcov,tcov` m8 file into the best hit
/// per query that clears the virulence-call gate (evalue AND query coverage).
/// Best = highest bit score, tie-broken by lowest E-value.
fn parse_m8_best(path: &Path) -> HashMap<String, M8Hit> {
    let mut best: HashMap<String, M8Hit> = HashMap::new();
    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return best,
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 7 {
            continue;
        }
        let query = p[0].to_string();
        let pident: f64 = p[2].parse().unwrap_or(0.0);
        let evalue: f64 = parse_evalue(p[3]);
        let bits: f32 = p[4].parse().unwrap_or(0.0);
        let qcov: f64 = p[5].parse().unwrap_or(0.0);
        let tcov: f64 = p.get(6).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        // mmseqs reports pident as a fraction (0..1) with `--format-output pident`;
        // normalise a percent (>1) to a fraction so the MIN_PIDENT (%) compare holds.
        let pident_pct = if pident <= 1.0 { pident * 100.0 } else { pident };
        if !(evalue <= MAX_EVALUE
            && qcov >= MIN_QCOV
            && tcov >= MIN_TCOV
            && pident_pct >= MIN_PIDENT)
        {
            continue;
        }
        let hit = M8Hit { target: p[1].to_string(), evalue, bits };
        match best.get(&query) {
            Some(cur) if !is_better(&hit, cur) => {}
            _ => {
                best.insert(query, hit);
            }
        }
    }
    best
}

/// True if `a` is a better hit than `b` (higher bits, then lower evalue).
fn is_better(a: &M8Hit, b: &M8Hit) -> bool {
    if a.bits != b.bits {
        return a.bits > b.bits;
    }
    a.evalue < b.evalue
}

/// Ensure the VFDB is gunzipped to `<dir>/VFDB_setA_pro.fas` and mmseqs-indexed at
/// `<dir>/vfdb_setA_db` (both cached). Returns the paths (plain FASTA, mmseqs DB).
///
/// NOTE: the cache name is `vfdb_setA_db`, deliberately DISTINCT from the legacy
/// setB index (`vfdb_db`). A stale setB `vfdb_db.dbtype` may still sit in the same
/// directory; a filename-only switch would silently keep searching that old index,
/// so the setA index gets its own name.
fn ensure_db(bin: &str, vfdb_dir: &Path) -> Result<(std::path::PathBuf, std::path::PathBuf), String> {
    let (gz, fasta, db) = vfdb_paths(vfdb_dir);

    if !fasta.exists() {
        if !gz.exists() {
            return Err(format!("neither {} nor {} found", fasta.display(), gz.display()));
        }
        // Decompress in-process (pure Rust) to the sibling plain FASTA, keeping .gz.
        crate::util_io::gunzip_file(&gz, &fasta)?;
        if !fasta.exists() {
            return Err(format!("gunzip did not produce {}", fasta.display()));
        }
    }
    if !db.with_extension("dbtype").exists() {
        crate::mmseqs::run_with(bin, &["createdb", &fasta.to_string_lossy(), &db.to_string_lossy(), "-v", "1"])?;
    }
    Ok((fasta, db))
}

/// Build the `target_id -> VfdbEntry` lookup by scanning the plain-FASTA headers.
fn load_entries(fasta: &Path) -> HashMap<String, VfdbEntry> {
    let mut out = HashMap::new();
    let file = match fs::File::open(fasta) {
        Ok(f) => f,
        Err(_) => return out,
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(entry) = parse_vfdb_header(rest) {
                out.insert(entry.target_id.clone(), entry);
            }
        }
    }
    out
}

/// Annotate CDS with VFDB virulence-factor evidence. `vfdb_dir` holds
/// `VFDB_setA_pro.fas.gz` (core set). Degrades gracefully: on any mmseqs / IO error it logs
/// to stderr and returns Ok(()) so the pipeline continues.
pub fn annotate(features: &mut [Feature], vfdb_dir: &str, threads: usize) -> Result<(), String> {
    // 1. Candidate CDS: every CDS with a protein translation.
    let mut query_index: HashMap<String, usize> = HashMap::new();
    let mut fasta = String::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        if let Some(aa) = &f.aa {
            let aa = aa.trim();
            if aa.is_empty() {
                continue;
            }
            fasta.push('>');
            fasta.push_str(&f.id);
            fasta.push('\n');
            fasta.push_str(aa);
            fasta.push('\n');
            query_index.insert(f.id.clone(), i);
        }
    }
    if query_index.is_empty() {
        return Ok(());
    }

    let vfdb_dir = Path::new(vfdb_dir);
    if !crate::mmseqs::available() {
        eprintln!("[vfdb] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return Ok(());
    }
    let bin = crate::mmseqs::bin();

    // 2. Ensure the cached plain FASTA + mmseqs DB, then load header descriptions.
    let (plain_fasta, vfdb_db) = match ensure_db(&bin, vfdb_dir) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("[vfdb] cannot prepare VFDB under {}: {e} — skipping", vfdb_dir.display());
            return Ok(());
        }
    };
    let entries = load_entries(&plain_fasta);
    if entries.is_empty() {
        eprintln!("[vfdb] no VFDB entries parsed from {} — skipping", plain_fasta.display());
        return Ok(());
    }

    // 3. Per-process temp workspace (wipe stale copy: mmseqs refuses existing DBs).
    let work = std::env::temp_dir().join(format!("bactars_vfdb_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[vfdb] cannot create temp dir {}: {e}", work.display());
        return Ok(());
    }
    let best = run_search(&bin, &work, &vfdb_db, &fasta, threads);
    let _ = fs::remove_dir_all(&work);

    let best = match best {
        Ok(b) => b,
        Err(e) => {
            eprintln!("[vfdb] mmseqs pipeline failed: {e} — continuing without VFDB calls");
            return Ok(());
        }
    };

    // 4. Apply notes + provenance inferences.
    let mut calls = 0usize;
    for (qid, hit) in best {
        let Some(&idx) = query_index.get(&qid) else { continue };
        let Some(entry) = entries.get(&hit.target) else { continue };
        apply_vf(&mut features[idx], entry);
        calls += 1;
    }
    if calls > 0 {
        eprintln!("[vfdb] annotated {calls} CDS with VFDB virulence-factor evidence");
    }
    Ok(())
}

/// Add the virulence `/note` + `/inference` to a CDS (de-duped by VFG id). Never
/// overwrites the product.
fn apply_vf(feat: &mut Feature, entry: &VfdbEntry) {
    let already = feat
        .func
        .inferences
        .iter()
        .any(|inf| inf.db == "VFDB" && inf.accession == entry.vfg_id);
    if already {
        return;
    }
    let desc = if entry.description.is_empty() {
        entry.vfg_id.clone()
    } else {
        entry.description.clone()
    };
    feat.func.note.push(format!("virulence factor: {desc}"));
    feat.func.inferences.push(Inference {
        category: None,
        kind: "similar to AA sequence".to_string(),
        same_species: false,
        db: "VFDB".to_string(),
        accession: entry.vfg_id.clone(),
    });
}

/// createdb(query) -> search -> convertalis -> best-hit-per-query.
fn run_search(
    bin: &str,
    work: &Path,
    vfdb_db: &Path,
    fasta: &str,
    threads: usize,
) -> Result<HashMap<String, M8Hit>, String> {
    let query_faa = work.join("query.faa");
    let query_db = work.join("queryDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");
    let nt = threads.clamp(1, 8).to_string();

    fs::write(&query_faa, fasta).map_err(|e| format!("write query.faa: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_faa, q_db, res_db, tmp_s, vfdb_s, m8_s) =
        (s(&query_faa), s(&query_db), s(&result_db), s(&tmp), s(vfdb_db), s(&m8));

    crate::mmseqs::run_with(bin, &["createdb", &q_faa, &q_db, "-v", "1"])?;
    crate::mmseqs::run_with(
        bin,
        &[
            "search", &q_db, &vfdb_s, &res_db, &tmp_s, "--threads", &nt, "-s", SEARCH_SENSITIVITY,
            "--max-seqs", MAX_SEQS, "-e", EVALUE_CUTOFF, "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        bin,
        &[
            "convertalis", &q_db, &vfdb_s, &res_db, &m8_s, "--format-output",
            "query,target,pident,evalue,bits,qcov,tcov", "--threads", &nt, "-v", "1",
        ],
    )?;
    Ok(parse_m8_best(&m8))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureKind, Functional};
    use std::io::Write as _;

    fn cds(id: &str, aa: &str) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c1".to_string(),
            id: id.to_string(),
            start: 1,
            end: 99,
            strand: 1,
            aa: Some(aa.to_string()),
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        }
    }

    #[test]
    fn vfdb_paths_use_seta_core_and_distinct_cache_name() {
        let (gz, fasta, db) = vfdb_paths(Path::new("/db/vfdb"));
        assert!(gz.ends_with("VFDB_setA_pro.fas.gz"), "{gz:?}");
        assert!(fasta.ends_with("VFDB_setA_pro.fas"), "{fasta:?}");
        // Distinct from the stale setB index name `vfdb_db`.
        assert!(db.ends_with("vfdb_setA_db"), "{db:?}");
        assert_ne!(db.file_name().unwrap(), "vfdb_db");
    }

    #[test]
    fn parse_header_full_grammar() {
        let h = ">VFG037170(gb|WP_001081754) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii 1656-2]";
        let e = parse_vfdb_header(h).expect("parses");
        assert_eq!(e.target_id, "VFG037170(gb|WP_001081754)");
        assert_eq!(e.vfg_id, "VFG037170");
        assert_eq!(e.gene.as_deref(), Some("plc1"));
        assert_eq!(e.description, "phospholipase C");
    }

    #[test]
    fn parse_header_without_gene_group() {
        // No second parenthesis (no gene symbol) — description follows the id paren.
        let h = "VFG000123(gb|ABC12345) some hypothetical virulence protein [Cat (VF0001)] [E. coli]";
        let e = parse_vfdb_header(h).expect("parses");
        assert_eq!(e.vfg_id, "VFG000123");
        assert_eq!(e.gene, None);
        assert_eq!(e.description, "some hypothetical virulence protein");
    }

    #[test]
    fn parse_header_rejects_non_vfg() {
        assert!(parse_vfdb_header(">sp|P12345|SOMETHING description").is_none());
        assert!(parse_vfdb_header(">").is_none());
        assert!(parse_vfdb_header("").is_none());
    }

    #[test]
    fn m8_best_hit_respects_evalue_and_qcov_gate() {
        let dir = std::env::temp_dir().join(format!("bactars_vfdb_m8_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        let m8 = dir.join("result.m8");
        let mut fh = fs::File::create(&m8).unwrap();
        // q1: two hits both passing the gate; higher bits wins.
        writeln!(fh, "q1\tVFG_A(gb|x)\t90.0\t1e-40\t200.0\t0.95\t0.98").unwrap();
        writeln!(fh, "q1\tVFG_B(gb|y)\t85.0\t1e-30\t150.0\t0.90\t0.92").unwrap();
        // q2: strong evalue but LOW query coverage -> rejected (partial/domain hit).
        writeln!(fh, "q2\tVFG_C(gb|z)\t99.0\t1e-50\t300.0\t0.30\t0.99").unwrap();
        // q3: good coverage but weak evalue -> rejected.
        writeln!(fh, "q3\tVFG_D(gb|w)\t40.0\t1e-5\t60.0\t0.90\t0.90").unwrap();
        drop(fh);

        let best = parse_m8_best(&m8);
        assert_eq!(best.get("q1").unwrap().target, "VFG_A(gb|x)");
        assert!(!best.contains_key("q2"));
        assert!(!best.contains_key("q3"));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn apply_vf_sets_note_and_inference_no_dupes_no_product_clobber() {
        let mut f = cds("cds1", "MKV");
        let e = VfdbEntry {
            target_id: "VFG037170(gb|WP_001081754)".to_string(),
            vfg_id: "VFG037170".to_string(),
            gene: Some("plc1".to_string()),
            description: "phospholipase C".to_string(),
        };
        apply_vf(&mut f, &e);
        apply_vf(&mut f, &e); // de-dupe: second call is a no-op.
        assert_eq!(f.func.note.len(), 1);
        assert_eq!(f.func.note[0], "virulence factor: phospholipase C");
        assert_eq!(f.func.inferences.len(), 1);
        assert_eq!(f.func.inferences[0].db, "VFDB");
        assert_eq!(f.func.inferences[0].accession, "VFG037170");
        // Product must be untouched.
        assert!(f.func.product.is_none());
    }
}
