//! PSC (Protein Sequence Cluster) naming: name unannotated CDS by protein
//! similarity to the native UniRef90 bacterial/archaeal/viral mmseqs DB.
//!
//! For each CDS with no curated/HMM product yet, run an mmseqs search against
//! `psc_db`, take the best hit's UniRef90 id, and set func.product from
//! `names.tsv` (UniRef90_ID -> product). This closes the ~25% unannotated-CDS gap
//! (Bakta's PSC tier). mmseqs must be on PATH (or full path); the DB + names.tsv
//! live under the psc bundle dir.
//!
//! Pipeline: write unnamed CDS proteins to a temp FASTA -> `mmseqs createdb` ->
//! `mmseqs search` (mmaps the precomputed psc_db.idx, ~194 GB) -> `convertalis`
//! to an m8 table -> best hit per query -> streaming lookup of the hit target
//! ids in names.tsv -> set func.product + provenance Annotation. Any mmseqs
//! failure is logged and swallowed (returns Ok) so the pipeline never crashes.

use crate::feature::{Annotation, Feature, FeatureKind};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

/// mmseqs search sensitivity (Bakta-style protein naming is fine at high-ish -s).
const SEARCH_SENSITIVITY: &str = "5.7";
/// Max target sequences returned per query.
const MAX_SEQS: &str = "25";
/// Reporting / naming E-value cutoff.
const EVALUE_CUTOFF: &str = "1e-6";
/// Accept a hit for naming only if percent-identity is at least this.
const MIN_PIDENT: f64 = 50.0;
/// ...AND the E-value is at least this good (both gates must hold; see
/// [`parse_m8_best`] for why the identity floor was dead when OR'd).
const MAX_EVALUE: f64 = 1e-6;

/// Emit a single warning the first time an m8 e-value column fails to parse, so a
/// malformed row is not silently coerced to +inf without any trace.
static WARNED_EVALUE: AtomicBool = AtomicBool::new(false);
fn parse_evalue(raw: &str) -> f64 {
    match raw.parse() {
        Ok(v) => v,
        Err(_) => {
            if !WARNED_EVALUE.swap(true, Ordering::Relaxed) {
                eprintln!(
                    "[psc] warning: malformed e-value '{raw}' in mmseqs m8 output — treating as +inf (further such rows silenced)"
                );
            }
            f64::INFINITY
        }
    }
}

/// A CDS needs PSC naming when it is a real protein-coding feature that is still
/// effectively unnamed: no resolved product (or only the "hypothetical protein"
/// placeholder) and no strong existing HMM/xref annotation that already yields a
/// real display name. `display_product()` folds in the best raw annotation, so a
/// CDS that any prior tier named will already report a non-placeholder product.
pub fn needs_psc(feat: &Feature) -> bool {
    if feat.kind != FeatureKind::Cds {
        return false;
    }
    // Must have a protein translation to search with.
    if feat.aa.as_deref().map(|s| s.trim().is_empty()).unwrap_or(true) {
        return false;
    }
    match feat.display_product() {
        None => true,
        Some(p) => is_placeholder_name(&p),
    }
}

/// True for names that carry no real functional information (treated as unnamed).
fn is_placeholder_name(name: &str) -> bool {
    let n = name.trim().to_ascii_lowercase();
    n.is_empty()
        || n == "hypothetical protein"
        || n == "uncharacterized protein"
        || n == "putative protein"
        || n == "protein of unknown function"
}

/// Normalise a raw UniRef90 product string into a display product. Uninformative
/// UniRef names collapse to the standard `hypothetical protein` placeholder so we
/// never advertise "Uncharacterized protein" as a real name.
pub fn clean_product(raw: &str) -> String {
    let p = raw.trim();
    if p.is_empty() || is_placeholder_name(p) {
        return "hypothetical protein".to_string();
    }
    p.to_string()
}

/// One parsed m8 alignment row (only the fields we requested).
#[derive(Clone, Debug)]
struct M8Hit {
    target: String,
    #[allow(dead_code)] // retained for provenance/debug; gating uses it inline
    pident: f64,
    evalue: f64,
    bits: f32,
    /// Query (this CDS/sORF protein) length in aa, from mmseqs `qlen`. `0` when
    /// the column is absent (old-format m8).
    qlen: usize,
    /// Target (UniRef90 reference protein) length in aa, from mmseqs `tlen`. `0`
    /// when the column is absent. Feeds `Annotation.ref_len` for truncation.
    tlen: usize,
    /// Fraction (0..1) of the query covered by the alignment (`qcov`). `0` when the
    /// column is absent. Used by the sORF PSC gate.
    qcov: f64,
    /// Fraction (0..1) of the TARGET covered by the alignment (`tcov`). `0` when the
    /// column is absent. The sORF PSC gate requires this to be high so a short ORF
    /// must match a similarly-short reference full-length, not a fragment of a large
    /// protein (the dominant sORF false-positive source).
    tcov: f64,
}

/// A resolved PSC hit for one query: the named UniRef90 target plus the scores and
/// the reference/query lengths. Returned by [`search`] so callers (PSC naming and
/// the sORF PSC path) can gate on coverage and record the reference length.
#[derive(Clone, Debug)]
pub struct PscHit {
    /// UniRef90 target id (the best hit).
    pub uniref_id: String,
    /// Cleaned display product from `names.tsv` (never the raw placeholder).
    pub product: String,
    /// Bit score of the best hit.
    pub bits: f32,
    /// E-value of the best hit.
    pub evalue: f64,
    /// Query protein length in aa (`0` if unknown).
    pub qlen: usize,
    /// Target (reference) protein length in aa (`0` if unknown).
    pub tlen: usize,
    /// Fraction (0..1) of the query covered by the alignment (`0` if unknown).
    pub qcov: f64,
    /// Fraction (0..1) of the target covered by the alignment (`0` if unknown).
    pub tcov: f64,
}

/// Parse a `query,target,pident,evalue,bits` m8 file into best-hit-per-query.
/// Best = highest bit score, tie-broken by lowest E-value. Only hits passing the
/// naming floor (pident >= MIN_PIDENT AND evalue <= MAX_EVALUE) are considered.
///
/// The gate is an AND, not an OR: the mmseqs search already ran with `-e
/// EVALUE_CUTOFF`, so every row here already satisfies the E-value arm and an OR
/// made the identity floor dead code. AND'ing makes MIN_PIDENT actually bite.
fn parse_m8_best(path: &Path) -> HashMap<String, M8Hit> {
    let mut best: HashMap<String, M8Hit> = HashMap::new();
    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return best,
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 5 {
            continue;
        }
        let query = p[0].to_string();
        let pident_raw: f64 = p[2].parse().unwrap_or(0.0);
        // mmseqs may report pident as a fraction (0..1); normalise to a percent so
        // the MIN_PIDENT (%) compare is consistent with the other stages (F8).
        let pident = if pident_raw <= 1.0 { pident_raw * 100.0 } else { pident_raw };
        let evalue: f64 = parse_evalue(p[3]);
        let bits: f32 = p[4].parse().unwrap_or(0.0);
        if !(pident >= MIN_PIDENT && evalue <= MAX_EVALUE) {
            continue;
        }
        // Appended length columns (qlen, tlen). Absent in old-format m8 -> 0.
        let qlen: usize = p.get(5).and_then(|s| s.parse().ok()).unwrap_or(0);
        let tlen: usize = p.get(6).and_then(|s| s.parse().ok()).unwrap_or(0);
        let qcov: f64 = p.get(7).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let tcov: f64 = p.get(8).and_then(|s| s.parse().ok()).unwrap_or(0.0);
        let hit = M8Hit { target: p[1].to_string(), pident, evalue, bits, qlen, tlen, qcov, tcov };
        match best.get(&query) {
            Some(cur) if !is_better(&hit, cur) => {}
            _ => {
                best.insert(query, hit);
            }
        }
    }
    best
}

/// True if `a` is a better naming hit than `b` (higher bits, then lower evalue).
fn is_better(a: &M8Hit, b: &M8Hit) -> bool {
    if a.bits != b.bits {
        return a.bits > b.bits;
    }
    a.evalue < b.evalue
}

/// Stream `names.tsv` (UniRef90_ID<TAB>product) once, keeping only the products
/// for ids in `wanted`. Never loads all 73M rows — bounded by |wanted|.
fn lookup_names(names_tsv: &Path, wanted: &HashSet<String>) -> HashMap<String, String> {
    let mut out: HashMap<String, String> = HashMap::new();
    if wanted.is_empty() {
        return out;
    }
    let file = match fs::File::open(names_tsv) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("[psc] cannot open names.tsv {}: {e}", names_tsv.display());
            return out;
        }
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        // Split only on the first tab: product text may itself contain tabs.
        if let Some((id, product)) = line.split_once('\t') {
            if wanted.contains(id) {
                out.insert(id.to_string(), product.to_string());
                if out.len() == wanted.len() {
                    break; // found them all; stop scanning the 3.7 GB file
                }
            }
        }
    }
    out
}

/// Name unannotated CDS features by UniRef90 protein similarity. `psc_dir` holds
/// `psc_db*` (mmseqs) + `names.tsv`. Degrades gracefully: on any mmseqs error it
/// logs to stderr and returns Ok(()) so the pipeline continues.
pub fn annotate(features: &mut [Feature], psc_dir: &str, threads: usize) -> Result<(), String> {
    // 1. Select the genuinely unnamed CDS and map query-id -> feature index.
    let mut query_index: HashMap<String, usize> = HashMap::new();
    let mut fasta = String::new();
    for (i, f) in features.iter().enumerate() {
        if needs_psc(f) {
            if let Some(aa) = &f.aa {
                fasta.push('>');
                fasta.push_str(&f.id);
                fasta.push('\n');
                fasta.push_str(aa.trim());
                fasta.push('\n');
                query_index.insert(f.id.clone(), i);
            }
        }
    }
    if query_index.is_empty() {
        return Ok(()); // nothing to name
    }

    let psc_dir = Path::new(psc_dir);
    let psc_db = psc_dir.join("psc_db");
    let names_tsv = psc_dir.join("names.tsv");
    if !psc_db.with_extension("dbtype").exists() {
        eprintln!("[psc] psc_db not found under {} — skipping PSC naming", psc_dir.display());
        return Ok(());
    }
    if !crate::mmseqs::available() {
        eprintln!("[psc] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return Ok(());
    }

    // 2. Per-process temp dir (wipe any stale copy: mmseqs refuses existing DBs).
    let work = std::env::temp_dir().join(format!("bactars_psc_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[psc] cannot create temp dir {}: {e}", work.display());
        return Ok(());
    }

    let result = run_psc_pipeline(&work, &psc_db, &names_tsv, &fasta, threads);

    // Always clean up the temp workspace.
    let _ = fs::remove_dir_all(&work);

    let names = match result {
        Ok(n) => n,
        Err(e) => {
            eprintln!("[psc] mmseqs pipeline failed: {e} — continuing without PSC names");
            return Ok(());
        }
    };

    // 6. Apply names + provenance. `names` maps query-id -> PscHit.
    let mut named = 0usize;
    for (qid, hit) in names {
        let Some(&idx) = query_index.get(&qid) else { continue };
        let feat = &mut features[idx];
        let cleaned = clean_product(&hit.product);
        // Record provenance regardless (the homology hit is real evidence). The
        // target length is carried as `ref_len` so the pseudogene truncation
        // signal can use it when the CDS has no ncbifams `hmm_length`.
        feat.annotations.push(Annotation {
            source: "psc:uniref90".to_string(),
            accession: hit.uniref_id,
            name: cleaned.clone(),
            score: hit.bits,
            evalue: Some(hit.evalue),
            ref_len: (hit.tlen > 0).then_some(hit.tlen),
        });
        // Only fill product if still empty (don't clobber a prior tier).
        if feat.func.product.is_none() {
            feat.func.product = Some(cleaned);
            named += 1;
        }
    }
    if named > 0 {
        eprintln!("[psc] named {named} previously-unannotated CDS via UniRef90");
    }
    Ok(())
}

/// Runs createdb -> search -> convertalis -> parse -> names lookup. Returns
/// query-id -> [`PscHit`] for hits that resolved a name.
fn run_psc_pipeline(
    work: &Path,
    psc_db: &Path,
    names_tsv: &Path,
    fasta: &str,
    threads: usize,
) -> Result<HashMap<String, PscHit>, String> {
    let bin = crate::mmseqs::bin();
    let query_faa = work.join("query.faa");
    let query_db = work.join("queryDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");
    let nt = threads.max(1).to_string();

    fs::write(&query_faa, fasta).map_err(|e| format!("write query.faa: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_faa, q_db, res_db, tmp_s, psc_s, m8_s) =
        (s(&query_faa), s(&query_db), s(&result_db), s(&tmp), s(psc_db), s(&m8));

    crate::mmseqs::run_with(&bin, &["createdb", &q_faa, &q_db, "-v", "1"])?;
    crate::mmseqs::run_with(
        &bin,
        &[
            "search", &q_db, &psc_s, &res_db, &tmp_s,
            "--threads", &nt,
            "-s", SEARCH_SENSITIVITY,
            "--max-seqs", MAX_SEQS,
            "-e", EVALUE_CUTOFF,
            "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        &bin,
        &[
            "convertalis", &q_db, &psc_s, &res_db, &m8_s,
            // Append qlen,tlen AFTER the original 5 columns so existing parsing of
            // the first 5 (and thus named-CDS behaviour) is byte-identical.
            "--format-output", "query,target,pident,evalue,bits,qlen,tlen,qcov,tcov",
            "--threads", &nt,
            "-v", "1",
        ],
    )?;

    // 4. Best hit per query.
    let best = parse_m8_best(&m8);
    if best.is_empty() {
        return Ok(HashMap::new());
    }

    // 5. Streaming names.tsv lookup for exactly the targets we hit.
    let wanted: HashSet<String> = best.values().map(|h| h.target.clone()).collect();
    let names = lookup_names(names_tsv, &wanted);

    let mut out: HashMap<String, PscHit> = HashMap::new();
    for (query, hit) in best {
        if let Some(product) = names.get(&hit.target) {
            out.insert(
                query,
                PscHit {
                    uniref_id: hit.target.clone(),
                    product: clean_product(product),
                    bits: hit.bits,
                    evalue: hit.evalue,
                    qlen: hit.qlen,
                    tlen: hit.tlen,
                    qcov: hit.qcov,
                    tcov: hit.tcov,
                },
            );
        }
    }
    Ok(out)
}

/// Reusable PSC search over an arbitrary set of proteins. Runs the SAME mmseqs
/// pipeline as [`annotate`] (createdb -> search -> convertalis) on every feature in
/// `features` that carries a protein translation, and returns `feature.id -> `
/// [`PscHit`] for the queries that resolved a UniRef90 name. Unlike [`annotate`] it
/// does not select by `needs_psc` and does not mutate the features — the caller
/// decides what to do with the hits (e.g. the sORF stage keeps only candidates with
/// good query coverage). Degrades gracefully: on any mmseqs/setup error it logs to
/// stderr and returns an empty map so the pipeline never crashes.
pub fn search(features: &[Feature], psc_dir: &str, threads: usize) -> HashMap<String, PscHit> {
    let mut fasta = String::new();
    let mut n_query = 0usize;
    for f in features {
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
            n_query += 1;
        }
    }
    if n_query == 0 {
        return HashMap::new();
    }

    let psc_dir = Path::new(psc_dir);
    let psc_db = psc_dir.join("psc_db");
    let names_tsv = psc_dir.join("names.tsv");
    if !psc_db.with_extension("dbtype").exists() {
        eprintln!("[psc] psc_db not found under {} — skipping PSC search", psc_dir.display());
        return HashMap::new();
    }
    if !crate::mmseqs::available() {
        eprintln!("[psc] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return HashMap::new();
    }

    // Per-process temp dir; suffix distinguishes it from the naming pass's dir so a
    // concurrent naming run never collides on the same path.
    let work = std::env::temp_dir().join(format!("bactars_pscsearch_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[psc] cannot create temp dir {}: {e}", work.display());
        return HashMap::new();
    }

    let result = run_psc_pipeline(&work, &psc_db, &names_tsv, &fasta, threads);
    let _ = fs::remove_dir_all(&work);

    match result {
        Ok(m) => m,
        Err(e) => {
            eprintln!("[psc] search pipeline failed: {e} — continuing without PSC hits");
            HashMap::new()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureKind, Functional};
    use std::io::Write as _;

    fn cds(id: &str, aa: Option<&str>, product: Option<&str>) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c1".to_string(),
            id: id.to_string(),
            start: 1,
            end: 99,
            strand: 1,
            aa: aa.map(|s| s.to_string()),
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional {
                product: product.map(|s| s.to_string()),
                ..Default::default()
            },
        }
    }

    #[test]
    fn selection_picks_only_genuinely_unnamed_cds() {
        // Unnamed CDS with a protein -> needs PSC.
        assert!(needs_psc(&cds("g1", Some("MKV"), None)));
        // Placeholder product -> still needs PSC.
        assert!(needs_psc(&cds("g2", Some("MKV"), Some("hypothetical protein"))));
        assert!(needs_psc(&cds("g3", Some("MKV"), Some("Uncharacterized protein"))));
        // Already has a real product -> skip.
        assert!(!needs_psc(&cds("g4", Some("MKV"), Some("DNA gyrase subunit A"))));
        // No protein translation -> cannot search, skip.
        assert!(!needs_psc(&cds("g5", None, None)));
        assert!(!needs_psc(&cds("g6", Some("   "), None)));
    }

    #[test]
    fn selection_skips_non_cds_and_hmm_named() {
        let mut trna = cds("t1", Some("MKV"), None);
        trna.kind = FeatureKind::Trna;
        assert!(!needs_psc(&trna));

        // A CDS with no func.product but a strong HMM annotation -> display_product
        // returns the HMM name, so PSC is not needed.
        let mut f = cds("g7", Some("MKV"), None);
        f.annotations.push(Annotation {
            source: "rustyhmmer:pfam".to_string(),
            accession: "PF00521".to_string(),
            name: "DNA gyrase subunit A".to_string(),
            score: 250.0,
            evalue: Some(1e-40),
            ref_len: None,
        });
        assert!(!needs_psc(&f));
    }

    #[test]
    fn clean_product_collapses_uninformative_names() {
        assert_eq!(clean_product("  LacI family regulator "), "LacI family regulator");
        assert_eq!(clean_product("Uncharacterized protein"), "hypothetical protein");
        assert_eq!(clean_product(""), "hypothetical protein");
        assert_eq!(clean_product("hypothetical protein"), "hypothetical protein");
    }

    #[test]
    fn names_lookup_streams_only_wanted_ids() {
        let dir = std::env::temp_dir().join(format!("bactars_psc_test_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        let tsv = dir.join("names.tsv");
        let mut fh = fs::File::create(&tsv).unwrap();
        // product text intentionally contains a tab to exercise split_once.
        writeln!(fh, "UniRef90_AAA\tLacI family DNA-binding transcriptional regulator").unwrap();
        writeln!(fh, "UniRef90_BBB\tUncharacterized protein").unwrap();
        writeln!(fh, "UniRef90_CCC\tM3 family\tmetallopeptidase").unwrap();
        writeln!(fh, "UniRef90_DDD\tShould not be picked").unwrap();
        drop(fh);

        let mut wanted = HashSet::new();
        wanted.insert("UniRef90_AAA".to_string());
        wanted.insert("UniRef90_CCC".to_string());
        wanted.insert("UniRef90_ZZZ".to_string()); // absent

        let got = lookup_names(&tsv, &wanted);
        assert_eq!(got.len(), 2);
        assert_eq!(
            got.get("UniRef90_AAA").map(String::as_str),
            Some("LacI family DNA-binding transcriptional regulator")
        );
        assert_eq!(
            got.get("UniRef90_CCC").map(String::as_str),
            Some("M3 family\tmetallopeptidase")
        );
        assert!(!got.contains_key("UniRef90_DDD"));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn m8_best_hit_by_bits_then_evalue() {
        let dir = std::env::temp_dir().join(format!("bactars_psc_m8_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        let m8 = dir.join("result.m8");
        let mut fh = fs::File::create(&m8).unwrap();
        // q1: two hits, both clear the AND gate (pident>=50 && evalue<=1e-6);
        // higher bits should win.
        writeln!(fh, "q1\tUniRef90_LOW\t60.0\t1e-20\t120.0").unwrap();
        writeln!(fh, "q1\tUniRef90_HIGH\t80.0\t1e-40\t250.0").unwrap();
        // q2: fails both arms (low pident, weak evalue) -> excluded.
        writeln!(fh, "q2\tUniRef90_BAD\t20.0\t1e-3\t30.0").unwrap();
        // q3: strong evalue but pident < MIN_PIDENT -> now excluded by the AND gate
        // (previously an OR let this through; the identity floor is real now, F4).
        writeln!(fh, "q3\tUniRef90_EOK\t35.0\t1e-10\t70.0").unwrap();
        // q4: passes on identity AND evalue -> kept.
        writeln!(fh, "q4\tUniRef90_OK\t72.0\t1e-15\t140.0").unwrap();
        drop(fh);

        let best = parse_m8_best(&m8);
        assert_eq!(best.get("q1").unwrap().target, "UniRef90_HIGH");
        assert!(!best.contains_key("q2"));
        assert!(!best.contains_key("q3"));
        assert_eq!(best.get("q4").unwrap().target, "UniRef90_OK");
        let _ = fs::remove_dir_all(&dir);
    }
}
