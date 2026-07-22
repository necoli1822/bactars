//! Integron detection (IntegronFinder-style), via rustyhmmer + infernox.
//!
//! An integron is a genetic element made of an **integron-integrase gene**
//! (`intI`, a specific sub-family of tyrosine site-specific recombinase) and an
//! adjacent array of **attC recombination sites** into which gene cassettes are
//! captured. IntegronFinder (Cury et al. 2016; github.com/gem-pasteur/
//! Integron_Finder) locates both signals independently and clusters them:
//!
//!   1. Scan CDS proteins with the intI HMM(s) (`IntI.hmm` /
//!      `integron_integrase.hmm`, the C-terminal intI motif) AND the generic
//!      tyrosine-recombinase HMM (`phage-int.hmm`, Pfam `PF00589`). A gene is an
//!      integron-integrase only when it hits BOTH — the intI motif discriminates
//!      the integron sub-family from ordinary phage integrases.
//!   2. Scan the CONTIG nucleotides with the attC covariance model
//!      (`attc_4.cm`) via infernox `cmsearch` -> attC sites.
//!   3. Cluster: an `intI` gene + >=1 nearby attC = a **complete** integron; an
//!      `intI` gene with no nearby attC = **In0** (integrase only); a cluster of
//!      >=2 attC sites with no nearby integrase = **CALIN** (clusters of attC-
//!      lacking integrase).
//!
//! Emits `Feature { kind: FeatureKind::Integron }` spanning the integrase + its
//! attC array. Model files come from the `db/integron/` dir
//! (IntegronFinder v2.0.2, GPL-3.0).

use std::collections::HashSet;
use std::io::Cursor;
use std::path::Path;

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use crate::hmm::{self, Cutoff};
use infernox::cm_file::cm_file_read_from_reader_opt;
use infernox::cm_search::{FaithfulConfig, FaithfulSearcher};

/// Max distance (bp) between an `intI` gene and an attC site (nearest edges) for
/// the site to be considered part of that integrase's array. IntegronFinder's
/// default cassette/attC search window around the integrase is ~4 kb.
const INTEGRASE_ATTC_WINDOW: i64 = 4_000;

/// Max gap (bp) between consecutive attC sites for them to belong to one array
/// (used both for complete-integron arrays and standalone CALIN clusters).
const ATTC_CLUSTER_GAP: i64 = 4_000;

/// Minimum number of clustered attC sites to call a standalone CALIN (no
/// integrase). A lone attC hit is too weak to report on its own.
const CALIN_MIN_ATTC: usize = 2;

/// attC cmsearch reporting E-value. attC sites are short and somewhat degenerate;
/// this is stricter than IntegronFinder's permissive `-E 1` because bactars calls
/// standalone CALINs (no integrase anchor) and favours precision.
const ATTC_EVALUE: f64 = 1e-3;

/// A genomic locus (1-based inclusive), used for both intI genes and attC sites.
#[derive(Clone, Debug)]
struct Locus {
    contig: String,
    start: i64,
    end: i64,
    strand: i8,
    /// Real infernox cmsearch E-value for an attC covariance-model site; `None`
    /// for an integrase locus (which comes from a CDS/HMM hit, not this search).
    evalue: Option<f64>,
}

impl Locus {
    /// Gap (bp) between two loci on the same contig: 0 if they overlap, else the
    /// distance between their nearest edges. `i64::MAX` if on different contigs.
    fn gap(&self, other: &Locus) -> i64 {
        if self.contig != other.contig {
            return i64::MAX;
        }
        if self.end < other.start {
            other.start - self.end
        } else if other.end < self.start {
            self.start - other.end
        } else {
            0
        }
    }
}

/// Integron classification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum IntegronClass {
    /// intI integrase + >=1 attC site.
    Complete,
    /// intI integrase, no nearby attC.
    In0,
    /// Cluster of attC sites, no nearby integrase.
    Calin,
}

impl IntegronClass {
    fn label(self) -> &'static str {
        match self {
            IntegronClass::Complete => "complete integron",
            IntegronClass::In0 => "In0 (integron-integrase, no attC)",
            IntegronClass::Calin => "CALIN (attC cluster, no integrase)",
        }
    }
    fn short(self) -> &'static str {
        match self {
            IntegronClass::Complete => "complete",
            IntegronClass::In0 => "In0",
            IntegronClass::Calin => "CALIN",
        }
    }
}

/// One assembled integron call (contig-local coordinates, 1-based inclusive).
#[derive(Clone, Debug)]
struct IntegronCall {
    contig: String,
    start: i64,
    end: i64,
    strand: i8,
    class: IntegronClass,
    n_attc: usize,
    /// Best (smallest) attC-site cmsearch E-value among the sites composing this
    /// integron; `None` for In0 (integrase only, no attC site).
    attc_evalue: Option<f64>,
}

/// Detect integrons across `contigs` (`(name, seq)`), using the CDS features in
/// `features` for the integrase search and the contig nucleotides for the attC
/// covariance search. `integron_dir` holds `IntI.hmm`, `integron_integrase.hmm`,
/// `attc_4.cm`, and (optionally) `phage-int.hmm`.
pub fn detect(
    contigs: &[(String, String)],
    features: &[Feature],
    integron_dir: &str,
) -> Result<Vec<Feature>, String> {
    let dir = integron_dir.trim_end_matches('/');

    // --- 1. Integron-integrase genes (intI) via rustyhmmer on CDS proteins. ---
    let intis_res = detect_integrases(features, dir);
    // --- 2. attC recombination sites via infernox cmsearch on the contigs. ---
    let attcs_res = detect_attc(contigs, dir);

    // The intI HMM and the attC covariance model are BOTH required models. If one
    // is missing but the other runs, degrade gracefully (integrase-only In0 calls,
    // or attC-only CALIN calls) with a warning. Only when NEITHER model could run
    // — the detector's DB is missing/unusable — is this a whole-stage failure that
    // is surfaced to the caller as an Err (consistent with trna/is/ncrna).
    let (intis, attcs) = match (intis_res, attcs_res) {
        (Ok(i), Ok(a)) => (i, a),
        (Ok(i), Err(e)) => {
            eprintln!("integron: attC search unavailable ({e}); reporting integrase-only calls");
            (i, Vec::new())
        }
        (Err(e), Ok(a)) => {
            eprintln!("integron: integrase search unavailable ({e}); reporting attC-only (CALIN) calls");
            (Vec::new(), a)
        }
        (Err(ei), Err(ea)) => {
            return Err(format!(
                "no usable integron models in {dir} (intI: {ei}; attC: {ea})"
            ));
        }
    };

    if intis.is_empty() && attcs.is_empty() {
        return Ok(Vec::new());
    }

    // --- 3. Cluster integrase + attC into integron calls. ---
    let calls = assemble_integrons(&intis, &attcs);

    // --- 4. Emit features (per-contig numbering, contig then start order). ---
    let mut calls = calls;
    calls.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));
    let mut out = Vec::new();
    let mut per_contig: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for c in calls {
        let n = per_contig.entry(c.contig.clone()).or_insert(0);
        *n += 1;
        out.push(Feature {
            kind: FeatureKind::Integron,
            contig: c.contig.clone(),
            id: format!("{}_integron{}", c.contig, n),
            start: c.start,
            end: c.end,
            strand: c.strand,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: format!("integron:{}", c.class.short()),
                accession: String::new(),
                name: format!("{} ({} attC)", c.class.label(), c.n_attc),
                score: c.n_attc as f32,
                // Best attC-site cmsearch E-value across the integron's sites;
                // `None` for In0 (integrase only, no attC search hit).
                evalue: c.attc_evalue,
                ref_len: None,
            }],
            func: Functional::default(),
        });
    }
    Ok(out)
}

/// Find integron-integrase genes: CDS that hit an intI motif HMM AND (when the
/// generic phage-integrase HMM is present) the `PF00589` tyrosine-recombinase
/// HMM. Returns each integrase gene as a genomic `Locus`.
fn detect_integrases(features: &[Feature], dir: &str) -> Result<Vec<Locus>, String> {
    // CDS proteins (id, aa) and an id -> feature index for coordinate lookup.
    let proteins: Vec<(String, String)> = features
        .iter()
        .filter(|f| f.kind == FeatureKind::Cds)
        .filter_map(|f| f.aa.as_ref().map(|aa| (f.id.clone(), aa.clone())))
        .collect();
    if proteins.is_empty() {
        return Ok(Vec::new());
    }

    // intI motif HMM(s): either/both of the shipped intI files.
    let mut inti_ids: HashSet<String> = HashSet::new();
    let mut any_inti_hmm = false;
    for name in ["IntI.hmm", "integron_integrase.hmm"] {
        let path = format!("{dir}/{name}");
        if !Path::new(&path).exists() {
            continue;
        }
        any_inti_hmm = true;
        // intI HMMs carry a curated GA cutoff (29.90) — use it (`--cut_ga`).
        for h in hmm::annotate(&proteins, &path, Cutoff::GatheringGa)? {
            inti_ids.insert(h.target_name);
        }
    }
    if !any_inti_hmm {
        return Err(format!("no intI HMM (IntI.hmm) in {dir}"));
    }

    // Generic tyrosine-recombinase HMM (PF00589). IntegronFinder requires an
    // integrase gene to hit BOTH the intI motif and PF00589; when the file is
    // present we intersect, otherwise the intI motif alone stands.
    let phage_path = format!("{dir}/phage-int.hmm");
    let integrase_ids: HashSet<String> = if Path::new(&phage_path).exists() {
        let mut phage_ids: HashSet<String> = HashSet::new();
        for h in hmm::annotate(&proteins, &phage_path, Cutoff::GatheringGa)? {
            phage_ids.insert(h.target_name);
        }
        inti_ids.intersection(&phage_ids).cloned().collect()
    } else {
        inti_ids
    };

    // Map the surviving CDS ids back to genomic loci.
    let mut out = Vec::new();
    for f in features.iter().filter(|f| f.kind == FeatureKind::Cds) {
        if integrase_ids.contains(&f.id) {
            let (start, end) = if f.start <= f.end {
                (f.start, f.end)
            } else {
                (f.end, f.start)
            };
            out.push(Locus {
                contig: f.contig.clone(),
                start,
                end,
                strand: f.strand,
                evalue: None,
            });
        }
    }
    Ok(out)
}

/// Search the contig nucleotides with `attc_4.cm` via infernox cmsearch and
/// return each attC site (E-value <= [`ATTC_EVALUE`]) as a `Locus`.
fn detect_attc(contigs: &[(String, String)], dir: &str) -> Result<Vec<Locus>, String> {
    let cm_path = format!("{dir}/attc_4.cm");
    if !Path::new(&cm_path).exists() {
        return Err(format!("no attc_4.cm in {dir}"));
    }
    let text = std::fs::read_to_string(&cm_path).map_err(|e| format!("{cm_path}: {e}"))?;
    // Read in GLOBAL config (do_localize=false); `FaithfulSearcher::new` expects
    // an already-global CM — see the note in `ncrna::detect_ncrna`.
    let cm = cm_file_read_from_reader_opt(Cursor::new(text.as_bytes()), false)
        .map_err(|e| format!("{cm_path}: {e:?}"))?;
    let searcher = FaithfulSearcher::new(cm)?;

    let seqs: Vec<&str> = contigs.iter().map(|(_, s)| s.as_str()).collect();
    let cfg = FaithfulConfig {
        e_report: ATTC_EVALUE,
        ..FaithfulConfig::default()
    };
    let hits = searcher.search(&seqs, &cfg);

    let mut out = Vec::new();
    for h in hits {
        if h.evalue > ATTC_EVALUE {
            continue;
        }
        // infernox reports start > stop for reverse-complement hits; normalise.
        let (start, end) = if h.start <= h.stop {
            (h.start, h.stop)
        } else {
            (h.stop, h.start)
        };
        out.push(Locus {
            contig: contigs[h.seq_idx].0.clone(),
            start,
            end,
            strand: if h.in_rc { -1 } else { 1 },
            evalue: Some(h.evalue),
        });
    }
    Ok(out)
}

/// Cluster integrase genes + attC sites into integron calls (pure geometry, no
/// IO — this is the unit-tested core).
///
///  * each intI integrase claims every attC within [`INTEGRASE_ATTC_WINDOW`] bp
///    -> a **complete** integron spanning integrase + its attC array (or **In0**
///    when it claims none);
///  * attC sites left unclaimed are grouped into runs with <=[`ATTC_CLUSTER_GAP`]
///    bp between neighbours; a run of >=[`CALIN_MIN_ATTC`] sites becomes a
///    **CALIN**.
/// Keep the smaller (more significant) of two optional E-values; `None` is absent.
fn min_opt(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.min(y)),
        (Some(x), None) | (None, Some(x)) => Some(x),
        (None, None) => None,
    }
}

fn assemble_integrons(intis: &[Locus], attcs: &[Locus]) -> Vec<IntegronCall> {
    // Sort attC sites by (contig, start) for deterministic clustering.
    let mut attcs: Vec<(usize, &Locus)> = attcs.iter().enumerate().collect();
    attcs.sort_by(|(_, a), (_, b)| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));

    let mut claimed: HashSet<usize> = HashSet::new();
    let mut calls = Vec::new();

    // --- Integrase-anchored: complete / In0. ---
    for inti in intis {
        let mut lo = inti.start;
        let mut hi = inti.end;
        let mut n_attc = 0;
        let mut attc_evalue: Option<f64> = None;
        for (idx, a) in &attcs {
            if claimed.contains(idx) {
                continue;
            }
            if inti.gap(a) <= INTEGRASE_ATTC_WINDOW {
                claimed.insert(*idx);
                lo = lo.min(a.start);
                hi = hi.max(a.end);
                n_attc += 1;
                attc_evalue = min_opt(attc_evalue, a.evalue);
            }
        }
        calls.push(IntegronCall {
            contig: inti.contig.clone(),
            start: lo,
            end: hi,
            strand: inti.strand,
            class: if n_attc > 0 {
                IntegronClass::Complete
            } else {
                IntegronClass::In0
            },
            n_attc,
            attc_evalue,
        });
    }

    // --- Standalone attC runs -> CALIN. ---
    let mut i = 0;
    while i < attcs.len() {
        let (idx, a) = attcs[i];
        if claimed.contains(&idx) {
            i += 1;
            continue;
        }
        // Extend the run over consecutive unclaimed attC within ATTC_CLUSTER_GAP.
        let mut lo = a.start;
        let mut hi = a.end;
        let mut count = 1;
        let mut attc_evalue = a.evalue;
        let mut j = i + 1;
        while j < attcs.len() {
            let (jdx, b) = attcs[j];
            if claimed.contains(&jdx) || b.contig != a.contig {
                break;
            }
            // Gap measured from the current run's right edge to the next site.
            if b.start - hi <= ATTC_CLUSTER_GAP {
                lo = lo.min(b.start);
                hi = hi.max(b.end);
                count += 1;
                attc_evalue = min_opt(attc_evalue, b.evalue);
                j += 1;
            } else {
                break;
            }
        }
        if count >= CALIN_MIN_ATTC {
            calls.push(IntegronCall {
                contig: a.contig.clone(),
                start: lo,
                end: hi,
                strand: 1,
                class: IntegronClass::Calin,
                n_attc: count,
                attc_evalue,
            });
        }
        i = j;
    }

    calls
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    fn locus(contig: &str, start: i64, end: i64) -> Locus {
        Locus {
            contig: contig.into(),
            start,
            end,
            strand: 1,
            evalue: None,
        }
    }

    #[test]
    fn gap_same_and_cross_contig() {
        let a = locus("c1", 100, 200);
        let b = locus("c1", 500, 600);
        assert_eq!(a.gap(&b), 300);
        assert_eq!(b.gap(&a), 300);
        // Overlapping -> 0.
        assert_eq!(a.gap(&locus("c1", 150, 250)), 0);
        // Different contig -> MAX.
        assert_eq!(a.gap(&locus("c2", 100, 200)), i64::MAX);
    }

    #[test]
    fn complete_integron_integrase_plus_attc() {
        let intis = vec![locus("c1", 1000, 2000)];
        // Two attC just downstream, within 4 kb.
        let attcs = vec![locus("c1", 2500, 2550), locus("c1", 3000, 3050)];
        let calls = assemble_integrons(&intis, &attcs);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].class, IntegronClass::Complete);
        assert_eq!(calls[0].n_attc, 2);
        assert_eq!(calls[0].start, 1000);
        assert_eq!(calls[0].end, 3050);
    }

    #[test]
    fn in0_integrase_no_attc() {
        let intis = vec![locus("c1", 1000, 2000)];
        // attC far away (>4 kb) -> not claimed, and only one -> no CALIN either.
        let attcs = vec![locus("c1", 50_000, 50_050)];
        let calls = assemble_integrons(&intis, &attcs);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].class, IntegronClass::In0);
        assert_eq!(calls[0].n_attc, 0);
        assert_eq!((calls[0].start, calls[0].end), (1000, 2000));
    }

    #[test]
    fn calin_cluster_without_integrase() {
        let intis: Vec<Locus> = Vec::new();
        let attcs = vec![
            locus("c1", 5000, 5050),
            locus("c1", 6000, 6050),
            locus("c1", 7000, 7050),
        ];
        let calls = assemble_integrons(&intis, &attcs);
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].class, IntegronClass::Calin);
        assert_eq!(calls[0].n_attc, 3);
        assert_eq!((calls[0].start, calls[0].end), (5000, 7050));
    }

    #[test]
    fn lone_attc_is_not_reported() {
        // A single attC with no integrase is below CALIN_MIN_ATTC -> nothing.
        let calls = assemble_integrons(&[], &[locus("c1", 5000, 5050)]);
        assert!(calls.is_empty());
    }

    #[test]
    fn attc_far_from_integrase_forms_separate_calin() {
        // Integrase with a nearby attC (complete), plus a distant 2-attC cluster
        // (CALIN) on the same contig.
        let intis = vec![locus("c1", 1000, 2000)];
        let attcs = vec![
            locus("c1", 2500, 2550),   // claimed by integrase
            locus("c1", 40_000, 40_050),
            locus("c1", 41_000, 41_050),
        ];
        let calls = assemble_integrons(&intis, &attcs);
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].class, IntegronClass::Complete);
        assert_eq!(calls[0].n_attc, 1);
        assert_eq!(calls[1].class, IntegronClass::Calin);
        assert_eq!(calls[1].n_attc, 2);
    }
}
