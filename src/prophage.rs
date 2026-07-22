//! Prophage detection by phage-hallmark gene clustering, via rustyhmmer.
//!
//! A prophage is an integrated (lysogenic) bacteriophage genome. It is recognised
//! not by any single gene but by a genomic REGION unusually dense in phage
//! "hallmark" functions — the virion morphogenesis + DNA-packaging machinery
//! (terminase, portal, major capsid, tail), plus integration and lysis genes —
//! that a bacterial genome does not otherwise carry clustered together. This is
//! the shared logic of PHASTER / VIBRANT / Prophage Hunter, reduced to a small,
//! license-clean Pfam hallmark set (`db/phage/phage_hallmark.hmm`, 25 Pfam
//! profiles, CC0).
//!
//! Method:
//!   1. Scan all CDS proteins with the hallmark HMMs (rustyhmmer, `--cut_ga`);
//!      each hit's model maps to a hallmark CLASS (terminase / portal / capsid /
//!      protease / tail / integrase / lysis).
//!   2. Cluster hallmark-hit CDS along each contig: consecutive hallmark genes
//!      within [`MAX_GENE_GAP`] bp join one region.
//!   3. Call a prophage when a region holds >=[`MIN_HALLMARK_CLASSES`] DISTINCT
//!      hallmark classes AND at least one CORE packaging/structural class
//!      (terminase / portal / capsid) — the strongest, most phage-exclusive
//!      signal. Favour precision over flooding.
//!
//! Emits `Feature { kind: FeatureKind::Prophage }` spanning the region, with the
//! hallmark inventory and a confidence tier ("intact" vs "questionable") in the
//! annotation.

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use crate::hmm::{self, Cutoff};
use std::collections::BTreeSet;
use std::path::Path;

/// Max distance (bp) between consecutive hallmark genes for them to belong to the
/// same prophage region. Prophages are ~10-50 kb and gene-dense; a gap larger
/// than this breaks the region (avoids chaining unrelated scattered hits).
const MAX_GENE_GAP: i64 = 10_000;

/// Minimum number of DISTINCT hallmark classes in a region to call a prophage.
/// Three independent phage functions clustered together (with at least one core
/// packaging class — see [`Hallmark::is_core`]) is a confident, precision-first
/// signal: on MG1655 it recovers the DLP12 and e14 cryptic prophages and nothing
/// spurious. Cryptic prophages too degraded to retain a core packaging gene
/// (Rac/Qin/CP4-57) are deliberately missed rather than risk false positives.
const MIN_HALLMARK_CLASSES: usize = 3;

/// Class count at/above which a region is tiered "intact" rather than
/// "questionable" (analogous to PHASTER's completeness tiers).
const INTACT_CLASSES: usize = 6;

/// A phage hallmark functional class.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Hallmark {
    Terminase,
    Portal,
    Capsid,
    Protease,
    Tail,
    Integrase,
    Lysis,
}

impl Hallmark {
    fn label(self) -> &'static str {
        match self {
            Hallmark::Terminase => "terminase",
            Hallmark::Portal => "portal",
            Hallmark::Capsid => "capsid",
            Hallmark::Protease => "protease",
            Hallmark::Tail => "tail",
            Hallmark::Integrase => "integrase",
            Hallmark::Lysis => "lysis",
        }
    }
    /// The DNA-packaging / capsid core: the most phage-exclusive signal. A region
    /// needs at least one of these to be called (a lone integrase+lysis pair is
    /// too generic — mobile elements and defence systems carry those too).
    fn is_core(self) -> bool {
        matches!(self, Hallmark::Terminase | Hallmark::Portal | Hallmark::Capsid)
    }
}

/// Map a hallmark HMM hit to its functional class.
///
/// Classification keys on the Pfam **ACCESSION** (the `ACC` line, e.g. `PF03354`)
/// rather than the model NAME. The accession is Pfam's stable identifier — a Pfam
/// release can rename a family (NAME drift) without changing its accession, so a
/// NAME-only match would silently stop classifying (yielding zero prophages) after
/// such a rename. Accessions carry a `.NN` version suffix (`PF03354.22`); we match
/// the bare accession so version bumps don't break the map.
///
/// The model NAME is retained as a compatibility fallback for profiles that carry
/// no `ACC` line (custom/non-Pfam hallmark HMMs). Accessions/names correspond to
/// the 25 Pfam profiles shipped in `db/phage/phage_hallmark.hmm`.
fn classify(model_name: &str, model_acc: &str) -> Option<Hallmark> {
    // Primary, stable key: bare Pfam accession.
    let bare_acc = model_acc.split('.').next().unwrap_or(model_acc);
    let by_acc = match bare_acc {
        // terminase large + small subunit
        "PF03237" | "PF03354" | "PF03592" | "PF04466" | "PF05876" => Some(Hallmark::Terminase),
        // portal
        "PF04860" | "PF05136" => Some(Hallmark::Portal),
        // major capsid / head
        "PF03864" | "PF04233" | "PF05065" => Some(Hallmark::Capsid),
        // prohead protease
        "PF04586" => Some(Hallmark::Protease),
        // tail (sheath / baseplate / tube / fibre / minor)
        "PF03906" | "PF04865" | "PF04984" | "PF04985" | "PF05100" | "PF06199" | "PF10145"
        | "PF10618" => Some(Hallmark::Tail),
        // integration
        "PF00589" | "PF07508" => Some(Hallmark::Integrase),
        // lysis
        "PF00959" | "PF05105" | "PF11860" | "PF16080" => Some(Hallmark::Lysis),
        _ => None,
    };
    if by_acc.is_some() {
        return by_acc;
    }
    // Fallback: model NAME (for profiles with no ACC line).
    match model_name {
        "TerL_ATPase" | "Terminase_3" | "Terminase_6N" | "GpA_ATPase" | "Terminase_2" => {
            Some(Hallmark::Terminase)
        }
        "Phage_portal" | "Phage_portal_2" => Some(Hallmark::Portal),
        "Phage_capsid" | "Phage_cap_E" | "Phage_Mu_F" => Some(Hallmark::Capsid),
        "Peptidase_S78" => Some(Hallmark::Protease),
        "Phage_sheath_1" | "Baseplate_J" | "Tail_tube" | "Phage_tail_2" | "Phage_tube"
        | "Phage_T7_tail" | "Phage_tail_L" | "PhageMin_Tail" => Some(Hallmark::Tail),
        "Phage_integrase" | "Recombinase" => Some(Hallmark::Integrase),
        "Phage_lysozyme" | "Muramidase" | "Phage_holin_4_1" | "Phage_holin_2_3" => {
            Some(Hallmark::Lysis)
        }
        _ => None,
    }
}

/// One CDS carrying a phage hallmark hit (genomic locus + class).
#[derive(Clone, Debug)]
struct HallmarkGene {
    contig: String,
    start: i64,
    end: i64,
    class: Hallmark,
    /// Real rustyhmmer sequence E-value of this CDS's hallmark HMM hit.
    evalue: Option<f64>,
}

/// A called prophage region.
#[derive(Clone, Debug)]
struct ProphageRegion {
    contig: String,
    start: i64,
    end: i64,
    classes: BTreeSet<Hallmark>,
    n_genes: usize,
    /// Best (smallest) hallmark-HMM E-value among the region's genes; `None` only
    /// if no gene carried one.
    evalue: Option<f64>,
}

/// Detect prophage regions across the genome. `features` supplies the CDS calls
/// (their proteins are scanned + their coordinates locate the region); `contigs`
/// is accepted for signature uniformity (region bounds come from the genes).
/// `phage_dir` holds `phage_hallmark.hmm`.
pub fn detect(
    _contigs: &[(String, String)],
    features: &[Feature],
    phage_dir: &str,
) -> Result<Vec<Feature>, String> {
    let hmm_path = format!("{}/phage_hallmark.hmm", phage_dir.trim_end_matches('/'));
    if !Path::new(&hmm_path).exists() {
        // The hallmark HMM set is the required model DB; a missing file is a
        // configuration failure surfaced to the caller, not a silent empty result.
        return Err(format!("no phage_hallmark.hmm in {phage_dir}"));
    }

    // CDS proteins (id, aa) + id -> feature index for coordinate lookup.
    let proteins: Vec<(String, String)> = features
        .iter()
        .filter(|f| f.kind == FeatureKind::Cds)
        .filter_map(|f| f.aa.as_ref().map(|aa| (f.id.clone(), aa.clone())))
        .collect();
    if proteins.is_empty() {
        return Ok(Vec::new());
    }

    // Hallmark HMMs carry curated GA cutoffs (Pfam) — search with `--cut_ga`.
    let hits = hmm::annotate(&proteins, &hmm_path, Cutoff::GatheringGa)
        .map_err(|e| format!("hallmark HMM search failed ({e})"))?;

    // Best (highest-scoring) hallmark class per CDS id.
    let n_hits = hits.len();
    let mut classified = 0usize;
    let mut best: std::collections::HashMap<String, (f32, Hallmark, f64)> =
        std::collections::HashMap::new();
    for h in hits {
        if let Some(class) = classify(&h.query_name, &h.query_acc) {
            classified += 1;
            let e = best.entry(h.target_name).or_insert((f32::MIN, class, f64::NAN));
            if h.seq_score > e.0 {
                *e = (h.seq_score, class, h.seq_evalue);
            }
        }
    }
    // Diagnose a silent-zero: hallmark HMMs DID hit CDS, but none mapped to a known
    // class. That almost always means the classification map has drifted out of
    // sync with the shipped `phage_hallmark.hmm` (renamed NAMEs / new accessions) —
    // make it loud rather than returning zero prophages without explanation.
    if n_hits > 0 && classified == 0 {
        eprintln!(
            "prophage: {n_hits} hallmark HMM hit(s) but NONE classified into a hallmark \
             class — the classify() accession/name map is likely out of sync with \
             {hmm_path}; reporting zero prophages"
        );
    }

    // Map each hallmark-hit CDS to a genomic locus.
    let mut genes: Vec<HallmarkGene> = Vec::new();
    for f in features.iter().filter(|f| f.kind == FeatureKind::Cds) {
        if let Some((_, class, seq_evalue)) = best.get(&f.id) {
            let (start, end) = if f.start <= f.end {
                (f.start, f.end)
            } else {
                (f.end, f.start)
            };
            genes.push(HallmarkGene {
                contig: f.contig.clone(),
                start,
                end,
                class: *class,
                evalue: (!seq_evalue.is_nan()).then_some(*seq_evalue),
            });
        }
    }

    let regions = cluster_regions(&mut genes);

    // Emit features (per-contig numbering; regions already contig/start sorted).
    let mut out = Vec::new();
    let mut per_contig: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    for r in regions {
        let n = per_contig.entry(r.contig.clone()).or_insert(0);
        *n += 1;
        let tier = if r.classes.len() >= INTACT_CLASSES {
            "intact"
        } else {
            "questionable"
        };
        let inventory: Vec<&str> = r.classes.iter().map(|c| c.label()).collect();
        out.push(Feature {
            kind: FeatureKind::Prophage,
            contig: r.contig.clone(),
            id: format!("{}_prophage{}", r.contig, n),
            start: r.start,
            end: r.end,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: format!("prophage:{tier}"),
                accession: String::new(),
                name: format!(
                    "prophage region [{}] ({} hallmark genes, {} classes: {})",
                    tier,
                    r.n_genes,
                    r.classes.len(),
                    inventory.join(",")
                ),
                score: r.classes.len() as f32,
                // Best hallmark-HMM sequence E-value across the region's genes.
                evalue: r.evalue,
                ref_len: None,
            }],
            func: Functional::default(),
        });
    }
    Ok(out)
}

/// Cluster hallmark genes into prophage regions (pure geometry — unit-tested).
/// Genes are grouped per contig into runs where consecutive genes are within
/// [`MAX_GENE_GAP`] bp; a run is kept only if it holds >=[`MIN_HALLMARK_CLASSES`]
/// distinct classes and >=1 core (terminase/portal/capsid) class. Output is in
/// (contig, start) order.
/// Keep the smaller (more significant) of two optional E-values; `None` is absent.
fn min_opt(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(x.min(y)),
        (Some(x), None) | (None, Some(x)) => Some(x),
        (None, None) => None,
    }
}

fn cluster_regions(genes: &mut [HallmarkGene]) -> Vec<ProphageRegion> {
    genes.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));

    let mut regions = Vec::new();
    let mut i = 0;
    while i < genes.len() {
        let g0 = &genes[i];
        let mut lo = g0.start;
        let mut hi = g0.end;
        let mut classes: BTreeSet<Hallmark> = BTreeSet::new();
        classes.insert(g0.class);
        let mut n_genes = 1;
        let mut evalue = g0.evalue;
        let mut j = i + 1;
        while j < genes.len() {
            let g = &genes[j];
            if g.contig != g0.contig || g.start - hi > MAX_GENE_GAP {
                break;
            }
            lo = lo.min(g.start);
            hi = hi.max(g.end);
            classes.insert(g.class);
            n_genes += 1;
            evalue = min_opt(evalue, g.evalue);
            j += 1;
        }
        // Keep the region only if it clears the density + core-signal bars.
        let has_core = classes.iter().any(|c| c.is_core());
        if classes.len() >= MIN_HALLMARK_CLASSES && has_core {
            regions.push(ProphageRegion {
                contig: g0.contig.clone(),
                start: lo,
                end: hi,
                classes,
                n_genes,
                evalue,
            });
        }
        i = j;
    }
    regions
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    fn gene(contig: &str, start: i64, class: Hallmark) -> HallmarkGene {
        HallmarkGene {
            contig: contig.into(),
            start,
            end: start + 900,
            class,
            evalue: None,
        }
    }

    #[test]
    fn classify_by_accession() {
        // Primary key: bare Pfam accession (from the shipped phage_hallmark.hmm).
        assert_eq!(classify("", "PF03354"), Some(Hallmark::Terminase)); // TerL_ATPase
        assert_eq!(classify("", "PF04860"), Some(Hallmark::Portal)); // Phage_portal
        assert_eq!(classify("", "PF05065"), Some(Hallmark::Capsid)); // Phage_capsid
        assert_eq!(classify("", "PF00589"), Some(Hallmark::Integrase)); // Phage_integrase
        assert_eq!(classify("", "PF00959"), Some(Hallmark::Lysis)); // Phage_lysozyme
        assert_eq!(classify("", "PF06199"), Some(Hallmark::Tail)); // Phage_tail_2
        assert_eq!(classify("", "PF04586"), Some(Hallmark::Protease)); // Peptidase_S78
        // Versioned accession still maps (bare-accession match).
        assert_eq!(classify("", "PF03354.22"), Some(Hallmark::Terminase));
        // Unknown accession/name -> None.
        assert_eq!(classify("not_a_phage_model", "PF99999"), None);
    }

    #[test]
    fn classify_survives_name_drift_via_accession() {
        // If a Pfam release RENAMES the family (NAME drift) but keeps the stable
        // accession, classification must still succeed off the accession.
        assert_eq!(
            classify("Terminase_large_ATPase_renamed", "PF03354.23"),
            Some(Hallmark::Terminase)
        );
    }

    #[test]
    fn classify_name_fallback_when_no_accession() {
        // Custom/non-Pfam profiles with no ACC line still classify by NAME.
        assert_eq!(classify("TerL_ATPase", ""), Some(Hallmark::Terminase));
        assert_eq!(classify("Phage_portal", ""), Some(Hallmark::Portal));
        assert_eq!(classify("Phage_lysozyme", ""), Some(Hallmark::Lysis));
        assert_eq!(classify("not_a_phage_model", ""), None);
    }

    #[test]
    fn dense_region_with_core_is_called_intact() {
        // Six distinct classes incl. core -> intact prophage.
        let mut genes = vec![
            gene("c1", 10_000, Hallmark::Integrase),
            gene("c1", 12_000, Hallmark::Terminase),
            gene("c1", 14_000, Hallmark::Portal),
            gene("c1", 16_000, Hallmark::Capsid),
            gene("c1", 18_000, Hallmark::Tail),
            gene("c1", 20_000, Hallmark::Lysis),
        ];
        let regions = cluster_regions(&mut genes);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].classes.len(), 6);
        assert!(regions[0].classes.len() >= INTACT_CLASSES);
        assert_eq!(regions[0].start, 10_000);
        assert_eq!(regions[0].end, 20_900);
    }

    #[test]
    fn region_without_core_is_rejected() {
        // Four classes but NO core (terminase/portal/capsid) -> not a prophage.
        let mut genes = vec![
            gene("c1", 10_000, Hallmark::Integrase),
            gene("c1", 11_000, Hallmark::Lysis),
            gene("c1", 12_000, Hallmark::Tail),
            gene("c1", 13_000, Hallmark::Protease),
        ];
        let regions = cluster_regions(&mut genes);
        assert!(regions.is_empty(), "no core packaging class -> reject");
    }

    #[test]
    fn too_few_classes_is_rejected() {
        // Core present but only 2 distinct classes -> below threshold.
        let mut genes = vec![
            gene("c1", 10_000, Hallmark::Terminase),
            gene("c1", 11_000, Hallmark::Capsid),
        ];
        let regions = cluster_regions(&mut genes);
        assert!(regions.is_empty());
    }

    #[test]
    fn distant_genes_do_not_chain() {
        // A qualifying cluster, then a far-away single gene (>10 kb gap) that must
        // not join it and is not itself a region.
        let mut genes = vec![
            gene("c1", 10_000, Hallmark::Terminase),
            gene("c1", 11_000, Hallmark::Portal),
            gene("c1", 12_000, Hallmark::Capsid),
            gene("c1", 13_000, Hallmark::Tail),
            gene("c1", 500_000, Hallmark::Capsid),
        ];
        let regions = cluster_regions(&mut genes);
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].n_genes, 4);
        assert!(regions[0].end < 500_000);
    }

    #[test]
    fn two_regions_on_two_contigs() {
        let mut genes = vec![
            gene("c1", 10_000, Hallmark::Terminase),
            gene("c1", 11_000, Hallmark::Portal),
            gene("c1", 12_000, Hallmark::Capsid),
            gene("c1", 13_000, Hallmark::Integrase),
            gene("c2", 5_000, Hallmark::Terminase),
            gene("c2", 6_000, Hallmark::Portal),
            gene("c2", 7_000, Hallmark::Capsid),
            gene("c2", 8_000, Hallmark::Tail),
        ];
        let regions = cluster_regions(&mut genes);
        assert_eq!(regions.len(), 2);
        assert_eq!(regions[0].contig, "c1");
        assert_eq!(regions[1].contig, "c2");
    }
}
