//! Conflict resolution over overlapping features.
//!
//! Subunits emit features independently, so the raw set can contain conflicts:
//! an RNA locus hit by several Rfam models (many overlapping ncRNA calls), or a
//! prodigal CDS lying across a tRNA/rRNA/ncRNA. Matching Bakta/Prokka, RNA
//! features mask overlapping CDS, and same-locus RNA calls are deduplicated to
//! the best-scoring one.
//!
//! Deliberately *not* resolved: CDS-vs-CDS overlaps (bacterial genes overlap in
//! operons — prodigal calls them on purpose), and IS elements / CRISPR arrays
//! (an IS element legitimately *contains* a transposase CDS), which pass through
//! untouched.

use crate::feature::{Feature, FeatureKind};

/// Two features conflict when they overlap by at least this fraction of the
/// shorter one. Below it the overlap is treated as incidental (e.g. a CDS that
/// merely abuts a tRNA) and both are kept.
const MIN_OVERLAP_FRAC: f64 = 0.5;

/// Higher wins a conflict. RNA outranks CDS; kinds excluded from resolution
/// (IS elements, CRISPR) are handled separately and never appear here.
fn priority(kind: FeatureKind) -> u8 {
    match kind {
        FeatureKind::Rrna => 4,
        FeatureKind::Trna | FeatureKind::Tmrna => 3,
        FeatureKind::Ncrna | FeatureKind::RegulatoryRegion => 2,
        FeatureKind::Cds => 1,
        // Non-resolvable overlay/structural kinds default to 0 (never reached
        // for these since `is_resolvable` excludes them, but the match is
        // exhaustive).
        FeatureKind::Crispr
        | FeatureKind::IsElement
        | FeatureKind::OriC
        | FeatureKind::AssemblyGap
        | FeatureKind::TandemRepeat
        | FeatureKind::Prophage
        | FeatureKind::Integron
        | FeatureKind::OriT => 0,
    }
}

/// Whether a feature participates in overlap resolution (RNA + CDS). IS elements,
/// CRISPR arrays, and structural/mobile overlays (oriC, gaps, repeats, prophage,
/// integron, oriT) are passed through untouched — they legitimately overlap CDS.
fn is_resolvable(kind: FeatureKind) -> bool {
    matches!(
        kind,
        FeatureKind::Rrna
            | FeatureKind::Trna
            | FeatureKind::Tmrna
            | FeatureKind::Ncrna
            | FeatureKind::Cds
    )
}

/// Best (highest) annotation bit score, or `f32::MIN` if unannotated — used to
/// order features of equal priority so the strongest call is kept.
fn best_score(f: &Feature) -> f32 {
    f.best_annotation().map(|a| a.score).unwrap_or(f32::MIN)
}

/// Overlap of two features as a fraction of the shorter feature's length.
/// `0.0` if they are on different contigs or do not overlap.
fn overlap_frac(a: &Feature, b: &Feature) -> f64 {
    if a.contig != b.contig {
        return 0.0;
    }
    let lo = a.start.max(b.start);
    let hi = a.end.min(b.end);
    if hi < lo {
        return 0.0;
    }
    let ov = (hi - lo + 1) as f64;
    let shorter = a.len_nt().min(b.len_nt()).max(1) as f64;
    ov / shorter
}

/// Resolve overlaps: RNA masks overlapping CDS, same-locus RNA calls collapse to
/// the best-scoring one, CDS-vs-CDS and IS/CRISPR are left alone. The returned
/// features are sorted by contig then start coordinate.
pub fn resolve(features: Vec<Feature>) -> Vec<Feature> {
    let (mut resolvable, passthrough): (Vec<Feature>, Vec<Feature>) =
        features.into_iter().partition(|f| is_resolvable(f.kind));

    // Strongest candidates first: higher priority, then higher score. A conflict
    // is always decided in favour of an already-kept feature, so ordering here is
    // what makes the winner deterministic.
    resolvable.sort_by(|a, b| {
        priority(b.kind)
            .cmp(&priority(a.kind))
            .then(best_score(b).partial_cmp(&best_score(a)).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut kept: Vec<Feature> = Vec::new();
    for cand in resolvable {
        let mut drop = false;
        for k in &kept {
            if overlap_frac(k, &cand) < MIN_OVERLAP_FRAC {
                continue;
            }
            // A strictly higher-priority feature occupies this locus (e.g. tRNA
            // over CDS), or an equal-priority RNA already claimed it (dedup of
            // multi-model ncRNA hits). CDS never dedups against CDS.
            if priority(k.kind) > priority(cand.kind)
                || (priority(k.kind) == priority(cand.kind) && cand.kind != FeatureKind::Cds)
            {
                drop = true;
                break;
            }
        }
        if !drop {
            kept.push(cand);
        }
    }

    kept.extend(passthrough);
    kept.sort_by(|a, b| {
        a.contig
            .cmp(&b.contig)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Annotation, Functional};

    fn feat(kind: FeatureKind, start: i64, end: i64, score: f32) -> Feature {
        Feature {
            kind,
            contig: "c1".into(),
            id: format!("{kind:?}_{start}_{end}"),
            start,
            end,
            strand: 1,
            aa: if kind == FeatureKind::Cds {
                Some("M".into())
            } else {
                None
            },
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: "t".into(),
                accession: "-".into(),
                name: "n".into(),
                score,
                evalue: Some(1.0),
                ref_len: None,
            }],
            func: Functional::default(),
        }
    }

    #[test]
    fn rna_masks_overlapping_cds() {
        // A tRNA fully inside a CDS: the CDS is dropped.
        let out = resolve(vec![
            feat(FeatureKind::Cds, 100, 400, 50.0),
            feat(FeatureKind::Trna, 150, 220, 80.0),
        ]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].kind, FeatureKind::Trna);
    }

    #[test]
    fn overlapping_cds_are_both_kept() {
        // Bacterial operon overlap: CDS never masks another CDS.
        let out = resolve(vec![
            feat(FeatureKind::Cds, 100, 400, 50.0),
            feat(FeatureKind::Cds, 380, 700, 60.0),
        ]);
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn multi_model_ncrna_dedups_to_best() {
        // Same locus hit by two ncRNA models: keep the higher-scoring one.
        let out = resolve(vec![
            feat(FeatureKind::Ncrna, 100, 200, 30.0),
            feat(FeatureKind::Ncrna, 102, 205, 90.0),
        ]);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].best_annotation().unwrap().score, 90.0);
    }

    #[test]
    fn is_element_passes_through_over_cds() {
        // An IS element legitimately contains a transposase CDS — both kept.
        let out = resolve(vec![
            feat(FeatureKind::IsElement, 100, 1300, 10.0),
            feat(FeatureKind::Cds, 200, 1100, 40.0),
        ]);
        assert_eq!(out.len(), 2);
    }
}
