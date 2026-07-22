//! Gene-model quality scoring via the learned `balrox` model (non-destructive).
//!
//! `balrox` is a pure-Rust dilated temporal CNN (architecture inspired by Balrog, own
//! weights on RefSeq CDS vs real intergenic/shadow/off-frame ORFs) that returns
//! *P(real protein-coding gene)* for a CDS's amino-acid sequence. bactars uses it as an
//! **additive quality signal**: every CDS gets a `gene_score`, and CDS whose score falls
//! below [`LOW_SCORE`] get a `/note` marking them possible spurious ORFs. **No CDS is
//! ever removed** — the faithful gene caller's calls are preserved; the score is advisory.
//!
//! On a clean complete genome this changes little (Prodigal rarely emits junk); its value
//! is on fragmented / draft / metagenomic assemblies, where contig-end and shadow ORFs
//! are common and score low. A separate hard-filter mode is intentionally not the default.

use crate::feature::{Feature, FeatureKind};

/// Below this balrox score a CDS is flagged (via `/note`) as a possible spurious ORF.
/// Calibrated on a fragmented-genome test (spurious ORFs scored a median ~0.7 there,
/// well-supported genes ~0.96); a conservative low cutoff flags only the clear tail.
const LOW_SCORE: f32 = 0.1;

/// Annotate a learned gene-model quality score on every CDS (non-destructive).
pub fn annotate(features: &mut [Feature]) {
    let model = balrox::Model::embedded();
    // Collect CDS amino-acid sequences, score them in one parallel batch (balrox's
    // score_many uses rayon internally — same speedup available to any caller), then
    // write results back. Numerically identical to per-sequence scoring.
    let idx: Vec<usize> = features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.kind == FeatureKind::Cds && f.aa.is_some())
        .map(|(i, _)| i)
        .collect();
    let seqs: Vec<&[u8]> = idx
        .iter()
        .map(|&i| features[i].aa.as_ref().unwrap().as_bytes())
        .collect();
    let scores = model.score_many(&seqs);
    for (&i, &s) in idx.iter().zip(scores.iter()) {
        features[i].func.gene_score = Some(s);
        if s < LOW_SCORE {
            features[i].func.note.push(format!(
                "low gene-model score, possible spurious ORF (balrox p={s:.2})"
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::Functional;

    fn cds(aa: &str) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c".into(),
            id: "c_1".into(),
            start: 1,
            end: (aa.len() as i64) * 3,
            strand: 1,
            aa: Some(aa.to_string()),
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        }
    }

    #[test]
    fn real_gene_scores_high_and_is_kept() {
        let mut feats = vec![cds(
            "MSKIVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKG",
        )];
        annotate(&mut feats);
        let s = feats[0].func.gene_score.expect("expected a gene_score");
        assert!(s > 0.5, "real gene should score high, got {s}");
        // never removed
        assert_eq!(feats.len(), 1);
    }

    #[test]
    fn non_cds_untouched() {
        let mut f = cds("MSKIVKIIGREIIDSRGNPTVEAEVHLEGG");
        f.kind = FeatureKind::Trna;
        let mut feats = vec![f];
        annotate(&mut feats);
        assert!(feats[0].func.gene_score.is_none());
    }
}
