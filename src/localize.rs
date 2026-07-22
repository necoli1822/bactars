//! Subcellular localization annotation via the learned `lotrs` model.
//!
//! `lotrs` is a pure-Rust bacterial subcellular-localization predictor (a sequence CNN;
//! architecture inspired by DeepLoc, own weights trained on UniProt/SwissProt). It reads
//! each CDS's protein sequence and predicts one of {cytoplasm, inner_membrane, periplasm,
//! outer_membrane, extracellular, cell_wall}. A confidence gate keeps only calls above
//! [`MIN_PROB`] — low-confidence predictions leave `localization` unset rather than
//! guessing. This is a new annotation bactars did not previously produce.

use crate::feature::{Feature, FeatureKind};

/// Minimum softmax probability for a localization call to be recorded. The model is
/// strongest on the well-populated classes; this gate suppresses the noisy tail.
const MIN_PROB: f32 = 0.7;

/// Annotate subcellular localization on CDS features in place.
pub fn annotate(features: &mut [Feature]) {
    let model = lotrs::Model::embedded();
    // Predict all CDS in one parallel batch (predict_many uses rayon internally).
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
    let preds = model.predict_many(&seqs);
    for (&i, p) in idx.iter().zip(preds.iter()) {
        if p.prob >= MIN_PROB {
            features[i].func.localization = Some(p.label.clone());
            features[i].func.note.push(format!(
                "predicted subcellular localization: {} (lotrs, p={:.2})",
                p.label, p.prob
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
    fn cytoplasmic_enzyme_localizes() {
        // an enolase-like cytoplasmic sequence should get a confident cytoplasm call
        let mut feats = vec![cds(
            "MSKIVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKGVTKAVAAVNGPI",
        )];
        annotate(&mut feats);
        assert_eq!(feats[0].func.localization.as_deref(), Some("cytoplasm"));
    }

    #[test]
    fn non_cds_untouched() {
        let mut f = cds("MSKIVKIIGREIIDSRGNPTVEAEVHLEGG");
        f.kind = FeatureKind::Trna;
        let mut feats = vec![f];
        annotate(&mut feats);
        assert!(feats[0].func.localization.is_none());
    }
}
