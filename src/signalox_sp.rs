//! Signal-peptide / cleavage / transmembrane annotation via the learned `signalox`
//! model (CNN + CRF), replacing the rule-based [`crate::signalpep`] detector.
//!
//! `signalox` is a pure-Rust bacterial signal-peptide predictor (architecture inspired
//! by DeepSig; own weights trained on UniProt/SwissProt). On held-out organism-disjoint
//! data it reaches SP-any precision ~0.95 / recall ~0.90 — far above the rule-based
//! heuristic (Sec precision 45.6%). It also yields the cleavage site (CRF) and predicted
//! transmembrane segments in one pass.
//!
//! Sets the same `feature.func.signal_peptide` attribute + `/note`s the rule-based path
//! did, so downstream output is unchanged in shape.

use crate::feature::{Feature, FeatureKind, SignalPeptide};

/// Annotate signal peptides + transmembrane topology on CDS features in place.
pub fn annotate(features: &mut [Feature]) {
    let model = signalox::Model::embedded();
    // Predict all CDS in one parallel batch (predict_many uses rayon internally),
    // then write results back. Numerically identical to per-sequence prediction.
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
        let f = &mut features[i];

        // Reusable transmembrane-protein annotation (learned TM head).
        if !p.tm_segments.is_empty() {
            let n = p.tm_segments.len();
            f.func.note.push(format!(
                "predicted transmembrane protein ({} TM {})",
                n,
                if n == 1 { "helix" } else { "helices" }
            ));
        }

        // Signal-peptide call. The N-terminal signal is only meaningful when the CDS
        // has a real start; skip 5'-partial genes for the SP call (TM note above still
        // applies as TM segments can be internal).
        if f.partial5 {
            continue;
        }
        let kind = match p.sp_type {
            signalox::SpType::Sec => "sec",
            signalox::SpType::Lipo => "lipoprotein",
            signalox::SpType::Tat => "tat",
            _ => continue, // NoSp / Tm -> no signal peptide
        };
        f.func.signal_peptide = Some(SignalPeptide {
            kind: kind.to_string(),
            cleavage: p.cleavage,
        });
        f.func.note.push(format!(
            "predicted {kind} signal peptide (signalox CNN+CRF, p={:.2})",
            p.sp_prob
        ));
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
    fn sec_signal_detected() {
        // classic Sec signal peptide N-terminus (OmpA-like)
        let mut feats = vec![cds(
            "MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGY",
        )];
        annotate(&mut feats);
        let sp = feats[0].func.signal_peptide.as_ref().expect("expected an SP call");
        assert_eq!(sp.kind, "sec");
    }

    #[test]
    fn cytoplasmic_not_called() {
        // an enolase-like cytoplasmic N-terminus should not get an SP call
        let mut feats = vec![cds(
            "MSKIVKIIGREIIDSRGNPTVEAEVHLEGGFVGMAAAPSGASTGSREALELRDGDKSRFLGKGVTKAVAAVNGPI",
        )];
        annotate(&mut feats);
        assert!(feats[0].func.signal_peptide.is_none());
    }
}
