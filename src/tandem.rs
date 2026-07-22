//! Tandem repeat detection — thin adapter over the standalone `tandemox` crate.
//!
//! The tandem-repeat scanner used to live in-tree here; it has been extracted to
//! the published `tandemox` crate (a TRF-comparable, pure-Rust, license-clean
//! reimplementation validated against C TRF). This module now just drives
//! `tandemox` per contig and maps each [`tandemox::TandemRepeat`] onto a
//! `Feature { kind: FeatureKind::TandemRepeat }` (→ GFF3 `repeat_region`),
//! preserving the previous feature shape (`<contig>_trf<n>` ids, the
//! `tandem repeat period=<p> copies=<c>` annotation name, and the TRF-comparable
//! score) so downstream output/resolve behaviour is unchanged.
//!
//! `tandemox::TandemRepeat` coordinates are already 1-based inclusive, matching
//! the bactars `Feature` convention, so no offset adjustment is applied.

use crate::fasta::Contig;
use crate::feature::{Annotation, Feature, FeatureKind, Functional};

/// Detect approximate tandem repeats in every contig via `tandemox` and return
/// them as `FeatureKind::TandemRepeat` features. DB-free; no external process.
pub fn detect(contigs: &[Contig]) -> Vec<Feature> {
    let mut features = Vec::new();

    for contig in contigs {
        // tandemox owns the period scan / maximal-segment extraction / harmonic
        // dedup. Deterministic, left-to-right order per contig.
        for (i, r) in tandemox::detect(&contig.seq).into_iter().enumerate() {
            let idx1 = i + 1;
            let name = format!("tandem repeat period={} copies={:.1}", r.period, r.copies);
            features.push(Feature {
                kind: FeatureKind::TandemRepeat,
                contig: contig.name.clone(),
                id: format!("{}_trf{}", contig.name, idx1),
                start: r.start as i64, // tandemox is already 1-based inclusive
                end: r.end as i64,
                strand: 1,
                aa: None,
                partial5: false,
                partial3: false,
                annotations: vec![Annotation {
                    source: "bactars:tandem".into(),
                    accession: String::new(),
                    name,
                    score: r.score,
                    // Score-only: tandemox reports a TRF alignment score, no E-value.
                    evalue: None,
                    ref_len: None,
                }],
                func: Functional::default(),
            });
        }
    }

    features
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(name: &str, seq: &[u8]) -> Contig {
        Contig {
            name: name.to_string(),
            seq: seq.to_vec(),
        }
    }

    /// Parse the `period=<n>` field out of a tandem-repeat annotation name.
    fn period_of(f: &Feature) -> usize {
        let ann = &f.annotations[0];
        ann.name
            .split_whitespace()
            .find_map(|tok| tok.strip_prefix("period="))
            .and_then(|v| v.parse::<usize>().ok())
            .expect("annotation must carry period=")
    }

    // Non-repetitive filler used to flank repeats.
    const FL_L: &[u8] = b"GACTTCAGGATCCGATACGTTAGCACTG";
    const FL_R: &[u8] = b"CATGACTGACAGTGCACTAGCTAGCTAA";

    /// The adapter surfaces a clean embedded period-4 repeat as a TandemRepeat
    /// feature with the expected id/annotation shape and 1-based coordinates.
    #[test]
    fn embedded_period4_repeat_is_reported() {
        let mut seq = Vec::new();
        seq.extend_from_slice(FL_L);
        let rep_start_0 = seq.len(); // 0-based start of the repeat
        for _ in 0..10 {
            seq.extend_from_slice(b"ATGC"); // 40 bp, period 4
        }
        seq.extend_from_slice(FL_R);

        let feats = detect(&[contig("c1", &seq)]);
        assert!(!feats.is_empty(), "expected a tandem repeat feature");

        let f = feats
            .iter()
            .find(|f| period_of(f) == 4)
            .expect("a period-4 feature");
        // Feature shape / provenance preserved by the adapter.
        assert_eq!(f.kind, FeatureKind::TandemRepeat);
        assert_eq!(f.annotations[0].source, "bactars:tandem");
        assert!(f.id.starts_with("c1_trf"));
        // 1-based inclusive coords that cover the embedded ATGC array.
        assert!(f.start >= 1);
        assert!(f.start as usize <= rep_start_0 + 2);
        assert!(f.end as usize >= rep_start_0 + 30);
    }

    /// Random / non-repetitive sequence yields no tandem features.
    #[test]
    fn non_repetitive_sequence_is_empty() {
        let seq = [FL_L, FL_R, FL_L].concat();
        let feats = detect(&[contig("c1", &seq)]);
        assert!(
            feats.iter().all(|f| period_of(f) >= 1),
            "any reported feature must carry a period"
        );
        // No long tandem structure here; the flanks were hand-picked to be clean.
        assert!(feats.is_empty(), "expected no tandem repeats, got {}", feats.len());
    }
}
