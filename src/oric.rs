//! oriC (replication origin) detection — thin adapter over the standalone
//! `oricle` crate.
//!
//! The GC-skew + DnaA-box oriC scanner used to live in-tree here; it has been
//! extracted to the published `oricle` crate (a pure-Rust, DB-free
//! reimplementation validated for parity). This module now drives `oricle` per
//! contig and maps each [`oricle::OriC`] onto a `Feature { kind: OriC }`
//! (→ GFF3 `origin_of_replication`), preserving the previous feature shape
//! (`<contig>_oriC` id, `replication origin` annotation name, the confidence
//! score) so downstream output/resolve behaviour is unchanged.
//!
//! `oricle::detect(seq, genes)` takes gene hints so a `dnaA` call can refine or
//! rescue the origin; we pass every already-named CDS on the contig as a
//! [`oricle::GeneHint`] (oricle matches `dnaA` by the name containing "dnaa",
//! which also catches the "…replication initiator protein DnaA" product). Both
//! `oricle` coordinates and the bactars `Feature` are 1-based inclusive, so no
//! offset adjustment is applied.

use crate::fasta::Contig;
use crate::feature::{Annotation, Feature, FeatureKind, Functional};

/// Detect oriC in every contig via `oricle` and return them as
/// `FeatureKind::OriC` features. DB-free; no external process. `cds` supplies
/// gene calls (already named) so oricle can use `dnaA` proximity.
pub fn detect(contigs: &[Contig], cds: &[Feature]) -> Vec<Feature> {
    let mut out = Vec::new();
    for contig in contigs {
        // Gene hints for this contig: oricle only needs coordinates + a name it
        // can test for "dnaa"; pass every CDS so the dnaA refinement can fire.
        let genes: Vec<oricle::GeneHint> = cds
            .iter()
            .filter(|f| f.contig == contig.name && f.kind == FeatureKind::Cds)
            .map(|f| oricle::GeneHint {
                start: f.start.max(1) as usize,
                end: f.end.max(f.start).max(1) as usize,
                name: gene_name(f),
            })
            .collect();

        // A contig usually yields a single oriC (`<contig>_oriC`), but oricle can
        // return more than one candidate trough; disambiguate their ids with a
        // 1-based index only in that case, so the single-oriC id stays unchanged.
        let ors = oricle::detect(&contig.seq, &genes);
        let multi = ors.len() > 1;
        for (i, o) in ors.into_iter().enumerate() {
            let id = if multi {
                format!("{}_oriC_{}", contig.name, i + 1)
            } else {
                format!("{}_oriC", contig.name)
            };
            out.push(Feature {
                kind: FeatureKind::OriC,
                contig: contig.name.clone(),
                id,
                start: o.start as i64, // oricle is already 1-based inclusive
                end: o.end as i64,
                strand: 1,
                aa: None,
                partial5: false,
                partial3: false,
                annotations: vec![Annotation {
                    source: "bactars:oric".into(),
                    accession: String::new(),
                    name: "replication origin".into(),
                    score: o.score,
                    // Score-only: oricle reports a GC-skew/DnaA-box score, no E-value.
                    evalue: None,
                    ref_len: None,
                }],
                func: Functional::default(),
            });
        }
    }
    out
}

/// The best available human name for a CDS, so oricle can spot `dnaA`. Prefer
/// the gene symbol, then the product, then the best annotation name.
fn gene_name(f: &Feature) -> String {
    if let Some(g) = &f.func.gene {
        if !g.is_empty() {
            return g.clone();
        }
    }
    if let Some(p) = &f.func.product {
        if !p.is_empty() {
            return p.clone();
        }
    }
    f.best_annotation()
        .map(|a| a.name.clone())
        .unwrap_or_default()
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

    /// A contig too short / with no skew signal yields no oriC (oricle gates on
    /// skew amplitude); the adapter simply surfaces the empty result.
    #[test]
    fn no_signal_contig_is_empty() {
        let seq = b"ACGT".repeat(50); // 200 bp, no cumulative-skew trough
        let feats = detect(&[contig("c1", &seq)], &[]);
        assert!(feats.is_empty());
    }

    /// gene_name prefers the gene symbol, then product, then annotation name.
    #[test]
    fn gene_name_prefers_symbol_then_product() {
        let mut f = Feature {
            kind: FeatureKind::Cds,
            contig: "c1".into(),
            id: "x".into(),
            start: 1,
            end: 9,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        };
        f.func.product = Some("Chromosomal replication initiator protein DnaA".into());
        assert!(gene_name(&f).to_lowercase().contains("dnaa"));
        f.func.gene = Some("dnaA".into());
        assert_eq!(gene_name(&f), "dnaA");
    }
}
