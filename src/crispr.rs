//! CRISPR-array detection — adapter over the published **`minced-rs`** crate
//! (byte-parity pure-Rust MinCED port).
//!
//! `minced-rs` is intentionally dependency-free and single-threaded. Each contig's
//! scan is independent, so we parallelise across contigs here with rayon (order
//! preserving → the flattened result reproduces MinCED's ordering, byte-parity
//! intact) and map `minced_rs::Feature` into the bactars [`Feature`] model.

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use rayon::prelude::*;

/// Detect CRISPR arrays in `contigs` (`(name, seq)` pairs), as `Feature{Crispr}`.
/// Contigs are scanned in parallel; results stay in contig then position order.
pub fn detect_crispr(contigs: &[(String, String)]) -> Result<Vec<Feature>, String> {
    let per_contig: Vec<Vec<Feature>> = contigs
        .par_iter()
        .map(|c| {
            // minced-rs numbers/labels arrays per input slice, so hand it one
            // contig at a time — identical to a whole-genome call, just parallel.
            let one = std::slice::from_ref(c);
            minced_rs::detect_crispr(one)
                .map(|fs| fs.into_iter().map(map_feature).collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>, String>>()?;
    // minced-rs numbers arrays per input slice, so each per-contig scan restarts at
    // "CRISPR1" — ids collide across contigs. Renumber globally (CRISPR1..N) in the
    // stable contig-then-position order the flatten preserves, so every array id is
    // unique for GFF3/GBFF.
    let mut out: Vec<Feature> = per_contig.into_iter().flatten().collect();
    for (i, f) in out.iter_mut().enumerate() {
        f.id = format!("CRISPR{}", i + 1);
    }
    Ok(out)
}

/// Serialise detected CRISPR arrays as MinCED's `-gff` output, byte-for-byte
/// (mirrors `minced_rs::to_minced_gff` on the bactars [`Feature`] model). Kept
/// here so the byte-parity test and `crispr_only` example exercise the full
/// bactars path; verified identical to MinCED's golden.
pub fn to_minced_gff(features: &[Feature]) -> String {
    let mut out = String::from("##gff-version 3\n");
    let mut n = 0usize;
    for f in features.iter().filter(|f| f.kind == FeatureKind::Crispr) {
        n += 1;
        // Guard against a CRISPR feature carrying no annotation (defensive: MinCED
        // always attaches the repeat record, but never index [0] blindly).
        let Some(a) = f.annotations.first() else {
            continue;
        };
        out.push_str(&format!(
            "{}\tminced:{}\trepeat_region\t{}\t{}\t{}\t.\t.\tID=CRISPR{};rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq={}\n",
            f.contig, minced_rs::MINCED_VERSION, f.start, f.end, a.score as i64, n, a.accession
        ));
    }
    out
}

/// Map a `minced_rs::Feature` into the bactars [`Feature`] (structurally identical
/// models; the CRISPR kind is fixed).
fn map_feature(f: minced_rs::Feature) -> Feature {
    Feature {
        kind: FeatureKind::Crispr,
        contig: f.contig,
        id: f.id,
        start: f.start,
        end: f.end,
        strand: f.strand,
        aa: f.aa,
        partial5: false,
        partial3: false,
        annotations: f
            .annotations
            .into_iter()
            .map(|a| Annotation {
                source: a.source,
                accession: a.accession,
                name: a.name,
                score: a.score,
                // minced-rs is a repeat-structure caller (no statistical search): its
                // annotation E-value is the NaN sentinel. Carry a real value only if
                // one is ever present; otherwise `None`.
                evalue: (!a.evalue.is_nan()).then_some(a.evalue),
                ref_len: None,
            })
            .collect(),
        func: Functional::default(),
    }
}
