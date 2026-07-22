//! tmRNA detection via oxagorn (Rust ARAGORN, in-process).
//!
//! tRNA is handled by trnascan-rs; here we run oxagorn in tmRNA-only mode
//! (`api::detect` with `SearchOptions{trna:false, tmrna:true}`) and emit
//! `Feature{Tmrna}`. oxagorn's gene calls are byte-identical to C ARAGORN
//! v1.2.41 (verified: 89/89 on MG1655 via the published `api::detect`).

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use oxagorn::api::{detect, SearchOptions};
use std::collections::HashMap;

/// Detect tmRNA genes in a genome FASTA, as `Feature{Tmrna}`.
pub fn detect_tmrna(genome_path: &str, threads: usize) -> Result<Vec<Feature>, String> {
    // tmRNA only (tRNA is trnascan-rs's job); bacterial defaults otherwise.
    let opts = SearchOptions {
        trna: false,
        ..Default::default()
    };
    let genes = detect(genome_path, &opts, threads)?;

    let mut per_contig: HashMap<String, usize> = HashMap::new();
    let mut feats = Vec::new();
    for g in genes {
        // TMRNA == 1; defensive filter in case tRNA slips through.
        if g.genetype != 1 {
            continue;
        }
        // oxagorn carries the FULL FASTA header as `seqname`; normalize to the
        // canonical first-token contig id (matching `fasta::read_genome`) so this
        // tmRNA's contig lines up with the infernox RF00023 ncRNA call over the
        // same locus and `resolve()` can dedup them (else both survive) — and so
        // the emitted id/contig aren't the mangled full header.
        let seqname = g
            .seqname
            .split([' ', '\t'])
            .next()
            .unwrap_or(&g.seqname)
            .to_string();
        let n = per_contig.entry(seqname.clone()).or_insert(0);
        *n += 1;
        let (start, end) = if g.start <= g.stop {
            (g.start, g.stop)
        } else {
            (g.stop, g.start)
        };
        feats.push(Feature {
            kind: FeatureKind::Tmrna,
            contig: seqname.clone(),
            id: format!("{}_tmrna{}", seqname, n),
            start,
            end,
            strand: if g.comp == 1 { -1 } else { 1 },
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: "oxagorn".into(),
                accession: String::new(),
                name: "tmRNA".into(),
                score: g.energy as f32,
                // Score-only: oxagorn (ARAGORN) reports a fold energy, no E-value.
                evalue: None,
                ref_len: None,
            }],
            func: Functional::default(),
        });
    }
    Ok(feats)
}
