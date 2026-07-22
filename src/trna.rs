//! tRNA detection via trnascan-rs (Rust tRNAscan-SE 2.0, embedded infernox).
//!
//! Wired in-process as a library. We drive the scanner with the exact
//! bacterial config that is byte-identical to C tRNAscan-SE 2.0 + Infernal
//! 1.1.5 (`-B -H --detail`): bacterial CM, default score cutoff 20.0, HMM-score
//! and detailed-isotype columns enabled. The scanner reads the genome FASTA
//! itself (via its own `SeqFileReader`), matching the verified CLI path exactly.

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use trnascan_rs::core::TrnaScanner;
use trnascan_rs::squid::SeqFileReader;
use trnascan_rs::trna::Strand as TStrand;

/// One detected tRNA, flattened from trnascan-rs's `TRna` into what bactars needs.
pub struct TrnaCall {
    pub contig: String,
    /// 1-based inclusive, `start <= end`.
    pub start: i64,
    pub end: i64,
    /// `+1` forward, `-1` reverse.
    pub strand: i8,
    /// Amino-acid isotype, e.g. `"Ala"` (`"Undet"` / `"Pseudo"` possible).
    pub isotype: String,
    /// Anticodon, e.g. `"TGC"`.
    pub anticodon: String,
    /// Infernal / covariance-model bit score.
    pub score: f64,
    /// HMM component score.
    pub hmm_score: f64,
    /// Secondary-structure component score.
    pub ss_score: f64,
}

/// Detect tRNAs in a genome FASTA using the bacterial (`-B`) faithful pipeline.
///
/// `models_dir` must point at trnascan-rs's `data/models` directory (the bundled
/// bacterial covariance + isotype models).
pub fn detect_trnas(genome_path: &str, models_dir: &str) -> Result<Vec<TrnaCall>, String> {
    // Bacterial mode, default cutoff 20.0 — the byte-parity configuration.
    let mut scanner = TrnaScanner::with_models_dir('B', 20.0, models_dir)
        .map_err(|e| format!("trnascan-rs init ({models_dir}): {e}"))?;
    scanner.set_quiet(true);
    // `-H`: populate the HMM / 2'-structure score breakdown.
    scanner.set_get_hmm_score(true);
    // `--detail`: populate isotype + note columns.
    scanner.set_detail(true);

    let mut reader =
        SeqFileReader::open(genome_path).map_err(|e| format!("{genome_path}: {e}"))?;
    loop {
        match reader.read_seq() {
            Ok(Some((seq, sqinfo))) => {
                scanner
                    .scan_sequence(&seq, &sqinfo)
                    .map_err(|e| format!("scan {}: {e}", sqinfo.name))?;
            }
            Ok(None) => break,
            Err(e) => return Err(format!("{genome_path}: {e}")),
        }
    }

    let calls = scanner
        .trna_results()
        .iter()
        .map(|t| TrnaCall {
            // Normalize to the canonical first-token contig id (matching
            // `fasta::read_genome` and tmrna.rs) so this tRNA's contig lines up
            // with the same locus elsewhere and `resolve()` overlap-dedup works.
            contig: t
                .seqname
                .split([' ', '\t'])
                .next()
                .unwrap_or(&t.seqname)
                .to_string(),
            start: t.start,
            end: t.end,
            strand: match t.strand {
                TStrand::Minus => -1,
                _ => 1,
            },
            isotype: t.isotype.clone(),
            anticodon: t.anticodon.clone(),
            score: t.score,
            hmm_score: t.hmm_score,
            ss_score: t.ss_score,
        })
        .collect();
    Ok(calls)
}

/// Turn detected tRNAs into `Feature`s, tagged with a `trnascan-rs` annotation
/// carrying the isotype/anticodon identity and the CM/HMM/2'-str score breakdown.
pub fn trna_features(calls: Vec<TrnaCall>) -> Vec<Feature> {
    let mut per_contig: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    calls
        .into_iter()
        .map(|c| {
            let n = per_contig.entry(c.contig.clone()).or_insert(0);
            *n += 1;
            let id = format!("{}_trna{}", c.contig, n);
            let name = format!("tRNA-{}", c.isotype);
            Feature {
                kind: FeatureKind::Trna,
                contig: c.contig,
                id,
                start: c.start,
                end: c.end,
                strand: c.strand,
                aa: None,
                partial5: false,
                partial3: false,
                annotations: vec![Annotation {
                    source: "trnascan-rs".into(),
                    accession: c.anticodon,
                    name,
                    score: c.score as f32,
                    // trnascan-rs reports covariance-model / HMM bit scores; its
                    // public `TRna` result carries no E-value, so there is none to
                    // plumb here (the internal cmsearch E-value is not exposed).
                    evalue: None,
                    ref_len: None,
                }],
                func: Functional::default(),
            }
        })
        .collect()
}
