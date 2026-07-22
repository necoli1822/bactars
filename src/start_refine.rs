//! Homology-guided start-codon refinement (PGAP-style, opt-in via `--refine-starts`).
//!
//! The ab-initio caller (rustygal / Prodigal) maximises the coding score and so
//! tends to pick CDS starts UPSTREAM of the biological start — over-extending the
//! 5' end. NCBI PGAP corrects this using protein homology: it trims the start to
//! the position where a conserved full-length protein family actually begins.
//!
//! We approximate that with the one homology signal rustyhmmer exposes: the HMM
//! alignment ENVELOPE on the protein. For a CDS whose best **ncbifams** hit (a
//! full-length equivalog family) has its envelope beginning `env_from` residues
//! into the protein (`env_from > 1`) yet still reaching the C-terminus (so the
//! model is matched end-to-end — an N-terminal over-extension, not an internal
//! domain), we trim the start to the in-frame start codon nearest that anchor.
//!
//! Strictly conservative: only ncbifams hits, only near-full-length matches, only
//! trims (never extends), bounded trim length, contig-edge partials untouched.
//! `new_protein = "M" + old_aa[k+1..]` — trimming k N-terminal codons just drops
//! the first k residues and re-marks the new start as Met, so no re-translation of
//! the body is needed; the contig is read only to locate the start codon.
//!
//! STATUS (2026-07-19): EXPERIMENTAL, off by default (`--refine-starts`). On a
//! 3-genome RefSeq-concordance A/B this does NOT improve start agreement — it is
//! slightly negative (start-match ~88% → ~87-88%). Root cause: rustyhmmer 0.1.3
//! exposes only the target-side alignment ENVELOPE (`env_from`), not the model
//! coordinate (`hmm_from`). HMMER envelopes routinely inset a few residues at the
//! N-terminus even for correctly-started genes, so `env_from - 1` overestimates the
//! true 5' over-extension and the trim overshoots RefSeq's start. A faithful PGAP
//! "align N-terminus to model start" refinement needs `hmm_from` — i.e. a rustyhmmer
//! change to record per-domain model bounds. This module is kept as the scaffold for
//! that upgrade; until then leave it OFF (default pipeline is byte-unchanged).

use crate::fasta::Contig;
use crate::feature::{Feature, FeatureKind};
use std::collections::HashMap;

/// Per-feature homology anchor captured during HMM annotation:
/// `(env_from, env_to, model_len, seq_score)` of the best ncbifams hit (1-based
/// envelope coords on the protein).
pub type Anchor = (i64, i64, usize, f32);

// --- tunables (conservative) ---
const MIN_OVERHANG_RES: i64 = 12; // require a LARGE N-terminal overhang (HMMER envelopes
                                  // inset a few residues even on correct starts, so small
                                  // env_from>1 is noise, not over-extension)
const MAX_TRIM_RES: i64 = 60; // never trim more than 60 aa (guard against a wrong model)
const CTERM_MARGIN_RES: i64 = 5; // envelope must reach within 5 aa of the protein C-terminus
const MIN_MODEL_COV: f64 = 0.70; // hit must cover >=70% of the model (full-length equivalog)
const MIN_SCORE: f32 = 25.0; // only trust reasonably strong hits
const MIN_RESULT_RES: i64 = 40; // keep the trimmed protein at least this long
// The true start sits a few residues UPSTREAM of the envelope start (envelope inset),
// so bias the start-codon search upstream of the anchor and never past it.
const SEARCH_BACK: i64 = 8; // scan start codons from anchor-8 ..
const SEARCH_FWD: i64 = -2; //             .. to anchor-2 (upstream-biased)

#[derive(Default, Debug)]
pub struct RefineStats {
    pub examined: usize,
    pub trimmed: usize,
    pub total_res_trimmed: i64,
}

/// Trim over-extended CDS starts using the ncbifams envelope anchors. Mutates
/// `features` in place (start/end/aa). Returns counts for reporting.
pub fn refine_starts(
    contigs: &[Contig],
    features: &mut [Feature],
    anchors: &HashMap<usize, Anchor>,
) -> RefineStats {
    let cmap: HashMap<&str, &Contig> = contigs.iter().map(|c| (c.name.as_str(), c)).collect();
    let mut st = RefineStats::default();
    for (&fi, &(env_from, env_to, model_len, score)) in anchors {
        let f = match features.get_mut(fi) {
            Some(f) => f,
            None => continue,
        };
        if f.kind != FeatureKind::Cds || f.partial5 {
            continue;
        }
        let aa = match &f.aa {
            Some(a) => a,
            None => continue,
        };
        let aalen = aa.len() as i64;
        st.examined += 1;
        let overhang = env_from - 1;
        if overhang < MIN_OVERHANG_RES || overhang > MAX_TRIM_RES {
            continue;
        }
        if score < MIN_SCORE || model_len == 0 {
            continue;
        }
        // model must be matched end-to-end (reaches the C-terminus) — otherwise
        // env_from>1 is a legitimate internal-domain start, not an over-extension.
        if env_to < aalen - CTERM_MARGIN_RES {
            continue;
        }
        let cov = (env_to - env_from + 1) as f64 / model_len as f64;
        if cov < MIN_MODEL_COV {
            continue;
        }
        let contig = match cmap.get(f.contig.as_str()) {
            Some(c) => *c,
            None => continue,
        };
        let five = if f.strand == -1 { f.end } else { f.start };
        // pick the in-frame start codon nearest the envelope anchor
        let k = match pick_trim_codon(contig, f.strand, five, overhang) {
            Some(k) => k,
            None => continue,
        };
        if k < 1 || k > MAX_TRIM_RES || aalen - k < MIN_RESULT_RES {
            continue;
        }
        // apply: move the 5' end downstream by k codons; body residues unchanged.
        if f.strand == -1 {
            f.end -= 3 * k;
        } else {
            f.start += 3 * k;
        }
        f.aa = Some(format!("M{}", &aa[(k as usize + 1)..]));
        st.trimmed += 1;
        st.total_res_trimmed += k;
    }
    st
}

/// Strand-oriented codon `k` counted from the 5' end (`five` is the 5' contig
/// coord, 1-based). Returns the (possibly reverse-complemented) triplet, or None
/// if out of bounds.
fn codon_at(contig: &Contig, strand: i8, five: i64, k: i64) -> Option<[u8; 3]> {
    let s = &contig.seq;
    if strand == -1 {
        // 5' at the high coord; codon k occupies forward [p0..p0+3], reverse-complemented.
        let p0 = (five - 1) - 3 * k - 2;
        if p0 < 0 || (p0 as usize + 3) > s.len() {
            return None;
        }
        let p = p0 as usize;
        Some([comp(s[p + 2]), comp(s[p + 1]), comp(s[p])])
    } else {
        let p0 = (five - 1) + 3 * k;
        if p0 < 0 || (p0 as usize + 3) > s.len() {
            return None;
        }
        let p = p0 as usize;
        Some([s[p], s[p + 1], s[p + 2]])
    }
}

/// Choose the codon index (>=1) to trim to: the in-frame start codon nearest the
/// homology anchor `overhang`, scanning `[overhang-SEARCH_BACK, overhang+SEARCH_FWD]`.
/// Ties break toward the anchor; among equals, prefer ATG over GTG/TTG.
fn pick_trim_codon(contig: &Contig, strand: i8, five: i64, overhang: i64) -> Option<i64> {
    let mut best: Option<(i64, u8, i64)> = None; // (dist_to_anchor, codon_pref, k)
    let lo = (overhang - SEARCH_BACK).max(1);
    let hi = overhang + SEARCH_FWD;
    let mut k = lo;
    while k <= hi {
        if let Some(c) = codon_at(contig, strand, five, k) {
            let pref = match &c {
                b"ATG" => 0u8,
                b"GTG" | b"TTG" => 1u8,
                _ => 2u8,
            };
            if pref < 2 {
                let dist = (k - overhang).abs();
                let cand = (dist, pref, k);
                if best.is_none() || cand < best.unwrap() {
                    best = Some(cand);
                }
            }
        }
        k += 1;
    }
    best.map(|(_, _, k)| k)
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}
