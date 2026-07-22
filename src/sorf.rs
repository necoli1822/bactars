//! Small ORF (sORF) detection: enumerate short ORFs below the gene caller's
//! floor, keep only those with a protein-homology hit (Bakta-style sORF).
//!
//! rustygal (the Prodigal port) has a minimum gene length (~90 nt / 30 codons),
//! so genuine short proteins are never called. This stage does a 6-frame scan
//! for short ORFs (start codon -> in-frame stop) in a window BELOW that floor,
//! then keeps only those that get a real protein-homology hit against `hmm_db`
//! via rustyhmmer. The homology filter is what separates credible sORFs from
//! the enormous background of random short ORFs.

use crate::fasta::Contig;
use crate::feature::{Annotation, Feature, FeatureKind, Functional};

/// Minimum protein length (aa, inclusive) for a candidate sORF. Includes the
/// start Met. Bakta's shortest homology-supported sORFs are ~7-8 aa (e.g. leader
/// peptides), so 8 is the practical floor; below this even a hit is not credible.
pub const MIN_AA: usize = 8;
/// Maximum protein length (aa, exclusive) for a candidate sORF. This sits below
/// / around the gene caller's ~30-codon floor with headroom, so we only look at
/// ORFs rustygal would have missed. Bakta's sORFs top out near 29 aa; keeping a
/// little slack (50) lets us also recover homology-supported 30-49 aa ORFs the
/// caller skipped without materially inflating the candidate set.
pub const MAX_AA: usize = 50;

/// Fallback E-value threshold for the homology filter. A candidate that clears
/// the trusted gathering cutoff of *any* model is kept outright; a candidate that
/// does not is still kept if it has an HMM hit at or below this E-value. Bakta
/// admits sORFs on protein-similarity support (UPS/IPS/PSC), not only trusted-
/// cutoff HMM hits, so this looser tier recovers the many real short proteins
/// whose family model has no GA cutoff (or a conservative one) yet still produce
/// an unambiguous alignment. 1e-3 keeps the false-positive rate low.
pub const SORF_EVALUE: f64 = 1e-3;

/// Minimum query coverage (fraction of the short ORF spanned by the alignment) for
/// a PSC/UniRef90 hit to rescue a candidate that failed the HMM filter. Bakta's
/// UPS/IPS/PSC tier accepts a short protein when a UniRef90 reference covers
/// essentially all of it; requiring most of the ORF (0.8) keeps partial/spurious
/// hits out. This is what lets the PSC path recover the many real sORFs whose only
/// support is protein similarity (no ncbifams/Pfam HMM).
pub const SORF_PSC_MIN_QCOV: f64 = 0.8;

/// Minimum TARGET coverage for the sORF PSC path. THIS is the decisive false-positive
/// gate: with query-coverage alone, a spurious short ORF that aligns full-length to a
/// small *fragment* of a large UniRef90 protein passes. Requiring the target to also
/// be mostly covered forces the match to be a similarly-short reference protein — a
/// genuine small gene, not a fragment of a big one. On MG1655, qcov-only (0.8) emitted
/// ~1198 short CDS (100 real / 1098 FP, precision ~8%); adding tcov collapses the FPs.
pub const SORF_PSC_MIN_TCOV: f64 = 0.7;

/// E-value cutoff for the sORF PSC path. Tightened well below the HMM tier because a
/// short (8–50 aa) query finds spurious UniRef90 hits at loose E-values in a 73 M-seq
/// DB; a real small protein clears a much stricter bar.
pub const SORF_PSC_MAX_EVALUE: f64 = 1e-6;

/// A PSC hit only supports a novel sORF if its product carries real functional
/// meaning. A short ORF matching only a *hypothetical / uncharacterized / DUF /
/// domain-containing* UniRef90 entry is weak evidence — that reference is itself an
/// unvalidated prediction, and MG1655's genome is full of short ORFs that hit such
/// entries by chance. On MG1655 this gate cut the sORF PSC false positives 482 -> 55
/// (precision 0.17 -> 0.62) while keeping 90 of 99 true short CDS.
fn is_informative_product(product: &str) -> bool {
    let p = product.trim().to_ascii_lowercase();
    if p.is_empty() {
        return false;
    }
    !(p.contains("hypothetical")
        || p.contains("uncharacter")
        || p.contains("domain-containing")
        || p.starts_with("duf")
        || p.contains(" duf"))
}

/// A candidate short ORF located by the 6-frame scan (pre homology filter).
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    /// Contig this ORF lives on.
    pub contig: String,
    /// 1-based inclusive start (left coordinate on the forward strand).
    pub start: i64,
    /// 1-based inclusive end (right coordinate; includes the stop codon).
    pub end: i64,
    /// `+1` forward, `-1` reverse.
    pub strand: i8,
    /// Protein translation (start codon rendered as `M`, stop excluded).
    pub protein: String,
}

/// Translate one DNA codon with the standard genetic code (table 11 differs
/// only in start codons, which are handled separately). `*` = stop, `X` =
/// codon containing an ambiguity/`N`.
fn translate_codon(c: &[u8]) -> u8 {
    match c {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}

/// Bacterial start codons (ATG plus the alternative GTG/TTG). Any of these, when
/// used as a start, translates to Met.
fn is_start(c: &[u8]) -> bool {
    matches!(c, b"ATG" | b"GTG" | b"TTG")
}

/// Reverse-complement an ASCII nucleotide slice (unknowns -> `N`).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// A raw ORF on a byte slice: 0-based inclusive `[start0, end0]` (end0 includes
/// the stop codon) plus its protein translation.
struct RawOrf {
    start0: usize,
    end0: usize,
    protein: String,
}

/// Scan the three forward frames of `seq` for short ORFs. An ORF runs from a
/// start codon to the first in-frame stop.
///
/// For each stop-delimited segment we normally take the most-upstream start (the
/// MAXIMAL ORF ending at that stop). BUT: if that maximal ORF overflows the window
/// (length >= [`MAX_AA`]), the maximal ORF alone would hide a genuine short ORF
/// beginning at a downstream in-frame ATG inside it. In that case we instead
/// recover the LONGEST nested ORF that fits the window — the most-upstream start
/// whose length is `< MAX_AA` — so a real short protein nested in a longer
/// (discarded) reading is not invisible.
///
/// At most ONE ORF is emitted per stop (the maximal if it fits, else the longest
/// nested that fits), which keeps the candidate set bounded (no Russian-doll
/// overlaps) and leans on the downstream homology filter + existing-CDS overlap
/// drop to control false positives. Only ORFs whose protein length is in
/// `[MIN_AA, MAX_AA)` are returned.
fn scan_frames(seq: &[u8]) -> Vec<RawOrf> {
    let mut out = Vec::new();
    for frame in 0..3 {
        // 0-based nt indices of start codons seen since the last in-frame stop, in
        // ascending (upstream-first) order.
        let mut starts: Vec<usize> = Vec::new();
        let mut i = frame;
        while i + 3 <= seq.len() {
            let codon = &seq[i..i + 3];
            if is_start(codon) {
                starts.push(i);
            } else if translate_codon(codon) == b'*' {
                emit_segment_orf(seq, &starts, i, &mut out);
                starts.clear();
            }
            i += 3;
        }
        // A trailing run with no in-frame stop is an unclosed ORF -> not emitted.
    }
    out
}

/// Choose and emit at most one ORF for the stop-delimited segment whose in-frame
/// start codons are `starts` (ascending 0-based nt indices) and whose terminating
/// stop codon begins at 0-based nt index `stop0`. See [`scan_frames`] for the
/// maximal-vs-nested selection rule.
fn emit_segment_orf(seq: &[u8], starts: &[usize], stop0: usize, out: &mut Vec<RawOrf>) {
    let Some(&first) = starts.first() else {
        return;
    };
    // Protein length (aa) of the ORF from start `s` to `stop0` = codons before the
    // stop, including the start Met.
    let aa_len = |s: usize| (stop0 - s) / 3;

    let chosen = if aa_len(first) < MAX_AA {
        // Maximal ORF already fits the window: original behaviour.
        Some(first)
    } else {
        // Maximal ORF overflows: recover the longest nested ORF that fits, i.e. the
        // most-upstream start whose length is < MAX_AA.
        starts.iter().copied().find(|&s| aa_len(s) < MAX_AA)
    };

    if let Some(s) = chosen {
        let alen = aa_len(s);
        if alen >= MIN_AA && alen < MAX_AA {
            let mut protein = String::with_capacity(alen);
            protein.push('M'); // any bacterial start => Met
            let mut j = s + 3;
            while j < stop0 {
                protein.push(translate_codon(&seq[j..j + 3]) as char);
                j += 3;
            }
            out.push(RawOrf {
                start0: s,
                end0: stop0 + 2, // include the stop codon
                protein,
            });
        }
    }
}

/// Enumerate all short ORFs on both strands of `contig`. Pure (no HMM); the
/// homology filter is applied later in [`detect`].
pub fn enumerate_short_orfs(contig: &Contig) -> Vec<Candidate> {
    let mut cands = Vec::new();
    let len = contig.seq.len();

    // Forward strand: 0-based [start0, end0] -> 1-based inclusive directly.
    for orf in scan_frames(&contig.seq) {
        cands.push(Candidate {
            contig: contig.name.clone(),
            start: orf.start0 as i64 + 1,
            end: orf.end0 as i64 + 1,
            strand: 1,
            protein: orf.protein,
        });
    }

    // Reverse strand: scan the reverse-complement, then map coordinates back.
    // A revcomp index p corresponds to forward index (len-1-p), so an ORF on
    // revcomp positions [a, b] covers forward positions [len-1-b, len-1-a].
    let rc = revcomp(&contig.seq);
    for orf in scan_frames(&rc) {
        let start = (len - 1 - orf.end0) as i64 + 1; // 1-based inclusive
        let end = (len - 1 - orf.start0) as i64 + 1;
        cands.push(Candidate {
            contig: contig.name.clone(),
            start,
            end,
            strand: -1,
            protein: orf.protein,
        });
    }

    cands
}

/// True if `cand` substantially overlaps (> 50% of its own length) any existing
/// CDS on the same contig and strand — i.e. it is already called by rustygal.
fn overlaps_existing(cand: &Candidate, existing: &[Feature]) -> bool {
    let cand_len = cand.end - cand.start + 1;
    for f in existing {
        if f.contig != cand.contig || f.strand != cand.strand {
            continue;
        }
        let ov = (cand.end.min(f.end) - cand.start.max(f.start) + 1).max(0);
        if ov * 2 > cand_len {
            return true;
        }
    }
    false
}

/// Detect sORFs: short ORFs (below the gene caller's floor) that carry a real
/// protein-homology hit in `hmm_db`. `existing` are already-called CDS features
/// used to drop ORFs that are simply re-discovered known genes.
///
/// `psc_dir` optionally enables a second homology tier: candidates that fail the
/// HMM filter are searched against the UniRef90 PSC DB (via [`crate::psc::search`]),
/// and any with a well-covered, sane-evalue hit are kept and named. This mirrors
/// Bakta's UPS/IPS/PSC sORF tier and recovers the many real short proteins that
/// have NO ncbifams/Pfam model (the structural gap in the HMM-only path). When
/// `psc_dir` is `None` the behaviour is exactly the original HMM-only path.
pub fn detect(
    contigs: &[Contig],
    existing: &[Feature],
    hmm_db: &str,
    psc_dir: Option<&str>,
    threads: usize,
) -> Vec<Feature> {
    // 1 + 2: enumerate short ORFs on both strands, drop those overlapping an
    // existing CDS.
    let mut cands: Vec<Candidate> = Vec::new();
    for contig in contigs {
        for c in enumerate_short_orfs(contig) {
            if !overlaps_existing(&c, existing) {
                cands.push(c);
            }
        }
    }
    if cands.is_empty() {
        return Vec::new();
    }

    // 3: homology filter. Each candidate gets a temporary unique id so we can
    // map hits back to it. A candidate is kept if it has EITHER a trusted
    // gathering-cutoff hit OR (fallback) a hit at/below `SORF_EVALUE`. Bakta
    // accepts sORFs on general protein-similarity support, not only trusted-
    // cutoff HMM hits, so the looser E-value tier is what lifts recall.
    use rustyhmmer::api::Cutoff;
    let proteins: Vec<(String, String)> = cands
        .iter()
        .enumerate()
        .map(|(i, c)| (format!("cand{i}"), c.protein.clone()))
        .collect();

    // Run the trusted-cutoff pass and the looser E-value pass. The E-value pass
    // is a near-superset, but the trusted pass still recovers GA hits whose
    // E-value is inflated past the threshold by the candidate-set size (the
    // annotator's E-value uses the number of query sequences as its Z), so we
    // union the two rather than relying on E-value alone.
    let ga_hits = match crate::hmm::annotate(&proteins, hmm_db, Cutoff::GatheringGa) {
        Ok(h) => h,
        Err(e) => {
            eprintln!("sorf: homology filter (gathering) failed ({e}); skipping sORF detection");
            return Vec::new();
        }
    };
    let ev_hits = match crate::hmm::annotate(&proteins, hmm_db, Cutoff::Evalue(SORF_EVALUE)) {
        Ok(h) => h,
        Err(e) => {
            eprintln!("sorf: homology filter (evalue) failed ({e}); skipping sORF detection");
            return Vec::new();
        }
    };

    // Best (highest-scoring) hit per candidate id across BOTH passes.
    use std::collections::HashMap;
    let mut best: HashMap<usize, &rustyhmmer::api::HmmHit> = HashMap::new();
    for h in ga_hits.iter().chain(ev_hits.iter()) {
        if let Some(idx) = h
            .target_name
            .strip_prefix("cand")
            .and_then(|s| s.parse::<usize>().ok())
        {
            match best.get(&idx) {
                Some(prev) if prev.seq_score >= h.seq_score => {}
                _ => {
                    best.insert(idx, h);
                }
            }
        }
    }

    // 3b: PSC rescue tier. Candidates that failed the HMM filter are searched
    // against the UniRef90 PSC DB; a candidate whose best hit is well-covered and
    // has a sane E-value is a credible short protein even with no HMM family. Only
    // runs when a psc_dir is configured — otherwise the path is the original one.
    let mut psc_named: HashMap<usize, crate::psc::PscHit> = HashMap::new();
    if let Some(psc_dir) = psc_dir {
        // Throwaway features for the un-HMM'd candidates only; ids echo `cand{i}`
        // so hits map straight back to the candidate index (no double-naming: a
        // candidate already kept by HMM is never searched here).
        let pending: Vec<Feature> = cands
            .iter()
            .enumerate()
            .filter(|(i, _)| !best.contains_key(i))
            .map(|(i, c)| Feature {
                kind: FeatureKind::Cds,
                contig: c.contig.clone(),
                id: format!("cand{i}"),
                start: c.start,
                end: c.end,
                strand: c.strand,
                aa: Some(c.protein.clone()),
                partial5: false,
                partial3: false,
                annotations: Vec::new(),
                func: Functional::default(),
            })
            .collect();
        if !pending.is_empty() {
            for (qid, hit) in crate::psc::search(&pending, psc_dir, threads) {
                let Some(idx) = qid.strip_prefix("cand").and_then(|s| s.parse::<usize>().ok()) else {
                    continue;
                };
                if hit.qcov >= SORF_PSC_MIN_QCOV
                    && hit.tcov >= SORF_PSC_MIN_TCOV
                    && hit.evalue <= SORF_PSC_MAX_EVALUE
                    && is_informative_product(&hit.product)
                {
                    psc_named.insert(idx, hit);
                }
            }
        }
    }

    // 4: emit surviving sORFs as (short) CDS features. A candidate is kept if it
    // has an HMM hit (preferred, unchanged) OR a qualifying PSC hit.
    let mut out = Vec::new();
    let mut counters: HashMap<String, usize> = HashMap::new();
    for (i, cand) in cands.iter().enumerate() {
        // Build the annotation + product from whichever tier supported this ORF.
        let (annotation, product) = if let Some(hit) = best.get(&i) {
            (
                Annotation {
                    source: "rustyhmmer:sorf".to_string(),
                    accession: hit.query_acc.clone(),
                    name: hit.query_name.clone(),
                    score: hit.seq_score,
                    evalue: Some(hit.seq_evalue),
                    ref_len: None,
                },
                hit.query_name.clone(),
            )
        } else if let Some(hit) = psc_named.get(&i) {
            (
                Annotation {
                    source: "psc:uniref90".to_string(),
                    accession: hit.uniref_id.clone(),
                    name: hit.product.clone(),
                    score: hit.bits,
                    evalue: Some(hit.evalue),
                    ref_len: (hit.tlen > 0).then_some(hit.tlen),
                },
                hit.product.clone(),
            )
        } else {
            continue;
        };

        let n = counters.entry(cand.contig.clone()).or_insert(0);
        *n += 1;
        let id = format!("{}_sorf{}", cand.contig, *n);

        out.push(Feature {
            kind: FeatureKind::Cds,
            contig: cand.contig.clone(),
            id,
            start: cand.start,
            end: cand.end,
            strand: cand.strand,
            aa: Some(cand.protein.clone()),
            partial5: false,
            partial3: false,
            annotations: vec![annotation],
            func: Functional {
                product: Some(product),
                ..Default::default()
            },
        });
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(name: &str, seq: &str) -> Contig {
        Contig {
            name: name.to_string(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    /// (a) A forward short ORF (ATG ... stop) is found with correct coords and
    /// translation.
    #[test]
    fn forward_short_orf() {
        // ATG + 14 x AAA(K) + TAA  => protein "M" + "K"*14 = 15 aa.
        let seq = format!("ATG{}TAA", "AAA".repeat(14));
        let c = contig("c1", &seq);
        let cands = enumerate_short_orfs(&c);
        let fwd: Vec<_> = cands.iter().filter(|x| x.strand == 1).collect();
        assert_eq!(fwd.len(), 1, "exactly one forward ORF");
        let orf = fwd[0];
        assert_eq!(orf.start, 1);
        assert_eq!(orf.end, seq.len() as i64); // includes the stop codon
        assert_eq!(orf.protein, format!("M{}", "K".repeat(14)));
        assert_eq!(orf.protein.len(), 15);
    }

    /// (b) The same ORF on the reverse strand is found with mapped coords.
    #[test]
    fn reverse_short_orf() {
        let fwd_seq = format!("ATG{}TAA", "AAA".repeat(14));
        // Reverse-complement so the ORF now reads on the minus strand.
        let rc = revcomp(fwd_seq.as_bytes());
        let c = Contig {
            name: "c2".to_string(),
            seq: rc.clone(),
        };
        let cands = enumerate_short_orfs(&c);
        let rev: Vec<_> = cands.iter().filter(|x| x.strand == -1).collect();
        assert_eq!(rev.len(), 1, "exactly one reverse ORF");
        let orf = rev[0];
        assert_eq!(orf.start, 1);
        assert_eq!(orf.end, rc.len() as i64);
        assert_eq!(orf.protein, format!("M{}", "K".repeat(14)));
    }

    /// (c) Alternative start codons GTG and TTG are accepted and both render as
    /// Met at the start.
    #[test]
    fn alternative_starts() {
        for start in ["GTG", "TTG"] {
            let seq = format!("{start}{}TAA", "AAA".repeat(12));
            let c = contig("c", &seq);
            let cands = enumerate_short_orfs(&c);
            let fwd: Vec<_> = cands.iter().filter(|x| x.strand == 1).collect();
            assert_eq!(fwd.len(), 1, "{start} start should give one ORF");
            assert!(
                fwd[0].protein.starts_with('M'),
                "{start} start must translate as M, got {}",
                fwd[0].protein
            );
            assert_eq!(fwd[0].protein, format!("M{}", "K".repeat(12)));
        }
    }

    /// (d) A candidate overlapping an existing CDS (same contig/strand) is
    /// dropped by the overlap filter.
    #[test]
    fn overlap_with_existing_dropped() {
        let cand = Candidate {
            contig: "c1".to_string(),
            start: 100,
            end: 148,
            strand: 1,
            protein: "M".to_string() + &"K".repeat(15),
        };
        // Existing CDS covering the whole candidate on the same strand.
        let existing = vec![Feature {
            kind: FeatureKind::Cds,
            contig: "c1".to_string(),
            id: "cds1".to_string(),
            start: 50,
            end: 400,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional::default(),
        }];
        assert!(overlaps_existing(&cand, &existing));

        // Opposite strand: not an overlap for our purposes.
        let mut opp = existing.clone();
        opp[0].strand = -1;
        assert!(!overlaps_existing(&cand, &opp));

        // Only a small tail overlap (< 50%): kept.
        let small = vec![Feature {
            strand: 1,
            start: 140,
            end: 500,
            ..existing[0].clone()
        }];
        assert!(!overlaps_existing(&cand, &small));
    }

    /// (e) ORFs longer than the window or shorter than the minimum are excluded;
    /// one inside the window is kept.
    #[test]
    fn length_window_filter() {
        // Too short: below MIN_AA (=8). M + (MIN_AA-2) K = MIN_AA-1 aa.
        let short = contig("s", &format!("ATG{}TAA", "AAA".repeat(MIN_AA - 2)));
        assert!(enumerate_short_orfs(&short)
            .iter()
            .all(|x| x.strand != 1));

        // Too long: 60 aa (M + 59 K). At/above MAX_AA.
        let long = contig("l", &format!("ATG{}TAA", "AAA".repeat(59)));
        assert!(enumerate_short_orfs(&long).iter().all(|x| x.strand != 1));

        // In window: 20 aa (M + 19 K).
        let ok = contig("o", &format!("ATG{}TAA", "AAA".repeat(19)));
        let fwd: Vec<_> = enumerate_short_orfs(&ok)
            .into_iter()
            .filter(|x| x.strand == 1)
            .collect();
        assert_eq!(fwd.len(), 1);
        assert_eq!(fwd[0].protein.len(), 20);
    }

    /// The widened lower bound admits genuinely-short ORFs (down to MIN_AA aa)
    /// that the old MIN_AA=10 floor excluded, while still rejecting anything one
    /// residue shorter. These are the Bakta-style leader-peptide-length sORFs.
    #[test]
    fn min_window_boundary() {
        // Exactly MIN_AA aa (M + (MIN_AA-1) K): kept by the enumerator.
        let at_min = contig("m", &format!("ATG{}TAA", "AAA".repeat(MIN_AA - 1)));
        let fwd: Vec<_> = enumerate_short_orfs(&at_min)
            .into_iter()
            .filter(|x| x.strand == 1)
            .collect();
        assert_eq!(fwd.len(), 1, "an ORF of exactly MIN_AA aa must be enumerated");
        assert_eq!(fwd[0].protein.len(), MIN_AA);

        // One residue shorter (MIN_AA-1 aa): still rejected.
        let below = contig("b", &format!("ATG{}TAA", "AAA".repeat(MIN_AA - 2)));
        assert!(enumerate_short_orfs(&below)
            .iter()
            .all(|x| x.strand != 1));
    }

    /// Nested-ORF recovery: when the maximal ORF for a stop overflows the window
    /// (>= MAX_AA, discarded), the longest nested in-frame ATG ORF that fits is
    /// emitted instead — so a genuine short protein hidden inside a longer reading
    /// is not invisible.
    #[test]
    fn nested_orf_recovered_when_maximal_overflows() {
        // Frame-0 segment: upstream ATG -> 60-aa ORF (> MAX_AA), with a nested
        // in-frame ATG 40 codons downstream -> a 20-aa ORF that fits the window.
        let mut seq = String::from("ATG");
        seq.push_str(&"AAA".repeat(39)); // codons 1..39
        seq.push_str("ATG"); // codon 40 (nt index 120): the nested start
        seq.push_str(&"AAA".repeat(19)); // codons 41..59
        seq.push_str("TAA"); // stop at codon 60 (nt index 180)

        let orfs = scan_frames(seq.as_bytes());
        // The maximal ORF (start0 == 0, 60 aa) is discarded...
        assert!(
            orfs.iter().all(|o| o.start0 != 0),
            "the overflowing maximal ORF must not be emitted"
        );
        // ...and the nested 20-aa ORF (start0 == 120) is recovered.
        let nested: Vec<_> = orfs.iter().filter(|o| o.start0 == 120).collect();
        assert_eq!(nested.len(), 1, "the nested in-window ORF must be recovered");
        assert_eq!(nested[0].protein, format!("M{}", "K".repeat(19)));
        assert_eq!(nested[0].end0, 182); // stop codon (180) + 2
    }

    /// Sanity: an ORF with no in-frame stop is not emitted (must be stop-closed).
    #[test]
    fn unclosed_orf_ignored() {
        let seq = format!("ATG{}", "AAA".repeat(20)); // no stop
        let c = contig("u", &seq);
        assert!(enumerate_short_orfs(&c).iter().all(|x| x.strand != 1));
    }

    #[test]
    fn informative_product_gate() {
        // Real functional names support a novel sORF.
        assert!(is_informative_product("Entericidin B membrane lipoprotein"));
        assert!(is_informative_product("Toxin CcdB"));
        // Placeholder / prediction-only names do not.
        assert!(!is_informative_product("hypothetical protein"));
        assert!(!is_informative_product("Uncharacterized protein YbgS"));
        assert!(!is_informative_product("DUF1382 family protein"));
        assert!(!is_informative_product("PIN domain-containing protein"));
        assert!(!is_informative_product(""));
    }
}
