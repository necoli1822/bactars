//! CDS prediction via rustygal (in-process): train → find_genes → translate.

use crate::fasta::Contig;
use rustygal::api::{find_genes, train_on_sequence, GeneFinderConfig};
use rustygal::{bitmap, gene};

/// A predicted coding sequence with its protein translation.
pub struct Cds {
    /// `{contig}_{index}` (prodigal-style id).
    pub id: String,
    pub contig: String,
    /// 1-based inclusive left coordinate on the contig.
    pub start: i64,
    /// 1-based inclusive right coordinate on the contig.
    pub end: i64,
    /// `+1` / `-1`.
    pub strand: i8,
    /// 5' end incomplete — no start codon (relative to `strand`). Parsed from the
    /// prodigal/rustygal header `partial=XY` attribute.
    pub partial5: bool,
    /// 3' end incomplete — no stop codon (relative to `strand`).
    pub partial3: bool,
    /// Protein translation (may end without a stop residue).
    pub aa: String,
}

/// Predict CDS across all contigs.
///
/// Two regimes, chosen by the LARGEST single contig:
///  * **Whole-genome (largest contig >= 20 kb):** train one gene model
///    (prodigal single-mode) on ALL contigs concatenated — matching prodigal
///    single-mode, which trains on the whole input rather than a single sequence —
///    then gene-call and translate every contig with the shared model.
///    Self-training is the most accurate choice when there is enough sequence to
///    fit a model.
///  * **Fragmented / short-contig (largest contig < 20 kb):** single-mode
///    self-training is unreliable on so little sequence, so instead we run
///    rustygal's metagenomic gene-finder (the 50 pretrained bins, byte-identical
///    to pyrodigal `-p meta`) per contig — the same `meta_api` engine rust-ise
///    uses for ISOSDB. Genes are called against the best-scoring pretrained bin.
///
/// `trans_table` selects the genetic code for the single-mode (whole-genome)
/// path: `Some(n)` forces NCBI translation table `n`; `None` **auto-detects** it,
/// mirroring prodigal's metagenomic bin-selection criterion — train a model under
/// each candidate table, and keep the table whose genome dynamic-programming **path
/// score** ([`genome_path_score`]) is highest. This is self-validating: a normal
/// bacterium always scores highest under its true table 11, while a Mollicute
/// (TGA=Trp) only wins table 4 when the readthrough model genuinely fits, so
/// auto-detection never regresses a table-11 genome. Auto-detection compares
/// **11** (standard) vs **4** (Mollicutes, TGA=Trp) — the only two the score can tell
/// apart. Table **25** (Gracilibacteria/SR1, TGA=Gly) is gene-structure-identical to
/// table 4 (both make TGA coding), so it always ties table 4 on score and can only be
/// requested explicitly via `Some(25)`; its coordinates match table 4 anyway. The
/// metagenomic bins carry their own per-bin tables (incl. table 4), so the
/// short-contig path is already code-aware and ignores `trans_table`.
/// Only genuinely empty input (no usable sequence) is a hard error.
pub fn predict_cds(contigs: &[Contig], trans_table: Option<i32>) -> Result<Vec<Cds>, String> {
    if contigs.is_empty() {
        return Err("no contigs to predict on".to_string());
    }
    let largest = contigs.iter().map(|c| c.seq.len()).max().unwrap_or(0);

    // Short-contig / fragmented regime: route through the real metagenomic mode
    // (its pretrained bins already include table-4 models and pick by score, so
    // the genetic code is auto-selected there regardless of `trans_table`).
    if largest < 20_000 {
        return predict_cds_meta(contigs);
    }

    // Pooled single-mode training input: all contigs concatenated with the same
    // 12-base "TTAATTAATTAA" spacer prodigal/rustygal inserts between records
    // during training (prevents genes from spanning contig boundaries). This is
    // the faithful rustygal 0.2.0 shape — `train_on_sequence` takes one sequence,
    // so we build the concatenation prodigal would have read.
    const SPACER: &[u8] = b"TTAATTAATTAA";
    let total: usize = contigs.iter().map(|c| c.seq.len()).sum::<usize>()
        + SPACER.len() * contigs.len().saturating_sub(1);
    let mut pooled: Vec<u8> = Vec::with_capacity(total);
    for c in contigs {
        if !pooled.is_empty() {
            pooled.extend_from_slice(SPACER);
        }
        pooled.extend_from_slice(&c.seq);
    }

    match trans_table {
        // Forced table: honour it exactly (no auto-detection). `Training` is a large
        // (~550 KB `mot_wt`/`gene_dc`) struct — box it so it lives on the heap and
        // does not inflate this frame (which would overflow a 2 MB worker stack).
        Some(t) => {
            let tinf = Box::new(train_on_sequence(&pooled, t, false)?);
            call_with_model(contigs, &tinf)
        }
        // Auto-detect: compare the prokaryotic candidate codes by the genome's
        // dynamic-programming path score (prodigal's own metagenomic bin-selection
        // metric) and keep the winning table. `>` (not `>=`) keeps table 11 on ties,
        // so the standard code is the default whenever an alternative doesn't clearly win.
        None => {
            // Candidates are 11 (standard) vs 4 (TGA=Trp, Mollicutes). Table 25
            // (TGA=Gly, Gracilibacteria/SR1) is intentionally NOT a candidate: it is
            // gene-structure-identical to table 4 (both make TGA a coding codon), so
            // its path score always TIES table 4 and it can never be distinguished by
            // gene-finding — telling 25 from 4 needs taxonomy (they differ only in the
            // amino acid at TGA, Gly vs Trp). A genome known to be Gracilibacteria can
            // still force it with `--translation-table 25`; its coordinates are correct
            // under table 4 regardless.
            const CANDIDATES: [i32; 2] = [11, 4];
            // `Training` is ~550 KB; box every model so the candidate we retain and the
            // per-iteration model stay on the heap (holding two by value would overflow
            // a 2 MB worker stack — the frame is reserved on entry regardless of input).
            let mut best: Option<(f64, i32, Box<rustygal::training::Training>)> = None;
            let mut scores: Vec<(i32, f64)> = Vec::with_capacity(CANDIDATES.len());
            for &t in &CANDIDATES {
                let tinf = Box::new(train_on_sequence(&pooled, t, false)?);
                let score = genome_path_score(contigs, &tinf);
                scores.push((t, score));
                if best.as_ref().map_or(true, |b| score > b.0) {
                    best = Some((score, t, tinf));
                }
            }
            let (_, chosen, tinf) = best.ok_or("gene-model training produced no candidate")?;
            if chosen != 11 {
                let s: Vec<String> = scores
                    .iter()
                    .map(|(t, sc)| format!("code{t}={sc:.0}"))
                    .collect();
                eprintln!(
                    "[bactars] genetic-code auto-detect: chose translation table {chosen} \
                     (non-standard, Mollicutes; {}). TGA is read as Trp here.",
                    s.join(" ")
                );
            }
            call_with_model(contigs, &tinf)
        }
    }
}

/// Gene-call every contig with a trained model and return the predicted CDS. This
/// is the OUTPUT pass (byte-identical to the previous single-model path); genetic-code
/// selection is done separately by [`genome_path_score`].
fn call_with_model(
    contigs: &[Contig],
    tinf: &rustygal::training::Training,
) -> Result<Vec<Cds>, String> {
    let cfg = GeneFinderConfig::default();
    let mut out = Vec::new();
    for c in contigs {
        let slen = c.seq.len() as i32;
        if slen < 1 {
            continue;
        }
        let result = find_genes(&c.seq, Some(tinf), &cfg)?;
        if result.num_genes == 0 {
            continue;
        }

        // Rebuild the digital sequence bitmaps (find_genes keeps them internal)
        // so we can translate; deterministic from the same input.
        let mut seq = vec![0u8; c.seq.len() / 4 + 1];
        let mut rseq = vec![0u8; c.seq.len() / 4 + 1];
        let mut useq = vec![0u8; c.seq.len() / 8 + 1];
        let mut gc = 0.0;
        bitmap::build_bitmaps(&c.seq, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

        let mut buf: Vec<u8> = Vec::new();
        gene::write_translations(
            &mut buf,
            &result.genes,
            result.num_genes,
            &result.nodes,
            &seq,
            &rseq,
            &useq,
            slen,
            tinf,
            1,
            &c.name,
        );
        parse_protein_fasta(&buf, &c.name, &mut out);
    }
    Ok(out)
}

/// Score a whole genome under one trained model with prodigal's genetic-code
/// selection metric: the sum over contigs of the best dynamic-programming **path
/// score** (`nodes[ipath].score`). This mirrors `find_genes`' internal node pipeline
/// (add → sort → score → record-overlaps → dprog) up to the traceback, and mirrors
/// how the metagenomic mode picks its winning bin — the path score already folds in
/// per-gene coding, start, and intergenic terms, so a wrong-code run that fragments
/// genes at spurious stops scores LOWER even though it emits more (shorter) genes.
/// That is exactly why the naive Σ-of-gene-scores does not work: it rewards fragment
/// count. Runs on the single-mode (≥20 kb) path only.
fn genome_path_score(contigs: &[Contig], tinf: &rustygal::training::Training) -> f64 {
    use rustygal::dprog::dprog;
    use rustygal::node::{
        add_nodes, compare_nodes, record_overlapping_starts, score_nodes, Node, STT_NOD,
    };
    use rustygal::sequence::{Mask, MAX_MASKS};

    // No sequence masking (nmask = 0), matching find_genes' default config.
    let mlist = vec![Mask { begin: 0, end: 0 }; MAX_MASKS];
    let mut total = 0.0f64;
    for c in contigs {
        let slen = c.seq.len() as i32;
        if slen < 1 {
            continue;
        }
        // Pad the 2-bit/4-bit packed buffers: the gene-finding node/score kernels read
        // a full codon (up to `slen+2`) via `base2`, one base past the sequence end, so
        // `slen/4 + 1` is off-by-one for some `slen`. `find_genes` sidesteps this by
        // allocating a fixed MAX_SEQ; here a small fixed margin (8 bytes = 32 bases) is
        // enough and cheap. (`call_with_model` needs no margin — its small buffers only
        // feed post-gene-finding translation, and the node kernels there run inside
        // `find_genes` on its own oversized buffers.)
        let mut seq = vec![0u8; c.seq.len() / 4 + 8];
        let mut rseq = vec![0u8; c.seq.len() / 4 + 8];
        let mut useq = vec![0u8; c.seq.len() / 8 + 8];
        let mut gc = 0.0;
        bitmap::build_bitmaps(&c.seq, slen, &mut seq, &mut rseq, &mut useq, &mut gc);

        let max_slen = if slen > STT_NOD as i32 * 8 {
            (slen / 8) as usize
        } else {
            STT_NOD
        };
        let mut nodes = vec![Node::default(); max_slen];
        // closed = 0, meta flag = 0 (single mode), start-overlap flag = 1 — the exact
        // arguments find_genes uses for the single-model gene call.
        let nn = add_nodes(&seq, &rseq, slen, &mut nodes, 0, &mlist, 0, tinf);
        if nn <= 0 {
            continue;
        }
        nodes[..nn as usize].sort_unstable_by(compare_nodes);
        score_nodes(&seq, &rseq, slen, &mut nodes, nn, tinf, 0, 0);
        record_overlapping_starts(&mut nodes, nn, tinf, 1);
        let ipath = dprog(&mut nodes, nn, tinf, 1);
        if ipath >= 0 {
            total += nodes[ipath as usize].score;
        }
    }
    total
}

/// Short-contig / fragmented regime: gene-call each contig with rustygal's
/// metagenomic mode (the 50 pretrained bins, `== pyrodigal -p meta`) instead of
/// single-mode self-training, which needs ~20 kb to fit a reliable model.
///
/// `meta_bins()` builds the read-only bin set once; `run_meta` picks the
/// best-scoring bin per contig and emits the exact `-a` protein FASTA the binary
/// would (`>{contig}_{i} # begin # end # strand # attrs`), so it feeds straight
/// into the same `parse_protein_fasta` mapper the single-mode path uses — identical
/// 1-based inclusive coords, strand, and translation conventions.
fn predict_cds_meta(contigs: &[Contig]) -> Result<Vec<Cds>, String> {
    if contigs.iter().all(|c| c.seq.is_empty()) {
        return Err("no usable sequence to predict genes on".to_string());
    }
    let largest = contigs.iter().map(|c| c.seq.len()).max().unwrap_or(0);
    eprintln!(
        "[bactars] short contigs (largest {largest} bp < 20000): using rustygal \
         metagenomic mode (pretrained bins, per contig)"
    );

    // Built once, shared read-only across all contigs (as rust-ise/ISOSDB does).
    let meta = rustygal::meta_api::meta_bins();

    let mut out = Vec::new();
    for (i, c) in contigs.iter().enumerate() {
        if c.seq.is_empty() {
            continue;
        }
        // seq_num is the 1-based contig index the binary would assign; only the
        // gff/ID tag uses it — the protein header id is `{short_header}_{n}`, and
        // short_header is the contig name's first token (== c.name here).
        let res = rustygal::meta_api::run_meta(i as i32 + 1, &c.name, &c.seq, &meta);
        parse_protein_fasta(&res.trans_faa, &c.name, &mut out);
    }
    Ok(out)
}

/// Parse the prodigal-format protein FASTA that `write_translations` emits.
/// Header: `>{contig}_{i} # {begin} # {end} # {strand} # {attrs}`.
fn parse_protein_fasta(buf: &[u8], contig: &str, out: &mut Vec<Cds>) {
    let text = String::from_utf8_lossy(buf);
    // (id, start, end, strand, partial5, partial3)
    type Pending = (String, i64, i64, i8, bool, bool);
    let mut pending: Option<Pending> = None;
    let mut aa = String::new();

    let flush = |p: &mut Option<Pending>, aa: &mut String, out: &mut Vec<Cds>| {
        if let Some((id, start, end, strand, partial5, partial3)) = p.take() {
            out.push(Cds {
                id,
                contig: contig.to_string(),
                start,
                end,
                strand,
                partial5,
                partial3,
                aa: std::mem::take(aa),
            });
        }
        aa.clear();
    };

    for line in text.lines() {
        if let Some(hdr) = line.strip_prefix('>') {
            flush(&mut pending, &mut aa, out);
            let mut parts = hdr.split(" # ");
            let id = parts.next().unwrap_or("").to_string();
            let start = parts.next().and_then(|s| s.trim().parse().ok()).unwrap_or(0);
            let end = parts.next().and_then(|s| s.trim().parse().ok()).unwrap_or(0);
            let strand = match parts.next().map(|s| s.trim()) {
                Some("-1") => -1,
                _ => 1,
            };
            // Remaining attrs (`ID=..;partial=XY;start_type=..;..`) carry prodigal's
            // partiality: `partial=XY` where X = the LOW (leftmost) coordinate is
            // incomplete, Y = the HIGH (rightmost) coordinate is incomplete — this
            // is the faithful source (rustygal 0.2.0 gene.rs emits it in both the
            // single-mode and metagenomic paths). Map the coordinate-based digits to
            // biological 5'/3' by strand: on `+` the 5' end is the low coordinate, on
            // `-` the 5' end is the high coordinate.
            let attrs = parts.next().unwrap_or("");
            let (lo_partial, hi_partial) = parse_partial(attrs);
            let (partial5, partial3) = if strand == -1 {
                (hi_partial, lo_partial)
            } else {
                (lo_partial, hi_partial)
            };
            pending = Some((id, start, end, strand, partial5, partial3));
        } else {
            aa.push_str(line.trim());
        }
    }
    flush(&mut pending, &mut aa, out);
}

/// Extract prodigal's `partial=XY` from a header attribute string, returning
/// `(low_coord_incomplete, high_coord_incomplete)`. X is the first digit (low /
/// leftmost coordinate), Y the second (high / rightmost coordinate); `1` means
/// incomplete. Absent/malformed → `(false, false)` (treat as complete).
fn parse_partial(attrs: &str) -> (bool, bool) {
    if let Some(rest) = attrs.split("partial=").nth(1) {
        let mut digits = rest.chars();
        let x = digits.next();
        let y = digits.next();
        return (x == Some('1'), y == Some('1'));
    }
    (false, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    // A real ~900 bp coding fragment (E. coli MG1655 thrA, from its ATG start).
    // Guaranteed to score as a gene, so the metagenomic path must call it.
    const THRA_900: &str = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGAAAAACTGCTGGCAGTGGGGCATTACCTCGAATCTACCGTCGATATTGCTGAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCGCCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGCTGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGTCAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCGCTAAAGTTCTTCACCCCCGCACCATTACCCCCATCGCCCAGTTCCAGATCCCTTGCCTGATTAAAAATACCGGAAATCCTCAAGCACCAGGTACGCTCATTGGTGCCAGCCGTGAT";

    // Build a short contig (< 20 kb) that contains an obvious ORF: AT-rich flanks
    // around the real coding fragment so the gene is internal and complete.
    fn short_contig(name: &str) -> Contig {
        let flank = "AATTTAAATTTAAATTTAAATTTAAATTTAAA";
        let mut seq = String::new();
        seq.push_str(flank);
        seq.push_str(THRA_900);
        seq.push_str("TAA"); // in-frame stop right after the fragment
        seq.push_str(flank);
        Contig {
            name: name.to_string(),
            seq: seq.into_bytes(),
        }
    }

    #[test]
    fn meta_mode_calls_genes_on_short_contigs() {
        // Three short contigs, each well under the 20 kb single-mode threshold.
        let contigs = vec![
            short_contig("ctg1"),
            short_contig("ctg2"),
            short_contig("ctg3"),
        ];
        let largest = contigs.iter().map(|c| c.seq.len()).max().unwrap();
        assert!(largest < 20_000, "test setup: contigs must be short");

        let cds = predict_cds(&contigs, Some(11)).expect("meta mode should not error");
        assert!(
            !cds.is_empty(),
            "metagenomic path called zero genes on short contigs with obvious ORFs"
        );

        // Every CDS must have sane 1-based inclusive coords, a valid strand, and a
        // non-empty translation — the same conventions the single-mode path yields.
        for g in &cds {
            let clen = contigs
                .iter()
                .find(|c| c.name == g.contig)
                .map(|c| c.seq.len())
                .expect("CDS references a known contig");
            assert!(g.start >= 1, "start must be 1-based: {}", g.start);
            assert!(g.end >= g.start, "end < start: {} < {}", g.end, g.start);
            assert!(g.end as usize <= clen, "coord runs past contig end");
            assert!(g.strand == 1 || g.strand == -1, "bad strand {}", g.strand);
            assert!(!g.aa.is_empty(), "empty protein translation");
            assert!(g.id.starts_with(&g.contig), "id must be prefixed by contig");
        }
    }

    #[test]
    fn parses_partial_and_maps_to_biological_ends() {
        // `partial=XY`: X = low (leftmost) coord incomplete, Y = high coord.
        // + strand: 5'=low(X), 3'=high(Y). − strand: 5'=high(Y), 3'=low(X).
        let faa = b">ctg_1 # 10 # 99 # 1 # ID=1;partial=10;start_type=Edge;\nMKV\n\
                    >ctg_2 # 200 # 300 # 1 # ID=2;partial=01;start_type=ATG;\nMKV\n\
                    >ctg_3 # 400 # 500 # -1 # ID=3;partial=10;start_type=ATG;\nMKV\n\
                    >ctg_4 # 600 # 700 # 1 # ID=4;partial=00;start_type=ATG;\nMKV\n";
        let mut out = Vec::new();
        parse_protein_fasta(faa, "ctg", &mut out);
        assert_eq!(out.len(), 4);
        // + strand, X=1 (low incomplete) → 5' partial.
        assert!(out[0].partial5 && !out[0].partial3);
        // + strand, Y=1 (high incomplete) → 3' partial.
        assert!(!out[1].partial5 && out[1].partial3);
        // − strand, X=1 (low incomplete) → 3' partial (biological ends swap).
        assert!(!out[2].partial5 && out[2].partial3);
        // Complete gene → neither.
        assert!(!out[3].partial5 && !out[3].partial3);
    }

    #[test]
    fn empty_input_errors() {
        assert!(predict_cds(&[], None).is_err());
        let empties = vec![
            Contig { name: "a".into(), seq: Vec::new() },
            Contig { name: "b".into(), seq: Vec::new() },
        ];
        assert!(predict_cds(&empties, None).is_err());
    }
}
