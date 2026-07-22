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
/// `trans_table` is the genetic code (11 for bacteria/archaea) used by the
/// single-mode path; the metagenomic bins carry their own per-bin translation
/// tables. Only genuinely empty input (no usable sequence) is a hard error.
pub fn predict_cds(contigs: &[Contig], trans_table: i32) -> Result<Vec<Cds>, String> {
    if contigs.is_empty() {
        return Err("no contigs to predict on".to_string());
    }
    let largest = contigs.iter().map(|c| c.seq.len()).max().unwrap_or(0);

    // Short-contig / fragmented regime: route through the real metagenomic mode.
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
    let tinf = train_on_sequence(&pooled, trans_table, false)?;
    let cfg = GeneFinderConfig::default();

    let mut out = Vec::new();
    for c in contigs {
        let slen = c.seq.len() as i32;
        if slen < 1 {
            continue;
        }
        let result = find_genes(&c.seq, Some(&tinf), &cfg)?;
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
            &tinf,
            1,
            &c.name,
        );
        parse_protein_fasta(&buf, &c.name, &mut out);
    }
    Ok(out)
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

        let cds = predict_cds(&contigs, 11).expect("meta mode should not error");
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
        assert!(predict_cds(&[], 11).is_err());
        let empties = vec![
            Contig { name: "a".into(), seq: Vec::new() },
            Contig { name: "b".into(), seq: Vec::new() },
        ];
        assert!(predict_cds(&empties, 11).is_err());
    }
}
