//! Assembly gap (run-of-N) annotation — DB-free, pure computation.
//!
//! Finds maximal runs of `N`/`n` of length >= MIN_GAP and emits
//! `Feature { kind: FeatureKind::AssemblyGap }` (→ `assembly_gap`).
//! Standard for draft assemblies.
//!
//! Note: only `N`/`n` counts as a gap residue — other IUPAC ambiguity codes
//! (R, Y, W, S, ...) are NOT treated as assembly gaps. `MIN_GAP` is set to 1
//! so no run is silently dropped here; callers that want to filter out short
//! runs (e.g. only report gaps >= 10 nt) can do so downstream.

use crate::fasta::Contig;
use crate::feature::{Feature, FeatureKind, Functional};

/// Minimum run length (nt) to call an assembly gap.
pub const MIN_GAP: usize = 1;

/// Detect assembly gaps (N-runs) across all contigs.
pub fn detect(contigs: &[Contig]) -> Vec<Feature> {
    let mut out = Vec::new();
    for contig in contigs {
        let mut idx1based = 0usize;
        let mut i = 0usize;
        let n = contig.seq.len();
        while i < n {
            if is_gap_residue(contig.seq[i]) {
                let start0 = i;
                let mut j = i + 1;
                while j < n && is_gap_residue(contig.seq[j]) {
                    j += 1;
                }
                let run_len = j - start0;
                if run_len >= MIN_GAP {
                    idx1based += 1;
                    out.push(Feature {
                        kind: FeatureKind::AssemblyGap,
                        contig: contig.name.clone(),
                        id: format!("{}_gap{}", contig.name, idx1based),
                        start: (start0 + 1) as i64,
                        end: j as i64,
                        strand: 1,
                        aa: None,
                        partial5: false,
                        partial3: false,
                        annotations: Vec::new(),
                        func: Functional::default(),
                    });
                }
                i = j;
            } else {
                i += 1;
            }
        }
    }
    out
}

#[inline]
fn is_gap_residue(b: u8) -> bool {
    matches!(b, b'N' | b'n')
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

    #[test]
    fn no_n_yields_no_gaps() {
        let contigs = vec![contig("c1", "ACGTACGTACGT")];
        let feats = detect(&contigs);
        assert!(feats.is_empty());
    }

    #[test]
    fn single_internal_run() {
        // 0-based: A C G N N N T -> N-run at indices 3..6 -> 1-based 4..6
        let contigs = vec![contig("c1", "ACGNNNT")];
        let feats = detect(&contigs);
        assert_eq!(feats.len(), 1);
        let f = &feats[0];
        assert_eq!(f.kind, FeatureKind::AssemblyGap);
        assert_eq!(f.contig, "c1");
        assert_eq!(f.id, "c1_gap1");
        assert_eq!(f.start, 4);
        assert_eq!(f.end, 6);
        assert_eq!(f.strand, 1);
        assert!(f.aa.is_none());
        assert!(f.annotations.is_empty());
    }

    #[test]
    fn two_separate_runs() {
        // indices: 0:A 1:N 2:N 3:C 4:G 5:N 6:T
        // run1: 1..3 (0-based) -> 1-based 2..3
        // run2: 5..6 (0-based) -> 1-based 6..6
        let contigs = vec![contig("c2", "ANNCGNT")];
        let feats = detect(&contigs);
        assert_eq!(feats.len(), 2);

        assert_eq!(feats[0].id, "c2_gap1");
        assert_eq!(feats[0].start, 2);
        assert_eq!(feats[0].end, 3);

        assert_eq!(feats[1].id, "c2_gap2");
        assert_eq!(feats[1].start, 6);
        assert_eq!(feats[1].end, 6);
    }

    #[test]
    fn run_at_start_and_end() {
        // NN ACGT NNN  (len 9): start run 0..2 (1-based 1..2), end run 6..9 (1-based 7..9)
        let contigs = vec![contig("c3", "NNACGTNNN")];
        let feats = detect(&contigs);
        assert_eq!(feats.len(), 2);

        assert_eq!(feats[0].start, 1);
        assert_eq!(feats[0].end, 2);

        assert_eq!(feats[1].start, 7);
        assert_eq!(feats[1].end, 9);
    }

    #[test]
    fn lowercase_n_is_detected() {
        let contigs = vec![contig("c4", "acgnnnacgt")];
        let feats = detect(&contigs);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].start, 4);
        assert_eq!(feats[0].end, 6);
    }

    #[test]
    fn multiple_contigs_reset_index() {
        let contigs = vec![contig("a", "NN"), contig("b", "AANN")];
        let feats = detect(&contigs);
        assert_eq!(feats.len(), 2);
        assert_eq!(feats[0].id, "a_gap1");
        assert_eq!(feats[0].contig, "a");
        assert_eq!(feats[1].id, "b_gap1");
        assert_eq!(feats[1].contig, "b");
    }
}
