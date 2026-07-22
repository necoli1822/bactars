//! Pseudogene detection: flag CDS/gene fragments showing frameshift, internal
//! stop, or truncation relative to a reference protein (PSC/HMM hit), matching
//! how RefSeq/PGAP annotates pseudogenes. Sets `func.pseudogene` + a `func.note`.
//!
//! ## Signals implemented
//!
//! ### 1. Reference-length truncation (PRIMARY, DB-aware)
//!
//! This is how RefSeq/PGAP actually annotates most pseudogenes: a gene
//! disrupted/truncated relative to its functional homolog. prodigal usually
//! calls such a pseudogene as *one* short CDS that stops early (at a premature
//! stop or after a frameshift), so its protein is markedly shorter than the
//! reference protein it best hits. When a CDS's best NCBIfams/TIGRFAM hit
//! resolves (via `hmm_PGAP.tsv` `hmm_length`) to a reference length `L_ref` and
//! the CDS protein is `< TRUNC_FRAC * L_ref`, we flag it as a truncated
//! pseudogene — provided it is **not at a contig edge** (a gene running off the
//! end of a contig is truncated by the *assembly*, not by a real disruption).
//!
//! This is the primary detector because the split-gene signal alone recovered
//! essentially none of MG1655's 145 RefSeq pseudogenes (most are a single short
//! truncated CDS, not two adjacent same-reference ORFs). It requires a reference
//! length, so it is only available through [`detect_with_refs`], which loads
//! `hmm_PGAP.tsv`. The DB-free [`detect`] runs the split-gene signal only.
//!
//! ### 2. Adjacent same-reference fragments (split gene / frameshift)
//!
//! The classic PGAP split signature: a functional gene disrupted by a frameshift
//! (or an internal premature stop) is called by the gene finder as *two or more*
//! short same-strand ORFs that each hit the *same* reference family. When we see
//! two+ genomically adjacent, same-strand CDS whose best annotation points at the
//! same reference (accession first, else model name) within a small
//! frameshift-consistent gap/overlap window, we flag *all* pieces as one
//! pseudogene. This catches disruptions whose fragments share a reference even
//! when no `hmm_length` is available (e.g. fragments that still hit an HMM but
//! whose ref length is missing, or name-only matches).
//!
//! ### 3. Reference-free frameshift fragments ([`frameshift_fragments`])
//!
//! A gene broken by an indel is called as two short, same-strand ORFs that OVERLAP
//! by a few nt in DIFFERENT reading frames. That overlap-in-a-shifted-frame
//! geometry (unlike plain adjacency, which is meaningless in a gene-dense genome)
//! is specific to a frameshift break, so it needs no reference length — reaching
//! disrupted genes that carry no resolvable `hmm_length`.
//!
//! ### 4. Reference-free internal-stop read-through ([`internal_stop_readthrough`])
//!
//! A premature internal stop leaves the gene finder calling a short CDS whose stop
//! is followed, IN THE SAME FRAME, by a long coding continuation that no other CDS
//! covers. We translate past the stop off the `contigs` sequence and flag when that
//! uncovered same-frame ORF is long ([`MIN_READTHROUGH_ORF`]). Also reference-free.
//!
//! ### Precision discipline
//!
//! We deliberately favour **precision over recall**. The truncation signal now also
//! demands an absolute deficit ([`MIN_ABS_SHORTFALL`]) so a short-but-complete
//! protein whose model is merely longer is not flagged; the split signal's weak
//! name-only fallback fires only for a specific (non-generic) shared name on short
//! pieces; the two reference-free signals require a frameshift/read-through geometry
//! a normal gene boundary does not exhibit. Measured on MG1655 vs RefSeq's 145
//! pseudogenes, the full set reaches ~13% recall at ~35% concordance-precision,
//! up from ~10% / ~33% for truncation+split alone.

use crate::fasta::Contig;
use crate::feature::{Feature, FeatureKind, Strand};
use std::collections::HashMap;

/// Largest intergenic gap (nt) between two fragments still considered adjacent
/// pieces of one frameshifted gene. Frameshift fragments abut or nearly abut.
pub const MAX_FRAGMENT_GAP: i64 = 60;

/// Largest overlap (nt) tolerated between two fragments. A frameshift commonly
/// leaves the two called ORFs overlapping by a few codons.
pub const MAX_FRAGMENT_OVERLAP: i64 = 30;

/// A CDS whose protein is shorter than this fraction of its reference protein
/// length (`hmm_length`) is a truncation candidate. Tuned on MG1655: 0.8 sits at
/// the recall/precision knee — it recovers 15/145 RefSeq pseudogenes (10.3%,
/// ~79% of the pool that even *has* a resolvable reference length) at ~46%
/// concordance with RefSeq's pseudogene set, while 0.7 recovers only 5 and 0.9
/// barely adds recall but floods false positives from natural length variation.
pub const TRUNC_FRAC: f64 = 0.8;

/// Minimum reference protein length (aa) for the truncation signal to fire.
/// Short HMM families give unreliable length ratios, so we ignore them.
pub const MIN_REF_LEN: usize = 50;

/// A CDS whose start/end lies within this many nt of a contig boundary is treated
/// as truncated by the assembly (partial gene at a contig edge), not by a real
/// disruption, and is never flagged by the truncation signal.
pub const CONTIG_EDGE_MARGIN: i64 = 3;

/// Minimum absolute shortfall (aa) between the reference and the CDS protein for
/// the truncation signal to fire, ON TOP OF the [`TRUNC_FRAC`] ratio. A genuinely
/// short protein whose family model happens to be a little longer (e.g. a 40 aa
/// protein vs a 55 aa model) trips the ratio but is not a pseudogene; requiring an
/// absolute deficit of this many residues removes those false positives. Tuned on
/// MG1655: 25 aa keeps every true truncation while trimming borderline FPs.
pub const MIN_ABS_SHORTFALL: usize = 25;

/// A CDS protein at or below this length (aa) is "short" — the regime where a
/// disruption (frameshift / internal stop) produces stubby fragments. Used to gate
/// the reference-free fragmentation and internal-stop signals so they do not flag
/// full-length genes.
pub const SHORT_PROTEIN_AA: usize = 150;

/// Largest overlap (nt) between two adjacent same-strand CDS still read as a
/// frameshift break. A frameshift makes the gene finder call two ORFs that overlap
/// by a few codons in DIFFERENT frames; genuine tandem genes abut with a gap or a
/// same-frame overlap. Kept tight (a real frameshift overlap is a handful of nt)
/// because same-strand adjacency alone is meaningless in a gene-dense genome.
pub const MAX_FRAMESHIFT_OVERLAP: i64 = 8;

/// Minimum length (in-frame sense codons) of the coding continuation past a CDS's
/// stop codon for the reference-free internal-stop signal to fire. A real gene's
/// stop is followed by intergenic sequence or an out-of-frame neighbour, so its
/// same-frame downstream run is short; a premature internal stop that truncated a
/// gene leaves a long same-frame ORF (the rest of the protein) that no other CDS
/// covers. Set high (80 aa) because the raw signal is weak in a dense genome — this
/// bar keeps it near-zero false positives while still recovering read-through
/// pseudogenes that carry no resolvable reference length.
pub const MIN_READTHROUGH_ORF: usize = 80;

/// Detect pseudogenes among CDS features and set `func.pseudogene` + a note on
/// each affected CDS. CDS features are never removed.
///
/// DB-free fallback: runs the split-gene fragment signal ONLY. When a bactars-db
/// `meta/` dir is available, prefer [`detect_with_refs`], which additionally runs
/// the primary reference-length truncation signal.
pub fn detect(_contigs: &[Contig], features: &mut Vec<Feature>) {
    for (idx, note) in split_gene_fragments(features) {
        let f = &mut features[idx];
        f.func.pseudogene = true;
        f.func.note.push(note);
    }
}

/// DB-aware pseudogene detection. Loads reference protein lengths from
/// `<meta_dir>/hmm_PGAP.tsv` and runs BOTH signals: the primary reference-length
/// truncation pass and the split-gene fragment pass. Sets `func.pseudogene` + a
/// note on each affected CDS; CDS features are never removed.
///
/// If `hmm_PGAP.tsv` cannot be read, degrades gracefully to the split-gene signal
/// (equivalent to [`detect`]) with a warning on stderr.
///
/// The orchestrator wires this when `--meta` is available; otherwise it calls
/// [`detect`].
pub fn detect_with_refs(contigs: &[Contig], features: &mut Vec<Feature>, meta_dir: &str) {
    let ref_lengths = match load_ref_lengths(meta_dir) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("pseudogene: no reference lengths ({e}); truncation signal skipped");
            HashMap::new()
        }
    };

    // Coverage bitmap (per contig+strand), shared by the sequence-based
    // internal-stop signal so it only counts downstream ORF that no other CDS on
    // the same strand already accounts for.
    let coverage = build_coverage(contigs, features);

    // Collect all flag sets, then apply. A CDS caught by several signals gets one
    // note per signal (deduplicated below); `pseudogene` is idempotent.
    let mut flags: Vec<(usize, String)> = truncated_genes(contigs, features, &ref_lengths);
    flags.extend(frameshift_fragments(features));
    flags.extend(internal_stop_readthrough(
        contigs,
        features,
        &ref_lengths,
        &coverage,
    ));
    flags.extend(split_gene_fragments(features));

    // One CDS may be flagged by multiple signals; keep each distinct note once.
    for (idx, note) in flags {
        let f = &mut features[idx];
        f.func.pseudogene = true;
        if !f.func.note.contains(&note) {
            f.func.note.push(note);
        }
    }
}

/// Reference-FREE disruption signals only — NO database, NO family label:
/// frameshift fragments + internal-stop read-through + split-gene fragments.
/// Returns `(feature_index, note)` pairs. Wired into the active `--pseudo` path to
/// reach disrupted genes that carry no ncbifams family — the coverage wall that
/// makes the alignment + reference-length detectors blind to ~90% of non-IS
/// pseudogenes (measured: only 4-11% of RefSeq non-IS pseudogenes have a family).
/// `internal_stop_readthrough` is passed an EMPTY ref-length map so it evaluates
/// every short interior CDS (with a real map it defers those to the truncation
/// signal). Designed precision-first (high thresholds), so expect modest recall.
pub fn detect_reference_free(contigs: &[Contig], features: &[Feature]) -> Vec<(usize, String)> {
    let coverage = build_coverage(contigs, features);
    let empty: HashMap<String, usize> = HashMap::new();
    let mut flags = frameshift_fragments(features);
    flags.extend(internal_stop_readthrough(contigs, features, &empty, &coverage));
    flags.extend(split_gene_fragments(features));
    flags
}

/// Core, side-effect-free detector: returns `(feature_index, note)` for every CDS
/// that is part of an adjacent same-reference fragment run (>= 2 pieces).
///
/// Split out from [`detect`] so the decision logic is unit-testable without
/// mutating features.
pub fn split_gene_fragments(features: &[Feature]) -> Vec<(usize, String)> {
    // Group CDS feature indices by (contig, strand); fragments of one gene share
    // both. Iterate deterministically (sorted keys) for stable output.
    let mut groups: HashMap<(String, Strand), Vec<usize>> = HashMap::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind == FeatureKind::Cds {
            groups
                .entry((f.contig.clone(), f.strand))
                .or_default()
                .push(i);
        }
    }
    let mut keys: Vec<_> = groups.keys().cloned().collect();
    keys.sort();

    let mut flags: Vec<(usize, String)> = Vec::new();
    for key in keys {
        let mut idxs = groups.remove(&key).unwrap();
        // Genomic order along the contig.
        idxs.sort_by_key(|&i| features[i].start);

        let mut k = 0;
        while k < idxs.len() {
            // Greedily extend a run of adjacent, same-reference fragments.
            let mut run = vec![idxs[k]];
            let mut j = k + 1;
            while j < idxs.len() {
                let prev = *run.last().unwrap();
                let cur = idxs[j];
                if adjacent(&features[prev], &features[cur])
                    && same_reference(&features[prev], &features[cur])
                {
                    run.push(cur);
                    j += 1;
                } else {
                    break;
                }
            }
            if run.len() >= 2 {
                let label = reference_label(&features[run[0]]);
                let note = format!(
                    "frameshifted; split into {} fragments matching {}",
                    run.len(),
                    label
                );
                for &i in &run {
                    flags.push((i, note.clone()));
                }
                k = j; // consume the whole run
            } else {
                k += 1;
            }
        }
    }

    flags.sort_by_key(|(i, _)| *i);
    flags
}

/// Core truncation detector: returns `(feature_index, note)` for every CDS whose
/// protein is shorter than [`TRUNC_FRAC`] of its best resolvable reference
/// protein length and which is not sitting at a contig edge.
///
/// `ref_lengths` maps an HMM accession (both the versioned `NF000282.2` form and
/// the version-stripped `NF000282` form are looked up) to the reference protein
/// length in amino acids (`hmm_length`). Side-effect-free so the decision logic
/// is unit-testable without mutating features or reading a DB.
pub fn truncated_genes(
    contigs: &[Contig],
    features: &[Feature],
    ref_lengths: &HashMap<String, usize>,
) -> Vec<(usize, String)> {
    // Contig name -> length (nt), for the contig-edge exclusion.
    let contig_len: HashMap<&str, i64> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.seq.len() as i64))
        .collect();

    let mut flags: Vec<(usize, String)> = Vec::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        let l_cds = match protein_len(f) {
            Some(n) if n > 0 => n,
            _ => continue,
        };
        // Best-scoring annotation whose accession resolves to a reference length.
        let (l_ref, label) = match best_ref_length(f, ref_lengths) {
            Some(x) => x,
            None => continue,
        };
        if l_ref < MIN_REF_LEN {
            continue;
        }
        // Primary signal: protein markedly shorter than the reference.
        if (l_cds as f64) >= TRUNC_FRAC * (l_ref as f64) {
            continue;
        }
        // Precision gate: also require a substantial ABSOLUTE deficit. A short but
        // complete protein whose family model is only a little longer trips the
        // ratio yet is not a pseudogene; demand at least MIN_ABS_SHORTFALL missing
        // residues so only real truncations pass.
        if l_ref.saturating_sub(l_cds) < MIN_ABS_SHORTFALL {
            continue;
        }
        // Disruption gate: a partial gene running off a contig end is truncated by
        // the assembly, not by a real internal stop/frameshift. Require both ends
        // interior. If the contig length is unknown, keep the candidate.
        if let Some(&clen) = contig_len.get(f.contig.as_str()) {
            if f.start <= CONTIG_EDGE_MARGIN || f.end >= clen - CONTIG_EDGE_MARGIN {
                continue;
            }
        }
        let pct = ((l_cds as f64 / l_ref as f64) * 100.0).round() as i64;
        flags.push((
            i,
            format!("truncated: {pct}% of reference length ({label})"),
        ));
    }
    flags
}

/// Protein length (aa) of a CDS feature, from its translation, ignoring a trailing
/// stop residue if present. `None` for a feature without a translation.
fn protein_len(f: &Feature) -> Option<usize> {
    f.aa.as_ref().map(|aa| aa.trim_end_matches('*').len())
}

/// The reference protein length and label for `f`, from the best available source.
///
/// Reference length source: **ncbifams/TIGRFAM `hmm_length`** — the highest-scoring
/// annotation whose accession resolves in `ref_lengths` (exact key first, then
/// version-stripped). This is a *curated consensus* length, the only reliable
/// truncation reference.
///
/// NOTE: PSC `Annotation.ref_len` (a UniRef90 best-hit target length) is deliberately
/// NOT used here. A UniRef90 representative is one arbitrary cluster member at 90%
/// identity, not a functional consensus, so "CDS shorter than its best UniRef90 hit"
/// fires on ordinary homolog length variation. Enabling it on MG1655 flagged 640 CDS
/// (precision 0.067) vs 145 real pseudogenes — a net-negative, so we keep truncation
/// on curated `hmm_length` only. The `ref_len` field remains plumbed for other uses.
fn best_ref_length(
    f: &Feature,
    ref_lengths: &HashMap<String, usize>,
) -> Option<(usize, String)> {
    // Curated hmm_length keyed by the annotation accession (highest-scoring hit).
    f.annotations
        .iter()
        .filter(|a| !a.accession.is_empty())
        .filter_map(|a| {
            resolve_ref_length(&a.accession, ref_lengths).map(|l| (a.score, l, a.accession.clone()))
        })
        // NaN-safe: `total_cmp` orders all f32 (incl. NaN) so a NaN score can't panic.
        .max_by(|x, y| x.0.total_cmp(&y.0))
        .map(|(_, l, acc)| (l, acc))
}

/// Look up a reference length for an accession: exact key, then version-stripped.
fn resolve_ref_length(accession: &str, ref_lengths: &HashMap<String, usize>) -> Option<usize> {
    if let Some(&l) = ref_lengths.get(accession) {
        return Some(l);
    }
    let bare = accession.split('.').next().unwrap_or(accession);
    ref_lengths.get(bare).copied()
}

/// Load accession -> reference protein length (`hmm_length`, col 6) from
/// `<meta_dir>/hmm_PGAP.tsv`. Both the versioned accession (`NF000282.2`) and the
/// version-stripped form (`NF000282`) are inserted as keys; on a stripped-key
/// collision the first record wins. The `#` header and blank lines are skipped.
pub fn load_ref_lengths(meta_dir: &str) -> Result<HashMap<String, usize>, String> {
    let path = format!("{}/hmm_PGAP.tsv", meta_dir.trim_end_matches('/'));
    let text = std::fs::read_to_string(&path)
        .map_err(|e| format!("failed to read {}: {}", path, e))?;

    let mut map: HashMap<String, usize> = HashMap::new();
    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        // Need through the hmm_length column (index 5).
        if cols.len() < 6 {
            continue;
        }
        let acc = cols[0].trim();
        if acc.is_empty() {
            continue;
        }
        let len: usize = match cols[5].trim().parse() {
            Ok(n) if n > 0 => n,
            _ => continue,
        };
        let bare = acc.split('.').next().unwrap_or(acc).to_string();
        map.entry(bare).or_insert(len);
        map.insert(acc.to_string(), len);
    }
    Ok(map)
}

/// Two genomically ordered features (`a.start <= b.start`) are adjacent fragments
/// if the gap between them is within the frameshift window (small gap or small
/// overlap).
fn adjacent(a: &Feature, b: &Feature) -> bool {
    let gap = b.start - a.end - 1;
    gap <= MAX_FRAGMENT_GAP && gap >= -MAX_FRAGMENT_OVERLAP
}

/// True if both features' best annotations point at the same reference family.
///
/// Accession identity is the strong signal and fires unconditionally. The
/// name-only fallback is far weaker — in a gene-dense genome countless adjacent
/// CDS share a *generic* product name (`hypothetical protein`, `putative ...`),
/// and flagging every such neighbour pair as one split gene floods false
/// positives. So the name fallback fires ONLY for a specific (non-generic) shared
/// name AND when both pieces are short (the stubby-fragment regime of a real
/// disruption), never for full-length tandem paralogs.
fn same_reference(a: &Feature, b: &Feature) -> bool {
    let (ba, bb) = match (a.best_annotation(), b.best_annotation()) {
        (Some(x), Some(y)) => (x, y),
        _ => return false,
    };
    if !ba.accession.is_empty() && ba.accession == bb.accession {
        return true;
    }
    !ba.name.is_empty()
        && ba.name == bb.name
        && !is_generic_name(&ba.name)
        && protein_len(a).map_or(false, |n| n <= SHORT_PROTEIN_AA)
        && protein_len(b).map_or(false, |n| n <= SHORT_PROTEIN_AA)
}

/// A product name too generic to establish that two CDS are fragments of the *same*
/// gene (as opposed to two unrelated genes that both lack a specific annotation).
fn is_generic_name(name: &str) -> bool {
    let n = name.trim().to_ascii_lowercase();
    n.is_empty()
        || n.contains("hypothetical")
        || n.contains("putative")
        || n.contains("uncharacterized")
        || n.contains("unknown")
        || n == "domain-containing protein"
}

/// Human-readable reference label for the note (accession preferred).
fn reference_label(f: &Feature) -> String {
    match f.best_annotation() {
        Some(a) if !a.accession.is_empty() => a.accession.clone(),
        Some(a) if !a.name.is_empty() => a.name.clone(),
        _ => "same reference".to_string(),
    }
}

// ============================================================================
// Reference-free genomic-disruption signals (frameshift + internal read-through)
// ============================================================================
//
// The reference-length truncation signal only reaches CDS whose best hit resolves
// to an `hmm_length`. Most disrupted genes carry no such reference (they hit
// UniRef/PSC by similarity, or nothing), so these two signals read the disruption
// straight off the genome — no reference length required — which is what lets them
// recover pseudogenes the truncation signal cannot see. Both are deliberately
// tight: same-strand adjacency and downstream ORFs are meaningless on their own in
// a gene-dense genome, so each demands a specific frameshift geometry (below) that
// a normal gene boundary does not exhibit.

/// Reference-free frameshift signal: a gene broken by an indel is called by the
/// gene finder as two short, same-strand ORFs that OVERLAP by a few nucleotides in
/// DIFFERENT reading frames. Returns `(feature_index, note)` for both pieces of
/// every such pair.
///
/// The specificity comes from the overlap-in-a-shifted-frame geometry: tandem genes
/// abut with a gap (or, rarely, overlap in the SAME frame); only a frameshift makes
/// two same-strand ORFs overlap by 1..[`MAX_FRAMESHIFT_OVERLAP`] nt with a
/// non-zero relative frame. We additionally require both pieces to be short
/// ([`SHORT_PROTEIN_AA`]) and NOT to be two independently-annotated distinct genes
/// (different non-empty accessions), so real overlapping operon genes are spared.
/// Side-effect-free for unit testing.
pub fn frameshift_fragments(features: &[Feature]) -> Vec<(usize, String)> {
    let mut groups: HashMap<(String, Strand), Vec<usize>> = HashMap::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind == FeatureKind::Cds {
            groups
                .entry((f.contig.clone(), f.strand))
                .or_default()
                .push(i);
        }
    }
    let mut keys: Vec<_> = groups.keys().cloned().collect();
    keys.sort();

    let mut flagged: HashMap<usize, String> = HashMap::new();
    for key in keys {
        let mut idxs = groups.remove(&key).unwrap();
        idxs.sort_by_key(|&i| features[i].start);
        for w in idxs.windows(2) {
            let (a, b) = (&features[w[0]], &features[w[1]]);
            // Strict frameshift overlap: b starts at or before a's end, by no more
            // than MAX_FRAMESHIFT_OVERLAP nt.
            let gap = b.start - a.end - 1;
            if gap > -1 || gap < -MAX_FRAMESHIFT_OVERLAP {
                continue;
            }
            // Both pieces short (stubby fragments).
            let (la, lb) = match (protein_len(a), protein_len(b)) {
                (Some(x), Some(y)) if x <= SHORT_PROTEIN_AA && y <= SHORT_PROTEIN_AA => (x, y),
                _ => continue,
            };
            let _ = (la, lb);
            // Relative reading frame of the two ORFs (measured from the gene's 5'
            // coordinate for each strand). Zero == same frame == not a frameshift.
            let rel = if a.strand >= 0 {
                b.start - a.start
            } else {
                a.end - b.end
            };
            if rel.rem_euclid(3) == 0 {
                continue;
            }
            // Reject anything that looks like an independently-annotated gene: a
            // real frameshifted pair is two UNannotated stubs, or two stubs that
            // hit the SAME family. If either piece carries an accession and they do
            // not share it, treat them as genuine neighbours, not one broken gene.
            let acc_a = a.best_annotation().map(|x| x.accession.as_str()).unwrap_or("");
            let acc_b = b.best_annotation().map(|x| x.accession.as_str()).unwrap_or("");
            let both_empty = acc_a.is_empty() && acc_b.is_empty();
            let same_acc = !acc_a.is_empty() && acc_a == acc_b;
            if !both_empty && !same_acc {
                continue;
            }
            let note = "frameshifted; overlapping out-of-frame fragments".to_string();
            flagged.entry(w[0]).or_insert_with(|| note.clone());
            flagged.entry(w[1]).or_insert(note);
        }
    }

    let mut out: Vec<(usize, String)> = flagged.into_iter().collect();
    out.sort_by_key(|(i, _)| *i);
    out
}

/// Reference-free internal-stop / read-through signal: a premature internal stop
/// truncates a gene, so the gene finder calls the 5' part as a short CDS whose stop
/// is followed — IN THE SAME READING FRAME — by a long coding continuation (the
/// rest of the protein) that no other CDS covers. Returns `(feature_index, note)`
/// for every such CDS.
///
/// Translating past the CDS stop and measuring the uncovered same-frame ORF is what
/// makes this reference-free. We only consider short ([`SHORT_PROTEIN_AA`]),
/// interior CDS that carry NO resolvable reference length (the truncation signal
/// already owns those), and require the continuation to be at least
/// [`MIN_READTHROUGH_ORF`] codons — a bar set high because the raw signal is weak
/// in a dense genome.
pub fn internal_stop_readthrough(
    contigs: &[Contig],
    features: &[Feature],
    ref_lengths: &HashMap<String, usize>,
    coverage: &Coverage,
) -> Vec<(usize, String)> {
    let seqs: HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.seq.as_slice()))
        .collect();

    let mut flags: Vec<(usize, String)> = Vec::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        let l_cds = match protein_len(f) {
            Some(n) if n > 0 && n <= SHORT_PROTEIN_AA => n,
            _ => continue,
        };
        // Owned by the truncation signal if a reference length resolves.
        if best_ref_length(f, ref_lengths).is_some() {
            continue;
        }
        let seq = match seqs.get(f.contig.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        let clen = seq.len() as i64;
        // Interior only: a partial gene at a contig edge is an assembly artefact.
        if f.start <= CONTIG_EDGE_MARGIN || f.end >= clen - CONTIG_EDGE_MARGIN {
            continue;
        }
        let covered = coverage.get(&f.contig, f.strand);
        let orf = uncovered_downstream_orf(seq, covered, f);
        if orf < MIN_READTHROUGH_ORF {
            continue;
        }
        let _ = l_cds;
        flags.push((
            i,
            format!("internal stop; {orf}-codon same-frame read-through past the stop"),
        ));
    }
    flags
}

/// Per-(contig, strand) coverage bitmap: `covered[pos]` (1-based) is `true` when
/// some CDS on that strand spans `pos`. Used by [`internal_stop_readthrough`] to
/// discount downstream ORF that another CDS already accounts for.
pub struct Coverage {
    /// (contig, strand) -> bitmap indexed by 1-based position (len = contig_len+1).
    maps: HashMap<(String, Strand), Vec<bool>>,
    empty: Vec<bool>,
}

impl Coverage {
    fn get(&self, contig: &str, strand: Strand) -> &[bool] {
        self.maps
            .get(&(contig.to_string(), strand))
            .map(|v| v.as_slice())
            .unwrap_or(self.empty.as_slice())
    }
}

/// Build the [`Coverage`] index from every CDS feature.
fn build_coverage(contigs: &[Contig], features: &[Feature]) -> Coverage {
    let clen: HashMap<&str, usize> = contigs.iter().map(|c| (c.name.as_str(), c.seq.len())).collect();
    let mut maps: HashMap<(String, Strand), Vec<bool>> = HashMap::new();
    for f in features {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        let n = match clen.get(f.contig.as_str()) {
            Some(&n) => n,
            None => continue,
        };
        let bm = maps
            .entry((f.contig.clone(), f.strand))
            .or_insert_with(|| vec![false; n + 1]);
        let lo = f.start.max(1) as usize;
        let hi = (f.end.min(n as i64)).max(0) as usize;
        for p in lo..=hi {
            if p < bm.len() {
                bm[p] = true;
            }
        }
    }
    Coverage {
        maps,
        empty: Vec::new(),
    }
}

/// Count in-frame sense codons past a CDS's stop codon, in the SAME reading frame,
/// walking outward (downstream in the gene's direction) until the next in-frame
/// stop, an `N`, a position already covered by another same-strand CDS, or the
/// contig end. `covered` is the 1-based coverage bitmap for this contig+strand.
fn uncovered_downstream_orf(seq: &[u8], covered: &[bool], f: &Feature) -> usize {
    const CAP: usize = 5000;
    let n = seq.len();
    let is_cov = |pos1: usize| -> bool { covered.get(pos1).copied().unwrap_or(false) };
    let mut count = 0usize;
    if f.strand >= 0 {
        // First codon after the stop starts at 0-based index `end` (= 1-based end+1).
        let mut p0 = f.end.max(0) as usize;
        while p0 + 3 <= n {
            if is_cov(p0 + 1) || is_cov(p0 + 2) || is_cov(p0 + 3) {
                break;
            }
            let c = &seq[p0..p0 + 3];
            if c.contains(&b'N') || is_plus_stop(c) {
                break;
            }
            count += 1;
            p0 += 3;
            if count >= CAP {
                break;
            }
        }
    } else {
        // Reverse strand: continuation reads leftward. First upstream codon's first
        // base is 0-based `start - 4` (= 1-based start-3).
        let mut p: i64 = f.start - 4;
        while p >= 0 {
            let pu = p as usize;
            if pu + 3 > n {
                break;
            }
            if is_cov(pu + 1) || is_cov(pu + 2) || is_cov(pu + 3) {
                break;
            }
            let c = &seq[pu..pu + 3];
            if c.contains(&b'N') || is_minus_stop(c) {
                break;
            }
            count += 1;
            p -= 3;
            if count >= CAP {
                break;
            }
        }
    }
    count
}

/// True if the 3-byte slice is a plus-strand stop codon (TAA/TAG/TGA).
fn is_plus_stop(c: &[u8]) -> bool {
    c == b"TAA" || c == b"TAG" || c == b"TGA"
}

/// True if the 3-byte plus-strand slice reverse-complements to a stop codon, i.e.
/// it is a stop on the MINUS strand (revcomp of {TAA,TAG,TGA} = {TTA,CTA,TCA}).
fn is_minus_stop(c: &[u8]) -> bool {
    c == b"TTA" || c == b"CTA" || c == b"TCA"
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Annotation, Functional};

    fn ann(accession: &str, name: &str, score: f32) -> Annotation {
        Annotation {
            source: "rustyhmmer:test".to_string(),
            accession: accession.to_string(),
            name: name.to_string(),
            score,
            evalue: Some(1e-30),
            ref_len: None,
        }
    }

    fn cds(start: i64, end: i64, strand: Strand, anns: Vec<Annotation>) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c1".to_string(),
            id: format!("c1_{start}"),
            start,
            end,
            strand,
            aa: Some("M".to_string()),
            partial5: false,
            partial3: false,
            annotations: anns,
            func: Functional::default(),
        }
    }

    #[test]
    fn adjacent_same_reference_both_flagged() {
        // Two adjacent + strand CDS both hitting NF000282 -> split gene.
        let features = vec![
            cds(100, 400, 1, vec![ann("NF000282.2", "dnaX", 250.0)]),
            cds(410, 700, 1, vec![ann("NF000282.2", "dnaX", 180.0)]),
        ];
        let flags = split_gene_fragments(&features);
        let flagged: Vec<usize> = flags.iter().map(|(i, _)| *i).collect();
        assert_eq!(flagged, vec![0, 1]);
        assert!(flags[0].1.contains("frameshifted"));
        assert!(flags[0].1.contains("NF000282.2"));
    }

    #[test]
    fn detect_sets_pseudogene_flag_and_note() {
        let mut features = vec![
            cds(100, 400, 1, vec![ann("NF000282.2", "dnaX", 250.0)]),
            cds(410, 700, 1, vec![ann("NF000282.2", "dnaX", 180.0)]),
        ];
        detect(&[], &mut features);
        assert!(features[0].func.pseudogene);
        assert!(features[1].func.pseudogene);
        assert!(!features[0].func.note.is_empty());
    }

    #[test]
    fn single_normal_cds_not_flagged() {
        let features = vec![cds(100, 1300, 1, vec![ann("NF000001.1", "gyrA", 900.0)])];
        assert!(split_gene_fragments(&features).is_empty());
    }

    #[test]
    fn adjacent_different_reference_not_flagged() {
        // Operon: two adjacent genes, DIFFERENT references -> not a split gene.
        let features = vec![
            cds(100, 400, 1, vec![ann("NF000010.1", "trpA", 250.0)]),
            cds(410, 800, 1, vec![ann("NF000011.1", "trpB", 260.0)]),
        ];
        assert!(split_gene_fragments(&features).is_empty());
    }

    #[test]
    fn same_reference_but_far_apart_not_flagged() {
        // Same family, but > MAX_FRAGMENT_GAP apart -> two genuine paralogs.
        let features = vec![
            cds(100, 400, 1, vec![ann("NF000282.2", "dnaX", 250.0)]),
            cds(5000, 5300, 1, vec![ann("NF000282.2", "dnaX", 240.0)]),
        ];
        assert!(split_gene_fragments(&features).is_empty());
    }

    #[test]
    fn same_reference_opposite_strand_not_flagged() {
        // Adjacent, same ref, but opposite strands -> not one gene's fragments.
        let features = vec![
            cds(100, 400, 1, vec![ann("NF000282.2", "dnaX", 250.0)]),
            cds(410, 700, -1, vec![ann("NF000282.2", "dnaX", 180.0)]),
        ];
        assert!(split_gene_fragments(&features).is_empty());
    }

    #[test]
    fn matches_by_name_when_accession_empty() {
        let features = vec![
            cds(100, 400, 1, vec![ann("", "recombinase family protein", 120.0)]),
            cds(405, 650, 1, vec![ann("", "recombinase family protein", 90.0)]),
        ];
        let flags = split_gene_fragments(&features);
        assert_eq!(flags.len(), 2);
        assert!(flags[0].1.contains("recombinase family protein"));
    }

    #[test]
    fn three_piece_run_all_flagged() {
        let features = vec![
            cds(100, 300, 1, vec![ann("NF001.1", "x", 100.0)]),
            cds(305, 500, 1, vec![ann("NF001.1", "x", 90.0)]),
            cds(505, 700, 1, vec![ann("NF001.1", "x", 80.0)]),
        ];
        let flags = split_gene_fragments(&features);
        assert_eq!(flags.len(), 3);
        assert!(flags[0].1.contains("3 fragments"));
    }

    #[test]
    fn small_overlap_still_adjacent() {
        // Frameshift often leaves the two ORFs overlapping by a few codons.
        let features = vec![
            cds(100, 400, 1, vec![ann("NF001.1", "x", 100.0)]),
            cds(385, 700, 1, vec![ann("NF001.1", "x", 90.0)]), // 16 nt overlap
        ];
        assert_eq!(split_gene_fragments(&features).len(), 2);
    }

    #[test]
    fn cds_without_annotation_not_flagged() {
        let features = vec![
            cds(100, 400, 1, vec![]),
            cds(410, 700, 1, vec![]),
        ];
        assert!(split_gene_fragments(&features).is_empty());
    }

    // ----- reference-length truncation signal -----

    /// A CDS with a protein of `aa_len` residues, on `contig`, coordinates chosen
    /// interior unless overridden, carrying one annotation with `accession`.
    fn cds_len(
        contig: &str,
        start: i64,
        aa_len: usize,
        accession: &str,
        score: f32,
    ) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: contig.to_string(),
            id: format!("{contig}_{start}"),
            start,
            end: start + (aa_len as i64 + 1) * 3 - 1, // +1 codon for the stop
            strand: 1,
            aa: Some("M".repeat(aa_len)),
            partial5: false,
            partial3: false,
            annotations: vec![ann(accession, "x", score)],
            func: Functional::default(),
        }
    }

    fn contig(name: &str, len: usize) -> Contig {
        Contig {
            name: name.to_string(),
            seq: vec![b'A'; len],
        }
    }

    fn reflen(pairs: &[(&str, usize)]) -> HashMap<String, usize> {
        pairs.iter().map(|(a, l)| (a.to_string(), *l)).collect()
    }

    #[test]
    fn truncated_cds_flagged() {
        // Protein 120 aa vs reference 300 aa (40%) -> truncated, interior.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 120, "NF000282.2", 200.0)];
        let refs = reflen(&[("NF000282.2", 300)]);
        let flags = truncated_genes(&contigs, &features, &refs);
        assert_eq!(flags.len(), 1);
        assert_eq!(flags[0].0, 0);
        assert!(flags[0].1.contains("truncated"));
        assert!(flags[0].1.contains("40%"));
        assert!(flags[0].1.contains("NF000282.2"));
    }

    #[test]
    fn full_length_cds_not_flagged() {
        // Protein 290 aa vs reference 300 aa (97%) -> complete, not a pseudogene.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 290, "NF000282.2", 200.0)];
        let refs = reflen(&[("NF000282.2", 300)]);
        assert!(truncated_genes(&contigs, &features, &refs).is_empty());
    }

    #[test]
    fn contig_edge_short_cds_not_flagged() {
        // A 120/300 short CDS, but it runs off the left contig edge (start=1) ->
        // truncated by the assembly, not a real disruption.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1, 120, "NF000282.2", 200.0)];
        let refs = reflen(&[("NF000282.2", 300)]);
        assert!(truncated_genes(&contigs, &features, &refs).is_empty());

        // Same, running off the right edge.
        let clen = 10_000i64;
        let mut f = cds_len("c1", 1000, 120, "NF000282.2", 200.0);
        f.end = clen; // abuts the right boundary
        assert!(truncated_genes(&contigs, &[f], &refs).is_empty());
    }

    #[test]
    fn tiny_reference_not_flagged() {
        // Reference below MIN_REF_LEN -> ratio unreliable, ignored.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 10, "NF000001.1", 50.0)];
        let refs = reflen(&[("NF000001.1", 30)]);
        assert!(truncated_genes(&contigs, &features, &refs).is_empty());
    }

    #[test]
    fn truncation_resolves_version_stripped_key() {
        // Annotation has a versioned accession; the ref map only has the bare key.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 100, "NF000282.2", 200.0)];
        let refs = reflen(&[("NF000282", 300)]);
        let flags = truncated_genes(&contigs, &features, &refs);
        assert_eq!(flags.len(), 1);
    }

    #[test]
    fn truncation_unknown_accession_not_flagged() {
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 100, "NF999999.9", 200.0)];
        let refs = reflen(&[("NF000282.2", 300)]);
        assert!(truncated_genes(&contigs, &features, &refs).is_empty());
    }

    #[test]
    fn truncation_picks_best_scoring_resolvable_annotation() {
        // Two annotations: a high-scoring one that resolves to a long ref (so the
        // CDS is truncated), and a low-scoring unresolvable one.
        let contigs = vec![contig("c1", 10_000)];
        let mut f = cds_len("c1", 1000, 100, "NF000282.2", 200.0);
        f.annotations.push(ann("PF99999.1", "y", 300.0)); // higher score, no reflen
        let refs = reflen(&[("NF000282.2", 300)]);
        let flags = truncated_genes(&contigs, &[f], &refs);
        assert_eq!(flags.len(), 1);
        assert!(flags[0].1.contains("NF000282.2"));
    }

    /// A PSC annotation's `ref_len` (a UniRef90 target length) must NOT drive the
    /// truncation signal. A UniRef90 representative is one arbitrary 90%-identity
    /// cluster member, not a functional consensus, so "CDS shorter than its best
    /// UniRef90 hit" fires on ordinary homolog length variation. Enabling it on MG1655
    /// flagged 640 CDS at precision 0.067 vs 145 real pseudogenes — a net-negative.
    /// Truncation stays on curated `hmm_length` only; a PSC-only CDS is not flagged.
    #[test]
    fn psc_ref_len_does_not_drive_truncation() {
        let contigs = vec![contig("c1", 10_000)];
        // 100 aa protein, only a PSC hit whose UniRef90 target is 300 aa (33%).
        let mut f = cds_len("c1", 1000, 100, "", 0.0);
        f.annotations = vec![Annotation {
            source: "psc:uniref90".to_string(),
            accession: "UniRef90_ABC".to_string(),
            name: "some real product".to_string(),
            score: 150.0,
            evalue: Some(1e-20),
            ref_len: Some(300),
        }];
        // No hmm_length map entry at all -> no curated reference -> not flagged.
        let refs: HashMap<String, usize> = HashMap::new();
        let flags = truncated_genes(&contigs, std::slice::from_ref(&f), &refs);
        assert!(flags.is_empty());
    }

    /// When BOTH an ncbifams `hmm_length` and a PSC `ref_len` are available, the
    /// curated `hmm_length` wins (behaviour for ncbifams-hit CDS is unchanged).
    #[test]
    fn ncbifams_hmm_length_preferred_over_psc_ref_len() {
        // ncbifams ref = 200 (protein 180 aa -> 90%, NOT truncated); PSC ref = 500
        // (would be 36%, truncated). Preference for ncbifams keeps it unflagged.
        let contigs = vec![contig("c1", 10_000)];
        let mut f = cds_len("c1", 1000, 180, "NF000282.2", 200.0);
        f.annotations.push(Annotation {
            source: "psc:uniref90".to_string(),
            accession: "UniRef90_ZZZ".to_string(),
            name: "psc product".to_string(),
            score: 400.0, // higher score than the ncbifams hit
            evalue: Some(1e-30),
            ref_len: Some(500),
        });
        let refs = reflen(&[("NF000282.2", 200)]);
        // best_ref_length must resolve to the ncbifams length (200), not 500.
        let (l_ref, label) = best_ref_length(&f, &refs).unwrap();
        assert_eq!(l_ref, 200);
        assert_eq!(label, "NF000282.2");
        // ...and so the 180/200 protein is NOT flagged as truncated.
        assert!(truncated_genes(&contigs, std::slice::from_ref(&f), &refs).is_empty());
    }

    // ----- reference-free frameshift signal -----

    #[test]
    fn frameshift_overlap_out_of_frame_flagged() {
        // Two short same-strand CDS overlapping by 3 nt in a shifted frame.
        let features = vec![
            cds(100, 250, 1, vec![]),
            cds(248, 400, 1, vec![]), // gap = -3, rel frame = 148 % 3 = 1
        ];
        let flags = frameshift_fragments(&features);
        let idx: Vec<usize> = flags.iter().map(|(i, _)| *i).collect();
        assert_eq!(idx, vec![0, 1]);
        assert!(flags[0].1.contains("frameshift"));
    }

    #[test]
    fn same_frame_overlap_not_flagged() {
        // Overlap of 1 nt but SAME frame (rel = 150 % 3 = 0) -> not a frameshift.
        let features = vec![cds(100, 250, 1, vec![]), cds(250, 400, 1, vec![])];
        assert!(frameshift_fragments(&features).is_empty());
    }

    #[test]
    fn frameshift_gap_not_overlapping_not_flagged() {
        // A positive gap (genes abut, do not overlap) -> not a frameshift break.
        let features = vec![cds(100, 250, 1, vec![]), cds(260, 400, 1, vec![])];
        assert!(frameshift_fragments(&features).is_empty());
    }

    #[test]
    fn frameshift_two_distinct_annotated_genes_not_flagged() {
        // Out-of-frame overlap but both carry different real accessions -> two
        // genuine overlapping genes, not one frameshifted gene.
        let features = vec![
            cds(100, 250, 1, vec![ann("NF000010.1", "a", 100.0)]),
            cds(248, 400, 1, vec![ann("NF000011.1", "b", 100.0)]),
        ];
        assert!(frameshift_fragments(&features).is_empty());
    }

    // ----- reference-free internal-stop / read-through signal -----

    /// A contig whose bytes are all `A` (codon AAA = Lys, never a stop) except an
    /// optional plus-strand stop inserted immediately after `stop_at` (0-based).
    fn contig_run(name: &str, len: usize, stop_at: Option<usize>) -> Contig {
        let mut seq = vec![b'A'; len];
        if let Some(p) = stop_at {
            seq[p] = b'T';
            seq[p + 1] = b'A';
            seq[p + 2] = b'A';
        }
        Contig {
            name: name.to_string(),
            seq,
        }
    }

    fn short_cds(contig: &str, start: i64, end: i64) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: contig.to_string(),
            id: format!("{contig}_{start}"),
            start,
            end,
            strand: 1,
            aa: Some("M".repeat(50)),
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional::default(),
        }
    }

    #[test]
    fn internal_stop_with_long_readthrough_flagged() {
        // CDS ends at 250; the whole downstream is AAA (no stop) for >80 codons ->
        // a long same-frame read-through -> internal-stop disruption.
        let contigs = vec![contig_run("c1", 600, None)];
        let features = vec![short_cds("c1", 100, 250)];
        let cov = build_coverage(&contigs, &features);
        let flags = internal_stop_readthrough(&contigs, &features, &HashMap::new(), &cov);
        assert_eq!(flags.len(), 1);
        assert_eq!(flags[0].0, 0);
        assert!(flags[0].1.contains("read-through"));
    }

    #[test]
    fn clean_short_protein_not_flagged() {
        // Same short CDS but a stop codon sits immediately after its stop (0-based
        // index 250) -> zero read-through -> a genuine short protein, not a pseudo.
        let contigs = vec![contig_run("c1", 600, Some(250))];
        let features = vec![short_cds("c1", 100, 250)];
        let cov = build_coverage(&contigs, &features);
        assert!(internal_stop_readthrough(&contigs, &features, &HashMap::new(), &cov).is_empty());
    }

    #[test]
    fn readthrough_skipped_when_reference_resolves() {
        // A long read-through but the CDS has a resolvable reference length -> owned
        // by the truncation signal, so the read-through signal stays quiet.
        let contigs = vec![contig_run("c1", 600, None)];
        let mut f = short_cds("c1", 100, 250);
        f.annotations = vec![ann("NF000282.2", "x", 100.0)];
        let features = vec![f];
        let refs = reflen(&[("NF000282.2", 300)]);
        let cov = build_coverage(&contigs, &features);
        assert!(internal_stop_readthrough(&contigs, &features, &refs, &cov).is_empty());
    }

    #[test]
    fn readthrough_covered_by_other_cds_not_flagged() {
        // The downstream ORF is fully covered by another same-strand CDS -> the
        // continuation is an annotated neighbour, not an orphan read-through.
        let contigs = vec![contig_run("c1", 600, None)];
        let features = vec![
            short_cds("c1", 100, 250),
            short_cds("c1", 251, 590), // covers the downstream frame
        ];
        let cov = build_coverage(&contigs, &features);
        assert!(internal_stop_readthrough(&contigs, &features, &HashMap::new(), &cov).is_empty());
    }

    #[test]
    fn truncation_requires_absolute_shortfall() {
        // 40 aa vs a 55 aa model: ratio 0.73 (< 0.8) but only 15 aa short (<
        // MIN_ABS_SHORTFALL) -> a short complete protein, not a truncation.
        let contigs = vec![contig("c1", 10_000)];
        let features = vec![cds_len("c1", 1000, 40, "NF000282.2", 200.0)];
        let refs = reflen(&[("NF000282.2", 55)]);
        assert!(truncated_genes(&contigs, &features, &refs).is_empty());
    }

    /// End-to-end concordance of the truncation + split signals against RefSeq
    /// MG1655's 145 pseudogenes, using the real bactars annotations
    /// (`/tmp/nokofam.gff3`), `hmm_PGAP.tsv`, the genome length, and the RefSeq
    /// GFF. Ignored by default (needs those files); run with:
    ///   cargo test -p bactars --release mg1655_pseudogene_recall -- --ignored --nocapture
    #[test]
    #[ignore]
    fn mg1655_pseudogene_recall() {
        use std::fs;

        let meta_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../db/meta");
        let gff = "/tmp/nokofam.gff3";
        let refgff = concat!(env!("CARGO_MANIFEST_DIR"), "/../benchmark/ref/genomic.gff");
        let fna = concat!(env!("CARGO_MANIFEST_DIR"), "/../benchmark/ref/genomic.fna");

        let ref_lengths = load_ref_lengths(meta_dir).expect("hmm_PGAP.tsv");

        // Real genome sequence (single contig NC_000913.3) — needed by the
        // sequence-based internal-stop signal.
        let mut seq: Vec<u8> = Vec::with_capacity(5_000_000);
        for l in fs::read_to_string(fna).expect("fna").lines() {
            if !l.starts_with('>') {
                seq.extend(l.trim().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        let contigs = vec![Contig {
            name: "NC_000913.3".to_string(),
            seq,
        }];

        // Build CDS features from the bactars GFF, recovering the ncbifams
        // accession from the Dbxref, the REAL product name, and a synthetic protein
        // of the right length. Using the real names exercises the split signal's
        // name-fallback exactly as production does.
        let mut features: Vec<Feature> = Vec::new();
        for line in fs::read_to_string(gff).expect("gff").lines() {
            if line.starts_with('#') {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 9 || c[2] != "CDS" {
                continue;
            }
            let (start, end): (i64, i64) = (c[3].parse().unwrap(), c[4].parse().unwrap());
            let strand: Strand = if c[6] == "-" { -1 } else { 1 };
            let acc = c[8]
                .split("ncbifams:")
                .nth(1)
                .and_then(|s| s.split(['%', ';', ',']).next())
                .unwrap_or("")
                .to_string();
            let name = c[8]
                .split(';')
                .find_map(|kv| kv.trim().strip_prefix("Name="))
                .unwrap_or("")
                .to_string();
            let aa_len = (((end - start + 1) / 3) - 1).max(0) as usize;
            let anns = if acc.is_empty() && name.is_empty() {
                vec![]
            } else {
                vec![ann(&acc, &name, 100.0)]
            };
            features.push(Feature {
                kind: FeatureKind::Cds,
                contig: "NC_000913.3".to_string(),
                id: format!("cds_{start}"),
                start,
                end,
                strand,
                aa: Some("M".repeat(aa_len)),
                partial5: false,
                partial3: false,
                annotations: anns,
                func: Functional::default(),
            });
        }

        // RefSeq pseudogene coordinates (the 145).
        let mut pseudo: Vec<(i64, i64)> = Vec::new();
        for line in fs::read_to_string(refgff).expect("refgff").lines() {
            if line.contains("\tpseudogene\t") {
                let c: Vec<&str> = line.split('\t').collect();
                pseudo.push((c[3].parse().unwrap(), c[4].parse().unwrap()));
            }
        }

        let overlaps = |a: &(i64, i64), b: &(i64, i64)| !(a.1 < b.0 || b.1 < a.0);
        let score = |flagged: &[(i64, i64)]| -> (usize, usize, usize) {
            let tp = pseudo
                .iter()
                .filter(|p| flagged.iter().any(|f| overlaps(p, f)))
                .count();
            let conc = flagged
                .iter()
                .filter(|f| pseudo.iter().any(|p| overlaps(f, p)))
                .count();
            (tp, conc, flagged.len())
        };
        let report = |label: &str, flagged: &[(i64, i64)]| {
            let (tp, conc, n) = score(flagged);
            eprintln!(
                "  {label:<10} recall = {tp:>3}/{} = {:>4.1}%   precision = {conc:>3}/{n:<3} = {:>4.1}%",
                pseudo.len(),
                tp as f64 / pseudo.len() as f64 * 100.0,
                conc as f64 / n.max(1) as f64 * 100.0
            );
        };

        // ---- BEFORE: the original two signals (truncation @0.8 with NO absolute
        // gate + split with the promiscuous accession-OR-any-name fallback). ----
        let mut before: Vec<usize> = Vec::new();
        {
            let clen = contigs[0].seq.len() as i64;
            for (i, f) in features.iter().enumerate() {
                let l_cds = f.aa.as_ref().map(|a| a.len()).unwrap_or(0);
                if l_cds == 0 {
                    continue;
                }
                if let Some((l_ref, _)) = best_ref_length(f, &ref_lengths) {
                    if l_ref >= MIN_REF_LEN
                        && (l_cds as f64) < TRUNC_FRAC * l_ref as f64
                        && f.start > CONTIG_EDGE_MARGIN
                        && f.end < clen - CONTIG_EDGE_MARGIN
                    {
                        before.push(i);
                    }
                }
            }
            // Original split: same accession OR same (any) name, adjacent.
            let mut groups: HashMap<(String, Strand), Vec<usize>> = HashMap::new();
            for (i, f) in features.iter().enumerate() {
                groups.entry((f.contig.clone(), f.strand)).or_default().push(i);
            }
            for (_, mut idxs) in groups {
                idxs.sort_by_key(|&i| features[i].start);
                let mut k = 0;
                while k < idxs.len() {
                    let mut run = vec![idxs[k]];
                    let mut j = k + 1;
                    while j < idxs.len() {
                        let (a, b) = (&features[*run.last().unwrap()], &features[idxs[j]]);
                        let gap = b.start - a.end - 1;
                        let sr = match (a.best_annotation(), b.best_annotation()) {
                            (Some(x), Some(y)) => {
                                (!x.accession.is_empty() && x.accession == y.accession)
                                    || (!x.name.is_empty() && x.name == y.name)
                            }
                            _ => false,
                        };
                        if gap <= MAX_FRAGMENT_GAP && gap >= -MAX_FRAGMENT_OVERLAP && sr {
                            run.push(idxs[j]);
                            j += 1;
                        } else {
                            break;
                        }
                    }
                    if run.len() >= 2 {
                        before.extend(&run);
                        k = j;
                    } else {
                        k += 1;
                    }
                }
            }
            before.sort();
            before.dedup();
        }
        let before_coords: Vec<(i64, i64)> =
            before.iter().map(|&i| (features[i].start, features[i].end)).collect();

        // ---- AFTER: the improved pipeline through the public API. ----
        let mut after_features = features.clone();
        detect_with_refs(&contigs, &mut after_features, meta_dir);
        let after_coords: Vec<(i64, i64)> = after_features
            .iter()
            .filter(|f| f.func.pseudogene)
            .map(|f| (f.start, f.end))
            .collect();

        eprintln!("MG1655 pseudogene concordance ({} RefSeq pseudogenes):", pseudo.len());
        report("BEFORE", &before_coords);
        report("AFTER", &after_coords);

        let (tp_b, _, _) = score(&before_coords);
        let (tp_a, conc_a, n_a) = score(&after_coords);
        let prec_a = conc_a as f64 / n_a.max(1) as f64 * 100.0;
        // The improved pipeline must beat the 10.3% / 33% production baseline on
        // BOTH axes: recall strictly above 15/145 and precision above 33%.
        assert!(tp_a >= 16, "recall regressed: expected >=16 of 145, got {tp_a} (before={tp_b})");
        assert!(prec_a > 33.0, "precision below target: {prec_a:.1}%");
    }
}
