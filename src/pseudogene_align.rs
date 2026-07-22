//! Alignment-based pseudogene detection (EXPERIMENTAL, off by default).
//!
//! PGAP-style disruption detection via an **mmseqs translated search**. Instead of
//! the length/geometry proxies in [`crate::pseudogene`], each candidate CDS's
//! genomic region (downstream-extended, kept as nucleotide) is searched — in all
//! six frames — against a small **consensus-protein reference** derived from the
//! ncbifams HMMs, and a disruption is read off the hit pattern:
//!
//!   * **internal premature stop** — the SAME reference family is hit by two query
//!     alignments that cover consecutive target ranges (N-terminal + C-terminal)
//!     IN THE SAME reading frame, with the query broken between them: a stop codon
//!     split one ORF into two in-frame pieces.
//!   * **frameshift** — the same split pattern but the two pieces are in
//!     DIFFERENT frames/strand: an indel shifted the frame mid-gene.
//!
//! Both are the classic PGAP pseudogene signatures. The two-piece requirement is
//! what makes this alignment-specific and high-precision — a naturally short but
//! intact family member yields a single full-or-prefix hit, never a split.
//!
//! ### Reference DB (LIGHT-tier)
//! The reference is a **majority-rule consensus protein per ncbifams family**
//! (argmax match-state emission per node), materialised once into a ~6 MB, ~19k
//! sequence FASTA by [`build_consensus_fasta`] streaming `ncbifams.hmm`. This is
//! derived only from ncbifams (already shipped in the LIGHT tier) — it does NOT
//! need the 194 GB UniRef90/PSC DB, so pseudogene detection stays available on
//! both LIGHT and FULL tiers.
//!
//! ### Why mmseqs
//! mmseqs is bactars' one sanctioned external binary (zero new dependency), gives
//! a k-mer prefilter + calibrated E-values + coverage, and auto-translates the
//! nucleotide query in six frames vs the protein reference — exactly the
//! blastx-style search this needs. Measured cost: 6–13 s per genome (16 threads,
//! <1 GB RAM) on top of the shared annotation run.
//!
//! Wired into the pipeline behind an off-by-default flag; the proxy baseline in
//! [`crate::pseudogene`] is untouched when the flag is off.

use crate::fasta::Contig;
use crate::feature::{Feature, FeatureKind};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// HMMER match-emission column order (amino, K=20).
const SYM: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// Nucleotides prepended upstream of each candidate CDS in the query region.
pub const REGION_UP_NT: i64 = 60;
/// Nucleotides appended downstream of each candidate CDS's stop, so a disrupted
/// gene's "missing" C-terminus is present for the aligner to recover.
pub const REGION_DOWN_NT: i64 = 450;
/// A CDS within this many nt of a contig boundary is a partial gene truncated by
/// the assembly, not a real disruption; never flagged.
pub const CONTIG_EDGE_MARGIN: i64 = 3;

/// A hit must be at least this significant (E-value) to establish family membership.
pub const EVALUE_MAX: f64 = 1e-10;
/// Each split piece must cover at least this many target (aa) residues.
pub const MIN_PIECE_AA: i64 = 15;
/// The two split pieces together must cover at least this fraction of the target
/// (lower bound: the recovered gene is nearly whole).
pub const MIN_COMBINED_TCOV: f64 = 0.50;
/// The SUM of the two pieces' aligned target lengths may cover at most this
/// fraction of the model. A genuine break splits ONE model between the pieces
/// (sum ~= 1x). Two pieces summing to ~2x the model are two FULL copies — a
/// tandem DUPLICATION of a complete gene (two real genes), the classic false
/// positive — and are rejected. (Geometry rule 2, model side.)
pub const MAX_COMBINED_TCOV: f64 = 1.25;
/// The two aligned pieces must lie within this many nt of each other in the
/// genome. A real disruption (a stop codon or a frameshifting indel) leaves the
/// fragments abutting; a large gap means the "pieces" are two separate genes.
/// (Geometry rule 1, raw gap.)
pub const MAX_Q_GAP_NT: i64 = 300;
/// The two pieces' combined genomic extent may span at most this multiple of one
/// model length (in nt; model_aa * 3). One broken gene occupies ~1 model length
/// of genome; two adjacent same-family genes read as a "split" occupy ~2x and are
/// rejected. This is the strongest discriminator against tandem paralogs.
/// (Geometry rule 1+2, genomic extent.)
pub const QSPAN_MODEL_FRAC_MAX: f64 = 1.5;
/// Target ranges of the two pieces may overlap by at most this many aa.
pub const T_OVERLAP_TOL: i64 = 20;
/// The C-terminal piece may start at most this many aa past the N-terminal piece's
/// target end (allows a small unaligned gap at the disruption point).
pub const T_GAP_TOL: i64 = 60;
/// Ignore reference families shorter than this (aa).
pub const MIN_TLEN: i64 = 60;

// --- truncation (single-piece) thresholds (B1/B0) ---
/// A single-piece hit covering at most this fraction of the model is a truncation
/// candidate (the gene is a short fragment of a longer family). Non-transposase.
pub const TRUNC_TCOV_MAX: f64 = 0.70;
/// Relaxed truncation ceiling for transposases / IS elements (B0): RefSeq/PGAP
/// pseudogenizes degraded IS copies liberally, so a milder deficit still counts.
pub const TRUNC_TCOV_MAX_IS: f64 = 0.85;
/// A truncation candidate must still be a real family member — cover at least this
/// fraction of the model (below this it is a spurious short hit, not a truncation).
pub const TRUNC_TCOV_MIN: f64 = 0.20;
/// One model end must be intact within this many aa (a real truncation keeps either
/// the N- or the C-terminus; a hit missing BOTH ends is an internal fragment/artefact).
pub const TRUNC_END_TOL: i64 = 10;

/// Stream `ncbifams.hmm` once and write one majority-rule consensus protein per
/// family to `out_path` (FASTA, id = accession). Returns the number of records.
///
/// Consensus residue per node = argmin match-emission score (HMMER stores
/// emissions as negative log-probabilities; the smallest score = highest prob).
/// Skips work and returns 0 if `out_path` already exists and is non-empty (cache).
pub fn build_consensus_fasta(hmm_path: &str, out_path: &str) -> Result<usize, String> {
    if let Ok(m) = std::fs::metadata(out_path) {
        if m.len() > 0 {
            return Ok(0); // cached
        }
    }
    let file = File::open(hmm_path).map_err(|e| format!("open {hmm_path}: {e}"))?;
    let reader = BufReader::new(file);
    let tmp = format!("{out_path}.tmp{}", std::process::id());
    let mut out = BufWriter::new(File::create(&tmp).map_err(|e| format!("create {tmp}: {e}"))?);

    let mut acc: Option<String> = None;
    let mut name = String::new();
    let mut cons = String::new();
    let mut in_body = false;
    let mut n = 0usize;

    for line in reader.lines() {
        let line = line.map_err(|e| format!("read {hmm_path}: {e}"))?;
        if let Some(r) = line.strip_prefix("ACC ") {
            acc = Some(r.trim().to_string());
        } else if let Some(r) = line.strip_prefix("NAME ") {
            name = r.trim().to_string();
        } else if line.starts_with("HMM ") {
            in_body = true;
            cons.clear();
        } else if line.starts_with("//") {
            if let Some(a) = acc.take() {
                if !cons.is_empty() {
                    writeln!(out, ">{a} {name}").map_err(|e| e.to_string())?;
                    for chunk in cons.as_bytes().chunks(60) {
                        out.write_all(chunk).map_err(|e| e.to_string())?;
                        out.write_all(b"\n").map_err(|e| e.to_string())?;
                    }
                    n += 1;
                }
            }
            name.clear();
            cons.clear();
            in_body = false;
        } else if in_body {
            let mut it = line.split_whitespace();
            let first = match it.next() {
                Some(t) => t,
                None => continue,
            };
            if first.parse::<usize>().is_err() {
                continue;
            }
            let mut best_i = usize::MAX;
            let mut best_v = f64::INFINITY;
            let mut k = 0usize;
            for tok in it.by_ref().take(20) {
                if tok != "*" {
                    if let Ok(v) = tok.parse::<f64>() {
                        if v < best_v {
                            best_v = v;
                            best_i = k;
                        }
                    }
                }
                k += 1;
            }
            if k == 20 && best_i < 20 {
                cons.push(SYM[best_i] as char);
            }
        }
    }
    out.flush().map_err(|e| e.to_string())?;
    drop(out);
    std::fs::rename(&tmp, out_path).map_err(|e| format!("rename {tmp}: {e}"))?;
    Ok(n)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

/// Write the candidate CDS regions (downstream-extended, strand-oriented
/// nucleotide) to a FASTA. The record id encodes the feature so a hit maps back:
/// `idx|start|end|strand|edge`. `idx` is the feature index in `features`.
fn write_regions(
    contigs: &[Contig],
    features: &[Feature],
    out_path: &str,
) -> Result<usize, String> {
    let seqs: HashMap<&str, &[u8]> = contigs
        .iter()
        .map(|c| (c.name.as_str(), c.seq.as_slice()))
        .collect();
    let mut out =
        BufWriter::new(File::create(out_path).map_err(|e| format!("create {out_path}: {e}"))?);
    let mut n = 0usize;
    for (i, f) in features.iter().enumerate() {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        let seq = match seqs.get(f.contig.as_str()) {
            Some(s) => *s,
            None => continue,
        };
        let clen = seq.len() as i64;
        let edge = if f.start <= CONTIG_EDGE_MARGIN || f.end >= clen - CONTIG_EDGE_MARGIN {
            1
        } else {
            0
        };
        let (s, e) = (f.start - 1, f.end - 1); // 0-based inclusive
        let sub: Vec<u8> = if f.strand >= 0 {
            let a = (s - REGION_UP_NT).max(0);
            let b = (e + REGION_DOWN_NT).min(clen - 1);
            seq[a as usize..=b as usize].to_vec()
        } else {
            let a = (s - REGION_DOWN_NT).max(0);
            let b = (e + REGION_UP_NT).min(clen - 1);
            revcomp(&seq[a as usize..=b as usize])
        };
        if sub.len() < 30 {
            continue;
        }
        writeln!(out, ">{i}|{}|{}|{}|{edge}", f.start, f.end, f.strand)
            .map_err(|e| e.to_string())?;
        for chunk in sub.chunks(70) {
            out.write_all(chunk).map_err(|e| e.to_string())?;
            out.write_all(b"\n").map_err(|e| e.to_string())?;
        }
        n += 1;
    }
    out.flush().map_err(|e| e.to_string())?;
    Ok(n)
}

/// One parsed m8 hit (translated search): nucleotide query coords, aa target coords.
#[derive(Clone, Copy)]
struct Hit {
    qstart: i64,
    qend: i64,
    tstart: i64,
    tend: i64,
    tlen: i64,
    evalue: f64,
}

/// (strand, phase) of an alignment from its nucleotide query coordinates.
fn frame_of(qs: i64, qe: i64) -> (i8, i64) {
    let strand = if qe >= qs { 1i8 } else { -1i8 };
    let lo = qs.min(qe);
    (strand, (lo - 1).rem_euclid(3))
}

/// Given all significant hits of one query to its PRIMARY family (`prim` target's
/// hits), decide whether they form a disruption split, returning the class note.
fn classify_split(hits: &[Hit], tlen: i64) -> Option<&'static str> {
    for a in hits {
        for b in hits {
            // a = N-terminal piece, b = C-terminal piece (in target/model coords).
            if !(a.tstart <= b.tstart && a.tend <= b.tend) {
                continue;
            }
            let alen = a.tend - a.tstart + 1;
            let blen = b.tend - b.tstart + 1;
            if alen < MIN_PIECE_AA || blen < MIN_PIECE_AA {
                continue;
            }
            // RULE 3: the pieces must cover CONSECUTIVE, largely NON-overlapping
            // model ranges (piece a ~ N-term, piece b ~ C-term). Two pieces over
            // the SAME model region are two paralogs, not one broken gene.
            let t_overlap = a.tend - b.tstart + 1; // >0 => model ranges overlap (aa)
            if t_overlap > T_OVERLAP_TOL {
                continue;
            }
            let t_gap = b.tstart - a.tend - 1; // unaligned model between pieces (aa)
            if t_gap > T_GAP_TOL {
                continue;
            }
            // RULE 2 (model side): combined aligned model length ~= ONE model.
            // Lower bound on the union span: the recovered gene is nearly whole.
            let union_cov = b.tend - a.tstart + 1;
            if (union_cov as f64) < MIN_COMBINED_TCOV * tlen as f64 {
                continue;
            }
            // Upper bound on the raw sum: reject ~2x = two full copies (tandem dup).
            let sum_cov = alen + blen;
            if (sum_cov as f64) > MAX_COMBINED_TCOV * tlen as f64 {
                continue;
            }
            // b must be a distinct query segment DOWNSTREAM of a (5'->3').
            let (a_lo, a_hi) = (a.qstart.min(a.qend), a.qstart.max(a.qend));
            let (b_lo, b_hi) = (b.qstart.min(b.qend), b.qstart.max(b.qend));
            if b_lo <= a_lo {
                continue;
            }
            // RULE 4: same contig + same strand. (The region is already strand-
            // oriented and single-contig per query; here we require both aligned
            // pieces to read off the SAME strand — a disruption keeps them co-linear.)
            let (sa, pa) = frame_of(a.qstart, a.qend);
            let (sb, pb) = frame_of(b.qstart, b.qend);
            if sa != sb {
                continue;
            }
            // RULE 1: the two pieces must be genomically CLOSE. A real disruption
            // leaves the fragments abutting; a whole gene's worth of sequence
            // between them means these are two separate genes (tandem paralogs).
            let q_gap = b_lo - a_hi - 1; // nt between the aligned pieces
            if q_gap > MAX_Q_GAP_NT {
                continue;
            }
            // ...and their combined genomic extent must span ~ONE model, not two.
            let q_span = b_hi - a_lo + 1; // total genomic extent of the pair (nt)
            if (q_span as f64) > QSPAN_MODEL_FRAC_MAX * (tlen as f64) * 3.0 {
                continue;
            }
            // Same strand already enforced; the split is a frameshift iff the two
            // pieces read in different phases, otherwise an internal premature stop.
            return Some(if pa != pb {
                "frameshift"
            } else {
                "internal_stop"
            });
        }
    }
    None
}

/// Per-CDS annotation context threaded into detection for the truncation signals.
#[derive(Clone, Copy, Default)]
struct FeatInfo {
    /// Protein length (aa) of the called CDS.
    aa_len: i64,
    /// Product/name marks this as a transposase / IS element (B0 relaxation).
    is_transpos: bool,
}

/// Single-piece TRUNCATION signal (B1/B0): the best hit to the primary family
/// covers only a fraction of the model while keeping ONE terminus intact, and the
/// called gene is itself short — a truncated fragment of a longer family. Runs only
/// when no split was found (splits are the stronger, already-handled signal).
fn classify_truncation(hits: &[Hit], tlen: i64, info: &FeatInfo) -> Option<String> {
    // Truncation-by-coverage is only reliable for TRANSPOSASES / IS elements (B0):
    // RefSeq pseudogenizes degraded IS copies liberally, and a truncated transposase
    // is almost always a dead IS. For general genes the signal is near-noise — the
    // per-family consensus is often longer than an individual member, so a normal short
    // gene covers only part of it (measured precision ~1-6%). So gate on is_transpos.
    if !info.is_transpos {
        return None;
    }
    // best (largest target span) hit to the primary family
    let best = hits
        .iter()
        .max_by_key(|h| h.tend - h.tstart)
        .copied()?;
    let cov_aa = best.tend - best.tstart + 1;
    if cov_aa < MIN_PIECE_AA {
        return None;
    }
    let tcov = cov_aa as f64 / tlen as f64;
    let ceil = if info.is_transpos {
        TRUNC_TCOV_MAX_IS
    } else {
        TRUNC_TCOV_MAX
    };
    if !(TRUNC_TCOV_MIN..=ceil).contains(&tcov) {
        return None;
    }
    // one model terminus must be intact (N-term near 1, or C-term near tlen).
    let n_intact = best.tstart <= 1 + TRUNC_END_TOL;
    let c_intact = best.tend >= tlen - TRUNC_END_TOL;
    if !(n_intact || c_intact) {
        return None;
    }
    // the CALLED gene must itself be short relative to the model (it is the truncated
    // fragment, not a full gene that merely aligns to part of a longer model).
    if info.aa_len >= (ceil * tlen as f64) as i64 {
        return None;
    }
    let end = if n_intact && !c_intact {
        "C-terminal"
    } else if c_intact && !n_intact {
        "N-terminal"
    } else {
        "internal"
    };
    Some(format!(
        "truncated ({end} loss, {:.0}% of family, aligned)",
        tcov * 100.0
    ))
}

/// Parse the mmseqs m8 and return `(feature_index, note)` for each disrupted CDS.
/// `info` supplies per-CDS context (protein length, transposase flag) for the
/// truncation signals; keyed by feature index.
///
/// m8 columns (as produced by [`detect_mmseqs`]):
/// `query target qstart qend tstart tend qlen tlen alnlen pident evalue bits`.
fn parse_and_detect(
    m8_path: &str,
    info: &HashMap<usize, FeatInfo>,
) -> Result<Vec<(usize, String)>, String> {
    let text = std::fs::read_to_string(m8_path).map_err(|e| format!("read {m8_path}: {e}"))?;
    // query id -> (target -> hits), plus best (lowest-evalue) target per query.
    let mut per_query: HashMap<String, HashMap<String, Vec<Hit>>> = HashMap::new();
    let mut best_target: HashMap<String, (String, f64, i64)> = HashMap::new();

    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 12 {
            continue;
        }
        let q = p[0].to_string();
        let t = p[1].to_string();
        let hit = Hit {
            qstart: p[2].parse().unwrap_or(0),
            qend: p[3].parse().unwrap_or(0),
            tstart: p[4].parse().unwrap_or(0),
            tend: p[5].parse().unwrap_or(0),
            tlen: p[7].parse().unwrap_or(0),
            evalue: p[10].parse().unwrap_or(f64::INFINITY),
        };
        if hit.evalue > EVALUE_MAX {
            continue;
        }
        match best_target.get(&q) {
            Some((_, ev, _)) if *ev <= hit.evalue => {}
            _ => {
                best_target.insert(q.clone(), (t.clone(), hit.evalue, hit.tlen));
            }
        }
        per_query
            .entry(q)
            .or_default()
            .entry(t)
            .or_default()
            .push(hit);
    }

    let mut flags: Vec<(usize, String)> = Vec::new();
    for (q, targets) in &per_query {
        // decode id: idx|start|end|strand|edge
        let mut it = q.split('|');
        let idx: usize = match it.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let _start = it.next();
        let _end = it.next();
        let _strand = it.next();
        let edge = it.next().unwrap_or("0");
        if edge == "1" {
            continue;
        }
        let (prim, _ev, tlen) = match best_target.get(q) {
            Some(v) => v,
            None => continue,
        };
        if *tlen < MIN_TLEN {
            continue;
        }
        let hits = match targets.get(prim) {
            Some(h) => h,
            None => continue,
        };
        if let Some(cls) = classify_split(hits, *tlen) {
            flags.push((
                idx,
                format!("aligned {cls} vs consensus family {prim} (mmseqs translated split)"),
            ));
        } else if let Some(note) =
            classify_truncation(hits, *tlen, info.get(&idx).unwrap_or(&FeatInfo::default()))
        {
            flags.push((idx, format!("{note} vs consensus family {prim}")));
        }
    }
    Ok(flags)
}

/// Full mmseqs-translated disruption detector. Writes candidate regions, builds
/// mmseqs DBs, runs the translated search + convertalis, parses the result. All
/// work happens under a per-process temp dir, cleaned up on return.
///
/// `consensus_fasta` must already exist (see [`build_consensus_fasta`]).
///
/// `prebuilt_target`: when `Some(db)`, search against the prebuilt mmseqs target DB
/// `db` (e.g. the `--psc` UniRef90 `psc_db`) instead of building one from
/// `consensus_fasta` — this gives BROAD coverage so the split / frameshift /
/// internal-stop signals reach non-IS disrupted genes the narrow ncbifams consensus
/// misses. `consensus_fasta` is then unused for the target. When `None`, the
/// ncbifams consensus is used (the default, DB-light path).
pub fn detect_mmseqs(
    contigs: &[Contig],
    features: &[Feature],
    consensus_fasta: &str,
    prebuilt_target: Option<&str>,
    threads: usize,
) -> Vec<(usize, String)> {
    if !crate::mmseqs::available() {
        eprintln!("pseudogene-align: mmseqs not found; alignment signal skipped");
        return Vec::new();
    }
    let bin = crate::mmseqs::bin();
    let work = std::env::temp_dir().join(format!("bactars_pseudoaln_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&work);
    if let Err(e) = std::fs::create_dir_all(&work) {
        eprintln!("pseudogene-align: cannot create temp dir {work:?}: {e}");
        return Vec::new();
    }
    let p = |n: &str| work.join(n).to_string_lossy().into_owned();
    let regions = p("regions.fna");
    let qdb = p("qdb");
    // Target DB: the prebuilt --psc UniRef90 DB (broad), or a consensus DB we build.
    let tdb = prebuilt_target
        .map(|s| s.to_string())
        .unwrap_or_else(|| p("tdb"));
    let res = p("res");
    let tmp = p("tmp");
    let m8 = p("hits.m8");
    let nthreads = if threads == 0 {
        num_cpus_str()
    } else {
        threads.to_string()
    };

    let n = match write_regions(contigs, features, &regions) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("pseudogene-align: region extraction failed: {e}");
            let _ = std::fs::remove_dir_all(&work);
            return Vec::new();
        }
    };
    if n == 0 {
        let _ = std::fs::remove_dir_all(&work);
        return Vec::new();
    }

    let steps: Result<(), String> = (|| {
        if prebuilt_target.is_none() {
            crate::mmseqs::run_with(&bin, &["createdb", consensus_fasta, &tdb, "-v", "0"])?;
        }
        crate::mmseqs::run_with(&bin, &["createdb", &regions, &qdb, "--dbtype", "2", "-v", "0"])?;
        crate::mmseqs::run_with(
            &bin,
            &[
                "search", &qdb, &tdb, &res, &tmp,
                "--search-type", "2",
                "-e", "1e-5",
                "--threads", &nthreads,
                "-v", "0",
            ],
        )?;
        crate::mmseqs::run_with(
            &bin,
            &[
                "convertalis", &qdb, &tdb, &res, &m8,
                "--format-output",
                "query,target,qstart,qend,tstart,tend,qlen,tlen,alnlen,pident,evalue,bits",
                "--threads", &nthreads,
                "-v", "0",
            ],
        )?;
        Ok(())
    })();

    // Per-CDS context for the truncation signals (protein length, transposase flag).
    let info: HashMap<usize, FeatInfo> = features
        .iter()
        .enumerate()
        .filter(|(_, f)| f.kind == FeatureKind::Cds)
        .map(|(i, f)| {
            (
                i,
                FeatInfo {
                    aa_len: f.aa.as_ref().map(|a| a.len() as i64).unwrap_or(0),
                    is_transpos: feat_is_transpos(f),
                },
            )
        })
        .collect();

    let flags = match steps {
        Ok(()) => parse_and_detect(&m8, &info).unwrap_or_else(|e| {
            eprintln!("pseudogene-align: parse failed: {e}");
            Vec::new()
        }),
        Err(e) => {
            eprintln!("pseudogene-align: mmseqs pipeline failed: {e}");
            Vec::new()
        }
    };
    let _ = std::fs::remove_dir_all(&work);
    flags
}

/// Minimum intergenic gap (nt) worth scanning for a missed/degraded gene. Below
/// this a real coding remnant is unlikely and the search is just noise.
pub const MIN_IG_LEN: i64 = 150;
/// An intergenic hit must cover at least this fraction of the reference family
/// protein to count as a degraded gene (vs a spurious short/domain match).
pub const IG_MIN_TCOV: f64 = 0.40;

/// Intergenic regions (contig_index, start1, end1; 1-based inclusive) — the gaps
/// between all called features on each contig, at least [`MIN_IG_LEN`] nt. A
/// pseudogene the gene finder missed ENTIRELY (e.g. degraded past the caller's
/// threshold) lives here, so these are the query set for the intergenic scan.
fn intergenic_regions(contigs: &[Contig], features: &[Feature]) -> Vec<(usize, i64, i64)> {
    let mut by_contig: HashMap<&str, Vec<(i64, i64)>> = HashMap::new();
    for f in features {
        by_contig
            .entry(f.contig.as_str())
            .or_default()
            .push((f.start, f.end));
    }
    let mut out = Vec::new();
    for (ci, c) in contigs.iter().enumerate() {
        let clen = c.seq.len() as i64;
        let mut iv = by_contig.get(c.name.as_str()).cloned().unwrap_or_default();
        iv.sort_unstable();
        // Walk the merged occupied intervals, emitting the gaps between them.
        let mut cursor = 1i64; // next free base (1-based)
        for (s, e) in iv {
            if s - cursor >= MIN_IG_LEN {
                out.push((ci, cursor, s - 1));
            }
            cursor = cursor.max(e + 1);
        }
        if clen - cursor + 1 >= MIN_IG_LEN {
            out.push((ci, cursor, clen));
        }
    }
    out
}

/// Write intergenic nucleotide regions as FASTA queries. Id = `IG|ci|start|end`
/// (contig index + 1-based inclusive contig coordinates). Forward strand only —
/// the translated mmseqs search (`--search-type 2`) checks all six frames.
fn write_intergenic(
    contigs: &[Contig],
    regions: &[(usize, i64, i64)],
    out_path: &str,
) -> Result<usize, String> {
    let mut out =
        BufWriter::new(File::create(out_path).map_err(|e| format!("create {out_path}: {e}"))?);
    let mut n = 0usize;
    for &(ci, s, e) in regions {
        let seq = &contigs[ci].seq;
        let sub = &seq[(s - 1) as usize..=(e - 1) as usize];
        writeln!(out, ">IG|{ci}|{s}|{e}").map_err(|err| err.to_string())?;
        for chunk in sub.chunks(70) {
            out.write_all(chunk).map_err(|err| err.to_string())?;
            out.write_all(b"\n").map_err(|err| err.to_string())?;
        }
        n += 1;
    }
    out.flush().map_err(|e| e.to_string())?;
    Ok(n)
}

/// Parse the intergenic m8 into NEW pseudogene features (one per intergenic region
/// with a confident reference hit). Keeps the single best (lowest-evalue) hit per
/// region; coordinates are mapped from the region-local nucleotide hit back onto
/// the contig, strand from the translated frame.
fn parse_intergenic(m8_path: &str, contigs: &[Contig]) -> Result<Vec<Feature>, String> {
    let text = std::fs::read_to_string(m8_path).map_err(|e| format!("read {m8_path}: {e}"))?;
    // region id -> best (evalue, feature)
    let mut best: HashMap<String, (f64, Feature)> = HashMap::new();
    for line in text.lines() {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 12 {
            continue;
        }
        let q = p[0];
        let target = p[1];
        let qstart: i64 = p[2].parse().unwrap_or(0);
        let qend: i64 = p[3].parse().unwrap_or(0);
        let tstart: i64 = p[4].parse().unwrap_or(0);
        let tend: i64 = p[5].parse().unwrap_or(0);
        let tlen: i64 = p[7].parse().unwrap_or(0);
        let evalue: f64 = p[10].parse().unwrap_or(f64::INFINITY);
        if evalue > EVALUE_MAX || tlen < MIN_TLEN {
            continue;
        }
        let tcov = ((tend - tstart).abs() + 1) as f64 / tlen as f64;
        if tcov < IG_MIN_TCOV {
            continue;
        }
        // decode id: IG|ci|rs|re
        let mut it = q.split('|');
        if it.next() != Some("IG") {
            continue;
        }
        let ci: usize = match it.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        let rs: i64 = match it.next().and_then(|s| s.parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if ci >= contigs.len() {
            continue;
        }
        let (strand, _) = frame_of(qstart, qend);
        let lo = rs + qstart.min(qend) - 1;
        let hi = rs + qstart.max(qend) - 1;
        let feat = Feature {
            kind: FeatureKind::Cds,
            contig: contigs[ci].name.clone(),
            id: format!("igpseudo_{}_{}", contigs[ci].name, lo),
            start: lo,
            end: hi,
            strand,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: crate::feature::Functional {
                pseudogene: true,
                note: vec![format!(
                    "intergenic homology to family {target} ({:.0}% of reference, translated) — degraded/uncalled gene",
                    tcov * 100.0
                )],
                ..Default::default()
            },
        };
        match best.get(q) {
            Some((ev, _)) if *ev <= evalue => {}
            _ => {
                best.insert(q.to_string(), (evalue, feat));
            }
        }
    }
    Ok(best.into_values().map(|(_, f)| f).collect())
}

/// Intergenic 6-frame disruption scan (Pseudofinder-style): translate every
/// inter-feature gap and search it against the reference family consensus, so a
/// gene the caller MISSED entirely (degraded past its threshold) is recovered as a
/// pseudogene feature. Reuses the translated-mmseqs pipeline of [`detect_mmseqs`];
/// returns NEW features (never mutates existing ones). `consensus_fasta` must exist.
///
/// `prebuilt_target`: as in [`detect_mmseqs`] — `Some(db)` searches the broad `--psc`
/// UniRef90 DB (recovers non-IS uncalled genes), `None` uses the ncbifams consensus.
pub fn detect_intergenic_mmseqs(
    contigs: &[Contig],
    features: &[Feature],
    consensus_fasta: &str,
    prebuilt_target: Option<&str>,
    threads: usize,
) -> Vec<Feature> {
    if !crate::mmseqs::available() {
        return Vec::new();
    }
    let regions = intergenic_regions(contigs, features);
    if regions.is_empty() {
        return Vec::new();
    }
    let bin = crate::mmseqs::bin();
    let work = std::env::temp_dir().join(format!("bactars_pseudoig_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&work);
    if std::fs::create_dir_all(&work).is_err() {
        return Vec::new();
    }
    let p = |n: &str| work.join(n).to_string_lossy().into_owned();
    let (regfa, qdb, res, tmp, m8) =
        (p("ig.fna"), p("qdb"), p("res"), p("tmp"), p("hits.m8"));
    let tdb = prebuilt_target
        .map(|s| s.to_string())
        .unwrap_or_else(|| p("tdb"));
    let nthreads = if threads == 0 {
        num_cpus_str()
    } else {
        threads.to_string()
    };
    let n = write_intergenic(contigs, &regions, &regfa).unwrap_or(0);
    if n == 0 {
        let _ = std::fs::remove_dir_all(&work);
        return Vec::new();
    }
    let steps: Result<(), String> = (|| {
        if prebuilt_target.is_none() {
            crate::mmseqs::run_with(&bin, &["createdb", consensus_fasta, &tdb, "-v", "0"])?;
        }
        crate::mmseqs::run_with(&bin, &["createdb", &regfa, &qdb, "--dbtype", "2", "-v", "0"])?;
        crate::mmseqs::run_with(
            &bin,
            &[
                "search", &qdb, &tdb, &res, &tmp, "--search-type", "2", "-e", "1e-5",
                "--threads", &nthreads, "-v", "0",
            ],
        )?;
        crate::mmseqs::run_with(
            &bin,
            &[
                "convertalis", &qdb, &tdb, &res, &m8, "--format-output",
                "query,target,qstart,qend,tstart,tend,qlen,tlen,alnlen,pident,evalue,bits",
                "--threads", &nthreads, "-v", "0",
            ],
        )?;
        Ok(())
    })();
    let feats = match steps {
        Ok(()) => parse_intergenic(&m8, contigs).unwrap_or_else(|e| {
            eprintln!("pseudogene-align: intergenic parse failed: {e}");
            Vec::new()
        }),
        Err(e) => {
            eprintln!("pseudogene-align: intergenic mmseqs pipeline failed: {e}");
            Vec::new()
        }
    };
    let _ = std::fs::remove_dir_all(&work);
    feats
}

/// Whether a CDS is a transposase / IS element, by product or HMM-hit name — drives
/// the B0 relaxed truncation ceiling (RefSeq pseudogenizes degraded IS copies freely).
pub fn feat_is_transpos(f: &Feature) -> bool {
    let hit = |s: &str| {
        let s = s.to_ascii_lowercase();
        s.contains("transpos") || s.contains("insertion sequence") || s.contains("is element")
    };
    if let Some(p) = &f.func.product {
        if hit(p) {
            return true;
        }
    }
    f.annotations.iter().any(|a| hit(&a.name))
}

fn num_cpus_str() -> String {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .to_string()
}

/// Default cache path for the consensus FASTA, alongside the ncbifams `.hmm`.
/// Falls back to the temp dir if the ncbifams directory is not writable.
pub fn consensus_cache_path(hmm_path: &str) -> String {
    let p = Path::new(hmm_path);
    if let Some(dir) = p.parent() {
        let cand = dir.join("ncbifams_consensus.faa");
        // Prefer alongside the HMM if the directory is writable.
        if dir.metadata().map(|m| !m.permissions().readonly()).unwrap_or(false) {
            return cand.to_string_lossy().into_owned();
        }
    }
    std::env::temp_dir()
        .join("bactars_ncbifams_consensus.faa")
        .to_string_lossy()
        .into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(qstart: i64, qend: i64, tstart: i64, tend: i64, tlen: i64) -> Hit {
        Hit { qstart, qend, tstart, tend, tlen, evalue: 1e-40 }
    }

    #[test]
    fn frame_of_distinguishes_phase() {
        assert_eq!(frame_of(1, 100).1, 0);
        assert_eq!(frame_of(2, 101).1, 1);
        assert_eq!(frame_of(3, 102).1, 2);
    }

    /// Two in-frame pieces covering consecutive target ranges = internal stop.
    #[test]
    fn classify_internal_stop() {
        let tlen = 200;
        // piece A: target 1..100, query nt 1..300 (phase 0)
        // piece B: target 101..200, query nt 310..610 (phase 0, 309%3==0)
        let hits = vec![h(1, 300, 1, 100, tlen), h(310, 610, 101, 200, tlen)];
        assert_eq!(classify_split(&hits, tlen), Some("internal_stop"));
    }

    /// Two pieces in DIFFERENT frames covering consecutive target ranges = frameshift.
    #[test]
    fn classify_frameshift() {
        let tlen = 200;
        // A: phase 0 (query 1..300); B: phase 1 (query 311..610 -> (311-1)%3==1)
        let hits = vec![h(1, 300, 1, 100, tlen), h(311, 610, 101, 200, tlen)];
        assert_eq!(classify_split(&hits, tlen), Some("frameshift"));
    }

    /// A single full-length hit (intact gene) is NOT a split.
    #[test]
    fn classify_single_hit_not_flagged() {
        let tlen = 200;
        let hits = vec![h(1, 600, 1, 200, tlen)];
        assert_eq!(classify_split(&hits, tlen), None);
    }

    /// Two pieces that overlap heavily in the target (repeat/paralog artefact),
    /// not consecutive, must NOT be called a split.
    #[test]
    fn classify_overlapping_pieces_not_flagged() {
        let tlen = 200;
        let hits = vec![h(1, 300, 1, 100, tlen), h(310, 610, 5, 105, tlen)];
        // b.tstart(5) < a.tend(100) - T_OVERLAP_TOL(20) => rejected
        assert_eq!(classify_split(&hits, tlen), None);
    }

    /// RULE 1: two consecutive-in-model pieces that are genomically FAR APART
    /// (a whole gene between them) are two separate genes (tandem paralogs), NOT
    /// one broken gene — rejected by the genomic gap / extent constraints.
    #[test]
    fn classify_distant_pieces_not_flagged() {
        let tlen = 200;
        // A: model 1..100, query nt 1..300; B: model 101..200 but query nt
        // 1601..1900 — ~1300 nt (a full gene) downstream of A's end.
        let hits = vec![h(1, 300, 1, 100, tlen), h(1601, 1900, 101, 200, tlen)];
        // q_gap = 1601-300-1 = 1300 > MAX_Q_GAP_NT(300); q_span also > 1.5*3*tlen.
        assert_eq!(classify_split(&hits, tlen), None);
    }

    /// RULE 2 (model side): two pieces that TOGETHER cover ~2x the model (two full
    /// copies of a short family = a tandem duplication) are rejected by the
    /// combined-coverage upper bound.
    #[test]
    fn classify_tandem_full_copies_not_flagged() {
        let tlen = 30;
        // A: model 1..25 (25 aa), B: model 6..30 (25 aa) — overlap 20 aa (== tol),
        // sum = 50 aa = 1.67x model > MAX_COMBINED_TCOV(1.25). Two ~full copies.
        let hits = vec![h(1, 75, 1, 25, tlen), h(76, 150, 6, 30, tlen)];
        assert_eq!(classify_split(&hits, tlen), None);
    }

    /// RULE 4: two pieces on OPPOSITE strands are not a co-linear disruption.
    #[test]
    fn classify_opposite_strand_not_flagged() {
        let tlen = 200;
        // B reads on the minus strand (qstart > qend) => different strand from A.
        let hits = vec![h(1, 300, 1, 100, tlen), h(610, 310, 101, 200, tlen)];
        assert_eq!(classify_split(&hits, tlen), None);
    }

    /// A genomically-close, single-model-extent, consecutive split is still called
    /// (regression guard that the tighter geometry did not kill the true signal).
    #[test]
    fn classify_close_split_still_flagged() {
        let tlen = 200;
        let hits = vec![h(1, 300, 1, 100, tlen), h(304, 610, 101, 200, tlen)];
        assert_eq!(classify_split(&hits, tlen), Some("internal_stop"));
    }

    #[test]
    fn consensus_fasta_from_synthetic_hmm() {
        // A tiny 3-node HMM: node1 min at col 0 (A), node2 min at col 2 (D),
        // node3 min at col 19 (Y). score_to_prob monotonic => argmin == our pick.
        let hmm = "\
HMMER3/f [test]
NAME  toy
ACC   NFZZZ.1
LENG  3
ALPH  amino
HMM      A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y
         m->m
  COMPO  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
         2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
         0 0 0 0 0 0 0
      1  0.1 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
         2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
         0 0 0 0 0 0 0
      2  9 9 0.1 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
         2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
         0 0 0 0 0 0 0
      3  9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 0.1
         2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
         0 0 0 0 0 0 0
//
";
        let dir = std::env::temp_dir().join(format!("bactars_cons_test_{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let hp = dir.join("toy.hmm");
        let op = dir.join("toy.faa");
        std::fs::write(&hp, hmm).unwrap();
        let n = build_consensus_fasta(hp.to_str().unwrap(), op.to_str().unwrap()).unwrap();
        assert_eq!(n, 1);
        let out = std::fs::read_to_string(&op).unwrap();
        assert!(out.contains(">NFZZZ.1"), "fasta: {out}");
        assert!(out.contains("ADY"), "consensus should be ADY; got: {out}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
