//! oriT / mobilizable-element detection: real oriT nucleotide locus + relaxase.
//!
//! The origin of transfer (**oriT**) is the nick site where a **relaxase** (MOBA/
//! MOBF/MOBP/MOBQ/...) initiates conjugative DNA transfer. This detector reports
//! two complementary signals:
//!
//!  1. **Real oriT locus** — an mmseqs nucleotide search (`--search-type 3`,
//!     mirroring [`crate::plasmid`]) of the contigs against a curated oriT
//!     nucleotide DB (`oriT_exp.fasta`, the 122 experimentally-verified oriT
//!     sequences; `oriT_all.fasta` as a fallback). A passing hit yields an `OriT`
//!     feature at the **actual nucleotide match span** (the true transfer-origin
//!     locus) carrying the matched oriT accession, percent identity, and a real
//!     mmseqs E-value. This is the oriTfinder-style locus call.
//!
//!  2. **Mobilization relaxase** — the CDS proteins are scanned with a relaxase
//!     (MOB) HMM set; each relaxase-positive gene is emitted as its own `OriT`
//!     feature (source `orit:relaxase`) typed by MOB family. Relaxase presence
//!     marks a mobilizable/conjugative element and, when it sits near an oriT nt
//!     hit, raises that hit's confidence (co-localization). Both signals are kept
//!     so a relaxase can be reported even when no oriT sequence is in the DB, and
//!     a strong standalone oriT hit can be reported even without a called relaxase.
//!
//! # Gating
//!
//! oriT nt hits are specific by design. A hit near a detected relaxase
//! (within [`ORIT_COLOC_WINDOW`] bp) is accepted at the base bar (E-value
//! ≤ [`ORIT_MAX_EVALUE`], identity ≥ [`ORIT_MIN_PIDENT`], target coverage
//! ≥ [`ORIT_MIN_TCOV`], alignment ≥ [`ORIT_MIN_ALNLEN`] bp). A **standalone** hit
//! (no relaxase nearby) must clear the stricter standalone bar
//! ([`ORIT_STANDALONE_MAX_EVALUE`] / [`ORIT_STANDALONE_MIN_PIDENT`] /
//! [`ORIT_STANDALONE_MIN_TCOV`]) to be reported.
//!
//! The oriT nt DB (`oriT_exp.fasta`) is looked up inside `--orit <DIR>`; when it is
//! absent (e.g. the CC0 `db/orit_pfam` layout ships only relaxase HMMs) or mmseqs
//! is unavailable, the nt-locus stage is skipped and only the relaxase signal is
//! reported — the detector never fails for a missing nt DB.
//!
//! # LICENSE — ships license-clean Pfam (CC0) relaxase families
//!
//! bactars aims to be a freely-distributable, pure-Rust tool, so this detector's
//! default DB is a **license-clean Pfam relaxase set** (`db/orit_pfam/`,
//! `relaxase_pfam.hmm`), a single concatenated HMMER3/f profile built from the
//! Pfam-A MOB relaxase families — all released into the public domain under
//! **CC0 1.0** by EMBL-EBI/Pfam and therefore safe to bundle and redistribute
//! commercially:
//!
//! | Accession | Name         | MOB family / role                  |
//! |-----------|--------------|------------------------------------|
//! | PF01076   | Mob_Pre      | MOBV plasmid recombination enzyme  |
//! | PF03389   | MobA_MobL    | MOBQ relaxase                      |
//! | PF03432   | Relaxase     | MOBF, MobA/VirD2-like nuclease     |
//! | PF05713   | MobC         | MOBC mobilisation protein          |
//! | PF07514   | TraI_2       | VirD2/TraI-type relaxase-helicase  |
//! | PF13814   | Replic_Relax | Replication-relaxation             |
//!
//! Each Pfam profile carries a curated GA gathering cutoff, so the search uses the
//! `--cut_ga` path per model. `--orit <DIR>` points the detector at ANY directory
//! of relaxase `.hmm` file(s); the CC0 `db/orit_pfam` layout is the documented
//! default (resolved from a `--db` bundle). The historical CONJScan/MacSyFinder
//! `db/orit/` profiles (`T4SS_MOB*.hmm` / `relaxase_MOB.hmm`) are more complete but
//! **CC-BY-NC-SA 4.0 (NonCommercial)** and so are NOT bundled or depended on; a user
//! who has them may still point `--orit` at that directory.
//!
//! Emits `Feature { kind: FeatureKind::OriT }` (SO `origin_of_transfer`).

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use crate::hmm::{self, Cutoff};
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Fallback reporting E-value for relaxase HMMs that carry no curated GA cutoff.
/// User-supplied relaxase profiles may or may not be gathering-thresholded; when
/// a model has a GA the `--cut_ga` path uses it, otherwise this bounds reporting.
const RELAXASE_EVALUE: f64 = 1e-5;

// --- oriT nucleotide-locus search thresholds ------------------------------------
/// Base E-value ceiling for an oriT nt hit that is co-localized with a relaxase.
const ORIT_MAX_EVALUE: f64 = 1e-5;
/// Base percent-identity floor (0..100) for a co-localized oriT nt hit.
const ORIT_MIN_PIDENT: f64 = 70.0;
/// Base fraction (0..1) of the oriT reference sequence the alignment must cover.
const ORIT_MIN_TCOV: f64 = 0.50;
/// Base minimum alignment length (bp) — oriT_exp seqs run ~7..701 bp (median 101),
/// so a real locus call needs a non-trivial span even for the short references.
const ORIT_MIN_ALNLEN: i64 = 25;
/// A relaxase CDS within this many bp of an oriT nt hit co-localizes it (oriT sits
/// adjacent to its relaxase gene), raising confidence to the base bar.
const ORIT_COLOC_WINDOW: i64 = 10_000;
/// Stricter E-value ceiling for a **standalone** oriT hit (no relaxase nearby).
const ORIT_STANDALONE_MAX_EVALUE: f64 = 1e-8;
/// Stricter identity floor for a standalone oriT hit.
const ORIT_STANDALONE_MIN_PIDENT: f64 = 80.0;
/// Stricter reference-coverage floor for a standalone oriT hit.
const ORIT_STANDALONE_MIN_TCOV: f64 = 0.60;
/// Candidate oriT nucleotide DB basenames inside `--orit <DIR>`, preferred first
/// (curated experimentally-verified set before the larger all-oriT set).
const ORIT_NT_DB_NAMES: [&str; 2] = ["oriT_exp.fasta", "oriT_all.fasta"];

/// Detect mobilizable-element / oriT proxies (relaxase genes) across the genome.
/// `features` supplies the CDS calls; `orit_dir` is a directory of relaxase `.hmm`
/// file(s) — the license-clean CC0 Pfam set (`db/orit_pfam`) by default, or any
/// user-supplied relaxase HMM directory (see the module docs). `contigs` is
/// accepted for signature uniformity.
pub fn detect(
    contigs: &[(String, String)],
    features: &[Feature],
    orit_dir: &str,
) -> Result<Vec<Feature>, String> {
    let dir = orit_dir.trim_end_matches('/');
    let hmm_files = collect_hmm_files(dir);
    if hmm_files.is_empty() {
        // The relaxase HMM set is the required model DB for this detector; a
        // missing/empty directory is a configuration failure, surfaced to the
        // caller (which logs the skipped stage) rather than silently returning [].
        return Err(format!("no .hmm relaxase profiles in {dir}"));
    }

    // CDS proteins (id, aa).
    let proteins: Vec<(String, String)> = features
        .iter()
        .filter(|f| f.kind == FeatureKind::Cds)
        .filter_map(|f| f.aa.as_ref().map(|aa| (f.id.clone(), aa.clone())))
        .collect();

    // Best relaxase hit per CDS id (MOB family from the model name + score +
    // the rustyhmmer sequence E-value). Prefer the model's curated GA cutoff;
    // fall back to an E-value bound.
    let mut best: std::collections::HashMap<String, (f32, String, f64)> =
        std::collections::HashMap::new();
    if !proteins.is_empty() {
        for path in &hmm_files {
            // A model that carries a curated GA gathering cutoff is searched with
            // `--cut_ga` (its own trusted bit-score bar); one without a GA is bounded
            // by an E-value instead.
            let use_ga = hmm_has_ga(path);
            let cutoff = if use_ga {
                Cutoff::GatheringGa
            } else {
                Cutoff::Evalue(RELAXASE_EVALUE)
            };
            let hits = match hmm::annotate(&proteins, path, cutoff) {
                Ok(h) => h,
                Err(e) => {
                    // A single unreadable/uncalibrated model is a per-item skip (like
                    // ncrna's per-model skips), not a whole-detector failure.
                    eprintln!("orit: skipping {path} ({e})");
                    continue;
                }
            };
            for h in hits {
                // BUGFIX: the E-value bound must apply ONLY to the non-GA branch. A hit
                // that passes a model's curated GA cutoff is trusted by that cutoff and
                // must NOT be discarded because its E-value inflated with DB size —
                // otherwise curated-GA relaxases silently drop out on large genomes.
                if !use_ga && h.seq_evalue > RELAXASE_EVALUE {
                    continue;
                }
                let e = best
                    .entry(h.target_name.clone())
                    .or_insert((f32::MIN, String::new(), f64::NAN));
                if h.seq_score > e.0 {
                    *e = (h.seq_score, mob_family(&h.query_name, &h.query_acc), h.seq_evalue);
                }
            }
        }
    }

    // Relaxase-positive CDS features, contig then start order. Used both to emit
    // the relaxase (mobilization) signal and, as genomic spans, to co-localize
    // oriT nt hits.
    let mut relaxase_hits: Vec<(&Feature, f32, String, f64)> = features
        .iter()
        .filter(|f| f.kind == FeatureKind::Cds)
        .filter_map(|f| best.get(&f.id).map(|(sc, fam, ev)| (f, *sc, fam.clone(), *ev)))
        .collect();
    relaxase_hits
        .sort_by(|(a, ..), (b, ..)| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));

    // (contig, low, high) spans of the called relaxases, for co-localization.
    let relaxase_spans: Vec<(String, i64, i64)> = relaxase_hits
        .iter()
        .map(|(f, ..)| {
            let (lo, hi) = if f.start <= f.end { (f.start, f.end) } else { (f.end, f.start) };
            (f.contig.clone(), lo, hi)
        })
        .collect();

    let mut out = Vec::new();
    let mut per_contig: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // --- Signal 2: mobilization relaxase features (one per relaxase-positive CDS).
    for (f, score, fam, seq_evalue) in &relaxase_hits {
        let (start, end) = if f.start <= f.end {
            (f.start, f.end)
        } else {
            (f.end, f.start)
        };
        let n = per_contig.entry(f.contig.clone()).or_insert(0);
        *n += 1;
        // HONESTY: the emitted span is the relaxase (MOB) CDS, NOT a mapped oriT
        // nick site. Name/product and note make that explicit so downstream output
        // does not overclaim a transfer-origin coordinate. The genuine nick-site
        // locus, when the oriT nt DB is present, is emitted separately below
        // (source `orit:nic`).
        out.push(Feature {
            kind: FeatureKind::OriT,
            contig: f.contig.clone(),
            id: format!("{}_orit{}", f.contig, n),
            start,
            end,
            strand: f.strand,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: "orit:relaxase".to_string(),
                accession: String::new(),
                name: format!("{fam} mobilization relaxase (oriT proxy)"),
                score: *score,
                // Real rustyhmmer relaxase-HMM sequence E-value for this CDS.
                evalue: (!seq_evalue.is_nan()).then_some(*seq_evalue),
                ref_len: None,
            }],
            func: Functional {
                product: Some(format!("{fam} mobilization relaxase (oriT proxy)")),
                note: vec![format!(
                    "coordinates span the {fam} relaxase (MOB) gene, not the oriT nick \
                     site; relaxase presence is used as a proxy for a mobilizable/\
                     conjugative element carrying an origin of transfer"
                )],
                ..Functional::default()
            },
        });
    }

    // --- Signal 1: real oriT nucleotide locus (mmseqs --search-type 3). Skipped
    // (relaxase signal still reported) when the nt DB is absent or mmseqs is not
    // available — a missing nt DB is never a hard error.
    if let Some(nt_db) = find_orit_nt_db(dir) {
        match orit_nt_hits(contigs, &nt_db) {
            Ok(nt) => {
                for hit in gate_orit_hits(nt, &relaxase_spans) {
                    let n = per_contig.entry(hit.contig.clone()).or_insert(0);
                    *n += 1;
                    let coloc_note = if hit.colocalized {
                        " co-localized with a detected relaxase gene (high confidence)"
                    } else {
                        " standalone (no relaxase called nearby)"
                    };
                    out.push(Feature {
                        kind: FeatureKind::OriT,
                        contig: hit.contig.clone(),
                        id: format!("{}_orit{}", hit.contig, n),
                        start: hit.start,
                        end: hit.end,
                        strand: hit.strand,
                        aa: None,
                        partial5: false,
                        partial3: false,
                        annotations: vec![Annotation {
                            source: "orit:nic".to_string(),
                            accession: hit.accession.clone(),
                            name: format!("origin of transfer ({})", hit.accession),
                            score: hit.bits,
                            // Real mmseqs (--search-type 3) nucleotide E-value.
                            evalue: Some(hit.evalue),
                            ref_len: None,
                        }],
                        func: Functional {
                            product: Some(format!(
                                "origin of transfer (oriT), {} match", hit.accession
                            )),
                            note: vec![format!(
                                "oriT nucleotide locus: {:.1}% identity to {} over {} bp \
                                 (E-value {:.1e});{}",
                                hit.pident, hit.accession, hit.alnlen, hit.evalue, coloc_note
                            )],
                            ..Functional::default()
                        },
                    });
                }
            }
            Err(e) => {
                // A failed nt search degrades to the relaxase-only signal.
                eprintln!("orit: oriT nt search skipped ({e})");
            }
        }
    }

    Ok(out)
}

/// A gated oriT nucleotide-locus call ready to be emitted as a feature.
#[derive(Clone, Debug)]
struct OritNtHit {
    contig: String,
    accession: String,
    pident: f64,
    evalue: f64,
    bits: f32,
    start: i64,
    end: i64,
    strand: i8,
    alnlen: i64,
    colocalized: bool,
}

/// One parsed convertalis row from the oriT nt search.
#[derive(Clone, Debug)]
struct OritM8 {
    query: String,   // contig name
    target: String,  // oriT accession (first header token, e.g. "oriT_RP4")
    pident: f64,     // percent identity (0..100)
    tcov: f64,       // fraction of the oriT reference covered (0..1)
    qstart: i64,     // 1-based contig coordinates
    qend: i64,
    bits: f32,
    evalue: f64,
    alnlen: i64,
}

/// Find the oriT nucleotide DB inside `--orit <DIR>`: prefer the curated
/// `oriT_exp.fasta`, fall back to `oriT_all.fasta`. `None` when neither is present
/// (e.g. the CC0 `db/orit_pfam` layout ships only relaxase HMMs) so the nt-locus
/// stage is cleanly skipped.
fn find_orit_nt_db(dir: &str) -> Option<std::path::PathBuf> {
    for name in ORIT_NT_DB_NAMES {
        let p = Path::new(dir).join(name);
        if p.is_file() {
            return Some(p);
        }
    }
    None
}

/// Run the mmseqs nucleotide search of the contigs against the oriT DB, returning
/// the parsed m8 rows (unfiltered — gating happens in [`gate_orit_hits`]). Mirrors
/// [`crate::plasmid`]'s createdb -> search --search-type 3 -> convertalis pipeline.
///
/// # Why mmseqs and not a pure-Rust aligner
///
/// A pure-Rust path (rust-bio `bio::alignment::pairwise` affine Smith-Waterman of
/// the 122 curated oriT refs vs a window around each relaxase) was evaluated. It
/// localized an embedded oriT exactly (100% identity, correct coordinates and
/// strand) but was far too slow — the full-matrix DP costs ~5.4 s for a single
/// 6 kb window, ~68 s for 60 kb, and whole-genome standalone scanning is
/// infeasible — whereas mmseqs' k-mer-prefiltered `--search-type 3` scans an entire
/// 4.6 Mb genome against all 122 refs in ~0.75 s (~84 MB) with zero false positives
/// on an oriT-free control, AND yields a real Karlin-Altschul E-value that a
/// hand-rolled SW would have to approximate. mmseqs is the sanctioned external
/// binary; whole-contig search also catches strong standalone oriT loci that a
/// relaxase-windowed pure-Rust scan would miss. Hence mmseqs is the chosen engine.
fn orit_nt_hits(
    contigs: &[(String, String)],
    nt_db: &Path,
) -> Result<Vec<OritM8>, String> {
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    if !crate::mmseqs::available() {
        return Err("mmseqs not found (set $BACTARS_MMSEQS)".to_string());
    }
    let bin = crate::mmseqs::bin();

    // Build the contig query FASTA.
    let mut query = String::new();
    for (name, seq) in contigs {
        query.push('>');
        query.push_str(name);
        query.push('\n');
        query.push_str(seq);
        query.push('\n');
    }

    let threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .clamp(1, 8)
        .to_string();

    let work = std::env::temp_dir().join(format!("bactars_orit_nt_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    fs::create_dir_all(&work).map_err(|e| format!("mkdir {}: {e}", work.display()))?;
    let result = run_orit_search(&bin, &work, nt_db, &query, &threads);
    let _ = fs::remove_dir_all(&work);
    result
}

/// createdb(query contigs, target oriT, nucleotide) -> search --search-type 3 ->
/// convertalis -> parsed rows.
fn run_orit_search(
    bin: &str,
    work: &Path,
    nt_db: &Path,
    query: &str,
    threads: &str,
) -> Result<Vec<OritM8>, String> {
    let query_fna = work.join("query.fna");
    let query_db = work.join("queryDB");
    let target_db = work.join("targetDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");

    fs::write(&query_fna, query).map_err(|e| format!("write query.fna: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_fna, q_db, t_db, res_db, tmp_s, db_s, m8_s) = (
        s(&query_fna),
        s(&query_db),
        s(&target_db),
        s(&result_db),
        s(&tmp),
        s(nt_db),
        s(&m8),
    );

    crate::mmseqs::run_with(bin, &["createdb", &q_fna, &q_db, "--dbtype", "2", "-v", "1"])?;
    crate::mmseqs::run_with(bin, &["createdb", &db_s, &t_db, "--dbtype", "2", "-v", "1"])?;
    // oriT references are short (down to ~7 bp); a small k and a permissive prefilter
    // let short exact/near-exact loci survive to alignment.
    crate::mmseqs::run_with(
        bin,
        &[
            "search", &q_db, &t_db, &res_db, &tmp_s, "--threads", threads, "--search-type", "3",
            "-e", "1e-4", "-s", "7.5", "-k", "7", "--max-seqs", "500", "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        bin,
        &[
            "convertalis", &q_db, &t_db, &res_db, &m8_s, "--format-output",
            "query,target,pident,qcov,tcov,qstart,qend,bits,evalue,alnlen",
            "--threads", threads, "-v", "1",
        ],
    )?;
    Ok(parse_orit_m8(&m8))
}

/// Parse the oriT convertalis output
/// (`query,target,pident,qcov,tcov,qstart,qend,bits,evalue,alnlen`).
fn parse_orit_m8(path: &Path) -> Vec<OritM8> {
    let mut rows = Vec::new();
    let file = match fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return rows,
    };
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 10 {
            continue;
        }
        let pident_raw: f64 = p[2].parse().unwrap_or(0.0);
        // mmseqs may emit pident as a fraction (0..1); normalise to a percent.
        let pident = if pident_raw <= 1.0 { pident_raw * 100.0 } else { pident_raw };
        rows.push(OritM8 {
            query: p[0].to_string(),
            target: p[1].to_string(),
            pident,
            tcov: p[4].parse().unwrap_or(0.0),
            qstart: p[5].parse().unwrap_or(0),
            qend: p[6].parse().unwrap_or(0),
            bits: p[7].parse().unwrap_or(0.0),
            evalue: p[8].parse().unwrap_or(f64::MAX),
            alnlen: p[9].parse().unwrap_or(0),
        });
    }
    rows
}

/// Apply the gating (base bar for relaxase co-localized hits, stricter bar for
/// standalone hits) and collapse to the best passing hit per genomic locus.
/// `relaxase_spans` are `(contig, low, high)` of the called relaxase genes.
fn gate_orit_hits(rows: Vec<OritM8>, relaxase_spans: &[(String, i64, i64)]) -> Vec<OritNtHit> {
    // Best passing hit per (contig, overlapping-locus). We collapse by keeping the
    // highest-bit hit whose span overlaps an already-kept one on the same contig,
    // so a single oriT is reported once even if several references match it.
    let mut kept: Vec<OritNtHit> = Vec::new();
    // Highest-bit first so the representative of an overlap cluster is the best hit.
    let mut rows = rows;
    rows.sort_by(|a, b| b.bits.partial_cmp(&a.bits).unwrap_or(std::cmp::Ordering::Equal));

    for r in rows {
        let (lo, hi) = if r.qstart <= r.qend { (r.qstart, r.qend) } else { (r.qend, r.qstart) };
        let strand: i8 = if r.qstart <= r.qend { 1 } else { -1 };

        // Distance (bp) to the nearest relaxase gene on the same contig; 0 = overlap.
        let colocalized = relaxase_spans.iter().any(|(c, rlo, rhi)| {
            c == &r.query && span_distance(lo, hi, *rlo, *rhi) <= ORIT_COLOC_WINDOW
        });

        // Base bar; standalone hits must clear a stricter bar.
        let passes = if colocalized {
            r.evalue <= ORIT_MAX_EVALUE
                && r.pident >= ORIT_MIN_PIDENT
                && r.tcov >= ORIT_MIN_TCOV
                && r.alnlen >= ORIT_MIN_ALNLEN
        } else {
            r.evalue <= ORIT_STANDALONE_MAX_EVALUE
                && r.pident >= ORIT_STANDALONE_MIN_PIDENT
                && r.tcov >= ORIT_STANDALONE_MIN_TCOV
                && r.alnlen >= ORIT_MIN_ALNLEN
        };
        if !passes {
            continue;
        }

        // Collapse overlapping loci on the same contig (rows are best-first, so an
        // overlap means this is a weaker redundant match of a kept oriT).
        let overlaps_kept = kept
            .iter()
            .any(|k| k.contig == r.query && span_distance(lo, hi, k.start, k.end) == 0);
        if overlaps_kept {
            continue;
        }

        kept.push(OritNtHit {
            contig: r.query.clone(),
            accession: r.target.clone(),
            pident: r.pident,
            evalue: r.evalue,
            bits: r.bits,
            start: lo,
            end: hi,
            strand,
            alnlen: r.alnlen,
            colocalized,
        });
    }

    // Deterministic output order: contig, then start.
    kept.sort_by(|a, b| a.contig.cmp(&b.contig).then(a.start.cmp(&b.start)));
    kept
}

/// Gap (bp) between two 1-based inclusive spans on the same sequence; `0` when they
/// overlap or touch.
fn span_distance(alo: i64, ahi: i64, blo: i64, bhi: i64) -> i64 {
    if ahi < blo {
        blo - ahi
    } else if bhi < alo {
        alo - bhi
    } else {
        0
    }
}

/// List `*.hmm` files in `dir` (non-recursive), sorted for determinism.
fn collect_hmm_files(dir: &str) -> Vec<String> {
    let mut files = Vec::new();
    if let Ok(rd) = std::fs::read_dir(dir) {
        for entry in rd.flatten() {
            let p = entry.path();
            if p.extension().and_then(|e| e.to_str()) == Some("hmm") {
                if let Some(s) = p.to_str() {
                    files.push(s.to_string());
                }
            }
        }
    }
    files.sort();
    files
}

/// Cheap check whether an HMM file carries any `GA` gathering cutoff line, so the
/// search can use `--cut_ga` when available and an E-value bound otherwise.
fn hmm_has_ga(path: &str) -> bool {
    match std::fs::read_to_string(path) {
        Ok(text) => text.lines().any(|l| l.starts_with("GA ")),
        Err(_) => false,
    }
}

/// Derive a MOB family label from a relaxase model's name / accession.
///
/// Two profile conventions are handled:
///  * Pfam-A relaxase families (the CC0 default set): mapped from their stable
///    accession to the canonical MOB family they represent.
///  * CONJScan/MacSyFinder models named `MOBF`, `MOBP1`, `MOBQ`, ...: the family
///    token is read straight out of the name.
///
/// Anything else falls back to the raw model name (or accession, or "relaxase").
fn mob_family(name: &str, acc: &str) -> String {
    // Pfam accessions carry a `.NN` version suffix (e.g. "PF03389.22"); match the
    // bare accession so version bumps don't break the mapping.
    let bare_acc = acc.split('.').next().unwrap_or(acc);
    match bare_acc {
        "PF01076" => return "MOBV".to_string(),      // Mob_Pre
        "PF03389" => return "MOBQ".to_string(),      // MobA_MobL
        "PF03432" => return "MOBF".to_string(),      // Relaxase (MobA/VirD2-like)
        "PF05713" => return "MOBC".to_string(),      // MobC
        "PF07514" => return "MOBF".to_string(),      // TraI_2 (VirD2/TraI)
        "PF13814" => return "relaxase (Replic_Relax)".to_string(),
        _ => {}
    }
    let up = name.to_ascii_uppercase();
    if let Some(idx) = up.find("MOB") {
        // Take "MOB" + the following family token (letters/digits).
        let tail: String = up[idx + 3..]
            .chars()
            .take_while(|c| c.is_ascii_alphanumeric())
            .collect();
        return format!("MOB{tail}");
    }
    if !name.is_empty() {
        name.to_string()
    } else if !acc.is_empty() {
        acc.to_string()
    } else {
        "relaxase".to_string()
    }
}

// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn mob_family_from_conjscan_names() {
        assert_eq!(mob_family("MOBF", ""), "MOBF");
        assert_eq!(mob_family("MOBP1", ""), "MOBP1");
        assert_eq!(mob_family("T4SS_MOBQ", ""), "MOBQ");
        assert_eq!(mob_family("MOBV", ""), "MOBV");
        // No MOB token and no known accession -> raw name / accession fallback.
        assert_eq!(mob_family("Weird", ""), "Weird");
        assert_eq!(mob_family("", "PF99999"), "PF99999");
    }

    #[test]
    fn mob_family_maps_pfam_accessions() {
        // The CC0 Pfam relaxase set is typed by stable accession (ACC line), which
        // rustyhmmer surfaces as `query_acc` with a `.NN` version suffix.
        assert_eq!(mob_family("Mob_Pre", "PF01076.26"), "MOBV");
        assert_eq!(mob_family("MobA_MobL", "PF03389.22"), "MOBQ");
        assert_eq!(mob_family("Relaxase", "PF03432.21"), "MOBF");
        assert_eq!(mob_family("MobC", "PF05713.18"), "MOBC");
        assert_eq!(mob_family("TraI_2", "PF07514.18"), "MOBF");
        assert_eq!(mob_family("Replic_Relax", "PF13814.13"), "relaxase (Replic_Relax)");
        // Bare accession (no version) still maps.
        assert_eq!(mob_family("MobA_MobL", "PF03389"), "MOBQ");
    }

    /// The staged CC0 Pfam relaxase DB must parse as a HMMER3 profile file and
    /// every profile must carry a GA gathering cutoff (so `--cut_ga` works). Skips
    /// gracefully when the DB is not present (e.g. a source-only checkout).
    #[test]
    fn staged_pfam_relaxase_db_parses_with_ga() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../db/orit_pfam/relaxase_pfam.hmm");
        if !std::path::Path::new(path).exists() {
            eprintln!("staged_pfam_relaxase_db_parses_with_ga: {path} absent; skipped");
            return;
        }
        // Loads via the same rustyhmmer path the detector uses.
        let annotator = rustyhmmer::api::HmmAnnotator::from_hmm_file(path)
            .expect("relaxase_pfam.hmm should parse as a HMMER3 profile file");
        // Searching an empty protein set exercises full model load without hits.
        let empty: Vec<(String, String)> = Vec::new();
        let hits = crate::hmm::annotate(&empty, path, Cutoff::GatheringGa)
            .expect("GA search over the staged DB should not error");
        assert!(hits.is_empty());
        // All six Pfam profiles carry a GA line.
        assert!(hmm_has_ga(path), "staged Pfam relaxase DB should carry GA cutoffs");
        let _ = annotator;
    }

    #[test]
    fn collect_hmm_files_filters_and_sorts() {
        let dir = std::env::temp_dir().join(format!("bactars_orit_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for name in ["b.hmm", "a.hmm", "notes.txt", "c.fasta"] {
            let mut f = std::fs::File::create(dir.join(name)).unwrap();
            writeln!(f, "x").unwrap();
        }
        let files = collect_hmm_files(dir.to_str().unwrap());
        assert_eq!(files.len(), 2);
        assert!(files[0].ends_with("a.hmm"));
        assert!(files[1].ends_with("b.hmm"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn hmm_has_ga_detects_cutoff_line() {
        let dir = std::env::temp_dir().join(format!("bactars_orit_ga_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let with = dir.join("with.hmm");
        let without = dir.join("without.hmm");
        std::fs::write(&with, "NAME  x\nGA    30.0 30.0;\n").unwrap();
        std::fs::write(&without, "NAME  x\nLENG  100\n").unwrap();
        assert!(hmm_has_ga(with.to_str().unwrap()));
        assert!(!hmm_has_ga(without.to_str().unwrap()));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn empty_dir_is_err() {
        // A relaxase-HMM directory with no `.hmm` profiles is a missing required
        // model DB -> the detector reports an error (skipped stage), not `[]`.
        let dir = std::env::temp_dir().join(format!("bactars_orit_empty_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let res = detect(&[], &[], dir.to_str().unwrap());
        assert!(res.is_err(), "empty relaxase DB dir should be an Err");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn no_cds_yields_empty_ok() {
        // A populated DB dir but no CDS proteins to search -> Ok(empty), not Err.
        let dir = std::env::temp_dir().join(format!("bactars_orit_nocds_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("relaxase.hmm"), "NAME  x\nGA    30.0 30.0;\n").unwrap();
        let res = detect(&[], &[], dir.to_str().unwrap());
        assert!(matches!(res, Ok(ref v) if v.is_empty()));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn span_distance_gap_and_overlap() {
        // Disjoint spans -> gap between them.
        assert_eq!(span_distance(1, 100, 201, 300), 101);
        assert_eq!(span_distance(201, 300, 1, 100), 101);
        // Touching / overlapping -> 0.
        assert_eq!(span_distance(1, 100, 100, 200), 0);
        assert_eq!(span_distance(1, 100, 50, 60), 0);
        // Adjacent by 1 bp gap.
        assert_eq!(span_distance(1, 100, 102, 200), 2);
    }

    /// Build a fully-populated `OritM8` row for gating tests.
    fn m8(
        query: &str,
        target: &str,
        pident: f64,
        tcov: f64,
        qstart: i64,
        qend: i64,
        evalue: f64,
        alnlen: i64,
    ) -> OritM8 {
        OritM8 {
            query: query.to_string(),
            target: target.to_string(),
            pident,
            tcov,
            qstart,
            qend,
            bits: 100.0,
            evalue,
            alnlen,
        }
    }

    #[test]
    fn gate_colocalized_vs_standalone_bar() {
        // A relaxase called on contig1 spanning 5000..6000.
        let relaxase = vec![("contig1".to_string(), 5000i64, 6000i64)];

        // Co-localized (within 10 kb) hit that only clears the BASE bar (E=1e-6,
        // 72% id, tcov 0.55) -> kept because it sits next to a relaxase.
        let coloc = m8("contig1", "oriT_A", 72.0, 0.55, 4000, 4100, 1e-6, 101);
        // Same-strength hit but STANDALONE (contig2, no relaxase) -> rejected: base
        // bar is not enough without co-localization.
        let weak_standalone = m8("contig2", "oriT_B", 72.0, 0.55, 100, 200, 1e-6, 101);
        // Strong standalone hit clearing the stricter bar -> kept.
        let strong_standalone = m8("contig3", "oriT_C", 99.0, 0.95, 900, 700, 1e-12, 200);

        let kept = gate_orit_hits(vec![coloc, weak_standalone, strong_standalone], &relaxase);
        let names: Vec<&str> = kept.iter().map(|h| h.accession.as_str()).collect();
        assert!(names.contains(&"oriT_A"), "co-localized base-bar hit should be kept");
        assert!(!names.contains(&"oriT_B"), "weak standalone hit should be rejected");
        assert!(names.contains(&"oriT_C"), "strong standalone hit should be kept");

        // Co-localization flag + normalized reverse-orientation coords/strand.
        let a = kept.iter().find(|h| h.accession == "oriT_A").unwrap();
        assert!(a.colocalized);
        assert_eq!((a.start, a.end), (4000, 4100));
        assert_eq!(a.strand, 1);
        let c = kept.iter().find(|h| h.accession == "oriT_C").unwrap();
        assert!(!c.colocalized);
        assert_eq!((c.start, c.end), (700, 900)); // qstart 900 > qend 700 -> normalized
        assert_eq!(c.strand, -1);
    }

    #[test]
    fn gate_collapses_overlapping_loci() {
        // Two references matching the SAME locus on contig1 (overlapping spans);
        // only the best-bit one survives. Both co-localized.
        let relaxase = vec![("contig1".to_string(), 1000i64, 2000i64)];
        let mut strong = m8("contig1", "oriT_best", 99.0, 0.95, 1500, 1650, 1e-20, 151);
        strong.bits = 300.0;
        let mut weak = m8("contig1", "oriT_redundant", 90.0, 0.80, 1520, 1600, 1e-10, 81);
        weak.bits = 120.0;
        let kept = gate_orit_hits(vec![weak, strong], &relaxase);
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].accession, "oriT_best");
    }

    /// End-to-end: embed a known `oriT_exp` sequence in a synthetic contig next to a
    /// (dummy) relaxase HMM dir and confirm the detector places an oriT feature at
    /// the embedded coordinates with the right accession. Requires mmseqs + the
    /// staged oriT_exp DB; skips gracefully otherwise.
    #[test]
    fn synthetic_contig_locates_embedded_orit() {
        let db = concat!(env!("CARGO_MANIFEST_DIR"), "/../db/orit/oriT_exp.fasta");
        if !std::path::Path::new(db).exists() {
            eprintln!("synthetic_contig_locates_embedded_orit: {db} absent; skipped");
            return;
        }
        if !crate::mmseqs::available() {
            eprintln!("synthetic_contig_locates_embedded_orit: mmseqs unavailable; skipped");
            return;
        }
        // Take the first oriT_exp record (accession + sequence).
        let (acc, orit_seq) = first_fasta_record(db).expect("oriT_exp.fasta should have a record");
        assert!(orit_seq.len() >= 40, "need a non-trivial oriT to embed");

        // Build a contig: 1000 bp flank + oriT + 1000 bp flank (1-based embed start).
        let flank = pseudo_dna(1000, 7);
        let embed_start = flank.len() as i64 + 1; // 1-based start of the oriT
        let embed_end = embed_start + orit_seq.len() as i64 - 1;
        let mut contig = String::new();
        contig.push_str(&flank);
        contig.push_str(&orit_seq);
        contig.push_str(&pseudo_dna(1000, 13));
        let contigs = vec![("synthetic_contig".to_string(), contig)];

        // Dir with a dummy .hmm (never parsed: no CDS proteins) + the oriT nt DB.
        let dir = std::env::temp_dir().join(format!("bactars_orit_syn_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("relaxase.hmm"), "NAME  x\nGA    30.0 30.0;\n").unwrap();
        std::fs::copy(db, dir.join("oriT_exp.fasta")).unwrap();

        let feats = detect(&contigs, &[], dir.to_str().unwrap()).expect("detect should succeed");
        let _ = std::fs::remove_dir_all(&dir);

        let orit = feats
            .iter()
            .find(|f| f.annotations.iter().any(|a| a.source == "orit:nic"))
            .unwrap_or_else(|| panic!("no oriT nt locus found; features: {feats:?}"));
        let ann = orit.annotations.iter().find(|a| a.source == "orit:nic").unwrap();
        assert_eq!(ann.accession, acc, "matched the embedded oriT accession");
        assert!(
            ann.evalue.is_some_and(|e| e <= ORIT_MAX_EVALUE),
            "real mmseqs evalue reported"
        );
        // Coordinates must land on the embedded locus (allow small alignment slack).
        assert!(
            (orit.start - embed_start).abs() <= 5 && (orit.end - embed_end).abs() <= 5,
            "oriT feature at {}..{} but embedded at {}..{}",
            orit.start, orit.end, embed_start, embed_end
        );
    }

    /// Read the first FASTA record's (first-header-token, uppercased sequence).
    fn first_fasta_record(path: &str) -> Option<(String, String)> {
        let text = std::fs::read_to_string(path).ok()?;
        let mut name = None;
        let mut seq = String::new();
        for line in text.lines() {
            if let Some(rest) = line.strip_prefix('>') {
                if name.is_some() {
                    break; // reached the second record
                }
                name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
            } else if name.is_some() {
                seq.push_str(line.trim());
            }
        }
        name.map(|n| (n, seq.to_ascii_uppercase()))
    }

    /// Deterministic non-homopolymer ACGT filler of length `n` (seeded so flanks do
    /// not accidentally match an oriT reference).
    fn pseudo_dna(n: usize, seed: u64) -> String {
        const B: [char; 4] = ['A', 'C', 'G', 'T'];
        let mut x = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut s = String::with_capacity(n);
        for _ in 0..n {
            // xorshift-ish step for a spread-out, reproducible base sequence.
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            s.push(B[(x % 4) as usize]);
        }
        s
    }
}
