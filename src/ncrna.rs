//! ncRNA detection via infernox (in-process Infernal cmsearch).
//!
//! Drives `infernox::cm_search::FaithfulSearcher` — the byte-parity-verified
//! faithful cmsearch pipeline (F1 MSV → F3 Forward → F3b bias → F4/F5 glocal env →
//! F6 banded CYK → F7 banded Inside + null3 + E-value) — over every covariance
//! model in a (possibly multi-model) Infernal `.cm` database, on both strands.
//!
//! The published `infernox-infernal` crate exposes a *single-model* reader, so a
//! multi-model DB (e.g. an Rfam bacterial subset) is split into per-model records
//! here and each is loaded + searched in turn. Each `.cm` record is a CM block
//! followed by its trailing p7 HMM filter block (both `//`-terminated), which
//! `cm_file_read_from_reader` consumes in one call.

use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use infernox::cm_file::cm_file_read_from_reader_opt;
use infernox::cm_search::{FaithfulConfig, FaithfulHit, FaithfulSearcher, ModelCutoff};
use std::collections::HashMap;
use rayon::prelude::*;
use std::io::Cursor;
use std::panic::{self, AssertUnwindSafe};

/// Detect ncRNAs by scanning `contigs` (`(name, seq)` pairs) against every CM in
/// the Infernal `.cm` database at `cm_db`. Emits `Feature{Ncrna|Rrna|Tmrna}`.
///
/// A full Rfam-scale DB (hundreds/thousands of models) against a whole genome is
/// inherently slow (cmsearch); each model's p7 HMM filter does the heavy lifting.
///
/// Each model is scanned inside `catch_unwind` as a defensive measure: a
/// malformed or uncalibrated model in a large Rfam DB is skipped with a warning
/// rather than aborting the whole scan. (Models must be read in GLOBAL config —
/// see `cm_file_read_from_reader_opt(.., false)` below — which is what
/// `FaithfulSearcher::new` expects; the localizing reader would make the CP9
/// psi sums trip an internal assertion.)
pub fn detect_ncrna(contigs: &[(String, String)], cm_db: &str) -> Result<Vec<Feature>, String> {
    let text = std::fs::read_to_string(cm_db).map_err(|e| format!("{cm_db}: {e}"))?;
    let models = split_models(&text);
    if models.is_empty() {
        return Err(format!("no INFERNAL1 models in {cm_db}"));
    }

    // All contigs are searched by every model in one `search` call; `seq_idx`
    // maps a hit back to its contig.
    let seqs: Vec<&str> = contigs.iter().map(|(_, s)| s.as_str()).collect();
    // Bakta's ncRNA-gene / ncRNA-region secondary E-value filter (features/nc_rna.py
    // HIT_EVALUE = 1E-4): discard --cut_tc survivors whose E-value exceeds this.
    const HIT_EVALUE: f64 = 1e-4;
    // Bakta's rRNA secondary COVERAGE filter (features/r_rna.py HIT_COVERAGE = 0.3):
    // discard --cut_tc rRNA survivors covering < 30% of the model consensus length.
    const HIT_COVERAGE: f64 = 0.3;

    // Threshold ncRNA hits the Bakta / Rfam-recommended way: use each model's
    // curated TC (trusted) bit-score cutoff (C `cmscan --cut_tc`) as the primary
    // reporting/inclusion threshold, then apply an E-value <= HIT_EVALUE secondary
    // filter (Bakta discards hits above 1e-4). `model_cutoff = Tc` resolves per
    // model to its own TC at search time; a model without a TC line falls back to
    // `e_report` (kept at HIT_EVALUE so those are held to the same E-value bar).
    // `--rfam` (C cm_pipeline.c:479-490): the Rfam-scan filter preset Bakta runs
    // (`cmscan --cut_tc --nohmmonly --rfam`). It raises the HMM filter thresholds
    // (F1=0.06, F2=F2b=0.02, F3..F5=0.0002, F6=0.0001) so far fewer windows reach
    // the expensive CM DP — the intended acceleration for large (Rfam-scale) DBs.
    // Matches Bakta's invocation exactly; byte-parity-verified against C `--rfam`.
    let cfg = FaithfulConfig {
        e_report: HIT_EVALUE,
        model_cutoff: Some(ModelCutoff::Tc),
        rfam: true,
        ..FaithfulConfig::default()
    };

    // NOTE: we deliberately do NOT install a process-global `panic::set_hook` here.
    // `panic::catch_unwind` (below) already isolates a per-model panic into a
    // skip-with-warning without any hook change; mutating the global hook is a
    // process-wide side effect that races with other concurrently-running stages
    // (and with any hook the host set), and could swallow/duplicate their output.
    // A caught model panic may still print the default backtrace line — an
    // acceptable, honest trade for not racing global state.

    // Parse + configure + scan each model INDEPENDENTLY, in parallel over the
    // rayon pool. infernox's per-model search is itself window-parallel, but on a
    // bacterial genome that inner parallelism is shallow (few windows) and each
    // model's serial setup (CM parse + p7 filter build) dominates — so processing
    // models one at a time pinned the run to ~3 cores. Fanning the model loop out
    // fills the pool. `par_iter().map(..).collect()` preserves model order, and the
    // per-contig hit numbering is assigned in a deterministic serial post-pass
    // below, so the emitted features are byte-identical to the sequential version.
    //
    // NOTE (task #44 investigation): a memory-bounded batched FLAT window pool was
    // implemented and byte-parity-verified (`FaithfulSearcher::search_many_batched`),
    // but MEASURED SLOWER than this per-model loop at 4/8/16 threads (t150 144-model
    // DB: +24/+38/+29 %; full 1104-model MG1655 sweep: +1 % wall, +0.8 GB RSS). The
    // flat pool must build every searcher in a batch BEFORE any window search starts
    // (windows can't be enumerated without a built searcher), so the slow giant-model
    // build (23S/16S QDB+CP9) gates the batch's search phase — whereas this loop
    // overlaps each model's build with other models' window searches (plus rayon's
    // nested window parallelism inside `.search()`), which wins as threads scale. The
    // flat-pool driver stays in the engine for FEW-model / pre-built-searcher callers
    // (where a single giant would otherwise gate the wall and the build barrier is
    // negligible); bactars' 1000+-model sweep keeps this per-model path.
    type ModelScan = Option<(String, String, FeatureKind, i32, Vec<FaithfulHit>)>;
    // LPT (longest-processing-time-first) scheduling: launch the largest models
    // first. A few giant CMs (23S/16S rRNA, tmRNA) parallelize shallowly on a
    // single-contig genome, so if they start late they gate the tail; starting
    // them first overlaps them with the long tail of small models. This is an
    // EXECUTION-ORDER change only — results are restored to original model order
    // below, so the per-contig numbering (and thus the byte output) is unchanged.
    // CM text length is a cheap proxy for model size / DP cost.
    let mut sched: Vec<usize> = (0..models.len()).collect();
    sched.sort_by_key(|&i| std::cmp::Reverse(models[i].len()));
    let scanned: Vec<(usize, ModelScan)> = sched
        .par_iter()
        .map(|&idx| {
            let chunk = &models[idx];
            let seqs_ref = &seqs;
            let cfg_ref = &cfg;
            let outcome = panic::catch_unwind(AssertUnwindSafe(
                || -> Result<(String, String, FeatureKind, i32, Vec<FaithfulHit>), String> {
                    // Read in GLOBAL config (do_localize=false); `FaithfulSearcher::new`
                    // expects an already-global CM (it derives `cm_global` for psi / QDB
                    // bands); the localizing reader mixes in local begin/end transitions,
                    // so split-state psi no longer sums to 1.0 and `cm_expected_state_
                    // occupancy` panics — see the note on `detect_ncrna`.
                    let cm = cm_file_read_from_reader_opt(Cursor::new(chunk.as_bytes()), false)
                        .map_err(|e| format!("{e:?}"))?;
                    let name = cm.name.clone();
                    let acc = cm.acc.clone().unwrap_or_else(|| "-".to_string());
                    let desc = cm.desc.clone().unwrap_or_default();
                    let kind = classify(&name, &desc);
                    let searcher = FaithfulSearcher::new(cm)?;
                    let clen = searcher.cm().clen; // model consensus length, for rRNA coverage
                    Ok((name, acc, kind, clen, searcher.search(seqs_ref, cfg_ref)))
                },
            ));
            let scan = match outcome {
                Ok(Ok(v)) => Some(v),
                // Uncalibrated model (no p7 filter) or parse error — expected, quiet.
                Ok(Err(_)) => None,
                Err(_) => {
                    eprintln!(
                        "ncrna: skipping model {} (infernox internal panic)",
                        model_label(chunk)
                    );
                    None
                }
            };
            (idx, scan)
        })
        .collect();

    // Restore original model order so the numbering post-pass is deterministic.
    let mut per_model: Vec<ModelScan> = (0..models.len()).map(|_| None).collect();
    for (idx, scan) in scanned {
        per_model[idx] = scan;
    }

    // Deterministic serial post-pass over the order-preserved results: assign the
    // per-contig ncRNA numbering exactly as the old sequential loop did.
    let mut per_contig: HashMap<String, usize> = HashMap::new();
    let mut feats = Vec::new();
    let mut skipped = 0usize;
    for entry in &per_model {
        let (name, acc, kind, clen, hits) = match entry {
            Some(v) => (&v.0, &v.1, v.2, v.3, &v.4),
            None => {
                skipped += 1;
                continue;
            }
        };

        for h in hits {
            // Bakta secondary filters on --cut_tc survivors, by class:
            //  - rRNA: coverage = hit length / model consensus length >= 0.3 (r_rna.py)
            //  - other ncRNA: E-value <= 1e-4 (nc_rna.py / nc_rna_region.py)
            if matches!(kind, FeatureKind::Rrna) {
                let hit_len = (h.start - h.stop).abs() + 1;
                let coverage = hit_len as f64 / (clen.max(1)) as f64;
                if coverage < HIT_COVERAGE {
                    continue;
                }
            } else if h.evalue > HIT_EVALUE {
                continue;
            }
            let contig = &contigs[h.seq_idx].0;
            // For reverse-complement hits infernox reports start > stop; store a
            // normalised (start <= end) span and record the strand separately.
            let (start, end) = if h.start <= h.stop {
                (h.start, h.stop)
            } else {
                (h.stop, h.start)
            };
            let n = per_contig.entry(contig.clone()).or_insert(0);
            *n += 1;
            feats.push(Feature {
                kind,
                contig: contig.clone(),
                id: format!("{}_ncrna{}", contig, n),
                start,
                end,
                strand: if h.in_rc { -1 } else { 1 },
                aa: None,
                partial5: false,
                partial3: false,
                annotations: vec![Annotation {
                    source: "infernox".to_string(),
                    accession: acc.clone(),
                    name: name.clone(),
                    score: h.score,
                    // Real infernox cmsearch E-value for this ncRNA/rRNA hit.
                    evalue: Some(h.evalue),
                    ref_len: None,
                }],
                func: Functional::default(),
            });
        }
    }

    if skipped > 0 {
        eprintln!("ncrna: scanned {} models, {skipped} skipped", models.len());
    }
    Ok(feats)
}

/// A short human label (`NAME` / `ACC`) for a model record, for warnings.
fn model_label(chunk: &str) -> String {
    let mut name = None;
    let mut acc = None;
    for line in chunk.lines().take(12) {
        let s = line.trim();
        if name.is_none() && s.starts_with("NAME") {
            name = s.split_whitespace().nth(1).map(String::from);
        } else if acc.is_none() && s.starts_with("ACC") {
            acc = s.split_whitespace().nth(1).map(String::from);
        }
    }
    match (name, acc) {
        (Some(n), Some(a)) => format!("{n} ({a})"),
        (Some(n), None) => n,
        (None, Some(a)) => a,
        (None, None) => "<unknown>".to_string(),
    }
}

/// Split a (possibly multi-model) Infernal `.cm` file into per-model text records.
/// Each record starts at an `INFERNAL1` header line and runs up to (not including)
/// the next one, so it contains one CM block plus its trailing p7 HMM filter.
fn split_models(text: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut cur = String::new();
    for line in text.lines() {
        if line.starts_with("INFERNAL1") && !cur.is_empty() {
            records.push(std::mem::take(&mut cur));
        }
        cur.push_str(line);
        cur.push('\n');
    }
    if !cur.trim().is_empty() {
        records.push(cur);
    }
    records
}

/// True if `needle` occurs in `hay` as a whole word — i.e. bounded on both sides
/// by a non-alphanumeric character (or a string edge) — rather than as an
/// incidental substring. Both arguments are expected already lowercased. Multi-word
/// needles (e.g. `"ribosomal rna"`) are matched with their internal spaces literal
/// and only their outer edges boundary-checked.
fn contains_word(hay: &str, needle: &str) -> bool {
    if needle.is_empty() {
        return false;
    }
    let hb = hay.as_bytes();
    let mut from = 0;
    while let Some(pos) = hay[from..].find(needle) {
        let start = from + pos;
        let end = start + needle.len();
        let before_ok = start == 0 || !hb[start - 1].is_ascii_alphanumeric();
        let after_ok = end >= hb.len() || !hb[end].is_ascii_alphanumeric();
        if before_ok && after_ok {
            return true;
        }
        from = start + 1;
    }
    false
}

/// Classify a covariance model by its name / description into a feature kind.
/// rRNA and tmRNA families get their own kind; everything else is generic ncRNA.
///
/// The Infernal `.cm` record exposes only free-text `NAME`/`ACC`/`DESC` — there is
/// no structured Rfam family-type (`TP`) or clan (`CL`) field to prefer — so we
/// classify from that free text. Matching is on WORD BOUNDARIES (not a raw
/// substring `contains`) so that a family which merely *mentions* rRNA/tmRNA in its
/// description (e.g. a snoRNA that "guides rRNA methylation") is NOT misfiled as the
/// rRNA/tmRNA itself. If a structured type/clan field becomes available on the CM
/// record, prefer it over this free-text heuristic.
fn classify(name: &str, desc: &str) -> FeatureKind {
    let hay = format!("{} {}", name, desc).to_lowercase();
    if contains_word(&hay, "rrna") || contains_word(&hay, "ribosomal rna") {
        FeatureKind::Rrna
    } else if contains_word(&hay, "tmrna") || contains_word(&hay, "transfer-messenger") {
        FeatureKind::Tmrna
    } else {
        FeatureKind::Ncrna
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn contains_word_boundaries() {
        assert!(contains_word("5s_rrna", "rrna"));
        assert!(contains_word("ssu_rrna_bacteria", "rrna"));
        assert!(contains_word("small subunit ribosomal rna", "ribosomal rna"));
        assert!(contains_word("bacterial tmrna", "tmrna"));
        // Incidental substrings inside a longer alnum run must NOT match.
        assert!(!contains_word("xrrnay", "rrna"));
        assert!(!contains_word("serrnase", "rrna"));
        assert!(!contains_word("attmrna1x", "tmrna"));
    }

    #[test]
    fn classify_rrna_tmrna_and_generic() {
        assert_eq!(classify("5S_rRNA", "5S ribosomal RNA"), FeatureKind::Rrna);
        assert_eq!(classify("LSU_rRNA_bacteria", ""), FeatureKind::Rrna);
        assert_eq!(classify("tmRNA", "transfer-messenger RNA"), FeatureKind::Tmrna);
        // A plain riboswitch/leader stays generic ncRNA.
        assert_eq!(
            classify("TPP", "TPP riboswitch (THI element)"),
            FeatureKind::Ncrna
        );
        // Regression: an incidental SUBSTRING match ("rrna" buried inside another
        // token) must NOT trip the rRNA classifier — the whole point of the
        // word-boundary fix. The old `contains("rrna")` misfiled this as rRNA.
        assert_eq!(
            classify("Xrrnay_leader", "riboswitch controlling the xrrnay operon"),
            FeatureKind::Ncrna
        );
    }
}
