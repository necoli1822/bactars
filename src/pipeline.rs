//! The orchestrator: scatter genome → subunits, gather into `Feature`s.
//!
//! First vertical slice: rustygal (CDS) → rustyhmmer (HMM annotation). Extensible
//! to oxagorn / trnascan-rs (RNA), infernox (ncRNA). AMR is done the Bakta way —
//! an AMR HMM DB on the same rustyhmmer path (see `HmmDb`), not a separate tool.

use crate::cds::predict_cds;
use crate::fasta::read_genome;
use crate::feature::{Annotation, Feature, FeatureKind, Functional};
use crate::hmm::{self, Cutoff};
use crate::plasmid::{self, PlasmidHit};
use crate::species::{self, SpeciesCall};
use crate::{
    amr_variant, crispr, gap, integron, is_elements, kofam, ncrna, oric, orit, prophage, psc,
    gene_score, localize, pseudogene_align, signalox_sp, sorf, tandem, tmrna, trna, vfdb, xref,
};
use std::collections::HashMap;

/// One HMM database to search CDS proteins against, with its cutoff and a label
/// recorded on each annotation (`source = "rustyhmmer:{label}"`).
pub struct HmmDb {
    pub label: String,
    pub path: String,
    pub cutoff: Cutoff,
}

/// Pipeline configuration.
pub struct Config {
    /// Genetic code for CDS translation. `None` = auto-detect (compare NCBI
    /// tables 11/4/25 by total gene score, prodigal-meta style — picks table 4 for
    /// Mollicutes, 25 for Gracilibacteria/SR1, 11 otherwise). `Some(n)` forces
    /// table `n`. See [`crate::cds::predict_cds`].
    pub trans_table: Option<i32>,
    /// HMM databases to annotate CDS with (Pfam / NCBIfams / AntiFam / ...).
    pub hmm_dbs: Vec<HmmDb>,
    /// trnascan-rs `data/models` directory. `Some` enables bacterial (`-B`)
    /// tRNA detection; `None` skips it.
    pub trna_models_dir: Option<String>,
    /// Enable tmRNA detection via oxagorn (pure-Rust, no external DB).
    pub detect_tmrna: bool,
    /// isscan IS-DB directory. `Some` enables rust-ise IS detection (needs
    /// `mmseqs` on PATH); `None` skips it.
    pub is_db_dir: Option<String>,
    /// Infernal `.cm` database (single- or multi-model, e.g. an Rfam bacterial
    /// subset). `Some` enables infernox ncRNA cmsearch; `None` skips it.
    pub ncrna_db: Option<String>,
    /// Enable CRISPR-array detection via the pure-Rust MinCED port (no DB).
    pub detect_crispr: bool,
    /// Enable oriC (replication origin) detection — GC-skew + DnaA-box (no DB).
    pub detect_oric: bool,
    /// Enable assembly-gap (run-of-N) annotation (no DB).
    pub detect_gaps: bool,
    /// Enable tandem-repeat detection (no DB).
    pub detect_tandem: bool,
    /// Enable signal-peptide / cleavage / TM annotation on CDS (signalox, no DB).
    pub annotate_signalpep: bool,
    /// Enable subcellular-localization annotation on CDS (lotrs, no DB).
    pub annotate_localize: bool,
    /// Enable learned gene-model quality scoring on CDS (balrox, no DB; non-destructive).
    pub annotate_gene_score: bool,
    /// bactars-db `meta/` dir (hmm_PGAP.tsv + pfam2go). `Some` enriches CDS with
    /// gene symbol / EC / GO / curated product from the HMM hit accession.
    pub meta_dir: Option<String>,
    /// KOfam bundle dir (`kofam_prok.hmm` + `ko_list` + `ko00001.keg`). `Some`
    /// assigns KEGG KO / EC / pathway / COG to each CDS.
    pub kofam_dir: Option<String>,
    /// Native PSC dir (`psc_db` mmseqs + `names.tsv`). `Some` names still-unnamed
    /// CDS by UniRef90 protein similarity (memory-heavy: mmaps a ~194GB index).
    pub psc_dir: Option<String>,
    /// AMRFinderPlus DB dir (AMRProt.fa + AMRProt-mutation.tsv). `Some` annotates
    /// AMR point-mutation (variant) resistance on CDS.
    pub amr_variant_dir: Option<String>,
    /// VFDB dir (`VFDB_setB_pro.fas.gz`). `Some` adds virulence-factor evidence
    /// (`/note` + `/inference`) to CDS by protein similarity. Additive only.
    pub vfdb_dir: Option<String>,
    /// PlasmidFinder dir (`plasmidfinder_db.tar.gz`). `Some` types plasmid
    /// replicons on the contigs (contig-level metadata, not a Feature).
    pub plasmid_dir: Option<String>,
    /// GTDB SSU dir (`bac120_ssu_reps.fna.gz` [+ `ar53_ssu_reps.fna.gz`]). `Some`
    /// identifies the genome's species from its 16S rRNA feature(s).
    pub species_dir: Option<String>,
    /// Enable pseudogene detection (truncation vs reference + split-gene).
    pub detect_pseudogenes: bool,
    /// `--full` profile active. Gates the heavy full-DB pseudogene detection
    /// (the UniRef90/PSC reference-length truncation, "lever A") to `--full` — it
    /// needs the large `--psc` DB and is a thorough-tier signal, so it does not run
    /// in the lighter `--fast` profile even when `--psc` is supplied.
    pub full_mode: bool,
    /// Enable sORF (small-ORF) detection (short ORFs + HMM homology filter).
    pub detect_sorf: bool,
    /// Integron DB dir (`IntI.hmm` + `attc_4.cm` [+ `phage-int.hmm`]). `Some`
    /// enables IntegronFinder-style integron detection (rustyhmmer + infernox).
    pub integron_dir: Option<String>,
    /// Prophage hallmark HMM dir (`phage_hallmark.hmm`). `Some` enables
    /// phage-hallmark region-clustering prophage detection (rustyhmmer).
    pub prophage_dir: Option<String>,
    /// User-supplied relaxase HMM dir. `Some` enables oriT / mobilizable-element
    /// (relaxase-presence) detection. Gated + license-sensitive — see `orit`.
    pub orit_dir: Option<String>,
    /// Worker threads for oxagorn / rust-ise (0 = all cores).
    pub threads: usize,
    /// Homology-guided start-codon refinement (PGAP-style): trim ab-initio
    /// over-extended CDS starts to the ncbifams full-length-family envelope anchor.
    /// Requires `--ncbifams`. Off by default (see [`crate::start_refine`]).
    pub refine_starts: bool,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            trans_table: None,
            hmm_dbs: Vec::new(),
            trna_models_dir: None,
            detect_tmrna: false,
            is_db_dir: None,
            ncrna_db: None,
            detect_crispr: false,
            detect_oric: false,
            detect_gaps: false,
            detect_tandem: false,
            annotate_signalpep: false,
            annotate_localize: false,
            annotate_gene_score: false,
            meta_dir: None,
            kofam_dir: None,
            psc_dir: None,
            amr_variant_dir: None,
            vfdb_dir: None,
            plasmid_dir: None,
            species_dir: None,
            detect_pseudogenes: false,
            full_mode: false,
            detect_sorf: false,
            integron_dir: None,
            prophage_dir: None,
            orit_dir: None,
            threads: 0,
            refine_starts: false,
        }
    }
}

/// Full pipeline output: the annotated feature list plus the genome/contig-level
/// metadata produced by the opt-in typing stages (plasmid replicons + species).
/// [`run`] returns only `features` for backwards compatibility; callers that want
/// the metadata call [`run_full`].
pub struct RunOutput {
    /// Every annotated feature (CDS + RNA + structural).
    pub features: Vec<Feature>,
    /// Plasmid replicon types detected on the contigs (empty unless `plasmid_dir`).
    pub plasmids: Vec<PlasmidHit>,
    /// Genome species assignment from 16S (None unless `species_dir` + a 16S hit).
    pub species: Option<SpeciesCall>,
    /// CheckM-lite genome-quality report (single-copy ribosomal marker
    /// completeness / contamination / disrupted-core), always computed.
    pub qc: crate::qc::QcReport,
}

/// Run the pipeline on a genome FASTA, returning annotated features. Thin wrapper
/// over [`run_full`] that discards the plasmid/species metadata — kept so existing
/// callers (and tests) compile unchanged.
pub fn run(genome_path: &str, config: &Config) -> Result<Vec<Feature>, String> {
    Ok(run_full(genome_path, config)?.features)
}

/// Run the pipeline on a genome FASTA, returning the annotated features together
/// with the opt-in plasmid-replicon + species metadata ([`RunOutput`]).
pub fn run_full(genome_path: &str, config: &Config) -> Result<RunOutput, String> {
    let contigs = read_genome(genome_path)?;
    if contigs.is_empty() {
        return Err(format!("no sequences in {genome_path}"));
    }

    // Optional-tier degradation tracker: every enrichment/typing stage that is
    // configured but fails is recorded here (short label + one-word reason) and
    // surfaced as a single concise end-of-run summary line, so a silently degraded
    // run (missing DB, absent mmseqs, ...) is visible instead of only scrolling by.
    let mut skipped_tiers: Vec<String> = Vec::new();

    // --- CDS prediction (rustygal) ---
    let cdss = predict_cds(&contigs, config.trans_table)?;
    let mut features: Vec<Feature> = cdss
        .iter()
        .map(|c| Feature {
            kind: FeatureKind::Cds,
            contig: c.contig.clone(),
            id: c.id.clone(),
            start: c.start,
            end: c.end,
            strand: c.strand,
            aa: Some(c.aa.clone()),
            partial5: c.partial5,
            partial3: c.partial3,
            annotations: Vec::new(),
            func: Functional::default(),
        })
        .collect();

    // --- HMM functional annotation (rustyhmmer) ---
    // feature index -> best ncbifams envelope anchor, for start-codon refinement.
    let mut start_anchors: HashMap<usize, crate::start_refine::Anchor> = HashMap::new();
    if !config.hmm_dbs.is_empty() {
        // id → feature index, for attaching hits back onto CDS features.
        let idx: HashMap<String, usize> = features
            .iter()
            .enumerate()
            .map(|(i, f)| (f.id.clone(), i))
            .collect();
        let proteins: Vec<(String, String)> = features
            .iter()
            .filter_map(|f| f.aa.as_ref().map(|aa| (f.id.clone(), aa.clone())))
            .collect();

        for db in &config.hmm_dbs {
            let hits = hmm::annotate(&proteins, &db.path, db.cutoff)?;
            let is_ncbifams = db.label == "ncbifams";
            for h in hits {
                if let Some(&fi) = idx.get(h.target_name.as_str()) {
                    // Capture the ncbifams (full-length equivalog) envelope as the
                    // best-scoring start-refinement anchor for this CDS, before the
                    // hit's String fields are moved into the Annotation.
                    if is_ncbifams {
                        let a = (h.env_from, h.env_to, h.model_len, h.seq_score);
                        start_anchors
                            .entry(fi)
                            .and_modify(|e| {
                                if a.3 > e.3 {
                                    *e = a;
                                }
                            })
                            .or_insert(a);
                    }
                    features[fi].annotations.push(Annotation {
                        source: format!("rustyhmmer:{}", db.label),
                        accession: h.query_acc,
                        name: h.query_name,
                        score: h.seq_score,
                        evalue: Some(h.seq_evalue),
                        ref_len: None,
                    });
                }
            }
        }
    }

    // --- homology-guided start-codon refinement (PGAP-style, opt-in) ---
    // Trim ab-initio over-extended 5' ends to the ncbifams envelope anchor. Runs
    // before every downstream tier so naming/pseudogene/sORF see refined starts.
    if config.refine_starts && !start_anchors.is_empty() {
        let st = crate::start_refine::refine_starts(&contigs, &mut features, &start_anchors);
        eprintln!(
            "start-refine: trimmed {}/{} CDS starts to ncbifams anchor ({} aa removed)",
            st.trimmed, st.examined, st.total_res_trimmed
        );
    }

    // --- xref enrichment: gene symbol / EC / GO / curated product from the HMM
    // hit accession (NCBIfams hmm_PGAP.tsv + pfam2go). Pure lookup, no new search.
    if let Some(meta_dir) = &config.meta_dir {
        match xref::XrefTable::load(meta_dir) {
            Ok(table) => xref::enrich(&mut features, &table),
            Err(e) => {
                eprintln!("xref: skipping enrichment ({e})");
                skipped_tiers.push(format!("xref({})", short_reason(&e)));
            }
        }
    }

    // --- KOfam KO / EC / pathway / COG (KofamScan method) ---
    if let Some(kofam_dir) = &config.kofam_dir {
        if let Err(e) = kofam::annotate(&mut features, kofam_dir) {
            eprintln!("kofam: skipping KO assignment ({e})");
            skipped_tiers.push(format!("kofam({})", short_reason(&e)));
        }
    }

    // --- PSC protein-similarity naming (native UniRef90 mmseqs). Runs LAST among
    // the CDS-naming tiers so it only names CDS still lacking a product. ---
    if let Some(psc_dir) = &config.psc_dir {
        if let Err(e) = psc::annotate(&mut features, psc_dir, config.threads) {
            eprintln!("psc: skipping UniRef naming ({e})");
            skipped_tiers.push(format!("psc({})", short_reason(&e)));
        }
    }

    // --- AMR point-mutation (variant) resistance on CDS (AMRFinderPlus catalog). ---
    if let Some(amr_dir) = &config.amr_variant_dir {
        if let Err(e) = amr_variant::annotate(&mut features, amr_dir, config.threads) {
            eprintln!("amr_variant: skipping ({e})");
            skipped_tiers.push(format!("amr_variant({})", short_reason(&e)));
        }
    }

    // --- VFDB virulence-factor evidence on CDS (protein similarity). Additive:
    // adds a /note + /inference, never overwrites the product. Runs after the
    // CDS-naming tiers so it only decorates already-named CDS. ---
    if let Some(vfdb_dir) = &config.vfdb_dir {
        if let Err(e) = vfdb::annotate(&mut features, vfdb_dir, config.threads) {
            eprintln!("vfdb: skipping ({e})");
            skipped_tiers.push(format!("vfdb({})", short_reason(&e)));
        }
    }

    // --- tRNA detection (trnascan-rs, bacterial -B) ---
    if let Some(models_dir) = &config.trna_models_dir {
        let calls = trna::detect_trnas(genome_path, models_dir)?;
        features.extend(trna::trna_features(calls));
    }

    // --- tmRNA detection (oxagorn, pure-Rust) ---
    if config.detect_tmrna {
        features.extend(tmrna::detect_tmrna(genome_path, config.threads)?);
    }

    // --- IS-element detection (rust-ise, needs mmseqs + IS DB) ---
    if let Some(is_db) = &config.is_db_dir {
        features.extend(is_elements::detect_is(genome_path, is_db, config.threads)?);
    }

    // ncRNA (infernox), CRISPR (MinCED port), and the mobile-element detectors
    // (integron / prophage / oriT) all work on in-memory contig sequences and/or
    // the CDS features gathered so far; build the `(name, seq)` view once if any
    // is enabled. These run AFTER CDS naming + IS detection (they need CDS).
    if config.ncrna_db.is_some()
        || config.detect_crispr
        || config.integron_dir.is_some()
        || config.prophage_dir.is_some()
        || config.orit_dir.is_some()
    {
        let seqs: Vec<(String, String)> = contigs
            .iter()
            .map(|c| (c.name.clone(), String::from_utf8_lossy(&c.seq).into_owned()))
            .collect();

        // --- ncRNA detection (infernox cmsearch vs an Rfam .cm DB) ---
        if let Some(cm_db) = &config.ncrna_db {
            features.extend(ncrna::detect_ncrna(&seqs, cm_db)?);
        }

        // --- CRISPR-array detection (pure-Rust MinCED port, no DB) ---
        if config.detect_crispr {
            features.extend(crispr::detect_crispr(&seqs)?);
        }

        // --- Integron detection (intI HMM on CDS + attC CM on contigs) ---
        // Each detector filters `features` to CDS internally, so the ncRNA/CRISPR
        // features appended above are ignored; the immutable borrow ends before
        // `extend`, so borrowing `&features` then mutating is fine.
        if let Some(dir) = &config.integron_dir {
            let new = match integron::detect(&seqs, &features, dir) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("integron: skipping ({e})");
                    skipped_tiers.push(format!("integron({})", short_reason(&e)));
                    Vec::new()
                }
            };
            features.extend(new);
        }

        // --- Prophage detection (phage-hallmark HMM region clustering) ---
        if let Some(dir) = &config.prophage_dir {
            let new = match prophage::detect(&seqs, &features, dir) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("prophage: skipping ({e})");
                    skipped_tiers.push(format!("prophage({})", short_reason(&e)));
                    Vec::new()
                }
            };
            features.extend(new);
        }

        // --- oriT / mobilizable-element detection (relaxase HMM; gated) ---
        if let Some(dir) = &config.orit_dir {
            let new = match orit::detect(&seqs, &features, dir) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("orit: skipping ({e})");
                    skipped_tiers.push(format!("orit({})", short_reason(&e)));
                    Vec::new()
                }
            };
            features.extend(new);
        }
    }

    // --- Rfam-type reclassification: split infernox ncRNA hits into rRNA / tRNA /
    // regulatory_region (riboswitch/leader) / ncRNA by their Rfam family type. ---
    if let Some(meta_dir) = &config.meta_dir {
        reclassify_ncrna(&mut features, meta_dir);
    }

    // --- sORF detection: short ORFs below the gene-caller floor kept when they
    // get an HMM homology hit OR (when a PSC DB is configured) a well-covered
    // UniRef90 protein-similarity hit. Uses the first HMM DB as the filter. ---
    if config.detect_sorf {
        if let Some(hmm_db) = config.hmm_dbs.first().map(|d| d.path.clone()) {
            let sorfs = sorf::detect(
                &contigs,
                &features,
                &hmm_db,
                config.psc_dir.as_deref(),
                config.threads,
            );
            features.extend(sorfs);
        } else {
            eprintln!("sorf: no HMM db (--ncbifams/--pfam) for the homology filter; skipped");
        }
    }

    // --- Signal peptide / cleavage / TM annotation via the learned signalox model ---
    // (replaces the rule-based signalpep::annotate; held-out SP-any P~0.95 vs 0.46).
    // Runs before the structural detectors so it only sees real CDS features.
    if config.annotate_signalpep {
        signalox_sp::annotate(&mut features);
    }

    // --- Subcellular localization annotation via the learned lotrs model ---
    if config.annotate_localize {
        localize::annotate(&mut features);
    }

    // --- Learned gene-model quality score via balrox (non-destructive: scores +
    //     flags low-scoring CDS as possible spurious ORFs; never removes a call) ---
    if config.annotate_gene_score {
        gene_score::annotate(&mut features);
    }

    // --- oriC (GC-skew + DnaA-box; uses CDS for dnaA proximity) ---
    // Compute against the CDS features gathered so far, then append.
    if config.detect_oric {
        let oric_feats = oric::detect(&contigs, &features);
        features.extend(oric_feats);
    }

    // --- Assembly gaps (run-of-N) and tandem repeats (DB-free) ---
    if config.detect_gaps {
        features.extend(gap::detect(&contigs));
    }
    if config.detect_tandem {
        features.extend(tandem::detect(&contigs));
    }

    // --- Pseudogene detection (truncation vs reference length + split-gene).
    // Runs after all CDS naming so each CDS's best annotation/accession is set. ---
    if config.detect_pseudogenes {
        // `--pseudo` IS the alignment-based detector (pseudofinder-style translated
        // search vs the ncbifams family consensus). It requires mmseqs on PATH + an
        // ncbifams HMM. When either is missing we WARN and SKIP — there is no
        // low-precision length/split proxy fallback (it was net-negative, ~7% precision,
        // and only ever presented near-noise as pseudogene calls; removed 2026-07-19).
        let has_ncbifams = config.hmm_dbs.iter().any(|d| d.label == "ncbifams");
        if has_ncbifams && crate::mmseqs::available() {
            let db = config
                .hmm_dbs
                .iter()
                .find(|d| d.label == "ncbifams")
                .expect("has_ncbifams checked");
            let cache = pseudogene_align::consensus_cache_path(&db.path);
            let t0 = std::time::Instant::now();
            match pseudogene_align::build_consensus_fasta(&db.path, &cache) {
                Ok(built) => {
                    let t_cons = t0.elapsed();
                    // Under --full with --psc, point the disruption + intergenic
                    // detectors at the BROAD UniRef90 psc_db instead of the narrow
                    // ncbifams consensus, so the split / frameshift / internal-stop and
                    // intergenic signals reach non-IS disrupted/uncalled genes (ncbifams
                    // covers only ~4-11% of non-IS pseudogenes). Falls back to the
                    // consensus when --psc is absent or not --full.
                    let psc_target: Option<String> = if config.full_mode {
                        config.psc_dir.as_ref().and_then(|d| {
                            let dbp = std::path::Path::new(d).join("psc_db");
                            dbp.with_extension("dbtype")
                                .exists()
                                .then(|| dbp.to_string_lossy().into_owned())
                        })
                    } else {
                        None
                    };
                    let ref_label = if psc_target.is_some() {
                        "UniRef90 psc_db (--full)"
                    } else {
                        "ncbifams consensus"
                    };
                    let t1 = std::time::Instant::now();
                    let flags = pseudogene_align::detect_mmseqs(
                        &contigs,
                        &features,
                        &cache,
                        psc_target.as_deref(),
                        config.threads,
                    );
                    // On the BROAD UniRef90 path, keep only IS/transposase disruption
                    // calls: the broad reference cleanly boosts IS pseudogenes (degraded
                    // transposases) but its non-IS disruption calls are near-noise
                    // (measured: +48 non-IS calls yielded +1 real RefSeq pseudogene) —
                    // non-IS truncation is already served precisely by lever A. The
                    // narrow ncbifams-consensus path keeps all calls (no such flood).
                    let broad = psc_target.is_some();
                    let mut n = 0usize;
                    for (idx, note) in flags {
                        if broad && !pseudogene_align::feat_is_transpos(&features[idx]) {
                            continue;
                        }
                        let f = &mut features[idx];
                        f.func.pseudogene = true;
                        if !f.func.note.contains(&note) {
                            f.func.note.push(note);
                        }
                        n += 1;
                    }
                    eprintln!(
                        "pseudogene: mmseqs-align detector — {} flagged (ref: {}; consensus {} [{} built] {:.1}s; detect {:.1}s)",
                        n, ref_label, cache, built, t_cons.as_secs_f64(), t1.elapsed().as_secs_f64()
                    );
                    // Intergenic 6-frame scan (Pseudofinder-style): recover genes
                    // the caller MISSED entirely — degraded past its threshold — by
                    // translating every inter-feature gap against the reference.
                    // Adds NEW pseudogene loci (the "0 ORFs called" disruption case).
                    // SKIPPED on the broad UniRef90 path: recovered gap ORFs have no
                    // product (→ all classed non-IS) and UniRef90's spurious weak
                    // homologies flood them (measured 10% precision, ~no real recall);
                    // it is only sound against the narrow curated ncbifams consensus.
                    if !broad {
                        let t2 = std::time::Instant::now();
                        let ig = pseudogene_align::detect_intergenic_mmseqs(
                            &contigs,
                            &features,
                            &cache,
                            psc_target.as_deref(),
                            config.threads,
                        );
                        let n_ig = ig.len();
                        features.extend(ig);
                        if n_ig > 0 {
                            eprintln!(
                                "pseudogene: +{n_ig} intergenic homology loci recovered (uncalled/degraded genes; {:.1}s)",
                                t2.elapsed().as_secs_f64()
                            );
                        }
                    }
                }
                Err(e) => {
                    eprintln!("warning: pseudogene consensus build failed ({e}); pseudogene detection SKIPPED");
                    skipped_tiers.push("pseudogene(consensus-build)".to_string());
                }
            }
        } else {
            let need = if !has_ncbifams { "an --ncbifams HMM" } else { "mmseqs on PATH" };
            eprintln!(
                "warning: --pseudo needs {need}; pseudogene detection SKIPPED (no low-precision \
                 proxy fallback)."
            );
            skipped_tiers.push(format!(
                "pseudogene(needs {})",
                if !has_ncbifams { "ncbifams" } else { "mmseqs" }
            ));
        }

        // Full-DB pseudogene detection (Pseudofinder-style truncation, "lever A"),
        // gated to `--full`. The PSC naming tier already searched each CDS against
        // UniRef90 and recorded the real homolog length in `Annotation.ref_len`; a CDS
        // much shorter than its real UniRef90 homolog is a truncated fragment. This is
        // the ONLY signal that reaches non-IS pseudogenes (measured ~32% non-IS
        // precision) because UniRef90 covers ~all proteins where the small ncbifams
        // family set does not — so it needs the heavy `--psc` DB and belongs in the
        // thorough `--full` profile, not `--fast`. Real protein length (not the
        // ncbifams consensus, which flags normal short members). Contig-edge partials
        // excluded. Fires only when BOTH `--full` and `--psc` (ref_len) are present.
        if config.full_mode {
        const TRUNC_REFLEN_FRAC: f64 = 0.55;
        let mut n_reflen = 0usize;
        for f in features.iter_mut() {
            if f.kind != FeatureKind::Cds || f.func.pseudogene || f.partial5 || f.partial3 {
                continue;
            }
            let aa = match &f.aa {
                Some(a) => a.len() as i64,
                None => continue,
            };
            if aa < 40 {
                continue;
            }
            let rref = f.annotations.iter().filter_map(|a| a.ref_len).max().unwrap_or(0) as i64;
            if rref < 60 {
                continue;
            }
            if (aa as f64) < TRUNC_REFLEN_FRAC * rref as f64 {
                let note = format!(
                    "truncated ({:.0}% of UniRef90 reference length)",
                    100.0 * aa as f64 / rref as f64
                );
                if !f.func.note.contains(&note) {
                    f.func.note.push(note);
                }
                f.func.pseudogene = true;
                n_reflen += 1;
            }
        }
        if n_reflen > 0 {
            eprintln!("pseudogene: +{n_reflen} ref-length truncations (PSC/UniRef90 real reference, --full)");
        }
        } // end config.full_mode (full-DB lever A)
        // NOTE: an embedded-ncbifams-length-DISTRIBUTION truncation lever was tried
        // here (distilled offline via GPU mmseqs vs UniRef90; see db/psc_distill/ +
        // benchmark/refseq_bench/pilot/sweep_trunc.py) to remove lever A's --psc
        // dependency for non-IS pseudogenes. It is a NEGATIVE result and was removed:
        // only 4-11% of RefSeq non-IS pseudogenes even carry an ncbifams family, so a
        // per-family length table is structurally blind to ~90% of them (a coverage
        // ceiling, not a threshold problem — the flagged set was disjoint from the
        // real pseudogenes at every cutoff). lever A works only because UniRef90
        // (--psc) covers ~all proteins. Distillation cannot replace per-CDS homolog
        // length. See memory bactars-start-refine-negative.

        // NOTE: the reference-free disruption signals (frameshift / internal-stop
        // read-through / split-gene fragments, in pseudogene::detect_reference_free)
        // were wired here and MEASURED — also a NEGATIVE result, removed. In dense
        // bacterial genomes chance adjacent same-strand / same-frame ORFs swamp the
        // signal: non-IS precision 0-17% (200+ flags/genome, e.g. 131 flags vs 18
        // real on GCF_009792355.1); split never fires. This is exactly the ~2-7%
        // precision that retired the old proxy detector. Confirms: usable non-IS
        // pseudogene detection needs a real per-CDS reference (UniRef90 via --psc,
        // lever A) — no DB-free / family-only / geometry-only signal reaches it.
        // See memory bactars-start-refine-negative.
    }

    // Defensive: normalize every feature's contig to its first whitespace token,
    // the canonical id used by `fasta::read_genome`. Subunits that carry the full
    // FASTA header (historically oxagorn tmRNA) would otherwise not match the
    // first-token contigs used elsewhere, breaking `resolve()` overlap dedup.
    for f in features.iter_mut() {
        if let Some(tok) = f.contig.split([' ', '\t']).next() {
            if tok.len() != f.contig.len() {
                f.contig = tok.to_string();
            }
        }
    }

    // --- Plasmid replicon typing (PlasmidFinder, nucleotide search on contigs).
    // Contig-level metadata, not a Feature. ---
    let plasmids = match &config.plasmid_dir {
        Some(dir) => match plasmid::detect(&contigs, dir, config.threads) {
            Ok(hits) => hits,
            Err(e) => {
                eprintln!("plasmid: skipping replicon typing ({e})");
                skipped_tiers.push(format!("plasmid({})", short_reason(&e)));
                Vec::new()
            }
        },
        None => Vec::new(),
    };

    // --- Species identification (GTDB SSU, 16S rRNA). Runs after ncRNA/rRNA
    // reclassification so the 16S features exist. Genome-level metadata. ---
    let species = match &config.species_dir {
        Some(dir) => match species::identify(&features, &contigs, dir, config.threads) {
            Ok(call) => call,
            Err(e) => {
                eprintln!("species: skipping ({e})");
                skipped_tiers.push(format!("species({})", short_reason(&e)));
                None
            }
        },
        None => None,
    };

    // Concise end-of-run degradation summary: one visible line naming every
    // configured tier that was skipped/failed, so a degraded run is diagnosable.
    if !skipped_tiers.is_empty() {
        eprintln!("[bactars] skipped tiers: {}", skipped_tiers.join(", "));
    }

    // CheckM-lite genome QC from the finalized annotations (single-copy ribosomal
    // markers already named by ncbifams). Cheap, always run, part of the main pipeline.
    let qc = crate::qc::compute(&features);

    Ok(RunOutput { features, plasmids, species, qc })
}

/// Reclassify infernox ncRNA features by Rfam family type, using
/// `<meta_dir>/rfam_type.tsv` (rfam_acc<TAB>rfam_type<TAB>so_type). A hit whose
/// SO type is rRNA/tRNA/regulatory_region is retyped; plain ncRNA is left as is.
/// Silently no-ops if the table is missing.
fn reclassify_ncrna(features: &mut [Feature], meta_dir: &str) {
    let path = std::path::Path::new(meta_dir).join("rfam_type.tsv");
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            // No table means every infernox ncRNA hit stays FeatureKind::Ncrna,
            // so rRNA/tRNA/regulatory families are mis-typed. Warn (don't fail) so
            // the mis-typing is diagnosable instead of a silent no-op.
            eprintln!(
                "[bactars] reclassify_ncrna: {} not readable ({e}); ncRNA hits will \
                 NOT be retyped to rRNA/tRNA/regulatory_region",
                path.display()
            );
            return;
        }
    };
    // RF##### -> so_type, and RF##### -> INSDC regulatory_class (regulatory only).
    let mut so: HashMap<String, FeatureKind> = HashMap::new();
    let mut regclass: HashMap<String, String> = HashMap::new();
    for line in text.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 || cols[0] == "rfam_acc" {
            continue;
        }
        let kind = match cols[2].trim() {
            "rRNA" => FeatureKind::Rrna,
            "tRNA" => FeatureKind::Trna,
            "regulatory_region" => FeatureKind::RegulatoryRegion,
            _ => continue, // "ncRNA" (default) — leave the feature as Ncrna
        };
        let acc = cols[0].trim().to_string();
        if kind == FeatureKind::RegulatoryRegion {
            regclass.insert(acc.clone(), regulatory_class_from_rfam_type(cols[1]));
        }
        so.insert(acc, kind);
    }
    for f in features.iter_mut() {
        if f.kind != FeatureKind::Ncrna {
            continue;
        }
        // Match the infernox ncRNA hit's Rfam accession (strip a `.N` version).
        let acc = f
            .annotations
            .iter()
            .find(|a| a.source == "infernox")
            .map(|a| a.accession.split('.').next().unwrap_or(&a.accession).to_string());
        if let Some(acc) = acc {
            if let Some(&kind) = so.get(&acc) {
                f.kind = kind;
                if kind == FeatureKind::RegulatoryRegion {
                    // NCBI requires a regulatory_class; default to `other` (+note)
                    // when the Rfam type maps to no INSDC controlled-vocab class.
                    f.func.regulatory_class = Some(
                        regclass
                            .get(&acc)
                            .cloned()
                            .unwrap_or_else(|| "other".to_string()),
                    );
                }
            }
        }
    }
}

/// Classify an optional-tier error into a one-word reason for the concise
/// end-of-run "skipped tiers" summary (the full error is already on its own
/// stderr line). Heuristic only — for a human-readable tag, not for control flow.
fn short_reason(e: &str) -> &'static str {
    let l = e.to_lowercase();
    if l.contains("mmseqs") {
        "mmseqs"
    } else if l.contains("no such file")
        || l.contains("not found")
        || l.contains("cannot find")
        || l.contains("missing")
        || l.contains("does not exist")
    {
        "db missing"
    } else {
        "failed"
    }
}

/// Map an Rfam `type_field` (e.g. `Cis-reg; riboswitch;`) to an INSDC
/// `regulatory_class` controlled-vocabulary token. Anything without a direct
/// class (thermoregulator, IRES, bare `Cis-reg;`) falls back to `other`, which
/// NCBI accepts only alongside a `/note` (emitted by the output writer).
fn regulatory_class_from_rfam_type(rfam_type: &str) -> String {
    let t = rfam_type.to_lowercase();
    if t.contains("riboswitch") {
        "riboswitch"
    } else if t.contains("leader") {
        // Cis-regulatory leaders (Thr_leader, His_leader, …) are transcriptional
        // attenuators.
        "attenuator"
    } else if t.contains("frameshift") {
        "recoding_stimulatory_region"
    } else {
        "other"
    }
    .to_string()
}
