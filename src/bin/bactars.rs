//! bactars CLI — first vertical slice: CDS + HMM annotation.
//!
//!   bactars <genome.fna> [--pfam <db.hmm>] [--ncbifams <db.hmm>] [--antifam <db.hmm>]
//!
//! Pfam/NCBIfams/AntiFam databases are searched with the gathering cutoff
//! (`--cut_ga`), matching Bakta. Emits a simple TSV of features to stdout.

use bactars::hmm::Cutoff;
use bactars::pipeline::{run_full, Config, HmmDb};
use bactars::{output, resolve};

fn main() {
    // Subcommand: `bactars setup-db --out DIR [--only ...] [--force]` fetches the
    // reference databases and exits before the annotation arg loop.
    let argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.first().map(String::as_str) == Some("setup-db") {
        run_setup_db(&argv[1..]);
        return;
    }

    let mut genome: Option<String> = None;
    let mut hmm_dbs: Vec<HmmDb> = Vec::new();
    let mut trna_models_dir: Option<String> = None;
    let mut detect_tmrna = false;
    let mut is_db_dir: Option<String> = None;
    let mut ncrna_db: Option<String> = None;
    let mut detect_crispr = false;
    let mut detect_oric = false;
    let mut detect_gaps = false;
    let mut detect_tandem = false;
    let mut annotate_signalpep = false;
    let mut annotate_localize = false;
    let mut annotate_gene_score = false;
    let mut refine_starts = false;
    let mut profile_fast = false;
    let mut profile_full = false;
    let mut meta_dir: Option<String> = None;
    let mut kofam_dir: Option<String> = None;
    let mut psc_dir: Option<String> = None;
    let mut amr_variant_dir: Option<String> = None;
    let mut vfdb_dir: Option<String> = None;
    let mut plasmid_dir: Option<String> = None;
    let mut species_dir: Option<String> = None;
    let mut db_bundle: Option<String> = None;
    let mut detect_pseudogenes = false;
    let mut detect_sorf = false;
    let mut integron_dir: Option<String> = None;
    let mut prophage_dir: Option<String> = None;
    let mut orit_dir: Option<String> = None;
    let mut gff3_out: Option<String> = None;
    let mut gbff_out: Option<String> = None;
    let mut threads: usize = 0;

    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--threads" | "-t" => {
                let v = require_arg(&arg, it.next());
                threads = v.parse().unwrap_or_else(|_| {
                    eprintln!("--threads needs a non-negative integer (0 = all cores)");
                    std::process::exit(2);
                });
            }
            "--pfam" => push_db(&mut hmm_dbs, "pfam", Cutoff::GatheringGa, it.next()),
            "--ncbifams" => push_db(&mut hmm_dbs, "ncbifams", Cutoff::GatheringGa, it.next()),
            "--antifam" => push_db(&mut hmm_dbs, "antifam", Cutoff::GatheringGa, it.next()),
            // AMR the Bakta way: an AMR HMM (NCBI AMRFinderPlus) searched on the
            // same rustyhmmer path. AMRFinderPlus HMMs carry trusted cutoffs (TC).
            "--amr" => push_db(&mut hmm_dbs, "amr", Cutoff::TrustedTc, it.next()),
            "--trna" => match it.next() {
                Some(d) => trna_models_dir = Some(d),
                None => {
                    eprintln!("--trna needs a trnascan-rs data/models directory");
                    std::process::exit(2);
                }
            },
            "--tmrna" => detect_tmrna = true,
            "--is-db" => match it.next() {
                Some(d) => is_db_dir = Some(d),
                None => {
                    eprintln!("--is-db needs an isscan IS-DB directory (needs mmseqs on PATH)");
                    std::process::exit(2);
                }
            },
            "--ncrna" => match it.next() {
                Some(d) => ncrna_db = Some(d),
                None => {
                    eprintln!("--ncrna needs an Infernal .cm database (Rfam subset)");
                    std::process::exit(2);
                }
            },
            "--crispr" => detect_crispr = true,
            "--oric" => detect_oric = true,
            "--gaps" => detect_gaps = true,
            "--tandem" => detect_tandem = true,
            "--signalpep" => annotate_signalpep = true,
            "--compartment" => annotate_localize = true,
            "--gene-score" => annotate_gene_score = true,
            "--refine-starts" => refine_starts = true,
            // Profiles: `--fast` = every auto stage EXCEPT the learned ML models
            // (signalpep/compartment/gene-score); `--full` = `--fast` plus those ML
            // models. DB-backed stages self-skip when their reference DB is absent
            // (supply it via --db or an explicit flag). `--all` is the old name for
            // `--full`, kept as a deprecated alias.
            "--fast" => profile_fast = true,
            "--full" => profile_full = true,
            "--all" => {
                eprintln!("note: --all was renamed to --full; treating as --full");
                profile_full = true;
            }
            "--meta" => meta_dir = Some(require_arg(&arg, it.next())),
            "--kofam" => kofam_dir = Some(require_arg(&arg, it.next())),
            "--psc" => psc_dir = Some(require_arg(&arg, it.next())),
            "--amr-variant" => amr_variant_dir = Some(require_arg(&arg, it.next())),
            "--vfdb" => vfdb_dir = Some(require_arg(&arg, it.next())),
            "--plasmid" => plasmid_dir = Some(require_arg(&arg, it.next())),
            "--species" => species_dir = Some(require_arg(&arg, it.next())),
            // Unified bundle: resolve every unset DB path by convention (see below).
            "--db" => db_bundle = Some(require_arg(&arg, it.next())),
            // --pseudo = the mmseqs alignment detector (needs mmseqs + --ncbifams;
            // warns + skips otherwise, never runs a low-precision proxy).
            "--pseudo" => detect_pseudogenes = true,
            "--sorf" => detect_sorf = true,
            "--integron" => integron_dir = Some(require_arg(&arg, it.next())),
            "--prophage" => prophage_dir = Some(require_arg(&arg, it.next())),
            "--orit" => orit_dir = Some(require_arg(&arg, it.next())),
            "--gff3" => gff3_out = Some(require_arg(&arg, it.next())),
            "--gbff" => gbff_out = Some(require_arg(&arg, it.next())),
            "-h" | "--help" => {
                eprintln!("usage: bactars <genome.fna> [--db BUNDLE] [--fast] [--full] [--pfam D] [--ncbifams D] [--antifam D] [--amr D] [--trna MODELS_DIR] [--tmrna] [--is-db DIR] [--ncrna CM_DB] [--crispr] [--oric] [--gaps] [--tandem] [--signalpep] [--compartment] [--gene-score] [--meta DIR] [--kofam DIR] [--psc DIR] [--amr-variant DIR] [--vfdb DIR] [--plasmid DIR] [--species DIR] [--sorf] [--pseudo] [--integron DIR] [--prophage DIR] [--orit DIR] [--threads N] [--gff3 FILE] [--gbff FILE]");
                eprintln!("       --fast enables every non-ML stage (oric/gaps/tandem/crispr/tmrna + sorf/pseudo); DB-backed stages self-skip when their DB is absent.");
                eprintln!("       --full = --fast plus the learned ML models (signalpep/compartment/gene-score); pure-Rust but conv-heavy (~3x wall time).");
                eprintln!("       --db BUNDLE resolves every unset DB path by the bactars-db bundle layout (explicit flags override).");
                eprintln!("       bactars setup-db --out DIR [--only pfam,ncbifams,antifam,amr,rfam,isdb,vfdb,plasmidfinder,gtdb_ssu,integron,phage,orit_pfam,meta] [--force]");
                return;
            }
            _ => {
                if genome.is_none() {
                    genome = Some(arg);
                } else {
                    eprintln!("unexpected argument: {arg}");
                    std::process::exit(2);
                }
            }
        }
    }

    let genome = match genome {
        Some(g) => g,
        None => {
            eprintln!("usage: bactars <genome.fna> [--pfam D] [--ncbifams D] [--antifam D]");
            std::process::exit(2);
        }
    };

    // Profiles. `--fast` (also implied by `--full`) turns on every auto stage EXCEPT
    // the learned ML models: the motif/structure detectors (oric/gaps/tandem/crispr/tmrna)
    // plus the DB-backed sorf/pseudo stages (which self-skip when their reference DB is
    // absent). `--full` additionally turns on the ML models (signalpep/compartment/
    // gene-score) — these are pure-Rust and DB-free but conv-heavy, so they roughly triple
    // wall time; that is why they are opt-in via --full rather than --fast. An explicit
    // per-stage flag and a profile are equivalent.
    if profile_fast || profile_full {
        detect_tmrna = true;
        detect_crispr = true;
        detect_oric = true;
        detect_gaps = true;
        detect_tandem = true;
        detect_sorf = true;
        detect_pseudogenes = true;
    }
    if profile_full {
        annotate_signalpep = true;
        annotate_localize = true;
        annotate_gene_score = true;
    }

    // Unified `--db <bundle>`: fill every DB path NOT already given explicitly, by
    // the `bactars-db` bundle layout convention. An explicit flag always wins; a
    // convention path is only used when present, so a LIGHT bundle (no kofam/psc)
    // simply leaves those stages disabled. See bactars-db/README.md.
    if let Some(bundle) = &db_bundle {
        let root = std::path::Path::new(bundle);
        let add_hmm = |dbs: &mut Vec<HmmDb>, label: &str, rel: &str, cutoff: Cutoff| {
            if !dbs.iter().any(|d| d.label == label) {
                let f = root.join(rel);
                if f.exists() {
                    dbs.push(HmmDb {
                        label: label.to_string(),
                        path: f.to_string_lossy().into_owned(),
                        cutoff,
                    });
                }
            }
        };
        add_hmm(&mut hmm_dbs, "ncbifams", "ncbifams/ncbifams.hmm", Cutoff::GatheringGa);
        add_hmm(&mut hmm_dbs, "antifam", "antifam/AntiFam.hmm", Cutoff::GatheringGa);
        add_hmm(&mut hmm_dbs, "amr", "amr/AMR.LIB", Cutoff::TrustedTc);
        add_hmm(&mut hmm_dbs, "pfam", "pfam/Pfam-A.hmm", Cutoff::GatheringGa);

        let set_dir = |slot: &mut Option<String>, rel: &str| {
            if slot.is_none() {
                let f = root.join(rel);
                if f.exists() {
                    *slot = Some(f.to_string_lossy().into_owned());
                }
            }
        };
        set_dir(&mut trna_models_dir, "trna_models");
        set_dir(&mut ncrna_db, "rfam/Rfam.cm");
        set_dir(&mut meta_dir, "meta");
        set_dir(&mut kofam_dir, "kofam");
        set_dir(&mut psc_dir, "psc");
        set_dir(&mut amr_variant_dir, "amr");
        set_dir(&mut vfdb_dir, "vfdb");
        set_dir(&mut plasmid_dir, "plasmidfinder");
        set_dir(&mut species_dir, "gtdb_ssu");
        set_dir(&mut is_db_dir, "isdb");
        set_dir(&mut integron_dir, "integron");
        set_dir(&mut prophage_dir, "phage");
        set_dir(&mut orit_dir, "orit_pfam");
    }

    // --fast/--full only flip stage switches; the database-backed stages they enable
    // (ncbifams functional HMM, ncRNA/Rfam, IS elements, mmseqs pseudogene, sORF HMM
    // filter) still need their reference DBs. A profile without any DB source silently
    // collapses to the DB-free subset (not even functional HMM naming), so warn loudly.
    if (profile_fast || profile_full) && db_bundle.is_none() && hmm_dbs.is_empty() {
        let p = if profile_full { "--full" } else { "--fast" };
        eprintln!("warning: {p} was given without --db (or explicit DB flags): the");
        eprintln!("         database-backed stages it enables (ncbifams functional HMM,");
        eprintln!("         ncrna, is, pseudogene, sorf) will be SKIPPED or degraded —");
        eprintln!("         only database-free detectors will run. Pass `--db <bundle>`");
        eprintln!("         (or explicit --ncbifams/--ncrna/... flags) for a full annotation.");
    }

    let config = Config {
        trans_table: 11,
        hmm_dbs,
        trna_models_dir,
        detect_tmrna,
        is_db_dir,
        ncrna_db,
        detect_crispr,
        detect_oric,
        detect_gaps,
        detect_tandem,
        annotate_signalpep,
        annotate_localize,
        annotate_gene_score,
        meta_dir,
        kofam_dir,
        psc_dir,
        amr_variant_dir,
        vfdb_dir,
        plasmid_dir,
        species_dir,
        detect_pseudogenes,
        full_mode: profile_full,
        detect_sorf,
        integron_dir,
        prophage_dir,
        orit_dir,
        threads,
        refine_starts,
    };

    // Fail fast: verify the genome file and every configured DB path exists and is
    // non-empty BEFORE running, so a missing/empty DB is a clear up-front error
    // rather than a silent mid-run degradation (or a subunit's opaque failure).
    if let Err(e) = validate_inputs(&genome, &config) {
        eprintln!("error: {e}");
        std::process::exit(2);
    }

    // Bound the rayon global pool (rustyhmmer/infernox parallelism) to the same
    // thread count. 0 leaves rayon at its default (all logical cores).
    if threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    }

    let run_out = match run_full(&genome, &config) {
        Ok(o) => o,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    let plasmids = run_out.plasmids;
    let species = run_out.species;

    // Resolve overlapping calls (RNA masks CDS, dedup multi-model ncRNA) and sort.
    let features = resolve::resolve(run_out.features);

    // Structured output files, if requested.
    if let Some(path) = &gff3_out {
        if let Err(e) = write_gff3_file(path, &genome, &features) {
            eprintln!("error writing {path}: {e}");
            std::process::exit(1);
        }
        eprintln!("wrote GFF3: {path}");
    }
    if let Some(path) = &gbff_out {
        if let Err(e) = write_gbff_file(path, &genome, &features) {
            eprintln!("error writing {path}: {e}");
            std::process::exit(1);
        }
        eprintln!("wrote GBFF: {path}");
    }

    let n_cds = features.iter().filter(|f| f.aa.is_some()).count();
    let n_annot = features.iter().filter(|f| !f.annotations.is_empty()).count();

    println!("# contig\tid\tstart\tend\tstrand\tkind\tbest_hit\tscore\tevalue");
    for f in &features {
        let (hit, score, evalue) = match f.best_annotation() {
            Some(a) => (format!("{} ({})", a.name, a.accession), a.score, a.evalue),
            None => ("-".to_string(), 0.0, None),
        };
        // `None` (a score-only feature, or no annotation) serialises as `.`; a real
        // search E-value keeps the scientific-notation column. Never emit `NaN`.
        let evalue = match evalue {
            Some(e) => format!("{e:.2e}"),
            None => ".".to_string(),
        };
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:?}\t{}\t{:.1}\t{}",
            f.contig, f.id, f.start, f.end, f.strand, f.kind, hit, score, evalue
        );
    }
    let n_total = features.len();
    eprintln!("{n_total} features ({n_cds} CDS, {n_annot} annotated)");

    // CheckM-lite genome quality (single-copy ribosomal markers), part of every run.
    let qc = &run_out.qc;
    eprintln!("{}", bactars::qc::summary_line(qc));
    if !qc.missing.is_empty() {
        eprintln!("  QC missing markers ({}): {}", qc.missing.len(), qc.missing.join(","));
    }
    if !qc.duplicated.is_empty() {
        eprintln!(
            "  QC duplicated markers ({}, contamination): {}",
            qc.duplicated.len(),
            qc.duplicated.join(",")
        );
    }

    // Genome/contig-level metadata from the opt-in typing stages.
    let n_vf = features
        .iter()
        .filter(|f| f.func.note.iter().any(|n| n.starts_with("virulence factor:")))
        .count();
    if n_vf > 0 {
        eprintln!("virulence factors: {n_vf} CDS with VFDB evidence");
    }
    if !plasmids.is_empty() {
        eprintln!("plasmid replicons: {} detected", plasmids.len());
        for p in &plasmids {
            eprintln!(
                "  {} : {} ({:.1}% id, {}..{})",
                p.contig, p.replicon_type, p.identity, p.start, p.end
            );
        }
    }
    if let Some(sp) = &species {
        eprintln!(
            "species: {} ({:.1}% 16S identity; {})",
            sp.species, sp.identity, sp.taxonomy
        );
    }
}

fn run_setup_db(args: &[String]) {
    use bactars::setup_db::{setup_db, SetupConfig};
    use std::path::PathBuf;

    let mut out_dir: Option<PathBuf> = None;
    let mut only: Vec<String> = Vec::new();
    let mut force = false;
    let mut yes = false;
    let mut build_psc = false;
    let mut build_kofam = false;
    let mut threads = 0usize;
    let mut isdb_url: Option<String> = std::env::var("BACTARS_ISDB_URL").ok();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--out" | "-o" => out_dir = it.next().map(PathBuf::from),
            "--only" => {
                if let Some(list) = it.next() {
                    only = list
                        .split(',')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();
                }
            }
            // Mode presets matching the annotation modes. `--fast` fetches the full
            // standard DB set (everything --fast/--full annotation needs); `--full`
            // adds the heavy FULL-tier UniRef90 PSC naming DB (built in-process).
            "--fast" => only = bactars::setup_db::mode_keys(),
            "--full" => {
                only = bactars::setup_db::mode_keys();
                build_psc = true;
                build_kofam = true;
            }
            "--threads" | "-t" => {
                threads = it.next().and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "--force" => force = true,
            "-y" | "--yes" => yes = true,
            // Source for the rust-ise IS profile-DB tarball (overrides env).
            // Preferred spelling `--is-db-url`; `--isdb-url` kept as a back-compat alias.
            "--is-db-url" | "--isdb-url" => isdb_url = it.next().cloned(),
            "-h" | "--help" => {
                eprintln!("usage: bactars setup-db --out DIR [--fast | --full] [--only pfam,ncbifams,antifam,amr,amr_variant,rfam,isdb,vfdb,plasmidfinder,gtdb_ssu,integron,phage,orit_pfam,meta,trna,psc,kofam] [--threads N] [--force] [--yes] [--isdb-url URL]");
                eprintln!("       --fast : all standard annotation DBs (light+standard tier, incl trna models)");
                eprintln!("       --full : --fast DBs + BUILD the FULL-tier UniRef90 PSC (~32 GB dl + ~194 GB index) + KEGG KOfam naming DBs");
                eprintln!("       --only psc / --only kofam build ONLY that FULL-tier DB");
                eprintln!("       (isdb also reads BACTARS_ISDB_URL from the environment)");
                return;
            }
            _ => {
                eprintln!("unexpected setup-db argument: {a}");
                std::process::exit(2);
            }
        }
    }
    // `psc` and `kofam` are BUILD targets, not download specs — pull them out of
    // `only` and flip their build flags (so `--only psc`/`--only kofam` and `--full`
    // all trigger the in-process build).
    if let Some(pos) = only.iter().position(|k| k == "psc") {
        only.remove(pos);
        build_psc = true;
    }
    if let Some(pos) = only.iter().position(|k| k == "kofam") {
        only.remove(pos);
        build_kofam = true;
    }
    if build_psc || build_kofam {
        eprintln!("note: --full / --only psc,kofam will BUILD the FULL-tier naming DBs in-process (heavy: PSC ~230 GB, KOfam ~5 GB).");
    }
    let out_dir = match out_dir {
        Some(d) => d,
        None => {
            eprintln!("setup-db needs --out DIR");
            std::process::exit(2);
        }
    };
    if let Err(e) = setup_db(&SetupConfig {
        out_dir,
        only,
        force,
        yes,
        isdb_url,
        build_psc,
        build_kofam,
        threads,
    }) {
        eprintln!("setup-db error: {e}");
        std::process::exit(1);
    }
}

/// Up-front input validation: the genome file must exist and be non-empty, and
/// every DB path the run was configured with (HMM db files + the `--is-db`/`--trna`/
/// `--vfdb`/`--plasmid`/… directories) must exist and be non-empty. Returns the
/// first problem found so the run fails fast with a clear message instead of a
/// silent mid-run skip or an opaque subunit error.
fn validate_inputs(genome: &str, config: &Config) -> Result<(), String> {
    let gp = std::path::Path::new(genome);
    if !gp.is_file() {
        return Err(format!("genome file not found: {genome}"));
    }
    if std::fs::metadata(gp).map(|m| m.len()).unwrap_or(0) == 0 {
        return Err(format!("genome file is empty: {genome}"));
    }

    // A configured DB path must exist and be non-empty (dir: has entries; file:
    // non-zero length). Anything else means the tier would silently no-op or crash.
    fn check(label: &str, path: &str) -> Result<(), String> {
        let p = std::path::Path::new(path);
        if !p.exists() {
            return Err(format!("{label} path does not exist: {path}"));
        }
        let empty = if p.is_dir() {
            std::fs::read_dir(p)
                .map(|mut it| it.next().is_none())
                .unwrap_or(true)
        } else {
            std::fs::metadata(p).map(|m| m.len() == 0).unwrap_or(true)
        };
        if empty {
            return Err(format!("{label} is empty: {path}"));
        }
        Ok(())
    }

    for db in &config.hmm_dbs {
        check(&format!("--{}", db.label), &db.path)?;
    }
    let dirs: [(&str, &Option<String>); 13] = [
        ("--trna", &config.trna_models_dir),
        ("--is-db", &config.is_db_dir),
        ("--ncrna", &config.ncrna_db),
        ("--meta", &config.meta_dir),
        ("--kofam", &config.kofam_dir),
        ("--psc", &config.psc_dir),
        ("--amr-variant", &config.amr_variant_dir),
        ("--vfdb", &config.vfdb_dir),
        ("--plasmid", &config.plasmid_dir),
        ("--species", &config.species_dir),
        ("--integron", &config.integron_dir),
        ("--prophage", &config.prophage_dir),
        ("--orit", &config.orit_dir),
    ];
    for (label, opt) in dirs {
        if let Some(p) = opt {
            check(label, p)?;
        }
    }
    Ok(())
}

fn require_arg(flag: &str, v: Option<String>) -> String {
    match v {
        Some(s) => s,
        None => {
            eprintln!("{flag} needs an argument");
            std::process::exit(2);
        }
    }
}

fn write_gff3_file(
    path: &str,
    genome: &str,
    features: &[bactars::Feature],
) -> std::io::Result<()> {
    // Submission-grade GFF3 embeds the sequence (##FASTA) + ##sequence-region, so
    // re-read the genome for its contigs (same as GBFF).
    let contigs: Vec<(String, Vec<u8>)> = bactars::fasta::read_genome(genome)
        .map_err(std::io::Error::other)?
        .into_iter()
        .map(|c| (c.name, c.seq))
        .collect();
    let w = std::io::BufWriter::new(std::fs::File::create(path)?);
    output::write_gff3_submission(w, features, &contigs, output::DEFAULT_LOCUS_PREFIX)
}

fn write_gbff_file(
    path: &str,
    genome: &str,
    features: &[bactars::Feature],
) -> std::io::Result<()> {
    // GBFF embeds the ORIGIN sequence, so re-read the genome for its contigs.
    let contigs: Vec<(String, Vec<u8>)> = bactars::fasta::read_genome(genome)
        .map_err(std::io::Error::other)?
        .into_iter()
        .map(|c| (c.name, c.seq))
        .collect();
    let w = std::io::BufWriter::new(std::fs::File::create(path)?);
    output::write_gbff(w, features, &contigs)
}

fn push_db(dbs: &mut Vec<HmmDb>, label: &str, cutoff: Cutoff, path: Option<String>) {
    match path {
        Some(p) => dbs.push(HmmDb {
            label: label.to_string(),
            path: p,
            cutoff,
        }),
        None => {
            eprintln!("--{label} needs a path");
            std::process::exit(2);
        }
    }
}
