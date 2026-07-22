//! `setup-db`: download and arrange the reference databases bactars needs.
//!
//! Fetches the HMM / CM / profile databases the pipeline consumes and lays them
//! out under a single directory, then prints the exact bactars flags that point
//! at them:
//!
//! ```text
//! <out>/pfam/Pfam-A.hmm             --pfam
//! <out>/ncbifams/ncbifams.hmm       --ncbifams
//! <out>/antifam/AntiFam.hmm         --antifam
//! <out>/amr/AMR.LIB                 --amr
//! <out>/amr/{AMRProt.fa,AMRProt-mutation.tsv} --amr-variant (point-mutation refs)
//! <out>/rfam/Rfam.cm                --ncrna
//! <out>/trna_models/{TRNA2-bact.cm,…} --trna      (tRNAscan-SE 2.0 CM bundle)
//! <out>/vfdb/VFDB_setA_pro.fas.gz   --vfdb       (light-tier auxiliary DBs)
//! <out>/plasmidfinder/plasmidfinder_db.tar.gz  --plasmid
//! <out>/gtdb_ssu/{bac120,ar53}_ssu_reps.fna.gz --species
//! <out>/integron/{IntI.hmm,attc_4.cm,…}         --integron
//! <out>/phage/phage_hallmark.hmm    --prophage
//! <out>/orit_pfam/{relaxase_pfam.hmm,oriT_exp.fasta} --orit
//! <out>/meta/{hmm_PGAP.tsv,pfam2go,ko00001.keg} --meta
//! <out>/isdb/…                      --is-db   (see note below)
//! ```
//!
//! The FULL-tier naming DBs are intentionally **not** fetched here: `kofam`
//! (2.5 GB) and especially `psc` (UniRef90, ~194 GB mmseqs index) are provisioned
//! separately — `psc` is *built* from UniRef90 by its own pipeline and is never
//! auto-downloaded. See `bactars-db/README.md` for the FULL-tier build.
//!
//! All I/O is pure-Rust and in-process: downloads use [`crate::util_io`]'s ureq
//! (rustls) client, and decompression uses flate2 (`gunzip`) / tar (safe
//! extraction). No host tools (`curl` / `wget` / `gunzip` / `tar` /
//! `sha256sum`) are spawned. Everything returns `Result<_, String>` and avoids
//! `unwrap` on any fallible path.
//!
//! # The IS database (`isdb`)
//!
//! rust-ise's union IS profile DB (ISOSDB ∪ ISfinder MMseqs2 profiles) is a
//! *versioned artifact calibrated together with the thresholds baked into the
//! matching rust-ise release* — the two are a locked pair, and there is no
//! stable canonical public URL to fetch it from. `setup-db` therefore fetches it
//! only when given an explicit source via `--isdb-url` (or the `BACTARS_ISDB_URL`
//! environment variable): it downloads the `.tar.gz`, extracts it under
//! `<out>/isdb/`, and locates the directory holding `manifest_union.tsv` to point
//! `--is-db` at. When no URL is supplied it prints guidance instead of failing
//! the whole run (see [`isdb_note`]).

use std::io::{BufRead, BufReader, IsTerminal, Write};
use std::path::{Path, PathBuf};

/// Compiled-in default source URL for the *fpc-carrying* rust-ise IS profile-DB
/// tarball (`rust-ise-isdb-fpc.tar.gz`, the one that ships `fpc/refset` so strict
/// FP-control actually filters host-recombinase false positives).
///
/// This is the published rust-ise release asset; `--isdb-url` and the
/// `BACTARS_ISDB_URL` environment variable override it. Because a non-empty
/// default is compiled in, IS-DB provisioning (both `setup-db --only isdb` and
/// the interactive prompt in `is_elements::detect_is`) is zero-config: when the
/// effective URL equals this default, the downloaded tarball is integrity-checked
/// against [`DEFAULT_ISDB_SHA256`] before extraction.
pub const DEFAULT_ISDB_URL: &str =
    "https://github.com/necoli1822/rust-ise/releases/download/v0.2.3/rust-ise-isdb-fpc.tar.gz";

/// Expected sha256 of the tarball at [`DEFAULT_ISDB_URL`]. Verified in-process
/// (pure-Rust sha2) *before* extraction whenever the effective URL is the
/// compiled-in default; on mismatch the tarball is rejected and not extracted.
/// A user-supplied `--isdb-url` override provides its own data, so this check is
/// skipped for overrides.
pub const DEFAULT_ISDB_SHA256: &str =
    "7fafe73dd4b795db61d921c9287b38dac95b698eb29af1649e8108bf712e5fb5";

/// Resolve the effective IS-DB tarball URL: an explicit override (from
/// `--isdb-url` / `BACTARS_ISDB_URL`) wins; otherwise the compiled-in
/// [`DEFAULT_ISDB_URL`] if it is non-empty. Returns `None` when neither is set
/// (the current default state), so callers can print guidance rather than fail.
pub fn effective_isdb_url(override_url: Option<&str>) -> Option<String> {
    match override_url {
        Some(u) if !u.is_empty() => Some(u.to_string()),
        _ if !DEFAULT_ISDB_URL.is_empty() => Some(DEFAULT_ISDB_URL.to_string()),
        _ => None,
    }
}

/// Download + extract the fpc-carrying IS-DB tarball from `url` into `dest_dir`
/// and return the directory to pass to `--is-db` (the one holding
/// `manifest_union.tsv`, next to `mmdb_union/` and `fpc/`).
///
/// This is the shared entry point reused by the interactive "download the
/// complete IS DB now?" prompt in [`crate::is_elements::detect_is`]; it wraps the
/// same downloader-detection + [`fetch_isdb`] pipeline `setup-db` uses. When
/// `force` is false an already-extracted DB under `dest_dir` is reused.
pub fn download_isdb(url: &str, dest_dir: &Path, force: bool) -> Result<PathBuf, String> {
    std::fs::create_dir_all(dest_dir)
        .map_err(|e| format!("cannot create {}: {e}", dest_dir.display()))?;
    fetch_isdb(url, dest_dir, force)
}

/// Configuration for [`setup_db`].
pub struct SetupConfig {
    /// Directory the databases are laid out under (created if absent).
    pub out_dir: PathBuf,
    /// Subset selector: any of `pfam`, `ncbifams`, `antifam`, `amr`, `amr_variant`,
    /// `rfam`, the light-tier auxiliary keys `vfdb`, `plasmidfinder`, `gtdb_ssu`,
    /// `integron`, `phage`, `orit_pfam`, `meta`, or `isdb`. Empty means "all".
    pub only: Vec<String>,
    /// Re-download even if the target file already exists.
    pub force: bool,
    /// Skip the interactive download-size confirmation prompt (assume "yes").
    /// Required for non-interactive stdin (piped / CI): without it, a run with no
    /// TTY refuses to download rather than silently pulling multiple GB.
    pub yes: bool,
    /// Optional source URL for the rust-ise IS profile-DB tarball (`.tar.gz`).
    /// When `isdb` is requested and this is set, the tarball is downloaded and
    /// extracted; when unset, [`isdb_note`] guidance is printed instead. The CLI
    /// populates it from `--isdb-url` or the `BACTARS_ISDB_URL` environment
    /// variable, since the DB is a release-locked artifact with no fixed URL.
    pub isdb_url: Option<String>,
    /// Build the FULL-tier UniRef90 PSC naming DB (pure-Rust port of build_psc.sh;
    /// the `--full` / `--only psc` path). Heavy (~32 GB download + ~194 GB mmseqs
    /// index). Off by default so a bare `setup-db` never triggers it.
    pub build_psc: bool,
    /// Build the FULL-tier KEGG KOfam naming DB (pure-Rust; the `--full` /
    /// `--only kofam` path). Downloads `ko_list.gz` + `profiles.tar.gz` from KEGG
    /// and concatenates the prokaryotic subset into `kofam_prok.hmm` (~2.4 GB
    /// download, ~5-6 GB extracted). Off by default so a bare `setup-db` never
    /// triggers it. No external binary required.
    pub build_kofam: bool,
    /// Thread count for the heavy mmseqs `createindex` in the PSC build (0 = 16).
    pub threads: usize,
}

/// How to turn a downloaded artifact into the final database file.
enum Fetch {
    /// Download straight to the final path (no decompression).
    Plain,
    /// Download a `.gz` then `gunzip` it to the final path.
    Gunzip,
    /// Download a `.tar.gz` then extract a single member to the final path.
    TarMember(&'static str),
}

/// One database bactars can consume.
struct DbSpec {
    /// Selector key (matches `SetupConfig::only`).
    key: &'static str,
    /// Subdirectory under `out_dir`.
    subdir: &'static str,
    /// Final on-disk file name inside `subdir`.
    filename: &'static str,
    /// Public source URL.
    url: &'static str,
    /// Post-download handling.
    fetch: Fetch,
    /// The bactars CLI flag that points at the finished file.
    flag: &'static str,
    /// Human-readable description for progress messages.
    desc: &'static str,
    /// Rough decompressed on-disk size in MB, used only for the pre-download
    /// size warning. A hard-coded estimate, not a live figure.
    approx_disk_mb: u64,
}

/// The downloadable databases, in setup order. URLs verified live 2026-07.
const SPECS: &[DbSpec] = &[
    DbSpec {
        key: "pfam",
        subdir: "pfam",
        filename: "Pfam-A.hmm",
        url: "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
        fetch: Fetch::Gunzip,
        flag: "--pfam",
        desc: "Pfam-A HMM library (EBI/InterPro)",
        approx_disk_mb: 1500,
    },
    DbSpec {
        key: "ncbifams",
        subdir: "ncbifams",
        filename: "ncbifams.hmm",
        // hmm_PGAP.LIB is the concatenated NCBIfam/TIGRFAM HMMER3 library PGAP
        // uses — a single, ready-to-search .hmm file (no decompression).
        url: "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB",
        fetch: Fetch::Plain,
        flag: "--ncbifams",
        desc: "NCBIfam/PGAP HMM library (NCBI)",
        approx_disk_mb: 2700,
    },
    DbSpec {
        key: "antifam",
        subdir: "antifam",
        filename: "AntiFam.hmm",
        url: "https://ftp.ebi.ac.uk/pub/databases/Pfam/AntiFam/current/Antifam.tar.gz",
        fetch: Fetch::TarMember("AntiFam.hmm"),
        flag: "--antifam",
        desc: "AntiFam spurious-ORF HMMs (EBI)",
        approx_disk_mb: 16,
    },
    DbSpec {
        key: "amr",
        subdir: "amr",
        filename: "AMR.LIB",
        url: "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/\
              AMRFinderPlus/database/latest/AMR.LIB",
        fetch: Fetch::Plain,
        flag: "--amr",
        desc: "AMRFinderPlus HMM library (NCBI)",
        approx_disk_mb: 102,
    },
    DbSpec {
        key: "rfam",
        subdir: "rfam",
        filename: "Rfam.cm",
        url: "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz",
        fetch: Fetch::Gunzip,
        flag: "--ncrna",
        desc: "Rfam covariance models (EBI)",
        approx_disk_mb: 1800,
    },
];

// ---------------------------------------------------------------------------
// Light-tier auxiliary databases
//
// Unlike [`SPECS`] (single HMM/CM model libraries with a per-file model count),
// these are FASTA / archive / table / concatenated-HMM DBs that the MGE and
// annotation-enrichment stages consume. Several ship more than one file, and the
// bactars flags that feed them (`--vfdb`, `--plasmid`, `--species`, `--integron`,
// `--prophage`, `--orit`, `--meta`) all take a *directory*, so `setup-db` lays
// each DB out under `<out>/<subdir>/` matching the `--db` bundle convention.
// ---------------------------------------------------------------------------

/// EBI InterPro HMM API endpoint prefix. `pfam_hmm_url` appends `/<ACC>?…`.
const PFAM_HMM_API: &str = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam";

/// Build the EBI InterPro API URL for a Pfam accession's (gzipped) HMM.
fn pfam_hmm_url(acc: &str) -> String {
    format!("{PFAM_HMM_API}/{acc}?annotation=hmm")
}

/// The 6 CC0 Pfam-A MOB relaxase families concatenated into `relaxase_pfam.hmm`
/// (license-clean oriT detector DB; see `db/orit_pfam/PROVENANCE.md`).
const ORIT_RELAXASE_PFAM: &[&str] =
    &["PF01076", "PF03389", "PF03432", "PF05713", "PF07514", "PF13814"];

/// oriTDB's curated `oriT_exp.fasta`: the 122 experimentally-validated
/// origin-of-transfer nucleotide sequences. Fetched verbatim into the SAME
/// `orit_pfam/` directory `--orit` resolves to, because [`crate::orit`] looks for
/// `oriT_exp.fasta` (fallback `oriT_all.fasta`) *inside* `--orit <DIR>` to run the
/// real oriT nucleotide-locus search (mmseqs `--search-type 3`). Without this file
/// the CC0 relaxase-HMM bundle ships only `relaxase_pfam.hmm`, so the nt-locus
/// stage is silently skipped and the detector falls back to the relaxase-presence
/// proxy alone.
///
/// Source: oriTDB2 (Shanghai Jiao Tong University, Ou lab); free for academic use.
/// Cite Li et al., "oriTDB: a database of the origin-of-transfer regions of
/// bacterial mobile genetic elements", Nucleic Acids Research 53(D1):D163 (2025).
/// The pinned sha256 in `bactars-db/manifest.tsv` (row `orit_nt`) matches this
/// asset byte-for-byte (verified against the shipped `db/orit/oriT_exp.fasta`).
const ORIT_NT_EXP_URL: &str =
    "https://bioinfo-mml.sjtu.edu.cn/oriTDB2/download/sequence/oriT_exp.fasta";

/// The 25 CC0 Pfam-A phage-hallmark families concatenated into
/// `phage_hallmark.hmm` (prophage detector DB; see `db/NEW_DETECTOR_DBS.md`).
/// Kept in ascending accession order so the concatenation reproduces the staged
/// `db/phage/phage_hallmark.hmm` byte-for-byte (the bundle was built that way).
const PHAGE_HALLMARK_PFAM: &[&str] = &[
    "PF00589", "PF00959", "PF03237", "PF03354", "PF03592", "PF03864", "PF03906",
    "PF04233", "PF04466", "PF04586", "PF04860", "PF04865", "PF04984", "PF04985",
    "PF05065", "PF05100", "PF05105", "PF05136", "PF05876", "PF06199", "PF07508",
    "PF10145", "PF10618", "PF11860", "PF16080",
];

/// Canonical tRNAscan-SE 2.0 distribution tarball (UCSC Lowe Lab). Its
/// `tRNAscan-SE-2.0/lib/models/` directory holds exactly the 99 covariance
/// models + cmpress indexes + isotype/mito signal files bactars' `--trna`
/// (trnascan-rs) consumes — verified byte-for-byte against `db/trna_models/`.
/// `setup-db` extracts ONLY that `models/` subdir, flattened into `trna_models/`.
const TRNA_MODELS_URL: &str = "https://trna.ucsc.edu/software/trnascan-se-2.0.12.tar.gz";

/// One downloadable artifact within an [`AuxSpec`].
enum AuxSource {
    /// Download `url` verbatim to `name` inside the DB's subdir.
    Plain { url: &'static str, name: &'static str },
    /// Concatenate gzipped Pfam-A HMMs (EBI InterPro API) into `name`. Each
    /// accession is fetched from [`pfam_hmm_url`], gunzipped, and appended in
    /// list order — reproducing the concatenated `.hmm` the detector scans.
    PfamConcat { accessions: &'static [&'static str], name: &'static str },
    /// Download a `.tar.gz` from `url`, extract it, and flatten every regular
    /// file that lives directly inside the archive's `models/` directory into
    /// the DB's subdir (reproducing the flat `db/trna_models/` layout `--trna`
    /// expects). `sentinel` names one representative model (e.g. `TRNA2-bact.cm`)
    /// used for the "already present" skip check and completion report.
    TarModels { url: &'static str, sentinel: &'static str },
}

/// A light-tier auxiliary reference DB laid out under `subdir`.
struct AuxSpec {
    /// Selector key (matches `SetupConfig::only`).
    key: &'static str,
    /// Subdirectory under `out_dir` (also the `--db` bundle convention path).
    subdir: &'static str,
    /// The bactars CLI flag (directory-valued) that points at the finished DB.
    flag: &'static str,
    /// Human-readable description for progress messages.
    desc: &'static str,
    /// Rough on-disk size in MB for the pre-download size warning.
    approx_disk_mb: u64,
    /// The artifacts that make up this DB, fetched in order.
    sources: &'static [AuxSource],
    /// Extra guidance printed once the DB is in place (e.g. derived files), or "".
    note: &'static str,
}

/// The light-tier auxiliary databases, in setup order. URLs verified live
/// 2026-07. Directory-valued flags — each DB resolves to its `<out>/<subdir>`.
const AUX_SPECS: &[AuxSpec] = &[
    AuxSpec {
        key: "amr_variant",
        // Shares the amr/ subdir with the `amr` model library (AMR.LIB): the
        // AMRFinderPlus point-mutation detector (crate::amr_variant) reads
        // AMRProt.fa + AMRProt-mutation.tsv from the SAME directory `--amr-variant`
        // resolves to, which the `--db` bundle maps to db/amr/ alongside AMR.LIB.
        // Without these two files the point-mutation stage silently skips.
        subdir: "amr",
        flag: "--amr-variant",
        desc: "AMRFinderPlus point-mutation reference (AMRProt.fa + mutation catalog)",
        approx_disk_mb: 5,
        // Both fetched verbatim (no decompression) from the same AMRFinderPlus
        // `database/latest/` base as AMR.LIB. amr_variant.rs opens exactly these two
        // files; AMRProt-susceptible.tsv / AMR_CDS.fa are not consumed → not fetched.
        // AMRProt.fa is listed first so aux_targets()[0] is the reference FASTA.
        sources: &[
            AuxSource::Plain {
                url: "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/\
                      AMRFinderPlus/database/latest/AMRProt.fa",
                name: "AMRProt.fa",
            },
            AuxSource::Plain {
                url: "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/\
                      AMRFinderPlus/database/latest/AMRProt-mutation.tsv",
                name: "AMRProt-mutation.tsv",
            },
        ],
        note: "amr_variant: AMRProt.fa + AMRProt-mutation.tsv land in the same amr/ \
               dir as AMR.LIB; pass that dir to --amr-variant for point-mutation calls.",
    },
    AuxSpec {
        key: "vfdb",
        subdir: "vfdb",
        flag: "--vfdb",
        desc: "VFDB setA core virulence factors (MGC/CAMS)",
        approx_disk_mb: 2,
        sources: &[AuxSource::Plain {
            url: "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz",
            name: "VFDB_setA_pro.fas.gz",
        }],
        note: "",
    },
    AuxSpec {
        key: "plasmidfinder",
        subdir: "plasmidfinder",
        flag: "--plasmid",
        desc: "PlasmidFinder replicon typing DB (DTU-CGE)",
        approx_disk_mb: 1,
        sources: &[AuxSource::Plain {
            url: "https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/master.tar.gz",
            name: "plasmidfinder_db.tar.gz",
        }],
        note: "",
    },
    AuxSpec {
        key: "gtdb_ssu",
        subdir: "gtdb_ssu",
        flag: "--species",
        desc: "GTDB 16S SSU species references (bac120 + ar53)",
        approx_disk_mb: 32,
        sources: &[
            AuxSource::Plain {
                url: "https://data.gtdb.ecogenomic.org/releases/latest/\
                      genomic_files_reps/bac120_ssu_reps.fna.gz",
                name: "bac120_ssu_reps.fna.gz",
            },
            AuxSource::Plain {
                url: "https://data.gtdb.ecogenomic.org/releases/latest/\
                      genomic_files_reps/ar53_ssu_reps.fna.gz",
                name: "ar53_ssu_reps.fna.gz",
            },
        ],
        note: "",
    },
    AuxSpec {
        key: "integron",
        subdir: "integron",
        flag: "--integron",
        desc: "IntegronFinder intI HMMs + attC CM (Institut Pasteur)",
        approx_disk_mb: 1,
        // All four models the integron stage reads (intI discriminators + the
        // phage-integrase exclusion HMM + the attC covariance model).
        sources: &[
            AuxSource::Plain {
                url: "https://raw.githubusercontent.com/gem-pasteur/Integron_Finder/\
                      master/integron_finder/data/Models/IntI.hmm",
                name: "IntI.hmm",
            },
            AuxSource::Plain {
                url: "https://raw.githubusercontent.com/gem-pasteur/Integron_Finder/\
                      master/integron_finder/data/Models/integron_integrase.hmm",
                name: "integron_integrase.hmm",
            },
            AuxSource::Plain {
                url: "https://raw.githubusercontent.com/gem-pasteur/Integron_Finder/\
                      master/integron_finder/data/Models/phage-int.hmm",
                name: "phage-int.hmm",
            },
            AuxSource::Plain {
                url: "https://raw.githubusercontent.com/gem-pasteur/Integron_Finder/\
                      master/integron_finder/data/Models/attc_4.cm",
                name: "attc_4.cm",
            },
        ],
        note: "",
    },
    AuxSpec {
        key: "phage",
        subdir: "phage",
        flag: "--prophage",
        desc: "prophage hallmark HMMs (25 CC0 Pfam families)",
        approx_disk_mb: 3,
        sources: &[AuxSource::PfamConcat {
            accessions: PHAGE_HALLMARK_PFAM,
            name: "phage_hallmark.hmm",
        }],
        note: "",
    },
    AuxSpec {
        key: "orit_pfam",
        subdir: "orit_pfam",
        flag: "--orit",
        desc: "MOB relaxase HMMs (6 CC0 Pfam families) + oriTDB oriT nt DB",
        approx_disk_mb: 1,
        // Two co-located assets, both landing in <out>/orit_pfam/ (the dir --orit
        // resolves to): the CC0 relaxase HMM library (relaxase-presence proxy) AND
        // oriTDB's curated oriT_exp.fasta, which crate::orit needs to run the real
        // oriT nucleotide-locus search. relaxase_pfam.hmm is listed first so
        // aux_targets()[0] stays the HMM library.
        sources: &[
            AuxSource::PfamConcat {
                accessions: ORIT_RELAXASE_PFAM,
                name: "relaxase_pfam.hmm",
            },
            AuxSource::Plain {
                url: ORIT_NT_EXP_URL,
                name: "oriT_exp.fasta",
            },
        ],
        note: "orit_pfam: oriT_exp.fasta is oriTDB's 122 experimentally-validated \
               oriT sequences (free for academic use; cite Li et al., NAR 2025 \
               53:D163). The nt-locus search additionally uses oriT_all.fasta as a \
               fallback if you drop it into this dir; only oriT_exp.fasta is fetched.",
    },
    AuxSpec {
        key: "meta",
        subdir: "meta",
        flag: "--meta",
        desc: "functional xref tables (PGAP names, pfam2go, KEGG brite)",
        approx_disk_mb: 23,
        sources: &[
            AuxSource::Plain {
                url: "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv",
                name: "hmm_PGAP.tsv",
            },
            AuxSource::Plain {
                url: "http://current.geneontology.org/ontology/external2go/pfam2go",
                name: "pfam2go",
            },
            AuxSource::Plain {
                url: "https://rest.kegg.jp/get/br:ko00001",
                name: "ko00001.keg",
            },
        ],
        // rfam_type.tsv, ko_cog.tsv and meta.sqlite are derived, not fetched.
        note: "meta: rfam_type.tsv, ko_cog.tsv and meta.sqlite are DERIVED artifacts \
               built from these files + Rfam/KEGG by the db/ builder scripts \
               (build_ko_cog.py, build_meta_sqlite.py); regenerate them there — \
               they are not downloaded by setup-db.",
    },
    AuxSpec {
        key: "trna",
        // The `--db` bundle resolver maps this to `trna_models/` (see
        // bin/bactars.rs `set_dir(&mut trna_models_dir, "trna_models")` and the
        // staged `db/trna_models/`), so the subdir is `trna_models`, NOT `trna`,
        // even though the selector key and flag are `trna` / `--trna`.
        subdir: "trna_models",
        flag: "--trna",
        desc: "tRNAscan-SE 2.0 covariance models (UCSC Lowe Lab)",
        approx_disk_mb: 18,
        // A single tarball whose `lib/models/` dir is flattened into trna_models/.
        // TRNA2-bact.cm is the representative sentinel for the skip check.
        sources: &[AuxSource::TarModels {
            url: TRNA_MODELS_URL,
            sentinel: "TRNA2-bact.cm",
        }],
        note: "trna: the 99 tRNAscan-SE 2.0 model files (covariance models + cmpress \
               indexes + isotype/mito signal files) are extracted flat into \
               trna_models/; pass that dir to --trna (trnascan-rs). Source: UCSC \
               Lowe Lab tRNAscan-SE 2.0 distribution (free for academic use).",
    },
];

/// All valid selector keys. The light-tier auxiliary keys sit between the model
/// libraries and the special-cased `isdb`.
const ALL_KEYS: &[&str] = &[
    "pfam",
    "ncbifams",
    "antifam",
    "amr",
    "amr_variant",
    "rfam",
    "vfdb",
    "plasmidfinder",
    "gtdb_ssu",
    "integron",
    "phage",
    "orit_pfam",
    "meta",
    "trna",
    "isdb",
];

/// Download and arrange the reference databases named by `cfg`.
///
/// On success the finished files sit under `cfg.out_dir` and the matching
/// bactars flags are printed to stdout. Databases whose target already exists
/// are skipped unless `cfg.force` is set. Returns an error only for a genuine
/// failure (bad selector, missing downloader, failed transfer/decompress); the
/// non-downloadable `isdb` selector is reported as guidance, not a failure.
pub fn setup_db(cfg: &SetupConfig) -> Result<(), String> {
    let selected = validate_selection(&cfg.only)?;

    // Warn about the on-disk footprint and confirm before downloading. Only the
    // downloadable specs contribute; `isdb` has no auto-download and no size.
    // Whether the IS DB will actually be downloaded (requested *and* a URL
    // resolvable — from --isdb-url / BACTARS_ISDB_URL, or the compiled-in default).
    let isdb_url = effective_isdb_url(cfg.isdb_url.as_deref());
    let will_fetch_isdb = selected.contains(&"isdb") && isdb_url.is_some();
    let mut total_mb: u64 = SPECS
        .iter()
        .filter(|s| selected.contains(&s.key))
        .map(|s| s.approx_disk_mb)
        .sum();
    let mut n_dbs = SPECS.iter().filter(|s| selected.contains(&s.key)).count();
    total_mb += AUX_SPECS
        .iter()
        .filter(|s| selected.contains(&s.key))
        .map(|s| s.approx_disk_mb)
        .sum::<u64>();
    n_dbs += AUX_SPECS.iter().filter(|s| selected.contains(&s.key)).count();
    if will_fetch_isdb {
        total_mb += 210; // rust-ise union DB, extracted (~204 MB).
        n_dbs += 1;
    }
    if cfg.build_psc {
        total_mb += 230_000; // UniRef90 gz (~32 GB) + filtered FASTA + ~194 GB mmseqs index.
        n_dbs += 1;
    }
    if cfg.build_kofam {
        total_mb += 6_000; // ko_list.gz + profiles.tar.gz (~2.4 GB) + extracted profiles/ + kofam_prok.hmm.
        n_dbs += 1;
    }
    if n_dbs > 0 {
        eprintln!(
            "About to download {n_dbs} database{} (~{:.1} GB on disk).",
            if n_dbs == 1 { "" } else { "s" },
            total_mb as f64 / 1024.0
        );
        if !confirm(cfg.yes)? {
            eprintln!("[abort] setup-db cancelled by user; nothing downloaded.");
            return Ok(());
        }
    }

    std::fs::create_dir_all(&cfg.out_dir)
        .map_err(|e| format!("cannot create {}: {e}", cfg.out_dir.display()))?;

    // (flag, absolute path) pairs to print once everything is in place.
    let mut ready: Vec<(String, PathBuf)> = Vec::new();
    // Per-DB completion rows: (key, size in bytes, model count).
    let mut summary: Vec<(&'static str, u64, u64)> = Vec::new();

    for spec in SPECS {
        if !selected.contains(&spec.key) {
            continue;
        }
        let dir = cfg.out_dir.join(spec.subdir);
        std::fs::create_dir_all(&dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
        let target = dir.join(spec.filename);

        if target.exists() && !cfg.force {
            eprintln!(
                "[skip] {} already present: {} (use --force to re-download)",
                spec.key,
                target.display()
            );
            record_completion(spec, &target, &mut summary);
            ready.push((spec.flag.to_string(), target));
            continue;
        }

        eprintln!("[fetch] {}: {}", spec.key, spec.desc);
        fetch_db(spec, &dir, &target)?;
        // Record the checksum for the freshly downloaded file.
        checksum_report(spec, &target);
        record_completion(spec, &target, &mut summary);
        ready.push((spec.flag.to_string(), target));
    }

    // Light-tier auxiliary DBs (FASTA / archive / table / concatenated HMM).
    for spec in AUX_SPECS {
        if !selected.contains(&spec.key) {
            continue;
        }
        let dir = cfg.out_dir.join(spec.subdir);
        std::fs::create_dir_all(&dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
        let targets = aux_targets(spec, &dir);

        if targets.iter().all(|t| t.exists()) && !cfg.force {
            eprintln!(
                "[skip] {} already present: {} (use --force to re-download)",
                spec.key,
                dir.display()
            );
            record_completion_aux(spec, &targets, &mut summary);
            if !spec.flag.is_empty() {
                ready.push((spec.flag.to_string(), dir));
            }
            continue;
        }

        eprintln!("[fetch] {}: {}", spec.key, spec.desc);
        fetch_aux(spec, &dir)?;
        for t in &targets {
            if let Ok(sha) = sha256_of(t) {
                let name = t.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                eprintln!("[checksum] {}: {sha}  {name}", spec.key);
            }
        }
        record_completion_aux(spec, &targets, &mut summary);
        if !spec.note.is_empty() {
            eprintln!("[note] {}", spec.note);
        }
        if !spec.flag.is_empty() {
            ready.push((spec.flag.to_string(), dir));
        }
    }

    // The IS DB has no fixed public URL: fetch+extract when one was supplied
    // (--isdb-url / BACTARS_ISDB_URL), otherwise fall back to guidance.
    if selected.contains(&"isdb") {
        match &isdb_url {
            Some(url) => {
                let dir = cfg.out_dir.join("isdb");
                std::fs::create_dir_all(&dir)
                    .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;
                eprintln!("[fetch] isdb: rust-ise IS profile DB (fpc-carrying)");
                let is_dir = fetch_isdb(url, &dir, cfg.force)?;
                // Confirm the extracted tree keeps `fpc/` next to mmdb_union — the
                // whole point of this tarball. Warn (don't fail) if it is missing,
                // because strict FP-control silently no-ops without it.
                if crate::is_elements::fpc_present(&is_dir) {
                    eprintln!(
                        "[isdb] IS DB ready at {} (fpc/refset present — strict FP-control active)",
                        is_dir.display()
                    );
                } else {
                    eprintln!(
                        "[isdb][warn] IS DB ready at {} but `fpc/refset.dbtype` is MISSING — strict\n\
                         [isdb][warn] FP-control will no-op (host-recombinase FPs NOT filtered). This\n\
                         [isdb][warn] tarball is not the fpc-carrying build; fetch rust-ise-isdb-fpc.tar.gz.",
                        is_dir.display()
                    );
                }
                ready.push(("--is-db".to_string(), is_dir));
            }
            None => eprintln!("{}", isdb_note(&cfg.out_dir)),
        }
    }

    // FULL-tier UniRef90 PSC naming DB (pure-Rust port of build_psc.sh; the only
    // external binary is mmseqs). Off unless requested via --full / --only psc.
    if cfg.build_psc {
        eprintln!(
            "[build] psc: UniRef90 PSC naming DB (FULL tier — heavy: ~32 GB download + ~194 GB index)"
        );
        let psc_dir = cfg.out_dir.join("psc");
        let uniref_gz = cfg.out_dir.join("uniref").join("uniref90.fasta.gz");
        crate::build_psc::build_psc(&crate::build_psc::PscBuildConfig {
            psc_dir: psc_dir.clone(),
            uniref_gz,
            threads: cfg.threads,
            force: cfg.force,
        })?;
        ready.push(("--psc".to_string(), psc_dir));
    }

    // FULL-tier KEGG KOfam naming DB (pure-Rust; no external binary). Off unless
    // requested via --full / --only kofam. Downloads ko_list.gz + profiles.tar.gz
    // from KEGG and concatenates the prokaryotic subset into kofam_prok.hmm.
    if cfg.build_kofam {
        eprintln!(
            "[build] kofam: KEGG KOfam naming DB (FULL tier — ~2.4 GB download, ~5-6 GB extracted)"
        );
        let kofam_dir = cfg.out_dir.join("kofam");
        crate::build_kofam::build_kofam(&crate::build_kofam::KofamBuildConfig {
            kofam_dir: kofam_dir.clone(),
            force: cfg.force,
        })?;
        ready.push(("--kofam".to_string(), kofam_dir));
    }

    print_flags(&ready);
    print_summary(&summary);
    Ok(())
}

/// Ask the user to confirm the download. Returns `Ok(true)` to proceed.
///
/// An explicit `yes` (from `--yes`) always proceeds. On an interactive terminal
/// with `yes` unset, we prompt and accept `y` / `yes` (case-insensitive).
/// Non-interactive stdin (piped / CI) with `yes` unset is NOT auto-confirmed:
/// silently downloading multiple GB in a script is unsafe, so we refuse with a
/// clear message and require an explicit `--yes` to proceed unattended.
fn confirm(yes: bool) -> Result<bool, String> {
    if yes {
        return Ok(true);
    }
    if !std::io::stdin().is_terminal() {
        eprintln!(
            "[abort] non-interactive (no TTY) and --yes not given: refusing to download \
             multi-GB databases unattended. Re-run with --yes to confirm."
        );
        return Ok(false);
    }
    eprint!("Continue? [y/N] ");
    std::io::stderr()
        .flush()
        .map_err(|e| format!("cannot flush prompt: {e}"))?;
    let mut line = String::new();
    std::io::stdin()
        .read_line(&mut line)
        .map_err(|e| format!("cannot read confirmation: {e}"))?;
    let ans = line.trim().to_ascii_lowercase();
    Ok(ans == "y" || ans == "yes")
}

/// Compute the final size + model count for `target` and push a summary row.
/// A failed count is logged and recorded as zero rather than aborting the run.
fn record_completion(spec: &DbSpec, target: &Path, summary: &mut Vec<(&'static str, u64, u64)>) {
    let size = std::fs::metadata(target).map(|m| m.len()).unwrap_or(0);
    let models = match count_models(spec, target) {
        Ok(n) => n,
        Err(e) => {
            eprintln!("[warn] {}: could not count models: {e}", spec.key);
            0
        }
    };
    eprintln!(
        "[done] {}: {} ({}, {} model{})",
        spec.key,
        target.display(),
        human_size(size),
        models,
        if models == 1 { "" } else { "s" }
    );
    summary.push((spec.key, size, models));
}

/// Count models in `target` by streaming lines and matching a per-format marker
/// header: `HMMER3/` for `.hmm` / AMR profile libraries, `INFERNAL1` for the
/// `.cm` covariance models. Reads via `BufReader` so multi-GB files never load
/// into memory at once.
fn count_models(spec: &DbSpec, target: &Path) -> Result<u64, String> {
    let marker = if spec.filename.ends_with(".cm") {
        "INFERNAL1"
    } else {
        "HMMER3/"
    };
    let file = std::fs::File::open(target)
        .map_err(|e| format!("cannot open {}: {e}", target.display()))?;
    let reader = BufReader::new(file);
    let mut n: u64 = 0;
    for line in reader.lines() {
        let line = line.map_err(|e| format!("read error on {}: {e}", target.display()))?;
        if line.starts_with(marker) {
            n += 1;
        }
    }
    Ok(n)
}

/// Compute and log the file's SHA-256 (pure-Rust, streamed). Best-effort:
/// a hashing error is logged and swallowed, never failing the run.
///
/// (The former best-effort sibling-`.md5` verification was dropped along with the
/// `curl` / `md5sum` shell-outs: no source in [`SPECS`] publishes a `.md5`, md5
/// cannot move onto the pure-Rust sha2 path, and adding an md5 crate is out of
/// scope — sha256 is recorded instead.)
fn checksum_report(spec: &DbSpec, target: &Path) {
    match sha256_of(target) {
        Ok(sha) => eprintln!("[checksum] {}: {sha}  {}", spec.key, spec.filename),
        Err(e) => eprintln!("[checksum] {}: sha256 unavailable ({e})", spec.key),
    }
}

/// Compute the file's SHA-256 hex digest in-process (pure-Rust sha2).
fn sha256_of(target: &Path) -> Result<String, String> {
    crate::util_io::sha256_file(target)
}

/// Format a byte count as a human-readable MB/GB string (binary units).
fn human_size(bytes: u64) -> String {
    const MB: f64 = 1024.0 * 1024.0;
    const GB: f64 = MB * 1024.0;
    let b = bytes as f64;
    if b >= GB {
        format!("{:.2} GB", b / GB)
    } else {
        format!("{:.1} MB", b / MB)
    }
}

/// Print the final per-DB completion table (key, human size, model count).
fn print_summary(summary: &[(&'static str, u64, u64)]) {
    if summary.is_empty() {
        return;
    }
    // Align columns to the widest entry for a clean table.
    let key_w = summary.iter().map(|(k, _, _)| k.len()).max().unwrap_or(3).max(3);
    let size_strs: Vec<String> = summary.iter().map(|(_, s, _)| human_size(*s)).collect();
    let size_w = size_strs.iter().map(|s| s.len()).max().unwrap_or(4).max(4);
    println!("\n# database summary:");
    println!("{:<key_w$}  {:>size_w$}  MODELS", "KEY", "SIZE");
    for ((key, _size, models), size_str) in summary.iter().zip(size_strs.iter()) {
        println!("{key:<key_w$}  {size_str:>size_w$}  {models}");
    }
}

/// Resolve `only` into the concrete set of keys to process, rejecting unknowns.
/// The downloadable DB keys a mode selects. `--fast` and `--full` both fetch the
/// full standard set ([`ALL_KEYS`]); `--full` additionally BUILDS the PSC DB (via
/// `SetupConfig::build_psc`) and points the user at kofam, neither of which is a
/// plain download. Returned as owned strings for the CLI's `only` list.
pub fn mode_keys() -> Vec<String> {
    ALL_KEYS.iter().map(|k| k.to_string()).collect()
}

fn validate_selection(only: &[String]) -> Result<Vec<&'static str>, String> {
    if only.is_empty() {
        return Ok(ALL_KEYS.to_vec());
    }
    let mut out: Vec<&'static str> = Vec::new();
    for want in only {
        match ALL_KEYS.iter().find(|k| **k == want.as_str()) {
            Some(k) => {
                if !out.contains(k) {
                    out.push(*k);
                }
            }
            None => {
                return Err(format!(
                    "unknown database '{want}'; valid: {}",
                    ALL_KEYS.join(", ")
                ))
            }
        }
    }
    Ok(out)
}

/// Download `spec` and produce `target`, handling decompression per `spec.fetch`.
fn fetch_db(spec: &DbSpec, dir: &Path, target: &Path) -> Result<(), String> {
    match &spec.fetch {
        Fetch::Plain => {
            download(spec.url, target)?;
        }
        Fetch::Gunzip => {
            // Download to "<target>.gz", then gunzip it to <target>.
            let gz = with_extra_ext(target, "gz");
            download(spec.url, &gz)?;
            gunzip(&gz)?;
            if !target.exists() {
                return Err(format!(
                    "gunzip did not produce expected file {}",
                    target.display()
                ));
            }
        }
        Fetch::TarMember(member) => {
            let tgz = dir.join("download.tar.gz");
            download(spec.url, &tgz)?;
            tar_extract_member(&tgz, member, dir)?;
            let extracted = dir.join(member);
            if extracted != *target {
                std::fs::rename(&extracted, target).map_err(|e| {
                    format!(
                        "cannot move {} to {}: {e}",
                        extracted.display(),
                        target.display()
                    )
                })?;
            }
            let _ = std::fs::remove_file(&tgz);
            if !target.exists() {
                return Err(format!(
                    "tar did not extract expected member '{member}' to {}",
                    target.display()
                ));
            }
        }
    }
    Ok(())
}

/// The final on-disk paths an [`AuxSpec`] produces under `dir`, one per source.
fn aux_targets(spec: &AuxSpec, dir: &Path) -> Vec<PathBuf> {
    spec.sources
        .iter()
        .map(|s| match s {
            AuxSource::Plain { name, .. } => dir.join(name),
            AuxSource::PfamConcat { name, .. } => dir.join(name),
            // The flattened model set is a directory of many files; its sentinel
            // (a representative `.cm`) stands in for the "already present" check.
            AuxSource::TarModels { sentinel, .. } => dir.join(sentinel),
        })
        .collect()
}

/// Download every source of an auxiliary DB into `dir`.
fn fetch_aux(spec: &AuxSpec, dir: &Path) -> Result<(), String> {
    for src in spec.sources {
        match src {
            AuxSource::Plain { url, name } => {
                let target = dir.join(name);
                download(url, &target)?;
                if !target.exists() {
                    return Err(format!("download did not produce {}", target.display()));
                }
            }
            AuxSource::PfamConcat { accessions, name } => {
                let target = dir.join(name);
                fetch_pfam_concat(accessions, dir, &target)?;
            }
            AuxSource::TarModels { url, sentinel } => {
                fetch_tar_models(url, dir, sentinel)?;
            }
        }
    }
    Ok(())
}

/// Download a tRNAscan-SE 2.0 distribution tarball and flatten its
/// `*/lib/models/*` files directly into `dir` — the covariance-model set that
/// `--trna` (trnascan-rs) consumes. Only mmseqs is an allowed external binary,
/// so extraction is pure-Rust (util_io tar/gzip); the whole tarball is unpacked
/// to a temp dir, the `lib/models` files copied up, and the temp cleaned.
fn fetch_tar_models(url: &str, dir: &Path, sentinel: &str) -> Result<(), String> {
    let tgz = dir.join(".trna_models.tar.gz");
    download(url, &tgz)?;
    let tmp = dir.join(".trna_extract");
    let _ = std::fs::remove_dir_all(&tmp);
    std::fs::create_dir_all(&tmp).map_err(|e| format!("cannot create {}: {e}", tmp.display()))?;
    crate::util_io::extract_tar_gz(&tgz, &tmp)?;
    let models = find_models_dir(&tmp)
        .ok_or_else(|| format!("no `*/lib/models` directory inside {url}"))?;
    let mut n = 0usize;
    for entry in std::fs::read_dir(&models).map_err(|e| format!("read {}: {e}", models.display()))? {
        let entry = entry.map_err(|e| e.to_string())?;
        let from = entry.path();
        if from.is_file() {
            let to = dir.join(entry.file_name());
            std::fs::copy(&from, &to)
                .map_err(|e| format!("copy model {}: {e}", from.display()))?;
            n += 1;
        }
    }
    let _ = std::fs::remove_dir_all(&tmp);
    let _ = std::fs::remove_file(&tgz);
    if !dir.join(sentinel).exists() {
        return Err(format!(
            "tar-models extract produced no sentinel {sentinel} ({n} files copied from {})",
            models.display()
        ));
    }
    Ok(())
}

/// Find a `<root>/*/lib/models` directory (the tRNAscan-SE layout), if present.
fn find_models_dir(root: &Path) -> Option<PathBuf> {
    for e in std::fs::read_dir(root).ok()?.flatten() {
        let cand = e.path().join("lib").join("models");
        if cand.is_dir() {
            return Some(cand);
        }
    }
    None
}

/// Build a concatenated `.hmm` from a list of Pfam-A accessions: fetch each
/// model's gzipped HMM from the EBI InterPro API, gunzip it, and append it to
/// `target` in list order. The per-accession temp files are cleaned up.
fn fetch_pfam_concat(
    accessions: &[&str],
    dir: &Path,
    target: &Path,
) -> Result<(), String> {
    use std::io::Write as _;
    let mut out = std::fs::File::create(target)
        .map_err(|e| format!("cannot create {}: {e}", target.display()))?;
    for acc in accessions {
        let url = pfam_hmm_url(acc);
        let gz = dir.join(format!("{acc}.hmm.gz"));
        download(&url, &gz)?;
        gunzip(&gz)?; // -> dir/<acc>.hmm
        let hmm = dir.join(format!("{acc}.hmm"));
        let data = std::fs::read(&hmm)
            .map_err(|e| format!("cannot read fetched {}: {e}", hmm.display()))?;
        out.write_all(&data)
            .map_err(|e| format!("cannot append {acc} to {}: {e}", target.display()))?;
        let _ = std::fs::remove_file(&hmm);
    }
    out.flush()
        .map_err(|e| format!("cannot flush {}: {e}", target.display()))?;
    Ok(())
}

/// Record a completion row for an auxiliary DB: total size across its files and
/// a "units" count (HMM/CM profiles for model files, FASTA records for uncompressed
/// FASTA, 0 for gzipped/tabular files it cannot cheaply count).
fn record_completion_aux(spec: &AuxSpec, targets: &[PathBuf], summary: &mut Vec<(&'static str, u64, u64)>) {
    let mut size = 0u64;
    let mut units = 0u64;
    for t in targets {
        size += std::fs::metadata(t).map(|m| m.len()).unwrap_or(0);
        units += count_aux_units(t);
    }
    eprintln!(
        "[done] {}: {} ({}, {} file{})",
        spec.key,
        spec.subdir,
        human_size(size),
        targets.len(),
        if targets.len() == 1 { "" } else { "s" }
    );
    summary.push((spec.key, size, units));
}

/// Best-effort count of "units" in an auxiliary file, chosen by extension:
/// HMMER3 profiles (`.hmm`), INFERNAL models (`.cm`), or FASTA records
/// (`.fas`/`.fasta`/`.fna`). Gzipped and tabular files return 0 (not counted).
fn count_aux_units(target: &Path) -> u64 {
    let name = target.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let marker = if name.ends_with(".hmm") {
        "HMMER3/"
    } else if name.ends_with(".cm") {
        "INFERNAL1"
    } else if name.ends_with(".fas") || name.ends_with(".fasta") || name.ends_with(".fna") {
        ">"
    } else {
        return 0; // .gz / .tar.gz / .tsv / .keg / pfam2go — not counted.
    };
    let file = match std::fs::File::open(target) {
        Ok(f) => f,
        Err(_) => return 0,
    };
    let mut n = 0u64;
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if line.starts_with(marker) {
            n += 1;
        }
    }
    n
}

/// Fetch `url` into `out` in-process over HTTP(S) (pure-Rust ureq/rustls,
/// following redirects). Fails on any transfer error, leaving no partial file.
fn download(url: &str, out: &Path) -> Result<(), String> {
    if let Err(e) = crate::util_io::http_download(url, out) {
        // Leave no truncated file behind to be mistaken for a good download.
        let _ = std::fs::remove_file(out);
        return Err(format!("download failed for {url}: {e}"));
    }
    Ok(())
}

/// Decompress `path` (a `.gz`) to its extension-stripped sibling and remove the
/// original `.gz`, mirroring `gunzip -f`. Pure-Rust (flate2, via `util_io`).
fn gunzip(path: &Path) -> Result<(), String> {
    // `foo.hmm.gz` -> `foo.hmm`, `Rfam.cm.gz` -> `Rfam.cm`, `x.fna.gz` -> `x.fna`.
    let dst = path.with_extension("");
    if dst == path {
        return Err(format!("not a .gz path: {}", path.display()));
    }
    crate::util_io::gunzip_file(path, &dst)?;
    let _ = std::fs::remove_file(path);
    Ok(())
}

/// Extract tar.gz `archive` into directory `dir` so `dir/<member>` lands in place.
/// Uses the pure-Rust, traversal-guarded, no-same-owner extractor; it unpacks the
/// whole (small) archive and the caller then picks out `member`.
fn tar_extract_member(archive: &Path, _member: &str, dir: &Path) -> Result<(), String> {
    crate::util_io::extract_tar_gz(archive, dir)?;
    Ok(())
}

/// Download and extract the rust-ise IS profile-DB tarball, returning the
/// directory to pass to `--is-db` (the one holding `manifest_union.tsv`).
///
/// The tarball unpacks to a versioned directory (`isscan_db_vX/`) containing the
/// MMseqs2 profile DB (`mmdb_union/…`) and `manifest_union.tsv`. We extract the
/// whole archive under `<out>/isdb/` and locate that directory so the returned
/// path is correct regardless of the exact version string inside the tarball.
fn fetch_isdb(url: &str, dir: &Path, force: bool) -> Result<PathBuf, String> {
    // If an already-extracted DB is present, reuse it unless --force.
    if !force {
        if let Some(found) = find_manifest_dir(dir) {
            eprintln!("[skip] isdb already present: {} (use --force to re-download)", found.display());
            return Ok(found);
        }
    }
    let tgz = dir.join("isscan_db.tar.gz");
    download(url, &tgz)?;
    // Integrity-check the tarball BEFORE extraction, but only when the effective
    // URL is the compiled-in default (whose sha256 we ship). A user `--isdb-url`
    // override supplies its own data, so we cannot pin its hash → skip the check.
    if url == DEFAULT_ISDB_URL && !DEFAULT_ISDB_SHA256.is_empty() {
        if let Err(e) = crate::util_io::verify_sha256(&tgz, DEFAULT_ISDB_SHA256) {
            let _ = std::fs::remove_file(&tgz);
            return Err(format!(
                "IS DB tarball integrity check failed (refusing to extract): {e}"
            ));
        }
        eprintln!("[isdb] sha256 verified against compiled-in DEFAULT_ISDB_SHA256");
    }
    crate::util_io::extract_tar_gz(&tgz, dir)?;
    let _ = std::fs::remove_file(&tgz);
    find_manifest_dir(dir).ok_or_else(|| {
        format!(
            "IS DB tarball extracted to {} but no directory containing \
             `manifest_union.tsv` was found (is this the isscan_db artifact?)",
            dir.display()
        )
    })
}

/// Find the directory under `root` (checking `root` itself, then its immediate
/// children, then one level deeper) that contains `manifest_union.tsv`.
fn find_manifest_dir(root: &Path) -> Option<PathBuf> {
    const MANIFEST: &str = "manifest_union.tsv";
    if root.join(MANIFEST).is_file() {
        return Some(root.to_path_buf());
    }
    // Immediate children, then grandchildren (the tarball nests one level).
    for depth1 in read_subdirs(root) {
        if depth1.join(MANIFEST).is_file() {
            return Some(depth1);
        }
        for depth2 in read_subdirs(&depth1) {
            if depth2.join(MANIFEST).is_file() {
                return Some(depth2);
            }
        }
    }
    None
}

/// Immediate subdirectories of `dir` (empty on any I/O error).
fn read_subdirs(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir() {
                out.push(p);
            }
        }
    }
    out
}

/// Append an extra extension to a path (`foo.hmm` + `gz` -> `foo.hmm.gz`).
fn with_extra_ext(path: &Path, ext: &str) -> PathBuf {
    let mut s = path.as_os_str().to_owned();
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

/// Guidance printed when `isdb` is requested but no source URL is resolvable
/// (neither `--isdb-url` / `BACTARS_ISDB_URL` nor the compiled-in default).
fn isdb_note(out_dir: &Path) -> String {
    let dir = out_dir.join("isdb");
    format!(
        "[isdb] The rust-ise IS profile DB cannot be auto-downloaded (no BACTARS_ISDB_URL\n\
         [isdb] set and no compiled-in default). It is a versioned artifact calibrated to a\n\
         [isdb] specific rust-ise release (ISOSDB ∪ ISfinder MMseqs2 profiles + matching\n\
         [isdb] thresholds) with no fixed public URL.\n\
         [isdb]\n\
         [isdb] IMPORTANT: use the *fpc-carrying* build (`rust-ise-isdb-fpc.tar.gz`). Strict\n\
         [isdb] precision mode (default) drops host-recombinase false positives ONLY when the\n\
         [isdb] DB ships an `fpc/refset` subdirectory; without it, strict FP-control silently\n\
         [isdb] no-ops. The fpc tarball extracts to `mmdb_union/`, `fpc/refset*` and\n\
         [isdb] `manifest_union.tsv` side by side.\n\
         [isdb]\n\
         [isdb] To provision it, host that tarball and either:\n\
         [isdb]   export BACTARS_ISDB_URL=https://<host>/rust-ise-isdb-fpc.tar.gz\n\
         [isdb]   bactars setup-db --only isdb --out <out>          # auto-fetch + extract\n\
         [isdb] or extract it by hand to {}, confirming it contains\n\
         [isdb] both `mmdb_union/profileDb*` and `fpc/refset.dbtype`.\n\
         [isdb] Then pass:  --is-db {}   (mmseqs must also be on PATH).",
        dir.display(),
        dir.display()
    )
}

/// Print the ready-to-use bactars flags for every database that is in place.
fn print_flags(ready: &[(String, PathBuf)]) {
    if ready.is_empty() {
        return;
    }
    println!("\n# databases ready — pass these flags to bactars:");
    // One flag per line for readability, plus a single copy-paste line.
    let mut joined = String::new();
    for (flag, path) in ready {
        println!("{flag} {}", path.display());
        if !joined.is_empty() {
            joined.push(' ');
        }
        joined.push_str(flag);
        joined.push(' ');
        joined.push_str(&path.display().to_string());
    }
    println!("\n# one line:\nbactars <genome.fna> {joined}");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pfam_hmm_url_is_interpro_api() {
        assert_eq!(
            pfam_hmm_url("PF01076"),
            "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF01076?annotation=hmm"
        );
        // Every phage/oriT accession yields a well-formed API URL.
        for acc in ORIT_RELAXASE_PFAM.iter().chain(PHAGE_HALLMARK_PFAM.iter()) {
            let u = pfam_hmm_url(acc);
            assert!(u.starts_with("https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/"));
            assert!(u.ends_with("?annotation=hmm"));
            assert!(u.contains(acc));
        }
    }

    #[test]
    fn pfam_accession_lists_match_provenance() {
        // db/orit_pfam/PROVENANCE.md: exactly these six CC0 relaxase families.
        assert_eq!(
            ORIT_RELAXASE_PFAM,
            &["PF01076", "PF03389", "PF03432", "PF05713", "PF07514", "PF13814"]
        );
        // db/NEW_DETECTOR_DBS.md: 25 CC0 phage-hallmark families.
        assert_eq!(PHAGE_HALLMARK_PFAM.len(), 25);
        // No duplicate accessions in either list, and each is a `PF#####` token.
        for list in [ORIT_RELAXASE_PFAM, PHAGE_HALLMARK_PFAM] {
            let mut seen = std::collections::HashSet::new();
            for acc in list {
                assert!(seen.insert(*acc), "duplicate accession {acc}");
                assert!(acc.starts_with("PF"), "bad accession {acc}");
                assert_eq!(acc.len(), 7, "bad accession length {acc}");
                assert!(acc[2..].chars().all(|c| c.is_ascii_digit()), "bad accession {acc}");
            }
        }
    }

    #[test]
    fn aux_keys_are_registered_and_unique() {
        // Every AUX_SPECS key is a valid selector, and there is no key collision
        // with the model-library SPECS or the special `isdb` selector.
        let spec_keys: Vec<&str> = SPECS.iter().map(|s| s.key).collect();
        let mut seen = std::collections::HashSet::new();
        for aux in AUX_SPECS {
            assert!(ALL_KEYS.contains(&aux.key), "{} missing from ALL_KEYS", aux.key);
            assert!(!spec_keys.contains(&aux.key), "{} collides with SPECS", aux.key);
            assert_ne!(aux.key, "isdb");
            assert!(seen.insert(aux.key), "duplicate aux key {}", aux.key);
            assert!(!aux.sources.is_empty(), "{} has no sources", aux.key);
        }
        // The light-tier keys we promised are all present.
        for want in ["vfdb", "plasmidfinder", "gtdb_ssu", "integron", "phage", "orit_pfam", "meta"] {
            assert!(AUX_SPECS.iter().any(|s| s.key == want), "missing aux DB {want}");
        }
        // kofam / psc are deliberately NOT auto-fetchable selectors.
        assert!(!ALL_KEYS.contains(&"kofam"));
        assert!(!ALL_KEYS.contains(&"psc"));
    }

    #[test]
    fn aux_flags_match_pipeline_bundle_convention() {
        // Directory-valued flags the `--db` bundle resolver maps by subdir.
        let expect = [
            ("vfdb", "vfdb", "--vfdb"),
            ("plasmidfinder", "plasmidfinder", "--plasmid"),
            ("gtdb_ssu", "gtdb_ssu", "--species"),
            ("integron", "integron", "--integron"),
            ("phage", "phage", "--prophage"),
            ("orit_pfam", "orit_pfam", "--orit"),
            ("meta", "meta", "--meta"),
        ];
        for (key, subdir, flag) in expect {
            let s = AUX_SPECS.iter().find(|s| s.key == key).expect("aux spec present");
            assert_eq!(s.subdir, subdir, "{key} subdir");
            assert_eq!(s.flag, flag, "{key} flag");
        }
    }

    #[test]
    fn aux_targets_names_match_on_disk_layout() {
        let dir = Path::new("/out/gtdb_ssu");
        let gtdb = AUX_SPECS.iter().find(|s| s.key == "gtdb_ssu").unwrap();
        let t = aux_targets(gtdb, dir);
        assert_eq!(t[0], dir.join("bac120_ssu_reps.fna.gz"));
        assert_eq!(t[1], dir.join("ar53_ssu_reps.fna.gz"));

        let vfdb = AUX_SPECS.iter().find(|s| s.key == "vfdb").unwrap();
        assert_eq!(
            aux_targets(vfdb, Path::new("/out/vfdb"))[0],
            Path::new("/out/vfdb/VFDB_setA_pro.fas.gz")
        );
        // orit_pfam ships two co-located assets in the --orit dir: the CC0
        // relaxase HMM library (first) and oriTDB's oriT nucleotide DB (second).
        let orit = AUX_SPECS.iter().find(|s| s.key == "orit_pfam").unwrap();
        let ot = aux_targets(orit, Path::new("/out/orit_pfam"));
        assert_eq!(ot[0], Path::new("/out/orit_pfam/relaxase_pfam.hmm"));
        assert_eq!(ot[1], Path::new("/out/orit_pfam/oriT_exp.fasta"));
    }

    #[test]
    fn amr_variant_provisions_point_mutation_files_into_amr_dir() {
        // crate::amr_variant reads AMRProt.fa + AMRProt-mutation.tsv from the dir
        // --amr-variant resolves to, which the --db bundle maps to db/amr/ (shared
        // with the AMR.LIB model library). setup-db must land BOTH files there or the
        // point-mutation stage silently skips. Guard the flag, the shared subdir, the
        // filenames, and that it rides along in the standard set.
        let amrv = AUX_SPECS
            .iter()
            .find(|s| s.key == "amr_variant")
            .expect("amr_variant aux spec present");
        assert_eq!(amrv.flag, "--amr-variant");
        assert_eq!(amrv.subdir, "amr"); // co-located with AMR.LIB
        assert!(
            ALL_KEYS.contains(&"amr_variant"),
            "amr_variant must be in the standard set (--fast/--full)"
        );
        // Both files land in the amr/ dir with their canonical names.
        let t = aux_targets(amrv, Path::new("/out/amr"));
        assert_eq!(t, vec![
            PathBuf::from("/out/amr/AMRProt.fa"),
            PathBuf::from("/out/amr/AMRProt-mutation.tsv"),
        ]);
        // Both are verbatim Plain fetches from the SAME AMRFinderPlus database base
        // as AMR.LIB (no decompression), so the URLs share that base path.
        for src in amrv.sources {
            match src {
                AuxSource::Plain { url, .. } => assert!(
                    url.contains(
                        "pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/"
                    ),
                    "unexpected amr_variant source base: {url}"
                ),
                _ => panic!("amr_variant sources must be Plain verbatim fetches"),
            }
        }
        // It is an AUX key, distinct from the `amr` model-library SPEC key.
        assert!(SPECS.iter().any(|s| s.key == "amr"));
        assert!(!SPECS.iter().any(|s| s.key == "amr_variant"));
    }

    #[test]
    fn orit_pfam_provisions_orit_nt_db_into_orit_dir() {
        // crate::orit looks for oriT_exp.fasta INSIDE the dir --orit resolves to.
        // setup-db must therefore land oriT_exp.fasta in <out>/orit_pfam/ (the
        // same subdir the --orit flag points at) or the real oriT nt-locus search
        // never runs. Guard both: the flag and the provisioned filename.
        let orit = AUX_SPECS.iter().find(|s| s.key == "orit_pfam").unwrap();
        assert_eq!(orit.flag, "--orit");
        assert_eq!(orit.subdir, "orit_pfam");
        // An oriT_exp.fasta source is present and it is a plain (verbatim) fetch
        // from the pinned oriTDB URL — not decompressed / concatenated.
        let has_orit_nt = orit.sources.iter().any(|s| matches!(
            s,
            AuxSource::Plain { url, name }
                if *name == "oriT_exp.fasta" && *url == ORIT_NT_EXP_URL
        ));
        assert!(has_orit_nt, "orit_pfam must fetch oriT_exp.fasta for the nt-locus search");
        // The relaxase HMM proxy DB is still provisioned alongside it.
        let has_relaxase = orit.sources.iter().any(|s| matches!(
            s,
            AuxSource::PfamConcat { name, .. } if *name == "relaxase_pfam.hmm"
        ));
        assert!(has_relaxase, "orit_pfam must still ship relaxase_pfam.hmm");
        // The pinned URL is the oriTDB2 experimentally-validated oriT FASTA.
        assert!(ORIT_NT_EXP_URL.ends_with("/oriT_exp.fasta"));
        assert!(ORIT_NT_EXP_URL.contains("oriTDB2"));
    }

    #[test]
    fn effective_isdb_url_prefers_override_then_default() {
        // An explicit non-empty override always wins.
        assert_eq!(
            effective_isdb_url(Some("https://h/db.tar.gz")).as_deref(),
            Some("https://h/db.tar.gz")
        );
        // Empty override is treated as unset -> falls through to the default.
        let empty = effective_isdb_url(Some(""));
        let none = effective_isdb_url(None);
        assert_eq!(empty, none, "empty override behaves like no override");
        // With the placeholder default empty, no override yields None (guidance path).
        if DEFAULT_ISDB_URL.is_empty() {
            assert!(none.is_none());
        } else {
            assert_eq!(none.as_deref(), Some(DEFAULT_ISDB_URL));
        }
    }

    #[test]
    fn validate_selection_accepts_new_keys_and_rejects_unknown() {
        let ok = validate_selection(&["vfdb".into(), "orit_pfam".into(), "meta".into()]).unwrap();
        assert_eq!(ok, vec!["vfdb", "orit_pfam", "meta"]);
        // De-duplicates repeated selectors.
        let dedup = validate_selection(&["phage".into(), "phage".into()]).unwrap();
        assert_eq!(dedup, vec!["phage"]);
        // Empty selection means "all".
        assert_eq!(validate_selection(&[]).unwrap(), ALL_KEYS.to_vec());
        // Unknown key errors and names the bad token.
        let err = validate_selection(&["kofam".into()]).unwrap_err();
        assert!(err.contains("kofam"));
    }
}
