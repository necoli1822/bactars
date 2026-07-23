# bactars

Rust-native, container-free prokaryotic genome annotation pipeline — the
orchestrator that wires the Rust subunits together as in-process libraries over a
shared `Feature` model.

**Everything is pure Rust except one external tool: `mmseqs2`.** Every annotation
engine — gene calling, HMM search, tRNA/tmRNA, ncRNA/rRNA (Infernal), CRISPR — is a
Rust reimplementation consumed in-process, each byte-for-byte faithful to its
canonical C tool. The single non-Rust dependency is `mmseqs2`, which the
protein-similarity stages (IS elements, UniRef90/PSC naming, VFDB, plasmid, 16S
species-ID, AMR-variant) shell out to. Everything else — including DB download,
gzip decompression, and tar extraction — runs in-process (`ureq`/`flate2`/`tar`/
`sha2`); no Docker, no conda, no Python, no `curl`/`wget`/`gunzip`/`tar` binaries,
no subprocess orchestration of our own beyond `mmseqs2`.

All subunits are consumed as **published crates.io libraries** (no path deps, no
subprocesses of our own) over a shared `Feature` model:

```
genome.fna
  ├─ rustygal 0.2.0    CDS prediction + translation        (in-process lib)
  │    └─ rustyhmmer 0.1.3   HMM annotation                 (in-process lib)
  │         ├─ Pfam / NCBIfams / AntiFam  (--cut_ga)
  │         ├─ AMR: NCBI AMRFinderPlus HMM (--cut_tc)       ← Bakta-style AMR
  │         └─ Feature { CDS + annotations }
  ├─ trnascan-rs 0.2.1  tRNA (bacterial -B, embeds infernox) (in-process lib)
  ├─ oxagorn 0.2.0      tmRNA (ARAGORN)                      (in-process lib)
  ├─ rust-ise 0.2.3     IS elements (ISEScan; needs mmseqs)  (in-process lib)
  ├─ infernox 0.2.0     ncRNA (Infernal cmsearch vs Rfam)    (in-process lib)
  └─ minced-rs          CRISPR arrays (MinCED port, no DB)   (in crate)
        ↓
  resolve (RNA masks CDS, dedup multi-model ncRNA)
        ↓
  output: GFF3 (noodles) / GBFF
```

`bactars setup-db --out DIR` provisions the reference databases and prints the
flags to point bactars at them. It is **pure-Rust — only mmseqs is ever spawned**:
downloads use an in-process rustls HTTP client, and even the FULL-tier UniRef90
PSC and KEGG KOfam DBs are **built in process** (no curl / pigz / gawk / python).
Modes mirror the annotation modes:

- `setup-db --fast` — every standard annotation DB (Pfam / NCBIfams / AntiFam /
  AMRFinderPlus HMM **+ point-mutation refs** / Rfam / **tRNA models** / VFDB /
  plasmidfinder / GTDB-SSU / integron / phage / oriT / meta).
- `setup-db --full` — `--fast` plus the FULL-tier naming DBs **built** in process:
  the UniRef90 PSC mmseqs DB (~32 GB download → ~194 GB index) and the KOfam HMM
  DB (`kofam_prok.hmm` from KEGG `prokaryote.hal`).
- `setup-db --only <keys>` — a subset; `--only psc` / `--only kofam` build just
  that FULL-tier DB, `--threads N` sizes the PSC mmseqs index build. The IS DB is
  release-locked (`--is-db-url` / `BACTARS_ISDB_URL`).

Every downloadable DB has a `db/` directory; `--db <bundle>` then resolves every
per-DB path from that one root by convention.

- **CDS**: `rustygal` (Rust prodigal). Byte-identical gene set to C prodigal
  (MG1655: 4319 CDS). The **genetic code is auto-detected**: bactars trains a gene
  model under NCBI table 11 (standard) and table 4 (Mollicutes, TGA=Trp) and keeps
  the one with the higher dynamic-programming path score — prodigal's own
  metagenomic bin-selection criterion, so a normal genome is never mis-flipped while
  a *Mycoplasma*/*Spiroplasma*/*Ureaplasma* genome (where table 11 fragments genes at
  every internal TGA) is correctly called under table 4. On the RefSeq benchmark this
  lifts Mollicute CDS precision from ~27 % (fragmented) to ~95–99 %. `--translation-table N`
  (alias `-g`) forces a specific code; `-g 25` selects Gracilibacteria/SR1 (TGA=Gly),
  which shares table 4's gene structure and so can only be requested explicitly.
- **HMM**: `rustyhmmer` (`api::HmmAnnotator`, byte-identical to C HMMER 3.4).
  Pfam / NCBIfams / AntiFam under the gathering cutoff, matching Bakta.
- **tRNA**: `trnascan-rs` (Rust tRNAscan-SE 2.0 + embedded infernox). Byte-identical
  to C tRNAscan-SE 2.0 + Infernal 1.1.5 (`-B`). Enabled with `--trna <data/models>`.
- **tmRNA**: `oxagorn` (`api::detect`, byte-identical to C ARAGORN v1.2.41).
  Pure-Rust, no external DB. Enabled with `--tmrna` (MG1655: 1 tmRNA = ssrA).
- **IS elements**: `rust-ise` (`run()`, byte-identical to the ISEScan-equivalent
  binary). The **only** subunit that uses an external tool — it shells out to
  `mmseqs2` for IS profile search. Needs an IS DB via `--is-db <dir>` (MG1655: 71 IS
  calls, all 49 curated ground-truth loci recovered). The IS DB is built **reproducibly
  from source in pure Rust** via `rust-ise build-db`: ISOSDB ∪ ISfinder ORFs (rustygal
  META gene-finding) → mmseqs clustering → in-process **poasta** POA multiple-sequence
  alignment (replaces clustalo) → mmseqs profiles. The build is **byte-identical across
  machines** (deterministic MSA + `msa2profile --threads 1`), and per-profile hit
  thresholds are **GTDB-calibrated** on 113k genomes (embedded in the crate), so a hit
  is only accepted above the profile's contaminant/promiscuous identity floor.
- **AMR**: done the **Bakta way** — an AMR HMM database (NCBI AMRFinderPlus, or the
  AMR subset of NCBIfams) searched on the same in-process `rustyhmmer` path with the
  trusted cutoff (`--cut_tc`). Enabled with `--amr <db.hmm>`; hits appear as
  `rustyhmmer:amr` annotations. No ARGenus / minimap2 / PAF orchestration.
- **ncRNA**: `infernox` (`infernal::faithful_search::FaithfulSearcher`, byte-parity
  Infernal cmsearch). Scans the genome against a (multi-model) Rfam `.cm` DB on both
  strands and reports significant hits (E ≤ 1e-3). Enabled with `--ncrna <cm_db>`.
  Models are read in **global** config (`cm_file_read_from_reader_opt(.., false)`), which
  is what `FaithfulSearcher::new` expects. Verified on MG1655: RNaseP (RF00010, E=3e-105),
  6S/ssrS (RF00013, E=4e-30), and small-SRP (RF00169, E=1e-14) at their known loci.
  **rRNA comes free**: point `--ncrna` at a full Rfam `.cm` (includes RF00177 16S / RF02541
  23S / RF00001 5S) and those hits classify as `Rrna` — no separate barrnap needed.

- **CRISPR**: `minced-rs` — a faithful pure-Rust port of MinCED (CRT-derived), no DB.
  Enabled with `--crispr`. On MG1655 it recovers both known K-12 CRISPR arrays near
  iap/ygcB (13- and 7-repeat).

**Conflict resolution** (`resolve`): RNA features mask overlapping CDS and same-locus
multi-model ncRNA hits collapse to the best-scoring call, matching Bakta/Prokka.
CDS-vs-CDS overlaps (operons) and IS/CRISPR pass through untouched.

**Output**: `--gff3 <file>` (written with **noodles**) and/or `--gbff <file>` (GenBank
flat file with the ORIGIN sequence). Without either flag a feature TSV goes to stdout.

Next: richer GFF3/GBFF qualifiers and per-genome locus-tag numbering.

## Usage

```sh
# one bundle from setup-db, then annotate by mode:
bactars setup-db --fast --out DB          # or --full for the PSC/KOfam naming tier
bactars <genome.fna> --db DB --fast       # every non-ML stage
bactars <genome.fna> --db DB --full       # --fast + the learned ML models (+ PSC/KOfam if built)

# or wire individual DBs explicitly:
cargo run --release -- <genome.fna> [--pfam D] [--ncbifams D] [--antifam D] [--amr D] \
    [--trna MODELS_DIR] [--tmrna] [--is-db DIR] [--ncrna CM_DB] [--gff3 FILE] [--gbff FILE]
```

`--fast` runs every non-ML stage; `--full` adds the learned ML models and, when
the PSC/KOfam DBs are present, the FULL-tier protein-similarity naming +
reference-length pseudogene signal. `--trna` takes trnascan-rs's `data/models`
directory and enables bacterial tRNA detection. Emits a TSV of features (contig,
id, coords, strand, kind, best hit) to stdout.

## Test

```sh
cargo test --release
```
