//! bactars — Rust-native, container-free prokaryotic genome annotation.
//!
//! The orchestrator: it wires the Rust subunits (rustygal CDS, rustyhmmer HMM
//! annotation + AMR, trnascan-rs tRNA, oxagorn tmRNA, rust-ise IS elements,
//! infernox ncRNA, minced-rs CRISPR) as in-process libraries over a shared
//! [`feature::Feature`] model, resolves overlaps ([`resolve`]), and serialises to
//! GFF3/GBFF ([`output`]). See `../plan.md`.

pub mod amr_variant;
pub mod cds;
pub mod crispr;
pub mod fasta;
pub mod feature;
pub mod gap;
pub mod gene_score;
pub mod hmm;
pub mod integron;
pub mod is_elements;
pub mod kofam;
pub mod localize;
pub mod mmseqs;
pub mod ncrna;
pub mod oric;
pub mod orit;
pub mod output;
pub mod plasmid;
pub mod prophage;
// TEMP (2026-07-19): the low-precision length/split pseudogene PROXY was retired from
// `--pseudo` (mmseqs alignment detector only). This module is retained ONLY to harvest
// its split-fragment / frameshift / internal-stop logic into the improved detector, then
// DELETE. `allow(dead_code)` until that harvest lands.
#[allow(dead_code)]
pub mod pseudogene;
pub mod pseudogene_align;
pub mod psc;
pub mod build_psc;
pub mod build_kofam;
pub mod qc;
pub mod pipeline;
pub mod resolve;
pub mod setup_db;
pub mod signalox_sp;
pub mod sorf;
pub mod start_refine;
pub mod species;
pub mod tandem;
pub mod tmrna;
pub mod trna;
pub mod util_io;
pub mod vfdb;
pub mod xref;

pub use feature::{Annotation, Feature, FeatureKind, Functional};
pub use pipeline::{run, run_full, Config, HmmDb, RunOutput};
pub use plasmid::PlasmidHit;
pub use species::SpeciesCall;
