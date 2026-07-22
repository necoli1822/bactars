//! The single internal feature model every subunit emits into.
//!
//! Each Rust subunit (rustygal CDS, rustyhmmer HMM annotation, later oxagorn /
//! trnascan-rs / infernox / ARGenus) produces `Feature`s over a shared genome,
//! which the orchestrator gathers, annotates, and finally serialises.

/// What kind of genomic feature this is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FeatureKind {
    Cds,
    Trna,
    Tmrna,
    Rrna,
    Ncrna,
    Crispr,
    IsElement,
    // --- Phase 0 additions (2026-07-16): new top-level feature kinds ---
    /// Replication origin (oriC) — `origin_of_replication`.
    OriC,
    /// Assembly gap (run of N) — `assembly_gap`.
    AssemblyGap,
    /// Tandem repeat region — `repeat_region`.
    TandemRepeat,
    /// Prophage region — `region` / `mobile_genetic_element`.
    Prophage,
    /// Integron — `mobile_genetic_element`.
    Integron,
    /// Origin of transfer (oriT) — `origin_of_transfer` (rep as misc/mobile).
    OriT,
    /// Cis-regulatory RNA element (riboswitch / leader / thermoregulator) —
    /// `regulatory_region`. Classified from the Rfam family type.
    RegulatoryRegion,
}

/// Strand: `+1` forward, `-1` reverse.
pub type Strand = i8;

/// One annotation attached to a feature by an expert subunit (e.g. an HMM hit).
///
/// This is the RAW hit record (many per feature). The resolved, enriched
/// functional annotation used for output lives in [`Functional`] on the feature.
#[derive(Clone, Debug)]
pub struct Annotation {
    /// Producing subunit / database, e.g. `"rustyhmmer:pfam"`.
    pub source: String,
    /// Database accession, e.g. `"PF00069.28"` / `"NF000282.2"`.
    pub accession: String,
    /// Human-readable model name.
    pub name: String,
    /// Bit score.
    pub score: f32,
    /// E-value from the producing statistical search — the HMM sequence/domain
    /// E-value (rustyhmmer), the covariance-model cmsearch E-value (infernox), or
    /// the mmseqs homology E-value — when the search yields one. `None` for
    /// score-only features whose producing method reports no E-value: tandem
    /// repeats (TRF score), oriC (GC-skew), tmRNA (ARAGORN fold energy), tRNA
    /// (tRNAscan-SE bit score) and CRISPR (repeat structure).
    pub evalue: Option<f64>,
    /// Reference (target) protein length in amino acids, when the producing tier
    /// knows it. Only the PSC/UniRep90 similarity tier sets this (from the mmseqs
    /// target length); HMM/xref/RNA tiers leave it `None`. Consumed by the
    /// pseudogene truncation signal as a reference-length source for CDS that have
    /// no resolvable ncbifams `hmm_length`.
    pub ref_len: Option<usize>,
}

/// A structured `/inference` evidence record (INSDC Feature Table v11.4 grammar):
/// `"[CATEGORY:]TYPE[ (same species)][:DB:ACCESSION]"`. Fields are colon-joined
/// at render time; multiple inferences are emitted as SEPARATE `/inference`
/// qualifiers (never combined). `kind` must be one of the fixed INSDC TYPEs
/// (e.g. "similar to AA sequence", "protein motif", "ab initio prediction").
#[derive(Clone, Debug)]
pub struct Inference {
    pub category: Option<String>, // e.g. "COORDINATES"
    pub kind: String,             // fixed INSDC TYPE
    pub same_species: bool,
    pub db: String,       // e.g. "UniProtKB", "AMRFinderPlus", "Prodigal"
    pub accession: String, // e.g. "P0AES4", "S83L", "2.6"
}

/// Signal peptide / lipoprotein prediction (rule-based: lipobox / Tat motif).
#[derive(Clone, Debug)]
pub struct SignalPeptide {
    /// "lipoprotein" (SPaseII lipobox) | "tat" (twin-arginine) | "sec".
    pub kind: String,
    /// Cleavage position (1-based AA index), if determined.
    pub cleavage: Option<usize>,
}

/// Ribosome binding site (Shine-Dalgarno) as predicted by the gene caller.
#[derive(Clone, Debug)]
pub struct Rbs {
    /// Motif string, e.g. "AGGAGG".
    pub motif: String,
    /// Spacer length (nt) between RBS and start codon.
    pub spacer: Option<i32>,
}

/// The resolved, enriched functional annotation of a feature — the layer the
/// output writers read. Default-constructible so every existing `Feature`
/// construction site just adds `func: Functional::default()`.
///
/// xref multi-values are kept as `Vec` in memory; they render to packed
/// `a;b;c` in `meta.sqlite` and to GFF3/GBFF per the AMR/inference spec.
#[derive(Clone, Debug, Default)]
pub struct Functional {
    /// Display product name (`hypothetical protein` if none resolved for a CDS).
    pub product: Option<String>,
    /// Gene symbol (dnaA, gyrA, ...).
    pub gene: Option<String>,
    /// EC numbers.
    pub ec: Vec<String>,
    /// GO terms.
    pub go: Vec<String>,
    /// COG functional-category bitmask (26 categories, one bit each).
    pub cog_cat: u32,
    /// KEGG Orthology id (K#####).
    pub ko: Option<String>,
    /// KEGG pathway/module ids.
    pub pathway: Vec<String>,
    /// Free-text `/note` lines (phenotype, mutation descriptions, ...).
    pub note: Vec<String>,
    /// Structured `/inference` evidence records.
    pub inferences: Vec<Inference>,
    /// Pseudogene flag (frameshift / internal stop / truncated vs a PSC hit).
    pub pseudogene: bool,
    /// INSDC `regulatory_class` for a `regulatory_region` feature (riboswitch,
    /// attenuator, recoding_stimulatory_region, other, …). Set during Rfam
    /// reclassification; `None` for non-regulatory features. NCBI requires this
    /// qualifier on every `regulatory_region` (class `other` additionally needs a
    /// `/note`).
    pub regulatory_class: Option<String>,
    /// Signal peptide / lipoprotein call, if any.
    pub signal_peptide: Option<SignalPeptide>,
    /// Ribosome binding site, if predicted.
    pub rbs: Option<Rbs>,
    /// Predicted subcellular localization (lotrs), e.g. "cytoplasm",
    /// "inner_membrane", "periplasm", "outer_membrane", "extracellular". `None`
    /// unless localization annotation ran and cleared its confidence threshold.
    pub localization: Option<String>,
    /// Learned gene-model quality score in [0,1] from balrox (P(real coding gene)).
    /// `None` unless gene-score annotation ran. Low values flag possible spurious ORFs;
    /// non-destructive (the CDS is kept regardless).
    pub gene_score: Option<f32>,
}

/// A genomic feature on one contig, with any expert annotations.
#[derive(Clone, Debug)]
pub struct Feature {
    pub kind: FeatureKind,
    /// Contig / sequence name this feature lives on.
    pub contig: String,
    /// Feature id, e.g. `"NC_000913.3_42"`.
    pub id: String,
    /// 1-based inclusive start (left coordinate on the contig).
    pub start: i64,
    /// 1-based inclusive end (right coordinate on the contig).
    pub end: i64,
    pub strand: Strand,
    /// Protein translation, for CDS features.
    pub aa: Option<String>,
    /// 5' end incomplete — the feature runs off the contig with no start codon
    /// (interpreted relative to the feature's own strand). Draft-genome /
    /// table2asn partiality: a `<`/`>` location marker in GBFF and
    /// `partial=true`/`start_range`/`end_range` in NCBI GFF3.
    pub partial5: bool,
    /// 3' end incomplete — the feature runs off the contig with no stop codon
    /// (interpreted relative to the feature's own strand).
    pub partial3: bool,
    /// Expert annotations (HMM hits, AMR, homology naming, ...).
    pub annotations: Vec<Annotation>,
    /// Resolved/enriched functional annotation (Phase 0, 2026-07-16).
    pub func: Functional,
}

impl Feature {
    /// Nucleotide length spanned on the contig.
    pub fn len_nt(&self) -> i64 {
        (self.end - self.start + 1).max(0)
    }

    /// Whether this feature has any incomplete end (5' or 3').
    pub fn is_partial(&self) -> bool {
        self.partial5 || self.partial3
    }

    /// Whether the LOW (leftmost / `start`) coordinate is the incomplete end.
    /// On the `+` strand the 5' end is the low coordinate; on the `-` strand the
    /// 3' end is the low coordinate. Drives the GBFF `<` prefix and the GFF3
    /// `start_range=.,<start>` attribute.
    pub fn low_coord_partial(&self) -> bool {
        if self.strand == -1 {
            self.partial3
        } else {
            self.partial5
        }
    }

    /// Whether the HIGH (rightmost / `end`) coordinate is the incomplete end.
    /// On the `+` strand the 3' end is the high coordinate; on the `-` strand the
    /// 5' end is the high coordinate. Drives the GBFF `>` prefix and the GFF3
    /// `end_range=<end>,.` attribute.
    pub fn high_coord_partial(&self) -> bool {
        if self.strand == -1 {
            self.partial5
        } else {
            self.partial3
        }
    }

    /// The best (highest-scoring) annotation, if any.
    pub fn best_annotation(&self) -> Option<&Annotation> {
        self.annotations
            .iter()
            // NaN-safe: a NaN score (e.g. from a degenerate tier) must not panic;
            // `total_cmp` gives a total order over all f32 including NaN.
            .max_by(|a, b| a.score.total_cmp(&b.score))
    }

    /// Display product for output: resolved `func.product` first, else the best
    /// raw annotation's name, else `None` (caller may default to hypothetical).
    pub fn display_product(&self) -> Option<String> {
        if let Some(p) = &self.func.product {
            return Some(p.clone());
        }
        self.best_annotation()
            .map(|a| a.name.clone())
            .filter(|n| !n.is_empty())
    }

    /// Product name for output, defaulting an unannotated CDS to the standard
    /// `hypothetical protein` (NCBI/Bakta convention) rather than nothing. Non-CDS
    /// features without a product return `None` (they carry their own type label).
    pub fn product_for_output(&self) -> Option<String> {
        if let Some(p) = self.display_product() {
            return Some(p);
        }
        if self.kind == FeatureKind::Cds {
            return Some("hypothetical protein".to_string());
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The new `ref_len` field is `None` for the non-PSC tiers and only carries a
    /// value when a tier explicitly records a reference length.
    #[test]
    fn annotation_ref_len_default_and_set() {
        let hmm = Annotation {
            source: "rustyhmmer:ncbifams".to_string(),
            accession: "NF000282.2".to_string(),
            name: "dnaX".to_string(),
            score: 250.0,
            evalue: Some(1e-40),
            ref_len: None,
        };
        assert_eq!(hmm.ref_len, None);

        let psc = Annotation {
            source: "psc:uniref90".to_string(),
            accession: "UniRef90_ABC".to_string(),
            name: "real product".to_string(),
            score: 150.0,
            evalue: Some(1e-20),
            ref_len: Some(300),
        };
        assert_eq!(psc.ref_len, Some(300));
    }
}
