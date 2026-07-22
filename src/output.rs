//! Serialisation of the resolved feature set.
//!
//! GFF3 is written with **noodles** (`noodles-gff`), per the plan's mandate that
//! output goes through noodles. GBFF (GenBank flat file) is emitted with a small
//! hand-rolled writer — noodles has no GenBank writer — including the ORIGIN
//! sequence so the record round-trips through Biopython/Artemis.
//!
//! # NCBI-submission-grade output
//!
//! Two output tiers are provided:
//!
//! * [`write_gff3`] / [`write_gbff`] — backward-compatible signatures (called by
//!   `bin/bactars.rs`). Both now emit the clean `locus_tag` scheme, sanitised
//!   product names, percent-encoded attribute values, and `/pseudo` handling.
//!   GBFF additionally emits the `gene`→feature block hierarchy in place.
//! * [`write_gff3_submission`] / [`write_gbff_submission`] — the *submission*
//!   variants. These take the contig sequences (for `##sequence-region` lengths
//!   and the trailing `##FASTA` block in GFF3) plus a configurable `locus_tag`
//!   prefix, and emit the full `gene`→`CDS`/RNA parent/child hierarchy.
//!
//! The orchestrator should call [`write_gff3_submission`] (self-contained GFF3
//! with embedded FASTA) and [`write_gbff`]/[`write_gbff_submission`] for
//! NCBI-ready deposits. See the function docs for exact signatures.

use std::io::{self, Write};

use noodles_core::Position;
use noodles_gff as gff;

use gff::feature::record::{Phase, Strand};
use gff::feature::record_buf::attributes::field::Value as AttrValue;
use gff::feature::record_buf::Attributes;
use gff::feature::RecordBuf;

use crate::feature::{Feature, FeatureKind, Inference};

/// Default `locus_tag` prefix when the caller does not supply one. NCBI issues a
/// registered prefix per BioProject; `BACT` is a sensible placeholder.
pub const DEFAULT_LOCUS_PREFIX: &str = "BACT";

/// GFF3 / GenBank feature key for a feature kind (Sequence-Ontology-ish terms
/// that both GFF3 and GenBank accept).
fn feature_type(kind: FeatureKind) -> &'static str {
    match kind {
        FeatureKind::Cds => "CDS",
        FeatureKind::Trna => "tRNA",
        FeatureKind::Tmrna => "tmRNA",
        FeatureKind::Rrna => "rRNA",
        FeatureKind::Ncrna => "ncRNA",
        FeatureKind::Crispr => "repeat_region",
        FeatureKind::IsElement => "mobile_genetic_element",
        FeatureKind::OriC => "origin_of_replication",
        FeatureKind::AssemblyGap => "assembly_gap",
        FeatureKind::TandemRepeat => "repeat_region",
        FeatureKind::Prophage => "region",
        FeatureKind::Integron => "mobile_genetic_element",
        FeatureKind::OriT => "origin_of_transfer",
        FeatureKind::RegulatoryRegion => "regulatory_region",
    }
}

/// Whether this feature kind is a *gene* (gets a parent `gene` feature carrying
/// the `locus_tag`, per NCBI/Prokka submission structure). Structural features
/// (repeats, mobile elements, origins, gaps, regulatory regions) stand alone.
fn is_gene_associated(kind: FeatureKind) -> bool {
    matches!(
        kind,
        FeatureKind::Cds
            | FeatureKind::Trna
            | FeatureKind::Tmrna
            | FeatureKind::Rrna
            | FeatureKind::Ncrna
    )
}

fn gff_strand(strand: i8) -> Strand {
    match strand {
        1 => Strand::Forward,
        -1 => Strand::Reverse,
        _ => Strand::None,
    }
}

fn pos(v: i64) -> Position {
    // Features are 1-based inclusive; clamp defensively so a stray 0/negative
    // never panics the writer.
    Position::try_from(v.max(1) as usize).expect("1-based position")
}

/// Light NCBI product-name hygiene: collapse internal whitespace, trim, and drop
/// a trailing period (NCBI rejects `/product` values ending in `.`). Conservative
/// on purpose — it does NOT rewrite "putative"/"homolog"/all-caps content, only
/// normalises whitespace and punctuation.
fn sanitize_product(p: &str) -> String {
    // Collapse any run of ASCII whitespace to a single space.
    let collapsed = p.split_whitespace().collect::<Vec<_>>().join(" ");
    // Drop trailing period(s) and any whitespace they leave behind.
    collapsed.trim_end_matches('.').trim_end().to_string()
}

/// Zero-pad width for `locus_tag` numbers, sized to the feature count but never
/// narrower than 5 digits (NCBI/Prokka convention, e.g. `BACT_00001`).
fn locus_tag_width(count: usize) -> usize {
    count.to_string().len().max(5)
}

/// Assign a clean, sequential `locus_tag` to each feature in the (already
/// coordinate-sorted) slice: `<PREFIX>_<zero-padded index>`, 1-based. The
/// returned vector is parallel to `features`.
pub fn assign_locus_tags(features: &[Feature], prefix: &str) -> Vec<String> {
    let width = locus_tag_width(features.len());
    (1..=features.len())
        .map(|n| format!("{prefix}_{n:0width$}"))
        .collect()
}

/// db_xref values valid as INSDC `/db_xref` (and GFF3 `Dbxref`). Only
/// **NCBI-registered** databases may appear here (https://www.ncbi.nlm.nih.gov/
/// genbank/collab/db_xref/), otherwise table2asn rejects every value with
/// `IllegalDbXref`. bactars currently populates exactly one registered database
/// this way: **RFAM** (the Rfam accession attached to ncRNA/rRNA/regulatory hits
/// by the infernox scanner). Everything else moves elsewhere:
///   * HMM / source-protein hits (`rustyhmmer:*`) → `/inference` ([`inference_values`]).
///   * EC numbers → `/EC_number` (GBFF) / `EC_number=` (GFF3).
///   * KEGG Orthology / anticodon → `/note` (KEGG is not a registered db_xref;
///     an anticodon is not a db_xref at all).
///   * GO terms → `Ontology_term` (GFF3) / `/db_xref="GO:..."` (GBFF), handled at
///     the call sites (GO *is* registered).
fn dbxref_values(f: &Feature) -> Vec<String> {
    let mut xrefs: Vec<String> = Vec::new();
    if let Some(a) = f.best_annotation() {
        // Rfam models come in under the `infernox` source with an `RF#####`
        // accession → the registered `RFAM` db_xref.
        if a.source == "infernox" && a.accession.starts_with("RF") && a.accession != "-" {
            xrefs.push(format!("RFAM:{}", a.accession));
        }
    }
    xrefs
}

/// INSDC `/inference` (GFF3 `inference=`) evidence strings for a feature. This is
/// where HMM / source-protein matches belong — NOT `/db_xref`. Combines the best
/// raw HMM hit (`rustyhmmer:*` → `protein motif:HMM:<DB>:<acc>`) with any
/// structured [`Inference`] records already recorded on the feature (VFDB, AMR…).
fn inference_values(f: &Feature) -> Vec<String> {
    let mut infs: Vec<String> = Vec::new();
    if let Some(a) = f.best_annotation() {
        if let Some(db) = hmm_inference_db(&a.source) {
            if !a.accession.is_empty() && a.accession != "-" {
                infs.push(format!("protein motif:HMM:{db}:{}", a.accession));
            }
        }
    }
    for inf in &f.func.inferences {
        infs.push(render_inference(inf));
    }
    infs
}

/// Evidence-database token for an HMM `/inference` from a `rustyhmmer:<db>`
/// source label. `None` for non-HMM sources (they are not protein-motif hits).
fn hmm_inference_db(source: &str) -> Option<String> {
    let db = source.strip_prefix("rustyhmmer:")?;
    Some(match db {
        "ncbifams" => "NCBIfam".to_string(),
        "pfam" => "Pfam".to_string(),
        other => {
            // Title-case as a reasonable evidence-db token.
            let mut c = other.chars();
            match c.next() {
                Some(first) => first.to_uppercase().collect::<String>() + c.as_str(),
                None => other.to_string(),
            }
        }
    })
}

/// Render a structured [`Inference`] into the INSDC feature-table grammar
/// `"[CATEGORY: ]TYPE[ (same species)][:DB:ACCESSION]"`.
fn render_inference(inf: &Inference) -> String {
    let mut s = String::new();
    if let Some(cat) = &inf.category {
        if !cat.is_empty() {
            s.push_str(cat);
            s.push_str(": ");
        }
    }
    s.push_str(&inf.kind);
    if inf.same_species {
        s.push_str(" (same species)");
    }
    if !inf.db.is_empty() {
        s.push(':');
        s.push_str(&inf.db);
        if !inf.accession.is_empty() {
            s.push(':');
            s.push_str(&inf.accession);
        }
    }
    s
}

/// Output-ready `/product` for a feature, with rRNA and tRNA model-name
/// normalisation applied (NCBI requires standard product strings). `None` when
/// the feature carries no product.
fn output_product(f: &Feature) -> Option<String> {
    let p = sanitize_product(&f.product_for_output()?);
    match f.kind {
        FeatureKind::Rrna => Some(normalize_rrna_product(&p)),
        FeatureKind::Trna => Some(normalize_trna_product(&p)),
        _ => Some(p),
    }
}

/// True when the (already-normalised) product is the generic hypothetical label.
/// A `hypothetical protein` CDS must NOT also carry a gene symbol.
fn is_hypothetical(product: Option<&str>) -> bool {
    matches!(product, Some(p) if p.eq_ignore_ascii_case("hypothetical protein"))
}

/// The gene symbol to emit for a feature: suppressed on `hypothetical protein`
/// CDS (NCBI flags a hypothetical CDS carrying a gene name).
fn output_gene<'a>(f: &'a Feature, product: Option<&str>) -> Option<&'a String> {
    if is_hypothetical(product) {
        return None;
    }
    f.func.gene.as_ref().filter(|g| !g.is_empty())
}

/// Map an Rfam rRNA model name to the standard NCBI product string. Unknown
/// names fall back to a coarse SSU/LSU/5S heuristic, else are kept verbatim.
fn normalize_rrna_product(name: &str) -> String {
    let key = name.trim();
    let mapped = match key {
        "SSU_rRNA_bacteria" | "SSU_rRNA_archaea" => Some("16S ribosomal RNA"),
        "LSU_rRNA_bacteria" | "LSU_rRNA_archaea" => Some("23S ribosomal RNA"),
        "5S_rRNA" => Some("5S ribosomal RNA"),
        "5_8S_rRNA" => Some("5.8S ribosomal RNA"),
        "SSU_rRNA_eukarya" => Some("18S ribosomal RNA"),
        "LSU_rRNA_eukarya" => Some("28S ribosomal RNA"),
        _ => None,
    };
    if let Some(m) = mapped {
        return m.to_string();
    }
    if key.starts_with("SSU_rRNA") {
        "16S ribosomal RNA".to_string()
    } else if key.starts_with("LSU_rRNA") {
        "23S ribosomal RNA".to_string()
    } else if key.starts_with("5_8S") {
        "5.8S ribosomal RNA".to_string()
    } else if key.starts_with("5S") {
        "5S ribosomal RNA".to_string()
    } else {
        key.to_string()
    }
}

/// Normalise a tRNA product so table2asn can always parse an amino acid. A blank,
/// bare `tRNA`, or undetermined/pseudo isotype becomes `tRNA-Xxx` (NCBI's
/// undetermined-amino-acid convention); a determined isotype is kept.
fn normalize_trna_product(p: &str) -> String {
    let iso = p
        .trim()
        .strip_prefix("tRNA-")
        .map(str::trim)
        .unwrap_or("");
    if iso.is_empty()
        || iso.eq_ignore_ascii_case("Undet")
        || iso.eq_ignore_ascii_case("Pseudo")
        || iso.eq_ignore_ascii_case("Sup")
        || iso.eq_ignore_ascii_case("Xxx")
    {
        return "tRNA-Xxx".to_string();
    }
    format!("tRNA-{iso}")
}

/// Free-text `/note` values: COG functional category and KEGG pathways are NOT
/// registered `/db_xref` databases (the letters/ids are not accessions), so they
/// go to `note` per INSDC rules, alongside any curated free-text notes.
fn note_values(f: &Feature) -> Vec<String> {
    let mut notes: Vec<String> = Vec::new();
    if f.func.cog_cat != 0 {
        notes.push(format!(
            "COG_category:{}",
            crate::kofam::cog_letters(f.func.cog_cat)
        ));
    }
    if !f.func.pathway.is_empty() {
        notes.push(format!("KEGG_pathway:{}", f.func.pathway.join(" ")));
    }
    // KEGG Orthology is not a registered /db_xref database → keep it as a note.
    if let Some(ko) = &f.func.ko {
        if !ko.is_empty() {
            notes.push(format!("KEGG:{ko}"));
        }
    }
    // tRNA anticodon provenance (previously mis-emitted as an illegal db_xref).
    if f.kind == FeatureKind::Trna {
        if let Some(a) = f.annotations.iter().find(|a| a.source == "trnascan-rs") {
            if !a.accession.is_empty() {
                notes.push(format!("anticodon:{}", a.accession));
            }
        }
    }
    // regulatory_region needs a note when regulatory_class is `other`; always
    // record the Rfam family so the requirement is met regardless of class.
    if f.kind == FeatureKind::RegulatoryRegion {
        if let Some(a) = f.best_annotation() {
            if !a.name.is_empty() {
                notes.push(format!("Rfam:{}", a.name));
            }
        }
    }
    // Predicted signal peptide (bactars signalpep). This was computed on the CDS but
    // never serialized — surface the class (+ cleavage) as a /note so `--signalpep`
    // actually produces output. Kept a note (not a sig_peptide feature) since it is a
    // rule-based prediction.
    if let Some(sp) = &f.func.signal_peptide {
        let cl = sp
            .cleavage
            .map(|c| format!(", cleavage after aa {c}"))
            .unwrap_or_default();
        notes.push(format!("predicted signal peptide ({}{})", sp.kind, cl));
    }
    for n in &f.func.note {
        if !n.is_empty() {
            notes.push(n.clone());
        }
    }
    notes
}

/// Build the attribute set shared by a feature (CDS/RNA/…) GFF3 record.
/// `locus_tag` is the clean scheme; `parent` links a child to its `gene` record.
/// The original feature id is preserved as the record `ID`.
fn feature_attributes(f: &Feature, locus_tag: &str, parent: Option<&str>) -> Attributes {
    let mut attributes = Attributes::default();
    let map = attributes.as_mut();

    // Keep the original (contig-embedded) id as the GFF3 ID; the clean scheme
    // lives in the dedicated locus_tag qualifier.
    map.insert("ID".into(), AttrValue::String(f.id.clone().into()));
    // An `assembly_gap` MUST NOT carry a /locus_tag — table2asn rejects it there.
    if f.kind != FeatureKind::AssemblyGap {
        map.insert("locus_tag".into(), AttrValue::String(locus_tag.into()));
    }
    if let Some(p) = parent {
        map.insert("Parent".into(), AttrValue::String(p.into()));
    }

    let product = output_product(f);
    if let Some(product) = &product {
        if !product.is_empty() {
            map.insert("Name".into(), AttrValue::String(product.clone().into()));
            map.insert("product".into(), AttrValue::String(product.clone().into()));
        }
    }
    if let Some(gene) = output_gene(f, product.as_deref()) {
        map.insert("gene".into(), AttrValue::String(gene.clone().into()));
    }

    // Genetic code 11 on every CDS: without it table2asn misreads GTG/TTG starts
    // as 5'-partial (RefSeq emits transl_table=11 on every CDS). /codon_start=1
    // pins the reading frame (features carry the full CDS from the first base).
    if f.kind == FeatureKind::Cds {
        map.insert("codon_start".into(), AttrValue::String("1".into()));
        map.insert("transl_table".into(), AttrValue::String("11".into()));
    }

    // assembly_gap: INSDC-required gap qualifiers (and never a /locus_tag, above).
    if f.kind == FeatureKind::AssemblyGap {
        map.insert(
            "estimated_length".into(),
            AttrValue::String(f.len_nt().to_string().into()),
        );
        map.insert(
            "gap_type".into(),
            AttrValue::String("within scaffold".into()),
        );
        map.insert(
            "linkage_evidence".into(),
            AttrValue::String("paired-ends".into()),
        );
    }

    // mobile_genetic_element requires a /mobile_element_type (controlled vocab).
    match f.kind {
        FeatureKind::IsElement => {
            map.insert(
                "mobile_element_type".into(),
                AttrValue::String("insertion sequence".into()),
            );
        }
        FeatureKind::Integron => {
            map.insert(
                "mobile_element_type".into(),
                AttrValue::String("integron".into()),
            );
        }
        _ => {}
    }

    // repeat_region requires an /rpt_type (INSDC controlled vocab); CRISPR arrays
    // additionally carry /rpt_family=CRISPR.
    match f.kind {
        FeatureKind::Crispr => {
            map.insert("rpt_type".into(), AttrValue::String("direct".into()));
            map.insert("rpt_family".into(), AttrValue::String("CRISPR".into()));
        }
        FeatureKind::TandemRepeat => {
            map.insert("rpt_type".into(), AttrValue::String("tandem".into()));
        }
        _ => {}
    }

    // regulatory_region requires a regulatory_class (INSDC controlled vocab).
    if f.kind == FeatureKind::RegulatoryRegion {
        let class = f
            .func
            .regulatory_class
            .clone()
            .unwrap_or_else(|| "other".to_string());
        map.insert("regulatory_class".into(), AttrValue::String(class.into()));
    }

    // Multi-value attributes MUST use Array so noodles emits real comma
    // separators and percent-encodes only WITHIN each element (a joined String
    // would have its separating commas encoded to %2C, corrupting the list).
    let xrefs = dbxref_values(f);
    if !xrefs.is_empty() {
        map.insert(
            "Dbxref".into(),
            AttrValue::Array(xrefs.into_iter().map(|s| s.into()).collect()),
        );
    }
    // EC numbers are NOT a db_xref database → their own EC_number attribute.
    let ecs: Vec<String> = f.func.ec.iter().filter(|e| !e.is_empty()).cloned().collect();
    if !ecs.is_empty() {
        map.insert(
            "EC_number".into(),
            AttrValue::Array(ecs.into_iter().map(|s| s.into()).collect()),
        );
    }
    // HMM / source-hit and structured evidence → inference (never db_xref).
    let infs = inference_values(f);
    if !infs.is_empty() {
        map.insert(
            "inference".into(),
            AttrValue::Array(infs.into_iter().map(|s| s.into()).collect()),
        );
    }
    if !f.func.go.is_empty() {
        map.insert(
            "Ontology_term".into(),
            AttrValue::Array(f.func.go.iter().map(|s| s.clone().into()).collect()),
        );
    }
    let notes = note_values(f);
    if !notes.is_empty() {
        map.insert(
            "Note".into(),
            AttrValue::Array(notes.into_iter().map(|s| s.into()).collect()),
        );
    }

    // Pseudogene: mark and (elsewhere) suppress the translation.
    if f.func.pseudogene {
        map.insert("pseudo".into(), AttrValue::String("true".into()));
    }

    // Draft-genome partiality (NCBI GFF3 convention): a feature running off a
    // contig end carries `partial=true` plus a coordinate-column range qualifier
    // — `start_range=.,<start>` when the LOW coordinate is incomplete,
    // `end_range=<end>,.` when the HIGH coordinate is incomplete. The 5'/3'→
    // coordinate mapping is strand-aware, so a minus-strand 5'-partial CDS emits
    // `end_range`. A complete feature adds nothing here (byte-identical output).
    if f.is_partial() {
        map.insert("partial".into(), AttrValue::String("true".into()));
        // `start_range`/`end_range` values carry a LITERAL comma (`.,N` / `N,.`);
        // an `Array` makes noodles emit a real comma separator instead of
        // percent-encoding it (`%2C`), which would corrupt the NCBI range grammar.
        if f.low_coord_partial() {
            map.insert(
                "start_range".into(),
                AttrValue::Array(vec![".".into(), f.start.to_string().into()]),
            );
        }
        if f.high_coord_partial() {
            map.insert(
                "end_range".into(),
                AttrValue::Array(vec![f.end.to_string().into(), ".".into()]),
            );
        }
    }

    attributes
}

/// Build the feature (CDS/RNA/…) GFF3 record.
fn build_feature_record(f: &Feature, locus_tag: &str, parent: Option<&str>) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(f.contig.clone())
        .set_source("bactars")
        .set_type(feature_type(f.kind))
        .set_start(pos(f.start))
        .set_end(pos(f.end.max(f.start)))
        .set_strand(gff_strand(f.strand))
        .set_attributes(feature_attributes(f, locus_tag, parent));

    if let Some(a) = f.best_annotation() {
        if a.score.is_finite() {
            builder = builder.set_score(a.score);
        }
    }
    // GFF3 phase is only meaningful for CDS; use 0 (features carry the full CDS).
    if f.kind == FeatureKind::Cds {
        builder = builder.set_phase(Phase::Zero);
    }
    builder.build()
}

/// Build the parent `gene` GFF3 record for a gene-associated feature.
fn build_gene_record(f: &Feature, locus_tag: &str, gene_id: &str) -> RecordBuf {
    let mut attributes = Attributes::default();
    {
        let map = attributes.as_mut();
        map.insert("ID".into(), AttrValue::String(gene_id.into()));
        map.insert("locus_tag".into(), AttrValue::String(locus_tag.into()));
        // Suppress the gene symbol on a hypothetical-protein CDS (NCBI flags a
        // hypothetical CDS that also carries a gene name).
        if let Some(gene) = output_gene(f, output_product(f).as_deref()) {
            map.insert("Name".into(), AttrValue::String(gene.clone().into()));
            map.insert("gene".into(), AttrValue::String(gene.clone().into()));
        }
        if f.func.pseudogene {
            map.insert("pseudo".into(), AttrValue::String("true".into()));
        }
        // Partiality mirrors the child feature (see `feature_attributes`); the
        // parent gene must carry the same partial/range flags for table2asn.
        if f.is_partial() {
            map.insert("partial".into(), AttrValue::String("true".into()));
            if f.low_coord_partial() {
                map.insert(
                    "start_range".into(),
                    AttrValue::Array(vec![".".into(), f.start.to_string().into()]),
                );
            }
            if f.high_coord_partial() {
                map.insert(
                    "end_range".into(),
                    AttrValue::Array(vec![f.end.to_string().into(), ".".into()]),
                );
            }
        }
    }
    RecordBuf::builder()
        .set_reference_sequence_name(f.contig.clone())
        .set_source("bactars")
        .set_type("gene")
        .set_start(pos(f.start))
        .set_end(pos(f.end.max(f.start)))
        .set_strand(gff_strand(f.strand))
        .set_attributes(attributes)
        .build()
}

/// Write the feature set as GFF3 via noodles (flat: one record per feature).
///
/// Backward-compatible signature. Emits the clean `locus_tag`, sanitised
/// products, percent-encoded values, and `pseudo=true` where applicable, but
/// **not** the `##sequence-region`/`##FASTA` directives or the `gene` parent
/// hierarchy — use [`write_gff3_submission`] for NCBI-ready output.
pub fn write_gff3<W: Write>(writer: W, features: &[Feature]) -> io::Result<()> {
    let mut w = gff::io::Writer::new(writer);
    // Version pragma (noodles writes records only; the header line is ours).
    w.get_mut().write_all(b"##gff-version 3\n")?;
    let tags = assign_locus_tags(features, DEFAULT_LOCUS_PREFIX);
    for (f, lt) in features.iter().zip(&tags) {
        w.write_record(&build_feature_record(f, lt, None))?;
    }
    Ok(())
}

/// Write a **submission-grade**, self-contained GFF3 (NCBI/Prokka convention):
///
/// * `##gff-version 3`
/// * one `##sequence-region <seqid> 1 <len>` per contig
/// * for each gene-associated feature, a `gene` record followed by the
///   `CDS`/`tRNA`/… child carrying `Parent=gene-<locus_tag>`; standalone
///   structural features are emitted directly
/// * a trailing `##FASTA` block with every contig sequence
///
/// `contigs` supplies `(name, nucleotides)` for the region lengths and FASTA.
/// `prefix` is the `locus_tag` prefix (pass [`DEFAULT_LOCUS_PREFIX`] for `BACT`).
pub fn write_gff3_submission<W: Write>(
    writer: W,
    features: &[Feature],
    contigs: &[(String, Vec<u8>)],
    prefix: &str,
) -> io::Result<()> {
    let mut w = gff::io::Writer::new(writer);
    w.get_mut().write_all(b"##gff-version 3\n")?;
    for (name, seq) in contigs {
        writeln!(w.get_mut(), "##sequence-region {} 1 {}", name, seq.len())?;
    }

    let tags = assign_locus_tags(features, prefix);
    for (f, lt) in features.iter().zip(&tags) {
        if is_gene_associated(f.kind) {
            let gene_id = format!("gene-{lt}");
            w.write_record(&build_gene_record(f, lt, &gene_id))?;
            w.write_record(&build_feature_record(f, lt, Some(&gene_id)))?;
        } else {
            w.write_record(&build_feature_record(f, lt, None))?;
        }
    }

    // Embedded FASTA makes the GFF3 self-contained (NCBI/Prokka `##FASTA`).
    let inner = w.get_mut();
    inner.write_all(b"##FASTA\n")?;
    for (name, seq) in contigs {
        write_fasta_record(inner, name, seq)?;
    }
    Ok(())
}

/// Write one FASTA record wrapped to 70 columns.
fn write_fasta_record<W: Write>(writer: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writeln!(writer, ">{name}")?;
    for chunk in seq.chunks(70) {
        writer.write_all(chunk)?;
        writeln!(writer)?;
    }
    Ok(())
}

/// Write the feature set as a GenBank flat file (GBFF), one LOCUS per contig,
/// including the ORIGIN sequence. `contigs` supplies each contig's nucleotides.
///
/// Backward-compatible signature; uses [`DEFAULT_LOCUS_PREFIX`]. Emits the
/// `gene`→feature block hierarchy, clean `locus_tag`s, sanitised products, and
/// `/pseudo` handling. For a configurable prefix use [`write_gbff_submission`].
pub fn write_gbff<W: Write>(
    writer: W,
    features: &[Feature],
    contigs: &[(String, Vec<u8>)],
) -> io::Result<()> {
    write_gbff_submission(writer, features, contigs, DEFAULT_LOCUS_PREFIX)
}

/// GBFF writer with a configurable `locus_tag` prefix. See [`write_gbff`].
/// GBFF feature-table location string for a feature, honouring draft-genome
/// partiality. A complete feature (both ends present) renders EXACTLY as before
/// (`start..end` / `complement(start..end)`); a partial end prefixes the LOW
/// coordinate with `<` and/or the HIGH coordinate with `>` (INSDC convention).
/// The 5'/3'→coordinate mapping is strand-aware (see `Feature::low_coord_partial`
/// / `high_coord_partial`): on a `complement(...)` location the 5' end is the
/// high coordinate, so a 5'-partial minus-strand CDS marks the `>` end.
fn feature_location(f: &Feature) -> String {
    let lo = if f.low_coord_partial() {
        format!("<{}", f.start)
    } else {
        f.start.to_string()
    };
    let hi = if f.high_coord_partial() {
        format!(">{}", f.end)
    } else {
        f.end.to_string()
    };
    if f.strand == -1 {
        format!("complement({lo}..{hi})")
    } else {
        format!("{lo}..{hi}")
    }
}

pub fn write_gbff_submission<W: Write>(
    mut writer: W,
    features: &[Feature],
    contigs: &[(String, Vec<u8>)],
    prefix: &str,
) -> io::Result<()> {
    let tags = assign_locus_tags(features, prefix);
    // No organism/config is threaded into the writer; default to a valid,
    // submission-safe placeholder rather than crashing or emitting nothing.
    let organism = "Bacteria";

    for (name, seq) in contigs {
        let len = seq.len();
        // The LOCUS name field is fixed-width; an over-long contig name must be
        // TRUNCATED (not just padded) so downstream fixed-column parsers align.
        let locus_name: String = name.chars().take(16).collect();
        writeln!(
            writer,
            "LOCUS       {:<16} {} bp    DNA     linear   BCT 01-JAN-2025",
            locus_name, len
        )?;
        writeln!(writer, "DEFINITION  bactars annotation of {name}.")?;
        writeln!(writer, "ACCESSION   {name}")?;
        writeln!(writer, "SOURCE      .")?;
        writeln!(writer, "FEATURES             Location/Qualifiers")?;
        writeln!(writer, "     source          1..{len}")?;
        writeln!(writer, "                     /organism=\"{organism}\"")?;
        writeln!(writer, "                     /mol_type=\"genomic DNA\"")?;

        for (i, f) in features
            .iter()
            .enumerate()
            .filter(|(_, f)| &f.contig == name)
        {
            let lt = &tags[i];
            let loc = feature_location(f);

            let product = output_product(f);
            let gene = output_gene(f, product.as_deref());

            // Parent gene block (NCBI/Prokka submission structure).
            if is_gene_associated(f.kind) {
                writeln!(writer, "     {:<15} {}", "gene", loc)?;
                writeln!(writer, "                     /locus_tag=\"{lt}\"")?;
                if let Some(gene) = gene {
                    writeln!(writer, "                     /gene=\"{gene}\"")?;
                }
                if f.func.pseudogene {
                    writeln!(writer, "                     /pseudo")?;
                }
            }

            // Feature block.
            writeln!(writer, "     {:<15} {}", feature_type(f.kind), loc)?;
            // An `assembly_gap` MUST NOT carry a /locus_tag (table2asn rejects it).
            if f.kind != FeatureKind::AssemblyGap {
                writeln!(writer, "                     /locus_tag=\"{lt}\"")?;
            }
            // assembly_gap: INSDC-required gap qualifiers.
            if f.kind == FeatureKind::AssemblyGap {
                writeln!(writer, "                     /estimated_length={}", f.len_nt())?;
                writeln!(writer, "                     /gap_type=\"within scaffold\"")?;
                writeln!(writer, "                     /linkage_evidence=\"paired-ends\"")?;
            }
            // mobile_genetic_element requires /mobile_element_type.
            match f.kind {
                FeatureKind::IsElement => writeln!(
                    writer,
                    "                     /mobile_element_type=\"insertion sequence\""
                )?,
                FeatureKind::Integron => writeln!(
                    writer,
                    "                     /mobile_element_type=\"integron\""
                )?,
                _ => {}
            }
            // repeat_region requires /rpt_type (controlled vocab); CRISPR arrays
            // additionally carry /rpt_family=CRISPR.
            match f.kind {
                FeatureKind::Crispr => {
                    writeln!(writer, "                     /rpt_type=direct")?;
                    writeln!(writer, "                     /rpt_family=\"CRISPR\"")?;
                }
                FeatureKind::TandemRepeat => {
                    writeln!(writer, "                     /rpt_type=tandem")?;
                }
                _ => {}
            }
            if let Some(gene) = gene {
                writeln!(writer, "                     /gene=\"{gene}\"")?;
            }
            if let Some(product) = &product {
                if !product.is_empty() {
                    writeln!(writer, "                     /product=\"{product}\"")?;
                }
            }
            // regulatory_region requires a regulatory_class (INSDC controlled vocab).
            if f.kind == FeatureKind::RegulatoryRegion {
                let class = f
                    .func
                    .regulatory_class
                    .clone()
                    .unwrap_or_else(|| "other".to_string());
                writeln!(writer, "                     /regulatory_class=\"{class}\"")?;
            }
            for ec in &f.func.ec {
                if !ec.is_empty() {
                    writeln!(writer, "                     /EC_number=\"{ec}\"")?;
                }
            }
            // HMM / source-hit + structured evidence → /inference (never db_xref).
            for inf in inference_values(f) {
                write_wrapped_qualifier(&mut writer, "inference", &inf)?;
            }
            // Only NCBI-registered databases in db_xref (RFAM) + GO.
            for x in dbxref_values(f) {
                writeln!(writer, "                     /db_xref=\"{x}\"")?;
            }
            for go in &f.func.go {
                if !go.is_empty() {
                    writeln!(writer, "                     /db_xref=\"{go}\"")?;
                }
            }
            let notes = note_values(f);
            if !notes.is_empty() {
                write_wrapped_qualifier(&mut writer, "note", &notes.join("; "))?;
            }
            // Genetic code 11 + fixed reading frame on every CDS (matches RefSeq;
            // GTG/TTG starts). /codon_start=1: prodigal/rustygal reports every gene
            // (including 5'-partial ones) trimmed to its own frame boundary, so the
            // translation always begins on the first base of the reported start
            // coordinate — the frame starts at 1 relative to that coordinate even
            // when 5'-partial. (rustygal does not expose a sub-codon 5' offset in the
            // header attrs, so a non-1 codon_start is not derivable here.)
            if f.kind == FeatureKind::Cds {
                writeln!(writer, "                     /codon_start=1")?;
                writeln!(writer, "                     /transl_table=11")?;
            }
            // Pseudogene: mark and suppress the translation (NCBI rule).
            if f.func.pseudogene {
                writeln!(writer, "                     /pseudo")?;
            } else if let Some(aa) = &f.aa {
                // A translation never includes the terminal stop codon.
                let aa = aa.strip_suffix('*').unwrap_or(aa);
                let aa = aa.strip_suffix('.').unwrap_or(aa);
                write_wrapped_qualifier(&mut writer, "translation", aa)?;
            }
        }

        writeln!(writer, "ORIGIN")?;
        write_origin(&mut writer, seq)?;
        writeln!(writer, "//")?;
    }
    Ok(())
}

/// Write a long qualifier value (e.g. /translation) wrapped to the GenBank
/// 79-column layout (58 value chars per continuation line).
fn write_wrapped_qualifier<W: Write>(writer: &mut W, key: &str, value: &str) -> io::Result<()> {
    const INDENT: &str = "                     ";
    let head = format!("/{key}=\"");

    // Split the value into per-line chunks on CHAR boundaries so a multi-byte
    // qualifier (α, β, ×, accented letters in product/note) is never sliced
    // through a UTF-8 boundary and silently dropped. The per-line budget counts
    // BYTES (as GenBank column width does), but a char is only added to the
    // current line if it fits whole — otherwise it starts the next line. For
    // pure-ASCII values this reproduces the previous byte-exact wrapping.
    let first_budget = 58usize.saturating_sub(head.len());
    let mut lines: Vec<String> = Vec::new();
    let mut current = String::new();
    let mut used = 0usize;
    let mut budget = first_budget;
    for ch in value.chars() {
        let w = ch.len_utf8();
        if used + w > budget && !current.is_empty() {
            lines.push(std::mem::take(&mut current));
            used = 0;
            budget = 58;
        }
        current.push(ch);
        used += w;
    }
    lines.push(current); // always emit at least one (possibly empty) line

    let last = lines.len() - 1;
    for (idx, chunk) in lines.iter().enumerate() {
        write!(writer, "{INDENT}")?;
        if idx == 0 {
            write!(writer, "{head}")?;
        }
        write!(writer, "{chunk}")?;
        if idx == last {
            write!(writer, "\"")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

/// Write the ORIGIN block: 60 bases per line in six space-separated groups of
/// ten, prefixed by a right-justified 1-based coordinate.
fn write_origin<W: Write>(writer: &mut W, seq: &[u8]) -> io::Result<()> {
    let lower: Vec<u8> = seq.iter().map(|b| b.to_ascii_lowercase()).collect();
    let mut i = 0;
    while i < lower.len() {
        write!(writer, "{:>9}", i + 1)?;
        let line_end = (i + 60).min(lower.len());
        let mut j = i;
        while j < line_end {
            let grp_end = (j + 10).min(line_end);
            write!(writer, " ")?;
            writer.write_all(&lower[j..grp_end])?;
            j = grp_end;
        }
        writeln!(writer)?;
        i = line_end;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Annotation, Feature, FeatureKind, Functional};

    fn cds(id: &str, start: i64, end: i64) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "contig1".to_string(),
            id: id.to_string(),
            start,
            end,
            strand: 1,
            aa: Some("MKV".to_string()),
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        }
    }

    #[test]
    fn locus_tags_zero_padded_and_sequential() {
        let feats = vec![cds("a", 1, 9), cds("b", 20, 29), cds("c", 40, 49)];
        let tags = assign_locus_tags(&feats, "BACT");
        assert_eq!(tags, vec!["BACT_00001", "BACT_00002", "BACT_00003"]);
    }

    #[test]
    fn locus_tag_width_scales_with_count() {
        // Fewer than 100k features → 5-digit minimum.
        assert_eq!(locus_tag_width(3), 5);
        assert_eq!(locus_tag_width(99_999), 5);
        // More → widen to fit.
        assert_eq!(locus_tag_width(1_000_000), 7);
    }

    #[test]
    fn product_commas_are_percent_encoded_in_gff3() {
        let mut f = cds("a", 1, 30);
        f.func.product = Some("Protein A, alpha; beta=1".to_string());
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        // Reserved chars inside the value must be encoded, not raw.
        assert!(out.contains("%2C"), "comma must encode to %2C: {out}");
        assert!(out.contains("%3B"), "semicolon must encode to %3B: {out}");
        assert!(out.contains("%3D"), "equals must encode to %3D: {out}");
        // And the raw product attribute must not leak an unescaped comma value.
        assert!(!out.contains("Protein A, alpha"));
    }

    fn ncrna(kind: FeatureKind, id: &str, source: &str, acc: &str, name: &str) -> Feature {
        Feature {
            kind,
            contig: "contig1".to_string(),
            id: id.to_string(),
            start: 1,
            end: 100,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: source.to_string(),
                accession: acc.to_string(),
                name: name.to_string(),
                score: 100.0,
                evalue: Some(1e-20),
                ref_len: None,
            }],
            func: Functional::default(),
        }
    }

    #[test]
    fn cds_hmm_hit_goes_to_inference_not_dbxref_and_ec_is_own_qualifier() {
        let mut f = cds("a", 1, 30);
        f.annotations.push(Annotation {
            source: "rustyhmmer:ncbifams".to_string(),
            accession: "NF006959.0".to_string(),
            name: "some enzyme".to_string(),
            score: 200.0,
            evalue: Some(1e-40),
            ref_len: None,
        });
        f.func.ec = vec!["1.1.1.3".to_string()];
        f.func.ko = Some("K00001".to_string());
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        // HMM hit → inference, NOT db_xref.
        assert!(
            out.contains("inference=protein motif:HMM:NCBIfam:NF006959.0"),
            "{out}"
        );
        assert!(!out.contains("rustyhmmer"), "no rustyhmmer db_xref: {out}");
        // EC → EC_number attribute, never a db_xref.
        assert!(out.contains("EC_number=1.1.1.3"), "{out}");
        assert!(!out.contains("Dbxref="), "no illegal db_xref emitted: {out}");
        // KEGG is not a registered db_xref → moved to a note.
        assert!(out.contains("KEGG:K00001"), "{out}");
        // transl_table=11 on every CDS.
        assert!(out.contains("transl_table=11"), "{out}");
    }

    #[test]
    fn rfam_hit_relabels_to_registered_rfam_dbxref() {
        let f = ncrna(FeatureKind::Ncrna, "n1", "infernox", "RF00506", "THI");
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("Dbxref=RFAM:RF00506"), "{out}");
        assert!(!out.contains("infernox"), "no infernox db_xref token: {out}");
    }

    #[test]
    fn rrna_products_map_to_standard_ncbi_strings() {
        assert_eq!(normalize_rrna_product("SSU_rRNA_bacteria"), "16S ribosomal RNA");
        assert_eq!(normalize_rrna_product("LSU_rRNA_bacteria"), "23S ribosomal RNA");
        assert_eq!(normalize_rrna_product("5S_rRNA"), "5S ribosomal RNA");
        assert_eq!(normalize_rrna_product("5_8S_rRNA"), "5.8S ribosomal RNA");
        assert_eq!(normalize_rrna_product("SSU_rRNA_archaea"), "16S ribosomal RNA");
        // Full-feature render path also applies the map.
        let f = ncrna(FeatureKind::Rrna, "r1", "infernox", "RF00177", "SSU_rRNA_bacteria");
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("product=16S ribosomal RNA"), "{out}");
        assert!(!out.contains("SSU_rRNA_bacteria"), "raw model name leaked: {out}");
    }

    #[test]
    fn trna_undetermined_never_blank_product() {
        assert_eq!(normalize_trna_product("tRNA-Undet"), "tRNA-Xxx");
        assert_eq!(normalize_trna_product("tRNA-"), "tRNA-Xxx");
        assert_eq!(normalize_trna_product("tRNA"), "tRNA-Xxx");
        assert_eq!(normalize_trna_product("tRNA-Pseudo"), "tRNA-Xxx");
        assert_eq!(normalize_trna_product("tRNA-Ile"), "tRNA-Ile");
    }

    #[test]
    fn regulatory_region_emits_class_and_note() {
        let mut f = ncrna(FeatureKind::RegulatoryRegion, "g1", "infernox", "RF00059", "TPP");
        f.func.regulatory_class = Some("riboswitch".to_string());
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("regulatory_class=riboswitch"), "{out}");
        // Rfam family recorded as a note (satisfies the `other`-needs-note rule).
        assert!(out.contains("Rfam:TPP"), "{out}");
        // A regulatory_region with no explicit class defaults to `other`.
        let f2 = ncrna(FeatureKind::RegulatoryRegion, "g2", "infernox", "RF99999", "Foo");
        let mut buf2 = Vec::new();
        write_gff3(&mut buf2, std::slice::from_ref(&f2)).unwrap();
        let out2 = String::from_utf8(buf2).unwrap();
        assert!(out2.contains("regulatory_class=other"), "{out2}");
    }

    #[test]
    fn hypothetical_cds_suppresses_gene_symbol() {
        let mut f = cds("a", 1, 30);
        f.func.product = Some("hypothetical protein".to_string());
        f.func.gene = Some("ybxG".to_string());
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGT".to_vec())];
        let mut gff = Vec::new();
        write_gff3_submission(&mut gff, std::slice::from_ref(&f), &contigs, "BACT").unwrap();
        let gff = String::from_utf8(gff).unwrap();
        assert!(gff.contains("hypothetical protein"), "{gff}");
        assert!(!gff.contains("ybxG"), "gene symbol must be suppressed: {gff}");
        let mut gbff = Vec::new();
        write_gbff(&mut gbff, std::slice::from_ref(&f), &contigs).unwrap();
        let gbff = String::from_utf8(gbff).unwrap();
        assert!(!gbff.contains("/gene=\"ybxG\""), "{gbff}");
    }

    #[test]
    fn gbff_translation_strips_trailing_stop() {
        let mut f = cds("a", 1, 30);
        f.aa = Some("MKVTCRE*".to_string());
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGT".to_vec())];
        let mut buf = Vec::new();
        write_gbff(&mut buf, std::slice::from_ref(&f), &contigs).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("/translation=\"MKVTCRE\""), "{out}");
        assert!(!out.contains("MKVTCRE*"), "trailing stop must be stripped: {out}");
        assert!(out.contains("/transl_table=11"), "{out}");
    }

    #[test]
    fn sanitize_product_cases() {
        assert_eq!(sanitize_product("foo."), "foo");
        assert_eq!(sanitize_product("  foo bar  "), "foo bar");
        assert_eq!(sanitize_product("foo   bar"), "foo bar");
        assert_eq!(sanitize_product("hypothetical protein"), "hypothetical protein");
        assert_eq!(sanitize_product("thing.."), "thing");
        assert_eq!(sanitize_product("DNA polymerase III."), "DNA polymerase III");
    }

    #[test]
    fn pseudogene_renders_pseudo_and_no_translation_gbff() {
        let mut f = cds("a", 1, 30);
        f.func.pseudogene = true;
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGT".to_vec())];
        let mut buf = Vec::new();
        write_gbff(&mut buf, std::slice::from_ref(&f), &contigs).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("/pseudo"), "expected /pseudo: {out}");
        assert!(
            !out.contains("/translation"),
            "pseudogene must not emit /translation: {out}"
        );
    }

    #[test]
    fn pseudogene_renders_pseudo_and_no_translation_gff3() {
        let mut f = cds("a", 1, 30);
        f.func.pseudogene = true;
        let mut buf = Vec::new();
        write_gff3(&mut buf, std::slice::from_ref(&f)).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("pseudo=true"), "{out}");
    }

    #[test]
    fn gene_parent_hierarchy_present_gff3() {
        let f = cds("a", 5, 40);
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGTACGTACGTACGT".to_vec())];
        let mut buf = Vec::new();
        write_gff3_submission(&mut buf, std::slice::from_ref(&f), &contigs, "BACT").unwrap();
        let out = String::from_utf8(buf).unwrap();
        // Header + region + FASTA directives.
        assert!(out.starts_with("##gff-version 3\n"), "{out}");
        assert!(out.contains("##sequence-region contig1 1 24"), "{out}");
        assert!(out.contains("##FASTA"), "{out}");
        // gene record + child Parent linkage.
        assert!(out.contains("\tgene\t"), "expected a gene line: {out}");
        assert!(out.contains("ID=gene-BACT_00001"), "{out}");
        assert!(out.contains("Parent=gene-BACT_00001"), "{out}");
        assert!(out.contains("locus_tag=BACT_00001"), "{out}");
    }

    #[test]
    fn gbff_gene_block_precedes_feature() {
        let f = cds("a", 1, 30);
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGT".to_vec())];
        let mut buf = Vec::new();
        write_gbff(&mut buf, std::slice::from_ref(&f), &contigs).unwrap();
        let out = String::from_utf8(buf).unwrap();
        let gene_at = out.find("     gene ").expect("gene block");
        let cds_at = out.find("     CDS ").expect("CDS block");
        assert!(gene_at < cds_at, "gene must precede CDS: {out}");
        assert!(out.contains("/locus_tag=\"BACT_00001\""), "{out}");
    }

    /// Build a bare non-CDS structural feature of the given kind.
    fn feat(kind: FeatureKind, start: i64, end: i64) -> Feature {
        Feature {
            kind,
            contig: "contig1".to_string(),
            id: "f1".to_string(),
            start,
            end,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        }
    }

    #[test]
    fn assembly_gap_no_locus_tag_and_gap_quals() {
        let f = feat(FeatureKind::AssemblyGap, 100, 149); // len 50
        let contigs = vec![("contig1".to_string(), vec![b'N'; 200])];
        // GFF3
        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(!g.contains("locus_tag"), "assembly_gap must not carry locus_tag: {g}");
        assert!(g.contains("estimated_length=50"), "{g}");
        assert!(g.contains("gap_type=within"), "{g}");
        assert!(g.contains("linkage_evidence=paired-ends"), "{g}");
        // GBFF
        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(!b.contains("/locus_tag"), "{b}");
        assert!(b.contains("/estimated_length=50"), "{b}");
        assert!(b.contains("/gap_type=\"within scaffold\""), "{b}");
        assert!(b.contains("/linkage_evidence=\"paired-ends\""), "{b}");
    }

    #[test]
    fn mobile_element_type_emitted() {
        let is = feat(FeatureKind::IsElement, 1, 100);
        let integ = feat(FeatureKind::Integron, 200, 400);
        let contigs = vec![("contig1".to_string(), vec![b'A'; 500])];
        let mut b = Vec::new();
        write_gbff(&mut b, &[is.clone(), integ.clone()], &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("/mobile_element_type=\"insertion sequence\""), "{b}");
        assert!(b.contains("/mobile_element_type=\"integron\""), "{b}");
        let mut g = Vec::new();
        write_gff3(&mut g, &[is, integ]).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("mobile_element_type=insertion"), "{g}");
        assert!(g.contains("mobile_element_type=integron"), "{g}");
    }

    #[test]
    fn repeat_region_rpt_type() {
        let crispr = feat(FeatureKind::Crispr, 1, 100);
        let tandem = feat(FeatureKind::TandemRepeat, 200, 250);
        let contigs = vec![("contig1".to_string(), vec![b'A'; 300])];
        let mut b = Vec::new();
        write_gbff(&mut b, &[crispr.clone(), tandem.clone()], &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("/rpt_type=direct"), "{b}");
        assert!(b.contains("/rpt_family=\"CRISPR\""), "{b}");
        assert!(b.contains("/rpt_type=tandem"), "{b}");
        let mut g = Vec::new();
        write_gff3(&mut g, &[crispr, tandem]).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("rpt_type=direct"), "{g}");
        assert!(g.contains("rpt_type=tandem"), "{g}");
    }

    #[test]
    fn gbff_source_has_organism_and_mol_type() {
        let contigs = vec![("contig1".to_string(), b"ACGT".to_vec())];
        let mut buf = Vec::new();
        write_gbff(&mut buf, &[], &contigs).unwrap();
        let out = String::from_utf8(buf).unwrap();
        assert!(out.contains("/organism="), "{out}");
        assert!(out.contains("/mol_type=\"genomic DNA\""), "{out}");
    }

    #[test]
    fn cds_has_codon_start() {
        let f = cds("a", 1, 30);
        let contigs = vec![("contig1".to_string(), b"ACGTACGTACGT".to_vec())];
        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("/codon_start=1"), "{b}");
        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("codon_start=1"), "{g}");
    }

    #[test]
    fn locus_name_truncated_to_16() {
        let long = "supercalifragilistic_contig_name"; // 32 chars
        let contigs = vec![(long.to_string(), b"ACGT".to_vec())];
        let mut buf = Vec::new();
        write_gbff(&mut buf, &[], &contigs).unwrap();
        let out = String::from_utf8(buf).unwrap();
        let locus_line = out.lines().find(|l| l.starts_with("LOCUS")).unwrap();
        assert!(locus_line.contains(&long[..16]), "{locus_line}");
        assert!(!locus_line.contains(long), "full name must be truncated: {locus_line}");
    }

    #[test]
    fn wrapped_qualifier_preserves_multibyte_chars() {
        // 100 multi-byte chars (2-byte each) → forces wrapping; a byte-index
        // slice would split a char and drop the chunk. Reconstruct and compare.
        let value: String = "αβγ×é".repeat(20);
        let mut buf = Vec::new();
        write_wrapped_qualifier(&mut buf, "note", &value).unwrap();
        let out = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = out.lines().collect();
        let last = lines.len() - 1;
        let mut recon = String::new();
        for (idx, line) in lines.iter().enumerate() {
            let mut s = line.trim_start();
            if idx == 0 {
                s = s.strip_prefix("/note=\"").expect("head prefix");
            }
            if idx == last {
                s = s.strip_suffix('"').unwrap_or(s);
            }
            recon.push_str(s);
        }
        assert_eq!(recon, value, "no chars dropped or split: {out}");
        assert!(lines.len() > 1, "expected wrapping: {out}");
    }

    // --- Draft-genome partiality (`<`/`>`) emission ------------------------

    /// A complete CDS (both ends present) must carry NO partiality markers in
    /// either GBFF or GFF3 — the existing byte layout is untouched.
    #[test]
    fn complete_cds_has_no_partial_markers() {
        let f = cds("f1", 10, 99); // strand +1, both ends complete
        let contigs = vec![("contig1".to_string(), vec![b'A'; 200])];
        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("     CDS             10..99"), "{b}");
        assert!(!b.contains('<') && !b.contains('>'), "no markers: {b}");

        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(!g.contains("partial=true"), "{g}");
        assert!(!g.contains("start_range") && !g.contains("end_range"), "{g}");
    }

    /// A 3'-partial + strand CDS (no stop, runs off the high coordinate) →
    /// `start..>end` in GBFF and `partial=true;end_range=<end>,.` in GFF3.
    #[test]
    fn partial3_plus_strand_marks_high_coord() {
        let mut f = cds("f1", 10, 99);
        f.partial3 = true; // 3' end incomplete; + strand → high coordinate
        let contigs = vec![("contig1".to_string(), vec![b'A'; 200])];

        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("     CDS             10..>99"), "{b}");
        assert!(!b.contains("<10"), "5' end must stay complete: {b}");

        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("partial=true"), "{g}");
        assert!(g.contains("end_range=99,."), "{g}");
        assert!(!g.contains("start_range"), "{g}");
    }

    /// A 5'-partial + strand CDS (no start, runs off the low coordinate) →
    /// `<start..end` in GBFF and `partial=true;start_range=.,<start>` in GFF3.
    #[test]
    fn partial5_plus_strand_marks_low_coord() {
        let mut f = cds("f1", 10, 99);
        f.partial5 = true; // 5' end incomplete; + strand → low coordinate
        let contigs = vec![("contig1".to_string(), vec![b'A'; 200])];

        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("     CDS             <10..99"), "{b}");

        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("partial=true"), "{g}");
        assert!(g.contains("start_range=.,10"), "{g}");
    }

    /// On the minus strand the 5'/3' ends swap relative to the coordinates: a
    /// 5'-partial − strand CDS is incomplete at its HIGH coordinate, so the GBFF
    /// `complement(...)` location marks the `>` end and GFF3 uses `end_range`.
    #[test]
    fn partial5_minus_strand_marks_high_coord() {
        let mut f = cds("f1", 10, 99);
        f.strand = -1;
        f.partial5 = true; // 5' end incomplete; − strand → high coordinate
        let contigs = vec![("contig1".to_string(), vec![b'A'; 200])];

        let mut b = Vec::new();
        write_gbff(&mut b, std::slice::from_ref(&f), &contigs).unwrap();
        let b = String::from_utf8(b).unwrap();
        assert!(b.contains("complement(10..>99)"), "{b}");

        let mut g = Vec::new();
        write_gff3(&mut g, std::slice::from_ref(&f)).unwrap();
        let g = String::from_utf8(g).unwrap();
        assert!(g.contains("partial=true"), "{g}");
        assert!(g.contains("end_range=99,."), "{g}");
        assert!(!g.contains("start_range"), "{g}");
    }
}
