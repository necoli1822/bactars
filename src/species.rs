//! 16S rRNA species identification against the GTDB SSU reference set.
//!
//! Given the 16S rRNA feature(s) already detected on the genome, we extract each
//! one's nucleotide sequence, mmseqs-search it (`--search-type 3`) against the
//! GTDB bacterial + archaeal SSU representative sequences, take the single best
//! hit, and read the GTDB taxonomy string off that reference's header. The rank
//! we report scales with identity: >=98.7% -> species; >=94.5% -> genus; else
//! family (standard SSU identity cut-offs). Genome-level result.
//!
//! Pipeline (mirrors psc.rs / vfdb.rs): gunzip the two GTDB `.fna.gz` files and
//! concatenate them to a cached refs FASTA + a cached `id -> taxonomy` TSV ONCE,
//! then `mmseqs createdb` the refs (cached). Per run: write the 16S query
//! sequences -> `createdb` -> `search --search-type 3` -> `convertalis` -> best
//! hit -> streaming taxonomy lookup. Any mmseqs / IO failure logs to stderr and
//! returns Ok(None) so the pipeline never crashes.

use crate::fasta::Contig;
use crate::feature::{Feature, FeatureKind};
use std::collections::HashMap;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// SSU identity (percent) at/above which a hit resolves to species rank.
const SPECIES_PIDENT: f64 = 98.7;
/// ...genus rank.
const GENUS_PIDENT: f64 = 94.5;

/// A genome-level species assignment from the best 16S vs GTDB match.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesCall {
    /// Full GTDB lineage string, e.g. `d__Bacteria;p__...;s__Escherichia coli`.
    pub taxonomy: String,
    /// The resolved assignment at the rank supported by the identity (species /
    /// genus / family), e.g. `Escherichia coli` or `Escherichia (genus)`.
    pub species: String,
    /// Best-hit percent identity (0..100).
    pub identity: f64,
    /// The GTDB reference id (accession) of the matched 16S sequence.
    pub matched_16s: String,
}

/// Rfam accessions for the small-subunit (16S/18S) rRNA families across domains.
const SSU_RFAM_ACCS: &[&str] = &["RF00177", "RF01959", "RF01960", "RF02542"];

/// True if a feature is a small-subunit (16S) rRNA. Matches on the product /
/// annotation name / annotation accession, since the producing tier may label it
/// "16S ribosomal RNA" (curated), "SSU_rRNA_bacteria" (Rfam model name), or only
/// carry the Rfam accession `RF00177`.
fn is_16s(f: &Feature) -> bool {
    if f.kind != FeatureKind::Rrna {
        return false;
    }
    // Assemble every label we might recognise, lowercased.
    let mut hay = f.display_product().unwrap_or_default().to_ascii_lowercase();
    for a in &f.annotations {
        hay.push(' ');
        hay.push_str(&a.name.to_ascii_lowercase());
    }
    if hay.contains("16s") || hay.contains("ssu") {
        return true;
    }
    // Fall back to the Rfam accession (strip any `.N` version suffix).
    f.annotations.iter().any(|a| {
        let acc = a.accession.split('.').next().unwrap_or(&a.accession);
        SSU_RFAM_ACCS.iter().any(|s| acc.eq_ignore_ascii_case(s))
    })
}

/// Reverse-complement a nucleotide byte slice (IUPAC-agnostic: only A/C/G/T
/// complemented, everything else passes through as N).
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            other => other,
        })
        .collect()
}

/// Extract the nucleotide sequence spanned by a feature (1-based inclusive,
/// reverse-complemented on the minus strand). Returns None if the contig is
/// missing or the coordinates are out of range.
fn feature_seq(feat: &Feature, by_contig: &HashMap<&str, &[u8]>) -> Option<String> {
    let seq = by_contig.get(feat.contig.as_str())?;
    if feat.start < 1 || feat.end < feat.start || feat.end as usize > seq.len() {
        return None;
    }
    let sub = &seq[(feat.start as usize - 1)..(feat.end as usize)];
    let bytes = if feat.strand < 0 { revcomp(sub) } else { sub.to_vec() };
    Some(String::from_utf8_lossy(&bytes).into_owned())
}

/// Parse a GTDB SSU header into (reference id, taxonomy lineage). The id is the
/// first whitespace token; the lineage runs from `d__` to the first ` [` metadata
/// bracket. Returns None if no lineage is present.
pub fn parse_gtdb_header(header: &str) -> Option<(String, String)> {
    let h = header.trim_start_matches('>').trim();
    let id = h.split_whitespace().next()?.to_string();
    let d = h.find("d__")?;
    let tax = &h[d..];
    let tax = match tax.find(" [") {
        Some(i) => &tax[..i],
        None => tax,
    };
    Some((id, tax.trim().to_string()))
}

/// Read a rank value (e.g. `s__`, `g__`, `f__`) out of a GTDB lineage string.
/// Returns the trimmed value without its `x__` prefix, or "" if absent/empty.
pub fn rank_value(taxonomy: &str, prefix: &str) -> String {
    for field in taxonomy.split(';') {
        let field = field.trim();
        if let Some(v) = field.strip_prefix(prefix) {
            return v.trim().to_string();
        }
    }
    String::new()
}

/// Resolve the reportable assignment string from a lineage + best-hit identity.
pub fn assign_rank(taxonomy: &str, identity: f64) -> String {
    if identity >= SPECIES_PIDENT {
        let s = rank_value(taxonomy, "s__");
        if !s.is_empty() {
            return s;
        }
    }
    if identity >= GENUS_PIDENT {
        let g = rank_value(taxonomy, "g__");
        if !g.is_empty() {
            return format!("{g} (genus)");
        }
    }
    let f = rank_value(taxonomy, "f__");
    if !f.is_empty() {
        return format!("{f} (family)");
    }
    // Fall back to the deepest named rank we can find.
    for prefix in ["g__", "o__", "c__", "p__", "d__"] {
        let v = rank_value(taxonomy, prefix);
        if !v.is_empty() {
            return format!("{v} (higher rank)");
        }
    }
    "unclassified".to_string()
}

/// One parsed convertalis row.
#[derive(Clone, Debug)]
struct M8Hit {
    target: String,
    pident: f64,
    bits: f32,
}

/// Parse `query,target,pident,evalue,bits` and return the single overall best hit
/// (highest bits) across all query 16S sequences.
fn parse_m8_best(path: &Path) -> Option<M8Hit> {
    let file = fs::File::open(path).ok()?;
    let mut best: Option<M8Hit> = None;
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() < 5 {
            continue;
        }
        let pident_raw: f64 = p[2].parse().unwrap_or(0.0);
        // mmseqs may report pident as a fraction (0..1); normalise to a percent so
        // the SSU identity cut-offs (percent) are consistent across all stages (F8).
        let pident = if pident_raw <= 1.0 { pident_raw * 100.0 } else { pident_raw };
        let hit = M8Hit {
            target: p[1].to_string(),
            pident,
            bits: p[4].parse().unwrap_or(0.0),
        };
        match &best {
            Some(b) if b.bits >= hit.bits => {}
            _ => best = Some(hit),
        }
    }
    best
}

/// gunzip a `.fna.gz` and append its records to `refs`, and its `id<TAB>taxonomy`
/// header rows to `tax`. Missing files are skipped (bac may exist without ar).
fn ingest_gz(gz: &Path, refs: &mut String, tax: &mut String) {
    if !gz.exists() {
        return;
    }
    // Decompress in-process (pure Rust): read the .gz then inflate to bytes.
    let raw = match fs::read(gz) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("[species] cannot read {}: {e}", gz.display());
            return;
        }
    };
    let out = match crate::util_io::gunzip_bytes(&raw) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("[species] gunzip failed for {}: {e}", gz.display());
            return;
        }
    };
    for line in String::from_utf8_lossy(&out).lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some((id, lineage)) = parse_gtdb_header(rest) {
                tax.push_str(&id);
                tax.push('\t');
                tax.push_str(&lineage);
                tax.push('\n');
            }
            // Re-emit a clean id-only header so the mmseqs target id = the id we
            // keyed the taxonomy TSV on.
            let id = line[1..].split_whitespace().next().unwrap_or("");
            refs.push('>');
            refs.push_str(id);
            refs.push('\n');
        } else {
            refs.push_str(line);
            refs.push('\n');
        }
    }
}

/// Ensure the GTDB refs are gunzipped/concatenated to `<dir>/gtdb_ssu_refs.fna`,
/// the taxonomy TSV built at `<dir>/gtdb_taxonomy.tsv`, and the mmseqs target DB
/// at `<dir>/gtdb_ssu_db` (all cached). Returns (refs, taxonomy_tsv, mmseqs_db).
fn ensure_db(
    bin: &str,
    species_dir: &Path,
) -> Result<(std::path::PathBuf, std::path::PathBuf, std::path::PathBuf), String> {
    let refs = species_dir.join("gtdb_ssu_refs.fna");
    let tax_tsv = species_dir.join("gtdb_taxonomy.tsv");
    let db = species_dir.join("gtdb_ssu_db");

    if !refs.exists() || !tax_tsv.exists() {
        let mut refs_txt = String::new();
        let mut tax_txt = String::new();
        ingest_gz(&species_dir.join("bac120_ssu_reps.fna.gz"), &mut refs_txt, &mut tax_txt);
        ingest_gz(&species_dir.join("ar53_ssu_reps.fna.gz"), &mut refs_txt, &mut tax_txt);
        if refs_txt.is_empty() {
            return Err("no GTDB SSU references (bac120/ar53 .fna.gz) found".to_string());
        }
        fs::write(&refs, refs_txt).map_err(|e| format!("write gtdb_ssu_refs.fna: {e}"))?;
        fs::write(&tax_tsv, tax_txt).map_err(|e| format!("write gtdb_taxonomy.tsv: {e}"))?;
    }
    if !db.with_extension("dbtype").exists() {
        crate::mmseqs::run_with(bin, &["createdb", &refs.to_string_lossy(), &db.to_string_lossy(), "--dbtype", "2", "-v", "1"])?;
    }
    Ok((refs, tax_tsv, db))
}

/// Look up a single reference id's taxonomy in the cached TSV (streaming).
fn lookup_taxonomy(tax_tsv: &Path, wanted: &str) -> Option<String> {
    let file = fs::File::open(tax_tsv).ok()?;
    for line in BufReader::new(file).lines().map_while(Result::ok) {
        if let Some((id, lineage)) = line.split_once('\t') {
            if id == wanted {
                return Some(lineage.to_string());
            }
        }
    }
    None
}

/// Identify the genome's species from its 16S rRNA feature(s) vs GTDB SSU.
/// `species_dir` holds `bac120_ssu_reps.fna.gz` (+ optional `ar53_ssu_reps.fna.gz`).
/// Returns Ok(None) when there is no 16S feature or no qualifying hit; degrades
/// gracefully (logs + Ok(None)) on any mmseqs / IO error.
pub fn identify(
    features: &[Feature],
    contigs: &[Contig],
    species_dir: &str,
    threads: usize,
) -> Result<Option<SpeciesCall>, String> {
    // 1. Collect 16S query sequences from the already-detected rRNA features.
    let by_contig: HashMap<&str, &[u8]> =
        contigs.iter().map(|c| (c.name.as_str(), c.seq.as_slice())).collect();
    let mut query = String::new();
    let mut n = 0usize;
    for (i, f) in features.iter().enumerate() {
        if !is_16s(f) {
            continue;
        }
        if let Some(seq) = feature_seq(f, &by_contig) {
            if seq.trim().is_empty() {
                continue;
            }
            query.push('>');
            query.push_str(&format!("q16s_{i}"));
            query.push('\n');
            query.push_str(&seq);
            query.push('\n');
            n += 1;
        }
    }
    if n == 0 {
        return Ok(None); // no 16S to classify with
    }

    let species_dir = Path::new(species_dir);
    if !crate::mmseqs::available() {
        eprintln!("[species] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return Ok(None);
    }
    let bin = crate::mmseqs::bin();
    let (_refs, tax_tsv, db) = match ensure_db(&bin, species_dir) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("[species] cannot prepare GTDB SSU DB under {}: {e} — skipping", species_dir.display());
            return Ok(None);
        }
    };

    let work = std::env::temp_dir().join(format!("bactars_species_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[species] cannot create temp dir {}: {e}", work.display());
        return Ok(None);
    }
    let result = run_search(&bin, &work, &db, &query, threads);
    let _ = fs::remove_dir_all(&work);

    let best = match result {
        Ok(Some(b)) => b,
        Ok(None) => return Ok(None),
        Err(e) => {
            eprintln!("[species] mmseqs pipeline failed: {e} — continuing without species call");
            return Ok(None);
        }
    };

    let Some(taxonomy) = lookup_taxonomy(&tax_tsv, &best.target) else {
        eprintln!("[species] best-hit id {} absent from taxonomy TSV — skipping", best.target);
        return Ok(None);
    };
    let species = assign_rank(&taxonomy, best.pident);
    Ok(Some(SpeciesCall {
        taxonomy,
        species,
        identity: best.pident,
        matched_16s: best.target,
    }))
}

/// createdb(query, refs nucleotide) -> search --search-type 3 -> convertalis ->
/// overall best hit.
fn run_search(
    bin: &str,
    work: &Path,
    db: &Path,
    query: &str,
    threads: usize,
) -> Result<Option<M8Hit>, String> {
    let query_fna = work.join("query.fna");
    let query_db = work.join("queryDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");
    let nt = threads.clamp(1, 8).to_string();

    fs::write(&query_fna, query).map_err(|e| format!("write query.fna: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_fna, q_db, res_db, tmp_s, db_s, m8_s) =
        (s(&query_fna), s(&query_db), s(&result_db), s(&tmp), s(db), s(&m8));

    crate::mmseqs::run_with(bin, &["createdb", &q_fna, &q_db, "--dbtype", "2", "-v", "1"])?;
    crate::mmseqs::run_with(
        bin,
        &[
            "search", &q_db, &db_s, &res_db, &tmp_s, "--threads", &nt, "--search-type", "3",
            "-e", "1e-20", "--max-seqs", "50", "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        bin,
        &[
            "convertalis", &q_db, &db_s, &res_db, &m8_s, "--format-output",
            "query,target,pident,evalue,bits", "--threads", &nt, "-v", "1",
        ],
    )?;
    Ok(parse_m8_best(&m8))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Annotation, Functional};
    use std::io::Write as _;

    #[test]
    fn parse_gtdb_header_extracts_id_and_lineage() {
        let h = ">RS_GCF_031457235.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [locus_tag=NZ_X] [location=5..1531]";
        let (id, tax) = parse_gtdb_header(h).expect("parses");
        assert_eq!(id, "RS_GCF_031457235.1");
        assert_eq!(
            tax,
            "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli"
        );
    }

    #[test]
    fn rank_value_reads_fields() {
        let tax = "d__Bacteria;p__X;g__Escherichia;s__Escherichia coli";
        assert_eq!(rank_value(tax, "s__"), "Escherichia coli");
        assert_eq!(rank_value(tax, "g__"), "Escherichia");
        assert_eq!(rank_value(tax, "f__"), "");
    }

    #[test]
    fn assign_rank_scales_with_identity() {
        let tax = "d__Bacteria;p__X;o__Y;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli";
        assert_eq!(assign_rank(tax, 99.5), "Escherichia coli");
        assert_eq!(assign_rank(tax, 96.0), "Escherichia (genus)");
        assert_eq!(assign_rank(tax, 90.0), "Enterobacteriaceae (family)");
    }

    #[test]
    fn is_16s_matches_product_or_annotation() {
        let mut f = Feature {
            kind: FeatureKind::Rrna,
            contig: "c1".to_string(),
            id: "r1".to_string(),
            start: 1,
            end: 1500,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional { product: Some("16S ribosomal RNA".to_string()), ..Default::default() },
        };
        assert!(is_16s(&f));
        f.func.product = Some("23S ribosomal RNA".to_string());
        assert!(!is_16s(&f));
        // Rfam model name "SSU_rRNA_bacteria" (no literal "16s") is recognised.
        f.func.product = None;
        f.annotations.push(Annotation {
            source: "infernox".to_string(),
            accession: "RF00177".to_string(),
            name: "SSU_rRNA_bacteria".to_string(),
            score: 100.0,
            evalue: Some(0.0),
            ref_len: None,
        });
        assert!(is_16s(&f));
        // ...and by Rfam accession alone (name deliberately uninformative).
        let mut g = f.clone();
        g.annotations[0].name = "x".to_string();
        g.func.product = None;
        assert!(is_16s(&g));
        // A large-subunit rRNA (RF02541) is NOT 16S.
        let mut lsu = f.clone();
        lsu.func.product = None;
        lsu.annotations[0].name = "LSU_rRNA_bacteria".to_string();
        lsu.annotations[0].accession = "RF02541".to_string();
        assert!(!is_16s(&lsu));
    }

    #[test]
    fn feature_seq_extracts_and_revcomps() {
        let seq = b"AAAACGTTTT".to_vec();
        let mut by = HashMap::new();
        by.insert("c1", seq.as_slice());
        let plus = Feature {
            kind: FeatureKind::Rrna,
            contig: "c1".to_string(),
            id: "r".to_string(),
            start: 5,
            end: 6,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional::default(),
        };
        // positions 5..6 (1-based) = "CG".
        assert_eq!(feature_seq(&plus, &by).as_deref(), Some("CG"));
        let minus = Feature { strand: -1, ..plus.clone() };
        // revcomp("CG") = "CG" -> complement reverse: C->G, G->C reversed = "CG".
        assert_eq!(feature_seq(&minus, &by).as_deref(), Some("CG"));
    }

    #[test]
    fn m8_best_hit_picks_highest_bits() {
        let dir = std::env::temp_dir().join(format!("bactars_species_m8_{}", std::process::id()));
        let _ = fs::create_dir_all(&dir);
        let m8 = dir.join("result.m8");
        let mut fh = fs::File::create(&m8).unwrap();
        writeln!(fh, "q16s_1\tRS_A\t99.1\t0.0\t1400").unwrap();
        writeln!(fh, "q16s_1\tRS_B\t97.0\t0.0\t1300").unwrap();
        writeln!(fh, "q16s_2\tRS_C\t95.0\t0.0\t1200").unwrap();
        drop(fh);
        let best = parse_m8_best(&m8).unwrap();
        assert_eq!(best.target, "RS_A");
        assert!((best.pident - 99.1).abs() < 1e-9);
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AAAC"), b"GTTT");
        assert_eq!(revcomp(b"NNAT"), b"ATNN");
    }
}
