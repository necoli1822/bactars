//! AMR point-mutation (variant) detection: for known resistance loci (gyrA S83L,
//! rpoB, parC, ...) align each candidate CDS protein to the specific
//! AMRFinderPlus reference protein and inspect the residue at each catalogued
//! resistance position, per `AMRProt-mutation.tsv`. When the query carries the
//! resistance (alt) allele we annotate the affected CDS in place with a `/note`
//! + a structured `/inference` — this does NOT create a new feature and it does
//! NOT overwrite the CDS product.
//!
//! Pipeline (mirrors psc.rs): parse the mutation catalog -> build a subset
//! reference FASTA of only the mutation-bearing AMR proteins -> `mmseqs createdb`
//! (query CDS + refs) -> `search` -> `convertalis` emitting `qaln,taln,tstart`
//! -> walk each alignment to map catalogued reference positions to the aligned
//! query residue -> compare against the wild-type / resistance alleles.
//!
//! Graceful: any mmseqs / IO / parse failure logs to stderr and returns Ok(())
//! so the surrounding pipeline continues.

use crate::feature::{Feature, FeatureKind, Inference};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

/// Minimum percent identity for a CDS<->AMR-reference alignment to be trusted for
/// residue-level mutation calling. Point-mutation detection requires we aligned
/// the *right* gene; a distant homolog would mis-map positions. AMRFinderPlus
/// uses a comparably strict identity floor for its mutation model.
const MIN_PIDENT: f64 = 90.0;

/// The catalogued effect of a mutation at a reference position. Substitutions,
/// deletions and nonsense (premature stop) are detectable at the protein-alignment
/// layer; insertions / frameshifts / duplications are not (see [`parse_mutation_spec`]).
#[derive(Clone, Debug, PartialEq)]
enum MutationEffect {
    /// Single-residue substitution to `alt`.
    Substitution(char),
    /// Deletion of `usize` reference residue(s) starting at the mutation position.
    Deletion(usize),
    /// Premature stop codon (nonsense, `Ter` / `*`) at the mutation position.
    Nonsense,
}

/// One catalogued mutation on a reference protein.
#[derive(Clone, Debug, PartialEq)]
struct PointMutation {
    /// Gene symbol from the mutation symbol (e.g. "gyrA").
    gene: String,
    /// 1-based residue position in the reference protein (start of the event).
    pos: usize,
    /// First wild-type (susceptible) reference residue at `pos`.
    ref_aa: char,
    /// What the mutation does (substitution / deletion / nonsense).
    effect: MutationEffect,
    /// Antibiotic class (col 6, e.g. "QUINOLONE").
    class: String,
    /// Full standard mutation symbol (e.g. "gyrA_S83L"); used as the inference
    /// accession.
    symbol: String,
}

/// Standard amino-acid one-letter codes we accept for ref/alt of a substitution.
fn is_std_aa(c: char) -> bool {
    matches!(
        c,
        'A' | 'R' | 'N' | 'D' | 'C' | 'Q' | 'E' | 'G' | 'H' | 'I' | 'L' | 'K' | 'M' | 'F' | 'P'
            | 'S' | 'T' | 'W' | 'Y' | 'V'
    )
}

/// Parse the `refPosEffect` tail of a standard mutation symbol into
/// `(first_ref_residue, position, effect)`.
///
/// Grammar handled (`<ref…><pos><tail>`):
///   - substitution  `S83L`      — one ref residue, tail = one std AA
///   - nonsense      `D233Ter`   — one ref residue, tail = `Ter` or `*`
///   - deletion      `WR65del`   — one or more ref residues, tail = `del`
///
/// Returns `None` for insertions (`E258EKE`), frameshifts (`S90YfsTer15ins2`),
/// duplications and any other complex event — these cannot be reliably called from
/// a protein-level alignment and are left unhandled (see module notes / F5).
fn parse_mutation_spec(spec: &str) -> Option<(char, usize, MutationEffect)> {
    let chars: Vec<char> = spec.chars().collect();
    // Leading run of reference residues (one for substitution/nonsense, one-or-more
    // for a deletion of consecutive residues).
    let mut i = 0;
    let mut ref_res: Vec<char> = Vec::new();
    while i < chars.len() && is_std_aa(chars[i]) {
        ref_res.push(chars[i]);
        i += 1;
    }
    if ref_res.is_empty() {
        return None;
    }
    // Position digits.
    let mut digits = String::new();
    while i < chars.len() && chars[i].is_ascii_digit() {
        digits.push(chars[i]);
        i += 1;
    }
    if digits.is_empty() {
        return None;
    }
    let pos: usize = digits.parse().ok()?;
    if pos == 0 {
        return None;
    }
    let tail: String = chars[i..].iter().collect();
    let first_ref = ref_res[0];

    // Substitution: exactly one ref residue and a single std-AA alt.
    if ref_res.len() == 1 && tail.chars().count() == 1 {
        let alt = tail.chars().next().unwrap();
        if is_std_aa(alt) {
            return Some((first_ref, pos, MutationEffect::Substitution(alt)));
        }
    }
    // Nonsense: single ref residue, `Ter` / `*` tail.
    if ref_res.len() == 1 && (tail == "Ter" || tail == "*") {
        return Some((first_ref, pos, MutationEffect::Nonsense));
    }
    // Deletion: `del` tail (one or more consecutive ref residues).
    if tail == "del" {
        return Some((first_ref, pos, MutationEffect::Deletion(ref_res.len())));
    }
    // Insertion / frameshift / duplication / complex — not handled.
    None
}

/// Thin wrapper: parse a *single amino-acid substitution* only, for callers/tests
/// that specifically want to reject every non-substitution event.
#[cfg(test)]
fn parse_substitution(spec: &str) -> Option<(char, usize, char)> {
    match parse_mutation_spec(spec)? {
        (r, p, MutationEffect::Substitution(a)) => Some((r, p, a)),
        _ => None,
    }
}

/// Parse `AMRProt-mutation.tsv` into: reference accession -> list of catalogued
/// mutations (substitution / deletion / nonsense). Insertions / frameshifts and
/// other complex events are skipped. Columns (tab-separated,
/// header line begins with `#`):
///   0 #taxgroup  1 accession_version  2 mutation_position
///   3 standard_mutation_symbol  4 reported_mutation_symbol
///   5 class  6 subclass  7 mutated_protein_name
fn parse_mutation_catalog(tsv: &str) -> HashMap<String, Vec<PointMutation>> {
    let mut out: HashMap<String, Vec<PointMutation>> = HashMap::new();
    for line in tsv.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }
        let accession = cols[1].trim();
        let symbol = cols[3].trim();
        let class = cols[5].trim();
        if accession.is_empty() || symbol.is_empty() {
            continue;
        }
        // Split gene_MUTSPEC on the last '_' (mutation spec never contains '_').
        let Some((gene, spec)) = symbol.rsplit_once('_') else {
            continue;
        };
        let Some((ref_aa, pos, effect)) = parse_mutation_spec(spec) else {
            continue;
        };
        out.entry(accession.to_string()).or_default().push(PointMutation {
            gene: gene.to_string(),
            pos,
            ref_aa,
            effect,
            class: class.to_string(),
            symbol: symbol.to_string(),
        });
    }
    out
}

/// The state of the query at a given reference position in a gapped alignment.
#[derive(Clone, Copy, Debug, PartialEq)]
enum QueryState {
    /// A query residue `c` is aligned to the reference position.
    Residue(char),
    /// The reference position is covered by the alignment but deleted in the query
    /// (a '-' in `qaln`).
    Deleted,
    /// The reference position falls outside the aligned span (before `tstart`, or
    /// beyond the last aligned reference column).
    Outside,
}

/// Walk a gapped alignment (mmseqs `qaln`/`taln`, `tstart` 1-based) and classify
/// the query at reference position `target_pos`.
///
/// Both aln strings share the same column count; a '-' in `taln` is an insertion
/// relative to the reference (does not advance the reference position), a '-' in
/// `qaln` is a deletion in the query (the reference position has no residue).
/// Distinguishing a *deletion* from *outside the span* is what makes deletion
/// calling possible (F5) — the old residue helper collapsed both to `None`.
fn query_state_at_ref_pos(
    qaln: &str,
    taln: &str,
    tstart: usize,
    target_pos: usize,
) -> QueryState {
    if target_pos < tstart {
        return QueryState::Outside;
    }
    let mut tpos = tstart; // 1-based reference position of the next non-gap taln char
    for (qc, tc) in qaln.chars().zip(taln.chars()) {
        if tc != '-' {
            if tpos == target_pos {
                return if qc == '-' {
                    QueryState::Deleted
                } else {
                    QueryState::Residue(qc)
                };
            }
            tpos += 1;
        }
    }
    QueryState::Outside
}

/// Convenience wrapper preserving the original semantics: the aligned query
/// residue, or `None` for a deletion / outside the span. Retained for the
/// substitution-mapping unit tests.
#[cfg(test)]
fn query_residue_at_ref_pos(
    qaln: &str,
    taln: &str,
    tstart: usize,
    target_pos: usize,
) -> Option<char> {
    match query_state_at_ref_pos(qaln, taln, tstart, target_pos) {
        QueryState::Residue(c) => Some(c),
        _ => None,
    }
}

/// Build a subset reference FASTA containing only AMRProt proteins whose accession
/// bears catalogued mutations. AMRProt.fa headers are pipe-delimited
/// (`>ACCESSION|...`); we re-emit a clean `>ACCESSION` header so the mmseqs target
/// id is exactly the accession we key mutations on. Returns the FASTA text.
fn build_ref_fasta(amrprot_fa: &str, wanted: &HashMap<String, Vec<PointMutation>>) -> String {
    let mut out = String::new();
    let mut keep = false;
    for line in amrprot_fa.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            let acc = rest.split('|').next().unwrap_or("").trim();
            keep = wanted.contains_key(acc);
            if keep {
                out.push('>');
                out.push_str(acc);
                out.push('\n');
            }
        } else if keep {
            out.push_str(line.trim());
            out.push('\n');
        }
    }
    out
}

/// A single alignment row from convertalis.
struct AlnHit {
    query: String,
    target: String,
    tstart: usize,
    qaln: String,
    taln: String,
    pident: f64,
}

/// Parse convertalis output with format
/// `query,target,qstart,qend,tstart,tend,qaln,taln,pident,evalue,bits`.
fn parse_alignments(m8: &str) -> Vec<AlnHit> {
    let mut hits = Vec::new();
    for line in m8.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 11 {
            continue;
        }
        let (Ok(tstart), Ok(pident_raw)) = (c[4].parse::<usize>(), c[8].parse::<f64>()) else {
            continue;
        };
        // mmseqs may report pident as a fraction (0..1); normalise to a percent so
        // the MIN_PIDENT (%) compare is consistent across all stages (F8).
        let pident = if pident_raw <= 1.0 { pident_raw * 100.0 } else { pident_raw };
        hits.push(AlnHit {
            query: c[0].to_string(),
            target: c[1].to_string(),
            tstart,
            qaln: c[6].to_string(),
            taln: c[7].to_string(),
            pident,
        });
    }
    hits
}

/// Run createdb -> search -> convertalis and return the parsed alignment rows.
fn run_alignment_pipeline(
    work: &Path,
    ref_fasta: &str,
    query_fasta: &str,
    threads: usize,
) -> Result<Vec<AlnHit>, String> {
    let bin = crate::mmseqs::bin();
    let query_faa = work.join("query.faa");
    let ref_faa = work.join("ref.faa");
    let query_db = work.join("queryDB");
    let ref_db = work.join("refDB");
    let result_db = work.join("resultDB");
    let tmp = work.join("tmp");
    let m8 = work.join("result.m8");
    let nt = threads.max(1).to_string();

    fs::write(&query_faa, query_fasta).map_err(|e| format!("write query.faa: {e}"))?;
    fs::write(&ref_faa, ref_fasta).map_err(|e| format!("write ref.faa: {e}"))?;
    fs::create_dir_all(&tmp).map_err(|e| format!("mkdir tmp: {e}"))?;

    let s = |p: &Path| p.to_string_lossy().to_string();
    let (q_faa, r_faa, q_db, r_db, res_db, tmp_s, m8_s) = (
        s(&query_faa),
        s(&ref_faa),
        s(&query_db),
        s(&ref_db),
        s(&result_db),
        s(&tmp),
        s(&m8),
    );

    crate::mmseqs::run_with(&bin, &["createdb", &q_faa, &q_db, "-v", "1"])?;
    crate::mmseqs::run_with(&bin, &["createdb", &r_faa, &r_db, "-v", "1"])?;
    crate::mmseqs::run_with(
        &bin,
        &[
            "search", &q_db, &r_db, &res_db, &tmp_s, "--threads", &nt, "-s", "6.0", "-e", "1e-10",
            "--max-seqs", "50", "-a", "-v", "1",
        ],
    )?;
    crate::mmseqs::run_with(
        &bin,
        &[
            "convertalis",
            &q_db,
            &r_db,
            &res_db,
            &m8_s,
            "--format-output",
            "query,target,qstart,qend,tstart,tend,qaln,taln,pident,evalue,bits",
            "--threads",
            &nt,
            "-v",
            "1",
        ],
    )?;

    let text = fs::read_to_string(&m8).map_err(|e| format!("read result.m8: {e}"))?;
    Ok(parse_alignments(&text))
}

/// Apply a catalogued mutation found in a CDS to its feature: append a `/note` and
/// a structured `/inference` (AMRFinderPlus). Does not touch the product. De-dupes
/// by mutation symbol so re-alignment against several references never
/// double-reports the same call on one CDS.
///
/// Returns `true` iff a NEW call was added (de-dupe returns `false`), so the caller
/// counts only calls that actually landed (F9).
fn apply_mutation(feat: &mut Feature, mu: &PointMutation) -> bool {
    let already = feat
        .func
        .inferences
        .iter()
        .any(|inf| inf.db == "AMRFinderPlus" && inf.accession == mu.symbol);
    if already {
        return false;
    }
    let descr = match &mu.effect {
        MutationEffect::Substitution(alt) => {
            format!("point mutation {} {}{}{}", mu.gene, mu.ref_aa, mu.pos, alt)
        }
        MutationEffect::Deletion(span) if *span == 1 => {
            format!("deletion {} {}{}del", mu.gene, mu.ref_aa, mu.pos)
        }
        MutationEffect::Deletion(span) => {
            format!("deletion {} {}{}del ({span} residues)", mu.gene, mu.ref_aa, mu.pos)
        }
        MutationEffect::Nonsense => {
            format!("nonsense {} {}{}Ter", mu.gene, mu.ref_aa, mu.pos)
        }
    };
    feat.func.note.push(format!("{} resistance: {descr}", mu.class));
    feat.func.inferences.push(Inference {
        category: None,
        kind: "protein motif".to_string(),
        same_species: false,
        db: "AMRFinderPlus".to_string(),
        accession: mu.symbol.clone(),
    });
    true
}

/// Detect AMR point-mutation (variant-type) resistance and annotate affected CDS
/// in place. `amr_dir` holds the AMRFinderPlus DB (`AMRProt.fa`,
/// `AMRProt-mutation.tsv`, ...). Graceful on every failure.
pub fn annotate(features: &mut [Feature], amr_dir: &str, threads: usize) -> Result<(), String> {
    let amr_dir = Path::new(amr_dir);
    let mut_tsv_path = amr_dir.join("AMRProt-mutation.tsv");
    let amrprot_path = amr_dir.join("AMRProt.fa");

    // 1. Parse the mutation catalog.
    let tsv = match fs::read_to_string(&mut_tsv_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!(
                "[amr-variant] cannot read {} : {e} — skipping point-mutation detection",
                mut_tsv_path.display()
            );
            return Ok(());
        }
    };
    let catalog = parse_mutation_catalog(&tsv);
    if catalog.is_empty() {
        eprintln!("[amr-variant] no single-substitution mutations parsed — skipping");
        return Ok(());
    }

    // 2. Candidate CDS: every CDS with a protein translation. mmseqs decides which
    //    actually align to a mutation-bearing AMR reference.
    let mut query_index: HashMap<String, usize> = HashMap::new();
    let mut query_fasta = String::new();
    for (i, f) in features.iter().enumerate() {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        if let Some(aa) = &f.aa {
            let aa = aa.trim();
            if aa.is_empty() {
                continue;
            }
            query_fasta.push('>');
            query_fasta.push_str(&f.id);
            query_fasta.push('\n');
            query_fasta.push_str(aa);
            query_fasta.push('\n');
            query_index.insert(f.id.clone(), i);
        }
    }
    if query_index.is_empty() {
        return Ok(());
    }
    if !crate::mmseqs::available() {
        eprintln!("[amr-variant] skipping — mmseqs not found (set $BACTARS_MMSEQS)");
        return Ok(());
    }

    // 3. Subset reference FASTA of only mutation-bearing AMR proteins.
    let amrprot = match fs::read_to_string(&amrprot_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!(
                "[amr-variant] cannot read {} : {e} — skipping",
                amrprot_path.display()
            );
            return Ok(());
        }
    };
    let ref_fasta = build_ref_fasta(&amrprot, &catalog);
    if ref_fasta.is_empty() {
        eprintln!("[amr-variant] no mutation reference proteins found in AMRProt.fa — skipping");
        return Ok(());
    }

    // 4. mmseqs align in a per-process temp workspace.
    let work: PathBuf =
        std::env::temp_dir().join(format!("bactars_amrvar_{}", std::process::id()));
    let _ = fs::remove_dir_all(&work);
    if let Err(e) = fs::create_dir_all(&work) {
        eprintln!("[amr-variant] cannot create temp dir {}: {e}", work.display());
        return Ok(());
    }
    let result = run_alignment_pipeline(&work, &ref_fasta, &query_fasta, threads);
    let _ = fs::remove_dir_all(&work);

    let hits = match result {
        Ok(h) => h,
        Err(e) => {
            eprintln!("[amr-variant] mmseqs pipeline failed: {e} — continuing without point-mutation calls");
            return Ok(());
        }
    };

    // 5. Map each catalogued reference position through every alignment; call
    //    resistance when the aligned query residue equals the alt (resistance) allele.
    let mut calls = 0usize;
    for hit in &hits {
        if hit.pident < MIN_PIDENT {
            continue;
        }
        let Some(muts) = catalog.get(&hit.target) else {
            continue;
        };
        let Some(&idx) = query_index.get(&hit.query) else {
            continue;
        };
        for mu in muts {
            let detected = match &mu.effect {
                // Substitution: the aligned query residue equals the resistance allele.
                MutationEffect::Substitution(alt) => {
                    matches!(
                        query_state_at_ref_pos(&hit.qaln, &hit.taln, hit.tstart, mu.pos),
                        QueryState::Residue(c) if c == *alt
                    )
                }
                // Deletion: every reference position in the deleted span is covered
                // by the alignment yet gapped in the query.
                MutationEffect::Deletion(span) => (0..*span).all(|k| {
                    matches!(
                        query_state_at_ref_pos(&hit.qaln, &hit.taln, hit.tstart, mu.pos + k),
                        QueryState::Deleted
                    )
                }),
                // Nonsense: a stop symbol (`*`) is aligned at the catalog position.
                // NOTE: gene callers usually truncate the ORF at a premature stop, so
                // in practice this fires only when the CDS translation carries an
                // explicit internal `*` — see F5 limitations in the module notes.
                MutationEffect::Nonsense => matches!(
                    query_state_at_ref_pos(&hit.qaln, &hit.taln, hit.tstart, mu.pos),
                    QueryState::Residue('*')
                ),
            };
            // Count only calls that actually landed (apply_mutation de-dupes) — F9.
            if detected && apply_mutation(&mut features[idx], mu) {
                calls += 1;
            }
        }
    }
    if calls > 0 {
        eprintln!("[amr-variant] detected {calls} AMR point-mutation resistance call(s)");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureKind, Functional};

    fn cds(id: &str, aa: &str) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c1".to_string(),
            id: id.to_string(),
            start: 1,
            end: 99,
            strand: 1,
            aa: Some(aa.to_string()),
            partial5: false,
            partial3: false,
            annotations: Vec::new(),
            func: Functional::default(),
        }
    }

    // Sample matching the real AMRProt-mutation.tsv column layout: a header, a
    // substitution, two more substitutions, a deletion + nonsense (now KEPT, F5),
    // and an insertion that must still be skipped.
    const SAMPLE_TSV: &str = "\
#taxgroup\taccession_version\tmutation_position\tstandard_mutation_symbol\treported_mutation_symbol\tclass\tsubclass\tmutated_protein_name
Escherichia\tWP_000019358.1\t12\tsoxS_A12S\tsoxS_A12S\tMULTIDRUG\tMULTIDRUG\tregulatory_protein_SoxS
Escherichia\tWP_001281243.1\t83\tgyrA_S83L\tgyrA_S83L\tQUINOLONE\tQUINOLONE\tquinolone_resistant_GyrA
Escherichia\tWP_001281243.1\t87\tgyrA_D87N\tgyrA_D87N\tQUINOLONE\tQUINOLONE\tquinolone_resistant_GyrA
Escherichia\tWP_XXX.1\t65\trplD_WR65del\trplD_WR65del\tMACROLIDE\tMACROLIDE\tribosomal
Escherichia\tWP_YYY.1\t233\tlpxA_D233Ter\tlpxA_D233Ter\tPOLYMYXIN\tPOLYMYXIN\tacyltransferase
Escherichia\tWP_ZZZ.1\t258\tftsI_E258EKE\tftsI_E258EKE\tBETA-LACTAM\tBETA-LACTAM\tPBP3";

    #[test]
    fn parse_catalog_keeps_substitutions_deletions_and_nonsense_skips_complex() {
        let cat = parse_mutation_catalog(SAMPLE_TSV);
        // soxS (sub) + gyrA (2 subs) + rplD (deletion) + lpxA (nonsense); the ftsI
        // insertion (WP_ZZZ.1) is skipped.
        assert_eq!(cat.len(), 4);
        let sox = &cat["WP_000019358.1"];
        assert_eq!(sox.len(), 1);
        assert_eq!(sox[0].gene, "soxS");
        assert_eq!(sox[0].pos, 12);
        assert_eq!(sox[0].ref_aa, 'A');
        assert_eq!(sox[0].effect, MutationEffect::Substitution('S'));
        assert_eq!(sox[0].class, "MULTIDRUG");
        assert_eq!(sox[0].symbol, "soxS_A12S");

        let gyra = &cat["WP_001281243.1"];
        assert_eq!(gyra.len(), 2);
        assert!(gyra
            .iter()
            .any(|m| m.pos == 83 && m.ref_aa == 'S' && m.effect == MutationEffect::Substitution('L')));
        assert!(gyra
            .iter()
            .any(|m| m.pos == 87 && m.ref_aa == 'D' && m.effect == MutationEffect::Substitution('N')));

        // Deletion + nonsense now parsed and KEPT.
        let rpld = &cat["WP_XXX.1"];
        assert_eq!(rpld.len(), 1);
        assert_eq!(rpld[0].pos, 65);
        assert_eq!(rpld[0].ref_aa, 'W');
        assert_eq!(rpld[0].effect, MutationEffect::Deletion(2));

        let lpxa = &cat["WP_YYY.1"];
        assert_eq!(lpxa.len(), 1);
        assert_eq!(lpxa[0].pos, 233);
        assert_eq!(lpxa[0].ref_aa, 'D');
        assert_eq!(lpxa[0].effect, MutationEffect::Nonsense);

        // The insertion event is still unhandled/skipped.
        assert!(!cat.contains_key("WP_ZZZ.1"));
    }

    #[test]
    fn parse_mutation_spec_classes() {
        assert_eq!(
            parse_mutation_spec("S83L"),
            Some(('S', 83, MutationEffect::Substitution('L')))
        );
        assert_eq!(
            parse_mutation_spec("D1287N"),
            Some(('D', 1287, MutationEffect::Substitution('N')))
        );
        // Deletion: single and multi-residue.
        assert_eq!(
            parse_mutation_spec("D149del"),
            Some(('D', 149, MutationEffect::Deletion(1)))
        );
        assert_eq!(
            parse_mutation_spec("WR65del"),
            Some(('W', 65, MutationEffect::Deletion(2)))
        );
        // Nonsense: `Ter` and `*`.
        assert_eq!(
            parse_mutation_spec("D233Ter"),
            Some(('D', 233, MutationEffect::Nonsense))
        );
        assert_eq!(
            parse_mutation_spec("Q42*"),
            Some(('Q', 42, MutationEffect::Nonsense))
        );
        // Unhandled complex events.
        assert_eq!(parse_mutation_spec("E258EKE"), None); // insertion
        assert_eq!(parse_mutation_spec("S90YfsTer15ins2"), None); // frameshift
        assert_eq!(parse_mutation_spec("S83"), None); // no tail
        assert_eq!(parse_mutation_spec("83L"), None); // no ref
        assert_eq!(parse_mutation_spec("X83L"), None); // non-standard ref
    }

    #[test]
    fn query_state_tri_state() {
        //   ref pos:  10 11 12 13
        //   taln:      A  B  C  D
        //   qaln:      A  -  C  D   (ref 11 deleted in query)
        let q = "A-CD";
        let t = "ABCD";
        assert_eq!(query_state_at_ref_pos(q, t, 10, 10), QueryState::Residue('A'));
        assert_eq!(query_state_at_ref_pos(q, t, 10, 11), QueryState::Deleted);
        assert_eq!(query_state_at_ref_pos(q, t, 10, 12), QueryState::Residue('C'));
        // Before/after the aligned span -> Outside (distinct from a deletion).
        assert_eq!(query_state_at_ref_pos(q, t, 10, 9), QueryState::Outside);
        assert_eq!(query_state_at_ref_pos(q, t, 10, 99), QueryState::Outside);
    }

    #[test]
    fn apply_mutation_formats_deletion_and_nonsense() {
        let mut f = cds("cds1", "MKCILF");
        let del = PointMutation {
            gene: "rplD".to_string(),
            pos: 65,
            ref_aa: 'W',
            effect: MutationEffect::Deletion(2),
            class: "MACROLIDE".to_string(),
            symbol: "rplD_WR65del".to_string(),
        };
        assert!(apply_mutation(&mut f, &del));
        assert!(!apply_mutation(&mut f, &del)); // de-dupe -> false
        assert_eq!(f.func.note.len(), 1);
        assert_eq!(
            f.func.note[0],
            "MACROLIDE resistance: deletion rplD W65del (2 residues)"
        );

        let non = PointMutation {
            gene: "lpxA".to_string(),
            pos: 233,
            ref_aa: 'D',
            effect: MutationEffect::Nonsense,
            class: "POLYMYXIN".to_string(),
            symbol: "lpxA_D233Ter".to_string(),
        };
        assert!(apply_mutation(&mut f, &non));
        assert_eq!(
            f.func.note[1],
            "POLYMYXIN resistance: nonsense lpxA D233Ter"
        );
    }

    #[test]
    fn parse_substitution_edge_cases() {
        assert_eq!(parse_substitution("S83L"), Some(('S', 83, 'L')));
        assert_eq!(parse_substitution("A12S"), Some(('A', 12, 'S')));
        assert_eq!(parse_substitution("D1287N"), Some(('D', 1287, 'N')));
        // Non-substitutions / malformed.
        assert_eq!(parse_substitution("WR65del"), None); // two ref residues + del
        assert_eq!(parse_substitution("D233Ter"), None); // nonsense
        assert_eq!(parse_substitution("E258EKE"), None); // insertion (multi alt)
        assert_eq!(parse_substitution("S83"), None); // no alt
        assert_eq!(parse_substitution("83L"), None); // no ref
        assert_eq!(parse_substitution("X83L"), None); // non-standard ref AA
        assert_eq!(parse_substitution("S83Z"), None); // non-standard alt AA
    }

    // ---- Alignment position-mapping: the bug-prone part. ----

    #[test]
    fn map_ungapped_alignment() {
        // Reference: M K C I L F  (positions 1..6), tstart=1, identical query.
        let q = "MKCILF";
        let t = "MKCILF";
        assert_eq!(query_residue_at_ref_pos(q, t, 1, 1), Some('M'));
        assert_eq!(query_residue_at_ref_pos(q, t, 3, 1), None); // signature: (qaln,taln,tstart,pos)
    }

    #[test]
    fn map_with_tstart_offset() {
        // Alignment begins at reference position 80. Query is the gyrA QRDR.
        //   ref pos:  80 81 82 83 84 85
        //   taln:      A  V  Y  S  T  I
        //   qaln:      A  V  Y  L  T  I   <- S83L present in query
        let q = "AVYLTI";
        let t = "AVYSTI";
        // reference position 83 -> query residue at that column.
        assert_eq!(query_residue_at_ref_pos(q, t, 80, 83), Some('L'));
        // reference position 82 -> 'Y' (wild-type match).
        assert_eq!(query_residue_at_ref_pos(q, t, 80, 82), Some('Y'));
        // position before the aligned span.
        assert_eq!(query_residue_at_ref_pos(q, t, 80, 79), None);
        // position after the aligned span.
        assert_eq!(query_residue_at_ref_pos(q, t, 80, 86), None);
    }

    #[test]
    fn map_with_gap_in_query_deletion() {
        //   ref pos:  10 11 12 13 14
        //   taln:      A  B  C  D  E   (all reference residues present)
        //   qaln:      A  -  C  D  E   (reference position 11 deleted in query)
        let q = "A-CDE";
        let t = "ABCDE";
        assert_eq!(query_residue_at_ref_pos(q, t, 10, 10), Some('A'));
        assert_eq!(query_residue_at_ref_pos(q, t, 10, 11), None); // deleted in query
        assert_eq!(query_residue_at_ref_pos(q, t, 10, 12), Some('C'));
        assert_eq!(query_residue_at_ref_pos(q, t, 10, 13), Some('D'));
    }

    #[test]
    fn map_with_gap_in_target_insertion() {
        // Query has an inserted residue (gap in the reference/taln). The insertion
        // must NOT advance the reference position counter.
        //   columns:   0  1  2  3  4  5
        //   qaln:      A  X  Y  Z  W  Q
        //   taln:      A  -  Y  Z  W  Q   <- 'X' inserted in query at column 1
        //   ref pos:   20    21 22 23 24
        let q = "AXYZWQ";
        let t = "A-YZWQ";
        assert_eq!(query_residue_at_ref_pos(q, t, 20, 20), Some('A'));
        // reference position 21 is the 'Y' column (insertion skipped), query='Y'.
        assert_eq!(query_residue_at_ref_pos(q, t, 20, 21), Some('Y'));
        assert_eq!(query_residue_at_ref_pos(q, t, 20, 22), Some('Z'));
        assert_eq!(query_residue_at_ref_pos(q, t, 20, 24), Some('Q'));
    }

    #[test]
    fn map_combined_gaps_both_strands() {
        //   columns:  0  1  2  3  4  5  6
        //   qaln:     M  -  A  R  N  Q  D
        //   taln:     M  K  A  -  N  Q  D
        //   ref pos:  50 51 52    53 54 55
        //     col0 ref50 q=M
        //     col1 ref51 q=- (deleted in query)
        //     col2 ref52 q=A
        //     col3 insertion in query (taln '-') -> no ref pos
        //     col4 ref53 q=N
        //     col5 ref54 q=Q
        //     col6 ref55 q=D
        let q = "M-ARNQD";
        let t = "MKA-NQD";
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 50), Some('M'));
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 51), None); // deletion
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 52), Some('A'));
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 53), Some('N'));
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 54), Some('Q'));
        assert_eq!(query_residue_at_ref_pos(q, t, 50, 55), Some('D'));
    }

    #[test]
    fn apply_mutation_sets_note_and_inference_no_dupes() {
        let mut f = cds("cds1", "MKCILF");
        let mu = PointMutation {
            gene: "gyrA".to_string(),
            pos: 83,
            ref_aa: 'S',
            effect: MutationEffect::Substitution('L'),
            class: "QUINOLONE".to_string(),
            symbol: "gyrA_S83L".to_string(),
        };
        assert!(apply_mutation(&mut f, &mu));
        assert!(!apply_mutation(&mut f, &mu)); // second call must be a no-op (de-dupe)
        assert_eq!(f.func.note.len(), 1);
        assert_eq!(
            f.func.note[0],
            "QUINOLONE resistance: point mutation gyrA S83L"
        );
        assert_eq!(f.func.inferences.len(), 1);
        let inf = &f.func.inferences[0];
        assert_eq!(inf.kind, "protein motif");
        assert_eq!(inf.db, "AMRFinderPlus");
        assert_eq!(inf.accession, "gyrA_S83L");
        assert!(!inf.same_species);
        assert!(inf.category.is_none());
        // Product must be untouched.
        assert!(f.func.product.is_none());
    }

    #[test]
    fn build_ref_fasta_subsets_and_cleans_headers() {
        let cat = parse_mutation_catalog(SAMPLE_TSV);
        let amrprot = "\
>WP_000019358.1|1|1|soxS|soxS|mutation|2|||regulatory_protein_SoxS
MKCILFAAAAA
>WP_999.9|1|1|other|other|core|1|||not_a_mutation_gene
QQQQQQ
>WP_001281243.1|1|1|gyrA|gyrA|mutation|2|||GyrA
AVYSTIAAAA";
        let out = build_ref_fasta(amrprot, &cat);
        // Clean accession-only headers, only the two mutation-bearing proteins.
        assert!(out.contains(">WP_000019358.1\n"));
        assert!(out.contains(">WP_001281243.1\n"));
        assert!(!out.contains("WP_999.9"));
        assert!(!out.contains('|'));
    }

    // ---- Full mmseqs integration (real DB). Run with:
    //      cargo test --release amr_variant -- --ignored
    // Synthetic E. coli gyrA fragment carrying the S83L resistance allele; the
    // pipeline must call QUINOLONE gyrA_S83L on it and leave a wild-type CDS clean.
    // The wild-type residue at position 83 is 'S' (verified against AMRProt.fa
    // WP_001281243.1). ----
    const GYRA_PREFIX_WT: &str = "MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGNDWNKAYKKSARVVGDVIGKYHPHGDSAVYDTIVRMAQPFSLRYMLVDGQGNFGSIDGDSAAAMRYTEIRLAKIAHELMADLEKETVDFVDNYDGTEKIPDVMPTKIPNLLVNGSSGIAVGMATNIPPHNLTEVINGCLAYIDD";

    #[test]
    #[ignore]
    fn integration_gyra_s83l_real_db() {
        let amr_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../db/amr");
        // Build the mutant (S83L) query and a wild-type negative control.
        let wt: Vec<char> = GYRA_PREFIX_WT.chars().collect();
        assert_eq!(wt[82], 'S', "position 83 of the WT fragment must be S");
        let mut mutant = wt.clone();
        mutant[82] = 'L';
        let mutant_seq: String = mutant.into_iter().collect();

        let mut feats = vec![
            cds("gyra_mut", &mutant_seq),
            cds("gyra_wt", GYRA_PREFIX_WT),
        ];
        annotate(&mut feats, amr_dir, 4).expect("annotate returns Ok");

        // Mutant CDS -> S83L call.
        let m = &feats[0];
        assert!(
            m.func.inferences.iter().any(|i| i.accession == "gyrA_S83L"),
            "expected gyrA_S83L inference on mutant, got {:?}",
            m.func.inferences
        );
        assert!(m.func.note.iter().any(|n| n.contains("gyrA S83L")));
        // Wild-type CDS -> no resistance call (negative control).
        assert!(
            feats[1].func.inferences.is_empty(),
            "wild-type CDS must have no resistance call, got {:?}",
            feats[1].func.inferences
        );
    }

    #[test]
    fn parse_alignments_reads_convertalis_columns() {
        // query target qstart qend tstart tend qaln taln pident evalue bits
        let m8 = "cds1\tWP_001281243.1\t1\t6\t80\t85\tAVYLTI\tAVYSTI\t98.5\t1e-30\t120\n";
        let hits = parse_alignments(m8);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits[0].query, "cds1");
        assert_eq!(hits[0].target, "WP_001281243.1");
        assert_eq!(hits[0].tstart, 80);
        assert_eq!(hits[0].qaln, "AVYLTI");
        assert_eq!(hits[0].taln, "AVYSTI");
        assert!((hits[0].pident - 98.5).abs() < 1e-9);
    }
}
