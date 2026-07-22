//! Metadata cross-reference (xref) lookup layer.
//!
//! Enriches already-called CDS features with curated functional metadata
//! (gene symbol / product name / EC numbers / GO terms) by looking up the
//! accession of each feature's best HMM hit in pre-downloaded curated tables.
//!
//! This is a PURE LOOKUP — no sequence search happens here. It consumes the
//! HMM hit accessions that expert subunits (rustyhmmer NCBIfams / Pfam, etc.)
//! already attached to features in [`Feature::annotations`].
//!
//! Data sources (bactars-db `meta/` dir):
//!  - `hmm_PGAP.tsv` — NCBIfams / TIGRFAM HMM accession → {product, gene, EC, GO}.
//!    Tab-separated, header starts `#ncbi_accession`. Columns used (1-based):
//!    1=ncbi_accession, 11=product_name, 12=gene_symbol, 14=ec_numbers,
//!    15=go_terms. Multi-values are comma/semicolon separated.
//!  - `pfam2go` — GO-project mapping file. Lines like
//!    `Pfam:PF00001 7tm_1 > GO:... ; GO:0004930`. Comment lines start `!`.
//!    Maps a Pfam accession (PFxxxxx, no version) → GO term(s).

use std::collections::HashMap;
use std::fs;

use crate::feature::{Feature, FeatureKind};

/// One resolved metadata record for an HMM accession.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct XrefHit {
    /// Gene symbol (e.g. `dnaA`), if known.
    pub gene: Option<String>,
    /// Product / description name, if known.
    pub product: Option<String>,
    /// EC numbers (e.g. `3.5.2.6`).
    pub ec: Vec<String>,
    /// GO terms (e.g. `GO:0008800`).
    pub go: Vec<String>,
}

impl XrefHit {
    fn is_empty(&self) -> bool {
        self.gene.is_none() && self.product.is_none() && self.ec.is_empty() && self.go.is_empty()
    }
}

/// Loaded metadata lookup tables.
///
/// `pgap` is keyed by the exact accession as it appears in `hmm_PGAP.tsv`
/// (e.g. `NF000008.1`, `TIGR00058.1`). `pgap_nover` is the same map keyed by
/// the accession with any trailing `.N` version stripped, so that a hit
/// accession whose version differs from the table's can still resolve.
/// `pfam_go` maps a bare Pfam accession (`PF00069`, no version) to its GO list.
pub struct XrefTable {
    /// accession (with version, as in file) -> record.
    pgap: HashMap<String, XrefHit>,
    /// accession (version stripped) -> record. First-writer-wins on collision.
    pgap_nover: HashMap<String, XrefHit>,
    /// bare Pfam accession (no version) -> GO terms.
    pfam_go: HashMap<String, Vec<String>>,
}

/// Strip a trailing `.N` version suffix from an accession, if present.
///
/// `NF000282.2` -> `NF000282`, `TIGR00001.1` -> `TIGR00001`,
/// `PF00069.28` -> `PF00069`. Only strips when the part after the last `.`
/// is all ASCII digits (so real accessions containing `.` in other ways are
/// left alone).
fn strip_version(acc: &str) -> &str {
    match acc.rfind('.') {
        Some(idx) => {
            let (head, tail) = acc.split_at(idx);
            let digits = &tail[1..];
            if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
                head
            } else {
                acc
            }
        }
        None => acc,
    }
}

/// Split a raw multi-value metadata field (EC or GO) into trimmed tokens.
/// Handles comma- and semicolon-separated lists, dropping empties.
fn split_multi(raw: &str) -> Vec<String> {
    raw.split(|c| c == ',' || c == ';')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn nonempty(s: &str) -> Option<String> {
    let t = s.trim();
    if t.is_empty() {
        None
    } else {
        Some(t.to_string())
    }
}

impl XrefTable {
    /// Load `hmm_PGAP.tsv` + `pfam2go` from a bactars-db `meta/` dir.
    pub fn load(meta_dir: &str) -> Result<XrefTable, String> {
        let pgap_path = format!("{}/hmm_PGAP.tsv", meta_dir.trim_end_matches('/'));
        let pfam_path = format!("{}/pfam2go", meta_dir.trim_end_matches('/'));

        let pgap_text = fs::read_to_string(&pgap_path)
            .map_err(|e| format!("failed to read {}: {}", pgap_path, e))?;
        let pfam_text = fs::read_to_string(&pfam_path)
            .map_err(|e| format!("failed to read {}: {}", pfam_path, e))?;

        let (pgap, pgap_nover) = parse_pgap(&pgap_text);
        let pfam_go = parse_pfam2go(&pfam_text);

        Ok(XrefTable {
            pgap,
            pgap_nover,
            pfam_go,
        })
    }

    /// Number of NCBIfams/TIGRFAM accessions loaded (versioned keys).
    pub fn pgap_len(&self) -> usize {
        self.pgap.len()
    }

    /// Number of distinct Pfam accessions with GO terms loaded.
    pub fn pfam_len(&self) -> usize {
        self.pfam_go.len()
    }

    /// Look up an HMM accession (NCBIfams `NFxxxxxx.x` / `TIGRxxxxx` / Pfam
    /// `PFxxxxx`). Returns any known fields, or `None` if nothing matches.
    ///
    /// Match order: exact PGAP key, then version-stripped PGAP key, then (for
    /// Pfam-style accessions) the pfam2go GO mapping. A Pfam hit that also has
    /// a PGAP record merges the pfam2go GO terms in.
    pub fn lookup(&self, accession: &str) -> Option<XrefHit> {
        let mut hit = XrefHit::default();
        let mut found = false;

        // 1) PGAP exact match, else version-stripped match.
        if let Some(rec) = self.pgap.get(accession) {
            hit = rec.clone();
            found = true;
        } else if let Some(rec) = self.pgap_nover.get(strip_version(accession)) {
            hit = rec.clone();
            found = true;
        }

        // 2) Pfam GO mapping (keyed without version). Merge GO terms in.
        let bare = strip_version(accession);
        if bare.starts_with("PF") {
            if let Some(gos) = self.pfam_go.get(bare) {
                for g in gos {
                    if !hit.go.contains(g) {
                        hit.go.push(g.clone());
                    }
                }
                found = true;
            }
        }

        if found && !hit.is_empty() {
            Some(hit)
        } else {
            None
        }
    }
}

/// Parse the `hmm_PGAP.tsv` text into (versioned-key map, version-stripped-key
/// map). The header line (starting `#`) and blank lines are skipped.
///
/// On a version-stripped-key collision the first record wins (records are read
/// top-to-bottom); the versioned map always keeps the exact key.
fn parse_pgap(text: &str) -> (HashMap<String, XrefHit>, HashMap<String, XrefHit>) {
    let mut by_ver: HashMap<String, XrefHit> = HashMap::new();
    let mut by_nover: HashMap<String, XrefHit> = HashMap::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        // Need at least through the go_terms column (index 14).
        if cols.len() < 15 {
            continue;
        }
        let acc = cols[0].trim();
        if acc.is_empty() {
            continue;
        }
        let rec = XrefHit {
            product: nonempty(cols[10]),
            gene: nonempty(cols[11]),
            ec: split_multi(cols[13]),
            go: split_multi(cols[14]),
        };
        if rec.is_empty() {
            // Still index it so lookups resolve (with no fields) but there is
            // nothing to store; skip to keep maps lean.
            continue;
        }
        let nover = strip_version(acc).to_string();
        by_nover.entry(nover).or_insert_with(|| rec.clone());
        by_ver.insert(acc.to_string(), rec);
    }

    (by_ver, by_nover)
}

/// Parse the `pfam2go` text into a map of bare Pfam accession -> GO terms.
///
/// Lines look like `Pfam:PF00001 7tm_1 > GO:description ; GO:0004930`.
/// Comment lines start `!`. The accession is the token after `Pfam:` up to the
/// first space; the GO id is the last whitespace-delimited `GO:<digits>` token
/// on the line (the descriptive `GO:...` after `>` is ignored).
fn parse_pfam2go(text: &str) -> HashMap<String, Vec<String>> {
    let mut map: HashMap<String, Vec<String>> = HashMap::new();

    for line in text.lines() {
        if line.is_empty() || line.starts_with('!') {
            continue;
        }
        // Accession: strip the "Pfam:" prefix, then take up to first space.
        let rest = match line.strip_prefix("Pfam:") {
            Some(r) => r,
            None => continue,
        };
        let acc = match rest.split_whitespace().next() {
            Some(a) if a.starts_with("PF") => a,
            _ => continue,
        };
        // GO id: the last token matching GO:<digits>.
        let go_id = line
            .split_whitespace()
            .rev()
            .find(|tok| {
                tok.strip_prefix("GO:")
                    .map(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
                    .unwrap_or(false)
            })
            .map(|t| t.to_string());
        if let Some(go) = go_id {
            let entry = map.entry(acc.to_string()).or_default();
            if !entry.contains(&go) {
                entry.push(go);
            }
        }
    }

    map
}

/// For each CDS feature, look up its HMM annotations' accessions in the table and
/// fill `func.gene` / `func.product` / `func.ec` / `func.go` where currently empty.
/// Existing values are never overwritten, and non-CDS features are left untouched.
///
/// Annotations are tried in DESCENDING score order and the FIRST one whose
/// accession resolves in the table is used. Consulting only the single top hit was
/// a bug (F6): a top-scoring accession with no curated metadata row would shadow a
/// slightly-lower hit that does resolve, dropping the enrichment entirely.
pub fn enrich(features: &mut [Feature], table: &XrefTable) {
    for feat in features.iter_mut() {
        if feat.kind != FeatureKind::Cds {
            continue;
        }
        // Annotations by descending score; take the first that both has an
        // accession and resolves in the table.
        let mut ranked: Vec<&crate::feature::Annotation> = feat
            .annotations
            .iter()
            .filter(|a| !a.accession.is_empty())
            .collect();
        ranked.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let hit = match ranked
            .iter()
            .find_map(|a| table.lookup(&a.accession))
        {
            Some(h) => h,
            None => continue,
        };

        if feat.func.gene.is_none() {
            if let Some(g) = hit.gene {
                feat.func.gene = Some(g);
            }
        }
        if feat.func.product.is_none() {
            if let Some(p) = hit.product {
                feat.func.product = Some(p);
            }
        }
        if feat.func.ec.is_empty() && !hit.ec.is_empty() {
            feat.func.ec = hit.ec;
        }
        if feat.func.go.is_empty() && !hit.go.is_empty() {
            feat.func.go = hit.go;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Annotation, Functional};

    const PGAP_SNIPPET: &str = "\
#ncbi_accession\tsource_identifier\tlabel\tsequence_cutoff\tdomain_cutoff\thmm_length\tfamily_type\tfor_structural_annotation\tfor_naming\tfor_AMRFinder\tproduct_name\tgene_symbol\tgene_synonyms\tec_numbers\tgo_terms\tpmids\ttaxonomic_range\ttaxonomic_range_name\ttaxonomic_rank_name\tn_refseq_protein_hits\tsource\tname_orig\thmm_name\tcomment
NF000459.2\t\tlabelX\t100\t100\t200\texception\tY\tY\tY\tclass B beta-lactamase\tblaSPG\t\t3.5.2.6\tGO:0008800\t\t\t\t\t0\tNCBIFAM\tn\tn\t
TIGR01663.1\t\tlabelY\t100\t100\t200\tequivalog\tY\tY\tN\tbifunctional enzyme\t\t\t2.7.1.78,3.1.3.32\tGO:0001,GO:0002\t\t\t\t\t0\tJCVI\tn\tn\t
NF999999.9\t\tempty\t1\t1\t1\thypoth\tN\tN\tN\t\t\t\t\t\t\t\t\t\t0\tNCBIFAM\tn\tn\t";

    const PFAM_SNIPPET: &str = "\
!version date: 2026/06/01
!description: mapping
Pfam:PF00069 Pkinase > GO:protein kinase activity ; GO:0004672
Pfam:PF00069 Pkinase > GO:protein phosphorylation ; GO:0006468
Pfam:PF00001 7tm_1 > GO:G protein-coupled receptor activity ; GO:0004930";

    #[test]
    fn parse_pgap_fields() {
        let (by_ver, by_nover) = parse_pgap(PGAP_SNIPPET);
        // Empty-field row is dropped.
        assert_eq!(by_ver.len(), 2);
        let rec = by_ver.get("NF000459.2").unwrap();
        assert_eq!(rec.gene.as_deref(), Some("blaSPG"));
        assert_eq!(rec.product.as_deref(), Some("class B beta-lactamase"));
        assert_eq!(rec.ec, vec!["3.5.2.6"]);
        assert_eq!(rec.go, vec!["GO:0008800"]);
        // Multi EC/GO split.
        let rec2 = by_nover.get("TIGR01663").unwrap();
        assert_eq!(rec2.ec, vec!["2.7.1.78", "3.1.3.32"]);
        assert_eq!(rec2.go, vec!["GO:0001", "GO:0002"]);
    }

    #[test]
    fn parse_pfam2go_multi() {
        let map = parse_pfam2go(PFAM_SNIPPET);
        assert_eq!(map.len(), 2);
        assert_eq!(map.get("PF00069").unwrap(), &vec!["GO:0004672", "GO:0006468"]);
        assert_eq!(map.get("PF00001").unwrap(), &vec!["GO:0004930"]);
    }

    #[test]
    fn strip_version_cases() {
        assert_eq!(strip_version("NF000282.2"), "NF000282");
        assert_eq!(strip_version("TIGR00001.1"), "TIGR00001");
        assert_eq!(strip_version("PF00069.28"), "PF00069");
        assert_eq!(strip_version("TIGR00058"), "TIGR00058");
        // Non-numeric tail is left intact.
        assert_eq!(strip_version("foo.bar"), "foo.bar");
    }

    fn build_table() -> XrefTable {
        let (pgap, pgap_nover) = parse_pgap(PGAP_SNIPPET);
        let pfam_go = parse_pfam2go(PFAM_SNIPPET);
        XrefTable {
            pgap,
            pgap_nover,
            pfam_go,
        }
    }

    #[test]
    fn lookup_exact_and_version_agnostic() {
        let t = build_table();
        // Exact match.
        let h = t.lookup("NF000459.2").unwrap();
        assert_eq!(h.gene.as_deref(), Some("blaSPG"));
        // Different version -> resolves via version-stripped map.
        let h2 = t.lookup("NF000459.5").unwrap();
        assert_eq!(h2.gene.as_deref(), Some("blaSPG"));
        assert_eq!(h2.ec, vec!["3.5.2.6"]);
        // Pfam-only accession -> GO from pfam2go.
        let h3 = t.lookup("PF00069.28").unwrap();
        assert!(h3.gene.is_none());
        assert_eq!(h3.go, vec!["GO:0004672", "GO:0006468"]);
        // Unknown -> None.
        assert!(t.lookup("NF123456.1").is_none());
    }

    fn mk_feature(kind: FeatureKind, ann_acc: &str, score: f32) -> Feature {
        Feature {
            kind,
            contig: "c1".to_string(),
            id: "c1_1".to_string(),
            start: 1,
            end: 99,
            strand: 1,
            aa: Some("M".to_string()),
            partial5: false,
            partial3: false,
            annotations: vec![Annotation {
                source: "rustyhmmer:ncbifams".to_string(),
                accession: ann_acc.to_string(),
                name: "raw name".to_string(),
                score,
                evalue: Some(1e-30),
                ref_len: None,
            }],
            func: Functional::default(),
        }
    }

    #[test]
    fn enrich_fills_cds_and_skips_noncds() {
        let t = build_table();
        let mut feats = vec![
            // CDS with a matching NCBIfams hit (differing version).
            mk_feature(FeatureKind::Cds, "NF000459.9", 250.0),
            // Non-CDS with the same accession — must stay untouched.
            mk_feature(FeatureKind::Trna, "NF000459.2", 250.0),
        ];
        enrich(&mut feats, &t);

        // CDS enriched.
        assert_eq!(feats[0].func.gene.as_deref(), Some("blaSPG"));
        assert_eq!(feats[0].func.product.as_deref(), Some("class B beta-lactamase"));
        assert_eq!(feats[0].func.ec, vec!["3.5.2.6"]);
        assert_eq!(feats[0].func.go, vec!["GO:0008800"]);

        // Non-CDS untouched.
        assert!(feats[1].func.gene.is_none());
        assert!(feats[1].func.ec.is_empty());
    }

    #[test]
    fn enrich_picks_best_score_and_preserves_existing() {
        let t = build_table();
        let mut f = mk_feature(FeatureKind::Cds, "NF123456.1", 10.0); // unknown, low
        // Add a higher-scoring known hit.
        f.annotations.push(Annotation {
            source: "rustyhmmer:ncbifams".to_string(),
            accession: "NF000459.2".to_string(),
            name: "n".to_string(),
            score: 500.0,
            evalue: Some(1e-40),
            ref_len: None,
        });
        // Pre-existing gene must be preserved.
        f.func.gene = Some("preset".to_string());
        let mut feats = vec![f];
        enrich(&mut feats, &t);
        assert_eq!(feats[0].func.gene.as_deref(), Some("preset")); // not overwritten
        assert_eq!(feats[0].func.product.as_deref(), Some("class B beta-lactamase")); // filled
        assert_eq!(feats[0].func.ec, vec!["3.5.2.6"]);
    }

    #[test]
    fn enrich_falls_through_to_lower_resolving_annotation() {
        // F6: the highest-scoring annotation has NO curated metadata row; a
        // lower-scoring one does. enrich must fall through to the resolving hit.
        let t = build_table();
        let mut f = mk_feature(FeatureKind::Cds, "NF999999.9", 900.0); // top score, empty row -> unresolved
        f.annotations.push(Annotation {
            source: "rustyhmmer:ncbifams".to_string(),
            accession: "NF000459.2".to_string(), // lower score but resolves
            name: "n".to_string(),
            score: 300.0,
            evalue: Some(1e-30),
            ref_len: None,
        });
        let mut feats = vec![f];
        enrich(&mut feats, &t);
        assert_eq!(feats[0].func.gene.as_deref(), Some("blaSPG"));
        assert_eq!(feats[0].func.product.as_deref(), Some("class B beta-lactamase"));
    }
}
