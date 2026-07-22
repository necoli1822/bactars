//! Genome-quality QC (CheckM-lite), computed in the MAIN pipeline from the
//! ncbifams annotations already produced — no extra database, no extra search.
//!
//! Uses the single-copy universal bacterial ribosomal proteins (the conserved
//! core of GTDB bac120 / CheckM's universal marker set) as markers. For a complete,
//! clean genome each marker is present exactly once; deviations quantify quality:
//!   * **completeness**  — fraction of expected markers present (assembly finished?)
//!   * **contamination** — markers present in >1 intact copy (mixed / contaminated)
//!   * **disrupted core** — markers present only as a pseudogene (degradation / a
//!     frameshift or internal stop in an essential gene — an assembly error or a
//!     genuinely reduced genome). This last signal is beyond CheckM1 and reuses the
//!     pseudogene detectors already in the pipeline.
//!
//! v1 marker set = the ribosomal proteins, which ncbifams names directly, so this
//! needs no new DB. A full bac120 HMM marker set is a precision follow-up.

use crate::feature::{Feature, FeatureKind};
use std::collections::HashMap;

/// Single-copy universal bacterial ribosomal-protein markers (L/S subunit ids).
/// Curated to the broadly single-copy set: the zinc/non-zinc paralog pairs
/// (L31, L33, L36, S14) and the non-universal ones (L25, S1, S21, L7/L12) are
/// omitted so a normal genome does not read as contaminated or incomplete. L34
/// (rpmH, ~46 aa) is also omitted: it is too short to be reliably gene-called /
/// HMM-named (missing on complete genomes in validation), so it would read as a
/// false incompleteness rather than a real quality signal.
pub const MARKERS: &[&str] = &[
    "L1", "L2", "L3", "L4", "L5", "L6", "L9", "L10", "L11", "L13", "L14", "L15",
    "L16", "L17", "L18", "L19", "L20", "L21", "L22", "L23", "L24", "L27", "L28",
    "L29", "L30", "L32", "L35", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
    "S9", "S10", "S11", "S12", "S13", "S15", "S16", "S17", "S18", "S19", "S20",
];

/// CheckM-lite genome-quality report.
#[derive(Debug, Clone, Default)]
pub struct QcReport {
    /// Number of expected markers ([`MARKERS`] length).
    pub expected: usize,
    /// Distinct markers present (intact or pseudogenised).
    pub present: usize,
    /// Markers found in >= 2 intact copies (contamination signal).
    pub duplicated: Vec<String>,
    /// Markers present ONLY as a pseudogene (disrupted essential gene).
    pub disrupted: Vec<String>,
    /// Expected markers with no copy at all.
    pub missing: Vec<String>,
    /// `100 * present / expected`.
    pub completeness: f64,
    /// `100 * duplicated / expected`.
    pub contamination: f64,
}

/// The ribosomal-marker key (`"L2"`, `"S7"`, …) for a CDS, from its product name,
/// or `None`. Matches e.g. `"50S ribosomal protein L2"` / `"30S ribosomal protein S7"`.
fn marker_key(f: &Feature) -> Option<String> {
    let product = f.func.product.as_deref()?;
    let lower = product.to_ascii_lowercase();
    let at = lower.find("ribosomal protein ")?;
    let rest = &product[at + "ribosomal protein ".len()..];
    // First whitespace/comma/paren-delimited token after the phrase.
    let tok = rest
        .split(|c: char| c.is_whitespace() || c == ',' || c == '(' || c == ';')
        .next()?;
    let up = tok.trim().to_ascii_uppercase();
    if MARKERS.contains(&up.as_str()) {
        Some(up)
    } else {
        None
    }
}

/// Compute the CheckM-lite report from the finalized feature set. Counts each
/// marker's intact vs pseudogenised copies; a marker present only as a pseudogene
/// counts toward `present` (the gene IS there) but is reported as `disrupted`.
pub fn compute(features: &[Feature]) -> QcReport {
    // marker key -> (intact copies, pseudogenised copies)
    let mut found: HashMap<String, (usize, usize)> = HashMap::new();
    for f in features {
        if f.kind != FeatureKind::Cds {
            continue;
        }
        if let Some(k) = marker_key(f) {
            let e = found.entry(k).or_insert((0, 0));
            if f.func.pseudogene {
                e.1 += 1;
            } else {
                e.0 += 1;
            }
        }
    }
    let expected = MARKERS.len();
    let present = found.len();
    let mut duplicated: Vec<String> = found
        .iter()
        .filter(|(_, (intact, _))| *intact >= 2)
        .map(|(k, _)| k.clone())
        .collect();
    let mut disrupted: Vec<String> = found
        .iter()
        .filter(|(_, (intact, pseudo))| *intact == 0 && *pseudo >= 1)
        .map(|(k, _)| k.clone())
        .collect();
    let mut missing: Vec<String> = MARKERS
        .iter()
        .filter(|m| !found.contains_key(**m))
        .map(|m| m.to_string())
        .collect();
    duplicated.sort();
    disrupted.sort();
    missing.sort();
    QcReport {
        expected,
        present,
        completeness: 100.0 * present as f64 / expected as f64,
        contamination: 100.0 * duplicated.len() as f64 / expected as f64,
        duplicated,
        disrupted,
        missing,
    }
}

/// One-line human summary (for stderr / logs).
pub fn summary_line(r: &QcReport) -> String {
    let disrupted = if r.disrupted.is_empty() {
        String::new()
    } else {
        format!(
            ", {} disrupted core gene(s): {}",
            r.disrupted.len(),
            r.disrupted.join(",")
        )
    };
    format!(
        "QC (CheckM-lite, {} ribosomal markers): completeness {:.1}%, contamination {:.1}%{}",
        r.expected, r.completeness, r.contamination, disrupted
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::{Feature, FeatureKind, Functional};

    fn ribo(product: &str, pseudo: bool) -> Feature {
        Feature {
            kind: FeatureKind::Cds,
            contig: "c".into(),
            id: "x".into(),
            start: 1,
            end: 300,
            strand: 1,
            aa: None,
            partial5: false,
            partial3: false,
            annotations: vec![],
            func: Functional {
                product: Some(product.to_string()),
                pseudogene: pseudo,
                ..Default::default()
            },
        }
    }

    #[test]
    fn marker_extraction() {
        assert_eq!(marker_key(&ribo("50S ribosomal protein L2", false)).as_deref(), Some("L2"));
        assert_eq!(marker_key(&ribo("30S ribosomal protein S7", false)).as_deref(), Some("S7"));
        // non-universal / paralog-prone markers are not in the set
        assert_eq!(marker_key(&ribo("30S ribosomal protein S1", false)), None);
        assert_eq!(marker_key(&ribo("hypothetical protein", false)), None);
    }

    #[test]
    fn complete_clean_genome() {
        let feats: Vec<Feature> = MARKERS
            .iter()
            .map(|m| {
                let kind = if m.starts_with('L') { "50S" } else { "30S" };
                ribo(&format!("{kind} ribosomal protein {m}"), false)
            })
            .collect();
        let r = compute(&feats);
        assert_eq!(r.present, MARKERS.len());
        assert!((r.completeness - 100.0).abs() < 1e-6);
        assert!(r.contamination.abs() < 1e-6);
        assert!(r.disrupted.is_empty() && r.missing.is_empty());
    }

    #[test]
    fn contamination_and_disruption() {
        let mut feats = vec![
            ribo("50S ribosomal protein L2", false),
            ribo("50S ribosomal protein L2", false), // duplicate -> contamination
            ribo("30S ribosomal protein S7", true),  // only pseudogene -> disrupted
        ];
        feats.push(ribo("50S ribosomal protein L3", false));
        let r = compute(&feats);
        assert_eq!(r.duplicated, vec!["L2"]);
        assert_eq!(r.disrupted, vec!["S7"]);
        assert!(r.present == 3); // L2, S7, L3
    }
}
