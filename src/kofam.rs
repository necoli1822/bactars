//! KOfam KO assignment + KEGG xref (KO → EC / pathway). Runs a KOfam HMM search
//! via rustyhmmer with per-KO adaptive thresholds (ko_list), assigns a KEGG
//! Orthology id to each CDS, and fills func.ko / func.ec / func.pathway.
//!
//! This is the KofamScan method: a single HMM search over the concatenated
//! prokaryotic KOfam profiles (`kofam_prok.hmm`) with a permissive reporting
//! cutoff, then a per-KO adaptive bit-score threshold from `ko_list` decides
//! which candidate hits are real KO assignments. EC numbers come from the
//! `ko_list` definition column (`[EC:...]`); pathways come from the KEGG BRITE
//! `ko00001.keg` hierarchy (a KO's enclosing `C` pathway line).

use crate::feature::{Feature, FeatureKind};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Which score of a hit the per-KO threshold applies to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScoreType {
    /// Compare against the full-sequence score (`seq_score`).
    Full,
    /// Compare against the best single-domain score (`best_dom_score`).
    Domain,
    /// No usable threshold (`-`) — the KO is uncallable, skip it.
    None,
}

/// Per-KO record distilled from `ko_list`.
#[derive(Clone, Debug)]
pub struct KoInfo {
    /// Adaptive bit-score threshold, or `None` when the row's threshold is `-`.
    pub threshold: Option<f32>,
    /// Which hit score the threshold is compared against.
    pub score_type: ScoreType,
    /// EC numbers parsed from the `[EC:...]` block in the definition.
    pub ec: Vec<String>,
    /// Free-text definition (product-ish description).
    pub definition: String,
}

/// Loaded KOfam reference tables (per-KO thresholds/EC + KO→pathway + KO→COG).
pub struct KofamTable {
    pub ko: HashMap<String, KoInfo>,
    pub pathways: HashMap<String, Vec<String>>,
    /// KO → COG functional-category bitmask (26 categories, one bit each; see
    /// [`COG_ORDER`] / [`cog_letters`]).
    pub cog: HashMap<String, u32>,
}

/// Fixed COG functional-category letter order; bit i = the letter at index i.
/// Standard COG single-letter categories (J..S).
pub const COG_ORDER: &[u8; 26] = b"JAKLBDYVTMNZWUOXCGEFHIPQRS";

/// Encode a set of COG category letters into a u32 bitmask (unknown letters
/// ignored).
fn cog_mask(letters: &str) -> u32 {
    let mut m = 0u32;
    for b in letters.bytes() {
        if let Some(i) = COG_ORDER.iter().position(|&c| c == b) {
            m |= 1 << i;
        }
    }
    m
}

/// Decode a COG category bitmask back to its letters in [`COG_ORDER`] order.
pub fn cog_letters(mask: u32) -> String {
    let mut s = String::new();
    for (i, &c) in COG_ORDER.iter().enumerate() {
        if mask & (1 << i) != 0 {
            s.push(c as char);
        }
    }
    s
}

impl KofamTable {
    /// Load `ko_list` (required) and `ko00001.keg` (optional, best-effort) for a
    /// KOfam bundle directory. The BRITE file is searched in the bundle dir and a
    /// few conventional sibling locations (e.g. `../meta/ko00001.keg`); if not
    /// found, pathway assignment is silently skipped.
    pub fn load(kofam_dir: &str) -> Result<Self, String> {
        let dir = Path::new(kofam_dir);
        let ko_list_path = dir.join("ko_list");
        let text = std::fs::read_to_string(&ko_list_path)
            .map_err(|e| format!("kofam: reading {}: {e}", ko_list_path.display()))?;
        let ko = parse_ko_list(&text);

        let pathways = match find_brite(dir) {
            Some(p) => match std::fs::read_to_string(&p) {
                Ok(t) => parse_brite(&t),
                Err(e) => {
                    eprintln!("kofam: reading {}: {e} (pathways skipped)", p.display());
                    HashMap::new()
                }
            },
            None => {
                eprintln!("kofam: ko00001.keg not found near {kofam_dir} (pathways skipped)");
                HashMap::new()
            }
        };

        // KO → COG-category bitmask (db/meta/ko_cog.tsv), best-effort like BRITE.
        let cog = match find_meta(dir, "ko_cog.tsv") {
            Some(p) => match std::fs::read_to_string(&p) {
                Ok(t) => parse_ko_cog(&t),
                Err(e) => {
                    eprintln!("kofam: reading {}: {e} (COG skipped)", p.display());
                    HashMap::new()
                }
            },
            None => {
                eprintln!("kofam: ko_cog.tsv not found near {kofam_dir} (COG skipped)");
                HashMap::new()
            }
        };

        Ok(Self { ko, pathways, cog })
    }
}

/// Parse `ko_cog.tsv` (`KO<TAB>COGletters`) into KO → category bitmask.
fn parse_ko_cog(text: &str) -> HashMap<String, u32> {
    let mut map = HashMap::new();
    for line in text.lines() {
        let mut it = line.split('\t');
        let (Some(ko), Some(letters)) = (it.next(), it.next()) else {
            continue;
        };
        let ko = ko.trim();
        if ko.is_empty() || !ko.starts_with('K') {
            continue;
        }
        let m = cog_mask(letters.trim());
        if m != 0 {
            map.insert(ko.to_string(), m);
        }
    }
    map
}

/// Locate a `db/meta/<name>` companion file for a KOfam bundle dir (same probe
/// order as [`find_brite`]).
fn find_meta(dir: &Path, name: &str) -> Option<PathBuf> {
    let mut cands: Vec<PathBuf> = vec![dir.join(name)];
    if let Some(parent) = dir.parent() {
        cands.push(parent.join("meta").join(name));
        cands.push(parent.join(name));
        if let Some(gp) = parent.parent() {
            cands.push(gp.join("meta").join(name));
            cands.push(gp.join("db/meta").join(name));
        }
    }
    cands.into_iter().find(|p| p.is_file())
}

/// Locate the KEGG BRITE `ko00001.keg` for a KOfam bundle dir. Checks the dir
/// itself and a handful of conventional neighbours (the resource actually lives
/// under `db/meta/` in this project, not in the kofam dir).
fn find_brite(dir: &Path) -> Option<PathBuf> {
    let mut cands: Vec<PathBuf> = vec![dir.join("ko00001.keg")];
    if let Some(parent) = dir.parent() {
        cands.push(parent.join("meta/ko00001.keg"));
        cands.push(parent.join("ko00001.keg"));
        if let Some(gp) = parent.parent() {
            cands.push(gp.join("meta/ko00001.keg"));
            cands.push(gp.join("db/meta/ko00001.keg"));
        }
    }
    cands.into_iter().find(|p| p.is_file())
}

/// Parse `[EC:1.1.1.1 2.7.2.4 1.1.1.-]` out of a definition string into a list
/// of EC numbers (regex-free). Returns an empty vec if there is no EC block.
fn parse_ec(def: &str) -> Vec<String> {
    let mut out = Vec::new();
    // Find "[EC:" then read until the closing ']'.
    let Some(start) = def.find("[EC:") else {
        return out;
    };
    let rest = &def[start + 4..];
    let end = rest.find(']').unwrap_or(rest.len());
    for tok in rest[..end].split_whitespace() {
        let tok = tok.trim();
        if !tok.is_empty() {
            out.push(tok.to_string());
        }
    }
    out
}

/// Parse the `ko_list` TSV into a KO → [`KoInfo`] map. Skips the header row and
/// any malformed lines. Columns: knum, threshold, score_type, profile_type,
/// F-measure, nseq, nseq_used, alen, mlen, eff_nseq, re/pos, definition.
pub fn parse_ko_list(text: &str) -> HashMap<String, KoInfo> {
    let mut map = HashMap::new();
    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        // Skip header.
        if line.starts_with("knum\t") || line.starts_with("knum ") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let knum = cols[0].trim();
        if knum.is_empty() {
            continue;
        }
        let threshold = match cols[1].trim() {
            "-" | "" => None,
            v => v.parse::<f32>().ok(),
        };
        let score_type = match cols[2].trim() {
            "full" => ScoreType::Full,
            "domain" => ScoreType::Domain,
            _ => ScoreType::None,
        };
        let definition = cols.last().copied().unwrap_or("").trim().to_string();
        let ec = parse_ec(&definition);
        map.insert(
            knum.to_string(),
            KoInfo {
                threshold,
                score_type,
                ec,
                definition,
            },
        );
    }
    map
}

/// Parse the KEGG BRITE `ko00001.keg` htext into a KO → pathway-ids map. A `D`
/// KO line's pathway is the `[PATH:ko#####]` id of its enclosing `C` line;
/// non-pathway `C` sections (`[BR:...]`) are ignored. A KO can appear under
/// several pathways.
pub fn parse_brite(text: &str) -> HashMap<String, Vec<String>> {
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let mut cur_pathway: Option<String> = None;
    for line in text.lines() {
        let mut chars = line.chars();
        let tag = chars.next();
        match tag {
            Some('C') => {
                // C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
                cur_pathway = extract_bracket_id(line, "[PATH:");
            }
            Some('D') => {
                if let Some(path) = &cur_pathway {
                    // D      K00844  HK; hexokinase [EC:2.7.1.1]
                    let rest = line[1..].trim_start();
                    let ko = rest.split_whitespace().next().unwrap_or("");
                    if ko.starts_with('K') && ko.len() >= 2 {
                        let e = map.entry(ko.to_string()).or_default();
                        if !e.contains(path) {
                            e.push(path.clone());
                        }
                    }
                }
            }
            // A / B lines and everything else: reset nothing except that a new
            // B section leaves cur_pathway until the next C line overwrites it.
            _ => {}
        }
    }
    map
}

/// Extract the id inside a `[PREFIX....]` marker (e.g. `[PATH:ko00010]` → `ko00010`).
fn extract_bracket_id(line: &str, prefix: &str) -> Option<String> {
    let start = line.find(prefix)?;
    let rest = &line[start + prefix.len()..];
    let end = rest.find(']')?;
    let id = rest[..end].trim();
    if id.is_empty() {
        None
    } else {
        Some(id.to_string())
    }
}

/// Annotate CDS features with KO / EC / pathway from a KOfam bundle dir
/// (`kofam_prok.hmm` + `ko_list` + `ko00001.keg`).
pub fn annotate(features: &mut [Feature], kofam_dir: &str) -> Result<(), String> {
    let table = KofamTable::load(kofam_dir)?;

    // Build the protein list from CDS features.
    let proteins: Vec<(String, String)> = features
        .iter()
        .filter(|f| f.kind == FeatureKind::Cds)
        .filter_map(|f| f.aa.as_ref().map(|aa| (f.id.clone(), aa.clone())))
        .collect();
    if proteins.is_empty() {
        return Ok(());
    }

    // Permissive reporting cutoff — we apply the real per-KO threshold ourselves.
    let hmm_db = Path::new(kofam_dir).join("kofam_prok.hmm");
    let hmm_db = hmm_db.to_string_lossy().to_string();
    let hits = crate::hmm::annotate(&proteins, &hmm_db, crate::hmm::Cutoff::Evalue(1e-2))?;

    // Per protein id, keep the best passing KO. We rank by margin over the KO's
    // own threshold, tie-broken by the raw applied score.
    struct Best {
        ko: String,
        margin: f32,
        score: f32,
    }
    let mut best: HashMap<String, Best> = HashMap::new();

    for h in &hits {
        let Some(info) = table.ko.get(&h.query_name) else {
            continue;
        };
        let Some(thr) = info.threshold else {
            continue; // uncallable KO ('-')
        };
        let applied = match info.score_type {
            ScoreType::Full => h.seq_score,
            ScoreType::Domain => h.best_dom_score,
            ScoreType::None => continue,
        };
        if applied < thr {
            continue; // below the KO's adaptive threshold
        }
        let margin = applied - thr;
        let cand = Best {
            ko: h.query_name.clone(),
            margin,
            score: applied,
        };
        match best.get(&h.target_name) {
            Some(cur)
                if (cur.margin, cur.score) >= (cand.margin, cand.score) => {}
            _ => {
                best.insert(h.target_name.clone(), cand);
            }
        }
    }

    // Fill the winning KO onto each CDS feature.
    let idx: HashMap<&str, usize> = features
        .iter()
        .enumerate()
        .map(|(i, f)| (f.id.as_str(), i))
        .collect();

    // Collect (feature index, KO) first to avoid borrow conflicts.
    let assignments: Vec<(usize, String)> = best
        .iter()
        .filter_map(|(target, b)| idx.get(target.as_str()).map(|&i| (i, b.ko.clone())))
        .collect();

    for (i, ko) in assignments {
        if features[i].kind != FeatureKind::Cds {
            continue;
        }
        let info = match table.ko.get(&ko) {
            Some(v) => v,
            None => continue,
        };
        let empty = Vec::new();
        let paths = table.pathways.get(&ko).unwrap_or(&empty);
        let f = &mut features[i];
        f.func.ko = Some(ko.clone());
        for ec in &info.ec {
            if !f.func.ec.contains(ec) {
                f.func.ec.push(ec.clone());
            }
        }
        for p in paths {
            if !f.func.pathway.contains(p) {
                f.func.pathway.push(p.clone());
            }
        }
        // COG functional category bitmask for this KO (union onto any existing).
        f.func.cog_cat |= table.cog.get(&ko).copied().unwrap_or(0);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const KO_LIST_SAMPLE: &str = "\
knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition
K00001\t210.97\tdomain\tall\t0.674093\t326\t256\t966\t374\t9.37\t0.590\talcohol dehydrogenase [EC:1.1.1.1]
K00004\t379.10\tfull\tall\t0.854087\t1918\t1605\t1270\t428\t6.05\t0.590\t(R,R)-butanediol dehydrogenase / diacetyl reductase [EC:1.1.1.4 1.1.1.- 1.1.1.303]
K00492\t-\t-\t-\t-\t1\t1\t396\t396\t0.39\t0.593\t1,3,7-trimethyluric acid 5-monooxygenase [EC:1.14.13.212]
K99999\t100.0\tfull\tall\t0.9\t1\t1\t1\t1\t1\t0.5\tsome protein with no ec";

    #[test]
    fn ko_list_threshold_and_score_type() {
        let m = parse_ko_list(KO_LIST_SAMPLE);
        assert_eq!(m.len(), 4);

        let k1 = &m["K00001"];
        assert_eq!(k1.threshold, Some(210.97));
        assert_eq!(k1.score_type, ScoreType::Domain);
        assert_eq!(k1.ec, vec!["1.1.1.1"]);

        let k4 = &m["K00004"];
        assert_eq!(k4.threshold, Some(379.10));
        assert_eq!(k4.score_type, ScoreType::Full);
        assert_eq!(k4.ec, vec!["1.1.1.4", "1.1.1.-", "1.1.1.303"]);

        // '-' threshold → uncallable
        let k492 = &m["K00492"];
        assert_eq!(k492.threshold, None);
        assert_eq!(k492.score_type, ScoreType::None);
        assert_eq!(k492.ec, vec!["1.14.13.212"]);

        // no EC block
        let k9 = &m["K99999"];
        assert!(k9.ec.is_empty());
        assert_eq!(k9.score_type, ScoreType::Full);
    }

    #[test]
    fn ec_parse_edge_cases() {
        assert!(parse_ec("plain product, no ec").is_empty());
        assert_eq!(parse_ec("thing [EC:1.2.3.4]"), vec!["1.2.3.4"]);
        assert_eq!(
            parse_ec("multi [EC:1.1.1.1 2.7.2.4]"),
            vec!["1.1.1.1", "2.7.2.4"]
        );
    }

    const BRITE_SAMPLE: &str = "\
+D\tKO
A09100 Metabolism
B
B  09101 Carbohydrate metabolism
C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
D      K00844  HK; hexokinase [EC:2.7.1.1]
D      K12407  GCK; glucokinase [EC:2.7.1.2]
C    00020 Citrate cycle (TCA cycle) [PATH:ko00020]
D      K00844  also here somehow
B  09180 Brite Hierarchies
C    01000 Enzymes [BR:ko01000]
D      K00001  E1.1.1.1; alcohol dehydrogenase";

    #[test]
    fn brite_ko_to_pathway() {
        let m = parse_brite(BRITE_SAMPLE);
        // K00844 appears under two pathways
        assert_eq!(m["K00844"], vec!["ko00010", "ko00020"]);
        assert_eq!(m["K12407"], vec!["ko00010"]);
        // K00001 sits under a [BR:...] section, not a pathway → no entry
        assert!(!m.contains_key("K00001"));
    }

    // End-to-end on a subset of the real MG1655 proteome + real KOfam HMM DB.
    // Slow (loads the 2.5GB HMM); run with `cargo test -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn real_search_mg1655_subset() {
        use crate::feature::{Feature, Functional};
        let faa = concat!(env!("CARGO_MANIFEST_DIR"), "/../prodigal/bench_results/faa/MG1655.rust.faa");
        let text = std::fs::read_to_string(faa).expect("read proteome");
        // Parse FASTA; keep the first N proteins to bound runtime.
        let mut feats: Vec<Feature> = Vec::new();
        let mut id = String::new();
        let mut aa = String::new();
        let push = |id: &str, aa: &str, feats: &mut Vec<Feature>| {
            if !id.is_empty() && !aa.is_empty() {
                feats.push(Feature {
                    kind: FeatureKind::Cds,
                    contig: "c".into(),
                    id: id.into(),
                    start: 1,
                    end: 1,
                    strand: 1,
                    aa: Some(aa.into()),
                    partial5: false,
                    partial3: false,
                    annotations: Vec::new(),
                    func: Functional::default(),
                });
            }
        };
        const N: usize = 200;
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                push(&id, &aa, &mut feats);
                if feats.len() >= N {
                    id.clear();
                    aa.clear();
                    break;
                }
                id = h.split_whitespace().next().unwrap_or("").to_string();
                aa.clear();
            } else {
                aa.push_str(line.trim());
            }
        }
        push(&id, &aa, &mut feats);
        feats.truncate(N);
        eprintln!("proteins: {}", feats.len());

        annotate(&mut feats, concat!(env!("CARGO_MANIFEST_DIR"), "/../db/kofam"))
            .expect("kofam annotate");

        let with_ko = feats.iter().filter(|f| f.func.ko.is_some()).count();
        eprintln!("CDS with KO: {}/{}", with_ko, feats.len());
        for f in feats.iter().filter(|f| f.func.ko.is_some()).take(8) {
            eprintln!(
                "  {} -> {} ec={:?} path={:?}",
                f.id,
                f.func.ko.as_deref().unwrap(),
                f.func.ec,
                f.func.pathway
            );
        }
        assert!(with_ko > 0, "expected some KO assignments");
    }

    // Loads the REAL reference tables — run with `cargo test -- --ignored`.
    #[test]
    #[ignore]
    fn load_real_tables() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../db/kofam");
        let t = KofamTable::load(dir).expect("load kofam tables");
        assert!(t.ko.len() > 20000, "ko_list rows: {}", t.ko.len());
        assert!(
            t.pathways.len() > 5000,
            "KO→pathway entries: {}",
            t.pathways.len()
        );
        let callable = t.ko.values().filter(|k| k.threshold.is_some()).count();
        eprintln!(
            "ko_list rows={} callable={} pathway_KOs={}",
            t.ko.len(),
            callable,
            t.pathways.len()
        );
    }
}
