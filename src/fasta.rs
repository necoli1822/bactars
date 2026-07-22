//! Genome input via needletail (plan §3: input = needletail, output = noodles).

use needletail::parse_fastx_file;
use std::collections::HashSet;

/// One contig: name (first header token) and its uppercased nucleotide bytes.
pub struct Contig {
    pub name: String,
    pub seq: Vec<u8>,
}

/// Read all contigs from a nucleotide FASTA/FASTQ file (plain or gzipped),
/// uppercasing residues. Contig name is the first whitespace-delimited header token.
pub fn read_genome(path: &str) -> Result<Vec<Contig>, String> {
    let mut reader = parse_fastx_file(path).map_err(|e| format!("{path}: {e}"))?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| format!("{path}: {e}"))?;
        let id = rec.id();
        let name = id
            .split(|&b| b == b' ' || b == b'\t')
            .next()
            .map(|t| String::from_utf8_lossy(t).into_owned())
            .unwrap_or_default();
        let seq = rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect();
        out.push(Contig { name, seq });
    }

    // Distinct contig ids matter: feature ids, GFF3 `##sequence-region` / seqid,
    // and GBFF LOCUS are all keyed on the contig name, and resolve()/dedup match
    // features by contig. Two records sharing a first-token id silently merge and
    // corrupt coordinates downstream — warn so it is diagnosable.
    let mut seen: HashSet<&str> = HashSet::with_capacity(out.len());
    for c in &out {
        if !seen.insert(c.name.as_str()) {
            eprintln!(
                "[bactars][warn] duplicate contig name '{}' in {path}: feature ids and \
                 GFF3/GBFF seqids will collide; use unique FASTA headers",
                c.name
            );
        }
    }
    Ok(out)
}
