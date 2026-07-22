//! Standalone CRISPR-only runner for benchmarking the Rust MinCED port against
//! the original MinCED (Java) in isolation (no CDS/HMM pipeline overhead).
//!
//!   cargo run --release --example crispr_only -- <genome.fna> [--gff]
//!
//! Default: one line per array (contig, start, end, n_repeats, repeat).
//! `--gff`: MinCED's exact `-gff` output, for byte-parity diffing.

use std::io::Write;

fn main() {
    let mut path: Option<String> = None;
    let mut gff = false;
    for a in std::env::args().skip(1) {
        match a.as_str() {
            "--gff" => gff = true,
            _ => path = Some(a),
        }
    }
    let path = match path {
        Some(p) => p,
        None => {
            eprintln!("usage: crispr_only <genome.fna> [--gff]");
            std::process::exit(2);
        }
    };
    let t0 = std::time::Instant::now();
    let contigs: Vec<(String, String)> = bactars::fasta::read_genome(&path)
        .unwrap_or_else(|e| {
            eprintln!("error: {e}");
            std::process::exit(1);
        })
        .into_iter()
        .map(|c| (c.name, String::from_utf8_lossy(&c.seq).into_owned()))
        .collect();
    let t_read = t0.elapsed();

    let t1 = std::time::Instant::now();
    let feats = bactars::crispr::detect_crispr(&contigs).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1);
    });
    let t_detect = t1.elapsed();
    eprintln!(
        "read={:.3}s detect={:.3}s",
        t_read.as_secs_f64(),
        t_detect.as_secs_f64()
    );

    let out = std::io::stdout();
    let mut w = out.lock();
    if gff {
        // MinCED-exact GFF (byte-parity target).
        w.write_all(bactars::crispr::to_minced_gff(&feats).as_bytes()).ok();
    } else {
        for f in &feats {
            let a = &f.annotations[0];
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                f.contig, f.start, f.end, a.score as i64, a.accession
            )
            .ok();
        }
    }
    eprintln!("{} CRISPR arrays", feats.len());
}
