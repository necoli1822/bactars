//! HMM functional annotation via rustyhmmer (in-process).

pub use rustyhmmer::api::Cutoff;
use rustyhmmer::api::{HmmAnnotator, HmmHit};
use rustyhmmer::seqio::digitize_amino;

/// Search `proteins` (`(id, aa)` pairs) against an HMM database under `cutoff`.
/// Returns every reported hit; `HmmHit::target_name` is the protein id.
pub fn annotate(
    proteins: &[(String, String)],
    hmm_db: &str,
    cutoff: Cutoff,
) -> Result<Vec<HmmHit>, String> {
    let annotator = HmmAnnotator::from_hmm_file(hmm_db)?.with_cutoff(cutoff);
    let seqs = digitize_amino(proteins.iter().map(|(n, a)| (n.as_str(), a.as_str())));
    Ok(annotator.search(&seqs))
}
