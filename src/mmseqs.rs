//! Shared `mmseqs` binary resolution + runner.
//!
//! `mmseqs` is the ONE sanctioned external binary in bactars (everything else is
//! pure Rust). The protein-similarity stages (psc, vfdb, plasmid, species,
//! amr_variant) previously each carried a copy of the resolver + runner plus a
//! hardcoded dev-machine path; this module centralises them so there is a single
//! resolution policy and one env override.
//!
//! Resolution order (first that works):
//!   1. `$BACTARS_MMSEQS` — explicit override (absolute path or a name on `PATH`),
//!   2. the pinned build at [`PINNED_PATH`] if that file exists (a convenience for
//!      the calibration machine; NOT required),
//!   3. bare `mmseqs` from `PATH`.

use std::path::Path;
use std::process::Command;

/// A pinned mmseqs build path used ONLY as a fallback convenience when present.
/// Not machine-required: resolution falls through to `$PATH` when it is absent.
/// Override with `$BACTARS_MMSEQS` on any other host.
pub const PINNED_PATH: &str = "/usr/local/bin/mmseqs/mmseqs_18-8cc5c/bin/mmseqs";

/// Resolve the mmseqs binary to invoke. See module docs for the order.
pub fn bin() -> String {
    if let Some(v) = std::env::var_os("BACTARS_MMSEQS") {
        let s = v.to_string_lossy().into_owned();
        if !s.is_empty() {
            return s;
        }
    }
    if Path::new(PINNED_PATH).is_file() {
        return PINNED_PATH.to_string();
    }
    "mmseqs".to_string()
}

/// Run an mmseqs sub-command with the resolved binary; `Err(stderr-tail)` on a
/// non-zero exit (last 6 stderr lines, matching the previous per-module runner).
pub fn run(args: &[&str]) -> Result<(), String> {
    run_with(&bin(), args)
}

/// Run mmseqs using an explicitly-provided binary (for callers that resolve once
/// and reuse). Equivalent behaviour to [`run`].
pub fn run_with(bin: &str, args: &[&str]) -> Result<(), String> {
    let out = Command::new(bin)
        .args(args)
        .output()
        .map_err(|e| format!("spawn mmseqs {}: {e}", args.first().copied().unwrap_or("")))?;
    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        let tail: String = stderr
            .lines()
            .rev()
            .take(6)
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .collect::<Vec<_>>()
            .join("\n");
        return Err(format!(
            "mmseqs {} exited {}: {tail}",
            args.first().copied().unwrap_or(""),
            out.status
        ));
    }
    Ok(())
}

/// True if an mmseqs binary appears runnable (`mmseqs version` exits 0). Used by
/// the protein-similarity stages to skip cleanly (with a log) when it is absent,
/// instead of failing mid-pipeline.
pub fn available() -> bool {
    Command::new(bin())
        .arg("version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    // `$BACTARS_MMSEQS` is process-global; `cargo test` runs tests as parallel
    // threads in one process, so these two env-mutating tests must not interleave
    // (else one's remove_var races the other's set_var → flaky assert). Serialize
    // them on a shared lock. (Poison is recovered so one failure does not cascade.)
    static ENV_LOCK: Mutex<()> = Mutex::new(());

    #[test]
    fn env_override_takes_precedence() {
        let _g = ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        // Set an explicit override; bin() must return it verbatim.
        std::env::set_var("BACTARS_MMSEQS", "/custom/path/to/mmseqs");
        assert_eq!(bin(), "/custom/path/to/mmseqs");
        std::env::remove_var("BACTARS_MMSEQS");
    }

    #[test]
    fn falls_back_to_path_name() {
        let _g = ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        std::env::remove_var("BACTARS_MMSEQS");
        // When neither the override nor the pinned path is usable, resolution is a
        // bare `mmseqs` (found via PATH at spawn time). We can't assert the pinned
        // branch portably, but the fallback must never be empty.
        assert!(!bin().is_empty());
    }
}
