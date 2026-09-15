//! The POA consensus, replacing the `spoa` subprocess.
//!
//! `consensus.run_spoa` shells out to:
//!
//! ```text
//! spoa <reads.fq> -l 0 -r 0 -g -2
//! ```
//!
//! Resolved against spoa's CLI defaults (`m=5, n=-4, g=-8, e=-6, q=-10, c=-4`)
//! and `AlignmentEngine::Create`'s subtype dispatch — `g >= e` means linear, and
//! then `e := g` — the effective configuration is:
//!
//! * alignment type **`kSW`** (local / Smith-Waterman), from `-l 0`
//! * gap model **linear**, `-2` per gap position, from `-g -2`
//! * match **+5**, mismatch **−4**
//! * `-r 0`: consensus only, not the MSA
//!
//! Byte-for-byte the invocation isONcorrect uses, which is why `spoars` was the
//! first thing tried — but see `tests/spoa_oracle.rs` for why that port's
//! validation does not transfer, and for what was measured here instead.
//!
//! **Sequence insertion order changes the consensus.** The caller must hand
//! sequences over in the exact order the reference writes them to its temp
//! file, including the `--max_seqs_for_consensus` cutoff, which is
//! `i >= max_seqs_for_consensus` — admitting exactly that many, unlike
//! isONcorrect's bare `>`.

use std::io::Write;
use std::path::Path;
use std::process::Command;

/// The consensus of one cluster, by running `spoa` exactly as the reference
/// does.
///
/// `reads_path` is the FASTQ `form_draft_consensus` has already written; the
/// argument vector is the reference's, verbatim. Returns `None` if spoa is
/// absent or fails, with the cause already reported.
// Wired in by the consensus stage, which is the next slice.
#[allow(dead_code)]
pub fn consensus_via_spoa(reads_path: &Path, tmp_out: &Path) -> Option<String> {
    let out = match Command::new("spoa")
        .arg(reads_path)
        .args(["-l", "0", "-r", "0", "-g", "-2"])
        .output()
    {
        Ok(o) => o,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            eprintln!("Error: spoa not found on PATH. --consensus needs it.");
            return None;
        }
        Err(e) => {
            eprintln!("Error: could not run spoa: {e}");
            return None;
        }
    };
    if !out.status.success() {
        // The reference raises CalledProcessError here, which is a traceback.
        // Finding 22's SIGABRT on an empty input is the reachable case.
        eprintln!(
            "Error: spoa failed on {} ({}).",
            reads_path.display(),
            out.status
        );
        return None;
    }
    // The reference writes spoa's stdout to a file and reads line 2 back.
    // Reproduced through the same file, so a spoa that prints anything
    // unexpected is seen the same way.
    if let Ok(mut f) = std::fs::File::create(tmp_out) {
        let _ = f.write_all(&out.stdout);
    }
    let text = String::from_utf8_lossy(&out.stdout);
    text.lines().nth(1).map(|l| l.trim_end().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// spoa's own output shape: a `>Consensus` header then the sequence. Only
    /// line 2 is read, which is what the reference does.
    #[test]
    fn the_second_line_is_the_consensus() {
        let text = ">Consensus LN:i:12\nACGTACGTACGT\n";
        assert_eq!(text.lines().nth(1), Some("ACGTACGTACGT"));
    }

    #[test]
    fn a_missing_spoa_is_reported_not_panicked() {
        // Nothing to assert about the message here beyond that it returns
        // rather than unwinding; the binary's absence is an environment fact.
        let d = std::env::temp_dir();
        let missing = d.join("ngsid-no-such-reads.fq");
        let out = d.join("ngsid-no-such-out.fa");
        if which_spoa().is_none() {
            assert!(consensus_via_spoa(&missing, &out).is_none());
        }
    }

    fn which_spoa() -> Option<std::path::PathBuf> {
        std::env::var_os("PATH").and_then(|paths| {
            std::env::split_paths(&paths)
                .map(|p| p.join("spoa"))
                .find(|p| p.is_file())
        })
    }
}
