//! Every fixed string the CLI emits, held verbatim.
//!
//! These are not hand-typed. Each one was extracted from a recorded golden
//! under `bench/golden/<corpus>/cli/` and lives in its own file next door, so
//! there is no escaping to get wrong and no transcription step to drift. The
//! test in `tests/text_matches_goldens.rs` re-derives them from the goldens and
//! fails if they have moved apart.
//!
//! Why constants rather than asking a library to reproduce them: argparse's
//! help layout and usage-line wrapping are argparse's, and no argument parser
//! in Rust produces them. The contract is the reference's exact bytes
//! (`bench/README.md`, "The 44 CLI cases are three classes"), so the bytes are
//! what is stored. See `cli.rs` for the four other axes on which argparse and
//! clap differ.

/// `--help` and `-h`, on **stdout**, exit 0. 8137 bytes of argparse layout,
/// including `ArgumentDefaultsHelpFormatter`'s `(default: X)` suffixes.
pub const HELP: &str = include_str!("text/help.txt");

/// `--version`, on **stdout**, exit 0. Deliberately not `CARGO_PKG_VERSION`.
pub const VERSION: &str = include_str!("text/version.txt");

/// `write_fastq --help`, on **stdout**, exit 0.
pub const WRITE_FASTQ_HELP: &str = include_str!("text/write_fastq_help.txt");

/// The 21-line usage block that precedes every argparse error, on **stderr**,
/// exit 2. Byte-identical across all 14 recorded exit-2 cases, which is why it
/// is one constant and not fourteen.
pub const USAGE: &str = include_str!("text/usage.txt");

/// `--ont --isoseq` together. On **stderr**, and it exits **0** — see
/// `cli::validate`. Note the trailing space before the newline: it is in the
/// reference's string literal and it is in the golden.
pub const PRESETS_EXCLUSIVE: &str = include_str!("text/presets_exclusive.txt");

/// `--w` outside `[--k, 100]`. On **stderr**, exit 1. One message for both
/// directions, and `100 < w` is checked before `w < k`.
pub const WINDOW_INVALID: &str = include_str!("text/window_invalid.txt");

/// The 2-line usage block for errors raised by the **`write_fastq` subparser**,
/// which has its own `prog` — `NGSpeciesID write_fastq` — and its own usage.
/// Measured: `write_fastq --N abc` reports
/// `NGSpeciesID write_fastq: error: argument --N: invalid int value: 'abc'`
/// under this block, not under the main one. Byte-identical to the first two
/// lines of `WRITE_FASTQ_HELP`, and derived from them.
pub const WRITE_FASTQ_USAGE: &str = include_str!("text/write_fastq_usage.txt");

/// An argparse error: the usage block, then `NGSpeciesID: error: <detail>`.
///
/// The program name is `NGSpeciesID` and not `argv[0]`'s basename. argparse
/// derives its `prog` from `sys.argv[0]`, so the reference would say
/// `NGSpeciesID` when installed and something else if the script were renamed;
/// the goldens were recorded through the repository path and say `NGSpeciesID`.
/// Hard-coding it is what makes the port's output independent of where its
/// binary sits, which is what a golden needs.
pub fn usage_error(detail: &str) -> String {
    format!("{USAGE}NGSpeciesID: error: {detail}\n")
}

/// The same, for an error the `write_fastq` subparser raises: its own usage
/// block and its own prog.
pub fn write_fastq_usage_error(detail: &str) -> String {
    format!("{WRITE_FASTQ_USAGE}NGSpeciesID write_fastq: error: {detail}\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constants_are_not_empty_and_end_in_a_newline() {
        // A truncated include_str! is silent, and a missing trailing newline is
        // the sort of thing that shows up as one failing golden much later.
        for (name, s) in [
            ("HELP", HELP),
            ("VERSION", VERSION),
            ("WRITE_FASTQ_HELP", WRITE_FASTQ_HELP),
            ("USAGE", USAGE),
            ("WRITE_FASTQ_USAGE", WRITE_FASTQ_USAGE),
            ("PRESETS_EXCLUSIVE", PRESETS_EXCLUSIVE),
            ("WINDOW_INVALID", WINDOW_INVALID),
        ] {
            assert!(!s.is_empty(), "{name} is empty");
            assert!(s.ends_with('\n'), "{name} does not end in a newline");
        }
    }

    #[test]
    fn usage_block_is_the_recorded_shape() {
        assert_eq!(USAGE.lines().count(), 21, "the usage block is 21 lines");
        assert!(USAGE.starts_with("usage: NGSpeciesID [-h] [--version] [--debug]\n"));
        assert!(USAGE.ends_with("{write_fastq} ...\n"));
    }

    #[test]
    fn version_is_the_reference_string_not_the_crate_version() {
        assert_eq!(VERSION, "NGSpeciesID 0.4.0\n");
    }

    #[test]
    fn presets_message_keeps_its_trailing_space() {
        // Removing it would be an improvement and a divergence. It is in the
        // reference's literal, so it is in the port.
        assert!(PRESETS_EXCLUSIVE.ends_with("--ont. \n"));
    }

    #[test]
    fn the_subparser_usage_is_the_head_of_its_help() {
        // Derived from WRITE_FASTQ_HELP rather than typed, so they cannot drift.
        assert_eq!(WRITE_FASTQ_USAGE.lines().count(), 2);
        assert!(WRITE_FASTQ_HELP.starts_with(WRITE_FASTQ_USAGE));
        assert!(WRITE_FASTQ_USAGE.starts_with("usage: NGSpeciesID write_fastq [-h]"));
    }

    #[test]
    fn subparser_errors_name_the_subparser() {
        let e = write_fastq_usage_error("argument --N: invalid int value: 'abc'");
        assert!(e.starts_with(WRITE_FASTQ_USAGE));
        assert_eq!(
            e.lines().last().unwrap(),
            "NGSpeciesID write_fastq: error: argument --N: invalid int value: 'abc'"
        );
        assert_eq!(e.lines().count(), 3);
    }

    #[test]
    fn usage_error_is_usage_then_one_line() {
        let e = usage_error("argument --k: invalid int value: 'abc'");
        assert!(e.starts_with(USAGE));
        assert_eq!(
            e.lines().last().unwrap(),
            "NGSpeciesID: error: argument --k: invalid int value: 'abc'"
        );
        assert_eq!(e.lines().count(), 22, "21 lines of usage plus one of error");
    }
}
