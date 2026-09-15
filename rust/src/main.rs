//! NGSpeciesID — Rust port.
//!
//! Clustering and consensus of long-read amplicon data. See PORTING.md for the
//! specification, which is byte-identity with the Python reference in this same
//! repository, and `bench/README.md` for the harness that checks it.
//!
//! **Only the CLI exists so far.** Everything past argument validation exits
//! `EXIT_NOT_IMPLEMENTED` with one line saying so. That code is deliberately
//! neither 0, 1 nor 2: those three are the reference's own exit codes, and a
//! placeholder that returned any of them would let cases pass for the wrong
//! reason — `bench/equivalence.sh` compares exit codes, and the 15 traceback
//! cases all want exit 1.

mod align;
mod blockalign;
mod cli;
mod cluster;
mod consensus;
mod fastq;
mod minimizers;
mod p_emp;
mod packed;
mod parallelize;
mod parasail;
#[cfg(feature = "parasail-ffi")]
mod parasail_ffi;
mod phred;
mod pipeline;
mod poa;
mod pyfloat;
mod pyrandom;
mod pyround;
mod sorting;
mod sweep;
mod text;

use std::io::Write;
use std::process::ExitCode;

/// Not 0, 1 or 2. See the module docs: a placeholder must not be mistakable for
/// the reference's own exit codes. 70 is `EX_SOFTWARE` from sysexits.h, which
/// is as close to "this program is incomplete" as the convention gets.
pub const EXIT_NOT_IMPLEMENTED: u8 = 70;

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();

    match cli::parse(&argv) {
        cli::Outcome::Stdout(s) => {
            // `--help`, `--version` and `write_fastq --help` are the only three
            // things the reference writes to stdout. Everything else goes to
            // stderr, because it goes through `logging`.
            print!("{s}");
            flush_stdout();
            ExitCode::SUCCESS
        }
        cli::Outcome::UsageError(detail) => {
            eprint!("{}", text::usage_error(&detail));
            ExitCode::from(2)
        }
        cli::Outcome::SubUsageError(detail) => {
            // The subparser has its own prog and its own usage block.
            eprint!("{}", text::write_fastq_usage_error(&detail));
            ExitCode::from(2)
        }
        cli::Outcome::Stderr(msg, code) => {
            eprint!("{msg}");
            ExitCode::from(code as u8)
        }
        cli::Outcome::Run(args) => {
            // The reference creates --outfolder here: after the --ont/--isoseq
            // check (which exits 0 without creating anything) and after the
            // window check. `cli::validate` has already run both, so reaching
            // this point means the directory should exist.
            //
            // Measured on the reference: a run that fails the window check
            // leaves NO directory, and one that fails later leaves an empty
            // one. Reproducing that means creating it exactly here.
            debug_assert!(cli::creates_outfolder_before_window_check());
            let outfolder = args
                .outfolder
                .as_deref()
                .expect("validate() rejects a missing --outfolder");
            if let Err(e) = std::fs::create_dir_all(outfolder) {
                // The reference lets os.makedirs raise, which is a traceback.
                // One line instead.
                eprintln!("Error: could not create --outfolder {outfolder}: {e}");
                return ExitCode::from(1);
            }
            match pipeline::run(&args) {
                pipeline::Outcome::Done => ExitCode::SUCCESS,
                pipeline::Outcome::Failed(code) => ExitCode::from(code as u8),
                pipeline::Outcome::Incomplete(stage) => not_implemented(stage),
            }
        }
        cli::Outcome::WriteFastq(wf) => {
            // Reading the payload here keeps it honest: the fields are the ones
            // the stage will need, and an unread struct is a struct nobody has
            // checked against the reference.
            if wf.clusters.is_none() {
                eprintln!("Error: write_fastq needs --clusters.");
                return ExitCode::from(1);
            }
            not_implemented("write_fastq")
        }
    }
}

fn not_implemented(what: &str) -> ExitCode {
    eprintln!(
        "NGSpeciesID (Rust port): {what} is not implemented yet. \
         Arguments parsed and validated successfully."
    );
    ExitCode::from(EXIT_NOT_IMPLEMENTED)
}

/// `print!` to a closed stdout is a silent no-op unless the flush is checked,
/// and `--help | head -1` closes it. The reference dies with a BrokenPipeError
/// traceback there; exiting quietly is the better behaviour and the difference
/// is not in any golden.
fn flush_stdout() {
    let _ = std::io::stdout().flush();
}

#[cfg(test)]
mod tests {
    #[test]
    fn the_not_implemented_code_is_not_one_the_reference_uses() {
        // 0, 1 and 2 are the reference's exit codes. A placeholder sharing one
        // of them would make traceback cases "pass" without the port having
        // implemented anything.
        assert!(!matches!(super::EXIT_NOT_IMPLEMENTED, 0..=2));
    }
}
