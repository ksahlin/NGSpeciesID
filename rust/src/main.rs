//! NGSpeciesID — Rust port.
//!
//! Clustering and consensus of long-read amplicon data. See PORTING.md for the
//! specification, which is byte-identity with the Python reference in this same
//! repository, and `bench/README.md` for the harness that checks it.
//!
//! **Every stage exists.** While the port was partial, anything past argument
//! validation exited `EXIT_NOT_IMPLEMENTED` (70) with one line saying so —
//! deliberately neither 0, 1 nor 2, because those three are the reference's own
//! exit codes and a placeholder returning one of them would let cases pass for
//! the wrong reason: `bench/equivalence.sh` compares exit codes, and the 15
//! traceback cases all want exit 1.
//!
//! That placeholder is gone, because nothing can reach it any more. Its absence
//! is the signal — as is `pipeline::Outcome` having no `Incomplete` variant. If
//! a future stage is stubbed out, bring the constant back rather than returning
//! 1: a stub that exits 1 is indistinguishable from a correctly reproduced
//! failure, and the harness would call it a pass.

mod align;
mod aligner;
mod blockalign;
mod cli;
mod cluster;
mod consensus;
mod edlib;
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
mod write_fastq;

use std::io::Write;
use std::process::ExitCode;

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
            }
        }
        cli::Outcome::WriteFastq(wf) => ExitCode::from(write_fastq::run(&wf) as u8),
    }
}

/// `print!` to a closed stdout is a silent no-op unless the flush is checked,
/// and `--help | head -1` closes it. The reference dies with a BrokenPipeError
/// traceback there; exiting quietly is the better behaviour and the difference
/// is not in any golden.
fn flush_stdout() {
    let _ = std::io::stdout().flush();
}
