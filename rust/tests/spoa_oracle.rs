//! **MEASURED AND REJECTED: `spoars` does not reproduce `spoa` on this data.**
//!
//! This file is the record of that decision, not a live gate. It is kept so the
//! next person does not repeat the experiment, and so a future POA candidate has
//! an oracle to be measured against. `src/poa.rs` shells out to `spoa`.
//!
//! # What was measured
//!
//! Ten invocations recorded from the reference by
//! `bench/dump_reference.py --stage spoa`, across both corpora:
//!
//! | attempt | result |
//! | --- | --- |
//! | `spoars` 0.1.4, weight 1 per base | **0 of 10** identical (e.g. 860 bp against 847) |
//! | ...with spoa's real CLI defaults for `e`/`q`/`c` | **0 of 10**, byte-for-byte the same failures |
//! | ...quality-weighted, as spoa actually does | **0 of 10**, but much closer: 848 against 847 |
//!
//! And the control, which is what makes the conclusion safe: the real `spoa`
//! binary, given the recorded input reconstructed as a FASTQ, reproduces the
//! recorded consensus in **6 of 6** cases. So the oracle captures everything
//! spoa needs -- sequences, qualities and insertion order -- and the remaining
//! 1-7 bp differences are `spoars` itself.
//!
//! # The quality trap, which cost the first two attempts
//!
//! `run_spoa` hands spoa a **FASTQ**, and spoa's CLI weights the graph by
//! per-base quality whenever the input has any:
//!
//! ```cpp
//! if (it->quality.empty()) graph.AddAlignment(alignment, it->data);
//! else                     graph.AddAlignment(alignment, it->data, it->quality);
//! ```
//!
//! Nothing in `run_spoa`'s argument list says so. Measured: the same 20
//! sequences give an **847 bp** consensus as FASTQ and **860 bp** as FASTA.
//! isONcorrect passes a FASTA, which is both why its 505/505 `spoars` result is
//! real and why it says nothing about this repository.
//!
//! # Why isONcorrect's answer did not transfer anyway
//!
//! | | isONcorrect | NGSpeciesID |
//! | --- | --- | --- |
//! | sequences per POA | up to 28 | up to **1 198** |
//! | sequence length | correction intervals | up to **1 600 bp** |
//! | input format | FASTA | **FASTQ** |
//!
//! # Re-recording the data
//!
//! ```text
//! bench/dump_reference.py --stage spoa --sorted-fastq OUT/sorted.fastq \
//!     --k 13 --w 20 --abundance_ratio 0.02 --out rust/tests/data/spoa_sup.tsv
//! ```
//!
//! `SPOA <n_seqs> <consensus>` then one `SEQ <i> <accession> <sequence>
//! <quality>` per sequence, in insertion order -- which is part of the input,
//! because the order sequences enter a POA graph changes the consensus.

use std::path::Path;
use std::time::Instant;

struct Invocation {
    seqs: Vec<String>,
    quals: Vec<String>,
    consensus: String,
}

fn load(path: &Path) -> Vec<Invocation> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut out: Vec<Invocation> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f.first() {
            Some(&"SPOA") => out.push(Invocation {
                seqs: Vec::new(),
                quals: Vec::new(),
                consensus: f.get(2).copied().unwrap_or("").to_string(),
            }),
            Some(&"SEQ") => {
                if let Some(inv) = out.last_mut() {
                    inv.seqs.push(f.get(3).copied().unwrap_or("").to_string());
                    // Field 5 is the quality string. spoa weights the graph by
                    // it; a dump without it compares two different problems.
                    inv.quals.push(f.get(4).copied().unwrap_or("").to_string());
                }
            }
            _ => {}
        }
    }
    out
}

/// Kept as a harness for the next POA candidate. It does not build against
/// `src/poa.rs`, which is now a subprocess wrapper with nothing to compare --
/// wire a candidate's `consensus(&seqs, &quals)` in here and run it.
#[test]
#[ignore = "the record of a rejected candidate; see the module docs"]
fn a_native_poa_candidate_would_be_measured_here() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mut total = 0usize;
    let mut agreed = 0usize;
    let mut failures: Vec<String> = Vec::new();

    // spoa_sup.tsv -- the 1 198-sequence invocations -- is NOT committed: it is
    // 4.4 MB, and this repository is in the middle of being taken from 500 MB to
    // single digits. Regenerate it with the command in the module docs when a
    // candidate needs testing at that scale; its measurements are in PORTING.md.
    for name in ["spoa_smoke.tsv", "spoa_sup_max20.tsv", "spoa_sup.tsv"] {
        let p = dir.join(name);
        let invocations = load(&p);
        if invocations.is_empty() {
            eprintln!("  SKIP {name}: not recorded");
            continue;
        }
        for (i, inv) in invocations.iter().enumerate() {
            let longest = inv.seqs.iter().map(|s| s.len()).max().unwrap_or(0);
            let t0 = Instant::now();
            // let got = candidate::consensus(&inv.seqs, &inv.quals);
            let got = inv.consensus.clone(); // no candidate wired in
            let dt = t0.elapsed();
            total += 1;
            if got == inv.consensus {
                agreed += 1;
                eprintln!(
                    "  ok   {name}[{i}] {:>5} seqs, longest {longest:>5} bp, \
                     consensus {:>5} bp, {:>8.2?}",
                    inv.seqs.len(),
                    got.len(),
                    dt
                );
            } else {
                let first = got
                    .chars()
                    .zip(inv.consensus.chars())
                    .position(|(a, b)| a != b);
                failures.push(format!(
                    "{name}[{i}]: {} seqs, got {} bp want {} bp, first difference at {:?}, {:.2?}",
                    inv.seqs.len(),
                    got.len(),
                    inv.consensus.len(),
                    first,
                    dt
                ));
                eprintln!("  FAIL {}", failures.last().unwrap());
            }
        }
    }

    eprintln!("\n  {agreed} of {total} invocations identical");
    assert!(
        total > 0,
        "no oracle data -- see the module docs for how to record it"
    );
    assert!(
        failures.is_empty(),
        "{} of {total} invocations differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
