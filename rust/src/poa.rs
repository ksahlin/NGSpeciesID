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

/// The POA consensus of `seqs`, weighted by `quals`, in insertion order.
///
/// **This is spoa's own C++ code**, vendored and linked by `spoa-sys`, not a
/// reimplementation that agrees with it. Exact by construction, which is what
/// byte-identity requires — see `tests/spoa_oracle.rs` for what happened when a
/// reimplementation was tried.
///
/// `quals` may be shorter than `seqs` or hold empty strings, in which case those
/// sequences weigh 1 per base — the FASTA path, which the reference never takes.
// Wired in by `form_draft_consensus`, which is the next slice.
#[allow(dead_code)]
pub fn consensus(seqs: &[String], quals: &[String]) -> String {
    if seqs.is_empty() {
        return String::new();
    }
    // `spoa <reads.fq> -l 0 -r 0 -g -2` resolves to kSW (local), m=5, n=-4,
    // g=-2 from the flag, and spoa's own defaults for e/q/c. Passing the
    // defaults explicitly rather than repeating -2 four times: spoa normalises
    // (g >= e means linear, then e := g) and the stored q/c are its own.
    let mut engine =
        crate::spoa::AlignmentEngine::new(crate::spoa::AlignmentType::kSW, 5, -4, -2, -6, -10, -4);
    let mut graph = crate::spoa::Graph::new();
    for (i, s) in seqs.iter().enumerate() {
        let b = s.as_bytes();
        let aln = engine.align(b, &graph);
        match quals.get(i) {
            // THE GRAPH IS QUALITY-WEIGHTED, and nothing in `run_spoa`'s
            // argument list says so: it hands spoa a FASTQ, and spoa's CLI
            // calls the quality overload whenever the input has qualities.
            // Measured: the same 20 sequences give 847 bp as FASTQ and 860 bp
            // as FASTA.
            Some(q) if q.len() == b.len() => graph.add_alignment_with_qual(&aln, b, q.as_bytes()),
            _ => graph.add_alignment(&aln, b, 1),
        }
    }
    String::from_utf8(graph.consensus()).expect("spoa's consensus is ASCII")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_sequences_gives_an_empty_consensus() {
        assert_eq!(consensus(&[], &[]), "");
    }

    #[test]
    fn one_sequence_is_its_own_consensus() {
        let s = "ACGTACGTACGTAAGGCCTTACGTACGT".to_string();
        assert_eq!(consensus(std::slice::from_ref(&s), &[]), s);
    }

    #[test]
    fn identical_sequences_give_that_sequence() {
        let s = "ACGTACGTACGTAAGGCCTTACGTACGT".to_string();
        assert_eq!(consensus(&[s.clone(), s.clone(), s.clone()], &[]), s);
    }

    /// A quality string of the wrong length falls back to weight 1 rather than
    /// panicking -- the binding asserts equal lengths, and a mismatched pair is
    /// a bug upstream of here, not a reason to abort a run.
    #[test]
    fn a_mismatched_quality_length_falls_back_to_unweighted() {
        let s = "ACGTACGTACGTAAGGCCTTACGTACGT".to_string();
        let short = "IIII".to_string();
        assert_eq!(consensus(std::slice::from_ref(&s), &[short]), s);
    }

    /// Quality weighting is not cosmetic: it changes the consensus. Two reads
    /// disagree at one position, and the higher-quality base wins even when it
    /// is in the minority.
    #[test]
    fn quality_weighting_changes_the_answer() {
        let a = "ACGTACGTACGTAAGGCCTTACGTACGTAC".to_string();
        let b = "ACGTACGTACGTAAGGCCTTACGTACGTAG".to_string();
        // Unweighted, two votes for ...AC beat one for ...AG.
        let unweighted = consensus(&[a.clone(), a.clone(), b.clone()], &[]);
        // Weighted, the single high-quality read outweighs two poor ones.
        let lowq: String = std::iter::repeat_n('!', a.len()).collect(); // phred 0
        let highq: String = std::iter::repeat_n('I', b.len()).collect(); // phred 40
        let weighted = consensus(
            &[a.clone(), a.clone(), b.clone()],
            &[lowq.clone(), lowq, highq],
        );
        assert_eq!(unweighted, a);
        assert_ne!(
            weighted, unweighted,
            "quality weights must reach the consensus"
        );
    }
}
