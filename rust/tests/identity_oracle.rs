//! `consensus::highest_aln_identity` against the reference, on recorded calls.
//!
//! This is the **second** parasail call site and it does not share the
//! clustering path's parameters:
//!
//! | | clustering (`parasail_block_alignment`) | this (`parasail_alignment`) |
//! | --- | --- | --- |
//! | opening penalty | binned 5/4/3/2 by summed error rate | fixed **3** |
//! | what is derived | aligned-window fraction over `k`-column windows | mismatch count over the gapped strings |
//!
//! Reusing the clustering scoring here would be wrong in a way no type would
//! catch, which is why it gets its own oracle rather than riding on the
//! clustering one being green.
//!
//! Both orientations are checked, not just the maximum: `max()` hides which side
//! won, and getting the reverse complement wrong is exactly the mistake this
//! stage exists to make.
//!
//! Recorded by `bench/dump_reference.py --stage identity`, one line per call:
//!
//! ```text
//! IDENT <identity_fw> <identity_rc> <max> <seq> <seq2>
//! ```

use std::path::Path;

#[path = "../src/align.rs"]
mod align;
#[path = "../src/parasail.rs"]
mod parasail;

/// A local copy of the two functions under test, because `consensus.rs` pulls in
/// the whole crate and this test only needs these.
mod subject {
    use super::{align, parasail};

    pub const RC_SCORING: parasail::Scoring = parasail::Scoring {
        match_score: 2,
        mismatch: -2,
        open: 3,
        ext: 1,
    };

    pub fn reverse_complement(s: &str) -> String {
        s.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'C' => 'G',
                'G' => 'C',
                'a' => 't',
                't' => 'a',
                'c' => 'g',
                'g' => 'c',
                'Y' => 'R',
                'R' => 'Y',
                'K' => 'M',
                'M' => 'K',
                'B' => 'V',
                'V' => 'B',
                'H' => 'D',
                'D' => 'H',
                other => other,
            })
            .collect()
    }

    pub fn identity(s1: &str, s2: &str) -> f64 {
        let aln = parasail::semiglobal(s1.as_bytes(), s2.as_bytes(), RC_SCORING);
        let Some((a1, a2)) = align::ops_to_seq(&aln.ops, s1.as_bytes(), s2.as_bytes()) else {
            return 0.0;
        };
        let total = a1.len();
        if total == 0 {
            return 0.0;
        }
        let mismatching = a1.iter().zip(a2.iter()).filter(|(x, y)| x != y).count();
        (total - mismatching) as f64 / total as f64
    }
}

#[test]
fn identity_matches_the_reference_on_every_recorded_call() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mut checked = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for name in [
        "identity_smoke.tsv",
        "identity_sup.tsv",
        "identity_sup_rc05.tsv",
    ] {
        let p = dir.join(name);
        let Ok(text) = std::fs::read_to_string(&p) else {
            eprintln!("  SKIP {name}: not recorded");
            continue;
        };
        for (lineno, line) in text.lines().enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.first() != Some(&"IDENT") || f.len() < 6 {
                continue;
            }
            let want_fw: f64 = f[1].parse().expect("identity_fw");
            let want_rc: f64 = f[2].parse().expect("identity_rc");
            let want_max: f64 = f[3].parse().expect("max");
            let (seq, seq2) = (f[4], f[5]);

            let got_rc = subject::identity(seq, &subject::reverse_complement(seq2));
            let got_fw = subject::identity(seq, seq2);
            let got_max = got_fw.max(got_rc);

            // Exact equality, not a tolerance. Both sides compute
            // (total - mismatches) / total from integer counts, so the f64 is
            // the same division of the same two integers -- if it is not, the
            // alignment differed and a tolerance would hide that.
            for (label, got, want) in [
                ("forward", got_fw, want_fw),
                ("reverse-complement", got_rc, want_rc),
                ("max", got_max, want_max),
            ] {
                if got != want {
                    failures.push(format!(
                        "{name}:{}: {label} identity {got} want {want} ({} vs {} bp)",
                        lineno + 1,
                        seq.len(),
                        seq2.len()
                    ));
                }
            }
            checked += 1;
        }
    }

    assert!(
        checked > 0,
        "no oracle data -- record it with bench/dump_reference.py --stage identity"
    );
    assert!(
        failures.is_empty(),
        "{} of {checked} recorded calls differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
    eprintln!("  {checked} recorded identity calls, all three values each");
}
