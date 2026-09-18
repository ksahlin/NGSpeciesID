//! `pyrandom` against CPython, on recorded draws.
//!
//! `pyrandom_cases.tsv` was produced by running `random.Random(seed).sample(
//! range(n), k)` in the reference interpreter; its own header says how. Each row
//! is `seed \t n \t k \t indices`, in **sample order** rather than sorted, so a
//! port that produces the right set in the wrong order still fails.
//!
//! Spot checks in `pyrandom.rs` cover the shapes by hand; this covers the axes:
//!
//! * both corpora's real sizes (274 and 3 000 surviving reads) at the matrix's
//!   `--sample_size` 100 and 200;
//! * seeds 0, 1, 7, 42, **-5** (the absolute-value path) and two above 2^32
//!   (the two-key-word path);
//! * the branch boundary, at `setsize(k) - 1`, `setsize(k)` and `setsize(k) + 1`
//!   for five values of `k`, so both algorithms are exercised either side of
//!   the switch;
//! * degenerate shapes: `k == 0`, `k == n`, `n == 1`.
//!
//! The first version of the hand-written spot checks had plausible-looking
//! wrong constants -- they were the `random()` stream rather than the
//! `getrandbits(32)` one, and every one of them was recalled rather than run.
//! This file exists so that cannot happen quietly again.

use std::path::Path;

#[path = "../src/pyrandom.rs"]
mod pyrandom;

#[test]
fn sample_matches_cpython_on_every_recorded_case() {
    let p = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/pyrandom_cases.tsv");
    let text = match std::fs::read_to_string(&p) {
        Ok(t) => t,
        Err(e) => panic!(
            "{} is missing ({e}); regenerate it, see its header",
            p.display()
        ),
    };

    let mut checked = 0usize;
    let mut pool_branch = 0usize;
    let mut set_branch = 0usize;
    for (lineno, line) in text.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            f.len(),
            4,
            "line {}: want 4 tab-separated fields",
            lineno + 1
        );
        let seed: i64 = f[0].parse().expect("seed");
        let n: usize = f[1].parse().expect("n");
        let k: usize = f[2].parse().expect("k");
        let want: Vec<usize> = if f[3].is_empty() {
            Vec::new()
        } else {
            f[3].split(',').map(|s| s.parse().expect("index")).collect()
        };

        let mut rng = pyrandom::PyRandom::seeded(seed);
        let got = pyrandom::sample_indices(&mut rng, n, k);
        assert_eq!(got, want, "line {}: seed={seed} n={n} k={k}", lineno + 1);

        if (n as f64) <= pyrandom::setsize_for(k) {
            pool_branch += 1;
        } else {
            set_branch += 1;
        }
        checked += 1;
    }

    assert!(
        checked >= 60,
        "only {checked} cases -- the file looks truncated"
    );
    // A pass that exercised one branch would say nothing about the other, and
    // the whole reason both corpora matter is that they straddle the switch.
    assert!(pool_branch >= 10, "only {pool_branch} pool-branch cases");
    assert!(set_branch >= 10, "only {set_branch} selection-set cases");
    eprintln!("  {checked} cases: {pool_branch} pool branch, {set_branch} selection set");
}
