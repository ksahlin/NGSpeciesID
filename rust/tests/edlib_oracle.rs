//! `edlib::align_hw` against the real edlib, on recorded calls.
//!
//! `barcode_trimmer.find_barcode_locations` calls
//! `edlib.align(primer, window, mode="HW", task="locations", k=primer_max_ed,
//! additionalEqualities=IUPAC_map)` and reads the edit distance and
//! `locations[0]`. This replays every such call the reference made.
//!
//! Recorded by `bench/dump_reference.py --stage barcode`, one line per call:
//!
//! ```text
//! EDLIB <primer> <target> <k> <edit_distance> <start,end;start,end;…>
//! ```
//!
//! # What the recorded calls cover, and what they do not
//!
//! | configuration | no hit (`ed = -1`) | one location | more than one |
//! | --- | --- | --- | --- |
//! | `--primer_file` | 26 | 6 | **0** |
//! | `--primer_file --primer_max_ed 0` | 32 | 0 | **0** |
//! | `--remove_universal_tails` | 24 | 8 | **0** |
//!
//! So the no-hit path is covered heavily, single hits are covered, `k = 0` is
//! covered — and **the multiple-location case is not covered at all**. Which
//! equally-optimal location edlib puts first is therefore unverified here, and
//! `align_hw`'s ordering is a guess that happens never to be tested. That is the
//! same shape as *Findings 19, 25 and 26*, and it is stated rather than papered
//! over: if a corpus ever produces two locations, this is the first thing to
//! check.
//!
//! The assertion below counts the classes it saw, so a future recording that
//! *does* produce multiple locations shows up as a changed count rather than as
//! silent extra coverage.

use std::path::Path;

#[path = "../src/edlib.rs"]
mod edlib;

#[test]
fn hw_alignment_matches_edlib_on_every_recorded_call() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let mut checked = 0usize;
    let mut no_hit = 0usize;
    let mut one_loc = 0usize;
    let mut many_loc = 0usize;
    let mut failures: Vec<String> = Vec::new();

    for name in [
        "edlib_primer.tsv",
        "edlib_primer_ed0.tsv",
        "edlib_tails.tsv",
    ] {
        let p = dir.join(name);
        let Ok(text) = std::fs::read_to_string(&p) else {
            eprintln!("  SKIP {name}: not recorded");
            continue;
        };
        for (lineno, line) in text.lines().enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.first() != Some(&"EDLIB") || f.len() < 6 {
                continue;
            }
            let (query, target) = (f[1], f[2]);
            let k: i64 = f[3].parse().unwrap_or(-1);
            let want_ed: i64 = f[4].parse().expect("edit distance");
            let want_locs: Vec<(usize, usize)> = if f[5].is_empty() {
                Vec::new()
            } else {
                f[5].split(';')
                    .filter_map(|p| {
                        let (a, b) = p.split_once(',')?;
                        Some((a.parse().ok()?, b.parse().ok()?))
                    })
                    .collect()
            };

            let got = edlib::align_hw(
                query.as_bytes(),
                target.as_bytes(),
                k,
                edlib::IUPAC_EQUALITIES,
            );

            if got.edit_distance != want_ed {
                failures.push(format!(
                    "{name}:{}: edit distance {} want {want_ed} (query {} bp, target {} bp, k={k})",
                    lineno + 1,
                    got.edit_distance,
                    query.len(),
                    target.len()
                ));
            }
            let got_locs: Vec<(usize, usize)> =
                got.locations.iter().map(|l| (l.start, l.end)).collect();
            if got_locs != want_locs {
                failures.push(format!(
                    "{name}:{}: locations {got_locs:?} want {want_locs:?}",
                    lineno + 1
                ));
            }

            match want_locs.len() {
                0 => no_hit += 1,
                1 => one_loc += 1,
                _ => many_loc += 1,
            }
            checked += 1;
        }
    }

    assert!(
        checked > 0,
        "no oracle data -- record it with bench/dump_reference.py --stage barcode"
    );
    assert!(
        failures.is_empty(),
        "{} of {checked} recorded calls differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
    eprintln!("  {checked} recorded edlib HW calls: {no_hit} no hit, {one_loc} one location, {many_loc} more than one");
    // Stated as a fact about the corpora, not as a requirement. If this ever
    // fires, the tie-break ordering in `trace_start` has become reachable and
    // needs verifying rather than assuming.
    assert_eq!(
        many_loc, 0,
        "a recorded call now returns multiple locations -- edlib's ordering is \
         no longer unexercised, and align_hw's guess needs checking"
    );
}
