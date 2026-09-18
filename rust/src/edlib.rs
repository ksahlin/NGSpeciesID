//! edlib's `HW` (infix) mode with `task="locations"`, as `barcode_trimmer` uses it.
//!
//! ```python
//! edlib.align(primer_seq, center, mode="HW", task="locations",
//!             k=primer_max_ed, additionalEqualities=IUPAC_map)
//! ```
//!
//! `HW` means the **query must align in full** but the target is free at both
//! ends — an infix search. The reference reads `result["editDistance"]` and
//! `result["locations"][0]`, and nothing else.
//!
//! # Why this is a reimplementation and not a binding
//!
//! isONcorrect's port measured both Rust edlib bindings and found them
//! unusable: `edlib_rs` only builds with `CMAKE_POLICY_VERSION_MINIMUM=3.5`
//! because its vendored CMakeLists requires compatibility CMake 4 removed, and
//! `rsedlib` compiles and then fails to link. Neither problem is worth taking on
//! for a function this small — the queries are ~20 bp primers and the targets
//! are `--trim_window` bases, 150 by default, so a plain O(nm) DP is
//! microseconds.
//!
//! # What the corpora do and do not exercise
//!
//! 96 calls recorded across three configurations (`bench/dump_reference.py
//! --stage barcode`):
//!
//! | | no hit (`ed = -1`) | one location | **more than one** |
//! | --- | --- | --- | --- |
//! | `--primer_file` | 26 | 6 | **0** |
//! | `--primer_file --primer_max_ed 0` | 32 | 0 | **0** |
//! | `--remove_universal_tails` | 24 | 8 | **0** |
//!
//! **No recorded call returns more than one location**, so edlib's choice of
//! which equally-optimal location comes first is never exercised. That is the
//! same shape as *Findings 19, 25 and 26*: a behaviour the goldens cannot see.
//! `locations` is returned in ascending end position here, which is edlib's own
//! order, and `first_location` is documented as the only one the reference
//! reads — but if a corpus ever produces two, this is the thing to check first.

// Wired in by `remove_barcodes`, which is the next slice.
#![allow(dead_code)]

/// One infix match: the target range `[start, end]`, inclusive at both ends,
/// exactly as edlib reports it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Location {
    pub start: usize,
    pub end: usize,
}

/// What `edlib.align(..., mode="HW", task="locations")` returns, reduced to the
/// two fields the reference reads.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HwResult {
    /// `-1` when no alignment within `k` exists, exactly as edlib reports it.
    pub edit_distance: i64,
    /// Empty when `edit_distance` is `-1`.
    pub locations: Vec<Location>,
}

impl HwResult {
    /// `result["locations"][0]` — the only one `find_barcode_locations` reads.
    pub fn first_location(&self) -> Option<Location> {
        self.locations.first().copied()
    }
}

/// The reference's `additionalEqualities`, verbatim.
///
/// **Ordered pairs, and deliberately asymmetric.** It contains `('M','A')` but
/// not `('A','M')`, so an ambiguity code in the **query** (the primer) matches a
/// concrete base in the **target** (the center) and not the other way round.
/// Given that centers are spoa consensus sequences over ACGT, that is the right
/// way round — but it is data, and a port must not "symmetrise" it.
///
/// Note it also maps `X` as a full wildcard, which is not IUPAC; it comes from
/// the `Bio.Data` table the reference's comment cites.
pub const IUPAC_EQUALITIES: &[(u8, u8)] = &[
    (b'A', b'A'),
    (b'C', b'C'),
    (b'G', b'G'),
    (b'T', b'T'),
    (b'M', b'A'),
    (b'M', b'C'),
    (b'R', b'A'),
    (b'R', b'G'),
    (b'W', b'A'),
    (b'W', b'T'),
    (b'S', b'C'),
    (b'S', b'G'),
    (b'Y', b'C'),
    (b'Y', b'T'),
    (b'K', b'G'),
    (b'K', b'T'),
    (b'V', b'A'),
    (b'V', b'C'),
    (b'V', b'G'),
    (b'H', b'A'),
    (b'H', b'C'),
    (b'H', b'T'),
    (b'D', b'A'),
    (b'D', b'G'),
    (b'D', b'T'),
    (b'B', b'C'),
    (b'B', b'G'),
    (b'B', b'T'),
    (b'X', b'G'),
    (b'X', b'A'),
    (b'X', b'T'),
    (b'X', b'C'),
    (b'N', b'G'),
    (b'N', b'A'),
    (b'N', b'T'),
    (b'N', b'C'),
];

/// A 256x256 equality table built once from `IUPAC_EQUALITIES`, plus identity.
///
/// edlib treats identical bytes as equal regardless of the extra pairs, so the
/// diagonal is always set.
struct Equality {
    table: Vec<bool>,
}

impl Equality {
    fn new(extra: &[(u8, u8)]) -> Equality {
        let mut table = vec![false; 256 * 256];
        for i in 0..256 {
            table[i * 256 + i] = true;
        }
        for (a, b) in extra {
            table[*a as usize * 256 + *b as usize] = true;
        }
        Equality { table }
    }
    #[inline]
    fn eq(&self, q: u8, t: u8) -> bool {
        self.table[q as usize * 256 + t as usize]
    }
}

/// `edlib.align(query, target, mode="HW", task="locations", k, additionalEqualities)`.
///
/// `k < 0` means unbounded, as edlib's `-1` does. The reference always passes a
/// non-negative `--primer_max_ed`.
pub fn align_hw(query: &[u8], target: &[u8], k: i64, extra: &[(u8, u8)]) -> HwResult {
    let eq = Equality::new(extra);
    let m = query.len();
    let n = target.len();

    // An empty query matches everywhere at distance 0. edlib reports the whole
    // target's positions; the reference never asks, since primers are non-empty.
    if m == 0 {
        return HwResult {
            edit_distance: 0,
            locations: Vec::new(),
        };
    }
    if n == 0 {
        // The full query must be deleted; that costs m, which may exceed k.
        return if k >= 0 && m as i64 > k {
            HwResult {
                edit_distance: -1,
                locations: Vec::new(),
            }
        } else {
            HwResult {
                edit_distance: m as i64,
                locations: Vec::new(),
            }
        };
    }

    // dp[j] = edit distance of the whole query against a suffix of target[..j].
    //
    // Row 0 is all zeros, which is what makes this infix rather than global: an
    // alignment may start anywhere in the target for free. The query still has
    // to be consumed in full, so column 0 counts up.
    let mut prev: Vec<usize> = vec![0; n + 1];
    let mut curr: Vec<usize> = vec![0; n + 1];
    // Traceback would need the whole matrix; the windows are ~150 bases and the
    // queries ~20, so keeping it is 3 KB and simpler than a reverse pass.
    let mut rows: Vec<Vec<usize>> = Vec::with_capacity(m + 1);
    rows.push(prev.clone());

    for i in 1..=m {
        curr[0] = i;
        for j in 1..=n {
            let cost = usize::from(!eq.eq(query[i - 1], target[j - 1]));
            curr[j] = (prev[j - 1] + cost)
                .min(prev[j] + 1) // deletion from the target
                .min(curr[j - 1] + 1); // insertion into the target
        }
        std::mem::swap(&mut prev, &mut curr);
        rows.push(prev.clone());
    }

    let last = &rows[m];
    let best = *last[1..=n].iter().min().unwrap_or(&usize::MAX);
    if k >= 0 && best as i64 > k {
        return HwResult {
            edit_distance: -1,
            locations: Vec::new(),
        };
    }

    // Every end position achieving `best`, ascending -- edlib's own order.
    let mut locations = Vec::new();
    for (j, d) in last.iter().enumerate().take(n + 1).skip(1) {
        if *d == best {
            let start = trace_start(&rows, query, target, &eq, j);
            locations.push(Location { start, end: j - 1 });
        }
    }
    HwResult {
        edit_distance: best as i64,
        locations,
    }
}

/// Walk the DP matrix back from `(m, end_j)` to row 0, and report the target
/// index the alignment began at.
///
/// The preference order on ties is diagonal, then up (a deletion from the
/// target), then left (an insertion). edlib's own tie-break is not specified by
/// its API, and **no recorded call here produces more than one location**, so
/// this ordering is unexercised rather than verified. See the module docs.
fn trace_start(
    rows: &[Vec<usize>],
    query: &[u8],
    target: &[u8],
    eq: &Equality,
    end_j: usize,
) -> usize {
    let mut i = rows.len() - 1;
    let mut j = end_j;
    while i > 0 {
        let here = rows[i][j];
        if j > 0 {
            let cost = usize::from(!eq.eq(query[i - 1], target[j - 1]));
            if here == rows[i - 1][j - 1] + cost {
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if here == rows[i - 1][j] + 1 {
            i -= 1;
            continue;
        }
        if j > 0 && here == rows[i][j - 1] + 1 {
            j -= 1;
            continue;
        }
        // Unreachable for a well-formed matrix; stop rather than loop.
        break;
    }
    j
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hw(q: &str, t: &str, k: i64) -> HwResult {
        align_hw(q.as_bytes(), t.as_bytes(), k, IUPAC_EQUALITIES)
    }

    #[test]
    fn an_exact_infix_is_distance_zero() {
        let r = hw("ACGT", "TTTACGTTTT", 2);
        assert_eq!(r.edit_distance, 0);
        assert_eq!(r.first_location(), Some(Location { start: 3, end: 6 }));
    }

    #[test]
    fn the_query_must_align_in_full_but_the_target_is_free_at_both_ends() {
        // A prefix match costs nothing at the target's ends...
        assert_eq!(hw("ACGT", "ACGTTTTTTT", 0).edit_distance, 0);
        assert_eq!(hw("ACGT", "TTTTTTACGT", 0).edit_distance, 0);
        // ...but the query is not free: a missing base costs 1.
        assert_eq!(hw("ACGTA", "TTTACGTTTT", 0).edit_distance, -1);
        assert_eq!(hw("ACGTA", "TTTACGTTTT", 1).edit_distance, 1);
    }

    #[test]
    fn k_bounds_the_answer_and_minus_one_means_no_hit() {
        // One substitution.
        let t = "TTTACCTTTT";
        assert_eq!(hw("ACGT", t, 0).edit_distance, -1);
        assert_eq!(hw("ACGT", t, 1).edit_distance, 1);
        assert!(hw("ACGT", t, 0).locations.is_empty());
    }

    #[test]
    fn k_zero_requires_an_exact_match() {
        // The --primer_max_ed 0 case, which is 32 of the 96 recorded calls and
        // finds nothing in all of them.
        assert_eq!(hw("ACGT", "TTTACGTTTT", 0).edit_distance, 0);
        assert_eq!(hw("ACGT", "TTTACCTTTT", 0).edit_distance, -1);
    }

    #[test]
    fn iupac_codes_in_the_query_match_concrete_bases() {
        // Y is C or T; R is A or G. This is the primer file's own alphabet.
        assert_eq!(
            hw("ACAAATCAYAARGAYATYGG", "ACAAATCACAAAGACATCGG", 0).edit_distance,
            0
        );
        assert_eq!(
            hw("ACAAATCAYAARGAYATYGG", "ACAAATCATAAGGATATTGG", 0).edit_distance,
            0
        );
        // ...and a base the code does not cover still costs.
        assert_eq!(
            hw("ACAAATCAYAARGAYATYGG", "ACAAATCAGAAAGACATCGG", 0).edit_distance,
            -1
        );
    }

    /// The equality table is ORDERED. An ambiguity code in the target does not
    /// match a concrete base in the query, which is the direction the reference
    /// does not use and must not acquire.
    #[test]
    fn the_equality_table_is_asymmetric() {
        let eq = Equality::new(IUPAC_EQUALITIES);
        assert!(eq.eq(b'Y', b'C'), "Y in the query matches C in the target");
        assert!(!eq.eq(b'C', b'Y'), "but not the other way round");
        assert!(eq.eq(b'C', b'C'), "identity always holds");
    }

    #[test]
    fn an_empty_target_costs_the_whole_query() {
        assert_eq!(hw("ACGT", "", 10).edit_distance, 4);
        assert_eq!(hw("ACGT", "", 2).edit_distance, -1);
    }

    #[test]
    fn locations_are_ascending_by_end_position() {
        // Two equally good matches; edlib returns both, ascending.
        let r = hw("ACGT", "ACGTAAACGT", 0);
        assert_eq!(r.edit_distance, 0);
        assert_eq!(r.locations.len(), 2);
        assert!(r.locations[0].end < r.locations[1].end);
        assert_eq!(r.first_location(), Some(Location { start: 0, end: 3 }));
    }
}
