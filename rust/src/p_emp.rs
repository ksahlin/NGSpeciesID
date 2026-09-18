//! The empirical minimizer-sharing probabilities.
//!
//! The reference keeps these in `modules/p_minimizers_shared.py`: a 1.79 MB
//! Python literal of 41 880 `(k, w, p, e1, e2)` rows, filtered at startup to the
//! rows where `k == args.k` and `abs(w - args.w) <= 2`, keyed by the rounded
//! error-rate pair, and inserted under both orderings.
//!
//! Here it is a packed binary blob (`p_emp_probs.bin`, 328 KB) generated from
//! that file by `tools/gen_p_emp_blob.py`, which also asserts the two facts
//! below and can re-verify the committed blob with `--check`.
//!
//! **This is NOT isONclust's blob**, and the difference is deliberate.
//! isONclust's table has 59 628 rows covering `k` 4..30, because its PR #13
//! added `k` 4..9; NGSpeciesID never took that PR and its table is exactly
//! isONclust's `k >= 10` subset -- same values, same order, checked. Shipping
//! the larger blob would silently make `--k 9` work, which the reference does
//! not: it builds an empty dict and dies with `KeyError`. That is a behaviour
//! change wearing a data change's clothing. 349 `(k, w)` pairs here against
//! isONclust's 465. See PORTING.md, *Finding 8*.
//!
//! Two facts measured from the table make the encoding safe:
//!
//! * **At most one `w` ever matches the ±2 filter.** The table's `w` values step
//!   by 5 for each `k`, so the filter selects exactly one row group -- asserted
//!   by the generator over every `(k, args_w)` combination. The reference's
//!   dict-overwrite behaviour therefore has nothing to overwrite, and there is
//!   no last-row-wins subtlety to reproduce.
//! * **Only 15 error rates are reachable.** `p_shared_minimizer_empirical`
//!   rounds to two decimals and clamps to `[0.01, 0.15]`, so a row outside that
//!   grid can never be looked up. What is left is a complete 15x15 grid, stored
//!   as its 120-entry upper triangle, for each of 349 `(k, w)` pairs.
//!
//! `tools/gen_p_emp_blob.py --check` re-derives the blob from the reference's
//! own file and fails if the two have drifted.

/// Error rates are `0.01 .. 0.15` in hundredths, indexed 0..15.
const N_E: usize = 15;
const TRI: usize = N_E * (N_E + 1) / 2; // 120
const HEADER: usize = 14 + 4;
const PAIR_BYTES: usize = 2 + TRI * 8;

static BLOB: &[u8] = include_bytes!("p_emp_probs.bin");

/// Upper-triangle index for an unordered pair of error-rate indices.
#[inline]
fn tri_index(i: usize, j: usize) -> usize {
    let (i, j) = if i <= j { (i, j) } else { (j, i) };
    i * N_E - i * (i.wrapping_sub(1)) / 2 + (j - i)
}

/// The 120 probabilities for one `(k, w)` pair, or `None` when the table has no
/// row group within ±2 of `w`.
///
/// `None` is not an error here: the reference builds an empty dict and then dies
/// with a `KeyError` on a tuple of two floats the first time it looks anything
/// up. **3 361 CLI-valid `(k, w)` combinations reach it** -- every `k` below 10
/// or above 30, plus twelve inside that range where `w` is 99 or 100. See
/// PORTING.md, *Finding 8*.
pub struct Table {
    values: &'static [u8],
}

impl Table {
    /// Select the row group for `(k, args_w)`, applying the reference's
    /// `abs(w - args.w) <= 2` filter.
    pub fn select(k: i64, args_w: i64) -> Option<Table> {
        if !(0..=255).contains(&k) {
            return None;
        }
        let n = u32::from_le_bytes(BLOB[14..18].try_into().expect("header")) as usize;
        for idx in 0..n {
            let off = HEADER + idx * PAIR_BYTES;
            let (rk, rw) = (BLOB[off] as i64, BLOB[off + 1] as i64);
            if rk == k && (rw - args_w).abs() <= 2 {
                return Some(Table {
                    values: &BLOB[off + 2..off + PAIR_BYTES],
                });
            }
        }
        None
    }

    /// `p_emp_probs[(e1, e2)]` for error-rate indices already rounded and
    /// clamped. Symmetric, as the reference's dict is.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let t = tri_index(i, j);
        let b = &self.values[t * 8..t * 8 + 8];
        f64::from_le_bytes(b.try_into().expect("8 bytes"))
    }

    /// Every stored entry, for the oracle.
    #[allow(dead_code)]
    pub fn entries(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        (0..N_E).flat_map(move |i| (i..N_E).map(move |j| (i, j, self.get(i, j))))
    }
}

/// Every `(k, w)` pair in the blob, for the oracle.
#[allow(dead_code)]
pub fn all_pairs() -> Vec<(i64, i64)> {
    let n = u32::from_le_bytes(BLOB[14..18].try_into().expect("header")) as usize;
    (0..n)
        .map(|idx| {
            let off = HEADER + idx * PAIR_BYTES;
            (BLOB[off] as i64, BLOB[off + 1] as i64)
        })
        .collect()
}

/// `p_shared_minimizer_empirical`'s rounding: two decimals, then clamped to
/// `[0.01, 0.15]`, expressed as an index into the error-rate axis.
///
/// The clamp is applied to the *rounded* value, exactly as the reference does,
/// so a raw error rate of 0.153 becomes 0.15 and 0.004 becomes 0.01.
pub fn error_rate_index(e: f64) -> usize {
    let r = crate::pyround::round2(e);
    // The reference compares the rounded float against 0.15 and 0.01.
    if r > 0.15 {
        return N_E - 1;
    }
    if r < 0.01 {
        return 0;
    }
    // r is a multiple of 0.01 in [0.01, 0.15]; recover the index without
    // trusting float arithmetic to land exactly.
    let hundredths = (r * 100.0).round() as i64;
    (hundredths.clamp(1, 15) - 1) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tri_index_covers_the_triangle_exactly_once() {
        let mut seen = vec![0u32; TRI];
        for i in 0..N_E {
            for j in i..N_E {
                seen[tri_index(i, j)] += 1;
                // symmetric
                assert_eq!(tri_index(i, j), tri_index(j, i));
            }
        }
        assert!(
            seen.iter().all(|c| *c == 1),
            "triangle index is not a bijection"
        );
    }

    #[test]
    fn the_blob_has_the_expected_shape() {
        assert_eq!(&BLOB[..13], b"NGSPECIESPEMP");
        assert_eq!(all_pairs().len(), 349, "k 10..30, not isONclust's 4..30");
    }

    /// The reference's ont and isoseq presets must both resolve.
    #[test]
    fn presets_select_a_table() {
        assert!(Table::select(13, 20).is_some(), "--ont");
        assert!(Table::select(15, 50).is_some(), "--isoseq");
    }

    /// Finding 8: 3 361 CLI-valid settings have no row group within +-2, and
    /// the reference dies with a KeyError on all of them.
    ///
    /// Two kinds, and both are in the goldens:
    ///
    /// * `k` outside 10..30 -- the whole range isONclust's PR #13 added and this
    ///   repository did not. `--k 9` and `--k 31` are cli/k_too_small and
    ///   cli/k_too_big.
    /// * `k` inside the range with `w` at the very top, where the row group
    ///   steps past 100. `--k 11 --w 100` is cli/kw_gap.
    #[test]
    fn the_settings_with_no_table_are_exactly_the_ones_measured() {
        for (k, w) in [(9, 20), (4, 10), (1, 1), (31, 50), (100, 100)] {
            assert!(
                Table::select(k, w).is_none(),
                "k={k} w={w} is outside the table's k range and must have no table"
            );
        }
        // Inside 10..30, the twelve gaps at the top of the w axis.
        for (k, w) in [
            (11, 99),
            (11, 100),
            (12, 100),
            (16, 99),
            (16, 100),
            (17, 100),
            (21, 99),
            (21, 100),
            (22, 100),
            (26, 99),
            (26, 100),
            (27, 100),
        ] {
            assert!(
                Table::select(k, w).is_none(),
                "k={k} w={w} should have no table (Finding 8)"
            );
        }
        // Count them all, so a change to the blob that quietly widens coverage
        // shows up here rather than as a behaviour change nobody asked for.
        let mut missing = 0usize;
        for k in 1..=100i64 {
            for w in k..=100 {
                if Table::select(k, w).is_none() {
                    missing += 1;
                }
            }
        }
        assert_eq!(
            missing, 3361,
            "the number of CLI-valid (k, w) pairs with no table"
        );
    }

    #[test]
    fn probabilities_are_symmetric_and_in_range() {
        let t = Table::select(13, 20).expect("--ont table");
        for i in 0..N_E {
            for j in 0..N_E {
                let p = t.get(i, j);
                assert_eq!(p, t.get(j, i), "asymmetric at ({i},{j})");
                assert!((0.0..=1.0).contains(&p), "p out of range at ({i},{j}): {p}");
            }
        }
    }

    #[test]
    fn error_rates_are_rounded_then_clamped() {
        assert_eq!(error_rate_index(0.0001), 0); // below the floor
        assert_eq!(error_rate_index(0.01), 0);
        assert_eq!(error_rate_index(0.104), 9); // rounds to 0.10
        assert_eq!(error_rate_index(0.15), 14);
        assert_eq!(error_rate_index(0.9), 14); // above the ceiling
    }
}
