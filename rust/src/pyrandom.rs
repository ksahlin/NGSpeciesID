//! CPython's Mersenne Twister and `random.sample`, reproduced exactly.
//!
//! `--sample_size` draws `random.Random(args.seed).sample(range(n), k)` and the
//! port has to produce the same subset. That is a harder contract than it looks
//! — the *set* of indices depends on the generator's exact bit stream and on
//! which of `sample`'s two selection algorithms runs — but it is entirely
//! specified, and it is the same class of work as `pyfloat.rs` and `pyround.rs`.
//!
//! Sources: CPython's `Modules/_randommodule.c` (`init_genrand`,
//! `init_by_array`, `genrand_uint32`, `getrandbits`) and `Lib/random.py`
//! (`_randbelow_with_getrandbits`, `sample`). Checked against the interpreter by
//! `tests/pyrandom_oracle.rs`, which replays recorded draws.
//!
//! # The two branches, and why both are covered
//!
//! `sample` picks its algorithm from a size heuristic:
//!
//! ```text
//! setsize = 21
//! if k > 5: setsize += 4 ** ceil(log(k * 3, 4))
//! if n <= setsize:  pool branch            (partial Fisher-Yates over a copy)
//! else:             selection-set branch   (draw and retry against a set)
//! ```
//!
//! At the case matrix's `--sample_size` 100 and 200 that threshold is 1045, and
//! the two committed corpora sit either side of it: `smoke` has 274 surviving
//! reads and takes the **pool** branch, `sup` has 3 000 and takes the
//! **selection-set** branch. So a port implementing only one fails on exactly
//! one corpus — which is luck worth banking, and worth not losing.
//!
//! # What is NOT reproduced
//!
//! `getrandbits(k)` for `k > 32`. CPython's fast path is one 32-bit word; above
//! that it assembles a big integer word by word. `_randbelow` only ever asks for
//! `n.bit_length()` bits where `n` is a read count, so 32 bits covers 4.29
//! billion reads. Asking for more is a programming error and panics rather than
//! returning something plausible.

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// CPython's `RandomObject`: the MT19937 state and its position.
pub struct PyRandom {
    mt: [u32; N],
    mti: usize,
}

impl PyRandom {
    /// `random.Random(seed)` for an integer seed.
    ///
    /// CPython takes the seed's **absolute value** and splits it into 32-bit
    /// words, little-endian, so `--seed -5` and `--seed 5` give the same stream.
    /// Measured against the interpreter, and worth knowing before anyone treats
    /// the sign as a second axis.
    ///
    /// A seed of 0 is not a special case in the generator: it produces a
    /// one-word key `[0]`, because CPython computes
    /// `keyused = bits == 0 ? 1 : (bits - 1) / 32 + 1`.
    pub fn seeded(seed: i64) -> PyRandom {
        // `unsigned_abs`, not `abs`: i64::MIN has no positive counterpart and
        // `abs` would panic in debug and wrap in release.
        let n = seed.unsigned_abs();
        let bits = 64 - n.leading_zeros();
        let keyused = if bits == 0 { 1 } else { (bits - 1) / 32 + 1 } as usize;
        let mut key = [0u32; 2];
        key[0] = (n & 0xffff_ffff) as u32;
        key[1] = (n >> 32) as u32;
        let mut r = PyRandom {
            mt: [0; N],
            mti: N + 1,
        };
        r.init_by_array(&key[..keyused]);
        r
    }

    fn init_genrand(&mut self, s: u32) {
        self.mt[0] = s;
        for i in 1..N {
            let prev = self.mt[i - 1];
            self.mt[i] = 1812433253u32
                .wrapping_mul(prev ^ (prev >> 30))
                .wrapping_add(i as u32);
        }
        self.mti = N;
    }

    fn init_by_array(&mut self, key: &[u32]) {
        self.init_genrand(19650218);
        let mut i = 1usize;
        let mut j = 0usize;
        let mut k = N.max(key.len());
        while k > 0 {
            let prev = self.mt[i - 1];
            self.mt[i] = (self.mt[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1664525))
                .wrapping_add(key[j])
                .wrapping_add(j as u32);
            i += 1;
            j += 1;
            if i >= N {
                self.mt[0] = self.mt[N - 1];
                i = 1;
            }
            if j >= key.len() {
                j = 0;
            }
            k -= 1;
        }
        k = N - 1;
        while k > 0 {
            let prev = self.mt[i - 1];
            self.mt[i] = (self.mt[i] ^ (prev ^ (prev >> 30)).wrapping_mul(1566083941))
                .wrapping_sub(i as u32);
            i += 1;
            if i >= N {
                self.mt[0] = self.mt[N - 1];
                i = 1;
            }
            k -= 1;
        }
        // MSB is 1, assuring a non-zero initial array.
        self.mt[0] = 0x8000_0000;
    }

    /// One 32-bit output, with MT19937's tempering.
    pub fn genrand_u32(&mut self) -> u32 {
        if self.mti >= N {
            for kk in 0..N - M {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            for kk in N - M..N - 1 {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] =
                    self.mt[kk + M - N] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            let y = (self.mt[N - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            self.mti = 0;
        }
        let mut y = self.mt[self.mti];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^ (y >> 18)
    }

    /// `getrandbits(k)` for `k <= 32`. CPython's fast path is one word shifted
    /// right, and `k == 0` returns 0 **without consuming a word**.
    pub fn getrandbits(&mut self, k: u32) -> u32 {
        assert!(
            k <= 32,
            "getrandbits({k}): only the <=32 fast path is ported"
        );
        if k == 0 {
            return 0;
        }
        self.genrand_u32() >> (32 - k)
    }

    /// `_randbelow_with_getrandbits(n)`: rejection sampling on `n.bit_length()`
    /// bits.
    ///
    /// The comment in CPython is worth keeping: the bit length is of `n` and not
    /// of `n - 1`, "because n can be 1". `n == 0` returns 0 without drawing.
    pub fn randbelow(&mut self, n: u32) -> u32 {
        if n == 0 {
            return 0;
        }
        let k = 32 - n.leading_zeros();
        loop {
            let r = self.getrandbits(k);
            if r < n {
                return r;
            }
        }
    }
}

/// `setsize` from `random.sample`, which decides the branch.
///
/// Computed in f64 exactly as CPython does: `math.log(k * 3, 4)` is
/// `log(k*3) / log(4)`, then `math.ceil`, then `4 ** that`. Reproducing the
/// float arithmetic rather than an integer equivalent is deliberate — the two
/// agree over the reachable range, but only one of them is what the reference
/// runs, and `tests/pyrandom_oracle.rs` checks it over every k the CLI accepts.
pub fn setsize_for(k: usize) -> f64 {
    let mut setsize = 21f64;
    if k > 5 {
        let e = ((k as f64) * 3.0).ln() / 4f64.ln();
        setsize += 4f64.powf(e.ceil());
    }
    setsize
}

/// `random.sample(range(n), k)` -- the indices, in the order `sample` returns
/// them.
///
/// The caller sorts them; the reference does
/// `sorted(random.sample(range(len(read_array)), args.sample_size))`. The order
/// is reproduced anyway, because a future caller that did not sort would
/// otherwise be silently wrong.
pub fn sample_indices(rng: &mut PyRandom, n: usize, k: usize) -> Vec<usize> {
    assert!(k <= n, "sample larger than population");
    let mut result = vec![0usize; k];
    if (n as f64) <= setsize_for(k) {
        // The pool branch: a partial Fisher-Yates over a copy, moving the
        // non-selected item at the end of the live range into the vacancy.
        let mut pool: Vec<usize> = (0..n).collect();
        for (i, slot) in result.iter_mut().enumerate() {
            let j = rng.randbelow((n - i) as u32) as usize;
            *slot = pool[j];
            pool[j] = pool[n - i - 1];
        }
    } else {
        // The selection-set branch: draw and retry until the index is new.
        let mut selected: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for slot in result.iter_mut() {
            let mut j = rng.randbelow(n as u32) as usize;
            while selected.contains(&j) {
                j = rng.randbelow(n as u32) as usize;
            }
            selected.insert(j);
            *slot = j;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first outputs of a freshly seeded generator, **recorded from
    /// CPython** rather than recalled:
    ///
    /// ```python
    /// r = random.Random(0); [r.getrandbits(32) for _ in range(4)]
    /// ```
    ///
    /// The first version of this test had plausible-looking wrong numbers in it
    /// -- they were the `random()` stream, not the `getrandbits(32)` one. Run
    /// the interpreter; do not remember it.
    #[test]
    fn the_bit_stream_matches_cpython() {
        for (seed, want) in [
            (0i64, [3626764237u32, 1654615998, 3255389356, 3823568514]),
            (1, [577090037, 2444712010, 3639700191, 3445702192]),
            (42, [2746317213, 478163327, 107420369, 3184935163]),
        ] {
            let mut r = PyRandom::seeded(seed);
            let got: Vec<u32> = (0..4).map(|_| r.getrandbits(32)).collect();
            assert_eq!(got, want.to_vec(), "seed {seed}");
        }
    }

    /// CPython seeds from the seed's ABSOLUTE value, so the sign is not a second
    /// axis. Surprising enough to pin.
    #[test]
    fn a_negative_seed_is_its_positive_twin() {
        let a: Vec<u32> = {
            let mut r = PyRandom::seeded(5);
            (0..4).map(|_| r.getrandbits(32)).collect()
        };
        let b: Vec<u32> = {
            let mut r = PyRandom::seeded(-5);
            (0..4).map(|_| r.getrandbits(32)).collect()
        };
        assert_eq!(a, b);
        // ...and i64::MIN must not panic on the way in.
        let _ = PyRandom::seeded(i64::MIN);
    }

    /// A seed needing two 32-bit key words takes a different path through
    /// `init_by_array` than a one-word seed.
    #[test]
    fn a_seed_above_2_32_uses_two_key_words() {
        let mut a = PyRandom::seeded(1);
        let mut b = PyRandom::seeded(1 + (1i64 << 32));
        assert_ne!(a.getrandbits(32), b.getrandbits(32));
    }

    #[test]
    fn getrandbits_zero_consumes_nothing() {
        let mut r = PyRandom::seeded(0);
        assert_eq!(r.getrandbits(0), 0);
        assert_eq!(r.getrandbits(32), 3626764237, "the stream did not advance");
    }

    #[test]
    fn setsize_matches_the_reference_at_the_matrix_sizes() {
        // Measured in Python: 4 ** ceil(log(k*3, 4)) + 21.
        assert_eq!(setsize_for(0), 21.0);
        assert_eq!(setsize_for(5), 21.0, "the k > 5 guard");
        assert_eq!(setsize_for(6), 21.0 + 64.0);
        assert_eq!(setsize_for(100), 1045.0);
        assert_eq!(setsize_for(200), 1045.0);
        assert_eq!(setsize_for(500), 4117.0);
    }

    /// Both corpora, at the matrix's sample sizes, and which branch each takes.
    /// This is the fact that makes the two committed corpora cover both.
    #[test]
    fn the_two_corpora_straddle_the_branch_threshold() {
        for k in [100usize, 200] {
            assert!(
                274f64 <= setsize_for(k),
                "smoke (274 reads) must take the pool branch at k={k}"
            );
            assert!(
                3000f64 > setsize_for(k),
                "sup (3000 reads) must take the selection-set branch at k={k}"
            );
        }
    }

    /// Recorded from CPython:
    /// `sorted(random.Random(SEED).sample(range(N), K))`.
    #[test]
    fn sample_matches_cpython_on_both_branches() {
        // pool branch: n <= setsize
        let mut r = PyRandom::seeded(0);
        let mut got = sample_indices(&mut r, 274, 100);
        got.sort_unstable();
        assert_eq!(got.len(), 100);
        assert_eq!(&got[..5], &[0, 3, 15, 16, 18], "pool branch, first five");
        assert_eq!(
            &got[95..],
            &[262, 263, 269, 271, 272],
            "pool branch, last five"
        );

        // selection-set branch: n > setsize
        let mut r = PyRandom::seeded(0);
        let mut got = sample_indices(&mut r, 3000, 100);
        got.sort_unstable();
        assert_eq!(got.len(), 100);
        assert_eq!(
            &got[..5],
            &[4, 57, 135, 165, 255],
            "selection-set branch, first five"
        );
    }

    #[test]
    fn sample_returns_distinct_indices_in_range() {
        for (n, k, seed) in [(274usize, 100usize, 0i64), (3000, 200, 7), (50, 50, 3)] {
            let mut r = PyRandom::seeded(seed);
            let got = sample_indices(&mut r, n, k);
            assert_eq!(got.len(), k);
            assert!(got.iter().all(|&i| i < n));
            let uniq: std::collections::HashSet<usize> = got.iter().copied().collect();
            assert_eq!(uniq.len(), k, "sample must not repeat an index");
        }
    }
}
