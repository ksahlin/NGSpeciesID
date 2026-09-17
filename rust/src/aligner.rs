//! Which parasail implementation the port actually calls.
//!
//! There are two, they are exact with respect to each other, and the choice is
//! a build-time feature rather than a correctness one:
//!
//! * `parasail_ffi` — parasail's own SIMD C library, through `libparasail-sys`.
//!   This IS the implementation the Python reference calls, so it is exact by
//!   construction.
//! * `parasail` — a scalar Rust reimplementation, carried across from the
//!   isONclust port where it was measured identical to the C on 18 633
//!   alignments, CIGAR and ratio both.
//!
//! # Why this module exists at all
//!
//! The dispatcher used to live privately in `blockalign.rs`, which meant it
//! covered the clustering call site and **only** the clustering call site.
//! `consensus::identity` — the reverse-complement detection path, and the second
//! of the two parasail call sites, with its own `opening_penalty=3` — called
//! `parasail::semiglobal` directly and therefore stayed scalar no matter how the
//! feature was set. Not a correctness bug, since the two agree, but the feature
//! did not mean what its name said. One dispatcher, both call sites.
//!
//! # Speed, measured, and why the default is what it is
//!
//! On a private long-read corpus (5 000 reads, 1 000-1 858 bp), `--ont --t 1`:
//!
//! | | wall |
//! | --- | --- |
//! | the Python reference (which calls the C) | 10.3 s |
//! | this port, scalar | 24.9 s |
//! | this port, FFI | 4.1 s |
//!
//! A Rust port being 2.4x slower than the Python it replaces is not a tradeoff
//! anyone chose; it is what happens when one side delegates to vectorised C and
//! the other does not. The amplicon corpora in this repository's own `test/`
//! directory hide it completely — at 816 bp both finish in under two seconds —
//! which is exactly why the original default was wrong. See PORTING.md,
//! *Performance*.

use crate::parasail::{Alignment, Scoring};

/// Semi-global affine alignment, `parasail.sg_trace_scan_16`'s contract.
///
/// Every call site in the port goes through here. Adding another that calls
/// `parasail::semiglobal` or `parasail_ffi::semiglobal` directly re-creates the
/// bug in this module's docs.
#[inline]
pub fn semiglobal(s1: &[u8], s2: &[u8], sc: Scoring) -> Alignment {
    #[cfg(feature = "parasail-ffi")]
    {
        crate::parasail_ffi::semiglobal(s1, s2, sc)
    }
    #[cfg(not(feature = "parasail-ffi"))]
    {
        crate::parasail::semiglobal(s1, s2, sc)
    }
}
