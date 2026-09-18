//! A minimal safe wrapper over `spoa-sys`, replacing the `spoa` convenience crate.
//!
//! # Why this exists
//!
//! `spoa 0.2.1` does not compile on **aarch64 Linux**. Its wrapper casts
//! sequence pointers with `as *const i8`, while the `spoa-sys` bridge correctly
//! declares them `*const c_char` — and `c_char` is `i8` on x86_64 and on all
//! Apple targets but **`u8` on aarch64 Linux**. Three `E0308: mismatched types`
//! in a dependency, on one of the four platforms this port claims to support.
//!
//! CI found it; nothing local could have. Every machine this was developed on
//! is one of the targets where `c_char == i8`, so the bug is invisible until
//! something builds it on ARM Linux — which is Graviton, most CI ARM runners,
//! and every Raspberry Pi.
//!
//! The upstream wrapper is about sixty lines of pointer casts over `spoa-sys`,
//! and this port uses six of its functions. Rather than pin a broken version,
//! fork it, or drop a platform, the six are written here with `c_char` casts
//! that are correct everywhere. `spoa-sys` — which is the actual C++ bridge, the
//! part that would be painful to reproduce — is used unchanged.
//!
//! Worth reporting upstream; it is a one-line fix in `spoa/src/lib.rs`.
//!
//! # What is NOT wrapped
//!
//! `consensus_with_min_coverage`, `multiple_sequence_alignment` and `clear`.
//! The reference calls none of them, and an unused wrapper is an unverified
//! wrapper.

use std::ffi::c_char;

pub use spoa_sys::ffi::AlignmentType;

/// An alignment of one sequence to a `Graph`.
pub struct Alignment(cxx::UniquePtr<spoa_sys::ffi::Alignment>);

/// Aligns sequences to a `Graph`.
pub struct AlignmentEngine(cxx::UniquePtr<spoa_sys::ffi::AlignmentEngine>);

impl AlignmentEngine {
    /// `spoa::AlignmentEngine::Create(type, m, n, g, e, q, c)`.
    ///
    /// The reference builds this as `kSW, 5, -4, -2, -6, -10, -4` — local
    /// alignment, match +5, mismatch -4, and a linear gap of -2 because spoa
    /// treats `g == e` as linear. See `poa.rs`.
    pub fn new(typ: AlignmentType, m: i8, n: i8, g: i8, e: i8, q: i8, c: i8) -> Self {
        AlignmentEngine(spoa_sys::ffi::create_alignment_engine(
            typ, m, n, g, e, q, c,
        ))
    }

    pub fn align(&mut self, sequence: &[u8], graph: &Graph) -> Alignment {
        let len = u32::try_from(sequence.len()).expect("a sequence longer than 4 GB");
        Alignment(unsafe {
            spoa_sys::ffi::align(
                self.0.pin_mut(),
                // `c_char`, NOT `i8`. This cast is the entire reason this module
                // exists; see the module docs.
                sequence.as_ptr() as *const c_char,
                len,
                self_graph(graph),
            )
        })
    }
}

/// The partial order alignment graph.
pub struct Graph(cxx::UniquePtr<spoa_sys::ffi::Graph>);

fn self_graph(g: &Graph) -> &spoa_sys::ffi::Graph {
    g.0.as_ref().expect("spoa returned a null graph")
}

impl Graph {
    pub fn new() -> Self {
        Graph(spoa_sys::ffi::create_graph())
    }

    /// Add a sequence at uniform `weight`. The reference uses weight 1 for FASTA
    /// input and for any record whose quality string is missing or the wrong
    /// length.
    pub fn add_alignment(&mut self, alignment: &Alignment, sequence: &[u8], weight: u32) {
        let len = u32::try_from(sequence.len()).expect("a sequence longer than 4 GB");
        unsafe {
            spoa_sys::ffi::add_alignment(
                self.0.pin_mut(),
                alignment
                    .0
                    .as_ref()
                    .expect("spoa returned a null alignment"),
                sequence.as_ptr() as *const c_char,
                len,
                weight,
            )
        }
    }

    /// Add a sequence weighted per base by its FASTQ quality.
    ///
    /// This is not an optimisation: spoa weights each base by `ord(q) - 33`, so
    /// a FASTQ input produces a different consensus from the same sequences
    /// without qualities. Measured at 847 bp against 860 bp on one cluster.
    pub fn add_alignment_with_qual(&mut self, alignment: &Alignment, sequence: &[u8], qual: &[u8]) {
        let len = u32::try_from(sequence.len()).expect("a sequence longer than 4 GB");
        let qlen = u32::try_from(qual.len()).expect("a quality string longer than 4 GB");
        assert_eq!(len, qlen, "spoa requires one quality character per base");
        unsafe {
            spoa_sys::ffi::add_alignment_with_qual(
                self.0.pin_mut(),
                alignment
                    .0
                    .as_ref()
                    .expect("spoa returned a null alignment"),
                sequence.as_ptr() as *const c_char,
                len,
                qual.as_ptr() as *const c_char,
                qlen,
            )
        }
    }

    pub fn consensus(&mut self) -> Vec<u8> {
        Vec::from(
            spoa_sys::ffi::generate_consensus(self.0.pin_mut())
                .as_ref()
                .expect("spoa returned a null consensus")
                .as_bytes(),
        )
    }
}

impl Default for Graph {
    fn default() -> Self {
        Self::new()
    }
}
