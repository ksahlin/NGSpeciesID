Changelog
=========

## 0.4.1

Packaging only. No change to behaviour or output; 0.4.1 is byte-identical to 0.4.0.

`setup.py` set `long_description` to the markdown README and never set
`long_description_content_type`, so PyPI would have parsed it as reStructuredText and mangled the
project page. The publish workflow's `twine check` caught it and refused to upload 0.4.0 — a version
cannot be replaced on PyPI once published, so this is released as 0.4.1 rather than retagged.

The licence classifier was also still pypa's sample-project MIT placeholder; the project is
GPL-3.0-or-later.

Linux binaries are now built with `cargo-zigbuild` against the glibc 2.17 ABI, with libc++ linked
statically. The 0.4.0 binaries required GLIBC_2.39 — Ubuntu 24.04 and nothing older — so they did not
run on CentOS 7, RHEL/Rocky 8, Ubuntu 20.04/22.04 or Debian 12. The workflow now asserts both the
glibc floor and the absence of a libstdc++ dependency, so a regression fails the release instead of
shipping.

## 0.4.0

### Rust implementation

NGSpeciesID is re-implemented in Rust. It produces **byte-identical output** to the Python
implementation — that is its specification, not an aspiration — and is 3-6x faster at clustering,
6.6x at `--t 8`.

Both implementations are in this repository and both are supported. `pip install NGSpeciesID`
installs the Python one; prebuilt Rust binaries are attached to this release.

How byte identity is checked, on every commit:

* 55 output cases x 6 corpora = 330 cases, plus 45 CLI cases x 6
* stage-level oracles replaying recorded reference calls: parasail alignments, edlib HW (96 calls),
  alignment identity (12 calls x 3 values), CPython's Mersenne Twister (64 draws), `str(float)`,
  `round(x, 2)`
* CI on Linux and macOS, x86_64 and arm64, including a full run against a conda-built reference

`spoa` is linked as a library rather than run as a subprocess, so the draft consensus needs no
external tools. `racon`, `minimap2` and `medaka` are still external programs.

### `--seed`, and `--sample_size` is reproducible

`--sample_size` drew its subsample from an unseeded global RNG, so two runs of the same command gave
different answers. Measured on the paper's own supplementary data: five runs, five results.

`--seed` (default 0) now seeds the draw. **This changes results** for existing `--sample_size` users:
the same command now returns the same subsample every time, where before it did not. `--top_reads`
takes the highest-scoring reads instead of a random subset and ignores `--seed`.

### `--consensus` without a polisher is draft-only

`--consensus` with neither `--medaka` nor `--racon` crashed with `UnboundLocalError` — after running
clustering, spoa and reverse-complement detection, and writing nothing. It now writes the spoa draft
consensus and exits 0. This is the only consensus mode that needs no external tools.

Conversely, `--medaka` or `--racon` **without** `--consensus` was silently accepted and did nothing.
It is now an error, so a run that polishes nothing says so.

### Installation

The published conda recipes did not resolve. All three failed on `osx-arm64` and two also on
`osx-64`, because of version pins (`medaka==2.0.1`, `openblas==0.3.3`) that exist on no macOS
platform, and because `pip install` without `--no-deps` reinstalls `parasail` from PyPI over the
working conda build — and PyPI has no ARM wheel, so it tries to compile it and fails.

The README recipe is fixed and verified end to end from a clean environment.

Python 3.12 or newer is recommended: before 3.12, `sum()` over a set of floats is order-dependent,
and repeated runs could differ in the last digits of reported read error rates.

### Known limitation

Given identical input and identical versions, `racon` produces different output on x86_64 Linux than
on arm64 macOS, so a polished consensus is reproducible on a platform but not across architectures.
This is racon's behaviour and was already true of earlier versions. Everything upstream of racon is
identical across all four platforms tested.

---

0.3.2 was prepared but never released; its changes are included above.
