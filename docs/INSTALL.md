Installation details
====================

The [README](../README.md) has the commands. This has the reasons, and the caveats that only matter
once something has gone wrong. Installing the **Python** implementation is a separate page:
[INSTALL-python.md](INSTALL-python.md).

Contents
--------

  * [Building the Rust implementation](#building-the-rust-implementation)
  * [What the Linux binaries require](#what-the-linux-binaries-require)
  * [macOS Gatekeeper](#macos-gatekeeper)
  * [What needs external tools, and what does not](#what-needs-external-tools-and-what-does-not)
  * [Reproducibility across machines](#reproducibility-across-machines)
  * [Reproducibility of --sample_size](#reproducibility-of---sample_size)

## Building the Rust implementation

```
cargo build --release --manifest-path rust/Cargo.toml
```

Needs **Rust 1.88 or newer**. Because it links parasail's C library through FFI it also needs:

| | Debian/Ubuntu | macOS |
| --- | --- | --- |
| cmake, libclang, pkg-config | `apt install cmake libclang-dev pkg-config` | included with the Xcode command line tools |

If you would rather not install those:

```
cargo build --release --manifest-path rust/Cargo.toml --no-default-features
```

That uses a pure-Rust aligner instead. It needs nothing but a Rust toolchain and produces the same
bytes; it is slower on reads much longer than ~1 kb, where the C library's vector instructions
matter.

The binaries attached to each release are built by CI on Linux and macOS, x86_64 and arm64, and are
byte-identity-checked against the Python on every commit.

## What the Linux binaries require

Built with [cargo-zigbuild](https://github.com/rust-cross/cargo-zigbuild) against an old glibc ABI,
with libc++ linked statically:

```
highest glibc symbol required: GLIBC_2.17
libstdc++:                     none
```

glibc 2.17 is CentOS 7 (2014), so the binaries run there and on everything newer. The macOS binaries
link only `libc++` and `libSystem` from the OS.

This is measured on every release, and **asserted** in the workflow: if a build ever needs a newer
glibc, or picks up a libstdc++ dependency, the release fails rather than shipping a binary much of
the intended audience cannot execute. The first 0.4.0 binaries needed GLIBC_2.39 — Ubuntu 24.04 and
nothing older — which is what prompted the change.

## macOS Gatekeeper

The macOS binaries carry only an ad-hoc signature, not an Apple Developer ID one, and are not
notarized. A file downloaded through a browser gets the `com.apple.quarantine` attribute, and
Gatekeeper then refuses to run unsigned code:

> Apple could not verify that "NGSpeciesID" is free of malware that may harm your Mac or compromise
> your privacy.

Clear the flag on the extracted folder:

```
xattr -dr com.apple.quarantine NGSpeciesID-*
```

`curl` and `gh release download` do not set the attribute, so a binary fetched that way runs without
this step — which is also why it is easy to miss when testing.

**The real fix is signing and notarizing**, which needs a paid Apple Developer account ($99/year) and
a `codesign --sign "Developer ID Application: ..."` plus `notarytool submit` step in the release
workflow. Until then this note is the workaround, and it belongs in the README rather than in a
support thread.

## What needs external tools, and what does not

`spoa` is **linked into** the Rust binary rather than run as a subprocess, so the draft consensus
needs nothing external. The polishers are still separate programs.

| what you run | external tools needed |
| --- | --- |
| clustering — `final_clusters.tsv`, `final_cluster_origins.tsv` | none |
| `write_fastq` | none |
| `--consensus` (draft) | none |
| `--consensus --racon` | `racon`, `minimap2` |
| `--consensus --medaka` | `medaka` |

## Reproducibility across machines

Given identical input and identical versions, **racon produces different output on x86_64 Linux than
on arm64 macOS**. So a polished consensus is reproducible on a given platform and is not guaranteed
to be identical if you move the same analysis to a different architecture.

This is a property of racon, not of NGSpeciesID, and it was already true of the Python. Everything
upstream of racon — clustering, the spoa draft consensus, primer trimming — is identical across all
four platforms tested.

## Reproducibility of --sample_size

Running the same command on the same input gives the same answer, down to the polished consensus
sequences. Two things are worth knowing about that:

* **`--sample_size` is seeded.** It draws a random subset of reads, and the draw comes from `--seed`
  (default 0), so it is reproducible. Pass a different `--seed` to draw a different subset — useful
  for checking how sensitive a consensus is to which reads went into it. This is the only randomness
  in the tool.

  Before v0.4.0 the draw was **not** seeded, so two runs of the same `--sample_size` command gave
  different clusters and different consensus sequences. If you are comparing against results
  produced by an older version, they will not match.

* **Use python 3.12 or newer.** On python 3.11 and earlier, `sum()` over a set of floats is
  order-dependent, so the read error rates NGSpeciesID reports can differ in their last digits
  between runs of the same command.

`--top_reads` is a different thing and still available: it takes the `--sample_size` highest-scoring
reads instead of a random subset, and ignores `--seed`.

