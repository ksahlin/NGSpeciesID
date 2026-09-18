Installation details
====================

The [README](../README.md) has the commands. This has the reasons, and the caveats that only matter
once something has gone wrong.

Contents
--------

  * [Why the conda command looks like that](#why-the-conda-command-looks-like-that)
  * [Building the Rust implementation](#building-the-rust-implementation)
  * [What needs external tools, and what does not](#what-needs-external-tools-and-what-does-not)
  * [Reproducibility across machines](#reproducibility-across-machines)
  * [Reproducibility of --sample_size](#reproducibility-of---sample_size)

## Why the conda command looks like that

Three details in the two lines above are deliberate, and getting any of them wrong is what makes the
install fail.

**`--no-deps` on the pip step.** NGSpeciesID needs two python libraries, `parasail` and `edlib`.
The conda command above already installs both — `medaka` depends on `parasail-python` and
`python-edlib`, and bioconda has prebuilt packages of each for linux-64, osx-64 and osx-arm64.
Without `--no-deps`, pip reads this package's `install_requires` and reinstalls `parasail` **from
PyPI, on top of the working conda build**. PyPI publishes no `parasail` wheel for any ARM platform,
so on Apple Silicon (and on ARM Linux) pip falls back to compiling it from source, which fails after
about two minutes with `RuntimeError: autoreconf -fi failed`. On x86_64 a PyPI wheel exists, the
reinstall succeeds, and you never notice — which is why this went unreported for so long.

**No version pins.** Earlier versions of these instructions pinned `medaka==2.0.1` (or `==0.11.5`)
and `openblas==0.3.3`. Those pins resolve on linux-64 and on neither macOS platform: `medaka 2.0.1`
has no macOS build at all, and conda-forge's earliest `openblas` for osx-arm64 is 0.3.11. Leave them
unpinned.

**`python=3.12`.** `medaka` 2.2.x requires it. It is also the interpreter to prefer for reproducible
results: on python 3.11 and earlier, `sum()` over a set of floats is order-dependent, and repeated
runs of NGSpeciesID on the same input can differ in the last digits of the reported read error rates.

If you would rather not use `medaka` at all, install the libraries directly and use `--racon`:

```
conda create -n NGSpeciesID -c conda-forge -c bioconda python=3.12 pip parasail-python python-edlib spoa racon minimap2
conda activate NGSpeciesID
pip install --no-deps NGSpeciesID
```


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

