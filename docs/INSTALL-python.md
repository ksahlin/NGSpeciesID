Installing the Python implementation
====================================

The Python implementation is the **reference**: the Rust one is checked against it, case by case, on
every commit. It is still supported and still what `pip install NGSpeciesID` gives you.

For most users the [Rust binary](../README.md#binaries) is the easier install — one download, no
dependency resolution. Use this if you want the reference implementation itself, or if you are
working on it.

```
conda create -n NGSpeciesID -c conda-forge -c bioconda python=3.12 pip medaka spoa racon minimap2 samtools
conda activate NGSpeciesID
pip install --no-deps NGSpeciesID
```

Then [test the installation](../README.md#testing-installation). You need `conda activate NGSpeciesID`
in each new shell.

## Which version you actually get

| where | version |
| --- | --- |
| PyPI (`pip install NGSpeciesID`) | **0.4.1** |
| bioconda (`conda install ngspeciesid`) | 0.3.1 until [the recipe update](https://github.com/bioconda/bioconda-recipes/pull/69371) merges |
| this repository, and the binaries | **0.4.1** |

0.4.0 exists as a git tag with working binaries, but it cannot be published to PyPI: its `setup.py`
omits `long_description_content_type`, so the project page would not render. 0.4.1 is that fix and is
byte-identical in behaviour.

## Why pip at all, when there is a conda package?

There is one — `bioconda::ngspeciesid` — and it is a version behind. The recipe above installs the
*dependencies* with conda and the *package* with pip, so you get the newest release without waiting
for the bioconda recipe to be updated. `conda install -c conda-forge -c bioconda ngspeciesid` works
and is simpler, if a slightly older version is fine.

## Why the commands look like that

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

