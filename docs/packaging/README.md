Packaging
=========

`bioconda-meta.yaml` is the updated bioconda recipe, ready to submit. It is kept here so the
checksum and the reasoning live with the release that produced them.

## Submitting it

```
git clone https://github.com/<you>/bioconda-recipes    # your fork
cd bioconda-recipes
git checkout -b ngspeciesid-0.4.0
cp <this repo>/docs/packaging/bioconda-meta.yaml recipes/ngspeciesid/meta.yaml
git commit -am "Update ngspeciesid to 0.4.0"
git push origin ngspeciesid-0.4.0
```

Then open a PR against `bioconda/bioconda-recipes`. Their CI builds it; a maintainer merges.

## What changed from the published 0.3.1 recipe, and why

Four lines.

| line | 0.3.1 | 0.4.0 | why |
| --- | --- | --- | --- |
| `version` | 0.3.1 | 0.4.0 | |
| `sha256` | `ddd378a6…` | `2fc90547…` | of `v0.4.0.tar.gz`, computed from the downloaded tag archive |
| `host: python` | >=3.10 | **>=3.12** | |
| `run: python` | >=3.10 | **>=3.12** | |

The python floor is the only judgement call. `medaka` 2.2.x requires 3.12 anyway, so in practice the
solve already pulled it — but more importantly, **before CPython 3.12 `sum()` over a set of floats is
order-dependent**, and this tool sums over a set of quality characters in four places. On 3.11 and
earlier, two runs on the same input can differ in the last digits of the reported read error rates.
Stating 3.12 makes the package reproducible rather than leaving it to the solver. See PORTING.md,
*Finding 2*.

`build: number` resets to 0, which is correct for a new version.

Nothing else moves: still `noarch: python`, same run dependencies, same tests.

## Verified before writing

* the tag archive downloads and its sha256 is the one above;
* the recipe's own build line — `pip install . -vv --no-deps --no-build-isolation` — succeeds on that
  archive;
* both of the recipe's test commands pass against the result: `NGSpeciesID --help`, and
  `NGSpeciesID --version` printing `NGSpeciesID 0.4.0`.

That last check matters more than it looks. The version string lives in four places in this
repository and `--version` is the one the recipe asserts; a bump that missed one would pass the build
and fail the test. CI's `version` job compares all four plus the recorded golden.

## What this recipe does NOT do

It packages the **Python** implementation, which is what the 0.3.1 recipe packaged. The Rust binary
is distributed through GitHub releases for now.

A recipe for the Rust binary would be a bigger change and a better package: it drops `python`,
`edlib`, `parasail-python` and `spoa` from the runtime dependencies entirely — spoa is linked in —
leaving only `racon`, `minimap2`, `medaka` and `samtools`, and those only for polishing. It stops
being `noarch`, so it needs per-platform builds; bioconda does linux-64 and osx-64 by default, with
linux-aarch64 and osx-arm64 opt-in. Worth doing after the Rust implementation has been in use for a
release or two.
