Packaging
=========

`bioconda-meta.yaml` is the updated bioconda recipe, ready to submit. It is kept here so the
checksum and the reasoning live with the release that produced them.

**It targets 0.4.1, not 0.4.0.** 0.4.0 is a valid git tag with valid binaries, but its `setup.py`
lacks `long_description_content_type`, so it cannot be published to PyPI and there is no reason to
put it on bioconda either. 0.4.1 is byte-identical in behaviour.

## Submitting it

```
git clone https://github.com/<you>/bioconda-recipes    # your fork
cd bioconda-recipes
git checkout -b ngspeciesid-0.4.1
cp <this repo>/docs/packaging/bioconda-meta.yaml recipes/ngspeciesid/meta.yaml
git commit -am "Update ngspeciesid to 0.4.1"
git push origin ngspeciesid-0.4.1
```

Then open a PR against `bioconda/bioconda-recipes`. Their CI builds it; a maintainer merges.

## What changed from the published 0.3.1 recipe, and why

Four lines.

| line | 0.3.1 | 0.4.0 | why |
| --- | --- | --- | --- |
| `version` | 0.3.1 | 0.4.1 | |
| `sha256` | `ddd378a6…` | `0800e5e3…` | of `v0.4.1.tar.gz`, computed from the downloaded tag archive |
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

Publishing to PyPI
------------------

`.github/workflows/publish-pypi.yml` builds, checks and uploads on a `v*` tag, using **trusted
publishing** — PyPI trusts this repository and workflow directly, and the upload is authenticated by
a short-lived token GitHub mints for the job. No API token to create, store, rotate or leak.

### One-time setup, which only you can do

1. Sign in to https://pypi.org as an owner of the `NGSpeciesID` project.
2. Go to **Manage → Publishing → Add a new publisher**, and choose GitHub.
3. Fill in exactly:

   | field | value |
   | --- | --- |
   | Owner | `ksahlin` |
   | Repository name | `NGSpeciesID` |
   | Workflow name | `publish-pypi.yml` |
   | Environment name | `pypi` |

4. Save.

### Publishing 0.4.0

The `v0.4.0` tag was pushed before this workflow existed, so it will not fire on its own. After the
setup above, run it by hand:

```
gh workflow run publish-pypi.yml -f tag=v0.4.0
```

Every later release publishes automatically when you push the tag.

### What the workflow checks before uploading

* `twine check` on both artefacts — a bad `long_description` render is only visible after publishing,
  and **a version cannot be re-uploaded to PyPI**;
* the built wheel is installed and `NGSpeciesID --version` is run from it, because `setup.py` ships
  `scripts=['NGSpeciesID']` and an entry point that does not land in `bin/` would still build fine.

### A packaging bug fixed on the way

`setup.py` set `long_description` to the markdown README and never set
`long_description_content_type`, so PyPI would have rendered it as reStructuredText and mangled the
project page. It also carried the sample-project's commented-out MIT classifier while the project is
GPL-3.0-or-later. Both fixed.

That fix is **not** in the `v0.4.0` tag archive, which was created earlier; it changes packaging
metadata only, no behaviour. The bioconda recipe builds from the tag and is unaffected. If you would
rather have them identical, publish this as 0.4.1 instead — nothing depends on 0.4.0 existing on PyPI
yet.
