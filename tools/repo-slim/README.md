# tools/repo-slim — shrinking the repository

The repo is 491.82 MiB to clone. **98.0% of that is committed test data** under `test/`, plus 1.13 MB
of committed CPython bytecode and two `.DS_Store` files. All source, scripts, docs and the fixtures
that are actually used come to 10.68 MB.

Measured on this history by `analyze.sh`:

| | |
| --- | --- |
| total distinct blob bytes in history | 0.531 GB |
| stripped by these tools | 0.520 GB (**98.0%**) |
| retained source, scripts, docs | 0.005 GB |
| retained fixtures (the KEEP list) | 5.815 MB |
| paths stripped | **15** |
| commits | 186 |
| tags | 4, and **all four** have trees containing stripped paths |

## The three steps

Run them in this order. Each stops before doing anything irreversible.

```bash
tools/repo-slim/archive_data.sh     # 0. get the data out of git, into an archive
tools/repo-slim/analyze.sh          # 1. compute + review the removal list
tools/repo-slim/slim.sh             # 2. rewrite history in a scratch clone, verify
```

Only after all three pass do you force-push, and `slim.sh` prints the exact command. **Nothing here
pushes or uploads on its own.**

| Script | Does | Never does |
| --- | --- | --- |
| `archive_data.sh` | Extracts every `test/` blob that is **not** in `HEAD` from **history**, reassembles the split `ccs.fastq.gz`, writes a checksummed manifest | Upload, unless you pass `--upload` |
| `analyze.sh` | Writes `removal-paths.txt` and `analysis.txt` | Modify the repository at all |
| `slim.sh` | Clones a scratch mirror, rewrites it, verifies it | Touch your working repo, or push |

`removal-paths.txt` is generated, reviewable, and the single input to the rewrite. Read it first. It
comes to 15 paths.

## Step 0 is probably already done, in another repository

NGSpeciesID was forked from isONclust and inherited its test data. **Every strippable blob here is
already archived**, byte-identical, at `isONclust/repo-slim-archive/raw/test/` from that repository's
slimming exercise. Verified by extracting each blob from this history and comparing sha256:

| path | archived? |
| --- | --- |
| `test/ccs.fastq.gz.part-{aa,ab,ac,ad}` (364 672 972 bytes reassembled) | **identical** to `raw/test/ccs.fastq.gz` |
| `test/ENS_100k.fastq.tar.gz` | **identical** |
| `test/old_sorted_ens_100k.fastq.tar.gz` | **identical** |
| `test/sample_alz_2k.fastq` | **identical** |
| `test/isonclust1.out` | **identical** |
| `test/terncytb_200.fastq` | **NOT archived** — 214 022 bytes, NGSpeciesID-only |

So `archive_data.sh` has exactly one file to rescue, and it is 214 KB. Check before spending 500 MB
on a second copy of the rest:

```bash
ls -la /Users/*/source/isONclust/repo-slim-archive/raw/test/
```

`test/terncytb_200.fastq` is worth a look rather than a blind archive. It is **200 real ONT reads of a
different amplicon** — Cytb/16S, median length 428 against `sample_h1`'s 632 and
`Supplementary_File1_reads.fastq`'s 816 — added in `a2128c8` and dropped in `11c516f` ("changed test
data set"). PORTING.md's *The corpora* names a thin corpus registry as a gap, and this is a third real
corpus for 214 KB. **Restoring it to `test/` after the rewrite is a live option and probably a better
one than archiving it.** It does not discriminate `--symmetric_map_align_thresholds` (2 clusters
either way), so it is not a substitute for `Supplementary_File1_reads.fastq`.

## What is different here from isONclust, whose scripts these are

Four things, and each needed the scripts changed rather than just re-run.

**1. `KEEP` is not empty.** isONclust had already deleted every `test/` file from `HEAD`, so it
protected nothing. NGSpeciesID has **five live fixtures** the README links by path on the default
branch: `sample_h1.fastq` (the install check, also run by `.travis.yml`),
`Supplementary_File1_reads.fastq` (the paper's supplementary data and, per PORTING.md *Finding 19*,
the corpus the port develops against), `Supplementary_File2_minibar.txt`,
`Supplementary_File3_primer.txt` (the only committed input that exercises `barcode_trimmer`'s IUPAC
path) and `consensus.sh`.

**2. Three of those five share a blob with a renamed historical path**, so `KEEP` names both
spellings of each — eight entries for five files. Removing either spelling removes the blob, and the
live file with it:

```
test/Supplementary_File1_reads.fastq == test/Supplementary_File2_reads.fastq
test/Supplementary_File2_minibar.txt == test/Supplementary_File3_minibar.txt
test/Supplementary_File3_primer.txt  == test/Supplementary_File4_primer.txt
```

**3. The verification could no longer assert "nothing under `test/` survives".** Eight `test/` paths
are meant to. `slim.sh` now compares the surviving path set against `removal-paths.txt` directly with
`comm -12`, for both the history walk and each tag's tree: a path on the removal list that is still
reachable is a failure, and nothing else is.

**4. `analyze.sh` had a `set -e` bug that only a non-empty `KEEP` could reach.** The line was

```bash
[[ ${#KEEP[@]} -eq 0 ]] && echo "  (none -- ...)"
```

which yields the condition's exit status when it is *false*, aborting the script under `set -e`
before the fixture verification ran. With `KEEP` empty — isONclust's case — the condition was always
true and the bug was invisible; it fired on the first run here, exiting 1 after printing a header and
no fixtures. Rewritten as `if/fi`. This is PORTING.md's method rule arriving on the tooling itself:
**a harness that has never failed has not been tested, it has only been run.**

## Requirements

`git-filter-repo` — `pipx install git-filter-repo` or `pip install git-filter-repo`. Point at an
unusual install with `GIT_FILTER_REPO=/path/to/git-filter-repo`.

macOS ships bash 3.2, which has no `mapfile` and errors on empty-array expansion under `set -u`.
These scripts are written for it. Do not add `mapfile`.

## The deduplication trap, and why it DOES bite here

`git rev-list --objects` emits each blob exactly *once*, with a single path, so any file whose content
is byte-identical to another is invisible in that listing. Building the removal list from it silently
missed a file in isONcorrect: the rewrite ran, reported success, and left an 18 MB file in place.

isONclust checked and found only three such blobs, none of which mattered. **Here there are seven,
and one of them matters:**

| blob | paths |
| --- | --- |
| `d78cdf85` | `test/Supplementary_File1_reads.fastq`, `test/Supplementary_File2_reads.fastq` |
| `5008ddfc` | `modules/.DS_Store`, **`scripts/.DS_Store`** |
| `94037afc` | `NGSpeciesID`, `isONclust` |
| `c964fb0d` | `cemetary/cluster_parallel.py`, `modules/cluster_parallel.py` |
| `d231d147` | `modules/compute_shared_minimizers_probabilities.py`, `scripts/compute_shared_minimizer_probabilities.py` |
| `1434c00b` | `test/Supplementary_File3_primer.txt`, `test/Supplementary_File4_primer.txt` |
| `ef1deb95` | `test/Supplementary_File2_minibar.txt`, `test/Supplementary_File3_minibar.txt` |

`scripts/.DS_Store` would **not** have appeared on a removal list built from the object listing, and
it is tracked in `HEAD`. The list is built from a tree walk (`git log --all --name-only`), because
that is the only way to know any of the above.

The same trap applies to `git rev-parse HEAD:<missing-path>`, which prints the unresolved string to
stdout *and* exits non-zero. Existence checks use `git cat-file -e`.

## Six files are tracked in `HEAD` and should be untracked first

```
modules/.DS_Store
scripts/.DS_Store
modules/__pycache__/__init__.cpython-36.pyc
modules/__pycache__/cluster.cpython-36.pyc
modules/__pycache__/get_sorted_fastq_for_cluster.cpython-36.pyc
modules/__pycache__/p_minimizers_shared.cpython-36.pyc
```

1.13 MB, and `.gitignore` already lists `*.pyc` and `modules/__pycache__/*` — which does nothing for
files that are already tracked. `git rm --cached` them in an **ordinary commit before** the rewrite,
so the untracking and the history rewrite are separately reviewable. The rewrite would remove them
either way; doing it first means the diff a reviewer sees is one commit removing six files, not a
force-push.

## What the rewrite costs

- **Every commit SHA changes.** All 16 forks diverge permanently and cannot be fast-forwarded; anyone
  holding a clone (68 stargazers, 16 forks) must re-clone.
- **All four tags move** to rewritten commits — `0.0.4`, `v0.1.2.1`, `v0.3.0`, `v0.3.1`. `slim.sh`
  asserts they all survive and that every tree is clean.
- **Commit-pinned links break**, including any in the paper (10.1002/ece3.7146) and in the
  field-protocol manuscript the README's EXAMPLE WORKFLOW describes.
- **The README's four `test/` links keep working**, because those paths are on the `KEEP` list and
  `slim.sh` asserts their blobs unchanged. Re-check them after the push anyway.
- **GitHub keeps old objects reachable for a while.** Stripped data may stay downloadable via old
  SHAs until GitHub garbage-collects; ask GitHub Support to run `gc` if it matters. Measured on
  isONclust immediately after its push: `git fetch` of an old SHA was refused, but the GitHub API
  still served the commit *and* the blob bytes, and the reported repository size did not change.
  Forks in the same network legitimately still hold the data on their live branches, and NGSpeciesID
  has 16. This is public test data behind a published paper, so it is a tidiness question rather
  than a disclosure one — and **the thing that actually matters is achieved either way: the
  repository clones in single-digit megabytes instead of 500.**

Mitigation for the backup: keep the old history **locally**, not as a tag on origin. A tag pointing
at old history keeps every stripped byte alive on the server and gets fetched by `git clone`.

## Verification

`slim.sh` refuses to declare success unless all of these hold:

- commits touching source preserved exactly (the total legitimately falls as data-only commits become
  empty and are pruned)
- 24 named source files, packaging files and live fixtures have byte-identical blobs before and after
- all six junk paths tracked in `HEAD` are absent from the rewritten `HEAD`
- **no** path from `removal-paths.txt` survives anywhere in history
- the tag set is preserved and every one of the four tags' trees carries no removal-list path

## Afterwards

Three things go back on top of the pushed history, and — unlike isONclust — **no fixture needs
restoring**, because every fixture the README and `.travis.yml` reference is on the `KEEP` list.

1. A `.gitignore` for the junk the rewrite removed, so it cannot come back:
   `printf '__pycache__/\n*.py[cod]\n.DS_Store\n' >> .gitignore`
2. `PORTING.md`, `bench/` and `tools/repo-slim/` itself, which is where the port starts.
3. Optionally `test/terncytb_200.fastq` — 214 KB, a third real corpus, see *Step 0* above.
