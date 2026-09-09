# bench/ — the equivalence harness

Byte-identity with the Python reference is the acceptance criterion for the Rust port. This is what
checks it. Nothing here is part of the shipped tool.

```bash
bench/setup_reference_env.sh          # build the pinned reference env
export PATH=~/miniforge3/envs/ngspeciesid-ref/bin:$PATH   # the external tools live here
bench/equivalence.sh env              # is it usable? which of the four tools are present?
bench/equivalence.sh seeds            # does the reference agree with itself?  <-- RUN FIRST
bench/equivalence.sh cli record       # capture the CLI contract
bench/equivalence.sh record           # record output goldens
bench/equivalence.sh stable           # record twice; the goldens must be identical
bench/equivalence.sh verify           # run the port, diff against the goldens
bench/equivalence.sh tools            # a missing binary must be named, not tracebacked
```

`CORPUS` picks the input (a path, or a name from `corpora.tsv`); `GOLDEN` picks where goldens live.
Goldens are per-corpus — including the CLI ones, which is not obvious and was checked rather than
assumed: recording the CLI contract against both corpora and diffing gives **16 of 44 cases
differing**, because every case that gets far enough to run prints its read and cluster counts
(`Starting Clustering: 274 reads` against `3000`). The other 28 are corpus-independent. Splitting
them would save about 180 KB and add a rule to remember, so they are simply recorded twice.

```bash
CORPUS=sup   GOLDEN=$PWD/bench/golden/sup   bench/equivalence.sh record
CORPUS=smoke GOLDEN=$PWD/bench/golden/smoke bench/equivalence.sh record
```

## What is here

| file | role |
| --- | --- |
| `equivalence.sh` | the harness. Every subcommand above |
| `cases.tsv` | 51 output cases. TAB-separated; see the warning at the top of the file |
| `corpora.tsv` | the corpus registry. Both entries are committed |
| `setup_reference_env.sh` | builds the reference environment. One conda line — see below |
| `dump_reference.py` | wraps the reference's own functions to dump a stage's inputs *and* outputs, for stages whose output never reaches a file |
| `diffsummary.py` | says which **column** moved and by how much. A line diff is useless on these files: `final_cluster_origins.tsv` carries the full read sequence and quality string in columns 3 and 4, so one wrong float in column 6 prints four kilobytes |
| `golden/<corpus>/manifest.tsv` | per-file sha256 for every case, plus the provenance the goldens are only valid under |
| `golden/<corpus>/cli/<case>/` | `exit`, `stdout`, `stderr` for 44 CLI cases |
| `golden/<corpus>/sample/` | the `default` case's output as **heads**, so there is something to read by eye without running anything. Heads and not whole files: `final_clusters.tsv` is one line per read, each carrying a whole ONT accession, which is 686 KB on the 3 000-read corpus. The hashes in `manifest.tsv` are the contract; this directory is a courtesy |
| `env/resolved-*.txt` | what the reference environment actually resolved to |

## The reference environment is one conda line

```
conda create -n ngspeciesid-ref -c conda-forge -c bioconda \
  python=3.12 pip parasail-python python-edlib spoa racon minimap2 samtools medaka
pip install --no-deps -e .
```

isONclust's equivalent script is three times longer, and almost all of it is there to build
`parasail` from source with a `$M4` export and two `glibtoolize` symlinks. None of that is needed:
bioconda ships `parasail-python` and `python-edlib` prebuilt for `osx-arm64`, and `medaka` depends on
both. `--no-deps` on the pip step is what stops pip reinstalling parasail from PyPI over the working
conda build — see PORTING.md, *Goal*.

## The goldens are hashes, not files

`manifest.tsv` records a sha256 per file per case. Recorded verbatim the 51 cases would be far larger
than this repository should carry, and hashes are enough to **fail** correctly. They cannot say
*what* moved, so on a mismatch `verify` re-runs the reference for that one case and diffs properly.

The manifest header records the corpus sha256, the interpreter version, whether `sum()` is
compensated, **and the version of every external binary**. That last part is not decoration: a `spoa`
upgrade changes every `consensus_reference_*.fasta` and a `racon` or `medaka` upgrade changes every
polished consensus, so the `--consensus` goldens are only valid for the four versions that produced
them.

## What is compared, and what is not

Every file a run writes, byte for byte — including the `--t > 1` merge intermediates and every
per-cluster consensus file. Three exceptions, each **measured** rather than assumed:

| excluded | why |
| --- | --- |
| `*_stderr_it_*.txt`, `mm2_stderr_it_*.txt`, `stdout.txt`, `stderr.txt` | racon's and minimap2's captured stderr carries timings. Three runs of `--consensus --racon` differed in exactly these and in nothing else |
| `calls_to_draft.bam`, `calls_to_draft.bam.bai` | minimap2 and samtools write their own command lines into the BAM `@PG` header, and those contain the **absolute path** of the output folder. Every case runs in a fresh temp dir, so the file can never match |
| — | `consensus_probs.hdf` is deliberately **not** excluded. It was measured stable. Excluding a file because it is the kind of file that might vary is how a contract quietly stops covering anything |

Their *existence* is still contract: `verify` compares the file count too, so a port that writes none
of them fails.

## The harness has been deliberately broken, and it noticed

A harness that has never failed has not been tested, it has only been run. Five "ports" were built,
each a shell wrapper around the reference, and run against the smoke goldens:

| the "port" | cases failed of 51 |
| --- | --- |
| the reference itself | **0** — and 0 on `sup` as well |
| the reference with `--min_shared` 5 → 6 | 15 |
| one digit changed in one `error_rate` field | 48 |
| the reference plus one extra output file | 51 |
| the reference minus `logfile.txt` | 48 |

`--min_shared 6` catching only 15 of 51 is itself a measurement, and it is about the corpus rather
than the harness: on 280 short reads most cases do not reach a shared-minimizer count where 5 and 6
differ. It is the same argument as PORTING.md's Finding 19 — develop against `sup`.

## How much the matrix actually discriminates

Comparing each case's **whole set of output files** — not `final_clusters.tsv` alone, since
`--consensus` does not change the clustering and sixteen cases share one on `sup`:

| corpus | cases | distinct output sets | collision groups |
| --- | --- | --- | --- |
| `smoke` | 51 | 30 | 5 |
| **`sup`** | 51 | **42** | 4 |

Two of `sup`'s four groups are **intended** and must stay — `default` == `ont` pins that the preset
resolves to `--k 13 --w 20`, and `t8` == `t8_total_nt` pins that `total_nt` is the default
`--batch_type`. The other two are real information: `--rc_identity_threshold 1.0`, `--primer_max_ed 0`
and `--trim_window 50` are the negative controls for the cases that do merge and do trim, and the
three `write_fastq` cases all crash before `--N` is read (PORTING.md, Finding 6).

**Four real bugs were found by those five runs, all of them in this harness or in the reference:**

1. `scrub()` ran `rm -f "$@".bak`, which with two arguments expands to
   `rm -f <first> <second>.bak` and **deletes the first file**. Every one of the 37 CLI cases then
   failed with `diff: .../o: No such file or directory`.
2. `grep -v` exits 1 on empty input, and under `set -o pipefail` that aborted the whole run — after
   the last case that wrote files and before the summary line, so it looked like a crash with no
   verdict. Reached by any case that legitimately writes nothing, which the `write_fastq` cases do.
3. `equivalence.sh stable` found that `cli/use_old_missing` was order-dependent: a failed
   `--use_old_sorted_file` run leaves an **empty** `sorted.fastq` behind, and the second run then
   takes the "use the existing sorted file" path and exits 0 having clustered zero reads. That is a
   reference defect (PORTING.md, Finding 23) that only a harness recording twice could have found.
4. `medaka_consensus --version` is not a supported invocation; it exits non-zero, and under
   `pipefail` that aborted `setup_reference_env.sh` silently, skipping medaka and the
   resolved-versions file.

## Stage oracles

End-to-end goldens say *that* the port is wrong, never *where* — and on the smoke corpus they cannot
even say *that* for eleven of the swept cases. `dump_reference.py` wraps the reference's own
functions and dumps each stage's inputs and outputs in a stable line format:

| stage | what it captures |
| --- | --- |
| `minimizers` | `get_kmer_minimizers` — the ordered `(position, minimizer)` list per read |
| `mapping` | `get_best_cluster` — the candidate ranking and the decision, replayed from the **live driver** |
| `parasail` | `parasail_block_alignment` — cigar and **both** alignment ratios (the second is what `--symmetric_map_align_thresholds` reads) |
| `spoa` | `form_draft_consensus` — the exact sequences handed to spoa, **in order**, and the consensus returned. Insertion order into a POA graph changes the consensus, so the order is contract |
| `identity` | `highest_aln_identity` — forward *and* reverse-complement identity per center pair, at the consensus path's opening penalty of 3 (not the clustering path's binned 2–5) |
| `barcode` | `find_barcode_locations` — every edlib HW call and its **full** locations list, not just the `locations[0]` the reference reads. Which equally-optimal location edlib puts first is not uniquely defined, and it is what a native reimplementation has to reproduce |

The reference is wrapped, never copied. Copying a stage into the harness is the one thing that
guarantees the harness stops measuring the thing it is named after.
