# NGSpeciesID — Rust rewrite

## Goal

Port NGSpeciesID from Python to Rust. **Identical CLI, byte-identical output, and an installation that
works.**

The Python implementation in `NGSpeciesID` and `modules/` is the **normative reference**: when Rust
and Python disagree, Python is right until a human decides otherwise. As in the isONclust port, the
specification is **byte-identity**, and there are **no deliberate divergences** until the port is
exact. Improvements go in *Deferred improvements* and land after exactness, each in its own commit.

Unlike the isONcorrect and isONform ports, **speed and memory are not the point here.** NGSpeciesID
is an amplicon tool; a whole barcode is a few thousand reads. Measured on the committed 3 000-read
fixture the reference clusters in 1.6 s at `--t 1` and 0.7 s at the default `--t 8`. There is nothing
to win. What there is to win is stated below and it is the reason this port is worth doing at all.

### The reason for this port is installation, and it is measured

**Both published ways of installing NGSpeciesID fail on Apple Silicon, and the last step of both
fails on every aarch64 machine.** They are fine on x86_64 Linux, which is presumably why nobody has
noticed. Solved with `conda create --dry-run` per subdir (with `CONDA_OVERRIDE_GLIBC` set so a
cross-platform solve is meaningful), and run verbatim on Darwin 25.5.0 / arm64 with miniforge3:

| Published instruction | linux-64 | osx-64 | osx-arm64 |
| --- | --- | --- | --- |
| README, "Recent update (2025-04-19)": `python=3.11 … medaka==2.0.1 openblas==0.3.3 spoa racon minimap2 samtools` | resolves | **`PackagesNotFoundError`** | **`PackagesNotFoundError`** |
| README, "Published installation instructions (2021-01-11)": `python=3.6 … medaka==0.11.5 openblas==0.3.3 spoa racon minimap2` | resolves | **fails**, nothing provides `libopenblas` | **fails**, Python 3.6 was never built for `osx-arm64` |
| `pip install NGSpeciesID` (both recipes end with it) | wheel | wheel | **fails after 1 m 45 s** building `parasail==1.2.4` |

Why each fails where it does:

- **medaka pins.** `medaka==2.0.1` exists for `linux-64` but for neither macOS subdir; on `osx-arm64`
  medaka exists only as 2.1.1, 2.2.0, 2.2.1, 2.2.2. `openblas==0.3.3` exists for `linux-64` and not
  for `osx-arm64`, whose earliest build is 0.3.11.
- **parasail from PyPI.** `setup.py` pins `parasail==1.2.4`, and **neither 1.2.4 nor 1.3.4 publishes
  an aarch64 wheel for any operating system.** 1.2.4 ships `macosx_10_9_x86_64`,
  `manylinux2010_i686`, `manylinux2010_x86_64`, `win32`, `win_amd64`; 1.3.4 ships
  `macosx_10_9_x86_64`, `manylinux_2_17_x86_64`, `musllinux_1_1_i686`, `win32`, `win_amd64`.
  `edlib` has no aarch64 wheel either (38 wheels, none). So Apple Silicon **and ARM Linux** fall back
  to a source build, and parasail's fails with `RuntimeError: autoreconf -fi failed` — an error
  naming neither the `m4` it rejected nor the `glibtoolize` it could not find. Measured: 1 m 45 s for
  1.2.4, 1 m 32 s for 1.3.4, both in a clean env with no cache.

#### And the fix is one flag, because conda already has the wheels

This is the part worth knowing, and it is the opposite of what the isONclust port concluded. That
port recorded parasail as "not on conda-forge or bioconda under any name (`parasail`,
`parasail-python`, `libparasail`: all absent)" and built it from source with a two-symlink
workaround. **That is now stale.** Measured today:

| package | osx-arm64 | note |
| --- | --- | --- |
| `parasail-python` 1.3.4 | **on bioconda**, py39/py312/py313 builds | newest build dated 2025-07-22 |
| `parasail` 2.6.2 (the C library) | **on bioconda** | |
| `python-edlib` 1.3.9.post1 | **on bioconda**, py38–py313 builds | |

And they arrive for free: **`medaka` 2.2.2 depends on both `parasail-python` and `python-edlib`.**
Confirmed by installing only `python=3.12 pip medaka spoa racon minimap2 samtools` into a fresh
environment and finding `parasail-python 1.3.4 py312h7d4aa78_6` and
`python-edlib 1.3.9.post1 py312h2ee7600_3` already present, importing, and correct
(`sg_trace_scan_16` → score 16, cigar `4=1X5=`, not saturated; `edlib.align(mode="HW",
task="locations")` → `[(2, 5)]`).

So the only broken step is the last one: `pip install NGSpeciesID` reads `install_requires` and
re-installs `parasail==1.2.4` **from PyPI, over the working conda build**, and that source build
fails on aarch64. On x86_64 the PyPI wheel exists, the override succeeds, and nothing looks wrong.

The recipe that works, on all three subdirs, with no compiler and no workaround:

```
conda create -n NGSpeciesID -c conda-forge -c bioconda python=3.12 pip medaka spoa racon minimap2 samtools
conda activate NGSpeciesID
pip install --no-deps NGSpeciesID
```

80 seconds for the solve on osx-arm64, yielding medaka 2.2.2, spoa 4.1.5, racon 1.5.0,
minimap2 2.31-r1302, samtools 1.24, python 3.12.14 — and `--no-deps` is the whole fix.

Two consequences beyond the README:

1. **`bench/setup_reference_env.sh` does not need the parasail source build at all**, which removes
   the most fragile part of the harness isONclust had to write.
2. `python=3.12` is not an arbitrary choice: medaka 2.2.2 requires `python >=3.12,<3.13`, and 3.12 is
   also the interpreter this port must pin for *Finding 2*. The two constraints agree.

### The clustering half is already done, and that is measured too

NGSpeciesID is a fork of isONclust with consensus, primer removal and polishing bolted on. The
isONclust Rust port is finished and byte-identical to *its* reference. The question that decides the
size of this project is whether it is also byte-identical to *this* one, and the answer is yes.

First, the two Python programs were compared directly. On the committed 3 000-read fixture at
`--k 13 --w 20 --t 1`, `isONclust` and `NGSpeciesID` produce:

| file | result |
| --- | --- |
| `logfile.txt` | **byte-identical** |
| `final_clusters.tsv`, `final_cluster_origins.tsv`, `sorted.fastq` | identical **after one substitution**: isONclust's `readfq` does `last[1:].replace(" ", "_")` and NGSpeciesID's does not, so every accession differs and nothing else does |

Normalised for that one difference, the cluster assignment diff is **0 lines over 3 000 reads**, and
both produce 49 clusters. Then the finished Rust binary was run against the Python reference of *this*
repository:

| | |
| --- | --- |
| configurations | 12: `k13w20`, `k15w50`, `k20w100`, `q8`, `min_shared 10`, `aligned_threshold 0.9` × 2 corpora |
| plus | `--t 1`, `--t 4`, `--t 8` on the 3 000-read corpus, including the per-iteration intermediates |
| `final_clusters.tsv` (accession-normalised) | **identical in all 12** |
| `final_cluster_origins.tsv` columns 1, 5, 6 (id, score, error rate) | **identical in all 12** |
| `logfile.txt` | **identical in all 12** |

And the empirical probability table matches: NGSpeciesID's 41 880 rows are **exactly** isONclust's
41 880 rows with `k >= 10`, same values, same order. isONclust's extra 17 748 rows are the `k` 4–9
support added by its PR #13, which this repository never took.

So the port does not start from zero. It starts from a finished, verified clustering engine, and the
new work is:

1. the accession handling (**delete** isONclust's `replace(" ", "_")` — it is the one line that differs);
2. `--m`/`--s` length filtering and `--sample_size`/`--top_reads` subsampling, which sit between the
   sort and the sweep and do not exist in isONclust;
3. `--symmetric_map_align_thresholds`, the one genuinely new piece of clustering logic;
4. the whole consensus stage: POA, reverse-complement detection, primer/tail trimming, and the
   medaka/racon drivers;
5. restricting the probability table to `k >= 10`, because shipping the full one would silently
   *fix* *Finding 8* and that is a behaviour change;
6. reproducing this repository's crash contract, which is not isONclust's — two of isONclust's
   upstream fixes were never taken here (*Finding 7*, *Finding 5*).

The entry point must keep its exact current name, flags and defaults:

- `NGSpeciesID` — clusters a fastq and optionally forms, trims and polishes consensus sequences
- `NGSpeciesID write_fastq` — a subcommand that splits a clustering into per-cluster fastq files.
  It is currently **unreachable and broken**; see *Finding 6*

## Branches

Two branches, because the installation fixes should not wait for the port.

| branch | base | contents | PR |
| --- | --- | --- | --- |
| `fix/installation` | `master` | the three non-port fixes: working install instructions, the unpinned parasail, the `--q` guard, and the untracking of six junk files | **now** — small, self-contained, and it fixes a tool that is currently uninstallable on Apple Silicon |
| `develop` | `master` | `PORTING.md`, `tools/repo-slim/`, `bench/`, and eventually `rust/` | **when the port is exact** |

Burying the install fixes behind a long-lived port PR would waste them: they are useful this week and
the port is not. `develop` gets `master` merged in after `fix/installation` lands.

`origin/develop` exists and is **fully merged into `master`** — zero unique commits, and `master` is
18 commits ahead of it. Unlike isONclust's, it holds nothing that would be lost, so re-creating it
from `master` costs nothing and needs no archival step.

The history rewrite is on neither branch. It is a force-push, it is reviewed on its own, and
`tools/repo-slim/` only carries the tooling and the analysis it produced.

## Layout

| Path | Role |
| --- | --- |
| `NGSpeciesID` | Reference: CLI, argparse, the `--m`/`--s`/`--sample_size` filter, the output writers, the `--consensus` driver, `write_fastq` |
| `modules/cluster.py` | Reference: minimizers, hit collection, mapping and alignment decisions, `reads_to_clusters`. **Diverges from isONclust's only by `--symmetric_map_align_thresholds`, logging, and the 8th tuple element that flag needs** |
| `modules/get_sorted_fastq_for_cluster.py` | Reference: quality scoring, filtering, the score sort, `sorted.fastq`. isONclust's minus the BAM path, and **without** isONclust's *Finding 9* fix |
| `modules/parallelize.py` | Reference: batching and the hierarchical merge driven by `--t`. isONclust's, but with `yield batch` unguarded — see *Finding 5* |
| `modules/consensus.py` | Reference: the spoa subprocess, parasail alignment, reverse-complement detection, medaka and racon drivers. **This is the new half of the port** |
| `modules/barcode_trimmer.py` | Reference: edlib HW-mode primer and universal-tail location and trimming, with an IUPAC equivalence map. **New; isONclust has no equivalent** |
| `modules/help_functions.py` | Reference: `readfq`, `cigar_to_seq`, `mkdir_p`. Note `readfq` here does **not** substitute spaces in accessions |
| `modules/p_minimizers_shared.py` | Reference: a 1.79 MB Python literal — 41 880 rows of empirical minimizer-sharing probabilities, `k` 10–30, `w` 10–100. Data, not code. Identical to isONclust's `k >= 10` subset |
| `rust/` | The port. Does not exist yet — this is where the next session starts |
| `bench/` | The equivalence harness, carried across from `isONclust/bench/` and adapted. `equivalence.sh`, `cases.tsv` (51 cases), `corpora.tsv`, `dump_reference.py` (6 stages), `diffsummary.py`, `setup_reference_env.sh`, `golden/{smoke,sup}/`, and its own `README.md`. **Runs; goldens recorded on both corpora** |
| `tools/repo-slim/` | The staged history-rewrite tooling, carried across from `isONclust/tools/repo-slim/` and adapted. `analyze.sh` has been **run**; `archive_data.sh` and `slim.sh` have not. See *Repo hygiene* |
| `test/sample_h1.fastq` | 280 ONT reads, 390 KB. The README's install check and the `.travis.yml` fixture. **Blind to five parameters** — *Finding 19* |
| `test/Supplementary_File1_reads.fastq` | 3 000 ONT reads from three fish species, 5.2 MB. The paper's supplementary data. **This is the corpus to develop against** — *Finding 19* |
| `test/Supplementary_File2_minibar.txt`, `test/Supplementary_File3_primer.txt` | The paper's demultiplexing index file and primer file. The primer file carries IUPAC codes and is what exercises `barcode_trimmer` |
| `scripts/` | Paper experiment scripts, referenced from nothing. **Not part of the port. Do not modify.** |
| `Dockerfile` | Pins `medaka==1.5.0` and biocontainer binaries. Not a reference; see *Repo hygiene* |
| `.travis.yml` | Dead CI: Travis, Python 3.6, `medaka=0.11.5`. Replace, do not migrate |

## Port status

Nothing is ported yet. This document is the reconnaissance, and it is the deliverable of the session
that wrote it.

| Stage | State | Verification |
| --- | --- | --- |
| reference environment | **done, and scripted** | `bench/setup_reference_env.sh`: one conda line, no parasail source build. Resolves to python 3.12.14, parasail-python 1.3.4, python-edlib 1.3.9.post1, spoa 4.1.5, racon 1.5.0, minimap2 2.31-r1302, samtools 1.24, medaka 2.2.2 |
| published install paths measured | **done, and fixed** | all three fail on `osx-arm64`, two also on `osx-64`, all three resolve on `linux-64`. *Goal*. Fixed on `fix/installation` (`76d3f88`), verified end to end from a clean env |
| reference runs on real data | **done** | 3 000 reads → 49 clusters at `--t 1`, 33 at `--t 8`, in 1.6 s / 0.7 s |
| determinism gate | **done, and it fails twice** | *Finding 1* (`--sample_size`) and *Finding 2* (Python ≤3.11) |
| interpreter decision | **taken: pin ≥3.12** | same decision as isONclust, same reason. The README recommending 3.11 is *Finding 2* |
| `--sample_size` decision | **not taken. This is the one blocking question in the document** | *Finding 1* |
| CLI contract captured | **done, 44 cases recorded** | `bench/golden/<corpus>/cli/` — exit code, stdout and stderr, scrubbed of paths, timings and traceback line numbers. *The exit-code contract* |
| output goldens recorded | **done, 51 cases on both corpora** | `bench/golden/<corpus>/manifest.tsv`. Includes every `--consensus` case, both polishers, and the `--t > 1` merge intermediates |
| goldens are reproducible | **done** | `equivalence.sh stable` records the whole matrix twice and diffs: 191 checks, identical. It found two real defects on the way — *Finding 23* and *Finding 24* |
| the harness itself is tested | **done** | five deliberately-broken "ports" run against the goldens; see *Has the harness got teeth?* |
| stage oracles | **written and exercised on both corpora; the replay half waits for the port** | `bench/dump_reference.py` covers six stages, three of them new here. *Finding 25* has the coverage counts |
| corpora | **done, 2 committed, both measured for discriminating power** | *The corpora*, *Finding 19* |
| case matrix swept on both corpora | **done** | 24 cases; `Supplementary_File1_reads.fastq` gives 19 distinct results and 1 unintended collision, `sample_h1.fastq` gives 12 and 8 |
| repository slimmed | **not started; analysed, and the tooling is in the tree and runs** | `tools/repo-slim/analyze.sh` reports 520.1 MB of 530.8 MB strippable (98.0%), 15 paths, and writes a reviewed `removal-paths.txt`. *Repo hygiene* |
| clustering engine | **already exists and is verified against this reference** | 12 configurations × 2 corpora, `--t` 1/4/8. *Goal* |
| `--symmetric_map_align_thresholds` | not started | new logic; visible on the 3 000-read corpus (85 clusters vs 49) and **invisible** on `sample_h1` |
| `--m`/`--s`/`--sample_size`/`--top_reads` | not started | blocked on *Finding 1* for the random path only; `--top_reads` and `--m`/`--s` are reproducible today |
| POA (`spoa`) | not started, **de-risked** | the invocation is byte-for-byte isONcorrect's, where `spoars` was measured identical on 505/505 cases. *spoa* |
| reverse-complement detection | not started | needs `parasail_alignment` at `opening_penalty=3` — a second call site with different defaults from the clustering path. The port already has parasail |
| primer / universal-tail trimming | not started | needs edlib **HW** mode with `task="locations"` and `additionalEqualities`. Different from the NW/CIGAR path isONcorrect reimplemented. *The aligners* |
| medaka / racon drivers | not started | stay subprocesses. Measured reproducible: 3 runs of `--consensus --racon` differ only in the captured minimap2 and racon **stderr logs**, which carry timings; every consensus fasta is identical |
| `write_fastq` | not started | broken in the reference (*Finding 6*). Decide before porting |

## The pipeline, in one pass

One fastq in, one clustering out, and optionally a consensus per abundant cluster.

Steps 1–7 are isONclust's, unchanged, and are specified in full in `isONclust/PORTING.md`. Only the
differences and the additions are written out here.

1. **Score and filter every read** (`get_sorted_fastq_for_cluster`). Homopolymer-compress; drop if
   `len(seq) < 2*k` or `len(hpol_compressed) < k`; drop if `10 * -log10(error_rate) <= --q`. Two phred
   tables, differing by one cap at `0.79433`: the rolling product uses the **capped** one, `error_rate`
   the **uncapped** one.

2. **Sort descending by score** and write `sorted.fastq`, with the score appended to every accession
   as `acc + "_" + str(score)` and parsed back out downstream as `float(acc.split("_")[-1])`. So
   `sorted.fastq` is an output file and the input to everything after it, and Python's `str(float)`
   formatting is part of the contract.

   **The accessions here contain spaces**, because this repository's `readfq` does not substitute
   them. Every ONT header does. That is the only textual difference from isONclust, it reaches every
   output file, and it is what breaks `write_fastq` (*Finding 6*).

3. **Filter by length and subsample — NGSpeciesID only.** Read `sorted.fastq` back into
   `read_array`, keeping reads with
   `--m - --s <= len(seq) <= --m + --s` when both are positive. Then, if `--top_reads`, take the
   first `--sample_size`; else if `0 < --sample_size < len(read_array)`, take
   `sorted(random.sample(range(len(read_array)), --sample_size))` — **an unseeded draw from the
   global `random` module**, which is *Finding 1*. Then
   `abundance_cutoff = int(--abundance_ratio * len(read_array))`, computed on the **post-subsample**
   count.

4. **Load the empirical probability table.** Keep rows where `k == --k` and `abs(w - --w) <= 2`, keyed
   by rounded error-rate pairs, inserted in **both** orders.

5. **Initialise every read as its own cluster and representative**, then **sweep in sorted order**
   (`reads_to_clusters`): homopolymer-compress, take minimizers (lexicographic, **first** occurrence
   on ties), compute the compressed error rate, collect hits, try to map, else try to align, else
   become a representative and add your minimizers to the database.

   **`--symmetric_map_align_thresholds` changes both decisions.** The representative tuple grows a
   7th element (the compressed error rate, as in isONclust) *and* an 8th (the representative's
   homopolymer-compressed sequence, which isONclust does not keep). With the flag set:

   | decision | default gate | symmetric gate |
   | --- | --- | --- |
   | mapping | `total_mapped / len(read_compressed) > --mapped_threshold` | `min(that, total_mapped / len(rep_compressed)) > --mapped_threshold` |
   | alignment | `sum(aligned_region) / len(s1) >= --aligned_threshold` | `min(that, sum(aligned_region) / len(s2)) >= --aligned_threshold` |

   and the value *returned* is the `min`, not the read-side ratio. The flag exists because a long,
   high-quality read containing an insertion the cluster does not share is not penalised by the
   default gate.

   The 8th element is why *Finding 18* exists: a read that reaches `continue` in this loop keeps a
   6-element tuple that the output writer unpacks as 8.

6. **Reassign**, one level deep — merge targets are always representatives and representatives are
   never merge sources.

7. **Write output, ordered by `(cluster size, representative score)` descending**, with fresh
   `output_cl_id` counting from 0, reads within a cluster sorted by score descending, and the score
   suffix stripped from accessions on the way out. Ties fall back to dict order, which after step 6
   is ascending internal id — so the port's key is `(-size, -score, internal_id)`.

8. **Optionally form, trim and polish consensus sequences (`--consensus`) — NGSpeciesID only, and
   this is the half that does not exist in the isONclust port.**

   a. **`form_draft_consensus`.** Walk clusters in the same `(size, score)` order. For each cluster
      of at least `abundance_cutoff` reads, write its reads to a temp fastq — **in cluster order, not
      score order**, and stopping at `--max_seqs_for_consensus` when that is `>= 0` — and run
      `spoa <reads> -l 0 -r 0 -g -2`, taking line 2 of the output as the consensus. Yields
      `[nr_reads, c_id, center, reads_path]`. **Sequence insertion order into the POA graph changes
      the consensus**, so the order reads are written in is contract.

   b. **Primer or universal-tail trimming, if requested** (`barcode_trimmer.remove_barcodes`), before
      reverse-complement detection. For each center, search a window of `--trim_window` bases at each
      end (halved to `len(center)//2` if `2*--trim_window > len(center)`) with
      `edlib.align(primer, window, mode="HW", task="locations", k=--primer_max_ed,
      additionalEqualities=IUPAC_map)`, taking `locations[0]` — the **first** location edlib returns.
      Cut start is the **latest** `stop` among beginning hits; cut end is derived from the
      **earliest** `start` among end hits. Mutates `centers[i][2]` in place and reports whether
      anything changed.

   c. **`detect_reverse_complements`.** For each center in order, align it against every later
      center's forward *and* reverse-complement forms with parasail at `opening_penalty=3` — not the
      clustering path's error-rate-binned 2–5 — take `max` of the two identities, and merge any
      center at or above `--rc_identity_threshold` into the earlier one. **This double-counts**;
      *Finding 10*.

   d. **`polish_sequences`.** Delete any `consensus_reference_*` and `medaka_cl_id_*`/`racon_cl_id_*`
      left from a previous run, write `consensus_reference_<c_id>.fasta` and
      `reads_to_consensus_<c_id>.fastq` (accessions truncated at the first space here, and only
      here), then run either `medaka_consensus -i ... -d ... -o ... -t 1` or
      `minimap2 -x map-ont` + `racon`, `--racon_iter` times. Reads line 2 of the polisher's output
      back into `centers[i][2]`.

   e. **If a primer file or `--remove_universal_tails` was given, trim again** on the polished
      centers, and if that changed anything, re-run (c) and (d). This is the loop commit `5463966`
      guards.

   Reaching `polish_sequences` with neither `--medaka` nor `--racon` is *Finding 4*.

### With `--t > 1` this is a different algorithm

`--t` does not parallelise the sweep; it replaces it (`parallelize.parallel_clustering`), exactly as
in isONclust. Reads are cut into `--t` batches by `--batch_type`, each batch is clustered in its own
process with its own minimizer database, survivors are pooled and re-clustered, and this repeats
until one batch remains.

**The default is `--t 8`, so the default invocation is the parallel path.** Measured on the
3 000-read corpus:

| `--t` | clusters | wall clock | extra directories written |
| --- | --- | --- | --- |
| 1 | 49 | 1.62 s | — |
| 2 | 42 | 1.07 s | `1/` |
| 4 | 35 | 0.80 s | `1/`, `2/` |
| 8 | **33** | 0.71 s | `1/`, `2/`, `3/` |

Each directory holds `pre_clusters.csv` and `cluster_origins.csv`. These are observable output and
the harness must diff them. Every `--t` value is its own equivalence case; a port that "parallelises
with rayon" and expects `--t 1` output is not a port.

Processes are spawned with `mp.set_start_method('spawn')`, so workers do not inherit parent state.

## Scope: what gets ported

### Ported — inside the equivalence contract

All 38 live flags. Clustering: `--fastq`, `--outfolder`, `--version`, `-h`/`--help`, `--debug`, `--k`,
`--w`, `--q`, `--t`, `--d`, `--ont`, `--isoseq`, `--min_shared`, `--mapped_threshold`,
`--aligned_threshold`, `--symmetric_map_align_thresholds`, `--min_fraction`, `--min_prob_no_hits`,
`--batch_type`, `--use_old_sorted_file`, `--m`, `--s`, `--sample_size`, `--top_reads`. Consensus:
`--consensus`, `--abundance_ratio`, `--rc_identity_threshold`, `--max_seqs_for_consensus`, `--medaka`,
`--racon`, `--medaka_model`, `--medaka_fastq`, `--racon_iter`, `--remove_universal_tails`,
`--primer_file`, `--primer_max_ed`, `--trim_window`. Plus the `write_fastq` subcommand with
`--clusters`, `--fastq`, `--outfolder`, `--N`.

`--ont` is exactly `--k 13 --w 20`; `--isoseq` is exactly `--k 15 --w 50`. Both need their own
equivalence case, because a preset silently resolving to the wrong numbers is invisible in output
that agrees for other reasons. Both are applied **after** explicit `--k`/`--w`, so `--ont --w 5`
silently becomes `--w 20` — measured, and contract.

**`--sample_size` without `--top_reads` cannot be inside a byte-identity contract as it stands.**
See *Finding 1*.

### Nothing is dropped

isONclust's port dropped `--consensus` on the grounds that consensus of a cluster is isONcorrect's
and isONform's job. That argument does not transfer: consensus **is** NGSpeciesID. Everything in the
tool is in scope.

The one thing that goes away for free is the BAM input path — `--ccs`/`--flnc` and the pysam
dependency were already deleted from this repository in commit `0c41cad`, so there is nothing to
defer.

### Diagnostic only

`--d` (`print_output`) controls a progress table emitted through `logging.debug` and does not change
any output file. Port it best-effort; matching the text is not part of the contract. `--d 0` is not a
flag to drop but an input to reproduce — see *Finding 11*.

`--debug` switches the log level. Everything the reference prints goes through `logging` with
`format='%(message)s'`, and much of it carries timings, so **log output is not in the byte-identity
contract**; the CLI goldens scrub timings and compare the rest. Note that at default level the tool
prints exactly three lines (`Starting Clustering`, `Finished Clustering`, and with `--consensus` two
more), which *is* worth pinning.

### External tools that stay external

`medaka`, `racon` and `minimap2` are invoked as subprocesses and their output files are theirs. The
port must invoke them with byte-identical argument vectors, in the same order, writing to the same
filenames — that is what is verifiable, and it was measured to be enough: three runs of
`--consensus --racon` produced identical `consensus.fasta` and differed only in the captured
`mm2_stderr_it_*.txt` and `racon_stderr_it_*.txt`, which contain timings.

Missing tools must **exit non-zero and name the tool**. Today they surface as a
`FileNotFoundError` traceback from `subprocess`, which is a non-zero exit with the right information
buried in it; the port should say `racon not found on PATH` and nothing else. This is the one place
where improving on the reference costs nothing measurable, and it should still be its own commit.

## spoa

`run_spoa` invokes:

```
spoa <reads.fq> -l 0 -r 0 -g -2
```

This is **byte-for-byte the invocation isONcorrect uses**, and isONcorrect's port measured `spoars`
against the real `spoa` binary on it: **505 of 505 unique cases produced an identical consensus**,
across 8 parameter combinations, with up to 28 sequences per POA. Resolved against spoa's CLI
defaults the configuration is `kSW` (local), **linear** gap `-2`, match `+5`, mismatch `-4`,
consensus only.

```rust
// spoa <fq> -l 0 -r 0 -g -2  ->  kSW (local), linear gap -2, m=+5, n=-4
let scoring = Scoring::new(5, -4, -2, -2, -2, -2)?;
let mut engine = SimdEngine::new(AlignmentType::Local, scoring);
```

Two differences from isONcorrect's use that must be checked rather than assumed:

1. **The input is a fastq, not a fasta.** spoa reads both, but confirm `spoars` is fed the same
   sequences in the same order — including the `--max_seqs_for_consensus` cutoff, which is
   `i >= max_seqs_for_consensus` (a `>=`, admitting exactly `max_seqs_for_consensus` sequences,
   unlike isONcorrect's `>`).
2. **The sequences are far longer and there are far more of them.** isONcorrect POAs correction
   intervals; this POAs whole amplicon reads — 300 to 1 600 bases, up to the full cluster. The
   oracle has to be re-recorded on this corpus, not inherited.

**Decision: adopt `spoars`, gated on its own differential oracle against the `spoa` binary recorded
from this repository's corpora.** If it fails, the fallbacks are `spoa`/`spoa-sys` bindings (identical
by construction, at the cost of a C++ toolchain) or keeping the subprocess (identical by
construction, at the cost of the installation win for `--consensus`).

## The aligners

Three call sites, three different problems.

| Call site | What it needs | Status |
| --- | --- | --- |
| `cluster.parasail_block_alignment` — the clustering decision | parasail `sg_trace_scan_16`, `match 2, mismatch -2, gap_ext 1`, gap opening binned 5/4/3/2 by summed error rate, falling back to `_32` on saturation | **done** in the isONclust port: 18 633 alignments identical, CIGAR *and* ratio |
| `consensus.parasail_alignment` — reverse-complement detection | the same parasail call at **`opening_penalty=3`**, and the identity is computed by zipping the two gapped strings and counting mismatches, so it charges leading and trailing gaps | port's `parasail.rs` covers the call; the **default is different** and the identity formula is new. Needs its own oracle |
| `barcode_trimmer.find_barcode_locations` — primer and tail trimming | edlib **HW** (infix) mode, `task="locations"`, `k=--primer_max_ed`, and a 36-pair `additionalEqualities` IUPAC map | **new.** Not what isONcorrect reimplemented |

The third is the one that needs a decision. isONcorrect's `align.rs` reimplements edlib's **NW**
traceback, and the tie-break had to be *measured* — only one of six preference orderings reproduces
edlib. HW mode with `task="locations"` is a different question and a smaller one: the return is a set
of end positions and a `[start, end]` pair per alignment, `locations[0]` is all that is read, and the
edit distance is uniquely defined. What is **not** uniquely defined is which location edlib puts
first when several achieve the optimum, and the reference reads exactly that one. Measure it against
a recorded corpus before writing it, the way isONcorrect did.

The IUPAC map is data and must be copied verbatim, including its asymmetry: it maps `('M','A')` but
not `('A','M')`, so an ambiguity code in the *primer* matches a concrete base in the *center* and not
the other way round.

`edlib_rs` and `rsedlib` are both C++ bindings and both were measured unusable in the isONcorrect
port — `edlib_rs` needs `CMAKE_POLICY_VERSION_MINIMUM=3.5`, `rsedlib` fails to link. Vendoring
edlib's single `.cpp` with `cc` is the fallback if a native HW implementation cannot be made exact,
and it is a better fallback here than there, because it is one file and no CMake.

## Determinism rules

The default code path is deterministic **on Python ≥3.12 only, and only without `--sample_size`**.
Read *Finding 1* and *Finding 2* before this list.

- **`--sample_size` without `--top_reads` is not reproducible at all.** *Finding 1*.
- **`sum()` over floats is compensated from CPython 3.12 and a naive left-fold before it.** The sites
  that matter sum over `set(qual)` and `set(qualcomp)`, whose iteration order is
  `PYTHONHASHSEED`-dependent. The port must reproduce the **exactly rounded** sum (Neumaier or
  equivalent), which then frees it to iterate in any order.
- **`reduce(mul, [p]*n, 1)` is not `p.powi(n)`.** It is a left-fold of `f64` multiplications starting
  from the integer `1`, and the port must fold the same way.
- **Minimizer ties go to the FIRST position in the window**, via `list(window_kmers).index(curr_min)`.
  isONcorrect takes the last. Do not carry that habit across.
- **`get_kmer_minimizers` can emit the same `(minimizer, position)` twice**, and `get_all_hits` counts
  per minimizer index, so duplicates inflate hit counts. Reproduce it.
- **`get_best_cluster`'s sort is a total order**, because its third key is the representative accession
  and accessions are unique — so the set-iteration order feeding it never reaches output. Checked in
  the isONclust port, not assumed.
- **Cluster output ids depend on dict order.** `sorted(..., reverse=True)` is stable, so ties in
  `(size, score)` resolve to insertion order, which is ascending internal cluster id. Sort by
  `(-size, -score, id)`.
- **`already_removed` is a Python `set` of cluster ids**, iterated never — only membership-tested. Safe.
- **`centers` is a list and its order is the `(size, score)` order from step 8a.** Both
  `detect_reverse_complements` and `polish_sequences` depend on it, and `detect_reverse_complements`
  is order-dependent in a way that matters (*Finding 10*).
- **`work_dir` is a `tempfile.mkdtemp()`**, so the per-cluster read files live at a path that changes
  every run. The path reaches no output file — but it reaches the argument vector handed to spoa, and
  `reads_path_name` is stored in `centers` and read back in `polish_sequences`. Reproduce the
  *structure*, not the path.
- **`write_fastq` iterates a `defaultdict` in insertion order**, which is the order cluster ids first
  appear in `final_clusters.tsv`. It decides only the order files are created in.
- Dict iteration in the reference is insertion-ordered (Python 3.7+). Where iteration order feeds
  output, use an order-preserving map in Rust, not `HashMap`.

## Verification

Byte-identity is the acceptance criterion, and it is checked, not assumed. `bench/` is carried across
from `isONclust/bench/` and adapted; it exists and runs. `bench/README.md` is its own documentation.

```bash
bench/setup_reference_env.sh                       # one conda line; no parasail source build
export PATH=~/miniforge3/envs/ngspeciesid-ref/bin:$PATH
bench/equivalence.sh env                           # usable? sum() compensated? which tools?
bench/equivalence.sh seeds                         # does the reference agree with itself?  <-- FIRST
CORPUS=sup GOLDEN=$PWD/bench/golden/sup bench/equivalence.sh cli record
CORPUS=sup GOLDEN=$PWD/bench/golden/sup bench/equivalence.sh record
CORPUS=sup GOLDEN=$PWD/bench/golden/sup bench/equivalence.sh stable
CORPUS=sup GOLDEN=$PWD/bench/golden/sup bench/equivalence.sh verify
bench/equivalence.sh tools                         # a missing binary must be named, not tracebacked
```

Where it stands, measured:

| check | result |
| --- | --- |
| `env` | 8 checks green; all four external tools present |
| `seeds`, python 3.12 | **30 green.** Three entries — `--ont --t 1`, `--ont --t 8`, `--ont --t 1 --consensus --racon` — stable across five `PYTHONHASHSEED` values, including all six `--t 8` merge intermediates and the whole racon output tree |
| `seeds`, python 3.11 | **fails, exit 1**, on `final_cluster_origins.tsv` and `logfile.txt`, and `diffsummary.py` names the column: `error_rate`. This is what proves the gate has teeth |
| `seeds`, `--sample_size` | asserts *Finding 1* explicitly: 5 identical runs, 5 distinct results. It **fails if someone fixes the defect** without updating the harness |
| `cli record` | **44 cases** |
| `record` | **51 cases**, both corpora |
| `stable` | **191 checks, identical across two recordings** |

**What counts as a difference.** Every file the tool writes, byte for byte:

| always | `final_clusters.tsv`, `final_cluster_origins.tsv`, `sorted.fastq`, `logfile.txt` |
| `--t > 1` | `<n>/pre_clusters.csv`, `<n>/cluster_origins.csv` for each merge iteration |
| `--consensus` | `consensus_reference_<c_id>.fasta`, `reads_to_consensus_<c_id>.fastq` |
| `--consensus --medaka` | `medaka_cl_id_<c_id>/consensus.fasta` (or `.fastq`), and the index and BAM files medaka leaves beside the draft |
| `--consensus --racon` | `racon_cl_id_<c_id>/consensus.fasta`, `racon_polished_it_<i>.fasta`, `read_alignments_it_<i>.paf` |
| never | `*_stderr_it_*.txt`, `mm2_stderr_it_*.txt`, `stdout.txt`, `stderr.txt` — measured to differ between runs of the reference itself, because they carry timings |

A file the port fails to write and a file it writes that the reference does not are both failures.
Both `record` and `seeds` must **discover** that list with `find` rather than holding their own copy
of it — isONclust's `seeds` held its own and therefore never checked the parallel-mode intermediates
for seed independence, and the file set here is larger and more conditional than isONclust's.

`sorted.fastq` is **not** an intermediate to be skipped: the score is formatted into every accession
and parsed back out downstream, so it is an output and an input at once.

`logfile.txt` is included because it is the only place the error-rate distribution is observable, and
`error_rate` is precisely where *Finding 2* surfaces.

**The goldens are a manifest of hashes, not the files**, and the manifest records the corpus sha256,
the interpreter version, whether `sum()` was compensated, **and the version of every external tool on
`PATH`** — spoa, racon, minimap2, medaka. The goldens are only valid for the environment that
produced them, and here that environment includes four binaries whose output is part of the contract.

### Has the harness got teeth? Measured, and it found four bugs

*A harness that has never failed has not been tested, it has only been run.* So five "ports" were
built — each a shell wrapper around the reference — and run against the smoke goldens:

| the "port" | cases failed of 51 |
| --- | --- |
| the reference itself, unmodified | **0** |
| the reference with `--min_shared` 5 → 6 | 15 |
| one digit changed in one `error_rate` field of one file | 48 |
| the reference plus one extra output file | 51 |
| the reference minus `logfile.txt` | 48 |

`--min_shared 6` catching only 15 of 51 is a measurement about the **corpus**, not the harness: on
280 short reads most cases never reach a shared-minimizer count where 5 and 6 differ. Same conclusion
as *Finding 19*.

Those five runs found **four real defects**, three of them in the harness and one in the reference:

| # | defect | how it presented |
| --- | --- | --- |
| 1 | `scrub()` ran `rm -f "$@".bak`, which with two arguments expands to `rm -f <first> <second>.bak` and **deletes the first file** | all 37 CLI cases failed with `diff: .../o: No such file or directory` |
| 2 | `grep -v` exits 1 on empty input, and under `set -o pipefail` that aborted the run | the run ended **after the last case that wrote files and before the summary line** — it looked like a crash with no verdict. Reached by any case that legitimately writes nothing, which the `write_fastq` cases do |
| 3 | `medaka_consensus --version` is not a supported invocation and exits non-zero; same `pipefail` shape | `setup_reference_env.sh` ended silently after samtools, skipping medaka *and* the resolved-versions file |
| 4 | **in the reference:** a failed `--use_old_sorted_file` run leaves an empty `sorted.fastq`, and the retry exits 0 on zero reads | *Finding 23*. Found by `stable`, which records twice in one process: exit 1, then exit 0 |

Defect 2 is the same shape as the one `tools/repo-slim/analyze.sh` had — a `set -e` interaction that
only fires on a path nobody had taken. Three instances of it in one session is enough to call it a
pattern: **in these scripts, put `|| true` inside the substitution, not after it.**

### The exit-code contract

**44 cases are recorded** in `bench/golden/<corpus>/cli/`, each with its exit code, stdout and
stderr. Every exit code below is measured, and several of them are wrong in an interesting way and
are contract regardless:

| Case | Exit | Note |
| --- | --- | --- |
| no arguments | **2** | argparse: the `--fastq`/`--use_old_sorted_file` group is `required=True`. The `len(sys.argv)==1: parser.print_help()` branch in `__main__` is **dead code** |
| unknown flag | **2** | and argparse reports the *missing required group*, not the unknown flag |
| `--version` | 0 | `NGSpeciesID 0.3.1` |
| `--ont --isoseq` together | **0** | prints "Arguments mutually exclusive" and exits 0 — a pipeline cannot detect this |
| `--medaka --racon` together | 2 | argparse's mutually-exclusive group |
| `--remove_universal_tails --primer_file` together | 2 | same |
| `--outfolder` omitted | 1 | `TypeError: expected str, bytes or os.PathLike object, not NoneType`. *Finding 9* |
| `--w` < `--k`, or `--w` > 100 | 1 | one shared message for both |
| `--d 0` | 1 | `ZeroDivisionError`. *Finding 11* |
| `--k 9` (or any `k` outside 10–30) | 1 | `KeyError` on an empty probability table. *Finding 8* |
| `--q 12` on the 3 000-read corpus | 1 | was `IndexError: list index out of range`; now two lines saying which flag filtered everything, same exit code. *Finding 7*, fixed in `ea7c209` |
| `--batch_type weighted` (documented) | 1 | `ValueError: min() iterable argument is empty`. *Finding 5* |
| `--use_old_sorted_file` with no `sorted.fastq` | 1 | `UnboundLocalError: read_array` — **but only the first time**: the failed run leaves an empty `sorted.fastq` and the retry exits 0. *Finding 23* |
| fastq with no trailing newline | 1 | `TypeError: 'NoneType' object is not iterable`. *Finding 12* |
| `--consensus` with neither `--medaka` nor `--racon` | 1 | `UnboundLocalError: polishing_pattern`, **after** forming and merging every consensus. *Finding 4* |
| `--consensus --max_seqs_for_consensus 0` | 1 | `CalledProcessError`: spoa died with `SIGABRT` on an empty input file. *Finding 22* |
| `--use_old_sorted_file --k 25` after sorting at `--k 13` | **1 on `smoke`, 0 on `sup`** | `ValueError: not enough values to unpack (expected 8, got 6)`. It needs a short read, and only the smaller corpus has one. *Finding 18* |
| `write_fastq` on a fastq with spaces in headers | 1 | `KeyError` on a truncated accession, after writing some files. *Finding 6* |

**Prefix matching is three behaviours, not two**, and all three are pinned:

| behaviour | case | result |
| --- | --- | --- |
| a prefix matching one flag is accepted | `--outfold` → `--outfolder` | exit 0 |
| an **exact** match wins over longer flags sharing it | `--m` is `--m` (`target_length`), not ambiguous with `--min_shared`/`--mapped_threshold`/`--medaka…`; `--s` is `--s`, not `--sample_size` | exit **1** — it parses, then dies on the missing `--outfolder` |
| a prefix matching several flags is rejected | `--me` → `--medaka`, `--medaka_model`, `--medaka_fastq`; `--min` → three; `--r` → four; `--prim` → two | exit 2, with the candidates listed |

The middle row is the one that surprises: `--m 5` is *not* an error. A port that treats `--m` as an
ambiguous prefix of the six `m` flags rejects a valid invocation.

And clap rewrites `field_name` to `--field-name`, so **all 20 multi-word flags need an explicit
`long = "..."`**: `--abundance_ratio`, `--aligned_threshold`, `--batch_type`, `--mapped_threshold`,
`--max_seqs_for_consensus`, `--medaka_fastq`, `--medaka_model`, `--min_fraction`,
`--min_prob_no_hits`, `--min_shared`, `--primer_file`, `--primer_max_ed`, `--racon_iter`,
`--rc_identity_threshold`, `--remove_universal_tails`, `--sample_size`,
`--symmetric_map_align_thresholds`, `--top_reads`, `--trim_window`, `--use_old_sorted_file`. The
**eight** double-dash single-letter options — `--N`, `--d`, `--k`, `--m`, `--q`, `--s`, `--t`, `--w` —
are not what clap's `short` produces either. Note `--m`/`--s`/`--t`/`--d`/`--q` also carry a `dest`
that differs from the flag (`target_length`, `target_deviation`, `nr_cores`, `print_output`,
`quality_threshold`).

### The case matrix

24 cases were swept on both committed corpora, at `--t 1`, and the question asked of each was: **do
any two cases produce the same `final_clusters.tsv`?** A case that duplicates another is not a test.

First on the 24 pre-consensus cases, comparing `final_clusters.tsv` alone:

| corpus | reads | cases | distinct | crashes | unintended collisions | wall clock |
| --- | --- | --- | --- | --- | --- | --- |
| `sample_h1.fastq` | 280 | 24 | 12 | 1 | **8** | 12 s |
| **`Supplementary_File1_reads.fastq`** | 3 000 | 24 | **19** | 1 | **1** | 38 s |

`sample_h1` gives one group of **11** identical results, one of 2, and ten singletons; the 3 000-read
corpus gives one group of 4, one of 2, and seventeen singletons.

Then on the full recorded matrix, comparing **each case's whole set of output files**. That is the
right comparison and `final_clusters.tsv` alone is not, because `--consensus` does not change the
clustering: sixteen cases share one `final_clusters.tsv` on the 3 000-read corpus and are
distinguished entirely by their consensus, polishing and trimming output.

| corpus | cases | **distinct output sets** | collision groups |
| --- | --- | --- | --- |
| `smoke` | 51 | **30** | 5 |
| **`sup`** | 51 | **42** | 4 |

`sup`'s four, all of them explicable:

| group | why |
| --- | --- |
| `default` == `ont` == `q0` == `top_huge` | `default`/`ont` is **intended** — it pins that the preset resolves to `--k 13 --w 20`. `q0` and `top_huge` are properties of the data: no read has mean quality between 0 and 7, and `--sample_size 999999` exceeds the read count so the subsample never fires |
| `t8` == `t8_total_nt` | **intended** — pins that `total_nt` is the default `--batch_type` |
| `cons_racon` == `cons_racon_rc1` == `cons_primer_ed0` == `cons_primer_tw50` | `--rc_identity_threshold 1.0` merges nothing, which on this data is what 0.9 already did; `--primer_max_ed 0` and `--trim_window 50` both find no primer, so both equal plain racon. Real information: they are the negative controls for the two cases that *do* trim |
| `wf_N0` == `wf_N2` == `wf_N10` | all three crash identically, before `--N` is ever read. *Finding 6* |

`smoke`'s five include everything above **plus** seven consensus cases collapsing into one (it has a
single abundant cluster and no primers are found in it) and `m750s50` == `m800s100`. Twelve of its
cases collapse onto `default`, `--symmetric_map_align_thresholds` among them.

Two collisions are intended on both and must stay: `default` == `k13w20` pins the defaults, and
`ont` == `k13w20` pins that the preset resolves to `--k 13 --w 20`. What the small corpus hides is
*Finding 19*, and it includes the flag this port has to write from scratch.

The matrix still needs extending in three directions that were not swept: `--t` at 2/4/8 (measured
separately — 42/35/33 clusters), all three `--batch_type` values, and every `--consensus` case. It
also needs `--q 8` and `--q 9`, which discriminate cleanly (31 and 14 clusters), rather than `--q 15`,
which crashes.

### The corpora

Only two exist and both are committed. That is a smaller registry than any of the other three ports
had, and it is a gap, not a virtue.

| corpus | reads | truth | what it is for |
| --- | --- | --- | --- |
| `sample_h1.fastq` | 280 | one haddock amplicon, `h1` in every header | committed, 390 KB. The README install check, `.travis.yml`, CI. Median length 632, and one read of length 14 that exercises the `2*k` filter |
| `Supplementary_File1_reads.fastq` | 3 000 | three species — cod, haddock, whiting — tagged `c1`/`h1`/`w1` in every header | committed, 5.2 MB. **Develop against this.** Median length 816. Discriminates 19 distinct results out of 24 cases, reaches the reverse-complement merge path (4 draft centers, one merged), and reaches primer trimming with the committed IUPAC primer file |

**One more is available for free.** `test/terncytb_200.fastq` — 200 real ONT reads of a Cytb/16S
amplicon, median length 428 — was committed in `a2128c8` and dropped in `11c516f`, and it is 214 KB.
It is the only strippable blob in this history that is **not** already archived from the isONclust
exercise (*Repo hygiene*), so it has to be handled either way, and restoring it as a third corpus
costs less than archiving it. It adds a length axis (428 / 632 / 816) and a third amplicon. It does
not discriminate `--symmetric_map_align_thresholds`.

**What is missing, in priority order:**

1. **A PacBio corpus.** `--isoseq` is a supported preset with no data behind it. It resolves to
   `--k 15 --w 50`, which the sweep does exercise on ONT reads, but nothing checks the path on data it
   was designed for.
2. **A corpus with a quality spread.** `--q 0` and the default `--q 7` produce identical output on
   both corpora, so the quality filter is only observable at `--q 8` and above, where it is already
   discarding most of the data.
3. **A depth axis.** 3 000 reads is one point. Method rule: *a conclusion at one depth is not a
   conclusion*. `--sample_size`, `--abundance_ratio` and `--t` all interact with depth and none of
   them has been measured at more than one.
4. **Reads with primers actually present at both ends in known positions**, so primer trimming can be
   checked against an answer rather than against the reference's own output.

`bench/corpora.tsv` should register these as they arrive; goldens are per-corpus and the manifest
records the corpus sha256, so adding one is cheap and changing one is not.

### Stage-level oracles are required, not optional

*Finding 19* is the argument in miniature — `sample_h1` cannot observe five parameters, one of which
is the only new clustering logic in the port — and the isONclust port's demonstration is the argument
in full: two plausible ways to get `get_kmer_minimizers` wrong were introduced deliberately, both
produced *exactly the same number of minimizers*, and only a full ordered-list diff caught them.

`bench/dump_reference.py` wraps the reference **without modifying it** and writes each stage's inputs
*and* outputs in a stable line format, to be replayed from Rust and diffed directly. Stages 1–9 are
isONclust's and are already verified there. The stages this port adds:

10. **the length filter and the subsample** — the surviving `read_array` indices, per `(--m, --s,
    --sample_size, --top_reads)`. Cheap, and it is where *Finding 1* lives
11. **`get_best_cluster` and `get_best_cluster_block_align` under `--symmetric_map_align_thresholds`**
    — the per-candidate pair of ratios, the `min`, and the winner. The isONclust oracle covers the
    default gate only
12. **`form_draft_consensus`** — the exact sequence list handed to spoa, in order, per cluster, and
    the consensus returned. This is the POA oracle and it is the one that decides whether `spoars`
    can be used
13. **`highest_aln_identity`** — both identities, forward and reverse-complement, per pair. Two
    parasail calls per pair at `opening_penalty=3`
14. **`detect_reverse_complements`** — the merge decisions and the resulting `centers` list. Stateful
    and order-dependent; diff, do not replay
15. **`find_barcode_locations`** — every `edlib.align` call's arguments and its full `locations` list,
    not just `locations[0]`, so the tie-break is measurable rather than guessable
16. **`remove_barcodes`** — the cut positions and the trimmed center, per center

End-to-end equivalence tells you *that* the port is wrong, never *where*. And **dump from the live
driver too, not only from a standalone dump binary**: oracles replay recorded reference inputs, so a
port following a different trajectory still passes them. The one real bug in the isONcorrect port
passed every oracle and was caught only by diffing a dump taken from the running driver.

## Findings in the reference

Ordered by how much they matter to the port.

### Finding 1 — `--sample_size` is not reproducible, and the README recommends it

`main` subsamples with

```python
read_array = [read_array[i] for i in sorted(random.sample(range(len(read_array)), args.sample_size))]
```

`random` is never seeded, there is no `--seed` flag, and CPython seeds the global Mersenne Twister
from OS entropy at import. So the draw is different every run.

Measured on `Supplementary_File1_reads.fastq` with `--sample_size 500 --t 1` at a **fixed**
`PYTHONHASHSEED=0`, five runs of the same command:

| run | `final_clusters.tsv` + `final_cluster_origins.tsv` sha256 (first 16) |
| --- | --- |
| 1 | `256de88f62b817f7` |
| 2 | `b6f6f7dbd7982da5` |
| 3 | `0b6f042529f2a1ea` |
| 4 | `1e9f091b0c86f252` |
| 5 | `82a33afd4c91af33` |

Five runs, five different answers. With `--top_reads` added, three runs gave one answer.

This is not a corner. The README's own worked example is
`NGSpeciesID --ont --consensus --sample_size 500 --m 800 --s 100 --medaka --primer_file primers.txt`,
`test/consensus.sh` ships that command in a loop over every fastq in a directory, and the protocol
manuscript the README describes is built around it. **The tool's recommended invocation cannot
reproduce its own results**, and every consensus sequence downstream of it inherits that.

**This is the one decision the port cannot take on its own.** Four options, with what each costs:

| option | byte-identity | user impact |
| --- | --- | --- |
| **(a)** add `--seed` to the Python, defaulting to a fixed value | achievable, goldens recorded once | changes today's output for every `--sample_size` user, and makes it reproducible. Recommended |
| **(b)** add `--seed` defaulting to entropy, i.e. opt-in reproducibility | achievable only with `--seed` set | no change to existing behaviour; nobody who does not read the changelog benefits |
| **(c)** reproduce Mersenne Twister and `random.sample` exactly in Rust, and add `--seed` | achievable | faithful, and a genuine amount of work for a defect |
| **(d)** exclude the random path from the contract; port `--top_reads` only | not achievable for that flag | leaves the recommended invocation unverifiable |

(a) is recommended, it is a fix to the *Python* under method point 8, it is one commit, and it is a
blocker: goldens for any `--sample_size` case are meaningless until it lands. Note that (a) and (c)
are not exclusive — (c) is what makes (a)'s goldens portable across implementations, and it can come
later.

Note also that `--top_reads` and `--m`/`--s` are reproducible today and can be ported and verified
immediately. Only the random draw is blocked.

### Finding 2 — the reference does not agree with itself on Python ≤3.11, and the README recommends 3.11

Four sites sum over a `set` of quality characters:

```python
poisson_mean = sum([qual.count(char_) * D_no_min[char_] for char_ in set(qual)])
```

Set iteration order is `PYTHONHASHSEED`-dependent, and floating-point addition is not associative.
From CPython 3.12 the interpreter's `sum()` is compensated and the result is order-independent;
before 3.12 it is a naive left-fold and it is not.

Measured, `PYTHONHASHSEED=random`, `--t 1`:

| interpreter | corpus | runs | distinct results |
| --- | --- | --- | --- |
| 3.12.14 | `sample_h1` | 5 | **1** |
| 3.12.14 | `sample_h1`, `--t 8` | 3 | **1** |
| 3.12.14 | `Supplementary_File1` | 3 | **1** |
| 3.11.16 | `sample_h1` | 6 | **6** |
| 3.11.16 | `Supplementary_File1` | 4 | **4** |

Which files move, on 3.11: `final_cluster_origins.tsv` and `logfile.txt`, both times, in the last few
digits of the error-rate column — `0.09274685314674982` against `0.09274685314674987`.
`final_clusters.tsv` and `sorted.fastq` happened to be stable on both corpora, and the cluster count
was 21 in all six `sample_h1` runs. **That is luck, not a property.** `error_rate` feeds
`p_shared_minimizer_empirical` through `round(e1, 2)`, so a last-ULP difference can cross a rounding
boundary and change a mapping decision; it did on some isONclust corpora.

**Decision taken: pin the reference to Python ≥3.12 and leave the Python alone**, the same decision
and the same reason as the isONclust port. The port targets 3.12 semantics — exactly rounded
summation — which frees it to iterate in any order, and the goldens are valid only for that
interpreter, which the manifest must record.

But this repository has something isONclust did not: **the README tells users to install Python 3.11**
("Recent update (2025-04-19)"). That instruction should say 3.12 or later, and saying so is repo
hygiene rather than port work. The one-line `math.fsum` fix that would make older interpreters agree
is written up under *Deferred improvements* and not applied.

### Finding 3 — the installation instructions do not work on aarch64, and do work on x86_64 Linux

Measured in *Goal*, and repeated here because it is a finding about the reference and not only a
motivation. **The scope matters and is easy to overstate:** both README recipes resolve on
`linux-64`, and `pip install NGSpeciesID` gets a wheel there. What fails is
`osx-arm64` (both recipes, plus parasail), `osx-64` (both recipes) and `linux-aarch64` (parasail).
The path that works everywhere is undocumented, and `parasail` — which `setup.py` pins at `1.2.4` —
publishes **no aarch64 wheel at any version for any OS**, so on ARM it always builds from source and
that build fails.

**Both are fixed**, on `fix/installation` commit `76d3f88`, because they are one-line problems in a
tool people are trying to install today and there is no reason to make them wait for a port:

1. The README now carries the recipe that works, plus `--no-deps`, plus the reasons — see *Goal*.
2. `setup.py` and `requirements.txt` no longer pin `parasail==1.2.4`. The dead 3.4–3.7 classifiers
   are gone and `python_requires` is `>=3.10`, the oldest interpreter bioconda still builds
   `parasail-python` and `python-edlib` for.

The README also now documents *Finding 1* (`--sample_size` is an unseeded draw, `--top_reads` is the
reproducible alternative), *Finding 2* (python ≤3.11 makes the reported error rates unstable) and
*Finding 20* (`--abundance_ratio` applies after subsampling). Documenting a defect is not fixing it —
all three stay open — but it stops users being surprised by it in the meantime.

### Finding 4 — `--consensus` without a polisher crashes, after doing all the work

`polish_sequences` sets `polishing_pattern` inside `if args.medaka: ... elif args.racon: ...` and then
reads it unconditionally. Neither flag is gated on `--consensus`, and neither is required by it. So:

```
NGSpeciesID --ont --fastq test/sample_h1.fastq --outfolder out --consensus
```

runs the clustering, runs spoa on every abundant cluster, runs reverse-complement detection, and then
dies with `UnboundLocalError: cannot access local variable 'polishing_pattern'`, exit 1. Everything
expensive has already happened and none of it is written out.

The mirror case is silent: `--medaka` or `--racon` **without** `--consensus` exits 0 and does nothing
at all — no consensus, no warning.

The port must reproduce the exit code. Whether it reproduces the traceback is a separate question, and
the fix — treat `--consensus` without a polisher as "draft consensus only", which is what a user
plainly means, and reject a polisher without `--consensus` — belongs in *Deferred improvements*.

### Finding 5 — `--batch_type weighted` is documented and crashes, and so does every typo

`--help` says: `how to split the reads into chunks "total_nt", "nr_reads", or "weighted"`. There is no
`weighted` branch in `batch_list`, and no `else`. An unrecognised `batch_type` makes the generator
yield nothing, and the first thing downstream of it is a `min()` over an empty list:

```
--batch_type weighted --t 4  ->  ValueError: min() iterable argument is empty, exit 1
--batch_type bogus   --t 4  ->  the same
```

isONclust has the identical defect in the identical function and it **does not crash there**, because
isONclust guards all three `yield batch` sites with `if batch:` and NGSpeciesID does not. The two
repositories fixed the same underlying problem in two different places: isONclust stopped yielding
empty batches; NGSpeciesID wrote `min(prev_b_indices or [1])` in `cluster.py`. NGSpeciesID's guard
covers one call site and not this one.

Contract: the port crashes here too, non-zero, until a commit says otherwise. The obvious fix —
validate `--batch_type` in the argument parser, which is where a bad enum belongs — is deferred, and
it is worth noting that it would turn an exit-1 traceback into an exit-2 argparse error, so it is a
contract change and needs its own case.

### Finding 6 — `write_fastq` is unreachable, and broken when reached

Two independent defects.

**It cannot be invoked.** The `--fastq`/`--use_old_sorted_file` mutually-exclusive group is
`required=True` on the *top-level* parser, so `NGSpeciesID write_fastq --clusters ... --fastq ...`
exits 2 asking for a top-level `--fastq`. The subcommand's own `--fastq` does not satisfy it. It only
runs if you pass a redundant top-level `--fastq` (or `--use_old_sorted_file`) as well.

**When it runs it fails on every ONT fastq.** `write_fastq` parses `final_clusters.tsv` with

```python
items = line.strip().split()
cl_id, acc = items[0], items[1]
```

`.split()` with no argument splits on **all** whitespace, and column 2 of `final_clusters.tsv` is a
full ONT accession containing spaces:

```
0	c948601e-1bd9-4039-94c5-3d8c741df65f runid=ed1de13037eb1bcfeefc20d46af38743bdac68e6 read=4709 ch=45 ...
```

So `acc` becomes the bare UUID, while `reads` is keyed by the whole header. Measured: it writes
`0.fastq` and then dies with `KeyError: 'c948601e-...'`, exit 1. It works only on a fastq whose
headers contain no spaces — which is why isONclust, whose `readfq` substitutes them, never sees this.

The fix is `line.rstrip("\n").split("\t", 1)`. It is one line, it makes a currently-dead feature
work, and it changes no output that anyone can currently obtain — so it is close to free. It is still
a behaviour change and needs its own commit and its own cases.

**Decide before porting:** fix it upstream first, or port it broken. Porting it broken means
reproducing a `KeyError` traceback, which is not a contract worth having.

### Finding 7 — `--q` above the corpus quality crashes with `IndexError`

```
NGSpeciesID --ont --fastq test/Supplementary_File1_reads.fastq --outfolder out --t 1 --q 12
  -> IndexError: list index out of range   at   min_e = error_rates[0]
```

Every read is filtered out, and the statistics block indexes an empty list. `--q 11` leaves one read
and works; `--q 12` is a perfectly reasonable thing to ask for.

**This is isONclust's *Finding 9*, and isONclust fixed it in its own Python** — a guard that writes an
explanatory line to the logfile, says "no reads passed the quality filter" on the log's error stream,
and exits 1, leaving the exit status where the traceback already had it.

**Taken.** `fix/installation` commit `ea7c209`. Verified a no-op across the whole 24-case parameter
matrix on `Supplementary_File1_reads.fastq`: 24 of 24 `final_clusters.tsv` byte-identical before and
after, and the one case that reaches the guard (`--q 12`) now prints two lines instead of a traceback
and still exits 1. Also checked on the `--t 8` parallel sort path, which has its own copy of the
filter loop.

### Finding 8 — `--k` outside 10–30 crashes on an empty probability table, and twelve settings inside the range do too

`p_minimizers_shared.L` covers `k` 10–30 and, for each `k`, `w` at `k + 5n` up to 100. `main` keeps
rows with `int(k) == args.k and abs(int(w) - args.w) <= 2`, and if nothing survives, `p_emp_probs` is
empty and the first lookup in `p_shared_minimizer_empirical` raises:

```
--k 9  --w 20  ->  KeyError: (0.08, 0.08), exit 1
--k 31 --w 50  ->  KeyError: (0.09, 0.08), exit 1
```

`--k` is validated for nothing but `--w >= --k` and `--w <= 100`, so **every `k` from 1 to 9 and from
31 up is CLI-valid and crashes** — 3 349 `(k, w)` pairs. Inside 10–30, twelve more do:

| `k` | `w` with an empty table |
| --- | --- |
| 11, 16, 21, 26 | 99, 100 |
| 12, 17, 22, 27 | 100 |

isONclust has the same defect (its *Finding 13*, fifteen settings) but a **larger table**: PR #13
added `k` 4–9, taking it to 59 628 rows. NGSpeciesID's 41 880 rows are exactly isONclust's `k >= 10`
subset — same values, same order — so the port can carry isONclust's frozen blob across and **must
restrict it to `k >= 10`**. Shipping the whole thing would silently make `--k 9` work, which is a
behaviour change wearing a data change's clothing.

The right fix is to validate `--k` and `--w` against the table at argument-parse time and say what
the valid range is. Deferred, and it changes exit 1 to exit 2.

### Finding 9 — `--outfolder` is optional in argparse and mandatory in fact

`--outfolder` has `default=None`. Omitting it reaches `os.path.join(args.outfolder, "sorted.fastq")`
and raises `TypeError: expected str, bytes or os.PathLike object, not NoneType`, exit 1. There is no
invocation of the main path in which it is optional.

Contract for now. `required=True` is the fix, and it moves the exit code from 1 to 2.

### Finding 10 — reverse-complement detection double-counts reads, and polishes one cluster twice

`detect_reverse_complements` skips a center that is in `already_removed` in the **outer** loop and
does not check it in the **inner** one. So a center merged into an earlier center can be merged again
into a later one.

Constructed and measured, with three centers where `identity(A,C) = 0.937`, `identity(B,C) = 0.937`,
`identity(A,B) = 0.873` and a threshold of 0.9:

| | in | out |
| --- | --- | --- |
| A | 10 reads | **13** reads, `[a.fq, c.fq]` |
| B | 5 reads | **8** reads, `[b.fq, c.fq]` |
| C | 3 reads | merged, twice |
| total | **18** | **21** |

Two consequences. The reported `total_supporting_reads` in the consensus fasta header is inflated —
it is written straight from `merged_nr_reads`. And `c.fq` appears in both centers' `all_reads`, so
cluster C's reads are handed to the polisher twice, pulling two different consensus sequences toward
the same data.

The condition needed is a genuine one — a center similar to two mutually dissimilar centers — so how
often it fires on real data is unmeasured and worth measuring. It did not fire on
`Supplementary_File1_reads.fastq`, where 4 draft centers produced one merge and 3 consensus.

Reproduce it. The fix is `if c_id2 in already_removed: continue` in the inner loop, and it belongs in
*Deferred improvements* with a measurement of how much it moves real output.

### Finding 11 — `--d 0` divides by zero

`if i % args.print_output == 0:` is evaluated before anything checks that `print_output` is non-zero,
and the guard that was meant to protect it (`if args.print_output:`) only wraps the header line
above. `--d 0` raises `ZeroDivisionError`, exit 1. Identical to isONclust's *Finding 2*. It has a
golden.

### Finding 12 — a fastq without a trailing newline crashes

`readfq` yields `(name, (seq, None))` when it hits EOF before reading enough quality characters, and
`expected_number_of_erroneous_kmers` iterates that `None`:

```
TypeError: 'NoneType' object is not iterable, exit 1
```

The same two-read file with a trailing newline works. Identical to isONclust's *Finding 11*, where the
port reproduces it and has a unit test for the path.

### Finding 13 — `--use_old_sorted_file` truncates the logfile

Measured: an output folder with a 171-byte `logfile.txt`, re-run with `--use_old_sorted_file`, comes
back with a **0-byte** `logfile.txt`. The clustering is identical; the log is destroyed. The file is
opened for writing before the branch that decides not to write anything to it. Identical to
isONclust's *Finding 8*.

This matters to the harness more than to users: `logfile.txt` is in the byte-identity contract, and a
`--use_old_sorted_file` case whose golden is an empty file is a case that cannot fail.

### Finding 14 — `--ont --isoseq` together exits 0, and `--ont` silently overrides `--k`/`--w`

```python
if args.ont and args.isoseq:
    logging.error("Arguments mutually exclusive, ...")
    sys.exit()
```

`sys.exit()` with no argument exits **0**. A pipeline cannot detect the mistake. The two flags are not
in an argparse mutually-exclusive group even though `--medaka`/`--racon` and
`--remove_universal_tails`/`--primer_file` both are.

Separately, the presets are applied *after* argument parsing and unconditionally overwrite `--k` and
`--w`, so `--ont --w 5` runs at `--w 20` — measured, 21 clusters, identical to plain `--ont`. Both
behaviours are contract and both need a case.

### Finding 15 — two pieces of dead code that look like safety checks

```python
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
```

Unreachable: the required `--fastq`/`--use_old_sorted_file` group makes argparse exit 2 first. It is
also placed *after* the `--ont`/`--isoseq` check, which would have to run first anyway.

```python
parasail_module = 'parasail'
if parasail_module not in sys.modules:
    logging.error('You have not imported the {0} module. ...')
    sys.exit(1)
```

Always false: `modules.consensus` imports `parasail` at module scope, so by the time this runs the
name is in `sys.modules`. It cannot fire, and if the import had failed the program would already have
died at line 16. The port needs neither. Both should be deleted from the Python — they are noise, not
behaviour — and that is a `master` commit.

### Finding 16 — `read_barcodes` keeps the fasta description in the key, and uppercases only half the primer

```python
barcodes = {acc + '_fw': seq.strip() for acc, (seq, _) in readfq(open(primer_file))}
for acc, seq in list(barcodes.items()):
    barcodes[acc[:-3] + '_rc'] = reverse_complement(seq.upper())
```

`readfq` returns the whole header after `>`, so the committed primer file — whose first header is
`>COIF-ALT ` with a trailing space — produces the keys `'COIF-ALT _fw'` and `'COIF-ALT _rc'`. Measured
and confirmed in the debug log. It reaches nothing but log output today, and it would reach output if
the barcode name were ever written into a filename or a fasta header.

More substantively: `seq.strip()` is stored **as written** for the forward primer while the reverse
complement is built from `seq.upper()`. `reverse_complement`'s table handles lowercase, so a lowercase
primer file yields a lowercase forward primer and an uppercase reverse complement — and edlib's
`additionalEqualities` IUPAC map is uppercase-only, so IUPAC codes in a lowercase primer stop
matching. The committed primer file is uppercase, so this is unexercised; it is still a real
asymmetry and the port must reproduce it.

### Finding 17 — the IUPAC equivalence map is one-directional

`additionalEqualities` contains `('M','A')` and `('M','C')` but not `('A','M')`. edlib treats the
pairs as ordered `(query, target)`, so an ambiguity code in the **primer** matches a concrete base in
the **center**, and a concrete base in the primer does not match an ambiguity code in the center.
Given that centers are spoa consensus sequences over ACGT that is the right way round — but it is
data, it must be copied verbatim, and the port must not "symmetrise" it.

The map also includes `('X','G')` through `('X','C')`, i.e. `X` as a full wildcard, which is not in
IUPAC and comes from the Bio.Data table the comment cites.

### Finding 18 — a skipped read keeps a 6-element tuple that the output writer unpacks as 8

`reads_to_clusters` skips reads whose homopolymer-compressed length is under `--k` with a bare
`continue`, before the block that grows `representatives[read_cl_id]` to eight elements. The read
stays its own cluster with a six-element tuple, and `main`'s output loop does

```python
read_cl_id, b_i, acc, c_seq, c_qual, score, error_rate, _ = representatives[c_id]
```

Within a single run this is unreachable, because the sort stage already drops those reads using the
same `--k`. Across two runs it is two ordinary commands:

```bash
NGSpeciesID --fastq test/sample_h1.fastq --outfolder out --t 1 --k 13 --w 20
NGSpeciesID --use_old_sorted_file        --outfolder out --t 1 --k 25 --w 50
#   -> ValueError: not enough values to unpack (expected 8, got 6), exit 1
```

Also reachable with a hand-written `sorted.fastq` containing a read of 60 bases in 5 homopolymer runs
— `len(seq) >= 2*k` passes, `len(hpol) < k` does not.

**And it is the one place the corpora invert.** The read has to survive the `--k 13` sort filter
(`len(seq) >= 26`, `len(hpol) >= 13`) and fail cluster.py's `--k 25` guard (`len(hpol) < 25`), which
means it has to be short. Counted:

| corpus | reads | reads that reach *Finding 18* | exit |
| --- | --- | --- | --- |
| `sample_h1` | 280 | **4** | 1, `ValueError: not enough values to unpack (expected 8, got 6)` |
| `Supplementary_File1_reads.fastq` | 3 000 | **0** | 0, 334 clusters |

So the corpus that is blind to five parameters (*Finding 19*) and finds no primers at all
(*Finding 25*) is the only one that reaches this crash — its shortest read is 14 bases against the
other's 183. **Keep both, and record this case's golden on both**: it is the clearest illustration in
the project that "the bigger corpus" is not the same thing as "the better corpus".

Contract. The underlying problem is that `--use_old_sorted_file` reuses a file sorted under different
parameters with no record of what they were; writing `--k`/`--w`/`--q` into the logfile and refusing
a mismatch is the fix, and it is deferred.

### Finding 19 — the small corpus is blind to five parameters, including the only new clustering logic

The 24-case sweep on `sample_h1.fastq` produced 12 distinct results. Eleven cases collapsed onto one
hash:

| collapsed onto `default` | why it hides something |
| --- | --- |
| `--mapped_threshold 0.3`, `--mapped_threshold 0.95` | the mapped fraction is far from both thresholds on 280 short reads |
| `--aligned_threshold 0.1` | the alignment gate never binds downward |
| `--min_fraction 0.5`, `--min_fraction 1.0` | ditto |
| `--min_prob_no_hits 0.01`, `--min_prob_no_hits 0.5` | ditto |
| **`--symmetric_map_align_thresholds`** | **the flag this port has to write from scratch produces identical output** |

On `Supplementary_File1_reads.fastq` every one of those is distinct, and
`--symmetric_map_align_thresholds` gives **85 clusters against the default's 49** — a large, obvious
effect that the smaller corpus cannot see at all — and combining it with `--aligned_threshold 0.9`
gives a third result distinct from both it and `aligned0.9` alone, so the flag's interaction with the
alignment gate is covered too. `--q 0` is the one remaining unintended collision
there, and it is a property of the data (no read has mean quality between 0 and 7) rather than of the
sweep. `--q 8` and `--q 9` do discriminate, at 31 and 14 clusters, and should replace `--q 0` and
`--q 15` in `cases.tsv`.

**Develop against `Supplementary_File1_reads.fastq`. Keep `sample_h1.fastq` as the smoke fixture and
the README check, and do not mistake it for coverage.** This is the isONclust port's *Finding 5*
arriving before any code was written, which is the only good time for it.

### Finding 20 — `abundance_cutoff` is computed after subsampling, and floors to zero

`abundance_cutoff = int(args.abundance_ratio * len(read_array))` runs on the **post-subsample**
`read_array`, so `--sample_size 500 --abundance_ratio 0.1` means "50 reads", not "10% of the file".
That is defensible and probably intended, but it means `--sample_size` and `--abundance_ratio`
interact, and it is not documented.

`int()` truncates, so at fewer than 10 reads the default `--abundance_ratio 0.1` gives a cutoff of 0.
Measured with `--sample_size 5 --top_reads`:
`Forming draft consensus with abundance_cutoff >= 0 (10.0% of 5 reads)`. `nr_reads_in_cluster >= 0`
is then always true, so every cluster including singletons goes to spoa and the
`elif nr_reads_in_cluster == 1` singleton counter becomes unreachable — which is why the same run
reports `0 singletons were discarded`.

### Finding 21 — seven blobs live at more than one path, so `git rev-list --objects` cannot build the removal list

`git rev-list --objects --all` lists each object **once**, with one arbitrary path. Building a removal
list from it looks correct and silently misses aliases. Measured here by walking every commit's tree
instead — 309 distinct `(blob, path)` pairs over 302 blobs, and **seven blobs at two paths each**:

| blob | paths |
| --- | --- |
| `d78cdf85` | `test/Supplementary_File1_reads.fastq`, `test/Supplementary_File2_reads.fastq` |
| `5008ddfc` | `modules/.DS_Store`, `scripts/.DS_Store` |
| `94037afc` | `NGSpeciesID`, `isONclust` |
| `c964fb0d` | `cemetary/cluster_parallel.py`, `modules/cluster_parallel.py` |
| `d231d147` | `modules/compute_shared_minimizers_probabilities.py`, `scripts/compute_shared_minimizer_probabilities.py` |
| `1434c00b` | `test/Supplementary_File3_primer.txt`, `test/Supplementary_File4_primer.txt` |
| `ef1deb95` | `test/Supplementary_File2_minibar.txt`, `test/Supplementary_File3_minibar.txt` |

`git rev-list --objects` reported only one path for each. This is exactly the trap isONclust's
`tools/repo-slim/analyze.sh` was rewritten to avoid, and it is live here: `scripts/.DS_Store` would
not have appeared on a naive removal list. **Build the list from a tree walk.** Checked, and no blob
straddles the boundary between a stripped path and a retained one, so the removal is clean — but that
is a measurement, not an assumption.

### Finding 22 — `--max_seqs_for_consensus 0` makes spoa abort

The cutoff is `if args.max_seqs_for_consensus >= 0 and i >= args.max_seqs_for_consensus: break`, with
`i` starting at 0. At `--max_seqs_for_consensus 0` the break fires immediately, the per-cluster fastq
is written empty, and spoa is handed it:

```
subprocess.CalledProcessError: Command '['spoa', '.../reads_c_id_17.fq', '-l', '0', '-r', '0', '-g', '-2']'
    died with <Signals.SIGABRT: 6>.        exit 1
```

Note the `>=` also means the flag admits exactly `max_seqs_for_consensus` sequences, not
`max_seqs_for_consensus + 1` — the opposite of isONcorrect's `--max_seqs_to_spoa`, whose cutoff is a
bare `>`. The default is `-1`, which disables the cutoff entirely, and `nr_reads_in_cluster` is
separately floored by `abundance_cutoff`, so only the literal `0` reaches this.

Contract: exit non-zero. A native POA has to decide what it does with an empty input, and "abort the
process" is not a behaviour worth reproducing faithfully — so this is one of the few places where the
port's `spoars` path will need an explicit guard *and* a matching exit code, and it needs its own
case. Rejecting `--max_seqs_for_consensus 0` in the parser is the fix and it is deferred.

### Finding 23 — a failed `--use_old_sorted_file` run poisons the folder, and the retry exits 0

`get_sorted_fastq_for_cluster.main` opens `sorted.fastq` for writing *before* the branch that would
have filled it, so a run that fails leaves an **empty** `sorted.fastq` behind. The next run in the
same folder then takes the "use the existing sorted file" path:

```
$ NGSpeciesID --use_old_sorted_file --outfolder out --t 1 --ont
UnboundLocalError: cannot access local variable 'read_array'        exit 1
$ ls out
logfile.txt  sorted.fastq        # both zero bytes

$ NGSpeciesID --use_old_sorted_file --outfolder out --t 1 --ont
Using already existing sorted file in specified directory, ...
Starting Clustering: 0 reads
Finished Clustering: 0 clusters formed                              exit 0
```

So a pipeline that retries after a failure gets a **clean exit and an empty `final_clusters.tsv`**,
which is worse than the failure it was retrying. The same shape reaches any workflow that reuses an
output folder.

**Found by the harness, not by reading the code.** `bench/equivalence.sh stable` records the whole
matrix twice in one process, and this case's exit code was 1 the first time and 0 the second — which
is precisely the class of defect that check exists for, and which no single recording could see.

Contract for the port. The fix is to open `sorted.fastq` only on the branch that writes it, and it is
one line; deferred.

### Finding 24 — medaka's BAM embeds absolute paths, so it cannot be a byte-identity target

`medaka_cl_id_<id>/calls_to_draft.bam` carries minimap2's and samtools' own command lines in its
`@PG` header, and those command lines contain the **absolute path** of the output folder:

```
@PG ID:minimap2  ... CL:minimap2 -x map-ont ... /tmp/xyz/consensus_reference_17.fasta /tmp/xyz/reads_to_consensus_17.fastq
@PG ID:samtools  ... CL:samtools view -@ 1 -T /tmp/xyz/consensus_reference_17.fasta -F 2308 -bS -
```

Every harness case runs in a fresh temporary directory, so the file can never match across two runs
and it is excluded from the contract along with its `.bai`. Its *existence* still is.

Two things worth noting rather than assuming:

1. `consensus_probs.hdf` is **not** excluded. It was measured stable across the same comparison, and
   excluding a file because it is the kind of file that might vary is how a contract quietly stops
   covering anything.
2. This was found by running `verify` against a "port" that **was** the reference — 48 of 51 cases
   passed and the three medaka cases failed on exactly these two files. Nothing in the reference or
   the goldens would have said so.

### Finding 25 — the smoke corpus finds no primers at all, so the barcode oracle is vacuous on it

*Finding 19* measured the case matrix. This measures the stage oracles, and the answer is worse.
`bench/dump_reference.py` was run on both committed corpora:

| stage | `sample_h1` (280 reads) | `Supplementary_File1_reads.fastq` (3 000) |
| --- | --- | --- |
| `minimizers` | 25 433 minimizers | 374 918 |
| `mapping` | 928 recorded decisions | 27 215 |
| `parasail` | 121 alignments | 960 |
| `spoa` | 255 lines, 2 POA calls | 2 766 |
| `identity` | **1** center pair (RC orientation won) | 6 pairs, RC won 4 |
| `barcode` | 16 edlib calls, **0 found a primer** | 32 calls, **6 found a primer** |

So on `sample_h1` the `barcode` oracle exercises the edlib call and **never the cut-position logic**
— `remove_barcodes` computes no cut at all — and `identity` has exactly one pair, which cannot
distinguish an ordering bug from a correct implementation. `equivalence.sh stage` prints these counts
for that reason: a pass on a corpus that reaches nothing must look different from coverage.

## The reference environment is not `pip install -r requirements.txt`

It is not `pip install NGSpeciesID` either, and it is not either of the README's conda recipes. See
*Goal* and *Finding 3*.

`bench/setup_reference_env.sh` is carried across from isONclust and is **simpler than isONclust's**,
because the hardest part of that script turned out to be unnecessary. What it has to do:

1. **One conda line.** `python=3.12 pip medaka spoa racon minimap2 samtools`, unpinned, from
   conda-forge and bioconda. That is the whole environment: `medaka` pulls in `parasail-python` and
   `python-edlib`, both of which bioconda builds for `osx-arm64`, so **there is no parasail source
   build and no `$M4`/`glibtoolize` workaround at all.** isONclust's script exists mostly to fight
   that build; delete that half rather than carrying it. Record what was resolved.
2. **Install the reference itself** with `pip install --no-deps -e .`, so pip does not reinstall
   parasail from PyPI over the working conda build.
3. **Build the 3.11 interpreter too.** `bench/setup_reference_env.sh 3.11 ngspeciesid-ref-311` is
   what makes *Finding 2* reproducible rather than a claim, and it is the only way to check that the
   determinism gate still fails when it should. Note bioconda has `parasail-python` and
   `python-edlib` builds for 3.10 through 3.13, so this works without a compiler as well; medaka is
   the only thing that insists on 3.12, and the determinism check does not need medaka.
4. **Refuse to record goldens for cases whose external tool is missing**, rather than recording a
   crash as a golden.

The environment this reconnaissance used, and which `bench/env/resolved-*.txt` should pin:

| | |
| --- | --- |
| clustering path | python 3.12.14, parasail 1.3.4, edlib 1.3.9.post1 |
| full pipeline | python 3.12.14, medaka 2.2.2, spoa 4.1.5, racon 1.5.0, minimap2 2.31-r1302, samtools 1.24, torch 2.9.1, numpy 2.5.3 |
| determinism counter-example | python 3.11.16, parasail 1.3.4, edlib 1.3.9.post1 |

`equivalence.sh env` must report which of the four binaries are present and refuse to record goldens
for cases that need a missing one, rather than recording a crash as a golden.

## Repo hygiene

**Not started.** Analysed, and the analysis is the same shape as isONclust's because this is the same
repository's history: NGSpeciesID was forked from isONclust and inherited its test data.

| | |
| --- | --- |
| `git clone` today | 491.82 MiB pack |
| commits | 186 |
| tags | 4: `0.0.4`, `v0.1.2.1`, `v0.3.0`, `v0.3.1` |
| distinct blobs in history | 302, 530.8 MB |
| **strippable** | **14 blobs, 520.1 MB — 98.0%** |
| retained source, scripts, docs, fixtures | 288 blobs, 10.68 MB |

The removal list, built from a tree walk (*Finding 21*):

| path | bytes |
| --- | --- |
| `test/ccs.fastq.gz.part-aa`, `-ab`, `-ac` | 94 371 840 each |
| `test/ccs.fastq.gz.part-ad` | 81 557 452 |
| `test/old_sorted_ens_100k.fastq.tar.gz` | 74 098 552 |
| `test/ENS_100k.fastq.tar.gz` | 73 094 875 |
| `test/sample_alz_2k.fastq` | 6 865 440 |
| `modules/__pycache__/p_minimizers_shared.cpython-36.pyc` | 1 125 957 |
| `test/terncytb_200.fastq` | 214 022 |
| `modules/__pycache__/cluster.cpython-36.pyc` | 16 099 |
| `modules/__pycache__/get_sorted_fastq_for_cluster.cpython-36.pyc` | 7 042 |
| `modules/.DS_Store` + `scripts/.DS_Store` | 6 148 (**one blob, two paths**) |
| `modules/__pycache__/__init__.cpython-36.pyc` | 146 |
| `test/isonclust1.out` | — |

Three things this repository needs that isONclust's tooling did not do:

1. **Six of those paths are still tracked in `HEAD`** — four `.pyc` files and two `.DS_Store` files,
   1.13 MB. `.gitignore` already lists `*.pyc` and `modules/__pycache__/*`, which does nothing for
   files that are already tracked. They need `git rm --cached` in an ordinary commit **before** any
   history rewrite, so that the rewrite and the untracking are separately reviewable.
2. **All four tags contain strippable paths** — 7 in `0.0.4`, 6 in each of the others. `filter-repo`
   rewrites tags, so this works, but two silent failure modes exist: a tag can vanish, or it can
   survive pointing at an *unrewritten* commit, which keeps every stripped byte reachable on the
   server **and** gets pulled by `git clone`, which fetches tags by default. `slim.sh` already asserts
   the tag set is preserved and each tag's tree is clean; it must be run with all four in scope.
3. **The four `ccs.fastq.gz` parts are a `split` of one file**, and the archive **already exists**.
   isONclust's `archive_data.sh` reassembles the parts and verifies with `gzip -t` before the
   originals are destroyed; NGSpeciesID inherited the same data from the same upstream history, so
   the blobs were compared rather than re-archived. Extracted from this history and checked by
   sha256 against `isONclust/repo-slim-archive/raw/test/`:

   | path | result |
   | --- | --- |
   | `test/ccs.fastq.gz.part-{aa,ab,ac,ad}`, reassembled to 364 672 972 bytes | **identical** |
   | `test/ENS_100k.fastq.tar.gz` | **identical** |
   | `test/old_sorted_ens_100k.fastq.tar.gz` | **identical** |
   | `test/sample_alz_2k.fastq` | **identical** |
   | `test/isonclust1.out` | **identical** |
   | `test/terncytb_200.fastq` | **not archived** — 214 022 bytes, NGSpeciesID-only |

   So `archive_data.sh` has exactly one 214 KB file to rescue, and step 0 costs nothing.
   `terncytb_200.fastq` is 200 real ONT reads of a **different amplicon** — Cytb/16S, median length
   428 against 632 and 816 — added in `a2128c8` and dropped in `11c516f`. Given that *The corpora*
   names a thin registry as this port's biggest measurement gap, **restoring it to `test/` after the
   rewrite is probably better than archiving it.** It does not discriminate
   `--symmetric_map_align_thresholds` (2 clusters either way), so it is an addition rather than a
   substitute.

4. **`analyze.sh` had a `set -e` bug that only a non-empty `KEEP` list could reach**, and it fired on
   the first run here. `[[ ${#KEEP[@]} -eq 0 ]] && echo ...` yields the condition's exit status when
   it is false, which under `set -e` aborted the script before the fixture verification ran — so it
   exited 1 after printing a header and no fixtures. Invisible in isONclust, where `KEEP` was empty
   and the condition was always true. Rewritten as `if/fi`, and it is the method rule landing on the
   tooling: **a harness that has never failed has not been tested, it has only been run.**

**The data will not become unreachable, and it is worth saying so up front.** isONclust measured this
after its own force-push: normal git access stops serving the old history, but the objects remain
addressable by SHA through the GitHub API until GitHub garbage-collects on its own schedule, and forks
in the same network legitimately still hold the data on their live branches. NGSpeciesID has forks.
What is actually achieved is the thing that matters: **the repository clones in single-digit megabytes
instead of 500.** Full removal was never available.

**Do this before any port work exists**, because it ends in a force-push and a port branch in flight
would make that worse.

Two things go back on top afterwards: a `.gitignore` covering the junk the rewrite removed, and
`bench/`, `tools/repo-slim/` and this file.

## Deferred improvements

Nothing here may land before the port is byte-identical, and each needs its own commit and its own
measurement. Ordered by how much they matter.

### Known bugs in the reference

| # | Fix | Effect |
| --- | --- | --- |
| *Finding 1* | `--seed`, defaulting to a fixed value | makes the recommended invocation reproducible. **This one is a blocker, not a deferral** — it has to be decided before goldens for any `--sample_size` case mean anything |
| *Finding 7* | take isONclust's guard for the empty-`error_rates` crash | turns a traceback into an explanation, same exit code. Nearly free; verify a no-op on every case that does not reach it |
| *Finding 6* | `split("\t", 1)` in `write_fastq`, and move the required-input group off the top-level parser | makes a dead feature work on ONT data |
| *Finding 4* | treat `--consensus` with no polisher as draft-only; reject a polisher with no `--consensus` | removes a crash after all the expensive work, and a silent no-op |
| *Finding 10* | `if c_id2 in already_removed: continue` in the inner loop | stops double-counting reads and polishing one cluster twice. **Measure how much it moves real output first** — it did not fire on either committed corpus |
| *Finding 5* | validate `--batch_type` as an enum in the parser | exit 1 traceback becomes exit 2 argparse error |
| *Finding 8* | validate `(--k, --w)` against the probability table at parse time | 3 361 CLI-valid settings stop crashing with a `KeyError` |
| *Finding 9* | `required=True` on `--outfolder` | exit 1 becomes exit 2 |
| *Finding 11* | move the `if args.print_output:` guard to cover the modulo | `--d 0` stops dividing by zero |
| *Finding 12* | reject a fastq whose last record has no quality string | traceback becomes a message |
| *Finding 13* | do not open `logfile.txt` for writing on the `--use_old_sorted_file` path | stops destroying the log |
| *Finding 14* | put `--ont`/`--isoseq` in a mutually-exclusive group; reject `--k`/`--w` alongside a preset | exit 0 becomes exit 2, and a silently-ignored `--w` becomes an error |
| *Finding 15* | delete both dead checks | noise removal, no behaviour |
| *Finding 18* | record `--k`/`--w`/`--q` in the logfile and refuse a `--use_old_sorted_file` mismatch | closes the only route to the 6-vs-8 tuple crash |
| *Finding 22* | reject `--max_seqs_for_consensus 0` in the parser | stops handing spoa an empty file and taking a `SIGABRT` |
| *Finding 20* | document that `--abundance_ratio` applies after subsampling, and floor the cutoff at 1 | stops singletons reaching spoa at small `--sample_size` |
| *Finding 23* | open `sorted.fastq` only on the branch that writes it | stops a failed run poisoning its output folder so the retry exits 0 on zero reads. One line, and the most user-visible of the small ones |
| *Finding 2* | `math.fsum` at the four `sum(...for...in set(...))` sites | makes Python ≤3.11 agree with ≥3.12. **Not applied**: the decision was to pin the interpreter instead. Written up so the option stays visible |

### Documentation and packaging, which are not port work but are the reason for it

- Replace both README conda recipes with the one that works (*Finding 3*).
- Say Python ≥3.12, not 3.11 (*Finding 2*).
- Unpin `parasail==1.2.4` in `setup.py` and `requirements.txt` (*Finding 3*).
- Retire `.travis.yml`. It runs Python 3.6 on Travis with `medaka=0.11.5`; none of those three things
  is obtainable. Replace with CI on Linux and macOS, x86_64 and arm64.
- The `Dockerfile` pins `python:3.6` and `medaka==1.5.0` and copies biocontainer binaries by digest.
  Once the port exists the Dockerfile is a much smaller thing — a Rust build stage plus, optionally,
  the two polishers.

### Performance and structure, once exact

Deliberately last, because performance is not this port's specification. The reference clusters
3 000 reads in 1.6 s; the interesting number is not there.

Where the time actually is, unprofiled but predictable from the other three ports: `spoa` on large
clusters, and `medaka`, which is a neural network and will dominate any run that uses it. Neither is
addressed by porting the Python. Profile before believing any of that — method point 5 — and note that
`--max_seqs_for_consensus` already exists as the knob for the first.

One structural improvement is worth naming now because the port gets it for free: the port can hold
`sorted.fastq` in memory rather than writing it and reading it back, **and must not**, because the
float formatting in the accessions is part of the contract and the file is an output.

## Method

Carried over from the isONcorrect, isONform and isONclust ports. The full versions, with the
measurements behind each point, are in those repositories' `PORTING.md`.

1. **CLI parity first**, locked by unit tests. Argument names, defaults, validation order, message
   text and exit codes. Twenty multi-word flags need explicit `long = "..."`; eight are double-dash
   single-letter; five of those carry a `dest` that differs from the flag; argparse prefix
   abbreviation is live.
2. **Differential oracles, not end-to-end tests.** Wrap the reference without modifying it; dump each
   stage's inputs *and* outputs in a stable line format; replay pure functions from Rust and diff
   stateful ones. End-to-end equivalence tells you *that* something is wrong, never *where*.
   *Finding 19* makes this mandatory here rather than advisable.
3. **Dump from the live driver too, not only from a standalone dump binary.** The one real bug in the
   isONcorrect port passed every oracle and was caught only by diffing a dump taken from the running
   driver.
4. **Build a real corpus before trusting anything.** Both corpora here are real, which is better than
   isONcorrect started with, and there are only two of them, which is worse than any of the three
   finished ports ended with. See *The corpora*.
5. **Profile before optimising, and re-profile after.** Bottlenecks move. Instrument at stage
   granularity, and *remove* sub-stage instrumentation after reading it.
6. **Measure, do not reason, about performance.** Recorded null results from the other ports: caching
   an edit-distance pattern was slower, reusing the POA engine was worth nothing, 4-bit DP cells were
   slower than 8-bit, hoisting a hash lookup out of a loop was worth nothing.
7. **Set up CI on day one, on Linux *and* macOS, x86_64 *and* arm64.** Three defects in isONcorrect
   existed only because every local check ran on one machine, and everything in this document was
   measured on one arm64 Mac. Note `.travis.yml` is dead, so CI is new work, not a migration — and
   note the entry point is `NGSpeciesID`, with capitals, which is the exact trap that cost isONform's
   CI a day on its first ext4 filesystem.
8. **Fix reference bugs upstream once measured, rather than reproducing them.** The largest accuracy
   win in the isONcorrect port was a one-character fix to the *Python*. Each such fix is its own
   commit, with goldens re-recorded. *Finding 1* is this port's first and it is a blocker; *Finding 7*
   is the second and it is already written in a sibling repository.

### Four more rules, earned in the isONform port

* **A conclusion at one depth is not a conclusion.** Sweep the depth before writing the sentence.
  Relevant immediately: everything here is measured at 280 and 3 000 reads.
* **Weight the corpus by its statistical power.** Do not pick a default from a corpus that cannot
  distinguish the options. *Finding 19* is this rule arriving before the port exists.
* **A mechanism confirmed on one case is not a prediction about the population.** *Finding 10* is a
  mechanism confirmed on a constructed case and unmeasured on real data. Do not write it up as an
  effect until it has been.
* **Measure one divergence at a time, against an exact baseline.** This port's "no deliberate
  divergences until exact" rule exists for that reason, and *Finding 1* is the awkward exception —
  the baseline cannot be exact for `--sample_size` until a seed exists.

### Two more, about the checks rather than the conclusions

* **A differential harness is only as good as the question it asks.** isONform's CI spent a run
  comparing two programs that were meant to differ and calling it a failure. Corollary earned in the
  isONclust port: **a harness that has never failed has not been tested, it has only been run.**
  Deliberately break the port and confirm the harness notices.
* **Do not edit a script while it is running.** bash reads a script incrementally, from a byte
  offset, so editing `bench/equivalence.sh` during a 20-minute recording made the running shell
  resume inside changed text: it died with `line 1140: d: unbound variable` and left a manifest that
  looked complete. The goldens were thrown away and re-recorded. Nothing in the output identified the
  cause, and the manifest's own provenance header could not have — it records the environment, not
  whether the harness changed underneath it.
* **A check that runs on one machine measures that machine.** Every number in this document came from
  one arm64 Mac. Two of the most important — that `medaka==2.0.1` and `python=3.6` are unavailable —
  are *platform-specific by construction*, and the linux-64 story is different: medaka 2.2.x is there
  too, and 2.0.1 may well be. Say "on osx-arm64" and mean it, and re-measure on Linux before writing
  the README.

## Working agreements

**Report results, not progress.** While anything is running, say nothing — no interim findings, no
partial tables, no "here is what I am about to do". One work stint gets one reply, in three parts and
no more:

1. **what was measured**, in a sentence or two;
2. **the numbers**, as a table;
3. **what to do next**, as a proposal.

Never restate a conclusion already given. Corrections are one sentence: what is now true, not a
retrospective on the error. Detail goes in this file, not in the reply — the reader will ask if they
want more.

- **Commit to a branch; never to `master`, and never force-push.** This differs from the isONcorrect,
  isONform and isONclust ports, where nothing was committed at all — here the author asked for
  branches and a PR when the port is ready (*Branches*). So: ordinary commits on `fix/installation`
  and `develop` are fine and expected; `git push`, anything touching `master` directly, and above all
  the force-push in `tools/repo-slim/slim.sh` remain human decisions. That last one is a rewrite over
  16 forks and a published paper's repository and is not something to run on anyone's behalf.
- Don't "improve" the algorithm while porting. A behaviour change and the port must not land in the
  same commit; an intentional divergence needs its own commit and a note here. **This port's
  specification is byte-identity, so there should be no intentional divergences at all** until it is
  exact — with the single, decided-in-advance exception of *Finding 1*.
- **"Behaviour" means observable output, not internal representation.** Different containers, dropping
  provably-dead entries, arena allocation, 2-bit-packed k-mers — none of that is a behaviour change if
  the emitted bytes are identical. What is not free: iteration order where it reaches results,
  tie-breaking, arithmetic and rounding. In this reference specifically: minimizer tie direction,
  summation order, the left-fold of products, float formatting, cluster-id assignment order, the
  order sequences enter the POA graph, and which of edlib's equally-optimal locations comes first.
- When you spot a possible improvement, write it into *Deferred improvements* and move on.

## First steps, in order

1. ~~Fix the installation instructions and the parasail pin~~ (*Finding 3*). **Done**, `76d3f88` on
   `fix/installation`, verified end to end from a clean env on osx-arm64. The fix turned out to be
   `--no-deps` plus dropping three version pins — conda already has ARM builds of both libraries.
2. ~~Take isONclust's *Finding 7* fix.~~ **Done**, `ea7c209`, verified a no-op on 24 of 24 cases.
3. ~~Untrack the six junk files in `HEAD`.~~ **Done**, `555f3e3`, before any history rewrite so the
   two are separately reviewable.
4. **Answer *Finding 1*.** Still the only question in this document that a human has to settle, and
   nothing that touches `--sample_size` can be recorded until it is. Recommended: add `--seed` with a
   fixed default. The README now warns about it and points at `--top_reads`, which buys time but is
   not the fix.
5. **Slim the repository.** 520.1 MB of 530.8 MB, 98.0%. Carry `tools/repo-slim/` across from
   isONclust, build the removal list from a tree walk (*Finding 21*), check the four tags, and check
   whether the archive already exists from the isONclust exercise before taking it again. Ends in a
   force-push; that is a separate human decision.
6. ~~Stand up `bench/`.~~ **Done.** `setup_reference_env.sh` (one conda line, no parasail source
   build), `equivalence.sh` with `env`/`seeds`/`cli`/`record`/`verify`/`stable`/`tools`/`stage`,
   `cases.tsv` with 51 cases including every `--consensus` combination, `corpora.tsv`,
   `dump_reference.py` with six stages, and `bench/README.md`.
7. ~~Record the CLI contract and the output goldens.~~ **Done**, on both corpora: 44 CLI cases and
   51 output cases, `stable` green over 191 checks, and the harness itself tested against five
   deliberately-broken ports (*Has the harness got teeth?*).
8. **Port the CLI**, locked by unit tests and the 44 differential cases. This is where the Rust work
   starts, and it is the next thing to do.
9. **Bring the isONclust Rust port across** and re-verify it against *this* reference, which has been
   measured to work at 12 configurations but not yet at the full case matrix. Delete the
   `replace(" ", "_")` in `readfq`. Restrict the probability table to `k >= 10` (*Finding 8*).
10. **Add the pieces isONclust does not have, in this order** — each is small and independently
    verifiable:
    a. `--m`/`--s` and `--top_reads` (reproducible today);
    b. `--sample_size`, once step 4 is settled;
    c. `--symmetric_map_align_thresholds`, with its own oracle, verified on
       `Supplementary_File1_reads.fastq` where it is visible and **not** on `sample_h1` where it is
       not.
11. **Then the consensus stage**, in dependency order: `form_draft_consensus` and its POA oracle
    (*spoa*), `highest_aln_identity`, `detect_reverse_complements`, `find_barcode_locations` and its
    tie-break measurement, `remove_barcodes`, then the medaka and racon drivers.
12. **CI on Linux and macOS, x86_64 and arm64, on day one** — method point 7, and doubly so here,
    because the central claim of this document is platform-specific and was measured on one machine.
13. **Then, and only then**, look at the deferred list. In particular *Finding 10*, which needs a
    measurement on real data before it needs a fix.

## Commits

The author asked for branches and a PR rather than a working tree left uncommitted, so unlike the
sibling ports this work is committed. See *Branches* for why there are two.

### `fix/installation`, off `master` — ready for a PR now

| sha | Message |
| --- | --- |
| `ea7c209` | `Explain, rather than traceback, when --q filters every read` |
| `76d3f88` | `Fix the install instructions, which failed on every ARM machine` |
| `555f3e3` | `Untrack committed bytecode and OS metadata` |

Every claim in `76d3f88`'s message is measured, including the per-subdir table and the end-to-end
verification from a clean environment. `ea7c209` carries its no-op evidence (24 of 24 cases). None of
the three depends on the port, and all three fix something a user hits today.

**Not in this branch, deliberately:** `--seed` (*Finding 1*), which is a user-visible behaviour change
and needs the author's decision, and the history rewrite, which is a force-push.

### `develop`, off `master` — PR when the port is exact

| sha | Contents | Message |
| --- | --- | --- |
| `981e1d3` | `PORTING.md` | `Add the Rust port plan, reconnaissance and findings` |
| `0b047c0` | `tools/repo-slim/` | `tools: add the staged history-rewrite tooling` |
| — | `bench/` | `bench: add the equivalence harness, corpora registry and goldens` |

`analyze.sh` has been run and its output — `removal-paths.txt`, `analysis.txt` — is committed;
`archive_data.sh` and `slim.sh` have not been run. `rust/` lands here when it is written.

`develop` should get `master` merged in once `fix/installation` lands, so the port is developed
against the fixed reference rather than the broken one. That matters for more than tidiness: the
`--q` guard changes what `--q 12` prints, and `cli/q_filters_all` is a recorded golden — so the
goldens in `bench/golden/` were deliberately recorded **with** the guard applied, and re-recording
after the merge should be a no-op. Check that it is.

## What is next, concretely

The reconnaissance is finished and so is the harness. Everything from here is Rust, in this order:

1. **`rust/` skeleton and the CLI.** 38 flags, 20 needing explicit `long`, 8 double-dash
   single-letter, 5 with a `dest` that differs from the flag, argparse prefix matching in all three
   of its behaviours, and the `write_fastq` subcommand behind a required top-level group. Locked by
   the 44 recorded CLI cases — `equivalence.sh cli verify` is the whole acceptance test, and 37 of
   the 44 are checkable before any clustering exists.
2. **Bring the isONclust engine across.** Measured to reproduce this reference at 12 configurations;
   re-verify against the full 51-case matrix. Delete `readfq`'s `replace(" ", "_")`. Restrict the
   probability table to `k >= 10`.
3. **`--m`/`--s`, `--top_reads`**, then `--sample_size` once *Finding 1* is settled.
4. **`--symmetric_map_align_thresholds`**, with `stage parasail` as its oracle — the dump already
   carries the second alignment ratio the flag reads — verified on `sup`, where it is visible.
5. **The consensus stage**, in dependency order, each against the oracle that already exists for it:
   `form_draft_consensus` + `spoars` (`stage spoa`), `highest_aln_identity` (`stage identity`),
   `detect_reverse_complements`, `find_barcode_locations` + the edlib HW tie-break measurement
   (`stage barcode`), `remove_barcodes`, then the medaka and racon drivers.
6. **CI on Linux and macOS, x86_64 and arm64.** Method point 7, and doubly so here: the central
   claim of this document is platform-specific and was measured on one machine.
