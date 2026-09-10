#!/usr/bin/env bash
#
# The equivalence harness. Byte-identity with the Python reference is the
# acceptance criterion for the port, and this is what checks it.
#
#   bench/equivalence.sh env      # is the reference environment usable? which tools are present?
#   bench/equivalence.sh seeds    # determinism gate -- run this BEFORE recording
#   bench/equivalence.sh cli      # CLI contract: exit codes, stdout, stderr
#   bench/equivalence.sh record   # record goldens from the reference
#   bench/equivalence.sh verify   # run the port, diff against the goldens
#   bench/equivalence.sh tools    # a missing external binary must be NAMED, not tracebacked
#   bench/equivalence.sh stable   # recording twice must give identical goldens
#   bench/equivalence.sh stage sort        # files a ported stage owns, across the matrix
#   bench/equivalence.sh stage minimizers  # a stage with no file output, via dumps
#   bench/equivalence.sh stage mapping     # get_best_cluster, replayed from the live driver
#   bench/equivalence.sh stage parasail    # the clustering aligner, on this tool's parameters
#   bench/equivalence.sh stage spoa        # the sequences handed to the POA, in order
#   bench/equivalence.sh stage identity    # the consensus path's parasail call, both orientations
#   bench/equivalence.sh stage barcode     # every edlib HW call and its FULL locations list
#   bench/equivalence.sh all      # everything
#
# Environment:
#   REF_PYTHON   interpreter that can import parasail and edlib
#                (default: the ngspeciesid-ref conda env; see setup_reference_env.sh)
#   PORT_BIN     the Rust binary under test (default: rust/target/release/NGSpeciesID)
#   CORPUS       input fastq path, or a name from bench/corpora.tsv (default: sup)
#   GOLDEN       where goldens live (default: bench/golden)
#
# WHAT COUNTS AS A DIFFERENCE
# ---------------------------
# Every file the tool writes, byte for byte. The list is longer and more
# conditional than isONclust's, because --consensus writes per-cluster files and
# shells out to three external programs:
#
#   always
#     final_clusters.tsv          the actual result
#     final_cluster_origins.tsv   representative per cluster, and its error_rate
#     sorted.fastq                the sorted input, with the score in each header
#     logfile.txt                 error-rate summary statistics
#   --t > 1
#     <n>/pre_clusters.csv        one dir per merge iteration
#     <n>/cluster_origins.csv
#   --consensus
#     consensus_reference_<id>.fasta      the spoa draft
#     reads_to_consensus_<id>.fastq       the reads handed to the polisher
#   --consensus --racon
#     racon_cl_id_<id>/consensus.fasta
#     racon_cl_id_<id>/racon_polished_it_<i>.fasta
#     racon_cl_id_<id>/read_alignments_it_<i>.paf
#   --consensus --medaka
#     medaka_cl_id_<id>/consensus.fasta   (or .fastq with --medaka_fastq)
#     plus the .fai/.mmi/.bam/.hdf files medaka leaves beside the draft
#
# NEVER compared -- measured to differ between runs of the REFERENCE itself,
# because they capture timings:
#   *_stderr_it_*.txt, mm2_stderr_it_*.txt, stdout.txt, stderr.txt
# Three runs of `--consensus --racon` differed in exactly those files and in
# nothing else; every consensus fasta was identical. So the polishers ARE
# reproducible and their output is contract -- their logs are not.
#
# sorted.fastq is NOT an intermediate to be skipped. The score is formatted into
# every read accession with Python's float repr and then parsed back out with
# float(), so it is both an output and an input, and Rust's default float
# formatting does not match Python's (Python writes 1234.0 and 1e-05 where Rust
# writes 1234 and 0.00001).
#
# logfile.txt is included because it is the only place the error-rate
# distribution is observable, and error_rate is where the reference's
# interpreter-dependent determinism defect shows up. See `seeds`.
#
# THE ACCESSIONS CONTAIN SPACES. Unlike isONclust, this reference's readfq does
# NOT substitute them, so column 2 of final_clusters.tsv is a whole ONT header.
# Anything in this harness that splits a line must split on TAB, not on
# whitespace -- that mistake is a live bug in the reference itself
# (PORTING.md, Finding 6).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

REF_PYTHON="${REF_PYTHON:-$HOME/miniforge3/envs/ngspeciesid-ref/bin/python}"
PORT_BIN="${PORT_BIN:-$ROOT/rust/target/release/NGSpeciesID}"
# CORPUS accepts a path, or a name from bench/corpora.tsv.
NGSPECIESID_DATA="${NGSPECIESID_DATA:-$HOME/data/amplicon}"
resolve_corpus() {
  local c="$1" p
  [[ -f "$c" ]] && { echo "$c"; return; }
  p="$(awk -F'\t' -v n="$c" '$1==n {print $2; exit}' "$ROOT/bench/corpora.tsv" 2>/dev/null)"
  [[ -z "$p" ]] && { echo "$c"; return; }          # not a known name: pass through and let it fail visibly
  case "$p" in
    /*)     echo "$p" ;;
    test/*) echo "$ROOT/$p" ;;
    *)      echo "$NGSPECIESID_DATA/$p" ;;
  esac
}
CORPUS="$(resolve_corpus "${CORPUS:-sup}")"
GOLDEN="${GOLDEN:-$ROOT/bench/golden}"
WORK="${WORK:-$(mktemp -d)}"

# NOT the equivalence comparison -- `record` hashes every file a run writes,
# found with `find`, so the numbered intermediate dirs are pinned too. This is
# only the floor for `seeds`: the list of files that must EXIST, so five runs
# that all crashed cannot agree with each other vacuously.
OUT_FILES=(final_clusters.tsv final_cluster_origins.tsv sorted.fastq logfile.txt)

PASS=0; FAIL=0
ok()   { printf '    ok    %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '    FAIL  %s\n' "$1"; FAIL=$((FAIL+1)); }
info() { printf '    info  %s\n' "$1"; }

# ---------------------------------------------------------------------------

cmd_env() {
  echo "==> reference environment"
  [[ -x "$REF_PYTHON" ]] || { bad "REF_PYTHON not executable: $REF_PYTHON"; return; }
  if "$REF_PYTHON" - <<'PY'
import sys
import parasail, edlib
m = parasail.matrix_create("ACGT", 2, -2)
r = parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 5, 1, m)
assert r.score == 16 and str(r.cigar.decode, "utf-8") == "4=1X5="
# The consensus path uses opening penalty 3, not the clustering path's binned
# 2..5. Separate call site, separate defaults; check both.
assert not parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 3, 1, m).saturated
# HW + task="locations" + additionalEqualities is what barcode_trimmer needs,
# and it is NOT the NW/CIGAR call isONcorrect's port reimplemented.
assert edlib.align("ACGT", "TTACGTTT", mode="HW", task="locations", k=1)["locations"] == [(2, 5)]
assert edlib.align("ACRT", "TTACGTTT", mode="HW", task="locations", k=0,
                   additionalEqualities=[("R", "A"), ("R", "G")])["editDistance"] == 0
xs = [0.1] * 10 + [1e17, -1e17]
compensated = sum(xs) == sum(reversed(xs))
print(f"    ok    python {sys.version.split()[0]}, parasail + edlib usable at both call sites")
print(f"    ok    sum() is {'COMPENSATED (>=3.12)' if compensated else 'NAIVE (<=3.11)'}")
if not compensated:
    print("    info  on this interpreter the reference is NOT deterministic run to run;")
    print("    info  `equivalence.sh seeds` will fail. See PORTING.md, Finding 2.")
PY
  then PASS=$((PASS+2)); else bad "reference env unusable"; fi
  [[ -f "$CORPUS" ]] && ok "corpus present: $CORPUS ($(( $(wc -l < "$CORPUS") / 4 )) reads)" \
                     || bad "corpus missing: $CORPUS"

  # The external binaries. isONclust's harness deliberately checked for NONE of
  # these, because --consensus was out of scope there. Here --consensus IS the
  # tool, so which of these is present decides which cases can be recorded at
  # all -- and recording a FileNotFoundError traceback as a golden would be
  # worse than skipping the case.
  local t v probe missing=()
  for t in spoa racon minimap2 medaka_consensus; do
    if command -v "$t" >/dev/null 2>&1; then
      # `|| true` INSIDE the substitution. `medaka_consensus --version` is not a
      # supported invocation -- it prints usage and exits non-zero -- and with
      # `set -o pipefail` that fails the assignment and, under `set -e`, aborts
      # the whole harness here. This is the third time this exact shape has bitten
      # in this repository; see also tools/repo-slim/analyze.sh.
      # medaka_consensus is a shell wrapper with no --version; ask `medaka`.
      probe="$t"; [[ "$t" == "medaka_consensus" ]] && probe="medaka"
      v="$( { "$probe" --version 2>&1 || true; } | head -1 | tr -d '\r')"
      ok "$t: ${v:-present}"
    else
      missing+=("$t")
    fi
  done
  if [[ ${#missing[@]} -eq 0 ]]; then
    ok "all four external tools present -- every case is recordable"
  else
    info "absent: ${missing[*]}"
    info "  cases needing an absent tool are SKIPPED by 'record', not recorded as crashes."
    info "  bench/setup_reference_env.sh installs all four; remember to put its bin on PATH."
  fi

  # The goldens are only valid for the tool versions that produced them, so the
  # manifest records them. Say so here too, because it is the thing people forget.
  info "external tool versions are recorded in the golden manifest -- upgrading any"
  info "  of spoa/racon/minimap2/medaka invalidates every --consensus golden."
}

# ---------------------------------------------------------------------------
# The determinism gate. Recording a golden from a non-deterministic reference
# records one of its several possible answers, so this runs first and refuses to
# be skipped.

cmd_seeds() {
  echo "==> determinism: does the reference agree with itself?"
  local seeds=(0 1 2 7 12345) s d first tag f
  # Entries chosen to cover the three code paths that could differ: the single
  # sweep, the multiprocessing sweep (which also writes the numbered merge
  # directories), and the consensus stage (which shells out to spoa and racon,
  # and whose captured stderr logs carry timings and are excluded below).
  for entry_args in "--ont --t 1" "--ont --t 8" "--ont --t 1 --consensus --racon"; do
    # A consensus entry with no spoa on PATH would leave five identical empty
    # directories and report five green "stable" lines -- a vacuous pass on no
    # data, which is exactly what the OUT_FILES floor below exists to prevent.
    # Skip it loudly instead.
    if [[ "$entry_args" == *--consensus* ]] && ! command -v spoa >/dev/null 2>&1; then
      info "skipping [$entry_args] -- spoa not on PATH"
      continue
    fi
    tag="$(echo "$entry_args" | tr -dc 'a-z0-9')"
    for s in "${seeds[@]}"; do
      d="$WORK/seeds_${s}_$tag"
      rm -rf "$d"; mkdir -p "$d"
      PYTHONHASHSEED="$s" "$REF_PYTHON" NGSpeciesID $entry_args \
        --fastq "$CORPUS" --outfolder "$d" >"$d.stdout" 2>"$d.stderr" || true
    done

    # The file list is DISCOVERED from the runs rather than taken from
    # OUT_FILES. isONclust's version held its own list and so never checked the
    # --t 8 merge intermediates for seed independence even though one of its
    # entries was --t 8. Here the gap would be much wider: the consensus entry
    # writes per-cluster drafts, read files and a racon directory, none of which
    # a fixed list would name.
    local listing="" got setdiffers=0
    for s in "${seeds[@]}"; do
      d="$WORK/seeds_${s}_$tag"
      got="$(cd "$d" && find . -type f ! -name '*stderr*' ! -name 'stdout.txt' ! -name 'stderr.txt' \
             | sed 's|^\./||' | LC_ALL=C sort | tr '\n' ' ')"
      if [[ -z "$listing" && "$s" == "${seeds[0]}" ]]; then
        listing="$got"
      elif [[ "$got" != "$listing" ]]; then
        # A seed-dependent iteration COUNT, or a seed-dependent set of clusters
        # passing --abundance_ratio, would otherwise read as agreement on
        # whichever files happen to exist under both seeds.
        bad "SEED-DEPENDENT FILE SET: [$entry_args]"
        info "  seed ${seeds[0]}: $listing"
        info "  seed $s: $got"
        setdiffers=1
      fi
    done
    [[ $setdiffers -eq 0 ]] && \
      ok "same file set across ${#seeds[@]} seeds: $(wc -w <<<"$listing" | tr -d ' ') files  [$entry_args]"

    # OUT_FILES is the floor, and it is why this stays here. Five runs that all
    # crashed leave five empty directories, and every file then agrees with
    # every other vacuously -- a green determinism gate on no data at all.
    for f in "${OUT_FILES[@]}"; do
      [[ " $listing " == *" $f "* ]] || bad "reference wrote no $f  [$entry_args]"
    done
    # And for a consensus entry, at least one draft must exist, for the same reason.
    if [[ "$entry_args" == *--consensus* ]]; then
      grep -q 'consensus_reference_' <<<"$listing" \
        || bad "reference formed no consensus at all  [$entry_args]"
    fi

    for f in $listing; do
      first=""; local differs=0
      for s in "${seeds[@]}"; do
        d="$WORK/seeds_${s}_$tag"
        [[ -f "$d/$f" ]] || continue
        local h; h="$(shasum -a 256 "$d/$f" | cut -d' ' -f1)"
        [[ -z "$first" ]] && first="$h" && continue
        [[ "$h" == "$first" ]] || differs=1
      done
      if [[ $differs -eq 0 ]]; then
        ok "stable across ${#seeds[@]} seeds: $f  [$entry_args]"
      else
        bad "SEED-DEPENDENT: $f  [$entry_args]"
        # Say where, not just that. This is the actionable half.
        info "  seed 0 vs seed 1:"
        "$REF_PYTHON" bench/diffsummary.py \
          "$WORK/seeds_0_$tag/$f" "$WORK/seeds_1_$tag/$f" "$(basename "$f")" 2>/dev/null || true
      fi
    done
  done

  cmd_seeds_sample_size
}

# --sample_size is seeded from --seed (default 0) since 0.3.2, so it belongs in
# the case matrix like any other flag. This function is what is left of the gate
# that used to assert the OPPOSITE, and it is kept rather than deleted because
# the property it checks is easy to lose again: seeding is three lines in one
# branch of one `elif`, and nothing else in the tool would notice if they went.
#
# It checks three things:
#
#   1. the same command twice gives the same answer;
#   2. a DIFFERENT --seed gives a different answer, so the flag is actually
#      reaching the sampler rather than being parsed and ignored;
#   3. --top_reads is reproducible and ignores --seed.
#
# (2) is the one worth having. A port that accepts --seed and then samples from
# an unseeded generator passes (1) and (3) and fails only (2).
#
# PYTHONHASHSEED is FIXED here on purpose. Varying it would confound this with
# Finding 2's interpreter-dependent defect, which is a different problem.
cmd_seeds_sample_size() {
  echo "==> --sample_size is seeded (Finding 1, fixed in 0.3.2)"
  # The sample size has to be SMALLER than the number of reads that survive
  # filtering, or the subsample never happens: the guard is
  # `0 < args.sample_size < len(read_array)`, so --sample_size 500 on the
  # 280-read smoke corpus silently takes every read. This check reported a false
  # alarm on its first run for exactly that reason. Derive the size from the
  # corpus.
  local nreads; nreads=$(( $(wc -l < "$CORPUS") / 4 ))
  local ss=$(( nreads / 2 ))
  if [[ "$ss" -lt 2 ]]; then
    info "corpus has $nreads reads -- too few to subsample; skipping"
    return 0
  fi
  info "corpus has $nreads reads; drawing --sample_size $ss"

  local i d hashes=() h
  for i in 1 2 3 4 5; do
    d="$WORK/ss_$i"
    rm -rf "$d"; mkdir -p "$d"
    PYTHONHASHSEED=0 "$REF_PYTHON" NGSpeciesID --ont --t 1 --sample_size "$ss" \
      --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
    if [[ ! -f "$d/final_clusters.tsv" ]]; then
      bad "--sample_size $ss produced no output at all -- cannot assess"
      return 0
    fi
    hashes+=("$(shasum -a 256 "$d/final_clusters.tsv" | cut -d' ' -f1)")
  done
  local distinct; distinct="$(printf '%s\n' "${hashes[@]}" | LC_ALL=C sort -u | wc -l | tr -d ' ')"
  if [[ "$distinct" == "1" ]]; then
    ok "reproducible: 5 identical runs gave 1 result"
  else
    bad "--sample_size is NOT reproducible: 5 identical runs gave $distinct results"
    info "  This was Finding 1 and it was fixed in 0.3.2 by seeding random.sample"
    info "  from --seed. If it is back, check that the elif branch in main() still"
    info "  builds a random.Random(args.seed) rather than calling random.sample."
    return 0
  fi

  # (2) A different seed must give a different subsample. Without this, a port
  # that parses --seed and ignores it passes every other check here.
  local other="$WORK/ss_seed7"
  rm -rf "$other"; mkdir -p "$other"
  PYTHONHASHSEED=0 "$REF_PYTHON" NGSpeciesID --ont --t 1 --sample_size "$ss" --seed 7 \
    --fastq "$CORPUS" --outfolder "$other" >/dev/null 2>&1 || true
  h="$(shasum -a 256 "$other/final_clusters.tsv" 2>/dev/null | cut -d' ' -f1)"
  if [[ -z "$h" ]]; then
    bad "--seed 7 produced no output"
  elif [[ "$h" == "${hashes[0]}" ]]; then
    bad "--seed 7 gives the SAME result as --seed 0 -- the flag is being ignored"
  else
    ok "--seed 7 gives a different subsample from the default --seed 0"
  fi

  # (3) --top_reads is the deterministic-by-construction path and must ignore
  # --seed. It is what the README recommends for "the best reads" rather than
  # "a reproducible random subset", so a regression here is a documentation bug
  # as well as a code one.
  local first="" tdiff=0
  for i in 1 2 3; do
    d="$WORK/tr_$i"
    rm -rf "$d"; mkdir -p "$d"
    PYTHONHASHSEED=0 "$REF_PYTHON" NGSpeciesID --ont --t 1 --sample_size "$ss" --top_reads \
      --seed "$i" --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
    h="$(shasum -a 256 "$d/final_clusters.tsv" 2>/dev/null | cut -d' ' -f1)"
    [[ -z "$h" ]] && { bad "--top_reads produced no output"; return 0; }
    [[ -z "$first" ]] && first="$h" && continue
    [[ "$h" == "$first" ]] || tdiff=1
  done
  [[ $tdiff -eq 0 ]] && ok "--top_reads is reproducible and ignores --seed (3 seeds, 1 result)" \
                     || bad "--top_reads varies with --seed -- it should not use the sampler at all"
  return 0
}

# ---------------------------------------------------------------------------
# The CLI contract. Argument names, defaults, validation order, the exact stderr
# and stdout text, and the exit code. Recorded from the reference, then replayed
# against the port.
#
# Several of these exit 0 on what is logically an error -- no arguments, and
# `--ont --isoseq` together, both print a message and `sys.exit()` with no
# argument. That is the contract, warts included; do not "fix" it in the port
# without giving the divergence its own commit and a note in PORTING.md.


# Two different reasons a CLI case cannot pass, kept apart on purpose. Calling
# them both "pending" would hide the fact that one of them will never resolve.
#
# PENDING: the invocation is valid, so the reference goes on to run the tool.
# These pass once the corresponding stages exist.
CLI_NEEDS_STAGES=" d_zero ont_over_k k_then_ont isoseq_over_w medaka_no_consensus abbrev_outf use_old_k_mismatch "

# DIVERGENT BY DESIGN. Nothing is dropped in this port -- see PORTING.md, Scope
# -- so unlike isONclust's harness this list is EMPTY, and the `dropped`
# subcommand it existed for is replaced by `tools`, which checks that a missing
# external binary is named rather than tracebacked.
CLI_DIVERGENT=" "

cli_case() { # cli_case <name> <args...>
  local name="$1"; shift
  local d="$GOLDEN/cli/$name"
  mkdir -p "$d"
  if [[ "$MODE" == "record" ]]; then
    set +e
    "$REF_PYTHON" NGSpeciesID "$@" >"$d/stdout" 2>"$d/stderr"; echo $? >"$d/exit"
    set -e
    # Paths, timings and tracebacks are not contract; scrub them so the golden
    # is portable. Traceback frames name the absolute path of the interpreter's
    # site-packages and of this checkout, so an unscrubbed golden is valid on
    # exactly one machine -- and a dozen of these cases exit through a traceback.
    scrub "$d/stdout" "$d/stderr"
    ok "recorded cli/$name (exit $(cat "$d/exit"))"
  else
    [[ -f "$d/exit" ]] || { bad "no golden for cli/$name -- re-run: equivalence.sh cli record"; return; }
    if [[ "$CLI_NEEDS_STAGES" == *" $name "* ]]; then
      info "pending cli/$name -- valid invocation, needs the clustering stages"
      return 0
    fi
    if [[ "$CLI_DIVERGENT" == *" $name "* ]]; then
      info "divergent by design cli/$name"
      return 0
    fi
    set +e
    "$PORT_BIN" "$@" >"$WORK/o" 2>"$WORK/e"; local rc=$?
    set -e
    scrub "$WORK/o" "$WORK/e"
    local want_rc; want_rc="$(cat "$d/exit")"
    local bad_parts=()
    [[ "$rc" == "$want_rc" ]] || bad_parts+=("exit $rc want $want_rc")
    diff -q "$WORK/e" "$d/stderr" >/dev/null || bad_parts+=("stderr")
    diff -q "$WORK/o" "$d/stdout" >/dev/null || bad_parts+=("stdout")
    if [[ ${#bad_parts[@]} -eq 0 ]]; then
      ok "cli/$name"
    else
      bad "cli/$name: ${bad_parts[*]}"
      diff "$d/stderr" "$WORK/e" 2>/dev/null | head -5 | sed 's/^/          stderr| /' || true
      diff "$d/stdout" "$WORK/o" 2>/dev/null | head -5 | sed 's/^/          stdout| /' || true
    fi
    # cli_case is called at top level under `set -e`. Without this, the exit
    # status of the last command above becomes the function's, and a failing
    # case aborts the whole run -- which is how isONclust's 28 cases came to
    # report 1.
    return 0
  fi
}

# Scrubbing is its own function because record and verify MUST apply exactly the
# same transformation; two copies of a long sed drift, and a drift here shows up
# as a permanent unexplainable diff.
scrub() {
  sed -i.bak -E \
    -e 's#(/private)?(/var/folders/[^ ]*|/tmp/[^ ]*)#<TMPDIR>#g' \
    -e 's#/[^ "]*/(NGSpeciesID|sample_h1\.fastq|Supplementary_File1_reads\.fastq|Supplementary_File3_primer\.txt)#<PATH>/\1#g' \
    -e 's#File "[^"]*/([^/"]+)", line [0-9]+#File "<PATH>/\1", line <LINE>#g' \
    -e 's/[0-9]+\.[0-9]{4,}/<TIME>/g' \
    "$@"
  # `rm -f "$@".bak` is WRONG and was here first: with two arguments it expands
  # to `rm -f <first> <second>.bak`, which DELETES the first file. Every CLI
  # case then failed with `diff: .../o: No such file or directory` -- 37 of 37,
  # and only because a "port" that was the reference itself was run against the
  # goldens. Suffix each element individually.
  local f
  for f in "$@"; do rm -f "$f.bak"; done
}

cmd_cli() {
  MODE="${1:-record}"
  echo "==> CLI contract ($MODE)"
  # Without this, verify mode returns from every case silently and the run
  # reports "0 passed, 0 failed", which reads like success.
  if [[ "$MODE" == "verify" && ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to verify yet"
    return
  fi

  # --- the easy half: things that exit before doing any work ---
  cli_case version      --version
  cli_case help         --help
  cli_case h_short      -h
  cli_case noargs
  cli_case wf_help      write_fastq --help

  # --- argparse's own error paths. All exit 2 with a fixed usage block plus one
  # --- distinguishing final line. A hand-written parser gets these wrong by
  # --- default, so they are contract, not decoration. Note `noargs` and
  # --- `unknown_flag` BOTH report the missing required --fastq group rather
  # --- than the thing that was actually wrong; that is the contract.
  cli_case unknown_flag --fastq "$CORPUS" --no-such-flag
  cli_case bad_int      --k abc --fastq "$CORPUS"
  cli_case bad_float    --q xyz --fastq "$CORPUS"
  cli_case bad_int_t    --t 1.5 --fastq "$CORPUS"
  cli_case missing_val  --fastq "$CORPUS" --k
  cli_case bad_subcmd   bogus_subcmd
  cli_case wf_stray_flag write_fastq --k 5
  # argparse accepts any unambiguous prefix; clap does not. Both halves of that
  # are contract, and they are three different behaviours, not two:
  #
  #   1. a prefix that matches one flag is ACCEPTED    (--outfold -> --outfolder)
  #   2. an EXACT match wins over longer flags sharing it, so --m is --m
  #      (target_length) and not ambiguous with --min_shared/--mapped_threshold/
  #      --medaka..., and --s is --s and not --sample_size
  #   3. a prefix matching several flags is REJECTED with exit 2 and a message
  #      listing the candidates
  #
  # Case 2 is the surprising one and is why `exact_m`/`exact_f` exit 1 (they
  # parse fine and then die on the missing --outfolder, Finding 9) rather than 2.
  cli_case abbrev_outf  --fastq "$CORPUS" --outfold "$WORK/ab" --t 1
  cli_case exact_m      --m 5 --fastq "$CORPUS"
  cli_case exact_f      --f "$CORPUS"
  cli_case ambiguous_me --me --fastq "$CORPUS"
  cli_case ambiguous_min --min 5 --fastq "$CORPUS"
  cli_case ambiguous_r  --r --fastq "$CORPUS"
  cli_case ambiguous_prim --prim x --fastq "$CORPUS"
  # --version and -h fire during parsing and win over arguments that would fail.
  cli_case version_wins --version --fastq /nope --k abc
  cli_case help_wins    --fastq /nope -h

  # --- mutually exclusive groups. Two of the three are real argparse groups and
  # --- exit 2; the --ont/--isoseq pair is hand-rolled and exits ZERO, which a
  # --- pipeline cannot detect. See PORTING.md, Finding 14.
  cli_case both_presets --ont --isoseq --fastq "$CORPUS"
  cli_case medaka_racon --fastq "$CORPUS" --outfolder "$WORK/mr" --t 1 --medaka --racon
  cli_case tails_primer --fastq "$CORPUS" --outfolder "$WORK/tp" --t 1 \
                        --remove_universal_tails --primer_file test/Supplementary_File3_primer.txt

  # --- validation that happens after parsing ---
  cli_case w_lt_k       --fastq "$CORPUS" --k 20 --w 15
  cli_case w_gt_100     --fastq "$CORPUS" --k 15 --w 101
  # --outfolder is optional in argparse and mandatory in fact: Finding 9.
  cli_case no_outfolder --ont --fastq "$CORPUS" --t 1

  # --- crashes that are contract until a commit says otherwise. Every one of
  # --- these was reproduced from a command; see PORTING.md's exit-code table.
  cli_case d_zero       --ont --fastq "$CORPUS" --outfolder "$WORK/dz" --t 1 --d 0
  cli_case k_too_small  --fastq "$CORPUS" --outfolder "$WORK/k9" --t 1 --k 9 --w 20
  cli_case k_too_big    --fastq "$CORPUS" --outfolder "$WORK/k31" --t 1 --k 31 --w 50
  cli_case kw_gap       --fastq "$CORPUS" --outfolder "$WORK/kw" --t 1 --k 11 --w 100
  cli_case q_filters_all --ont --fastq "$CORPUS" --outfolder "$WORK/qa" --t 1 --q 12
  cli_case batch_weighted --ont --fastq "$CORPUS" --outfolder "$WORK/bw" --t 4 --batch_type weighted
  cli_case batch_bogus  --ont --fastq "$CORPUS" --outfolder "$WORK/bb" --t 4 --batch_type nosuchtype
  # A FRESH directory, and removed first. Not fussiness: a failed run of this
  # invocation leaves an EMPTY sorted.fastq behind -- get_sorted_fastq_for_cluster
  # opens it for writing before the branch that would have filled it -- and a
  # second run in the same folder then takes the "use the existing sorted file"
  # path and exits 0 having clustered zero reads. So this case's exit code
  # depends on whether it has run before, and `equivalence.sh stable` caught it
  # by recording twice in one process: exit 1 the first time, 0 the second.
  # See PORTING.md, Finding 23.
  rm -rf "$WORK/uo_missing"
  cli_case use_old_missing --use_old_sorted_file --outfolder "$WORK/uo_missing" --t 1 --ont
  cli_case consensus_no_polisher --ont --fastq "$CORPUS" --outfolder "$WORK/cnp" --t 1 --consensus
  cli_case max_seqs_zero --ont --fastq "$CORPUS" --outfolder "$WORK/msz" --t 1 \
                         --consensus --racon --max_seqs_for_consensus 0

  # --- a fastq whose last record has no trailing newline. Finding 12. Built
  # --- here rather than committed, because a file whose whole point is a
  # --- missing final byte does not survive an editor or a git checkout.
  local nonl="$WORK/no_trailing_newline.fastq"
  head -8 "$CORPUS" > "$nonl.tmp"
  printf '%s' "$(cat "$nonl.tmp")" > "$nonl"
  cli_case no_trailing_nl --ont --fastq "$nonl" --outfolder "$WORK/nn" --t 1

  # --- silently doing nothing is also contract: a polisher without --consensus
  # --- exits 0 and produces no consensus at all. Finding 4's mirror.
  cli_case medaka_no_consensus --ont --fastq "$CORPUS" --outfolder "$WORK/mnc" --t 1 --medaka

  # --- the presets OVERWRITE an explicit --k/--w regardless of order, so
  # --- --ont --k 99 is k=13 w=20 and runs fine rather than failing w<k.
  cli_case ont_over_k   --ont --k 99 --fastq "$CORPUS" --outfolder "$WORK/ok1" --t 1
  cli_case k_then_ont   --k 99 --ont --fastq "$CORPUS" --outfolder "$WORK/ok2" --t 1
  cli_case isoseq_over_w --isoseq --w 7 --fastq "$CORPUS" --outfolder "$WORK/ok3" --t 1

  # --- write_fastq. Unreachable without a redundant top-level --fastq, and then
  # --- broken on any header containing a space -- which is every ONT header.
  # --- Finding 6. Both halves are pinned: the unreachability above
  # --- (wf_stray_flag) and the KeyError here.
  # Needs a real clustering. Produce one here rather than depending on `record`
  # having run first in the same $WORK -- a CLI case that silently skips
  # depending on invocation order is a case that can never fail.
  local clusters="$WORK/wf_input_clusters.tsv"
  if [[ ! -f "$clusters" ]]; then
    local cd_="$WORK/wf_cluster_src"
    rm -rf "$cd_"; mkdir -p "$cd_"
    PYTHONHASHSEED=0 "$REF_PYTHON" NGSpeciesID --ont --t 1 \
      --fastq "$CORPUS" --outfolder "$cd_" >/dev/null 2>&1 || true
    cp "$cd_/final_clusters.tsv" "$clusters" 2>/dev/null || true
  fi
  if [[ -f "$clusters" ]]; then
    cli_case wf_spaces --fastq "$CORPUS" write_fastq \
             --clusters "$clusters" --fastq "$CORPUS" --outfolder "$WORK/wfs"
  else
    bad "cli/wf_spaces: could not produce a clustering to feed write_fastq"
  fi

  # --- --use_old_sorted_file with a sorted.fastq made at a different --k, which
  # --- reaches the 6-vs-8 tuple unpack. Finding 18. Needs two invocations, so
  # --- it is built here rather than being a single-command case.
  local uo="$WORK/uo_kmismatch"
  rm -rf "$uo"; mkdir -p "$uo"
  "$REF_PYTHON" NGSpeciesID --fastq "$CORPUS" --outfolder "$uo" --t 1 --k 13 --w 20 \
    >/dev/null 2>&1 || true
  if [[ -s "$uo/sorted.fastq" ]]; then
    cli_case use_old_k_mismatch --use_old_sorted_file --outfolder "$uo" --t 1 --k 25 --w 50
  else
    info "skipping cli/use_old_k_mismatch -- could not produce a sorted.fastq"
  fi
}

# ---------------------------------------------------------------------------

# Which external binary a case needs, derived from its args rather than listed
# separately -- a hand-maintained second list would drift the first time someone
# adds a case. Returns the name of the missing tool, or nothing.
case_missing_tool() { # case_missing_tool <args>
  local args="$1" t
  [[ "$args" != *--consensus* ]] && return 0
  # spoa is needed by every --consensus case; the polisher depends on the flag.
  for t in spoa; do
    command -v "$t" >/dev/null 2>&1 || { echo "$t"; return 0; }
  done
  if [[ "$args" == *--racon* ]]; then
    for t in racon minimap2; do
      command -v "$t" >/dev/null 2>&1 || { echo "$t"; return 0; }
    done
  fi
  if [[ "$args" == *--medaka* ]]; then
    command -v medaka_consensus >/dev/null 2>&1 || { echo "medaka_consensus"; return 0; }
  fi
  return 0
}

run_case() { # run_case <name> <entry> <args> <runner> <outdir>  -> echoes exit code or SKIP
  local name="$1" entry="$2" args="$3" runner="$4" outdir="$5"
  rm -rf "$outdir"; mkdir -p "$outdir"

  # A case whose external tool is absent must be SKIPPED, not run. Running it
  # records a FileNotFoundError traceback as the golden, which then "passes"
  # forever on any machine that also lacks the tool -- a case that can never
  # fail. isONclust's harness had no such cases; here more than a fifth of the
  # matrix shells out to something.
  local miss; miss="$(case_missing_tool "$args")"
  [[ -n "$miss" ]] && { echo "SKIP:$miss"; return; }

  set +e
  if [[ "$entry" == "write_fastq" ]]; then
    # write_fastq consumes a clustering, so it needs one to exist first.
    # The FULL clustering, not the head -- write_fastq consumes it. Kept in
    # $WORK rather than $GOLDEN so it is not committed: it is an input to three
    # cases, not a golden, and on the 3000-read corpus it is 686 KB.
    #
    # PRODUCED HERE ON DEMAND, from the reference, rather than relying on the
    # `default` case having run first. Depending on that made `verify` fail all
    # three write_fastq cases: `record` populates $WORK, `verify` gets a fresh
    # $WORK, and the cases then reported a missing input as a mismatched golden.
    # A case whose result depends on what ran before it in the same $WORK is a
    # case that says nothing.
    local clusters="$WORK/wf_input_clusters.tsv"
    if [[ ! -f "$clusters" ]]; then
      local wfsrc="$WORK/wf_cluster_src"
      rm -rf "$wfsrc"; mkdir -p "$wfsrc"
      PYTHONHASHSEED=0 "$REF_PYTHON" NGSpeciesID --ont --t 1 \
        --fastq "$CORPUS" --outfolder "$wfsrc" >/dev/null 2>&1 || true
      cp "$wfsrc/final_clusters.tsv" "$clusters" 2>/dev/null || true
    fi
    [[ -f "$clusters" ]] || { echo "SKIP:clustering"; return; }
    # NOTE the redundant top-level --fastq. It is not a mistake and it is not
    # optional: the --fastq/--use_old_sorted_file mutually exclusive group is
    # `required=True` on the TOP-LEVEL parser, so `NGSpeciesID write_fastq
    # --fastq X ...` exits 2 asking for a top-level --fastq that the subcommand's
    # own --fastq does not satisfy. PORTING.md, Finding 6. The port must
    # reproduce that, so the harness has to invoke it the way a user is forced to.
    $runner --fastq "$CORPUS" write_fastq --clusters "$clusters" --fastq "$CORPUS" \
            --outfolder "$outdir" $args >"$outdir.stdout" 2>"$outdir.stderr"
  else
    $runner $args --fastq "$CORPUS" --outfolder "$outdir" \
            >"$outdir.stdout" 2>"$outdir.stderr"
  fi
  local rc=$?
  set -e
  echo "$rc"
}

# Files whose CONTENT is not contract because it carries timings. Measured, not
# assumed: three runs of `--consensus --racon` against the same input differed
# in exactly these and in nothing else -- every consensus fasta was identical.
# So the polishers are reproducible and their output IS contract; their captured
# logs are not.
#
# Their EXISTENCE is still contract, and `verify` still checks the file count,
# so a port that writes none of them fails.
#
# The medaka BAM and its index are here for a different reason, and it is not
# timings: minimap2 and samtools write their own command lines into the BAM's
# @PG header, and those command lines contain the ABSOLUTE PATH of the output
# folder --
#
#   @PG ID:minimap2 ... CL:minimap2 ... /tmp/xyz/consensus_reference_17.fasta ...
#
# -- so the file can never match across two runs in different directories, and
# every case in this harness runs in a fresh temp dir. Found by `verify` against
# a "port" that WAS the reference: 48 of 51 cases passed and the three medaka
# cases failed on exactly these two files.
#
# consensus_probs.hdf is deliberately NOT excluded. It was measured stable
# across the same comparison, so it stays in the contract; excluding a file
# because it is the kind of file that might vary is how a contract quietly
# stops covering anything.
NOT_CONTRACT_RE='(_stderr_it_[0-9]+\.txt|mm2_stderr_it_[0-9]+\.txt|/stdout\.txt|/stderr\.txt|calls_to_draft\.bam(\.bai)?)$'
is_contract() { ! [[ "$1" =~ $NOT_CONTRACT_RE ]]; }


# A case line with no tab silently degrades into "run with no arguments", which
# then records a plausible-looking golden for the wrong invocation. Refuse.
check_cases() {
  local bad_lines
  bad_lines="$(grep -vE '^#|^$' bench/cases.tsv | grep -vcE $'^[^\t]+\t[^\t]+\t' || true)"
  if [[ "$bad_lines" != "0" ]]; then
    bad "bench/cases.tsv has $bad_lines line(s) that are not TAB-separated into 3 fields"
    info "an editor probably expanded tabs to spaces; see the comment in that file"
    exit 1
  fi
  info "$(grep -vcE '^#|^$' bench/cases.tsv) cases, all tab-separated"
}

cmd_record() {
  echo "==> recording goldens from the reference"
  check_cases
  mkdir -p "$GOLDEN/out"

  # Goldens are a MANIFEST OF HASHES, not the files themselves. Recorded
  # verbatim, the 27 cases come to 318 MB: final_cluster_origins.tsv carries
  # every representative's full sequence and quality string (3.5 MB a case),
  # sorted.fastq is the whole input again (6.5 MB a case), and `wf_N0` alone is
  # 1240 files. None of that belongs in a repository this exercise just shrank
  # from 492 MB to 1 MB.
  #
  # Hashes are enough to FAIL correctly. They are not enough to say what broke,
  # so `verify` re-runs the reference for a failing case and diffs against that.
  # The reference is 2 seconds a case; the storage is not worth it.
  {
    echo "# NGSpeciesID reference goldens -- per-file sha256"
    echo "# No timestamp on purpose: it made this file differ on every re-record,"
    echo "# which leaves git permanently dirty and hides real changes in the noise."
    echo "# git records when it was committed; what matters for validity is below."
    echo "# corpus:   $(basename "$CORPUS")  sha256 $(shasum -a 256 "$CORPUS" | cut -d' ' -f1)"
    "$REF_PYTHON" -c "import sys,parasail,edlib; print(f'# reference: python {sys.version.split()[0]}, parasail + edlib')"
    # The external tool versions are part of golden validity, not decoration: a
    # spoa upgrade changes every consensus_reference_*.fasta, and a racon or
    # medaka upgrade changes every polished consensus. Record them, or the
    # --consensus goldens are unfalsifiable.
    for t in spoa racon minimap2 medaka; do
      if command -v "$t" >/dev/null 2>&1; then
        echo "# $t: $( { "$t" --version 2>&1 || true; } | head -1 | tr -d '\r')"
      else
        echo "# $t: ABSENT -- cases needing it were skipped"
      fi
    done
    "$REF_PYTHON" -c "xs=[0.1]*10+[1e17,-1e17]; print('# sum():    ' + ('compensated (>=3.12)' if sum(xs)==sum(reversed(xs)) else 'NAIVE (<=3.11) -- these goldens are NOT reproducible'))"
    echo "# PYTHONHASHSEED=0 for every case"
    echo "#"
    echo "# case	exit	relpath	sha256	bytes"
  } > "$GOLDEN/manifest.tsv"

  local n=0
  while IFS=$'\t' read -r name entry args; do
    [[ "$name" =~ ^# ]] && continue
    [[ -z "${name// }" ]] && continue
    local d="$WORK/rec/$name"
    local rc; rc="$(PYTHONHASHSEED=0 run_case "$name" "$entry" "$args" "$REF_PYTHON NGSpeciesID" "$d")"
    if [[ "$rc" == SKIP:* ]]; then
      # A skip is RECORDED, as a meta row with no hash rows, so `verify` can
      # tell "this case was never recorded" from "this case legitimately has no
      # golden here". Without the row, verify reports "no golden" and moves on,
      # which is neither a pass nor a fail -- the port could do anything.
      printf '%s\t%s\t%s\t%s\t%s\n' "$name" "-" "__skip__" "${rc#SKIP:}" "0" >> "$GOLDEN/manifest.tsv"
      info "skipped $name (needs ${rc#SKIP:})"
      continue
    fi
    local nf=0
    while IFS= read -r rel; do
      # Hash only the files whose content is contract. The polishers' captured
      # stderr carries timings, so hashing it makes every case permanently fail
      # -- which trains you to ignore the harness. Their existence is still
      # counted, below, via files=.
      is_contract "$rel" || continue
      printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$rc" "$rel" \
        "$(shasum -a 256 "$d/$rel" | cut -d' ' -f1)" "$(wc -c < "$d/$rel" | tr -d ' ')" \
        >> "$GOLDEN/manifest.tsv"
      nf=$((nf+1))
    done < <(cd "$d" && find . -type f | sed 's|^\./||' | LC_ALL=C sort)
    # One meta row per case, ALWAYS. A case that legitimately writes no files
    # (wf_N10 on a corpus with no cluster of 10+ reads) otherwise contributes no
    # rows at all, and `verify` then finds no expected exit code, reports "no
    # golden" and silently skips it -- neither pass nor fail. The port could do
    # anything there and nothing would say so.
    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$rc" "__meta__" "files=$nf" "0" >> "$GOLDEN/manifest.tsv"
    # Keep ONE case's small files verbatim, so there is something to read by eye
    # without running anything. sorted.fastq and final_cluster_origins.tsv are
    # excluded by size; their hashes are in the manifest like everything else.
    if [[ "$name" == "default" ]]; then
      mkdir -p "$GOLDEN/sample"
      # HEADS, not whole files. On the 3000-read corpus final_clusters.tsv is
      # 686 KB -- one line per read, each carrying a whole ONT accession -- and
      # this directory exists only so a human can read something without
      # running anything. The hashes in manifest.tsv are the contract; these are
      # a courtesy, and a courtesy has no business being the largest thing in
      # the repository after the fixtures.
      head -40 "$d/final_clusters.tsv" > "$GOLDEN/sample/final_clusters.tsv.head" 2>/dev/null || true
      cp "$d/logfile.txt" "$GOLDEN/sample/" 2>/dev/null || true
      head -8 "$d/sorted.fastq" > "$GOLDEN/sample/sorted.fastq.head" 2>/dev/null || true
      cut -f1,5,6 "$d/final_cluster_origins.tsv" \
        > "$GOLDEN/sample/final_cluster_origins.id_score_errorrate.tsv" 2>/dev/null || true
    fi
    n=$((n+1))
    ok "recorded $name (exit $rc, $nf files)"
  done < bench/cases.tsv
  info "$n cases -> $GOLDEN/manifest.tsv ($(wc -c < "$GOLDEN/manifest.tsv" | tr -d ' ') bytes)"
}

cmd_verify() {
  echo "==> verifying the port against the goldens"
  check_cases
  [[ -f "$GOLDEN/manifest.tsv" ]] || { bad "no manifest -- run: bench/equivalence.sh record"; return; }
  if [[ ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to verify yet"
    info "expected until rust/ exists; build with"
    info "  cargo build --release --manifest-path rust/Cargo.toml"
    return
  fi
  while IFS=$'\t' read -r name entry args; do
    [[ "$name" =~ ^# ]] && continue
    [[ -z "${name// }" ]] && continue
    local d="$WORK/port/$name"
    local skipped; skipped="$(awk -F'\t' -v n="$name" '$1==n && $3=="__skip__" {print $4; exit}' "$GOLDEN/manifest.tsv")"
    if [[ -n "$skipped" ]]; then
      info "skipped $name -- golden was not recorded ($skipped absent when recording)"
      continue
    fi
    local want_exit; want_exit="$(awk -F'\t' -v n="$name" '$1==n && $3=="__meta__" {print $2; exit}' "$GOLDEN/manifest.tsv")"
    [[ -n "$want_exit" ]] || { bad "no golden for $name -- re-run: equivalence.sh record"; continue; }
    # If the port needs a tool this machine lacks, say so rather than failing it.
    local miss; miss="$(case_missing_tool "$args")"
    if [[ -n "$miss" ]]; then
      info "cannot verify $name here -- $miss not on PATH (golden exists)"
      continue
    fi
    local want_files; want_files="$(awk -F'\t' -v n="$name" '$1==n && $3=="__meta__" {sub(/^files=/,"",$4); print $4; exit}' "$GOLDEN/manifest.tsv")"
    local rc; rc="$(run_case "$name" "$entry" "$args" "$PORT_BIN" "$d")"

    local mismatched=() missing=()
    while IFS=$'\t' read -r rel want_sha want_bytes; do
      if [[ ! -f "$d/$rel" ]]; then missing+=("$rel"); continue; fi
      local got; got="$(shasum -a 256 "$d/$rel" | cut -d' ' -f1)"
      [[ "$got" == "$want_sha" ]] || mismatched+=("$rel")
    done < <(awk -F'\t' -v n="$name" '$1==n && $3!="__meta__" {print $3"\t"$4"\t"$5}' "$GOLDEN/manifest.tsv")

    # `grep -v` exits 1 when its input is EMPTY, and under `set -o pipefail`
    # that fails the substitution and, under `set -e`, aborts the whole run --
    # after the last case that wrote files and before the summary line, so the
    # output looked like a crash with no verdict. It happens whenever a case
    # legitimately writes nothing, which the write_fastq cases do. `|| true`
    # goes INSIDE the substitution, around the grep.
    local got_files; got_files="$(cd "$d" 2>/dev/null && find . -type f | sed 's|^\./||' \
      | { grep -vE "$NOT_CONTRACT_RE" || true; } | wc -l | tr -d ' ')"
    [[ "${got_files:-0}" == "$want_files" ]] || mismatched+=("file count: got ${got_files:-0}, want $want_files")

    # A file the port writes that the reference does not is also a failure.
    local extra=()
    while IFS= read -r rel; do
      is_contract "$rel" || continue
      awk -F'\t' -v n="$name" -v r="$rel" '$1==n && $3!="__meta__" && $3==r {found=1} END {exit !found}' \
        "$GOLDEN/manifest.tsv" || extra+=("$rel")
    done < <(cd "$d" 2>/dev/null && find . -type f | sed 's|^\./||' | LC_ALL=C sort)

    if [[ ${#mismatched[@]} -eq 0 && ${#missing[@]} -eq 0 && ${#extra[@]} -eq 0 && "$rc" == "$want_exit" ]]; then
      ok "$name"
      continue
    fi
    bad "$name (exit $rc, want $want_exit)"
    [[ ${#missing[@]}  -gt 0 ]] && info "  not written by the port: ${missing[*]}"
    [[ ${#extra[@]}    -gt 0 ]] && info "  written by the port only: ${extra[*]}"
    # Hashes cannot say what moved, so re-derive the reference output for this
    # one case and diff properly. Two seconds beats 318 MB in git.
    if [[ ${#mismatched[@]} -gt 0 ]]; then
      info "  differing: ${mismatched[*]}"
      local r="$WORK/refre/$name"
      PYTHONHASHSEED=0 run_case "$name" "$entry" "$args" "$REF_PYTHON NGSpeciesID" "$r" >/dev/null
      for rel in "${mismatched[@]}"; do
        [[ -f "$r/$rel" && -f "$d/$rel" ]] || continue
        info "  --- $rel ---"
        "$REF_PYTHON" bench/diffsummary.py "$r/$rel" "$d/$rel" "$(basename "$rel")" || true
      done
    fi
  done < bench/cases.tsv
}


# ---------------------------------------------------------------------------
# Missing external binaries must be NAMED.
#
# isONclust's harness had a `dropped` subcommand here, asserting that the port
# refuses four out-of-scope flags. Nothing is dropped in this port -- consensus
# IS NGSpeciesID -- so that check has no analogue. This is its replacement, and
# it checks the one place where the port is allowed to improve on the reference
# for free.
#
# Today, a missing tool surfaces as a Python traceback:
#
#   FileNotFoundError: [Errno 2] No such file or directory: 'racon'
#
# which is a non-zero exit with the right information buried in twelve lines of
# stack. The port must exit non-zero and NAME THE TOOL, and nothing else. Two
# things must hold: non-zero exit, and the tool named in the output.
#
# This CANNOT be a recorded golden, because the reference does not do it. It is
# the port's own contract, so it is asserted directly -- and it is a deliberate
# divergence with its own commit and its own note in PORTING.md.
#
# The check works by running the port with a PATH that has been emptied of the
# tool in question, so it does not require the tool to be genuinely absent from
# the machine.
cmd_tools() {
  echo "==> a missing external tool must be refused, by name"
  if [[ ! -x "$PORT_BIN" ]]; then
    bad "no port binary at $PORT_BIN -- nothing to check yet"
    return 0
  fi
  # A shim directory that shadows the real tool with nothing. Putting an empty
  # dir FIRST on PATH does not shadow anything, so each iteration builds a PATH
  # containing only the other tools.
  local shim="$WORK/toolshim"
  local tool flags rc out
  for spec in "spoa:--consensus --racon" \
              "racon:--consensus --racon" \
              "minimap2:--consensus --racon" \
              "medaka_consensus:--consensus --medaka"; do
    tool="${spec%%:*}"; flags="${spec#*:}"
    rm -rf "$shim"; mkdir -p "$shim"
    # Symlink every tool EXCEPT this one into the shim, then use only the shim
    # plus the system directories the port itself needs.
    local t src
    for t in spoa racon minimap2 medaka_consensus medaka samtools; do
      [[ "$t" == "$tool" ]] && continue
      src="$(command -v "$t" 2>/dev/null || true)"
      [[ -n "$src" ]] && ln -sf "$src" "$shim/$t"
    done
    set +e
    out="$(PATH="$shim:/usr/bin:/bin" "$PORT_BIN" $flags \
           --fastq "$CORPUS" --outfolder "$WORK/tools_$tool" --t 1 --ont 2>&1)"
    rc=$?
    set -e
    if [[ $rc -eq 0 ]]; then
      bad "$tool missing but the port exited 0 -- must fail loudly"
    elif ! grep -qF -- "$tool" <<<"$out"; then
      bad "$tool missing, port exited $rc, but the message does not name it"
      sed 's/^/          /' <<<"$out" | head -3
    elif [[ "$(wc -l <<<"$out" | tr -d ' ')" -gt 6 ]]; then
      # A traceback is not a message. Six lines is generous.
      bad "$tool named, but the port printed $(wc -l <<<"$out" | tr -d ' ') lines -- looks like a traceback"
      sed 's/^/          /' <<<"$out" | head -3
    else
      ok "$tool missing: exit $rc, named, $(wc -l <<<"$out" | tr -d ' ') line(s)"
    fi
  done
  return 0
}

# ---------------------------------------------------------------------------
# Recording twice must give the same goldens. A golden containing a timestamp, a
# temp path, a PID or a duration can never be matched by anything -- including
# the reference itself -- so it is not a check, it is a permanent failure that
# trains you to ignore the harness. Two got through: a `tempfile.mkdtemp()` path
# in the --consensus stdout, and the manifest's own "recorded:" line.

cmd_stable() {
  echo "==> recording twice must be byte-identical"
  local a="$WORK/stable_a" b="$WORK/stable_b"
  GOLDEN="$a" cmd_record >/dev/null 2>&1
  GOLDEN="$a" cmd_cli record >/dev/null 2>&1
  GOLDEN="$b" cmd_record >/dev/null 2>&1
  GOLDEN="$b" cmd_cli record >/dev/null 2>&1
  local unstable
  unstable="$(diff -rq "$a" "$b" 2>&1 || true)"
  if [[ -z "$unstable" ]]; then
    ok "goldens are reproducible across two recordings"
  else
    bad "goldens are NOT reproducible -- these contain run-varying data:"
    sed 's/^/          /' <<<"$unstable" | head -8
    while read -r _ f1 _ _; do
      [[ -f "$f1" ]] || continue
      diff "$f1" "${f1/$a/$b}" 2>/dev/null | grep -E '^[<>]' | head -2 | sed 's/^/            /'
    done <<<"$unstable"
  fi
}

# ---------------------------------------------------------------------------
# Per-stage verification. The port cannot produce final_clusters.tsv yet, so
# `verify` fails every output case -- which says nothing about the stages that
# ARE done. This checks the files the ported stages actually own, across the
# same case matrix.
#
# NGSPECIESID_STAGE is an environment variable, not a flag, on purpose: the CLI is
# a byte-for-byte contract and must not grow options the reference lacks.


# `get_kmer_minimizers` never writes a file, so it is compared through a dump on
# both sides: bench/dump_reference.py against the port's NGSPECIESID_STAGE=
# minimizers, which emit the same format. PORTING.md Finding 5 is why this stage
# in particular cannot be trusted to end-to-end goldens.
#
# Both sides consume the REFERENCE's sorted.fastq, so a difference here is the
# minimizer selection and not the sort. (`stage sort` already proves the sort.)
cmd_stage_minimizers() {
  local pairs=("13 20" "15 50" "9 25" "20 100" "15 15" "4 10")
  local d="$WORK/mz"
  rm -rf "$d"; mkdir -p "$d"
  PYTHONHASHSEED=0 $REF_PYTHON NGSpeciesID --k 15 --w 50 --t 1 \
    --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
  if [[ ! -s "$d/sorted.fastq" ]]; then
    bad "could not produce a sorted.fastq from $CORPUS"
    return
  fi
  local kw k w rl pl
  for kw in "${pairs[@]}"; do
    k="${kw% *}"; w="${kw#* }"
    $REF_PYTHON bench/dump_reference.py --stage minimizers \
      --sorted-fastq "$d/sorted.fastq" --k "$k" --w "$w" --out "$d/ref.tsv" 2>/dev/null
    NGSPECIESID_STAGE=minimizers "$PORT_BIN" --k "$k" --w "$w" \
      --fastq "$d/sorted.fastq" --outfolder "$d/o" > "$d/port.tsv" 2>/dev/null
    rl=$(wc -l < "$d/ref.tsv" | tr -d ' '); pl=$(wc -l < "$d/port.tsv" | tr -d ' ')
    if cmp -s "$d/ref.tsv" "$d/port.tsv"; then
      # Report what the case actually exercised, so a vacuous pass is visible.
      local empties subk
      empties=$(awk -F'\t' '$2>=0 && $3=="" ' "$d/ref.tsv" | wc -l | tr -d ' ')
      subk=$(awk -F'\t' -v k="$k" '$2>=0 && $3!="" && length($3)<k' "$d/ref.tsv" | wc -l | tr -d ' ')
      ok "k=$k w=$w: $rl minimizers identical (empty: $empties, sub-k: $subk)"
    else
      bad "k=$k w=$w: ref $rl lines, port $pl lines"
      diff "$d/ref.tsv" "$d/port.tsv" | head -6 | sed 's/^/          /' || true
    fi
  done
  # See cli_case: without this the last command's status becomes the function's
  # and `set -e` aborts before the summary line. Third time this shape has bitten.
  return 0
}


# `get_best_cluster` is stateful -- the minimizer database grows as reads become
# representatives -- so its calls are captured by wrapping the LIVE driver and
# replayed here. That is method point 3: an oracle fed only recorded inputs
# cannot catch a port that follows a different trajectory, but this records what
# the reference actually did on the way through.
#
# NOTE the corpus matters more here than anywhere else. Both simulated corpora
# decide ZERO reads by mapping -- every assignment goes through the alignment
# fallback -- so running this on `smoke` proves nothing. The counts are printed
# so a vacuous pass is visible.
cmd_stage_mapping() {
  local settings=("13 20" "15 50")
  local d="$WORK/mp"
  rm -rf "$d"; mkdir -p "$d"
  PYTHONHASHSEED=0 $REF_PYTHON NGSpeciesID --k 15 --w 50 --t 1 \
    --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
  if [[ ! -s "$d/sorted.fastq" ]]; then
    bad "could not produce a sorted.fastq from $CORPUS"
    return
  fi
  local kw k w calls assigned
  for kw in "${settings[@]}"; do
    k="${kw% *}"; w="${kw#* }"
    $REF_PYTHON bench/dump_reference.py --stage mapping \
      --sorted-fastq "$d/sorted.fastq" --k "$k" --w "$w" --out "$d/dump.tsv" 2>/dev/null
    grep '^RES' "$d/dump.tsv" > "$d/ref.tsv" || true
    NGSPECIESID_MAPPING_DUMP="$d/dump.tsv" NGSPECIESID_STAGE=mapping "$PORT_BIN" \
      --k "$k" --w "$w" --fastq "$d/sorted.fastq" --outfolder "$d/o" \
      > "$d/port.tsv" 2>/dev/null
    calls=$(wc -l < "$d/ref.tsv" | tr -d ' ')
    assigned=$(awk -F'\t' '$2>=0' "$d/ref.tsv" | wc -l | tr -d ' ')
    if cmp -s "$d/ref.tsv" "$d/port.tsv"; then
      ok "k=$k w=$w: $calls decisions identical ($assigned assigned a cluster)"
      [[ "$assigned" == "0" ]] && info "  WARNING: no read mapped -- this corpus does not exercise the stage"
    else
      bad "k=$k w=$w: $(diff "$d/ref.tsv" "$d/port.tsv" | grep -c '^<') of $calls decisions differ"
      diff "$d/ref.tsv" "$d/port.tsv" | head -6 | sed 's/^/          /' || true
    fi
  done
  # See cli_case: without this the last command's status becomes the function's
  # and `set -e` aborts before the summary line. Third time this shape has bitten.
  return 0
}


# The aligner, on isONclust's OWN parameters. parasail.rs is carried across from
# isONform, which verified it against isONcorrect's scoring (match 4, mismatch
# -8, open 12). isONclust uses match 2, mismatch -2 and an opening penalty of
# 2..5 chosen per comparison, which can reach different tie-breaking paths -- so
# reusing a verified module is not the same as having verified it here.
#
# Both the CIGAR and the derived alignment ratio are compared: the ratio is what
# the clustering decision actually reads, and two different optimal paths can
# score the same while giving different ratios.
cmd_stage_parasail() {
  local d="$WORK/pa"
  rm -rf "$d"; mkdir -p "$d"
  PYTHONHASHSEED=0 $REF_PYTHON NGSpeciesID --k 13 --w 20 --t 1 \
    --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
  if [[ ! -s "$d/sorted.fastq" ]]; then
    bad "could not produce a sorted.fastq from $CORPUS"
    return 0
  fi
  $REF_PYTHON bench/dump_reference.py --stage parasail \
    --sorted-fastq "$d/sorted.fastq" --k 13 --w 20 --out "$d/ref.tsv" 2>/dev/null
  local n; n=$(wc -l < "$d/ref.tsv" | tr -d ' ')
  if [[ "$n" == "0" ]]; then
    ok "no alignments recorded"
    info "  WARNING: this corpus never reaches the alignment path"
    return 0
  fi
  NGSPECIESID_PARASAIL_DUMP="$d/ref.tsv" NGSPECIESID_STAGE=parasail "$PORT_BIN" \
    --k 13 --w 20 --fastq "$d/sorted.fastq" --outfolder "$d/o" > "$d/port.tsv" 2>/dev/null
  local opens; opens=$(awk -F'\t' '{print $2}' "$d/ref.tsv" | sort -u | paste -sd, -)
  if cmp -s "$d/ref.tsv" "$d/port.tsv"; then
    ok "$n alignments identical, cigar and ratio (opening penalties: $opens)"
  else
    bad "$(diff "$d/ref.tsv" "$d/port.tsv" | grep -c '^<') of $n alignments differ"
    diff "$d/ref.tsv" "$d/port.tsv" | head -4 | cut -c1-140 | sed 's/^/          /' || true
  fi
  return 0
}

# The three consensus-stage oracles share one shape: dump from the reference,
# replay from the port, diff. They are written as one function rather than three
# because the differences between them are arguments, and three near-copies
# drift -- which is how a stage ends up being checked against the harness's own
# idea of the reference rather than against the reference.
#
# Each prints WHAT IT ACTUALLY EXERCISED, because a vacuous pass must be
# visible. Measured on the two committed corpora:
#
#   stage        sample_h1                Supplementary_File1_reads.fastq
#   spoa         255 lines                2 766 lines
#   identity     1 call  (RC won 1)       6 calls (RC won 4)
#   barcode      16 calls, ZERO hits      32 calls, 6 hits
#
# `barcode` on the smoke corpus finds no primer at all, so it exercises the
# edlib call and never the cut-position logic. Do not read a green `barcode` on
# `smoke` as coverage of remove_barcodes.
cmd_stage_dump() { # cmd_stage_dump <stage> <extra dump args...>
  local stage="$1"; shift
  local d="$WORK/dump_$stage"
  rm -rf "$d"; mkdir -p "$d"
  PYTHONHASHSEED=0 $REF_PYTHON NGSpeciesID --ont --t 1 \
    --fastq "$CORPUS" --outfolder "$d" >/dev/null 2>&1 || true
  if [[ ! -s "$d/sorted.fastq" ]]; then
    bad "could not produce a sorted.fastq from $CORPUS"
    return 0
  fi
  $REF_PYTHON bench/dump_reference.py --stage "$stage" \
    --sorted-fastq "$d/sorted.fastq" --k 13 --w 20 "$@" --out "$d/ref.tsv" 2>/dev/null || true
  local n; n=$(wc -l < "$d/ref.tsv" 2>/dev/null | tr -d ' ')
  if [[ "${n:-0}" == "0" ]]; then
    bad "$stage: the reference recorded nothing -- this corpus does not reach the stage"
    return 0
  fi
  NGSPECIESID_STAGE="$stage" NGSPECIESID_DUMP="$d/ref.tsv" "$PORT_BIN" --ont --t 1 \
    --fastq "$CORPUS" --outfolder "$d/o" > "$d/port.tsv" 2>/dev/null || true
  if cmp -s "$d/ref.tsv" "$d/port.tsv"; then
    ok "$stage: $n lines identical  $(stage_coverage "$stage" "$d/ref.tsv")"
  else
    bad "$stage: $(diff "$d/ref.tsv" "$d/port.tsv" 2>/dev/null | grep -c '^<') of $n lines differ"
    diff "$d/ref.tsv" "$d/port.tsv" 2>/dev/null | head -4 | cut -c1-140 | sed 's/^/          /' || true
  fi
  return 0
}

# What a stage actually exercised, so a pass on a corpus that reaches nothing
# says so rather than looking like coverage.
stage_coverage() {
  case "$1" in
    spoa)     echo "($(grep -c '^SPOA' "$2" | tr -d ' ') POA calls)" ;;
    identity) echo "($(wc -l < "$2" | tr -d ' ') pairs, RC orientation won $(awk -F'\t' '$5>$4' "$2" | wc -l | tr -d ' '))" ;;
    barcode)  echo "($(wc -l < "$2" | tr -d ' ') edlib calls, $(awk -F'\t' '$4!=-1' "$2" | wc -l | tr -d ' ') found a primer)" ;;
    *)        echo "" ;;
  esac
}

cmd_stage() {
  local which="${1:-sort}"
  if [[ ! -x "$PORT_BIN" ]]; then
    echo "==> stage '$which'"
    bad "no port binary at $PORT_BIN"
    return 0
  fi
  # The dump-based stages sweep their own (k, w) settings rather than the case
  # matrix, so they get their own header and skip check_cases.
  if [[ "$which" == "minimizers" ]]; then
    echo "==> stage 'minimizers': (minimizer, position) lists, via dumps"
    cmd_stage_minimizers
    return 0
  fi
  if [[ "$which" == "mapping" ]]; then
    echo "==> stage 'mapping': get_best_cluster decisions, replayed from the live driver"
    cmd_stage_mapping
    return 0
  fi
  if [[ "$which" == "parasail" ]]; then
    echo "==> stage 'parasail': alignments replayed from the live driver"
    cmd_stage_parasail
    return 0
  fi
  if [[ "$which" == "spoa" ]]; then
    echo "==> stage 'spoa': the sequences handed to the POA, in order, and the consensus"
    command -v spoa >/dev/null 2>&1 || { bad "spoa not on PATH"; return 0; }
    cmd_stage_dump spoa --abundance_ratio 0.02
    return 0
  fi
  if [[ "$which" == "identity" ]]; then
    echo "==> stage 'identity': forward and reverse-complement identity per center pair"
    command -v spoa >/dev/null 2>&1 || { bad "spoa not on PATH"; return 0; }
    cmd_stage_dump identity --abundance_ratio 0.02
    return 0
  fi
  if [[ "$which" == "barcode" ]]; then
    echo "==> stage 'barcode': every edlib HW call and its FULL locations list"
    command -v spoa >/dev/null 2>&1 || { bad "spoa not on PATH"; return 0; }
    cmd_stage_dump barcode --abundance_ratio 0.02 \
      --primer_file test/Supplementary_File3_primer.txt
    return 0
  fi
  echo "==> stage '$which': the files this stage owns, across the case matrix"
  check_cases
  local files
  case "$which" in
    sort) files="sorted.fastq logfile.txt" ;;
    *) bad "unknown stage '$which'"; return ;;
  esac

  while IFS=$'\t' read -r name entry args; do
    [[ "$name" =~ ^# ]] && continue
    [[ -z "${name// }" ]] && continue
    # write_fastq does not run the sorting stage at all
    [[ "$entry" == "write_fastq" ]] && continue
    local r="$WORK/stage_ref/$name" p="$WORK/stage_port/$name"
    rm -rf "$r" "$p"; mkdir -p "$r" "$p"
    PYTHONHASHSEED=0 $REF_PYTHON NGSpeciesID $args --fastq "$CORPUS" --outfolder "$r" \
      >/dev/null 2>&1 || true
    NGSPECIESID_STAGE="$which" "$PORT_BIN" $args --fastq "$CORPUS" --outfolder "$p" \
      >/dev/null 2>&1 || true
    local bad_files=()
    for f in $files; do
      if [[ ! -f "$r/$f" && ! -f "$p/$f" ]]; then continue; fi
      cmp -s "$r/$f" "$p/$f" 2>/dev/null || bad_files+=("$f")
    done
    if [[ ${#bad_files[@]} -eq 0 ]]; then
      ok "$name"
    else
      bad "$name: ${bad_files[*]}"
      for f in "${bad_files[@]}"; do
        [[ -f "$r/$f" && -f "$p/$f" ]] || { info "  one side missing $f"; continue; }
        info "  --- $f ---"
        diff "$r/$f" "$p/$f" | head -4 | sed 's/^/          /' || true
      done
    fi
  done < bench/cases.tsv
}
# ---------------------------------------------------------------------------

case "${1:-all}" in
  env)     cmd_env ;;
  seeds)   cmd_seeds ;;
  cli)     cmd_cli "${2:-record}" ;;
  record)  cmd_record ;;
  verify)  cmd_verify ;;
  tools)   cmd_tools ;;
  stable)  cmd_stable ;;
  stage)   cmd_stage "${2:-sort}" ;;
  all)     cmd_env; cmd_seeds; cmd_cli record; cmd_record; cmd_stable
           for st in sort minimizers mapping parasail spoa identity barcode; do cmd_stage "$st"; done
           cmd_verify; cmd_tools ;;
  *) echo "usage: $0 {env|seeds|cli [record|verify]|record|verify|tools|stable|stage [sort|minimizers|mapping|parasail|spoa|identity|barcode]|all}" >&2; exit 2 ;;
esac

echo
echo "==> $PASS passed, $FAIL failed"
[[ $FAIL -eq 0 ]] || exit 1
