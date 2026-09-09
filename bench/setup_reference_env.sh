#!/usr/bin/env bash
#
# Build a pinned conda environment that can run the NGSpeciesID Python reference.
#
#   bench/setup_reference_env.sh                      # default: python 3.12, env ngspeciesid-ref
#   bench/setup_reference_env.sh 3.11 ngspeciesid-ref-311   # the determinism counter-example
#
# This is far simpler than isONclust's version of the same script, and the
# reason is worth recording because that script's complexity was almost
# entirely wasted effort by the time this one was written.
#
# isONclust's script exists to build `parasail` from source, because PyPI
# publishes no parasail wheel for macOS arm64 and parasail's setup.py fights
# back twice on Darwin (it prepends /usr/bin to PATH so macOS's m4 1.4.6 always
# wins the version probe, then downloads and fails to build m4-1.4.17; and it
# demands the Homebrew spellings `glibtoolize`/`glibtool`). Two symlinks and a
# $M4 export.
#
# None of that is needed. bioconda ships BOTH python libraries prebuilt for
# osx-arm64 -- parasail-python 1.3.4 and python-edlib 1.3.9.post1 -- and
# `medaka` depends on both, so one conda line installs everything. Measured:
# 24 s for the solve, and the resulting parasail passes the same sanity check
# the source build was there to reach.
#
# The interpreter version is not a detail. It decides whether the reference is
# deterministic: `sum()` over floats is compensated from CPython 3.12 and a
# naive left-fold before it, and four sites in this reference sum over a SET of
# quality characters. See PORTING.md, Finding 2. Build both if you want to
# reproduce that measurement -- the 3.11 environment is what proves the
# determinism gate still fails when it should.
#
# WHAT THIS ENVIRONMENT CONTAINS, AND WHY EACH PIECE
#   parasail-python  the clustering aligner, and consensus's RC detection
#   python-edlib     primer and universal-tail location (HW mode + IUPAC map)
#   spoa             the draft consensus. --consensus shells out to this binary
#   racon, minimap2  the --racon polishing path
#   medaka           the --medaka polishing path, and the reason the other two
#                    python libraries arrive for free
#   samtools         medaka's own dependency; listed so a missing one is a
#                    solve error here rather than a subprocess failure later
set -euo pipefail

PYVER="${1:-3.12}"
ENVNAME="${2:-ngspeciesid-ref}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

CONDA="${CONDA_EXE:-$(command -v conda || true)}"
if [[ -z "$CONDA" ]]; then
  for c in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda"; do
    [[ -x "$c" ]] && CONDA="$c" && break
  done
fi
[[ -n "$CONDA" ]] || { echo "error: conda not found. Install miniforge." >&2; exit 1; }

CONDA_ROOT="$(dirname "$(dirname "$CONDA")")"
ENVDIR="$CONDA_ROOT/envs/$ENVNAME"

# medaka requires python >=3.12,<3.13, so the 3.11 environment cannot have it.
# That is fine: the 3.11 env exists only to reproduce Finding 2, and the
# determinism defect is in the sorting stage, which needs neither medaka nor any
# other external binary. Asking for medaka on 3.11 would make the solve fail and
# look like a broken script.
PKGS=(pip parasail-python python-edlib spoa racon minimap2 samtools)
case "$PYVER" in
  3.12|3.13) PKGS+=(medaka) ;;
  *) echo "note: python $PYVER cannot have medaka (it needs >=3.12); building without it" ;;
esac

echo "==> creating env '$ENVNAME' (python $PYVER)"
echo "    NOTHING IS PINNED, deliberately. Every version pin in the README's"
echo "    historical recipes -- medaka==2.0.1, medaka==0.11.5, openblas==0.3.3,"
echo "    python=3.6 -- is unsatisfiable on osx-arm64. See PORTING.md, Finding 3."
"$CONDA" create -y -n "$ENVNAME" -c conda-forge -c bioconda \
  "python=$PYVER" "${PKGS[@]}" 2>&1 | tail -3

echo "==> installing the reference itself, --no-deps"
echo "    Without --no-deps, pip reads install_requires and reinstalls parasail"
echo "    FROM PyPI over the working conda build -- which on any ARM machine"
echo "    means a source build, and that build fails."
"$ENVDIR/bin/pip" install -q --no-deps -e "$ROOT"

echo "==> verifying"
"$ENVDIR/bin/python" - <<'PY'
import sys
import parasail, edlib

m = parasail.matrix_create("ACGT", 2, -2)
r = parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 5, 1, m)
assert r.score == 16 and str(r.cigar.decode, "utf-8") == "4=1X5=", "parasail behaves unexpectedly"
assert not r.saturated

# The consensus path calls parasail with a DIFFERENT opening penalty (3) from
# the clustering path's error-rate-binned 2..5, so check that call too -- it is
# a separate call site with separate defaults and the port needs both.
r3 = parasail.sg_trace_scan_16("ACGTACGTAA", "ACGTTCGTAA", 3, 1, m)
assert not r3.saturated

# HW mode with task="locations" and an IUPAC additionalEqualities map is what
# barcode_trimmer uses, and it is NOT the NW/CIGAR call isONcorrect verified.
loc = edlib.align("ACGT", "TTACGTTT", mode="HW", task="locations", k=1)
assert loc["locations"] == [(2, 5)], loc
iu = edlib.align("ACRT", "TTACGTTT", mode="HW", task="locations", k=0,
                 additionalEqualities=[("R", "A"), ("R", "G")])
assert iu["editDistance"] == 0, iu

print("    python  ", sys.version.split()[0])
print("    parasail  sg_trace_scan_16 OK at opening penalty 5 and 3")
print("    edlib     HW/locations OK, and IUPAC additionalEqualities OK")

# The interpreter's summation behaviour is part of the pinned contract.
xs = [0.1] * 10 + [1e17, -1e17]
print("    sum() is", "COMPENSATED (>=3.12)" if sum(xs) == sum(reversed(xs)) else "NAIVE (<=3.11)")
PY

echo "==> external tools"
for t in spoa racon minimap2 samtools medaka_consensus; do
  probe=""
  if [[ -x "$ENVDIR/bin/$t" ]]; then
    # `|| true` inside the substitution, not after it. `medaka_consensus
    # --version` is not a supported invocation: it prints usage and exits
    # non-zero, and with `set -o pipefail` that fails the whole pipeline, fails
    # the assignment, and -- under `set -e` -- aborted this script silently
    # right here, skipping medaka AND the resolved-versions file below. Measured,
    # not theorised: the first run of this script ended after samtools.
    # medaka_consensus is a shell wrapper with no --version; ask `medaka`.
    probe="$t"; [[ "$t" == "medaka_consensus" ]] && probe="medaka"
    v="$( { PATH="$ENVDIR/bin:$PATH" "$ENVDIR/bin/$probe" --version 2>&1 || true; } | head -1 | tr -d '\r')"
    printf '    %-18s %s\n' "$t" "${v:-present}"
  else
    printf '    %-18s ABSENT -- cases needing it cannot be recorded\n' "$t"
  fi
done

echo "==> writing $HERE/env/resolved-$(uname -s)-$(uname -m)-py$PYVER.txt"
mkdir -p "$HERE/env"
{
  echo "# $ENVNAME  on $(uname -s) $(uname -r) $(uname -m)"
  "$ENVDIR/bin/python" -V
  "$CONDA" list -n "$ENVNAME" 2>/dev/null \
    | grep -E '^(python|parasail-python|python-edlib|spoa|racon|minimap2|samtools|medaka|numpy|pytorch|torch) '
} > "$HERE/env/resolved-$(uname -s)-$(uname -m)-py$PYVER.txt"
cat "$HERE/env/resolved-$(uname -s)-$(uname -m)-py$PYVER.txt"

echo
echo "Use it with:"
echo "  REF_PYTHON=$ENVDIR/bin/python bench/equivalence.sh env"
echo
echo "Note the external tools must be on PATH for --consensus cases:"
echo "  export PATH=$ENVDIR/bin:\$PATH"
