#!/usr/bin/env python3
"""Generate `rust/src/p_emp_probs.bin` from `modules/p_minimizers_shared.py`.

The reference keeps the empirical minimizer-sharing probabilities as a 1.79 MB
Python literal: 41 880 `(k, w, p, e1, e2)` rows, filtered at startup to the rows
where `k == args.k and abs(w - args.w) <= 2`, keyed by the rounded error-rate
pair and inserted under both orderings.

Parsing 1.79 MB of Python source at every start is not something a Rust binary
should do, so the port carries a packed blob instead. This is what makes it, and
it reads the reference's own file so the two cannot drift.

    tools/gen_p_emp_blob.py                 # write rust/src/p_emp_probs.bin
    tools/gen_p_emp_blob.py --check         # verify the committed blob, write nothing

WHY THE ENCODING IS SAFE, measured rather than assumed
------------------------------------------------------
* **At most one `w` ever matches the +-2 filter.** The table's `w` values step by
  5 within each `k`, so the filter selects exactly one row group. The
  reference's dict-overwrite behaviour therefore has nothing to overwrite and
  there is no last-row-wins subtlety to reproduce. Asserted below over every
  (k, args_w) combination.
* **Only 15 error rates are reachable.** `p_shared_minimizer_empirical` rounds
  to two decimals and clamps to [0.01, 0.15], so any row outside that grid can
  never be looked up. Asserted below.

THIS IS NOT isONclust's BLOB
----------------------------
isONclust's table has 59 628 rows and covers k 4..30, because its PR #13 added
k 4..9. NGSpeciesID never took that PR: its table is exactly isONclust's k >= 10
subset -- same values, same order, checked. Shipping isONclust's blob here would
silently make `--k 9` work, which the reference does not: it builds an empty
dict and dies with `KeyError`. That is a behaviour change wearing a data
change's clothing, so the blob is generated from THIS repository's table and
covers k 10..30 only. See PORTING.md, Finding 8.

FORMAT
------
    magic   13 bytes  "NGSPECIESPEMP"
    version  1 byte   1
    count    4 bytes  u32 little-endian, number of (k, w) pairs
    then `count` records of:
        k        1 byte
        w        1 byte
        values   120 * 8 bytes, f64 little-endian, the upper triangle of the
                 15x15 error-rate grid in the order produced by
                 `for i in range(15): for j in range(i, 15)`
"""
import argparse
import ast
import os
import struct
import sys

MAGIC = b"NGSPECIESPEMP"
VERSION = 1
N_E = 15
TRI = N_E * (N_E + 1) // 2  # 120

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SOURCE = os.path.join(ROOT, "modules", "p_minimizers_shared.py")
TARGET = os.path.join(ROOT, "rust", "src", "p_emp_probs.bin")


def tri_index(i, j):
    """Must match `p_emp.rs::tri_index` exactly."""
    if i > j:
        i, j = j, i
    return i * N_E - i * (i - 1) // 2 + (j - i)


def load_rows():
    with open(SOURCE) as fh:
        src = fh.read()
    literal = src.split("=", 1)[1].split("\ndef ")[0].strip()
    return ast.literal_eval(literal)


def e_index(e):
    """The reference's rounding: two decimals, clamped to [0.01, 0.15]."""
    r = round(e, 2)
    if r > 0.15:
        return N_E - 1
    if r < 0.01:
        return 0
    return int(round(r * 100)) - 1


def build():
    rows = load_rows()

    ks = sorted({int(r[0]) for r in rows})
    assert ks == list(range(10, 31)), (
        "expected k 10..30, got %r. If this repository's table has changed, the "
        "blob and PORTING.md's Finding 8 both need revisiting." % (ks,)
    )

    # Assertion 1: the +-2 filter selects at most one row group, for every
    # (k, args_w) the CLI can produce.
    ws_by_k = {}
    for k, w, _p, _e1, _e2 in rows:
        ws_by_k.setdefault(int(k), set()).add(int(w))
    for k, ws in ws_by_k.items():
        for args_w in range(1, 101):
            hits = [w for w in ws if abs(w - args_w) <= 2]
            assert len(hits) <= 1, (
                "k=%d args_w=%d matches %r -- the encoding assumes at most one, "
                "and the reference's dict would overwrite" % (k, args_w, hits)
            )

    # Assertion 2: every row whose error rates are inside the reachable grid
    # lands on a distinct cell, and the grid is complete.
    groups = {}
    for k, w, p, e1, e2 in rows:
        k, w = int(k), int(w)
        i, j = e_index(float(e1)), e_index(float(e2))
        # Rows outside [0.01, 0.15] cannot be looked up. They exist in the
        # table (0.002 and 0.005 in isONclust's); skip them rather than letting
        # them clamp on top of a real cell.
        if round(float(e1), 2) < 0.01 or round(float(e2), 2) < 0.01:
            continue
        cell = groups.setdefault((k, w), {})
        t = tri_index(i, j)
        if t in cell:
            assert cell[t] == float(p), (
                "k=%d w=%d cell %d has two different values: %r and %r"
                % (k, w, t, cell[t], p)
            )
        cell[t] = float(p)

    incomplete = {kw: len(c) for kw, c in groups.items() if len(c) != TRI}
    assert not incomplete, "incomplete (k, w) groups: %r" % (incomplete,)

    pairs = sorted(groups)
    out = bytearray()
    out += MAGIC
    out += bytes([VERSION])
    out += struct.pack("<I", len(pairs))
    for k, w in pairs:
        out += bytes([k, w])
        cell = groups[(k, w)]
        for t in range(TRI):
            out += struct.pack("<d", cell[t])
    return bytes(out), pairs, len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify the committed blob matches the reference table; write nothing")
    args = ap.parse_args()

    blob, pairs, n_rows = build()
    print("  source     %s" % SOURCE)
    print("  rows       %d" % n_rows)
    print("  (k, w)     %d pairs, k %d..%d" % (len(pairs), pairs[0][0], pairs[-1][0]))
    print("  blob       %d bytes" % len(blob))

    if args.check:
        if not os.path.exists(TARGET):
            print("  FAIL       %s does not exist" % TARGET)
            return 1
        with open(TARGET, "rb") as fh:
            have = fh.read()
        if have == blob:
            print("  ok         the committed blob matches the reference table")
            return 0
        print("  FAIL       the committed blob does NOT match; re-run without --check")
        print("             committed %d bytes, generated %d bytes" % (len(have), len(blob)))
        return 1

    with open(TARGET, "wb") as fh:
        fh.write(blob)
    print("  wrote      %s" % TARGET)
    return 0


if __name__ == "__main__":
    sys.exit(main())
