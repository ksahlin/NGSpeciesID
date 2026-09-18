#!/usr/bin/env python3
"""Dump a stage's inputs and outputs from the reference, without modifying it.

End-to-end equivalence says *that* the port is wrong, never *where*. Stages
whose output reaches a file can be diffed directly (`equivalence.sh stage sort`);
stages whose output stays in memory cannot, and this is how those get checked.

PORTING.md's Finding 19 is the argument that this file is mandatory rather than
advisable here: the committed smoke corpus produces byte-identical output for
ELEVEN of the 24 swept cases, and one of the eleven is
`--symmetric_map_align_thresholds` -- the only piece of clustering logic this
port writes from scratch. A port that got it completely wrong would pass every
end-to-end case on that corpus.

The reference is imported and its functions are WRAPPED, never copied or
edited. If the reference changes, this changes with it. Copying a stage into the
harness is the one thing that guarantees the harness stops measuring the thing
it is named after.

    bench/dump_reference.py --stage minimizers --sorted-fastq OUT/sorted.fastq \
        --k 13 --w 20 > minimizers.tsv

Stages, in dependency order:

    minimizers  get_kmer_minimizers -- (position, minimizer) per read
    mapping     get_best_cluster -- the candidate ranking and the decision
    parasail    parasail_block_alignment -- cigar, and BOTH alignment ratios
    spoa        form_draft_consensus -- the exact sequences handed to spoa, in
                order, and the consensus returned
    identity    highest_aln_identity -- forward and reverse-complement identity
                per center pair, at the consensus path's opening penalty of 3
    barcode     find_barcode_locations -- every edlib HW call and its FULL
                locations list, not just locations[0]

Output format, one line per minimizer, in the order the reference produces them:

    <read index>\t<position>\t<minimizer>

Read index is the 0-based position in sorted.fastq, which is the order
`reads_to_clusters` iterates. Reads the reference skips (homopolymer-compressed
length < k) emit a single line with position -1 and an empty minimizer, so the
skip itself is part of the contract rather than an absence.
"""
import argparse
import itertools
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from modules import help_functions  # noqa: E402
from modules import cluster  # noqa: E402


def dump_minimizers(sorted_fastq, k, w, out):
    for idx, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(sorted_fastq))):
        # Exactly what reads_to_clusters does before calling the function.
        seq_hpol_comp = "".join(ch for ch, _ in itertools.groupby(seq))
        if len(seq_hpol_comp) < k:
            out.write("{0}\t-1\t\n".format(idx))
            continue
        for m, pos in cluster.get_kmer_minimizers(seq_hpol_comp, k, w):
            out.write("{0}\t{1}\t{2}\n".format(idx, pos, m))


def dump_mapping(args, out):
    """Record every get_best_cluster call the real driver makes.

    `cluster.get_best_cluster` is wrapped, not reimplemented: the sweep is
    stateful -- the minimizer database grows as reads become representatives --
    so the only faithful way to capture its inputs is to let the reference run
    and watch. This is PORTING.md method point 3, dumping from the live driver.

    Format, per call:

        CALL  <read_cl_id> <compressed_seq_len> <n_minimizers> <error_rate_read>
        CAND  <cl_id> <error_rate> <indices,...> <positions,...> <acc>
        ...one CAND per candidate, in the reference's dict order...
        RES   <best_cluster_id> <nr_shared_kmers> <mapped_ratio>

    Floats are written with repr() so the replay reads back the identical double.
    """
    import argparse as _argparse
    from modules import p_minimizers_shared

    real = cluster.get_best_cluster
    state = {"n": 0, "stop": False}

    def wrapper(read_cl_id, compressed_seq_len, hit_clusters_ids,
                hit_clusters_hit_positions, minimizers, nummber_of_minimizers,
                hit_clusters_hit_index, representatives, p_emp_probs, a):
        res = real(read_cl_id, compressed_seq_len, hit_clusters_ids,
                   hit_clusters_hit_positions, minimizers, nummber_of_minimizers,
                   hit_clusters_hit_index, representatives, p_emp_probs, a)
        if not state["stop"]:
            out.write("CALL\t{0}\t{1}\t{2}\t{3!r}\n".format(
                read_cl_id, compressed_seq_len, nummber_of_minimizers,
                # index 6 is the homopolymer-compressed error rate. NGSpeciesID's
                # representative tuple has EIGHT elements, not isONclust's seven
                # -- index 7 is the compressed sequence, added for
                # --symmetric_map_align_thresholds -- but the error rate is at 6
                # in both.
                representatives[read_cl_id][6]))
            for cl_id in hit_clusters_hit_positions:
                out.write("CAND\t{0}\t{1!r}\t{2}\t{3}\t{4}\n".format(
                    cl_id, representatives[cl_id][6],
                    ",".join(str(x) for x in hit_clusters_hit_index[cl_id]),
                    ",".join(str(x) for x in hit_clusters_hit_positions[cl_id]),
                    representatives[cl_id][2]))
            out.write("RES\t{0}\t{1}\t{2!r}\n".format(res[0], res[1], res[2]))
            state["n"] += 1
            if args.max_calls and state["n"] >= args.max_calls:
                state["stop"] = True
        return res

    cluster.get_best_cluster = wrapper
    try:
        p_min_shared = p_minimizers_shared.read_empirical_p()
        p_emp_probs = {}
        for k, w, p, e1, e2 in p_min_shared:
            if int(k) == args.k and abs(int(w) - args.w) <= 2:
                p_emp_probs[(float(e1), float(e2))] = float(p)
                p_emp_probs[(float(e2), float(e1))] = float(p)

        read_array = [(i, 0, acc, seq, qual, float(acc.split("_")[-1]))
                      for i, (acc, (seq, qual))
                      in enumerate(help_functions.readfq(open(args.sorted_fastq)))]
        clusters, representatives = {}, {}
        for i, b_i, acc, seq, qual, score in read_array:
            clusters[i] = [acc]
            representatives[i] = (i, b_i, acc, seq, qual, score)

        a = _driver_namespace(args)
        cluster.reads_to_clusters(clusters, representatives, read_array,
                                  p_emp_probs, {}, 1, a)
    finally:
        cluster.get_best_cluster = real


def dump_parasail(args, out):
    """Record every parasail alignment the real driver performs.

    `cluster.parasail_block_alignment` is wrapped, so what is captured is what
    the clustering path actually asks parasail for -- match 2, mismatch -2,
    gap_ext 1, and an opening penalty chosen per comparison from the two reads'
    summed error rates. isONform verified its parasail port against
    isONcorrect's parameters (match 4, mismatch -8, open 12); different
    penalties can reach different tie-breaking paths, so this checks this call
    site rather than trusting a module that was verified elsewhere.

    NOTE this is only ONE of the two parasail call sites in NGSpeciesID. The
    other is `consensus.parasail_alignment`, at a FIXED opening penalty of 3,
    and it is covered by `--stage identity`.

    Format, one record per call:

        PARA\t<opening_penalty>\t<k>\t<match_id>\t<s1>\t<s2>\t<cigar>\t<alignment_ratio>\t<target_alignment_ratio>
    """
    import argparse as _argparse
    from modules import p_minimizers_shared

    real = cluster.parasail_block_alignment
    state = {"n": 0}

    def wrapper(s1, s2, k, match_id, match_score=2, mismatch_penalty=-2,
                opening_penalty=5, gap_ext=1):
        # Recompute the cigar the same way the reference does, so the dump
        # carries it without changing what the driver receives.
        import parasail as _p
        m = _p.matrix_create("ACGT", match_score, mismatch_penalty)
        r = _p.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, m)
        if r.saturated:
            r = _p.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, m)
        cig = str(r.cigar.decode, "utf-8")
        res = real(s1, s2, k, match_id, match_score, mismatch_penalty,
                   opening_penalty, gap_ext)
        # NGSpeciesID's inner tuple has FOUR elements where isONclust's has
        # three: (s1_alignment, s2_alignment, alignment_ratio,
        # target_alignment_ratio). The fourth is what
        # --symmetric_map_align_thresholds reads, so it MUST be dumped -- a port
        # that computes it over the wrong sequence would otherwise pass this
        # oracle and fail only in the one case the smoke corpus cannot see.
        ratio, target_ratio = res[2][2], res[2][3]
        if not (args.max_calls and state["n"] >= args.max_calls):
            out.write("PARA\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6!r}\t{7!r}\n".format(
                opening_penalty, k, match_id, s1, s2, cig, ratio, target_ratio))
            state["n"] += 1
        return res

    cluster.parasail_block_alignment = wrapper
    try:
        p_min_shared = p_minimizers_shared.read_empirical_p()
        p_emp_probs = {}
        for k, w, p, e1, e2 in p_min_shared:
            if int(k) == args.k and abs(int(w) - args.w) <= 2:
                p_emp_probs[(float(e1), float(e2))] = float(p)
                p_emp_probs[(float(e2), float(e1))] = float(p)
        read_array = [(i, 0, acc, seq, qual, float(acc.split("_")[-1]))
                      for i, (acc, (seq, qual))
                      in enumerate(help_functions.readfq(open(args.sorted_fastq)))]
        clusters, representatives = {}, {}
        for i, b_i, acc, seq, qual, score in read_array:
            clusters[i] = [acc]
            representatives[i] = (i, b_i, acc, seq, qual, score)
        a = _driver_namespace(args)
        cluster.reads_to_clusters(clusters, representatives, read_array,
                                  p_emp_probs, {}, 1, a)
    finally:
        cluster.parasail_block_alignment = real



def _driver_namespace(args):
    """The argparse Namespace `reads_to_clusters` reads, and nothing more.

    Built here rather than inline in each stage so the two driver-wrapping
    stages cannot drift apart -- and because it has to grow every time the
    reference reads a new attribute. `symmetric_map_align_thresholds` is the one
    NGSpeciesID added and isONclust does not have; leaving it out raises
    AttributeError from inside get_best_cluster, which reads like a harness bug
    rather than a missing field.
    """
    import argparse as _argparse
    return _argparse.Namespace(
        k=args.k, w=args.w, min_shared=args.min_shared,
        mapped_threshold=args.mapped_threshold,
        aligned_threshold=args.aligned_threshold,
        min_fraction=args.min_fraction,
        min_prob_no_hits=args.min_prob_no_hits,
        symmetric_map_align_thresholds=args.symmetric_map_align_thresholds,
        # print_output is a modulo divisor and 0 raises ZeroDivisionError
        # (PORTING.md, Finding 11). A huge value is how you turn the progress
        # table off without reaching that.
        print_output=10 ** 9)


def dump_spoa(args, out):
    """Record every spoa invocation form_draft_consensus makes.

    This is the oracle that decides whether `spoars` can replace the spoa
    binary. What matters is not only the consensus that comes back but the
    EXACT list of sequences handed in, in order: sequence insertion order into
    a POA graph changes the consensus, so a port that picks the same reads in a
    different order produces different output and is not a port.

    `consensus.run_spoa` is wrapped rather than reimplemented, so the reads file
    it is about to hand to the subprocess is read back and dumped verbatim. That
    also captures the --max_seqs_for_consensus cutoff, which is `i >=
    max_seqs_for_consensus` -- admitting exactly that many sequences, unlike
    isONcorrect's bare `>`.

    THE QUALITY STRING IS PART OF THE INPUT. `run_spoa` hands spoa a FASTQ, and
    spoa's CLI weights the graph by per-base quality whenever the input has any:

        if (it->quality.empty()) graph.AddAlignment(alignment, it->data);
        else                     graph.AddAlignment(alignment, it->data, it->quality);

    with weight = ord(q) - 33. Nothing in run_spoa's argument list says so, and
    an earlier version of this dump recorded only the sequences -- which made
    every comparison against a Rust POA fail for a reason that had nothing to do
    with the POA. Measured: the same 20 sequences give an 847 bp consensus as
    FASTQ and 860 bp as FASTA.

    isONcorrect passes a FASTA, so its spoars validation genuinely does not
    cover this.

    Format, per invocation:

        SPOA\t<n_seqs>\t<consensus>
        SEQ\t<i>\t<accession>\t<sequence>\t<quality>
        ...one SEQ per sequence, in the order written to the temp file...
    """
    from modules import consensus as _cons

    real = _cons.run_spoa
    state = {"n": 0}

    def wrapper(reads, spoa_out_file, spoa_path):
        res = real(reads, spoa_out_file, spoa_path)
        if not (args.max_calls and state["n"] >= args.max_calls):
            seqs = list(help_functions.readfq(open(reads)))
            out.write("SPOA\t{0}\t{1}\n".format(len(seqs), res))
            for i, (acc, (seq, qual)) in enumerate(seqs):
                out.write("SEQ\t{0}\t{1}\t{2}\t{3}\n".format(
                    i, acc, seq, qual if qual is not None else ""))
            state["n"] += 1
        return res

    _cons.run_spoa = wrapper
    try:
        _run_consensus(args)
    finally:
        _cons.run_spoa = real


def dump_identity(args, out):
    """Record every highest_aln_identity call detect_reverse_complements makes.

    This is the SECOND parasail call site and it does not share the clustering
    path's parameters: `consensus.parasail_alignment` defaults to
    opening_penalty=3 where the clustering path bins 2..5 by error rate. It also
    computes identity differently -- by zipping the two gapped strings and
    counting mismatches -- which charges leading and trailing gaps, so a
    semi-global alignment that shifts one base changes the number.

    Both orientations are recorded, not just the max, because the port has to
    get the reverse complement right and `max()` hides which side won.

    THE SEQUENCES ARE RECORDED TOO. An earlier version wrote only the lengths
    and the three identities, which is enough to read but not enough to REPLAY:
    a port cannot be checked against a number whose input it does not have. The
    centers are a few hundred bytes each and there are a handful of pairs, so
    the cost is nothing.

    Format, per call:

        IDENT\t<identity_fw>\t<identity_rc>\t<max>\t<seq>\t<seq2>
    """
    from modules import consensus as _cons

    real_pa = _cons.parasail_alignment
    real_hi = _cons.highest_aln_identity
    state = {"n": 0}

    def wrapper(seq, seq2):
        res = real_hi(seq, seq2)
        if not (args.max_calls and state["n"] >= args.max_calls):
            # Recompute both halves rather than parsing them out of the log, so
            # the dump does not depend on logging configuration.
            rc = _cons.reverse_complement(seq2)
            a1, a2, _c, _t, _s = real_pa(seq, rc)
            id_rc = (len(a1) - sum(1 for x, y in zip(a1, a2) if x != y)) / float(len(a1))
            b1, b2, _c, _t, _s = real_pa(seq, seq2)
            id_fw = (len(b1) - sum(1 for x, y in zip(b1, b2) if x != y)) / float(len(b1))
            out.write("IDENT\t{0!r}\t{1!r}\t{2!r}\t{3}\t{4}\n".format(
                id_fw, id_rc, res, seq, seq2))
            state["n"] += 1
        return res

    _cons.highest_aln_identity = wrapper
    try:
        _run_consensus(args)
    finally:
        _cons.highest_aln_identity = real_hi


def dump_barcode(args, out):
    """Record every edlib HW call find_barcode_locations makes.

    The FULL locations list is dumped, not just `locations[0]` which is all the
    reference reads. That is deliberate: edlib's choice of which equally-optimal
    location comes first is not uniquely defined, and it is exactly the thing a
    native reimplementation has to reproduce. Dumping only the one the reference
    used would make a port that gets the ordering wrong look correct whenever
    the first two happen to coincide.

    THE TARGET IS RECORDED, not just its length. The earlier version wrote
    `len(target)`, which is enough to read and not enough to replay -- the same
    mistake the identity dump made. The windows are `--trim_window` bases each,
    so the cost is small.

    Format, per call:

        EDLIB\t<primer_seq>\t<target>\t<k>\t<edit_distance>\t<loc0_start,loc0_end;...>
    """
    from modules import barcode_trimmer as _bt

    real_align = None
    state = {"n": 0}
    import edlib as _edlib
    real_align = _edlib.align

    def wrapper(query, target, **kw):
        res = real_align(query, target, **kw)
        if kw.get("mode") == "HW" and not (args.max_calls and state["n"] >= args.max_calls):
            locs = res.get("locations") or []
            out.write("EDLIB\t{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                query, target, kw.get("k"), res.get("editDistance"),
                ";".join("{0},{1}".format(a, b) for a, b in locs)))
            state["n"] += 1
        return res

    _edlib.align = wrapper
    try:
        _run_consensus(args)
    finally:
        _edlib.align = real_align


def _run_consensus(args):
    """Drive the reference far enough to reach the consensus stage.

    Deliberately calls the reference's own functions rather than re-deriving
    the clusters, because the ORDER form_draft_consensus walks clusters in --
    (size, representative score) descending -- is part of what is being
    verified, and reimplementing it here would verify this file instead of the
    reference.
    """
    import tempfile
    from modules import consensus as _cons
    from modules import p_minimizers_shared
    from modules import barcode_trimmer as _bt

    p_min_shared = p_minimizers_shared.read_empirical_p()
    p_emp_probs = {}
    for k, w, p, e1, e2 in p_min_shared:
        if int(k) == args.k and abs(int(w) - args.w) <= 2:
            p_emp_probs[(float(e1), float(e2))] = float(p)
            p_emp_probs[(float(e2), float(e1))] = float(p)

    read_array = [(i, 0, acc, seq, qual, float(acc.split("_")[-1]))
                  for i, (acc, (seq, qual))
                  in enumerate(help_functions.readfq(open(args.sorted_fastq)))]
    clusters, representatives = {}, {}
    for i, b_i, acc, seq, qual, score in read_array:
        clusters[i] = [acc]
        representatives[i] = (i, b_i, acc, seq, qual, score)

    a = _driver_namespace(args)
    result = cluster.reads_to_clusters(clusters, representatives, read_array,
                                       p_emp_probs, {}, 1, a)
    clusters, representatives, _, _ = list(result.values())[0]

    # Reassign, exactly as NGSpeciesID's single_clustering path leaves it.
    work_dir = tempfile.mkdtemp()
    a.max_seqs_for_consensus = args.max_seqs_for_consensus
    a.abundance_ratio = args.abundance_ratio
    a.outfolder = work_dir
    a.trim_window = args.trim_window
    a.primer_max_ed = args.primer_max_ed
    abundance_cutoff = int(args.abundance_ratio * len(read_array))
    centers = _cons.form_draft_consensus(clusters, representatives,
                                         args.sorted_fastq, work_dir,
                                         abundance_cutoff, a)
    if args.primer_file:
        barcodes = _bt.read_barcodes(args.primer_file)
        _bt.remove_barcodes(centers, barcodes, a)
    elif args.remove_universal_tails:
        barcodes = _bt.get_universal_tails()
        _bt.remove_barcodes(centers, barcodes, a)
    _cons.detect_reverse_complements(centers, args.rc_identity_threshold)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True,
                    choices=["minimizers", "mapping", "parasail",
                             "spoa", "identity", "barcode"])
    ap.add_argument("--sorted-fastq", required=True,
                    help="the reference's sorted.fastq -- NOT the raw input")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--w", type=int, required=True)
    ap.add_argument("--q", type=float, default=7.0)
    ap.add_argument("--min_shared", type=int, default=5)
    ap.add_argument("--mapped_threshold", type=float, default=0.7)
    ap.add_argument("--aligned_threshold", type=float, default=0.4)
    ap.add_argument("--min_fraction", type=float, default=0.8)
    ap.add_argument("--min_prob_no_hits", type=float, default=0.1)
    ap.add_argument("--symmetric_map_align_thresholds", action="store_true")
    ap.add_argument("--abundance_ratio", type=float, default=0.1)
    ap.add_argument("--rc_identity_threshold", type=float, default=0.9)
    ap.add_argument("--max_seqs_for_consensus", type=int, default=-1)
    ap.add_argument("--trim_window", type=int, default=150)
    ap.add_argument("--primer_max_ed", type=int, default=2)
    ap.add_argument("--primer_file", default="")
    ap.add_argument("--remove_universal_tails", action="store_true")
    ap.add_argument("--max-calls", type=int, default=0,
                    help="stop after this many recorded calls (0 = no limit)")
    ap.add_argument("--out", default=None,
                    help="write the dump here instead of stdout. USE THIS: the "
                         "reference prints progress to stdout from inside "
                         "reads_to_clusters, so a stdout dump is interleaved "
                         "with 5 stray lines like 'Saved: 0 iterations.'")
    args = ap.parse_args()
    out = open(args.out, "w") if args.out else sys.stdout
    try:
        if args.stage == "minimizers":
            dump_minimizers(args.sorted_fastq, args.k, args.w, out)
        elif args.stage == "mapping":
            dump_mapping(args, out)
        elif args.stage == "parasail":
            dump_parasail(args, out)
        elif args.stage == "spoa":
            dump_spoa(args, out)
        elif args.stage == "identity":
            dump_identity(args, out)
        elif args.stage == "barcode":
            dump_barcode(args, out)
    finally:
        if args.out:
            out.close()


if __name__ == "__main__":
    main()
