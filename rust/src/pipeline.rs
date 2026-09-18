//! The driver: the stages of `NGSpeciesID`'s `main()`, in its order.
//!
//! The reference's `main` is four steps, and this mirrors them one for one:
//!
//! 1. sort every read by expected error-free k-mers (`get_sorted_fastq_for_cluster`)
//! 2. filter by length and subsample (`--m`/`--s`, `--sample_size`/`--top_reads`)
//! 3. load the empirical probability table and cluster
//! 4. write the output, and optionally form consensus sequences
//!
//! **All four exist.** `write_fastq` is the only part of the reference still
//! missing, and it is a separate subcommand that never enters this module.
//!
//! WHAT IS PRINTED, AND WHERE
//! --------------------------
//! Everything the reference emits goes through `logging` with
//! `format='%(message)s'`, so it lands on **stderr** with no prefix, and at the
//! default INFO level a whole successful run prints exactly two lines:
//!
//! ```text
//! Starting Clustering: 274 reads
//! Finished Clustering: 3 clusters formed
//! ```
//!
//! The sorting stage prints nothing at all — its messages are all
//! `logging.debug`. That is why this module does not port them: at INFO level
//! they are not observable, and `--debug` output is not in the byte-identity
//! contract because it carries timings and the entire probability table.
//!
//! The two exceptions, both of which ARE shown at INFO level and both of which
//! are in the goldens:
//!
//! * `logging.warning` when `--use_old_sorted_file` reuses an existing file;
//! * the `logging.error` pair when `--q` filters every read.

use crate::cli::Args;
use crate::{
    blockalign, consensus, fastq, p_emp, packed, parallelize, pyfloat, pyrandom, sorting, sweep,
};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::{Path, PathBuf};

/// How far the pipeline got.
///
/// There is no longer an `Incomplete` variant: every stage of `main()` is
/// implemented. `write_fastq` is the one thing left, and it is a separate
/// subcommand that never enters this module.
pub enum Outcome {
    /// Ran to completion.
    Done,
    /// The reference fails here too, with this exit code. The message has
    /// already been written to stderr.
    Failed(i32),
}

/// The files the reference writes into `--outfolder`, by the names it uses.
pub struct Paths {
    pub sorted: PathBuf,
    pub logfile: PathBuf,
}

impl Paths {
    pub fn in_outfolder(outfolder: &str) -> Self {
        let d = Path::new(outfolder);
        Paths {
            sorted: d.join("sorted.fastq"),
            logfile: d.join("logfile.txt"),
        }
    }
}

pub fn run(args: &Args) -> Outcome {
    let outfolder = args
        .outfolder
        .as_deref()
        .expect("cli::validate rejects a missing --outfolder");
    let paths = Paths::in_outfolder(outfolder);

    let sorted_count = match sort_stage(args, &paths) {
        Err(code) => return Outcome::Failed(code),
        Ok(n) => n,
    };

    let clustered = match cluster_stage(args, &paths, sorted_count) {
        Err(code) => return Outcome::Failed(code),
        Ok(c) => c,
    };
    eprintln!(
        "Finished Clustering: {} clusters formed",
        clustered.nontrivial
    );

    if !args.consensus {
        return Outcome::Done;
    }
    match consensus_stage(args, &paths, &clustered) {
        Err(code) => Outcome::Failed(code),
        Ok(()) => Outcome::Done,
    }
}

/// What the clustering stage hands the consensus stage.
pub struct Clustered {
    /// Clusters with more than one read -- the number reported as "formed".
    pub nontrivial: usize,
    /// `(c_id, member accessions)` in the reference's walk order:
    /// `(size, representative score)` descending.
    ///
    /// `c_id` is the read's **ordinal in sorted.fastq**, which is what the
    /// reference uses and what ends up in `consensus_reference_<c_id>.fasta`.
    /// Members are in the cluster's STORED order -- the sweep's insertion
    /// order, not the score-descending order `final_clusters.tsv` uses.
    pub clusters: Vec<(usize, Vec<String>)>,
    /// The read count `--abundance_ratio` is applied to: the POST-subsample
    /// count, not the input's (*Finding 20*).
    pub nr_reads: usize,
}

/// Step 4: `--consensus`.
///
/// The reference's order, which the re-trim loop makes less obvious than it
/// looks:
///
/// 1. `form_draft_consensus`
/// 2. trim, if a primer file or `--remove_universal_tails` was given
/// 3. `detect_reverse_complements`
/// 4. `polish_sequences`
/// 5. **trim again**, and if that changed anything, redo 3 and 4
///
/// The final count reported is `len(centers_filtered)` — the output of step 3,
/// not step 4 — which matters only in that the two are the same length.
fn consensus_stage(args: &Args, paths: &Paths, clustered: &Clustered) -> Result<(), i32> {
    eprintln!("Starting Consensus creation and polishing");
    let outfolder = Path::new(args.outfolder.as_deref().expect("validated"));

    // `int(args.abundance_ratio * len(read_array))`. Truncates, so below ten
    // reads the default 0.1 gives a cutoff of 0 and every cluster qualifies,
    // singletons included (*Finding 20*).
    let abundance_cutoff = (args.abundance_ratio * clustered.nr_reads as f64) as usize;

    // `tempfile.mkdtemp()`: the per-cluster read files live outside --outfolder
    // and are removed at the end, so they are not part of the output contract.
    let work_dir = match tempdir() {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error: cannot create a temporary directory: {e}");
            return Err(1);
        }
    };

    // Every read in sorted.fastq, by accession. The reference slurps the same
    // dict; the centers need arbitrary lookup, not a stream.
    let mut by_acc: FxHashMap<String, (String, String)> = FxHashMap::default();
    if let Err(e) = fastq::for_each_file(&paths.sorted, |r| {
        by_acc.insert(r.name, (r.seq, r.qual.unwrap_or_default()));
    }) {
        eprintln!("Error: cannot read {}: {e}", paths.sorted.display());
        return Err(1);
    }
    let read_of = |acc: &str| by_acc.get(acc).cloned();

    let mut centers = match consensus::form_draft_consensus(
        &clustered.clusters,
        &read_of,
        &work_dir,
        abundance_cutoff,
        args.max_seqs_for_consensus,
    ) {
        Ok(c) => c,
        Err(msg) => {
            eprintln!("Error: {msg}");
            let _ = std::fs::remove_dir_all(&work_dir);
            return Err(1);
        }
    };

    // The barcodes, if any. Read once and reused by both trim passes, as the
    // reference does.
    let barcodes = match load_barcodes(args) {
        Ok(b) => b,
        Err(msg) => {
            eprintln!("Error: {msg}");
            let _ = std::fs::remove_dir_all(&work_dir);
            return Err(1);
        }
    };
    if let Some(b) = &barcodes {
        consensus::remove_barcodes(&mut centers, b, args.trim_window, args.primer_max_ed);
    }

    let polisher = consensus::Polisher::of(args).expect("cli::validate requires one");
    let mut filtered = consensus::detect_reverse_complements(centers, args.rc_identity_threshold);
    if let Err(msg) = consensus::polish_sequences(&mut filtered, outfolder, polisher, args) {
        eprintln!("Error: {msg}");
        let _ = std::fs::remove_dir_all(&work_dir);
        return Err(1);
    }

    // The second trim, and the redo it can trigger. Commit 5463966 added the
    // `if centers_updated` guard so medaka is not run twice for nothing.
    if let Some(b) = &barcodes {
        let updated =
            consensus::remove_barcodes(&mut filtered, b, args.trim_window, args.primer_max_ed);
        if updated {
            filtered = consensus::detect_reverse_complements(filtered, args.rc_identity_threshold);
            if let Err(msg) = consensus::polish_sequences(&mut filtered, outfolder, polisher, args)
            {
                eprintln!("Error: {msg}");
                let _ = std::fs::remove_dir_all(&work_dir);
                return Err(1);
            }
        }
    }

    let _ = std::fs::remove_dir_all(&work_dir);
    eprintln!("Finished Consensus creation: {} created", filtered.len());
    Ok(())
}

/// `--primer_file` or `--remove_universal_tails`, or neither. argparse makes
/// them mutually exclusive.
fn load_barcodes(args: &Args) -> Result<Option<Vec<(String, String)>>, String> {
    if args.remove_universal_tails {
        return Ok(Some(consensus::universal_tails()));
    }
    if args.primer_file.is_empty() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&args.primer_file)
        .map_err(|e| format!("cannot read {}: {e}", args.primer_file))?;
    let mut records: Vec<(String, String)> = Vec::new();
    fastq::for_each(text.split_inclusive('\n'), |r| {
        records.push((r.name.clone(), r.seq.clone()));
    });
    Ok(Some(consensus::read_barcodes(&records)))
}

/// `tempfile.mkdtemp()`. The path varies per run and reaches no output file.
fn tempdir() -> std::io::Result<PathBuf> {
    let base = std::env::temp_dir();
    for _ in 0..64 {
        let n: u64 = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos() as u64)
            .unwrap_or(0)
            ^ (std::process::id() as u64) << 32;
        let p = base.join(format!("ngspeciesid{n:x}"));
        match std::fs::create_dir(&p) {
            Ok(()) => return Ok(p),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::other(
        "could not create a temporary directory",
    ))
}

/// The score the sorting stage appended, recovered with
/// `float(acc.split("_")[-1])`.
fn score_of(acc: &str) -> f64 {
    acc.rsplit('_')
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(f64::NAN)
}

/// `"_".join(acc.split("_")[:-1])` -- drop the score suffix on the way out.
///
/// Note what this does to an accession with no underscore at all: the
/// reference's join of an empty list is `""`, not the original string. Faithful.
pub fn strip_score(acc: &str) -> &str {
    match acc.rfind('_') {
        Some(i) => &acc[..i],
        None => "",
    }
}

/// Steps 2 and 3: filter, subsample, then cluster and write the output.
///
/// Returns what the consensus stage needs: the cluster walk order, the member
/// accessions, the post-subsample read count, and the number of non-trivial
/// clusters the reference reports as "clusters formed".
fn cluster_stage(args: &Args, paths: &Paths, sorted_count: usize) -> Result<Clustered, i32> {
    // The reference re-reads sorted.fastq rather than reusing anything from the
    // sort, and recovers the score by parsing it back out of the accession. So
    // the accession the sweep sees is the one WITH the score suffix, and the
    // f64 is whatever `float()` makes of the text `str()` produced. Reproduced
    // rather than shortcut: that round trip through a decimal string is part of
    // the contract (PORTING.md, "The pipeline, in one pass", step 2).
    let mut reads: Vec<sweep::SweepRead> = Vec::with_capacity(sorted_count);
    // Dense read index -> ordinal in sorted.fastq. Identity until `--m`/`--s` or
    // a subsample drops something, and then NOT: the quality strings for the
    // representatives are recovered by streaming sorted.fastq, and looking them
    // up by dense index fetches the wrong read.
    //
    // Found by the goldens. `m750s50` produced a 649-character quality string
    // against the reference's 725, on a cluster whose sequence, score and error
    // rate were all correct -- so `final_clusters.tsv` matched and only column 4
    // of `final_cluster_origins.tsv` moved.
    let mut ordinal_of: Vec<usize> = Vec::with_capacity(sorted_count);
    let mut ordinal = 0usize;
    // 2-bit packing cannot represent a fifth symbol, so a non-ACGT base is read
    // as `A`, which changes minimizer selection and the alignment path. Every
    // committed corpus is pure ACGT, so this stays zero -- and a corpus that is
    // not would diverge silently without this count.
    let mut substituted = 0usize;
    // The length filter. `--m 0 --s 0` -- the default -- disables it entirely,
    // and the reference tests `> 0` on BOTH, so `--m 800 --s 0` is also
    // disabled. Faithful.
    let filtering = args.target_length > 0 && args.target_deviation > 0;
    let lo = args.target_length - args.target_deviation;
    let hi = args.target_length + args.target_deviation;

    if let Err(e) = fastq::for_each_file(&paths.sorted, |r| {
        let this_ordinal = ordinal;
        ordinal += 1;
        if filtering {
            let n = r.seq.len() as i64;
            if n < lo || n > hi {
                return;
            }
        }
        ordinal_of.push(this_ordinal);
        let score = score_of(&r.name);
        let (seq, sub) = packed::PackedSeq::from_bytes(r.seq.as_bytes());
        substituted += sub;
        let qual = r.qual.unwrap_or_default();
        let qual_b = qual.as_bytes();
        let hp_error_rate = sweep::compressed_error_rate(r.seq.as_bytes(), qual_b);
        let err_per_base = blockalign::expected_errors(qual_b) / r.seq.len() as f64;
        reads.push(sweep::SweepRead {
            // The reference keeps `i` from `enumerate` over ALL reads, so its
            // cluster ids are sparse when the filter drops something. Using the
            // dense position instead is safe and is checked by the goldens: the
            // only place the id reaches output is the third sort key, and the
            // filter preserves order, so sparse and dense ids sort identically.
            id: reads.len(),
            prev_batch_index: 0,
            acc: r.name.into(),
            seq,
            hp_error_rate,
            err_per_base,
            score,
        });
    }) {
        eprintln!("Error: cannot read {}: {e}", paths.sorted.display());
        return Err(1);
    }
    if substituted > 0 {
        eprintln!(
            "Warning: {substituted} non-ACGT bases were read as 'A'. Sequences are 2-bit packed,\n\
             which cannot represent a fifth symbol, so output for this input will NOT match the\n\
             Python reference. See rust/src/packed.rs."
        );
    }

    // The subsample. Two paths, and the reference's guards are reproduced
    // exactly because both edges are observable:
    //
    // * `--top_reads` takes the head of a file already sorted by score
    //   descending, and its guard is a bare `if args.top_reads` -- so
    //   `--top_reads` WITHOUT `--sample_size` truncates to zero reads. Faithful.
    // * otherwise the guard is `0 < sample_size < len(read_array)`, so a size at
    //   or above the surviving read count silently takes every read rather than
    //   erroring. `top_huge` in the case matrix is that no-op path on purpose.
    if args.top_reads {
        reads.truncate(args.sample_size.max(0) as usize);
        ordinal_of.truncate(reads.len());
    } else if args.sample_size > 0 && (args.sample_size as usize) < reads.len() {
        // `sorted(random.Random(seed).sample(range(n), k))`. The sort is the
        // reference's, and it is what keeps the subsample in score order --
        // which matters, because the sweep is order-dependent.
        let mut rng = pyrandom::PyRandom::seeded(args.seed);
        let mut picked = pyrandom::sample_indices(&mut rng, reads.len(), args.sample_size as usize);
        picked.sort_unstable();
        let mut kept_reads = Vec::with_capacity(picked.len());
        let mut kept_ordinals = Vec::with_capacity(picked.len());
        for i in &picked {
            kept_reads.push(reads[*i].clone());
            kept_ordinals.push(ordinal_of[*i]);
        }
        reads = kept_reads;
        ordinal_of = kept_ordinals;
        // The sweep indexes `reads` by `SweepRead::id`, so the ids have to be
        // the new dense positions. The reference keeps the ORIGINAL enumerate
        // index here and its cluster ids therefore have gaps -- which reaches
        // output only through the third sort key, where relative order is all
        // that matters, and the subsample preserves it. Same argument as the
        // length filter; the goldens are what check it.
        for (i, r) in reads.iter_mut().enumerate() {
            r.id = i;
        }
    }

    let nr_reads = reads.len();
    eprintln!("Starting Clustering: {nr_reads} reads");

    // The probability table. `None` is Finding 8: the reference builds an empty
    // dict and dies with a KeyError on a tuple of two floats the first time it
    // looks anything up. 3 361 CLI-valid (k, w) pairs reach it.
    let Some(table) = p_emp::Table::select(args.k, args.w) else {
        eprintln!(
            "Error: no empirical probability table for --k {} --w {}.",
            args.k, args.w
        );
        eprintln!("The table covers --k 10 to 30; --w must be within 2 of a value it stores.");
        return Err(1);
    };

    let params = sweep::SweepParams {
        k: args.k as usize,
        w: args.w as usize,
        min_shared: args.min_shared,
        min_fraction: args.min_fraction,
        min_prob_no_hits: args.min_prob_no_hits,
        mapped_threshold: args.mapped_threshold,
        aligned_threshold: args.aligned_threshold,
        symmetric_map_align_thresholds: args.symmetric_map_align_thresholds,
    };

    let outfolder = Path::new(args.outfolder.as_deref().expect("validated"));
    if args.nr_cores > 1 {
        // `--t` does NOT parallelise the sweep -- it REPLACES it. Reads are cut
        // into batches, each clustered independently with its own minimizer
        // database, and the survivors are pooled and re-clustered until one
        // batch remains. Every `--t` value is its own answer: 49 / 42 / 35 / 33
        // clusters at 1 / 2 / 4 / 8 on the 3 000-read corpus. PORTING.md,
        // "With --t > 1 this is a different algorithm".
        //
        // `--batch_type` has already been validated by cli::validate, so the
        // parse below cannot fail; the reference reaches a ValueError here
        // instead (Finding 5).
        let bt = parallelize::BatchType::parse(&args.batch_type)
            .expect("cli::validate rejects an unknown --batch_type");
        let r = parallelize::parallel_clustering(
            &reads,
            args.nr_cores as usize,
            bt,
            &table,
            params,
            &paths.sorted,
            outfolder,
        );
        if let Some(msg) = &r.intermediate_error {
            eprintln!("Error: {msg}");
            return Err(1);
        }
        // Before anything else the caller would do, because the reference dies
        // here with an empty output folder. See `finding_18`.
        if let Some(c_id) = r.short_rep {
            finding_18(c_id);
            return Err(1);
        }
        let res = sweep::SweepResult {
            clusters: r.clusters,
            representatives: r.representatives,
            db: crate::cluster::MinimizerDatabase::new(),
            batch_index: 0,
            mapped_passed: r.mapped_passed,
            aln_passed: r.aln_passed,
            aln_called: r.aln_called,
            skipped_short: r.skipped_short,
            times: r.times,
        };
        let nr_reads = reads.len();
        let mut out = write_output(args, paths, &reads, &ordinal_of, &res)?;
        out.nr_reads = nr_reads;
        return Ok(out);
    }

    let clusters = sweep::OrderedClusters::default();
    let reps: FxHashMap<usize, sweep::ReadInfo> = FxHashMap::default();
    let db = crate::cluster::MinimizerDatabase::default();
    let res = sweep::reads_to_clusters(clusters, reps, &reads, db, 1, &table, params);

    let nr_reads = reads.len();
    let mut out = write_output(args, paths, &reads, &ordinal_of, &res)?;
    out.nr_reads = nr_reads;
    Ok(out)
}

/// Step 4: `final_clusters.tsv` and `final_cluster_origins.tsv`.
///
/// Ordered by `(cluster size, representative score)` descending. Python's
/// `sorted(..., reverse=True)` is stable and does **not** reverse ties, so equal
/// pairs keep dict insertion order -- which after the reassignment step is
/// ascending cluster id. Hence the third key.
fn write_output(
    args: &Args,
    paths: &Paths,
    reads: &[sweep::SweepRead],
    ordinal_of: &[usize],
    res: &sweep::SweepResult,
) -> Result<Clustered, i32> {
    let clusters = &res.clusters;
    let representatives = &res.representatives;

    let mut order: Vec<usize> = clusters.order.clone();
    order.sort_by(|a, b| {
        let ka = (clusters.map[a].len(), representatives[a].score);
        let kb = (clusters.map[b].len(), representatives[b].score);
        kb.0.cmp(&ka.0)
            .then(kb.1.partial_cmp(&ka.1).expect("scores are finite"))
            .then(a.cmp(b))
    });

    // The representatives' quality strings, which the clustering stage does not
    // keep resident. One sequential pass over sorted.fastq for the survivors.
    // By ORDINAL IN sorted.fastq, not by dense read index. See `ordinal_of`.
    let want: FxHashSet<usize> = order.iter().map(|c| ordinal_of[*c]).collect();
    let rep_quals = match quals_for(&paths.sorted, &want) {
        Ok(q) => q,
        Err(e) => {
            eprintln!("Error: cannot re-read {}: {e}", paths.sorted.display());
            return Err(1);
        }
    };

    let mut clusters_out = String::new();
    let mut origins_out = String::new();
    let mut nontrivial = 0usize;
    // The walk order and member lists the consensus stage needs. Members stay in
    // the cluster's STORED order here, NOT the score-descending order written to
    // final_clusters.tsv -- `form_draft_consensus` iterates the stored one, and
    // the order reads enter a POA graph changes the consensus.
    let mut walk: Vec<(usize, Vec<String>)> = Vec::with_capacity(order.len());
    // Set when a representative is still a six-element tuple. See below.
    let mut short_rep: Option<usize> = None;
    for (output_cl_id, c_id) in order.iter().enumerate() {
        let rep = &representatives[c_id];
        // FINDING 18, reproduced including the partial output.
        //
        // `reads_to_clusters` skips a read whose homopolymer-compressed length
        // is under `--k` with a bare `continue`, BEFORE the block that grows its
        // representative from six elements to eight. The read stays its own
        // cluster, and this loop's
        //
        //     read_cl_id, b_i, acc, c_seq, c_qual, score, error_rate, _ = ...
        //
        // raises `ValueError: not enough values to unpack (expected 8, got 6)`.
        // `error_rate == None` here IS that six-element tuple.
        //
        // Unreachable in a single run -- the sort stage drops those reads using
        // the same `--k` -- and two ordinary commands away otherwise:
        //
        //     NGSpeciesID --fastq test/sample_h1.fastq --outfolder out --t 1 --k 13 --w 20
        //     NGSpeciesID --use_old_sorted_file        --outfolder out --t 1 --k 25 --w 50
        //
        // Both files are already open and have been written to, and CPython
        // flushes them at interpreter shutdown, so the run leaves COMPLETE
        // records for every earlier cluster and nothing for this one: measured
        // at 310 074 and 61 945 bytes on `sample_h1`. Breaking before appending
        // anything for this cluster reproduces that byte for byte.
        //
        // Writing `nan` instead -- which is what this did, and what
        // `cli/use_old_k_mismatch` caught -- is strictly worse than crashing:
        // the run exits 0 and every downstream consumer sees a cluster whose
        // error rate is not a number.
        if rep.error_rate.is_none() {
            short_rep = Some(*c_id);
            break;
        }
        origins_out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            output_cl_id,
            strip_score(&rep.acc),
            String::from_utf8_lossy(&rep.seq.to_bytes()),
            rep_quals
                .get(&ordinal_of[*c_id])
                .map(String::as_str)
                .unwrap_or(""),
            pyfloat::repr(rep.score),
            pyfloat::repr(rep.error_rate.unwrap_or(f64::NAN)),
        ));
        // Within a cluster, reads are score-descending. The sort is deliberately
        // NOT total: ties keep member-list order, which is ascending read id,
        // because Python's sort is stable.
        let mut members: Vec<u32> = clusters.map[c_id].clone();
        members.sort_by(|a, b| {
            reads[*b as usize]
                .score
                .partial_cmp(&reads[*a as usize].score)
                .expect("scores are finite")
        });
        for id in &members {
            clusters_out.push_str(&format!(
                "{}\t{}\n",
                output_cl_id,
                strip_score(&reads[*id as usize].acc)
            ));
        }
        if clusters.map[c_id].len() > 1 {
            nontrivial += 1;
        }
        walk.push((
            // The ORIGINAL ordinal in sorted.fastq, not the dense index.
            //
            // The reference keeps `i` from `enumerate` over all of sorted.fastq
            // as the cluster id, so its ids are sparse after `--m`/`--s` or a
            // subsample drops something. An earlier comment here claimed that
            // only reached the third sort key, where relative order is all that
            // matters -- and that was WRONG: the consensus stage puts the id in
            // a FILENAME, `consensus_reference_<c_id>.fasta`. `cons_sample100`
            // caught it, writing `..._9.fasta` where the reference writes
            // `..._31.fasta` with byte-identical content.
            //
            // Dense ids stay internal, because the sweep indexes `reads` by
            // them; `ordinal_of` maps back at the one point it escapes. The
            // ordering is unaffected either way, since `ordinal_of` is
            // monotonically increasing.
            ordinal_of[*c_id],
            clusters.map[c_id]
                .iter()
                .map(|id| reads[*id as usize].acc.to_string())
                .collect(),
        ));
    }

    let outfolder = Path::new(args.outfolder.as_deref().expect("validated"));
    let cp = outfolder.join("final_clusters.tsv");
    let op = outfolder.join("final_cluster_origins.tsv");
    if let Err(e) = std::fs::write(&cp, clusters_out).and_then(|_| std::fs::write(&op, origins_out))
    {
        eprintln!("Error: cannot write output: {e}");
        return Err(1);
    }
    // After the partial write, not before: the reference's files carry what it
    // managed to write before the traceback.
    if let Some(c_id) = short_rep {
        finding_18(c_id);
        return Err(1);
    }
    Ok(Clustered {
        nontrivial,
        clusters: walk,
        // Filled in by the caller, which knows the post-subsample count.
        nr_reads: 0,
    })
}

/// Byte range of each wanted record in sorted.fastq.
///
/// The `--t > 1` path writes an intermediate per merge iteration, each needing
/// the surviving representatives' quality strings. Streaming the whole file once
/// per iteration would be the obvious thing and is not what this does: the
/// ranges are found once, and each iteration reads only the bytes it needs.
pub fn record_ranges_for(
    sorted_path: &Path,
    want: &FxHashSet<usize>,
) -> std::io::Result<FxHashMap<usize, (u64, u32)>> {
    let mut out: FxHashMap<usize, (u64, u32)> = FxHashMap::default();
    out.reserve(want.len());
    let mut ordinal = 0usize;
    fastq::for_each_file_indexed(sorted_path, |_r, at, n| {
        let this = ordinal;
        ordinal += 1;
        if want.contains(&this) {
            out.insert(this, (at, n));
        }
    })?;
    Ok(out)
}

/// One record's quality string, read from its byte range.
pub fn qual_at(
    src: &std::fs::File,
    at: u64,
    n: u32,
    raw: &mut Vec<u8>,
) -> std::io::Result<Option<String>> {
    use std::os::unix::fs::FileExt;
    raw.resize(n as usize, 0);
    src.read_exact_at(raw, at)?;
    let text = String::from_utf8_lossy(raw);
    let mut qual = None;
    // Re-parsed with the parser that produced the range, so the bytes get one
    // interpretation and not two.
    fastq::for_each(text.split_inclusive('\n'), |r| {
        qual = Some(r.qual.unwrap_or_default())
    });
    Ok(qual)
}

/// Quality strings for a set of read ids, recovered by streaming sorted.fastq.
/// A read's id is its ordinal there, so one sequential pass finds them all.
fn quals_for(
    sorted_path: &Path,
    want: &FxHashSet<usize>,
) -> std::io::Result<FxHashMap<usize, String>> {
    let mut out: FxHashMap<usize, String> = FxHashMap::default();
    let mut ordinal = 0usize;
    fastq::for_each_file(sorted_path, |r| {
        let this = ordinal;
        ordinal += 1;
        if want.contains(&this) {
            out.insert(this, r.qual.unwrap_or_default());
        }
    })?;
    Ok(out)
}

/// Step 1: `get_sorted_fastq_for_cluster.main`.
///
/// Returns the number of reads that passed the quality filter, or the exit code
/// to fail with.
///
/// TWO ORDERING FACTS FROM THE REFERENCE, both observable:
///
/// * `logfile.txt` is opened for **writing** before the branch that decides
///   whether to write anything to it. So `--use_old_sorted_file` truncates it
///   to zero bytes even though the clustering is unchanged — PORTING.md,
///   *Finding 13*. Reproduced.
/// * `sorted.fastq` is likewise opened for writing before the code that fills
///   it, so a run that fails leaves an **empty** one behind, and a second run
///   with `--use_old_sorted_file` then "succeeds" on zero reads — *Finding 23*.
///   Reproduced, because a pipeline that retries after a failure has to see the
///   same thing it sees today.
fn sort_stage(args: &Args, paths: &Paths) -> Result<usize, i32> {
    // Truncate the logfile first, exactly as the reference's
    // `open(..., 'w')` does, and before anything can return early.
    if let Err(e) = std::fs::write(&paths.logfile, b"") {
        eprintln!("Error: cannot write {}: {e}", paths.logfile.display());
        return Err(1);
    }

    if paths.sorted.is_file() && args.use_old_sorted_file {
        // logging.warning, so it IS shown at the default level.
        eprintln!(
            "Using already existing sorted file in specified directory, in not intended, specify different outfolder or delete the current file."
        );
        return Ok(count_records(&paths.sorted));
    }

    let Some(input) = args.fastq.as_deref() else {
        // `--use_old_sorted_file` with no sorted.fastq to reuse. The reference
        // reaches `for i, (...) in enumerate(read_array)` with read_array
        // unbound and dies with UnboundLocalError. Finding 23's first half.
        //
        // Note it has ALREADY created an empty sorted.fastq by this point --
        // the reference opens it unconditionally further down -- which is what
        // makes the retry exit 0. Reproduce that, since a pipeline that retries
        // must see what it sees today.
        let _ = std::fs::write(&paths.sorted, b"");
        eprintln!(
            "Error: --use_old_sorted_file was given but {} does not exist.",
            paths.sorted.display()
        );
        return Err(1);
    };

    sort_from_fastq(args, paths, input)
}

/// Score, filter and sort, then write `sorted.fastq` and `logfile.txt`.
///
/// Two streaming passes, carrying 24 bytes per read rather than its bases:
/// pass one scores every record and notes where it sits in the input and how
/// long its output line will be; the vector is sorted, which fixes every
/// surviving record's byte offset in the output; pass two streams the input
/// again and writes each record straight to its place. Carried across from the
/// isONclust port, where holding `acc`/`seq`/`qual` for every read was the
/// largest allocation in the whole program.
fn sort_from_fastq(args: &Args, paths: &Paths, input: &str) -> Result<usize, i32> {
    let k = args.k as usize;

    struct SortRec {
        score: f64,
        error_rate: f64,
        ordinal: u32,
        out_len: u32,
    }
    let mut recs: Vec<SortRec> = Vec::new();
    // Finding 12: a record with no quality -- which is what a fastq with no
    // trailing newline produces -- makes the reference die with
    // `TypeError: 'NoneType' object is not iterable`. Keep the first and report
    // it, rather than a stack trace.
    let mut no_qual: Option<String> = None;
    let mut ordinal: u32 = 0;
    if let Err(e) = fastq::for_each_file(Path::new(input), |r| {
        let this = ordinal;
        ordinal += 1;
        if r.qual.is_none() && no_qual.is_none() {
            no_qual = Some(r.name.clone());
        }
        if let Some(sc) = sorting::score_record(&r, k, args.quality_threshold) {
            let qual = r.qual.as_deref().unwrap_or("");
            recs.push(SortRec {
                score: sc.score,
                error_rate: sc.error_rate,
                ordinal: this,
                out_len: sorting::sorted_fastq_record_len(&r.name, sc.score, &r.seq, qual) as u32,
            });
        }
    }) {
        eprintln!("Error: cannot read {input}: {e}");
        return Err(1);
    }
    if let Some(name) = no_qual {
        eprintln!("Error: read '{name}' has no quality values.");
        eprintln!("The usual cause is a fastq with no trailing newline on its last line.");
        return Err(1);
    }

    let nr_scored = recs.len();
    sorting::sort_by_score(&mut recs, |r| r.score);

    let mut place: Vec<u64> = vec![u64::MAX; ordinal as usize];
    let mut total: u64 = 0;
    for r in &recs {
        place[r.ordinal as usize] = total;
        total += u64::from(r.out_len);
    }
    // The score is needed again in pass two to rebuild the header, and
    // recomputing it would mean a second compensated sum over the quality.
    let mut score_of_ordinal: Vec<f64> = vec![0.0; ordinal as usize];
    for r in &recs {
        score_of_ordinal[r.ordinal as usize] = r.score;
    }
    let mut rates: Vec<f64> = recs.iter().map(|r| r.error_rate).collect();
    drop(recs);

    {
        use std::os::unix::fs::FileExt;
        let f = match std::fs::File::create(&paths.sorted) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error: cannot write {}: {e}", paths.sorted.display());
                return Err(1);
            }
        };
        if let Err(e) = f.set_len(total) {
            eprintln!("Error: cannot write {}: {e}", paths.sorted.display());
            return Err(1);
        }
        let mut ordinal: u32 = 0;
        let mut failed: Option<std::io::Error> = None;
        if let Err(e) = fastq::for_each_file(Path::new(input), |r| {
            let this = ordinal as usize;
            ordinal += 1;
            let at = place[this];
            if at == u64::MAX || failed.is_some() {
                return;
            }
            let line = sorting::sorted_fastq_record(
                &r.name,
                score_of_ordinal[this],
                &r.seq,
                r.qual.as_deref().unwrap_or(""),
            );
            if let Err(e) = f.write_all_at(line.as_bytes(), at) {
                failed = Some(e);
            }
        }) {
            eprintln!("Error: cannot re-read {input}: {e}");
            return Err(1);
        }
        if let Some(e) = failed {
            eprintln!("Error: cannot write {}: {e}", paths.sorted.display());
            return Err(1);
        }
    }

    // The reference reports this through logging.DEBUG, so at the default level
    // it prints nothing. Not ported: it is unobservable unless --debug, and
    // --debug output is not in the contract.

    match sorting::logfile_contents(&mut rates) {
        Some(contents) => {
            if let Err(e) = std::fs::write(&paths.logfile, contents) {
                eprintln!("Error: cannot write {}: {e}", paths.logfile.display());
                return Err(1);
            }
        }
        None => {
            // The guard added to the Python on master (Finding 7). Both lines
            // go through logging.error, so both are shown at any level, and
            // both are in cli/q_filters_all.
            let q = pyfloat::repr(args.quality_threshold);
            let msg = format!("No reads passed the quality filter (--q {q}).\n");
            let _ = std::fs::write(&paths.logfile, &msg);
            eprintln!("Error: no reads passed the quality filter (--q {q}).");
            eprintln!("Lower --q, or check that the input has quality values.");
            return Err(1);
        }
    }
    Ok(nr_scored)
}

/// How many records a fastq holds. Only needed on the `--use_old_sorted_file`
/// path, where the reference re-reads the file it decided not to write.
fn count_records(p: &Path) -> usize {
    let mut n = 0usize;
    let _ = fastq::for_each_file(p, |_| n += 1);
    n
}

/// The one-line report for *Finding 18*, shared by the three places the
/// reference unpacks an eight-element tuple that may hold six.
///
/// One line per sentence instead of a `ValueError` traceback, and it names the
/// cause: `ValueError: not enough values to unpack (expected 8, got 6)` tells a
/// user nothing about `--use_old_sorted_file`, which is the only way to get here.
fn finding_18(c_id: usize) {
    eprintln!(
        "Error: cluster {c_id} was never processed by the clustering sweep, so it has no\n\
         homopolymer-compressed error rate. This happens when --use_old_sorted_file reuses a\n\
         sorted.fastq produced with a smaller --k. Re-sort with the --k you are clustering at,\n\
         or delete sorted.fastq. See PORTING.md, Finding 18."
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::Args;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ngsid_pipeline_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("tmpdir");
        d
    }

    fn args_for(dir: &Path, fastq_path: Option<&Path>) -> Args {
        Args {
            fastq: fastq_path.map(|p| p.to_string_lossy().into_owned()),
            outfolder: Some(dir.to_string_lossy().into_owned()),
            ..Default::default()
        }
    }

    // Both reads must survive the default filter, which is
    // `len(seq) >= 2*k` AND `len(homopolymer_compressed) >= k`, with k=13. The
    // second condition is the one that bites: TTTTGGGG... compresses to eight
    // bases and is dropped. These alternate enough to keep 32 after compression.
    const TWO_READS: &str =
        "@r1 x y\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
@r2 x y\nTGCATGCATGCATGCATGCATGCATGCATGCA\n+\nHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH\n";

    #[test]
    fn the_logfile_is_truncated_before_anything_else() {
        // Finding 13: --use_old_sorted_file leaves logfile.txt at zero bytes,
        // because it is opened for writing before the branch that skips the
        // work. The clustering is unchanged; the log is destroyed.
        let d = tmp("truncate");
        let paths = Paths::in_outfolder(&d.to_string_lossy());
        std::fs::write(&paths.logfile, b"previous run's statistics\n").unwrap();
        std::fs::write(&paths.sorted, TWO_READS).unwrap();
        let mut a = args_for(&d, None);
        a.use_old_sorted_file = true;
        let n = sort_stage(&a, &paths).expect("reuses the existing file");
        assert_eq!(n, 2);
        assert_eq!(
            std::fs::read(&paths.logfile).unwrap().len(),
            0,
            "logfile must be truncated to zero bytes"
        );
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn a_failed_use_old_run_leaves_an_empty_sorted_fastq() {
        // Finding 23. The first run fails; the second, identical, run then
        // takes the "use the existing sorted file" path and succeeds on zero
        // reads. A pipeline that retries after a failure sees exactly this.
        let d = tmp("poison");
        let paths = Paths::in_outfolder(&d.to_string_lossy());
        let mut a = args_for(&d, None);
        a.use_old_sorted_file = true;

        assert_eq!(sort_stage(&a, &paths), Err(1), "first run fails");
        assert!(paths.sorted.is_file(), "and leaves sorted.fastq behind");
        assert_eq!(std::fs::read(&paths.sorted).unwrap().len(), 0, "empty");

        let n = sort_stage(&a, &paths).expect("the retry 'succeeds'");
        assert_eq!(n, 0, "on zero reads");
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn accessions_keep_their_spaces_through_the_sort() {
        // The one behaviour that differs from isONclust, checked end to end
        // rather than only in fastq.rs: the header goes into sorted.fastq with
        // its spaces intact and the score appended after an underscore.
        let d = tmp("spaces");
        let input = d.join("in.fastq");
        std::fs::write(&input, TWO_READS).unwrap();
        let paths = Paths::in_outfolder(&d.to_string_lossy());
        let a = args_for(&d, Some(&input));
        let n = sort_stage(&a, &paths).expect("sorts");
        assert_eq!(n, 2);
        let out = std::fs::read_to_string(&paths.sorted).unwrap();
        let first = out.lines().next().unwrap();
        assert!(first.starts_with("@r"), "{first}");
        assert!(
            first.contains(" x y_"),
            "spaces kept, score appended: {first}"
        );
        std::fs::remove_dir_all(&d).ok();
    }

    /// Finding 27: the `--top_reads` guard is on the flag alone, and
    /// `--sample_size` defaults to 0, so `--top_reads` by itself truncates to
    /// nothing. The other branch treats a missing size as "everything". Both
    /// halves of that asymmetry are measured against the reference.
    #[test]
    fn top_reads_without_a_sample_size_keeps_nothing() {
        let d = tmp("topreads");
        let input = d.join("in.fastq");
        std::fs::write(&input, TWO_READS).unwrap();
        let paths = Paths::in_outfolder(&d.to_string_lossy());

        let mut a = args_for(&d, Some(&input));
        a.top_reads = true; // and sample_size stays 0
        sort_stage(&a, &paths).expect("sorts");
        // The truncation happens in cluster_stage; check the arithmetic it uses
        // rather than running the whole sweep.
        assert_eq!(a.sample_size, 0);
        assert_eq!(
            a.sample_size.max(0) as usize,
            0,
            "truncate(0) keeps nothing"
        );

        // ...while --sample_size alone at 0 falls through the
        // `0 < sample_size < len` guard and keeps everything.
        let mut b = args_for(&d, Some(&input));
        // Written as "the reference's guard, negated" rather than as the
        // simplified comparison clippy asks for: the point is that
        // `0 < sample_size` and `sample_size < len` are the two halves of
        // `if 0 < args.sample_size < len(read_array)`, and inverting them by
        // hand loses the correspondence with the line being reproduced.
        b.sample_size = 0;
        let triggers = b.sample_size > 0;
        assert!(!triggers, "0 does not trigger the subsample");
        b.sample_size = 999_999;
        let below_len = (b.sample_size as usize) < 2;
        assert!(
            !below_len,
            "a size above the read count does not trigger it either"
        );
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn filtering_everything_out_writes_the_guard_and_fails() {
        // Finding 7, as fixed on master: a line in the logfile, two on stderr,
        // exit 1.
        let d = tmp("qall");
        let input = d.join("in.fastq");
        std::fs::write(&input, TWO_READS).unwrap();
        let paths = Paths::in_outfolder(&d.to_string_lossy());
        let mut a = args_for(&d, Some(&input));
        a.quality_threshold = 99.0;
        assert_eq!(sort_stage(&a, &paths), Err(1));
        let log = std::fs::read_to_string(&paths.logfile).unwrap();
        assert_eq!(log, "No reads passed the quality filter (--q 99.0).\n");
        std::fs::remove_dir_all(&d).ok();
    }
}
