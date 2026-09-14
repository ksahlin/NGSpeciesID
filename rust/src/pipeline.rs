//! The driver: the stages of `NGSpeciesID`'s `main()`, in its order.
//!
//! The reference's `main` is four steps, and this mirrors them one for one:
//!
//! 1. sort every read by expected error-free k-mers (`get_sorted_fastq_for_cluster`)
//! 2. filter by length and subsample (`--m`/`--s`, `--sample_size`/`--top_reads`)
//! 3. load the empirical probability table and cluster
//! 4. write the output, and optionally form consensus sequences
//!
//! **Only step 1 exists so far.** The rest returns `Incomplete`, which `main`
//! turns into exit 70 — deliberately not one of the reference's own codes, so a
//! stage that is missing cannot be mistaken for a stage that agrees.
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
use crate::{blockalign, fastq, p_emp, packed, parallelize, pyfloat, sorting, sweep};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::{Path, PathBuf};

/// How far the pipeline got.
pub enum Outcome {
    /// Ran to completion. Nothing constructs this yet -- the last stage is not
    /// written -- but `main` already maps it to exit 0, so the day it is
    /// constructed nothing else has to change.
    #[allow(dead_code)]
    Done,
    /// A stage that is not written yet. `main` maps this to exit 70.
    Incomplete(&'static str),
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

    match cluster_stage(args, &paths, sorted_count) {
        Err(code) => Outcome::Failed(code),
        Ok(nontrivial) => {
            if args.consensus {
                return Outcome::Incomplete("consensus");
            }
            // The second and last line a default run prints.
            eprintln!("Finished Clustering: {nontrivial} clusters formed");
            Outcome::Done
        }
    }
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
/// Returns the number of clusters with more than one read, which is the number
/// the reference reports as "clusters formed".
fn cluster_stage(args: &Args, paths: &Paths, sorted_count: usize) -> Result<usize, i32> {
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

    // --top_reads: the highest-scoring `--sample_size` reads, which is the head
    // of a file already sorted by score descending. `--sample_size` WITHOUT
    // --top_reads is the seeded random draw and is not implemented yet.
    if args.top_reads {
        reads.truncate(args.sample_size.max(0) as usize);
        ordinal_of.truncate(reads.len());
    } else if args.sample_size > 0 && (args.sample_size as usize) < reads.len() {
        return Err(not_implemented_stage("--sample_size without --top_reads"));
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
        return write_output(args, paths, &reads, &ordinal_of, &res);
    }

    let clusters = sweep::OrderedClusters::default();
    let reps: FxHashMap<usize, sweep::ReadInfo> = FxHashMap::default();
    let db = crate::cluster::MinimizerDatabase::default();
    let res = sweep::reads_to_clusters(clusters, reps, &reads, db, 1, &table, params);

    write_output(args, paths, &reads, &ordinal_of, &res)
}

/// A stage that is not written yet, reported from inside `cluster_stage` where
/// the error type is an exit code. Keeps `EXIT_NOT_IMPLEMENTED` in one place.
fn not_implemented_stage(what: &str) -> i32 {
    eprintln!(
        "NGSpeciesID (Rust port): {what} is not implemented yet. \
         Arguments parsed and validated successfully."
    );
    crate::EXIT_NOT_IMPLEMENTED as i32
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
) -> Result<usize, i32> {
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
    for (output_cl_id, c_id) in order.iter().enumerate() {
        let rep = &representatives[c_id];
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
    }

    let outfolder = Path::new(args.outfolder.as_deref().expect("validated"));
    let cp = outfolder.join("final_clusters.tsv");
    let op = outfolder.join("final_cluster_origins.tsv");
    if let Err(e) = std::fs::write(&cp, clusters_out).and_then(|_| std::fs::write(&op, origins_out))
    {
        eprintln!("Error: cannot write output: {e}");
        return Err(1);
    }
    Ok(nontrivial)
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
