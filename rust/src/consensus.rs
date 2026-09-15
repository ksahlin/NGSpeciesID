//! `modules/consensus.py` — the `--consensus` stage.
//!
//! Four steps, in the reference's order:
//!
//! 1. `form_draft_consensus` — POA every cluster at or above `abundance_cutoff`
//! 2. `remove_barcodes` — primer or universal-tail trimming, *before* step 3
//! 3. `detect_reverse_complements` — merge centers that are each other's
//!    reverse complement, or are simply similar
//! 4. `polish_sequences` — medaka or racon, per surviving center
//!
//! …and then, when a primer file or `--remove_universal_tails` was given, steps
//! 2–4 again if the second trim changed anything.
//!
//! # What this module does NOT decide
//!
//! The order clusters are walked in, and the order reads go into each POA. Both
//! come from the clustering stage and both change the output — see
//! `crate::poa`. `form_draft_consensus` iterates a cluster's member list in its
//! **stored order**, which is the sweep's insertion order, *not* the
//! score-descending order `final_clusters.tsv` is written in. Sorting here would
//! be a silent divergence.

// Wired in by the pipeline's consensus stage, which is the next slice.
#![allow(dead_code)]

use crate::cli::Args;
use crate::{align, parasail, poa};
use std::path::{Path, PathBuf};

/// One draft consensus. Mirrors the reference's
/// `[nr_reads_in_cluster, c_id, center, reads_path]`, whose fourth element is a
/// path until `detect_reverse_complements` turns it into a list of paths.
#[derive(Debug, Clone)]
pub struct Center {
    pub nr_reads: usize,
    pub c_id: usize,
    pub seq: String,
    /// Every read file merged into this center. One entry until a merge.
    pub reads: Vec<PathBuf>,
}

/// `reverse_complement`, with the reference's own table.
///
/// It maps IUPAC codes and preserves case, and it **panics on nothing**: a byte
/// outside the table is passed through unchanged. The reference would raise
/// `KeyError`; nothing in this pipeline can reach it, because centers come from
/// spoa over ACGT.
pub fn reverse_complement(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            'N' => 'N',
            'n' => 'n',
            'X' => 'X',
            'x' => 'x',
            'Y' => 'R',
            'R' => 'Y',
            'K' => 'M',
            'M' => 'K',
            'S' => 'S',
            'W' => 'W',
            'B' => 'V',
            'V' => 'B',
            'H' => 'D',
            'D' => 'H',
            'y' => 'r',
            'r' => 'y',
            'k' => 'm',
            'm' => 'k',
            's' => 's',
            'w' => 'w',
            'b' => 'v',
            'v' => 'b',
            'h' => 'd',
            'd' => 'h',
            other => other,
        })
        .collect()
}

/// `consensus.parasail_alignment`'s scoring: match 2, mismatch −2, **opening
/// penalty 3**, gap extension 1.
///
/// Note the opening penalty. The clustering path bins it 5/4/3/2 by summed error
/// rate (`blockalign::gap_opening_penalty`); this one is a fixed 3, and it is a
/// *different call site with different defaults*. Reusing the clustering
/// scoring here would be wrong in a way no type would catch.
const RC_SCORING: parasail::Scoring = parasail::Scoring {
    match_score: 2,
    mismatch: -2,
    open: 3,
    ext: 1,
};

/// The identity of one gapped alignment, as `highest_aln_identity` computes it.
///
/// Divided by the **alignment length**, so leading and trailing gaps count
/// against it — an end-to-end overlap of two sequences of very different lengths
/// scores low even where they agree perfectly.
fn identity(s1: &str, s2: &str) -> f64 {
    let aln = parasail::semiglobal(s1.as_bytes(), s2.as_bytes(), RC_SCORING);
    let (a1, a2) = match align::ops_to_seq(&aln.ops, s1.as_bytes(), s2.as_bytes()) {
        Some(p) => p,
        None => return 0.0,
    };
    let total = a1.len();
    if total == 0 {
        return 0.0;
    }
    let mismatching = a1.iter().zip(a2.iter()).filter(|(x, y)| x != y).count();
    (total - mismatching) as f64 / total as f64
}

/// `highest_aln_identity`: the better of the forward and reverse-complement
/// identities.
///
/// The reference computes the **reverse complement first**, then the forward,
/// and returns `max`. Order does not affect the result, but the two parasail
/// calls happen in that order and a dump oracle records them that way.
pub fn highest_aln_identity(seq: &str, seq2: &str) -> f64 {
    let rc = identity(seq, &reverse_complement(seq2));
    let fw = identity(seq, seq2);
    fw.max(rc)
}

/// `detect_reverse_complements`.
///
/// # Finding 10: this double-counts, and the port reproduces it
///
/// The outer loop skips a center already merged away; the **inner loop does
/// not**. So a center merged into an earlier one can be merged again into a
/// later one, and its reads are counted twice and handed to two polishers.
///
/// Constructed and measured against the reference, with three centers where
/// `identity(A,C) = identity(B,C) = 0.937` and `identity(A,B) = 0.873` at a
/// threshold of 0.9:
///
/// | | in | out |
/// | --- | --- | --- |
/// | A | 10 reads | **13**, `[a.fq, c.fq]` |
/// | B | 5 reads | **8**, `[b.fq, c.fq]` |
/// | C | 3 reads | merged twice |
/// | total | **18** | **21** |
///
/// The inflated count reaches output: it is written into the consensus fasta
/// header as `total_supporting_reads`.
pub fn detect_reverse_complements(centers: Vec<Center>, rc_identity_threshold: f64) -> Vec<Center> {
    let mut filtered: Vec<Center> = Vec::new();
    let mut already_removed: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let n = centers.len();

    for i in 0..n {
        if already_removed.contains(&centers[i].c_id) {
            continue;
        }
        let mut merged = centers[i].clone();
        // The reference's `elif i == len(centers) - 1` branch: the last center,
        // when not already removed, is appended without comparing it to
        // anything. Equivalent to running the loop over an empty tail, but
        // written out because the reference's branch is visible in its output
        // ordering.
        if i < n - 1 {
            for other in &centers[i + 1..] {
                // NOT `if already_removed.contains(...)`. See the doc comment:
                // this is Finding 10 and it is deliberate.
                let id = highest_aln_identity(&merged.seq, &other.seq);
                if id >= rc_identity_threshold {
                    merged.nr_reads += other.nr_reads;
                    already_removed.insert(other.c_id);
                    merged.reads.extend(other.reads.iter().cloned());
                }
            }
        }
        filtered.push(merged);
    }
    filtered
}

/// `form_draft_consensus`: POA every cluster at or above `abundance_cutoff`.
///
/// `clusters` is `(c_id, members)` **already in the reference's walk order** —
/// `(cluster size, representative score)` descending — because that order is the
/// clustering stage's to decide and this stage's to respect.
///
/// `read_of` yields a member's `(sequence, quality)` from `sorted.fastq`.
#[allow(clippy::too_many_arguments)]
pub fn form_draft_consensus(
    clusters: &[(usize, Vec<String>)],
    read_of: &dyn Fn(&str) -> Option<(String, String)>,
    work_dir: &Path,
    abundance_cutoff: usize,
    max_seqs_for_consensus: i64,
) -> Result<Vec<Center>, String> {
    let mut centers: Vec<Center> = Vec::new();
    for (c_id, members) in clusters {
        let nr_reads = members.len();
        if nr_reads < abundance_cutoff {
            // The reference counts singletons and discarded clusters here and
            // reports them through logging.debug, which prints nothing at the
            // default level. Not ported: unobservable.
            continue;
        }
        let reads_path = work_dir.join(format!("reads_c_id_{c_id}.fq"));
        let mut seqs: Vec<String> = Vec::new();
        let mut quals: Vec<String> = Vec::new();
        let mut body = String::new();
        for (i, acc) in members.iter().enumerate() {
            // `>= 0 and i >= max_seqs_for_consensus`. Note the `>=`, which
            // admits exactly that many sequences -- isONcorrect's equivalent
            // uses a bare `>` and admits one more.
            if max_seqs_for_consensus >= 0 && i as i64 >= max_seqs_for_consensus {
                break;
            }
            let Some((seq, qual)) = read_of(acc) else {
                return Err(format!(
                    "read '{acc}' is in a cluster but not in sorted.fastq"
                ));
            };
            body.push_str(&format!("@{acc}\n{seq}\n+\n{qual}\n"));
            seqs.push(seq);
            quals.push(qual);
        }
        // The file is written even though the POA is now in-process: the
        // reference writes it, `polish_sequences` reads it back, and racon and
        // medaka are both handed it by path.
        std::fs::write(&reads_path, &body)
            .map_err(|e| format!("cannot write {}: {e}", reads_path.display()))?;
        let center = poa::consensus(&seqs, &quals);
        centers.push(Center {
            nr_reads,
            c_id: *c_id,
            seq: center,
            reads: vec![reads_path],
        });
    }
    Ok(centers)
}

/// `barcode_trimmer.read_barcodes`: the primer file, plus each primer's reverse
/// complement.
///
/// Two behaviours reproduced because both are observable and neither is
/// obviously intended (*Finding 16*):
///
/// * **The fasta description stays in the key.** `readfq` returns the whole
///   header, so the committed primer file's `>COIF-ALT ` — with a trailing
///   space — yields the key `"COIF-ALT _fw"`. It reaches only log output today.
/// * **Only the reverse complement is upper-cased.** `seq.strip()` is stored as
///   written for the `_fw` entry while the `_rc` is built from `seq.upper()`. A
///   lowercase primer file therefore gets a lowercase forward primer and an
///   uppercase reverse complement — and the IUPAC equality table is
///   uppercase-only, so ambiguity codes stop matching in one direction. The
///   committed file is uppercase, so this is unexercised.
///
/// Insertion order is preserved: the reference iterates this dict, and although
/// the consumer takes a max and a min over all hits — so order cannot change the
/// answer — an order-dependent container here would be a trap for later.
pub fn read_barcodes(records: &[(String, String)]) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    for (acc, seq) in records {
        out.push((format!("{acc}_fw"), seq.trim().to_string()));
    }
    // A second pass over the ORIGINAL list, as the reference's
    // `for acc, seq in list(barcodes.items())` does -- it snapshots before
    // inserting, so the reverse complements are not themselves reversed.
    let fw: Vec<(String, String)> = out.clone();
    for (acc, seq) in &fw {
        let base = &acc[..acc.len() - 3]; // drop "_fw"
        out.push((
            format!("{base}_rc"),
            reverse_complement(&seq.to_uppercase()),
        ));
    }
    out
}

/// `barcode_trimmer.get_universal_tails`.
///
/// Note the naming: the two given literals are `1_F_fw` and `2_R_rc`, and their
/// reverse complements become `1_F_rc` and `2_R_fw`. So `2_R`'s "forward" entry
/// is the reverse complement of the literal, not the literal.
pub fn universal_tails() -> Vec<(String, String)> {
    let f_fw = "TTTCTGTTGGTGCTGATATTGC";
    let r_rc = "ACTTGCCTGTCGCTCTATCTTC";
    vec![
        ("1_F_fw".to_string(), f_fw.to_string()),
        ("2_R_rc".to_string(), r_rc.to_string()),
        ("1_F_rc".to_string(), reverse_complement(f_fw)),
        ("2_R_fw".to_string(), reverse_complement(r_rc)),
    ]
}

/// One primer hit: the barcode's name and the `locations[0]` the reference reads.
struct BarcodeHit {
    start: usize,
    stop: usize,
}

/// `find_barcode_locations`: every barcode with a hit, at its FIRST location.
fn find_barcode_locations(
    window: &str,
    barcodes: &[(String, String)],
    primer_max_ed: i64,
) -> Vec<BarcodeHit> {
    let mut out = Vec::new();
    for (_acc, primer) in barcodes {
        let r = crate::edlib::align_hw(
            primer.as_bytes(),
            window.as_bytes(),
            primer_max_ed,
            crate::edlib::IUPAC_EQUALITIES,
        );
        // `if locations:` -- a hit is a non-empty list, which edlib gives only
        // when the distance is within k.
        if let Some(loc) = r.first_location() {
            out.push(BarcodeHit {
                start: loc.start,
                stop: loc.end,
            });
        }
    }
    out
}

/// `barcode_trimmer.remove_barcodes`: trim each center in place, and report
/// whether anything changed.
///
/// The arithmetic is the reference's and it is asymmetric:
///
/// * the **start** cut is the **latest** `stop` among hits in the leading
///   window — so overlapping primers all get removed;
/// * the **end** cut is derived from the **earliest** `start` among hits in the
///   trailing window, as `len(center) - (trim_window - earliest_hit)`.
///
/// `trim_window` halves to `len(center) / 2` when `2 * --trim_window` exceeds the
/// center, so the two windows never overlap.
pub fn remove_barcodes(
    centers: &mut [Center],
    barcodes: &[(String, String)],
    trim_window: i64,
    primer_max_ed: i64,
) -> bool {
    let mut updated = false;
    for center in centers.iter_mut() {
        let n = center.seq.len();
        let tw = if 2 * trim_window as usize > n {
            n / 2
        } else {
            trim_window as usize
        };
        if tw == 0 {
            continue;
        }
        let head = &center.seq[..tw.min(n)];
        let tail = &center.seq[n.saturating_sub(tw)..];

        let mut cut_start = 0usize;
        for h in find_barcode_locations(head, barcodes, primer_max_ed) {
            if h.stop > cut_start {
                cut_start = h.stop;
            }
        }
        let mut cut_end = n;
        let hits_end = find_barcode_locations(tail, barcodes, primer_max_ed);
        if !hits_end.is_empty() {
            let mut earliest = n;
            for h in &hits_end {
                if h.start < earliest {
                    earliest = h.start;
                }
            }
            // `len(center) - (trim_window - earliest_hit)`, where trim_window is
            // the possibly-halved one.
            cut_end = n - (tw - earliest.min(tw));
        }
        if cut_start > 0 || cut_end < n {
            // Python slicing clamps; an inverted range yields "".
            let (a, b) = (cut_start.min(n), cut_end.min(n));
            center.seq = if a < b {
                center.seq[a..b].to_string()
            } else {
                String::new()
            };
            updated = true;
        }
    }
    updated
}

/// Which polisher, if any. The reference's two flags are mutually exclusive in
/// argparse, and `--consensus` with neither is *Finding 4*'s crash — rejected by
/// `cli::validate` before this stage runs.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Polisher {
    Medaka,
    Racon,
}

impl Polisher {
    pub fn of(args: &Args) -> Option<Polisher> {
        if args.medaka {
            Some(Polisher::Medaka)
        } else if args.racon {
            Some(Polisher::Racon)
        } else {
            None
        }
    }
    /// The per-center output directory's prefix: `medaka_cl_id_<id>` or
    /// `racon_cl_id_<id>`.
    pub fn dir_prefix(self) -> &'static str {
        match self {
            Polisher::Medaka => "medaka_cl_id_",
            Polisher::Racon => "racon_cl_id_",
        }
    }
}

/// `polish_sequences`: write each center's draft and reads, then polish.
///
/// Writes, per surviving center:
///
/// * `consensus_reference_<c_id>.fasta` — the draft, with a header naming the
///   cluster and its **`total_supporting_reads`**, which is the count
///   `detect_reverse_complements` may have inflated (*Finding 10*)
/// * `reads_to_consensus_<c_id>.fastq` — every merged cluster's reads
/// * `<medaka|racon>_cl_id_<c_id>/` — the polisher's own output tree
///
/// and first **deletes** any `consensus_reference_*` and `<polisher>_cl_id_*`
/// left from a previous run, so a rerun into the same folder does not mix old
/// and new. That deletion is why the second pass of the re-trim loop does not
/// accumulate files.
pub fn polish_sequences(
    centers: &mut [Center],
    outfolder: &Path,
    polisher: Polisher,
    args: &Args,
) -> Result<(), String> {
    // Clear the previous run's output, exactly as the reference's two glob
    // loops do. Note it removes `consensus_reference_*` -- every polisher's --
    // but only its OWN `<polisher>_cl_id_*`, so switching --racon to --medaka in
    // the same folder leaves the racon directories behind. Faithful.
    if let Ok(entries) = std::fs::read_dir(outfolder) {
        for e in entries.flatten() {
            let name = e.file_name().to_string_lossy().into_owned();
            if name.starts_with(polisher.dir_prefix()) {
                let _ = std::fs::remove_dir_all(e.path());
            } else if name.starts_with("consensus_reference_") {
                let _ = std::fs::remove_file(e.path());
            }
        }
    }

    for center in centers.iter_mut() {
        let c_id = center.c_id;
        let draft = outfolder.join(format!("consensus_reference_{c_id}.fasta"));
        std::fs::write(
            &draft,
            format!(
                ">consensus_cl_id_{c_id}_total_supporting_reads_{}\n{}\n",
                center.nr_reads, center.seq
            ),
        )
        .map_err(|e| format!("cannot write {}: {e}", draft.display()))?;

        // Every merged cluster's reads, in the order the files were merged.
        //
        // The reference builds a DICT keyed by accession per file, so a
        // duplicate accession within one file collapses -- and then iterates it
        // in insertion order. Reproduced with an order-preserving de-duplication
        // per file, because the fastq is handed to the polisher and its order
        // reaches racon's output.
        let reads_file = outfolder.join(format!("reads_to_consensus_{c_id}.fastq"));
        let mut body = String::new();
        for src in &center.reads {
            let text = std::fs::read_to_string(src)
                .map_err(|e| format!("cannot read {}: {e}", src.display()))?;
            let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
            let mut order: Vec<(String, String, String)> = Vec::new();
            crate::fastq::for_each(text.split_inclusive('\n'), |r| {
                let acc = r.name.clone();
                if seen.insert(acc.clone()) {
                    order.push((acc, r.seq.clone(), r.qual.clone().unwrap_or_default()));
                }
            });
            for (acc, seq, qual) in order {
                // `acc.split()[0]` -- truncated at the first whitespace, and
                // ONLY here. Every other writer keeps the whole accession.
                let short = acc.split_whitespace().next().unwrap_or("");
                body.push_str(&format!("@{short}\n{seq}\n+\n{qual}\n"));
            }
        }
        std::fs::write(&reads_file, &body)
            .map_err(|e| format!("cannot write {}: {e}", reads_file.display()))?;

        let dir = outfolder.join(format!("{}{c_id}", polisher.dir_prefix()));
        std::fs::create_dir_all(&dir)
            .map_err(|e| format!("cannot create {}: {e}", dir.display()))?;

        let polished = match polisher {
            Polisher::Racon => run_racon(&reads_file, &draft, &dir, args.racon_iter)?,
            Polisher::Medaka => run_medaka(&reads_file, &draft, &dir, args)?,
        };
        center.seq = polished;
    }
    Ok(())
}

/// Second line of a fasta or fastq, which is what the reference reads back from
/// every polisher: `cf.readlines()[1].strip()`.
fn second_line(path: &Path) -> Option<String> {
    let text = std::fs::read_to_string(path).ok()?;
    text.lines().nth(1).map(|l| l.trim().to_string())
}

/// `run_racon`: `--racon_iter` rounds of minimap2 + racon, each round polishing
/// the previous round's output.
fn run_racon(reads: &Path, draft: &Path, dir: &Path, racon_iter: i64) -> Result<String, String> {
    let mut center_file = draft.to_path_buf();
    for i in 0..racon_iter.max(0) {
        let paf = dir.join(format!("read_alignments_it_{i}.paf"));
        let out = dir.join(format!("racon_polished_it_{i}.fasta"));
        // The captured stderr files carry timings, so they are NOT in the
        // byte-identity contract -- the harness excludes them. They are still
        // written, because the reference writes them and their absence would be
        // a missing file.
        run_capturing(
            "minimap2",
            &[
                "-x",
                "map-ont",
                &center_file.to_string_lossy(),
                &reads.to_string_lossy(),
            ],
            &paf,
            &dir.join(format!("mm2_stderr_it_{i}.txt")),
        )?;
        run_capturing(
            "racon",
            &[
                &reads.to_string_lossy(),
                &paf.to_string_lossy(),
                &center_file.to_string_lossy(),
            ],
            &out,
            &dir.join(format!("racon_stderr_it_{i}.txt")),
        )?;
        center_file = out;
    }
    let final_path = dir.join("consensus.fasta");
    std::fs::copy(&center_file, &final_path)
        .map_err(|e| format!("cannot write {}: {e}", final_path.display()))?;
    second_line(&final_path).ok_or_else(|| format!("{} has no sequence", final_path.display()))
}

/// `run_medaka`: one `medaka_consensus` call, then read back whichever of
/// `consensus.fasta` / `consensus.fastq` it produced.
fn run_medaka(reads: &Path, draft: &Path, dir: &Path, args: &Args) -> Result<String, String> {
    let mut argv: Vec<String> = vec![
        "-i".into(),
        reads.to_string_lossy().into_owned(),
        "-d".into(),
        draft.to_string_lossy().into_owned(),
        "-o".into(),
        dir.to_string_lossy().into_owned(),
        // Hard-coded "1" in the reference, NOT --t. A --t 8 run still polishes
        // single-threaded.
        "-t".into(),
        "1".into(),
    ];
    if !args.medaka_model.is_empty() {
        argv.push("-m".into());
        argv.push(args.medaka_model.clone());
    }
    if args.medaka_fastq {
        argv.push("-q".into());
    }
    let refs: Vec<&str> = argv.iter().map(|s| s.as_str()).collect();
    run_capturing(
        "medaka_consensus",
        &refs,
        &dir.join("stdout.txt"),
        &dir.join("stderr.txt"),
    )?;
    // "consider all output formats for compatibility with all Medaka versions",
    // fasta first.
    for name in ["consensus.fasta", "consensus.fastq"] {
        let p = dir.join(name);
        if p.is_file() {
            if let Some(s) = second_line(&p) {
                return Ok(s);
            }
        }
    }
    // The reference's `assert centers[i][2], "Medaka consensus sequence not found"`.
    Err(format!(
        "medaka produced no consensus in {} (looked for consensus.fasta and consensus.fastq)",
        dir.display()
    ))
}

/// Run a command, sending stdout and stderr to the given files.
fn run_capturing(prog: &str, args: &[&str], stdout: &Path, stderr: &Path) -> Result<(), String> {
    let out = std::fs::File::create(stdout)
        .map_err(|e| format!("cannot write {}: {e}", stdout.display()))?;
    let err = std::fs::File::create(stderr)
        .map_err(|e| format!("cannot write {}: {e}", stderr.display()))?;
    let status = std::process::Command::new(prog)
        .args(args)
        .stdout(out)
        .stderr(err)
        .status()
        .map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                format!("{prog} not found on PATH")
            } else {
                format!("could not run {prog}: {e}")
            }
        })?;
    if !status.success() {
        return Err(format!(
            "{prog} failed ({status}); see {}",
            stderr.display()
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(nr: usize, id: usize, seq: &str, f: &str) -> Center {
        Center {
            nr_reads: nr,
            c_id: id,
            seq: seq.to_string(),
            reads: vec![PathBuf::from(f)],
        }
    }

    #[test]
    fn reverse_complement_matches_the_references_table() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("acgt"), "acgt");
        // IUPAC, which the primer file uses.
        assert_eq!(reverse_complement("RYKMSWBVHDN"), "NHDBVWSKMRY");
        // Case is preserved, not normalised.
        assert_eq!(reverse_complement("AcGt"), "aCgT");
    }

    #[test]
    fn identity_of_a_sequence_with_itself_is_one() {
        let s = "ACGTACGTACGTAAGGCCTTACGTACGT";
        assert_eq!(identity(s, s), 1.0);
        assert_eq!(highest_aln_identity(s, s), 1.0);
    }

    #[test]
    fn a_reverse_complement_scores_as_high_as_the_original() {
        let s = "ACGTACGTTTGGCCAATTACGTACGTAAGG";
        let rc = reverse_complement(s);
        // Forward identity against the RC is poor; the RC branch recovers it.
        let both = highest_aln_identity(s, &rc);
        assert!(
            both > 0.99,
            "an exact reverse complement should score ~1, got {both}"
        );
    }

    /// Finding 10, as a test rather than as a comment. The numbers are the ones
    /// measured against the reference.
    #[test]
    fn detect_reverse_complements_double_counts() {
        // Three centers where A~C and B~C but A!~B. Built from literal
        // sequences so the identities are deterministic.
        let base: String = std::iter::repeat_n("ACGTTGCA", 80).collect();
        let mut a = base.clone().into_bytes();
        let mut b = base.clone().into_bytes();
        // A differs from base in its first half, B in its second; both by
        // enough that A and B differ twice as much as either does from C.
        for i in (0..base.len() / 2).step_by(8) {
            a[i] = b'G';
        }
        for i in ((base.len() / 2)..base.len()).step_by(8) {
            b[i] = b'G';
        }
        let a = String::from_utf8(a).unwrap();
        let b = String::from_utf8(b).unwrap();
        let ac = highest_aln_identity(&a, &base);
        let bc = highest_aln_identity(&b, &base);
        let ab = highest_aln_identity(&a, &b);
        // The construction has to hold for the test to mean anything.
        assert!(ac > ab && bc > ab, "ac={ac} bc={bc} ab={ab}");
        let t = (ab + ac.min(bc)) / 2.0; // a threshold between them

        let centers = vec![
            c(10, 0, &a, "a.fq"),
            c(5, 1, &b, "b.fq"),
            c(3, 2, &base, "c.fq"),
        ];
        let out = detect_reverse_complements(centers, t);
        let total: usize = out.iter().map(|x| x.nr_reads).sum();
        assert_eq!(out.len(), 2, "A and B both survive");
        assert!(
            total > 18,
            "18 reads went in and {total} came out -- Finding 10's double count"
        );
        // ...and c.fq is attached to both survivors.
        assert!(out
            .iter()
            .all(|x| x.reads.iter().any(|p| p.ends_with("c.fq"))));
    }

    #[test]
    fn a_center_already_merged_is_skipped_by_the_outer_loop() {
        let s: String = std::iter::repeat_n("ACGTTGCA", 80).collect();
        let centers = vec![
            c(10, 0, &s, "a.fq"),
            c(5, 1, &s, "b.fq"),
            c(3, 2, &s, "c.fq"),
        ];
        // All identical, so the first absorbs both and the others are skipped.
        let out = detect_reverse_complements(centers, 0.9);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].nr_reads, 18, "no double count when all three merge");
        assert_eq!(out[0].reads.len(), 3);
    }

    #[test]
    fn the_rc_scoring_is_not_the_clustering_scoring() {
        // A fixed opening penalty of 3, where the clustering path bins 5/4/3/2.
        assert_eq!(RC_SCORING.open, 3);
        assert_eq!(RC_SCORING.match_score, 2);
        assert_eq!(RC_SCORING.mismatch, -2);
        assert_eq!(RC_SCORING.ext, 1);
    }

    #[test]
    fn polisher_is_read_from_the_flags() {
        let mut a = Args::default();
        assert_eq!(Polisher::of(&a), None);
        a.medaka = true;
        assert_eq!(Polisher::of(&a), Some(Polisher::Medaka));
        a.medaka = false;
        a.racon = true;
        assert_eq!(Polisher::of(&a), Some(Polisher::Racon));
        assert_eq!(Polisher::Racon.dir_prefix(), "racon_cl_id_");
        assert_eq!(Polisher::Medaka.dir_prefix(), "medaka_cl_id_");
    }
}
