//! The `write_fastq` subcommand: split a clustering into per-cluster fastq files.
//!
//! ```python
//! clusters = defaultdict(list)
//! with open(args.clusters) as f:
//!     for line in f:
//!         items = line.strip().split()
//!         cl_id, acc = items[0], items[1]
//!         clusters[cl_id].append(acc)
//! reads = { acc : (seq, qual) for acc, (seq, qual) in readfq(open(args.fastq)) }
//! for cl_id in clusters:
//!     r = clusters[cl_id]
//!     if len(r) >= args.N:
//!         curr_file = open(os.path.join(args.outfolder, str(cl_id) + ".fastq"), "w")
//!         for acc in r:
//!             seq, qual = reads[acc]
//!             ...
//! ```
//!
//! # This is broken, and the port reproduces the breakage exactly
//!
//! `line.strip().split()` splits on **all** whitespace, and column 2 of
//! `final_clusters.tsv` is a full ONT accession *containing spaces*:
//!
//! ```text
//! 0	c948601e-1bd9-4039-94c5-3d8c741df65f runid=ed1de13 read=4709 ch=45 …
//! ```
//!
//! So `acc` becomes the bare UUID while `reads` is keyed by the whole header,
//! and the very first lookup raises `KeyError`. Measured on both corpora: the
//! run creates `0.fastq`, writes **nothing** to it, and exits 1. That zero-byte
//! file is in the goldens.
//!
//! It only fails on a fastq whose headers contain spaces — which is every
//! ONT-basecalled fastq, and therefore NGSpeciesID's entire target input.
//! isONclust never sees it because its `readfq` substitutes spaces for
//! underscores. PORTING.md, *Finding 6*, and the one-line fix
//! (`split("\t", 1)`) is in *Deferred improvements*.
//!
//! # Three details that are contract rather than accident
//!
//! * **Cluster ids stay strings.** The reference uses `str(cl_id)` unparsed as
//!   the filename, so a file holding `007` would produce `007.fastq`. Parsing to
//!   an integer would quietly rename it.
//! * **Output order is the order ids first appear** in the clusters file,
//!   because the reference iterates a `defaultdict`. It decides only which file
//!   is created first — which is exactly what decides which one exists when the
//!   `KeyError` lands.
//! * **Accessions are looked up in the ORIGINAL `--fastq`**, not `sorted.fastq`,
//!   so they carry no score suffix.

use crate::cli::WriteFastqArgs;
use crate::fastq;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Run the subcommand. Returns the process exit code.
pub fn run(args: &WriteFastqArgs) -> i32 {
    let Some(clusters_path) = args.clusters.as_deref() else {
        eprintln!("Error: write_fastq needs --clusters.");
        return 1;
    };
    let Some(fastq_path) = args.fastq.as_deref() else {
        eprintln!("Error: write_fastq needs --fastq.");
        return 1;
    };
    let Some(outfolder) = args.outfolder.as_deref() else {
        // The reference reaches os.path.join(None, ...) and raises TypeError,
        // the same shape as Finding 9.
        eprintln!("Error: write_fastq needs --outfolder.");
        return 1;
    };

    let text = match std::fs::read_to_string(clusters_path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error: cannot read {clusters_path}: {e}");
            return 1;
        }
    };

    // An order-preserving multimap, because the reference iterates a
    // defaultdict and that order decides which file is created first.
    let mut order: Vec<String> = Vec::new();
    let mut members: HashMap<String, Vec<String>> = HashMap::new();
    for line in text.lines() {
        // `line.strip().split()` -- ALL whitespace, which is the bug. Taking
        // items[0] and items[1] means an accession containing spaces is
        // truncated at the first one.
        let mut it = line.split_whitespace();
        let (Some(cl_id), Some(acc)) = (it.next(), it.next()) else {
            // `items[0], items[1]` on a short line raises IndexError. A blank
            // trailing line is the common case and the reference dies on it
            // too -- but `final_clusters.tsv` has no trailing blank line, so
            // this is unreachable from the tool's own output.
            if line.trim().is_empty() {
                continue;
            }
            eprintln!("Error: {clusters_path} has a line with fewer than two columns.");
            return 1;
        };
        let e = members.entry(cl_id.to_string()).or_insert_with(|| {
            order.push(cl_id.to_string());
            Vec::new()
        });
        e.push(acc.to_string());
    }

    // Keyed by the FULL accession, which is what makes the truncated lookup
    // above fail. A dict, so a duplicate accession keeps the last record.
    let mut reads: HashMap<String, (String, String)> = HashMap::new();
    if let Err(e) = fastq::for_each_file(Path::new(fastq_path), |r| {
        reads.insert(r.name, (r.seq, r.qual.unwrap_or_default()));
    }) {
        eprintln!("Error: cannot read {fastq_path}: {e}");
        return 1;
    }

    if let Err(e) = std::fs::create_dir_all(outfolder) {
        eprintln!("Error: cannot create {outfolder}: {e}");
        return 1;
    }

    for cl_id in &order {
        let r = &members[cl_id];
        if (r.len() as i64) < args.n {
            continue;
        }
        // Created BEFORE the reads are looked up, which is why a failure leaves
        // a zero-byte file behind. The goldens record exactly that.
        let path = Path::new(outfolder).join(format!("{cl_id}.fastq"));
        let mut f = match std::fs::File::create(&path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error: cannot write {}: {e}", path.display());
                return 1;
            }
        };
        for acc in r {
            let Some((seq, qual)) = reads.get(acc) else {
                // The reference's KeyError. One line instead of a traceback,
                // with the cause named -- this is reachable on any ONT fastq
                // and "KeyError: '<uuid>'" explains nothing.
                eprintln!("Error: read '{acc}' is in {clusters_path} but not in {fastq_path}.");
                eprintln!(
                    "The cluster file's accessions are split on whitespace, so a header \
                     containing spaces is truncated. See PORTING.md, Finding 6."
                );
                return 1;
            };
            if let Err(e) = write!(f, "@{acc}\n{seq}\n+\n{qual}\n") {
                eprintln!("Error: cannot write {}: {e}", path.display());
                return 1;
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("ngsid_wf_{name}_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("tmpdir");
        d
    }

    fn args(dir: &Path, clusters: &Path, fq: &Path, n: i64) -> WriteFastqArgs {
        WriteFastqArgs {
            clusters: Some(clusters.to_string_lossy().into_owned()),
            fastq: Some(fq.to_string_lossy().into_owned()),
            outfolder: Some(dir.join("out").to_string_lossy().into_owned()),
            n,
        }
    }

    /// Space-free accessions work, which is the case isONclust always has and
    /// this reference almost never does.
    #[test]
    fn space_free_accessions_round_trip() {
        let d = tmp("plain");
        let fq = d.join("in.fastq");
        std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n").unwrap();
        let cl = d.join("clusters.tsv");
        std::fs::write(&cl, "0\tr1\n0\tr2\n").unwrap();
        let a = args(&d, &cl, &fq, 0);
        assert_eq!(run(&a), 0);
        let out = std::fs::read_to_string(d.join("out/0.fastq")).unwrap();
        assert_eq!(out, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n");
        std::fs::remove_dir_all(&d).ok();
    }

    /// Finding 6: an accession with spaces is truncated at the first one, the
    /// lookup fails, and the file is left at zero bytes. Both halves matter --
    /// the exit code AND the empty file, which is what the goldens record.
    #[test]
    fn an_accession_with_spaces_fails_leaving_an_empty_file() {
        let d = tmp("spaces");
        let fq = d.join("in.fastq");
        std::fs::write(&fq, "@r1 runid=x ch=9\nACGT\n+\nIIII\n").unwrap();
        let cl = d.join("clusters.tsv");
        std::fs::write(&cl, "0\tr1 runid=x ch=9\n").unwrap();
        let a = args(&d, &cl, &fq, 0);
        assert_eq!(run(&a), 1, "the reference exits 1 here");
        let out = d.join("out/0.fastq");
        assert!(out.is_file(), "the file is created before the lookup");
        assert_eq!(std::fs::read(&out).unwrap().len(), 0, "and left empty");
        std::fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn n_filters_by_cluster_size() {
        let d = tmp("nfilter");
        let fq = d.join("in.fastq");
        std::fs::write(
            &fq,
            "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n",
        )
        .unwrap();
        let cl = d.join("clusters.tsv");
        std::fs::write(&cl, "0\tr1\n0\tr2\n1\tr3\n").unwrap();

        // --N 2 keeps cluster 0 (two reads) and drops cluster 1 (one).
        let a = args(&d, &cl, &fq, 2);
        assert_eq!(run(&a), 0);
        assert!(d.join("out/0.fastq").is_file());
        assert!(!d.join("out/1.fastq").is_file(), "cluster 1 is too small");
        std::fs::remove_dir_all(&d).ok();
    }

    /// Cluster ids are strings, not integers. `007` names `007.fastq`.
    #[test]
    fn cluster_ids_are_not_parsed() {
        let d = tmp("strid");
        let fq = d.join("in.fastq");
        std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n").unwrap();
        let cl = d.join("clusters.tsv");
        std::fs::write(&cl, "007\tr1\n").unwrap();
        let a = args(&d, &cl, &fq, 0);
        assert_eq!(run(&a), 0);
        assert!(d.join("out/007.fastq").is_file(), "not 7.fastq");
        std::fs::remove_dir_all(&d).ok();
    }

    /// Files are created in the order cluster ids FIRST APPEAR, not sorted.
    /// It decides which file exists when a KeyError lands.
    #[test]
    fn output_order_is_first_appearance() {
        let d = tmp("order");
        let fq = d.join("in.fastq");
        std::fs::write(&fq, "@r1\nACGT\n+\nIIII\n@r2 x\nTTTT\n+\nJJJJ\n").unwrap();
        let cl = d.join("clusters.tsv");
        // Cluster 5 appears first; cluster 1's read has a space and will fail.
        std::fs::write(&cl, "5\tr1\n1\tr2 x\n").unwrap();
        let a = args(&d, &cl, &fq, 0);
        assert_eq!(run(&a), 1);
        assert!(
            d.join("out/5.fastq").is_file(),
            "5 was written before 1 failed"
        );
        // "@r1\n" + "ACGT\n" + "+\n" + "IIII\n" = 4 + 5 + 2 + 5.
        assert_eq!(
            std::fs::read(d.join("out/5.fastq")).unwrap().len(),
            16,
            "and written in full"
        );
        assert!(d.join("out/1.fastq").is_file());
        assert_eq!(std::fs::read(d.join("out/1.fastq")).unwrap().len(), 0);
        std::fs::remove_dir_all(&d).ok();
    }
}
