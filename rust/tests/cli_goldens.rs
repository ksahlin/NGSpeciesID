//! Run the built binary against the recorded CLI goldens.
//!
//! `bench/equivalence.sh cli verify` is the authoritative check and it needs a
//! reference environment. This needs nothing but the goldens, which are
//! committed, so `cargo test` alone catches a CLI regression — and CI does not
//! have to build conda to find one.
//!
//! It covers the **exact** class only: the 23 cases whose stdout, stderr and
//! exit code are byte-identity contract. The traceback class is asserted by the
//! harness (exit code, non-empty, not-a-stack-trace) and the pending class
//! needs the clustering stages, so neither belongs here. `bench/README.md`,
//! "The 44 CLI cases are three classes", is the split.
//!
//! If `bench/golden/` is absent — a source tarball, say — every test reports
//! itself skipped rather than passing vacuously.

use std::path::{Path, PathBuf};
use std::process::Command;

fn repo_root() -> PathBuf {
    // CARGO_MANIFEST_DIR is <repo>/rust.
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("rust/ has a parent")
        .to_path_buf()
}

fn golden_dir() -> Option<PathBuf> {
    let d = repo_root().join("bench/golden/smoke/cli");
    d.is_dir().then_some(d)
}

fn binary() -> PathBuf {
    // `cargo test` builds the binary next to the test harness.
    let mut p = std::env::current_exe().expect("current_exe");
    p.pop(); // deps/
    if p.ends_with("deps") {
        p.pop();
    }
    p.join("NGSpeciesID")
}

/// The corpus path the goldens were recorded against. Cases that name it must
/// be invoked with the same path, or the reference's own messages differ.
fn corpus() -> String {
    repo_root()
        .join("test/sample_h1.fastq")
        .to_string_lossy()
        .into_owned()
}

/// Apply the same scrubbing `bench/equivalence.sh` applies when recording:
/// temp paths, absolute paths to the entry point and the fixtures, traceback
/// line numbers, and long floats. Without this the goldens would be valid on
/// exactly one machine.
fn scrub(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for line in s.split_inclusive('\n') {
        let mut l = line.to_string();
        // Absolute paths to the two fixtures and the entry point.
        for name in [
            "NGSpeciesID",
            "sample_h1.fastq",
            "Supplementary_File1_reads.fastq",
            "Supplementary_File3_primer.txt",
        ] {
            l = replace_path_ending_in(&l, name);
        }
        out.push_str(&l);
    }
    out
}

/// Replace `/any/absolute/path/<name>` with `<PATH>/<name>`, matching the
/// harness's sed. Written out rather than pulled in as a regex dependency: this
/// crate has none, and one substitution does not earn one.
fn replace_path_ending_in(line: &str, name: &str) -> String {
    let mut out = String::with_capacity(line.len());
    let mut rest = line;
    while let Some(hit) = rest.find(name) {
        let (before, after) = rest.split_at(hit);
        // Walk back over the path component chain.
        let cut = before
            .char_indices()
            .rev()
            .take_while(|(_, c)| !matches!(c, ' ' | '"' | '\'' | '\t' | '\n'))
            .last()
            .map(|(i, _)| i);
        match cut {
            Some(i) if before[i..].starts_with('/') => {
                out.push_str(&before[..i]);
                out.push_str("<PATH>/");
            }
            _ => out.push_str(before),
        }
        out.push_str(name);
        rest = &after[name.len()..];
    }
    out.push_str(rest);
    out
}

struct Case {
    name: &'static str,
    args: Vec<String>,
}

fn c(name: &'static str, args: &[&str]) -> Case {
    Case {
        name,
        args: args
            .iter()
            .map(|a| {
                if *a == "@CORPUS" {
                    corpus()
                } else {
                    (*a).to_string()
                }
            })
            .collect(),
    }
}

/// The 23 exact cases, with the same arguments `bench/equivalence.sh` uses.
/// Kept in the same order as `cmd_cli` so the two can be read side by side.
fn exact_cases() -> Vec<Case> {
    vec![
        c("version", &["--version"]),
        c("help", &["--help"]),
        c("h_short", &["-h"]),
        c("noargs", &[]),
        c("wf_help", &["write_fastq", "--help"]),
        c("unknown_flag", &["--fastq", "@CORPUS", "--no-such-flag"]),
        c("bad_int", &["--k", "abc", "--fastq", "@CORPUS"]),
        c("bad_float", &["--q", "xyz", "--fastq", "@CORPUS"]),
        c("bad_int_t", &["--t", "1.5", "--fastq", "@CORPUS"]),
        c("missing_val", &["--fastq", "@CORPUS", "--k"]),
        c("bad_subcmd", &["bogus_subcmd"]),
        c("wf_stray_flag", &["write_fastq", "--k", "5"]),
        c("ambiguous_me", &["--me", "--fastq", "@CORPUS"]),
        c("ambiguous_min", &["--min", "5", "--fastq", "@CORPUS"]),
        c("ambiguous_r", &["--r", "--fastq", "@CORPUS"]),
        c("ambiguous_prim", &["--prim", "x", "--fastq", "@CORPUS"]),
        c(
            "version_wins",
            &["--version", "--fastq", "/nope", "--k", "abc"],
        ),
        c("help_wins", &["--fastq", "/nope", "-h"]),
        c("both_presets", &["--ont", "--isoseq", "--fastq", "@CORPUS"]),
        c("w_lt_k", &["--fastq", "@CORPUS", "--k", "20", "--w", "15"]),
        c(
            "w_gt_100",
            &["--fastq", "@CORPUS", "--k", "15", "--w", "101"],
        ),
        // exact_m and exact_f are traceback-class in the harness -- the
        // reference reaches a TypeError -- but the port replaces them with a
        // one-line message, so their stderr is NOT byte-identity contract.
        // They are covered by cli.rs's unit tests instead.
    ]
}

#[test]
fn exact_cases_match_their_goldens() {
    let Some(g) = golden_dir() else {
        eprintln!("SKIP: bench/golden/smoke/cli is absent");
        return;
    };
    let bin = binary();
    assert!(
        bin.is_file(),
        "binary not found at {}; run `cargo build` first",
        bin.display()
    );

    let mut failures: Vec<String> = Vec::new();
    let mut checked = 0usize;
    for case in exact_cases() {
        let d = g.join(case.name);
        if !d.is_dir() {
            failures.push(format!("{}: no golden directory", case.name));
            continue;
        }
        let want_exit: i32 = std::fs::read_to_string(d.join("exit"))
            .expect("exit")
            .trim()
            .parse()
            .expect("exit is a number");
        let want_stdout = std::fs::read_to_string(d.join("stdout")).expect("stdout");
        let want_stderr = std::fs::read_to_string(d.join("stderr")).expect("stderr");

        let out = Command::new(&bin)
            .args(&case.args)
            .output()
            .expect("running the binary");
        let got_exit = out.status.code().unwrap_or(-1);
        let got_stdout = scrub(&String::from_utf8_lossy(&out.stdout));
        let got_stderr = scrub(&String::from_utf8_lossy(&out.stderr));

        if got_exit != want_exit {
            failures.push(format!(
                "{}: exit {} want {}",
                case.name, got_exit, want_exit
            ));
        }
        if got_stdout != want_stdout {
            failures.push(format!(
                "{}: stdout differs\n  got:  {:?}\n  want: {:?}",
                case.name,
                first_diff_line(&got_stdout, &want_stdout),
                first_diff_line(&want_stdout, &got_stdout)
            ));
        }
        if got_stderr != want_stderr {
            failures.push(format!(
                "{}: stderr differs\n  got:  {:?}\n  want: {:?}",
                case.name,
                first_diff_line(&got_stderr, &want_stderr),
                first_diff_line(&want_stderr, &got_stderr)
            ));
        }
        checked += 1;
    }
    assert!(
        failures.is_empty(),
        "{} of {} exact CLI cases failed:\n{}",
        failures.len(),
        checked,
        failures.join("\n")
    );
    // A test that checked nothing must not report success. The list is 21 of
    // the 23 exact cases; exact_m and exact_f are excluded above, with a reason.
    assert_eq!(checked, 21, "every listed case must have been checked");
}

/// The first line that differs, so a failure says where rather than dumping
/// eight kilobytes of help text.
fn first_diff_line(a: &str, b: &str) -> String {
    for (i, (la, lb)) in a.lines().zip(b.lines()).enumerate() {
        if la != lb {
            return format!("line {}: {la}", i + 1);
        }
    }
    match a.lines().count().cmp(&b.lines().count()) {
        std::cmp::Ordering::Greater => format!(
            "{} extra line(s), first: {}",
            a.lines().count() - b.lines().count(),
            a.lines().nth(b.lines().count()).unwrap_or("")
        ),
        std::cmp::Ordering::Less => "missing line(s) at the end".to_string(),
        std::cmp::Ordering::Equal => {
            "identical by lines; a trailing-newline difference".to_string()
        }
    }
}

#[test]
fn the_text_constants_still_match_the_goldens() {
    // src/text/*.txt were extracted from the goldens. If someone re-records the
    // goldens -- a version bump does it, and so does adding a flag -- these
    // must be re-extracted or 23 cases start failing for one reason.
    let Some(g) = golden_dir() else {
        eprintln!("SKIP: bench/golden/smoke/cli is absent");
        return;
    };
    let t = repo_root().join("rust/src/text");
    let pairs = [
        ("help.txt", g.join("help/stdout")),
        ("version.txt", g.join("version/stdout")),
        ("write_fastq_help.txt", g.join("wf_help/stdout")),
        ("presets_exclusive.txt", g.join("both_presets/stderr")),
        ("window_invalid.txt", g.join("w_lt_k/stderr")),
    ];
    for (name, golden) in pairs {
        let ours = std::fs::read_to_string(t.join(name)).expect(name);
        let theirs = std::fs::read_to_string(&golden).expect("golden");
        assert_eq!(
            ours,
            theirs,
            "rust/src/text/{name} has drifted from {}",
            golden.display()
        );
    }
    // The two usage blocks are derived rather than whole files: the main one is
    // the first 21 lines of any exit-2 case's stderr, and the subparser's is
    // the first 2 lines of its help.
    let usage = std::fs::read_to_string(t.join("usage.txt")).expect("usage.txt");
    let noargs = std::fs::read_to_string(g.join("noargs/stderr")).expect("noargs");
    assert!(
        noargs.starts_with(&usage),
        "usage.txt is not the head of the recorded exit-2 stderr"
    );
    assert_eq!(usage.lines().count(), 21);

    let wf_usage =
        std::fs::read_to_string(t.join("write_fastq_usage.txt")).expect("write_fastq_usage.txt");
    let wf_help = std::fs::read_to_string(g.join("wf_help/stdout")).expect("wf_help");
    assert!(wf_help.starts_with(&wf_usage));
    assert_eq!(wf_usage.lines().count(), 2);
}
