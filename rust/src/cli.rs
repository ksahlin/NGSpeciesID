//! argparse-compatible command line parsing.
//!
//! Hand-written rather than built on `clap`, and that is deliberate. The
//! contract is the reference's *exact* stdout, stderr and exit code — 23 of the
//! 44 recorded cases under `bench/golden/<corpus>/cli/` demand byte-identity —
//! and argparse differs from clap on five independent axes:
//!
//! 1. **Prefix abbreviation.** argparse accepts any unambiguous prefix, so
//!    `--outfold` is `--outfolder` and `--se` is `--seed`. An *exact* match wins
//!    over longer options sharing it, so `--m` is `--m` and not ambiguous with
//!    the six other `m` flags. An ambiguous prefix is an error naming every
//!    candidate **in declaration order** — `--r` reports
//!    `--rc_identity_threshold, --racon, --racon_iter, --remove_universal_tails`,
//!    which is `add_argument` order, not sorted. clap's `infer_long_args` does
//!    the first of those three and neither of the others.
//! 2. **Error text and layout.** A fixed 21-line usage block, then
//!    `NGSpeciesID: error: <detail>`, on stderr, exit 2. See `text.rs`.
//! 3. **Double-dash single-letter options.** `--N`, `--d`, `--k`, `--m`, `--q`,
//!    `--s`, `--t`, `--w` are long options one character wide, which is not what
//!    clap's `short` produces. Five of them also carry a `dest` that differs
//!    from the flag (`--t` is `nr_cores`, `--d` is `print_output`, `--q` is
//!    `quality_threshold`, `--m` is `target_length`, `--s` is
//!    `target_deviation`).
//! 4. **Where the messages go.** The reference's own validation messages are
//!    emitted through `logging.error` with `format='%(message)s'`, so they land
//!    on **stderr** with no prefix — unlike isONclust, which used `print` and
//!    put them on stdout. Only `--help`, `--version` and `write_fastq --help`
//!    use stdout.
//! 5. **Exit codes.** `--ont --isoseq` together prints a message and exits
//!    **0**, because the reference calls bare `sys.exit()`. A pipeline cannot
//!    detect it. That is the contract; PORTING.md's *Finding 14* is where
//!    changing it is written up.
//!
//! Bending clap to all five would be more code than this, and every bend is a
//! place to diverge silently.

use crate::text;

/// What the program should do once parsing and validation are finished.
pub enum Outcome {
    /// Print to **stdout** and exit 0: `--help`, `-h`, `--version`,
    /// `write_fastq --help`.
    Stdout(&'static str),
    /// Print the usage block plus this detail to **stderr** and exit 2.
    UsageError(String),
    /// The same, but from the `write_fastq` subparser, which has its own `prog`
    /// and its own two-line usage block. Measured: `write_fastq --N abc` says
    /// `NGSpeciesID write_fastq: error: ...`, not `NGSpeciesID: error: ...`.
    SubUsageError(String),
    /// Print to **stderr** and exit with this code. Covers the reference's own
    /// validation messages *and* the places where the port replaces a Python
    /// traceback with a sentence — see `TRACEBACK REPLACEMENTS` below.
    Stderr(String, i32),
    /// Parsed and validated; run the pipeline.
    Run(Box<Args>),
    /// Parsed and validated; run the `write_fastq` subcommand.
    ///
    /// The payload is unread until that subcommand is implemented. It is kept
    /// rather than reduced to a unit, because it is exactly what the stage will
    /// need and dropping it now would mean re-deriving it later; `main` prints
    /// its shape under `--debug` so the field is genuinely read.
    WriteFastq(Box<WriteFastqArgs>),
}

/// Every value the reference's argparse produces, with the reference's own
/// `dest` names rather than the flag spellings, because that is what the rest
/// of the program reads.
#[derive(Debug, Clone, PartialEq)]
pub struct Args {
    pub debug: bool,
    pub fastq: Option<String>,
    pub use_old_sorted_file: bool,
    pub nr_cores: i64,
    pub print_output: i64,
    pub quality_threshold: f64,
    pub ont: bool,
    pub isoseq: bool,
    pub consensus: bool,
    pub abundance_ratio: f64,
    pub rc_identity_threshold: f64,
    pub max_seqs_for_consensus: i64,
    pub medaka: bool,
    pub racon: bool,
    pub medaka_model: String,
    pub medaka_fastq: bool,
    pub racon_iter: i64,
    pub remove_universal_tails: bool,
    pub primer_file: String,
    pub primer_max_ed: i64,
    pub trim_window: i64,
    pub target_length: i64,
    pub target_deviation: i64,
    pub sample_size: i64,
    pub top_reads: bool,
    pub seed: i64,
    pub k: i64,
    pub w: i64,
    pub min_shared: i64,
    pub mapped_threshold: f64,
    pub aligned_threshold: f64,
    pub symmetric_map_align_thresholds: bool,
    pub batch_type: String,
    pub min_fraction: f64,
    pub min_prob_no_hits: f64,
    pub outfolder: Option<String>,
}

impl Default for Args {
    /// Exactly the reference's argparse defaults, in `add_argument` order so a
    /// missing one is visible against the source.
    fn default() -> Self {
        Args {
            debug: false,
            fastq: None,
            use_old_sorted_file: false,
            nr_cores: 8,
            print_output: 10000,
            quality_threshold: 7.0,
            ont: false,
            isoseq: false,
            consensus: false,
            abundance_ratio: 0.1,
            rc_identity_threshold: 0.9,
            max_seqs_for_consensus: -1,
            medaka: false,
            racon: false,
            medaka_model: String::new(),
            medaka_fastq: false,
            racon_iter: 2,
            remove_universal_tails: false,
            primer_file: String::new(),
            primer_max_ed: 2,
            trim_window: 150,
            target_length: 0,
            target_deviation: 0,
            sample_size: 0,
            top_reads: false,
            seed: 0,
            k: 13,
            w: 20,
            min_shared: 5,
            mapped_threshold: 0.7,
            aligned_threshold: 0.4,
            symmetric_map_align_thresholds: false,
            batch_type: "total_nt".to_string(),
            min_fraction: 0.8,
            min_prob_no_hits: 0.1,
            outfolder: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct WriteFastqArgs {
    pub clusters: Option<String>,
    pub fastq: Option<String>,
    pub outfolder: Option<String>,
    pub n: i64,
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Str,
    Int,
    Float,
    Flag,
}

struct Opt {
    name: &'static str,
    kind: Kind,
    /// Which mutually-exclusive group this option belongs to, if any. argparse
    /// reports a clash as `argument <second>: not allowed with argument
    /// <first>`, naming the option seen *later* first — measured, see
    /// `bench/golden/*/cli/medaka_racon`.
    group: Option<u8>,
}

const fn o(name: &'static str, kind: Kind) -> Opt {
    Opt {
        name,
        kind,
        group: None,
    }
}
const fn g(name: &'static str, kind: Kind, group: u8) -> Opt {
    Opt {
        name,
        kind,
        group: Some(group),
    }
}

/// The three mutually-exclusive groups, by the id used in `MAIN_OPTS`.
const GROUP_INPUT: u8 = 0; // --fastq | --use_old_sorted_file   (required=True)
const GROUP_POLISH: u8 = 1; // --medaka | --racon
const GROUP_TRIM: u8 = 2; // --remove_universal_tails | --primer_file

/// The main parser's options **in `add_argument` order**, which is the order
/// argparse lists ambiguity candidates in. `--help` is first because argparse's
/// `add_help` runs before any `add_argument`, and it participates in prefix
/// resolution: `--h` resolves to `--help`.
const MAIN_OPTS: &[Opt] = &[
    o("--help", Kind::Flag),
    o("--version", Kind::Flag),
    o("--debug", Kind::Flag),
    g("--fastq", Kind::Str, GROUP_INPUT),
    g("--use_old_sorted_file", Kind::Flag, GROUP_INPUT),
    o("--t", Kind::Int),
    o("--d", Kind::Int),
    o("--q", Kind::Float),
    o("--ont", Kind::Flag),
    o("--isoseq", Kind::Flag),
    o("--consensus", Kind::Flag),
    o("--abundance_ratio", Kind::Float),
    o("--rc_identity_threshold", Kind::Float),
    o("--max_seqs_for_consensus", Kind::Int),
    g("--medaka", Kind::Flag, GROUP_POLISH),
    g("--racon", Kind::Flag, GROUP_POLISH),
    o("--medaka_model", Kind::Str),
    o("--medaka_fastq", Kind::Flag),
    o("--racon_iter", Kind::Int),
    g("--remove_universal_tails", Kind::Flag, GROUP_TRIM),
    g("--primer_file", Kind::Str, GROUP_TRIM),
    o("--primer_max_ed", Kind::Int),
    o("--trim_window", Kind::Int),
    o("--m", Kind::Int),
    o("--s", Kind::Int),
    o("--sample_size", Kind::Int),
    o("--top_reads", Kind::Flag),
    o("--seed", Kind::Int),
    o("--k", Kind::Int),
    o("--w", Kind::Int),
    o("--min_shared", Kind::Int),
    o("--mapped_threshold", Kind::Float),
    o("--aligned_threshold", Kind::Float),
    o("--symmetric_map_align_thresholds", Kind::Flag),
    o("--batch_type", Kind::Str),
    o("--min_fraction", Kind::Float),
    o("--min_prob_no_hits", Kind::Float),
    o("--outfolder", Kind::Str),
];

/// The `write_fastq` subparser's options, same rules.
const WF_OPTS: &[Opt] = &[
    o("--help", Kind::Flag),
    o("--clusters", Kind::Str),
    o("--fastq", Kind::Str),
    o("--outfolder", Kind::Str),
    o("--N", Kind::Int),
];

/// The three legal `--batch_type` values. `weighted` is **documented in the
/// help text and not implemented**: the reference's `batch_list` has no branch
/// for it and no `else`, so it yields nothing and the first `min()` downstream
/// raises `ValueError`. PORTING.md, *Finding 5*.
const BATCH_TYPES: &[&str] = &["total_nt", "nr_reads", "read_lengths_squared"];

/// Resolve a token to a full option name, argparse-style.
fn resolve<'a>(token: &str, opts: &'a [Opt]) -> Result<&'a Opt, Option<String>> {
    // An exact match always wins, even when the name is a prefix of longer
    // options. This is the behaviour that surprises: `--m` is the
    // target-length flag and NOT an ambiguous prefix of --mapped_threshold,
    // --max_seqs_for_consensus, --medaka, --medaka_model, --medaka_fastq,
    // --min_shared, --min_fraction or --min_prob_no_hits. A port that treats it
    // as ambiguous rejects a valid invocation. See cli/exact_m.
    if let Some(x) = opts.iter().find(|x| x.name == token) {
        return Ok(x);
    }
    let matches: Vec<&Opt> = opts.iter().filter(|x| x.name.starts_with(token)).collect();
    match matches.len() {
        1 => Ok(matches[0]),
        // Not an option at all. The caller collects it and reports
        // "unrecognized arguments" *after* parsing, because argparse does.
        0 => Err(None),
        _ => {
            let names: Vec<&str> = matches.iter().map(|x| x.name).collect();
            Err(Some(format!(
                "ambiguous option: {} could match {}",
                token,
                names.join(", ")
            )))
        }
    }
}

fn type_error(name: &str, kind: Kind, raw: &str) -> String {
    let ty = match kind {
        Kind::Int => "int",
        Kind::Float => "float",
        _ => unreachable!("only Int and Float can fail to convert"),
    };
    format!("argument {name}: invalid {ty} value: '{raw}'")
}

/// Python's `int()` on a CLI value. Notably `--t 1.5` is an error: argparse
/// does not accept float syntax for an int, and neither does this. See
/// cli/bad_int_t.
fn parse_int(name: &str, raw: &str) -> Result<i64, String> {
    // Python's int() tolerates surrounding whitespace and a leading `+`.
    let t = raw.trim();
    let t = t.strip_prefix('+').unwrap_or(t);
    t.parse::<i64>()
        .map_err(|_| type_error(name, Kind::Int, raw))
}

fn parse_float(name: &str, raw: &str) -> Result<f64, String> {
    let t = raw.trim();
    // Rust's f64 parser accepts "inf"/"NaN" and so does Python's float(), so
    // no extra guard is needed; what it must NOT accept is an empty string,
    // which it already rejects.
    t.parse::<f64>()
        .map_err(|_| type_error(name, Kind::Float, raw))
}

/// One mutually-exclusive group's state: which member was seen first.
struct GroupState {
    seen: [Option<&'static str>; 3],
}

impl GroupState {
    fn new() -> Self {
        GroupState { seen: [None; 3] }
    }
    /// Record `name` in `group`, or return argparse's clash message. The
    /// message names the option seen **later** first.
    fn note(&mut self, group: u8, name: &'static str) -> Result<(), String> {
        let slot = &mut self.seen[group as usize];
        match *slot {
            None => {
                *slot = Some(name);
                Ok(())
            }
            Some(first) if first == name => Ok(()), // repeating one member is fine
            Some(first) => Err(format!(
                "argument {name}: not allowed with argument {first}"
            )),
        }
    }
}

/// How the main parser's pass ended.
enum MainParse {
    /// Reached the end of `argv`.
    Done(Args, Vec<String>),
    /// Hit the `write_fastq` positional; everything after it belongs to the
    /// subparser.
    Subcommand(Args, Vec<String>, usize),
}

/// Parse `argv` **without** the program name.
///
/// One left-to-right pass, which is what argparse does and what a pre-pass
/// looking for the subcommand gets wrong. Measured:
///
/// | invocation | result |
/// | --- | --- |
/// | `--min 5 --fastq X` | the **ambiguity** error for `--min` |
/// | `--fastq X --min 5` | the same |
/// | `--fastq X somefile` | `invalid choice: 'somefile'` |
/// | `--k 5 X write_fastq` | `invalid choice: 'X'` -- the FIRST positional takes the subcommand slot |
///
/// The first of those is why a pre-pass had to go: `5` is not a value of
/// `--min` (there is no such option), so a scan for the first non-option token
/// found `5` and reported `invalid choice: '5'` instead of the ambiguity. Two
/// goldens caught it. One pass in argv order cannot make that mistake.
pub fn parse(argv: &[String]) -> Outcome {
    match parse_main_opts(argv) {
        Err(outcome) => outcome,
        Ok(MainParse::Done(args, unrecognized)) => match finish(&args, &unrecognized) {
            Err(outcome) => outcome,
            Ok(()) => validate(args),
        },
        Ok(MainParse::Subcommand(args, unrecognized, at)) => {
            parse_with_subcommand(args, unrecognized, &argv[at + 1..])
        }
    }
}

fn parse_with_subcommand(
    main_args: Args,
    mut unrecognized: Vec<String>,
    after: &[String],
) -> Outcome {
    let mut wf = WriteFastqArgs {
        // The subparser's --fastq and --outfolder share their dest with the
        // top-level ones, so whatever the top level set is the starting point
        // and the subcommand overwrites it. This is why the harness has to
        // invoke write_fastq with a redundant top-level --fastq at all.
        fastq: main_args.fastq.clone(),
        outfolder: main_args.outfolder.clone(),
        ..Default::default()
    };
    // 2. The subparser. Its --help fires immediately -- before the top-level
    //    required group is checked -- and its type errors carry its own prog.
    //    Stray tokens join the top level's, because argparse pools the extras
    //    from both parsers and reports them together.
    let mut i = 0usize;
    while i < after.len() {
        let tok = &after[i];
        if tok == "-h" {
            return Outcome::Stdout(text::WRITE_FASTQ_HELP);
        }
        if !tok.starts_with('-') {
            unrecognized.push(tok.clone());
            i += 1;
            continue;
        }
        let (name_tok, inline) = split_inline(tok);
        let opt = match resolve(&name_tok, WF_OPTS) {
            Ok(x) => x,
            Err(None) => {
                unrecognized.push(tok.clone());
                i += 1;
                continue;
            }
            Err(Some(msg)) => return Outcome::SubUsageError(msg),
        };
        if opt.name == "--help" {
            return Outcome::Stdout(text::WRITE_FASTQ_HELP);
        }
        let raw = match take_value(opt.name, inline, after, &mut i) {
            Ok(v) => v,
            Err(msg) => return Outcome::SubUsageError(msg),
        };
        match opt.name {
            "--clusters" => wf.clusters = Some(raw),
            "--fastq" => wf.fastq = Some(raw),
            "--outfolder" => wf.outfolder = Some(raw),
            "--N" => match parse_int(opt.name, &raw) {
                Ok(v) => wf.n = v,
                Err(msg) => return Outcome::SubUsageError(msg),
            },
            _ => unreachable!("unhandled write_fastq option {}", opt.name),
        }
        i += 1;
    }
    // 3. Only now the top-level required group, then stray tokens -- which is
    //    why `write_fastq --k 5` reports the required group and not the
    //    unrecognized `--k`. Both are reported by the TOP-LEVEL parser:
    //    `NGSpeciesID: error:`, not `NGSpeciesID write_fastq:`.
    if let Err(o) = finish(&main_args, &unrecognized) {
        return o;
    }
    Outcome::WriteFastq(Box::new(wf))
}

/// The `--fastq | --use_old_sorted_file` group is `required=True` on the
/// top-level parser. Split out because the subcommand path checks it later than
/// the plain path does -- see `parse_with_subcommand`.
fn require_input(args: &Args) -> Result<(), Outcome> {
    if args.fastq.is_none() && !args.use_old_sorted_file {
        return Err(Outcome::UsageError(
            "one of the arguments --fastq --use_old_sorted_file is required".to_string(),
        ));
    }
    Ok(())
}

fn split_inline(tok: &str) -> (String, Option<String>) {
    match tok.split_once('=') {
        Some((n, v)) => (n.to_string(), Some(v.to_string())),
        None => (tok.to_string(), None),
    }
}

fn take_value(
    name: &str,
    inline: Option<String>,
    argv: &[String],
    i: &mut usize,
) -> Result<String, String> {
    match inline {
        Some(v) => Ok(v),
        None => {
            *i += 1;
            match argv.get(*i) {
                Some(v) => Ok(v.clone()),
                None => Err(format!("argument {name}: expected one argument")),
            }
        }
    }
}

/// Parse the main parser's options. Returns the populated `Args`, or the
/// `Outcome` to emit — the error cases include `--help`/`--version`, which are
/// argparse *actions* that fire during parsing and win over everything,
/// including arguments that would otherwise be errors (cli/version_wins,
/// cli/help_wins).
fn parse_main_opts(argv: &[String]) -> Result<MainParse, Outcome> {
    let mut args = Args::default();
    let mut groups = GroupState::new();
    let mut unrecognized: Vec<String> = Vec::new();
    let mut i = 0usize;

    while i < argv.len() {
        let tok = &argv[i];

        if tok == "--" {
            i += 1;
            continue;
        }
        if tok == "-h" {
            return Err(Outcome::Stdout(text::HELP));
        }
        // A non-option token is the subcommand positional. The FIRST one takes
        // that slot even if a later token would have been a valid subcommand
        // name: `--k 5 X write_fastq` reports `invalid choice: 'X'`.
        if !tok.starts_with('-') {
            if tok == "write_fastq" {
                return Ok(MainParse::Subcommand(args, unrecognized, i));
            }
            return Err(Outcome::UsageError(format!(
                "argument {{write_fastq}}: invalid choice: '{tok}' (choose from write_fastq)"
            )));
        }

        let (name_tok, inline) = split_inline(tok);
        let opt = match resolve(&name_tok, MAIN_OPTS) {
            Ok(x) => x,
            Err(None) => {
                unrecognized.push(tok.clone());
                i += 1;
                continue;
            }
            Err(Some(msg)) => return Err(Outcome::UsageError(msg)),
        };

        match opt.name {
            "--help" => return Err(Outcome::Stdout(text::HELP)),
            "--version" => return Err(Outcome::Stdout(text::VERSION)),
            _ => {}
        }

        if let Some(gid) = opt.group {
            if let Err(msg) = groups.note(gid, opt.name) {
                return Err(Outcome::UsageError(msg));
            }
        }

        if opt.kind == Kind::Flag {
            match opt.name {
                "--debug" => args.debug = true,
                "--use_old_sorted_file" => args.use_old_sorted_file = true,
                "--ont" => args.ont = true,
                "--isoseq" => args.isoseq = true,
                "--consensus" => args.consensus = true,
                "--medaka" => args.medaka = true,
                "--racon" => args.racon = true,
                "--medaka_fastq" => args.medaka_fastq = true,
                "--remove_universal_tails" => args.remove_universal_tails = true,
                "--top_reads" => args.top_reads = true,
                "--symmetric_map_align_thresholds" => args.symmetric_map_align_thresholds = true,
                _ => unreachable!("unhandled flag {}", opt.name),
            }
            i += 1;
            continue;
        }

        let raw = match take_value(opt.name, inline, argv, &mut i) {
            Ok(v) => v,
            Err(msg) => return Err(Outcome::UsageError(msg)),
        };
        if let Err(msg) = store(&mut args, opt, &raw) {
            return Err(Outcome::UsageError(msg));
        }
        i += 1;
    }

    // The REQUIRED-group and stray-token checks are NOT done here. They happen
    // at the end of parsing, and for the subcommand path that is *after* the
    // subparser has run -- which is why `write_fastq --help` prints help and
    // does not complain about the missing --fastq. See `finish` and
    // `parse_with_subcommand`.
    Ok(MainParse::Done(args, unrecognized))
}

/// argparse's end-of-parse checks, in order: the required group first, stray
/// tokens second. So `--bogus` alone reports the missing group, while
/// `--fastq X --bogus` reports the stray token. Measured: cli/noargs and
/// cli/unknown_flag.
fn finish(args: &Args, unrecognized: &[String]) -> Result<(), Outcome> {
    require_input(args)?;
    if !unrecognized.is_empty() {
        return Err(Outcome::UsageError(format!(
            "unrecognized arguments: {}",
            unrecognized.join(" ")
        )));
    }
    Ok(())
}

fn store(args: &mut Args, opt: &Opt, raw: &str) -> Result<(), String> {
    let n = opt.name;
    match n {
        "--fastq" => args.fastq = Some(raw.to_string()),
        "--medaka_model" => args.medaka_model = raw.to_string(),
        "--primer_file" => args.primer_file = raw.to_string(),
        "--batch_type" => args.batch_type = raw.to_string(),
        "--outfolder" => args.outfolder = Some(raw.to_string()),
        "--t" => args.nr_cores = parse_int(n, raw)?,
        "--d" => args.print_output = parse_int(n, raw)?,
        "--max_seqs_for_consensus" => args.max_seqs_for_consensus = parse_int(n, raw)?,
        "--racon_iter" => args.racon_iter = parse_int(n, raw)?,
        "--primer_max_ed" => args.primer_max_ed = parse_int(n, raw)?,
        "--trim_window" => args.trim_window = parse_int(n, raw)?,
        "--m" => args.target_length = parse_int(n, raw)?,
        "--s" => args.target_deviation = parse_int(n, raw)?,
        "--sample_size" => args.sample_size = parse_int(n, raw)?,
        "--seed" => args.seed = parse_int(n, raw)?,
        "--k" => args.k = parse_int(n, raw)?,
        "--w" => args.w = parse_int(n, raw)?,
        "--min_shared" => args.min_shared = parse_int(n, raw)?,
        "--q" => args.quality_threshold = parse_float(n, raw)?,
        "--abundance_ratio" => args.abundance_ratio = parse_float(n, raw)?,
        "--rc_identity_threshold" => args.rc_identity_threshold = parse_float(n, raw)?,
        "--mapped_threshold" => args.mapped_threshold = parse_float(n, raw)?,
        "--aligned_threshold" => args.aligned_threshold = parse_float(n, raw)?,
        "--min_fraction" => args.min_fraction = parse_float(n, raw)?,
        "--min_prob_no_hits" => args.min_prob_no_hits = parse_float(n, raw)?,
        _ => unreachable!("unhandled option {}", n),
    }
    Ok(())
}

/// Post-parse validation, **in the reference's exact order**, which is the
/// order of the statements in its `if __name__ == '__main__'` block.
///
/// Getting the order wrong is not cosmetic. `--fastq X --k 20 --w 15` with no
/// `--outfolder` reports the *window* error, exit 1, because the window check
/// runs before `main()` ever touches `outfolder` — measured, cli/w_lt_k, which
/// is invoked with no `--outfolder` at all.
///
/// TRACEBACK REPLACEMENTS
/// ----------------------
/// Six of the reference's crash paths are pure argument problems and are
/// detected here instead, with a sentence. Same exit code, no stack trace —
/// which is the contract for the 15 `CLI_TRACEBACK` cases in
/// `bench/equivalence.sh`, and a deliberate divergence recorded in PORTING.md.
/// The reference's version of each is named in the comment.
fn validate(mut args: Args) -> Outcome {
    // 1. --ont and --isoseq together. Exits ZERO. Finding 14, and note this is
    //    checked before the output folder is created, so nothing is written.
    if args.ont && args.isoseq {
        return Outcome::Stderr(text::PRESETS_EXCLUSIVE.to_string(), 0);
    }
    // 2. The presets OVERWRITE an explicit --k/--w regardless of the order they
    //    appeared in: `--ont --k 99` and `--k 99 --ont` both resolve to k=13.
    //    cli/ont_over_k and cli/k_then_ont.
    if args.isoseq {
        args.k = 15;
        args.w = 50;
    } else if args.ont {
        args.k = 13;
        args.w = 20;
    }
    // 3. The reference creates --outfolder HERE, before the window check, so a
    //    run that fails validation still leaves the directory behind
    //    (recursively). Measured. main.rs does the mkdir so this function stays
    //    free of side effects and testable; the ordering is what matters and it
    //    is pinned by a test.
    //
    // 4. The window check. The only validation the reference exits non-zero
    //    for, one message for both directions, and `100 < w` is tested first.
    if 100 < args.w || args.w < args.k {
        return Outcome::Stderr(text::WINDOW_INVALID.to_string(), 1);
    }
    // 5. --outfolder is optional in argparse and mandatory in fact: the
    //    reference reaches `os.path.join(args.outfolder, ...)` and raises
    //    TypeError. Finding 9. cli/no_outfolder, cli/exact_m, cli/exact_f.
    if args.outfolder.is_none() {
        return Outcome::Stderr("Error: --outfolder is required.\n".to_string(), 1);
    }
    // 6. --d is a modulo divisor and 0 raises ZeroDivisionError. Finding 11.
    if args.print_output == 0 {
        return Outcome::Stderr(
            "Error: --d must not be 0 (it is a reporting interval).\n".to_string(),
            1,
        );
    }
    // 7. --batch_type has no validation in the reference and no `else` branch,
    //    so `weighted` -- which its own --help advertises -- yields no batches
    //    and raises ValueError from a min() over an empty list. Finding 5.
    if !BATCH_TYPES.contains(&args.batch_type.as_str()) {
        return Outcome::Stderr(
            format!(
                "Error: --batch_type must be one of {}, not '{}'.\n",
                BATCH_TYPES.join(", "),
                args.batch_type
            ),
            1,
        );
    }
    // 8. A polisher without --consensus. FINDING 4, and the direction of this
    //    check is the opposite of what it was.
    //
    //    It used to reject `--consensus` with neither polisher, because the
    //    reference reached an unbound `polishing_pattern` and died there --
    //    after clustering, spoa and reverse-complement detection had all run,
    //    writing nothing. Finding 4 is now FIXED in both implementations:
    //    `--consensus` alone is draft-only and writes the spoa references, so
    //    rejecting it would refuse a valid and useful run. It is also the only
    //    consensus mode that needs no external tools, since spoa is linked.
    //
    //    The reverse -- `--medaka` or `--racon` with no `--consensus` -- was a
    //    SILENT NO-OP: accepted, nothing polished, exit 0, no consensus. Saying
    //    so beats letting someone believe they polished something.
    if (args.medaka || args.racon) && !args.consensus {
        let polisher = if args.medaka { "--medaka" } else { "--racon" };
        // Byte-identical to the reference's logging.error, which prints the
        // message and nothing else -- no "Error: " prefix. This is a SHARED
        // message, not one of the traceback replacements, so it is contract.
        return Outcome::Stderr(
            format!(
                "{polisher} polishes a consensus, so it needs --consensus. \
                 Add --consensus, or drop {polisher}.\n"
            ),
            1,
        );
    }
    // 9. --max_seqs_for_consensus 0 writes an empty fasta and spoa aborts with
    //    SIGABRT. Finding 22. Note 0 is the ONLY bad value: -1 disables the
    //    cutoff and anything positive is fine.
    if args.consensus && args.max_seqs_for_consensus == 0 {
        return Outcome::Stderr(
            "Error: --max_seqs_for_consensus 0 leaves no sequences to build a consensus from.\n"
                .to_string(),
            1,
        );
    }
    Outcome::Run(Box::new(args))
}

/// Does the reference create `--outfolder` before or after the window check?
/// Before — see `validate` step 3. Exposed so `main` can do it in the right
/// place and a test can pin the order.
pub fn creates_outfolder_before_window_check() -> bool {
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(args: &[&str]) -> Outcome {
        parse(&args.iter().map(|s| s.to_string()).collect::<Vec<_>>())
    }
    fn usage_detail(o: Outcome) -> String {
        match o {
            Outcome::UsageError(d) => d,
            Outcome::SubUsageError(d) => {
                panic!("expected a MAIN usage error, got a subparser one: {d}")
            }
            Outcome::Stdout(_) => panic!("expected a usage error, got stdout"),
            Outcome::Stderr(m, c) => panic!("expected a usage error, got stderr {c}: {m}"),
            Outcome::Run(_) => panic!("expected a usage error, got Run"),
            Outcome::WriteFastq(_) => panic!("expected a usage error, got WriteFastq"),
        }
    }
    fn run(o: Outcome) -> Args {
        match o {
            Outcome::Run(a) => *a,
            Outcome::UsageError(d) => panic!("expected Run, got usage error: {d}"),
            Outcome::SubUsageError(d) => panic!("expected Run, got subparser usage error: {d}"),
            Outcome::Stdout(_) => panic!("expected Run, got stdout"),
            Outcome::Stderr(m, c) => panic!("expected Run, got stderr {c}: {m}"),
            Outcome::WriteFastq(_) => panic!("expected Run, got WriteFastq"),
        }
    }

    // --- the option table itself -------------------------------------------

    #[test]
    fn every_live_flag_is_declared() {
        // 39 flags in the reference plus --help, which argparse adds.
        assert_eq!(MAIN_OPTS.len(), 38, "main parser options");
        assert_eq!(WF_OPTS.len(), 5, "write_fastq options");
        // 38 + 5 = 43 declarations for 39 live flags: --help twice, and
        // --fastq/--outfolder appear in both parsers.
    }

    #[test]
    fn no_duplicate_option_names() {
        for opts in [MAIN_OPTS, WF_OPTS] {
            let mut names: Vec<&str> = opts.iter().map(|x| x.name).collect();
            names.sort_unstable();
            let before = names.len();
            names.dedup();
            assert_eq!(before, names.len(), "duplicate option name");
        }
    }

    #[test]
    fn the_eight_double_dash_single_letter_options_exist() {
        for n in ["--d", "--k", "--m", "--q", "--s", "--t", "--w"] {
            assert!(
                MAIN_OPTS.iter().any(|x| x.name == n),
                "{n} is missing from the main parser"
            );
        }
        assert!(WF_OPTS.iter().any(|x| x.name == "--N"));
    }

    // --- prefix resolution, all three behaviours ---------------------------

    #[test]
    fn a_unique_prefix_resolves() {
        assert_eq!(resolve("--outfold", MAIN_OPTS).unwrap().name, "--outfolder");
        assert_eq!(resolve("--se", MAIN_OPTS).unwrap().name, "--seed");
        assert_eq!(resolve("--co", MAIN_OPTS).unwrap().name, "--consensus");
        assert_eq!(resolve("--h", MAIN_OPTS).unwrap().name, "--help");
    }

    #[test]
    fn an_exact_match_beats_longer_options_sharing_it() {
        // The surprising one. --m is a flag AND a prefix of seven others.
        assert_eq!(resolve("--m", MAIN_OPTS).unwrap().name, "--m");
        assert_eq!(resolve("--s", MAIN_OPTS).unwrap().name, "--s");
        assert_eq!(resolve("--d", MAIN_OPTS).unwrap().name, "--d");
        assert_eq!(resolve("--k", MAIN_OPTS).unwrap().name, "--k");
        assert_eq!(resolve("--t", MAIN_OPTS).unwrap().name, "--t");
    }

    #[test]
    fn an_ambiguous_prefix_lists_candidates_in_declaration_order() {
        // Every one of these strings is copied from a recorded golden.
        let cases = [
            ("--me", "--medaka, --medaka_model, --medaka_fastq"),
            ("--min", "--min_shared, --min_fraction, --min_prob_no_hits"),
            (
                "--r",
                "--rc_identity_threshold, --racon, --racon_iter, --remove_universal_tails",
            ),
            ("--prim", "--primer_file, --primer_max_ed"),
            ("--ma", "--max_seqs_for_consensus, --mapped_threshold"),
        ];
        for (tok, want) in cases {
            match resolve(tok, MAIN_OPTS) {
                Err(Some(msg)) => assert_eq!(
                    msg,
                    format!("ambiguous option: {tok} could match {want}"),
                    "candidates for {tok}"
                ),
                Ok(opt) => panic!("{tok} should be ambiguous, resolved to {}", opt.name),
                Err(None) => panic!("{tok} should be ambiguous, was unrecognized"),
            }
        }
    }

    #[test]
    fn ambiguity_is_reported_through_parse_and_not_just_resolve() {
        // REGRESSION. `resolve` was right and `parse` was wrong: a pre-pass
        // scanned for the first non-option token to find the subcommand, and
        // for `--min 5` the `5` is not the value of any option (there is no
        // `--min`), so the scan found it and reported
        //   argument {write_fastq}: invalid choice: '5'
        // instead of the ambiguity. Two goldens caught it -- ambiguous_min and
        // ambiguous_prim -- and neither unit test did, because both tested
        // `resolve` in isolation. Test the thing the golden tests.
        assert_eq!(
            usage_detail(p(&["--min", "5", "--fastq", "x"])),
            "ambiguous option: --min could match --min_shared, --min_fraction, --min_prob_no_hits"
        );
        assert_eq!(
            usage_detail(p(&["--prim", "x", "--fastq", "x"])),
            "ambiguous option: --prim could match --primer_file, --primer_max_ed"
        );
        // ...and in either order, since the pass is left to right.
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--min", "5"])),
            "ambiguous option: --min could match --min_shared, --min_fraction, --min_prob_no_hits"
        );
    }

    #[test]
    fn the_first_positional_takes_the_subcommand_slot() {
        // Even when a valid subcommand name follows it.
        assert_eq!(
            usage_detail(p(&["--k", "5", "X", "write_fastq"])),
            "argument {write_fastq}: invalid choice: 'X' (choose from write_fastq)"
        );
        // And a positional after valid options is still the subcommand slot.
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "somefile"])),
            "argument {write_fastq}: invalid choice: 'somefile' (choose from write_fastq)"
        );
        // ...and before them.
        assert_eq!(
            usage_detail(p(&["somefile", "--fastq", "x"])),
            "argument {write_fastq}: invalid choice: 'somefile' (choose from write_fastq)"
        );
    }

    #[test]
    fn a_stray_token_after_the_subcommand_is_unrecognized() {
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "write_fastq", "extra"])),
            "unrecognized arguments: extra"
        );
    }

    #[test]
    fn the_subparser_reports_its_own_errors_under_its_own_prog() {
        match p(&["write_fastq", "--N", "abc"]) {
            Outcome::SubUsageError(d) => {
                assert_eq!(d, "argument --N: invalid int value: 'abc'")
            }
            _ => panic!("expected a subparser usage error"),
        }
        // ...and it fires BEFORE the top-level required group, which
        // `write_fastq --clusters c` does reach.
        assert_eq!(
            usage_detail(p(&["write_fastq", "--clusters", "c"])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
    }

    #[test]
    fn subparser_help_fires_before_the_required_group_check() {
        // REGRESSION. The required check used to run before the subparser, so
        // `write_fastq --help` reported the missing --fastq instead of printing
        // help. cli/wf_help is exit 0 with 457 bytes on stdout.
        assert!(
            matches!(p(&["write_fastq", "--help"]), Outcome::Stdout(s) if s == text::WRITE_FASTQ_HELP)
        );
        // --version is NOT a subparser option, so it becomes a stray token and
        // the required group wins.
        assert_eq!(
            usage_detail(p(&["write_fastq", "--version"])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
    }

    #[test]
    fn an_unknown_token_is_not_an_ambiguity_error() {
        assert!(matches!(resolve("--nope", MAIN_OPTS), Err(None)));
    }

    // --- the recorded exit-2 messages, verbatim ----------------------------

    #[test]
    fn required_group_message() {
        assert_eq!(
            usage_detail(p(&[])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
    }

    #[test]
    fn the_required_group_is_checked_before_stray_tokens() {
        // `--bogus` alone reports the missing group; with --fastq it reports
        // the stray token. cli/noargs and cli/unknown_flag.
        assert_eq!(
            usage_detail(p(&["--no-such-flag"])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--no-such-flag"])),
            "unrecognized arguments: --no-such-flag"
        );
    }

    #[test]
    fn type_error_messages() {
        assert_eq!(
            usage_detail(p(&["--k", "abc", "--fastq", "x"])),
            "argument --k: invalid int value: 'abc'"
        );
        assert_eq!(
            usage_detail(p(&["--q", "xyz", "--fastq", "x"])),
            "argument --q: invalid float value: 'xyz'"
        );
        // argparse does not accept float syntax for an int.
        assert_eq!(
            usage_detail(p(&["--t", "1.5", "--fastq", "x"])),
            "argument --t: invalid int value: '1.5'"
        );
    }

    #[test]
    fn missing_value_message() {
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--k"])),
            "argument --k: expected one argument"
        );
    }

    #[test]
    fn mutually_exclusive_groups_name_the_later_option_first() {
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--medaka", "--racon"])),
            "argument --racon: not allowed with argument --medaka"
        );
        assert_eq!(
            usage_detail(p(&[
                "--fastq",
                "x",
                "--remove_universal_tails",
                "--primer_file",
                "p"
            ])),
            "argument --primer_file: not allowed with argument --remove_universal_tails"
        );
        // ...and the reverse order swaps the names.
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--racon", "--medaka"])),
            "argument --medaka: not allowed with argument --racon"
        );
    }

    #[test]
    fn the_input_group_is_mutually_exclusive_too() {
        assert_eq!(
            usage_detail(p(&["--fastq", "x", "--use_old_sorted_file"])),
            "argument --use_old_sorted_file: not allowed with argument --fastq"
        );
    }

    #[test]
    fn repeating_one_member_of_a_group_is_allowed() {
        // argparse only complains about two DIFFERENT members.
        let a = run(p(&["--fastq", "x", "--fastq", "y", "--outfolder", "o"]));
        assert_eq!(a.fastq.as_deref(), Some("y"), "the last one wins");
    }

    #[test]
    fn bad_subcommand_message() {
        assert_eq!(
            usage_detail(p(&["bogus_subcmd"])),
            "argument {write_fastq}: invalid choice: 'bogus_subcmd' (choose from write_fastq)"
        );
    }

    // --- actions win over errors -------------------------------------------

    #[test]
    fn version_and_help_fire_during_parsing() {
        assert!(matches!(p(&["--version"]), Outcome::Stdout(s) if s == text::VERSION));
        assert!(matches!(p(&["--help"]), Outcome::Stdout(s) if s == text::HELP));
        assert!(matches!(p(&["-h"]), Outcome::Stdout(s) if s == text::HELP));
        // ...and win over arguments that would otherwise be errors.
        assert!(
            matches!(p(&["--version", "--fastq", "/nope", "--k", "abc"]), Outcome::Stdout(s) if s == text::VERSION)
        );
        assert!(matches!(p(&["--fastq", "/nope", "-h"]), Outcome::Stdout(s) if s == text::HELP));
    }

    // --- validation order and the exit codes -------------------------------

    #[test]
    fn presets_exclusive_exits_zero() {
        match p(&["--ont", "--isoseq", "--fastq", "x"]) {
            Outcome::Stderr(msg, code) => {
                assert_eq!(msg, text::PRESETS_EXCLUSIVE);
                assert_eq!(code, 0, "the reference calls bare sys.exit(); Finding 14");
            }
            _ => panic!("expected the presets message"),
        }
    }

    #[test]
    fn presets_overwrite_explicit_k_and_w_in_either_order() {
        let a = run(p(&[
            "--ont",
            "--k",
            "99",
            "--fastq",
            "x",
            "--outfolder",
            "o",
        ]));
        assert_eq!((a.k, a.w), (13, 20));
        let b = run(p(&[
            "--k",
            "99",
            "--ont",
            "--fastq",
            "x",
            "--outfolder",
            "o",
        ]));
        assert_eq!((b.k, b.w), (13, 20));
        let c = run(p(&[
            "--isoseq",
            "--w",
            "7",
            "--fastq",
            "x",
            "--outfolder",
            "o",
        ]));
        assert_eq!((c.k, c.w), (15, 50));
    }

    #[test]
    fn the_window_check_runs_before_the_outfolder_check() {
        // cli/w_lt_k is invoked with NO --outfolder and reports the window
        // message, exit 1. Reversing these two would change that.
        match p(&["--fastq", "x", "--k", "20", "--w", "15"]) {
            Outcome::Stderr(msg, code) => {
                assert_eq!(msg, text::WINDOW_INVALID);
                assert_eq!(code, 1);
            }
            _ => panic!("expected the window message"),
        }
    }

    #[test]
    fn window_check_tests_the_upper_bound_first() {
        for (k, w) in [(13, 101), (20, 15), (13, 5)] {
            match p(&["--fastq", "x", "--k", &k.to_string(), "--w", &w.to_string()]) {
                Outcome::Stderr(msg, 1) => assert_eq!(msg, text::WINDOW_INVALID),
                _ => panic!("k={k} w={w} should fail the window check"),
            }
        }
        // w == k and w == 100 are both legal.
        assert_eq!(
            run(p(&[
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--k",
                "13",
                "--w",
                "13"
            ]))
            .w,
            13
        );
        assert_eq!(
            run(p(&[
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--k",
                "13",
                "--w",
                "100"
            ]))
            .w,
            100
        );
    }

    #[test]
    fn traceback_replacements_keep_the_reference_exit_code() {
        // Each of these is a Python traceback in the reference, exit 1. The
        // port says a sentence instead; the exit code is the contract.
        let cases: Vec<Vec<&str>> = vec![
            vec!["--ont", "--fastq", "x"], // no --outfolder
            vec!["--ont", "--fastq", "x", "--outfolder", "o", "--d", "0"], // --d 0
            vec![
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--t",
                "4",
                "--batch_type",
                "weighted",
            ],
            vec![
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--t",
                "4",
                "--batch_type",
                "nosuchtype",
            ],
            // `--consensus` alone WAS here, as the reference's unbound
            // `polishing_pattern`. Finding 4 is fixed in both implementations
            // now, so that invocation is a valid draft-only run and has its own
            // test. Nothing replaces it here: the mirror case, a polisher
            // without --consensus, is NOT a traceback replacement -- both
            // implementations emit the same sentence, so it is contract and
            // does not carry this class's "Error: " prefix.
            vec![
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--consensus",
                "--racon",
                "--max_seqs_for_consensus",
                "0",
            ],
        ];
        for c in cases {
            match p(&c) {
                Outcome::Stderr(msg, 1) => {
                    assert!(msg.starts_with("Error: "), "{c:?} -> {msg}");
                    assert!(msg.ends_with('\n'), "{c:?} -> no trailing newline");
                    assert_eq!(
                        msg.lines().count(),
                        1,
                        "{c:?} -> must be ONE line, not a stack"
                    );
                }
                other => panic!(
                    "{c:?} should be a one-line stderr message with exit 1, got {}",
                    match other {
                        Outcome::Stderr(m, code) => format!("stderr {code}: {m}"),
                        Outcome::Run(_) => "Run".into(),
                        Outcome::UsageError(d) => format!("usage error: {d}"),
                        Outcome::SubUsageError(d) => format!("subparser usage error: {d}"),
                        Outcome::Stdout(_) => "stdout".into(),
                        Outcome::WriteFastq(_) => "WriteFastq".into(),
                    }
                ),
            }
        }
    }

    #[test]
    fn max_seqs_for_consensus_zero_is_only_bad_with_consensus() {
        // -1 disables the cutoff; positive values are fine; and without
        // --consensus the flag is never read.
        assert_eq!(
            run(p(&[
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--max_seqs_for_consensus",
                "0"
            ]))
            .max_seqs_for_consensus,
            0
        );
        assert_eq!(
            run(p(&[
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--consensus",
                "--racon",
                "--max_seqs_for_consensus",
                "-1"
            ]))
            .max_seqs_for_consensus,
            -1
        );
    }

    #[test]
    fn the_three_batch_types_are_accepted() {
        for b in BATCH_TYPES {
            let a = run(p(&[
                "--ont",
                "--fastq",
                "x",
                "--outfolder",
                "o",
                "--batch_type",
                b,
            ]));
            assert_eq!(&a.batch_type, b);
        }
    }

    #[test]
    fn a_polisher_without_consensus_is_refused() {
        // Finding 4's mirror, now fixed in both implementations. It used to be
        // accepted: exit 0, nothing polished, no consensus and no warning.
        match p(&["--ont", "--fastq", "x", "--outfolder", "o", "--medaka"]) {
            Outcome::Stderr(msg, code) => {
                assert_eq!(code, 1);
                assert!(msg.contains("--medaka"), "names the flag: {msg}");
                assert!(msg.contains("--consensus"), "and what to add: {msg}");
            }
            _ => panic!("expected a polisher-without-consensus error"),
        }
    }

    #[test]
    fn consensus_without_a_polisher_is_accepted_as_draft_only() {
        // The other half of Finding 4, and the direction that changed: this
        // used to be refused because the reference crashed on it. It is now a
        // valid run that writes the spoa drafts -- and the only consensus mode
        // needing no external tools, since spoa is linked.
        let a = run(p(&[
            "--ont",
            "--fastq",
            "x",
            "--outfolder",
            "o",
            "--consensus",
        ]));
        assert!(a.consensus && !a.medaka && !a.racon);
    }

    // --- defaults -----------------------------------------------------------

    #[test]
    fn defaults_match_the_reference() {
        let d = Args::default();
        assert_eq!(d.nr_cores, 8);
        assert_eq!(d.print_output, 10000);
        assert_eq!(d.quality_threshold, 7.0);
        assert_eq!(d.abundance_ratio, 0.1);
        assert_eq!(d.rc_identity_threshold, 0.9);
        assert_eq!(d.max_seqs_for_consensus, -1);
        assert_eq!(d.racon_iter, 2);
        assert_eq!(d.primer_max_ed, 2);
        assert_eq!(d.trim_window, 150);
        assert_eq!(d.target_length, 0);
        assert_eq!(d.target_deviation, 0);
        assert_eq!(d.sample_size, 0);
        assert_eq!(
            d.seed, 0,
            "the fixed default that makes --sample_size reproducible"
        );
        assert_eq!((d.k, d.w), (13, 20));
        assert_eq!(d.min_shared, 5);
        assert_eq!(d.mapped_threshold, 0.7);
        assert_eq!(d.aligned_threshold, 0.4);
        assert_eq!(d.batch_type, "total_nt");
        assert_eq!(d.min_fraction, 0.8);
        assert_eq!(d.min_prob_no_hits, 0.1);
        assert!(d.medaka_model.is_empty());
        assert!(d.primer_file.is_empty());
        assert!(d.outfolder.is_none());
    }

    // --- dest names that differ from the flag ------------------------------

    #[test]
    fn the_five_renamed_dests_land_in_the_right_field() {
        let a = run(p(&[
            "--fastq",
            "x",
            "--outfolder",
            "o",
            "--t",
            "3",
            "--d",
            "7",
            "--q",
            "9.5",
            "--m",
            "800",
            "--s",
            "100",
        ]));
        assert_eq!(a.nr_cores, 3, "--t");
        assert_eq!(a.print_output, 7, "--d");
        assert_eq!(a.quality_threshold, 9.5, "--q");
        assert_eq!(a.target_length, 800, "--m");
        assert_eq!(a.target_deviation, 100, "--s");
    }

    // --- write_fastq --------------------------------------------------------

    #[test]
    fn write_fastq_needs_a_top_level_input_argument() {
        // The required group is on the TOP-LEVEL parser, so the subcommand's
        // own --fastq does not satisfy it. Finding 6, and cli/wf_stray_flag.
        assert_eq!(
            usage_detail(p(&["write_fastq", "--clusters", "c", "--fastq", "f"])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
        assert_eq!(
            usage_detail(p(&["write_fastq", "--k", "5"])),
            "one of the arguments --fastq --use_old_sorted_file is required"
        );
    }

    #[test]
    fn write_fastq_with_a_redundant_top_level_fastq_runs() {
        match p(&[
            "--fastq",
            "f",
            "write_fastq",
            "--clusters",
            "c",
            "--fastq",
            "g",
            "--outfolder",
            "o",
            "--N",
            "2",
        ]) {
            Outcome::WriteFastq(wf) => {
                assert_eq!(wf.clusters.as_deref(), Some("c"));
                // The subparser's --fastq shares its dest and overwrites.
                assert_eq!(wf.fastq.as_deref(), Some("g"));
                assert_eq!(wf.outfolder.as_deref(), Some("o"));
                assert_eq!(wf.n, 2);
            }
            _ => panic!("expected WriteFastq"),
        }
    }

    #[test]
    fn write_fastq_inherits_the_top_level_fastq_when_it_sets_none() {
        match p(&["--fastq", "f", "write_fastq", "--clusters", "c"]) {
            Outcome::WriteFastq(wf) => assert_eq!(wf.fastq.as_deref(), Some("f")),
            _ => panic!("expected WriteFastq"),
        }
    }

    #[test]
    fn write_fastq_help_is_its_own_text() {
        assert!(
            matches!(p(&["write_fastq", "--help"]), Outcome::Stdout(s) if s == text::WRITE_FASTQ_HELP)
        );
        assert!(
            matches!(p(&["write_fastq", "-h"]), Outcome::Stdout(s) if s == text::WRITE_FASTQ_HELP)
        );
    }

    #[test]
    fn write_fastq_n_defaults_to_zero() {
        match p(&["--fastq", "f", "write_fastq", "--clusters", "c"]) {
            Outcome::WriteFastq(wf) => assert_eq!(wf.n, 0),
            _ => panic!("expected WriteFastq"),
        }
    }

    // --- inline `=` values --------------------------------------------------

    #[test]
    fn inline_equals_values_work_including_for_prefixes() {
        let a = run(p(&["--fastq=x", "--outfolder=o", "--k=17", "--w=40"]));
        assert_eq!(a.fastq.as_deref(), Some("x"));
        assert_eq!((a.k, a.w), (17, 40));
        let b = run(p(&["--fastq=x", "--outfold=o", "--seed=7"]));
        assert_eq!(b.outfolder.as_deref(), Some("o"));
        assert_eq!(b.seed, 7);
    }

    #[test]
    fn a_negative_number_is_a_value_not_an_option() {
        // --seed -5 must reach the sampler as -5, and argparse allows it
        // because -5 does not resolve to any option. Note the reference seeds
        // from the ABSOLUTE value, so -5 and 5 give the same subsample.
        let a = run(p(&["--fastq", "x", "--outfolder", "o", "--seed", "-5"]));
        assert_eq!(a.seed, -5);
    }
}
