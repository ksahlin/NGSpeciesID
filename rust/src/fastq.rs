//! `help_functions.readfq`, ported.
//!
//! Carried across from the isONclust port with ONE deliberate change, and it is
//! the single line on which the two references differ:
//!
//! **Spaces in the header are kept.** isONclust's `readfq` does
//! `last[1:].replace(" ", "_")`; NGSpeciesID's does not. Measured on 3 000 real
//! ONT reads, that substitution is the *only* textual difference between the
//! two programs' output — normalised for it, their cluster assignments are
//! byte-identical. So deleting it here is what makes the engine reproduce this
//! reference rather than the other one.
//!
//! It is load-bearing in both directions. Accessions are split on `_` and the
//! score is appended with `_`, so `float(acc.split("_")[-1])` still works: an
//! ONT header ends in a token like `h1`, and the score becomes `h1_1234.5`.
//! What it costs is that column 2 of `final_clusters.tsv` now contains spaces,
//! which is exactly what breaks `write_fastq` (PORTING.md, *Finding 6*).
//!
//! The second behaviour is where the two references differ AGAIN, in the
//! opposite direction, and this one is easy to inherit by accident:
//!
//! 2. **Lines are truncated with `l[:-1]`, not chomped.** NGSpeciesID's
//!    `readfq` removes the final character of every line *whatever it is*. On a
//!    line that ends in `\n` that is a chomp; on the last line of a file with no
//!    trailing newline it eats a real base or quality character. When that line
//!    is a quality string the record comes back with **no quality at all** and
//!    the caller dies with `TypeError: 'NoneType' object is not iterable`
//!    (PORTING.md, *Finding 12*).
//!
//!    isONclust FIXED this -- its `help_functions.py` has a `_chomp` that only
//!    removes `\n` -- and the module this file was copied from matches the fixed
//!    behaviour. NGSpeciesID did not take that fix, so the fix had to be undone
//!    here. Measured on `@r1\nACGT\n+\nIIII` with no trailing newline:
//!
//!    | | acc | seq | qual |
//!    | --- | --- | --- | --- |
//!    | NGSpeciesID | `r1` | `ACGT` | **None** |
//!    | isONclust | `r1` | `ACGT` | `IIII` |
//!
//!    The harness is what caught it: `cli/no_trailing_nl` wants exit 1 and the
//!    port was exiting 70, because with the fixed chomp the file parsed
//!    perfectly and the run carried on.
//!
//!    A CRLF file keeps its `\r` in both, since there the last character *is*
//!    the `\n`.

/// One record. `qual` is `None` for fasta input.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    pub name: String,
    pub seq: String,
    pub qual: Option<String>,
}

/// The reference's `l[:-1]`: remove the final character, **whatever it is**.
///
/// Not `strip_suffix('\n')`, and not `trim_end`. On a line ending in `\n` this
/// is a chomp; on a file's last line with no trailing newline it removes a real
/// base or quality character, which is what makes Finding 12's crash reachable.
/// isONclust fixed this in its Python and NGSpeciesID did not, so this is the
/// unfixed behaviour on purpose. See the module docs.
///
/// An empty line yields an empty string, matching Python's `""[:-1]`.
fn chop(line: &str) -> &str {
    let mut it = line.chars();
    it.next_back();
    it.as_str()
}

/// Parse fastq/fasta exactly as the reference's generator does.
// Used by the tests here and by the stage oracles; the streaming `for_each_file`
// is what the pipeline uses, because holding a whole corpus costs gigabytes.
#[allow(dead_code)]
pub fn read(text: &str) -> Vec<Record> {
    let mut out = Vec::new();
    for_each(text.split_inclusive('\n'), |r| out.push(r));
    out
}

/// The parser proper: hands each record to `f` instead of collecting them.
///
/// Lines arrive **inclusive of their trailing `\n`**, exactly as
/// `split_inclusive('\n')` yields them, because `chop` is the only thing allowed
/// to remove it -- the CRLF behaviour documented above depends on nothing else
/// being stripped.
///
/// This exists so `sorted.fastq` can be turned into reads without holding the
/// file and the parsed records in memory at the same time. Measured on
/// droso_100k, the slurp alone was 138 MB of a 839 MB heap peak; see PORTING.md.
pub fn for_each<S, I, F>(lines: I, mut f: F)
where
    S: AsRef<str>,
    I: IntoIterator<Item = S>,
    F: FnMut(Record),
{
    for_each_indexed(
        lines.into_iter().map(|l| {
            let n = l.as_ref().len();
            (l, n)
        }),
        |r, _, _| f(r),
    );
}

/// `for_each`, plus each record's byte range in the input.
///
/// Lines arrive as `(line, byte length)`, and `f` is called with the record, the
/// byte offset its header line starts at, and how many bytes the whole record
/// occupies. `write_fastq` uses this to index a fastq by accession and then read
/// records back one at a time instead of holding them all -- the offsets have to
/// come from this parser rather than a scan for `@`, because a quality line can
/// begin with `@` and only the state machine knows which is which.
pub fn for_each_indexed<S, I, F>(lines: I, mut f: F)
where
    S: AsRef<str>,
    I: IntoIterator<Item = (S, usize)>,
    F: FnMut(Record, u64, u32),
{
    let mut it = lines.into_iter();
    // The lookahead, with the offset its line started at.
    let mut last: Option<(String, u64)> = None;
    let mut pos: u64 = 0;

    loop {
        if last.is_none() {
            // Look for the next header.
            for (l, n) in it.by_ref() {
                let start = pos;
                pos += n as u64;
                let l = l.as_ref();
                if l.starts_with('>') || l.starts_with('@') {
                    last = Some((chop(l).to_string(), start));
                    break;
                }
            }
        }
        let (header, record_start) = match last.take() {
            Some(h) => h,
            None => break,
        };

        // NOT `.replace(' ', "_")`. See the module docs: isONclust substitutes
        // spaces here and NGSpeciesID does not, and that one line is the whole
        // textual difference between the two programs.
        let name = header[1..].to_string();
        let mut seq = String::new();
        let mut next_header: Option<(String, u64)> = None;
        for (l, n) in it.by_ref() {
            let start = pos;
            pos += n as u64;
            let l = l.as_ref();
            if l.starts_with('@') || l.starts_with('+') || l.starts_with('>') {
                next_header = Some((chop(l).to_string(), start));
                break;
            }
            seq.push_str(chop(l));
        }

        let is_fastq = matches!(&next_header, Some((h, _)) if h.starts_with('+'));
        if !is_fastq {
            // fasta record: it ends where the next header begins, or at EOF.
            let end = next_header.as_ref().map_or(pos, |(_, o)| *o);
            f(
                Record {
                    name,
                    seq,
                    qual: None,
                },
                record_start,
                (end - record_start) as u32,
            );
            match next_header {
                Some(h) => last = Some(h),
                None => break,
            }
            continue;
        }

        let seq_chars = seq.chars().count();
        let mut quals = String::new();
        let mut leng = 0usize;
        let mut completed = false;
        for (l, n) in it.by_ref() {
            pos += n as u64;
            let q = chop(l.as_ref());
            quals.push_str(q);
            leng += q.chars().count();
            if leng >= seq_chars {
                last = None;
                f(
                    Record {
                        // Cloned because the compiler cannot see that the
                        // `!completed` branch below is unreachable once this has
                        // run.
                        name: name.clone(),
                        seq: std::mem::take(&mut seq),
                        qual: Some(std::mem::take(&mut quals)),
                    },
                    record_start,
                    (pos - record_start) as u32,
                );
                completed = true;
                break;
            }
        }
        if !completed {
            // EOF before enough quality: the reference yields a fasta record
            // and stops entirely.
            f(
                Record {
                    name,
                    seq,
                    qual: None,
                },
                record_start,
                (pos - record_start) as u32,
            );
            break;
        }
    }
}

/// Stream a fastq/fasta file, handing each record to `f`.
///
/// Reads a line at a time rather than the whole file, and validates UTF-8 per
/// line -- equivalent to `read_to_string`, because `\n` cannot appear inside a
/// multi-byte UTF-8 sequence.
/// `for_each_file`, plus each record's byte range; see `for_each_indexed`.
// For the clustering stage, which needs each record's ordinal to key the
// minimizer database. Unused until `sweep` lands.
#[allow(dead_code)]
pub fn for_each_file_indexed<F>(path: &std::path::Path, f: F) -> std::io::Result<()>
where
    F: FnMut(Record, u64, u32),
{
    use std::io::BufRead;
    let file = std::fs::File::open(path)?;
    let mut rdr = std::io::BufReader::with_capacity(1 << 20, file);
    let mut err: Option<std::io::Error> = None;
    {
        let mut buf = Vec::new();
        let lines = std::iter::from_fn(|| {
            buf.clear();
            match rdr.read_until(b'\n', &mut buf) {
                Ok(0) => None,
                Ok(n) => match std::str::from_utf8(&buf) {
                    Ok(s) => Some((s.to_string(), n)),
                    Err(_) => {
                        err = Some(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "stream did not contain valid UTF-8",
                        ));
                        None
                    }
                },
                Err(e) => {
                    err = Some(e);
                    None
                }
            }
        });
        for_each_indexed(lines, f);
    }
    match err {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

pub fn for_each_file<F>(path: &std::path::Path, f: F) -> std::io::Result<()>
where
    F: FnMut(Record),
{
    use std::io::BufRead;
    let file = std::fs::File::open(path)?;
    let mut rdr = std::io::BufReader::with_capacity(1 << 20, file);
    let mut err: Option<std::io::Error> = None;
    {
        let mut buf = Vec::new();
        let lines = std::iter::from_fn(|| {
            buf.clear();
            match rdr.read_until(b'\n', &mut buf) {
                Ok(0) => None,
                Ok(_) => match std::str::from_utf8(&buf) {
                    Ok(s) => Some(s.to_string()),
                    Err(_) => {
                        err = Some(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "stream did not contain valid UTF-8",
                        ));
                        None
                    }
                },
                Err(e) => {
                    err = Some(e);
                    None
                }
            }
        });
        for_each(lines, f);
    }
    match err {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A record's reported byte range must be exactly the bytes of that record,
    /// so that reading `len` bytes at `start` and re-parsing yields it again.
    /// `write_fastq` indexes a fastq this way and then reads records back one at
    /// a time, so an off-by-one here would silently emit the wrong read.
    #[test]
    fn reported_byte_ranges_round_trip_through_the_parser() {
        let cases = [
            "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n",
            // a quality line that begins with '@', which a naive scan would take
            // for a header
            "@r1\nACGT\n+\n@@@@\n@r2\nTTTT\n+\nJJJJ\n",
            // multi-line sequence and quality
            "@r1\nAC\nGT\n+\nII\nII\n@r2\nTT\n+\nJJ\n",
            // fasta
            ">r1\nACGT\n>r2\nTTTT\n",
            // no trailing newline
            "@r1\nACGT\n+\nIIII",
            // junk before the first header
            "noise\n@r1\nACGT\n+\nIIII\n",
        ];
        for (i, text) in cases.iter().enumerate() {
            let mut got = Vec::new();
            for_each_indexed(
                text.split_inclusive('\n').map(|l| (l, l.len())),
                |r, at, n| got.push((r, at, n)),
            );
            let want = read(text);
            assert_eq!(got.len(), want.len(), "case {i}");
            for ((rec, at, n), w) in got.iter().zip(&want) {
                assert_eq!(rec, w, "case {i}: record");
                let slice = &text.as_bytes()[*at as usize..(*at as usize + *n as usize)];
                let reparsed = read(std::str::from_utf8(slice).expect("utf8"));
                assert_eq!(
                    reparsed.len(),
                    1,
                    "case {i}: slice {:?} should hold exactly one record",
                    std::str::from_utf8(slice).unwrap()
                );
                assert_eq!(&reparsed[0], w, "case {i}: re-parsed slice");
            }
        }
    }

    /// `for_each_file` must parse byte-for-byte what `read` parses from the
    /// whole file. The streaming path exists only to save memory; any
    /// difference between the two is a silent divergence in every run.
    ///
    /// Run against both committed fixtures, and against `$NGSPECIESID_DATA`
    /// corpora when they are present. The committed ones matter most here: real
    /// ONT headers carry eight spaces each, which is exactly the behaviour this
    /// port changed.
    #[test]
    fn streaming_and_slurping_parse_identically() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("repo root");
        let mut checked = 0usize;
        let mut paths: Vec<std::path::PathBuf> = [
            "test/sample_h1.fastq",
            "test/Supplementary_File1_reads.fastq",
        ]
        .iter()
        .map(|r| root.join(r))
        .filter(|p| p.exists())
        .collect();
        if let Ok(data) = std::env::var("NGSPECIESID_DATA") {
            // A list of one today. Kept as a list because the corpus registry
            // is the port's largest measurement gap (PORTING.md, "The corpora")
            // and the next entry goes here.
            let p = std::path::Path::new(&data).join("sirv/SIRV_real_10k.fastq");
            if p.exists() {
                paths.push(p);
            }
        }
        for path in paths {
            let text = std::fs::read_to_string(&path).expect("fixture readable");
            let want = read(&text);
            let mut got = Vec::new();
            for_each_file(&path, |r| got.push(r)).expect("streams");
            assert_eq!(got.len(), want.len(), "record count for {}", path.display());
            assert_eq!(got, want, "records differ for {}", path.display());
            checked += 1;
        }
        assert!(checked >= 1);
    }

    /// The cases the whole-file parser is subtle about, driven through the
    /// streaming path via a temporary file so both see the same bytes.
    #[test]
    fn streaming_matches_on_the_awkward_shapes() {
        let cases = [
            // no trailing newline (Finding 12)
            "@r1\nACGT\n+\nIIII",
            // multi-line sequence and quality
            "@r1\nAC\nGT\n+\nII\nII\n@r2\nTTTT\n+\nJJJJ\n",
            // fasta, no quality at all
            ">r1\nACGT\n>r2\nTTTT\n",
            // CRLF, whose \r must survive into name, seq and qual
            "@r1\r\nACGT\r\n+\r\nIIII\r\n",
            // quality shorter than the sequence: yields qual: None
            "@r1\nACGTACGT\n+\nII\n",
            // junk before the first header
            "noise\n@r1\nACGT\n+\nIIII\n",
            // empty input
            "",
        ];
        let dir = std::env::temp_dir().join("isonclust-fastq-stream-test");
        std::fs::create_dir_all(&dir).expect("tempdir");
        for (i, text) in cases.iter().enumerate() {
            let path = dir.join(format!("case{i}.fastq"));
            std::fs::write(&path, text).expect("write");
            let want = read(text);
            let mut got = Vec::new();
            for_each_file(&path, |r| got.push(r)).expect("streams");
            assert_eq!(got, want, "case {i}: {text:?}");
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn reads_a_basic_fastq() {
        let r = read("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].name, "r1");
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual.as_deref(), Some("IIII"));
        assert_eq!(r[1].name, "r2");
    }

    /// Spaces are KEPT. This is the one line on which NGSpeciesID's readfq
    /// differs from isONclust's, and getting it wrong changes every accession
    /// in every output file.
    #[test]
    fn spaces_in_the_header_are_kept() {
        let r = read("@read 1 strand=+\nACGT\n+\nIIII\n");
        assert_eq!(r[0].name, "read 1 strand=+");
    }

    /// A real ONT header, which is what makes the difference matter: it has
    /// eight spaces and four underscores, and the score is appended after the
    /// last token.
    #[test]
    fn a_real_ont_header_survives_intact() {
        let h = "@c948601e-1bd9-4039-94c5-3d8c741df65f runid=ed1de13 read=4709 ch=45 \
start_time=2020-02-10T16:11:44Z flow_cell_id=ACE547 protocol_group_id=barcode_test \
sample_id=barcodes_fish H-(1,9),H+(6,11) h1";
        let r = read(&format!("{h}\nACGT\n+\nIIII\n"));
        assert_eq!(r[0].name, &h[1..]);
        // And the round trip the rest of the pipeline depends on: append the
        // score with '_', split it back off with rsplit('_').
        let with_score = format!("{}_{}", r[0].name, 1234.5);
        assert_eq!(
            with_score
                .rsplit('_')
                .next()
                .unwrap()
                .parse::<f64>()
                .unwrap(),
            1234.5
        );
        let stripped: Vec<&str> = with_score.split('_').collect();
        assert_eq!(stripped[..stripped.len() - 1].join("_"), r[0].name);
    }

    /// Finding 12: a file whose last line has no newline loses its final
    /// character, and when that line is a quality string the record comes back
    /// with NO quality and the caller dies.
    ///
    /// This is the inverse of the isONclust port's test of the same shape. That
    /// repository fixed `readfq`; this one did not, and inheriting the fix made
    /// `cli/no_trailing_nl` fail with exit 70 instead of 1.
    ///
    /// The values below were measured against the reference, not reasoned out.
    #[test]
    fn a_missing_final_newline_eats_the_last_character() {
        let with = read("@r1\nACGT\n+\nIIII\n");
        assert_eq!(with[0].qual.as_deref(), Some("IIII"));

        // Without the newline the quality line is truncated to "III", which is
        // shorter than the sequence, so readfq yields a record with no quality
        // at all rather than a short one.
        let without = read("@r1\nACGT\n+\nIIII");
        assert_eq!(without[0].seq, "ACGT");
        assert_eq!(without[0].qual, None, "Finding 12: the caller dies on this");
        assert_ne!(with, without);
    }

    #[test]
    fn a_missing_final_newline_keeps_the_last_read_but_not_its_quality() {
        let r = read("@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].qual.as_deref(), Some("IIII"), "the first is intact");
        assert_eq!(r[1].name, "r2");
        assert_eq!(r[1].qual, None, "the last one loses its quality");
    }

    /// The same truncation on a SEQUENCE line, where it silently shortens the
    /// read rather than failing. A fasta with no trailing newline loses a base.
    #[test]
    fn a_missing_final_newline_shortens_a_fasta_sequence() {
        assert_eq!(read(">r1\nACGT\n")[0].seq, "ACGT");
        assert_eq!(read(">r1\nACGT")[0].seq, "ACG", "the T is eaten");
    }

    /// CRLF is deliberately NOT handled: the reference carries the `\r` into
    /// the name, sequence and quality, so this does too.
    #[test]
    fn crlf_carries_the_carriage_return_through() {
        let r = read("@r1\r\nACGT\r\n+\r\nIIII\r\n");
        assert_eq!(r[0].name, "r1\r");
        assert_eq!(r[0].seq, "ACGT\r");
        assert_eq!(r[0].qual.as_deref(), Some("IIII\r"));
    }

    /// Quality genuinely shorter than the sequence still yields `None`, with or
    /// without a trailing newline. That path is untouched by the fix and still
    /// crashes the reference, so the port reports it.
    #[test]
    fn quality_shorter_than_sequence_yields_none() {
        let r = read("@r1\nACGTA\n+\nIII\n");
        assert_eq!(r[0].seq, "ACGTA");
        assert_eq!(r[0].qual, None);
    }

    #[test]
    fn multiline_sequence_and_quality_are_joined() {
        let r = read("@r1\nAC\nGT\n+\nII\nII\n");
        assert_eq!(r[0].seq, "ACGT");
        assert_eq!(r[0].qual.as_deref(), Some("IIII"));
    }

    #[test]
    fn fasta_records_have_no_quality() {
        let r = read(">r1\nACGT\n>r2\nTTTT\n");
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].qual, None);
        assert_eq!(r[1].seq, "TTTT");
    }

    #[test]
    fn empty_input_yields_nothing() {
        assert!(read("").is_empty());
        assert!(read("\n\n").is_empty());
    }
}
