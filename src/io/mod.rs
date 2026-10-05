//! Input readers (`docs/AGENT_PLAN.md` sections 2.2 and 3).
//!
//! Two entry points:
//!
//! * [`portable`] reads the tab-separated fixture and annotation formats.
//! * [`methylome`] reads the six methylation file formats upstream's
//!   `read_methylome()` accepts, and applies its three post-steps in order.

pub mod methylome;
pub mod portable;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use flate2::read::MultiGzDecoder;

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Read buffer size. Large enough that a 268 717-row TFBS file and a genome-scale
/// methylome both stream with very few syscalls, small enough to be irrelevant
/// next to the data itself.
pub const BUF_CAPACITY: usize = 1 << 20;

/// Open a file, transparently decompressing it when it starts with the gzip
/// magic bytes `1f 8b`. Upstream relies on `data.table::fread` doing the same
/// thing, and on nothing else.
pub fn open_maybe_gzipped(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|e| Error::io(path, e))?;
    let mut reader = BufReader::with_capacity(BUF_CAPACITY, file);
    // `fill_buf` peeks without consuming, so the magic bytes stay in the stream
    // for whichever decoder we hand the buffer to.
    let gzipped = reader
        .fill_buf()
        .map_err(|e| Error::io(path, e))
        .is_ok_and(|head| head.starts_with(&GZIP_MAGIC));
    if gzipped {
        Ok(Box::new(BufReader::with_capacity(
            BUF_CAPACITY,
            MultiGzDecoder::new(reader),
        )))
    } else {
        Ok(Box::new(reader))
    }
}

/// Buffered line reader over an already-opened stream.
///
/// Reuses one `Vec<u8>` for every line, so parsing a file allocates nothing per
/// record. Returns `&str` slices borrowed from that buffer, which is why the
/// borrow lives on `&mut self`: only one line is valid at a time.
pub struct LineReader<R: BufRead> {
    inner: R,
    buf: Vec<u8>,
    done: bool,
    /// 1-based number of the line most recently returned, for errors.
    line_no: usize,
}

impl<R: BufRead> LineReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::with_capacity(256),
            done: false,
            line_no: 0,
        }
    }

    /// Number of the line most recently returned.
    pub fn line_no(&self) -> usize {
        self.line_no
    }

    /// The next line without its terminator, or `None` at end of file.
    ///
    /// Returns the 1-based line number alongside the text so that callers can
    /// build errors without borrowing `self` while holding the line.
    ///
    /// `\r\n` is trimmed, because the fixtures are checked out on every
    /// platform and a stray `\r` would otherwise be parsed into the last column.
    pub fn next_line(&mut self) -> Result<Option<(usize, &str)>> {
        if self.done {
            return Ok(None);
        }
        self.buf.clear();
        let n = self
            .inner
            .read_until(b'\n', &mut self.buf)
            .map_err(|e| Error::io("<stream>", e))?;
        if n == 0 {
            self.done = true;
            return Ok(None);
        }
        self.line_no += 1;
        if self.buf.last() == Some(&b'\n') {
            self.buf.pop();
        }
        if self.buf.last() == Some(&b'\r') {
            self.buf.pop();
        }
        // Fixtures are plain tab-separated text with no quoting, so a UTF-8
        // check here costs a scan and saves a per-field one.
        let no = self.line_no;
        match std::str::from_utf8(&self.buf) {
            Ok(s) => Ok(Some((no, s))),
            Err(_) => Err(Error::parse(
                "<stream>",
                self.line_no,
                "line is not valid UTF-8",
            )),
        }
    }
}

/// Convenience: `LineReader` over a path, gzip-aware.
pub fn read_lines(path: &Path) -> Result<LineReader<Box<dyn BufRead>>> {
    Ok(LineReader::new(open_maybe_gzipped(path)?))
}

/// The path a reader should blame in errors, kept next to the reader so call
/// sites do not have to pass it twice.
#[derive(Clone, Debug)]
pub struct Reader {
    pub path: PathBuf,
}

impl Reader {
    pub fn new(path: impl AsRef<Path>) -> Self {
        Self {
            path: path.as_ref().to_path_buf(),
        }
    }
}

/// Split a tab-separated line into fields, borrowed in place.
///
/// `split('\t')` on the whole line keeps every field a `&str` slice of the
/// reader's buffer, so the hot parsers allocate nothing per record.
#[inline]
pub fn fields(line: &str) -> impl Iterator<Item = &str> {
    line.split('\t')
}

/// How many tab-separated fields a line has.
#[inline]
pub fn count_fields(line: &str) -> usize {
    line.split('\t').count()
}

/// Parse a whole field as `i64`, reporting the file and line on failure.
pub fn parse_i64(value: &str, path: &Path, line: usize, column: &str) -> Result<i64> {
    value.parse::<i64>().map_err(|_| {
        Error::parse(
            path,
            line,
            format!("column {column}: '{value}' is not an integer"),
        )
    })
}

/// Parse a whole field as `f64`, reporting the file and line on failure.
pub fn parse_f64(value: &str, path: &Path, line: usize, column: &str) -> Result<f64> {
    // `f64::from_str` accepts "inf", "NaN" and "1e5", which is what fread does
    // too, so this is deliberately the permissive parser.
    value.parse::<f64>().map_err(|_| {
        Error::parse(
            path,
            line,
            format!("column {column}: '{value}' is not a number"),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn lines_of(text: &str) -> Vec<String> {
        let mut r = LineReader::new(Cursor::new(text.as_bytes().to_vec()));
        let mut out = Vec::new();
        while let Some((_, l)) = r.next_line().unwrap() {
            out.push(l.to_string());
        }
        out
    }

    #[test]
    fn strips_both_line_endings() {
        assert_eq!(lines_of("a\nb\r\nc"), vec!["a", "b", "c"]);
        assert_eq!(lines_of(""), Vec::<String>::new());
    }

    #[test]
    fn keeps_empty_lines() {
        assert_eq!(lines_of("a\n\nb\n"), vec!["a", "", "b"]);
    }

    #[test]
    fn field_access() {
        let l = "a\tb\tc";
        let got: Vec<&str> = fields(l).collect();
        assert_eq!(got, vec!["a", "b", "c"]);
        assert_eq!(count_fields(l), 3);
        assert_eq!(fields("a").count(), 1);
    }

    #[test]
    fn number_errors_name_the_line() {
        let e = parse_i64("x", Path::new("f.tsv"), 3, "start").unwrap_err();
        assert_eq!(
            e.to_string(),
            "f.tsv:3: column start: 'x' is not an integer"
        );
        assert!(parse_f64("y", Path::new("f.tsv"), 1, "score").is_err());
        assert_eq!(
            parse_f64("1e3", Path::new("f.tsv"), 1, "score").unwrap(),
            1000.0
        );
    }
}
