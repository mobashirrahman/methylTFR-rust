//! The six methylation file formats upstream's `read_methylome()` accepts
//! (`docs/AGENT_PLAN.md` section 2.2).
//!
//! The row each format implements, with columns 1-based:
//!
//! | type | header | chr | start | end | strand | score | coverage | min cols |
//! |---|---|---|---|---|---|---|---|---|
//! | `epp` | none | c1 | c2 | c3 | c6 | `c5 / 1000` | number after `/` in c4 | 6 |
//! | `bissnp` | line 1 always skipped | c1 | c2 | c3 | c6 | `c4 / 100` | c5 | 6 |
//! | `allc` | none | c1 | c2 | c2 | c3 | `c5 / c6` | c6 | 6 |
//! | `bismarkcytosine` | none | c1 | c2 | c2 | c3 | `c4 / (c4 + c5)` | `c4 + c5` | 5 |
//! | `bismarkcov` | none | c1 | c2 | c3 | always `*` | `c4 / 100` | `c5 + c6` | 6 |
//! | `encode` | auto-detected | c1 | c2 | c3 | c6 | `c11 / 100` | c10 | 11 |
//!
//! Then, in this order: drop records whose score is NaN, drop coverage below the
//! threshold, and reject any surviving score outside `[0, 1]`.
//!
//! Scores go through [`round6`] on the way, because every upstream format rounds
//! to six decimals and the golden comparison there is bit-exact
//! (`docs/divergences.md` D2).

use std::path::Path;

use crate::error::{Error, Result};
use crate::io::{LineReader, count_fields, parse_f64, parse_i64, read_lines};
use crate::model::{ChromTable, Methylome, Site, Strand};
use crate::rmath::round6;

/// The six accepted `type` strings, in the order upstream lists them.
pub const TYPES: [&str; 6] = [
    "bissnp",
    "epp",
    "allc",
    "bismarkcytosine",
    "bismarkcov",
    "encode",
];

/// A record as each parser produces it, before the shared post-steps run.
///
/// `strand` keeps an `Invalid` variant so that the strand check happens after
/// the whole file has been read, which is where upstream's `GRanges` call does
/// it: a file with a bad strand on line 2 and a bad number on line 9 must
/// complain about the number first.
#[derive(Debug, Clone, Copy)]
struct Raw {
    chr: u32,
    start: i64,
    end: i64,
    strand: RawStrand,
    score: f64,
    coverage: f64,
    line: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RawStrand {
    Plus,
    Minus,
    Star,
    Invalid,
}

impl RawStrand {
    fn parse(s: &str) -> Self {
        match s {
            "+" => RawStrand::Plus,
            "-" => RawStrand::Minus,
            "*" => RawStrand::Star,
            _ => RawStrand::Invalid,
        }
    }

    fn resolve(self, path: &Path, line: usize) -> Result<Strand> {
        match self {
            RawStrand::Plus => Ok(Strand::Plus),
            RawStrand::Minus => Ok(Strand::Minus),
            RawStrand::Star => Ok(Strand::Star),
            RawStrand::Invalid => Err(Error::parse(
                path,
                line,
                "strand values must be in '+' '-' '*'",
            )),
        }
    }
}

/// `read_methylome(path, type, cov_threshold)`.
///
/// `type` is matched case-insensitively, as `read_methylome` does with
/// `tolower()`. `chroms` is the run's shared chromosome table, so the ids in the
/// returned methylome are comparable with the annotation's.
pub fn read_methylome(
    path: &Path,
    ty: &str,
    cov_threshold: f64,
    chroms: &mut ChromTable,
) -> Result<Methylome> {
    let lower = ty.to_ascii_lowercase();
    if !TYPES.contains(&lower.as_str()) {
        return Err(Error::run(format!("{lower} is not a valid file type!")));
    }
    if cov_threshold.is_nan() || cov_threshold < 0.0 {
        return Err(Error::run(format!(
            "{cov_threshold} is not a valid coverage threshold!"
        )));
    }

    let mut lines = read_lines(path)?;
    let mut raw: Vec<Raw> = Vec::new();
    match lower.as_str() {
        "epp" => read_epp(&mut lines, path, chroms, &mut raw)?,
        "bissnp" => read_bissnp(&mut lines, path, chroms, &mut raw)?,
        "allc" => read_allc(&mut lines, path, chroms, &mut raw)?,
        "bismarkcytosine" => read_bismark_cytosine(&mut lines, path, chroms, &mut raw)?,
        "bismarkcov" => read_bismark_cov(&mut lines, path, chroms, &mut raw)?,
        "encode" => read_encode(&mut lines, path, chroms, &mut raw)?,
        _ => unreachable!("type list and dispatch agree"),
    }

    // Post-step 1: drop records whose score is NaN, which happens when coverage
    // is zero.
    // Post-step 2: drop records below the coverage threshold.
    let mut sites = Vec::with_capacity(raw.len());
    for r in raw {
        if r.score.is_nan() {
            continue;
        }
        if r.coverage < cov_threshold {
            continue;
        }
        sites.push(Site {
            chr: r.chr,
            start: r.start,
            end: r.end,
            strand: r.strand.resolve(path, r.line)?,
            score: r.score,
            coverage: r.coverage,
        });
    }

    // Post-step 3: what survives must be a fraction. Upstream tests
    // `is.finite(min) && (min < 0 || max > 1)`, so a -Inf minimum passes; keep
    // that rather than "fixing" it.
    if !sites.is_empty() {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for s in &sites {
            if s.score < lo {
                lo = s.score;
            }
            if s.score > hi {
                hi = s.score;
            }
        }
        if lo.is_finite() && (lo < 0.0 || hi > 1.0) {
            return Err(Error::file(
                path,
                format!(
                    "Methylation scores parsed from {} span {} to {}, which is not a \
                     fraction. Check that '{lower}' is the right format for this file.",
                    path.display(),
                    signif(lo),
                    signif(hi)
                ),
            ));
        }
    }

    Ok(Methylome { sites })
}

/// R's `signif(x, 3)`, used only inside the "not a fraction" message.
fn signif(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let mag = x.abs().log10().floor() as i32;
    format!("{:.*}", (2 - mag).max(0) as usize, x)
}

/// Upstream raises the "not enough columns" error from the parser, after
/// `fread` has already read the file, so it has no line number either.
fn require_cols(path: &Path, found: usize, min: usize, ty: &str) -> Result<()> {
    if found < min {
        return Err(Error::file(
            path,
            format!("{ty} file must contain at least {min} columns!"),
        ));
    }
    Ok(())
}

/// EPP: `chr start end mc/cov score_per_mille strand`.
///
/// Coverage is the number after the `/` in column 4 with any `'` stripped.
/// Upstream splits the whole column on `/` and keeps every even element, which
/// agrees with taking the last field for a well-formed row.
fn read_epp<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        require_cols(path, f.len(), 6, "epp")?;
        let cov = epp_coverage(f[3], path, no)?;
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: RawStrand::parse(f[5]),
            score: round6(parse_f64(f[4], path, no, "score")? / 1000.0),
            coverage: cov,
            line: no,
        });
    }
    Ok(())
}

fn epp_coverage(field: &str, path: &Path, line: usize) -> Result<f64> {
    let after = field.rsplit('/').next().unwrap_or(field);
    parse_f64(after.replace('\'', "").trim(), path, line, "coverage")
}

/// BisSNP: one header line is always skipped, whether or not one is present.
fn read_bissnp<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    lines.next_line()?;
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        require_cols(path, f.len(), 6, "bissnp")?;
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: RawStrand::parse(f[5]),
            score: round6(parse_f64(f[3], path, no, "score")? / 100.0),
            coverage: parse_f64(f[4], path, no, "coverage")?,
            line: no,
        });
    }
    Ok(())
}

/// allc: `chr pos strand context meth cov`.
///
/// Column 4 is the methylation context and upstream does *not* filter on it.
/// Every context is kept; that is faithful, see `docs/divergences.md`.
fn read_allc<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        require_cols(path, f.len(), 6, "allc")?;
        let cov = parse_f64(f[5], path, no, "coverage")?;
        let meth = parse_f64(f[4], path, no, "methylated")?;
        let pos = parse_i64(f[1], path, no, "start")?;
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: pos,
            end: pos,
            strand: RawStrand::parse(f[2]),
            score: round6(meth / cov),
            coverage: cov,
            line: no,
        });
    }
    Ok(())
}

/// bismarkCytosine: `chr pos strand meth unmeth total context`.
///
/// Coverage is the sum, so a `0/0` row divides to NaN and is dropped by the
/// first post-step.
fn read_bismark_cytosine<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        require_cols(path, f.len(), 5, "bismarkCytosine")?;
        let meth = parse_f64(f[3], path, no, "methylated")?;
        let unmeth = parse_f64(f[4], path, no, "unmethylated")?;
        let cov = meth + unmeth;
        let pos = parse_i64(f[1], path, no, "start")?;
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: pos,
            end: pos,
            strand: RawStrand::parse(f[2]),
            score: round6(meth / cov),
            coverage: cov,
            line: no,
        });
    }
    Ok(())
}

/// bismarkCov: `chr start end percent meth unmeth`. Strand is always `*`.
fn read_bismark_cov<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        require_cols(path, f.len(), 6, "bismarkCov")?;
        let meth = parse_f64(f[4], path, no, "methylated")?;
        let unmeth = parse_f64(f[5], path, no, "unmethylated")?;
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: RawStrand::Star,
            score: round6(parse_f64(f[3], path, no, "percent")? / 100.0),
            coverage: meth + unmeth,
            line: no,
        });
    }
    Ok(())
}

/// ENCODE bedMethyl: 11 columns, coverage in c10, percent methylated in c11.
///
/// The header is auto-detected. Upstream uses `fread(header = "auto")`, whose
/// inference inspects the whole column; this crate uses the single-field rule
/// recorded as `docs/divergences.md` D1.
fn read_encode<R: std::io::BufRead>(
    lines: &mut LineReader<R>,
    path: &Path,
    chroms: &mut ChromTable,
    out: &mut Vec<Raw>,
) -> Result<()> {
    let mut first = true;
    while let Some((no, line)) = lines.next_line()? {
        let f: Vec<&str> = line.split('\t').collect();
        if first {
            first = false;
            if count_fields(line) >= 3 && looks_like_header(line) {
                continue;
            }
        }
        require_cols(path, f.len(), 11, "encode")?;
        let cov = parse_f64(f[9], path, no, "coverage")?;
        let pct = parse_f64(f[10], path, no, "percentMeth")?;
        if !(0.0..=100.0).contains(&pct) {
            return Err(Error::file(
                path,
                format!(
                    "Column 11 of {} is outside 0-100, so it is not a methylation \
                     percentage. Check that this is a bedMethyl file.",
                    path.display()
                ),
            ));
        }
        out.push(Raw {
            chr: chroms.intern(f[0]),
            start: parse_i64(f[1], path, no, "start")?,
            end: parse_i64(f[2], path, no, "end")?,
            strand: RawStrand::parse(f[5]),
            score: round6(pct / 100.0),
            coverage: cov,
            line: no,
        });
    }
    Ok(())
}

/// Divergence D1: the first line is a header iff its 2nd or 3rd field is not an
/// integer.
fn looks_like_header(line: &str) -> bool {
    let f: Vec<&str> = line.split('\t').collect();
    let bad = |v: &&str| v.parse::<i64>().is_err();
    f.get(1).is_some_and(bad) || f.get(2).is_some_and(bad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_header_detection() {
        assert!(looks_like_header("chrom\tstart\tend\tname"));
        assert!(!looks_like_header("chr1\t100\t101\tname"));
        assert!(!looks_like_header("chr1\t100\t101"));
        assert!(looks_like_header("chr1\t100.5\t101\tname"));
    }

    #[test]
    fn epp_coverage_takes_the_denominator() {
        assert_eq!(epp_coverage("27/27", Path::new("f"), 1).unwrap(), 27.0);
        assert_eq!(
            epp_coverage("1'000/2'000", Path::new("f"), 1).unwrap(),
            2000.0
        );
    }

    #[test]
    fn signif_formats_three_significant_digits() {
        assert_eq!(signif(0.0), "0");
        assert_eq!(signif(1.23456), "1.23");
        assert_eq!(signif(-12.3456), "-12.3");
        assert_eq!(signif(1.0), "1.00");
    }

    #[test]
    fn unknown_type_is_rejected() {
        let mut chroms = ChromTable::new();
        let e = read_methylome(Path::new("nope"), "wiggle", 1.0, &mut chroms).unwrap_err();
        assert!(e.to_string().contains("is not a valid file type"));
        let e = read_methylome(Path::new("nope"), "epp", -1.0, &mut chroms).unwrap_err();
        assert!(e.to_string().contains("not a valid coverage threshold"));
    }

    #[test]
    fn type_matching_is_case_insensitive() {
        let mut chroms = ChromTable::new();
        // Reaches the file open, so the dispatch did not reject the name.
        let e = read_methylome(Path::new("/nonexistent.tsv"), "EPP", 1.0, &mut chroms).unwrap_err();
        assert!(e.to_string().contains("nonexistent"));
    }
}
