//! Shared helpers for the integration tests.
//!
//! The comparison rule comes from `AGENT_PLAN.md` section 3:
//!
//! * Everything compared against R uses `|rust - r| <= 1e-10` **and**
//!   `|rust - r| <= 1e-10 * |r|`. Both, not either: a value near zero has to be
//!   reproduced to the bit, which is why the plan can print a golden column to 7
//!   significant digits and still know it is right.
//! * Parser scores, parser coverage and the `rmath` tables are compared
//!   **bit-exactly** instead, because `docs/divergences.md` D2 makes those the
//!   places where a float-parser difference would be real but invisible at
//!   1e-10.
//!
//! [`Errors::finish`] reports the worst absolute and relative difference it saw
//! even when it passes, so a green run still says how much headroom there was.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::sync::Arc;

use methyltfr::io::portable::{read_gc_windows, read_gcfreq, read_msites, read_ranges};
use methyltfr::model::{Annotation, ChromTable, MotifAnnotation, Site, Strand};

pub const TOL_ABS: f64 = 1e-10;
pub const TOL_REL: f64 = 1e-10;

pub fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

pub fn fixture(rel: &str) -> PathBuf {
    repo_root().join("tests/fixtures").join(rel)
}

pub fn batf(rel: &str) -> PathBuf {
    fixture(&format!("batf/{rel}"))
}

pub fn batf_1d99721(rel: &str) -> PathBuf {
    fixture(&format!("batf_1d99721/{rel}"))
}

/// Running worst-case tracker implementing the section 3 comparison rule.
#[derive(Debug, Default)]
pub struct Errors {
    pub count: usize,
    pub worst_abs: f64,
    pub worst_rel: f64,
    /// Where the worst absolute difference was seen.
    pub worst_where: String,
    /// Descriptions of every comparison that violated the rule, capped.
    pub violations: Vec<String>,
}

impl Errors {
    pub fn new() -> Self {
        Self::default()
    }

    /// Compare one value. `NaN` matches `NaN`, as R's `is.nan` does.
    pub fn add(&mut self, label: &str, rust: f64, r: f64) {
        self.count += 1;
        if rust.is_nan() && r.is_nan() {
            return;
        }
        let abs = (rust - r).abs();
        let rel = if r == 0.0 { 0.0 } else { abs / r.abs() };
        if abs > self.worst_abs {
            self.worst_abs = abs;
            self.worst_where = format!("{label}: rust={rust:?} r={r:?}");
        }
        if rel > self.worst_rel {
            self.worst_rel = rel;
        }
        // The plan's rule, literally: both bounds must hold.
        if !(abs <= TOL_ABS && abs <= TOL_REL * r.abs()) && self.violations.len() < 20 {
            self.violations.push(format!(
                "{label}: rust={rust:?} r={r:?} abs={abs:.3e} rel={rel:.3e}"
            ));
        }
    }

    /// Compare with bit equality, which is what parser scores and the `rmath`
    /// tables require.
    pub fn add_bits(&mut self, label: &str, rust: f64, r: f64) {
        self.count += 1;
        if rust.to_bits() == r.to_bits() {
            return;
        }
        let abs = (rust - r).abs();
        let rel = if r == 0.0 { 0.0 } else { abs / r.abs() };
        self.worst_abs = self.worst_abs.max(abs);
        self.worst_rel = self.worst_rel.max(rel);
        self.worst_where = format!("{label}: rust={rust:?} r={r:?} (not bit-equal)");
        if self.violations.len() < 20 {
            self.violations.push(format!(
                "{label}: rust={rust:?} ({:?}) r={r:?} ({:?})",
                rust.to_bits(),
                r.to_bits()
            ));
        }
    }

    /// Report and panic on any violation. Returns the worst absolute and relative
    /// difference so callers can print them.
    pub fn finish(&self, what: &str) -> (f64, f64) {
        assert!(
            self.violations.is_empty(),
            "{what}: {} of {} values outside the plan's tolerance:\n{}",
            self.violations.len(),
            self.count,
            self.violations.join("\n")
        );
        eprintln!(
            "{what}: {} values, worst abs {:.3e}, worst rel {:.3e} ({})",
            self.count, self.worst_abs, self.worst_rel, self.worst_where
        );
        (self.worst_abs, self.worst_rel)
    }
}

/// Read a two-column TSV as `(f64, f64)` pairs.
pub fn read_pairs(path: &Path) -> Vec<(f64, f64)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let mut f = l.split('\t');
            let a: f64 = f.next().expect("first column").parse().expect("number");
            let b: f64 = f.next().expect("second column").parse().expect("number");
            (a, b)
        })
        .collect()
}

/// Read a tab-separated file with a header into `(header, rows)` of strings.
pub fn read_columns(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split('\t')
        .map(|s| s.to_string())
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect();
    (header, rows)
}

pub fn col<'a>(header: &[String], rows: &'a [Vec<String>], name: &str) -> Vec<&'a str> {
    let i = header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("no column '{name}' in {header:?}"));
    rows.iter().map(|r| r[i].as_str()).collect()
}

/// Load one BATF fixture revision into an annotation plus its sites.
pub fn load_batf(dir: &Path) -> (Annotation, Vec<Site>) {
    let mut chroms = ChromTable::new();
    let sites = read_msites(&dir.join("msites.tsv"), &mut chroms)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    let windows = read_gc_windows(&dir.join("gc_windows.tsv"), &mut chroms)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    let tfbs = read_ranges(&dir.join("tfbs.tsv.gz"), &mut chroms)
        .unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    let gcfreq =
        read_gcfreq(&dir.join("gcfreq.tsv")).unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
    let annotation = Annotation {
        gc_windows: windows,
        motifs: vec![MotifAnnotation {
            name: "BATF".to_string(),
            tfbs: Arc::new(tfbs),
            gcfreq: Arc::new(gcfreq),
        }],
        enhancer: None,
        ignore_strand: true,
    };
    (annotation, sites)
}

/// A small xorshift generator, so the fuzz tests are deterministic without adding
/// a dependency for one `rand::rng()` call.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    /// Uniform in `[0, n)`.
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n.max(1)
    }

    pub fn range_i64(&mut self, lo: i64, hi: i64) -> i64 {
        lo + self.below((hi - lo + 1) as u64) as i64
    }

    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }
}

/// A strand chosen at random.
pub fn any_strand(rng: &mut Rng) -> Strand {
    match rng.below(3) {
        0 => Strand::Plus,
        1 => Strand::Minus,
        _ => Strand::Star,
    }
}

/// A random range on `chr`.
pub fn any_range(rng: &mut Rng, chr: u32, lo: i64, hi: i64) -> methyltfr::Range {
    let start = rng.range_i64(lo, hi);
    let end = start + rng.below(60) as i64;
    methyltfr::Range {
        chr,
        start,
        end,
        strand: any_strand(rng),
    }
}

/// A random GC window on `chr`, with a bin in `1..=5`.
pub fn any_window(rng: &mut Rng, chr: u32, lo: i64, hi: i64) -> methyltfr::GcWindow {
    let r = any_range(rng, chr, lo, hi);
    methyltfr::GcWindow {
        chr: r.chr,
        start: r.start,
        end: r.end,
        strand: r.strand,
        gc_bin: (rng.below(5) + 1) as u8,
    }
}

// --------------------------------------------------------------- parsers ---

use std::io::Write;

/// One record of a `*.expected.tsv` file.
#[derive(Clone, Debug, PartialEq)]
pub struct ExpectedSite {
    pub chr: String,
    pub start: i64,
    pub end: i64,
    pub strand: String,
    pub score: f64,
    pub coverage: f64,
}

/// Read one of the Tier A `*.expected.tsv` files written by
/// `scripts/gen_parser_expectations.R`.
pub fn expected_sites(ty: &str, suffix: &str) -> Vec<ExpectedSite> {
    let path = fixture(&format!("parsers/{ty}.expected{suffix}.tsv"));
    let (header, rows) = read_columns(&path);
    let chr = col(&header, &rows, "chr");
    let start = col(&header, &rows, "start");
    let end = col(&header, &rows, "end");
    let strand = col(&header, &rows, "strand");
    let score = col(&header, &rows, "score");
    let coverage = col(&header, &rows, "coverage");
    (0..rows.len())
        .map(|i| ExpectedSite {
            chr: chr[i].to_string(),
            start: start[i].parse().unwrap(),
            end: end[i].parse().unwrap(),
            strand: strand[i].to_string(),
            score: score[i].parse().unwrap(),
            coverage: coverage[i].parse().unwrap(),
        })
        .collect()
}

/// The bundled example file for a `type`, whose upstream file name does not always
/// match the type.
pub fn parser_example(ty: &str) -> PathBuf {
    let file = match ty {
        "bismarkcytosine" => "bismarkCytosine.tsv.gz",
        "bismarkcov" => "bismarkCov.tsv.gz",
        other => &format!("{other}.tsv.gz"),
    };
    fixture(&format!("parsers/{file}"))
}

/// Compare a methylome against the Tier A records: positions and strand exactly,
/// score and coverage **bit** exactly (`docs/divergences.md` D2).
pub fn compare_to_oracle(ty: &str, suffix: &str, sites: &[methyltfr::Site]) {
    let want = expected_sites(ty, suffix);
    assert_eq!(
        sites.len(),
        want.len(),
        "{ty}: {} records vs {} in the oracle",
        sites.len(),
        want.len()
    );
    for (i, (got, exp)) in sites.iter().zip(want.iter()).enumerate() {
        assert_eq!(got.start, exp.start, "{ty}: record {i} start");
        assert_eq!(got.end, exp.end, "{ty}: record {i} end");
        assert_eq!(
            got.strand.as_char().to_string(),
            exp.strand,
            "{ty}: record {i} strand"
        );
        assert_eq!(
            got.score.to_bits(),
            exp.score.to_bits(),
            "{ty}: record {i} score {:?} vs {:?}",
            got.score,
            exp.score
        );
        assert_eq!(
            got.coverage.to_bits(),
            exp.coverage.to_bits(),
            "{ty}: record {i} coverage {:?} vs {:?}",
            got.coverage,
            exp.coverage
        );
    }
    eprintln!(
        "{ty}{suffix}: {} records, positions and strand equal, score and coverage bit-identical",
        sites.len()
    );
}

/// A scratch directory unique to one test process.
pub fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("methyltfr-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("scratch dir");
    d
}

/// Write a file, gzipping it when the name ends in `.gz`.
pub fn write_file(dir: &Path, name: &str, body: &str) -> PathBuf {
    let p = dir.join(name);
    if name.ends_with(".gz") {
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(&p).expect("create"),
            flate2::Compression::default(),
        );
        enc.write_all(body.as_bytes()).expect("write gzip");
        enc.finish().expect("finish gzip");
    } else {
        std::fs::write(&p, body).expect("write");
    }
    p
}

/// Read one sample through `read_methylome`.
pub fn read(path: &Path, ty: &str, cov_threshold: f64) -> methyltfr::Result<Vec<methyltfr::Site>> {
    let mut chroms = ChromTable::new();
    methyltfr::io::methylome::read_methylome(path, ty, cov_threshold, &mut chroms).map(|m| m.sites)
}

/// The four error cases every parser shares, from `AGENT_PLAN.md` T30..T35.
pub fn assert_shared_errors(dir: &Path, ty: &str, min_cols: usize) {
    // Too few columns. The last row is short, so a parser that reads greedily
    // would silently accept it.
    // Keep the first row intact and truncate the second one.
    let short = truncate_last_row(&sample_rows(ty, min_cols), min_cols - 1);
    let p = write_file(dir, &format!("{ty}_short.tsv"), &short);
    let e = read(&p, ty, 1.0).unwrap_err();
    assert!(
        e.to_string()
            .contains(&format!("at least {min_cols} columns")),
        "{ty}: unexpected message for a short file: {e}"
    );

    // A non-numeric value where a number belongs, reported with file and line.
    // It goes in the *second* row on purpose: `bissnp` skips the first line and
    // `encode` may read the first line as a header, so a fault on line 1 would
    // never reach the parser.
    let broken = corrupt_last_row(&sample_rows(ty, min_cols), 1, "nope");
    let p = write_file(dir, &format!("{ty}_nonnumeric.tsv"), &broken);
    let e = read(&p, ty, 1.0).unwrap_err();
    assert_eq!(e.line(), Some(2), "{ty}: {e}");
    assert!(
        e.to_string().contains("is not a number") || e.to_string().contains("is not an integer"),
        "{ty}: unexpected message: {e}"
    );

    // A strand outside {+, -, *}, which is what IRanges rejects. Also on the second
    // row, for the same reason.
    let strand_col = strand_column(ty);
    let bad = corrupt_last_row(&sample_rows(ty, min_cols), strand_col, ".");
    let p = write_file(dir, &format!("{ty}_strand.tsv"), &bad);
    let e = read(&p, ty, 1.0).unwrap_err();
    assert!(
        e.to_string().contains("strand values must be in"),
        "{ty}: unexpected message for a bad strand: {e}"
    );

    // Plain and gzipped input must produce the same records.
    let plain_path = write_file(dir, &format!("{ty}_plain.tsv"), &sample_rows(ty, min_cols));
    let gz_path = write_file(
        dir,
        &format!("{ty}_plain.tsv.gz"),
        &sample_rows(ty, min_cols),
    );
    let a = read(&plain_path, ty, 1.0).expect("plain");
    let b = read(&gz_path, ty, 1.0).expect("gzip");
    assert_eq!(a.len(), b.len(), "{ty}: gzip and plain disagree on length");
    for (x, y) in a.iter().zip(b.iter()) {
        assert_eq!(x.start, y.start);
        assert_eq!(x.end, y.end);
        assert_eq!(x.score.to_bits(), y.score.to_bits());
        assert_eq!(x.coverage.to_bits(), y.coverage.to_bits());
        assert_eq!(x.strand, y.strand);
    }
}

/// Two well-formed rows for `ty`, enough to reach every column.
pub fn sample_rows(ty: &str, min_cols: usize) -> String {
    let rows: Vec<Vec<&str>> = match ty {
        "epp" => vec![
            vec!["chr1", "10", "11", "1/2", "500", "+"],
            vec!["chr1", "12", "13", "2/4", "1000", "-"],
        ],
        "bissnp" => vec![
            vec!["chr1", "10", "11", "50", "4", "+"],
            vec!["chr1", "12", "13", "100", "8", "-"],
        ],
        "allc" => vec![
            vec!["chr1", "10", "+", "CGT", "1", "2"],
            vec!["chr1", "11", "-", "CGA", "2", "4"],
        ],
        "bismarkcytosine" => vec![
            vec!["chr1", "10", "+", "1", "3", "CG", "CGG"],
            vec!["chr1", "11", "-", "2", "2", "CG", "CGG"],
        ],
        "bismarkcov" => vec![
            vec!["chr1", "10", "11", "50", "1", "1"],
            vec!["chr1", "12", "13", "0", "2", "2"],
        ],
        "encode" => vec![
            vec![
                "chr1", "10", "11", "x", "0", "+", "10", "11", "0,0,0", "5", "20",
            ],
            vec![
                "chr1", "12", "13", "y", "0", "-", "12", "13", "0,0,0", "7", "100",
            ],
        ],
        other => panic!("no sample rows for {other}"),
    };
    let mut out = String::new();
    for r in rows {
        assert!(r.len() >= min_cols, "{ty} sample row is too short");
        out.push_str(&r.join("\t"));
        out.push('\n');
    }
    out
}

/// Which column carries the strand for each format.
pub fn strand_column(ty: &str) -> usize {
    match ty {
        "epp" | "bissnp" | "encode" => 5,
        // bismarkcov synthesises its strand and has no strand column.
        _ => 2,
    }
}

/// Replace one field of the **last** row, leaving every other byte alone.
pub fn corrupt_last_row(body: &str, field: usize, with: &str) -> String {
    let mut lines: Vec<String> = body.lines().map(|l| l.to_string()).collect();
    let last = lines.len() - 1;
    let mut cells: Vec<String> = lines[last].split('\t').map(|s| s.to_string()).collect();
    cells[field] = with.to_string();
    lines[last] = cells.join("\t");
    format!("{}\n", lines.join("\n"))
}

/// Truncate `body` so its last line has `keep` fields.
pub fn truncate_last_row(body: &str, keep: usize) -> String {
    let mut lines: Vec<String> = body.lines().map(|l| l.to_string()).collect();
    let last = lines.len() - 1;
    let fields: Vec<&str> = lines[last].split('\t').collect();
    lines[last] = fields[..keep].join("\t");
    format!("{}\n", lines.join("\n"))
}
