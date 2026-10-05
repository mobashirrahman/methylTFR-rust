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
