//! The data model (`docs/AGENT_PLAN.md` section 3, T11).
//!
//! Two decisions here are load bearing and easy to get wrong:
//!
//! * Coordinates are signed `i64`, 1-based, closed `[start, end]`. Upstream
//!   copies BED-style `start`/`end` straight into `IRanges` without converting,
//!   so a BED record has width 2 and an `allc` / `bismarkCytosine` record has
//!   width 1. Nothing in this crate subtracts or adds 1 "to convert BED".
//! * Chromosome names are interned to `u32` through one shared [`ChromTable`].
//!   Comparing names as strings is what upstream does and is what makes
//!   `chr1` and `1` different chromosomes; interning preserves that, and makes
//!   the comparison a single integer compare in the hot loops.

use crate::error::{Error, Result};
use std::collections::HashMap;
use std::sync::Arc;

/// Chromosome name interning. One instance is shared by every reader in a run,
/// so the `u32` ids in [`Site`], [`Range`] and [`GcWindow`] are comparable
/// across a sample, the GC windows, the binding sites and the enhancer regions.
#[derive(Debug, Default)]
pub struct ChromTable {
    names: Vec<Box<str>>,
    ids: HashMap<Box<str>, u32>,
}

impl ChromTable {
    pub fn new() -> Self {
        Self::default()
    }

    /// Intern `name`, returning its id. The same string always gets the same
    /// id; different strings always get different ids.
    pub fn intern(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.ids.get(name) {
            return id;
        }
        let id = self.names.len() as u32;
        let boxed: Box<str> = name.into();
        self.names.push(boxed.clone());
        self.ids.insert(boxed, id);
        id
    }

    /// The name behind an id, for messages.
    pub fn name(&self, id: u32) -> Option<&str> {
        self.names.get(id as usize).map(|s| &**s)
    }

    /// Number of distinct chromosomes interned so far.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Resolve a user-supplied chromosome name to an id, for callers that look
    /// a chromosome up by name (used by the differential tests).
    pub fn lookup(&self, name: &str) -> Option<u32> {
        self.ids.get(name).copied()
    }
}

/// Strand, restricted to the three values `GRanges` accepts. Upstream does not
/// validate the column itself; `IRanges` does, and it rejects anything else --
/// including `.` -- with "strand values must be in '+' '-' '*'" (confirmed
/// against methylTFR 0.99.9, see `docs/reference-version.md`).
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum Strand {
    #[default]
    Star,
    Plus,
    Minus,
}

impl Strand {
    pub fn parse(s: &str) -> Result<Strand> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            "*" => Ok(Strand::Star),
            _ => Err(Error::run("strand values must be in '+' '-' '*'")),
        }
    }

    pub fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Star => '*',
        }
    }

    /// Strand compatibility (`AGENT_PLAN.md` section 2.5).
    ///
    /// `ignore_strand = true` is the upstream default and makes every pair
    /// compatible. Otherwise two strands match when they are equal or either
    /// side is `*`.
    #[inline]
    pub fn compatible(self, other: Strand, ignore_strand: bool) -> bool {
        ignore_strand || self == other || self == Strand::Star || other == Strand::Star
    }
}

/// A methylation site: one row of a methylome file after parsing and filtering.
#[derive(Clone, Copy, Debug)]
pub struct Site {
    pub chr: u32,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub score: f64,
    /// Kept as `f64` because upstream compares `coverage >= cov_threshold`
    /// without ever converting it to an integer.
    pub coverage: f64,
}

impl Site {
    #[inline]
    pub fn width(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// A plain 1-based closed range: a TFBS, an enhancer region, or a GC window
/// without its bin.
#[derive(Clone, Copy, Debug)]
pub struct Range {
    pub chr: u32,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
}

impl Range {
    #[inline]
    pub fn width(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// A genome-wide GC window, which carries the bin it belongs to.
#[derive(Clone, Copy, Debug)]
pub struct GcWindow {
    pub chr: u32,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub gc_bin: u8,
}

impl GcWindow {
    #[inline]
    pub fn width(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// A `gcfreq` matrix: one row per GC bin, `cols` columns. Row-major, so the
/// inner product of section 2.7 walks one row contiguously.
#[derive(Clone, Debug)]
pub struct GcFreq {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl GcFreq {
    pub fn new(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self> {
        if data.len() != rows * cols {
            return Err(Error::run(format!(
                "gcfreq has {} values but {rows}x{cols} = {} were expected",
                data.len(),
                rows * cols
            )));
        }
        Ok(Self { rows, cols, data })
    }

    #[inline]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }
}

/// One sample's methylation sites, in input order.
#[derive(Clone, Debug, Default)]
pub struct Methylome {
    pub sites: Vec<Site>,
}

impl Methylome {
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }
}

/// One motif's annotation: its binding sites in file order (original widths,
/// never resized) and its GC frequency matrix.
#[derive(Clone, Debug)]
pub struct MotifAnnotation {
    pub name: String,
    /// Shared, not cloned, when two manifest rows point at the same TFBS file:
    /// a 268 717-row BATF binding-site list is not something to copy per motif.
    pub tfbs: Arc<Vec<Range>>,
    pub gcfreq: Arc<GcFreq>,
}

/// Everything a run needs besides the samples themselves.
#[derive(Clone, Debug)]
pub struct Annotation {
    pub gc_windows: Vec<GcWindow>,
    /// Manifest row order is motif order (`AGENT_PLAN.md` section 2.8).
    pub motifs: Vec<MotifAnnotation>,
    pub enhancer: Option<Vec<Range>>,
    pub ignore_strand: bool,
}

/// Parse a chromosome name in a context that can report where it came from.
///
/// `GRanges` accepts only `+`, `-` and `*`; this keeps the same message so the
/// two implementations fail the same way on the same file.
pub(crate) fn parse_strand(
    value: &str,
    path: impl AsRef<std::path::Path>,
    line: usize,
) -> Result<Strand> {
    match value {
        "+" => Ok(Strand::Plus),
        "-" => Ok(Strand::Minus),
        "*" => Ok(Strand::Star),
        other => Err(Error::parse(
            path,
            line,
            format!("strand values must be in '+' '-' '*', found '{other}'"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interning_is_stable_and_distinct() {
        let mut t = ChromTable::new();
        let a = t.intern("chr1");
        let b = t.intern("1");
        assert_ne!(a, b, "chr1 and 1 are different chromosomes upstream too");
        assert_eq!(t.intern("chr1"), a);
        assert_eq!(t.name(a).unwrap(), "chr1");
        assert_eq!(t.name(b).unwrap(), "1");
        assert_eq!(t.len(), 2);
    }

    #[test]
    fn strand_rules() {
        assert!(Strand::Plus.compatible(Strand::Minus, true));
        assert!(!Strand::Plus.compatible(Strand::Minus, false));
        assert!(Strand::Plus.compatible(Strand::Star, false));
        assert!(Strand::Star.compatible(Strand::Star, false));
        assert!(Strand::Plus.compatible(Strand::Plus, false));
    }

    #[test]
    fn strand_parse_rejects_dot() {
        let e = Strand::parse(".").unwrap_err();
        assert!(e.to_string().contains("strand values must be in"));
        assert_eq!(Strand::parse("*").unwrap(), Strand::Star);
    }

    #[test]
    fn gcfreq_rejects_wrong_length() {
        assert!(GcFreq::new(2, 3, vec![0.0; 5]).is_err());
        assert!(GcFreq::new(2, 3, vec![0.0; 6]).is_ok());
    }
}
