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

/// The methylome as the hot loops want it: grouped by chromosome, sorted by
/// start, in struct-of-arrays form.
///
/// Two reasons, both from `AGENT_PLAN.md` M7 T72:
///
/// * **Sorted once.** A methylation file is usually already in coordinate order,
///   but it is not guaranteed to be, and every consumer of the sample would
///   otherwise have to establish that itself. Sorting here is done once per
///   sample.
/// * **Struct of arrays.** `starts` is the array a sweep walks; `scores` is the
///   array a hit reads. Keeping them apart means the sweep's inner loop touches
///   one contiguous `i64` array rather than striding over 32-byte structs, and it
///   drops `coverage`, which nothing below this point still needs.
///
/// Consequence to be explicit about: hits are now visited in sorted order rather
/// than input order, so a sum over many hits can differ from the same sum taken in
/// input order in its last bit. Every committed fixture is already sorted, so
/// their output is unchanged; see `docs/divergences.md` D4 for the measured size
/// of the difference on unsorted input.
#[derive(Clone, Debug)]
pub struct SortedMethylome {
    starts: Vec<i64>,
    ends: Vec<i64>,
    strands: Vec<u8>,
    scores: Vec<f64>,
    /// One entry per chromosome, in ascending chromosome id, each a half-open
    /// range into the arrays above.
    groups: Vec<ChromGroup>,
    len: usize,
}

/// The sites of one chromosome, borrowed.
#[derive(Clone, Copy, Debug)]
pub struct ChromGroup {
    pub chr: u32,
    pub begin: usize,
    pub end: usize,
}

/// Borrowed struct-of-arrays view of one chromosome's sites.
#[derive(Clone, Copy, Debug)]
pub struct ChromSites<'a> {
    pub chr: u32,
    pub starts: &'a [i64],
    pub ends: &'a [i64],
    pub strands: &'a [u8],
    pub scores: &'a [f64],
}

impl ChromSites<'_> {
    #[inline]
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    #[inline]
    pub fn strand(&self, i: usize) -> Strand {
        match self.strands[i] {
            1 => Strand::Plus,
            2 => Strand::Minus,
            _ => Strand::Star,
        }
    }
}

impl SortedMethylome {
    /// Group and sort a methylome. Stable, so equal (start, end) pairs keep their
    /// input order and a pre-sorted input comes out byte-identical.
    pub fn new(sites: &[Site]) -> Self {
        let mut order: Vec<u32> = (0..sites.len()).map(|i| i as u32).collect();
        order.sort_by_key(|&i| {
            let s = &sites[i as usize];
            (s.chr, s.start, s.end)
        });

        let len = sites.len();
        let mut starts = Vec::with_capacity(len);
        let mut ends = Vec::with_capacity(len);
        let mut strands = Vec::with_capacity(len);
        let mut scores = Vec::with_capacity(len);
        for &i in &order {
            let s = &sites[i as usize];
            starts.push(s.start);
            ends.push(s.end);
            strands.push(match s.strand {
                Strand::Plus => 1u8,
                Strand::Minus => 2,
                Strand::Star => 0,
            });
            scores.push(s.score);
        }

        // One group per run of equal chromosome ids in the sorted order.
        let mut groups: Vec<ChromGroup> = Vec::new();
        let mut run_start = 0usize;
        // From the *sorted* order, not from `sites[0]`: the input may start on any
        // chromosome.
        let mut run_chr = order.first().map(|&i| sites[i as usize].chr);
        for (i, &idx) in order.iter().enumerate() {
            let c = sites[idx as usize].chr;
            if Some(c) != run_chr {
                groups.push(ChromGroup {
                    chr: run_chr.unwrap_or(c),
                    begin: run_start,
                    end: i,
                });
                run_start = i;
                run_chr = Some(c);
            }
        }
        if len > 0 {
            groups.push(ChromGroup {
                chr: run_chr.unwrap_or_default(),
                begin: run_start,
                end: len,
            });
        }

        Self {
            starts,
            ends,
            strands,
            scores,
            groups,
            len,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// One entry per chromosome, ascending by chromosome id.
    pub fn groups(&self) -> &[ChromGroup] {
        &self.groups
    }

    /// The group for `chr`, if the methylome has any site on it.
    pub fn group(&self, chr: u32) -> Option<&ChromGroup> {
        self.groups
            .binary_search_by_key(&chr, |g| g.chr)
            .ok()
            .map(|i| &self.groups[i])
    }

    /// Borrow one chromosome's sites.
    pub fn chrom(&self, g: &ChromGroup) -> ChromSites<'_> {
        ChromSites {
            chr: g.chr,
            starts: &self.starts[g.begin..g.end],
            ends: &self.ends[g.begin..g.end],
            strands: &self.strands[g.begin..g.end],
            scores: &self.scores[g.begin..g.end],
        }
    }

    /// The sites of `chr`, or an empty view when there are none.
    pub fn chrom_by_id(&self, chr: u32) -> ChromSites<'_> {
        match self.group(chr) {
            Some(g) => self.chrom(g),
            None => ChromSites {
                chr,
                starts: &[],
                ends: &[],
                strands: &[],
                scores: &[],
            },
        }
    }

    /// The group order's chromosome ids, for iterating the annotation side too.
    pub fn chromosome_ids(&self) -> Vec<u32> {
        self.groups.iter().map(|g| g.chr).collect()
    }
}

#[cfg(test)]
mod sorted_tests {
    use super::*;

    fn site(chr: u32, start: i64, end: i64, score: f64) -> Site {
        Site {
            chr,
            start,
            end,
            strand: Strand::Star,
            score,
            coverage: 1.0,
        }
    }

    #[test]
    fn groups_and_sorts_by_chromosome_then_start() {
        let sites = vec![
            site(1, 50, 51, 0.5),
            site(0, 10, 11, 0.1),
            site(1, 10, 11, 0.2),
            site(0, 20, 21, 0.3),
        ];
        let m = SortedMethylome::new(&sites);
        assert_eq!(m.len(), 4);
        assert_eq!(
            m.groups().iter().map(|g| g.chr).collect::<Vec<_>>(),
            vec![0, 1]
        );
        let c0 = m.chrom_by_id(0);
        assert_eq!(c0.len(), 2);
        assert_eq!(c0.starts, &[10, 20]);
        assert_eq!(c0.scores, &[0.1, 0.3]);
        let c1 = m.chrom_by_id(1);
        assert_eq!(c1.starts, &[10, 50]);
        assert_eq!(c1.scores, &[0.2, 0.5]);
    }

    #[test]
    fn a_missing_chromosome_is_an_empty_view_not_a_panic() {
        let m = SortedMethylome::new(&[site(0, 1, 2, 1.0)]);
        assert!(m.chrom_by_id(7).is_empty());
        assert!(m.group(7).is_none());
        assert_eq!(m.group(0).map(|g| g.chr), Some(0));
    }

    #[test]
    fn an_empty_methylome_has_no_groups() {
        let m = SortedMethylome::new(&[]);
        assert!(m.is_empty());
        assert!(m.groups().is_empty());
        assert_eq!(m.chromosome_ids(), Vec::<u32>::new());
    }

    #[test]
    fn an_already_sorted_input_keeps_its_order() {
        let sites: Vec<Site> = (0..100)
            .map(|i| site(0, i * 10, i * 10 + 1, i as f64))
            .collect();
        let m = SortedMethylome::new(&sites);
        let chrom = m.chrom_by_id(0);
        assert_eq!(chrom.len(), sites.len());
        for (i, s) in sites.iter().enumerate() {
            assert_eq!(chrom.starts[i], s.start);
            assert_eq!(chrom.scores[i].to_bits(), s.score.to_bits());
        }
    }

    #[test]
    fn equal_coordinates_keep_input_order() {
        // Three identical coordinates with different scores: the order matters for
        // the f64 sum, so stability is a requirement, not a nicety.
        let sites = vec![site(0, 5, 6, 0.1), site(0, 5, 6, 0.2), site(0, 5, 6, 0.3)];
        let m = SortedMethylome::new(&sites);
        let c = m.chrom_by_id(0);
        assert_eq!(c.scores, &[0.1, 0.2, 0.3]);
    }

    #[test]
    fn strands_round_trip_through_the_byte_encoding() {
        let sites = vec![
            Site {
                chr: 0,
                start: 1,
                end: 1,
                strand: Strand::Plus,
                score: 0.0,
                coverage: 0.0,
            },
            Site {
                chr: 0,
                start: 2,
                end: 2,
                strand: Strand::Minus,
                score: 0.0,
                coverage: 0.0,
            },
            Site {
                chr: 0,
                start: 3,
                end: 3,
                strand: Strand::Star,
                score: 0.0,
                coverage: 0.0,
            },
        ];
        let m = SortedMethylome::new(&sites);
        let c = m.chrom_by_id(0);
        assert_eq!(c.strand(0), Strand::Plus);
        assert_eq!(c.strand(1), Strand::Minus);
        assert_eq!(c.strand(2), Strand::Star);
    }
}
