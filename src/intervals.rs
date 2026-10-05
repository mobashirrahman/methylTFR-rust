//! Interval queries (`docs/AGENT_PLAN.md` section 2.4, T12).
//!
//! Three queries are needed, and all three reduce to binary searches over a
//! per-chromosome start-sorted index -- no interval tree:
//!
//! * [`StartIndex::any_overlaps`] -- `findOverlaps(type = "any")`, used for sites
//!   against GC windows.
//! * [`StartIndex::within_uniform`] -- `findOverlaps(type = "within")` against
//!   binding sites that all share one width after the centre resize.
//! * [`StartIndex::overlaps_any`] -- the boolean form, for enhancer subsetting.
//!
//! The key observation for the first two: the set of candidate ranges is a
//! contiguous slice of the start-sorted list.
//!
//! * A range overlaps `[s, e]` when `range.start <= e && range.end >= s`. The
//!   upper bound is a plain `start <= e`. The lower bound needs the range's
//!   width, and since every range in a chromosome is at most `max_width` wide,
//!   `range.start >= s - max_width + 1` is a valid (never too narrow) bound.
//! * Every resized TFBS has the same width `W`, so containing `[s, e]` is
//!   `e - W + 1 <= range.start <= s` exactly, again a contiguous slice.

use crate::model::{GcWindow, Range, Site, Strand};

/// The three fields every interval query needs, so the index can be built
/// generically over [`GcWindow`], [`Range`] and [`Site`] without a closure per
/// element.
pub trait Interval {
    fn chr_id(&self) -> u32;
    fn start_pos(&self) -> i64;
    fn end_pos(&self) -> i64;
    fn strand_of(&self) -> Strand;
    /// `end - start + 1`, i.e. the closed-interval width.
    fn width_pos(&self) -> i64 {
        self.end_pos() - self.start_pos() + 1
    }
}

impl Interval for GcWindow {
    fn chr_id(&self) -> u32 {
        self.chr
    }
    fn start_pos(&self) -> i64 {
        self.start
    }
    fn end_pos(&self) -> i64 {
        self.end
    }
    fn strand_of(&self) -> Strand {
        self.strand
    }
}

impl Interval for Range {
    fn chr_id(&self) -> u32 {
        self.chr
    }
    fn start_pos(&self) -> i64 {
        self.start
    }
    fn end_pos(&self) -> i64 {
        self.end
    }
    fn strand_of(&self) -> Strand {
        self.strand
    }
}

/// Plain tuples are the ergonomic query type for callers that already have the
/// four fields (the enhancer regions in `AGENT_PLAN.md` section 2.6, and the
/// differential tests).
impl Interval for (u32, i64, i64, Strand) {
    fn chr_id(&self) -> u32 {
        self.0
    }
    fn start_pos(&self) -> i64 {
        self.1
    }
    fn end_pos(&self) -> i64 {
        self.2
    }
    fn strand_of(&self) -> Strand {
        self.3
    }
}

impl Interval for Site {
    fn chr_id(&self) -> u32 {
        self.chr
    }
    fn start_pos(&self) -> i64 {
        self.start
    }
    fn end_pos(&self) -> i64 {
        self.end
    }
    fn strand_of(&self) -> Strand {
        self.strand
    }
}

/// A query range together with the strand mode it should be judged under.
///
/// Bundling these five values into one argument is not cosmetic: the three query
/// entry points take them as a group, and passing them separately is how a
/// `within_uniform` call ends up with the site end in the start slot.
#[derive(Clone, Copy, Debug)]
pub struct Query {
    pub chr: u32,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub ignore_strand: bool,
}

impl Query {
    pub fn new(chr: u32, start: i64, end: i64, strand: Strand, ignore_strand: bool) -> Self {
        Self {
            chr,
            start,
            end,
            strand,
            ignore_strand,
        }
    }

    pub fn from_site(site: &Site, ignore_strand: bool) -> Self {
        Self {
            chr: site.chr,
            start: site.start,
            end: site.end,
            strand: site.strand,
            ignore_strand,
        }
    }

    pub fn from_tuple(t: (u32, i64, i64, Strand), ignore_strand: bool) -> Self {
        Self::new(t.0, t.1, t.2, t.3, ignore_strand)
    }
}

impl Interval for Query {
    fn chr_id(&self) -> u32 {
        self.chr
    }
    fn start_pos(&self) -> i64 {
        self.start
    }
    fn end_pos(&self) -> i64 {
        self.end
    }
    fn strand_of(&self) -> Strand {
        self.strand
    }
}

/// Per-chromosome ranges sorted by start, remembering each one's original index,
/// its start, and the maximum width in its chromosome.
///
/// Ties on `start` keep the original relative order, which is what makes the hit
/// order of `AGENT_PLAN.md` section 2.3 ("sites in input order, then windows by
/// ascending start") reproducible.
#[derive(Debug, Default)]
pub struct StartIndex {
    /// `order[c][k]` is the original index of the k-th range on chromosome `c`,
    /// by ascending start.
    order: Vec<Vec<u32>>,
    /// `starts[c][k]` is that range's start. Kept alongside `order` so that the
    /// binary searches never have to gather it.
    starts: Vec<Vec<i64>>,
    /// `offsets[c]..offsets[c + 1]` is the flat range covering `order[c]`.
    offsets: Vec<usize>,
    /// Widest range on each chromosome.
    max_width: Vec<i64>,
    /// The chromosome id behind each group, in group order.
    group_chr: Vec<u32>,
    /// Chromosome id -> group index, or `NO_CHR` when the chromosome has no
    /// ranges at all. Groups are created in first-appearance order, which is not
    /// the same as ascending chromosome id, so this lookup is what stops a query
    /// on chr 1 from being answered with chr 0's ranges.
    slot_of: Vec<u32>,
}

/// Sentinel in [`StartIndex::slot_of`] for a chromosome with no ranges.
const NO_CHR: u32 = u32::MAX;

impl StartIndex {
    /// Build an index over `ranges`, bucketed by chromosome id.
    ///
    /// Bucketing rather than splitting on contiguous runs matters: a GC window
    /// file interleaves chromosomes freely, and a run-based split would then put
    /// most windows in groups nothing can reach. Each bucket is then sorted by
    /// start, stably, so equal starts keep the input order the plan's hit order
    /// depends on.
    pub fn build<T: Interval>(ranges: &[T]) -> Self {
        let max_chr = ranges.iter().map(|r| r.chr_id()).max().unwrap_or(0);
        let n_chr = max_chr as usize + 1;
        let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); n_chr];
        for (i, r) in ranges.iter().enumerate() {
            buckets[r.chr_id() as usize].push(i as u32);
        }
        buckets.retain(|b| !b.is_empty());
        let max_width: Vec<i64> = buckets
            .iter()
            .map(|b| {
                b.iter()
                    .map(|&i| ranges[i as usize].width_pos())
                    .max()
                    .unwrap_or(0)
            })
            .collect();
        let mut starts: Vec<Vec<i64>> = Vec::with_capacity(buckets.len());
        for b in buckets.iter_mut() {
            let mut keyed: Vec<(i64, u32)> = b
                .iter()
                .map(|&i| (ranges[i as usize].start_pos(), i))
                .collect();
            // `sort_by_key` is stable, so equal starts keep input order.
            keyed.sort_by_key(|(s, _)| *s);
            *b = keyed.iter().map(|&(_, i)| i).collect();
            starts.push(keyed.iter().map(|&(s, _)| s).collect());
        }

        let mut offsets = Vec::with_capacity(buckets.len() + 1);
        offsets.push(0usize);
        for b in &buckets {
            offsets.push(offsets.last().copied().unwrap() + b.len());
        }
        let mut slot_of = vec![NO_CHR; n_chr];
        let mut group_chr = Vec::with_capacity(buckets.len());
        for (slot, b) in buckets.iter().enumerate() {
            // The bucket's chromosome id is the first entry's, which is exact
            // because a bucket only holds ranges of one chromosome.
            let c = ranges[b[0] as usize].chr_id();
            slot_of[c as usize] = slot as u32;
            group_chr.push(c);
        }
        Self {
            order: buckets,
            starts,
            offsets,
            max_width,
            group_chr,
            slot_of,
        }
    }

    /// Number of chromosomes with at least one indexed range.
    pub fn chromosomes(&self) -> usize {
        self.order.len()
    }

    /// Total number of indexed ranges.
    pub fn len(&self) -> usize {
        self.offsets.last().copied().unwrap_or(0)
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The group holding `chr`, or `None` when it has no ranges.
    ///
    /// Chromosome ids and group indices are different things: groups are built in
    /// first-appearance order, so a file whose first range is on chr2 and second
    /// on chr1 would otherwise answer a chr1 query from chr2's ranges.
    #[inline]
    fn chr_slot(&self, chr: u32) -> Option<usize> {
        let slot = *self.slot_of.get(chr as usize)?;
        if slot == NO_CHR {
            None
        } else {
            Some(slot as usize)
        }
    }

    /// The chromosome id of a group, for tests and diagnostics.
    pub fn group_chr(&self, slot: usize) -> Option<u32> {
        self.group_chr.get(slot).copied()
    }

    /// `findOverlaps(site, ranges, type = "any")`: every range overlapping
    /// `[start, end]`, in ascending start order.
    pub fn any_overlaps<'a, T: Interval>(
        &'a self,
        ranges: &'a [T],
        q: Query,
    ) -> OverlapIter<'a, T> {
        let Some(c) = self.chr_slot(q.chr) else {
            return OverlapIter::empty(ranges);
        };
        let starts = &self.starts[c];
        let lo = lo_ge(starts, q.start, self.max_width[c]);
        let hi = hi_le(starts, q.end);
        OverlapIter {
            cursor: lo,
            hi,
            start_limit: q.end,
            site_start: q.start,
            site_chr: q.chr,
            site_strand: q.strand,
            ignore_strand: q.ignore_strand,
            slice: &self.order[c],
            ranges,
        }
    }

    /// True when at least one indexed range overlaps `[start, end]`.
    ///
    /// A separate method rather than `any_overlaps(..).next().is_some()` because
    /// the boolean form never builds an iterator and can stop at the first hit,
    /// which is what reducing a genome-scale GC window set needs.
    pub fn overlaps_any<T: Interval>(&self, ranges: &[T], q: Query) -> bool {
        let Some(c) = self.chr_slot(q.chr) else {
            return false;
        };
        let starts = &self.starts[c];
        let lo = lo_ge(starts, q.start, self.max_width[c]);
        let hi = hi_le(starts, q.end);
        (lo..hi).any(|k| {
            let r = &ranges[self.order[c][k] as usize];
            r.chr_id() == q.chr
                && r.end_pos() >= q.start
                && r.strand_of().compatible(q.strand, q.ignore_strand)
        })
    }

    /// `findOverlaps(site, tfbs, type = "within")` when every range has width
    /// `width`: the ranges containing `[start, end]` are exactly those with
    /// `end - width + 1 <= range.start <= start`.
    pub fn within_uniform<'a, T: Interval>(
        &'a self,
        ranges: &'a [T],
        q: Query,
        width: i64,
    ) -> WithinIter<'a, T> {
        let Some(c) = self.chr_slot(q.chr) else {
            return WithinIter::empty(ranges);
        };
        let starts = &self.starts[c];
        let lo = lo_ge(starts, q.end, width);
        let hi = hi_le(starts, q.start);
        WithinIter {
            cursor: lo,
            hi,
            site_start: q.start,
            site_end: q.end,
            site_chr: q.chr,
            site_strand: q.strand,
            ignore_strand: q.ignore_strand,
            slice: &self.order[c],
            ranges,
        }
    }
}

/// The first position in an ascending slice whose value is `>= needle`.
#[inline]
fn lower_bound(starts: &[i64], needle: i64) -> usize {
    starts.partition_point(|&s| s < needle)
}

/// Upper end of the slice satisfying `range.start <= end`: every entry up to the
/// first one starting after `end`.
///
/// `end + 1` is the comparison threshold and is not representable at `i64::MAX`.
/// Saturating there would drop a range that starts exactly at `i64::MAX` and
/// therefore does overlap, so overflow widens the bound to the whole slice
/// instead. The `range.end >= start` test culls the extras.
#[inline]
fn hi_le(starts: &[i64], end: i64) -> usize {
    match end.checked_add(1) {
        Some(needle) => lower_bound(starts, needle),
        None => starts.len(),
    }
}

/// Lower end of the slice that can satisfy `range.end >= start` given that no
/// range is wider than `max_width`: every entry with
/// `range.start >= start - max_width + 1`.
///
/// `i64::MIN` coordinates make `start - max_width + 1` unrepresentable. Widening
/// the bound there is always safe: it only adds candidates, and each is filtered
/// by the real `range.end >= start` test.
#[inline]
fn lo_ge(starts: &[i64], start: i64, max_width: i64) -> usize {
    let needle = start
        .checked_sub(max_width)
        .and_then(|v| v.checked_add(1))
        .unwrap_or(i64::MIN);
    lower_bound(starts, needle)
}

/// Lazy iterator over the ranges overlapping a query range, in ascending start
/// order. Allocating nothing per query.
pub struct OverlapIter<'a, T> {
    cursor: usize,
    hi: usize,
    /// A range starting after this cannot overlap, and starts ascend, so once
    /// one does the iteration is over.
    start_limit: i64,
    site_start: i64,
    site_chr: u32,
    site_strand: Strand,
    ignore_strand: bool,
    slice: &'a [u32],
    ranges: &'a [T],
}

impl<'a, T: Interval> OverlapIter<'a, T> {
    fn empty(ranges: &'a [T]) -> Self {
        Self {
            cursor: 0,
            hi: 0,
            start_limit: i64::MIN,
            site_start: i64::MAX,
            site_chr: u32::MAX,
            site_strand: Strand::Star,
            ignore_strand: true,
            slice: &[],
            ranges,
        }
    }
}

impl<T: Interval> Iterator for OverlapIter<'_, T> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.cursor < self.hi {
            let idx = self.slice[self.cursor] as usize;
            self.cursor += 1;
            let r = &self.ranges[idx];
            if r.start_pos() > self.start_limit {
                return None;
            }
            if r.chr_id() != self.site_chr || r.end_pos() < self.site_start {
                continue;
            }
            if !r
                .strand_of()
                .compatible(self.site_strand, self.ignore_strand)
            {
                continue;
            }
            return Some(idx);
        }
        None
    }
}

/// Lazy iterator over the uniform-width ranges containing a query range.
pub struct WithinIter<'a, T> {
    cursor: usize,
    hi: usize,
    site_start: i64,
    site_end: i64,
    site_chr: u32,
    site_strand: Strand,
    ignore_strand: bool,
    slice: &'a [u32],
    ranges: &'a [T],
}

impl<'a, T: Interval> WithinIter<'a, T> {
    fn empty(ranges: &'a [T]) -> Self {
        Self {
            cursor: 0,
            hi: 0,
            site_start: i64::MAX,
            site_end: i64::MIN,
            site_chr: u32::MAX,
            site_strand: Strand::Star,
            ignore_strand: true,
            slice: &[],
            ranges,
        }
    }
}

impl<T: Interval> Iterator for WithinIter<'_, T> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.cursor < self.hi {
            let idx = self.slice[self.cursor] as usize;
            self.cursor += 1;
            let r = &self.ranges[idx];
            if r.chr_id() != self.site_chr || r.start_pos() > self.site_start {
                continue;
            }
            if r.end_pos() < self.site_end {
                continue;
            }
            if !r
                .strand_of()
                .compatible(self.site_strand, self.ignore_strand)
            {
                continue;
            }
            return Some(idx);
        }
        None
    }
}

/// `subsetByOverlaps(query, ranges, type = "any")` as a boolean mask **over
/// `ranges`**: entry `i` is true when `ranges[i]` overlaps at least one query
/// region.
///
/// The mask is indexed by `ranges`, not by `query`. That direction matters: the
/// enhancer set is usually far smaller than the set being reduced, so it is the
/// one worth indexing, and a mask indexed by `query` would be the wrong length.
///
/// This is what `AGENT_PLAN.md` section 2.6 needs, twice: to reduce the GC
/// windows and to filter the resized TFBS.
pub fn subset_by_overlaps<T: Interval, Q: Interval>(
    ranges: &[T],
    query: &[Q],
    ignore_strand: bool,
) -> Vec<bool> {
    // Index the query set -- the enhancer regions -- and probe it once per
    // range being reduced.
    let index = StartIndex::build(query);
    ranges
        .iter()
        .map(|r| {
            index
                .any_overlaps(
                    query,
                    Query::new(
                        r.chr_id(),
                        r.start_pos(),
                        r.end_pos(),
                        r.strand_of(),
                        ignore_strand,
                    ),
                )
                .next()
                .is_some()
        })
        .collect()
}

/// Resize every range to `width` about its centre.
///
/// `AGENT_PLAN.md` section 2.4, step 2, with `width_of_range` the range's own
/// width: `new_start = start + (width_of_range - width) %/% 2` and
/// `new_end = new_start + width - 1`. `new_start` may be zero or negative and is
/// kept, because upstream does not filter on it.
///
/// The halving is of the *difference*, not of each width separately:
/// `(100 - 11) %/% 2` is 44, while `100 %/% 2 - 11 %/% 2` is 45. They differ
/// whenever the two widths have different parity, so the difference is halved
/// whole.
pub fn resize_center(ranges: &[Range], width: i64) -> Vec<Range> {
    debug_assert!(width > 0, "resize width must be positive");
    ranges
        .iter()
        .map(|r| {
            let new_start = r.start + (r.width() - width).div_euclid(2);
            Range {
                chr: r.chr,
                start: new_start,
                end: new_start + width - 1,
                strand: r.strand,
            }
        })
        .collect()
}

/// The midpoint of an already-resized range:
/// `round(end + (start - end) / 2)`, computed in `f64` and rounded half-to-even.
#[inline]
pub fn midpoint(start: i64, end: i64) -> i64 {
    crate::rmath::midpoint(start, end)
}

/// `W` in `AGENT_PLAN.md` section 2.4, step 1: the width of the **first** TFBS,
/// plus 130, taken before any filtering.
pub fn motif_width(first: &Range) -> i64 {
    first.width() + 130
}

/// [`StartIndex::any_overlaps`] for a [`Site`] rather than its four fields, which
/// is the shape the GC assignment in `src/gc.rs` wants.
pub fn site_windows<'a, T: Interval>(
    index: &'a StartIndex,
    ranges: &'a [T],
    site: &Site,
    ignore_strand: bool,
) -> OverlapIter<'a, T> {
    index.any_overlaps(ranges, Query::from_site(site, ignore_strand))
}

/// Brute-force `type = "any"` overlap set, for the differential test in
/// `tests/interval_semantics.rs`. Deliberately quadratic.
pub fn brute_force_overlaps<T: Interval>(
    ranges: &[T],
    query: &(u32, i64, i64, Strand),
    ignore_strand: bool,
) -> Vec<usize> {
    let (chr, s, e, strand) = *query;
    let mut out = Vec::new();
    for (i, r) in ranges.iter().enumerate() {
        if r.chr_id() == chr
            && r.start_pos() <= e
            && r.end_pos() >= s
            && r.strand_of().compatible(strand, ignore_strand)
        {
            out.push(i);
        }
    }
    out
}

/// Brute-force `type = "within"` overlap set for uniform-width ranges.
pub fn brute_force_within<T: Interval>(
    ranges: &[T],
    query: &(u32, i64, i64, Strand),
    ignore_strand: bool,
) -> Vec<usize> {
    let (chr, s, e, strand) = *query;
    let mut out = Vec::new();
    for (i, r) in ranges.iter().enumerate() {
        if r.chr_id() == chr
            && r.start_pos() <= s
            && r.end_pos() >= e
            && r.strand_of().compatible(strand, ignore_strand)
        {
            out.push(i);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn w(chr: u32, start: i64, end: i64, strand: Strand, bin: u8) -> GcWindow {
        GcWindow {
            chr,
            start,
            end,
            strand,
            gc_bin: bin,
        }
    }

    fn r(chr: u32, start: i64, end: i64, strand: Strand) -> Range {
        Range {
            chr,
            start,
            end,
            strand,
        }
    }

    fn hits(ranges: &[GcWindow], chr: u32, s: i64, e: i64, st: Strand, ig: bool) -> Vec<usize> {
        let idx = StartIndex::build(ranges);
        idx.any_overlaps(ranges, Query::new(chr, s, e, st, ig))
            .collect()
    }

    #[test]
    fn touching_ends_overlap() {
        let g = vec![w(0, 10, 20, Strand::Star, 1)];
        assert_eq!(hits(&g, 0, 20, 25, Strand::Plus, true), vec![0]);
    }

    #[test]
    fn off_by_one_misses_on_both_sides() {
        let g = vec![w(0, 10, 20, Strand::Star, 1)];
        assert!(hits(&g, 0, 21, 25, Strand::Plus, true).is_empty());
        assert!(hits(&g, 0, 5, 9, Strand::Plus, true).is_empty());
    }

    #[test]
    fn width_two_site_straddling_two_windows_gives_two_hits() {
        let g = vec![w(0, 10, 20, Strand::Star, 1), w(0, 21, 30, Strand::Star, 2)];
        assert_eq!(hits(&g, 0, 20, 21, Strand::Plus, true), vec![0, 1]);
    }

    #[test]
    fn wide_window_is_not_missed_by_the_max_width_bound() {
        // A narrow window next to a wide one: the wide one must stay reachable
        // from a site that is far from its start.
        let g = vec![
            w(0, 1000, 1099, Strand::Star, 1),
            w(0, 10, 20, Strand::Star, 1),
        ];
        let idx = StartIndex::build(&g);
        // Query at 1000..1001: the wide window is found even though its start is
        // far below `start - max_width + 1` only if max_width is 100.
        assert_eq!(
            idx.any_overlaps(&g, Query::new(0, 1000, 1001, Strand::Plus, true))
                .collect::<Vec<_>>(),
            vec![0]
        );
    }

    #[test]
    fn different_chromosomes_never_overlap() {
        let g = vec![w(0, 10, 20, Strand::Star, 1), w(1, 10, 20, Strand::Star, 1)];
        assert_eq!(hits(&g, 1, 15, 16, Strand::Plus, true).len(), 1);
        assert!(hits(&g, 7, 15, 16, Strand::Plus, true).is_empty());
        // The chromosome-id -> group lookup, not the group index, picks the bucket.
    }

    #[test]
    fn hits_come_back_in_ascending_start_order() {
        let g = vec![
            w(0, 100, 120, Strand::Star, 1),
            w(0, 10, 200, Strand::Star, 1),
            w(0, 50, 60, Strand::Star, 1),
        ];
        // The index order is 1 (start 10), 2 (start 50), 0 (start 100).
        assert_eq!(hits(&g, 0, 110, 111, Strand::Plus, true), vec![1, 0]);
    }

    #[test]
    fn resize_center_matches_upstream() {
        let tfbs = vec![r(0, 47430, 47840, Strand::Plus)];
        let resized = resize_center(&tfbs, motif_width(&tfbs[0]));
        assert_eq!(resized[0].start, 47365);
        assert_eq!(resized[0].end, 47905);
        assert_eq!(resized[0].width(), 541);
        assert_eq!(midpoint(resized[0].start, resized[0].end), 47635);
    }

    #[test]
    fn resize_center_handles_mixed_and_negative() {
        // Mixed widths, odd and even, so `(width - target) %/% 2` has both
        // parities of difference.
        let tfbs = vec![r(0, 100, 199, Strand::Plus), r(0, 10, 12, Strand::Minus)];
        let resized = resize_center(&tfbs, 11);
        assert_eq!(resized[0].width(), 11);
        assert_eq!(resized[0].start, 100 + (100i64 - 11).div_euclid(2));
        assert_eq!(resized[0].start, 144, "(100-11)%/%2 is 44, not 45");
        assert_eq!(resized[1].start, 10 + (3i64 - 11).div_euclid(2));
        assert_eq!(resized[1].start, 6);
        let tfbs = vec![r(0, 5, 10, Strand::Star)];
        let resized = resize_center(&tfbs, 100);
        assert!(resized[0].start <= 0, "negative new_start must be kept");
        assert_eq!(resized[0].width(), 100);
    }

    #[test]
    fn resize_center_matches_upstream_exactly_for_411_to_541() {
        // (411 - 541) %/% 2 == -65 exactly, which is the BATF fixture case.
        let tfbs = vec![r(0, 47430, 47840, Strand::Plus)];
        let width = motif_width(&tfbs[0]);
        assert_eq!(width, 541);
        assert_eq!(resize_center(&tfbs, width)[0].start, 47430 - 65);
    }

    #[test]
    fn within_uniform_containment() {
        let g = vec![r(0, 100, 140, Strand::Star), r(0, 200, 240, Strand::Star)];
        let idx = StartIndex::build(&g);
        let n = |s: i64, e: i64| {
            idx.within_uniform(&g, Query::new(0, s, e, Strand::Plus, true), 41)
                .count()
        };
        assert_eq!(n(110, 111), 1);
        assert_eq!(n(100, 140), 1, "equal to the TFBS is contained");
        assert_eq!(n(99, 111), 0, "sticks out by one on the left");
        assert_eq!(n(110, 141), 0, "sticks out by one on the right");
    }

    #[test]
    fn strand_modes() {
        let g = vec![
            w(0, 10, 20, Strand::Plus, 1),
            w(0, 30, 40, Strand::Minus, 1),
            w(0, 50, 60, Strand::Star, 1),
        ];
        let idx = StartIndex::build(&g);
        // A query spanning 11..60 reaches all three windows.
        let n = |site: Strand, ignore: bool| {
            idx.any_overlaps(&g, Query::new(0, 11, 60, site, ignore))
                .count()
        };
        assert_eq!(n(Strand::Star, false), 3);
        assert_eq!(n(Strand::Plus, false), 2);
        assert_eq!(n(Strand::Minus, false), 2);
        assert_eq!(n(Strand::Plus, true), 3);
    }

    #[test]
    fn subset_mask_is_indexed_by_the_reduced_ranges() {
        let g = vec![
            w(0, 10, 20, Strand::Star, 1),
            w(0, 100, 200, Strand::Star, 2),
        ];
        let q: Vec<(u32, i64, i64, Strand)> = vec![(0, 105, 106, Strand::Star)];
        // One query, two ranges: the mask must have one entry per range, and the
        // second range is the one that overlaps.
        assert_eq!(subset_by_overlaps(&g, &q, true), vec![false, true]);
    }

    #[test]
    fn subset_mask_matches_brute_force() {
        let g = vec![
            w(0, 10, 20, Strand::Star, 1),
            w(0, 100, 200, Strand::Star, 2),
            w(1, 105, 106, Strand::Star, 3),
        ];
        let q: Vec<(u32, i64, i64, Strand)> = vec![
            (0, 15, 150, Strand::Star),
            (1, 105, 106, Strand::Plus),
            (0, 500, 501, Strand::Star),
        ];
        let mask = subset_by_overlaps(&g, &q, true);
        let brute: Vec<bool> = (0..g.len())
            .map(|i| {
                let r = (g[i].chr, g[i].start, g[i].end, g[i].strand);
                q.iter()
                    .any(|x| !brute_force_overlaps(std::slice::from_ref(&r), x, true).is_empty())
            })
            .collect();
        assert_eq!(mask, brute);
        assert_eq!(mask, vec![true, true, true]);
    }

    #[test]
    fn empty_index_is_harmless() {
        let g: Vec<GcWindow> = Vec::new();
        let idx = StartIndex::build(&g);
        assert!(idx.is_empty());
        assert_eq!(idx.chromosomes(), 0);
        assert!(
            idx.any_overlaps(&g, Query::new(0, 1, 2, Strand::Star, true))
                .next()
                .is_none()
        );
        assert!(!idx.overlaps_any(&g, Query::new(0, 1, 2, Strand::Star, true)));
    }
}
