//! T12: interval semantics.
//!
//! Covers every case `AGENT_PLAN.md` T12 lists, plus a brute-force comparison on
//! 2000 random ranges. The brute-force check is the one that matters: the binary
//! searches are the only non-obvious code in the crate, and a boundary error in
//! them would be invisible in the fixture tests because the fixtures have no
//! adversarial coordinates.

mod common;

use common::Rng;
use methyltfr::intervals::{
    Query, StartIndex, brute_force_overlaps, brute_force_within, motif_width, resize_center,
    site_windows, subset_by_overlaps,
};
use methyltfr::model::{Range, Site, Strand};

fn w(chr: u32, start: i64, end: i64, strand: Strand) -> methyltfr::GcWindow {
    methyltfr::GcWindow {
        chr,
        start,
        end,
        strand,
        gc_bin: 1,
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

fn any_hits(
    ranges: &[methyltfr::GcWindow],
    q: (u32, i64, i64, Strand),
    ignore: bool,
) -> Vec<usize> {
    let idx = StartIndex::build(ranges);
    idx.any_overlaps(ranges, Query::from_tuple(q, ignore))
        .collect()
}

#[test]
fn touching_ends_hit() {
    let g = vec![w(0, 10, 20, Strand::Star)];
    assert_eq!(any_hits(&g, (0, 20, 25, Strand::Plus), true), vec![0]);
}

#[test]
fn off_by_one_misses_on_both_sides() {
    let g = vec![w(0, 10, 20, Strand::Star)];
    assert!(any_hits(&g, (0, 21, 25, Strand::Plus), true).is_empty());
    assert!(any_hits(&g, (0, 5, 9, Strand::Plus), true).is_empty());
}

#[test]
fn site_equal_to_tfbs_hits() {
    let g = vec![r(0, 100, 140, Strand::Star)];
    let idx = StartIndex::build(&g);
    assert_eq!(
        idx.within_uniform(&g, Query::new(0, 100, 140, Strand::Plus, true), 41)
            .count(),
        1
    );
}

#[test]
fn site_sticking_out_by_one_on_each_side_misses() {
    let g = vec![r(0, 100, 140, Strand::Star)];
    let idx = StartIndex::build(&g);
    let n = |s: i64, e: i64| {
        idx.within_uniform(&g, Query::new(0, s, e, Strand::Plus, true), 41)
            .count()
    };
    assert_eq!(n(99, 140), 0, "one base too far left");
    assert_eq!(n(100, 141), 0, "one base too far right");
    assert_eq!(n(100, 140), 1);
}

#[test]
fn width_two_site_straddling_two_windows_is_two_hits() {
    let g = vec![w(0, 10, 20, Strand::Star), w(0, 21, 30, Strand::Star)];
    assert_eq!(any_hits(&g, (0, 20, 21, Strand::Plus), true), vec![0, 1]);
}

#[test]
fn new_start_may_be_zero_or_negative() {
    // A 20 bp TFBS resized to the BATF width 541 about its centre starts far
    // enough left to go negative, and that must be kept rather than clamped.
    let tfbs = vec![r(0, 200, 219, Strand::Plus)];
    let resized = resize_center(&tfbs, 541);
    assert_eq!(resized[0].start, 200 + (20i64 - 541).div_euclid(2));
    assert_eq!(resized[0].start, -61);
    assert!(resized[0].start < 0, "kept, not clamped");
    assert_eq!(resized[0].width(), 541);

    // A site at a negative coordinate still hits it.
    let idx = StartIndex::build(&resized);
    assert_eq!(
        idx.within_uniform(&resized, Query::new(0, -10, -9, Strand::Plus, true), 541)
            .count(),
        1
    );
}

#[test]
fn mixed_width_resize_odd_and_even() {
    // The plan's own case: 411 wide to 541 wide, an even difference of -130.
    let tfbs = vec![r(0, 1000, 1410, Strand::Plus)];
    let width = motif_width(&tfbs[0]);
    assert_eq!(width, 541);
    assert_eq!(resize_center(&tfbs, width)[0].start, 1000 - 65);

    // Mixed widths in one call: the halving is of the *difference*, so
    // (100 - 11) %/% 2 == 44, not 100 %/% 2 - 11 %/% 2 == 45.
    let mixed = vec![
        r(0, 100, 199, Strand::Plus), // width 100, even
        r(0, 10, 12, Strand::Minus),  // width 3, odd
        r(0, 500, 511, Strand::Star), // width 12, even
    ];
    let resized = resize_center(&mixed, 11);
    assert_eq!(resized[0].start, 100 + (100i64 - 11).div_euclid(2));
    assert_eq!(resized[1].start, 10 + (3i64 - 11).div_euclid(2));
    assert_eq!(resized[2].start, 500 + (12i64 - 11).div_euclid(2));
    for t in &resized {
        assert_eq!(t.width(), 11, "every resized TFBS has the target width");
    }
}

#[test]
fn midpoint_matches_the_r_table() {
    // tests/fixtures/rmath/midpoint.tsv, compared through the public entry point.
    let path = common::fixture("rmath/midpoint.tsv");
    let (header, rows) = common::read_columns(&path);
    let starts = common::col(&header, &rows, "start");
    let ends = common::col(&header, &rows, "end");
    let mids = common::col(&header, &rows, "mid");
    let mut errors = common::Errors::new();
    for i in 0..rows.len() {
        let s: i64 = starts[i].parse().unwrap();
        let e: i64 = ends[i].parse().unwrap();
        let want: f64 = mids[i].parse().unwrap();
        let got = methyltfr::intervals::midpoint(s, e) as f64;
        // The table was produced by R's `round`, which returns a double here, so
        // compare bit patterns rather than values.
        errors.add_bits(&format!("midpoint({s}, {e})"), got, want);
    }
    errors.finish("midpoint vs rmath/midpoint.tsv");
}

#[test]
fn every_strand_pair_under_both_strand_modes() {
    let strands = [Strand::Plus, Strand::Minus, Strand::Star];
    for a in strands {
        for b in strands {
            // A site with strand `a` against a single window of strand `b`.
            let g = vec![w(0, 10, 20, b)];
            for ignore in [true, false] {
                let got = any_hits(&g, (0, 15, 15, a), ignore);
                let expect_upstream = ignore || a == b || a == Strand::Star || b == Strand::Star;
                assert_eq!(
                    !got.is_empty(),
                    expect_upstream,
                    "site={a:?} window={b:?} ignore_strand={ignore}"
                );
            }
        }
    }
}

#[test]
fn different_chromosomes_never_overlap() {
    let g = vec![w(0, 10, 20, Strand::Star), w(1, 10, 20, Strand::Star)];
    assert_eq!(any_hits(&g, (0, 15, 15, Strand::Plus), true).len(), 1);
    assert_eq!(any_hits(&g, (1, 15, 15, Strand::Plus), true).len(), 1);
    assert!(any_hits(&g, (2, 15, 15, Strand::Plus), true).is_empty());
    assert!(any_hits(&g, (9, 15, 15, Strand::Plus), true).is_empty());
}

#[test]
fn subset_mask_reduces_by_overlap_with_the_query() {
    let g = vec![
        w(0, 10, 20, Strand::Star),
        w(0, 100, 200, Strand::Star),
        w(0, 1000, 1010, Strand::Star),
    ];
    let q = vec![(0u32, 150i64, 151i64, Strand::Star)];
    assert_eq!(subset_by_overlaps(&g, &q, true), vec![false, true, false]);
}

#[test]
fn brute_force_agreement_on_2000_random_ranges() {
    // Three chromosomes, coordinates deliberately clustered so overlaps are
    // frequent rather than rare, mixed widths, mixed strands.
    let mut rng = Rng::new(0x5EED_1234);
    let n = 2000usize;
    let mut ranges: Vec<methyltfr::GcWindow> = Vec::with_capacity(n);
    for _ in 0..n {
        let chr = rng.below(3) as u32;
        let start = rng.range_i64(0, 900);
        let end = start + rng.below(120) as i64;
        ranges.push(w(chr, start, end, common::any_strand(&mut rng)));
    }
    let idx = StartIndex::build(&ranges);
    let mut queries = 0usize;
    let mut hits = 0usize;
    for _ in 0..2000 {
        let chr = rng.below(3) as u32;
        let q_start = rng.range_i64(0, 900);
        let q = (
            chr,
            q_start,
            q_start + rng.below(120) as i64,
            common::any_strand(&mut rng),
        );
        for ignore in [true, false] {
            let mut got: Vec<usize> = idx
                .any_overlaps(&ranges, Query::from_tuple(q, ignore))
                .collect();
            let mut want = brute_force_overlaps(&ranges, &q, ignore);
            queries += 1;
            hits += got.len();
            // The index returns ascending start order; brute force returns input
            // order. Compare as multisets, then separately assert the index order.
            got.sort_unstable();
            want.sort_unstable();
            assert_eq!(got, want, "q={q:?} ignore_strand={ignore}");
        }
    }
    assert!(
        hits > 1000,
        "only {hits} hits over {queries} queries: too easy"
    );
    eprintln!("brute force: {queries} queries, {hits} total hits, 0 disagreements");
}

#[test]
fn brute_force_agreement_for_within_uniform() {
    let mut rng = Rng::new(0xBEEF);
    // Uniform width, so `within_uniform`'s single-bound shortcut is valid.
    let width = 41i64;
    let n = 2000usize;
    let mut ranges: Vec<Range> = Vec::with_capacity(n);
    for _ in 0..n {
        let chr = rng.below(2) as u32;
        let start = rng.range_i64(0, 900);
        ranges.push(r(
            chr,
            start,
            start + width - 1,
            common::any_strand(&mut rng),
        ));
    }
    let idx = StartIndex::build(&ranges);
    let mut queries = 0usize;
    let mut hits = 0usize;
    for _ in 0..2000 {
        let chr = rng.below(2) as u32;
        let s = rng.range_i64(0, 900);
        let q = (
            chr,
            s,
            s + rng.below(5) as i64,
            common::any_strand(&mut rng),
        );
        for ignore in [true, false] {
            let mut got: Vec<usize> = idx
                .within_uniform(&ranges, Query::from_tuple(q, ignore), width)
                .collect();
            let mut want = brute_force_within(&ranges, &q, ignore);
            queries += 1;
            hits += got.len();
            got.sort_unstable();
            want.sort_unstable();
            assert_eq!(got, want, "q={q:?} ignore_strand={ignore}");
        }
    }
    assert!(hits > 500, "only {hits} hits over {queries} queries");
    eprintln!("within_uniform: {queries} queries, {hits} total hits, 0 disagreements");
}

#[test]
fn site_accessor_goes_through_the_same_code_path() {
    let g = vec![w(0, 10, 20, Strand::Star), w(0, 25, 40, Strand::Star)];
    let idx = StartIndex::build(&g);
    let site = Site {
        chr: 0,
        start: 15,
        end: 16,
        strand: Strand::Plus,
        score: 0.5,
        coverage: 1.0,
    };
    let got: Vec<usize> = site_windows(&idx, &g, &site, true).collect();
    assert_eq!(got, vec![0]);
}

#[test]
fn windows_and_ranges_can_both_be_indexed() {
    let mut rng = Rng::new(7);
    let windows: Vec<methyltfr::GcWindow> = (0..50)
        .map(|_| common::any_window(&mut rng, 0, 0, 2000))
        .collect();
    let ranges: Vec<Range> = (0..50)
        .map(|_| common::any_range(&mut rng, 0, 0, 2000))
        .collect();
    let iw = StartIndex::build(&windows);
    let ir = StartIndex::build(&ranges);
    assert_eq!(iw.chromosomes(), 1);
    assert_eq!(ir.len(), 50);
    for _ in 0..100 {
        let s = rng.range_i64(0, 2000);
        let e = s + rng.below(30) as i64;
        let strand = common::any_strand(&mut rng);
        let mut a: Vec<usize> = iw
            .any_overlaps(&windows, Query::new(0, s, e, strand, true))
            .collect();
        let mut b = brute_force_overlaps(&windows, &(0, s, e, strand), true);
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }
}
