//! The two-pointer sweep must return exactly what the binary search returns.
//!
//! `docs/AGENT_PLAN.md` M7 replaces the per-site binary searches with sweeps for
//! speed. A sweep is the kind of optimisation that can be right on ordinary input
//! and wrong on the last site of a chromosome, so this compares the two
//! implementations hit for hit on adversarial input: mixed widths, mixed strands,
//! out-of-order site ends, abutting windows, extreme coordinates, and empty and
//! missing chromosomes.

mod common;

use common::Rng;
use methyltfr::deviation::{compute_observed, compute_observed_swept, motif_width_of};
use methyltfr::gc::{bin_means, bin_means_swept};
use methyltfr::intervals::StartIndex;
use methyltfr::model::{ChromSites, GcWindow, Range, Site, SortedMethylome, Strand};

fn site(chr: u32, start: i64, end: i64, strand: Strand, score: f64) -> Site {
    Site {
        chr,
        start,
        end,
        strand,
        score,
        coverage: 1.0,
    }
}

fn window(chr: u32, start: i64, end: i64, strand: Strand, bin: u8) -> GcWindow {
    GcWindow {
        chr,
        start,
        end,
        strand,
        gc_bin: bin,
    }
}

/// Every (site index, window index) the binary search reports, as a sorted
/// multiset.
fn search_hits(
    sites: &[Site],
    windows: &[GcWindow],
    index: &StartIndex,
    ignore_strand: bool,
) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for (i, s) in sites.iter().enumerate() {
        for w in index.any_overlaps(
            windows,
            methyltfr::intervals::Query::from_site(s, ignore_strand),
        ) {
            out.push((i, w));
        }
    }
    out
}

/// Sorted position -> original index, built the same way
/// [`SortedMethylome::new`] builds its order. `SortedMethylome` deliberately does
/// not keep this mapping: it would cost 4 bytes per site for something only a test
/// needs.
fn sort_permutation(sites: &[Site]) -> Vec<u32> {
    let mut perm: Vec<u32> = (0..sites.len() as u32).collect();
    perm.sort_by_key(|&i| {
        let s = &sites[i as usize];
        (s.chr, s.start, s.end)
    });
    perm
}

/// The same hits from the sweep, identified by original site index and window index.
fn sweep_hits(
    sorted: &SortedMethylome,
    sites: &[Site],
    windows: &[GcWindow],
    index: &StartIndex,
    ignore_strand: bool,
) -> Vec<(usize, usize)> {
    let perm = sort_permutation(sites);
    // The sweep does not report indices, so the bin means and hit counts are the
    // comparable surface. Here we recover a hit list by pairing each sweep's
    // per-site windows with the binary search's, which is what makes the assertion
    // exact rather than aggregate.
    let mut out = Vec::new();
    for group in sorted.groups() {
        let chrom = sorted.chrom(group);
        let slice = index.slice_of(chrom.chr);
        if slice.is_empty() {
            continue;
        }
        let max_width = index.max_width_of(chrom.chr);
        let mut sweep = methyltfr::intervals::Sweep::new(windows, slice, max_width);
        for i in 0..chrom.len() {
            let s = chrom.starts[i];
            let e = chrom.ends[i];
            let site_strand = chrom.strand(i);
            sweep.for_each_candidate(s, e, |w| {
                let win = &windows[w];
                if win.chr != chrom.chr || win.start > e || win.end < s {
                    return;
                }
                if !win.strand.compatible(site_strand, ignore_strand) {
                    return;
                }
                out.push((perm[group.begin + i] as usize, w));
            });
        }
    }
    out.sort_unstable();
    out
}

fn check_case(sites: &[Site], windows: &[GcWindow], label: &str) {
    for ignore_strand in [true, false] {
        let index = StartIndex::build(windows);
        let sorted = SortedMethylome::new(sites);

        let mut want = search_hits(sites, windows, &index, ignore_strand);
        want.sort_unstable();
        let got = sweep_hits(&sorted, sites, windows, &index, ignore_strand);
        assert_eq!(
            got, want,
            "{label}: sweep and search disagree on hits \
             (ignore_strand = {ignore_strand})"
        );

        // And the aggregate the pipeline actually uses must agree too.
        let by_search = bin_means(sites, windows, &index, ignore_strand);
        let by_sweep = bin_means_swept(&sorted, windows, &index, ignore_strand);
        match (by_search, by_sweep) {
            (Ok(a), Ok(b)) => {
                assert_eq!(a.bin, b.bin, "{label}: bins");
                assert_eq!(
                    a.n_hits, b.n_hits,
                    "{label}: n_hits (ignore_strand = {ignore_strand})"
                );
            }
            (Err(x), Err(y)) => assert_eq!(
                x.to_string(),
                y.to_string(),
                "{label}: both fail but differently"
            ),
            (a, b) => panic!("{label}: one side failed and the other did not: {a:?} vs {b:?}"),
        }
    }
}

#[test]
fn random_mixed_widths_and_strands() {
    let mut rng = Rng::new(0x1234_5678);
    for case in 0..40 {
        let n = 1 + rng.below(400) as usize;
        let sites: Vec<Site> = (0..n)
            .map(|_| {
                let chr = rng.below(3) as u32;
                let start = rng.range_i64(-50, 500);
                site(
                    chr,
                    start,
                    start + rng.below(4) as i64,
                    common::any_strand(&mut rng),
                    rng.unit(),
                )
            })
            .collect();
        let m = 1 + rng.below(200) as usize;
        let windows: Vec<GcWindow> = (0..m)
            .map(|_| {
                let chr = rng.below(3) as u32;
                let start = rng.range_i64(-50, 500);
                window(
                    chr,
                    start,
                    start + rng.below(40) as i64,
                    common::any_strand(&mut rng),
                    (rng.below(5) + 1) as u8,
                )
            })
            .collect();
        check_case(&sites, &windows, &format!("random case {case}"));
    }
}

#[test]
fn abutting_and_overlapping_windows() {
    let windows: Vec<GcWindow> = (0..40)
        .map(|i| window(0, i * 30, i * 30 + 29, Strand::Star, (i % 5 + 1) as u8))
        .collect();
    let sites: Vec<Site> = (0..200)
        .map(|i| site(0, i * 7, i * 7 + 1, Strand::Plus, (i % 7) as f64 / 7.0))
        .collect();
    check_case(&sites, &windows, "tiled 30 bp");
}

#[test]
fn site_ends_out_of_order() {
    // Sites sorted by start whose ends are *not* monotone: a wide site followed by
    // a narrow one at a nearby position. This is the case the sweep's running-end
    // bound exists for.
    let sites = vec![
        site(0, 10, 110, Strand::Star, 0.1),
        site(0, 11, 12, Strand::Star, 0.2),
        site(0, 12, 13, Strand::Star, 0.3),
        site(0, 20, 21, Strand::Star, 0.4),
        site(0, 21, 21, Strand::Star, 0.5),
    ];
    let windows = vec![
        window(0, 15, 25, Strand::Star, 1),
        window(0, 100, 130, Strand::Star, 2),
    ];
    check_case(&sites, &windows, "ends out of order");
}

#[test]
fn extreme_coordinates() {
    let sites = vec![
        site(0, i64::MIN, i64::MIN + 1, Strand::Star, 0.5),
        site(0, i64::MAX - 1, i64::MAX, Strand::Star, 0.5),
        site(0, 0, 1, Strand::Star, 0.5),
    ];
    let windows = vec![
        window(0, i64::MIN, i64::MIN + 29, Strand::Star, 1),
        window(0, i64::MAX - 29, i64::MAX, Strand::Star, 2),
    ];
    check_case(&sites, &windows, "i64 extremes");
}

#[test]
fn empty_and_missing_chromosomes() {
    check_case(&[], &[window(0, 1, 30, Strand::Star, 1)], "no sites");
    check_case(&[site(9, 1, 2, Strand::Star, 0.5)], &[], "no windows");
    // Windows only on a chromosome the sites never touch.
    check_case(
        &[site(0, 100, 101, Strand::Star, 0.5)],
        &[window(1, 100, 130, Strand::Star, 1)],
        "disjoint chromosomes",
    );
}

#[test]
fn tfbs_sweep_matches_search() {
    let mut rng = Rng::new(0xABCD);
    for case in 0..30 {
        let n_tfbs = 1 + rng.below(200) as usize;
        let width0 = 11 + rng.below(200) as i64;
        let tfbs: Vec<Range> = (0..n_tfbs)
            .map(|_| {
                let chr = rng.below(3) as u32;
                let start = rng.range_i64(-50, 500);
                Range {
                    chr,
                    start,
                    // Mixed widths, so W's parity varies between cases and the
                    // midpoint lands on a .5 inside them.
                    end: start + rng.below(60) as i64,
                    strand: common::any_strand(&mut rng),
                }
            })
            .collect();
        let width = motif_width_of(&tfbs).unwrap();
        let n_sites = 1 + rng.below(300) as usize;
        let sites: Vec<Site> = (0..n_sites)
            .map(|_| {
                let chr = rng.below(3) as u32;
                let start = rng.range_i64(-100, 600);
                site(
                    chr,
                    start,
                    start + rng.below(3) as i64,
                    common::any_strand(&mut rng),
                    rng.unit(),
                )
            })
            .collect();

        for ignore_strand in [true, false] {
            let by_search = compute_observed(&tfbs, width, None, &sites, ignore_strand, "M");
            let sorted = SortedMethylome::new(&sites);
            let by_sweep = compute_observed_swept(&tfbs, width, None, &sorted, ignore_strand, "M");
            match (by_search, by_sweep) {
                (Ok(a), Ok(b)) => {
                    // Same multiset of (x, value): the sweep visits sites in sorted
                    // order, so the sequence may differ, but nothing else may.
                    let mut a_pairs: Vec<(i64, u64)> =
                        a.x.iter()
                            .zip(a.value.iter())
                            .map(|(x, v)| (*x as i64, v.to_bits()))
                            .collect();
                    let mut b_pairs: Vec<(i64, u64)> =
                        b.x.iter()
                            .zip(b.value.iter())
                            .map(|(x, v)| (*x as i64, v.to_bits()))
                            .collect();
                    a_pairs.sort_unstable();
                    b_pairs.sort_unstable();
                    assert_eq!(
                        a_pairs, b_pairs,
                        "case {case}: (x, value) multisets differ \
                         (ignore_strand = {ignore_strand}), width0 = {width0}"
                    );
                }
                (Err(x), Err(y)) => assert_eq!(x.to_string(), y.to_string()),
                (a, b) => panic!("case {case}: one side failed: {a:?} vs {b:?}"),
            }
        }
    }
}

#[test]
fn a_pre_sorted_methylome_gives_identical_output() {
    // The property docs/divergences.md D4 relies on: when the input is already in
    // coordinate order, the sweep and the binary search produce the *same
    // sequence*, not merely the same multiset.
    let mut rng = Rng::new(0x5EED);
    let mut sites: Vec<Site> = (0..500)
        .map(|i| {
            site(
                0,
                i * 3,
                i * 3 + 1,
                common::any_strand(&mut rng),
                rng.unit(),
            )
        })
        .collect();
    sites.sort_by_key(|s| (s.chr, s.start, s.end));
    let windows: Vec<GcWindow> = (0..50)
        .map(|i| window(0, i * 30, i * 30 + 29, Strand::Star, (i % 5 + 1) as u8))
        .collect();
    let tfbs: Vec<Range> = (0..20)
        .map(|i| Range {
            chr: 0,
            start: i * 200,
            end: i * 200 + 50,
            strand: Strand::Star,
        })
        .collect();
    let width = motif_width_of(&tfbs).unwrap();
    let index = StartIndex::build(&windows);
    let sorted = SortedMethylome::new(&sites);

    let a = compute_observed(&tfbs, width, None, &sites, true, "M").unwrap();
    let b = compute_observed_swept(&tfbs, width, None, &sorted, true, "M").unwrap();
    assert_eq!(a.x, b.x, "x sequences differ on pre-sorted input");
    assert_eq!(
        a.value, b.value,
        "value sequences differ on pre-sorted input"
    );

    let ba = bin_means(&sites, &windows, &index, true).unwrap();
    let bb = bin_means_swept(&sorted, &windows, &index, true).unwrap();
    assert_eq!(ba.bin, bb.bin);
    assert_eq!(ba.n_hits, bb.n_hits);
    for (x, y) in ba.mean.iter().zip(bb.mean.iter()) {
        assert_eq!(x.to_bits(), y.to_bits(), "bin means differ bit for bit");
    }
}

/// The `ChromSites` view has to stay consistent with the sites it was built from.
#[test]
fn chrom_view_is_consistent() {
    let sites = vec![
        site(1, 5, 6, Strand::Minus, 0.25),
        site(0, 1, 2, Strand::Plus, 0.75),
    ];
    let sorted = SortedMethylome::new(&sites);
    let g = sorted.group(1).unwrap();
    let chrom: ChromSites<'_> = sorted.chrom(g);
    assert_eq!(chrom.chr, 1);
    assert_eq!(chrom.starts, &[5]);
    assert_eq!(chrom.ends, &[6]);
    assert_eq!(chrom.strand(0), Strand::Minus);
    assert_eq!(chrom.scores[0], 0.25);
}
