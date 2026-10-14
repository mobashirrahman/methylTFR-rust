//! T20: GC bin means against the Tier A / Tier B oracle outputs.
//!
//! `expected_bins.tsv` has three columns: `gc_bin mean n_hits`. The means are
//! compared at the section 3 tolerance; `n_hits` is an integer count and is
//! compared exactly, because that is what the whole "mean over hits, not over
//! sites" rule is about.

mod common;

use common::{Errors, batf, batf_1d99721};
use methyltfr::gc::bin_means_indexed;
use methyltfr::model::Strand;

fn check(dir: &std::path::Path, label: &str) {
    let (annotation, sites) = common::load_batf(dir);
    let bins = bin_means_indexed(&sites, &annotation.gc_windows, annotation.ignore_strand)
        .unwrap_or_else(|e| panic!("{label}: {e}"));

    let path = dir.join("expected_tierB/expected_bins.tsv");
    let (header, rows) = common::read_columns(&path);
    let want_bin = common::col(&header, &rows, "gc_bin");
    let want_mean = common::col(&header, &rows, "mean");
    let want_hits = common::col(&header, &rows, "n_hits");

    assert_eq!(
        bins.len(),
        want_bin.len(),
        "{label}: {} bins vs {} in the oracle",
        bins.len(),
        want_bin.len()
    );

    let mut errors = Errors::new();
    for i in 0..rows.len() {
        assert_eq!(
            bins.bin[i].to_string(),
            want_bin[i],
            "{label}: bin label at position {i}"
        );
        errors.add(
            &format!("{label} bin {}", want_bin[i]),
            bins.mean[i],
            want_mean[i].parse::<f64>().unwrap(),
        );
        assert_eq!(
            bins.n_hits[i].to_string(),
            want_hits[i],
            "{label}: n_hits for bin {}",
            want_bin[i]
        );
    }
    errors.finish(&format!("gc_bins {label}"));

    // The fixture facts AGENT_PLAN.md section 4 states for the pinned SHA.
    if label == "batf" {
        assert_eq!(
            bins.n_hits,
            vec![8, 59, 96, 152, 729],
            "plan section 4 GC hits per bin"
        );
        assert_eq!(bins.n_hits.iter().sum::<u64>(), 1044);
        // 1000 sites, 1044 hits: 44 sites straddle two abutting windows.
        assert_eq!(sites.len(), 1000);
        assert_eq!(bins.n_hits.iter().sum::<u64>() - sites.len() as u64, 44);
        assert_eq!(
            bins.bin,
            vec![1, 2, 3, 4, 5],
            "all five bins are populated here"
        );
    }
}

#[test]
fn matches_the_oracle_for_the_pinned_sha() {
    check(&batf(""), "batf");
}

#[test]
fn matches_the_oracle_for_the_published_revision() {
    check(&batf_1d99721(""), "batf_1d99721");
}

#[test]
fn no_hits_is_the_upstream_error() {
    let (_annotation, sites) = common::load_batf(&batf(""));
    let empty: Vec<methyltfr::GcWindow> = Vec::new();
    let e = bin_means_indexed(&sites, &empty, true).unwrap_err();
    assert_eq!(
        e.to_string(),
        "No methylation sites found in the GC distribution"
    );

    // Windows on another chromosome are as good as no windows at all.
    let (annotation, _) = common::load_batf(&batf(""));
    let elsewhere: Vec<_> = annotation
        .gc_windows
        .iter()
        .map(|w| methyltfr::GcWindow { chr: 999, ..*w })
        .collect();
    let e = bin_means_indexed(&sites, &elsewhere, true).unwrap_err();
    assert_eq!(
        e.to_string(),
        "No methylation sites found in the GC distribution"
    );
}

#[test]
fn an_unpopulated_bin_is_absent_not_zero() {
    // Windows in bins 1 and 3, one site in bin 1. Bin 3 must not appear, and the
    // bin list must be 1, not 1 and 3 with a zero.
    let w = |start: i64, end: i64, bin: u8| methyltfr::GcWindow {
        chr: 0,
        start,
        end,
        strand: Strand::Star,
        gc_bin: bin,
    };
    let site = |start: i64, score: f64| methyltfr::Site {
        chr: 0,
        start,
        end: start + 1,
        strand: Strand::Star,
        score,
        coverage: 1.0,
    };
    let bins = bin_means_indexed(
        &[site(15, 0.4), site(215, 0.8)],
        &[w(10, 20, 1), w(210, 220, 3)],
        true,
    )
    .unwrap();
    assert_eq!(bins.bin, vec![1, 3]);
    assert_eq!(bins.mean, vec![0.4, 0.8]);
    assert_eq!(bins.position_of(1), Some(0));
    assert_eq!(bins.position_of(2), None);
    assert_eq!(bins.position_of(3), Some(1));
    assert_eq!(bins.position_of(5), None);
}

#[test]
fn bin_means_survive_extreme_coordinates() {
    // A site at i64::MAX must not overflow the arithmetic; the overlap test is
    // on `window.start <= e`, which is a comparison, not an addition.
    let w = methyltfr::GcWindow {
        chr: 0,
        start: i64::MAX - 10,
        end: i64::MAX,
        strand: Strand::Star,
        gc_bin: 4,
    };
    let site = methyltfr::Site {
        chr: 0,
        start: i64::MAX - 5,
        end: i64::MAX - 4,
        strand: Strand::Plus,
        score: 0.5,
        coverage: 1.0,
    };
    let bins = bin_means_indexed(&[site], &[w], true).unwrap();
    assert_eq!(bins.bin, vec![4]);
    assert_eq!(bins.mean, vec![0.5]);

    // And at the bottom of the range, where a resize would go negative.
    let w = methyltfr::GcWindow {
        chr: 0,
        start: i64::MIN,
        end: i64::MIN + 29,
        strand: Strand::Star,
        gc_bin: 1,
    };
    let site = methyltfr::Site {
        chr: 0,
        start: i64::MIN,
        end: i64::MIN + 1,
        strand: Strand::Plus,
        score: 0.25,
        coverage: 1.0,
    };
    let bins = bin_means_indexed(&[site], &[w], true).unwrap();
    assert_eq!(bins.mean, vec![0.25]);
}
