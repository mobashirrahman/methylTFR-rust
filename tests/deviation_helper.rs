//! T22: `dev_helper`, the observed profile, and the two together.

mod common;

use common::{Errors, batf};
use methyltfr::deviation::{compute_observed, dev_helper, motif_width_of};
use methyltfr::expected::Profile;
use methyltfr::model::{Range, Site, Strand};

fn profile(x: &[f64], v: &[f64]) -> Profile {
    Profile {
        x: x.to_vec(),
        value: v.to_vec(),
    }
}

/// One value inside each of the five cut intervals.
fn one_per_interval(values: [f64; 5]) -> Profile {
    profile(&[-225.0, -100.0, 0.0, 100.0, 225.0], &values)
}

#[test]
fn k_0_is_na() {
    // Nothing lands inside the cut range.
    let p = profile(&[-1e9, -251.0, -300.0, 300.0], &[1.0, 1.0, 1.0, 1.0]);
    assert_eq!(dev_helper(&p), None);
    // A profile with no values at all.
    assert_eq!(dev_helper(&profile(&[], &[])), None);
}

#[test]
fn k_1_gives_one() {
    // v = 0 is deliberately absent: it makes the denominator zero, which the
    // zero-denominator test below covers.
    for v in [0.5, 1.0, -2.0, 1e-9] {
        let p = profile(&[-225.0], &[v]);
        let d = dev_helper(&p).expect("k = 1 is not NA");
        assert_eq!(d, 1.0, "v = {v}: D = mean / ((mean + mean) / 2) = 1");
    }
}

#[test]
fn k_2_takes_the_first_mean() {
    let p = profile(&[-225.0, 225.0], &[0.5, 1.5]);
    assert_eq!(dev_helper(&p), Some(0.5 / ((0.5 + 1.5) / 2.0)));
    // R's `mean[(k + 1) %/% 2]` with k = 2 is mean[1], not the midpoint.
    let p = profile(&[-225.0, 225.0], &[2.0, 6.0]);
    assert_eq!(dev_helper(&p), Some(2.0 / ((2.0 + 6.0) / 2.0)));
}

#[test]
fn k_3_takes_the_second_mean() {
    // Only intervals 0, 1 and 4 are populated, so k = 3 and the means are
    // 0.2, 0.4, 0.8; the numerator is the second of the three.
    let p = profile(&[-225.0, -100.0, 225.0], &[0.2, 0.4, 0.8]);
    assert_eq!(dev_helper(&p), Some(0.4 / ((0.2 + 0.8) / 2.0)));
}

#[test]
fn k_4_takes_the_second_mean() {
    // Four populated intervals; the numerator is still the second.
    let p = profile(&[-225.0, -100.0, 100.0, 225.0], &[0.2, 0.4, 0.8, 0.6]);
    assert_eq!(dev_helper(&p), Some(0.4 / ((0.2 + 0.6) / 2.0)));
}

#[test]
fn k_5_takes_the_middle_mean() {
    let p = one_per_interval([0.2, 0.4, 0.6, 0.8, 1.0]);
    assert_eq!(dev_helper(&p), Some(0.6 / ((0.2 + 1.0) / 2.0)));
}

#[test]
fn every_break_value_lands_in_the_interval_above() {
    // -250 is NA; -200 in interval 0; -25 in 1; 25 in 2; 200 in 3; 250 in 4.
    // Each of the five populated, so the numerator is the third: 3.0.
    let p = profile(
        &[-200.0, -25.0, 25.0, 200.0, 250.0],
        &[1.0, 2.0, 3.0, 4.0, 5.0],
    );
    assert_eq!(dev_helper(&p), Some(3.0 / ((1.0 + 5.0) / 2.0)));

    // Shifting -25 to -24 must move it into interval 2 and change k.
    let p = profile(
        &[-200.0, -24.0, 25.0, 200.0, 250.0],
        &[1.0, 2.0, 3.0, 4.0, 5.0],
    );
    // Now interval 1 is empty and interval 2 holds 2.0 and 3.0 -> mean 2.5.
    assert_eq!(dev_helper(&p), Some(2.5 / ((1.0 + 5.0) / 2.0)));
}

#[test]
fn break_edges_just_outside_are_na() {
    for x in [-250.0, -251.0, 250.5, 1e9] {
        assert_eq!(
            dev_helper(&profile(&[x], &[1.0])),
            None,
            "x = {x} is outside the cut range"
        );
    }
    for x in [-249.999, 250.0] {
        assert!(dev_helper(&profile(&[x], &[1.0])).is_some(), "x = {x}");
    }
}

#[test]
fn a_zero_denominator_is_passed_through_not_caught() {
    // Only interval 0 populated with mean 0 -> 0 / 0.
    assert!(dev_helper(&profile(&[-225.0], &[0.0])).unwrap().is_nan());
    // Only interval 0 with a negative mean -> -Inf.
    assert_eq!(
        dev_helper(&profile(&[-225.0], &[-1.0])).unwrap(),
        1.0,
        "a negative mean over itself is still 1"
    );
    // Intervals 0 and 4 with means summing to zero -> infinite result.
    let d = dev_helper(&profile(&[-225.0, 225.0], &[1.0, -1.0])).unwrap();
    assert!(d.is_infinite(), "got {d}");
}

#[test]
fn the_observed_profile_matches_the_oracle_as_a_multiset() {
    let (annotation, sites) = common::load_batf(&batf(""));
    let motif = &annotation.motifs[0];
    let width = motif_width_of(&motif.tfbs).unwrap();
    let observed = compute_observed(
        &motif.tfbs,
        width,
        None,
        &sites,
        annotation.ignore_strand,
        &motif.name,
    )
    .unwrap();

    let path = batf("expected_tierB/observed_profile.tsv");
    let (header, rows) = common::read_columns(&path);
    let want_x = common::col(&header, &rows, "x");
    let want_value = common::col(&header, &rows, "value");
    assert_eq!(observed.len(), want_x.len(), "hit count");

    let mut got: Vec<(i64, String)> = observed
        .x
        .iter()
        .zip(observed.value.iter())
        .map(|(x, v)| (*x as i64, format!("{v:?}")))
        .collect();
    let mut want: Vec<(i64, String)> = want_x
        .iter()
        .zip(want_value.iter())
        .map(|(x, v)| {
            (
                x.parse::<i64>().unwrap(),
                format!("{:?}", v.parse::<f64>().unwrap()),
            )
        })
        .collect();
    got.sort();
    want.sort();
    assert_eq!(got, want, "observed profile as (x, value) multisets");

    // Plan section 4 for the pinned SHA: 29 hits, interval counts 3, 9, 1, 9, 5.
    assert_eq!(observed.len(), 29);
    let mut counts = [0usize; 5];
    let mut outside = 0usize;
    for x in &observed.x {
        match methyltfr::rmath::cut_index(*x) {
            Some(i) => counts[i] += 1,
            None => outside += 1,
        }
    }
    assert_eq!(counts, [3, 9, 1, 9, 5], "plan section 4 observed intervals");
    assert_eq!(
        outside, 2,
        "two hits fall outside +/-250 and dev_helper skips them"
    );
}

#[test]
fn observed_deviation_matches_the_oracle() {
    let (annotation, sites) = common::load_batf(&batf(""));
    let motif = &annotation.motifs[0];
    let width = motif_width_of(&motif.tfbs).unwrap();
    let observed = compute_observed(
        &motif.tfbs,
        width,
        None,
        &sites,
        annotation.ignore_strand,
        &motif.name,
    )
    .unwrap();
    let d = dev_helper(&observed).unwrap();

    let (header, rows) = common::read_columns(&batf("expected_tierB/expected_dev.tsv"));
    let want = common::col(&header, &rows, "obs_d")[0]
        .parse::<f64>()
        .unwrap();
    let mut errors = Errors::new();
    errors.add("obs_d", d, want);
    errors.finish("observed deviation");
    // Plan section 4: obs_d = 2.727272727.
    assert!((d - 2.727272727272727).abs() < 1e-9, "got {d}");
}

#[test]
fn x_uses_site_start_never_site_end() {
    let tfbs = vec![Range {
        chr: 0,
        start: 1000,
        end: 1410,
        strand: Strand::Star,
    }];
    let width = motif_width_of(&tfbs).unwrap();
    // Two sites one base apart, so using site.end instead would shift every x.
    let sites: Vec<Site> = [1200i64, 1201]
        .iter()
        .map(|s| Site {
            chr: 0,
            start: *s,
            end: *s + 1,
            strand: Strand::Star,
            score: 1.0,
            coverage: 1.0,
        })
        .collect();
    let obs = compute_observed(&tfbs, width, None, &sites, true, "T").unwrap();
    let mut xs = obs.x.clone();
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(xs, vec![-5.0, -4.0]);
}

#[test]
fn a_site_in_k_tfbs_gives_k_hits() {
    let tfbs = vec![
        Range {
            chr: 0,
            start: 100,
            end: 510,
            strand: Strand::Star,
        },
        Range {
            chr: 0,
            start: 300,
            end: 710,
            strand: Strand::Star,
        },
        Range {
            chr: 0,
            start: 900,
            end: 1310,
            strand: Strand::Star,
        },
    ];
    let width = motif_width_of(&tfbs).unwrap();
    let sites = vec![Site {
        chr: 0,
        start: 400,
        end: 401,
        strand: Strand::Star,
        score: 0.5,
        coverage: 1.0,
    }];
    let obs = compute_observed(&tfbs, width, None, &sites, true, "T").unwrap();
    assert_eq!(obs.len(), 2, "the third TFBS does not contain the site");
}

#[test]
fn no_hits_names_the_motif() {
    let tfbs = vec![Range {
        chr: 0,
        start: 1_000_000,
        end: 1_000_410,
        strand: Strand::Star,
    }];
    let width = motif_width_of(&tfbs).unwrap();
    let sites = vec![Site {
        chr: 0,
        start: 5,
        end: 6,
        strand: Strand::Star,
        score: 1.0,
        coverage: 1.0,
    }];
    let e = compute_observed(&tfbs, width, None, &sites, true, "BATF").unwrap_err();
    assert_eq!(
        e.to_string(),
        "No methylation sites found in the BATF binding sites"
    );
}
