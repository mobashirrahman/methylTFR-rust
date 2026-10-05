//! T21: the expected profile against the oracle outputs, for both fixture sets.

mod common;

use common::{Errors, batf, batf_1d99721};
use methyltfr::expected::compute_expectations;
use methyltfr::gc::bin_means_indexed;
use methyltfr::model::GcFreq;

fn check(dir: &std::path::Path, label: &str) {
    let (annotation, sites) = common::load_batf(dir);
    let bins = bin_means_indexed(&sites, &annotation.gc_windows, annotation.ignore_strand).unwrap();
    let gcfreq: &GcFreq = &annotation.motifs[0].gcfreq;
    let profile = compute_expectations(gcfreq, &bins).unwrap_or_else(|e| panic!("{label}: {e}"));

    let path = dir.join("expected_tierB/expected_profile.tsv");
    let (header, rows) = common::read_columns(&path);
    let want_x = common::col(&header, &rows, "x");
    let want_value = common::col(&header, &rows, "value");

    assert_eq!(profile.len(), gcfreq.cols);
    assert_eq!(profile.len(), want_x.len(), "{label}: profile length");

    let mut errors = Errors::new();
    for i in 0..rows.len() {
        errors.add(
            &format!("{label} x[{i}]"),
            profile.x[i],
            want_x[i].parse::<f64>().unwrap(),
        );
        errors.add(
            &format!("{label} value[{i}]"),
            profile.value[i],
            want_value[i].parse::<f64>().unwrap(),
        );
    }
    errors.finish(&format!("expected_profile {label}"));

    // AGENT_PLAN.md section 2.7: even L skips zero.
    if gcfreq.cols % 2 == 0 {
        assert!(!profile.x.contains(&0.0), "{label}: even L must skip zero");
        assert_eq!(
            profile.x.first().copied(),
            Some(-(gcfreq.cols as f64 / 2.0))
        );
        assert_eq!(profile.x.last().copied(), Some(gcfreq.cols as f64 / 2.0));
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

fn bins(labels: &[u8], means: &[f64]) -> methyltfr::GcBins {
    methyltfr::GcBins {
        bin: labels.to_vec(),
        mean: means.to_vec(),
        n_hits: means.iter().map(|_| 1u64).collect(),
    }
}

#[test]
fn grid_lengths_1_2_4_5_and_odd_and_even() {
    for l in [1usize, 2, 3, 4, 5, 101, 512, 601] {
        let f = GcFreq::new(1, l, vec![1.0 / l as f64; l]).unwrap();
        let p = compute_expectations(&f, &bins(&[1], &[0.5])).unwrap();
        assert_eq!(p.len(), l);
        let h = (l / 2) as f64;
        if l % 2 == 1 {
            assert_eq!(p.x[0], -h);
            assert_eq!(p.x[l - 1], h);
            assert!(p.x.contains(&0.0), "odd L must include zero");
        } else {
            assert_eq!(p.x[0], -h);
            assert_eq!(p.x[l - 1], h);
            assert!(!p.x.contains(&0.0), "even L must skip zero");
        }
        // Every value is gcfreq[0][j] * mean[0].
        for v in &p.value {
            assert!((v - (1.0 / l as f64) * 0.5).abs() < 1e-15);
        }
    }
}

#[test]
fn a_missing_bin_is_an_error() {
    let f = GcFreq::new(5, 8, vec![0.125; 40]).unwrap();
    for populated in [
        vec![1u8, 2, 3, 4],
        vec![1, 2, 3],
        vec![1],
        vec![1, 2, 3, 4, 5, 1],
    ] {
        let means = vec![0.5; populated.len()];
        let e = compute_expectations(&f, &bins(&populated, &means)).unwrap_err();
        assert!(
            e.to_string().contains("non-conformable") || e.to_string().contains("must be 1..=5"),
            "unexpected message for {populated:?}: {e}"
        );
    }
}

#[test]
fn bins_must_start_at_one() {
    let f = GcFreq::new(3, 4, vec![0.25; 12]).unwrap();
    let e = compute_expectations(&f, &bins(&[2, 3, 4], &[0.5; 3])).unwrap_err();
    assert!(e.to_string().contains("must be 1..=3"), "{e}");
}

#[test]
fn columns_of_gcfreq_sum_to_one() {
    // Fixture fact from AGENT_PLAN.md section 4; if this ever fails the expected
    // profile is not a distribution and the golden values mean nothing.
    for dir in [batf(""), batf_1d99721("")] {
        let (annotation, _) = common::load_batf(&dir);
        let g = &annotation.motifs[0].gcfreq;
        for j in 0..g.cols {
            let sum: f64 = (0..g.rows).map(|i| g.row(i)[j]).sum();
            assert!(
                (sum - 1.0).abs() < 1e-12,
                "{}: column {j} sums to {sum}",
                dir.display()
            );
        }
    }
}
