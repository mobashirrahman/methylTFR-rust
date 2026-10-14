//! T23: the end-to-end golden values.
//!
//! Two assertions, deliberately separate:
//!
//! * Against `expected_dev.tsv` from the oracles, at the section 3 tolerance of
//!   1e-10. This is the parity check.
//! * Against the published pkgdown column for `batf_1d99721` -- `dev = 1.743674`
//!   and `exp_dev = 0.9835985` -- at 5e-7, which is what 7 printed significant
//!   digits means. This is the check that the fixtures are the right ones.

mod common;

use common::{Errors, batf, batf_1d99721};
use methyltfr::deviation::{compute_observed, deviations, motif_width_of};
use methyltfr::expected::compute_expectations;
use methyltfr::gc::bin_means_indexed;
use methyltfr::pipeline::{InputFormat, Options, Sample};

/// Run one fixture through the whole pipeline and return
/// `(deviation, expected_deviation)`.
fn run_fixture(dir: &std::path::Path, label: &str) -> (f64, f64) {
    let (annotation, sites) = common::load_batf(dir);
    let bins = bin_means_indexed(&sites, &annotation.gc_windows, annotation.ignore_strand)
        .unwrap_or_else(|e| panic!("{label}: GC bins: {e}"));
    let motif = &annotation.motifs[0];
    let width = motif_width_of(&motif.tfbs).unwrap_or_else(|e| panic!("{label}: motif width: {e}"));
    let observed = compute_observed(
        &motif.tfbs,
        width,
        None,
        &sites,
        annotation.ignore_strand,
        &motif.name,
    )
    .unwrap_or_else(|e| panic!("{label}: observed profile: {e}"));
    let expected = compute_expectations(&motif.gcfreq, &bins)
        .unwrap_or_else(|e| panic!("{label}: expected profile: {e}"));
    let d = deviations(&observed, &expected);
    (
        d.deviation()
            .unwrap_or_else(|| panic!("{label}: deviation is NA")),
        d.expected
            .unwrap_or_else(|| panic!("{label}: exp_dev is NA")),
    )
}

fn oracle_dev(dir: &std::path::Path) -> (f64, f64) {
    let (header, rows) = common::read_columns(&dir.join("expected_tierB/expected_dev.tsv"));
    let dev = common::col(&header, &rows, "dev")[0]
        .parse::<f64>()
        .unwrap();
    let exp = common::col(&header, &rows, "exp_d")[0]
        .parse::<f64>()
        .unwrap();
    (dev, exp)
}

#[test]
fn golden_batf_matches_the_oracle_at_1e_10() {
    let dir = batf("");
    let (got_dev, got_exp) = run_fixture(&dir, "batf");
    let (want_dev, want_exp) = oracle_dev(&dir);
    let mut errors = Errors::new();
    errors.add("dev", got_dev, want_dev);
    errors.add("exp_dev", got_exp, want_exp);
    errors.finish("golden_batf");

    // Plan section 4's provisional column, now confirmed against Tier A.
    assert!((got_dev - 1.7474268).abs() < 5e-7, "dev = {got_dev}");
    assert!((got_exp - 0.9798459268).abs() < 5e-7, "exp_dev = {got_exp}");
}

#[test]
fn golden_batf_1d99721_reproduces_the_published_column() {
    let dir = batf_1d99721("");
    let (got_dev, got_exp) = run_fixture(&dir, "batf_1d99721");
    let (want_dev, want_exp) = oracle_dev(&dir);
    let mut errors = Errors::new();
    errors.add("dev", got_dev, want_dev);
    errors.add("exp_dev", got_exp, want_exp);
    errors.finish("golden_batf_1d99721");

    // The published pkgdown values, to the 7 significant digits they are printed
    // with. 5e-7 is half of the last printed digit.
    assert!((got_dev - 1.743674).abs() < 5e-7, "dev = {got_dev}");
    assert!((got_exp - 0.9835985).abs() < 5e-7, "exp_dev = {got_exp}");
}

#[test]
fn the_two_fixtures_disagree_exactly_where_the_plan_says() {
    // Bins 3, 4 and 5 agree between the two fixture revisions; bins 1 and 2 do
    // not, because `gcdist_subset.rda` changed between them while
    // `BATF_tf_bindsites.rda` did not. The observed deviation is identical.
    let (_, a) = common::load_batf(&batf(""));
    let (_, b) = common::load_batf(&batf_1d99721(""));
    let bins_a = bin_means_indexed(&a, &a_windows(), true).unwrap();
    let bins_b = bin_means_indexed(&b, &b_windows(), true).unwrap();
    assert_ne!(bins_a.mean[0], bins_b.mean[0], "bin 1 differs");
    assert_ne!(bins_a.mean[1], bins_b.mean[1], "bin 2 differs");
    for i in 2..5 {
        assert!(
            (bins_a.mean[i] - bins_b.mean[i]).abs() < 1e-15,
            "bin {} agrees: {} vs {}",
            i + 1,
            bins_a.mean[i],
            bins_b.mean[i]
        );
    }
    // The observed hits are the same 29 in both, so obs_d is the same.
    let dir_a = batf("");
    let dir_b = batf_1d99721("");
    assert!((oracle_dev(&dir_a).1 - oracle_dev(&dir_b).1).abs() > 1e-4);
}

fn a_windows() -> Vec<methyltfr::GcWindow> {
    common::load_batf(&batf("")).0.gc_windows.clone()
}

fn b_windows() -> Vec<methyltfr::GcWindow> {
    common::load_batf(&batf_1d99721("")).0.gc_windows.clone()
}

#[test]
fn the_run_level_pipeline_produces_the_same_row() {
    // The same thing through `pipeline::run`, with the portable sample format, so
    // the library entry point and the inline computation are cross-checked.
    let dir = batf("");
    let samples = vec![Sample {
        path: dir.join("msites.tsv"),
        format: InputFormat::Portable,
    }];
    let rows = methyltfr::pipeline::run(&samples, &dir, &Options::default()).unwrap();
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].sample, "msites.tsv");
    assert_eq!(rows[0].motif, "BATF");
    let (want_dev, want_exp) = oracle_dev(&dir);
    let mut errors = Errors::new();
    errors.add("dev", rows[0].deviation.unwrap(), want_dev);
    errors.add("exp_dev", rows[0].expected_deviation.unwrap(), want_exp);
    errors.finish("pipeline::run golden_batf");
}

#[test]
fn an_empty_motif_list_is_fatal() {
    let dir = batf("");
    let samples = vec![Sample {
        path: dir.join("msites.tsv"),
        format: InputFormat::Portable,
    }];
    let mut annotation_dir = std::env::temp_dir();
    annotation_dir.push(format!("methyltfr-empty-motifs-{}", std::process::id()));
    std::fs::create_dir_all(&annotation_dir).unwrap();
    // A manifest whose only motif has no binding sites.
    std::fs::write(
        annotation_dir.join("gc_windows.tsv"),
        "chr\tstart\tend\tstrand\tgc_bin\nchr1\t1\t30\t*\t1\n",
    )
    .unwrap();
    std::fs::write(
        annotation_dir.join("empty_tfbs.tsv"),
        "chr\tstart\tend\tstrand\n",
    )
    .unwrap();
    std::fs::write(annotation_dir.join("g.tsv"), "1\n").unwrap();
    std::fs::write(
        annotation_dir.join("motifs.tsv"),
        "motif\ttfbs_path\tgcfreq_path\nEMPTY\tempty_tfbs.tsv\tg.tsv\n",
    )
    .unwrap();
    let e = methyltfr::pipeline::run(&samples, &annotation_dir, &Options::default()).unwrap_err();
    assert_eq!(e.to_string(), "No valid motifs remaining after validation.");
    std::fs::remove_dir_all(&annotation_dir).ok();
}
