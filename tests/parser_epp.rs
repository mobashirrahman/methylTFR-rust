//! T31: the EPP reader.
//!
//! EPP is the format upstream's own vignette reads, so it is the one where a
//! regression is most visible.
//!
//! Row: `chr start end mc/cov score_per_mille strand`. Score is `c5 / 1000`
//! rounded to six decimals, coverage is the denominator of `c4` with any `'`
//! removed, and strand is `c6`.

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("epp"), "epp", 1.0).expect("read epp");
    compare_to_oracle("epp", "", &sites);
    // Plan section 5: the bundled EPP example has 6 records.
    assert_eq!(sites.len(), 6);
}

#[test]
fn coordinates_are_copied_verbatim() {
    // BED-style start/end means width 2. Upstream does not convert, so neither
    // do we; `chr1 10 11` stays 10..11 and not 9..10.
    let sites = read(&parser_example("epp"), "epp", 1.0).unwrap();
    assert_eq!(sites[0].start, 3010957);
    assert_eq!(sites[0].end, 3010958);
    assert_eq!(sites[0].end - sites[0].start + 1, 2);
}

#[test]
fn scores_are_per_mille_over_one_thousand() {
    let sites = read(&parser_example("epp"), "epp", 1.0).unwrap();
    // 1000/1000, 500/1000, 1000/1000, 500/1000, 814/1000, 500/1000.
    let want = [1.0f64, 0.5, 1.0, 0.5, 0.814, 0.5];
    for (i, w) in want.iter().enumerate() {
        assert_eq!(sites[i].score.to_bits(), (*w).to_bits(), "record {i}");
    }
    // Coverage is the denominator of column 4.
    let cov = [27.0f64, 7.0, 20.0, 20.0, 70.0, 100.0];
    for (i, w) in cov.iter().enumerate() {
        assert_eq!(sites[i].coverage.to_bits(), (*w).to_bits(), "record {i}");
    }
}

#[test]
fn thousands_separators_in_the_coverage_are_removed() {
    // Upstream runs `str_replace(mcov, "'", "")` over the split column, so
    // `1'000/2'000` gives coverage 2000.
    let dir = scratch("parser-epp-quote");
    let p = common::write_file(&dir, "epp.tsv", "chr1\t10\t11\t1'000/2'000\t1000\t+\n");
    let sites = read(&p, "epp", 1.0).unwrap();
    assert_eq!(sites.len(), 1);
    assert_eq!(sites[0].coverage, 2000.0);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn coverage_threshold_filters_and_thresholds_may_be_fractional() {
    let dir = scratch("parser-epp-threshold");
    let body = "chr1\t10\t11\t1/2\t500\t+\n\
                chr1\t12\t13\t3/4\t1000\t-\n\
                chr1\t14\t15\t9/10\t1000\t+\n";
    let p = common::write_file(&dir, "epp.tsv", body);
    // Coverages are the denominators: 2, 4 and 10.
    assert_eq!(
        read(&p, "epp", 0.0).unwrap().len(),
        3,
        "threshold 0 keeps all"
    );
    assert_eq!(
        read(&p, "epp", 1.0).unwrap().len(),
        3,
        "default threshold 1"
    );
    assert_eq!(read(&p, "epp", 2.0).unwrap().len(), 3, "2 >= 2");
    assert_eq!(read(&p, "epp", 2.5).unwrap().len(), 2, "2 < 2.5 drops 1/2");
    assert_eq!(read(&p, "epp", 5.0).unwrap().len(), 1, "5 > 4 drops 3/4");
    assert_eq!(read(&p, "epp", 10.0).unwrap().len(), 1);
    assert_eq!(read(&p, "epp", 11.0).unwrap().len(), 0);
    // The threshold is an f64 comparison, never an integer one.
    assert_eq!(read(&p, "epp", 2.000_001).unwrap().len(), 2);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_negative_coverage_threshold_is_rejected() {
    let e = read(&parser_example("epp"), "epp", -1.0).unwrap_err();
    assert!(
        e.to_string().contains("not a valid coverage threshold"),
        "{e}"
    );
}

#[test]
fn an_unknown_type_is_rejected() {
    let e = read(&parser_example("epp"), "wiggle", 1.0).unwrap_err();
    assert!(
        e.to_string().contains("wiggle is not a valid file type"),
        "{e}"
    );
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-epp-errors");
    assert_shared_errors(&dir, "epp", 6);
    std::fs::remove_dir_all(&dir).ok();
}
