//! T34: the bismarkCytosine reader.
//!
//! Row: `chr pos strand meth_unmeth total context`. Coverage is `c4 + c5` and the
//! score is `c4 / coverage` rounded to six decimals, so a `0/0` row scores NaN
//! and is dropped by the first post-step -- not by the coverage threshold.
//!
//! Coordinates are width 1: `parse_bismarkcytosine` uses column 2 for both.

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("bismarkcytosine"), "bismarkcytosine", 1.0)
        .expect("read bismarkcytosine");
    compare_to_oracle("bismarkcytosine", "", &sites);
}

#[test]
fn the_zero_over_zero_row_is_dropped() {
    // Plan section 5: the bundled example has one 0/0 row that must be dropped.
    // It holds 7 data rows and 6 survive.
    let sites = read(&parser_example("bismarkcytosine"), "bismarkcytosine", 1.0).unwrap();
    assert_eq!(sites.len(), 6);
    assert!(sites.iter().all(|s| s.coverage > 0.0));
}

#[test]
fn the_nan_drop_is_independent_of_the_threshold() {
    // The row is dropped because its score is NaN, so lowering the threshold to 0
    // must not bring it back. `expected.thr0` is the Tier A file for exactly this.
    let sites = read(&parser_example("bismarkcytosine"), "bismarkcytosine", 0.0).unwrap();
    compare_to_oracle("bismarkcytosine", ".thr0", &sites);
    assert_eq!(sites.len(), 6);
}

#[test]
fn coverage_is_the_sum_of_the_two_count_columns() {
    let dir = scratch("parser-cytosine-cov");
    let body = "chr1\t10\t+\t1\t3\tCG\tCGG\n\
                chr1\t11\t-\t2\t2\tCG\tCGG\n";
    let p = common::write_file(&dir, "bismark.tsv", body);
    let sites = read(&p, "bismarkcytosine", 1.0).unwrap();
    assert_eq!(sites[0].coverage, 4.0);
    assert_eq!(sites[0].score.to_bits(), 0.25f64.to_bits());
    assert_eq!(sites[1].coverage, 4.0);
    assert_eq!(sites[1].score.to_bits(), 0.5f64.to_bits());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-cytosine-errors");
    assert_shared_errors(&dir, "bismarkcytosine", 5);
    std::fs::remove_dir_all(&dir).ok();
}
