//! T33: the allc reader.
//!
//! Row: `chr pos strand context meth cov`. **The context column is not a filter**:
//! upstream keeps every context, and so does this port (see
//! `docs/divergences.md`, "Not divergences").
//!
//! Coordinates are width 1 because `parse_allc` uses column 2 for both `start` and
//! `end`. Coverage 0 divides to NaN, which the first post-step drops.

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("allc"), "allc", 1.0).expect("read allc");
    compare_to_oracle("allc", "", &sites);
    // Plan section 5: the allc example gives 3 records.
    assert_eq!(sites.len(), 3);
}

#[test]
fn every_context_is_kept() {
    // The bundled example mixes CGA and CGT contexts. Filtering on the context
    // would drop two of the three records.
    let sites = read(&parser_example("allc"), "allc", 1.0).unwrap();
    assert_eq!(sites.len(), 3);
    let strands: Vec<char> = sites.iter().map(|s| s.strand.as_char()).collect();
    assert_eq!(strands, vec!['+', '-', '+']);
}

#[test]
fn sites_are_width_one() {
    let sites = read(&parser_example("allc"), "allc", 1.0).unwrap();
    for s in &sites {
        assert_eq!(s.start, s.end, "allc width is 1");
    }
    assert_eq!(sites[0].start, 18283342);
}

#[test]
fn zero_coverage_is_dropped_as_nan() {
    let dir = scratch("parser-allc-nan");
    let body = "chr1\t10\t+\tCGT\t1\t2\n\
                chr1\t11\t-\tCGA\t0\t0\n\
                chr1\t12\t+\tCGG\t3\t4\n";
    let p = common::write_file(&dir, "allc.tsv", body);
    let sites = read(&p, "allc", 0.0).unwrap();
    assert_eq!(sites.len(), 2, "the 0/0 row is NaN and is dropped");
    assert_eq!(sites[0].start, 10);
    assert_eq!(sites[1].start, 12);
    // And 1/2 and 3/4 to six decimals.
    assert_eq!(sites[0].score.to_bits(), 0.5f64.to_bits());
    assert_eq!(sites[1].score.to_bits(), 0.75f64.to_bits());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn threshold_zero_matches_the_oracle() {
    let sites = read(&parser_example("allc"), "allc", 0.0).unwrap();
    compare_to_oracle("allc", ".thr0", &sites);
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-allc-errors");
    assert_shared_errors(&dir, "allc", 6);
    std::fs::remove_dir_all(&dir).ok();
}
