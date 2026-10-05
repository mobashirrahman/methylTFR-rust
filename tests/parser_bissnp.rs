//! T32: the BisSNP reader.
//!
//! Row: `chr start end percent_methylated coverage strand ...`. Score is
//! `c4 / 100` rounded to six decimals, coverage is `c5`, strand is `c6`.
//!
//! The quirk worth stating: `parse_bissnp` passes `skip = 1` to `fread`, so the
//! **first line is always skipped whether or not it is a header**. A headerless
//! BisSNP file therefore loses its first record. That is upstream's behaviour and
//! this port reproduces it.

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("bissnp"), "bissnp", 1.0).expect("read bissnp");
    compare_to_oracle("bissnp", "", &sites);
    // Plan section 5 does not name a count for this one; the example has 4 data
    // rows after the header.
    assert_eq!(sites.len(), 4);
}

#[test]
fn the_first_line_is_always_skipped() {
    let dir = scratch("parser-bissnp-skip");
    // No header at all: upstream still drops line 1.
    let headerless = "chr1\t10\t11\t50\t4\t+\n\
         chr1\t12\t13\t100\t8\t-\n\
         chr1\t14\t15\t0\t2\t+\n";
    let p = common::write_file(&dir, "headerless.tsv", headerless);
    let sites = read(&p, "bissnp", 1.0).unwrap();
    assert_eq!(
        sites.len(),
        2,
        "the first data row is dropped as if it were a header"
    );
    assert_eq!(sites[0].start, 12);
    assert_eq!(sites[1].start, 14);

    // With a header, all three data rows survive.
    let with_header = "chr\tstart\tend\tscore\tcoverage\tstrand\n\
         chr1\t10\t11\t50\t4\t+\n\
         chr1\t12\t13\t100\t8\t-\n\
         chr1\t14\t15\t0\t2\t+\n";
    let p = common::write_file(&dir, "headed.tsv", with_header);
    assert_eq!(read(&p, "bissnp", 1.0).unwrap().len(), 3);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn score_is_percent_over_one_hundred() {
    let sites = read(&parser_example("bissnp"), "bissnp", 1.0).unwrap();
    // 79.69, 90.62, 58.7, 50 over 100.
    let want = [0.7969f64, 0.9062, 0.587, 0.5];
    for (i, w) in want.iter().enumerate() {
        assert_eq!(sites[i].score.to_bits(), (*w).to_bits(), "record {i}");
    }
    let cov = [64.0f64, 64.0, 46.0, 4.0];
    for (i, w) in cov.iter().enumerate() {
        assert_eq!(sites[i].coverage.to_bits(), (*w).to_bits(), "record {i}");
    }
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-bissnp-errors");
    assert_shared_errors(&dir, "bissnp", 6);
    std::fs::remove_dir_all(&dir).ok();
}
