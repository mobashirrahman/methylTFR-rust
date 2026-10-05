//! T35: the ENCODE bedMethyl reader.
//!
//! Row: bedMethyl, 11 columns, with coverage in `c10` and percent methylated in
//! `c11`. Score is `c11 / 100` rounded to six decimals, strand is `c6`.
//!
//! Two ENCODE-specific behaviours:
//!
//! * `c11` outside `0..=100` is an error, checked before anything else.
//! * The header is auto-detected. Upstream uses `fread(header = "auto")`; this
//!   port uses the single-field rule recorded as `docs/divergences.md` D1, namely
//!   "the first line is a header iff its 2nd or 3rd field is not an integer".

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("encode"), "encode", 1.0).expect("read encode");
    compare_to_oracle("encode", "", &sites);
    // Plan section 5: 6 records, these scores, these coverages.
    assert_eq!(sites.len(), 6);
    let scores = [0.06f64, 0.03, 0.0, 1.0, 0.55, 1.0];
    for (i, w) in scores.iter().enumerate() {
        assert_eq!(sites[i].score.to_bits(), (*w).to_bits(), "record {i} score");
    }
    let cov = [62.0f64, 62.0, 31.0, 5.0, 31.0, 10.0];
    for (i, w) in cov.iter().enumerate() {
        assert_eq!(
            sites[i].coverage.to_bits(),
            (*w).to_bits(),
            "record {i} coverage"
        );
    }
}

#[test]
fn a_coverage_threshold_of_20_keeps_four_records() {
    // Plan section 5: 4 records at cov_threshold = 20. The Tier A file for this
    // is `encode.expected.thr20.tsv`.
    let sites = read(&parser_example("encode"), "encode", 20.0).unwrap();
    compare_to_oracle("encode", ".thr20", &sites);
    assert_eq!(sites.len(), 4);
    assert_eq!(sites[3].start, 1000199, "the coverage-5 record is gone");
}

#[test]
fn the_header_is_auto_detected() {
    let dir = scratch("parser-encode-header");
    let row = |s: i64, e: i64, cov: i64, pct: i64| {
        format!("chr1\t{s}\t{e}\tname\t0\t+\t{s}\t{e}\t255,255,0\t{cov}\t{pct}\n")
    };
    // A header line is detected and skipped, leaving all three data rows.
    let with_header = format!(
        "chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tcoverage\tpercentMeth\n{}{}{}",
        row(10, 11, 5, 20),
        row(12, 13, 7, 100),
        row(14, 15, 3, 0)
    );
    let p = common::write_file(&dir, "headed.bedMethyl", &with_header);
    assert_eq!(read(&p, "encode", 1.0).unwrap().len(), 3);

    // The same file without the header line is still read as three records: the
    // first data row has integer start/end, so D1 does not call it a header.
    let headerless = format!(
        "{}{}{}",
        row(10, 11, 5, 20),
        row(12, 13, 7, 100),
        row(14, 15, 3, 0)
    );
    let p = common::write_file(&dir, "headerless.bedMethyl", &headerless);
    assert_eq!(read(&p, "encode", 1.0).unwrap().len(), 3);

    // D1's residual risk, made explicit: a first row with a non-integer coordinate
    // in a headerless file is treated as a header and dropped.
    let risky = format!(
        "chr1\t10.5\t11\tname\t0\t+\t10.5\t11\t255,255,0\t5\t20\n{}",
        row(12, 13, 7, 100)
    );
    let p = common::write_file(&dir, "risky.bedMethyl", &risky);
    let sites = read(&p, "encode", 1.0).unwrap();
    assert_eq!(sites.len(), 1, "the first row was read as a header");
    assert_eq!(sites[0].start, 12);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_percentage_outside_zero_to_hundred_is_an_error() {
    let dir = scratch("parser-encode-pct");
    let mk = |pct: &str| format!("chr1\t10\t11\tname\t0\t+\t10\t11\t255,255,0\t5\t{pct}\n");
    for pct in ["-1", "101", "1000"] {
        let p = common::write_file(&dir, "bad.bedMethyl", &mk(pct));
        let e = read(&p, "encode", 1.0).unwrap_err();
        assert!(e.to_string().contains("outside 0-100"), "pct = {pct}: {e}");
    }
    // 0 and 100 are both fine.
    for pct in ["0", "100"] {
        let p = common::write_file(&dir, "ok.bedMethyl", &mk(pct));
        assert_eq!(read(&p, "encode", 1.0).unwrap().len(), 1);
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-encode-errors");
    assert_shared_errors(&dir, "encode", 11);
    std::fs::remove_dir_all(&dir).ok();
}
