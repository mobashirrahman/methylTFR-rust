//! T30: the bismarkCov reader.
//!
//! Row: `chr start end percent_methylated meth unmeth`. Strand is **always** `*`:
//! `parse_bismarkcov` overwrites the column with `rep("*", nrow)`, so a file with
//! a strand column is not one, and there is nothing to validate.
//!
//! Score is `c4 / 100` rounded to six decimals, coverage is `c5 + c6`.

mod common;

use common::{assert_shared_errors, compare_to_oracle, parser_example, read, scratch};

#[test]
fn matches_the_tier_a_oracle() {
    let sites = read(&parser_example("bismarkcov"), "bismarkcov", 1.0).expect("read bismarkcov");
    compare_to_oracle("bismarkcov", "", &sites);
    assert_eq!(sites.len(), 6);
}

#[test]
fn strand_is_always_star() {
    let dir = scratch("parser-cov-strand");
    // A seven-column file whose 7th column is nonsense must still come out with
    // strand `*`, because upstream never reads a strand column for this format.
    let body = "chr1\t10\t11\t50\t1\t1\textra\n\
                chr1\t12\t13\t0\t2\t2\textra\n";
    let p = common::write_file(&dir, "bismarkCov.tsv", body);
    let sites = read(&p, "bismarkcov", 1.0).unwrap();
    assert_eq!(sites.len(), 2);
    for s in &sites {
        assert_eq!(s.strand, methyltfr::Strand::Star);
    }
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn score_and_coverage() {
    let sites = read(&parser_example("bismarkcov"), "bismarkcov", 1.0).unwrap();
    let want = [1.0f64, 0.0, 1.0, 0.0, 1.0, 1.0];
    let cov = [1.0f64, 1.0, 1.0, 1.0, 1.0, 1.0];
    for i in 0..6 {
        assert_eq!(sites[i].score.to_bits(), want[i].to_bits(), "record {i}");
        assert_eq!(sites[i].coverage.to_bits(), cov[i].to_bits(), "record {i}");
        assert_eq!(sites[i].end - sites[i].start + 1, 2, "BED width is kept");
    }
}

#[test]
fn shared_error_cases() {
    let dir = scratch("parser-cov-errors");
    // bismarkcov's strand column is synthesized, so the bad-strand case does not
    // apply; every other shared case does.
    let path = parser_example("bismarkcov");
    assert!(read(&path, "bismarkcov", 1.0).is_ok());
    let _ = assert_shared_errors;
    // Too few columns.
    let short = common::truncate_last_row(&common::sample_rows("bismarkcov", 6), 5);
    let p = common::write_file(&dir, "short.tsv", &short);
    let e = read(&p, "bismarkcov", 1.0).unwrap_err();
    assert!(e.to_string().contains("at least 6 columns"), "{e}");
    // Non-numeric.
    let broken = common::sample_rows("bismarkcov", 6).replacen('\t', "\tnope", 1);
    let p = common::write_file(&dir, "nonnumeric.tsv", &broken);
    let e = read(&p, "bismarkcov", 1.0).unwrap_err();
    assert_eq!(e.line(), Some(1));
    assert!(
        e.to_string().contains("is not a number") || e.to_string().contains("is not an integer"),
        "{e}"
    );
    // gzip and plain agree.
    let plain = common::write_file(&dir, "p.tsv", &common::sample_rows("bismarkcov", 6));
    let gz = common::write_file(&dir, "p.tsv.gz", &common::sample_rows("bismarkcov", 6));
    let a = read(&plain, "bismarkcov", 1.0).unwrap();
    let b = read(&gz, "bismarkcov", 1.0).unwrap();
    assert_eq!(a.len(), b.len());
    assert_eq!(a[0].score.to_bits(), b[0].score.to_bits());
    std::fs::remove_dir_all(&dir).ok();
}
