//! T51: every committed differential case, against the real methylTFR outputs.
//!
//! `tests/fixtures/differential/case_NNN/` holds a complete portable annotation,
//! one methylome, the case's settings, and either:
//!
//! * `expected/expected_bins.tsv`, `expected/observed_profile.tsv`,
//!   `expected/expected_profile.tsv` and `expected/expected_dev.tsv` -- produced by
//!   methylTFR 0.99.9 itself through `scripts/gen_differential_cases.R`; or
//! * `error.txt` -- the message upstream raised, for the cases that make it raise.
//!
//! Every case is checked:
//!
//! * the GC bin labels and hit counts must match **exactly**, since the mean is
//!   over hits and a miscounted hit is a wrong answer;
//! * the profiles are compared at the section 3 tolerance. The observed profile is
//!   compared as a multiset of `(x, value)` pairs, because its order is a
//!   consequence of the hit order, not part of the result;
//! * the three scalars must match, and `NA` must appear exactly where R's does;
//! * a case with `error.txt` must fail, with a message that agrees.
//!
//! The report at the end names the worst case and the worst absolute and relative
//! difference, which is what T51's Done-when asks for.

mod common;

use std::path::{Path, PathBuf};

use common::{Errors, read_columns, read_text, resolve_case_file};
use methyltfr::deviation::dev_helper;
use methyltfr::pipeline::{InputFormat, Options, Sample, load_annotation, run_case};

fn case_root() -> PathBuf {
    common::fixture("differential")
}

/// Every committed case directory, in name order.
fn cases() -> Vec<PathBuf> {
    let root = case_root();
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(&root) {
        for e in entries.flatten() {
            let p = e.path();
            if p.is_dir()
                && p.file_name()
                    .and_then(|s| s.to_str())
                    .is_some_and(|n| n.starts_with("case_"))
            {
                out.push(p);
            }
        }
    }
    out.sort();
    out
}

struct CaseResult {
    abs: f64,
    rel: f64,
    compared: usize,
}

fn run_case_dir(dir: &Path) -> CaseResult {
    let name = dir.file_name().unwrap().to_string_lossy().into_owned();
    let ignore_strand = read_ignore_strand(dir);

    // `options.tsv` carries ignoreStrand; enhancer presence is the enhancer.tsv
    // file itself.
    // An absent enhancer file is the signal that the case has none.
    let enhancer = {
        let plain = dir.join("enhancer.tsv");
        let gz = dir.join("enhancer.tsv.gz");
        if plain.exists() {
            Some(plain)
        } else if gz.exists() {
            Some(gz)
        } else {
            None
        }
    };
    let options = Options {
        ignore_strand,
        cov_threshold: 0.0,
        keep_going: false,
        enhancer,
    };

    let (annotation, mut chroms) =
        load_annotation(dir, &options).unwrap_or_else(|e| panic!("{name}: annotation: {e}"));
    assert_eq!(
        annotation.motifs.len(),
        1,
        "{name}: the generator writes one motif per case"
    );
    let sample = Sample {
        path: resolve_case_file(dir, "msites.tsv"),
        format: InputFormat::Portable,
    };

    let expected_error: Option<String> = {
        let p = dir.join("error.txt");
        if p.exists() {
            Some(read_text(&p).trim().to_string())
        } else {
            None
        }
    };

    let got = run_case(&sample, &annotation, &options, &mut chroms);

    match (expected_error, &got) {
        (Some(want), Ok(c)) => panic!(
            "{name}: upstream raises \"{want}\" but the port succeeded with \
             dev = {:?}",
            c.deviation.deviation()
        ),
        (None, Err(e)) => panic!("{name}: upstream succeeds but the port fails: {e}"),
        (Some(want), Err(e)) => {
            let message = e.to_string();
            // Both fail, and for the same reason. Upstream's message is compared
            // verbatim except for the bin-count mismatch, where upstream only
            // surfaces R's opaque `%*%` failure and this port names the cause; that
            // is recorded in docs/divergences.md under "Not divergences".
            if want == "non-conformable arguments" {
                assert!(
                    message.contains("non-conformable"),
                    "{name}: upstream says \"{want}\" but the port says \"{message}\""
                );
            } else {
                assert_eq!(
                    message, want,
                    "{name}: both fail, but with different messages"
                );
            }
            CaseResult {
                abs: 0.0,
                rel: 0.0,
                compared: 0,
            }
        }
        (None, Ok(out)) => compare_case(dir, &name, out),
    }
}

fn compare_case(dir: &Path, name: &str, out: &methyltfr::CaseOutput) -> CaseResult {
    let expected = dir.join("expected");
    assert!(
        expected.is_dir(),
        "{name}: neither expected/ nor error.txt is present"
    );
    let mut errors = Errors::new();

    // --- GC bins: labels and hit counts exact, means at the tolerance ---
    let (h, rows) = read_columns(&resolve_case_file(&expected, "expected_bins.tsv"));
    let want_bin = common::col(&h, &rows, "gc_bin");
    let want_mean = common::col(&h, &rows, "mean");
    let want_hits = common::col(&h, &rows, "n_hits");
    assert_eq!(
        out.bins.len(),
        want_bin.len(),
        "{name}: {} populated bins vs {} in the oracle",
        out.bins.len(),
        want_bin.len()
    );
    for i in 0..rows.len() {
        assert_eq!(
            out.bins.bin[i].to_string(),
            want_bin[i],
            "{name}: bin label at position {i}"
        );
        assert_eq!(
            out.bins.n_hits[i].to_string(),
            want_hits[i],
            "{name}: n_hits for bin {}",
            want_bin[i]
        );
        errors.add(
            &format!("{name} bin {}", want_bin[i]),
            out.bins.mean[i],
            want_mean[i].parse().unwrap(),
        );
    }

    // --- observed profile: a multiset of (x, value) ---
    let (h, rows) = read_columns(&resolve_case_file(&expected, "observed_profile.tsv"));
    let want_x = common::col(&h, &rows, "x");
    let want_value = common::col(&h, &rows, "value");
    assert_eq!(
        out.observed.len(),
        want_x.len(),
        "{name}: {} observed hits vs {} in the oracle",
        out.observed.len(),
        want_x.len()
    );
    let mut got_pairs: Vec<(i64, u64)> = out
        .observed
        .x
        .iter()
        .zip(out.observed.value.iter())
        .map(|(x, v)| (*x as i64, v.to_bits()))
        .collect();
    let mut want_pairs: Vec<(i64, u64)> = want_x
        .iter()
        .zip(want_value.iter())
        .map(|(x, v)| {
            (
                x.parse::<i64>().unwrap(),
                v.parse::<f64>().unwrap().to_bits(),
            )
        })
        .collect();
    got_pairs.sort_unstable();
    want_pairs.sort_unstable();
    assert_eq!(
        got_pairs, want_pairs,
        "{name}: observed profiles differ as (x, value) multisets"
    );

    // --- expected profile: x exactly, value at the tolerance ---
    let (h, rows) = read_columns(&resolve_case_file(&expected, "expected_profile.tsv"));
    let want_x = common::col(&h, &rows, "x");
    let want_value = common::col(&h, &rows, "value");
    assert_eq!(
        out.expected.len(),
        want_x.len(),
        "{name}: expected profile length"
    );
    for i in 0..rows.len() {
        errors.add(
            &format!("{name} expected x[{i}]"),
            out.expected.x[i],
            want_x[i].parse().unwrap(),
        );
        errors.add(
            &format!("{name} expected value[{i}]"),
            out.expected.value[i],
            want_value[i].parse().unwrap(),
        );
    }

    // --- the three scalars ---
    let (h, rows) = read_columns(&resolve_case_file(&expected, "expected_dev.tsv"));
    let want_obs = common::col(&h, &rows, "obs_d")[0].to_string();
    let want_exp = common::col(&h, &rows, "exp_d")[0].to_string();
    let want_dev = common::col(&h, &rows, "dev")[0].to_string();
    assert_eq!(
        out.motif,
        common::col(&h, &rows, "motif")[0],
        "{name}: motif"
    );

    for (label, got, want) in [
        ("obs_d", out.deviation.observed, want_obs.as_str()),
        ("exp_d", out.deviation.expected, want_exp.as_str()),
        ("dev", out.deviation.deviation(), want_dev.as_str()),
    ] {
        match (got, want == "NA") {
            (None, true) => {}
            (None, false) => panic!("{name}: {label} is NA but the oracle says {want}"),
            (Some(g), true) => panic!("{name}: {label} is {g} but the oracle says NA"),
            (Some(g), false) => errors.add(&format!("{name} {label}"), g, want.parse().unwrap()),
        }
    }

    // The two deviations must be internally consistent too: dev = obs - exp.
    if let (Some(o), Some(e), Some(d)) = (
        out.deviation.observed,
        out.deviation.expected,
        out.deviation.deviation(),
    ) {
        errors.add(&format!("{name} dev == obs_d - exp_d"), d, o - e);
        // And dev_helper on the observed profile must be what we reported.
        errors.add(
            &format!("{name} obs_d == dev_helper(observed)"),
            dev_helper(&out.observed).unwrap_or(f64::NAN),
            o,
        );
    }

    (errors.worst_abs, errors.worst_rel, errors.count).into()
}

impl From<(f64, f64, usize)> for CaseResult {
    fn from(v: (f64, f64, usize)) -> Self {
        CaseResult {
            abs: v.0,
            rel: v.1,
            compared: v.2,
        }
    }
}

fn read_ignore_strand(dir: &Path) -> bool {
    let p = dir.join("options.tsv");
    if !p.exists() {
        return true;
    }
    let text = read_text(&p);
    let lines: Vec<&str> = text
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .collect();
    let value = if lines.first() == Some(&"ignore_strand") {
        lines.get(1)
    } else {
        lines.first()
    };
    match value {
        Some(v) => v.eq_ignore_ascii_case("true") || *v == "1" || *v == "TRUE",
        None => true,
    }
}

#[test]
fn every_committed_case_matches_r() {
    let dirs = cases();
    assert!(
        dirs.len() >= 200,
        "only {} cases committed, T50 asks for 200",
        dirs.len()
    );
    let mut worst = CaseResult {
        abs: 0.0,
        rel: 0.0,
        compared: 0,
    };
    let mut worst_case = String::new();
    let mut n_errors = 0usize;
    let mut n_enhancer = 0usize;
    let mut n_strand_aware = 0usize;
    let mut total = 0usize;

    for dir in &dirs {
        let r = run_case_dir(dir);
        total += r.compared;
        if dir.join("error.txt").exists() {
            n_errors += 1;
        }
        if dir.join("enhancer.tsv").exists() || dir.join("enhancer.tsv.gz").exists() {
            n_enhancer += 1;
        }
        if !read_ignore_strand(dir) {
            n_strand_aware += 1;
        }
        if r.abs > worst.abs {
            worst = CaseResult {
                abs: r.abs,
                rel: r.rel,
                compared: r.compared,
            };
            worst_case = dir.file_name().unwrap().to_string_lossy().into_owned();
        }
    }

    // Coverage of the paths that matter: a differential suite where every case
    // succeeds tests far less than one where the error cases are present.
    assert!(n_errors > 0, "no case exercises an upstream error");
    assert!(n_enhancer > 0, "no case has an enhancer");
    assert!(
        n_strand_aware > 0 && n_strand_aware < dirs.len(),
        "ignoreStrand must be false in some cases and true in others"
    );

    eprintln!(
        "{} cases: {n_errors} raise an upstream error, {n_enhancer} have an \
         enhancer, {n_strand_aware} are strand-aware",
        dirs.len()
    );
    eprintln!("{total} values compared");
    eprintln!(
        "worst case {worst_case}: worst abs {:.3e}, worst rel {:.3e}, tolerance 1e-10",
        worst.abs, worst.rel
    );
}
