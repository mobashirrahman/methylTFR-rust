//! T41: the CLI, end to end, with `assert_cmd`.
//!
//! The BATF golden is driven through the CLI the way upstream's own example does:
//! `msites.tsv` is written back out as a `bismarkCov` file and read with
//! `--format bismarkcov`. That exercises the parser, the pipeline and the CSV
//! writer in one go, and it is the check that matters most for "drop-in
//! replacement": the same command produces the same deviation as R.
//!
//! Two samples x two motifs covers the ordering rule, because that is a
//! contract, not an implementation detail: rows come out by sample (argument
//! order) then motif (manifest order), whatever order the motifs were named in.

mod common;

use std::path::{Path, PathBuf};
use std::process::Command;

use assert_cmd::prelude::*;

/// `cargo test` puts the binary next to the test executable's directory.
fn binary() -> PathBuf {
    let mut p = std::env::current_exe().expect("test executable path");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    p.join("methyltfr")
}

fn cmd() -> Command {
    Command::cargo_bin("methyltfr").expect("methyltfr binary")
}

/// `tests/fixtures/batf/msites.tsv` as a `bismarkCov` file, which is how
/// upstream's own example feeds a methylome back in.
fn write_bismarkcov_sample(dir: &Path, name: &str) -> PathBuf {
    let (annotation, sites) = common::load_batf(&common::batf(""));
    let body: String = sites
        .iter()
        .map(|s| {
            // bismarkCov is `chr start end percent meth unmeth`, and the score is
            // `percent / 100` rounded to six decimals. The fixture scores are all
            // exact multiples of 1/1000, so percent = score * 1000 divides back.
            // percent = score * 100, written with full precision so that the
            // parser's `round(percent / 100, 6)` lands on exactly the score the
            // fixture holds, which is what makes the CLI's output comparable with
            // the oracle bit for bit.
            let pct = s.score * 100.0;
            let meth = s.score * s.coverage;
            let unmeth = s.coverage - meth;
            format!(
                "chr1\t{}\t{}\t{pct:?}\t{meth:.6}\t{unmeth:.6}",
                s.start, s.end
            )
        })
        .collect::<Vec<_>>()
        .join("\n");
    let p = dir.join(name);
    std::fs::write(&p, format!("{body}\n")).unwrap();
    let _ = &annotation;
    p
}

/// An annotation directory holding the BATF resources, optionally with a second
/// motif pointing at the same files under a different name.
fn write_annotation(dir: &Path, motifs: &[&str]) -> PathBuf {
    let src = common::batf("");
    for name in ["gc_windows.tsv", "tfbs.tsv.gz", "gcfreq.tsv"] {
        std::fs::copy(src.join(name), dir.join(name)).unwrap();
    }
    let manifest: String = motifs
        .iter()
        .map(|m| format!("{m}\ttfbs.tsv.gz\tgcfreq.tsv"))
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(
        dir.join("motifs.tsv"),
        format!("motif\ttfbs_path\tgcfreq_path\n{manifest}\n"),
    )
    .unwrap();
    dir.to_path_buf()
}

#[test]
fn batf_golden_through_the_cli() {
    let dir = common::scratch("cli-golden");
    let sample = write_bismarkcov_sample(&dir, "sample1.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);

    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "exit {:?}, stderr: {}",
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let want_dev = oracle("batf", "dev");
    let want_exp = oracle("batf", "exp_d");

    let mut lines = stdout.lines();
    assert_eq!(
        lines.next(),
        Some("sample,motif,deviation,expected_deviation")
    );
    let fields: Vec<&str> = lines.next().expect("one row").split(',').collect();
    assert!(lines.next().is_none(), "exactly one row");

    // The sample id is the file's basename, which is upstream's `sample_ids`.
    assert_eq!(fields[0], "sample1.bismarkCov");
    assert_eq!(fields[1], "BATF");

    // `deviation` is byte-identical to the oracle's. `expected_deviation` differs
    // in the last one or two ulp and cannot be otherwise: R's `%*%` goes through
    // BLAS and its `mean()` accumulates in long double, while this port sums in
    // f64. Section 3's tolerance is 1e-10, and the measured gap is ~2 ulp. See
    // docs/reference-version.md, "T04 oracle agreement".
    assert_eq!(
        fields[2],
        want_dev.as_str(),
        "deviation must be byte-identical to the oracle"
    );
    let mut errors = common::Errors::new();
    errors.add(
        "expected_deviation",
        fields[3].parse().unwrap(),
        want_exp.parse().unwrap(),
    );
    errors.finish("cli expected_deviation");
    let want_bits = want_exp.parse::<f64>().unwrap().to_bits();
    let got_bits = fields[3].parse::<f64>().unwrap().to_bits();
    let ulps = want_bits.abs_diff(got_bits);
    assert!(
        ulps <= 4,
        "expected_deviation differs by {ulps} ulp: {} vs {}",
        fields[3],
        want_exp
    );

    // And it must equal the golden column to the printed precision.
    let dev: f64 = want_dev.parse().unwrap();
    assert!((dev - 1.7474268).abs() < 5e-7, "dev = {dev}");
    let exp: f64 = fields[3].parse().unwrap();
    assert!((exp - 0.9798459268).abs() < 5e-7, "exp_dev = {exp}");
    std::fs::remove_dir_all(&dir).ok();
}

trait TapCreate {
    fn tap_create(&self) -> PathBuf;
}

impl TapCreate for Path {
    fn tap_create(&self) -> PathBuf {
        std::fs::create_dir_all(self).unwrap();
        self.to_path_buf()
    }
}

fn oracle(fixture: &str, column: &str) -> String {
    let path = if fixture == "batf" {
        common::batf("expected_tierB/expected_dev.tsv")
    } else {
        common::batf_1d99721("expected_tierB/expected_dev.tsv")
    };
    let (header, rows) = common::read_columns(&path);
    common::col(&header, &rows, column)[0].to_string()
}

#[test]
fn two_samples_and_two_motifs_are_ordered_by_sample_then_manifest() {
    let dir = common::scratch("cli-order");
    let a = write_bismarkcov_sample(&dir, "sample_a.bismarkCov");
    let b = write_bismarkcov_sample(&dir, "sample_b.bismarkCov");
    // Manifest order is BATF then BATF2; --motif is asked for in the other order.
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF", "BATF2"]);

    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg("--motif")
        .arg("BATF2")
        .arg("--motif")
        .arg("BATF")
        .arg(&a)
        .arg(&b)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let stdout = String::from_utf8(out).unwrap();
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines[0], "sample,motif,deviation,expected_deviation");
    assert_eq!(lines.len(), 5, "2 samples x 2 motifs plus the header");
    let names: Vec<&str> = lines[1..]
        .iter()
        .map(|l| l.split(',').next().unwrap())
        .collect();
    assert_eq!(
        names,
        vec![
            "sample_a.bismarkCov",
            "sample_a.bismarkCov",
            "sample_b.bismarkCov",
            "sample_b.bismarkCov",
        ]
    );
    let motifs: Vec<&str> = lines[1..]
        .iter()
        .map(|l| l.split(',').nth(1).unwrap())
        .collect();
    assert_eq!(
        motifs,
        vec!["BATF", "BATF2", "BATF", "BATF2"],
        "manifest order, not --motif argument order"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn two_runs_produce_byte_identical_output() {
    let dir = common::scratch("cli-stable");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    let run = || {
        cmd()
            .args(["run", "--format", "bismarkcov", "--annotation"])
            .arg(&ann)
            .arg(&sample)
            .output()
            .unwrap()
            .stdout
    };
    assert_eq!(run(), run());
    // And identical to the library path, which is the same code without the
    // process boundary.
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn output_can_go_to_a_file() {
    let dir = common::scratch("cli-out");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    let out_path = dir.join("out.csv");
    cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg("-o")
        .arg(&out_path)
        .arg(&sample)
        .assert()
        .success()
        .stdout("");
    let written = std::fs::read_to_string(&out_path).unwrap();
    assert!(written.starts_with("sample,motif,deviation,expected_deviation\n"));
    assert_eq!(written.lines().count(), 2);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_motif_with_no_hits_exits_non_zero() {
    let dir = common::scratch("cli-nohits");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann_dir = dir.join("ann").tap_create();
    // A second motif whose binding sites are nowhere near any site: the BATF TFBS
    // files with every coordinate shifted far past the end of the sample.
    write_annotation(&ann_dir, &["BATF"]);
    let mut tfbs = String::from("chr\tstart\tend\tstrand\n");
    let (ann, _) = common::load_batf(&common::batf(""));
    for r in ann.motifs[0].tfbs.iter().take(50) {
        tfbs.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            "chr1",
            r.start + 100_000_000,
            r.end + 100_000_000,
            r.strand.as_char()
        ));
    }
    std::fs::write(ann_dir.join("far_tfbs.tsv"), &tfbs).unwrap();
    std::fs::write(
        ann_dir.join("motifs.tsv"),
        "motif\ttfbs_path\tgcfreq_path\nFAR\tfar_tfbs.tsv\tgcfreq.tsv\n",
    )
    .unwrap();

    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann_dir)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(!out.status.success(), "must not succeed");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("No methylation sites found in the FAR binding sites"),
        "stderr: {stderr}"
    );
    assert!(out.stdout.is_empty(), "stdout must stay clean");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn keep_going_writes_na_instead_of_failing() {
    let dir = common::scratch("cli-keepgoing");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann_dir = dir.join("ann").tap_create();
    write_annotation(&ann_dir, &["BATF"]);
    let (far_src, _) = common::load_batf(&common::batf(""));
    let mut tfbs = String::from("chr\tstart\tend\tstrand\n");
    for r in far_src.motifs[0].tfbs.iter().take(50) {
        tfbs.push_str(&format!(
            "{}\t{}\t{}\t{}\n",
            "chr1",
            r.start + 100_000_000,
            r.end + 100_000_000,
            r.strand.as_char()
        ));
    }
    std::fs::write(ann_dir.join("far_tfbs.tsv"), &tfbs).unwrap();
    std::fs::write(
        ann_dir.join("motifs.tsv"),
        "motif\ttfbs_path\tgcfreq_path\nFAR\tfar_tfbs.tsv\tgcfreq.tsv\n",
    )
    .unwrap();

    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann_dir)
        .arg("--keep-going")
        .arg(&sample)
        .output()
        .unwrap();
    assert!(out.status.success(), "--keep-going must succeed");
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert_eq!(
        stdout,
        "sample,motif,deviation,expected_deviation\ns.bismarkCov,FAR,NA,NA\n"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("keep-going"), "stderr: {stderr}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_missing_sample_is_a_clean_error() {
    let dir = common::scratch("cli-missing");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg(dir.join("nope.bismarkCov"))
        .assert()
        .failure()
        .stderr(predicates::str::contains("nope.bismarkCov"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn an_unknown_format_is_rejected() {
    let dir = common::scratch("cli-format");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    cmd()
        .args(["run", "--format", "wiggle", "--annotation"])
        .arg(&ann)
        .arg(&sample)
        .assert()
        .failure()
        .stderr(predicates::str::contains("wiggle is not a valid file type"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn an_unknown_motif_is_rejected() {
    let dir = common::scratch("cli-motif");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg("--motif")
        .arg("NOSUCH")
        .arg(&sample)
        .assert()
        .failure()
        .stderr(predicates::str::contains("is not in the annotation"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn logs_never_reach_stdout() {
    let dir = common::scratch("cli-logs");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg(&sample)
        .output()
        .unwrap();
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert_eq!(stdout.lines().count(), 2, "stdout is CSV and nothing else");
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(
        stderr.contains("Processing s.bismarkCov"),
        "stderr: {stderr}"
    );
    assert!(stderr.contains("Finished processing s.bismarkCov"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn help_and_version_work() {
    cmd()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicates::str::contains("methyltfr"));
    cmd().arg("--version").assert().success();
}

#[test]
fn a_strand_aware_run_still_produces_the_golden_value() {
    // The fixture's sites and binding sites agree in strand often enough that the
    // BATF result should be close but not identical; what matters here is that the
    // flag is wired through and does not crash.
    let dir = common::scratch("cli-strand");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann = write_annotation(&dir.join("ann").tap_create(), &["BATF"]);
    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann)
        .arg("--strand-aware")
        .arg(&sample)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(String::from_utf8(out.stdout).unwrap().lines().count(), 2);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_gzipped_annotation_is_read() {
    // `gc_windows.tsv[.gz]` is allowed by the portable format spec.
    let dir = common::scratch("cli-gz");
    let sample = write_bismarkcov_sample(&dir, "s.bismarkCov");
    let ann_dir = dir.join("ann").tap_create();
    let src = common::batf("");
    std::fs::copy(src.join("tfbs.tsv.gz"), ann_dir.join("tfbs.tsv.gz")).unwrap();
    std::fs::copy(src.join("gcfreq.tsv"), ann_dir.join("gcfreq.tsv")).unwrap();
    let windows = std::fs::read(src.join("gc_windows.tsv")).unwrap();
    {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(ann_dir.join("gc_windows.tsv.gz")).unwrap(),
            flate2::Compression::default(),
        );
        enc.write_all(&windows).unwrap();
        enc.finish().unwrap();
    }
    std::fs::write(
        ann_dir.join("motifs.tsv"),
        "motif\ttfbs_path\tgcfreq_path\nBATF\ttfbs.tsv.gz\tgcfreq.tsv\n",
    )
    .unwrap();
    let out = cmd()
        .args(["run", "--format", "bismarkcov", "--annotation"])
        .arg(&ann_dir)
        .arg(&sample)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(stdout.contains("BATF"), "{stdout}");
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn the_binary_is_where_the_tests_expect_it() {
    assert!(binary().exists(), "{} not built", binary().display());
}
