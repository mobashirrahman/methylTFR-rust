//! The command line interface (`docs/AGENT_PLAN.md` T40).
//!
//! ```text
//! methyltfr run --format <type> --annotation <dir> [--motif NAME]...
//!              [--enhancer FILE] [--strand-aware] [--cov-threshold 1]
//!              [--keep-going] [-o out.csv] SAMPLE...
//! ```
//!
//! `<dir>` holds `gc_windows.tsv[.gz]` and `motifs.tsv`. Output is CSV on stdout
//! or to `-o`, with the header `sample,motif,deviation,expected_deviation`, rows
//! ordered by sample (argument order) then motif (manifest order), floats in
//! Rust's shortest round-trip form and missing values as `NA`. Log lines go to
//! stderr only, so stdout is always a clean CSV.

use std::io::Write;
use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand};

use crate::error::{Error, Result};
use crate::pipeline::{InputFormat, Options, Row, Sample, run};

/// Top-level command line.
#[derive(Parser, Debug)]
#[command(
    name = "methyltfr",
    version,
    about = "Quantification of DNA methylation signatures in TFBS",
    long_about = "Rust port of the Bioconductor package methylTFR 0.99.9. \
                  A drop-in replacement for the methylation-signature part of that \
                  package: the same inputs, the same numbers."
)]
pub struct Cli {
    /// The subcommand to run.
    #[command(subcommand)]
    pub command: Command,
}

/// Subcommands.
#[derive(Subcommand, Debug)]
pub enum Command {
    /// Compute deviations for every sample against one annotation directory.
    Run(RunArgs),
}

/// Arguments of `methyltfr run`.
#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Methylation file format: epp, bissnp, allc, bismarkcytosine,
    /// bismarkcov or encode.
    #[arg(long, value_name = "type")]
    pub format: String,

    /// Directory holding `gc_windows.tsv[.gz]` and `motifs.tsv`.
    #[arg(long, value_name = "dir")]
    pub annotation: PathBuf,

    /// Use only these motifs, in the manifest's order rather than the argument
    /// order. May be repeated; unknown names are an error.
    #[arg(long = "motif", value_name = "NAME")]
    pub motifs: Vec<String>,

    /// Restrict the analysis to sites in these regions.
    #[arg(long, value_name = "FILE")]
    pub enhancer: Option<PathBuf>,

    /// Set `ignoreStrand = false`: require compatible strands on every overlap.
    #[arg(long)]
    pub strand_aware: bool,

    /// Minimum coverage per site (`cov_threshold`, default 1).
    #[arg(long, default_value_t = 1.0, value_name = "N")]
    pub cov_threshold: f64,

    /// Write `NA` and carry on when one cell fails, instead of aborting the run.
    /// A Rust extension; off by default, see docs/divergences.md D3.
    #[arg(long)]
    pub keep_going: bool,

    /// Write the CSV here instead of stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Sample files, in the order their rows should appear.
    #[arg(value_name = "SAMPLE", required = true)]
    pub samples: Vec<PathBuf>,
}

/// The CSV header, exactly as T40 specifies.
pub const HEADER: &str = "sample,motif,deviation,expected_deviation";

impl RunArgs {
    /// Validate the arguments and build the run options.
    pub fn options(&self) -> Result<Options> {
        if !self.cov_threshold.is_finite() || self.cov_threshold < 0.0 {
            return Err(Error::run(format!(
                "{} is not a valid coverage threshold!",
                self.cov_threshold
            )));
        }
        for p in self.samples.iter().chain(self.enhancer.iter()) {
            if !p.exists() {
                return Err(Error::Missing(p.to_string_lossy().into_owned()));
            }
        }
        if !self.annotation.is_dir() {
            return Err(Error::Missing(format!(
                "{} (annotation directory)",
                self.annotation.to_string_lossy()
            )));
        }
        Ok(Options {
            ignore_strand: !self.strand_aware,
            cov_threshold: self.cov_threshold,
            keep_going: self.keep_going,
            enhancer: self.enhancer.clone(),
        })
    }

    /// The samples, in argument order, each tagged with the requested format.
    pub fn samples(&self) -> Vec<Sample> {
        self.samples
            .iter()
            .map(|path| Sample {
                path: path.clone(),
                format: InputFormat::Methylome(self.format.to_ascii_lowercase()),
            })
            .collect()
    }
}

/// Render the CSV body. Floats use Rust's shortest round-trip form, which is
/// exact: the shortest decimal that parses back to the same `f64`.
pub fn to_csv(rows: &[Row]) -> String {
    let mut out = String::with_capacity(64 + rows.len() * 48);
    out.push_str(HEADER);
    out.push('\n');
    for r in rows {
        out.push_str(&r.sample);
        out.push(',');
        out.push_str(&r.motif);
        out.push(',');
        out.push_str(&fmt_value(r.deviation));
        out.push(',');
        out.push_str(&fmt_value(r.expected_deviation));
        out.push('\n');
    }
    out
}

/// A `None` deviation is `NA`; a `Some` is written with `{:?}`, which is Rust's
/// shortest round-trip representation. `Inf` and `NaN` are written the way R
/// spells them, so the file stays readable next to upstream's output.
fn fmt_value(v: Option<f64>) -> String {
    match v {
        None => "NA".to_string(),
        Some(x) if x.is_nan() => "NaN".to_string(),
        Some(x) if x.is_infinite() => {
            if x > 0.0 {
                "Inf".to_string()
            } else {
                "-Inf".to_string()
            }
        }
        Some(x) => format!("{x:?}"),
    }
}

/// Run the CLI and return the CSV. Kept separate from `main` so the integration
/// tests can drive it without a subprocess.
pub fn run_cli(args: &RunArgs) -> Result<String> {
    let options = args.options()?;
    let rows = run(&args.samples(), &args.annotation, &options)?;
    Ok(to_csv(&select_motifs(
        rows,
        &args.motifs,
        &args.annotation,
    )?))
}

/// Apply `--motif NAME...`, keeping manifest order rather than argument order.
fn select_motifs(rows: Vec<Row>, wanted: &[String], _annotation: &Path) -> Result<Vec<Row>> {
    if wanted.is_empty() {
        return Ok(rows);
    }
    let mut seen: Vec<&str> = Vec::new();
    for w in wanted {
        if !rows.iter().any(|r| r.motif == *w) {
            return Err(Error::run(format!("motif '{w}' is not in the annotation")));
        }
        if !seen.contains(&w.as_str()) {
            seen.push(w);
        }
    }
    Ok(rows
        .into_iter()
        .filter(|r| seen.contains(&r.motif.as_str()))
        .collect())
}

/// Entry point used by `main`: run and write the CSV.
pub fn main_with(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Run(args) => {
            let csv = run_cli(&args)?;
            match &args.output {
                Some(path) => {
                    let mut f = std::fs::File::create(path).map_err(|e| Error::io(path, e))?;
                    f.write_all(csv.as_bytes())
                        .map_err(|e| Error::io(path, e))?;
                }
                None => {
                    let mut out = std::io::stdout().lock();
                    out.write_all(csv.as_bytes())
                        .map_err(|e| Error::io("<stdout>", e))?;
                    out.flush().map_err(|e| Error::io("<stdout>", e))?;
                }
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(sample: &str, motif: &str, d: Option<f64>, e: Option<f64>) -> Row {
        Row {
            sample: sample.to_string(),
            motif: motif.to_string(),
            deviation: d,
            expected_deviation: e,
        }
    }

    #[test]
    fn csv_shape_and_order() {
        let rows = vec![
            row("s1.bed", "A", Some(1.0), Some(0.5)),
            row("s1.bed", "B", None, Some(2.0)),
            row("s2.bed", "A", Some(3.5), Some(1.25)),
        ];
        assert_eq!(
            to_csv(&rows),
            "sample,motif,deviation,expected_deviation\n\
             s1.bed,A,1.0,0.5\n\
             s1.bed,B,NA,2.0\n\
             s2.bed,A,3.5,1.25\n"
        );
    }

    #[test]
    fn missing_values_are_na() {
        assert_eq!(fmt_value(None), "NA");
        assert_eq!(fmt_value(Some(f64::NAN)), "NaN");
        assert_eq!(fmt_value(Some(f64::INFINITY)), "Inf");
        assert_eq!(fmt_value(Some(f64::NEG_INFINITY)), "-Inf");
    }

    #[test]
    fn floats_round_trip() {
        let v = 0.9798459267795764_f64;
        let s = fmt_value(Some(v));
        assert_eq!(s.parse::<f64>().unwrap().to_bits(), v.to_bits());
    }

    #[test]
    fn motif_filter_keeps_manifest_order() {
        let rows = vec![
            row("s", "A", Some(1.0), Some(1.0)),
            row("s", "B", Some(2.0), Some(2.0)),
            row("s", "C", Some(3.0), Some(3.0)),
        ];
        let out = select_motifs(rows, &["C".into(), "A".into()], Path::new(".")).unwrap();
        assert_eq!(
            out.iter().map(|r| r.motif.as_str()).collect::<Vec<_>>(),
            vec!["A", "C"]
        );
        let e = select_motifs(vec![], &["Z".into()], Path::new(".")).unwrap_err();
        assert!(e.to_string().contains("is not in the annotation"));
    }
}
