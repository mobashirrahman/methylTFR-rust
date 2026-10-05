//! The run (`docs/AGENT_PLAN.md` section 2.8 run-level rules, T23).
//!
//! One pass per sample, then one pass per motif, sequentially, in manifest order.
//! Upstream splits the motifs into chunks and hands each chunk to
//! `BiocParallel`; with `threads = 1` that is a `SerialParam`, so the sequential
//! order here *is* the order upstream produces by default.
//!
//! Two rules that are easy to miss:
//!
//! * When an enhancer is supplied it reduces the **GC windows** once, before any
//!   bin mean is computed, *and* filters the resized TFBS. Both happen, so the
//!   bin means themselves change with the enhancer.
//! * Upstream aborts the whole run on a sample with no GC hits, a motif with no
//!   hits inside its binding sites, or a bin-count mismatch. The default here
//!   does the same. `--keep-going` is a Rust extension, off by default, that
//!   writes `NA` for that cell and logs why; see `docs/divergences.md` D3.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use rayon::prelude::*;

use crate::deviation::{
    Deviation, compute_observed, compute_observed_swept, deviations, motif_width_of,
};
use crate::error::{Error, Result};
use crate::expected::{Profile, compute_expectations};
use crate::gc::{GcBins, bin_means_swept, reduce_windows_by_enhancer};
use crate::intervals::StartIndex;
use crate::io::methylome::read_methylome;
use crate::io::portable::{read_annotation, read_msites};
use crate::model::{Annotation, ChromTable, Methylome, MotifAnnotation, Site, SortedMethylome};

/// How one sample file should be read.
#[derive(Clone, Debug)]
pub enum InputFormat {
    /// `msites.tsv`, the portable form used by the fixtures and the
    /// differential cases.
    Portable,
    /// One of the six methylation formats `read_methylome()` accepts.
    Methylome(String),
}

/// Run-level knobs.
#[derive(Clone, Debug)]
pub struct Options {
    /// `ignoreStrand`, default `true` upstream and here.
    pub ignore_strand: bool,
    /// `cov_threshold`, default 1.
    pub cov_threshold: f64,
    /// Rust extension: write `NA` and keep going instead of aborting.
    pub keep_going: bool,
    pub enhancer: Option<PathBuf>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            ignore_strand: true,
            cov_threshold: 1.0,
            keep_going: false,
            enhancer: None,
        }
    }
}

/// One output row. `None` in a numeric column is written as `NA`.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub sample: String,
    pub motif: String,
    pub deviation: Option<f64>,
    pub expected_deviation: Option<f64>,
}

/// One sample to process. Its identifier is the file's basename, which is what
/// upstream uses (`sample_ids = basename(files_list)`).
#[derive(Clone, Debug)]
pub struct Sample {
    pub path: PathBuf,
    pub format: InputFormat,
}

impl Sample {
    pub fn id(&self) -> String {
        sample_id(&self.path)
    }
}

/// `basename()` as R defines it: the last path component.
pub fn sample_id(path: &Path) -> String {
    path.file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// Where log lines go: stderr only, never stdout (T40).
pub fn log(msg: impl AsRef<str>) {
    eprintln!("{}", msg.as_ref());
}

/// Load the annotation, then apply the enhancer reduction to the GC windows
/// exactly once (`AGENT_PLAN.md` section 2.6).
///
/// The chromosome table is returned alongside the annotation and must be the
/// **same** table every sample is read into. Chromosome identity in this crate is
/// a `u32` from that one table, so two tables would silently compare unrelated
/// ids as if they were chromosomes.
pub fn load_annotation(dir: &Path, options: &Options) -> Result<(Annotation, ChromTable)> {
    let mut chroms = ChromTable::new();
    let mut annotation = read_annotation(
        dir,
        options.enhancer.as_deref(),
        options.ignore_strand,
        &mut chroms,
    )?;
    if let Some(enhancer) = &annotation.enhancer {
        let reduced =
            reduce_windows_by_enhancer(&annotation.gc_windows, enhancer, options.ignore_strand);
        log(format!(
            "enhancer: {} of {} GC windows retained",
            reduced.len(),
            annotation.gc_windows.len()
        ));
        annotation.gc_windows = reduced;
    }
    Ok((annotation, chroms))
}

/// Drop motifs with no binding sites or no matrix, as `valid_core_motifs` does.
///
/// A manifest row always carries a matrix on our side, so only the empty-TFBS
/// case can fire here, and an annotation with nothing left is fatal.
pub fn valid_motifs(annotation: &Annotation) -> Result<Vec<usize>> {
    let valid: Vec<usize> = (0..annotation.motifs.len())
        .filter(|&i| !annotation.motifs[i].tfbs.is_empty())
        .collect();
    if valid.len() < annotation.motifs.len() {
        log(format!(
            "Discarding {} motifs due to empty TFBS or missing matrix.",
            annotation.motifs.len() - valid.len()
        ));
    }
    if valid.is_empty() {
        return Err(Error::run("No valid motifs remaining after validation."));
    }
    Ok(valid)
}

/// Everything one sample against one annotation produces, in the shape the
/// reference oracles write. The CLI needs only [`Row`]; the differential tests
/// need the intermediate profiles to compare against `expected_bins.tsv`,
/// `observed_profile.tsv` and `expected_profile.tsv`.
#[derive(Clone, Debug)]
pub struct CaseOutput {
    pub sample: String,
    pub motif: String,
    pub bins: GcBins,
    pub observed: Profile,
    pub expected: Profile,
    pub deviation: Deviation,
}

/// Run one sample against a single-motif annotation and return every intermediate,
/// or upstream's error.
pub fn run_case(
    sample: &Sample,
    annotation: &Annotation,
    options: &Options,
    chroms: &mut ChromTable,
) -> Result<CaseOutput> {
    let methylome = load_sample(sample, options, chroms)?;
    // Sort once per sample, here, rather than in each consumer.
    let sorted = SortedMethylome::new(&methylome.sites);
    let bins = gc_bins(&methylome, annotation, options)?;
    let motif = &annotation.motifs[0];
    let width = motif_width_of(&motif.tfbs)?;
    let observed = compute_observed_swept(
        &motif.tfbs,
        width,
        annotation.enhancer.as_deref(),
        &sorted,
        options.ignore_strand,
        &motif.name,
    )?;
    let expected = compute_expectations(&motif.gcfreq, &bins)?;
    let deviation = deviations(&observed, &expected);
    Ok(CaseOutput {
        sample: sample.id(),
        motif: motif.name.clone(),
        bins,
        observed,
        expected,
        deviation,
    })
}

/// Run every sample against the annotation, sequentially.
///
/// Rows come out ordered by sample (argument order) then motif (manifest order).
pub fn run(samples: &[Sample], annotation_dir: &Path, options: &Options) -> Result<Vec<Row>> {
    let (annotation, mut chroms) = load_annotation(annotation_dir, options)?;
    let motifs = valid_motifs(&annotation)?;
    let mut rows = Vec::with_capacity(samples.len() * motifs.len());

    for sample in samples {
        log(format!("Processing {}", sample.id()));
        let methylome = load_sample(sample, options, &mut chroms)?;
        let bins = match gc_bins(&methylome, &annotation, options) {
            Ok(b) => b,
            Err(e) if options.keep_going => {
                log(format!("keep-going: {}: {e}", sample.id()));
                for &i in &motifs {
                    rows.push(Row {
                        sample: sample.id(),
                        motif: annotation.motifs[i].name.clone(),
                        deviation: None,
                        expected_deviation: None,
                    });
                }
                continue;
            }
            Err(e) => return Err(e),
        };

        // Over motifs, never inside one (M7 T73). Two properties make this safe:
        //
        // * `compute_cell` reads only `methylome`, `annotation` and `bins`, and
        //   writes to nothing shared, so the cells are independent.
        // * `par_iter().map(..).collect()` preserves the input order, so the rows
        //   come out in manifest order regardless of thread count. Each motif's
        //   internal f64 reductions are untouched, because nothing inside a cell
        //   is parallel -- so the numbers are identical for 1 thread and for 16,
        //   not merely close.
        //
        // The error case is the awkward one: the sequential build aborts the whole
        // run on the first failure, so the error is collected here and raised
        // after the fact rather than inside the parallel closure, which cannot
        // return a `Result` without paying for one per motif.
        let sample_id = sample.id();
        let outcomes: Vec<Result<Row>> = motifs
            .par_iter()
            .map(|&i| {
                let motif = &annotation.motifs[i];
                compute_cell(&methylome, &annotation, motif, &bins, options).map(|mut row| {
                    row.sample = sample_id.clone();
                    row
                })
            })
            .collect();
        for (i, outcome) in outcomes.into_iter().enumerate() {
            let motif = &annotation.motifs[motifs[i]];
            match outcome {
                Ok(row) => rows.push(row),
                Err(e) if options.keep_going => {
                    log(format!("keep-going: {} {}: {e}", sample_id, motif.name));
                    rows.push(Row {
                        sample: sample_id.clone(),
                        motif: motif.name.clone(),
                        deviation: None,
                        expected_deviation: None,
                    });
                }
                Err(e) => return Err(e),
            }
        }
        log(format!("Finished processing {}", sample_id));
    }
    Ok(rows)
}

/// Read one sample, interning its chromosome names into the run's table.
pub fn load_sample(
    sample: &Sample,
    options: &Options,
    chroms: &mut ChromTable,
) -> Result<Methylome> {
    match &sample.format {
        InputFormat::Portable => {
            let sites = read_msites(&sample.path, chroms)?;
            // The portable form has already been through a methylome file, so
            // only the threshold applies here -- but the NaN drop is kept so the
            // two entry points agree on what a site is.
            Ok(Methylome {
                sites: sites
                    .into_iter()
                    .filter(|s| !s.score.is_nan() && s.coverage >= options.cov_threshold)
                    .collect(),
            })
        }
        InputFormat::Methylome(ty) => {
            read_methylome(&sample.path, ty, options.cov_threshold, chroms)
        }
    }
}

/// Convenience for the differential tests: an in-memory annotation plus sites.
pub fn run_in_memory(
    sample_id: &str,
    sites: &[Site],
    annotation: &Annotation,
    options: &Options,
) -> Result<Vec<Row>> {
    let motifs = valid_motifs(annotation)?;
    let methylome = Methylome {
        sites: sites
            .iter()
            .copied()
            .filter(|s| !s.score.is_nan() && s.coverage >= options.cov_threshold)
            .collect(),
    };
    let bins = gc_bins(&methylome, annotation, options)?;
    let mut rows = Vec::with_capacity(motifs.len());
    for &i in &motifs {
        match compute_cell(
            &methylome,
            annotation,
            &annotation.motifs[i],
            &bins,
            options,
        ) {
            Ok(mut r) => {
                r.sample = sample_id.to_string();
                rows.push(r);
            }
            Err(e) if options.keep_going => {
                rows.push(Row {
                    sample: sample_id.to_string(),
                    motif: annotation.motifs[i].name.clone(),
                    deviation: None,
                    expected_deviation: None,
                });
                log(format!("keep-going: {sample_id}: {e}"));
            }
            Err(e) => return Err(e),
        }
    }
    Ok(rows)
}

/// GC bin means for one sample, using the (already enhancer-reduced) windows.
///
/// The window index is built once and reused for every motif, which is the point:
/// at genome scale it is millions of entries and building it per motif would cost
/// more than the whole rest of the run.
pub fn gc_bins(
    methylome: &Methylome,
    annotation: &Annotation,
    options: &Options,
) -> Result<GcBins> {
    let index = StartIndex::build(&annotation.gc_windows);
    let sorted = SortedMethylome::new(&methylome.sites);
    bin_means_swept(
        &sorted,
        &annotation.gc_windows,
        &index,
        options.ignore_strand,
    )
}

/// One (sample, motif) cell: the observed profile, the expected profile, and
/// their deviations.
///
/// `Deviation::observed` is `dev_helper(observed profile)` -- upstream's
/// `obs_dev` -- and `Deviation::expected` is `dev_helper(expected profile)`,
/// upstream's `exp_dev`.
fn compute_cell(
    methylome: &Methylome,
    annotation: &Annotation,
    motif: &MotifAnnotation,
    bins: &GcBins,
    options: &Options,
) -> Result<Row> {
    let width = motif_width_of(&motif.tfbs)?;
    let observed = compute_observed(
        &motif.tfbs,
        width,
        annotation.enhancer.as_deref(),
        &methylome.sites,
        options.ignore_strand,
        &motif.name,
    )?;
    let expected = compute_expectations(&motif.gcfreq, bins)?;
    let d = deviations(&observed, &expected);
    Ok(Row {
        sample: String::new(),
        motif: motif.name.clone(),
        deviation: d.deviation(),
        expected_deviation: d.expected,
    })
}

/// Wrap an already-loaded annotation for the differential tests.
pub fn annotation_with(
    gc_windows: Vec<crate::model::GcWindow>,
    motifs: Vec<(String, Vec<crate::model::Range>, crate::model::GcFreq)>,
    ignore_strand: bool,
) -> Annotation {
    Annotation {
        gc_windows,
        motifs: motifs
            .into_iter()
            .map(|(name, tfbs, gcfreq)| MotifAnnotation {
                name,
                tfbs: Arc::new(tfbs),
                gcfreq: Arc::new(gcfreq),
            })
            .collect(),
        enhancer: None,
        ignore_strand,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{GcFreq, Range, Strand};

    fn motif(name: &str, n: usize) -> MotifAnnotation {
        MotifAnnotation {
            name: name.to_string(),
            tfbs: Arc::new(
                (0..n)
                    .map(|i| Range {
                        chr: 0,
                        start: i as i64 * 1000,
                        end: i as i64 * 1000 + 410,
                        strand: Strand::Star,
                    })
                    .collect(),
            ),
            gcfreq: Arc::new(GcFreq::new(1, 2, vec![0.5, 0.5]).unwrap()),
        }
    }

    #[test]
    fn sample_id_is_the_basename() {
        assert_eq!(sample_id(Path::new("/a/b/c.bed")), "c.bed");
        assert_eq!(sample_id(Path::new("c.bed")), "c.bed");
    }

    #[test]
    fn empty_tfbs_motif_is_dropped() {
        let annotation = Annotation {
            gc_windows: Vec::new(),
            motifs: vec![motif("A", 1), motif("EMPTY", 0), motif("B", 2)],
            enhancer: None,
            ignore_strand: true,
        };
        assert_eq!(valid_motifs(&annotation).unwrap(), vec![0, 2]);

        let all_empty = Annotation {
            gc_windows: Vec::new(),
            motifs: vec![motif("EMPTY", 0)],
            enhancer: None,
            ignore_strand: true,
        };
        assert_eq!(
            valid_motifs(&all_empty).unwrap_err().to_string(),
            "No valid motifs remaining after validation."
        );
    }
}
