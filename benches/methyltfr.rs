//! T70: the four benchmarks `docs/AGENT_PLAN.md` M7 asks for.
//!
//! The data is built once, in memory, before the measured closure, so criterion
//! measures the computation rather than the disk. Scale is chosen so one iteration
//! is tens of milliseconds: a benchmark that takes a second per sample measures
//! the scheduler, and one that takes ten minutes is not something anyone reruns.
//! `docs/benchmarks.md` records the end-to-end numbers at genome scale, which is
//! where the reported speed-up comes from.
//!
//! The scale here is 400 000 sites, 100 000 GC windows and 4 motifs of 20 000
//! binding sites across 4 chromosomes.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use methyltfr::deviation::{compute_observed, compute_observed_swept, deviations, motif_width_of};
use methyltfr::expected::compute_expectations;
use methyltfr::gc::{bin_means, bin_means_swept};
use methyltfr::intervals::{Query, StartIndex, site_windows};
use methyltfr::model::{GcFreq, GcWindow, Range, Site, SortedMethylome, Strand};

const CHROMOSOMES: u32 = 4;
const SPAN: i64 = 60_000_000;
const SITES: usize = 400_000;
const WINDOWS: usize = 100_000;
const MOTIFS: usize = 4;
const TFBS_PER_MOTIF: usize = 20_000;
const WINDOW_WIDTH: i64 = 30;
const TFBS_WIDTH: i64 = 411;

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n.max(1)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A methylome: width-2 sites, strands on all three values, sorted by chromosome
/// then start, which is what a real methylation file looks like.
fn build_sites() -> Vec<Site> {
    let mut rng = Rng::new(0xBEEF_1234);
    let mut out = Vec::with_capacity(SITES);
    for _ in 0..SITES {
        let chr = rng.below(CHROMOSOMES as u64) as u32;
        let start = (rng.below(SPAN as u64)) as i64 + 1;
        let strand = match rng.below(3) {
            0 => Strand::Plus,
            1 => Strand::Minus,
            _ => Strand::Star,
        };
        out.push(Site {
            chr,
            start,
            end: start + 1,
            strand,
            score: (rng.unit() * 1000.0).round() / 1000.0,
            coverage: 1.0 + rng.below(20) as f64,
        });
    }
    out.sort_by_key(|s| (s.chr, s.start, s.end));
    out
}

/// GC windows: 30 bp, tiled, bins cycling 1..5.
fn build_windows() -> Vec<GcWindow> {
    let per_chr = WINDOWS / CHROMOSOMES as usize;
    let tile = SPAN / per_chr as i64;
    let mut out = Vec::with_capacity(WINDOWS);
    for c in 0..CHROMOSOMES {
        for k in 0..per_chr {
            let start = k as i64 * tile + 1;
            out.push(GcWindow {
                chr: c,
                start,
                end: start + WINDOW_WIDTH - 1,
                strand: Strand::Star,
                gc_bin: (k % 5) as u8 + 1,
            });
        }
    }
    out
}

/// Binding sites for one motif: uniform width, sorted by start.
fn build_tfbs(motif: usize) -> Vec<Range> {
    let mut rng = Rng::new(0xC0FFEE + motif as u64);
    let mut out = Vec::with_capacity(TFBS_PER_MOTIF);
    for _ in 0..TFBS_PER_MOTIF {
        let chr = rng.below(CHROMOSOMES as u64) as u32;
        let start = (rng.below(SPAN as u64)) as i64 + 1;
        let strand = match rng.below(2) {
            0 => Strand::Plus,
            _ => Strand::Minus,
        };
        out.push(Range {
            chr,
            start,
            end: start + TFBS_WIDTH - 1,
            strand,
        });
    }
    out.sort_by_key(|r| (r.chr, r.start, r.end));
    out
}

fn build_gcfreq() -> GcFreq {
    let mut rng = Rng::new(0x5A5A);
    let cols = 512;
    let mut data = vec![0.0; 5 * cols];
    for j in 0..cols {
        let mut total = 0.0;
        for i in 0..5 {
            let v = 0.05 + rng.unit();
            data[i * cols + j] = v;
            total += v;
        }
        for i in 0..5 {
            data[i * cols + j] /= total;
        }
    }
    GcFreq::new(5, cols, data).unwrap()
}

/// The interval queries on their own, with the hits summed so the optimiser cannot
/// delete the work.
fn bench_interval_overlap(c: &mut Criterion) {
    let sites = build_sites();
    let windows = build_windows();
    let tfbs = build_tfbs(0);
    let width = motif_width_of(&tfbs).unwrap();
    let index = StartIndex::build(&windows);
    let tfbs_index = StartIndex::build(&tfbs);

    c.bench_function("interval_overlap", |b| {
        b.iter(|| {
            let mut hits = 0usize;
            for s in &sites {
                hits += site_windows(&index, &windows, s, true).count();
            }
            for s in sites.iter().take(SITES / 4) {
                hits += tfbs_index
                    .within_uniform(&tfbs, Query::from_site(s, true), width)
                    .count();
            }
            hits
        })
    });
}

/// GC bin assignment, the section 2.3 hot loop, both implementations.
///
/// The `_swept` variant is what the pipeline uses; the other is kept so the gain
/// from M7 T72 stays measurable instead of being a claim.
fn bench_gc_assignment(c: &mut Criterion) {
    let sites = build_sites();
    let windows = build_windows();
    let index = StartIndex::build(&windows);
    let sorted = SortedMethylome::new(&sites);

    c.bench_function("gc_assignment", |b| {
        b.iter(|| bin_means(&sites, &windows, &index, true).unwrap())
    });
    c.bench_function("gc_assignment_swept", |b| {
        b.iter(|| bin_means_swept(&sorted, &windows, &index, true).unwrap())
    });
}

/// One motif end to end: resize, hits, both profiles, both deviations.
fn bench_motif_deviation(c: &mut Criterion) {
    let sites = build_sites();
    let windows = build_windows();
    let index = StartIndex::build(&windows);
    let bins = bin_means(&sites, &windows, &index, true).unwrap();
    let gcfreq = build_gcfreq();
    let mut rng = Rng::new(77);
    let motifs: Vec<Vec<Range>> = (0..MOTIFS).map(build_tfbs).collect();

    let sorted = SortedMethylome::new(&sites);

    c.bench_function("motif_deviation", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for tfbs in &motifs {
                let width = motif_width_of(tfbs).unwrap();
                let observed = compute_observed(tfbs, width, None, &sites, true, "M").unwrap();
                let expected = compute_expectations(&gcfreq, &bins).unwrap();
                let d = deviations(&observed, &expected);
                acc += d.deviation().unwrap_or(f64::NAN);
                acc += rng.unit().min(1.0);
            }
            acc
        })
    });

    c.bench_function("motif_deviation_swept", |b| {
        b.iter(|| {
            let mut acc = 0.0f64;
            for tfbs in &motifs {
                let width = motif_width_of(tfbs).unwrap();
                let observed =
                    compute_observed_swept(tfbs, width, None, &sorted, true, "M").unwrap();
                let expected = compute_expectations(&gcfreq, &bins).unwrap();
                let d = deviations(&observed, &expected);
                acc += d.deviation().unwrap_or(f64::NAN);
                acc += rng.unit().min(1.0);
            }
            acc
        })
    });
}

/// Reading a methylome plus the whole run, i.e. what the CLI does per sample.
fn bench_end_to_end(c: &mut Criterion) {
    let sites = build_sites();
    let windows = build_windows();
    let index = StartIndex::build(&windows);
    let bins = bin_means(&sites, &windows, &index, true).unwrap();
    let gcfreq = build_gcfreq();
    let motifs: Vec<Vec<Range>> = (0..MOTIFS).map(build_tfbs).collect();
    let motif_width = motif_width_of(&motifs[0]).unwrap();
    let sorted = SortedMethylome::new(&sites);

    c.bench_function("end_to_end", |b| {
        b.iter(|| {
            let mut rows = Vec::with_capacity(motifs.len());
            for (i, tfbs) in motifs.iter().enumerate() {
                let observed =
                    compute_observed_swept(tfbs, motif_width, None, &sorted, true, "M").unwrap();
                let expected = compute_expectations(&gcfreq, &bins).unwrap();
                let d = deviations(&observed, &expected);
                rows.push((i, d.deviation(), d.expected));
            }
            rows
        })
    });
}

criterion_group!(
    benches,
    bench_interval_overlap,
    bench_gc_assignment,
    bench_motif_deviation,
    bench_end_to_end
);
criterion_main!(benches);

/// Criterion needs a warm-up allowance on the first run of a criterion target;
/// this keeps the default 3 s measurement period from being spent entirely on a
/// cold page cache.
#[allow(dead_code)]
fn warmup_hint() -> Duration {
    Duration::from_secs(3)
}
