//! T70: a deterministic genome-scale data generator.
//!
//! ```text
//! cargo run --release --example gen_synthetic -- --out DIR [options]
//! ```
//!
//! Writes a complete annotation directory plus one `bismarkcov` sample, in the
//! portable formats of `AGENT_PLAN.md` section 3, so the R baselines in
//! `scripts/bench_r.R` and the Rust benchmarks in `benches/methyltfr.rs` measure
//! the same work.
//!
//! The scale follows T70 -- about 28 M sites, 7.7 M GC windows and 500 motifs of
//! 250 000 binding sites -- but every number is an option, and the exact command
//! that produced a given benchmark is recorded in `docs/benchmarks.md`.
//!
//! Everything is generated from a single seed with an xorshift generator, so two
//! runs with the same options produce byte-identical files. That is what lets the
//! R baseline and the Rust benchmark be compared on identical input rather than on
//! two statistically similar ones.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use methyltfr::rmath::round6;

/// xorshift64*, so the generator needs no dependency and is reproducible from the
/// seed alone.
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

struct Options {
    out: PathBuf,
    sites: u64,
    windows: u64,
    motifs: usize,
    tfbs_per_motif: u64,
    chromosomes: usize,
    span: u64,
    seed: u64,
    cells: usize,
    gz: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            out: PathBuf::from("/scratch/mdra00001/tmp/opencode/synth"),
            // T70's numbers.
            sites: 28_000_000,
            windows: 7_700_000,
            motifs: 500,
            tfbs_per_motif: 250_000,
            chromosomes: 25,
            span: 250_000_000,
            seed: 0x5EED_1234_5678_9ABC,
            cells: 101,
            gz: false,
        }
    }
}

fn parse_args() -> Result<Options, String> {
    let mut o = Options::default();
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let key = args[i].as_str();
        let value = || -> Result<&str, String> {
            args.get(i + 1)
                .map(|s| s.as_str())
                .ok_or_else(|| format!("{key} needs a value"))
        };
        match key {
            "--out" => {
                o.out = PathBuf::from(value()?);
                i += 2;
            }
            "--sites" => {
                o.sites = value()?.parse().map_err(|e| format!("--sites: {e}"))?;
                i += 2;
            }
            "--windows" => {
                o.windows = value()?.parse().map_err(|e| format!("--windows: {e}"))?;
                i += 2;
            }
            "--motifs" => {
                o.motifs = value()?.parse().map_err(|e| format!("--motifs: {e}"))?;
                i += 2;
            }
            "--tfbs" => {
                o.tfbs_per_motif = value()?.parse().map_err(|e| format!("--tfbs: {e}"))?;
                i += 2;
            }
            "--chromosomes" => {
                o.chromosomes = value()?
                    .parse()
                    .map_err(|e| format!("--chromosomes: {e}"))?;
                i += 2;
            }
            "--span" => {
                o.span = value()?.parse().map_err(|e| format!("--span: {e}"))?;
                i += 2;
            }
            "--seed" => {
                o.seed = value()?.parse().map_err(|e| format!("--seed: {e}"))?;
                i += 2;
            }
            "--cells" => {
                o.cells = value()?.parse().map_err(|e| format!("--cells: {e}"))?;
                i += 2;
            }
            "--gzip" => {
                o.gz = true;
                i += 1;
            }
            other => return Err(format!("unknown option {other}")),
        }
    }
    Ok(o)
}

fn open_out(path: &Path, gz: bool) -> BufWriter<Box<dyn Write>> {
    let file = File::create(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let sink: Box<dyn Write> = if gz {
        Box::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::fast(),
        ))
    } else {
        Box::new(file)
    };
    BufWriter::with_capacity(1 << 20, sink)
}

fn main() {
    let o = match parse_args() {
        Ok(o) => o,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };
    std::fs::create_dir_all(&o.out).expect("output directory");
    let ann = o.out.join("annotation");
    std::fs::create_dir_all(&ann).expect("annotation directory");

    let chr_name = |c: usize| format!("chr{}", c + 1);
    let mut rng = Rng::new(o.seed);

    // ---------------------------------------------------------- GC windows
    // 30 bp, tiled within a chromosome, so windows abut and a width-2 site can
    // straddle two of them -- the "mean over hits" case from section 2.3.
    let gz_path = if o.gz {
        ann.join("gc_windows.tsv.gz")
    } else {
        ann.join("gc_windows.tsv")
    };
    let mut w = open_out(&gz_path, o.gz);
    let _ = writeln!(w, "chr\tstart\tend\tstrand\tgc_bin");
    let per_chr_w = o.windows.div_ceil(o.chromosomes as u64);
    let tile = o.span / per_chr_w.max(1);
    for c in 0..o.chromosomes {
        let name = chr_name(c);
        for k in 0..per_chr_w {
            let start = k * tile + 1;
            // Bins cycle 1..5 so every bin is populated, as a real annotation has.
            let bin = (k % 5) as u32 + 1;
            let _ = writeln!(w, "{name}\t{}\t{}\t*\t{bin}", start, start + tile - 1);
        }
    }
    w.flush().expect("flush windows");
    println!(
        "gc_windows.tsv: {} windows",
        per_chr_w * o.chromosomes as u64
    );

    // ---------------------------------------------------------------- gcfreq
    // One matrix per motif, all rows the same length. `cells` stands in for the
    // motif profile length L; upstream uses 512, and 101 keeps the benchmark
    // comparable without letting the matrix dominate the runtime.
    // One matrix per motif, as upstream's `gcfreqs[[motif]]` is. Five rows (one
    // per GC bin) by `cells` columns; the columns sum to 1 so the expected
    // profile is a distribution.
    let mut rng_g = Rng::new(o.seed ^ 0xA5A5_A5A5);
    for m in 0..o.motifs {
        let mut gf = open_out(&ann.join(format!("gcfreq_{m:04}.tsv")), false);
        for _ in 0..5 {
            let cols: Vec<f64> = (0..o.cells).map(|_| 0.05 + rng_g.unit()).collect();
            let total: f64 = cols.iter().sum();
            let row: String = cols
                .iter()
                .map(|c| format!("{:?}", round6(c / total)))
                .collect::<Vec<_>>()
                .join("\t");
            let _ = writeln!(gf, "{row}");
        }
    }

    // ----------------------------------------------------------------- TFBS
    // 500 separate files, one per motif, as upstream's annotation is one set per
    // motif. Widths are uniform within a file because W comes from the first one.
    let width = 411i64;
    let per_chr_t = o.tfbs_per_motif.div_ceil(o.chromosomes as u64);
    let tspan = o.span / per_chr_t.max(1);
    for m in 0..o.motifs {
        let path = ann.join(format!("tfbs_{m:04}.tsv"));
        let mut f = open_out(&path, false);
        let _ = writeln!(f, "chr\tstart\tend\tstrand");
        for c in 0..o.chromosomes {
            let name = chr_name(c);
            for k in 0..per_chr_t {
                let start = (k * tspan + 1) as i64;
                let strand = if (k as usize + m) % 2 == 0 { '+' } else { '-' };
                let _ = writeln!(f, "{name}\t{start}\t{}\t{strand}", start + width - 1);
            }
        }
        f.flush().expect("flush tfbs");
    }
    println!(
        "tfbs: {} motifs x {} binding sites",
        o.motifs,
        per_chr_t * o.chromosomes as u64
    );

    // ------------------------------------------------------------- manifest
    let mut mf = open_out(&ann.join("motifs.tsv"), false);
    let _ = writeln!(mf, "motif\ttfbs_path\tgcfreq_path");
    for m in 0..o.motifs {
        let _ = writeln!(mf, "M{m}\ttfbs_{m:04}.tsv\tgcfreq_{m:04}.tsv");
    }
    drop(mf);

    // --------------------------------------------------------------- sample
    // One bismarkcov sample: `chr start end percent meth unmeth`, BED-style so
    // every site has width 2.
    let sample = o.out.join("sample0.bismarkCov");
    let mut f = open_out(&sample, false);
    let per_chr_s = o.sites.div_ceil(o.chromosomes as u64);
    let sspan = o.span / per_chr_s.max(1);
    for c in 0..o.chromosomes {
        let name = chr_name(c);
        for k in 0..per_chr_s {
            let start = k * sspan + 1;
            let meth = 1 + rng.below(20);
            let unmeth = rng.below(20);
            let pct = round6(100.0 * meth as f64 / (meth + unmeth) as f64);
            let _ = writeln!(f, "{name}\t{start}\t{}\t{pct}\t{meth}\t{unmeth}", start + 1);
        }
    }
    f.flush().expect("flush sample");
    println!(
        "sample0.bismarkCov: {} sites",
        per_chr_s * o.chromosomes as u64
    );

    // Sanity: the files must be readable by the same reader the pipeline uses.
    {
        use methyltfr::io::portable::{read_gc_windows, read_gcfreq};
        use methyltfr::model::ChromTable;
        let mut chroms = ChromTable::new();
        let windows = read_gc_windows(&gz_path, &mut chroms).expect("re-read windows");
        let freq = read_gcfreq(&ann.join("gcfreq_0000.tsv")).expect("re-read gcfreq");
        assert_eq!(windows.len() as u64, per_chr_w * o.chromosomes as u64);
        assert_eq!(freq.rows, 5);
        assert_eq!(freq.cols, o.cells);
        use methyltfr::io::methylome::read_methylome;
        let m = read_methylome(&sample, "bismarkcov", 1.0, &mut chroms).expect("re-read sample");
        assert_eq!(m.sites.len() as u64, per_chr_s * o.chromosomes as u64);
    }
    println!("wrote {}", o.out.display());
}
