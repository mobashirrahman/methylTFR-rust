# methylTFR-rs — agent implementation plan (v0.1)

This file is the working spec for implementing agents. It replaces the need to
re-derive behaviour from the R source: every rule below was read from upstream
at the pinned commit, and the ones marked **[verified]** were also checked by
running R 4.3.3. Rules marked **[confirm-A]** are from source reading only and
must be confirmed by the Tier A oracle (task T04) before v0.1 is declared done.

- Upstream: `https://github.com/EpigenomeInformatics/methylTFR`
- Pinned SHA: `8aeab03cb469a9c213b7cd7fd2cbe6eeb8db5856` (2026-09-29, package version 0.99.9)
- Upstream licence: MIT, `YEAR: 2025`, `COPYRIGHT HOLDER: Epigenome Informatics`

## 0. Rules for every agent

1. Do exactly one task per session. Do not start the next task.
2. Before coding, read this file's section 2 (semantics) for your module. Do not
   "improve" any rule in it, even if it looks like a bug. Odd behaviour is listed
   on purpose.
3. A task is done only when its **Done when** command passes. Paste the command
   output in your final message. If it fails, say so; do not weaken the test.
4. Never change a tolerance, a fixture file, or an expected value to make a test
   pass. If a fixture looks wrong, stop and report.
5. Never use `f64::round()` for anything that mirrors R. Use `rmath` (T10).
6. All coordinates are signed `i64`, 1-based, closed `[start, end]`. Never add or
   subtract 1 to "convert BED". The parsers do not convert, so neither do you.
7. No `rayon`, no `unsafe`, no new dependency beyond section 3 without asking.
8. After each task: `cargo fmt --check && cargo clippy --all-targets -- -D warnings && cargo test`.
9. If the R source and this file disagree, stop and report the line numbers.
   Do not pick one silently.

## 1. What changed from the original plan

| # | Original plan | Correction |
|---|---|---|
| 1 | Golden values `1.743674 / 0.9835985`, bins `0.5315789, 0.6466688, …` | These are from the pkgdown vignette rendered 2026-08-28. The fixtures were changed afterwards (commit `db2f00f`). They do **not** match the pinned SHA. See section 4. |
| 2 | Enhancer filtering applies to TFBS | It also subsets the GC windows (`methyltfr_core`), so GC-bin means change too. |
| 3 | "GC-bin mean = mean(score)" | It is the mean over *overlap hits*, not sites. Windows are abutting 30 bp tiles and BED-style sites have width 2, so a site straddling two windows counts twice (44 of 1000 sites in the fixture). |
| 4 | "Preserve six-decimal rounding" | R's `round(x, 6)` is not `(x*1e6).round()/1e6`. Naive rounding differs on real inputs, e.g. `1/128 → 0.007812` in R, `0.007813` naive. |
| 5 | "Test midpoint rounding" | R `round()` is half-to-even. Midpoints and expected positions both depend on it. |
| 6 | "Reproduce position generation" | For even profile length the positions skip `0` (L=512 gives −256…−1, 1…256). |
| 7 | "Evaluate noodles / Rust-Bio" | Decided: use neither. None of the six formats is a noodles format, and the two overlap queries reduce to binary searches (section 2.4). |
| 8 | Differential tests run R | R on the dev machine has no Bioconductor. Split into Tier A (real package, offline, committed outputs) and Tier B (base-R oracle, always runnable). CI is Rust-only. |
| 9 | No `rmath` module | Added. It is the highest-risk unit and is built and tested first. |
| 10 | Missing bins/intervals "behave like upstream" | Upstream aborts the whole run (`stop()`). Made explicit in section 2.8. |

## 2. Reference semantics

### 2.1 R numeric shim (`src/rmath.rs`)

`round_half_even(x: f64) -> f64` — R `round(x)`: nearest integer, ties to even.
`round(2.5)=2`, `round(0.5)=0`, `round(-0.5)=-0`, `round(3.5)=4`. Use
`f64::round_ties_even`. **[verified]**

`round_digits(x: f64, d: i32) -> f64` — R ≥ 4.0 `round(x, d)`, for `d > 0`,
all arithmetic in `f64`: **[verified: 0 mismatches on 7.5 M values, d = 6]**

```
if x is NaN or infinite or 0: return x
if x < 0: return -round_digits(-x, d)
p   = 10f64.powi(d)
x10 = p * x
i10 = x10.floor()
xd  = i10 / p
xu  = x10.ceil() / p
dd  = x - xd
du  = xu - x
if dd < du || (dd == du && i10 % 2.0 == 0.0) { xd } else { xu }
```

`seq_len_out(from: f64, to: f64, n: usize) -> Vec<f64>` — R `seq(from, to, length.out = n)`:

```
n == 0 -> []          n == 1 -> [from]          n == 2 -> [from, to]
n > 2 and from == to -> [from; n]
n > 2 -> by = (to - from) / (n - 1) as f64
         [from, from + 1.0*by, from + 2.0*by, …, from + (n-2)*by, to]   // last is exactly `to`
```

`cut_index(x: f64) -> Option<usize>` — R `cut(x, c(-250,-200,-25,25,200,250))`.
Intervals are open on the left, closed on the right: **[verified]**

| index | interval |
|---|---|
| 0 | (−250, −200] |
| 1 | (−200, −25] |
| 2 | (−25, 25] |
| 3 | (25, 200] |
| 4 | (200, 250] |

`x = -250` → `None`. `x = 250` → `Some(4)`. `x < -250` or `x > 250` or NaN → `None`.

`floor_div2(a: i64) -> i64` — R `a %/% 2L`, floor division. `-131 → -66`. Use `a.div_euclid(2)`.

### 2.2 Methylome parsing (`read_methylome`)

All files are tab-separated, optionally gzip (detect by magic bytes `1f 8b`).
Columns are 1-indexed below. Rounding is `round_digits(·, 6)`.

| type (case-insensitive) | header | chr | start | end | strand | score | coverage | min cols |
|---|---|---|---|---|---|---|---|---|
| `epp` | none | c1 | c2 | c3 | c6 | `c5 / 1000` | number after `/` in c4, with any `'` removed | 6 |
| `bissnp` | line 1 always skipped | c1 | c2 | c3 | c6 | `c4 / 100` | c5 | 6 |
| `allc` | none | c1 | c2 | c2 | c3 | `c5 / c6` | c6 | 6 |
| `bismarkcytosine` | none | c1 | c2 | c2 | c3 | `c4 / (c4 + c5)` | `c4 + c5` | 5 |
| `bismarkcov` | none | c1 | c2 | c3 | always `*` | `c4 / 100` | `c5 + c6` | 6 |
| `encode` | auto-detected | c1 | c2 | c3 | c6 | `c11 / 100` | c10 | 11 |

Then, in this order:

1. Drop records whose score is NaN (happens when coverage is 0 in `allc` / `bismarkcytosine`).
2. Keep records with `coverage >= cov_threshold` (default 1; `f64` comparison).
3. If any record remains and `min(score) < 0` or `max(score) > 1`: error.

Further rules:

- `start`/`end` are copied verbatim into 1-based closed ranges. BED-style inputs
  therefore have width 2, `allc` / `bismarkcytosine` width 1. This is upstream behaviour.
- `allc` is **not** filtered by context column c4. All contexts are kept.
- `encode`: error if any c11 is `< 0` or `> 100`. Header rule for Rust: the first
  line is a header iff c2 or c3 of that line does not parse as an integer. This
  approximates `fread(header = "auto")` and is recorded as a divergence.
- Strand must be `+`, `-` or `*`. Anything else (including `.`) is an error **[confirm-A]**.
- A non-numeric value in a numeric column is an error with file name and line number.
- Fewer than "min cols" columns: error.
- Parse numbers with `str::parse::<f64>` / `::<i64>`. `fread`'s float parser can
  differ by 1 ulp before rounding; recorded as a divergence, not chased.

### 2.3 GC-bin means (`addGCBintoMethylome`)

Inputs: sites, GC windows (`chr, start, end, strand, gc_bin` with `gc_bin` an integer 1..5).

- A site hits a window when same `chr` (exact string match, no `chr` prefix
  normalisation), `window.start <= site.end && window.end >= site.start`, and the
  strands are compatible (2.5). This is `findOverlaps(type = "any")`.
- One site may hit several windows; each hit contributes separately.
- Zero hits in total: error `No methylation sites found in the GC distribution`.
- `mean[bin] = sum(score over hits in bin) / n_hits[bin]`, plain sequential `f64`
  sum in hit order (sites in input order, then windows by ascending start).
- Result is ordered by bin and contains only populated bins.

### 2.4 Per-motif observed profile (`computeDeviation`)

Given the motif's TFBS list (`chr, start, end, strand`) in file order:

1. `W = (end[0] - start[0] + 1) + 130`, using the **first** TFBS before any filtering.
2. Resize every TFBS to width `W` about its centre. Strand does not affect this **[confirm-A]**:
   `new_start = start + floor_div2((end - start + 1) - W)`, `new_end = new_start + W - 1`.
   `new_start` may be ≤ 0. Keep it.
3. If enhancer regions are given: keep resized TFBS that overlap (type "any",
   strand rule 2.5) at least one enhancer region.
4. A site is *within* a TFBS when same `chr`, `tfbs.new_start <= site.start &&
   tfbs.new_end >= site.end`, strands compatible. Every (site, TFBS) pair is a hit;
   a site inside k TFBS yields k hits.
5. Zero hits: error `No methylation sites found in the <motif> binding sites`.
6. Midpoint, computed literally in `f64` then converted to `i64`:
   `mid = round_half_even(new_end as f64 + ((new_start - new_end) as f64) / 2.0)`.
   When `W` is even this lands on `.5` and rounds to the even neighbour — do not
   replace it with integer arithmetic.
7. For each hit: `x = site.start - mid` (never `site.end`; never strand-flipped),
   `value = site.score`.

Implementation note (allowed because all resized TFBS share width `W`): per
chromosome, sort TFBS by `new_start`; the TFBS containing site `[s, e]` are
exactly those with `e - W + 1 <= new_start <= s`, i.e. one contiguous slice found
by two binary searches. For 2.3 use the windows sorted by start plus the maximum
window width `M`: candidates have `s - M + 1 <= window.start <= e`, then filter
`window.end >= s`.

### 2.5 Strand compatibility

- `ignore_strand = true` (default): always compatible.
- `ignore_strand = false`: compatible iff `a == b || a == '*' || b == '*'`.

Applies to every overlap in 2.3, 2.4 and 2.6.

### 2.6 Enhancer option

When enhancer regions are supplied, **both** of these happen:

- GC windows are reduced to those overlapping any enhancer region (type "any",
  strand rule 2.5), once, before any sample is processed. Bin means are then
  computed from the reduced windows.
- Resized TFBS are filtered as in 2.4 step 3.

### 2.7 Expected profile (`computeExpectations`)

`gcfreq` is a matrix with one row per GC bin (5) and `L` columns (512 for BATF).

- If `number of populated bins != gcfreq row count`: error (R fails with
  "non-conformable arguments"). Rust also errors if the populated bins are not
  exactly `1..=rows`.
- `expected[j] = Σ_i gcfreq[i][j] * mean[i]`, summed in ascending `i`.
- `h = floor(L / 2)`; `x[j] = round_half_even(seq_len_out(-h, h, L)[j])`.
  Odd `L` gives −h…h. Even `L` skips 0. **[verified]**

### 2.8 `dev_helper` and final values

For a profile `(x[], value[])`:

1. `idx = cut_index(x)`; skip entries with `None`.
2. Per interval: `n`, and `mean = sum / n` (sequential `f64`, input order).
3. Keep intervals with `n > 0` and non-NaN mean, in ascending index. Let `k` be the count.
4. `k == 0` → NA.
5. `D = mean[(k + 1) / 2 - 1] / ((mean[0] + mean[k - 1]) / 2.0)` (integer division).
   So k=5→3rd, k=4→2nd, k=3→2nd, k=2→1st, k=1→1st (D = 1). Do not special-case.
   Division by zero yields ±inf or NaN and is passed through.

`expected_deviation = D(expected profile)`, `deviation = D(observed) − D(expected)`.
NA in either makes `deviation` NA.

Run-level behaviour:

- Motif order = order of the annotation manifest. Motifs with no TFBS or no
  matrix are dropped with a log line. No motif left: error.
- Sample id = file basename.
- Upstream aborts the whole run on the errors in 2.3 / 2.4 step 5 / 2.7. The CLI
  default does the same (non-zero exit). `--keep-going` is a Rust extension that
  writes `NA` for that cell and logs the reason; it is off by default.

### 2.9 Out of scope for v0.1

Plotting, Bioconductor classes, AnnotationHub, RnBeads, differential tests,
variability, the `z` assay (row z-scores), HDF5 sinks, chunking, RDS/RDA parsing
in Rust, R bindings, parallelism.

## 3. Repository layout and dependencies

```
Cargo.toml  LICENSE  NOTICE  README.md  .github/workflows/ci.yml
docs/      AGENT_PLAN.md  reference-version.md  divergences.md  benchmarks.md
scripts/   export_reference_data.R  reference_base.R  run_reference.R
           gen_rmath_tables.R  gen_differential_cases.R
src/       lib.rs  error.rs  model.rs  rmath.rs  intervals.rs  gc.rs
           expected.rs  deviation.rs  pipeline.rs  cli.rs  main.rs
           io/{mod,portable,epp,bissnp,allc,bismark_cytosine,bismark_cov,encode}.rs
tests/     rmath.rs  interval_semantics.rs  gc_bins.rs  expected_profile.rs
           deviation_helper.rs  parser_{epp,bissnp,allc,bismark_cytosine,bismark_cov,encode}.rs
           golden_batf.rs  differential_vs_r.rs  cli.rs
tests/fixtures/  rmath/  parsers/  batf/  batf_1d99721/  differential/
benches/   (M7 only)
reference/ (gitignored clone of upstream at the pinned SHA)
```

Dependencies: `clap` (derive), `flate2`, `thiserror`; `anyhow` in `main.rs` only.
Dev: `tempfile`, `assert_cmd`. `criterion` and `rayon` only in M7. No `csv`,
`serde`, `noodles`, `bio`: fields are plain tab-separated and never quoted.

### Portable file formats (all TSV, 1-based closed coordinates, header line unless stated)

| file | columns |
|---|---|
| `msites.tsv` | `chr start end strand score coverage` |
| `tfbs.tsv[.gz]` | `chr start end strand` — original, **not** resized, original order |
| `gc_windows.tsv[.gz]` | `chr start end strand gc_bin` |
| `gcfreq.tsv` | no header; one row per bin, `L` columns |
| `enhancer.tsv` | `chr start end strand` |
| `motifs.tsv` (manifest) | `motif tfbs_path gcfreq_path`, paths relative to the manifest; row order = motif order |

Floats are written by R with `sprintf("%.17g", x)` so they round-trip exactly.

Tolerance for every comparison against R: `|rust − r| <= 1e-10` **and**
`|rust − r| <= 1e-10 * |r|`, except parser scores and `rmath` tables, which must
be bit-identical. Tests print max abs and max rel error.

## 4. Golden values

The numbers in the original plan come from fixture revision `1d99721` and are the
only values upstream has published. Both revisions are used as golden cases.

| quantity | `batf_1d99721` (published vignette) | `batf` (pinned SHA, provisional) |
|---|---|---|
| bin 1 mean | 0.5315789 | 0.8125000000 |
| bin 2 mean | 0.6466688 | 0.5871203390 |
| bin 3 mean | 0.7217031 | 0.7217031250 |
| bin 4 mean | 0.7091566 | 0.7091565789 |
| bin 5 mean | 0.7838198 | 0.7838197531 |
| GC hits per bin | — | 8, 59, 96, 152, 729 (1044 total) |
| observed hits | — | 29, intervals n = 3, 9, 1, 9, 5 |
| observed D | — | 2.727272727 |
| `exp_dev` | 0.9835985 | 0.9798459268 |
| `dev` | 1.743674 | 1.7474268 |

"Provisional" means: computed by an independent base-R reimplementation of
section 2 (not by the methylTFR package, which is not installed on the dev
machine). That reimplementation reproduces the published column to all 7 printed
digits, which is the evidence that section 2 is right. T04 replaces the
provisional column with full-precision Tier A output. If Tier A disagrees with
the provisional column beyond 1e-9, stop and report — section 2 is then wrong.

Fixture facts (pinned SHA): 1000 sites, all `chr1`, width 2, strands `+`/`-`,
sorted; 268 717 TFBS, all width 411 (so `W = 541`, odd), both strands; GC windows
567 × 30 bp, strand `*`, 234 abutting pairs; `gcfreq` 5 × 512, columns sum to 1.
`tf_bindsites` is a plain R `list` of `GRanges`, not a `GRangesList`.

## 5. Tasks

Difficulty: **S** = small model fine; **M** = needs care, review the diff; **H** =
give to the strongest available model or a human.

### M0 — bootstrap

**T00 Repo skeleton (S).** `git init`; `cargo init --lib --name methyltfr`; add
`src/main.rs` stub, MIT `LICENSE` for this project, `NOTICE` naming upstream
(authors from `DESCRIPTION`, MIT, "2025 Epigenome Informatics") and stating that
`tests/fixtures/batf*` are derived from upstream `inst/extdata`; `.gitignore`
with `target/`, `reference/`, `methylTFR_tmp/`; CI workflow running fmt, clippy
`-D warnings`, test on `ubuntu-latest` stable.
Done when: the three cargo commands pass on an empty lib.

**T01 Pin reference (S).** Clone upstream into `reference/methylTFR`, checkout the
pinned SHA, write `docs/reference-version.md` (SHA, commit date, package version,
licence, list of the R files section 2 was derived from, `R --version`). Create
`docs/divergences.md` with the three divergences already named in section 2.2
and 2.8 (`encode` header rule, float parsing, `--keep-going`).
Done when: `git -C reference/methylTFR rev-parse HEAD` prints the pinned SHA.

**T02 R environment for Tier A (H, needs human approval before installing).**
Install into a user library (`R_LIBS_USER`), or a container: Bioconductor 3.18
for R 4.3 with `GenomicRanges`, `IRanges`, `S4Vectors`, `SummarizedExperiment`,
`BiocParallel`, `DelayedArray`, `HDF5Array`, plus CRAN `data.table`, `R.utils`,
`logger`, `matrixStats`, `stringr`, `ggplot2`; then `R CMD INSTALL reference/methylTFR`.
Append `sessionInfo()` to `docs/reference-version.md`.
Done when: `Rscript -e 'library(methylTFR); packageVersion("methylTFR")'` prints `0.99.9`.
This task does not block T03, T05 or M1–M4.

**T03 Fixture exporter (M).** `scripts/export_reference_data.R <extdata_dir> <out_dir>`.
Must run under **base R only**: `load()` the `.rda` and read S4 slots with `attr()`:

```r
rle <- function(r) { v <- attr(r, "values"); if (is.factor(v)) v <- as.character(v); rep(v, attr(r, "lengths")) }
gr  <- function(o) {
  rg <- attr(o, "ranges")
  d  <- data.frame(chr = rle(attr(o, "seqnames")), start = attr(rg, "start"),
                   end = attr(rg, "start") + attr(rg, "width") - 1L,
                   strand = rle(attr(o, "strand")), stringsAsFactors = FALSE)
  md <- attr(attr(o, "elementMetadata"), "listData")
  for (n in names(md)) d[[n]] <- md[[n]]
  d
}
```

Write `msites.tsv`, `tfbs.tsv.gz`, `gc_windows.tsv`, `gcfreq.tsv`, `motifs.tsv`
(one row, `BATF`) in the formats of section 3. Run it twice: on the pinned
checkout → `tests/fixtures/batf/`, and on the four files extracted with
`git show 1d99721:inst/extdata/<file>` → `tests/fixtures/batf_1d99721/`. Also copy
the six `*.tsv.gz` parser examples to `tests/fixtures/parsers/`. Write
`SHA256SUMS` in each fixture directory.
Done when: row counts match section 4 (1000 / 268717 / 567 / 5×512) and
`sha256sum -c SHA256SUMS` passes.

**T04 Reference outputs (M).** Two scripts with identical output files:
`expected_bins.tsv` (`gc_bin mean n_hits`), `observed_profile.tsv` (`x value`),
`expected_profile.tsv` (`x value`), `expected_dev.tsv` (`motif obs_d exp_d dev`).

- `scripts/reference_base.R` (Tier B): base-R implementation of section 2,
  reading the portable files. Loops are fine; clarity over speed.
- `scripts/run_reference.R` (Tier A, needs T02): calls `addGCBintoMethylome`,
  `computeDeviation` and the internals `methylTFR:::computeExpectations`,
  `methylTFR:::dev_helper` on GRanges rebuilt from the same portable files.

Commit Tier B outputs as `expected_tierB/`, Tier A as `expected/`.
Done when: Tier B on `batf_1d99721` prints `dev` = 1.743674 and `exp_dev` =
0.9835985 at 7 significant digits; and, once T02 exists, Tier A and Tier B agree
within 1e-10 on both fixture sets.

**T05 rmath tables (S).** `scripts/gen_rmath_tables.R` (base R) writes to
`tests/fixtures/rmath/`: `round6.tsv` (`x`, `round(x, 6)`; at least 200 000 rows:
all `mc/cov` for `cov ≤ 300`, `k/2^n` for `n ≤ 20`, random percentages `/100`,
`0:1000/1000`), `round0.tsv` (ties `±k.5`, random), `seq.tsv` (for `L` in 1..1200:
`L`, comma-joined `round(seq(-floor(L/2), floor(L/2), length.out = L))`),
`cut.tsv` (integers −260..260 plus `±0.5` offsets → interval index or `NA`),
`midpoint.tsv` (random `start`, `end` with odd and even widths →
`round(end + (start - end)/2)`). All floats `%.17g`.
Done when: files exist and `round6.tsv` contains the row `0.0078125 → 0.007812`.

### M1 — primitives

**T10 `rmath.rs` (M).** Implement section 2.1 exactly. `tests/rmath.rs` reads every
T05 table and requires bit equality (`to_bits`).
Done when: `cargo test --test rmath` passes with 0 mismatches.

**T11 `model.rs`, `error.rs`, `io/portable.rs` (S).** Types: `Strand {Plus, Minus, Star}`,
`Site {chr, start, end, strand, score, coverage}`, `Range {chr, start, end, strand}`,
`GcWindow`, `GcFreq {rows, cols, data}`, `Methylome`. Chromosomes interned to `u32`
ids through one shared table. Readers for every portable format in section 3
(gzip by magic bytes). Errors via `thiserror`, always with path and line number.
Done when: unit tests load `tests/fixtures/batf/*` and assert the section 4 row counts.

**T12 `intervals.rs` (M).** `StartIndex` (per-chromosome ranges sorted by start, with
original indices and max width) and three functions: `any_overlaps(site, index,
ignore_strand) -> iter of window idx`, `within_uniform(site, index, W,
ignore_strand) -> iter of tfbs idx`, `overlaps_any(range, index, ignore_strand) -> bool`.
Also `resize_center(ranges) -> (Vec<Range>, W)` and `midpoint(start, end) -> i64`.
`tests/interval_semantics.rs` must cover: touching ends (hit), off-by-one misses
on both sides, site equal to TFBS, site sticking out by 1 on each side, width-2
site straddling two abutting windows (2 hits), every strand pair under both strand
modes, different chromosomes, `new_start <= 0`, mixed-width input to
`resize_center` with odd and even differences, `midpoint` against `midpoint.tsv`,
and a brute-force O(n·m) comparison on 2 000 random ranges.
Done when: `cargo test --test interval_semantics` passes.

### M2 — core

**T20 `gc.rs` (S).** Section 2.3. Test against `batf/expected_tierB/expected_bins.tsv`
(means and `n_hits`), plus: no hits → error; a bin with no hits is absent.
Done when: `cargo test --test gc_bins` passes and prints max abs/rel error.

**T21 `expected.rs` (S).** Section 2.7. Test against `expected_profile.tsv` for both
fixture sets; unit tests for `L` = 1, 2, 4, 5; missing bin → error.
Done when: `cargo test --test expected_profile` passes.

**T22 `deviation.rs` (M).** Section 2.8 `dev_helper`, and the observed profile of 2.4
steps 4–7. `tests/deviation_helper.rs`: hand-built profiles for k = 0..5,
boundary values −250, −200, −25, 25, 200, 250, zero denominator, and the fixture
`observed_profile.tsv` (compare as multisets of `(x, value)`).
Done when: `cargo test --test deviation_helper` passes.

**T23 `pipeline.rs` + golden test (M).** `run(samples, annotation, options) -> Vec<Row>`
where `Row {sample, motif, deviation: Option<f64>, expected_deviation: Option<f64>}`,
implementing 2.6 and the run-level rules of 2.8, sequentially.
`tests/golden_batf.rs` asserts `obs_d`, `exp_d`, `dev` for `batf` and
`batf_1d99721` against `expected_dev.tsv` at 1e-10, and separately asserts the
`batf_1d99721` values equal the published `1.743674` / `0.9835985` at 5e-7.
Done when: `cargo test --test golden_batf` passes. **Nothing in M3+ starts before this.**

### M3 — parsers (one task each, S; do them in this order)

**T30 `bismark_cov`, T31 `epp`, T32 `bissnp`, T33 `allc`, T34 `bismark_cytosine`,
T35 `encode`.** Each: implement the row of the section 2.2 table and one
`tests/parser_<type>.rs` that compares against `tests/fixtures/parsers/<type>.expected.tsv`
(positions and strand equal, score and coverage bit-identical). The expected
files come from `read_methylome()` under Tier A; until T02 exists, generate them
with a base-R transcription of the table (`read.delim`, `round(…, 6)`) and mark
them `tierB` in the filename. Each test also covers: too few columns, a
non-numeric field, a bad strand, gzip and plain input.
Known expectations from upstream's own tests: `allc` example → 3 records;
`encode` example → 6 records, scores `0.06 0.03 0 1 0.55 1`, coverage
`62 62 31 5 31 10`, and 4 records at `cov_threshold = 20`; `bismarkCytosine`
example has one `0/0` row that must be dropped.
**T36 `io/mod.rs` (S).** `read_methylome(path, type, cov_threshold)`: dispatch on
lowercase type, then the three post-steps of 2.2 in order. Tests: threshold 0, 1,
20; score out of range → error; unknown type → error.
Done when (each): `cargo test --test parser_<type>` passes.

### M4 — CLI

**T40 `cli.rs` / `main.rs` (S).**

```
methyltfr run --format <type> --annotation <dir> [--motif NAME]... [--enhancer FILE]
              [--strand-aware] [--cov-threshold 1] [--keep-going] [-o out.csv] SAMPLE...
```

`<dir>` holds `gc_windows.tsv[.gz]` and `motifs.tsv`. Output CSV, header
`sample,motif,deviation,expected_deviation`, rows ordered by sample (argument
order) then motif (manifest order), floats in Rust's shortest round-trip form,
missing values as `NA`. Logs go to stderr only.
**T41 integration (S).** `tests/cli.rs` with `assert_cmd`: the BATF golden through
the CLI (write `msites.tsv` out as a `bismarkcov` file first, as upstream's own
example does); two samples × two motifs (duplicate BATF under a second name)
checking order; byte-identical output across two runs; non-zero exit on a motif
with no hits; `NA` with `--keep-going`.
Done when: `cargo test --test cli` passes.

### M5 — differential tests against R

**T50 Case generator (H, needs T02).** `scripts/gen_differential_cases.R <seed> <n>`
writes `tests/fixtures/differential/case_NNN/` with portable inputs and Tier A
outputs. Each case draws at random: 1–3 chromosomes, 200–5 000 sites, site width
1 or 2, strands including `*`, coverage including 0-adjacent values; 50–2 000
TFBS with **mixed** widths (odd and even, so even `W` and half-integer midpoints
occur); window width 30 tiled or random and overlapping; `L` odd and even, 101..601;
enhancer present in a third of cases; `ignoreStrand` both ways. Cases where R
raises an error are kept, with the message in `error.txt`.
Commit 200 cases (compressed). CI does not run R.
**T51 `tests/differential_vs_r.rs` (S).** Run every case, compare at the section 3
tolerance, require an error wherever `error.txt` exists, print the worst case id
and max abs/rel error.
Done when: all 200 pass. Any **[confirm-A]** rule that fails here is fixed in
section 2 of this file first, then in code, with a regression test.

### M6 — release gate

**T60 (S).** `README.md` (usage, scope, parity statement, attribution),
`docs/divergences.md` complete, `docs/benchmarks.md` describing the M7 method
only (no numbers), CI green. Checklist: fmt, clippy, all tests, six parsers
tested, both BATF goldens, 200 differential cases, pinned SHA documented, no
**[confirm-A]** marker left unresolved in this file.

### M7 — performance (after v0.1 is tagged; not part of the gate)

**T70** criterion benches (`interval_overlap`, `gc_assignment`, `motif_deviation`,
`end_to_end`) on a synthetic genome-scale generator (≈28 M sites, ≈7.7 M GC
windows, 500 motifs × 250 k TFBS). **T71** R baselines with `threads = 1` and
`threads = N` on the same inputs, wall time and peak RSS via `/usr/bin/time -v`,
three runs each, medians. **T72** profile, then optimise in this order: load and
sort the methylome once per sample; struct-of-arrays per chromosome; replace
per-site binary searches by a two-pointer sweep; fixed `[sum; 5], [n; 5]`
accumulators with no per-hit allocation. **T73** `rayon` over motifs with a fixed
reduction order so output stays byte-identical to the sequential build; measure
1/2/4/8/16 threads. Every optimisation must leave `golden_batf` and
`differential_vs_r` passing unchanged. No speed-up is reported without the
measurement and the command that produced it.

## 6. Order and parallelism

```
T00 → T01 → T03 → T04(Tier B) → T05 → T10 → T11 → T12 → T20 → T21 → T22 → T23
                                                   T23 → T30…T36 → T40 → T41
T02 (any time, human-approved) → T04(Tier A) → T50 → T51 → T60 → M7
```

T30–T35 are independent of each other once T11 and T10 exist and may be given to
separate agents. Everything else is strictly sequential.
