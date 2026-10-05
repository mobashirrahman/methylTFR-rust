# Pinned reference

Every rule in `AGENT_PLAN.md` section 2 was read from the source below. Nothing
in this port is modelled on a later revision.

## Upstream

| field | value |
|---|---|
| repository | `https://github.com/EpigenomeInformatics/methylTFR` |
| pinned commit | `8aeab03cb469a9c213b7cd7fd2cbe6eeb8db5856` |
| commit subject | `Merge pull request #35 from EpigenomeInformatics/devel` |
| commit date | 2026-09-29T17:45:32+02:00 (author and committer) |
| `DESCRIPTION` `Version:` | 0.99.9 |
| `DESCRIPTION` `Date:` | 2026-09-29 |
| position relative to tags | 10 commits after tag `v1.0.0` (`git describe` → `v1.0.0-10-g8aeab03`); the tag itself is *not* what we pin, and the package is still versioned `0.99.9` (pre-1.0 Bioconductor development version) |
| licence | MIT (`DESCRIPTION`: `License: MIT + file LICENSE`) |
| licence file | `YEAR: 2025`, `COPYRIGHT HOLDER: Epigenome Informatics` |
| authors | Irem B. Gündüz (aut, cre), Sarath Kumar Murugan (aut), Fabian Mueller (aut) |
| funders | ERA-NET Transcan-Neu III - EPILUNAR (Grant 01KT2409); Saarland University (NanoBioMed Young Investigator Grant) |

Local clone: `reference/methylTFR` (gitignored). Recreate and verify with

```sh
git clone https://github.com/EpigenomeInformatics/methylTFR reference/methylTFR
git -C reference/methylTFR checkout 8aeab03cb469a9c213b7cd7fd2cbe6eeb8db5856
git -C reference/methylTFR rev-parse HEAD   # must print the pinned SHA
```

### Second revision, used for fixtures only

The published pkgdown vignette values in `AGENT_PLAN.md` section 4 were
rendered from an earlier fixture revision, not from the pinned commit. That
revision is pinned too, and is the source of the `batf_1d99721` fixtures only:

| field | value |
|---|---|
| commit | `1d99721b308d0e404270ff2fb06ca8ef5544c9a3` |
| subject | `update` |
| date | 2025-11-26 |

## R sources section 2 was derived from

| plan section | behaviour | upstream file (lines at the pinned SHA) |
|---|---|---|
| 2.2 | `read_methylome`, the six `parse_*` readers, `granges_helper`, the NaN / coverage / range post-steps | `R/data_reader.R` (`read_methylome` 26–80, `parse_epp` 83–91, `parse_bissnp` 94–103, `parse_allc` 106–116, `parse_bismarkcytosine` 119–131, `parse_bismarkcov` 134–145, `parse_encode` 148–174, `granges_helper` 186–199) |
| 2.3 | `addGCBintoMethylome`: `findOverlaps` per site/window pair, mean over *hits*, ordered by bin | `R/expected_deviations.R` 33–66 |
| 2.4 | `computeDeviation`: `resize(width(tfbs)[1] + 130, fix = "center")`, enhancer `subsetByOverlaps`, `findOverlaps(type = "within")`, `mid_point`, `x = start(msites) - mid` | `R/compute_deviations.R` 89–132 (resize 99–101, enhancer 103–107, hits 108–111, midpoint 120–121, profile 122–125) |
| 2.5 | strand handling: `ignore.strand = ignoreStrand` on every `findOverlaps` / `subsetByOverlaps` | `R/expected_deviations.R` 51–53, `R/compute_deviations.R` 104–110, `R/methyltfr_core.R` 344–346 |
| 2.6 | enhancer reduces the GC windows once per run, before bin means are computed | `R/methyltfr_core.R` 343–346, passed to `addGCBintoMethylome` at 217 and to `computeDeviation` at 229–230 |
| 2.7 | `computeExpectations`: `t(gcfreq) %*% binMsites[, 2]`, `round(seq(-floor(L/2), floor(L/2), length.out = L))` | `R/expected_deviations.R` 82–100 |
| 2.8 | `dev_helper`: `cut(x, c(-250,-200,-25,25,200,250))`, `n`/mean per interval, `mean[(k+1) %/% 2] / ((mean[1] + mean[k]) / 2)`, `NA` when `k == 0` | `R/compute_deviations.R` 139–152 |
| 2.8 (run level) | `valid_core_motifs` drops motifs with empty TFBS or missing matrix and aborts when none remain; sample id is the file basename; per-sample per-motif loop | `R/methyltfr_core.R` 138–160 (`valid_core_motifs`), 203–251 (`process_core_sample`); `R/run_methyltfr.R` 202 (`sample_ids = basename(files_list)`) |
| 2.1 | R numeric semantics (`round`, `round(x, 6)`, `seq`, `cut`, `%/%`) | not upstream code — R 4.3.3 base semantics, verified by running R; see below |

Files deliberately **not** read for section 2, because they are out of scope for
v0.1 (plan 2.9): `R/plots.R`, `R/plot_helpers.R`, `R/differential_analysis.R`,
`R/variability.R`, `R/rnbeads_interface.R`, `R/class.R`, `R/memory_helpers.R`,
`R/zzz.R`.

One reading worth recording, because the plan's wording and the source differ in
emphasis: `valid_core_motifs` takes the motif order from `names(gcfreqs)`
(`R/methyltfr_core.R` 139), not from `names(tf_bindsites)`. In the portable
format of plan section 3 both names come from the same `motifs.tsv` manifest
row, so "manifest order" (plan 2.8) is the same order.

## R used for the Tier A / Tier B oracles

Two R installations are used, and it matters which script needs which.

**Tier B** (`scripts/reference_base.R`, `scripts/gen_rmath_tables.R`,
`scripts/gen_differential_cases.R` helpers) is base R only, so it runs on the
system R with no extra packages:

```
R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
platform: x86_64-pc-linux-gnu
```

**Tier A** (`scripts/run_reference.R`, `scripts/export_reference_data.R` on the
`.rda` route, and the T02/T04/T50 oracles) needs Bioconductor 3.18, so it runs
in a micromamba environment (`T02`) rather than in the system library: the box
has no `sudo`, and its system R is missing the development headers
(`libcurl4-openssl-dev`, `libxml2-dev`, `libhdf5-dev`) that a from-source
Bioconductor build needs.

- micromamba 2.9.0, `MAMBA_ROOT_PREFIX=/scratch/mdra00001/micromamba`
- environment `methyltfr-ref`, channel `conda-forge` only
- outside the repository, so nothing in `.gitignore` has to cover it

Recreate it:

```sh
curl -sSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=/scratch/mdra00001/micromamba
micromamba create -y -n methyltfr-ref -c conda-forge --strict-channel-priority \
  r-base=4.3 \
  bioconductor-genomicranges bioconductor-iranges bioconductor-s4vectors \
  bioconductor-summarizedexperiment bioconductor-biocparallel \
  bioconductor-delayedarray bioconductor-hdf5array \
  r-data.table r-r.utils r-logger r-matrixstats r-stringr r-ggplot2
micromamba run -n methyltfr-ref R CMD INSTALL reference/methylTFR
```

Resolved versions of the packages the plan names, which together identify
Bioconductor 3.18 (the release series for R 4.3):

| plan | package | version |
|---|---|---|
| — | `r-base` | 4.3.3 |
| `GenomicRanges` | GenomicRanges | 1.54.1 |
| `IRanges` | IRanges | 2.36.0 |
| `S4Vectors` | S4Vectors | 0.40.2 |
| `SummarizedExperiment` | SummarizedExperiment | 1.32.0 |
| `BiocParallel` | BiocParallel | 1.36.0 |
| `DelayedArray` | DelayedArray | 0.28.0 |
| `HDF5Array` | HDF5Array | 1.30.0 |
| — | BiocGenerics | 0.48.1 |
| — | GenomeInfoDb | 1.38.1 |
| — | rhdf5 | 2.46.1 |
| `data.table` | data.table | 1.17.8 |
| `R.utils` | R.utils | 2.13.0 |
| `logger` | logger | 0.4.0 |
| `matrixStats` | matrixStats | 1.5.0 |
| `stringr` | stringr | 1.5.2 |
| `ggplot2` | ggplot2 | 3.5.2 |
| — | methylTFR | 0.99.9 (installed from the pinned clone) |

`data.table` resolves to 1.17.8 here, not the CRAN head, because `r-base=4.3`
pins the whole conda-forge dependency closure; upstream only requires
`>= 1.14.0`, so the pin is not a source of divergence.

### `sessionInfo()` of the Tier A environment

```
R version 4.3.3 (2024-02-29)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Ubuntu 24.04.5 LTS

Matrix products: default
BLAS/LAPACK: /scratch/mdra00001/micromamba/envs/methyltfr-ref/lib/libopenblasp-r0.3.34.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] methylTFR_0.99.9            SummarizedExperiment_1.32.0
 [3] Biobase_2.62.0              GenomicRanges_1.54.1
 [5] GenomeInfoDb_1.38.1         IRanges_2.36.0
 [7] S4Vectors_0.40.2            BiocGenerics_0.48.1
 [9] MatrixGenerics_1.14.0       matrixStats_1.5.0
[11] data.table_1.17.8

loaded via a namespace (and not attached):
 [1] SparseArray_1.2.2       bitops_1.0-9            stringi_1.8.7
 [4] lattice_0.22-7          magrittr_2.0.3          grid_4.3.3
 [7] RColorBrewer_1.1-3      R.oo_1.27.1             Matrix_1.6-5
[10] R.utils_2.13.0          scales_1.4.0            HDF5Array_1.30.0
[13] codetools_0.2-20        abind_1.4-5             cli_3.6.5
[16] rlang_1.1.6             crayon_1.5.3            XVector_0.42.0
[19] R.methodsS3_1.8.2       DelayedArray_0.28.0     S4Arrays_1.2.0
[22] tools_4.3.3             parallel_4.3.3          BiocParallel_1.36.0
[25] Rhdf5lib_1.24.0         ggplot2_3.5.2           GenomeInfoDbData_1.2.11
[28] vctrs_0.6.5             logger_0.4.0            R6_2.6.1
[31] rhdf5_2.46.1            lifecycle_1.0.4         zlibbioc_1.48.0
[34] stringr_1.5.2           pkgconfig_2.0.3         pillar_1.11.0
[37] gtable_0.3.6            glue_1.8.0              tibble_3.3.0
[40] rhdf5filters_1.14.1     farver_2.1.2            compiler_4.3.3
[43] RCurl_1.98-1.17
```

CI remains Rust-only and never starts R.

## Attestation

The `[verified]` markers in `AGENT_PLAN.md` section 2.1 and 2.7 were confirmed
on R 4.3.3. The `[confirm-A]` markers (strand validation in `read_methylome`,
and the resize/midpoint rules being strand-independent) are source-reading
conclusions only; the Tier A environment documented above now exists, so they
close as soon as T04 and the differential cases (T50/T51) run against it.

### Closed `[confirm-A]` markers

Both markers `AGENT_PLAN.md` section 2 left for the Tier A oracle have been
closed by running methylTFR 0.99.9 itself; the plan text stays as written, but
the behaviour below is now observed rather than read.

| marker | how it was closed | observed behaviour |
|---|---|---|
| 2.2, strand validation | `read_methylome()` on a one-line EPP file whose 6th column is `.`, then `x`, then `+` | `strand values must be in '+' '-' '*'` for both `.` and `x`; `+` succeeds. Strand `*` is accepted, so a Rust reader that rejects anything outside `{+, -, *}` is faithful. |
| 2.4, resize and midpoint are strand-independent | `identical()` between `resize(tfbs, W, fix = "center")` and the plan's formula `new_start = start + floor((width - W) / 2)`, `new_end = new_start + W - 1` over all 268 717 BATF TFBS, then `findOverlaps(..., type = "within", ignore.strand = FALSE)` | starts identical for every TFBS, ends differ from the plan's formula by 0 for every TFBS (all resized widths are exactly 541), and `findOverlaps` with `ignore.strand = FALSE` returns 29 hits -- the same count as with strand ignored, because the fixture's sites and resized TFBS agree in this subset. The resize never reads `strand`. |

## T04 oracle agreement

`scripts/compare_oracles.R` compares the two committed oracle outputs:

```
$ Rscript scripts/compare_oracles.R tests/fixtures/batf/expected tests/fixtures/batf/expected_tierB
expected_bins.tsv         5 rows  max abs 2.220e-16  max rel 2.833e-16  ok
observed_profile.tsv     29 rows  max abs 0.000e+00  max rel 0.000e+00  ok
expected_profile.tsv   512 rows  max abs 2.220e-16  max rel 3.208e-16  ok
expected_dev.tsv         1 rows  max abs 2.220e-16  max rel 2.266e-16  ok
tests/fixtures/batf/expected vs tests/fixtures/batf/expected_tierB: worst abs 2.220e-16 (expected_bins.tsv), worst rel 3.208e-16, tolerance 1.0e-10 -- PASS

$ Rscript scripts/compare_oracles.R tests/fixtures/batf_1d99721/expected tests/fixtures/batf_1d99721/expected_tierB
expected_bins.tsv         5 rows  max abs 2.220e-16  max rel 2.833e-16  ok
observed_profile.tsv     29 rows  max abs 0.000e+00  max rel 0.000e+00  ok
expected_profile.tsv   512 rows  max abs 3.331e-16  max rel 4.854e-16  ok
expected_dev.tsv         1 rows  max abs 2.220e-16  max rel 2.257e-16  ok
tests/fixtures/batf_1d99721/expected vs tests/fixtures/batf_1d99721/expected_tierB: worst abs 3.331e-16 (expected_profile.tsv), worst rel 4.854e-16, tolerance 1.0e-10 -- PASS
```

The two oracles are not bit-identical and cannot be: R's `mean()` accumulates in
long double and `%*%` runs through BLAS, while the Tier B loops sum in `f64`.
`observed_profile.tsv` -- the only file with no accumulation in it, because every
value is copied straight out of a site -- is byte-identical. The worst
disagreement anywhere is 3.3e-16, six orders of magnitude below the 1e-10 the
plan requires.

### Golden values (AGENT_PLAN.md section 4), recomputed

`tests/fixtures/batf_1d99721/expected_tierB/expected_dev.tsv` reproduces the
published pkgdown column to all 7 printed digits:

```
motif	obs_d	exp_d	dev
BATF	2.7272727272727271	0.98359851852227742	1.7436742087504498
```

which is `dev = 1.743674` and `exp_dev = 0.9835985`, and the Tier A oracle
agrees. `tests/fixtures/batf/expected/expected_dev.tsv` reproduces the
provisional pinned-SHA column, `obs_d = 2.7272727`, `exp_dev = 0.9798459`,
`dev = 1.7474268`, again with Tier A agreeing, so the provisional column is
confirmed rather than merely unrefuted.
