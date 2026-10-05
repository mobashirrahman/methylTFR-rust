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

```
R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
platform: x86_64-pc-linux-gnu
```

`sessionInfo()` of the R installation that ran the oracle scripts will be
appended below by task T02, once Bioconductor is installed into `R_LIBS_USER`
and `methylTFR` itself is installed from `reference/methylTFR`. Until that
exists, only the Tier B (base-R) oracle can be run; CI is Rust-only and never
starts R.

## Attestation

The `[verified]` markers in `AGENT_PLAN.md` section 2.1 and 2.7 were confirmed
on this R 4.3.3. The `[confirm-A]` markers (strand validation in
`read_methylome`, and the resize/midpoint rules being strand-independent) are
source-reading conclusions only and remain open until Tier A (T04) and the
differential cases (T50/T51) run.
