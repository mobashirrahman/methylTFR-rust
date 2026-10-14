#!/usr/bin/env Rscript

# run_reference.R -- T04 **Tier A** reference oracle, command-line entry point.
#
#   usage: run_reference.R <fixture_dir> <out_dir>
#
# Writes the four reference output files for one fixture directory by calling the
# real methylTFR 0.99.9 package. The work is in reference_common.R, which
# gen_differential_cases.R also uses, so the BATF fixtures and the differential
# cases go through one code path.
#
#   expected_bins.tsv      gc_bin mean n_hits
#   observed_profile.tsv   x value
#   expected_profile.tsv   x value
#   expected_dev.tsv       motif obs_d exp_d dev
#
# Compare against reference_base.R with scripts/compare_oracles.R.

source("scripts/reference_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: run_reference.R <fixture_dir> <out_dir>")
}
fixture_dir <- args[[1]]
out_dir <- args[[2]]

run_reference_case(fixture_dir, out_dir)

# run_reference_case() already logged the summary line; nothing else to do here.
