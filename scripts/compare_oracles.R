#!/usr/bin/env Rscript

# compare_oracles.R -- T04 acceptance check.
#
#   usage: compare_oracles.R <dir_a> <dir_b> [tolerance]
#
# Compares two reference-output directories (Tier A `expected/` against Tier B
# `expected_tierB/`) and prints the largest absolute and relative difference per
# file, then exits non-zero if anything exceeds the tolerance.
#
# The tolerance default is the AGENT_PLAN.md section 3 value, 1e-10.  The two
# oracles do not agree bit-for-bit and cannot: R's `mean()` accumulates in
# long double and `%*%` goes through BLAS, while the Tier B loops sum in f64.
# The observed differences are one or two ulp, which is why the plan compares
# with a tolerance instead of `diff`.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: compare_oracles.R <dir_a> <dir_b> [tolerance]")
}
dir_a <- args[[1]]
dir_b <- args[[2]]
tol <- if (length(args) >= 3L) as.numeric(args[[3]]) else 1e-10

FILES <- c(
    "expected_bins.tsv", "observed_profile.tsv",
    "expected_profile.tsv", "expected_dev.tsv"
)

worst_abs <- 0
worst_rel <- 0
worst_file <- NA_character_
failed <- FALSE

for (f in FILES) {
    a <- utils::read.delim(file.path(dir_a, f), colClasses = "character")
    b <- utils::read.delim(file.path(dir_b, f), colClasses = "character")
    if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
        cat(sprintf("%-22s SHAPE MISMATCH\n", f))
        failed <- TRUE
        next
    }
    bad <- 0L
    file_abs <- 0
    file_rel <- 0
    for (col in names(a)) {
        ca <- suppressWarnings(as.numeric(a[[col]]))
        cb <- suppressWarnings(as.numeric(b[[col]]))
        if (any(is.na(ca) != is.na(cb))) {
            cat(sprintf("%-22s %-12s NA placement differs\n", f, col))
            failed <- TRUE
            next
        }
        ok <- !is.na(ca)
        if (!any(ok)) {
            next
        }
        # Text columns (motif name) must match exactly.
        if (all(is.na(suppressWarnings(as.numeric(a[[col]][!ok]))) &
                is.na(suppressWarnings(as.numeric(b[[col]][!ok])))) &&
                any(!ok)) {
            if (!identical(a[[col]], b[[col]])) {
                cat(sprintf("%-22s %-12s text column differs\n", f, col))
                failed <- TRUE
            }
            next
        }
        d <- abs(ca[ok] - cb[ok])
        r <- ifelse(cb[ok] == 0, 0, d / abs(cb[ok]))
        file_abs <- max(file_abs, max(d))
        file_rel <- max(file_rel, max(r))
        bad <- bad + sum(!(d <= tol | r <= tol))
    }
    if (file_abs > worst_abs) {
        worst_abs <- file_abs
        worst_file <- f
    }
    worst_rel <- max(worst_rel, file_rel)
    if (bad > 0L) {
        failed <- TRUE
    }
    cat(sprintf(
        "%-22s %4d rows  max abs %.3e  max rel %.3e  %s\n",
        f, nrow(a), file_abs, file_rel,
        if (bad > 0L) "FAIL" else "ok"
    ))
}

cat(sprintf(
    "%s vs %s: worst abs %.3e (%s), worst rel %.3e, tolerance %.1e -- %s\n",
    dir_a, dir_b, worst_abs, worst_file, worst_rel, tol,
    if (failed) "FAIL" else "PASS"
))
if (failed) {
    quit(status = 1L)
}