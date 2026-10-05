#!/usr/bin/env Rscript

# compare_csv.R -- T71's output check.
#
#   usage: compare_csv.R <r_csv> <rust_csv> [tolerance]
#
# The two implementations do the same work on the same files, so their outputs
# must agree. They are not expected to be byte-identical: R's `mean()` accumulates
# in long double and its `%*%` goes through BLAS, while this port sums in f64, so
# the last digit or two can differ. The tolerance is the section 3 one.
#
# Exits non-zero if the two disagree by more than that, which is what keeps the
# benchmark from ever reporting a speed-up between two programs that compute
# different numbers.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: compare_csv.R <r_csv> <rust_csv> [tolerance]")
}
tol <- if (length(args) >= 3L) as.numeric(args[[3]]) else 1e-10

a <- read.delim(args[[1]], sep = ",", colClasses = "character", check.names = FALSE)
b <- read.delim(args[[2]], sep = ",", colClasses = "character", check.names = FALSE)
want <- c("sample", "motif", "deviation", "expected_deviation")
stopifnot(
    identical(as.character(names(a)), want),
    identical(as.character(names(b)), want),
    identical(a$sample, b$sample),
    identical(a$motif, b$motif),
    nrow(a) == nrow(b)
)

if (identical(a, b)) {
    cat(sprintf("byte-identical over %d rows\n", nrow(a)))
    quit(status = 0)
}

worst <- 0
worst_col <- ""
worst_row <- 0L
for (col in c("deviation", "expected_deviation")) {
    x <- as.numeric(a[[col]])
    y <- as.numeric(b[[col]])
    if (!identical(is.na(x), is.na(y))) {
        cat(sprintf("FAIL: NA placement differs in %s\\n", col))
        quit(status = 1)
    }
    d <- abs(x - y)
    d[is.na(d)] <- 0
    if (max(d) > worst) {
        worst <- max(d)
        worst_col <- col
        worst_row <- which.max(d)
    }
}
cat(sprintf(
    "not byte-identical: worst absolute difference %.3e in %s (row %d, motif %s), \
tolerance %.1e\n",
    worst, worst_col, worst_row, a$motif[[worst_row]], tol
))
if (!(worst <= tol)) {
    cat("FAIL: outside the section 3 tolerance\n")
    quit(status = 1)
}
