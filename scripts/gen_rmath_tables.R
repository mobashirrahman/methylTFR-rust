#!/usr/bin/env Rscript

# gen_rmath_tables.R -- T05, generates the R numeric oracle tables that
# tests/rmath.rs is checked against.
#
#   usage: gen_rmath_tables.R [out_dir]
#
# Base R only, and deterministic: the same seed always produces the same files,
# so a regeneration that changes a byte means the semantics changed.
#
# Files written to tests/fixtures/rmath/:
#
#   round6.tsv     x, round(x, 6)   -- parser scores and the 1/128 tie
#   round0.tsv     x, round(x)      -- half-to-even ties
#   seq.tsv        L, comma-joined round(seq(-floor(L/2), floor(L/2), length.out = L))
#   cut.tsv        x, interval index (1-based, as R factor levels) or NA
#   midpoint.tsv   start, end, round(end + (start - end) / 2)
#
# Every float is written with sprintf("%.17g", x) so that it round-trips.

set.seed(20240929)

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1L) {
    args[[1]]
} else {
    "tests/fixtures/rmath"
}
if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("could not create ", out_dir)
}

fmt <- function(x) sprintf("%.17g", x)
fmt_int <- function(x) sprintf("%d", x)

write_table <- function(lines, name) {
    path <- file.path(out_dir, name)
    writeLines(lines, path)
    cat(sprintf("%-16s %8d rows\n", name, length(lines) - 1L))
}

# ------------------------------------------------------------------ round6 --

x <- numeric(0)

# Every methylation fraction a parser can produce: m/cov with cov <= 300.
# This is the exact value set of round(msites$V5 / 1000, 6) and friends once
# the divisons are reduced, so it is the set that matters for parity.
for (cov in 1:300) {
    x <- c(x, (0:cov) / cov)
}

# Every dyadic rational with denominator 2^n, n <= 20, in its reduced form.
# round(x, 6) is at a tie exactly when x * 1e6 is an odd half-integer, which is
# how 1/128 becomes 0.007812 instead of 0.007813; the whole dyadic family is
# the cheapest way to cover that boundary densely.
for (n in 1:20) {
    x <- c(x, seq(1, 2^n - 1, by = 2) / 2^n)
}

# Random percentages, the shape of a bismark / encode score.
x <- c(x, sample.int(100L, 5000L, replace = TRUE) / 100)

# The full grid of 0:1000/1000.
x <- c(x, (0:1000) / 1000)

x <- unique(x)
write_table(c("x\trounded", paste(fmt(x), fmt(round(x, 6)), sep = "\t")),
            "round6.tsv")

# ------------------------------------------------------------------ round0 --

# Half-to-even ties on both sides, plus random doubles so that the non-tie
# branches are covered too.
ties <- c(rep(seq(-2000, 2000, by = 1) + 0.5, 1L), seq(-2000, 2000) - 0.5)
rnd <- runif(20000L, -5000, 5000)
all0 <- c(ties, rnd, c(0, -0, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5,
                       0.49999999999999994, 4503599627370496.5))
write_table(
    c("x\trounded", paste(fmt(all0), fmt(round(all0)), sep = "\t")),
    "round0.tsv"
)

# --------------------------------------------------------------------- seq --

# seq() is only ever asked for the profile grid, i.e. from = -floor(L/2) to
# +floor(L/2), for every L from 1 to 1200.  Even L skips 0; that falls out of
# the table rather than being asserted separately.
Ls <- 1:1200
rows <- vapply(Ls, function(L) {
    h <- floor(L / 2)
    v <- round(seq(-h, h, length.out = L))
    paste0(
        L, "\t",
        paste(if (L == 0L) "" else fmt(v), collapse = ",")
    )
}, character(1))
write_table(c("L\tpositions", rows), "seq.tsv")

# --------------------------------------------------------------------- cut --

# cut(x, c(-250, -200, -25, 25, 200, 250)) is right-closed, so the breaks
# themselves matter; include them, the values just inside and outside, and the
# halves.
cut_x <- c(
    (-260:260),
    (-260:260) + 0.5,
    (-260:260) - 0.5,
    c(-250.0000001, -249.9999999, 250.0000001, 249.9999999)
)
cut_x <- cut_x[order(cut_x)]
f <- cut(cut_x, c(-250, -200, -25, 25, 200, 250))
lvl <- as.integer(f)
lvl[is.na(lvl)] <- NA_integer_
write_table(
    c("x\tinterval", paste(fmt(cut_x), fmt_int(lvl), sep = "\t")),
    "cut.tsv"
)

# --------------------------------------------------------------- midpoint --

# mid_point = round(end + (start - end) / 2), computed in doubles and only then
# rounded to an integer.  For even widths the exact result is a half-integer,
# which is where half-to-even shows up.
n <- 4000L
w <- sample(1:400, n, replace = TRUE)
s <- sample.int(10^7, n, replace = TRUE) - 5 * 10^6
e <- s + w - 1L
# Make sure both parities of (end - start) are represented, including the ones
# that put the exact half on a .5 boundary.
e2 <- e + sample(c(-1L, 0L, 1L), n, replace = TRUE)
s_all <- c(s, s)
e_all <- c(e, e2)
mid <- round(e_all + ((s_all - e_all) / 2))
write_table(
    c("start\tend\tmid",
      paste(fmt_int(s_all), fmt_int(e_all), fmt(mid), sep = "\t")),
    "midpoint.tsv"
)

cat("done:", normalizePath(out_dir), "\n")