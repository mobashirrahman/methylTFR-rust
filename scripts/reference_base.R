#!/usr/bin/env Rscript

# reference_base.R -- T04 **Tier B** reference oracle.
#
#   usage: reference_base.R <fixture_dir> <out_dir>
#
# An independent implementation of AGENT_PLAN.md section 2 in **base R only**:
# it reads the portable fixture files and recomputes the four reference output
# files from the rules in the plan.  It never loads methylTFR, so agreement
# between this script and run_reference.R (Tier A) is evidence about section 2
# rather than about a single implementation.
#
# Output files (AGENT_PLAN.md section 5, T04):
#
#   expected_bins.tsv      gc_bin mean n_hits
#   observed_profile.tsv   x value
#   expected_profile.tsv   x value
#   expected_dev.tsv       motif obs_d exp_d dev
#
# Every float is written with sprintf("%.17g", x) so that it round-trips
# exactly.  The loops are deliberately quadratic; this script is the oracle,
# not a competitor.

# ---------------------------------------------------------------- section 2.1

# R's round(): nearest integer, ties to even, for finite x.  NAs pass through.
round_half_even <- function(x) {
    na <- is.na(x)
    x <- x[!na]
    fl <- floor(x)
    frac <- x - fl
    out <- fl
    out[frac > 0.5] <- out[frac > 0.5] + 1
    half <- frac == 0.5
    # Exact tie: the even neighbour wins, so 0.5 -> 0 and 1.5 -> 2.
    out[half] <- ifelse(fl[half] %% 2 == 0, fl[half], fl[half] + 1)
    res <- rep(NA_real_, sum(na) + length(x))
    res[seq_along(x)] <- out
    res
}

# R's seq(from, to, length.out = n), as used for the profile grid.
seq_len_out <- function(from, to, n) {
    if (n == 0L) {
        return(numeric(0))
    }
    if (n == 1L) {
        return(from)
    }
    if (n > 2L && from == to) {
        return(rep(from, n))
    }
    by <- (to - from) / (n - 1L)
    out <- from + (seq_len(n) - 1L) * by
    out[n] <- to
    out
}

# R's cut(x, c(-250, -200, -25, 25, 200, 250)): open on the left, closed on the
# right, so a break belongs to the interval above it.  Returns 1-based indices,
# or NA when x falls outside the closed range of the breaks.
CUT_BREAKS <- c(-250, -200, -25, 25, 200, 250)

cut_index <- function(x) {
    idx <- integer(length(x))
    # The interval index is the number of breaks strictly below x.  That puts
    # each break itself in the interval above it, which is what cut()'s
    # right-closed labels say.
    for (b in CUT_BREAKS) {
        below <- !is.na(x) & x > b
        idx[below] <- idx[below] + 1L
    }
    idx[idx < 1L | idx > length(CUT_BREAKS) - 1L] <- NA_integer_
    idx
}

# ------------------------------------------------------------------- readers

fmt_num <- function(x) {
    ifelse(
        is.na(x) & !is.nan(x),
        "NA",
        ifelse(
            is.finite(x) & x == trunc(x) & abs(x) < 1e15,
            sprintf("%.0f", x),
            sprintf("%.17g", x)
        )
    )
}

read_portable <- function(path, header = TRUE) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con))
    read.delim(
        con,
        header = header,
        stringsAsFactors = FALSE,
        colClasses = "character",
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )
}

write_lines_to <- function(lines, path) {
    con <- if (grepl("\\.gz$", path)) {
        gzfile(path, "wt")
    } else {
        file(path, "wt")
    }
    on.exit(close(con))
    writeLines(lines, con)
}

# Strand compatibility, AGENT_PLAN.md section 2.5.
strand_ok <- function(a, b, ignore_strand) {
    if (ignore_strand) {
        return(rep(TRUE, length(a)))
    }
    a == b | a == "*" | b == "*"
}

# Brute-force "keep rows of x overlapping any row of q", i.e. IRanges'
# subsetByOverlaps(q, x) / findOverlaps type "any".
reduce <- function(x, q, ignore_strand = TRUE) {
    keep <- logical(nrow(x))
    for (i in seq_len(nrow(q))) {
        hit <- which(
            x$chr == q$chr[i] &
                x$start <= q$end[i] &
                x$end >= q$start[i]
        )
        if (length(hit) > 0L) {
            hit <- hit[strand_ok(
                x$strand[hit], rep(q$strand[i], length(hit)),
                ignore_strand
            )]
            keep[hit] <- TRUE
        }
    }
    keep
}

# ----------------------------------------------------------------- section 2.3

# addGCBintoMethylome: the mean is over *overlap hits*, not over sites, so a
# site straddling two abutting windows contributes twice.
add_gc_bins <- function(msites, gcdist, ignore_strand = TRUE) {
    bins <- integer(0)
    scores <- numeric(0)
    for (i in seq_len(nrow(msites))) {
        s <- msites$start[i]
        e <- msites$end[i]
        hit <- which(
            gcdist$chr == msites$chr[i] &
                gcdist$start <= e & gcdist$end >= s
        )
        # hit order: windows by ascending start, ties by original index
        hit <- hit[order(gcdist$start[hit], hit)]
        if (length(hit) > 0L) {
            hit <- hit[strand_ok(
                rep(msites$strand[i], length(hit)),
                gcdist$strand[hit],
                ignore_strand
            )]
        }
        if (length(hit) > 0L) {
            bins <- c(bins, gcdist$gc_bin[hit])
            scores <- c(scores, rep(msites$score[i], length(hit)))
        }
    }
    if (length(bins) == 0L) {
        stop("No methylation sites found in the GC distribution")
    }
    ub <- sort(unique(bins))
    list(
        bin = ub,
        mean = vapply(ub, function(b) mean(scores[bins == b]), numeric(1)),
        n_hits = vapply(ub, function(b) sum(bins == b), integer(1))
    )
}

# ----------------------------------------------------------------- section 2.4

# computeDeviation: resize every TFBS to width W about its centre, keep the
# sites *within* a TFBS, and turn each hit into (x = site.start - midpoint,
# value = site.score).
compute_observed <- function(msites, tfbs, ignore_strand = TRUE, motif) {
    w <- (tfbs$end[1] - tfbs$start[1] + 1L) + 130L
    new_start <- tfbs$start + floor((tfbs$end - tfbs$start + 1L - w) / 2L)
    new_end <- new_start + w - 1L
    mid <- round_half_even(
        as.numeric(new_end) + (as.numeric(new_start) - as.numeric(new_end)) / 2
    )

    xs <- numeric(0)
    vals <- numeric(0)
    for (i in seq_len(nrow(msites))) {
        s <- msites$start[i]
        e <- msites$end[i]
        hit <- which(
            tfbs$chr == msites$chr[i] & new_start <= s & new_end >= e
        )
        hit <- hit[order(new_start[hit], hit)]
        if (length(hit) > 0L) {
            hit <- hit[strand_ok(
                rep(msites$strand[i], length(hit)),
                tfbs$strand[hit],
                ignore_strand
            )]
        }
        if (length(hit) > 0L) {
            xs <- c(xs, s - mid[hit])
            vals <- c(vals, rep(msites$score[i], length(hit)))
        }
    }
    if (length(xs) == 0L) {
        stop("No methylation sites found in the ", motif, " binding sites")
    }
    list(x = xs, value = vals)
}

# ----------------------------------------------------------------- section 2.7

compute_expectations <- function(gcfreq, bin_mean) {
    L <- ncol(gcfreq)
    if (nrow(gcfreq) != length(bin_mean)) {
        stop("non-conformable arguments")
    }
    expected <- numeric(L)
    for (j in seq_len(L)) {
        acc <- 0
        for (i in seq_along(bin_mean)) {
            acc <- acc + gcfreq[i, j] * bin_mean[i]
        }
        expected[j] <- acc
    }
    h <- floor(L / 2)
    list(x = round_half_even(seq_len_out(-h, h, L)), value = expected)
}

# ----------------------------------------------------------------- section 2.8

dev_helper <- function(x, value) {
    idx <- cut_index(x)
    keep <- !is.na(idx)
    idx <- idx[keep]
    value <- value[keep]
    if (length(idx) == 0L) {
        return(NA_real_)
    }
    ub <- sort(unique(idx))
    means <- vapply(ub, function(b) mean(value[idx == b]), numeric(1))
    k <- length(ub)
    if (k == 0L) {
        return(NA_real_)
    }
    # mean[(k + 1) %/% 2] / ((mean[1] + mean[k]) / 2), integer division.
    means[(k + 1L) %/% 2L] / ((means[1] + means[k]) / 2)
}

# ------------------------------------------------------------------- the run

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: reference_base.R <fixture_dir> <out_dir>")
}
fixture_dir <- args[[1]]
out_dir <- args[[2]]

if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("could not create ", out_dir)
}

read_fixture <- function(name) {
    d <- read_portable(file.path(fixture_dir, name))
    for (col in c("start", "end")) {
        if (col %in% names(d)) {
            d[[col]] <- as.integer(d[[col]])
        }
    }
    d
}

msites <- read_fixture("msites.tsv")
msites$score <- as.numeric(msites$score)
gcdist <- read_fixture("gc_windows.tsv")
gcdist$gc_bin <- as.integer(gcdist$gc_bin)
motifs <- read_fixture("motifs.tsv")
stopifnot(nrow(motifs) >= 1L)
motif <- motifs$motif[[1]]

manifest_path <- function(p) {
    if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(fixture_dir, p)
}
tfbs <- read_portable(manifest_path(motifs$tfbs_path[[1]]))
tfbs$start <- as.integer(tfbs$start)
tfbs$end <- as.integer(tfbs$end)
gcfreq <- as.matrix(
    read_portable(manifest_path(motifs$gcfreq_path[[1]]), header = FALSE)
)
storage.mode(gcfreq) <- "double"

ignore_strand <- TRUE
opt_path <- file.path(fixture_dir, "options.tsv")
if (file.exists(opt_path)) {
    opts <- read_portable(opt_path)
    if ("ignore_strand" %in% names(opts)) {
        ignore_strand <- as.logical(opts$ignore_strand[[1]])
    }
}

# AGENT_PLAN.md section 2.6: the enhancer reduces the GC windows once, before
# any bin mean is computed, and filters the resized TFBS.
enhancer_path <- file.path(fixture_dir, "enhancer.tsv")
if (file.exists(enhancer_path)) {
    enhancer <- read_fixture("enhancer.tsv")
    gcdist <- gcdist[reduce(gcdist, enhancer, ignore_strand), , drop = FALSE]
    tfbs <- tfbs[reduce(tfbs, enhancer, ignore_strand), , drop = FALSE]
}

bins <- add_gc_bins(msites, gcdist, ignore_strand)
observed <- compute_observed(msites, tfbs, ignore_strand, motif)
expected <- compute_expectations(gcfreq, bins$mean)

obs_d <- dev_helper(observed$x, observed$value)
exp_d <- dev_helper(expected$x, expected$value)

write_lines_to(c(
    "gc_bin\tmean\tn_hits",
    paste(
        fmt_num(bins$bin), fmt_num(bins$mean),
        fmt_num(bins$n_hits), sep = "\t"
    )
), file.path(out_dir, "expected_bins.tsv"))

write_lines_to(c(
    "x\tvalue",
    paste(fmt_num(observed$x), fmt_num(observed$value), sep = "\t")
), file.path(out_dir, "observed_profile.tsv"))

write_lines_to(c(
    "x\tvalue",
    paste(fmt_num(expected$x), fmt_num(expected$value), sep = "\t")
), file.path(out_dir, "expected_profile.tsv"))

write_lines_to(c(
    "motif\tobs_d\texp_d\tdev",
    paste(
        motif, fmt_num(obs_d), fmt_num(exp_d), fmt_num(obs_d - exp_d),
        sep = "\t"
    )
), file.path(out_dir, "expected_dev.tsv"))

cat(sprintf(
    "%s: motif %s obs_d %.7g exp_dev %.7g dev %.7g (%d bins, %d observed hits)\n",
    out_dir, motif, obs_d, exp_d, obs_d - exp_d,
    length(bins$bin), length(observed$x)
))