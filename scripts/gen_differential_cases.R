#!/usr/bin/env Rscript

# gen_differential_cases.R -- T50, the differential case generator.
#
#   usage: gen_differential_cases.R <seed> <n> [out_dir]
#
# Writes `tests/fixtures/differential/case_NNN/` for N = 1..n, each a complete
# portable annotation plus the Tier A oracle outputs for it (or `error.txt` when
# upstream raises). `tests/differential_vs_r.rs` then runs every case through the
# Rust port and compares at the AGENT_PLAN.md section 3 tolerance, requiring an
# error wherever `error.txt` exists.
#
# Needs the Tier A environment (docs/reference-version.md). CI never runs this; it
# runs the committed cases, which is why they are committed.
#
# Each case draws at random (AGENT_PLAN.md T50):
#
#   * 1-3 chromosomes
#   * 200-5000 sites, width 1 or 2, strands including `*`, coverage near 0
#   * GC windows 30 bp wide, tiled abutting or randomly overlapping
#   * 1-5 GC bins, sometimes with one bin left empty so the bin-count mismatch of
#     section 2.7 fires
#   * 50-2000 TFBS with **mixed** widths, so even W and half-integer midpoints occur
#   * L odd and even, 101..601
#   * an enhancer in a third of the cases
#   * ignoreStrand both ways
#
# Every case is reproducible from its index: case N uses `set.seed(seed + N)`, so
# regenerating any single case gives the same bytes.

source("scripts/reference_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: gen_differential_cases.R <seed> <n> [out_dir]")
}
seed <- as.integer(args[[1]])
n <- as.integer(args[[2]])
out_root <- if (length(args) >= 3L) args[[3]] else "tests/fixtures/differential"

# `%.17g` for the doubles, integers as integers: the same rule the portable format
# uses, so the Rust readers round-trip exactly.
fmt <- function(x) sprintf("%.17g", x)
# Every generated file is gzipped: 200 cases at the sizes the plan asks for are
# around 60 MB of text, and around 15 MB after zlib. `resolve_case_file` in
# reference_common.R finds them again.
# `apply(df, 1, paste, ...)` would go through `format()` and pad the numeric
# columns to a common width with leading spaces, so every field is converted
# explicitly instead.
write_tsv <- function(df, path, header = TRUE) {
    con <- gzfile(path, "wt")
    on.exit(close(con))
    if (header) {
        writeLines(paste(names(df), collapse = "\t"), con)
    }
    cols <- lapply(df, as.character)
    for (i in seq_len(nrow(df))) {
        fields <- vapply(cols, function(c) c[[i]], character(1))
        writeLines(paste(fields, collapse = "\t"), con)
    }
}

strands <- c("+", "-", "*")

draw_case <- function(case_id) {
    set.seed(seed + case_id)
    dir <- file.path(out_root, sprintf("case_%03d", case_id))
    unlink(dir, recursive = TRUE)
    dir.create(dir, recursive = TRUE)

    n_chr <- sample.int(3, 1)
    chr_names <- sprintf("chr%d", seq_len(n_chr))
    # A genome region short enough that the site and annotation sets overlap a lot,
    # but long enough that they do not always overlap.
    span <- sample(c(2000L, 20000L, 200000L), 1)

    # ------------------------------------------------------------------ sites
    n_sites <- sample(c(200L, 800L, 5000L), 1, replace = TRUE)
    site_width <- sample(1:2, n_sites, replace = TRUE)
    site_chr <- sample(chr_names, n_sites, replace = TRUE)
    site_start <- as.integer(sample.int(span, n_sites, replace = TRUE))
    site_end <- site_start + site_width - 1L
    # `*` is included on purpose: strand matching must not drop star sites.
    site_strand <- sample(strands, n_sites, replace = TRUE, prob = c(0.45, 0.45, 0.1))
    # Coverage values adjacent to zero, so the threshold boundary is exercised.
    site_cov <- sample(c(0L, 1L, 1L, 2L, 3L, sample.int(30L, 20L, TRUE)), n_sites, TRUE)
    # Scores are already-rounded six-decimal fractions, as read_methylome leaves them.
    site_score <- round(runif(n_sites, 0, 1), 6)
    sites <- data.frame(
        chr = site_chr, start = site_start, end = site_end,
        strand = site_strand, score = fmt(site_score), coverage = site_cov,
        stringsAsFactors = FALSE
    )
    sites <- sites[order(sites$chr, sites$start, sites$end), ]
    write_tsv(sites, file.path(dir, "msites.tsv.gz"))

    # ----------------------------------------------------------- GC windows
    n_bins <- sample(1:5, 1)
    win_width <- 30L
    if (runif(1) < 0.5) {
        # Tiled: abutting 30 bp windows, which is what makes a width-2 site land
        # in two bins and turns "mean over hits" into something other than "mean
        # over sites".
        n_windows <- max(10L, span %/% win_width)
        win_chr <- sample(chr_names, 1)
        win_start <- seq.int(1L, by = win_width, length.out = n_windows)
        win_end <- win_start + win_width - 1L
        # Some gaps are punched into the tiling so the "hits no window" path runs as
        # well as the "straddles two windows" path.
        gap <- runif(n_windows) < 0.3
        win_start[gap] <- win_start[gap] + win_width * sample.int(4L, sum(gap), TRUE)
        win_end[gap] <- win_start[gap] + win_width - 1L
        windows <- data.frame(
            chr = rep(win_chr, n_windows), start = as.integer(win_start),
            end = as.integer(win_end), strand = "*",
            gc_bin = sample.int(n_bins, n_windows, TRUE),
            stringsAsFactors = FALSE
        )
    } else {
        n_windows <- sample(c(20L, 200L, 2000L), 1)
        w_start <- as.integer(sample.int(span, n_windows, TRUE))
        w_width <- sample(c(30L, 30L, 30L, 45L, 60L), n_windows, TRUE)
        windows <- data.frame(
            chr = sample(chr_names, n_windows, TRUE), start = w_start,
            end = w_start + w_width - 1L, strand = sample(strands, n_windows, TRUE),
            gc_bin = sample.int(n_bins, n_windows, TRUE),
            stringsAsFactors = FALSE
        )
        windows <- windows[order(windows$chr, windows$start, windows$end), ]
    }
    write_tsv(windows, file.path(dir, "gc_windows.tsv.gz"))

    # ---------------------------------------------------------------- TFBS
    n_tfbs <- sample(c(50L, 500L, 2000L), 1)
    tfb_chr <- sample(chr_names, n_tfbs, TRUE)
    tfb_start <- as.integer(sample.int(span, n_tfbs, TRUE))
    # Mixed widths, odd and even, so W = width(first) + 130 lands on both parities
    # and the midpoint lands on a .5 that has to round half-to-even.
    tfb_width <- sample(c(11L, 12L, 21L, 22L, 101L, 102L, 411L), n_tfbs, TRUE)
    tfbs <- data.frame(
        chr = tfb_chr, start = tfb_start, end = tfb_start + tfb_width - 1L,
        strand = sample(strands, n_tfbs, TRUE, prob = c(0.45, 0.45, 0.1)),
        stringsAsFactors = FALSE
    )
    write_tsv(tfbs, file.path(dir, "tfbs.tsv.gz"))

    # ---------------------------------------------------------------- gcfreq
    # L odd and even in 101..601. Columns sum to 1, as they do in the real
    # annotation, so the expected profile is a distribution.
    # Odd and even values inside the plan's 101..601 range.
    L <- sample(c(101L, 102L, 151L, 152L, 301L, 302L, 512L, 601L), 1)
    gcfreq <- matrix(runif(n_bins * L, 0.05, 1), nrow = n_bins, ncol = L)
    gcfreq <- gcfreq / colSums(gcfreq)
    con <- gzfile(file.path(dir, "gcfreq.tsv.gz"), "wt")
    writeLines(apply(gcfreq, 1, paste, collapse = "\t"), con)
    close(con)

    # -------------------------------------------------------------- manifest
    write_tsv(
        data.frame(
            motif = "M1", tfbs_path = "tfbs.tsv.gz", gcfreq_path = "gcfreq.tsv.gz",
            stringsAsFactors = FALSE
        ),
        file.path(dir, "motifs.tsv.gz")
    )

    # ------------------------------------------------------------- enhancer
    if (runif(1) < 1 / 3) {
        n_enh <- sample(c(1L, 5L, 50L), 1)
        e_start <- as.integer(sample.int(span, n_enh, TRUE))
        e_end <- e_start + sample(c(100L, 5000L, span %/% 4L), n_enh, TRUE)
        enhancer <- data.frame(
            chr = sample(chr_names, n_enh, TRUE), start = e_start,
            end = e_end, strand = sample(strands, n_enh, TRUE),
            stringsAsFactors = FALSE
        )
        write_tsv(enhancer, file.path(dir, "enhancer.tsv.gz"))
    }

    # --------------------------------------------------------------- options
    ignore_strand <- runif(1) < 0.5
    writeLines(
        c(
            "ignore_strand",
            if (ignore_strand) "TRUE" else "FALSE"
        ),
        file.path(dir, "options.tsv")
    )

    # --------------------------------------------------------------- oracle
    out_dir <- file.path(dir, "expected")
    message <- tryCatch(
        {
            run_reference_case(dir, out_dir)
            NULL
        },
        error = function(e) conditionMessage(e)
    )
    if (!is.null(message)) {
        unlink(out_dir, recursive = TRUE)
        writeLines(message, file.path(dir, "error.txt"))
    }
    invisible(list(dir = dir, error = message))
}

if (!dir.exists(out_root) && !dir.create(out_root, recursive = TRUE)) {
    stop("could not create ", out_root)
}

n_errors <- 0L
for (i in seq_len(n)) {
    r <- draw_case(i)
    if (is.null(r$error)) {
        cat(sprintf("case_%03d ok\n", i))
    } else {
        n_errors <- n_errors + 1L
        cat(sprintf("case_%03d error: %s\n", i, r$error))
    }
}
cat(sprintf(
    "%d cases in %s (%d raise an error upstream)\n", n, out_root, n_errors
))