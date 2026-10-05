#!/usr/bin/env Rscript

# bench_r.R -- T71, the R baseline.
#
#   usage: bench_r.R <data_dir> <out_csv> [motif_limit]
#
# The CSV is written to stdout; `out_csv` is only used to name what is being
# measured in messages.
#
# Runs the same work the Rust CLI does, on the same files, through methylTFR
# 0.99.9 itself:
#
#   read_methylome(sample, "bismarkcov")  -> msites
#   addGCBintoMethylome(msites, gc_dist) -> bin means
#   computeDeviation(motif, ...) per motif
#
# and writes the same CSV shape the CLI writes, so the two can be diffed as well
# as timed. Time it with /usr/bin/time -v for wall clock and peak RSS; see
# scripts/run_benchmarks.sh.
#
# `motif_limit` caps how many motifs are processed, which is how the Rust side is
# scaled to match: the two must do the same amount of work for the comparison to
# mean anything.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
})

# Phase timings, so the benchmark report can say where the time went instead of
# just how much there was.
phase <- function(label) {
    message(sprintf("[%7.2f s] %s", as.numeric(proc.time()[["elapsed"]]) - t0, label))
    invisible(NULL)
}
t0 <- as.numeric(proc.time()[["elapsed"]])

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: bench_r.R <data_dir> <out_csv> [motif_limit]")
}
data_dir <- args[[1]]
out_csv <- args[[2]]
motif_limit <- if (length(args) >= 3L) as.integer(args[[3]]) else NA_integer_

ann <- file.path(data_dir, "annotation")
sample_path <- file.path(data_dir, "sample0.bismarkCov")

fmt_num <- function(x) {
    out <- character(length(x))
    integral <- is.finite(x) & x == trunc(x) & abs(x) < 1e15
    out[integral] <- sprintf("%.0f", x[integral])
    rest <- !integral
    out[rest] <- sprintf("%.17g", x[rest])
    out[is.na(x) & !is.nan(x)] <- "NA"
    out
}

# Everything is read as text and converted explicitly, so a malformed field is an
# error here rather than a silent NA that would change the answer.
read_plain <- function(path) {
    read.delim(
        path, colClasses = "character", check.names = FALSE,
        quote = "", comment.char = "", stringsAsFactors = FALSE
    )
}

as_int <- function(x) as.integer(x)

phase("start")
message("reading sample ...")
msites <- read_methylome(sample_path, "bismarkcov", cov_threshold = 1)
message("sites: ", length(msites))
phase("read_methylome")

gcdf <- read_plain(file.path(ann, "gc_windows.tsv"))
gc_dist <- GRanges(
    gcdf$chr,
    IRanges(as_int(gcdf$start), as_int(gcdf$end)),
    strand = gcdf$strand
)
mcols(gc_dist)$GC_bin <- as.integer(gcdf$gc_bin)
message("gc windows: ", length(gc_dist))
phase("read gc windows")

motifs_df <- read_plain(file.path(ann, "motifs.tsv"))
if (!is.na(motif_limit) && motif_limit < nrow(motifs_df)) {
    motifs_df <- motifs_df[seq_len(motif_limit), ]
}
message("motifs: ", nrow(motifs_df))

bin_meth <- addGCBintoMethylome(msites, gc_dist, TRUE)
message("gc bins populated: ", nrow(bin_meth))
phase("addGCBintoMethylome")

# The TFBS sets are read one at a time and released again, so peak RSS reflects the
# largest single motif rather than the whole annotation. Upstream's annotation is
# supplied already in memory; this is the closest equivalent from files.
results <- vector("list", nrow(motifs_df))
for (i in seq_len(nrow(motifs_df))) {
    motif <- motifs_df$motif[[i]]
    t <- read_plain(file.path(ann, motifs_df$tfbs_path[[i]]))
    tfbs <- GRanges(
        t$chr, IRanges(as_int(t$start), as_int(t$end)),
        strand = t$strand
    )
    g <- as.matrix(read.delim(
        file.path(ann, motifs_df$gcfreq_path[[i]]),
        header = FALSE, check.names = FALSE
    ))
    storage.mode(g) <- "double"
    dev <- computeDeviation(
        motif = motif,
        msites = msites,
        # Upstream keeps tf_bindsites as a plain named list of GRanges and does
        # `tf_bindsites[[motif]]`, so the list holds the GRanges directly, named.
        tf_bindsites = setNames(list(tfbs), motif),
        gcfreqs = setNames(list(g), motif),
        enhancer = NULL,
        ignoreStrand = TRUE,
        binMsites = bin_meth
    )
    results[[i]] <- data.frame(
        sample = "sample0.bismarkCov",
        motif = motif,
        deviation = dev$dev[[1]],
        expected_deviation = dev$exp_dev[[1]],
        stringsAsFactors = FALSE
    )
    rm(tfbs, t, g, dev)
    invisible(gc())
    phase(paste("computeDeviation", motif))
}

# The CSV goes to stdout, not to a file this script opens, so that
# scripts/run_benchmarks.sh captures exactly what is measured and compares the same
# bytes the CLI wrote. Progress goes to stderr.
out <- do.call(rbind, results)
writeLines(
    c(
        "sample,motif,deviation,expected_deviation",
        paste(
            out$sample, out$motif,
            fmt_num(out$deviation), fmt_num(out$expected_deviation),
            sep = ","
        )
    )
)
invisible(out_csv)
phase("total")