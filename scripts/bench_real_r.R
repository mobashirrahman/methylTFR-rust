#!/usr/bin/env Rscript

# bench_real_r.R -- the R side of the real-data comparison.
#
#   usage: bench_real_r.R <rds_dir> <sample.bed.gz> <out_csv> [motif_list|"all"]
#
# This runs methylTFR itself, from the published hg38 annotation objects, on a
# published ENCODE methylome. It is the canonical upstream path: the annotation
# arrives in memory as .rds, exactly as it does for the authors, rather than
# being re-read from the TSV directory the Rust CLI uses. The exported
# directory is verified to round-trip to the same doubles
# (scripts/verify_real_annotation.R), so the two sides are solving the same
# problem; only the input decoding differs, and the phase timings below report
# that separately instead of hiding it inside the total.
#
# Per-motif wall time goes to <out_csv>.timings.tsv because the interesting
# quantity on real data is not the single total but how cost scales with the
# number of TFBS in a motif: the motifs here span 23,061 to 5,192,071 sites, so
# R's per-motif cost is not a constant and a single benchmark number would hide
# that.
#
# `motif_list` caps the work; the default is all 632 motifs, which on the
# ENCODE T-cell sample is a multi-hour single-threaded run. Progress and
# timings go to stderr, the CSV goes to stdout.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
})

phase <- function(label) {
    elapsed <- as.numeric(proc.time()[["elapsed"]]) - t0
    message(sprintf("[%8.1f s] %s", elapsed, label))
    invisible(elapsed)
}
t0 <- as.numeric(proc.time()[["elapsed"]])

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("usage: bench_real_r.R <rds_dir> <sample.bed.gz> <out_csv> [motif_list]")
}
rds_dir <- args[[1]]
sample_path <- args[[2]]
out_csv <- args[[3]]
motif_arg <- if (length(args) >= 4L) args[[4]] else "all"

COV <- 5  # the paper's methods: "CpGs with coverage of at least 5"

fmt_num <- function(x) {
    out <- character(length(x))
    integral <- is.finite(x) & x == trunc(x) & abs(x) < 1e15
    out[integral] <- sprintf("%.0f", x[integral])
    rest <- !integral
    out[rest] <- sprintf("%.17g", x[rest])
    out[is.na(x) & !is.nan(x)] <- "NA"
    out
}

# --- annotation ------------------------------------------------------------
gcfreqs <- readRDS(file.path(rds_dir, "jaspar2020_motif_gcfreq.rds"))
phase("read gcfreqs")
tfbs_all <- readRDS(file.path(rds_dir, "jaspar2020_tf_bindsites.rds"))
phase("read tf_bindsites")
gc_dist <- readRDS(file.path(rds_dir, "genomewide_GC_hg38.rds"))
phase("read genomewide GC")
message("gc windows: ", length(gc_dist))
message("tfbs total: ", sum(lengths(tfbs_all)))

if (motif_arg == "all") {
    wanted <- names(gcfreqs)
} else {
    wanted <- trimws(readLines(motif_arg, warn = FALSE))
    wanted <- wanted[nzchar(wanted) & !startsWith(wanted, "#")]
    unknown <- setdiff(wanted, names(gcfreqs))
    if (length(unknown)) {
        stop("unknown motifs requested: ", paste(unknown, collapse = ", "))
    }
    wanted <- names(gcfreqs)[names(gcfreqs) %in% wanted]
}
message("motifs to process: ", length(wanted))

# --- methylome -------------------------------------------------------------
msites <- read_methylome(sample_path, "encode", cov_threshold = COV)
message("sites: ", length(msites))
phase("read_methylome")

bin_meth <- addGCBintoMethylome(msites, gc_dist, TRUE)
message("gc bins populated: ", nrow(bin_meth))
phase("addGCBintoMethylome")

# --- per motif -------------------------------------------------------------
results <- vector("list", length(wanted))
timings <- vector("list", length(wanted))
last <- as.numeric(proc.time()[["elapsed"]]) - t0

for (i in seq_along(wanted)) {
    motif <- wanted[[i]]
    tfbs <- tfbs_all[[motif]]
    dev <- computeDeviation(
        motif = motif,
        msites = msites,
        # Upstream indexes tf_bindsites[[motif]], so the container is a named
        # list of GRanges. A one-element list keeps that contract identical to
        # the full-annotation case while bounding peak memory to one motif.
        tf_bindsites = setNames(list(tfbs), motif),
        gcfreqs = setNames(list(gcfreqs[[motif]]), motif),
        enhancer = NULL,
        ignoreStrand = TRUE,
        binMsites = bin_meth
    )
    results[[i]] <- data.frame(
        sample = basename(sample_path),
        motif = motif,
        deviation = dev$dev[[1]],
        expected_deviation = dev$exp_dev[[1]],
        stringsAsFactors = FALSE
    )
    # `seconds` is this motif alone, `elapsed_s` is the clock since the start of
    # the run. Keeping both matters on real data, where per-motif cost tracks
    # motif size over a 200-fold range rather than being constant.
    now <- as.numeric(proc.time()[["elapsed"]]) - t0
    timings[[i]] <- data.frame(
        motif = motif,
        n_tfbs = length(tfbs),
        seconds = now - last,
        elapsed_s = now,
        stringsAsFactors = FALSE
    )
    last <- now
    rm(tfbs, dev)
    invisible(gc())
    message(sprintf("[%3d/%3d] %-24s sites %9d  dev % .12f  exp % .6f",
                    i, length(wanted), motif, timings[[i]]$n_tfbs,
                    results[[i]]$deviation, results[[i]]$expected_deviation))
}

out <- do.call(rbind, results)
writeLines(
    c(
        "sample,motif,deviation,expected_deviation",
        paste(out$sample, out$motif,
              fmt_num(out$deviation), fmt_num(out$expected_deviation),
              sep = ",")
    )
)
# The loop's own wall time, excluding annotation load and methylome parsing.
tim <- do.call(rbind, timings)
write.table(
    tim, paste0(out_csv, ".timings.tsv"), sep = "\t",
    row.names = FALSE, quote = FALSE
)
phase("total")
