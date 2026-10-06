#!/usr/bin/env Rscript

# footprint_check.R -- does the paper's own sample reproduce the paper's own
# published footprint numbers?
#
#   usage: footprint_check.R <rds_dir> <sample.bed.gz> <out_tsv> [motif...]
#
# The paper reports, for the aggregate methylation profile around motif
# occurrences (Figure 2E), a single number per motif per cell type. For T cells
# it gives CEBPB 0.01 and SPI1 -0.09, describing the former as "unchanged" and
# the latter as mildly depleted.
#
# Those summary values come from the authors' analysis code, which is not in the
# package, so the exact statistic they used cannot be read off. This script
# therefore computes a plainly defined one and says so:
#
#     depletion = mean(middle bins) / mean(first, last bins) - 1
#
# using the same five position intervals methylTFR's own dev_helper uses,
# WITHOUT the GC correction. That puts it on the same scale as the paper's
# numbers (a value below zero means the centre is less methylated than the
# flanks) without claiming to be the identical computation.
#
# Two further caveats, both stated in the output:
#   * The paper restricts all analyses to gene-distal regions. The published hg38
#     annotation has no distal flag, so these are whole-genome numbers.
#   * The paper's "T cells" is a BLUEPRINT population; this is one ENCODE
#     primary T cell. Same lineage, not the same sample.
#
# This is a directional check, not a reproduction. A large disagreement would
# mean the data or annotation is wrong; agreement in sign and rough size means
# they are right.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
    stop("usage: footprint_check.R <rds_dir> <sample.bed.gz> <out_tsv> [motif...]")
}
rds_dir <- args[[1]]
sample_path <- args[[2]]
out_tsv <- args[[3]]
motifs <- if (length(args) >= 4L) args[-(1:3)] else c("CEBPB", "SPI1")

gcfreqs <- readRDS(file.path(rds_dir, "jaspar2020_motif_gcfreq.rds"))
tfbs_all <- readRDS(file.path(rds_dir, "jaspar2020_tf_bindsites.rds"))
msites <- read_methylome(sample_path, "encode", cov_threshold = 5)
message("sites: ", length(msites))

# The paper's published T-cell values, for comparison only.
published <- c(CEBPB = 0.01, SPI1 = -0.09)

out <- list()
for (motif in motifs) {
    if (!motif %in% names(tfbs_all)) {
        message("motif not in annotation, skipping: ", motif)
        next
    }
    # computeFootprint is internal (not in NAMESPACE), so it is reached through
    # the pinned reference build in reference/methylTFR rather than being
    # reimplemented here. The point is to ask the authors' own code the
    # question, and a reimplementation would answer a different one.
    fp <- methylTFR:::computeFootprint(
        motif_name = motif,
        tf_bindsites = setNames(list(tfbs_all[[motif]]), motif),
        msites = msites,
        enhancer = NULL
    )
    # Same five intervals as dev_helper: the middle one against the outer two.
    # Sites further than 250 bp from the centre land outside every interval, so
    # the NA group is dropped before anything is compared. Note that x runs
    # -270..269 rather than -257..257, because mid_point is round()ed
    # half-to-even and lands a base off centre for even motif widths; the same
    # offset is present in methylTFR's own deviation arithmetic, so the
    # intervals are kept identical on purpose.
    fp[, grp := cut(x, c(-250, -200, -25, 25, 200, 250), labels = FALSE)]
    means <- fp[!is.na(grp), .(centre_methyl = mean(avg_methyl)), by = grp]
    means <- means[order(grp)]
    stopifnot(identical(means$grp, 1:5))
    # `(grp == 1) | (grp == 5)`, not `c(grp == 1, grp == 5)`: the latter builds
    # a ten-element index for a five-row table and indexes off the end.
    middle <- means$centre_methyl[means$grp == 3]
    flank <- mean(means$centre_methyl[(means$grp == 1) | (means$grp == 5)])
    depletion <- middle / flank - 1

    out[[length(out) + 1L]] <- data.table(
        motif = motif,
        n_tfbs = length(tfbs_all[[motif]]),
        centre = middle,
        flank = flank,
        depletion = depletion,
        paper_t_cell = if (motif %in% names(published)) published[[motif]] else NA_real_,
        x_min = min(fp$x),
        x_max = max(fp$x),
        covered_bp = sum(fp$n)
    )
    message(sprintf(
        "%-12s centre %.4f  flank %.4f  depletion %+.3f%s",
        motif, middle, flank, depletion,
        if (motif %in% names(published))
            sprintf("   (paper T cells %+.2f)", published[[motif]]) else ""
    ))
    rm(fp)
    invisible(gc())
}

res <- rbindlist(out)
fwrite(res, out_tsv, sep = "\t")
cat("\nwrote ", out_tsv, "\n", sep = "")
cat("whole-genome, no distal restriction, single ENCODE T cell; ",
    "compare direction and order of magnitude only\n", sep = "")
