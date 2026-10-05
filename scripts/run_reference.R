#!/usr/bin/env Rscript

# run_reference.R -- T04 **Tier A** reference oracle.
#
#   usage: run_reference.R <fixture_dir> <out_dir>
#
# Same four output files as reference_base.R, but computed by the *real*
# methylTFR 0.99.9 package (Bioconductor 3.18 / R 4.3, see
# docs/reference-version.md).  The portable fixture files are read back into
# GRanges and handed to:
#
#   methylTFR::addGCBintoMethylome        (AGENT_PLAN.md section 2.3)
#   methylTFR::computeDeviation           (section 2.4)
#   methylTFR:::computeExpectations       (section 2.7)
#   methylTFR:::dev_helper                (section 2.8)
#
# Because the output files must be byte-identical to the Tier B ones, the hit
# lists are re-ordered into the canonical order the plan defines -- sites in
# input order, then TFBS by ascending resized start -- before being written.
# `dev_helper` is order-insensitive, so this does not change any result; it
# only makes the two oracles comparable with `diff`.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
    library(data.table)
})

# ------------------------------------------------------------------- helpers

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
    con <- file(path, "wt")
    on.exit(close(con))
    writeLines(lines, con)
}

read_fixture <- function(dir, name) {
    d <- read_portable(file.path(dir, name))
    for (col in c("start", "end")) {
        if (col %in% names(d)) {
            d[[col]] <- as.integer(d[[col]])
        }
    }
    d
}

to_granges <- function(d, extra = NULL) {
    r <- GRanges(
        seqnames = d$chr,
        ranges = IRanges(start = d$start, end = d$end),
        strand = d$strand
    )
    if (!is.null(extra)) {
        for (n in names(extra)) {
            mcols(r)[[n]] <- extra[[n]]
        }
    }
    r
}

# ------------------------------------------------------------------- the run

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: run_reference.R <fixture_dir> <out_dir>")
}
fixture_dir <- args[[1]]
out_dir <- args[[2]]

if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("could not create ", out_dir)
}

msites_df <- read_fixture(fixture_dir, "msites.tsv")
gcdist_df <- read_fixture(fixture_dir, "gc_windows.tsv")
motifs_df <- read_fixture(fixture_dir, "motifs.tsv")
stopifnot(nrow(motifs_df) >= 1L)
motif <- motifs_df$motif[[1]]

manifest_path <- function(p) {
    if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(fixture_dir, p)
}
tfbs_df <- read_portable(manifest_path(motifs_df$tfbs_path[[1]]))
tfbs_df$start <- as.integer(tfbs_df$start)
tfbs_df$end <- as.integer(tfbs_df$end)
gcfreq <- as.matrix(
    read_portable(manifest_path(motifs_df$gcfreq_path[[1]]), header = FALSE)
)
storage.mode(gcfreq) <- "double"

msites <- to_granges(msites_df, list(
    score = as.numeric(msites_df$score),
    coverage = as.integer(msites_df$coverage)
))
gcdist <- to_granges(gcdist_df, list(GC_bin = as.integer(gcdist_df$gc_bin)))
# Upstream keeps tf_bindsites as a plain list of GRanges, not a GRangesList.
tf_bindsites <- list()
tf_bindsites[[motif]] <- to_granges(tfbs_df)
gcfreqs <- list()
gcfreqs[[motif]] <- gcfreq

ignore_strand <- TRUE
opt_path <- file.path(fixture_dir, "options.tsv")
if (file.exists(opt_path)) {
    opts <- read_portable(opt_path)
    if ("ignore_strand" %in% names(opts)) {
        ignore_strand <- as.logical(opts$ignore_strand[[1]])
    }
}

# AGENT_PLAN.md section 2.6: the enhancer reduces the GC windows once, before
# any bin mean is computed.
enhancer <- NULL
enhancer_path <- file.path(fixture_dir, "enhancer.tsv")
if (file.exists(enhancer_path)) {
    enhancer_df <- read_fixture(fixture_dir, "enhancer.tsv")
    enhancer <- to_granges(enhancer_df)
    gcdist <- subsetByOverlaps(
        gcdist, enhancer,
        ignore.strand = ignore_strand
    )
}

bin_meth <- addGCBintoMethylome(msites, gcdist, ignore_strand)

# n_hits is not returned by addGCBintoMethylome, so the same findOverlaps call
# it makes internally is repeated here.
gc_hits <- findOverlaps(msites, gcdist, ignore.strand = ignore_strand)
n_hits <- table(factor(
    gcdist[gc_hits@to]$GC_bin,
    levels = as.character(bin_meth[, 1])
))

devs <- computeDeviation(
    motif = motif,
    msites = msites,
    tf_bindsites = tf_bindsites,
    gcfreqs = gcfreqs,
    enhancer = enhancer,
    ignoreStrand = ignore_strand,
    binMsites = bin_meth
)
# The observed profile, taken from the same internals computeDeviation uses.
tfbs_r <- resize(
    tf_bindsites[[motif]],
    width(tf_bindsites[[motif]])[1] + 130,
    fix = "center"
)
if (!is.null(enhancer)) {
    tfbs_r <- subsetByOverlaps(
        tfbs_r, enhancer,
        ignore.strand = ignore_strand
    )
}
mid <- round(end(tfbs_r) + ((start(tfbs_r) - end(tfbs_r)) / 2))
tf_hits <- findOverlaps(
    msites, tfbs_r,
    type = "within", ignore.strand = ignore_strand
)
# Canonical hit order: sites in input order, TFBS by ascending resized start.
ord <- order(tf_hits@from, start(tfbs_r)[tf_hits@to])
obs_x <- start(msites[tf_hits@from[ord]]) - mid[tf_hits@to[ord]]
obs_value <- msites[tf_hits@from[ord]]$score

# computeDeviation returns data.table(dev, exp_dev): `dev` is the
# bias-corrected deviation and `exp_dev` the expected deviation.  The observed
# deviation is not returned, so take it from dev_helper on the observed
# profile and assert the two agree.
dev_value <- as.numeric(devs$dev[[1]])
exp_d <- as.numeric(devs$exp_dev[[1]])
obs_d <- as.numeric(methylTFR:::dev_helper(data.table::data.table(
    x = as.numeric(obs_x),
    avg_methyl = as.numeric(obs_value)
)))
residual <- dev_value - (obs_d - exp_d)
if (!(abs(residual) <= 1e-12)) {
    stop("deviation is not obs_d - exp_d: ", dev_value, " vs ",
         obs_d - exp_d)
}

exp_profile <- methylTFR:::computeExpectations(bin_meth, gcfreq)

write_lines_to(c(
    "gc_bin\tmean\tn_hits",
    paste(
        fmt_num(bin_meth[, 1]), fmt_num(bin_meth[, 2]),
        fmt_num(as.integer(n_hits)), sep = "\t"
    )
), file.path(out_dir, "expected_bins.tsv"))

write_lines_to(c(
    "x\tvalue",
    paste(fmt_num(as.integer(obs_x)), fmt_num(obs_value), sep = "\t")
), file.path(out_dir, "observed_profile.tsv"))

write_lines_to(c(
    "x\tvalue",
    paste(
        fmt_num(as.integer(exp_profile$x)),
        fmt_num(exp_profile$avg_methyl), sep = "\t"
    )
), file.path(out_dir, "expected_profile.tsv"))

write_lines_to(c(
    "motif\tobs_d\texp_d\tdev",
    paste(
        motif, fmt_num(obs_d), fmt_num(exp_d), fmt_num(dev_value),
        sep = "\t"
    )
), file.path(out_dir, "expected_dev.tsv"))

cat(sprintf(
    "%s: motif %s obs_d %.7g exp_dev %.7g dev %.7g (%d bins, %d observed hits)\n",
    out_dir, motif, obs_d, exp_d, dev_value,
    nrow(bin_meth), length(obs_x)
))