#!/usr/bin/env Rscript

# reference_common.R -- the Tier A oracle, as a function.
#
#   run_reference_case(fixture_dir, out_dir) -> invisible(TRUE), or stops with
#                                               upstream's own error message
#
# Computes the four reference output files from the portable fixture files using
# the *real* methylTFR 0.99.9 package (Bioconductor 3.18 / R 4.3, see
# docs/reference-version.md):
#
#   methylTFR::addGCBintoMethylome        (AGENT_PLAN.md section 2.3)
#   methylTFR::computeDeviation           (section 2.4)
#   methylTFR:::computeExpectations       (section 2.7)
#   methylTFR:::dev_helper                (section 2.8)
#
# It writes the same four files as reference_base.R, so the two oracles can be
# compared file by file (scripts/compare_oracles.R does that).
#
# The hit lists are re-ordered into the canonical order the plan defines -- sites
# in input order, then TFBS by ascending resized start -- before being written.
# `dev_helper` is order-insensitive, so this does not change any result; it only
# makes the two oracles directly diffable.
#
# Sourced by run_reference.R and by gen_differential_cases.R, so the BATF fixtures
# and the 200 differential cases are produced by exactly the same code path.

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
    # A `.gz` name selects gzip on the way out, matching what `resolve_case_file`
    # expects on the way back in.
    con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
    on.exit(close(con))
    writeLines(lines, con)
}

read_fixture <- function(path) {
    d <- read_portable(path)
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

# ------------------------------------------------------------ the case reader

# `options.tsv` carries `ignoreStrand`. It is accepted with or without a header
# line, because the differential cases are the only writer and a one-line file
# reads better without one.
# A case's files are gzipped (`msites.tsv.gz` and friends) so that 200 of them fit
# in a repository at the sizes the plan asks for. Every reader here sniffs the
# content, not the name, so resolution is by trying the plain name first.
resolve_case_file <- function(dir, name) {
    plain <- file.path(dir, name)
    if (file.exists(plain)) {
        return(plain)
    }
    gz <- paste0(plain, ".gz")
    if (file.exists(gz)) {
        return(gz)
    }
    stop("case is missing ", plain, " or ", gz)
}

# The enhancer is optional, and its file may be plain or gzipped, so its presence
# cannot be checked with a bare `file.exists` on the uncompressed name.
case_file_or_null <- function(dir, name) {
    plain <- file.path(dir, name)
    if (file.exists(plain)) {
        return(plain)
    }
    gz <- paste0(plain, ".gz")
    if (file.exists(gz)) {
        return(gz)
    }
    NULL
}

read_ignore_strand <- function(path) {
    if (!file.exists(path)) {
        return(TRUE)
    }
    con <- file(path, "rt")
    on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0L) {
        return(TRUE)
    }
    if (lines[[1]] == "ignore_strand") {
        lines <- lines[-1]
    }
    if (length(lines) == 0L) {
        return(TRUE)
    }
    as.logical(lines[[1]])
}

# ------------------------------------------------------------------ the oracle

# Run one fixture directory. Stops, with upstream's message, when upstream would
# stop; the caller decides whether to record that as `error.txt` or let it
# propagate.
run_reference_case <- function(fixture_dir, out_dir) {
    if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
        stop("could not create ", out_dir)
    }
    msites_df <- read_fixture(resolve_case_file(fixture_dir, "msites.tsv"))
    gcdist_df <- read_fixture(resolve_case_file(fixture_dir, "gc_windows.tsv"))
    motifs_df <- read_fixture(resolve_case_file(fixture_dir, "motifs.tsv"))
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

    ignore_strand <- read_ignore_strand(file.path(fixture_dir, "options.tsv"))

    # AGENT_PLAN.md section 2.6: the enhancer reduces the GC windows once, before
    # any bin mean is computed.
    enhancer <- NULL
    enhancer_path <- case_file_or_null(fixture_dir, "enhancer.tsv")
    if (!is.null(enhancer_path)) {
        enhancer_df <- read_fixture(enhancer_path)
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
    ), file.path(out_dir, "expected_bins.tsv.gz"))

    write_lines_to(c(
        "x\tvalue",
        paste(fmt_num(as.integer(obs_x)), fmt_num(obs_value), sep = "\t")
    ), file.path(out_dir, "observed_profile.tsv.gz"))

    write_lines_to(c(
        "x\tvalue",
        paste(
            fmt_num(as.integer(exp_profile$x)),
            fmt_num(exp_profile$avg_methyl), sep = "\t"
        )
    ), file.path(out_dir, "expected_profile.tsv.gz"))

    write_lines_to(c(
        "motif\tobs_d\texp_d\tdev",
        paste(
            motif, fmt_num(obs_d), fmt_num(exp_d), fmt_num(dev_value),
            sep = "\t"
        )
    ), file.path(out_dir, "expected_dev.tsv.gz"))

    cat(sprintf(
        "%s: motif %s obs_d %.7g exp_dev %.7g dev %.7g (%d bins, %d observed hits)\n",
        out_dir, motif, obs_d, exp_d, dev_value,
        nrow(bin_meth), length(obs_x)
    ))
    invisible(TRUE)
}
