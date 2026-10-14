#!/usr/bin/env Rscript

# verify_real_annotation.R -- prove the exported TSVs are the published objects
# again, not an approximation of them.
#
#   usage: verify_real_annotation.R <rds_dir> <ann_dir>
#
# This is the step that makes the real-data comparison meaningful. The Rust CLI
# reads the exported directory and the R reference reads the .rds files, so if
# the export lost a bit anywhere the two sides would be solving slightly
# different problems and any agreement afterwards would be luck.
#
# The check is on the parsed f64 values, not on the text: `%.17g` should give a
# bit-identical double, and comparing doubles is what the callers care about.

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: verify_real_annotation.R <rds_dir> <ann_dir>")
}
rds_dir <- args[[1]]
ann_dir <- args[[2]]

fail <- 0L
check <- function(ok, what) {
    if (ok) {
        cat(sprintf("  ok    %s\n", what))
    } else {
        cat(sprintf("  FAIL  %s\n", what))
        fail <<- fail + 1L
    }
}

gcfreqs <- readRDS(file.path(rds_dir, "jaspar2020_motif_gcfreq.rds"))
gc_src <- readRDS(file.path(rds_dir, "genomewide_GC_hg38.rds"))
manifest <- fread(file.path(ann_dir, "motifs.tsv"), header = TRUE)
cat("manifest rows: ", nrow(manifest), "\n\n")

cat("motif names and order\n")
# The export may be a subset, but it must be a subset in published order: a
# reordered manifest would silently reorder the output rows.
want_names <- names(gcfreqs)[names(gcfreqs) %in% manifest$motif]
check(identical(manifest$motif, want_names),
      sprintf("motifs.tsv is %d motifs in published order",
              length(want_names)))
check(all(manifest$motif %in% names(gcfreqs)), "every motif is a real motif")

cat("\nGC frequency matrices (parsed f64, compared bitwise)\n")
for (i in seq_len(nrow(manifest))) {
    m <- gcfreqs[[manifest$motif[[i]]]]
    got <- as.matrix(fread(
        file.path(ann_dir, manifest$gcfreq_path[[i]]),
        header = FALSE, sep = "\t"
    ))
    storage.mode(got) <- "double"
    same_dim <- identical(dim(got), dim(m))
    check(same_dim, sprintf("%s: dimensions %s", manifest$motif[[i]],
                            paste(dim(m), collapse = "x")))
    if (same_dim) {
        # 17 significant digits is the shortest width that names every double
        # uniquely, so agreement there is bitwise agreement.
        want <- sprintf("%.17g", as.vector(m))
        have <- sprintf("%.17g", as.vector(got))
        bad <- which(want != have)
        check(length(bad) == 0L,
              sprintf("%s: all %d values round-trip exactly",
                      manifest$motif[[i]], length(want)))
    }
}

cat("\nGC windows\n")
gc_got <- fread(file.path(ann_dir, "gc_windows.tsv.gz"), header = TRUE, sep = "\t")
check(nrow(gc_got) == length(gc_src), sprintf("row count %d", length(gc_src)))
check(identical(gc_got$chr, as.character(seqnames(gc_src))), "chr column")
check(identical(gc_got$start, as.integer(start(gc_src))), "start column")
check(identical(gc_got$end, as.integer(end(gc_src))), "end column")
check(identical(gc_got$strand, as.character(strand(gc_src))), "strand column")
check(identical(gc_got$gc_bin, as.integer(mcols(gc_src)$GC_bin)), "gc_bin column")

cat("\nTFBS (spot-checked: first, middle and last 1000 rows of each motif)\n")
tfbs_all <- readRDS(file.path(rds_dir, "jaspar2020_tf_bindsites.rds"))
for (i in seq_len(nrow(manifest))) {
    motif <- manifest$motif[[i]]
    g <- tfbs_all[[motif]]
    got <- fread(file.path(ann_dir, manifest$tfbs_path[[i]]),
                 header = TRUE, sep = "\t")
    idx <- unique(c(1:1000, seq(nrow(got) %/% 2, nrow(got) %/% 2 + 999),
                     (nrow(got) - 999):nrow(got)))
    idx <- idx[idx >= 1 & idx <= nrow(got)]
    ok <- nrow(got) == length(g) &&
        all(got$chr[idx] == as.character(seqnames(g))[idx]) &&
        all(got$start[idx] == as.integer(start(g))[idx]) &&
        all(got$end[idx] == as.integer(end(g))[idx]) &&
        all(got$strand[idx] == as.character(strand(g))[idx])
    check(ok, sprintf("%s: %d sites, %d rows spot-checked",
                      motif, length(g), length(idx)))
}

cat(sprintf("\n%d failure(s)\n", fail))
quit(status = if (fail > 0L) 1L else 0L)
