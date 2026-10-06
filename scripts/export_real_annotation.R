#!/usr/bin/env Rscript

# export_real_annotation.R -- turn the published methylTFRAnnotationHg38 objects
# into the portable annotation directory the Rust CLI reads.
#
#   usage: export_real_annotation.R <rds_dir> <out_dir> [motif_list] ["all"]
#
# The upstream package ships its hg38 annotation as three .rds objects
# (Zenodo record 22206980, CC-BY 4.0), which is the form the paper's own code
# consumes. The Rust CLI reads a directory of TSVs instead. This script is the
# bridge, and it exists so the real-data comparison is against the published
# annotation rather than a reconstruction of it.
#
# Everything is written with the precision needed to round-trip: %.17g for the
# GC frequency matrices, integers for coordinates. A GC frequency rounded to R's
# default 7 significant digits would still parse, but the deviations computed
# from it would differ in the last bits, which is the thing this whole exercise
# is measuring. So no lossy formatting anywhere.
#
# `motif_list` is a file of motif names, one per line, `#` for comments. The
# default, "all", exports all 632 motifs and needs about 14 GiB of headroom for
# the exported TFBS (263,425,993 sites). A subset is a few minutes and a few GiB.

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop("usage: export_real_annotation.R <rds_dir> <out_dir> [motif_list] [all]")
}
rds_dir <- args[[1]]
out_dir <- args[[2]]
motif_arg <- if (length(args) >= 3L) args[[3]] else "all"

phase <- function(label) {
    message(sprintf("[%7.1f s] %s", as.numeric(proc.time()[["elapsed"]]) - t0, label))
    invisible(NULL)
}
t0 <- as.numeric(proc.time()[["elapsed"]])

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tfbs"), showWarnings = FALSE)
dir.create(file.path(out_dir, "gcfreq"), showWarnings = FALSE)

gcfreqs <- readRDS(file.path(rds_dir, "jaspar2020_motif_gcfreq.rds"))
phase("read jaspar2020_motif_gcfreq.rds")

# The two objects are published separately but describe the same 632 motifs in
# the same order; check that rather than assuming it, because a silent mismatch
# would pair every motif with another motif's expected profile.
tfbs_all <- readRDS(file.path(rds_dir, "jaspar2020_tf_bindsites.rds"))
phase("read jaspar2020_tf_bindsites.rds")

stopifnot(length(gcfreqs) == length(tfbs_all))
stopifnot(setequal(names(gcfreqs), names(tfbs_all)))
if (!identical(unname(names(gcfreqs)), unname(names(tfbs_all)))) {
    stop("motif order differs between the two annotation objects")
}

if (motif_arg == "all") {
    wanted <- names(gcfreqs)
} else {
    wanted <- trimws(readLines(motif_arg, warn = FALSE))
    wanted <- wanted[nzchar(wanted) & !startsWith(wanted, "#")]
    unknown <- setdiff(wanted, names(gcfreqs))
    if (length(unknown)) {
        stop("unknown motifs requested: ", paste(unknown, collapse = ", "))
    }
    # Manifest order, not request order, so the output is reproducible.
    wanted <- names(gcfreqs)[names(gcfreqs) %in% wanted]
}
message("exporting ", length(wanted), " motifs")

# Motif names carry characters that are legal in a JASPAR name but not in a
# path ("MAX::MYC", "MZF1(var.2)"). The index prefix keeps the mapping unique
# whatever the name contains, and PROVENANCE.tsv records it.
safe <- function(x) {
    y <- gsub("[^A-Za-z0-9._-]", "_", x)
    substr(y, 1, 60)
}

# fwrite has no "%.17g" mode, and R's own as.character() gives 15 significant
# digits, so both would quietly lose the last two bits of every expected
# profile value. sprintf is the only formatter here that round-trips an f64.
write_gcfreq_precise <- function(m, path) {
    rows <- lapply(seq_len(nrow(m)), function(i) {
        paste(sprintf("%.17g", as.numeric(m[i, ])), collapse = "\t")
    })
    con <- gzfile(path, "wt")
    on.exit(close(con))
    writeLines(unlist(rows, use.names = FALSE), con)
    invisible(path)
}

# --- GC windows -------------------------------------------------------------
gc_dist <- readRDS(file.path(rds_dir, "genomewide_GC_hg38.rds"))
phase("read genomewide_GC_hg38.rds")

gc_df <- data.table(
    chr = as.character(seqnames(gc_dist)),
    start = as.integer(start(gc_dist)),
    end = as.integer(end(gc_dist)),
    strand = as.character(strand(gc_dist)),
    gc_bin = as.integer(mcols(gc_dist)$GC_bin)
)
# This data.table's fwrite has no gzip option, and the file is ~2.5 GB of text,
# so write it plain and compress it in one pass rather than buffering it.
gc_plain <- file.path(out_dir, "gc_windows.tsv")
fwrite(gc_df, gc_plain, sep = "\t", quote = FALSE)
gc_path <- paste0(gc_plain, ".gz")
status <- system2("gzip", c("-1", "-f", shQuote(gc_plain)))
if (!identical(as.integer(status), 0L)) {
    stop("gzip failed on ", gc_plain)
}
message("gc windows: ", nrow(gc_df), " -> ", gc_path)
gc_window_count <- nrow(gc_df)
rm(gc_df, gc_dist)
invisible(gc())
phase("write gc_windows.tsv.gz")

# --- per-motif TFBS and GC frequency ---------------------------------------
manifest <- vector("list", length(wanted))
prov <- vector("list", length(wanted))

for (i in seq_along(wanted)) {
    motif <- wanted[[i]]
    stem <- sprintf("%03d_%s", i, safe(motif))

    g <- tfbs_all[[motif]]
    tfbs_df <- data.table(
        chr = as.character(seqnames(g)),
        start = as.integer(start(g)),
        end = as.integer(end(g)),
        strand = as.character(strand(g))
    )
    tfbs_rel <- file.path("tfbs", paste0(stem, ".tsv.gz"))
    tfbs_plain <- file.path(out_dir, "tfbs", paste0(stem, ".tsv"))
    fwrite(tfbs_df, tfbs_plain, sep = "\t", quote = FALSE)
    if (!identical(as.integer(system2("gzip", c("-1", "-f", shQuote(tfbs_plain)))), 0L)) {
        stop("gzip failed on ", tfbs_plain)
    }

    # 5 GC bins x motif-width columns, no header.
    m <- gcfreqs[[motif]]
    gcf_rel <- file.path("gcfreq", paste0(stem, ".tsv"))
    write_gcfreq_precise(m, file.path(out_dir, gcf_rel))

    manifest[[i]] <- data.table(
        motif = motif,
        tfbs_path = tfbs_rel,
        gcfreq_path = gcf_rel
    )
    prov[[i]] <- data.table(
        motif = motif,
        n_tfbs = nrow(tfbs_df),
        gcfreq_rows = nrow(m),
        gcfreq_cols = ncol(m),
        tfbs_file = tfbs_rel,
        gcfreq_file = gcf_rel
    )
    message(sprintf("[%3d/%3d] %-22s sites %9d  gcfreq %dx%d",
                    i, length(wanted), motif, nrow(tfbs_df), nrow(m), ncol(m)))
    rm(tfbs_df, g)
    if (i %% 25L == 0L) invisible(gc())
}

fwrite(rbindlist(manifest), file.path(out_dir, "motifs.tsv"),
       sep = "\t", quote = FALSE)
fwrite(rbindlist(prov), file.path(out_dir, "PROVENANCE.tsv"),
       sep = "\t", quote = FALSE)
phase("write motifs.tsv and PROVENANCE.tsv")

cat(sprintf(
    "exported %d motifs, %d GC windows to %s\n",
    length(wanted), gc_window_count, out_dir
))
phase("done")
