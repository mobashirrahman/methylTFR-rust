#!/usr/bin/env Rscript

# export_reference_data.R -- T03 fixture exporter.
#
#   usage: export_reference_data.R <extdata_dir> <out_dir> [parsers_out_dir]
#
# Reads upstream methylTFR's `inst/extdata` and writes the portable fixture
# files described in AGENT_PLAN.md section 3.  It runs under **base R only**:
# the `.rda` files are loaded with `load()` and the `GRanges` / `IRanges` /
# `DFrame` slots are read with `attr()`, so no Bioconductor package is needed.
#
# Which object feeds which output file (see AGENT_PLAN.md section 4, "Fixture
# facts"):
#
#   msites.tsv     <- example_data.rda       `$msites`   (1000 sites, chr1)
#   gc_windows.tsv <- gcdist_subset.rda      `$gcdist`   (567 windows, 30 bp)
#   tfbs.tsv.gz    <- BATF_tf_bindsites.rda  `$tf_bindsites$BATF` (268717)
#   gcfreq.tsv     <- BATF_gcfreqs.rda       `$gcfreqs$BATF` (5 x 512)
#   motifs.tsv     <- written by hand, one row: BATF
#
# Run it twice -- once on the pinned checkout and once on the `1d99721`
# checkout -- to produce tests/fixtures/batf and tests/fixtures/batf_1d99721.

# ---------------------------------------------------------------- helpers ---

# The three S4 slots a `GRanges` carries are all `Rle` / `IRanges` / `DFrame`
# objects, and every one of them records its payload in attributes.  Reading
# them that way is what keeps this script base-R only.
rle_values <- function(r) {
    v <- attr(r, "values")
    if (is.factor(v)) {
        v <- as.character(v)
    }
    rep(v, attr(r, "lengths"))
}

# Coerce a `GRanges` to a plain data.frame of chr / start / end / strand plus
# every elementMetadata column.  `end` is reconstructed as start + width - 1,
# i.e. upstream's 1-based closed convention, with no off-by-one conversion.
granges_to_df <- function(o) {
    rg <- attr(o, "ranges")
    start <- attr(rg, "start")
    d <- data.frame(
        chr = rle_values(attr(o, "seqnames")),
        start = start,
        end = start + attr(rg, "width") - 1L,
        strand = rle_values(attr(o, "strand")),
        stringsAsFactors = FALSE
    )
    md <- attr(attr(o, "elementMetadata"), "listData")
    for (n in names(md)) {
        d[[n]] <- md[[n]]
    }
    d
}

# AGENT_PLAN.md section 3: floats are written by R with sprintf("%.17g", x) so
# that they round-trip exactly.  Integers go through as.character() and keep no
# decimal point.
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

# A connection that transparently gzips when the path ends in `.gz`.  R's
# `gzfile()` writes mtime 0, so repeated runs are byte-identical and
# SHA256SUMS stays valid.
open_out <- function(path) {
    if (grepl("\\.gz$", path)) {
        gzfile(path, "wt")
    } else {
        file(path, "wt")
    }
}

write_table <- function(lines, path) {
    con <- open_out(path)
    on.exit(close(con))
    writeLines(lines, con)
}

tab <- function(...) paste(..., sep = "\t")

# `SHA256SUMS` is what `sha256sum -c` consumes.  base R can only do md5, so
# shell out; failing loudly beats silently writing a file that cannot be
# verified.
write_sha256sums <- function(dir, names) {
    if (!nzchar(Sys.which("sha256sum"))) {
        stop("sha256sum not found on PATH; cannot write SHA256SUMS")
    }
    tmp <- tempfile()
    on.exit(unlink(tmp))
    status <- system2(
        "sha256sum", shQuote(file.path(dir, names)),
        stdout = tmp, stderr = FALSE
    )
    if (!identical(status, 0L) || !file.exists(tmp)) {
        stop("sha256sum failed in ", dir)
    }
    sums <- readLines(tmp, warn = FALSE)
    # Drop the directory prefix sha256sum echoes back, keeping the file names
    # relative so that `sha256sum -c SHA256SUMS` works from inside `dir`.
    sums <- sub(paste0(dir, "/"), "", sums, fixed = TRUE)
    writeLines(sums, file.path(dir, "SHA256SUMS"))
    invisible(sums)
}

load_one <- function(path, object) {
    e <- new.env()
    loaded <- load(path, envir = e)
    if (!object %in% loaded) {
        stop(basename(path), " does not contain an object named ", object)
    }
    get(object, envir = e)
}

# ------------------------------------------------------------------- main ---

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
    stop(
        "usage: export_reference_data.R <extdata_dir> <out_dir> ",
        "[parsers_out_dir]"
    )
}
extdata_dir <- args[[1]]
out_dir <- args[[2]]
parsers_out_dir <- if (length(args) >= 3L) args[[3]] else NULL

if (!dir.exists(extdata_dir)) {
    stop(extdata_dir, " is not a directory")
}
if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) {
    stop("could not create ", out_dir)
}

msites <- granges_to_df(load_one(
    file.path(extdata_dir, "example_data.rda"), "msites"
))
gcdist <- granges_to_df(load_one(
    file.path(extdata_dir, "gcdist_subset.rda"), "gcdist"
))
tf_bindsites <- load_one(
    file.path(extdata_dir, "BATF_tf_bindsites.rda"), "tf_bindsites"
)
gcfreqs <- load_one(file.path(extdata_dir, "BATF_gcfreqs.rda"), "gcfreqs")

if (!is.list(tf_bindsites)) {
    stop("tf_bindsites is not a list; AGENT_PLAN.md section 4 expects one")
}
motif <- "BATF"
tfbs <- granges_to_df(tf_bindsites[[motif]])
gcfreq <- gcfreqs[[motif]]

if (!is.matrix(gcfreq)) {
    stop("gcfreqs$", motif, " is not a matrix")
}

data_files <- c("msites.tsv", "gc_windows.tsv", "tfbs.tsv.gz", "gcfreq.tsv", "motifs.tsv")

# msites.tsv -- the methylome sample itself.
write_table(c(
    tab("chr", "start", "end", "strand", "score", "coverage"),
    tab(
        msites$chr, msites$start, msites$end, msites$strand,
        fmt_num(msites$score), msites$coverage
    )
), file.path(out_dir, "msites.tsv"))

# gc_windows.tsv -- the genome-wide GC distribution subset.
write_table(c(
    tab("chr", "start", "end", "strand", "gc_bin"),
    tab(gcdist$chr, gcdist$start, gcdist$end, gcdist$strand, gcdist$GC_bin)
), file.path(out_dir, "gc_windows.tsv"))

# tfbs.tsv.gz -- original widths, original row order, never resized here.
write_table(c(
    tab("chr", "start", "end", "strand"),
    tab(tfbs$chr, tfbs$start, tfbs$end, tfbs$strand)
), file.path(out_dir, "tfbs.tsv.gz"))

# gcfreq.tsv -- no header, one row per GC bin, L columns.
gcfreq_rows <- lapply(
    seq_len(nrow(gcfreq)),
    function(i) paste(fmt_num(gcfreq[i, ]), collapse = "\t")
)
write_table(unlist(gcfreq_rows), file.path(out_dir, "gcfreq.tsv"))

# motifs.tsv -- the annotation manifest.  Its row order *is* the motif order
# (AGENT_PLAN.md section 2.8).
write_table(c(
    tab("motif", "tfbs_path", "gcfreq_path"),
    tab(motif, "tfbs.tsv.gz", "gcfreq.tsv")
), file.path(out_dir, "motifs.tsv"))

write_sha256sums(out_dir, data_files)

# Optional: the six upstream parser examples, copied verbatim.
if (!is.null(parsers_out_dir)) {
    if (!dir.exists(parsers_out_dir) &&
        !dir.create(parsers_out_dir, recursive = TRUE)) {
        stop("could not create ", parsers_out_dir)
    }
    parser_files <- c(
        "allc.tsv.gz", "bismarkCov.tsv.gz", "bismarkCytosine.tsv.gz",
        "bissnp.tsv.gz", "encode.tsv.gz", "epp.tsv.gz"
    )
    missing <- parser_files[
        !file.exists(file.path(extdata_dir, parser_files))
    ]
    if (length(missing) > 0L) {
        stop("missing parser examples: ", paste(missing, collapse = ", "))
    }
    dest <- file.path(parsers_out_dir, parser_files)
    ok <- file.copy(
        file.path(extdata_dir, parser_files), parsers_out_dir,
        overwrite = TRUE
    )
    if (!all(ok)) {
        stop("could not copy parser examples to ", parsers_out_dir)
    }
    write_sha256sums(parsers_out_dir, parser_files)
}

cat(
    sprintf(
        "%s: %d sites, %d gc windows, %d tfbs, gcfreq %dx%d, motif %s\n",
        out_dir, nrow(msites), nrow(gcdist), nrow(tfbs),
        nrow(gcfreq), ncol(gcfreq), motif
    )
)