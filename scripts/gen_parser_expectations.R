#!/usr/bin/env Rscript

# gen_parser_expectations.R -- T30..T35, Tier A parser oracle.
#
#   usage: gen_parser_expectations.R <parsers_dir>
#
# Runs the six examples bundled in `tests/fixtures/parsers/` through the *real*
# methylTFR 0.99.9 `read_methylome()` and writes one `<type>.expected.tsv` per
# format, plus the two coverage-threshold variants the plan names. This is what
# `tests/parser_*.rs` compares against, bit for bit on score and coverage --
# `docs/divergences.md` D2 is why a tolerance is not allowed there.
#
# Needs the Tier A environment (docs/reference-version.md); CI never runs this.
#
# Output columns, tab separated, one header line:
#
#   chr  start  end  strand  score  coverage
#
# with floats written by sprintf("%.17g", x) except integral values, which are
# written as integers. That is the same rule the portable format uses, so the two
# are byte-comparable.

suppressPackageStartupMessages({
    library(methylTFR)
    library(GenomicRanges)
})

# R's file name for each `type` argument, which is not always the type itself.
file_for_type <- c(
    epp = "epp.tsv.gz",
    bissnp = "bissnp.tsv.gz",
    allc = "allc.tsv.gz",
    bismarkcytosine = "bismarkCytosine.tsv.gz",
    bismarkcov = "bismarkCov.tsv.gz",
    encode = "encode.tsv.gz"
)

fmt_num <- function(x) {
    out <- character(length(x))
    integral <- is.finite(x) & x == trunc(x) & abs(x) < 1e15
    out[integral] <- sprintf("%.0f", x[integral])
    rest <- !integral
    out[rest] <- sprintf("%.17g", x[rest])
    na <- is.na(x) & !is.nan(x)
    out[na] <- "NA"
    out
}

dump_one <- function(path, type, out, thr = 1) {
    g <- read_methylome(path, type, cov_threshold = thr)
    d <- data.frame(
        chr = as.character(seqnames(g)),
        start = start(g),
        end = end(g),
        strand = as.character(strand(g)),
        score = g$score,
        coverage = g$coverage,
        stringsAsFactors = FALSE
    )
    con <- file(out, "wt")
    on.exit(close(con))
    writeLines(c(
        "chr\tstart\tend\tstrand\tscore\tcoverage",
        paste(
            d$chr, d$start, d$end, d$strand,
            fmt_num(d$score), fmt_num(d$coverage),
            sep = "\t"
        )
    ), con)
    cat(sprintf(
        "%-22s thr=%-4g %3d records -> %s\n",
        basename(path), thr, nrow(d), basename(out)
    ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
    stop("usage: gen_parser_expectations.R <parsers_dir>")
}
dir <- args[[1]]
if (!dir.exists(dir)) {
    stop(dir, " is not a directory")
}

for (ty in names(file_for_type)) {
    f <- file.path(dir, file_for_type[[ty]])
    if (!file.exists(f)) {
        stop("missing example: ", f)
    }
    dump_one(f, ty, file.path(dir, paste0(ty, ".expected.tsv")))
}

# The two variants AGENT_PLAN.md section 5 names explicitly: `encode` at
# cov_threshold = 20 keeps 4 of 6 records, and `bismarkCytosine` at
# cov_threshold = 0 still drops the 0/0 row, because that row is dropped by the
# NaN filter and not by the threshold.
dump_one(
    file.path(dir, file_for_type[["encode"]]), "encode",
    file.path(dir, "encode.expected.thr20.tsv"), thr = 20
)
dump_one(
    file.path(dir, file_for_type[["bismarkcytosine"]]), "bismarkcytosine",
    file.path(dir, "bismarkcytosine.expected.thr0.tsv"), thr = 0
)
dump_one(
    file.path(dir, file_for_type[["allc"]]), "allc",
    file.path(dir, "allc.expected.thr0.tsv"), thr = 0
)

# Refresh SHA256SUMS so the expected files are covered too. The six example
# files come first, in the order export_reference_data.R wrote them, then the
# expected files in sorted order.
sums_file <- file.path(dir, "SHA256SUMS")
examples <- unname(file_for_type)
expected <- sort(setdiff(
    list.files(dir, pattern = "\\.tsv$"),
    c("msites.tsv", "gc_windows.tsv", "tfbs.tsv", "gcfreq.tsv", "motifs.tsv",
      "enhancer.tsv")
))
tmp <- tempfile()
on.exit(unlink(tmp))
status <- system2(
    "sha256sum", shQuote(file.path(dir, c(examples, expected))),
    stdout = tmp
)
if (!identical(status, 0L)) {
    stop("sha256sum failed")
}
sums <- readLines(tmp, warn = FALSE)
sums <- sub(paste0(dir, "/"), "", sums, fixed = TRUE)
writeLines(sums, sums_file)
cat("SHA256SUMS:", length(sums), "files\n")