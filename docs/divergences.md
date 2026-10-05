# Divergences from upstream methylTFR

Places where methylTFR-rs deliberately does not behave like
`methylTFR` 0.99.9 at commit `8aeab03cb469a9c213b7cd7fd2cbe6eeb8db5856`
(see `reference-version.md`). Every other difference in behaviour is a bug in
this port and must be fixed, not documented here.

Each entry is anchored to the upstream line it comes from, so a reviewer can
check the claim against the pinned clone. New entries must be added here before
the corresponding code lands, and `T60` closes this file.

---

## D1 — `encode` header detection

**Plan:** 2.2, table row `encode` ("header: auto-detected").

**Upstream** (`R/data_reader.R` 153–156):

```r
msites <- data.table::fread(filename, header = "auto", showProgress = FALSE)
```

`header = "auto"` makes `fread` decide by inspecting the data: it compares the
types it infers for the first line against the rest of the file and treats the
first line as a header when it cannot be reconciled. That inference uses the
whole column, so it is not reproducible from a single field.

**methylTFR-rs:** the first line is treated as a header iff its 2nd or 3rd
field does not parse as an integer.

**Why:** the Rust reader is single-pass and column-addressed, and the rule above
reproduces the intended behaviour for bedMethyl files (whose `chromStart` /
`chromEnd` are integers, and whose first line is a `#`-prefixed or textual
header) without buffering the file.

**Residual risk:** a headerless bedMethyl file whose first record has
non-integer `chromStart` or `chromEnd` is read as a header row and its data is
dropped. Upstream's inference would keep it. Accepted for v0.1; revisit if a
differential case (T50) exhibits it.

---

## D2 — floating-point parsing of input fields

**Plan:** 2.2, last bullet.

**Upstream** reads every file with `data.table::fread` (`R/data_reader.R` 84,
95, 107, 120, 135, 153). `fread` has its own C-level numeric parser, which is
not specified to produce the same `double` as a correctly-rounded
`strtod`.

**methylTFR-rs** parses numbers with `str::parse::<f64>()` / `str::parse::<i64>()`.

**Why:** `fread`'s parser can differ from correctly-rounded parsing by up to
1 ulp *before* the `round(·, 6)` of plan section 2.1 is applied, and a 1-ulp
difference can in principle flip the 6th decimal at a rounding boundary.

**Consequence:** parser scores and coverage values are compared **bit-exactly**
against R (`T30`–`T36`), not with the 1e-10 tolerance of plan section 3. A
mismatch here is a real divergence to report, not to paper over by loosening
the comparison.

**Status:** not chased. Recorded so that a future failure is understood.

---

## D3 — `--keep-going` instead of aborting the run

**Plan:** 2.8, run-level behaviour; CLI flag in T40.

**Upstream** aborts the entire run with `stop()` when a sample has no GC hits
(`R/expected_deviations.R` 54–56), a motif has no hits inside its binding sites
(`R/compute_deviations.R` 112–117), or the GC-bin count does not match the
`gcfreq` row count (`R/expected_deviations.R` 83–88).

**methylTFR-rs:** the default behaviour is identical — same conditions, same
messages, non-zero exit. `--keep-going` is an **extension**: that one cell is
written as `NA` and the reason is logged to stderr, and the run continues.

**Why:** useful when one motif in a large annotation is empty, which upstream
treats as fatal.

**Scope:** off by default, so every result produced without it is a faithful
port. Results produced *with* it are not comparable with upstream output for
that cell.

---

## Not divergences

Recorded because they look like they might be, but are in fact faithful:

- **Upstream has no `--keep-going` equivalent, and no error for a missing
  populated GC bin either** — the `computeExpectations` guard only checks that
  `binMsites` and `gcfreq` are matrices; the row-count mismatch surfaces as
  R's "non-conformable arguments" from `%*%`. Plan section 2.7 reproduces the
  check as an explicit error.
- **Motif order comes from `names(gcfreqs)`** (`R/methyltfr_core.R` 139), not
  from the TFBS list. In the portable format both come from one manifest row,
  so plan 2.8's "order of the annotation manifest" matches.
- **Coordinates are not converted.** Upstream copies BED-style `start`/`end`
  straight into `IRanges` (`R/data_reader.R` 189), so a BED record has width 2
  in `GRanges` terms and `allc` / `bismarkCytosine` records have width 1. This
  port keeps that verbatim, per plan rule 6.
- **`allc` is not filtered by its context column** (`R/data_reader.R` 106–116),
  even though the file's 4th column is the methylation context. Faithful.
