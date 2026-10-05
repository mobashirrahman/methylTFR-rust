//! methyltfr — Rust port of the Bioconductor package `methylTFR`.
//!
//! Quantification of DNA methylation signatures in transcription factor binding
//! sites (TFBS). The reference implementation is `methylTFR` 0.99.9 at the
//! commit pinned in `docs/reference-version.md`; the port is behaviour
//! compatible, including its quirks, as specified in `docs/AGENT_PLAN.md`.

#![forbid(unsafe_code)]
