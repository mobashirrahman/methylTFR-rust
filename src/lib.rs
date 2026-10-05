//! methyltfr -- a Rust port of the Bioconductor package `methylTFR`.
//!
//! Quantification of DNA methylation signatures in transcription factor binding
//! sites. The reference implementation is `methylTFR` 0.99.9 at the commit
//! pinned in `docs/reference-version.md`; this port reproduces its behaviour
//! including its quirks, as specified in `docs/AGENT_PLAN.md` section 2, and
//! differs from it only where `docs/divergences.md` says so.
//!
//! The two rules that catch everyone out:
//!
//! * Coordinates are 1-based and closed, and are copied from the input verbatim.
//!   A BED-style record has width 2 here, exactly as it does upstream. Nothing
//!   converts BED to half-open or 0-based.
//! * The GC-bin mean is a mean over *overlap hits*, not over sites, so a site
//!   straddling two abutting windows is counted twice. See [`gc`].

#![forbid(unsafe_code)]
#![deny(clippy::all)]

pub mod cli;
pub mod deviation;
pub mod error;
pub mod expected;
pub mod gc;
pub mod intervals;
pub mod io;
pub mod model;
pub mod pipeline;
pub mod rmath;

pub use error::{Error, Result};
pub use gc::GcBins;
pub use model::{
    Annotation, ChromTable, GcFreq, GcWindow, Methylome, MotifAnnotation, Range, Site, Strand,
};
pub use pipeline::{CaseOutput, InputFormat, Options, Row, Sample, run};
