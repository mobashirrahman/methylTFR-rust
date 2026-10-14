//! Error type for the whole crate (`docs/AGENT_PLAN.md` T11).
//!
//! Every error that can be attributed to an input file carries the path and,
//! when the file has line numbers, the 1-based line number. The strings that
//! upstream's `stop()` calls produce are reproduced verbatim where they are
//! reachable from the CLI, because a drop-in replacement should say the same
//! thing when it refuses the same input.

use std::path::{Path, PathBuf};

/// Result alias used everywhere in the crate.
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A malformed record in a file. Always carries path and line.
    #[error("{path}:{line}: {message}")]
    Parse {
        path: PathBuf,
        line: usize,
        message: String,
    },

    /// A file could not be opened or read.
    #[error("{path}: {message}")]
    Io { path: PathBuf, message: String },

    /// A whole-file condition: the file parsed, but it is not usable. Upstream
    /// reaches these through `stop()` inside `read_methylome`, which runs after
    /// the whole file has been read, so there is no line number to report.
    #[error("{path}: {message}")]
    File { path: PathBuf, message: String },

    /// A run-level condition from `AGENT_PLAN.md` section 2.3, 2.4 or 2.7. These
    /// are not tied to a file: the same motif fails for every sample.
    #[error("{message}")]
    Run { message: String },

    /// A path named on the command line does not exist.
    #[error("{0}")]
    Missing(String),
}

impl Error {
    pub fn parse(path: impl AsRef<Path>, line: usize, message: impl Into<String>) -> Self {
        Error::Parse {
            path: path.as_ref().to_path_buf(),
            line,
            message: message.into(),
        }
    }

    pub fn io(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Error::Io {
            path: path.as_ref().to_path_buf(),
            message: source.to_string(),
        }
    }

    pub fn file(path: impl AsRef<Path>, message: impl Into<String>) -> Self {
        Error::File {
            path: path.as_ref().to_path_buf(),
            message: message.into(),
        }
    }

    pub fn run(message: impl Into<String>) -> Self {
        Error::Run {
            message: message.into(),
        }
    }

    /// The path this error is about, when it has one.
    pub fn path(&self) -> Option<&Path> {
        match self {
            Error::Parse { path, .. } | Error::Io { path, .. } | Error::File { path, .. } => {
                Some(path.as_path())
            }
            Error::Run { .. } | Error::Missing(_) => None,
        }
    }

    /// The 1-based line number this error is about, when it has one.
    pub fn line(&self) -> Option<usize> {
        match self {
            Error::Parse { line, .. } => Some(*line),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_errors_report_path_and_line() {
        let e = Error::parse("a/b.tsv", 7, "not a number");
        assert_eq!(e.to_string(), "a/b.tsv:7: not a number");
        assert_eq!(e.line(), Some(7));
        assert_eq!(e.path().unwrap().to_str().unwrap(), "a/b.tsv");
    }

    #[test]
    fn run_errors_have_no_path() {
        let e = Error::run("no hits");
        assert_eq!(e.to_string(), "no hits");
        assert!(e.path().is_none());
        assert!(e.line().is_none());
    }
}
