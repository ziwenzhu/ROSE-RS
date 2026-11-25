//! Error types for genomic file parsing, binary track access, and atomic writes.

use std::io;
use std::path::PathBuf;

use rose_core::RoseCoreError;
use thiserror::Error;

/// Unified error returned by `rose-io` operations.
#[must_use]
#[non_exhaustive]
#[derive(Debug, Error)]
pub enum RoseIoError {
    /// Underlying filesystem or stream I/O error.
    #[error("I/O error on '{path}': {source}")]
    Io {
        /// File path associated with the error.
        path: PathBuf,
        /// Underlying OS I/O error.
        #[source]
        source: io::Error,
    },

    /// Domain validation error from `rose-core`.
    #[error(transparent)]
    Core(#[from] RoseCoreError),

    /// Malformed tabular or text line in BED, GFF, GTF, refGene, `.seg`, or TSV input.
    #[error("parse error in '{path}' at line {line}: {message}")]
    Parse {
        /// Source file path.
        path: PathBuf,
        /// 1-based line number.
        line: usize,
        /// Diagnostic message.
        message: String,
    },

    /// Missing `.bai`, `.csi`, or `.crai` alignment index required for indexed access.
    #[error(
        "missing alignment index for '{path}'; sort and index the file (e.g. `{path}.bai` or `{path}.csi`)"
    )]
    MissingIndex {
        /// BAM or CRAM file path.
        path: PathBuf,
    },

    /// Malformed binary BAM, CRAM, or bigWig stream.
    #[error("invalid binary format in '{path}': {message}")]
    BinaryFormat {
        /// Source file path.
        path: PathBuf,
        /// Diagnostic message.
        message: String,
    },

    /// JSON serialization error when writing `run.json`.
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
}

impl RoseIoError {
    /// Wraps an [`io::Error`] with file path context.
    pub fn io(path: impl Into<PathBuf>, source: io::Error) -> Self {
        Self::Io {
            path: path.into(),
            source,
        }
    }
}
