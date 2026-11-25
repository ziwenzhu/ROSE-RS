//! Domain error types for `rose-core`.

use thiserror::Error;

/// Errors returned by pure algorithms in `rose-core`.
#[must_use]
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Error)]
pub enum RoseCoreError {
    /// Start coordinate is greater than or equal to end coordinate in half-open `[start, end)`.
    #[error("invalid half-open genomic interval [{start}, {end}) on '{chrom}'")]
    InvalidInterval {
        /// Chromosome or contig identifier.
        chrom: String,
        /// 0-based start coordinate.
        start: u64,
        /// 0-based exclusive end coordinate.
        end: u64,
    },

    /// Duplicate constituent enhancer identifier encountered in input peaks.
    #[error("duplicate constituent region identifier '{id}'")]
    DuplicateRegionId {
        /// Duplicated peak identifier.
        id: String,
    },

    /// Input regions failed chromosome-name reconciliation against the alignment/track header.
    #[error(
        "only {matched} of {total} regions ({pct:.1}%) matched reference contigs (minimum required: {min_pct:.1}%)"
    )]
    ContigMismatch {
        /// Number of regions whose contig matched after `chr` reconciliation.
        matched: usize,
        /// Total regions checked.
        total: usize,
        /// Percentage matched (`0..=100`).
        pct: f64,
        /// Required percentage threshold (`0..=100`).
        min_pct: f64,
    },

    /// Control track is too shallow for reliable bin-based copy-number estimation.
    #[error(
        "control coverage too shallow for --cn-correct: autosomal median bin count is {median_reads:.1} (minimum required: {min_required:.1}); supply --cn-segments instead"
    )]
    InsufficientControlDepth {
        /// Observed median reads/fragments per mappable autosomal bin.
        median_reads: f64,
        /// Minimum required median reads per bin.
        min_required: f64,
    },

    /// Invalid algorithm hyperparameter or empty input where non-empty data is required.
    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
}
