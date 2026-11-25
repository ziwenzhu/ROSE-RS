//! Pure, I/O-free algorithms for Rank Ordering of Super-Enhancers (`rose-rs`).
//!
//! `rose-core` normalises all genomic spans to 0-based half-open [`Interval`] values and provides:
//! - [`tss`]: $O(\log N)$ sorted transcript and TSS indexing with primary-contig preference.
//! - [`stitch`]: Blacklist removal, TSS `contained`/`subtract` handling, 12.5 kb stitching,
//!   `--max-tss-genes` disaggregation, and multi-sample consensus peak merging.
//! - [`signal`]: Exact overlap and reusable difference-array coverage, multi-tier `xxh3`
//!   downsampling, width-weighted bigWig integration, and contig-name reconciliation.
//! - [`cutoff`]: Exact $O(N)$ chord-tangent cutoff and R Brent `optimize` compatibility solver.
//! - [`genes`]: Binary-search `Enhancer -> Gene` and `Gene -> Enhancer` assignment.
//! - [`stability`]: Seeded SplitMix64 region bootstrap and read-depth saturation curves.
//! - [`cn`]: Copy-number estimation (`100 kb` bins) and external `.seg` segment normalisation.
//! - [`multi`]: Multi-sample signal matrices, DESeq2 median-of-ratios scaling, and group consensus.

pub mod cn;
pub mod cutoff;
pub mod error;
pub mod genes;
pub mod multi;
pub mod signal;
pub mod stability;
pub mod stitch;
pub mod stitch_sweep;
pub mod tss;
pub mod types;

pub use cn::{CopyNumberConfig, CopyNumberProfile, CopyNumberSegment};
pub use cutoff::{
    call_super_enhancers, compute_cutoff_on_slice, CallConfig, CallResult, CutoffSummary,
    RankedEnhancer, RegionSignal,
};
pub use error::RoseCoreError;
pub use genes::{map_enhancers_to_genes, EnhancerGeneRow, GeneEnhancerRow, GeneMapConfig};
pub use multi::{
    compute_median_ratio_size_factors, evaluate_multi_cohort, DiffLabel, MultiCohortResult,
    MultiRegionSummary, SampleSpec,
};
pub use signal::{
    integrate_bigwig_values, normalize_coverage, reconcile_chrom_names, ReadSpan, SignalScratch,
};
pub use stability::{
    evaluate_stability, DepthSaturationPoint, RegionStability, StabilityReport,
    DEFAULT_DOWNSAMPLE_FRACTIONS,
};
pub use stitch::{
    merge_consensus_peaks, merge_intervals, stitch_constituents, StitchConfig, StitchStats,
};
pub use stitch_sweep::{
    estimate_auto_stitch_distance, AutoStitchConfig, AutoStitchReport, ElementaryTrackValues,
    MaterializedStitch, PiecewiseStitchIndex, StitchParam, StitchSensitivityRow,
    CANONICAL_STITCH_BP, DEFAULT_SENSITIVITY_DISTANCES,
};
pub use tss::{is_primary_contig, AnnotationIndex, Transcript};
pub use types::{
    AssayPreset, CompatMode, Constituent, ControlSubtraction, CutoffMethod, FragmentMode,
    GeneBiotypeFilter, Interval, NormDenom, PromoterMode, RegionClass, RegionOrigin,
    StitchedRegion, Strand,
};
