//! Genomic format readers, BAM/CRAM/bigWig adapters, and atomic table writers for `rose-rs`.

pub mod annotation;
pub mod atomic;
pub mod bam;
pub mod bigwig;
pub mod error;
pub mod peaks;
pub mod tables;

pub use annotation::read_annotation_index;
pub use atomic::{sha256_file, write_atomic_bytes, AtomicWriter};
pub use bam::{
    compute_bam_signals, inspect_alignment_header, AlignmentHeaderInfo, BamSignalConfig,
    BamSignalOutput, DEFAULT_EXCLUDE_FLAGS,
};
pub use bigwig::{write_bigwig, BigWigFile};
pub use error::RoseIoError;
pub use peaks::{
    read_blacklist_bed, read_cn_segments, read_constituents, read_sample_sheet, RegionFileFormat,
    SampleSheetRow,
};
pub use tables::{
    read_ranked_enhancers_tsv, read_region_map_tsv, write_enhancer_beds,
    write_enhancer_call_tables, write_gene_tables, write_multi_outputs, write_points_and_manifest,
    write_region_map_tsv, write_stitch_sensitivity_tsv, write_stitched_gff, RunManifest,
};
