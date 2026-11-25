//! Command-line argument definitions (`clap` derive) matching `PLAN.md` Section 4.

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use rose_core::{
    AssayPreset, CompatMode, ControlSubtraction, CutoffMethod, FragmentMode, GeneBiotypeFilter,
    NormDenom, PromoterMode, StitchParam,
};

/// Fast, self-contained Rust reimplementation of ROSE (Rank Ordering of Super-Enhancers).
#[derive(Debug, Parser)]
#[command(name = "rose", version, about)]
pub struct Cli {
    /// Worker thread count (`0` uses all available logical cores).
    #[arg(long, global = true, default_value_t = 0)]
    pub threads: usize,

    /// Subcommand to execute.
    #[command(subcommand)]
    pub command: Commands,
}

/// Top-level `rose` subcommands.
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Run the complete stitch -> signal -> call -> genes pipeline in memory.
    Run(RunArgs),
    /// Stitch constituent peaks into candidate enhancers (with optional TSS handling).
    Stitch(StitchCmdArgs),
    /// Quantify BAM/CRAM or bigWig signal across regions (`*_ENHANCER_REGION_MAP.txt`).
    Signal(SignalCmdArgs),
    /// Compute the super-enhancer cutoff, ranked tables, BED tracks, and hockey-stick plot.
    Call(CallCmdArgs),
    /// Assign overlapping, proximal, and closest genes to ranked enhancers.
    Genes(GenesCmdArgs),
    /// Score multiple samples and replicates on a shared consensus region set.
    Multi(MultiCmdArgs),
    /// Download or install UCSC RefSeq annotation tables into a local directory.
    FetchAnnotation(FetchAnnotationArgs),
}

/// CLI enum for `--assay`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum AssayArg {
    /// ChIP-seq defaults (`floor = 1`, subtract control when `-c` is given).
    Chip,
    /// CUT&RUN defaults (`floor = 0`, paired-end fragments, no IgG subtraction unless requested).
    Cutrun,
}

impl From<AssayArg> for AssayPreset {
    fn from(v: AssayArg) -> Self {
        match v {
            AssayArg::Chip => Self::Chip,
            AssayArg::Cutrun => Self::CutRun,
        }
    }
}

/// CLI enum for `--promoter-mode`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum PromoterModeArg {
    /// Drop constituent peaks wholly contained inside the TSS window (ROSE).
    Contained,
    /// Subtract the TSS window out of overlapping peaks (LILY).
    Subtract,
}

impl From<PromoterModeArg> for PromoterMode {
    fn from(v: PromoterModeArg) -> Self {
        match v {
            PromoterModeArg::Contained => Self::Contained,
            PromoterModeArg::Subtract => Self::Subtract,
        }
    }
}

/// CLI enum for `--signal`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum SignalSourceArg {
    /// Read alignments from BAM or CRAM.
    Bam,
    /// Read continuous signal density from UCSC bigWig.
    Bigwig,
}

/// CLI enum for `--fragment`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum FragmentArg {
    /// Auto-detect paired-end (`TLEN`) vs single-end (`CIGAR + extend`).
    Auto,
    /// Single-read reference-consuming CIGAR span + `--extend`.
    Cigar,
    /// Single-read `SEQ` length + `--extend` (legacy ROSE).
    Seq,
    /// Paired-end template length `[POS, POS + TLEN)`.
    Tlen,
}

impl From<FragmentArg> for FragmentMode {
    fn from(v: FragmentArg) -> Self {
        match v {
            FragmentArg::Auto => Self::Auto,
            FragmentArg::Cigar => Self::Cigar,
            FragmentArg::Seq => Self::Seq,
            FragmentArg::Tlen => Self::Tlen,
        }
    }
}

/// CLI enum for `--norm-denom`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum NormDenomArg {
    /// Post-filter fragment / read count.
    Filtered,
    /// Mapped alignment count from `.bai` / `.csi` index (`idxstats`).
    Index,
}

impl From<NormDenomArg> for NormDenom {
    fn from(v: NormDenomArg) -> Self {
        match v {
            NormDenomArg::Filtered => Self::Filtered,
            NormDenomArg::Index => Self::Index,
        }
    }
}

/// CLI enum for `--subtract-control`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum ToggleArg {
    /// Enable feature.
    On,
    /// Disable feature.
    Off,
}

impl From<ToggleArg> for ControlSubtraction {
    fn from(v: ToggleArg) -> Self {
        match v {
            ToggleArg::On => Self::On,
            ToggleArg::Off => Self::Off,
        }
    }
}

/// CLI enum for `--cutoff`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CutoffArg {
    /// Exact single-pass chord-tangent scan.
    Exact,
    /// Legacy R `stats::optimize` Brent step-function solver.
    RoseCompat,
}

impl From<CutoffArg> for CutoffMethod {
    fn from(v: CutoffArg) -> Self {
        match v {
            CutoffArg::Exact => Self::Exact,
            CutoffArg::RoseCompat => Self::RoseCompat,
        }
    }
}

/// CLI enum for `--compat`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CompatArg {
    /// Clean 0-based half-open modern behaviour.
    Clean,
    /// Enable all legacy Young Lab ROSE 0.1 quirks.
    Rose,
}

impl From<CompatArg> for CompatMode {
    fn from(v: CompatArg) -> Self {
        match v {
            CompatArg::Clean => Self::Clean,
            CompatArg::Rose => Self::RoseLegacy,
        }
    }
}

/// CLI enum for `--gene-biotype`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum BiotypeArg {
    /// Keep all annotated transcripts.
    All,
    /// Restrict to protein-coding transcripts (`protein_coding` / `NM_`).
    ProteinCoding,
}

impl From<BiotypeArg> for GeneBiotypeFilter {
    fn from(v: BiotypeArg) -> Self {
        match v {
            BiotypeArg::All => Self::All,
            BiotypeArg::ProteinCoding => Self::ProteinCoding,
        }
    }
}

/// Shared stitching and annotation options.
#[derive(Debug, Clone, Args)]
pub struct SharedStitchOpts {
    /// Maximum stitching distance in bp (`12500` default) or `auto` (signal-blind gap excess).
    #[arg(short = 's', long = "stitch", default_value_t = StitchParam::default())]
    pub stitch: StitchParam,

    /// Disable generating `*_Stitch_sensitivity.tsv` across `0..=25 kb` on `run` / `multi`.
    #[arg(long)]
    pub no_stitch_report: bool,

    /// TSS promoter exclusion half-window in bp (`0` disables promoter exclusion; `2500` recommended).
    #[arg(short = 't', long = "tss", default_value_t = 0)]
    pub tss: u64,

    /// Promoter handling mode (`contained` drop vs `subtract` window cut).
    #[arg(long, value_enum, default_value_t = PromoterModeArg::Contained)]
    pub promoter_mode: PromoterModeArg,

    /// Un-stitch regions overlapping strictly more distinct gene TSSs than this threshold.
    #[arg(long, default_value_t = 2)]
    pub max_tss_genes: usize,

    /// Biotype filter applied when loading GTF/refGene annotations.
    #[arg(long, value_enum, default_value_t = BiotypeArg::All)]
    pub gene_biotype: BiotypeArg,

    /// Optional BED file of blacklisted artifact regions removed before stitching.
    #[arg(long)]
    pub blacklist: Option<PathBuf>,

    /// Chromosomes excluded from stitching (comma-separated; pass `""` to keep all chromosomes).
    #[arg(long, default_value = "chrY,Y")]
    pub exclude_chrom: String,
}

/// Shared signal extraction options.
#[derive(Debug, Clone, Args)]
pub struct SharedSignalOpts {
    /// Experimental assay preset (`chip` or `cutrun`).
    #[arg(long, value_enum, default_value_t = AssayArg::Chip)]
    pub assay: AssayArg,

    /// Signal track format (`bam` for BAM/CRAM, `bigwig` for continuous density tracks).
    #[arg(long, value_enum)]
    pub signal: Option<SignalSourceArg>,

    /// Single-end 5'->3' read extension in bp.
    #[arg(long, default_value_t = 200)]
    pub extend: u64,

    /// Per-base coverage floor (defaults to `1` for `chip`, `0` for `cutrun`).
    #[arg(long)]
    pub floor: Option<u32>,

    /// Minimum mapping quality (`MAPQ`).
    #[arg(long, default_value_t = 0)]
    pub min_mapq: u8,

    /// Reference genome FASTA (indexed with `samtools faidx`) used to decode CRAM input.
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Excluded SAM bitwise flags (defaults to `3844` in clean mode, `4` in `--compat rose`).
    #[arg(long)]
    pub exclude_flags: Option<u16>,

    /// Fragment / read span derivation mode.
    #[arg(long, value_enum, default_value_t = FragmentArg::Auto)]
    pub fragment: FragmentArg,

    /// Minimum paired-end fragment length in bp.
    #[arg(long, default_value_t = 0)]
    pub min_fragment: u64,

    /// Maximum paired-end fragment length in bp.
    #[arg(long, default_value_t = 1_000)]
    pub max_fragment: u64,

    /// Scaling denominator (`filtered` post-QC units vs `index` `idxstats` alignments).
    #[arg(long, value_enum)]
    pub norm_denom: Option<NormDenomArg>,

    /// Explicit scale factor for the ranking track (replaces RPM scaling).
    #[arg(long)]
    pub scale_factor: Option<f64>,

    /// Explicit scale factor for the control track.
    #[arg(long)]
    pub control_scale_factor: Option<f64>,
}

/// Shared calling, stability, and copy-number options.
#[derive(Debug, Clone, Args)]
pub struct SharedCallOpts {
    /// Subtract control track when `-c` is supplied (`on` default for `chip`, `off` for `cutrun`).
    #[arg(long, value_enum)]
    pub subtract_control: Option<ToggleArg>,

    /// Run calling both with and without control subtraction and report call concordance.
    #[arg(long)]
    pub compare_control: bool,

    /// Tangent cutoff solver (`exact` or `rose-compat`).
    #[arg(long, value_enum)]
    pub cutoff: Option<CutoffArg>,

    /// Exclude top `N` extreme outlier signals when calculating the cutoff slope.
    #[arg(long, default_value_t = 0)]
    pub drop_top: usize,

    /// Enforce at least `N` super-enhancers (`0` = off; LILY uses `600`).
    #[arg(long, default_value_t = 0)]
    pub min_se: usize,

    /// Number of bootstrap resamples for cutoff stability & depth saturation (`0` = disabled).
    #[arg(long, default_value_t = 0)]
    pub stability: usize,

    /// Deterministic PRNG and QNAME-downsampling seed.
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    /// Estimate local copy-number ratios from fixed bins on the control BAM.
    #[arg(long)]
    pub cn_correct: bool,

    /// Apply externally computed copy-number segments (`.seg` or 4-col BED).
    #[arg(long)]
    pub cn_segments: Option<PathBuf>,

    /// Fixed bin size in bp for `--cn-correct`.
    #[arg(long, default_value_t = 100_000)]
    pub cn_bin: u64,

    /// Legacy compatibility preset (`--compat rose` enables all ROSE 0.1 quirks).
    #[arg(long, value_enum, default_value_t = CompatArg::Clean)]
    pub compat: CompatArg,
}

/// Arguments for `rose run`.
#[derive(Debug, Clone, Args)]
pub struct RunArgs {
    /// Constituent peaks file (BED, narrowPeak, broadPeak, or GFF).
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,

    /// Ranking BAM, CRAM, or bigWig track.
    #[arg(short = 'r', long = "rankby")]
    pub rankby: PathBuf,

    /// Optional control (Input / IgG) BAM, CRAM, or bigWig track.
    #[arg(short = 'c', long = "control")]
    pub control: Option<PathBuf>,

    /// Transcript annotation file (UCSC refGene `.ucsc` or `.gtf`).
    #[arg(short = 'a', long = "annotation")]
    pub annotation: Option<PathBuf>,

    /// Optional genome build shortcut (`hg38`, `hg19`, `mm10`, ...) resolved in local `annotations/`.
    #[arg(short = 'g', long = "genome")]
    pub genome: Option<String>,

    /// Output directory.
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,

    /// Proximal gene search window in bp.
    #[arg(short = 'w', long = "window", default_value_t = 50_000)]
    pub gene_window: u64,

    #[command(flatten)]
    pub stitch_opts: SharedStitchOpts,

    #[command(flatten)]
    pub signal_opts: SharedSignalOpts,

    #[command(flatten)]
    pub call_opts: SharedCallOpts,
}

/// Arguments for `rose stitch`.
#[derive(Debug, Clone, Args)]
pub struct StitchCmdArgs {
    /// Constituent peaks file (BED, narrowPeak, broadPeak, or GFF).
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,
    /// Optional annotation file (required when `--tss > 0`).
    #[arg(short = 'a', long = "annotation")]
    pub annotation: Option<PathBuf>,
    /// Output GFF/BED path (or directory).
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Legacy compatibility mode.
    #[arg(long, value_enum, default_value_t = CompatArg::Clean)]
    pub compat: CompatArg,
    #[command(flatten)]
    pub stitch_opts: SharedStitchOpts,
}

/// Arguments for `rose signal`.
#[derive(Debug, Clone, Args)]
pub struct SignalCmdArgs {
    /// Stitched or constituent regions file (BED or GFF).
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,
    /// Ranking BAM, CRAM, or bigWig track.
    #[arg(short = 'r', long = "rankby")]
    pub rankby: PathBuf,
    /// Optional control track.
    #[arg(short = 'c', long = "control")]
    pub control: Option<PathBuf>,
    /// Output `region_map.tsv` path.
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Legacy compatibility mode.
    #[arg(long, value_enum, default_value_t = CompatArg::Clean)]
    pub compat: CompatArg,
    #[command(flatten)]
    pub signal_opts: SharedSignalOpts,
}

/// Arguments for `rose call`.
#[derive(Debug, Clone, Args)]
pub struct CallCmdArgs {
    /// Input `region_map.tsv` / `*_ENHANCER_REGION_MAP.txt`.
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,
    /// Output directory.
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Sample / prefix name for output files.
    #[arg(short = 'n', long = "name")]
    pub name: Option<String>,
    #[command(flatten)]
    pub call_opts: SharedCallOpts,
}

/// Arguments for `rose genes`.
#[derive(Debug, Clone, Args)]
pub struct GenesCmdArgs {
    /// Ranked enhancer table (`*_AllEnhancers.table.txt` or `*_SuperEnhancers.table.txt`).
    #[arg(short = 'i', long = "input")]
    pub input: PathBuf,
    /// Transcript annotation file (UCSC refGene `.ucsc` or `.gtf`).
    #[arg(short = 'a', long = "annotation")]
    pub annotation: PathBuf,
    /// Output directory.
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Proximal TSS search window in bp.
    #[arg(short = 'w', long = "window", default_value_t = 50_000)]
    pub window: u64,
    /// Biotype filter.
    #[arg(long, value_enum, default_value_t = BiotypeArg::All)]
    pub gene_biotype: BiotypeArg,
    /// Legacy compatibility mode.
    #[arg(long, value_enum, default_value_t = CompatArg::Clean)]
    pub compat: CompatArg,
}

/// Arguments for `rose multi`.
#[derive(Debug, Clone, Args)]
pub struct MultiCmdArgs {
    /// Sample sheet TSV (`sample`, `group`, `peaks`, `bam`, optional `control`, `scale_factor`).
    #[arg(long = "samples")]
    pub samples: PathBuf,
    /// Transcript annotation file (UCSC refGene `.ucsc` or `.gtf`).
    #[arg(short = 'a', long = "annotation")]
    pub annotation: Option<PathBuf>,
    /// Output directory.
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Minimum number of samples supporting a merged constituent peak before stitching.
    #[arg(long, default_value_t = 1)]
    pub min_peak_samples: usize,
    /// Minimum replicates within a group required to call a consensus super-enhancer (default: all).
    #[arg(long)]
    pub min_replicates: Option<usize>,
    #[command(flatten)]
    pub stitch_opts: SharedStitchOpts,
    #[command(flatten)]
    pub signal_opts: SharedSignalOpts,
    #[command(flatten)]
    pub call_opts: SharedCallOpts,
}

/// Arguments for `rose fetch-annotation`.
#[derive(Debug, Clone, Args)]
pub struct FetchAnnotationArgs {
    /// Genome build (`hg18`, `hg19`, `hg38`, `mm8`, `mm9`, `mm10`).
    #[arg(short = 'g', long = "genome")]
    pub genome: String,
    /// Target directory to store `<genome>_refseq.ucsc`.
    #[arg(short = 'o', long = "out", default_value = "annotations")]
    pub out: PathBuf,
    /// Optional local directory of pre-downloaded UCSC files (offline HPC install).
    #[arg(long)]
    pub from_dir: Option<PathBuf>,
}
