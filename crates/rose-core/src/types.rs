//! Canonical 0-based half-open genomic intervals, strand types, and configuration enums.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::error::RoseCoreError;

/// Canonical 0-based half-open `[start, end)` genomic coordinate span (`start < end`).
///
/// # Examples
///
/// ```
/// use rose_core::Interval;
///
/// # fn main() -> Result<(), rose_core::RoseCoreError> {
/// let a = Interval::try_new("chr1", 1_000, 2_000)?;
/// let b = Interval::try_new("chr1", 1_500, 2_500)?;
/// assert_eq!(a.width(), 1_000);
/// assert_eq!(a.overlap_len(b), 500);
/// # Ok(())
/// # }
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Interval {
    start: u64,
    end: u64,
}

impl Interval {
    /// Constructs a validated 0-based half-open `[start, end)` interval.
    ///
    /// # Errors
    ///
    /// Returns [`RoseCoreError::InvalidInterval`] if `start >= end`.
    pub fn try_new(chrom: &str, start: u64, end: u64) -> Result<Self, RoseCoreError> {
        if start >= end {
            return Err(RoseCoreError::InvalidInterval {
                chrom: chrom.to_owned(),
                start,
                end,
            });
        }
        Ok(Self { start, end })
    }

    /// Constructs a 1-bp half-open `[pos, pos + 1)` point interval (e.g. a TSS coordinate).
    #[must_use]
    #[inline]
    pub const fn point(pos: u64) -> Self {
        Self {
            start: pos,
            end: pos.saturating_add(1),
        }
    }

    /// Inclusive 0-based start coordinate.
    #[must_use]
    #[inline]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Exclusive 0-based end coordinate.
    #[must_use]
    #[inline]
    pub const fn end(self) -> u64 {
        self.end
    }

    /// Length in base pairs (`end - start >= 1`).
    #[must_use]
    #[inline]
    pub const fn width(self) -> u64 {
        self.end - self.start
    }

    /// Midpoint `(start + end) / 2` in 0-based coordinates.
    #[must_use]
    #[inline]
    pub const fn center(self) -> u64 {
        self.start + (self.end - self.start) / 2
    }

    /// Returns `true` when `self` and `other` share at least 1 base pair.
    #[must_use]
    #[inline]
    pub const fn overlaps(self, other: Self) -> bool {
        self.start < other.end && other.start < self.end
    }

    /// Returns `true` when every base pair of `other` lies inside `self`.
    #[must_use]
    #[inline]
    pub const fn contains(self, other: Self) -> bool {
        self.start <= other.start && other.end <= self.end
    }

    /// Shared overlap width in base pairs (`0` if disjoint).
    #[must_use]
    #[inline]
    pub fn overlap_len(self, other: Self) -> u64 {
        let lo = self.start.max(other.start);
        let hi = self.end.min(other.end);
        hi.saturating_sub(lo)
    }

    /// Expands `[start, end)` symmetrically by `pad` bp, saturating at `0`.
    #[must_use]
    #[inline]
    pub const fn expanded(self, pad: u64) -> Self {
        let lo = self.start.saturating_sub(pad);
        let hi = self.end.saturating_add(pad);
        Self {
            start: lo,
            end: if hi <= lo { lo + 1 } else { hi },
        }
    }

    /// Subtracts `cut` from `self`, returning `0`, `1`, or `2` remaining sub-intervals.
    #[must_use]
    pub fn subtract(self, cut: Self) -> Vec<Self> {
        if !self.overlaps(cut) {
            return vec![self];
        }
        let mut rem = Vec::with_capacity(2);
        if self.start < cut.start {
            rem.push(Self {
                start: self.start,
                end: cut.start,
            });
        }
        if cut.end < self.end {
            rem.push(Self {
                start: cut.end,
                end: self.end,
            });
        }
        rem
    }
}

/// Genomic strand orientation (`+`, `-`, or `.`).
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default, Serialize, Deserialize,
)]
pub enum Strand {
    /// Plus / Watson strand (`+`).
    Forward,
    /// Minus / Crick strand (`-`).
    Reverse,
    /// Unstranded or ambidextrous (`.`).
    #[default]
    Unstranded,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ch = match self {
            Self::Forward => '+',
            Self::Reverse => '-',
            Self::Unstranded => '.',
        };
        write!(f, "{ch}")
    }
}

impl FromStr for Strand {
    type Err = RoseCoreError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim() {
            "+" | "1" | "sense" => Ok(Self::Forward),
            "-" | "-1" | "antisense" => Ok(Self::Reverse),
            "." | "" | "*" | "both" => Ok(Self::Unstranded),
            other => Err(RoseCoreError::InvalidConfig(format!(
                "invalid strand token '{other}'"
            ))),
        }
    }
}

/// Input constituent enhancer or peak region.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Constituent {
    /// Unique constituent identifier.
    pub id: String,
    /// Chromosome / contig name.
    pub chrom: String,
    /// 0-based half-open coordinate span.
    pub interval: Interval,
    /// Strand annotation.
    pub strand: Strand,
    /// Optional input peak score (e.g. MACS3 signalValue).
    pub score: Option<f64>,
}

/// Biological region origin tracked for LILY-style classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
pub enum RegionOrigin {
    /// Distal or stitched enhancer candidate.
    #[default]
    Enhancer,
    /// Promoter-proximal element retained for annotation.
    Promoter,
}

/// Stitched (or disaggregated) candidate regulatory region produced by [`crate::stitch`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StitchedRegion {
    /// Stitched region identifier (`<count>_<leftmost_id>_lociStitched` or original ID).
    pub id: String,
    /// Chromosome / contig name.
    pub chrom: String,
    /// 0-based half-open span of the stitched locus.
    pub interval: Interval,
    /// Strand (`.` when stitched across constituents).
    pub strand: Strand,
    /// Number of constituent peaks merged into this region (`NUM_LOCI`).
    pub num_loci: usize,
    /// Sum of widths of original constituents overlapping this region (`CONSTITUENT_SIZE`).
    pub constituent_size: u64,
    /// Ordered IDs of merged constituents.
    pub constituent_ids: Vec<String>,
    /// Distal enhancer vs promoter origin.
    pub origin: RegionOrigin,
}

/// Final LILY/ROSE regulatory class assigned after cutoff calling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RegionClass {
    /// Super-enhancer (`signal > cutoff`).
    SuperEnhancer,
    /// Typical enhancer (`signal <= cutoff`).
    Enhancer,
    /// Promoter-overlapping element.
    Promoter,
}

impl fmt::Display for RegionClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SuperEnhancer => write!(f, "super_enhancer"),
            Self::Enhancer => write!(f, "enhancer"),
            Self::Promoter => write!(f, "promoter"),
        }
    }
}

/// High-level experimental assay preset (`--assay`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AssayPreset {
    /// Chromatin immunoprecipitation sequencing (`ChIP-seq`).
    #[default]
    Chip,
    /// Cleavage Under Targets and Release Using Nuclease (`CUT&RUN`).
    CutRun,
}

/// Promoter TSS exclusion policy (`--promoter-mode`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum PromoterMode {
    /// Drop constituent peaks wholly contained inside `±tss_window` (ROSE default).
    #[default]
    Contained,
    /// Subtract `±tss_window` intervals out of overlapping peaks (LILY mode).
    Subtract,
}

/// Clean modern semantics vs strict legacy ROSE 0.1 compatibility (`--compat rose`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CompatMode {
    /// Clean 0-based half-open semantics, full `f64` precision, and exact chord cutoff.
    #[default]
    Clean,
    /// Legacy Young Lab ROSE 0.1 quirks (1-bp GFF `end - start`, `±50 bp` TSS un-stitching,
    /// step depth floor, 4-decimal MMR/density rounding, Brent `optimize` cutoff).
    RoseLegacy,
}

/// Tangent cutoff algorithm (`--cutoff`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum CutoffMethod {
    /// Exact single-pass chord-tangent scan with rightmost tie-breaking.
    #[default]
    Exact,
    /// Port of R `stats::optimize` (Brent's method) on `numPts_below_line`.
    RoseCompat,
}

/// Library-size scaling denominator source (`--norm-denom`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum NormDenom {
    /// Divide by post-filter unique fragments (paired-end) or reads (single-end).
    #[default]
    Filtered,
    /// Divide by `.bai`/`.csi` `idxstats` mapped alignment records (`flagstat` compatible).
    Index,
}

/// Read/fragment span derivation policy (`--fragment`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum FragmentMode {
    /// Paired-end proper pairs use `[POS, POS + TLEN)`; single-end uses CIGAR + `--extend`.
    #[default]
    Auto,
    /// Single-read reference-consuming CIGAR span + `--extend`.
    Cigar,
    /// Single-read `SEQ` length + `--extend` (legacy ROSE `readsToLoci` behaviour).
    Seq,
    /// Paired-end template span `[POS, POS + TLEN)` only.
    Tlen,
}

/// Transcript biotype filter applied to GTF annotations (`--gene-biotype`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GeneBiotypeFilter {
    /// Retain all transcript biotypes.
    #[default]
    All,
    /// Retain only protein-coding transcripts (`protein_coding` / `NM_`).
    ProteinCoding,
}

/// Control background subtraction toggle (`--subtract-control`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ControlSubtraction {
    /// Subtract control signal from ranking signal and clamp negatives to `0.0`.
    #[default]
    On,
    /// Ignore control signal during ranking (default for `--assay cutrun`).
    Off,
}
