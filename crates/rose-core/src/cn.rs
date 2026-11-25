//! Built-in copy-number estimation (`--cn-correct` from 100 kb control bins) and external
//! segment (`--cn-segments`) lookup.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error::RoseCoreError;
use crate::types::{Interval, StitchedRegion};

/// Parameters controlling bin-based copy-number estimation from a control BAM.
#[derive(Debug, Clone, PartialEq)]
pub struct CopyNumberConfig {
    /// Fixed genomic bin width in bp (`100_000` by default, matching HMCan).
    pub bin_size: u64,
    /// Running-median smoothing window across neighbouring mappable bins (`5` by default).
    pub smooth_window_bins: usize,
    /// Minimum autosomal median reads per bin required to permit estimation (`50.0` by default).
    pub min_median_reads_per_bin: f64,
    /// Fraction of raw positive-bin median below which a bin is marked unmappable (`0.10`).
    pub unmappable_floor_fraction: f64,
    /// Minimum clamped copy-number ratio (`0.25`).
    pub min_cn_ratio: f64,
    /// Maximum clamped copy-number ratio (`20.0`).
    pub max_cn_ratio: f64,
}

impl Default for CopyNumberConfig {
    fn default() -> Self {
        Self {
            bin_size: 100_000,
            smooth_window_bins: 5,
            min_median_reads_per_bin: 50.0,
            unmappable_floor_fraction: 0.10,
            min_cn_ratio: 0.25,
            max_cn_ratio: 20.0,
        }
    }
}

/// External or estimated copy-number segment (`chrom`, `interval`, `cn_ratio`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CopyNumberSegment {
    /// Chromosome / contig name.
    pub chrom: String,
    /// 0-based half-open segment span.
    pub interval: Interval,
    /// Relative copy-number ratio (`1.0` = diploid / baseline).
    pub cn_ratio: f64,
    /// `true` when mappable and estimated; `false` when unmappable/blacklisted (`cn_ratio = 1.0`).
    pub is_mappable: bool,
}

/// Indexed copy-number profile supporting fast overlap-weighted region ratio queries.
#[derive(Debug, Clone, Default)]
pub struct CopyNumberProfile {
    by_chrom: BTreeMap<String, Vec<CopyNumberSegment>>,
    autosomal_median_reads: f64,
}

impl CopyNumberProfile {
    /// Constructs a [`CopyNumberProfile`] from external `.seg` / BED copy-number segments.
    #[must_use]
    pub fn from_segments(segments: Vec<CopyNumberSegment>, min_ratio: f64, max_ratio: f64) -> Self {
        let mut by_chrom: BTreeMap<String, Vec<CopyNumberSegment>> = BTreeMap::new();
        for mut seg in segments {
            if seg.is_mappable && seg.cn_ratio.is_finite() && seg.cn_ratio > 0.0 {
                seg.cn_ratio = seg.cn_ratio.clamp(min_ratio, max_ratio);
            } else {
                seg.cn_ratio = 1.0;
                seg.is_mappable = false;
            }
            by_chrom.entry(seg.chrom.clone()).or_default().push(seg);
        }
        for list in by_chrom.values_mut() {
            list.sort_by_key(|s| (s.interval.start(), s.interval.end()));
        }
        Self {
            by_chrom,
            autosomal_median_reads: 0.0,
        }
    }

    /// Estimates a smoothed [`CopyNumberProfile`] from per-chromosome 100 kb bin fragment counts.
    ///
    /// # Errors
    ///
    /// Returns [`RoseCoreError::InsufficientControlDepth`] when the control BAM is too shallow
    /// (e.g. sparse CUT&RUN IgG) to support stable bin ratios.
    pub fn from_bin_counts(
        bin_counts: &BTreeMap<String, Vec<u64>>,
        contig_lengths: &BTreeMap<String, u64>,
        blacklist: &BTreeMap<String, Vec<Interval>>,
        cfg: &CopyNumberConfig,
    ) -> Result<Self, RoseCoreError> {
        let mut auto_positive: Vec<f64> = Vec::new();
        for (chrom, counts) in bin_counts {
            if is_sex_or_mito(chrom) {
                continue;
            }
            let bl = blacklist.get(chrom).map_or(&[][..], Vec::as_slice);
            for (b_idx, &cnt) in counts.iter().enumerate() {
                let start = (b_idx as u64) * cfg.bin_size;
                let len = contig_lengths
                    .get(chrom)
                    .copied()
                    .unwrap_or((counts.len() as u64) * cfg.bin_size);
                let end = (start + cfg.bin_size).min(len).max(start + 1);
                if let Ok(iv) = Interval::try_new(chrom, start, end) {
                    if !overlaps_blacklist(bl, iv) && cnt > 0 {
                        auto_positive.push(cnt as f64);
                    }
                }
            }
        }

        let raw_pos_median = median_slice(&mut auto_positive);
        let mappable_floor = raw_pos_median * cfg.unmappable_floor_fraction;
        let mut auto_mappable: Vec<f64> = auto_positive
            .into_iter()
            .filter(|&v| v >= mappable_floor)
            .collect();
        let auto_median = median_slice(&mut auto_mappable);

        if auto_median < cfg.min_median_reads_per_bin {
            return Err(RoseCoreError::InsufficientControlDepth {
                median_reads: auto_median,
                min_required: cfg.min_median_reads_per_bin,
            });
        }

        let mut by_chrom: BTreeMap<String, Vec<CopyNumberSegment>> = BTreeMap::new();
        for (chrom, counts) in bin_counts {
            let bl = blacklist.get(chrom).map_or(&[][..], Vec::as_slice);
            let chr_len = contig_lengths
                .get(chrom)
                .copied()
                .unwrap_or((counts.len() as u64) * cfg.bin_size);

            // Compute chromosome-specific baseline for sex chromosomes; autosomal median otherwise.
            let baseline = if is_sex_or_mito(chrom) {
                let mut sex_vals: Vec<f64> = counts
                    .iter()
                    .copied()
                    .filter(|&c| (c as f64) >= mappable_floor * 0.5 && c > 0)
                    .map(|c| c as f64)
                    .collect();
                let sm = median_slice(&mut sex_vals);
                if sm >= cfg.min_median_reads_per_bin * 0.25 {
                    sm
                } else {
                    auto_median
                }
            } else {
                auto_median
            };

            let mut raw_ratios: Vec<Option<f64>> = Vec::with_capacity(counts.len());
            let mut intervals: Vec<Interval> = Vec::with_capacity(counts.len());
            for (b_idx, &cnt) in counts.iter().enumerate() {
                let start = (b_idx as u64) * cfg.bin_size;
                let end = (start + cfg.bin_size).min(chr_len).max(start + 1);
                let iv =
                    Interval::try_new(chrom, start, end).unwrap_or_else(|_| Interval::point(start));
                intervals.push(iv);
                if overlaps_blacklist(bl, iv) || (cnt as f64) < mappable_floor {
                    raw_ratios.push(None);
                } else {
                    raw_ratios.push(Some((cnt as f64) / baseline));
                }
            }

            let half_w = cfg.smooth_window_bins / 2;
            let mut segs = Vec::with_capacity(counts.len());
            for idx in 0..counts.len() {
                if raw_ratios[idx].is_none() {
                    segs.push(CopyNumberSegment {
                        chrom: chrom.clone(),
                        interval: intervals[idx],
                        cn_ratio: 1.0,
                        is_mappable: false,
                    });
                    continue;
                }
                let lo = idx.saturating_sub(half_w);
                let hi = (idx + half_w + 1).min(counts.len());
                let mut win: Vec<f64> = raw_ratios[lo..hi].iter().filter_map(|&r| r).collect();
                let smoothed = median_slice(&mut win).clamp(cfg.min_cn_ratio, cfg.max_cn_ratio);
                segs.push(CopyNumberSegment {
                    chrom: chrom.clone(),
                    interval: intervals[idx],
                    cn_ratio: smoothed,
                    is_mappable: true,
                });
            }
            by_chrom.insert(chrom.clone(), segs);
        }

        Ok(Self {
            by_chrom,
            autosomal_median_reads: auto_median,
        })
    }

    /// Autosomal median bin count observed (`0.0` when loaded from external segments).
    #[must_use]
    #[inline]
    pub const fn autosomal_median_reads(&self) -> f64 {
        self.autosomal_median_reads
    }

    /// Computes overlap-weighted local copy-number ratio for `region` (`1.0` when unmappable or
    /// absent).
    #[must_use]
    pub fn ratio_for_region(&self, region: &StitchedRegion) -> f64 {
        let Some(segs) = self.by_chrom.get(&region.chrom) else {
            return 1.0;
        };
        let q = region.interval;
        let upper = segs.partition_point(|s| s.interval.start() < q.end());
        let mut weighted_sum = 0.0_f64;
        let mut weight_bp = 0_u64;
        for seg in segs[..upper].iter().rev() {
            if seg.interval.end() <= q.start() {
                // Segments are sorted by `start`, though widths may vary for `.seg` files, so don't
                // break early on `.seg` unless non-overlapping; continue scanning overlapping range.
                continue;
            }
            if seg.is_mappable {
                let ov = seg.interval.overlap_len(q);
                if ov > 0 {
                    weighted_sum += (ov as f64) * seg.cn_ratio;
                    weight_bp += ov;
                }
            }
        }
        if weight_bp == 0 {
            1.0
        } else {
            weighted_sum / (weight_bp as f64)
        }
    }
}

fn is_sex_or_mito(chrom: &str) -> bool {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    matches!(bare, "X" | "Y" | "M" | "MT")
}

fn overlaps_blacklist(sorted_bl: &[Interval], q: Interval) -> bool {
    let upper = sorted_bl.partition_point(|iv| iv.start() < q.end());
    sorted_bl[..upper].iter().any(|iv| iv.overlaps(q))
}

fn median_slice(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.sort_by(f64::total_cmp);
    let mid = vals.len() / 2;
    if vals.len() % 2 == 1 {
        vals[mid]
    } else {
        0.5 * (vals[mid - 1] + vals[mid])
    }
}
