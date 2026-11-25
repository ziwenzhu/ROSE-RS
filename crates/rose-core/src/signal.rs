//! Coverage accumulation, depth-floor difference arrays, nested `xxh3` downsampling,
//! bigWig weighted sums, and chromosome-name reconciliation.

use std::collections::{HashMap, HashSet};

use xxhash_rust::xxh3::xxh3_64_with_seed;

use crate::error::RoseCoreError;
use crate::types::{CompatMode, Interval};

/// Single decoded read or paired-end fragment span in 0-based half-open coordinates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadSpan {
    /// 0-based half-open genomic span `[start, end)`.
    pub interval: Interval,
    /// Deterministic `xxh3_64` hash of `(seed, qname)` mapped into `[0.0, 1.0)` via `u64::MAX`.
    pub subsample_hash: u64,
}

impl ReadSpan {
    /// Computes deterministic `xxh3_64` hash from `qname` bytes and `seed`.
    #[must_use]
    #[inline]
    pub fn hash_qname(qname: &[u8], seed: u64) -> u64 {
        xxh3_64_with_seed(qname, seed)
    }

    /// Returns `true` when this read/fragment belongs to nested downsampling `fraction` in `(0.0, 1.0]`.
    #[must_use]
    #[inline]
    pub fn in_fraction(&self, fraction: f64) -> bool {
        if fraction >= 1.0 {
            return true;
        }
        if fraction <= 0.0 {
            return false;
        }
        let threshold = (fraction * (u64::MAX as f64)) as u64;
        self.subsample_hash <= threshold
    }
}

/// Reusable per-thread scratch workspace for region coverage with depth flooring (`floor >= 1`)
/// and optional nested downsampling tiers.
#[derive(Debug, Default, Clone)]
pub struct SignalScratch {
    diff: Vec<i32>,
    tier_diffs: Vec<Vec<i32>>,
}

impl SignalScratch {
    /// Creates an empty reusable scratch workspace.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Computes base-pair coverage total over `region` from overlapping `reads` given `floor`
    /// and `compat` mode.
    #[must_use]
    pub fn compute_region_bp_sum(
        &mut self,
        region: Interval,
        reads: &[ReadSpan],
        floor: u32,
        compat: CompatMode,
    ) -> u64 {
        if reads.is_empty() {
            return 0;
        }
        let Some(region) = counted_region(region, compat) else {
            return 0;
        };
        if floor == 0 {
            return reads.iter().map(|r| r.interval.overlap_len(region)).sum();
        }

        let width = usize::try_from(region.width()).unwrap_or(0);
        if width == 0 {
            return 0;
        }
        self.diff.clear();
        self.diff.resize(width + 1, 0);

        let reg_start = region.start();
        let reg_end = region.end();
        for r in reads {
            if !r.interval.overlaps(region) {
                continue;
            }
            let lo = (r.interval.start().max(reg_start) - reg_start) as usize;
            let hi = (r.interval.end().min(reg_end) - reg_start) as usize;
            self.diff[lo] += 1;
            self.diff[hi] -= 1;
        }

        let floor_i32 = floor as i32;
        let mut depth = 0_i32;
        let mut sum = 0_u64;
        for &delta in &self.diff[..width] {
            depth += delta;
            // Both modes follow `ROSE_bamToGFF.py:113`: bases with `depth > floor` contribute
            // their full depth (not `depth - floor`).
            if depth > floor_i32 {
                sum += depth as u64;
            }
        }
        sum
    }

    /// Simultaneously computes base-pair coverage sums across multiple nested downsampling
    /// `fractions` (e.g. `[0.25, 0.50, 0.75, 0.90, 1.0]`) in a single pass over `reads`.
    #[must_use]
    pub fn compute_multi_tier_bp_sums(
        &mut self,
        region: Interval,
        reads: &[ReadSpan],
        floor: u32,
        fractions: &[f64],
        compat: CompatMode,
    ) -> Vec<u64> {
        let k = fractions.len();
        if k == 0 {
            return Vec::new();
        }
        if reads.is_empty() {
            return vec![0; k];
        }
        let Some(region) = counted_region(region, compat) else {
            return vec![0; k];
        };
        if floor == 0 {
            let mut sums = vec![0_u64; k];
            for r in reads {
                let ov = r.interval.overlap_len(region);
                if ov > 0 {
                    for (t_idx, &frac) in fractions.iter().enumerate() {
                        if r.in_fraction(frac) {
                            sums[t_idx] += ov;
                        }
                    }
                }
            }
            return sums;
        }

        let width = usize::try_from(region.width()).unwrap_or(0);
        if width == 0 {
            return vec![0; k];
        }
        if self.tier_diffs.len() < k {
            self.tier_diffs.resize_with(k, Vec::new);
        }
        for td in &mut self.tier_diffs[..k] {
            td.clear();
            td.resize(width + 1, 0);
        }

        let reg_start = region.start();
        let reg_end = region.end();
        for r in reads {
            if !r.interval.overlaps(region) {
                continue;
            }
            let lo = (r.interval.start().max(reg_start) - reg_start) as usize;
            let hi = (r.interval.end().min(reg_end) - reg_start) as usize;
            for (t_idx, &frac) in fractions.iter().enumerate() {
                if r.in_fraction(frac) {
                    self.tier_diffs[t_idx][lo] += 1;
                    self.tier_diffs[t_idx][hi] -= 1;
                }
            }
        }

        let floor_i32 = floor as i32;
        let mut out = Vec::with_capacity(k);
        for td in &self.tier_diffs[..k] {
            let mut depth = 0_i32;
            let mut sum = 0_u64;
            for &delta in &td[..width] {
                depth += delta;
                if depth > floor_i32 {
                    sum += depth as u64;
                }
            }
            out.push(sum);
        }
        out
    }
}

/// Bases that are counted for `region`. Legacy ROSE (`ROSE_bamToGFF.py:115`) keeps only bases
/// strictly inside the GFF locus (`start < x < end`), which on the legacy number line is
/// `[start + 1, end)`. The density divisor stays `end - start`.
#[must_use]
#[inline]
pub fn counted_region(region: Interval, compat: CompatMode) -> Option<Interval> {
    match compat {
        CompatMode::Clean => Some(region),
        // A 1 bp locus has no interior bases; ROSE reports `NA` and later treats it as no signal.
        CompatMode::RoseLegacy => Interval::try_new("", region.start() + 1, region.end()).ok(),
    }
}

/// Converts raw base-pair coverage sum `bp_sum` over `width` into normalized density (`rpm/bp`)
/// and total integrated signal (`rpm * bp`).
#[must_use]
pub fn normalize_coverage(
    bp_sum: u64,
    width: u64,
    denom_millions: f64,
    scale_factor: Option<f64>,
    compat: CompatMode,
) -> (f64, f64) {
    if width == 0 {
        return (0.0, 0.0);
    }
    let w = width as f64;
    let raw_density = (bp_sum as f64) / w;
    let norm_density = if let Some(sf) = scale_factor {
        raw_density * sf
    } else {
        let mmr = match compat {
            CompatMode::Clean => denom_millions.max(1e-12),
            CompatMode::RoseLegacy => round_4(denom_millions).max(1e-4),
        };
        let d = raw_density / mmr;
        match compat {
            CompatMode::Clean => d,
            CompatMode::RoseLegacy => round_4(d),
        }
    };
    let total_signal = norm_density * w;
    (norm_density, total_signal)
}

#[inline]
fn round_4(val: f64) -> f64 {
    (val * 10_000.0).round() / 10_000.0
}

/// Integrates width-weighted bigWig intervals `(iv, value)` overlapping `region`, returning
/// `(mean_density, total_signal)`.
#[must_use]
pub fn integrate_bigwig_values(region: Interval, values: &[(Interval, f64)]) -> (f64, f64) {
    let w = region.width();
    if w == 0 {
        return (0.0, 0.0);
    }
    let total: f64 = values
        .iter()
        .map(|&(iv, val)| (iv.overlap_len(region) as f64) * val)
        .sum();
    let total_clamped = if total.is_finite() {
        total.max(0.0)
    } else {
        0.0
    };
    (total_clamped / (w as f64), total_clamped)
}

/// Reconciles region chromosome names (`chr1` vs `1`, `chrM` vs `MT`) against `header_contigs`.
///
/// # Errors
///
/// Returns [`RoseCoreError::ContigMismatch`] if fewer than `min_match_share` (e.g. `0.5`) of
/// `region_chroms` match any contig in `header_contigs`.
pub fn reconcile_chrom_names(
    region_chroms: &[String],
    header_contigs: &[String],
    min_match_share: f64,
) -> Result<HashMap<String, String>, RoseCoreError> {
    if region_chroms.is_empty() {
        return Ok(HashMap::new());
    }
    let contig_set: HashSet<&str> = header_contigs.iter().map(String::as_str).collect();
    let mut mapping: HashMap<String, String> = HashMap::new();
    let mut matched = 0_usize;

    for rc in region_chroms {
        if let Some(resolved) = mapping.get(rc) {
            if contig_set.contains(resolved.as_str()) {
                matched += 1;
            }
            continue;
        }
        let candidate = if contig_set.contains(rc.as_str()) {
            Some(rc.clone())
        } else if let Some(stripped) = rc.strip_prefix("chr") {
            if contig_set.contains(stripped) {
                Some(stripped.to_owned())
            } else if stripped == "M" && contig_set.contains("MT") {
                Some("MT".to_owned())
            } else {
                None
            }
        } else {
            let with_chr = if rc == "MT" {
                "chrM".to_owned()
            } else {
                format!("chr{rc}")
            };
            contig_set.contains(with_chr.as_str()).then_some(with_chr)
        };

        if let Some(target) = candidate {
            matched += 1;
            mapping.insert(rc.clone(), target);
        } else {
            mapping.insert(rc.clone(), rc.clone());
        }
    }

    let total = region_chroms.len();
    let share = (matched as f64) / (total as f64);
    if share < min_match_share {
        return Err(RoseCoreError::ContigMismatch {
            matched,
            total,
            pct: share * 100.0,
            min_pct: min_match_share * 100.0,
        });
    }
    Ok(mapping)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(s: u64, e: u64) -> Interval {
        Interval::try_new("t", s, e).expect("valid interval")
    }

    fn span(s: u64, e: u64) -> ReadSpan {
        ReadSpan {
            interval: iv(s, e),
            subsample_hash: 0,
        }
    }

    #[test]
    fn legacy_counts_only_strict_interior_bases() {
        // ROSE keeps bases `start < x < end`; a read covering the whole region adds width - 1.
        let region = iv(100, 110);
        let reads = [span(90, 120)];
        let mut s = SignalScratch::new();
        assert_eq!(
            s.compute_region_bp_sum(region, &reads, 0, CompatMode::Clean),
            10
        );
        assert_eq!(
            s.compute_region_bp_sum(region, &reads, 0, CompatMode::RoseLegacy),
            9
        );
        // With floor 1 the same rule applies to bases with depth > 1, which add full depth.
        let two = [span(90, 120), span(90, 120)];
        assert_eq!(
            s.compute_region_bp_sum(region, &two, 1, CompatMode::Clean),
            20
        );
        assert_eq!(
            s.compute_region_bp_sum(region, &two, 1, CompatMode::RoseLegacy),
            18
        );
    }

    #[test]
    fn legacy_one_bp_region_has_no_interior() {
        let mut s = SignalScratch::new();
        let reads = [span(0, 1000)];
        assert_eq!(
            s.compute_region_bp_sum(iv(100, 101), &reads, 0, CompatMode::RoseLegacy),
            0
        );
        assert_eq!(
            s.compute_multi_tier_bp_sums(
                iv(100, 101),
                &reads,
                1,
                &[0.5, 1.0],
                CompatMode::RoseLegacy
            ),
            vec![0, 0]
        );
    }
}
