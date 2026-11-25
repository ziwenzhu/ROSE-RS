//! Cutoff stability bootstrap (`--stability N`) and nested read-downsampling saturation curves.

use serde::{Deserialize, Serialize};

use crate::cutoff::compute_cutoff_on_slice;
use crate::types::{CompatMode, CutoffMethod};

/// Default nested downsampling read fractions evaluated in a single BAM pass when `--stability` is active.
pub const DEFAULT_DOWNSAMPLE_FRACTIONS: [f64; 5] = [0.25, 0.50, 0.75, 0.90, 1.00];

/// Per-region bootstrap call stability metrics.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionStability {
    /// Share of `n_bootstrap` resamples (`0.0..=1.0`) in which this region exceeded the resampled cutoff.
    pub call_frequency: f64,
    /// `true` when `call_frequency` lies in `[0.10, 0.90]` (fragile / borderline call).
    pub is_borderline: bool,
}

/// Single point on the sequencing-depth saturation curve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DepthSaturationPoint {
    /// Subsampled read/fragment fraction (`0.0..=1.0`).
    pub read_fraction: f64,
    /// Cutoff threshold computed at this depth.
    pub cutoff_signal: f64,
    /// Number of super-enhancers called at this depth.
    pub super_enhancer_count: usize,
    /// Jaccard similarity of called super-enhancers vs the 100%-depth call set.
    pub jaccard_vs_full: f64,
}

/// Aggregate stability report produced by [`evaluate_stability`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StabilityReport {
    /// Number of region bootstrap resamples executed.
    pub n_bootstrap: usize,
    /// Deterministic PRNG seed used.
    pub seed: u64,
    /// Mean bootstrap cutoff.
    pub cutoff_mean: f64,
    /// 2.5th percentile bootstrap cutoff.
    pub cutoff_ci_low: f64,
    /// 97.5th percentile bootstrap cutoff.
    pub cutoff_ci_high: f64,
    /// Mean bootstrap super-enhancer count.
    pub se_count_mean: f64,
    /// Number of borderline regions (`0.10 <= call_frequency <= 0.90`).
    pub num_borderline_regions: usize,
    /// Per-region stability aligned 1:1 with the input `rank_signals` slice.
    pub per_region: Vec<RegionStability>,
    /// Optional sequencing-depth saturation curve (present when BAM read downsampling was collected).
    pub saturation_curve: Vec<DepthSaturationPoint>,
}

/// Deterministic 64-bit SplitMix64 generator for cross-platform reproducible bootstrapping.
#[derive(Debug, Clone, Copy)]
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn next_index(&mut self, len: usize) -> usize {
        (self.next_u64() % (len as u64)) as usize
    }
}

/// Runs `n_bootstrap` region resamples with replacement and evaluates optional read-downsampling
/// tiers (`tier_signals[tier_idx][region_idx]`).
#[must_use]
pub fn evaluate_stability(
    rank_signals: &[f64],
    n_bootstrap: usize,
    seed: u64,
    method: CutoffMethod,
    drop_top: usize,
    compat: CompatMode,
    downsampled_tiers: Option<&[(f64, Vec<f64>)]>,
) -> StabilityReport {
    let n = rank_signals.len();
    if n == 0 || n_bootstrap == 0 {
        return StabilityReport {
            n_bootstrap,
            seed,
            cutoff_mean: 0.0,
            cutoff_ci_low: 0.0,
            cutoff_ci_high: 0.0,
            se_count_mean: 0.0,
            num_borderline_regions: 0,
            per_region: vec![
                RegionStability {
                    call_frequency: 0.0,
                    is_borderline: false,
                };
                n
            ],
            saturation_curve: Vec::new(),
        };
    }

    let mut rng = SplitMix64::new(seed);
    let mut boot_cutoffs = Vec::with_capacity(n_bootstrap);
    let mut boot_counts = Vec::with_capacity(n_bootstrap);
    let mut sample_buf = vec![0.0_f64; n];
    let mut called_hits = vec![0_usize; n];

    for _ in 0..n_bootstrap {
        for slot in &mut sample_buf {
            *slot = rank_signals[rng.next_index(n)];
        }
        let sum = compute_cutoff_on_slice(&sample_buf, method, drop_top, compat);
        boot_cutoffs.push(sum.cutoff_signal);
        let mut se_cnt = 0_usize;
        for (idx, &sig) in rank_signals.iter().enumerate() {
            if sig > sum.cutoff_signal {
                called_hits[idx] += 1;
                se_cnt += 1;
            }
        }
        boot_counts.push(se_cnt as f64);
    }

    boot_cutoffs.sort_by(f64::total_cmp);
    let cutoff_mean = boot_cutoffs.iter().sum::<f64>() / (n_bootstrap as f64);
    let se_count_mean = boot_counts.iter().sum::<f64>() / (n_bootstrap as f64);
    let cutoff_ci_low = percentile_sorted(&boot_cutoffs, 0.025);
    let cutoff_ci_high = percentile_sorted(&boot_cutoffs, 0.975);

    let mut num_borderline = 0_usize;
    let per_region: Vec<RegionStability> = called_hits
        .into_iter()
        .map(|hits| {
            let freq = (hits as f64) / (n_bootstrap as f64);
            let is_borderline = (0.10..=0.90).contains(&freq);
            if is_borderline {
                num_borderline += 1;
            }
            RegionStability {
                call_frequency: freq,
                is_borderline,
            }
        })
        .collect();

    let full_sum = compute_cutoff_on_slice(rank_signals, method, drop_top, compat);
    let full_super_mask: Vec<bool> = rank_signals
        .iter()
        .map(|&s| s > full_sum.cutoff_signal)
        .collect();

    let mut saturation_curve = Vec::new();
    if let Some(tiers) = downsampled_tiers {
        for (frac, sigs) in tiers {
            let t_sum = compute_cutoff_on_slice(sigs, method, drop_top, compat);
            let mut inter = 0_usize;
            let mut union = 0_usize;
            let mut t_cnt = 0_usize;
            for (&s, &full_is_se) in sigs.iter().zip(full_super_mask.iter()) {
                let cur_is_se = s > t_sum.cutoff_signal;
                if cur_is_se {
                    t_cnt += 1;
                }
                if cur_is_se && full_is_se {
                    inter += 1;
                }
                if cur_is_se || full_is_se {
                    union += 1;
                }
            }
            let jaccard = if union == 0 {
                1.0
            } else {
                (inter as f64) / (union as f64)
            };
            saturation_curve.push(DepthSaturationPoint {
                read_fraction: *frac,
                cutoff_signal: t_sum.cutoff_signal,
                super_enhancer_count: t_cnt,
                jaccard_vs_full: jaccard,
            });
        }
    }

    StabilityReport {
        n_bootstrap,
        seed,
        cutoff_mean,
        cutoff_ci_low,
        cutoff_ci_high,
        se_count_mean,
        num_borderline_regions: num_borderline,
        per_region,
        saturation_curve,
    }
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let idx = p.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    let g = idx - (lo as f64);
    (1.0 - g) * sorted[lo] + g * sorted[hi]
}
