//! Super-enhancer tangent cutoff calling (`exact` chord scan and `rose-compat` Brent method),
//! LILY outlier options (`--drop-top`, `--min-se`), and ranking table assembly.

use serde::{Deserialize, Serialize};

use crate::types::{
    CompatMode, ControlSubtraction, CutoffMethod, RegionClass, RegionOrigin, StitchedRegion,
};

/// Configuration for super-enhancer cutoff calling (`rose call` / `rose run`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CallConfig {
    /// Tangent cutoff solver (`Exact` vs `RoseCompat`).
    pub method: CutoffMethod,
    /// Control background subtraction policy (`On` vs `Off`).
    pub subtract_control: ControlSubtraction,
    /// Exclude top `drop_top` extreme outlier regions when computing the slope/cutoff (`0` = off).
    pub drop_top: usize,
    /// Minimum number of super-enhancers enforced after calling (`0` = off, LILY uses `600`).
    pub min_super_enhancers: usize,
    /// Compatibility policy (`Clean` vs `RoseLegacy`).
    pub compat: CompatMode,
}

impl Default for CallConfig {
    fn default() -> Self {
        Self {
            method: CutoffMethod::Exact,
            subtract_control: ControlSubtraction::On,
            drop_top: 0,
            min_super_enhancers: 0,
            compat: CompatMode::Clean,
        }
    }
}

/// Measured signal entry for one stitched region entering [`call_super_enhancers`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RegionSignal {
    /// Stitched region metadata.
    pub region: StitchedRegion,
    /// Ranking track mean density (`rpm/bp`).
    pub treat_density: f64,
    /// Ranking track integrated signal (`rpm * bp`).
    pub treat_signal: f64,
    /// Optional control track mean density (`rpm/bp`).
    pub ctrl_density: Option<f64>,
    /// Optional control track integrated signal (`rpm * bp`).
    pub ctrl_signal: Option<f64>,
    /// Optional local copy-number ratio applied (`1.0` when uncorrected).
    pub cn_ratio: Option<f64>,
}

/// Ranked enhancer call row written to `*_AllEnhancers.table.txt` and `*_SuperEnhancers.table.txt`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RankedEnhancer {
    /// Underlying region and raw signals.
    pub item: RegionSignal,
    /// Raw background-subtracted ranking signal before copy-number division (`>= 0.0`).
    pub raw_net_signal: f64,
    /// Final ranking signal (`raw_net_signal / cn_ratio` when CN correction is enabled).
    pub rank_signal: f64,
    /// 1-based rank (`1` = strongest enhancer, `N` = weakest).
    pub enhancer_rank: usize,
    /// `true` if classified as a super-enhancer (`isSuper == 1`).
    pub is_super: bool,
    /// Regulatory class (`super_enhancer`, `enhancer`, `promoter`).
    pub class: RegionClass,
    /// BED column-5 score (`0..=1000` in `Clean` mode, raw `enhancer_rank` in `RoseLegacy`).
    pub bed_score: u32,
}

/// Summary of a cutoff call across a region set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CutoffSummary {
    /// Signal threshold separating super-enhancers (`signal > cutoff`) from typical enhancers.
    pub cutoff_signal: f64,
    /// 1-based ascending index `x_pt` (`1..=n`) of the tangent region on the sorted curve.
    pub tangent_index_1based: usize,
    /// Slope `s` of the sliding diagonal line.
    pub slope: f64,
    /// Ratio `cutoff_signal / median(signals)`.
    pub fold_over_median: f64,
    /// Ratio `cutoff_signal / mean(signals)`.
    pub fold_over_mean: f64,
    /// Total regions evaluated.
    pub num_regions: usize,
    /// Number of called super-enhancers.
    pub num_super_enhancers: usize,
    /// Number of regions with exact `0.0` ranking signal.
    pub num_zero_regions: usize,
}

/// Complete output of [`call_super_enhancers`].
#[derive(Debug, Clone, PartialEq)]
pub struct CallResult {
    /// Cutoff summary statistics.
    pub summary: CutoffSummary,
    /// Regions ordered by ascending `enhancer_rank` (`1` = top super-enhancer).
    pub ranked: Vec<RankedEnhancer>,
    /// Ascending sorted curve values `y[0..n]` used for hockey-stick plotting.
    pub sorted_signals_asc: Vec<f64>,
}

/// Computes the cutoff threshold on an arbitrary unsorted slice of non-negative signals.
#[must_use]
pub fn compute_cutoff_on_slice(
    raw_signals: &[f64],
    method: CutoffMethod,
    drop_top: usize,
    compat: CompatMode,
) -> CutoffSummary {
    if raw_signals.is_empty() {
        return CutoffSummary {
            cutoff_signal: 0.0,
            tangent_index_1based: 0,
            slope: 0.0,
            fold_over_median: 0.0,
            fold_over_mean: 0.0,
            num_regions: 0,
            num_super_enhancers: 0,
            num_zero_regions: 0,
        };
    }

    let mut sorted: Vec<f64> = raw_signals
        .iter()
        .map(|&v| if v.is_finite() && v > 0.0 { v } else { 0.0 })
        .collect();
    sorted.sort_by(f64::total_cmp);

    let n_total = sorted.len();
    let num_zero = sorted.iter().filter(|&&v| v == 0.0).count();
    let effective_len = n_total.saturating_sub(drop_top).max(1);
    let curve = &sorted[..effective_len];

    let (x_pt_1based, slope) = if curve.len() <= 1 || curve[curve.len() - 1] <= curve[0] {
        (curve.len(), 0.0)
    } else {
        match method {
            CutoffMethod::Exact => exact_tangent_point(curve, compat),
            CutoffMethod::RoseCompat => brent_rose_tangent_point(curve),
        }
    };

    let cutoff_signal = curve[x_pt_1based - 1];
    let num_super = sorted.iter().filter(|&&v| v > cutoff_signal).count();
    let mean = sorted.iter().sum::<f64>() / (n_total as f64);
    let median = if n_total % 2 == 1 {
        sorted[n_total / 2]
    } else {
        0.5 * (sorted[n_total / 2 - 1] + sorted[n_total / 2])
    };

    CutoffSummary {
        cutoff_signal,
        tangent_index_1based: x_pt_1based,
        slope,
        fold_over_median: if median > 0.0 {
            cutoff_signal / median
        } else {
            f64::INFINITY
        },
        fold_over_mean: if mean > 0.0 {
            cutoff_signal / mean
        } else {
            0.0
        },
        num_regions: n_total,
        num_super_enhancers: num_super,
        num_zero_regions: num_zero,
    }
}

/// Exact $O(N)$ linear scan for `argmin_{i in 1..=n} (y[i] - s * i)` with rightmost tie-breaking.
fn exact_tangent_point(sorted_asc: &[f64], compat: CompatMode) -> (usize, f64) {
    let n = sorted_asc.len();
    let denom = match compat {
        CompatMode::Clean => (n - 1).max(1) as f64,
        CompatMode::RoseLegacy => n as f64,
    };
    let slope = (sorted_asc[n - 1] - sorted_asc[0]) / denom;
    let mut best_i1 = 1_usize;
    let mut min_intercept = f64::INFINITY;
    let eps = 1e-12 * sorted_asc[n - 1].abs().max(1.0);

    for (idx0, &y) in sorted_asc.iter().enumerate() {
        let i1 = idx0 + 1;
        let b = y - slope * (i1 as f64);
        // Break ties (`<= min_intercept + eps`) toward the rightmost index so zero-signal
        // plateaus select the exact lift-off point.
        if b <= min_intercept + eps {
            min_intercept = min_intercept.min(b);
            best_i1 = i1;
        }
    }
    (best_i1, slope)
}

/// Port of R `ROSE_callSuper.R::calculate_cutoff` using Brent's `stats::optimize` on `[1, n]`.
fn brent_rose_tangent_point(sorted_asc: &[f64]) -> (usize, f64) {
    let n = sorted_asc.len();
    let slope = (sorted_asc[n - 1] - sorted_asc[0]) / (n as f64);

    let num_pts_below = |x: f64| -> f64 {
        let idx1 = (x.trunc() as usize).clamp(1, n);
        let y_pt = sorted_asc[idx1 - 1];
        // Note: R's `numPts_below_line` uses continuous `b <- yPt - (slope * x)`!
        let b = y_pt - slope * x;
        let mut count = 0_usize;
        for (k0, &yk) in sorted_asc.iter().enumerate() {
            let xk = (k0 + 1) as f64;
            if yk <= xk * slope + b {
                count += 1;
            }
        }
        count as f64
    };

    let mut a = 1.0_f64;
    let mut b = n as f64;
    let c_gold = 0.5 * (3.0 - 5.0_f64.sqrt());
    let eps = f64::EPSILON.sqrt();
    let tol = f64::EPSILON.powf(0.25); // R default `tol = .Machine$double.eps^0.25`

    let mut v = a + c_gold * (b - a);
    let mut w = v;
    let mut x = v;
    let mut d = 0.0_f64;
    let mut e = 0.0_f64;
    let mut fx = num_pts_below(x);
    let mut fv = fx;
    let mut fw = fx;

    loop {
        let xm = 0.5 * (a + b);
        let tol1 = eps * x.abs() + tol / 3.0;
        let tol2 = 2.0 * tol1;
        if (x - xm).abs() <= tol2 - 0.5 * (b - a) {
            break;
        }
        let mut use_golden = true;
        if e.abs() > tol1 {
            let mut r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            } else {
                q = -q;
            }
            r = e;
            e = d;
            if p.abs() < (0.5 * q * r).abs() && p > q * (a - x) && p < q * (b - x) {
                d = p / q;
                let u_try = x + d;
                if (u_try - a) < tol2 || (b - u_try) < tol2 {
                    d = if x < xm { tol1 } else { -tol1 };
                }
                use_golden = false;
            }
        }
        if use_golden {
            e = if x >= xm { a - x } else { b - x };
            d = c_gold * e;
        }
        let u = if d.abs() >= tol1 {
            x + d
        } else if d > 0.0 {
            x + tol1
        } else {
            x - tol1
        };
        let fu = num_pts_below(u);
        if fu <= fx {
            if u < x {
                b = x;
            } else {
                a = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || (w - x).abs() == 0.0 {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || (v - x).abs() == 0.0 || (v - w).abs() == 0.0 {
                v = u;
                fv = fu;
            }
        }
    }

    let x_pt = (x.floor() as usize).clamp(1, n);
    (x_pt, slope)
}

/// Ranks regions, computes the super-enhancer cutoff, applies optional `--min-se`, and builds
/// [`CallResult`].
#[must_use]
pub fn call_super_enhancers(items: &[RegionSignal], cfg: &CallConfig) -> CallResult {
    let n = items.len();
    let mut raw_nets = Vec::with_capacity(n);
    let mut rank_signals = Vec::with_capacity(n);

    for it in items {
        let ctrl = match cfg.subtract_control {
            ControlSubtraction::On => it.ctrl_signal.unwrap_or(0.0),
            ControlSubtraction::Off => 0.0,
        };
        let raw_net = (it.treat_signal - ctrl).max(0.0);
        let cn = it.cn_ratio.unwrap_or(1.0).max(1e-6);
        let corrected = raw_net / cn;
        raw_nets.push(raw_net);
        rank_signals.push(corrected);
    }

    let mut summary = compute_cutoff_on_slice(&rank_signals, cfg.method, cfg.drop_top, cfg.compat);

    // Match R's `nrow - rank(rankBy_vector, ties.method = "first") + 1`:
    // Ascending sort breaking ties by original row index `orig_idx` ascending -> `asc_rank` in `1..=n`,
    // then `enhancer_rank = n - asc_rank + 1`!
    let mut order_asc: Vec<usize> = (0..n).collect();
    order_asc.sort_by(|&i, &j| {
        rank_signals[i]
            .total_cmp(&rank_signals[j])
            .then_with(|| i.cmp(&j))
    });

    let sorted_signals_asc: Vec<f64> = order_asc.iter().map(|&idx| rank_signals[idx]).collect();

    if cfg.min_super_enhancers > 0 && n > 0 {
        let target_se = cfg.min_super_enhancers.min(n);
        if summary.num_super_enhancers < target_se {
            let cut_pos = n - target_se;
            summary.cutoff_signal = if cut_pos == 0 {
                -1.0
            } else {
                sorted_signals_asc[cut_pos - 1]
            };
            summary.num_super_enhancers = target_se;
        }
    }

    let mut asc_rank_of = vec![0_usize; n];
    for (pos0, &orig_idx) in order_asc.iter().enumerate() {
        asc_rank_of[orig_idx] = pos0 + 1;
    }

    let mut ranked: Vec<RankedEnhancer> = items
        .iter()
        .enumerate()
        .map(|(idx, item)| {
            let enh_rank = n - asc_rank_of[idx] + 1;
            let sig = rank_signals[idx];
            let is_super = if cfg.min_super_enhancers > 0
                && summary.num_super_enhancers == cfg.min_super_enhancers.min(n)
            {
                enh_rank <= summary.num_super_enhancers
            } else {
                sig > summary.cutoff_signal
            };
            let class = if item.region.origin == RegionOrigin::Promoter {
                RegionClass::Promoter
            } else if is_super {
                RegionClass::SuperEnhancer
            } else {
                RegionClass::Enhancer
            };
            let bed_score = match cfg.compat {
                CompatMode::RoseLegacy => u32::try_from(enh_rank).unwrap_or(u32::MAX),
                CompatMode::Clean => {
                    if n <= 1 {
                        1_000
                    } else {
                        let frac = ((n - enh_rank) as f64) / ((n - 1) as f64);
                        (frac * 1_000.0).round().clamp(0.0, 1_000.0) as u32
                    }
                }
            };
            RankedEnhancer {
                item: item.clone(),
                raw_net_signal: raw_nets[idx],
                rank_signal: sig,
                enhancer_rank: enh_rank,
                is_super,
                class,
                bed_score,
            }
        })
        .collect();

    ranked.sort_by_key(|r| r.enhancer_rank);

    CallResult {
        summary,
        ranked,
        sorted_signals_asc,
    }
}
