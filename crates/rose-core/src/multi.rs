//! Multi-sample and replicate cohort engine (`rose multi`): DESeq2-style median-of-ratios
//! normalisation, per-sample calling, consensus super-enhancers, and two-group fold-change labels.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::cutoff::{call_super_enhancers, CallConfig, CallResult, RegionSignal};
use crate::error::RoseCoreError;
use crate::types::StitchedRegion;

/// Metadata for one sample in a `rose multi` cohort sheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SampleSpec {
    /// Unique sample identifier.
    pub sample_id: String,
    /// Biological condition / group label (e.g. `"Treated"`, `"DMSO"`).
    pub group: String,
    /// Optional explicit spike-in or external scale factor overriding median-of-ratios.
    pub scale_factor: Option<f64>,
}

/// Consensus status for a shared region across groups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DiffLabel {
    /// Consensus super-enhancer in group A (`target`) and NOT in group B (`reference`).
    Gained,
    /// Consensus super-enhancer in group B (`reference`) and NOT in group A (`target`).
    Lost,
    /// Consensus super-enhancer in both groups.
    SharedSuper,
    /// Not a consensus super-enhancer in either group.
    Typical,
}

/// Group-level consensus and comparison record for one shared stitched region.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MultiRegionSummary {
    /// Shared stitched region.
    pub region: StitchedRegion,
    /// Normalized ranking signal per sample (`len == samples.len()`).
    pub sample_signals: Vec<f64>,
    /// Per-sample super-enhancer call (`true` = called in sample `j`).
    pub sample_calls: Vec<bool>,
    /// Per-group count of replicates calling this region as a super-enhancer.
    pub group_call_counts: BTreeMap<String, usize>,
    /// Per-group consensus super-enhancer boolean (`calls >= min_replicates`).
    pub group_consensus: BTreeMap<String, bool>,
    /// Optional two-group comparison label (`group_b` vs `group_a` when exactly 2 groups exist).
    pub diff_label: Option<DiffLabel>,
    /// Optional `log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))`.
    pub log2_fold_change: Option<f64>,
}

/// Full result of [`evaluate_multi_cohort`].
#[derive(Debug, Clone, PartialEq)]
pub struct MultiCohortResult {
    /// Sample descriptors.
    pub samples: Vec<SampleSpec>,
    /// Inter-sample size factors applied (`1.0` when unchanged).
    pub size_factors: Vec<f64>,
    /// Per-sample [`CallResult`] on the shared region set.
    pub per_sample_calls: Vec<CallResult>,
    /// Shared region matrix, consensus calls, and optional group comparison.
    pub regions: Vec<MultiRegionSummary>,
    /// Ordered unique group names (`[group1, group2, ...]`).
    pub group_order: Vec<String>,
    /// Pseudocount added when computing `log2_fold_change` (`1.0` signal unit).
    pub pseudocount: f64,
}

/// Computes inter-sample DESeq2 median-of-ratios size factors over regions with strictly positive
/// signal in at least `min_positive_fraction` (e.g. `0.80`) of samples.
#[must_use]
pub fn compute_median_ratio_size_factors(
    raw_matrix: &[Vec<f64>],
    n_samples: usize,
    min_positive_fraction: f64,
) -> Option<Vec<f64>> {
    if raw_matrix.is_empty() || n_samples == 0 {
        return None;
    }
    let min_pos_count = ((n_samples as f64) * min_positive_fraction).ceil() as usize;
    let min_pos_count = min_pos_count.clamp(1, n_samples);

    let mut pseudo_refs: Vec<(usize, f64)> = Vec::new();
    for (r_idx, row) in raw_matrix.iter().enumerate() {
        let mut log_sum = 0.0_f64;
        let mut pos_cnt = 0_usize;
        for &val in row {
            if val.is_finite() && val > 0.0 {
                log_sum += val.ln();
                pos_cnt += 1;
            }
        }
        if pos_cnt >= min_pos_count {
            let geom_mean = (log_sum / (pos_cnt as f64)).exp();
            if geom_mean > 0.0 && geom_mean.is_finite() {
                pseudo_refs.push((r_idx, geom_mean));
            }
        }
    }

    // Require at least 10 shared positive regions for stable median-of-ratios scaling.
    if pseudo_refs.len() < 10 {
        return None;
    }

    let mut factors = Vec::with_capacity(n_samples);
    for s_idx in 0..n_samples {
        let mut ratios: Vec<f64> = pseudo_refs
            .iter()
            .filter_map(|&(r_idx, g_mean)| {
                let v = *raw_matrix[r_idx].get(s_idx)?;
                (v.is_finite() && v > 0.0).then_some(v / g_mean)
            })
            .collect();
        if ratios.is_empty() {
            return None;
        }
        ratios.sort_by(f64::total_cmp);
        let mid = ratios.len() / 2;
        let med = if ratios.len() % 2 == 1 {
            ratios[mid]
        } else {
            0.5 * (ratios[mid - 1] + ratios[mid])
        };
        if med <= 0.0 || !med.is_finite() {
            return None;
        }
        factors.push(med);
    }
    Some(factors)
}

/// Evaluates a multi-sample cohort across shared `regions`, applying inter-sample normalisation,
/// per-sample super-enhancer calling, group replicate consensus, and 2-group comparison.
///
/// # Errors
///
/// Returns [`RoseCoreError::InvalidConfig`] on dimension mismatches between `samples` and
/// `per_sample_signals`.
pub fn evaluate_multi_cohort(
    regions: &[StitchedRegion],
    samples: &[SampleSpec],
    per_sample_signals: &[Vec<RegionSignal>],
    call_cfg: &CallConfig,
    min_replicates: Option<usize>,
) -> Result<MultiCohortResult, RoseCoreError> {
    let n_samples = samples.len();
    let n_regions = regions.len();
    if n_samples == 0 || per_sample_signals.len() != n_samples {
        return Err(RoseCoreError::InvalidConfig(
            "samples count must be >= 1 and match per_sample_signals length".to_owned(),
        ));
    }
    if per_sample_signals.iter().any(|v| v.len() != n_regions) {
        return Err(RoseCoreError::InvalidConfig(
            "every sample signal vector must have length equal to shared regions".to_owned(),
        ));
    }

    // Build initial raw net signal matrix `[region_idx][sample_idx]`
    let mut raw_matrix = vec![vec![0.0_f64; n_samples]; n_regions];
    for s_idx in 0..n_samples {
        for r_idx in 0..n_regions {
            let it = &per_sample_signals[s_idx][r_idx];
            let ctrl = it.ctrl_signal.unwrap_or(0.0);
            raw_matrix[r_idx][s_idx] = (it.treat_signal - ctrl).max(0.0);
        }
    }

    let all_have_explicit_sf = samples.iter().all(|s| s.scale_factor.is_some());
    let size_factors: Vec<f64> = if all_have_explicit_sf {
        vec![1.0; n_samples]
    } else {
        compute_median_ratio_size_factors(&raw_matrix, n_samples, 0.80)
            .unwrap_or_else(|| vec![1.0; n_samples])
    };

    let mut per_sample_calls = Vec::with_capacity(n_samples);
    let mut scaled_matrix = vec![vec![0.0_f64; n_samples]; n_regions];
    let mut call_matrix = vec![vec![false; n_samples]; n_regions];

    for s_idx in 0..n_samples {
        let sf = size_factors[s_idx];
        let adjusted: Vec<RegionSignal> = per_sample_signals[s_idx]
            .iter()
            .map(|orig| {
                let mut item = orig.clone();
                if samples[s_idx].scale_factor.is_none() && (sf - 1.0).abs() > f64::EPSILON {
                    item.treat_density /= sf;
                    item.treat_signal /= sf;
                    if let Some(cd) = &mut item.ctrl_density {
                        *cd /= sf;
                    }
                    if let Some(cs) = &mut item.ctrl_signal {
                        *cs /= sf;
                    }
                }
                item
            })
            .collect();

        let res = call_super_enhancers(&adjusted, call_cfg);
        // Map back to shared `regions` index by region ID
        let id_to_pos: BTreeMap<&str, usize> = regions
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.id.as_str(), idx))
            .collect();
        for rk in &res.ranked {
            if let Some(&r_idx) = id_to_pos.get(rk.item.region.id.as_str()) {
                scaled_matrix[r_idx][s_idx] = rk.rank_signal;
                call_matrix[r_idx][s_idx] = rk.is_super;
            }
        }
        per_sample_calls.push(res);
    }

    let mut group_order: Vec<String> = Vec::new();
    let mut group_members: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (s_idx, s) in samples.iter().enumerate() {
        if !group_members.contains_key(&s.group) {
            group_order.push(s.group.clone());
        }
        group_members
            .entry(s.group.clone())
            .or_default()
            .push(s_idx);
    }

    let pseudocount = 1.0_f64;
    let mut out_regions = Vec::with_capacity(n_regions);
    for (r_idx, reg) in regions.iter().enumerate() {
        let mut group_call_counts = BTreeMap::new();
        let mut group_consensus = BTreeMap::new();
        for g_name in &group_order {
            let members = &group_members[g_name];
            let c_cnt = members
                .iter()
                .filter(|&&s_idx| call_matrix[r_idx][s_idx])
                .count();
            let req = min_replicates.map_or(members.len(), |m| m.min(members.len()).max(1));
            group_call_counts.insert(g_name.clone(), c_cnt);
            group_consensus.insert(g_name.clone(), c_cnt >= req);
        }

        let (diff_label, log2_fc) = if group_order.len() == 2 {
            let g_ref = &group_order[0];
            let g_tgt = &group_order[1];
            let is_ref = group_consensus[g_ref];
            let is_tgt = group_consensus[g_tgt];
            let lbl = match (is_ref, is_tgt) {
                (false, true) => DiffLabel::Gained,
                (true, false) => DiffLabel::Lost,
                (true, true) => DiffLabel::SharedSuper,
                (false, false) => DiffLabel::Typical,
            };
            let mean_ref = mean_for_indices(&scaled_matrix[r_idx], &group_members[g_ref]);
            let mean_tgt = mean_for_indices(&scaled_matrix[r_idx], &group_members[g_tgt]);
            let fc = ((mean_tgt + pseudocount) / (mean_ref + pseudocount)).log2();
            (Some(lbl), Some(fc))
        } else {
            (None, None)
        };

        out_regions.push(MultiRegionSummary {
            region: reg.clone(),
            sample_signals: scaled_matrix[r_idx].clone(),
            sample_calls: call_matrix[r_idx].clone(),
            group_call_counts,
            group_consensus,
            diff_label,
            log2_fold_change: log2_fc,
        });
    }

    Ok(MultiCohortResult {
        samples: samples.to_vec(),
        size_factors,
        per_sample_calls,
        regions: out_regions,
        group_order,
        pseudocount,
    })
}

fn mean_for_indices(row: &[f64], indices: &[usize]) -> f64 {
    if indices.is_empty() {
        return 0.0;
    }
    let sum: f64 = indices.iter().map(|&i| row[i]).sum();
    sum / (indices.len() as f64)
}
