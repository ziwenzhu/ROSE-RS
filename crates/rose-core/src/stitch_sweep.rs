//! Zero-I/O piecewise stitch-distance scanning, sensitivity reporting (`0..=25 kb`), and
//! signal-blind inter-peak gap excess estimation (`--stitch auto`).

use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::cn::CopyNumberProfile;
use crate::cutoff::{call_super_enhancers, CallConfig, RegionSignal};
use crate::error::RoseCoreError;
use crate::signal::normalize_coverage;
use crate::stitch::{
    merge_intervals, prepare_constituents, stitch_prepared, sweep_stitch_chrom,
    PreparedConstituents, StitchConfig, StitchStats,
};
use crate::tss::AnnotationIndex;
use crate::types::{CompatMode, Constituent, Interval, RegionOrigin, StitchedRegion, Strand};

/// Default stitch distances (in bp) evaluated for the sensitivity report (`0..=25 kb`).
pub const DEFAULT_SENSITIVITY_DISTANCES: [u64; 10] = [
    0, 1_000, 2_500, 5_000, 7_500, 10_000, 12_500, 15_000, 20_000, 25_000,
];

/// Canonical reference stitch distance (`12,500 bp`) from Whyte et al. (2013).
pub const CANONICAL_STITCH_BP: u64 = 12_500;

/// User-selected stitch distance specification (`--stitch <bp|auto>`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StitchParam {
    /// Fixed distance in base pairs (default `12_500`).
    Fixed(u64),
    /// Opt-in geometric inter-peak gap excess selector (`--stitch auto`).
    Auto,
}

impl Default for StitchParam {
    fn default() -> Self {
        Self::Fixed(CANONICAL_STITCH_BP)
    }
}

impl fmt::Display for StitchParam {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fixed(bp) => write!(f, "{bp}"),
            Self::Auto => write!(f, "auto"),
        }
    }
}

impl FromStr for StitchParam {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if trimmed.eq_ignore_ascii_case("auto") {
            return Ok(Self::Auto);
        }
        let cleaned = trimmed.replace('_', "");
        cleaned.parse::<u64>().map(Self::Fixed).map_err(|_| {
            format!(
                "invalid --stitch value '{trimmed}'; expected non-negative bp integer or 'auto'"
            )
        })
    }
}

/// Guard-rail configuration for `--stitch auto`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AutoStitchConfig {
    /// Fallback distance when no clear clustering scale exists (`12_500`).
    pub fallback_bp: u64,
    /// Minimum allowed auto-selected distance (`2_000`).
    pub min_clamp_bp: u64,
    /// Maximum allowed auto-selected distance (`25_000`).
    pub max_clamp_bp: u64,
    /// Upper gap cutoff for active euchromatin (`1_000_000` excludes centromeric deserts).
    pub max_euchromatic_gap_bp: u64,
    /// Minimum euchromatic inter-peak gaps required before estimating $\hat{\lambda}_{\text{bg}}$ (`100`).
    pub min_gaps_required: usize,
    /// Minimum short-gap ECDF excess above Poisson baseline required (`0.05` = 5% of gaps).
    pub min_excess_fraction: f64,
}

impl Default for AutoStitchConfig {
    fn default() -> Self {
        Self {
            fallback_bp: CANONICAL_STITCH_BP,
            min_clamp_bp: 2_000,
            max_clamp_bp: 25_000,
            max_euchromatic_gap_bp: 1_000_000,
            min_gaps_required: 100,
            min_excess_fraction: 0.05,
        }
    }
}

/// Diagnostic summary produced by [`estimate_auto_stitch_distance`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AutoStitchReport {
    /// Resolved stitch distance in bp (`[min_clamp_bp, max_clamp_bp]` or `fallback_bp`).
    pub selected_stitch_bp: u64,
    /// Raw unclamped ECDF-excess argmax in bp (if estimated).
    pub raw_excess_argmax_bp: Option<u64>,
    /// `true` when guard rails triggered fallback to `fallback_bp`.
    pub used_fallback: bool,
    /// Human-readable explanation when `used_fallback` is `true` (or clamping note).
    pub note: String,
    /// Number of euchromatic inter-peak gaps (`0 < gap <= max_euchromatic_gap_bp`) analysed.
    pub num_euchromatic_gaps: usize,
    /// Estimated Poisson background gap rate $\hat{\lambda}_{\text{bg}}$ (gaps per bp).
    pub bg_rate_per_bp: Option<f64>,
    /// Maximum observed cumulative short-gap excess $\max_d [\hat{F}(d) - (1 - e^{-\hat{\lambda} d})]$.
    pub max_short_gap_excess: f64,
}

/// Estimates an opt-in stitch distance purely from inter-peak gap geometry (never inspecting
/// ChIP/CUT&RUN/bigWig signal, preventing circular ranking optimization).
///
/// # Errors
///
/// Returns [`RoseCoreError::DuplicateRegionId`] if `constituents` contain duplicate IDs.
pub fn estimate_auto_stitch_distance(
    constituents: &[Constituent],
    annotation: Option<&AnnotationIndex>,
    blacklist: &BTreeMap<String, Vec<Interval>>,
    stitch_cfg: &StitchConfig,
    auto_cfg: &AutoStitchConfig,
) -> Result<AutoStitchReport, RoseCoreError> {
    let prepared = prepare_constituents(constituents, annotation, blacklist, stitch_cfg)?;
    Ok(estimate_auto_from_prepared(
        &prepared,
        stitch_cfg.compat,
        auto_cfg,
    ))
}

fn estimate_auto_from_prepared(
    prepared: &PreparedConstituents,
    compat: CompatMode,
    auto_cfg: &AutoStitchConfig,
) -> AutoStitchReport {
    let mut gaps: Vec<u64> = Vec::new();
    for peaks in prepared.filtered_by_chrom.values() {
        let blocks = sweep_stitch_chrom(peaks, 0, compat);
        let mut prev_end: Option<u64> = None;
        for cl in blocks {
            let b_start = cl[0].interval.start();
            let b_end = cl
                .iter()
                .map(|c| c.interval.end())
                .max()
                .unwrap_or(b_start + 1);
            if let Some(pe) = prev_end {
                let gap = b_start.saturating_sub(pe);
                if gap > 0 && gap <= auto_cfg.max_euchromatic_gap_bp {
                    gaps.push(gap);
                }
            }
            prev_end = Some(b_end);
        }
    }

    gaps.sort_unstable();
    let m = gaps.len();
    if m < auto_cfg.min_gaps_required {
        return AutoStitchReport {
            selected_stitch_bp: auto_cfg.fallback_bp,
            raw_excess_argmax_bp: None,
            used_fallback: true,
            note: format!(
                "insufficient euchromatic inter-peak gaps ({m} < {}); fell back to {} bp",
                auto_cfg.min_gaps_required, auto_cfg.fallback_bp
            ),
            num_euchromatic_gaps: m,
            bg_rate_per_bp: None,
            max_short_gap_excess: 0.0,
        };
    }

    // Fit log-survival slope `ln S(g) = ln(1 - pi_clust) - lambda_bg * g` strictly in the
    // unclustered euchromatic tail `[lo_tail, hi_tail]` so even heavy short-gap clustering
    // (e.g. 75% clustered peaks) and upper-tail deserts cancel out cleanly.
    let hi_idx = ((m - 1) * 96) / 100;
    let hi_tail = gaps[hi_idx];
    let q65 = gaps[((m - 1) * 65) / 100];
    let q70 = gaps[((m - 1) * 70) / 100];
    let lo_tail = if hi_tail > 50_000 {
        25_000_u64.max(q65)
    } else {
        q70
    };

    let lo_pos = gaps.partition_point(|&g| g < lo_tail);
    let hi_pos = gaps.partition_point(|&g| g <= hi_tail).min(m - 1);

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let step = ((hi_pos.saturating_sub(lo_pos)) / 32).max(1);
    let mut idx = lo_pos;
    while idx <= hi_pos && idx < m {
        let surv = ((m - idx) as f64) / (m as f64);
        if surv > 1e-6 {
            xs.push(gaps[idx] as f64);
            ys.push(surv.ln());
        }
        idx += step;
    }

    if xs.len() < 4 {
        return AutoStitchReport {
            selected_stitch_bp: auto_cfg.fallback_bp,
            raw_excess_argmax_bp: None,
            used_fallback: true,
            note: format!(
                "insufficient tail span for background gap-rate fit; fell back to {} bp",
                auto_cfg.fallback_bp
            ),
            num_euchromatic_gaps: m,
            bg_rate_per_bp: None,
            max_short_gap_excess: 0.0,
        };
    }

    let n_pts = xs.len() as f64;
    let mean_x = xs.iter().sum::<f64>() / n_pts;
    let mean_y = ys.iter().sum::<f64>() / n_pts;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mean_x;
        num += dx * (y - mean_y);
        den += dx * dx;
    }

    let lambda_bg = if den > 0.0 {
        (-num / den).max(0.0)
    } else {
        0.0
    };
    if !lambda_bg.is_finite() || lambda_bg <= 0.0 {
        return AutoStitchReport {
            selected_stitch_bp: auto_cfg.fallback_bp,
            raw_excess_argmax_bp: None,
            used_fallback: true,
            note: format!(
                "degenerate tail gap variance; fell back to {} bp",
                auto_cfg.fallback_bp
            ),
            num_euchromatic_gaps: m,
            bg_rate_per_bp: None,
            max_short_gap_excess: 0.0,
        };
    }

    // Scan `d` in 100 bp increments up to `max(35_000, max_clamp_bp)` to find where empirical
    // short-gap density drops back to the Poisson null (`argmax [F_obs(d) - (1 - exp(-lambda * d))]`).
    let scan_hi = auto_cfg.max_clamp_bp.max(35_000);
    let mut best_d = auto_cfg.fallback_bp;
    let mut best_excess = f64::NEG_INFINITY;
    let m_f64 = m as f64;

    let mut d = 100_u64;
    while d <= scan_hi {
        let le_cnt = gaps.partition_point(|&g| g <= d);
        let f_obs = (le_cnt as f64) / m_f64;
        let f_null = 1.0 - (-lambda_bg * (d as f64)).exp();
        let excess = f_obs - f_null;
        if excess >= best_excess - 1e-12 {
            best_excess = excess;
            best_d = d;
        }
        d += 100;
    }

    // Guard rail: require both KS `p < 0.01` (`1.63 / sqrt(M)`) and practical excess `>= min_excess_fraction`.
    let ks_thresh = (1.63 / m_f64.sqrt()).max(auto_cfg.min_excess_fraction);
    if best_excess < ks_thresh {
        return AutoStitchReport {
            selected_stitch_bp: auto_cfg.fallback_bp,
            raw_excess_argmax_bp: Some(best_d),
            used_fallback: true,
            note: format!(
                "short-gap excess {:.3} below threshold {:.3} (random-like spacing); fell back to {} bp",
                best_excess, ks_thresh, auto_cfg.fallback_bp
            ),
            num_euchromatic_gaps: m,
            bg_rate_per_bp: Some(lambda_bg),
            max_short_gap_excess: best_excess.max(0.0),
        };
    }

    let clamped = best_d.clamp(auto_cfg.min_clamp_bp, auto_cfg.max_clamp_bp);
    let note = if clamped != best_d {
        format!(
            "selected {clamped} bp (clamped from raw gap-excess elbow {best_d} bp; excess={best_excess:.3})"
        )
    } else {
        format!("selected {clamped} bp from inter-peak gap excess (excess={best_excess:.3})")
    };

    AutoStitchReport {
        selected_stitch_bp: clamped,
        raw_excess_argmax_bp: Some(best_d),
        used_fallback: false,
        note,
        num_euchromatic_gaps: m,
        bg_rate_per_bp: Some(lambda_bg),
        max_short_gap_excess: best_excess,
    }
}

/// Elementary per-piece signal values captured during a single BAM/CRAM or BigWig scan.
#[derive(Debug, Clone)]
pub enum ElementaryTrackValues {
    /// Integer base-pair coverage sums from BAM/CRAM (`bp_sums` aligned 1:1 with `elementary_regions`).
    Bam {
        /// Per-piece integer base-pair coverage sum.
        bp_sums: Vec<u64>,
        /// Effective million-mapped-reads denominator.
        eff_mmr: f64,
        /// Optional explicit user scale factor overriding RPM.
        scale_factor: Option<f64>,
        /// Legacy vs clean rounding policy.
        compat: CompatMode,
        /// Optional downsampling fractions (`--stability`).
        tier_fractions: Vec<f64>,
        /// Optional `tier_fractions.len() x elementary_regions.len()` integer `bp_sums`.
        tier_bp_sums: Vec<Vec<u64>>,
    },
    /// Continuous width-weighted integrated signals (`rpm * bp`) from BigWig.
    BigWig {
        /// Per-piece integrated signal aligned 1:1 with `elementary_regions`.
        total_signals: Vec<f64>,
    },
}

/// Exact regions, signals, and coverage fractions materialized from [`PiecewiseStitchIndex`] at one
/// stitch distance.
#[derive(Debug, Clone)]
pub struct MaterializedStitch {
    /// Stitched (and multi-TSS-disaggregated) enhancer regions ordered by width descending.
    pub regions: Vec<StitchedRegion>,
    /// Filtering and disaggregation statistics at this stitch distance.
    pub stats: StitchStats,
    /// Quantified region signals aligned 1:1 with `regions`.
    pub signals: Vec<RegionSignal>,
    /// Optional per-tier downsampled total signals `(fraction, Vec<f64>)` (`--stability`).
    pub downsampled_totals: Vec<(f64, Vec<f64>)>,
    /// Share `(0.0..=1.0]` of each region's span covered by merged constituent peaks.
    pub peak_coverage_fracs: Vec<f64>,
}

/// One row of the stitch-distance sensitivity report (`*_Stitch_sensitivity.tsv`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StitchSensitivityRow {
    /// Evaluated stitch distance in bp.
    pub stitch_bp: u64,
    /// `true` if this row matches the active stitch distance used for primary outputs.
    pub is_active: bool,
    /// Total stitched (and multi-TSS-disaggregated) enhancer regions at `stitch_bp`.
    pub num_regions: usize,
    /// Number of super-enhancers called at `stitch_bp`.
    pub num_super_enhancers: usize,
    /// Tangent cutoff signal at `stitch_bp`.
    pub cutoff_signal: f64,
    /// Median width (in bp) of called super-enhancers (`0` if none).
    pub median_se_width_bp: u64,
    /// 90th-percentile width (in bp) of called super-enhancers (`0` if none).
    pub p90_se_width_bp: u64,
    /// Median fraction `(0.0..=1.0]` of each super-enhancer's span covered by constituent peaks.
    pub median_se_peak_coverage: f64,
    /// Median fraction `(0.0..=1.0]` of all regions' spans covered by constituent peaks.
    pub median_all_peak_coverage: f64,
    /// Exact base-pair Jaccard overlap against super-enhancers called at `12,500 bp`.
    pub base_jaccard_vs_12_5kb: f64,
    /// Fraction of `12,500 bp` super-enhancers overlapped by at least one super-enhancer at `stitch_bp`.
    pub se_overlap_frac_vs_12_5kb: f64,
}

#[derive(Debug, Clone)]
struct ChromPieceLayout {
    block_intervals: Vec<Interval>,
    block_elem_indices: Vec<usize>,
    gap_elem_indices: Vec<Option<usize>>,
    block_width_prefix: Vec<u64>,
}

/// Pre-decomposed peak + inter-peak gap index making any stitch distance `<= max_sweep_bp`
/// evaluable in milliseconds from a single BAM/CRAM/BigWig quantification pass.
#[derive(Debug, Clone)]
pub struct PiecewiseStitchIndex {
    prepared: PreparedConstituents,
    base_cfg: StitchConfig,
    max_sweep_bp: u64,
    elementary_regions: Vec<StitchedRegion>,
    chrom_layouts: BTreeMap<String, ChromPieceLayout>,
    single_constituent_elem: HashMap<(String, String), usize>,
    auto_report: Option<AutoStitchReport>,
}

impl PiecewiseStitchIndex {
    /// Decomposes `constituents` into strictly disjoint `d = 0` peak blocks and inter-peak gaps
    /// `<= max_sweep_bp` (plus any rare overlapping individual constituents for exact disaggregation),
    /// optionally resolving `--stitch auto` at the same time.
    ///
    /// # Errors
    ///
    /// Returns [`RoseCoreError`] if constituent IDs are duplicated.
    pub fn build(
        constituents: &[Constituent],
        annotation: Option<&AnnotationIndex>,
        blacklist: &BTreeMap<String, Vec<Interval>>,
        stitch_param: StitchParam,
        base_cfg: &StitchConfig,
        min_max_sweep_bp: u64,
    ) -> Result<(Self, u64), RoseCoreError> {
        let prepared = prepare_constituents(constituents, annotation, blacklist, base_cfg)?;
        let (active_stitch_bp, auto_report) = match stitch_param {
            StitchParam::Fixed(bp) => (bp, None),
            StitchParam::Auto => {
                let rep = estimate_auto_from_prepared(
                    &prepared,
                    base_cfg.compat,
                    &AutoStitchConfig::default(),
                );
                (rep.selected_stitch_bp, Some(rep))
            }
        };

        let max_sweep_bp = min_max_sweep_bp
            .max(active_stitch_bp)
            .max(CANONICAL_STITCH_BP);

        let mut elementary_regions: Vec<StitchedRegion> = Vec::new();
        let mut chrom_layouts: BTreeMap<String, ChromPieceLayout> = BTreeMap::new();
        let mut single_constituent_elem: HashMap<(String, String), usize> = HashMap::new();

        for (chrom, peaks) in &prepared.filtered_by_chrom {
            let clusters0 = sweep_stitch_chrom(peaks, 0, base_cfg.compat);
            let k = clusters0.len();
            let mut block_intervals = Vec::with_capacity(k);
            let mut block_elem_indices = Vec::with_capacity(k);
            let mut block_width_prefix = Vec::with_capacity(k + 1);
            block_width_prefix.push(0);

            for (b_i, cl) in clusters0.iter().enumerate() {
                let s = cl[0].interval.start();
                let e = cl.iter().map(|c| c.interval.end()).max().unwrap_or(s + 1);
                let iv = Interval::try_new(chrom, s, e)?;
                let elem_idx = elementary_regions.len();
                elementary_regions.push(StitchedRegion {
                    id: format!("__ELEM_PK_{chrom}_{b_i}"),
                    chrom: chrom.clone(),
                    interval: iv,
                    strand: Strand::Unstranded,
                    num_loci: 1,
                    constituent_size: iv.width(),
                    constituent_ids: Vec::new(),
                    origin: RegionOrigin::Enhancer,
                });
                block_intervals.push(iv);
                block_elem_indices.push(elem_idx);
                let prev_w = *block_width_prefix.last().unwrap_or(&0);
                block_width_prefix.push(prev_w + iv.width());

                if cl.len() == 1 {
                    single_constituent_elem.insert((chrom.clone(), cl[0].id.clone()), elem_idx);
                } else {
                    for single in cl {
                        let aux_idx = elementary_regions.len();
                        elementary_regions.push(StitchedRegion {
                            id: format!("__ELEM_SINGLE_{}_{}", chrom, single.id),
                            chrom: chrom.clone(),
                            interval: single.interval,
                            strand: single.strand,
                            num_loci: 1,
                            constituent_size: single.interval.width(),
                            constituent_ids: Vec::new(),
                            origin: RegionOrigin::Enhancer,
                        });
                        single_constituent_elem.insert((chrom.clone(), single.id.clone()), aux_idx);
                    }
                }
            }

            let mut gap_elem_indices = Vec::with_capacity(k.saturating_sub(1));
            for b_i in 0..k.saturating_sub(1) {
                let g_start = block_intervals[b_i].end();
                let g_end = block_intervals[b_i + 1].start();
                let gap_w = g_end.saturating_sub(g_start);
                if gap_w > 0 && gap_w <= max_sweep_bp {
                    let g_iv = Interval::try_new(chrom, g_start, g_end)?;
                    let g_idx = elementary_regions.len();
                    elementary_regions.push(StitchedRegion {
                        id: format!("__ELEM_GAP_{chrom}_{b_i}"),
                        chrom: chrom.clone(),
                        interval: g_iv,
                        strand: Strand::Unstranded,
                        num_loci: 0,
                        constituent_size: 0,
                        constituent_ids: Vec::new(),
                        origin: RegionOrigin::Enhancer,
                    });
                    gap_elem_indices.push(Some(g_idx));
                } else {
                    gap_elem_indices.push(None);
                }
            }

            chrom_layouts.insert(
                chrom.clone(),
                ChromPieceLayout {
                    block_intervals,
                    block_elem_indices,
                    gap_elem_indices,
                    block_width_prefix,
                },
            );
        }

        Ok((
            Self {
                prepared,
                base_cfg: base_cfg.clone(),
                max_sweep_bp,
                elementary_regions,
                chrom_layouts,
                single_constituent_elem,
                auto_report,
            },
            active_stitch_bp,
        ))
    }

    /// Elementary peak + gap regions to pass ONCE into BAM/CRAM/BigWig quantification.
    #[must_use]
    pub fn elementary_regions(&self) -> &[StitchedRegion] {
        &self.elementary_regions
    }

    /// Optional `--stitch auto` diagnostic report.
    #[must_use]
    pub const fn auto_report(&self) -> Option<&AutoStitchReport> {
        self.auto_report.as_ref()
    }

    /// Maximum stitch distance pre-indexed for zero-I/O scanning.
    #[must_use]
    pub const fn max_sweep_bp(&self) -> u64 {
        self.max_sweep_bp
    }

    /// Materializes exact [`StitchedRegion`]s, [`RegionSignal`]s, and optional downsampled tier
    /// totals at any `stitch_distance <= self.max_sweep_bp()` in milliseconds without touching I/O.
    ///
    /// # Errors
    ///
    /// Returns [`RoseCoreError`] if `stitch_distance > self.max_sweep_bp()`.
    pub fn materialize_at_distance(
        &self,
        stitch_distance: u64,
        annotation: Option<&AnnotationIndex>,
        treat_elem: &ElementaryTrackValues,
        ctrl_elem: Option<&ElementaryTrackValues>,
        cn_profile: Option<&CopyNumberProfile>,
    ) -> Result<MaterializedStitch, RoseCoreError> {
        if stitch_distance > self.max_sweep_bp {
            return Err(RoseCoreError::InvalidConfig(format!(
                "requested stitch_distance {stitch_distance} exceeds pre-indexed max_sweep_bp {}",
                self.max_sweep_bp
            )));
        }
        let mut cfg = self.base_cfg.clone();
        cfg.stitch_distance = stitch_distance;
        let (regions, stats) = stitch_prepared(&self.prepared, annotation, stitch_distance, &cfg)?;

        let mut signals = Vec::with_capacity(regions.len());
        let mut peak_cov_fracs = Vec::with_capacity(regions.len());

        let n_tiers = match treat_elem {
            ElementaryTrackValues::Bam { tier_fractions, .. } => tier_fractions.len(),
            ElementaryTrackValues::BigWig { .. } => 0,
        };
        let mut downsampled_vecs: Vec<Vec<f64>> = vec![Vec::with_capacity(regions.len()); n_tiers];

        for reg in &regions {
            let pieces = self.resolve_region_pieces(reg);
            peak_cov_fracs.push(pieces.peak_coverage_frac);

            let (t_den, t_sig) =
                eval_piece_sum(treat_elem, &pieces.elem_indices, reg.interval.width());
            let ctrl_pair =
                ctrl_elem.map(|ce| eval_piece_sum(ce, &pieces.elem_indices, reg.interval.width()));

            if let ElementaryTrackValues::Bam {
                eff_mmr,
                scale_factor,
                compat,
                tier_fractions,
                tier_bp_sums,
                ..
            } = treat_elem
            {
                for (t_i, &frac) in tier_fractions.iter().enumerate() {
                    let bp: u64 = pieces
                        .elem_indices
                        .iter()
                        .map(|&ei| tier_bp_sums[t_i][ei])
                        .sum();
                    let tier_mmr = (*eff_mmr * frac).max(1e-12);
                    let tier_sf = scale_factor.map(|sf| sf / frac.max(1e-6));
                    let sig =
                        normalize_coverage(bp, reg.interval.width(), tier_mmr, tier_sf, *compat).1;
                    downsampled_vecs[t_i].push(sig);
                }
            }

            signals.push(RegionSignal {
                region: reg.clone(),
                treat_density: t_den,
                treat_signal: t_sig,
                ctrl_density: ctrl_pair.map(|p| p.0),
                ctrl_signal: ctrl_pair.map(|p| p.1),
                cn_ratio: cn_profile.map(|p| p.ratio_for_region(reg)),
            });
        }

        let downsampled = match treat_elem {
            ElementaryTrackValues::Bam { tier_fractions, .. } => tier_fractions
                .iter()
                .copied()
                .zip(downsampled_vecs)
                .collect(),
            ElementaryTrackValues::BigWig { .. } => Vec::new(),
        };

        Ok(MaterializedStitch {
            regions,
            stats,
            signals,
            downsampled_totals: downsampled,
            peak_coverage_fracs: peak_cov_fracs,
        })
    }

    /// Evaluates the stitch-distance sensitivity sweep across `distances` (plus `active_distance`
    /// and `12,500 bp` reference) with zero additional BAM/BigWig I/O.
    ///
    /// # Errors
    ///
    /// Returns [`RoseCoreError`] if any distance exceeds `self.max_sweep_bp()`.
    #[expect(
        clippy::too_many_arguments,
        reason = "sensitivity sweep combines distance grid, optional annotation, optional control/CN tracks, and call config"
    )]
    pub fn evaluate_sensitivity_sweep(
        &self,
        distances: &[u64],
        active_distance: u64,
        annotation: Option<&AnnotationIndex>,
        treat_elem: &ElementaryTrackValues,
        ctrl_elem: Option<&ElementaryTrackValues>,
        cn_profile: Option<&CopyNumberProfile>,
        call_cfg: &CallConfig,
    ) -> Result<Vec<StitchSensitivityRow>, RoseCoreError> {
        let mut grid: Vec<u64> = distances
            .iter()
            .copied()
            .filter(|&d| d <= self.max_sweep_bp)
            .collect();
        if active_distance <= self.max_sweep_bp {
            grid.push(active_distance);
        }
        if CANONICAL_STITCH_BP <= self.max_sweep_bp {
            grid.push(CANONICAL_STITCH_BP);
        }
        grid.sort_unstable();
        grid.dedup();

        // Build `12.5 kb` reference super-enhancer intervals first.
        let ref_mat = self.materialize_at_distance(
            CANONICAL_STITCH_BP,
            annotation,
            treat_elem,
            ctrl_elem,
            cn_profile,
        )?;
        let ref_call = call_super_enhancers(&ref_mat.signals, call_cfg);
        let ref_se_raw: Vec<(String, Interval)> = ref_call
            .ranked
            .iter()
            .filter(|r| r.is_super)
            .map(|r| (r.item.region.chrom.clone(), r.item.region.interval))
            .collect();
        let ref_se_merged = group_and_merge_by_chrom(&ref_se_raw);

        let mut out = Vec::with_capacity(grid.len());
        for d in grid {
            let mut cur_mat =
                self.materialize_at_distance(d, annotation, treat_elem, ctrl_elem, cn_profile)?;
            let call = call_super_enhancers(&cur_mat.signals, call_cfg);

            let id_to_cov: HashMap<&str, f64> = cur_mat
                .regions
                .iter()
                .zip(cur_mat.peak_coverage_fracs.iter())
                .map(|(r, &cf)| (r.id.as_str(), cf))
                .collect();

            let mut se_widths: Vec<u64> = Vec::new();
            let mut se_covs: Vec<f64> = Vec::new();
            let mut cur_se_raw: Vec<(String, Interval)> = Vec::new();

            for rk in &call.ranked {
                if rk.is_super {
                    let iv = rk.item.region.interval;
                    se_widths.push(iv.width());
                    if let Some(&cf) = id_to_cov.get(rk.item.region.id.as_str()) {
                        se_covs.push(cf);
                    }
                    cur_se_raw.push((rk.item.region.chrom.clone(), iv));
                }
            }

            let cur_se_merged = group_and_merge_by_chrom(&cur_se_raw);
            let base_jaccard = bp_jaccard_chrom_maps(&cur_se_merged, &ref_se_merged);
            let se_recall = if ref_se_raw.is_empty() {
                if cur_se_raw.is_empty() {
                    1.0
                } else {
                    0.0
                }
            } else {
                let hits = ref_se_raw
                    .iter()
                    .filter(|(chr, iv)| {
                        cur_se_merged.get(chr).is_some_and(|list| {
                            let up = list.partition_point(|x| x.start() < iv.end());
                            list[..up].iter().any(|x| x.overlaps(*iv))
                        })
                    })
                    .count();
                (hits as f64) / (ref_se_raw.len() as f64)
            };

            out.push(StitchSensitivityRow {
                stitch_bp: d,
                is_active: d == active_distance,
                num_regions: call.summary.num_regions,
                num_super_enhancers: call.summary.num_super_enhancers,
                cutoff_signal: call.summary.cutoff_signal,
                median_se_width_bp: quantile_u64(&mut se_widths, 0.50),
                p90_se_width_bp: quantile_u64(&mut se_widths, 0.90),
                median_se_peak_coverage: quantile_f64(&mut se_covs, 0.50),
                median_all_peak_coverage: quantile_f64(&mut cur_mat.peak_coverage_fracs, 0.50),
                base_jaccard_vs_12_5kb: base_jaccard,
                se_overlap_frac_vs_12_5kb: se_recall,
            });
        }
        Ok(out)
    }

    fn resolve_region_pieces(&self, reg: &StitchedRegion) -> RegionPieces {
        // Disaggregated single constituent (`id` does not end with `_lociStitched`):
        if reg.num_loci == 1 && !reg.id.ends_with("_lociStitched") {
            if let Some(&single_idx) = self
                .single_constituent_elem
                .get(&(reg.chrom.clone(), reg.id.clone()))
            {
                return RegionPieces {
                    elem_indices: vec![single_idx],
                    peak_coverage_frac: 1.0,
                };
            }
        }

        let Some(layout) = self.chrom_layouts.get(&reg.chrom) else {
            return RegionPieces {
                elem_indices: Vec::new(),
                peak_coverage_frac: 0.0,
            };
        };

        let lo = layout
            .block_intervals
            .partition_point(|b| b.start() < reg.interval.start());
        let hi_excl = layout
            .block_intervals
            .partition_point(|b| b.end() <= reg.interval.end());
        if lo >= hi_excl || hi_excl > layout.block_intervals.len() {
            return RegionPieces {
                elem_indices: Vec::new(),
                peak_coverage_frac: 0.0,
            };
        }

        let mut elem_indices = Vec::with_capacity((hi_excl - lo) * 2);
        for k in lo..hi_excl {
            elem_indices.push(layout.block_elem_indices[k]);
            if k + 1 < hi_excl {
                if let Some(g_idx) = layout.gap_elem_indices[k] {
                    elem_indices.push(g_idx);
                }
            }
        }

        let peak_bp = layout.block_width_prefix[hi_excl] - layout.block_width_prefix[lo];
        let frac = ((peak_bp as f64) / (reg.interval.width().max(1) as f64)).clamp(0.0, 1.0);
        RegionPieces {
            elem_indices,
            peak_coverage_frac: frac,
        }
    }
}

struct RegionPieces {
    elem_indices: Vec<usize>,
    peak_coverage_frac: f64,
}

fn eval_piece_sum(
    track: &ElementaryTrackValues,
    indices: &[usize],
    region_width: u64,
) -> (f64, f64) {
    match track {
        ElementaryTrackValues::Bam {
            bp_sums,
            eff_mmr,
            scale_factor,
            compat,
            ..
        } => {
            let bp: u64 = indices.iter().map(|&i| bp_sums[i]).sum();
            normalize_coverage(bp, region_width, *eff_mmr, *scale_factor, *compat)
        }
        ElementaryTrackValues::BigWig { total_signals } => {
            let sig: f64 = indices.iter().map(|&i| total_signals[i]).sum();
            let den = sig / (region_width.max(1) as f64);
            (den, sig)
        }
    }
}

fn group_and_merge_by_chrom(items: &[(String, Interval)]) -> BTreeMap<String, Vec<Interval>> {
    let mut map: BTreeMap<String, Vec<Interval>> = BTreeMap::new();
    for (chr, iv) in items {
        map.entry(chr.clone()).or_default().push(*iv);
    }
    for list in map.values_mut() {
        *list = merge_intervals(list);
    }
    map
}

fn bp_jaccard_chrom_maps(
    a: &BTreeMap<String, Vec<Interval>>,
    b: &BTreeMap<String, Vec<Interval>>,
) -> f64 {
    let total_a: u64 = a.values().flat_map(|v| v.iter()).map(|iv| iv.width()).sum();
    let total_b: u64 = b.values().flat_map(|v| v.iter()).map(|iv| iv.width()).sum();
    if total_a == 0 && total_b == 0 {
        return 1.0;
    }

    let mut inter_bp = 0_u64;
    for (chr, list_a) in a {
        if let Some(list_b) = b.get(chr) {
            let mut i = 0;
            let mut j = 0;
            while i < list_a.len() && j < list_b.len() {
                inter_bp += list_a[i].overlap_len(list_b[j]);
                if list_a[i].end() < list_b[j].end() {
                    i += 1;
                } else {
                    j += 1;
                }
            }
        }
    }
    let union_bp = (total_a + total_b).saturating_sub(inter_bp);
    if union_bp == 0 {
        1.0
    } else {
        (inter_bp as f64) / (union_bp as f64)
    }
}

fn quantile_u64(vals: &mut [u64], q: f64) -> u64 {
    if vals.is_empty() {
        return 0;
    }
    vals.sort_unstable();
    let idx = (((vals.len() - 1) as f64) * q.clamp(0.0, 1.0)).round() as usize;
    vals[idx.min(vals.len() - 1)]
}

fn quantile_f64(vals: &mut [f64], q: f64) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.sort_by(f64::total_cmp);
    let idx = (((vals.len() - 1) as f64) * q.clamp(0.0, 1.0)).round() as usize;
    vals[idx.min(vals.len() - 1)]
}
