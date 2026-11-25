//! Constituent enhancer filtering, TSS promoter handling, 12.5 kb stitching,
//! `--max-tss-genes` disaggregation, and multi-sample consensus peak merging.

use std::collections::{BTreeMap, HashSet};

use crate::error::RoseCoreError;
use crate::tss::AnnotationIndex;
use crate::types::{
    CompatMode, Constituent, Interval, PromoterMode, RegionOrigin, StitchedRegion, Strand,
};

/// Parameters controlling enhancer stitching and TSS exclusion (`rose stitch` / `rose run`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StitchConfig {
    /// Maximum gap in bp bridged between adjacent constituents (`12_500` by default).
    pub stitch_distance: u64,
    /// Symmetric half-window around TSSs excluded before stitching (`0` = disabled, `2_500` typical).
    pub tss_window: u64,
    /// Promoter handling mode (`Contained` drop vs `Subtract` interval carving).
    pub promoter_mode: PromoterMode,
    /// Maximum distinct gene symbols allowed inside a stitched enhancer before disaggregation (`2`).
    pub max_tss_genes: usize,
    /// Chromosomes excluded from output (`vec!["chrY".into()]` by default).
    pub exclude_chroms: Vec<String>,
    /// Compatibility policy (`Clean` vs `RoseLegacy`).
    pub compat: CompatMode,
}

impl Default for StitchConfig {
    fn default() -> Self {
        Self {
            stitch_distance: 12_500,
            tss_window: 0,
            promoter_mode: PromoterMode::Contained,
            max_tss_genes: 2,
            exclude_chroms: vec!["chrY".to_owned(), "Y".to_owned()],
            compat: CompatMode::Clean,
        }
    }
}

/// Summary statistics recorded during stitching (`run.json` & diagnostics).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StitchStats {
    /// Input constituent count.
    pub input_constituents: usize,
    /// Constituents removed by `--blacklist`.
    pub removed_blacklisted: usize,
    /// Constituents removed or carved by `--tss` promoter handling.
    pub removed_tss_contained: usize,
    /// Stitched clusters disaggregated due to spanning `> max_tss_genes` distinct gene TSSs.
    pub disaggregated_clusters: usize,
    /// Constituents restored by disaggregation.
    pub restored_constituents: usize,
    /// Regions removed by `--exclude-chrom`.
    pub removed_excluded_chroms: usize,
}

/// Pre-filtered, per-chromosome sorted constituents shared between single-distance stitching
/// and multi-distance `PiecewiseStitchIndex` sweeps.
#[derive(Debug, Clone)]
pub(crate) struct PreparedConstituents {
    pub stats: StitchStats,
    pub ref_by_chrom: BTreeMap<String, Vec<Interval>>,
    pub filtered_by_chrom: BTreeMap<String, Vec<Constituent>>,
}

pub(crate) fn prepare_constituents(
    constituents: &[Constituent],
    annotation: Option<&AnnotationIndex>,
    blacklist: &BTreeMap<String, Vec<Interval>>,
    cfg: &StitchConfig,
) -> Result<PreparedConstituents, RoseCoreError> {
    let mut seen_ids = HashSet::with_capacity(constituents.len());
    for c in constituents {
        if !seen_ids.insert(c.id.as_str()) {
            return Err(RoseCoreError::DuplicateRegionId { id: c.id.clone() });
        }
    }

    let mut stats = StitchStats {
        input_constituents: constituents.len(),
        ..StitchStats::default()
    };

    let excluded_set: HashSet<&str> = cfg.exclude_chroms.iter().map(String::as_str).collect();

    // Index original input constituents per chromosome for exact `CONSTITUENT_SIZE` computation
    // matching `ROSE_main.py:169-172` (`referenceCollection.getOverlap`).
    let mut ref_by_chrom: BTreeMap<String, Vec<Interval>> = BTreeMap::new();
    let mut filtered_by_chrom: BTreeMap<String, Vec<Constituent>> = BTreeMap::new();

    for c in constituents {
        if excluded_set.contains(c.chrom.as_str()) {
            stats.removed_excluded_chroms += 1;
            continue;
        }
        if let Some(bl_list) = blacklist.get(&c.chrom) {
            if overlaps_any_sorted(bl_list, c.interval) {
                stats.removed_blacklisted += 1;
                continue;
            }
        }
        ref_by_chrom
            .entry(c.chrom.clone())
            .or_default()
            .push(c.interval);

        if cfg.tss_window > 0 {
            if let Some(annot) = annotation {
                let promo_wins = collect_promoter_windows(annot, &c.chrom, c.interval, cfg);
                match cfg.promoter_mode {
                    PromoterMode::Contained => {
                        if promo_wins.iter().any(|w| w.contains(c.interval)) {
                            stats.removed_tss_contained += 1;
                            continue;
                        }
                    }
                    PromoterMode::Subtract => {
                        let pieces = subtract_many(c.interval, &promo_wins);
                        if pieces.is_empty() {
                            stats.removed_tss_contained += 1;
                            continue;
                        }
                        if pieces.len() == 1 && pieces[0] == c.interval {
                            filtered_by_chrom
                                .entry(c.chrom.clone())
                                .or_default()
                                .push(c.clone());
                        } else {
                            for (p_idx, piece) in pieces.into_iter().enumerate() {
                                filtered_by_chrom.entry(c.chrom.clone()).or_default().push(
                                    Constituent {
                                        id: format!("{}_part{}", c.id, p_idx + 1),
                                        chrom: c.chrom.clone(),
                                        interval: piece,
                                        strand: c.strand,
                                        score: c.score,
                                    },
                                );
                            }
                        }
                        continue;
                    }
                }
            }
        }

        filtered_by_chrom
            .entry(c.chrom.clone())
            .or_default()
            .push(c.clone());
    }

    for list in ref_by_chrom.values_mut() {
        list.sort();
    }
    for peaks in filtered_by_chrom.values_mut() {
        peaks.sort_by(|a, b| {
            a.interval
                .start()
                .cmp(&b.interval.start())
                .then_with(|| a.interval.end().cmp(&b.interval.end()))
                .then_with(|| a.id.cmp(&b.id))
        });
    }

    Ok(PreparedConstituents {
        stats,
        ref_by_chrom,
        filtered_by_chrom,
    })
}

pub(crate) fn stitch_prepared(
    prepared: &PreparedConstituents,
    annotation: Option<&AnnotationIndex>,
    stitch_distance: u64,
    cfg: &StitchConfig,
) -> Result<(Vec<StitchedRegion>, StitchStats), RoseCoreError> {
    let mut stats = prepared.stats.clone();
    let mut out_regions: Vec<StitchedRegion> = Vec::new();
    for (chrom, peaks) in &prepared.filtered_by_chrom {
        let clusters = sweep_stitch_chrom(peaks, stitch_distance, cfg.compat);
        let ref_slice = prepared
            .ref_by_chrom
            .get(chrom)
            .map_or(&[][..], Vec::as_slice);

        for cluster in clusters {
            let merged_start = cluster[0].interval.start();
            let merged_end = cluster
                .iter()
                .map(|c| c.interval.end())
                .max()
                .unwrap_or(merged_start + 1);
            let merged_iv = Interval::try_new(chrom, merged_start, merged_end)?;

            let should_disaggregate = if cfg.tss_window > 0 {
                if let Some(annot) = annotation {
                    let (lo, hi) = match cfg.compat {
                        CompatMode::Clean => (merged_iv.start(), merged_iv.end().saturating_sub(1)),
                        // Legacy `ROSE_main.py:101` uses `makeTSSLocus(..., 50, 50)` against
                        // `[g_start, g_end]`.
                        CompatMode::RoseLegacy => (
                            merged_iv.start().saturating_sub(50),
                            merged_iv.end().saturating_add(50),
                        ),
                    };
                    let tss_hits = annot.transcripts_with_tss_in(chrom, lo, hi);
                    let mut distinct_genes: HashSet<&str> = HashSet::new();
                    for tx in tss_hits {
                        distinct_genes.insert(tx.gene_name.as_str());
                    }
                    distinct_genes.len() > cfg.max_tss_genes && cluster.len() > 1
                } else {
                    false
                }
            } else {
                false
            };

            if should_disaggregate {
                stats.disaggregated_clusters += 1;
                stats.restored_constituents += cluster.len();
                for single in cluster {
                    let c_size = sum_overlapping_widths(ref_slice, single.interval, cfg.compat);
                    out_regions.push(StitchedRegion {
                        id: single.id.clone(),
                        chrom: chrom.clone(),
                        interval: single.interval,
                        strand: single.strand,
                        num_loci: 1,
                        constituent_size: c_size,
                        constituent_ids: vec![single.id],
                        origin: RegionOrigin::Enhancer,
                    });
                }
            } else {
                let count = cluster.len();
                let leftmost_id = &cluster[0].id;
                let stitched_strand = if count == 1 {
                    cluster[0].strand
                } else {
                    Strand::Unstranded
                };
                let stitched_id = format!("{count}_{leftmost_id}_lociStitched");
                let c_size = sum_overlapping_widths(ref_slice, merged_iv, cfg.compat);
                let ids = cluster.into_iter().map(|c| c.id).collect();
                out_regions.push(StitchedRegion {
                    id: stitched_id,
                    chrom: chrom.clone(),
                    interval: merged_iv,
                    strand: stitched_strand,
                    num_loci: count,
                    constituent_size: c_size,
                    constituent_ids: ids,
                    origin: RegionOrigin::Enhancer,
                });
            }
        }
    }

    // Legacy `ROSE_main.py:159` orders `ENHANCER_REGION_MAP.txt` by region length descending.
    out_regions.sort_by(|a, b| {
        b.interval
            .width()
            .cmp(&a.interval.width())
            .then_with(|| a.chrom.cmp(&b.chrom))
            .then_with(|| a.interval.start().cmp(&b.interval.start()))
            .then_with(|| a.id.cmp(&b.id))
    });

    Ok((out_regions, stats))
}

/// Validates unique constituent IDs, removes blacklisted peaks, applies TSS promoter exclusion,
/// stitches constituents within `cfg.stitch_distance`, and disaggregates multi-gene spans.
///
/// # Errors
///
/// Returns [`RoseCoreError::DuplicateRegionId`] if two input constituents share the same `id`.
pub fn stitch_constituents(
    constituents: &[Constituent],
    annotation: Option<&AnnotationIndex>,
    blacklist: &BTreeMap<String, Vec<Interval>>,
    cfg: &StitchConfig,
) -> Result<(Vec<StitchedRegion>, StitchStats), RoseCoreError> {
    let prepared = prepare_constituents(constituents, annotation, blacklist, cfg)?;
    stitch_prepared(&prepared, annotation, cfg.stitch_distance, cfg)
}

pub(crate) fn sweep_stitch_chrom(
    peaks: &[Constituent],
    stitch_distance: u64,
    _compat: CompatMode,
) -> Vec<Vec<Constituent>> {
    let mut clusters: Vec<Vec<Constituent>> = Vec::new();
    for peak in peaks {
        if let Some(cur) = clusters.last_mut() {
            let cur_end = cur.iter().map(|c| c.interval.end()).max().unwrap_or(0);
            if peak.interval.start() <= cur_end.saturating_add(stitch_distance) {
                cur.push(peak.clone());
                continue;
            }
        }
        clusters.push(vec![peak.clone()]);
    }
    clusters
}

fn collect_promoter_windows(
    annot: &AnnotationIndex,
    chrom: &str,
    interval: Interval,
    cfg: &StitchConfig,
) -> Vec<Interval> {
    let lo = interval.start().saturating_sub(cfg.tss_window + 2);
    let hi = interval.end().saturating_add(cfg.tss_window + 2);
    let wins: Vec<Interval> = annot
        .transcripts_with_tss_in(chrom, lo, hi)
        .into_iter()
        .map(|tx| tx.promoter_window(cfg.tss_window, cfg.compat))
        .collect();
    merge_intervals(&wins)
}

/// Merges overlapping or directly abutting 0-based half-open intervals, sorting internally if
/// needed so callers cannot accidentally pass non-monotonic slices.
#[must_use]
pub fn merge_intervals(intervals: &[Interval]) -> Vec<Interval> {
    if intervals.len() <= 1 {
        return intervals.to_vec();
    }
    let mut sorted = intervals.to_vec();
    if !sorted.is_sorted() {
        sorted.sort_unstable();
    }
    let mut merged: Vec<Interval> = Vec::with_capacity(sorted.len());
    for iv in sorted {
        if let Some(last) = merged.last_mut() {
            if iv.start() <= last.end() {
                if let Ok(unioned) = Interval::try_new("", last.start(), last.end().max(iv.end())) {
                    *last = unioned;
                }
                continue;
            }
        }
        merged.push(iv);
    }
    merged
}

fn subtract_many(seed: Interval, cuts: &[Interval]) -> Vec<Interval> {
    let mut current = vec![seed];
    for &cut in cuts {
        current = current
            .into_iter()
            .flat_map(|iv| iv.subtract(cut))
            .collect();
        if current.is_empty() {
            break;
        }
    }
    current
}

fn overlaps_any_sorted(sorted_bl: &[Interval], query: Interval) -> bool {
    let upper = sorted_bl.partition_point(|iv| iv.start() < query.end());
    sorted_bl[..upper].iter().any(|iv| iv.overlaps(query))
}

fn sum_overlapping_widths(ref_sorted: &[Interval], query: Interval, compat: CompatMode) -> u64 {
    let upper = ref_sorted.partition_point(|iv| iv.start() < query.end());
    ref_sorted[..upper]
        .iter()
        .filter(|iv| iv.overlaps(query))
        .map(|iv| match compat {
            CompatMode::Clean => iv.width(),
            CompatMode::RoseLegacy => iv.width() + 1,
        })
        .sum()
}

/// Pools constituent peaks across multiple samples (`rose multi`) and collapses overlapping intervals
/// supported by at least `min_samples` distinct samples into canonical consensus constituents.
///
/// # Errors
///
/// Propagates [`RoseCoreError::InvalidConfig`] if `min_samples == 0`.
pub fn merge_consensus_peaks(
    sample_peaks: &[Vec<Constituent>],
    min_samples: usize,
) -> Result<Vec<Constituent>, RoseCoreError> {
    if min_samples == 0 {
        return Err(RoseCoreError::InvalidConfig(
            "min_peak_samples must be at least 1".to_owned(),
        ));
    }
    let mut by_chrom: BTreeMap<String, Vec<(Interval, usize)>> = BTreeMap::new();
    for (s_idx, peaks) in sample_peaks.iter().enumerate() {
        for p in peaks {
            by_chrom
                .entry(p.chrom.clone())
                .or_default()
                .push((p.interval, s_idx));
        }
    }

    let mut consensus = Vec::new();
    let mut seq_id = 1_usize;
    for (chrom, mut items) in by_chrom {
        items.sort_by_key(|&(iv, s_idx)| (iv.start(), iv.end(), s_idx));
        let mut idx = 0;
        while idx < items.len() {
            let start = items[idx].0.start();
            let mut end = items[idx].0.end();
            let mut supporters = HashSet::new();
            supporters.insert(items[idx].1);
            let mut next = idx + 1;
            while next < items.len() && items[next].0.start() < end {
                end = end.max(items[next].0.end());
                supporters.insert(items[next].1);
                next += 1;
            }
            if supporters.len() >= min_samples {
                let iv = Interval::try_new(&chrom, start, end)?;
                consensus.push(Constituent {
                    id: format!("CONSENSUS_PEAK_{seq_id}"),
                    chrom: chrom.clone(),
                    interval: iv,
                    strand: Strand::Unstranded,
                    score: Some(supporters.len() as f64),
                });
                seq_id += 1;
            }
            idx = next;
        }
    }
    Ok(consensus)
}
