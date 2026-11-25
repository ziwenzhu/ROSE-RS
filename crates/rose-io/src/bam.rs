//! BAM and CRAM alignment reader, read/fragment filter kernel, and parallel `rayon` signal
//! engine (`rose signal` / `rose run` / `rose multi`).

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::path::Path;

use noodles::bam;
use noodles::core::{Position, Region};
use noodles::cram;
use noodles::csi::{self, BinningIndex as _};
use noodles::sam::{
    self,
    alignment::record::{cigar::op::Kind, Cigar as _, Flags},
};
use rayon::prelude::*;
use rose_core::{
    merge_intervals, normalize_coverage, CompatMode, FragmentMode, Interval, NormDenom, ReadSpan,
    SignalScratch, StitchedRegion,
};

use crate::error::RoseIoError;

/// Default bitwise SAM flags excluded in clean modes (`0x4 | 0x100 | 0x200 | 0x400 | 0x800` = `3844`).
pub const DEFAULT_EXCLUDE_FLAGS: u16 = 0x0004 | 0x0100 | 0x0200 | 0x0400 | 0x0800;

/// Parameters controlling BAM/CRAM read decoding, filtering, and region coverage accumulation.
#[derive(Debug, Clone)]
pub struct BamSignalConfig {
    /// Single-end 5'->3' read extension in bp (`200` by default).
    pub extend: u64,
    /// Per-base coverage floor (`1` in `chip`, `0` in `cutrun`).
    pub floor: u32,
    /// Minimum mapping quality (`0` by default).
    pub min_mapq: u8,
    /// Excluded bitwise SAM flags (`3844` by default; `4` under `--compat rose`).
    pub exclude_flags: u16,
    /// Fragment / span derivation mode (`Auto`, `Cigar`, `Seq`, `Tlen`).
    pub fragment_mode: FragmentMode,
    /// Minimum paired-end template length (`0` by default).
    pub min_fragment: u64,
    /// Maximum paired-end template length (`1_000` by default).
    pub max_fragment: u64,
    /// Scaling denominator mode (`Filtered` vs `Index`).
    pub norm_denom: NormDenom,
    /// Optional explicit scale factor overriding RPM (`--scale-factor`).
    pub scale_factor: Option<f64>,
    /// Compatibility policy (`Clean` vs `RoseLegacy`).
    pub compat: CompatMode,
    /// Optional nested downsampling fractions (`[0.25, 0.50, 0.75, 0.90, 1.0]` when `--stability` is active).
    pub downsample_fractions: Vec<f64>,
    /// Seed for `xxh3_64` QNAME hashing.
    pub seed: u64,
    /// Optional fixed bin width (`100_000`) to simultaneously collect bin counts (`--cn-correct`).
    pub collect_cn_bin_size: Option<u64>,
}

impl Default for BamSignalConfig {
    fn default() -> Self {
        Self {
            extend: 200,
            floor: 1,
            min_mapq: 0,
            exclude_flags: DEFAULT_EXCLUDE_FLAGS,
            fragment_mode: FragmentMode::Auto,
            min_fragment: 0,
            max_fragment: 1_000,
            norm_denom: NormDenom::Filtered,
            scale_factor: None,
            compat: CompatMode::Clean,
            downsample_fractions: Vec::new(),
            seed: 42,
            collect_cn_bin_size: None,
        }
    }
}

/// Metadata and computed signals returned from one BAM/CRAM evaluation pass.
#[derive(Debug, Clone)]
pub struct BamSignalOutput {
    /// Per-region `(density_rpm_per_bp, total_signal_rpm_bp)` aligned 1:1 with input `regions`.
    pub region_signals: Vec<(f64, f64)>,
    /// Raw integer base-pair coverage sums aligned 1:1 with input `regions` (used for exact
    /// piecewise stitch-distance summation).
    pub raw_bp_sums: Vec<u64>,
    /// Raw integer base-pair coverage sums per downsampling tier (`tier x region`).
    pub raw_tier_bp_sums: Vec<Vec<u64>>,
    /// Optional per-tier total signals `(fraction, Vec<f64>)` when `downsample_fractions` is non-empty.
    pub downsampled_totals: Vec<(f64, Vec<f64>)>,
    /// Optional per-chromosome fixed-bin fragment counts collected during streaming (`--cn-correct`).
    pub cn_bin_counts: BTreeMap<String, Vec<u64>>,
    /// Raw mapped record count from `.bai`/`.csi` index metadata (if available).
    pub index_mapped_count: Option<u64>,
    /// Post-filter fragment/read count observed during streaming (if streamed).
    pub filtered_unit_count: Option<u64>,
    /// Effective million-reads denominator used (`MMR`).
    pub effective_mmr: f64,
    /// `true` if `@PG` headers mention a duplicate-marking tool (`markdup`, `MarkDuplicates`, etc.).
    pub has_duplicate_marker_pg: bool,
    /// `true` if single-end mapped reads were encountered.
    pub saw_single_end_reads: bool,
}

/// Header summary extracted before running region queries.
#[derive(Debug, Clone)]
pub struct AlignmentHeaderInfo {
    /// Reference contigs `(name, length_bp)`.
    pub contigs: Vec<(String, u64)>,
    /// Whether any `@PG` record indicates duplicate marking.
    pub has_duplicate_marker_pg: bool,
    /// True if file is CRAM rather than BAM.
    pub is_cram: bool,
}

/// Reads the SAM header of a BAM or CRAM file without scanning alignment records.
///
/// # Errors
///
/// Returns [`RoseIoError`] if the file cannot be opened or its header is malformed.
pub fn inspect_alignment_header(
    path: impl AsRef<Path>,
) -> Result<AlignmentHeaderInfo, RoseIoError> {
    let p = path.as_ref();
    let is_cram = p
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("cram"));

    let header = if is_cram {
        let mut reader = File::open(p)
            .map(cram::io::Reader::new)
            .map_err(|e| RoseIoError::io(p, e))?;
        reader.read_header().map_err(|e| RoseIoError::io(p, e))?
    } else {
        let mut reader = File::open(p)
            .map(bam::io::Reader::new)
            .map_err(|e| RoseIoError::io(p, e))?;
        reader.read_header().map_err(|e| RoseIoError::io(p, e))?
    };

    let contigs: Vec<(String, u64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            (
                String::from_utf8_lossy(name.as_ref()).into_owned(),
                usize::from(map.length()) as u64,
            )
        })
        .collect();

    let header_text = format!("{header:?}").to_ascii_lowercase();
    let has_dup_pg = [
        "markdup",
        "markduplicates",
        "picard",
        "sambamba",
        "samblaster",
        "dedup",
        "umi_tools",
        "bamsormadup",
    ]
    .iter()
    .any(|kw| header_text.contains(kw));

    Ok(AlignmentHeaderInfo {
        contigs,
        has_duplicate_marker_pg: has_dup_pg,
        is_cram,
    })
}

fn read_bam_index_with_fallback(bam_path: &Path) -> Result<bam::Index, RoseIoError> {
    if let Ok(idx) = bam::fs::read_associated_index(bam_path) {
        return Ok(idx);
    }
    let alt_bai = bam_path.with_extension("bai");
    if let Ok(idx) = bam::bai::fs::read(&alt_bai) {
        return Ok(bam::Index::Bai(idx));
    }
    let alt_csi = bam_path.with_extension("csi");
    if let Ok(idx) = csi::fs::read(&alt_csi) {
        return Ok(bam::Index::Csi(idx));
    }
    Err(RoseIoError::MissingIndex {
        path: bam_path.to_path_buf(),
    })
}

/// Evaluates region coverage across `regions` from a BAM or CRAM file using `rayon` worker threads.
///
/// Automatically chooses between merged indexed window queries (`Bam::query`) and single-pass
/// streaming based on `cfg.norm_denom`, `cfg.collect_cn_bin_size`, format (`BAM` vs `CRAM`), and
/// genomic coverage fraction.
///
/// # Errors
///
/// Returns [`RoseIoError`] on missing indexes (when indexed mode is required) or corrupted records.
pub fn compute_bam_signals(
    path: impl AsRef<Path>,
    regions: &[StitchedRegion],
    chrom_map: &HashMap<String, String>,
    cfg: &BamSignalConfig,
) -> Result<BamSignalOutput, RoseIoError> {
    let p = path.as_ref().to_path_buf();
    let info = inspect_alignment_header(&p)?;
    let contig_lens: BTreeMap<String, u64> = info.contigs.iter().cloned().collect();
    let total_genome_bp: u64 = info.contigs.iter().map(|(_, l)| *l).sum::<u64>().max(1);

    // Group target region indices by reconciled BAM contig name.
    let mut by_contig: BTreeMap<String, Vec<(usize, Interval)>> = BTreeMap::new();
    let mut region_bp_sum = 0_u64;
    for (r_idx, reg) in regions.iter().enumerate() {
        let target_chr = chrom_map
            .get(&reg.chrom)
            .cloned()
            .unwrap_or_else(|| reg.chrom.clone());
        if contig_lens.contains_key(&target_chr) {
            region_bp_sum += reg.interval.width();
            by_contig
                .entry(target_chr)
                .or_default()
                .push((r_idx, reg.interval));
        }
    }

    let bam_index_opt = if info.is_cram {
        None
    } else {
        read_bam_index_with_fallback(&p).ok()
    };

    let index_mapped_count: Option<u64> = bam_index_opt.as_ref().and_then(|idx| {
        let mut total = 0_u64;
        let mut any_meta = false;
        for rs in idx.reference_sequences() {
            if let Some(meta) = rs.metadata() {
                total += meta.mapped_record_count();
                any_meta = true;
            }
        }
        any_meta.then_some(total)
    });

    let upstream_pad = match cfg.compat {
        CompatMode::RoseLegacy => cfg.extend,
        CompatMode::Clean => cfg.extend.max(cfg.max_fragment),
    };
    let downstream_pad = cfg.extend;

    // Estimate coalesced BGZF fetch-window count and span (64 KB BGZF blocks span ~65 kb+).
    // When > 512 disjoint BGZF seeks are required across the genome, per-chromosome contiguous
    // BAI streaming (`run_streaming_pass`, still using `.bai` `idxstats` for `NormDenom::Index`)
    // eliminates redundant block re-decompression and is strictly faster.
    let mut coalesced_win_count = 0_usize;
    let coalesced_query_bp: u64 = by_contig
        .iter()
        .map(|(contig, reg_list)| {
            let chr_len = contig_lens.get(contig).copied().unwrap_or(u64::MAX);
            let raw_wins: Vec<Interval> = reg_list
                .iter()
                .filter_map(|&(_, iv)| {
                    let s = iv.start().saturating_sub(upstream_pad);
                    let e = iv
                        .end()
                        .saturating_add(downstream_pad)
                        .min(chr_len)
                        .max(s + 1);
                    Interval::try_new(contig, s, e).ok()
                })
                .collect();
            let merged = merge_nearby_windows(&raw_wins, 65_536);
            coalesced_win_count += merged.len();
            merged.iter().map(|w| w.width()).sum::<u64>()
        })
        .sum();

    let can_use_indexed = !info.is_cram
        && cfg.collect_cn_bin_size.is_none()
        && !(cfg.scale_factor.is_none() && cfg.norm_denom == NormDenom::Filtered)
        && !(cfg.scale_factor.is_none()
            && cfg.norm_denom == NormDenom::Index
            && index_mapped_count.is_none()
            && bam_index_opt.is_some());

    if can_use_indexed {
        let Some(bam_index) = bam_index_opt.as_ref() else {
            return Err(RoseIoError::MissingIndex { path: p });
        };
        if coalesced_win_count <= 512
            && (coalesced_query_bp as f64) / (total_genome_bp as f64) < 0.25
            && (region_bp_sum as f64) / (total_genome_bp as f64) < 0.25
        {
            return run_indexed_bam_queries(
                &p,
                bam_index,
                regions,
                &by_contig,
                &contig_lens,
                index_mapped_count.unwrap_or(1),
                info.has_duplicate_marker_pg,
                cfg,
            );
        }
    }

    run_streaming_pass(
        &p,
        info.is_cram,
        bam_index_opt.as_ref(),
        regions,
        &by_contig,
        &info.contigs,
        index_mapped_count,
        info.has_duplicate_marker_pg,
        cfg,
    )
}

#[expect(
    clippy::too_many_arguments,
    reason = "internal execution pass passes pre-resolved header and index metadata"
)]
fn run_indexed_bam_queries(
    bam_path: &Path,
    bam_index: &bam::Index,
    regions: &[StitchedRegion],
    by_contig: &BTreeMap<String, Vec<(usize, Interval)>>,
    contig_lens: &BTreeMap<String, u64>,
    idx_mapped: u64,
    has_dup_pg: bool,
    cfg: &BamSignalConfig,
) -> Result<BamSignalOutput, RoseIoError> {
    let upstream_pad = match cfg.compat {
        CompatMode::RoseLegacy => cfg.extend,
        CompatMode::Clean => cfg.extend.max(cfg.max_fragment),
    };
    let downstream_pad = cfg.extend;

    let contig_tasks: Vec<(&String, &Vec<(usize, Interval)>)> = by_contig.iter().collect();
    let per_contig_res: Result<Vec<_>, RoseIoError> = contig_tasks
        .par_iter()
        .map(|&(contig, reg_list)| {
            let mut file = File::open(bam_path).map_err(|e| RoseIoError::io(bam_path, e))?;
            let mut reader = bam::io::Reader::new(&mut file);
            let header = reader
                .read_header()
                .map_err(|e| RoseIoError::io(bam_path, e))?;
            let chr_len = contig_lens.get(contig).copied().unwrap_or(u64::MAX);

            // Merge nearby widened query windows (gap <= 65_536 bp) so enhancers sharing a 64 KB BGZF
            // block are decoded in a single seek, then dispatch reads to regions via binary search.
            let mut sorted_regs = reg_list.clone();
            sorted_regs.sort_by_key(|&(_, iv)| (iv.start(), iv.end()));

            let raw_wins: Vec<Interval> = sorted_regs
                .iter()
                .filter_map(|&(_, iv)| {
                    let s = iv.start().saturating_sub(upstream_pad);
                    let e = iv
                        .end()
                        .saturating_add(downstream_pad)
                        .min(chr_len)
                        .max(s + 1);
                    Interval::try_new(contig, s, e).ok()
                })
                .collect();
            let merged_wins = merge_nearby_windows(&raw_wins, 65_536);

            let mut scratch = SignalScratch::new();
            let mut local_bp: Vec<(usize, u64, Vec<u64>)> = Vec::with_capacity(sorted_regs.len());
            let mut saw_se = false;
            let mut reg_cursor = 0_usize;

            for win in merged_wins {
                let mut win_end_cursor = reg_cursor;
                while win_end_cursor < sorted_regs.len() {
                    let (_, reg_iv) = sorted_regs[win_end_cursor];
                    if !win.contains(reg_iv) && reg_iv.start() >= win.end() {
                        break;
                    }
                    win_end_cursor += 1;
                }
                let win_regs = &sorted_regs[reg_cursor..win_end_cursor];
                if win_regs.is_empty() {
                    continue;
                }

                let start_1b =
                    Position::try_from((win.start() as usize) + 1).unwrap_or(Position::MIN);
                let end_1b = Position::try_from(win.end() as usize).unwrap_or(Position::MIN);
                let query_reg = Region::new(contig.as_str(), start_1b..=end_1b);
                let query = reader
                    .query(&header, bam_index, &query_reg)
                    .map_err(|e| RoseIoError::io(bam_path, e))?;

                let mut buckets: Vec<Vec<ReadSpan>> = vec![Vec::new(); win_regs.len()];
                for rec_res in query.records() {
                    let rec = rec_res.map_err(|e| RoseIoError::io(bam_path, e))?;
                    if let Some((span, is_se)) = decode_bam_record(&rec, cfg) {
                        saw_se = saw_se || is_se;
                        let upper =
                            win_regs.partition_point(|&(_, iv)| iv.start() < span.interval.end());
                        for (b_idx, &(_, reg_iv)) in win_regs[..upper].iter().enumerate().rev() {
                            if reg_iv.overlaps(span.interval) {
                                buckets[b_idx].push(span.clone());
                            } else if reg_iv.end() + 50_000 < span.interval.start() {
                                break;
                            }
                        }
                    }
                }

                for (&(orig_idx, reg_iv), spans) in win_regs.iter().zip(buckets.iter()) {
                    let bp = scratch.compute_region_bp_sum(reg_iv, spans, cfg.floor, cfg.compat);
                    let tiers = if cfg.downsample_fractions.is_empty() {
                        Vec::new()
                    } else {
                        scratch.compute_multi_tier_bp_sums(
                            reg_iv,
                            spans,
                            cfg.floor,
                            &cfg.downsample_fractions,
                            cfg.compat,
                        )
                    };
                    local_bp.push((orig_idx, bp, tiers));
                }
                reg_cursor = win_end_cursor;
            }
            Ok((local_bp, saw_se))
        })
        .collect();

    let parts = per_contig_res?;
    let n_reg = regions.len();
    let mut bp_sums = vec![0_u64; n_reg];
    let n_tiers = cfg.downsample_fractions.len();
    let mut tier_bp = vec![vec![0_u64; n_reg]; n_tiers];
    let mut saw_se = false;

    for (local_bp, se) in parts {
        saw_se = saw_se || se;
        for (orig_idx, bp, tiers) in local_bp {
            bp_sums[orig_idx] = bp;
            for (t_i, &tb) in tiers.iter().enumerate() {
                tier_bp[t_i][orig_idx] = tb;
            }
        }
    }

    let mmr = (idx_mapped as f64) / 1_000_000.0;
    Ok(finalize_signal_output(
        regions,
        &bp_sums,
        &tier_bp,
        BTreeMap::new(),
        Some(idx_mapped),
        None,
        mmr,
        has_dup_pg,
        saw_se,
        cfg,
    ))
}

#[expect(
    clippy::too_many_arguments,
    reason = "internal streaming pass aggregates coverage, filtered counts, and CN bins"
)]
fn run_streaming_pass(
    path: &Path,
    is_cram: bool,
    bam_index: Option<&bam::Index>,
    regions: &[StitchedRegion],
    by_contig: &BTreeMap<String, Vec<(usize, Interval)>>,
    contigs: &[(String, u64)],
    index_mapped_count: Option<u64>,
    has_dup_pg: bool,
    cfg: &BamSignalConfig,
) -> Result<BamSignalOutput, RoseIoError> {
    let n_reg = regions.len();
    let n_tiers = cfg.downsample_fractions.len();

    // If BAM has an index, parallelize streaming across chromosomes via `rayon`!
    if !is_cram {
        if let Some(idx) = bam_index {
            let par_res: Result<Vec<_>, RoseIoError> = contigs
                .par_iter()
                .map(|(contig, chr_len)| {
                    stream_single_bam_contig(
                        path,
                        idx,
                        contig,
                        *chr_len,
                        by_contig.get(contig).map_or(&[][..], Vec::as_slice),
                        cfg,
                    )
                })
                .collect();

            let mut bp_sums = vec![0_u64; n_reg];
            let mut tier_bp = vec![vec![0_u64; n_reg]; n_tiers];
            let mut cn_bins = BTreeMap::new();
            let mut filtered_total = 0_u64;
            let mut raw_mapped_total = 0_u64;
            let mut saw_se = false;

            for part in par_res? {
                filtered_total += part.filtered_cnt;
                raw_mapped_total += part.raw_mapped_cnt;
                saw_se = saw_se || part.saw_se;
                if !part.bins.is_empty() {
                    cn_bins.insert(part.contig, part.bins);
                }
                for (orig_idx, bp, tiers) in part.local_bp {
                    bp_sums[orig_idx] = bp;
                    for (t_i, &tb) in tiers.iter().enumerate() {
                        tier_bp[t_i][orig_idx] = tb;
                    }
                }
            }

            let denom_units = match cfg.norm_denom {
                NormDenom::Filtered => filtered_total,
                NormDenom::Index => index_mapped_count.unwrap_or(raw_mapped_total),
            };
            let mmr = (denom_units as f64) / 1_000_000.0;
            return Ok(finalize_signal_output(
                regions,
                &bp_sums,
                &tier_bp,
                cn_bins,
                index_mapped_count.or(Some(raw_mapped_total)),
                Some(filtered_total),
                mmr,
                has_dup_pg,
                saw_se,
                cfg,
            ));
        }
    }

    // Sequential fallback for unindexed BAM or CRAM
    stream_sequential_file(path, is_cram, regions, by_contig, contigs, has_dup_pg, cfg)
}

struct ContigStreamPart {
    contig: String,
    local_bp: Vec<(usize, u64, Vec<u64>)>,
    bins: Vec<u64>,
    filtered_cnt: u64,
    raw_mapped_cnt: u64,
    saw_se: bool,
}

fn stream_single_bam_contig(
    bam_path: &Path,
    idx: &bam::Index,
    contig: &str,
    chr_len: u64,
    contig_regs: &[(usize, Interval)],
    cfg: &BamSignalConfig,
) -> Result<ContigStreamPart, RoseIoError> {
    let mut file = File::open(bam_path).map_err(|e| RoseIoError::io(bam_path, e))?;
    let mut reader = bam::io::Reader::new(&mut file);
    let header = reader
        .read_header()
        .map_err(|e| RoseIoError::io(bam_path, e))?;

    let end_pos = Position::try_from((chr_len as usize).max(1)).unwrap_or(Position::MIN);
    let query_reg = Region::new(contig, Position::MIN..=end_pos);
    let query = reader
        .query(&header, idx, &query_reg)
        .map_err(|e| RoseIoError::io(bam_path, e))?;

    let mut sorted_regs = contig_regs.to_vec();
    sorted_regs.sort_by_key(|&(_, iv)| (iv.start(), iv.end()));

    let mut bins = if let Some(bs) = cfg.collect_cn_bin_size {
        let n_b = chr_len.div_ceil(bs) as usize;
        vec![0_u64; n_b]
    } else {
        Vec::new()
    };

    let mut filtered_cnt = 0_u64;
    let mut raw_mapped_cnt = 0_u64;
    let mut saw_se = false;
    let mut overlapping_spans: Vec<Vec<ReadSpan>> = vec![Vec::new(); sorted_regs.len()];

    for rec_res in query.records() {
        let rec = rec_res.map_err(|e| RoseIoError::io(bam_path, e))?;
        if !rec.flags().is_unmapped() {
            raw_mapped_cnt += 1;
        }
        if let Some((span, is_se)) = decode_bam_record(&rec, cfg) {
            saw_se = saw_se || is_se;
            filtered_cnt += 1;
            if let Some(bs) = cfg.collect_cn_bin_size {
                let mid = span.interval.center();
                let b_idx = (mid / bs) as usize;
                if let Some(slot) = bins.get_mut(b_idx) {
                    *slot += 1;
                }
            }
            // Binary-search regions whose `start < span.end`
            let upper = sorted_regs.partition_point(|&(_, iv)| iv.start() < span.interval.end());
            for (slot_idx, &(_, reg_iv)) in sorted_regs[..upper].iter().enumerate().rev() {
                if reg_iv.overlaps(span.interval) {
                    overlapping_spans[slot_idx].push(span.clone());
                } else if reg_iv.end() + 50_000 < span.interval.start() {
                    break;
                }
            }
        }
    }

    let mut scratch = SignalScratch::new();
    let mut local_bp = Vec::with_capacity(sorted_regs.len());
    for ((orig_idx, reg_iv), spans) in sorted_regs.into_iter().zip(overlapping_spans) {
        let bp = scratch.compute_region_bp_sum(reg_iv, &spans, cfg.floor, cfg.compat);
        let tiers = if cfg.downsample_fractions.is_empty() {
            Vec::new()
        } else {
            scratch.compute_multi_tier_bp_sums(
                reg_iv,
                &spans,
                cfg.floor,
                &cfg.downsample_fractions,
                cfg.compat,
            )
        };
        local_bp.push((orig_idx, bp, tiers));
    }

    Ok(ContigStreamPart {
        contig: contig.to_owned(),
        local_bp,
        bins,
        filtered_cnt,
        raw_mapped_cnt,
        saw_se,
    })
}

fn stream_sequential_file(
    path: &Path,
    is_cram: bool,
    regions: &[StitchedRegion],
    by_contig: &BTreeMap<String, Vec<(usize, Interval)>>,
    contigs: &[(String, u64)],
    has_dup_pg: bool,
    cfg: &BamSignalConfig,
) -> Result<BamSignalOutput, RoseIoError> {
    let mut sorted_by_ref: Vec<Vec<(usize, Interval)>> = Vec::with_capacity(contigs.len());
    let mut spans_by_ref: Vec<Vec<Vec<ReadSpan>>> = Vec::with_capacity(contigs.len());
    let mut cn_bins: BTreeMap<String, Vec<u64>> = BTreeMap::new();

    for (name, len) in contigs {
        let mut r_list = by_contig.get(name).cloned().unwrap_or_default();
        r_list.sort_by_key(|&(_, iv)| (iv.start(), iv.end()));
        spans_by_ref.push(vec![Vec::new(); r_list.len()]);
        sorted_by_ref.push(r_list);
        if let Some(bs) = cfg.collect_cn_bin_size {
            cn_bins.insert(name.clone(), vec![0_u64; len.div_ceil(bs) as usize]);
        }
    }

    let mut filtered_cnt = 0_u64;
    let mut raw_mapped_cnt = 0_u64;
    let mut saw_se = false;

    let mut push_span = |ref_id: usize, span: ReadSpan, is_se: bool| {
        saw_se = saw_se || is_se;
        filtered_cnt += 1;
        if let Some(bs) = cfg.collect_cn_bin_size {
            if let Some((c_name, _)) = contigs.get(ref_id) {
                if let Some(b_vec) = cn_bins.get_mut(c_name) {
                    let b_idx = (span.interval.center() / bs) as usize;
                    if let Some(slot) = b_vec.get_mut(b_idx) {
                        *slot += 1;
                    }
                }
            }
        }
        if let Some(r_list) = sorted_by_ref.get(ref_id) {
            let upper = r_list.partition_point(|&(_, iv)| iv.start() < span.interval.end());
            for (s_i, &(_, reg_iv)) in r_list[..upper].iter().enumerate().rev() {
                if reg_iv.overlaps(span.interval) {
                    spans_by_ref[ref_id][s_i].push(span.clone());
                } else if reg_iv.end() + 50_000 < span.interval.start() {
                    break;
                }
            }
        }
    };

    if is_cram {
        let mut reader = File::open(path)
            .map(cram::io::Reader::new)
            .map_err(|e| RoseIoError::io(path, e))?;
        let header = reader.read_header().map_err(|e| RoseIoError::io(path, e))?;
        for rec_res in reader.records(&header) {
            let rec = rec_res.map_err(|e| RoseIoError::io(path, e))?;
            if let Ok(flags) = sam::alignment::Record::flags(&rec) {
                if !flags.is_unmapped() {
                    raw_mapped_cnt += 1;
                }
            }
            if let Some(Ok(ref_id)) = sam::alignment::Record::reference_sequence_id(&rec, &header) {
                if let Some((span, is_se)) = decode_generic_record(&rec, cfg) {
                    push_span(ref_id, span, is_se);
                }
            }
        }
    } else {
        let mut reader = File::open(path)
            .map(bam::io::Reader::new)
            .map_err(|e| RoseIoError::io(path, e))?;
        let _header = reader.read_header().map_err(|e| RoseIoError::io(path, e))?;
        for rec_res in reader.records() {
            let rec = rec_res.map_err(|e| RoseIoError::io(path, e))?;
            if !rec.flags().is_unmapped() {
                raw_mapped_cnt += 1;
            }
            if let Some(Ok(ref_id)) = rec.reference_sequence_id() {
                if let Some((span, is_se)) = decode_bam_record(&rec, cfg) {
                    push_span(ref_id, span, is_se);
                }
            }
        }
    }

    let n_reg = regions.len();
    let n_tiers = cfg.downsample_fractions.len();
    let mut bp_sums = vec![0_u64; n_reg];
    let mut tier_bp = vec![vec![0_u64; n_reg]; n_tiers];

    let mut scratch = SignalScratch::new();
    for (r_list, s_vecs) in sorted_by_ref.into_iter().zip(spans_by_ref) {
        for ((orig_idx, reg_iv), spans) in r_list.into_iter().zip(s_vecs) {
            bp_sums[orig_idx] =
                scratch.compute_region_bp_sum(reg_iv, &spans, cfg.floor, cfg.compat);
            if n_tiers > 0 {
                let tiers = scratch.compute_multi_tier_bp_sums(
                    reg_iv,
                    &spans,
                    cfg.floor,
                    &cfg.downsample_fractions,
                    cfg.compat,
                );
                for (t_i, tb) in tiers.into_iter().enumerate() {
                    tier_bp[t_i][orig_idx] = tb;
                }
            }
        }
    }

    let denom_units = match cfg.norm_denom {
        NormDenom::Filtered => filtered_cnt,
        NormDenom::Index => raw_mapped_cnt,
    };
    let mmr = (denom_units as f64) / 1_000_000.0;
    Ok(finalize_signal_output(
        regions,
        &bp_sums,
        &tier_bp,
        cn_bins,
        Some(raw_mapped_cnt),
        Some(filtered_cnt),
        mmr,
        has_dup_pg,
        saw_se,
        cfg,
    ))
}

fn decode_bam_record(rec: &bam::Record, cfg: &BamSignalConfig) -> Option<(ReadSpan, bool)> {
    decode_generic_record(rec, cfg)
}

fn decode_generic_record<R: sam::alignment::Record>(
    rec: &R,
    cfg: &BamSignalConfig,
) -> Option<(ReadSpan, bool)> {
    let flags: Flags = rec.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }

    let raw_bits = u16::from(flags);
    let eff_exclude = match cfg.compat {
        CompatMode::RoseLegacy => 0x0004_u16,
        CompatMode::Clean => cfg.exclude_flags | 0x0004_u16,
    };
    if (raw_bits & eff_exclude) != 0 {
        return None;
    }

    if cfg.compat == CompatMode::Clean && cfg.min_mapq > 0 {
        let mq = rec
            .mapping_quality()
            .and_then(Result::ok)
            .map_or(255, u8::from);
        if mq < cfg.min_mapq {
            return None;
        }
    }

    let pos_1b = usize::from(rec.alignment_start()?.ok()?) as u64;
    let pos0 = pos_1b.saturating_sub(1);

    let cigar = rec.cigar();
    if cfg.compat == CompatMode::RoseLegacy {
        for op in cigar.iter().flatten() {
            if op.kind() == Kind::Skip {
                return None;
            }
        }
    }

    let use_paired_tlen = cfg.compat == CompatMode::Clean
        && matches!(cfg.fragment_mode, FragmentMode::Auto | FragmentMode::Tlen)
        && flags.is_segmented();

    let qname_bytes: &[u8] = rec.name().map_or(b"*", |b| b.as_ref());
    let sub_hash = if cfg.downsample_fractions.is_empty() {
        0
    } else {
        ReadSpan::hash_qname(qname_bytes, cfg.seed)
    };

    if use_paired_tlen {
        if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
            return None;
        }
        let tlen = rec.template_length().ok()?;
        let mate_pos0 = rec
            .mate_alignment_start()
            .and_then(Result::ok)
            .map(|p| (usize::from(p) as u64).saturating_sub(1));

        let is_canonical_left = if tlen > 0 {
            if mate_pos0 == Some(pos0) {
                flags.is_first_segment()
            } else {
                true
            }
        } else if tlen == 0 && mate_pos0 == Some(pos0) {
            flags.is_first_segment()
        } else {
            false
        };

        if !is_canonical_left {
            return None;
        }
        let frag_len = (tlen.unsigned_abs() as u64).max(1);
        if frag_len < cfg.min_fragment || frag_len > cfg.max_fragment {
            return None;
        }
        let iv = Interval::try_new("", pos0, pos0 + frag_len).ok()?;
        return Some((
            ReadSpan {
                interval: iv,
                subsample_hash: sub_hash,
            },
            false,
        ));
    }

    if cfg.fragment_mode == FragmentMode::Tlen {
        return None;
    }

    let read_len = match (cfg.compat, cfg.fragment_mode) {
        (CompatMode::RoseLegacy, _) | (_, FragmentMode::Seq) => rec.sequence().len() as u64,
        (_, FragmentMode::Auto | FragmentMode::Cigar) => {
            cigar.alignment_span().ok().unwrap_or(0) as u64
        }
        (_, FragmentMode::Tlen) => return None,
    }
    .max(1);

    let end0 = pos0 + read_len;
    let (ext_start, ext_end) = if flags.is_reverse_complemented() {
        (pos0.saturating_sub(cfg.extend), end0)
    } else {
        (pos0, end0.saturating_add(cfg.extend))
    };
    let iv = Interval::try_new("", ext_start, ext_end.max(ext_start + 1)).ok()?;
    Some((
        ReadSpan {
            interval: iv,
            subsample_hash: sub_hash,
        },
        !flags.is_segmented(),
    ))
}

fn merge_nearby_windows(wins: &[Interval], max_gap: u64) -> Vec<Interval> {
    if wins.is_empty() {
        return Vec::new();
    }
    let mut sorted = wins.to_vec();
    sorted.sort();
    let mut out: Vec<Interval> = Vec::with_capacity(sorted.len());
    for iv in sorted {
        if let Some(last) = out.last_mut() {
            if iv.start() <= last.end().saturating_add(max_gap) {
                if let Ok(u) = Interval::try_new("", last.start(), last.end().max(iv.end())) {
                    *last = u;
                }
                continue;
            }
        }
        out.push(iv);
    }
    merge_intervals(&out)
}

#[expect(
    clippy::too_many_arguments,
    reason = "internal output normaliser packs final signal vectors and provenance counts"
)]
fn finalize_signal_output(
    regions: &[StitchedRegion],
    bp_sums: &[u64],
    tier_bp: &[Vec<u64>],
    cn_bin_counts: BTreeMap<String, Vec<u64>>,
    index_mapped_count: Option<u64>,
    filtered_unit_count: Option<u64>,
    mmr: f64,
    has_duplicate_marker_pg: bool,
    saw_single_end_reads: bool,
    cfg: &BamSignalConfig,
) -> BamSignalOutput {
    let eff_mmr = mmr.max(1e-12);
    let region_signals: Vec<(f64, f64)> = regions
        .iter()
        .zip(bp_sums.iter())
        .map(|(reg, &bp)| {
            normalize_coverage(
                bp,
                reg.interval.width(),
                eff_mmr,
                cfg.scale_factor,
                cfg.compat,
            )
        })
        .collect();

    let downsampled_totals: Vec<(f64, Vec<f64>)> = cfg
        .downsample_fractions
        .iter()
        .zip(tier_bp.iter())
        .map(|(&frac, t_vec)| {
            let tier_mmr = (eff_mmr * frac).max(1e-12);
            let tier_sf = cfg.scale_factor.map(|sf| sf / frac.max(1e-6));
            let sigs: Vec<f64> = regions
                .iter()
                .zip(t_vec.iter())
                .map(|(reg, &bp)| {
                    normalize_coverage(bp, reg.interval.width(), tier_mmr, tier_sf, cfg.compat).1
                })
                .collect();
            (frac, sigs)
        })
        .collect();

    BamSignalOutput {
        region_signals,
        raw_bp_sums: bp_sums.to_vec(),
        raw_tier_bp_sums: tier_bp.to_vec(),
        downsampled_totals,
        cn_bin_counts,
        index_mapped_count,
        filtered_unit_count,
        effective_mmr: eff_mmr,
        has_duplicate_marker_pg,
        saw_single_end_reads,
    }
}
