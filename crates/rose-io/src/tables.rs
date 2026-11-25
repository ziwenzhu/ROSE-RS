//! Atomic writers and readers for ROSE-compatible tables, BED tracks, multi-sample matrices,
//! stability reports, and `run.json` provenance manifests.

use std::io::Write as _;
use std::path::Path;

use rose_core::{
    AutoStitchReport, CallResult, CompatMode, EnhancerGeneRow, GeneEnhancerRow, Interval,
    MultiCohortResult, PromoterMode, RankedEnhancer, RegionOrigin, RegionSignal, StabilityReport,
    StitchSensitivityRow, StitchedRegion, Strand,
};
use serde::{Deserialize, Serialize};

use crate::atomic::{write_atomic_bytes, AtomicWriter};
use crate::error::RoseIoError;

/// Provenance record written to `run.json`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    /// `rose-rs` semantic version.
    pub version: String,
    /// Subcommand executed.
    pub command: String,
    /// Resolved CLI arguments / presets.
    pub arguments: Vec<String>,
    /// SHA-256 checksums of input files `(role, path, sha256)`.
    pub input_checksums: Vec<(String, String, String)>,
    /// Ranking track mapped / filtered counts `(index_mapped, filtered_units, effective_mmr)`.
    pub ranking_read_stats: Option<(Option<u64>, Option<u64>, f64)>,
    /// Control track mapped / filtered counts `(index_mapped, filtered_units, effective_mmr)`.
    pub control_read_stats: Option<(Option<u64>, Option<u64>, f64)>,
    /// Computed super-enhancer cutoff signal.
    pub cutoff_signal: Option<f64>,
    /// Total stitched regions evaluated.
    pub num_regions: usize,
    /// Number of super-enhancers called.
    pub num_super_enhancers: Option<usize>,
    /// Optional `--compare-control` summary.
    pub control_comparison_note: Option<String>,
    /// Resolved active stitch distance in bp (`rose run` / `rose multi`).
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub active_stitch_bp: Option<u64>,
    /// Optional `--stitch auto` gap-excess estimation report.
    #[serde(skip_serializing_if = "Option::is_none", default)]
    pub auto_stitch_report: Option<AutoStitchReport>,
}

/// Writes stitched regions to a 9-column GFF (`*_STITCHED.gff`) and/or 6-column BED (`rose stitch`).
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_stitched_gff(
    path: impl AsRef<Path>,
    regions: &[StitchedRegion],
    compat: CompatMode,
) -> Result<(), RoseIoError> {
    let p = path.as_ref();
    let mut aw = AtomicWriter::create(p)?;
    for r in regions {
        let start = match compat {
            CompatMode::Clean => r.interval.start() + 1,
            CompatMode::RoseLegacy => r.interval.start(),
        };
        let end = r.interval.end();
        writeln!(
            aw,
            "{}\t{}\t\t{}\t{}\t\t{}\t\t{}",
            r.chrom, r.id, start, end, r.strand, r.id
        )
        .map_err(|e| RoseIoError::io(p, e))?;
    }
    aw.finish()
}

/// Writes `*_ENHANCER_REGION_MAP.txt` (`rose signal` / `rose run`), recording both ROSE's
/// expected integrated signal columns (`7` and `8`) and explicit `DENSITY_*` columns in clean mode.
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_region_map_tsv(
    path: impl AsRef<Path>,
    signals: &[RegionSignal],
    treat_name: &str,
    ctrl_name: Option<&str>,
    compat: CompatMode,
) -> Result<(), RoseIoError> {
    let p = path.as_ref();
    let mut aw = AtomicWriter::create(p)?;
    write!(
        aw,
        "REGION_ID\tCHROM\tSTART\tSTOP\tNUM_LOCI\tCONSTITUENT_SIZE\t{treat_name}"
    )
    .map_err(|e| RoseIoError::io(p, e))?;
    if let Some(cn) = ctrl_name {
        write!(aw, "\t{cn}").map_err(|e| RoseIoError::io(p, e))?;
    }
    if compat == CompatMode::Clean {
        write!(aw, "\tTREAT_DENSITY\tCTRL_DENSITY\tCN_RATIO").map_err(|e| RoseIoError::io(p, e))?;
    }
    writeln!(aw).map_err(|e| RoseIoError::io(p, e))?;

    for s in signals {
        let r = &s.region;
        write!(
            aw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
            r.id,
            r.chrom,
            r.interval.start(),
            r.interval.end(),
            r.num_loci,
            r.constituent_size,
            s.treat_signal
        )
        .map_err(|e| RoseIoError::io(p, e))?;
        if ctrl_name.is_some() {
            write!(aw, "\t{:.4}", s.ctrl_signal.unwrap_or(0.0))
                .map_err(|e| RoseIoError::io(p, e))?;
        }
        if compat == CompatMode::Clean {
            write!(
                aw,
                "\t{:.6}\t{:.6}\t{:.4}",
                s.treat_density,
                s.ctrl_density.unwrap_or(0.0),
                s.cn_ratio.unwrap_or(1.0)
            )
            .map_err(|e| RoseIoError::io(p, e))?;
        }
        writeln!(aw).map_err(|e| RoseIoError::io(p, e))?;
    }
    aw.finish()
}

/// Reads a `region_map.tsv` / `*_ENHANCER_REGION_MAP.txt` back into [`RegionSignal`] records
/// (`rose call` standalone subcommand).
///
/// # Errors
///
/// Returns [`RoseIoError`] if the table is unreadable or malformed.
pub fn read_region_map_tsv(
    path: impl AsRef<Path>,
) -> Result<(Vec<RegionSignal>, String, Option<String>), RoseIoError> {
    let p = path.as_ref();
    let content = std::fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let mut lines = content
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty() && !l.starts_with('#'));

    let Some((hdr_idx, hdr_line)) = lines.next() else {
        return Err(RoseIoError::Parse {
            path: p.to_path_buf(),
            line: 1,
            message: "empty region_map TSV".to_owned(),
        });
    };
    let headers: Vec<&str> = hdr_line.trim_end_matches('\r').split('\t').collect();
    if headers.len() < 7 {
        return Err(RoseIoError::Parse {
            path: p.to_path_buf(),
            line: hdr_idx + 1,
            message: "region_map TSV requires at least 7 columns".to_owned(),
        });
    }
    let treat_name = headers[6].to_owned();
    let has_ctrl = headers
        .get(7)
        .is_some_and(|&h| h != "TREAT_DENSITY" && h != "enhancerRank");
    let ctrl_name = has_ctrl.then(|| headers[7].to_owned());
    let cn_col_idx = headers.iter().position(|&h| h == "CN_RATIO");

    let mut out = Vec::new();
    for (idx0, line) in lines {
        let cols: Vec<&str> = line.trim_end_matches('\r').split('\t').collect();
        if cols.len() < 7 {
            continue;
        }
        let id = cols[0].to_owned();
        let chrom = cols[1].to_owned();
        let start = cols[2].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("bad START: {e}"),
        })?;
        let end = cols[3].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("bad STOP: {e}"),
        })?;
        let num_loci = cols[4].parse::<usize>().unwrap_or(1);
        let constituent_size = cols[5].parse::<u64>().unwrap_or(end.saturating_sub(start));
        let treat_signal = cols[6].parse::<f64>().unwrap_or(0.0);
        let ctrl_signal = if has_ctrl {
            cols.get(7).and_then(|s| s.parse::<f64>().ok())
        } else {
            None
        };
        let cn_ratio = cn_col_idx
            .and_then(|ci| cols.get(ci))
            .and_then(|s| s.parse::<f64>().ok());
        let iv = Interval::try_new(&chrom, start, end.max(start + 1))?;
        let w = iv.width() as f64;
        out.push(RegionSignal {
            region: StitchedRegion {
                id,
                chrom,
                interval: iv,
                strand: Strand::Unstranded,
                num_loci,
                constituent_size,
                constituent_ids: Vec::new(),
                origin: RegionOrigin::Enhancer,
            },
            treat_density: treat_signal / w,
            treat_signal,
            ctrl_density: ctrl_signal.map(|c| c / w),
            ctrl_signal,
            cn_ratio,
        });
    }
    Ok((out, treat_name, ctrl_name))
}

/// Writes `*_AllEnhancers.table.txt` and `*_SuperEnhancers.table.txt` (`ROSE_callSuper.R` layout,
/// plus optional `CLASS`, `RAW_SIGNAL`, `CN_RATIO`, `STABILITY_FREQ`, `BORDERLINE` in clean/extended modes).
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
#[expect(
    clippy::too_many_arguments,
    reason = "table writer emits optional CN and stability columns alongside core ROSE schema"
)]
pub fn write_enhancer_call_tables(
    all_path: impl AsRef<Path>,
    super_path: impl AsRef<Path>,
    sample_name: &str,
    source_file: &str,
    treat_name: &str,
    ctrl_name: Option<&str>,
    call: &CallResult,
    stability: Option<&StabilityReport>,
    has_cn: bool,
    compat: CompatMode,
) -> Result<(), RoseIoError> {
    let write_subset = |path: &Path, only_super: bool| -> Result<(), RoseIoError> {
        let mut aw = AtomicWriter::create(path)?;
        writeln!(aw, "#{sample_name} Enhancers").map_err(|e| RoseIoError::io(path, e))?;
        writeln!(aw, "#Created from {source_file}").map_err(|e| RoseIoError::io(path, e))?;
        writeln!(aw, "#Ranked by {treat_name}").map_err(|e| RoseIoError::io(path, e))?;
        writeln!(
            aw,
            "#Using cutoff of {:.4} for Super-Enhancers",
            call.summary.cutoff_signal
        )
        .map_err(|e| RoseIoError::io(path, e))?;

        write!(
            aw,
            "REGION_ID\tCHROM\tSTART\tSTOP\tNUM_LOCI\tCONSTITUENT_SIZE\t{treat_name}"
        )
        .map_err(|e| RoseIoError::io(path, e))?;
        if let Some(cn) = ctrl_name {
            write!(aw, "\t{cn}").map_err(|e| RoseIoError::io(path, e))?;
        }
        write!(aw, "\tenhancerRank\tisSuper").map_err(|e| RoseIoError::io(path, e))?;
        if compat == CompatMode::Clean {
            write!(aw, "\tCLASS").map_err(|e| RoseIoError::io(path, e))?;
        }
        if has_cn {
            write!(aw, "\tRAW_NET_SIGNAL\tCORRECTED_SIGNAL\tCN_RATIO")
                .map_err(|e| RoseIoError::io(path, e))?;
        }
        if stability.is_some() {
            write!(aw, "\tBOOTSTRAP_FREQ\tBORDERLINE").map_err(|e| RoseIoError::io(path, e))?;
        }
        writeln!(aw).map_err(|e| RoseIoError::io(path, e))?;

        for r in &call.ranked {
            if only_super && !r.is_super {
                continue;
            }
            let reg = &r.item.region;
            write!(
                aw,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}",
                reg.id,
                reg.chrom,
                reg.interval.start(),
                reg.interval.end(),
                reg.num_loci,
                reg.constituent_size,
                r.item.treat_signal
            )
            .map_err(|e| RoseIoError::io(path, e))?;
            if ctrl_name.is_some() {
                write!(aw, "\t{:.4}", r.item.ctrl_signal.unwrap_or(0.0))
                    .map_err(|e| RoseIoError::io(path, e))?;
            }
            write!(aw, "\t{}\t{}", r.enhancer_rank, u8::from(r.is_super))
                .map_err(|e| RoseIoError::io(path, e))?;
            if compat == CompatMode::Clean {
                write!(aw, "\t{}", r.class).map_err(|e| RoseIoError::io(path, e))?;
            }
            if has_cn {
                write!(
                    aw,
                    "\t{:.4}\t{:.4}\t{:.4}",
                    r.raw_net_signal,
                    r.rank_signal,
                    r.item.cn_ratio.unwrap_or(1.0)
                )
                .map_err(|e| RoseIoError::io(path, e))?;
            }
            if let Some(stab) = stability {
                // `call.ranked` is ordered by `enhancer_rank` (`1..=N`); map by `enhancer_rank - 1`
                if let Some(sr) = stab.per_region.get(r.enhancer_rank - 1) {
                    write!(
                        aw,
                        "\t{:.4}\t{}",
                        sr.call_frequency,
                        u8::from(sr.is_borderline)
                    )
                    .map_err(|e| RoseIoError::io(path, e))?;
                }
            }
            writeln!(aw).map_err(|e| RoseIoError::io(path, e))?;
        }
        aw.finish()
    };

    write_subset(all_path.as_ref(), false)?;
    write_subset(super_path.as_ref(), true)?;
    Ok(())
}

/// Reads ranked enhancers from `*_AllEnhancers.table.txt` or `*_SuperEnhancers.table.txt` (`rose genes`).
///
/// # Errors
///
/// Returns [`RoseIoError`] on malformed enhancer tables.
pub fn read_ranked_enhancers_tsv(
    path: impl AsRef<Path>,
) -> Result<Vec<RankedEnhancer>, RoseIoError> {
    let p = path.as_ref();
    let content = std::fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let mut data_lines = content
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty() && !l.starts_with('#'));

    let Some((_, hdr)) = data_lines.next() else {
        return Ok(Vec::new());
    };
    let hcols: Vec<&str> = hdr.trim_end_matches('\r').split('\t').collect();
    let i_rank = hcols
        .iter()
        .position(|&c| c == "enhancerRank")
        .unwrap_or(hcols.len().saturating_sub(2));
    let i_super = hcols
        .iter()
        .position(|&c| c == "isSuper")
        .unwrap_or(hcols.len().saturating_sub(1));

    let mut out = Vec::new();
    for (idx0, line) in data_lines {
        let cols: Vec<&str> = line.trim_end_matches('\r').split('\t').collect();
        if cols.len() <= i_super || cols.len() < 7 {
            continue;
        }
        let chrom = cols[1].to_owned();
        let start = cols[2].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("invalid START: {e}"),
        })?;
        let end = cols[3].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("invalid STOP: {e}"),
        })?;
        let iv = Interval::try_new(&chrom, start, end.max(start + 1))?;
        let treat_sig = cols[6].parse::<f64>().unwrap_or(0.0);
        let enh_rank = cols[i_rank].parse::<usize>().unwrap_or(idx0);
        let is_super = cols[i_super].trim() == "1" || cols[i_super].eq_ignore_ascii_case("true");
        let class = if is_super {
            rose_core::RegionClass::SuperEnhancer
        } else {
            rose_core::RegionClass::Enhancer
        };
        out.push(RankedEnhancer {
            item: RegionSignal {
                region: StitchedRegion {
                    id: cols[0].to_owned(),
                    chrom,
                    interval: iv,
                    strand: Strand::Unstranded,
                    num_loci: cols[4].parse().unwrap_or(1),
                    constituent_size: cols[5].parse().unwrap_or(iv.width()),
                    constituent_ids: Vec::new(),
                    origin: RegionOrigin::Enhancer,
                },
                treat_density: treat_sig / (iv.width() as f64),
                treat_signal: treat_sig,
                ctrl_density: None,
                ctrl_signal: None,
                cn_ratio: None,
            },
            raw_net_signal: treat_sig,
            rank_signal: treat_sig,
            enhancer_rank: enh_rank,
            is_super,
            class,
            bed_score: 1_000,
        });
    }
    Ok(out)
}

/// Writes `*_Enhancers_withSuper.bed`, `*_Gateway_Enhancers.bed`, and `*_Gateway_SuperEnhancers.bed`
/// atomically (never appending to stale files!).
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_enhancer_beds(
    out_dir: &Path,
    prefix_name: &str,
    call: &CallResult,
) -> Result<(), RoseIoError> {
    let with_super = out_dir.join(format!("{prefix_name}_Enhancers_withSuper.bed"));
    let mut aw = AtomicWriter::create(&with_super)?;
    writeln!(
        aw,
        "track name=\"{prefix_name}_Enhancers\" description=\"{prefix_name} Enhancers\" itemRGB=On color=0,0,0"
    )
    .map_err(|e| RoseIoError::io(&with_super, e))?;
    for r in &call.ranked {
        let reg = &r.item.region;
        writeln!(
            aw,
            "{}\t{}\t{}\t{}\t{}",
            reg.chrom,
            reg.interval.start(),
            reg.interval.end(),
            reg.id,
            r.bed_score
        )
        .map_err(|e| RoseIoError::io(&with_super, e))?;
    }
    writeln!(
        aw,
        "\ntrack name=\"Super_{prefix_name}_Enhancers\" description=\"Super {prefix_name} Enhancers\" itemRGB=On color=255,0,0"
    )
    .map_err(|e| RoseIoError::io(&with_super, e))?;
    for r in call.ranked.iter().filter(|r| r.is_super) {
        let reg = &r.item.region;
        writeln!(
            aw,
            "{}\t{}\t{}\t{}\t{}",
            reg.chrom,
            reg.interval.start(),
            reg.interval.end(),
            reg.id,
            r.bed_score
        )
        .map_err(|e| RoseIoError::io(&with_super, e))?;
    }
    aw.finish()?;

    for (suffix, only_super) in [
        ("_Gateway_Enhancers.bed", false),
        ("_Gateway_SuperEnhancers.bed", true),
    ] {
        let gw_path = out_dir.join(format!("{prefix_name}{suffix}"));
        let mut gw = AtomicWriter::create(&gw_path)?;
        for r in &call.ranked {
            if only_super && !r.is_super {
                continue;
            }
            let reg = &r.item.region;
            writeln!(
                gw,
                "{}\t{}\t{}\t{}\t{}\t.",
                reg.chrom,
                reg.interval.start(),
                reg.interval.end(),
                reg.id,
                r.bed_score
            )
            .map_err(|e| RoseIoError::io(&gw_path, e))?;
        }
        gw.finish()?;
    }
    Ok(())
}

/// Writes `*_ENHANCER_TO_GENE.txt` and `*_GENE_TO_ENHANCER.txt` with clean typed headers.
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_gene_tables(
    enh_to_gene_path: impl AsRef<Path>,
    gene_to_enh_path: impl AsRef<Path>,
    enh_rows: &[EnhancerGeneRow],
    gene_rows: &[GeneEnhancerRow],
    compat: CompatMode,
) -> Result<(), RoseIoError> {
    let p1 = enh_to_gene_path.as_ref();
    let mut aw1 = AtomicWriter::create(p1)?;
    write!(
        aw1,
        "REGION_ID\tCHROM\tSTART\tSTOP\tNUM_LOCI\tCONSTITUENT_SIZE\tSIGNAL\tOVERLAP_GENES\tPROXIMAL_GENES\tCLOSEST_GENE\tenhancerRank\tisSuper"
    )
    .map_err(|e| RoseIoError::io(p1, e))?;
    if compat == CompatMode::Clean {
        write!(aw1, "\tCLASS").map_err(|e| RoseIoError::io(p1, e))?;
    }
    writeln!(aw1).map_err(|e| RoseIoError::io(p1, e))?;

    for row in enh_rows {
        let e = &row.enhancer;
        let reg = &e.item.region;
        write!(
            aw1,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
            reg.id,
            reg.chrom,
            reg.interval.start(),
            reg.interval.end(),
            reg.num_loci,
            reg.constituent_size,
            e.rank_signal,
            row.overlap_genes.join(","),
            row.proximal_genes.join(","),
            row.closest_gene.as_deref().unwrap_or(""),
            e.enhancer_rank,
            u8::from(e.is_super)
        )
        .map_err(|e| RoseIoError::io(p1, e))?;
        if compat == CompatMode::Clean {
            write!(aw1, "\t{}", e.class).map_err(|e| RoseIoError::io(p1, e))?;
        }
        writeln!(aw1).map_err(|e| RoseIoError::io(p1, e))?;
    }
    aw1.finish()?;

    let p2 = gene_to_enh_path.as_ref();
    let mut aw2 = AtomicWriter::create(p2)?;
    writeln!(
        aw2,
        "GENE_NAME\tREFSEQ_ID\tPROXIMAL_ENHANCERS\tENHANCER_RANKS\tIS_SUPER"
    )
    .map_err(|e| RoseIoError::io(p2, e))?;
    for gr in gene_rows {
        let ranks_str: Vec<String> = gr.enhancer_ranks.iter().map(usize::to_string).collect();
        writeln!(
            aw2,
            "{}\t{}\t{}\t{}\t{}",
            gr.gene_name,
            gr.transcript_id,
            gr.proximal_enhancers.join(","),
            ranks_str.join(","),
            u8::from(gr.is_super)
        )
        .map_err(|e| RoseIoError::io(p2, e))?;
    }
    aw2.finish()
}

/// Writes hockey-stick plot points TSV, `stability_saturation.tsv` (if present), and `run.json`.
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_points_and_manifest(
    points_tsv: &Path,
    call: &CallResult,
    saturation_tsv: Option<(&Path, &StabilityReport)>,
    manifest_path: &Path,
    manifest: &RunManifest,
) -> Result<(), RoseIoError> {
    let mut aw = AtomicWriter::create(points_tsv)?;
    writeln!(aw, "x_index\tenhancer_rank\tsignal\tis_super")
        .map_err(|e| RoseIoError::io(points_tsv, e))?;
    let n = call.sorted_signals_asc.len();
    for (idx0, &sig) in call.sorted_signals_asc.iter().enumerate() {
        let x1 = idx0 + 1;
        let rank = n - idx0;
        let is_sup = u8::from(sig > call.summary.cutoff_signal);
        writeln!(aw, "{x1}\t{rank}\t{sig:.6}\t{is_sup}")
            .map_err(|e| RoseIoError::io(points_tsv, e))?;
    }
    aw.finish()?;

    if let Some((sat_path, stab)) = saturation_tsv {
        let mut sw = AtomicWriter::create(sat_path)?;
        writeln!(
            sw,
            "read_fraction\tcutoff_signal\tsuper_enhancer_count\tjaccard_vs_full"
        )
        .map_err(|e| RoseIoError::io(sat_path, e))?;
        for pt in &stab.saturation_curve {
            writeln!(
                sw,
                "{:.2}\t{:.4}\t{}\t{:.4}",
                pt.read_fraction, pt.cutoff_signal, pt.super_enhancer_count, pt.jaccard_vs_full
            )
            .map_err(|e| RoseIoError::io(sat_path, e))?;
        }
        sw.finish()?;
    }

    let json_bytes = serde_json::to_vec_pretty(manifest)?;
    write_atomic_bytes(manifest_path, &json_bytes)?;
    let _ = PromoterMode::Contained;
    Ok(())
}

/// Writes `rose multi` cohort outputs (`signal_matrix.tsv`, `consensus_calls.tsv`, and
/// `group_comparison.tsv`).
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_multi_outputs(out_dir: &Path, res: &MultiCohortResult) -> Result<(), RoseIoError> {
    let mat_path = out_dir.join("multi_signal_matrix.tsv");
    let mut mw = AtomicWriter::create(&mat_path)?;
    write!(mw, "REGION_ID\tCHROM\tSTART\tSTOP").map_err(|e| RoseIoError::io(&mat_path, e))?;
    for s in &res.samples {
        write!(mw, "\t{}", s.sample_id).map_err(|e| RoseIoError::io(&mat_path, e))?;
    }
    writeln!(mw).map_err(|e| RoseIoError::io(&mat_path, e))?;
    for row in &res.regions {
        let r = &row.region;
        write!(
            mw,
            "{}\t{}\t{}\t{}",
            r.id,
            r.chrom,
            r.interval.start(),
            r.interval.end()
        )
        .map_err(|e| RoseIoError::io(&mat_path, e))?;
        for &sig in &row.sample_signals {
            write!(mw, "\t{sig:.4}").map_err(|e| RoseIoError::io(&mat_path, e))?;
        }
        writeln!(mw).map_err(|e| RoseIoError::io(&mat_path, e))?;
    }
    mw.finish()?;

    let cons_path = out_dir.join("multi_consensus_calls.tsv");
    let mut cw = AtomicWriter::create(&cons_path)?;
    write!(cw, "REGION_ID\tCHROM\tSTART\tSTOP").map_err(|e| RoseIoError::io(&cons_path, e))?;
    for s in &res.samples {
        write!(cw, "\t{}_isSuper", s.sample_id).map_err(|e| RoseIoError::io(&cons_path, e))?;
    }
    for g in &res.group_order {
        write!(cw, "\t{g}_replicates\t{g}_consensus")
            .map_err(|e| RoseIoError::io(&cons_path, e))?;
    }
    if res.group_order.len() == 2 {
        write!(cw, "\tDIFF_STATUS\tLOG2_FC").map_err(|e| RoseIoError::io(&cons_path, e))?;
    }
    writeln!(cw).map_err(|e| RoseIoError::io(&cons_path, e))?;

    for row in &res.regions {
        let r = &row.region;
        write!(
            cw,
            "{}\t{}\t{}\t{}",
            r.id,
            r.chrom,
            r.interval.start(),
            r.interval.end()
        )
        .map_err(|e| RoseIoError::io(&cons_path, e))?;
        for &called in &row.sample_calls {
            write!(cw, "\t{}", u8::from(called)).map_err(|e| RoseIoError::io(&cons_path, e))?;
        }
        for g in &res.group_order {
            write!(
                cw,
                "\t{}\t{}",
                row.group_call_counts[g],
                u8::from(row.group_consensus[g])
            )
            .map_err(|e| RoseIoError::io(&cons_path, e))?;
        }
        if let (Some(lbl), Some(fc)) = (row.diff_label, row.log2_fold_change) {
            let s_lbl = match lbl {
                rose_core::DiffLabel::Gained => "gained",
                rose_core::DiffLabel::Lost => "lost",
                rose_core::DiffLabel::SharedSuper => "shared",
                rose_core::DiffLabel::Typical => "typical",
            };
            write!(cw, "\t{s_lbl}\t{fc:.4}").map_err(|e| RoseIoError::io(&cons_path, e))?;
        }
        writeln!(cw).map_err(|e| RoseIoError::io(&cons_path, e))?;
    }
    cw.finish()
}

/// Writes the multi-distance stitch sensitivity report (`*_Stitch_sensitivity.tsv`).
///
/// # Errors
///
/// Returns [`RoseIoError`] on write failure.
pub fn write_stitch_sensitivity_tsv(
    path: impl AsRef<Path>,
    rows: &[StitchSensitivityRow],
) -> Result<(), RoseIoError> {
    let p = path.as_ref();
    let mut aw = AtomicWriter::create(p)?;
    writeln!(
        aw,
        "STITCH_BP\tIS_ACTIVE\tNUM_REGIONS\tNUM_SUPER\tCUTOFF_SIGNAL\tMEDIAN_SE_WIDTH_BP\tP90_SE_WIDTH_BP\tMEDIAN_SE_PEAK_COVERAGE\tMEDIAN_ALL_PEAK_COVERAGE\tBASE_JACCARD_VS_12_5KB\tSE_OVERLAP_FRAC_VS_12_5KB"
    )
    .map_err(|e| RoseIoError::io(p, e))?;

    for row in rows {
        writeln!(
            aw,
            "{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            row.stitch_bp,
            u8::from(row.is_active),
            row.num_regions,
            row.num_super_enhancers,
            row.cutoff_signal,
            row.median_se_width_bp,
            row.p90_se_width_bp,
            row.median_se_peak_coverage,
            row.median_all_peak_coverage,
            row.base_jaccard_vs_12_5kb,
            row.se_overlap_frac_vs_12_5kb
        )
        .map_err(|e| RoseIoError::io(p, e))?;
    }
    aw.finish()
}
