//! Parsers for constituent peak regions (BED, narrowPeak, broadPeak, GFF), blacklists,
//! copy-number `.seg` files, and `rose multi` sample sheets.

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use rose_core::{
    merge_intervals, CompatMode, Constituent, CopyNumberSegment, Interval, RoseCoreError, Strand,
};

use crate::error::RoseIoError;

/// Detected region file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionFileFormat {
    /// 0-based half-open BED / narrowPeak / broadPeak.
    Bed,
    /// 1-based closed GFF / GTF.
    Gff,
}

impl RegionFileFormat {
    /// Infers region format from file extension (`.gff`/`.gtf` -> `Gff`, otherwise `Bed` unless
    /// first data row has 9 columns with numeric cols 4 & 5).
    #[must_use]
    pub fn detect(path: &Path, first_data_cols: &[&str]) -> Self {
        if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
            let lower = ext.to_ascii_lowercase();
            if lower == "gff" || lower == "gff3" || lower == "gtf" {
                return Self::Gff;
            }
            if matches!(
                lower.as_str(),
                "bed" | "narrowpeak" | "broadpeak" | "gappedpeak"
            ) {
                return Self::Bed;
            }
        }
        if first_data_cols.len() >= 7
            && first_data_cols[1].parse::<u64>().is_err()
            && first_data_cols[3].parse::<u64>().is_ok()
            && first_data_cols[4].parse::<u64>().is_ok()
        {
            Self::Gff
        } else {
            Self::Bed
        }
    }
}

/// Reads constituent peaks from a BED, narrowPeak, broadPeak, or GFF file, normalizing coordinates
/// to 0-based half-open [`Interval`]s and verifying unique constituent IDs.
///
/// # Errors
///
/// Returns [`RoseIoError`] on unreadable files, malformed coordinates, or duplicate IDs.
pub fn read_constituents(
    path: impl AsRef<Path>,
    compat: CompatMode,
) -> Result<Vec<Constituent>, RoseIoError> {
    let p = path.as_ref();
    let raw_content = fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let content = raw_content.trim_start_matches('\u{feff}');

    let mut format_opt: Option<RegionFileFormat> = None;
    let mut out = Vec::new();
    let mut seen_ids = HashSet::new();
    let mut auto_idx = 1_usize;

    for (idx0, raw_line) in content.lines().enumerate() {
        let line_no = idx0 + 1;
        let trimmed = raw_line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track ")
            || trimmed.starts_with("browser ")
        {
            continue;
        }
        let cols: Vec<&str> = raw_line
            .trim_end_matches(&['\r', '\n'][..])
            .split('\t')
            .collect();
        let fmt = *format_opt.get_or_insert_with(|| RegionFileFormat::detect(p, &cols));

        let constituent = match fmt {
            RegionFileFormat::Bed => parse_bed_row(p, line_no, &cols, &mut auto_idx)?,
            RegionFileFormat::Gff => parse_gff_row(p, line_no, &cols, compat)?,
        };

        if !seen_ids.insert(constituent.id.clone()) {
            return Err(RoseIoError::Core(RoseCoreError::DuplicateRegionId {
                id: constituent.id,
            }));
        }
        out.push(constituent);
    }

    Ok(out)
}

fn parse_bed_row(
    path: &Path,
    line_no: usize,
    cols: &[&str],
    auto_idx: &mut usize,
) -> Result<Constituent, RoseIoError> {
    if cols.len() < 3 {
        return Err(RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("expected at least 3 BED columns, found {}", cols.len()),
        });
    }
    let chrom = cols[0].trim().to_owned();
    let start = cols[1]
        .trim()
        .parse::<u64>()
        .map_err(|e| RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("invalid BED start '{}': {e}", cols[1]),
        })?;
    let end = cols[2]
        .trim()
        .parse::<u64>()
        .map_err(|e| RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("invalid BED end '{}': {e}", cols[2]),
        })?;
    let interval = Interval::try_new(&chrom, start, end)?;
    let id = match cols.get(3).map(|s| s.trim()) {
        Some(name) if !name.is_empty() && name != "." => name.to_owned(),
        _ => {
            let gen_id = format!("PEAK_{auto_idx}");
            *auto_idx += 1;
            gen_id
        }
    };
    let score = cols
        .get(6) // narrowPeak/broadPeak `signalValue` if present, else col 4 score
        .or_else(|| cols.get(4))
        .and_then(|s| s.trim().parse::<f64>().ok());
    let strand = match cols.get(5) {
        Some(s) => Strand::from_str(s)?,
        None => Strand::Unstranded,
    };
    Ok(Constituent {
        id,
        chrom,
        interval,
        strand,
        score,
    })
}

fn parse_gff_row(
    path: &Path,
    line_no: usize,
    cols: &[&str],
    compat: CompatMode,
) -> Result<Constituent, RoseIoError> {
    if cols.len() < 7 {
        return Err(RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("expected at least 7 GFF columns, found {}", cols.len()),
        });
    }
    let chrom = cols[0].trim().to_owned();
    let raw_start = cols[3]
        .trim()
        .parse::<u64>()
        .map_err(|e| RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("invalid GFF start '{}': {e}", cols[3]),
        })?;
    let raw_end = cols[4]
        .trim()
        .parse::<u64>()
        .map_err(|e| RoseIoError::Parse {
            path: path.to_path_buf(),
            line: line_no,
            message: format!("invalid GFF end '{}': {e}", cols[4]),
        })?;
    let (lo, hi) = if raw_start <= raw_end {
        (raw_start, raw_end)
    } else {
        (raw_end, raw_start)
    };
    let interval = match compat {
        // Standard GFF is 1-based closed `[lo, hi]` -> 0-based half-open `[lo - 1, hi)`.
        CompatMode::Clean => Interval::try_new(&chrom, lo.saturating_sub(1), hi)?,
        // Note: In legacy ROSE (`ROSE_bamToGFF.py:122` & `ROSE_main.py:205`), GFF coordinates
        // `start` and `end` use width `end - start` (treating `[start, end)` directly when built
        // from `bedToGFF`)! Storing `[lo, hi)` when `lo < hi` preserves exact `end - start` width!
        CompatMode::RoseLegacy => {
            if lo < hi {
                Interval::try_new(&chrom, lo, hi)?
            } else {
                Interval::point(lo)
            }
        }
    };
    let strand = Strand::from_str(cols[6])?;
    let col2 = cols[1].trim();
    let col9 = cols.get(8).map_or("", |s| s.trim());
    let id = if !col2.is_empty() {
        col2.to_owned()
    } else if !col9.is_empty() {
        col9.to_owned()
    } else {
        format!("{chrom}:{strand}:{lo}-{hi}")
    };
    let score = cols.get(5).and_then(|s| s.trim().parse::<f64>().ok());
    Ok(Constituent {
        id,
        chrom,
        interval,
        strand,
        score,
    })
}

/// Reads a BED blacklist into sorted, merged per-chromosome intervals.
///
/// # Errors
///
/// Returns [`RoseIoError`] if the file cannot be read or contains invalid coordinates.
pub fn read_blacklist_bed(
    path: impl AsRef<Path>,
) -> Result<BTreeMap<String, Vec<Interval>>, RoseIoError> {
    let p = path.as_ref();
    let raw_content = fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let content = raw_content.trim_start_matches('\u{feff}');
    let mut raw: BTreeMap<String, Vec<Interval>> = BTreeMap::new();
    for (idx0, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || trimmed.starts_with('#')
            || trimmed.starts_with("track ")
            || trimmed.starts_with("browser ")
        {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() < 3 {
            return Err(RoseIoError::Parse {
                path: p.to_path_buf(),
                line: idx0 + 1,
                message: "blacklist BED requires at least 3 columns".to_owned(),
            });
        }
        let chrom = cols[0].to_owned();
        let start = cols[1].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("invalid start: {e}"),
        })?;
        let end = cols[2].parse::<u64>().map_err(|e| RoseIoError::Parse {
            path: p.to_path_buf(),
            line: idx0 + 1,
            message: format!("invalid end: {e}"),
        })?;
        let iv = Interval::try_new(&chrom, start, end)?;
        raw.entry(chrom).or_default().push(iv);
    }
    let mut merged = BTreeMap::new();
    for (chrom, mut list) in raw {
        list.sort();
        merged.insert(chrom, merge_intervals(&list));
    }
    Ok(merged)
}

/// Reads external copy-number segments from either a 4-column BED (`chrom`, `start`, `end`,
/// `linear_cn_ratio`) or a 6-column CBS `.seg` file (`ID`, `chrom`, `loc.start`, `loc.end`,
/// `num.mark`, `seg.mean` where `seg.mean = log2(CN_ratio)`).
///
/// # Errors
///
/// Returns [`RoseIoError`] on malformed rows.
pub fn read_cn_segments(path: impl AsRef<Path>) -> Result<Vec<CopyNumberSegment>, RoseIoError> {
    let p = path.as_ref();
    let raw_content = fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let content = raw_content.trim_start_matches('\u{feff}');
    let is_seg = p
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("seg"));

    let mut out = Vec::new();
    for (idx0, line) in content.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        // Skip header row if coordinate columns are non-numeric.
        if idx0 == 0
            && (cols.get(1).and_then(|c| c.parse::<u64>().ok()).is_none()
                && cols.get(2).and_then(|c| c.parse::<u64>().ok()).is_none())
        {
            continue;
        }
        if is_seg || cols.len() >= 6 {
            let chrom = cols[1].trim().to_owned();
            let start = cols[2]
                .trim()
                .parse::<u64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid seg start: {e}"),
                })?;
            let end = cols[3]
                .trim()
                .parse::<u64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid seg end: {e}"),
                })?;
            let seg_mean = cols[5]
                .trim()
                .parse::<f64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid seg.mean: {e}"),
                })?;
            let ratio = 2.0_f64.powf(seg_mean);
            let iv = Interval::try_new(&chrom, start.saturating_sub(1), end.max(start))?;
            out.push(CopyNumberSegment {
                chrom,
                interval: iv,
                cn_ratio: ratio,
                is_mappable: ratio.is_finite() && ratio > 0.0,
            });
        } else if cols.len() >= 4 {
            let chrom = cols[0].trim().to_owned();
            let start = cols[1]
                .trim()
                .parse::<u64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid BED start: {e}"),
                })?;
            let end = cols[2]
                .trim()
                .parse::<u64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid BED end: {e}"),
                })?;
            let ratio = cols[3]
                .trim()
                .parse::<f64>()
                .map_err(|e| RoseIoError::Parse {
                    path: p.to_path_buf(),
                    line: idx0 + 1,
                    message: format!("invalid CN ratio: {e}"),
                })?;
            let iv = Interval::try_new(&chrom, start, end)?;
            out.push(CopyNumberSegment {
                chrom,
                interval: iv,
                cn_ratio: ratio,
                is_mappable: ratio.is_finite() && ratio > 0.0,
            });
        } else {
            return Err(RoseIoError::Parse {
                path: p.to_path_buf(),
                line: idx0 + 1,
                message: "CN segment file requires 4 (BED) or 6 (.seg) columns".to_owned(),
            });
        }
    }
    Ok(out)
}

/// Parsed row from a `rose multi --samples sheet.tsv` file.
#[derive(Debug, Clone)]
pub struct SampleSheetRow {
    /// Sample identifier (`sample`).
    pub sample: String,
    /// Condition / group name (`group`).
    pub group: String,
    /// Constituent peak file path (`peaks`).
    pub peaks: PathBuf,
    /// Ranking BAM/CRAM/bigWig path (`bam`).
    pub bam: PathBuf,
    /// Optional control BAM/CRAM/bigWig path (`control`).
    pub control: Option<PathBuf>,
    /// Optional explicit scale factor (`scale_factor`).
    pub scale_factor: Option<f64>,
}

/// Parses a tab-delimited sample sheet (`sample`, `group`, `peaks`, `bam`, optional `control`,
/// optional `scale_factor`). Relative paths are resolved relative to the sheet's parent directory.
///
/// # Errors
///
/// Returns [`RoseIoError`] if required columns are missing or malformed.
pub fn read_sample_sheet(path: impl AsRef<Path>) -> Result<Vec<SampleSheetRow>, RoseIoError> {
    let p = path.as_ref();
    let base_dir = p.parent().unwrap_or_else(|| Path::new("."));
    let raw_content = fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let content = raw_content.trim_start_matches('\u{feff}');
    let mut lines = content
        .lines()
        .enumerate()
        .filter(|(_, l)| !l.trim().is_empty() && !l.trim_start().starts_with('#'));

    let Some((hdr_idx, hdr_line)) = lines.next() else {
        return Err(RoseIoError::Parse {
            path: p.to_path_buf(),
            line: 1,
            message: "empty sample sheet".to_owned(),
        });
    };
    let headers: Vec<String> = hdr_line
        .trim_end_matches('\r')
        .split('\t')
        .map(|h| h.trim().to_ascii_lowercase())
        .collect();

    let find_col = |name: &str| -> Result<usize, RoseIoError> {
        headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| RoseIoError::Parse {
                path: p.to_path_buf(),
                line: hdr_idx + 1,
                message: format!("missing required column '{name}' in sample sheet"),
            })
    };
    let i_sample = find_col("sample")?;
    let i_group = find_col("group")?;
    let i_peaks = find_col("peaks")?;
    let i_bam = find_col("bam")?;
    let i_ctrl = headers.iter().position(|h| h == "control");
    let i_sf = headers.iter().position(|h| h == "scale_factor");

    let resolve = |raw: &str| -> PathBuf {
        let pb = PathBuf::from(raw);
        if pb.is_absolute() {
            pb
        } else {
            base_dir.join(pb)
        }
    };

    let mut rows = Vec::new();
    for (idx0, line) in lines {
        let cols: Vec<&str> = line.trim_end_matches('\r').split('\t').collect();
        let max_req = i_sample.max(i_group).max(i_peaks).max(i_bam);
        if cols.len() <= max_req {
            return Err(RoseIoError::Parse {
                path: p.to_path_buf(),
                line: idx0 + 1,
                message: "fewer columns than header".to_owned(),
            });
        }
        let control = i_ctrl
            .and_then(|idx| cols.get(idx))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("none") && *s != ".")
            .map(resolve);
        let scale_factor = match i_sf
            .and_then(|idx| cols.get(idx))
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && *s != "." && !s.eq_ignore_ascii_case("na"))
        {
            Some(v) => Some(v.parse::<f64>().map_err(|e| RoseIoError::Parse {
                path: p.to_path_buf(),
                line: idx0 + 1,
                message: format!("invalid scale_factor '{v}': {e}"),
            })?),
            None => None,
        };
        rows.push(SampleSheetRow {
            sample: cols[i_sample].trim().to_owned(),
            group: cols[i_group].trim().to_owned(),
            peaks: resolve(cols[i_peaks].trim()),
            bam: resolve(cols[i_bam].trim()),
            control,
            scale_factor,
        });
    }
    Ok(rows)
}
