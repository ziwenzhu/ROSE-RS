//! UCSC `refGene` (`*.ucsc` / `refGene.txt`) and GTF / GFF3 transcript annotation parsers.

use std::collections::btree_map::Entry;
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::str::FromStr;

use rose_core::{AnnotationIndex, CompatMode, GeneBiotypeFilter, Interval, Strand, Transcript};

use crate::error::RoseIoError;

/// Loads an [`AnnotationIndex`] from either a UCSC `refGene` table (`*.ucsc` / `*.txt`) or a
/// GTF / GFF3 annotation file (`*.gtf` / `*.gff` / `*.gff3`), auto-detected from extension and
/// column structure.
///
/// # Errors
///
/// Returns [`RoseIoError`] if the file cannot be read or contains malformed transcript rows.
pub fn read_annotation_index(
    path: impl AsRef<Path>,
    biotype: GeneBiotypeFilter,
    compat: CompatMode,
) -> Result<AnnotationIndex, RoseIoError> {
    let p = path.as_ref();
    let raw_content = fs::read_to_string(p).map_err(|e| RoseIoError::io(p, e))?;
    let content = raw_content.trim_start_matches('\u{feff}');
    let ext = p
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    let is_gtf = ext == "gtf"
        || ext == "gff"
        || ext == "gff3"
        || content
            .lines()
            .find(|l| !l.trim().is_empty() && !l.starts_with('#'))
            .is_some_and(|first| {
                let cols: Vec<&str> = first.split('\t').collect();
                cols.len() == 9 && (cols[8].contains("gene_id") || cols[8].contains("ID="))
            });

    let raw = if is_gtf {
        parse_gtf_transcripts(p, content)?
    } else {
        parse_ucsc_refgene(p, content)?
    };
    Ok(AnnotationIndex::build(raw, biotype, compat))
}

fn parse_ucsc_refgene(path: &Path, content: &str) -> Result<Vec<Transcript>, RoseIoError> {
    let mut out = Vec::new();
    for (idx0, raw_line) in content.lines().enumerate() {
        let trimmed = raw_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw_line
            .trim_end_matches(&['\r', '\n'][..])
            .split('\t')
            .collect();
        // Skip unprefixed header row like `bin\tname\tchrom...`
        if idx0 == 0 && (cols[0] == "bin" || cols[0] == "name") {
            continue;
        }
        // UCSC refGene with `bin` has strand at col 3 (16 cols); without `bin` has strand at col 2 (15 cols).
        let offset = if cols.len() >= 13 && matches!(cols[3].trim(), "+" | "-") {
            1_usize
        } else if cols.len() >= 12 && matches!(cols[2].trim(), "+" | "-") {
            0_usize
        } else {
            return Err(RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("unrecognized UCSC refGene row with {} columns", cols.len()),
            });
        };

        let transcript_id = cols[offset].trim().to_owned();
        let chrom = cols[offset + 1].trim().to_owned();
        let strand = Strand::from_str(cols[offset + 2])?;
        let tx_start = cols[offset + 3]
            .trim()
            .parse::<u64>()
            .map_err(|e| RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("invalid txStart: {e}"),
            })?;
        let tx_end = cols[offset + 4]
            .trim()
            .parse::<u64>()
            .map_err(|e| RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("invalid txEnd: {e}"),
            })?;
        let cds_start = cols[offset + 5].trim().parse::<u64>().unwrap_or(tx_start);
        let cds_end = cols[offset + 6].trim().parse::<u64>().unwrap_or(tx_start);
        let gene_name = cols
            .get(offset + 11)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .unwrap_or(transcript_id.as_str())
            .to_owned();

        let interval = Interval::try_new(&chrom, tx_start, tx_end.max(tx_start + 1))?;
        let is_coding = cds_start < cds_end || transcript_id.starts_with("NM_");
        out.push(Transcript {
            transcript_id,
            gene_name,
            chrom,
            strand,
            interval,
            is_protein_coding: is_coding,
        });
    }
    Ok(out)
}

#[derive(Debug, Clone)]
struct ExonAcc {
    gene_name: String,
    chrom: String,
    strand: Strand,
    min_start: u64,
    max_end: u64,
    is_coding: bool,
}

fn parse_gtf_transcripts(path: &Path, content: &str) -> Result<Vec<Transcript>, RoseIoError> {
    let mut explicit_txs: Vec<Transcript> = Vec::new();
    let mut exon_groups: BTreeMap<String, ExonAcc> = BTreeMap::new();
    let mut order_keys: Vec<String> = Vec::new();

    for (idx0, raw_line) in content.lines().enumerate() {
        let trimmed = raw_line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = raw_line
            .trim_end_matches(&['\r', '\n'][..])
            .split('\t')
            .collect();
        if cols.len() < 9 {
            return Err(RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("GTF/GFF row requires 9 columns, found {}", cols.len()),
            });
        }
        let feature = cols[2].trim();
        let is_tx_feature = matches!(
            feature,
            "transcript" | "mRNA" | "lnc_RNA" | "ncRNA" | "primary_transcript"
        );
        let is_exon_or_cds = feature == "exon" || feature == "CDS";
        if !is_tx_feature && !is_exon_or_cds {
            continue;
        }

        let chrom = cols[0].trim().to_owned();
        let start_1b = cols[3]
            .trim()
            .parse::<u64>()
            .map_err(|e| RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("invalid GTF start: {e}"),
            })?;
        let end_1b = cols[4]
            .trim()
            .parse::<u64>()
            .map_err(|e| RoseIoError::Parse {
                path: path.to_path_buf(),
                line: idx0 + 1,
                message: format!("invalid GTF end: {e}"),
            })?;
        let start_0b = start_1b.saturating_sub(1);
        let end_0b = end_1b.max(start_0b + 1);
        let strand = Strand::from_str(cols[6])?;
        let attrs = parse_gtf_attrs(cols[8]);

        let tx_id = attrs
            .get("transcript_id")
            .or_else(|| attrs.get("ID"))
            .or_else(|| attrs.get("Parent"))
            .or_else(|| attrs.get("gene_id"))
            .cloned()
            .unwrap_or_else(|| format!("TX_{}", idx0 + 1));
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene"))
            .or_else(|| attrs.get("Name"))
            .or_else(|| attrs.get("gene_id"))
            .cloned()
            .unwrap_or_else(|| tx_id.clone());
        let biotype_str = attrs
            .get("transcript_biotype")
            .or_else(|| attrs.get("transcript_type"))
            .or_else(|| attrs.get("gene_biotype"))
            .or_else(|| attrs.get("gene_type"))
            .map_or("", String::as_str);
        let is_coding = feature == "CDS"
            || feature == "mRNA"
            || biotype_str == "protein_coding"
            || (biotype_str.is_empty() && tx_id.starts_with("NM_"));

        if is_tx_feature {
            let interval = Interval::try_new(&chrom, start_0b, end_0b)?;
            explicit_txs.push(Transcript {
                transcript_id: tx_id,
                gene_name,
                chrom,
                strand,
                interval,
                is_protein_coding: is_coding || biotype_str.is_empty(),
            });
        } else {
            match exon_groups.entry(tx_id.clone()) {
                Entry::Vacant(vac) => {
                    order_keys.push(tx_id);
                    vac.insert(ExonAcc {
                        gene_name,
                        chrom,
                        strand,
                        min_start: start_0b,
                        max_end: end_0b,
                        is_coding,
                    });
                }
                Entry::Occupied(mut occ) => {
                    let acc = occ.get_mut();
                    acc.min_start = acc.min_start.min(start_0b);
                    acc.max_end = acc.max_end.max(end_0b);
                    acc.is_coding = acc.is_coding || is_coding;
                }
            }
        }
    }

    if !explicit_txs.is_empty() {
        return Ok(explicit_txs);
    }

    // Fallback for UCSC GTFs that contain only `exon` / `CDS` rows without `transcript` rows!
    let mut reconstructed = Vec::with_capacity(order_keys.len());
    for tx_id in order_keys {
        if let Some(acc) = exon_groups.remove(&tx_id) {
            let interval = Interval::try_new(&acc.chrom, acc.min_start, acc.max_end)?;
            reconstructed.push(Transcript {
                transcript_id: tx_id,
                gene_name: acc.gene_name,
                chrom: acc.chrom,
                strand: acc.strand,
                interval,
                is_protein_coding: acc.is_coding,
            });
        }
    }
    Ok(reconstructed)
}

fn parse_gtf_attrs(raw: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    for part in raw.split(';') {
        let p = part.trim();
        if p.is_empty() {
            continue;
        }
        if let Some((k, v)) = p.split_once('=') {
            map.insert(k.trim().to_owned(), v.trim().trim_matches('"').to_owned());
        } else if let Some((k, v)) = p.split_once(char::is_whitespace) {
            map.insert(k.trim().to_owned(), v.trim().trim_matches('"').to_owned());
        }
    }
    map
}
