//! Sorted per-chromosome transcript and Transcription Start Site (TSS) annotation index.

use std::collections::{BTreeMap, HashMap};

use serde::{Deserialize, Serialize};

use crate::types::{CompatMode, GeneBiotypeFilter, Interval, Strand};

/// Annotated transcript record normalized to 0-based half-open `[start, end)` coordinates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Transcript {
    /// Transcript accession (e.g. `NM_002467` or `ENST00000621592`).
    pub transcript_id: String,
    /// Gene symbol (e.g. `MYC`).
    pub gene_name: String,
    /// Chromosome / contig name.
    pub chrom: String,
    /// Transcript strand (`+` or `-`).
    pub strand: Strand,
    /// 0-based half-open transcript body `[tx_start, tx_end)`.
    pub interval: Interval,
    /// True when protein-coding (`cmpl` CDS in refGene, `NM_` prefix, or `protein_coding` in GTF).
    pub is_protein_coding: bool,
}

impl Transcript {
    /// Returns the 0-based 1-bp TSS coordinate (`tx_start` on `+`, `tx_end - 1` on `-`).
    #[must_use]
    #[inline]
    pub const fn tss_pos(&self) -> u64 {
        match self.strand {
            Strand::Reverse => self.interval.end().saturating_sub(1),
            Strand::Forward | Strand::Unstranded => self.interval.start(),
        }
    }

    /// Returns the legacy ROSE TSS coordinate (`txStart` on `+`, `txEnd` on `-` from UCSC refGene).
    #[must_use]
    #[inline]
    pub const fn legacy_rose_tss_coord(&self) -> u64 {
        match self.strand {
            Strand::Reverse => self.interval.end(),
            Strand::Forward | Strand::Unstranded => self.interval.start(),
        }
    }

    /// Builds a promoter window around the TSS clamped to `>= 0`.
    #[must_use]
    pub fn promoter_window(&self, half_window: u64, compat: CompatMode) -> Interval {
        let (lo, hi) = match compat {
            CompatMode::Clean => {
                let center = self.tss_pos();
                (
                    center.saturating_sub(half_window),
                    center.saturating_add(half_window).saturating_add(1),
                )
            }
            // Legacy ROSE builds `Locus(center - w, center + w)` against GFF `[g_start, g_end]`
            // stored as `[g_start, g_end)`, so `contains` checks `center - w <= g_start` and
            // `g_end <= center + w`.
            CompatMode::RoseLegacy => {
                let center = self.legacy_rose_tss_coord();
                (
                    center.saturating_sub(half_window),
                    center.saturating_add(half_window),
                )
            }
        };
        Interval::try_new(&self.chrom, lo, hi.max(lo + 1)).unwrap_or_else(|_| Interval::point(lo))
    }
}

/// Returns `true` when `chrom` is a primary assembly chromosome rather than `_alt`, `_random`,
/// `chrUn`, or decoy/fix patch contigs.
#[must_use]
pub fn is_primary_contig(chrom: &str) -> bool {
    !chrom.contains('_')
        && !chrom.starts_with("chrUn")
        && !chrom.starts_with("HLA-")
        && !chrom.starts_with("GL")
        && !chrom.starts_with("KI")
}

/// Indexed, chromosome-partitioned collection of transcripts for $O(\log N)$ TSS and overlap queries.
#[derive(Debug, Clone, Default)]
pub struct AnnotationIndex {
    /// Deduplicated transcripts in deterministic order.
    transcripts: Vec<Transcript>,
    /// Per-chromosome indices into `transcripts`, sorted ascending by `tss_pos`.
    by_chrom_tss: BTreeMap<String, Vec<(u64, usize)>>,
    /// Per-chromosome indices into `transcripts`, sorted ascending by `interval.start()`,
    /// paired with prefix-maximum `interval.end()` for fast overlap pruning.
    by_chrom_body: BTreeMap<String, Vec<(u64, u64, usize)>>,
}

impl AnnotationIndex {
    /// Builds a deduplicated, sorted [`AnnotationIndex`] from raw parsed transcripts.
    #[must_use]
    pub fn build(raw: Vec<Transcript>, biotype: GeneBiotypeFilter, compat: CompatMode) -> Self {
        let mut seen_pos: HashMap<String, usize> = HashMap::with_capacity(raw.len());
        let mut transcripts: Vec<Transcript> = Vec::with_capacity(raw.len());

        for tx in raw {
            if biotype == GeneBiotypeFilter::ProteinCoding && !tx.is_protein_coding {
                continue;
            }
            if let Some(&existing_idx) = seen_pos.get(&tx.transcript_id) {
                if compat == CompatMode::Clean
                    && !is_primary_contig(&transcripts[existing_idx].chrom)
                    && is_primary_contig(&tx.chrom)
                {
                    transcripts[existing_idx] = tx;
                }
            } else {
                let idx = transcripts.len();
                seen_pos.insert(tx.transcript_id.clone(), idx);
                transcripts.push(tx);
            }
        }

        let mut by_chrom_tss: BTreeMap<String, Vec<(u64, usize)>> = BTreeMap::new();
        let mut by_chrom_body_raw: BTreeMap<String, Vec<(u64, u64, usize)>> = BTreeMap::new();

        for (idx, tx) in transcripts.iter().enumerate() {
            let tss_coord = match compat {
                CompatMode::Clean => tx.tss_pos(),
                CompatMode::RoseLegacy => tx.legacy_rose_tss_coord(),
            };
            by_chrom_tss
                .entry(tx.chrom.clone())
                .or_default()
                .push((tss_coord, idx));
            by_chrom_body_raw
                .entry(tx.chrom.clone())
                .or_default()
                .push((tx.interval.start(), tx.interval.end(), idx));
        }

        for list in by_chrom_tss.values_mut() {
            list.sort_by_key(|&(pos, idx)| (pos, idx));
        }

        let mut by_chrom_body: BTreeMap<String, Vec<(u64, u64, usize)>> = BTreeMap::new();
        for (chrom, mut list) in by_chrom_body_raw {
            list.sort_by_key(|&(start, end, idx)| (start, end, idx));
            let mut prefix_max_end = 0_u64;
            for entry in &mut list {
                prefix_max_end = prefix_max_end.max(entry.1);
                entry.1 = prefix_max_end;
            }
            by_chrom_body.insert(chrom, list);
        }

        Self {
            transcripts,
            by_chrom_tss,
            by_chrom_body,
        }
    }

    /// Number of indexed transcripts.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.transcripts.len()
    }

    /// Returns `true` if the annotation index holds zero transcripts.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.transcripts.is_empty()
    }

    /// Borrows all indexed transcripts.
    #[must_use]
    #[inline]
    pub fn transcripts(&self) -> &[Transcript] {
        &self.transcripts
    }

    /// Rewrites chromosome names (`chr1` <-> `1`) in place using a canonical mapper.
    #[must_use]
    pub fn with_reconciled_chroms<F>(&self, mut map_chrom: F, compat: CompatMode) -> Self
    where
        F: FnMut(&str) -> String,
    {
        let mapped: Vec<Transcript> = self
            .transcripts
            .iter()
            .map(|tx| {
                let mut cloned = tx.clone();
                cloned.chrom = map_chrom(&tx.chrom);
                cloned
            })
            .collect();
        Self::build(mapped, GeneBiotypeFilter::All, compat)
    }

    /// Returns transcripts on `chrom` whose TSS coordinate lies in `[lo_inclusive, hi_inclusive]`.
    #[must_use]
    pub fn transcripts_with_tss_in(
        &self,
        chrom: &str,
        lo_inclusive: u64,
        hi_inclusive: u64,
    ) -> Vec<&Transcript> {
        let Some(entries) = self.by_chrom_tss.get(chrom) else {
            return Vec::new();
        };
        let first = entries.partition_point(|&(pos, _)| pos < lo_inclusive);
        let last = entries.partition_point(|&(pos, _)| pos <= hi_inclusive);
        entries[first..last]
            .iter()
            .map(|&(_, idx)| &self.transcripts[idx])
            .collect()
    }

    /// Returns transcripts on `chrom` whose full body `[tx_start, tx_end)` overlaps `query`.
    #[must_use]
    pub fn transcripts_overlapping(&self, chrom: &str, query: Interval) -> Vec<&Transcript> {
        let Some(entries) = self.by_chrom_body.get(chrom) else {
            return Vec::new();
        };
        // Only transcripts with `start < query.end()` can overlap `query`.
        let upper = entries.partition_point(|&(start, _, _)| start < query.end());
        // Skip prefix entries whose prefix-maximum `end <= query.start()`.
        let lower = entries[..upper].partition_point(|&(_, max_end, _)| max_end <= query.start());
        let mut hits: Vec<(usize, &Transcript)> = entries[lower..upper]
            .iter()
            .filter_map(|&(_, _, idx)| {
                let tx = &self.transcripts[idx];
                tx.interval.overlaps(query).then_some((idx, tx))
            })
            .collect();
        hits.sort_by_key(|&(idx, _)| idx);
        hits.into_iter().map(|(_, tx)| tx).collect()
    }

    /// Finds the transcript on `chrom` with minimum `|center - tss_pos|` (optionally bounded by
    /// `max_dist` for `--compat rose`).
    #[must_use]
    pub fn closest_tss(
        &self,
        chrom: &str,
        center: u64,
        max_dist: Option<u64>,
    ) -> Option<&Transcript> {
        let entries = self.by_chrom_tss.get(chrom)?;
        if entries.is_empty() {
            return None;
        }
        let pos_idx = entries.partition_point(|&(pos, _)| pos < center);
        let mut best: Option<(u64, usize)> = None;

        for cand_slot in [pos_idx.checked_sub(1), Some(pos_idx)] {
            if let Some(mut k) = cand_slot {
                while k < entries.len() {
                    let (pos, tx_idx) = entries[k];
                    let dist = pos.abs_diff(center);
                    if let Some((best_d, _)) = best {
                        if dist > best_d {
                            break;
                        }
                    }
                    if max_dist.is_none_or(|limit| dist <= limit) {
                        match best {
                            None => best = Some((dist, tx_idx)),
                            Some((best_d, best_idx))
                                if dist < best_d || (dist == best_d && tx_idx < best_idx) =>
                            {
                                best = Some((dist, tx_idx));
                            }
                            _ => {}
                        }
                    }
                    if cand_slot == pos_idx.checked_sub(1) {
                        if k == 0 {
                            break;
                        }
                        k -= 1;
                    } else {
                        k += 1;
                    }
                }
            }
        }
        best.map(|(_, tx_idx)| &self.transcripts[tx_idx])
    }
}
