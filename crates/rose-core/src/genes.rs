//! Binary-search `Enhancer -> Gene` and `Gene -> Enhancer` mapper (`rose genes`).

use std::collections::{HashMap, HashSet};

use serde::{Deserialize, Serialize};

use crate::cutoff::RankedEnhancer;
use crate::tss::AnnotationIndex;
use crate::types::CompatMode;

/// Configuration for [`map_enhancers_to_genes`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneMapConfig {
    /// Proximal TSS search window in bp (`50_000` by default).
    pub search_window: u64,
    /// Deduplicate `GENE_TO_ENHANCER` rows by gene symbol (`true` in ROSE).
    pub unique_genes: bool,
    /// Compatibility mode (`Clean` searches chromosome-wide for closest TSS; `RoseLegacy` bounds to `1 Mb`).
    pub compat: CompatMode,
}

impl Default for GeneMapConfig {
    fn default() -> Self {
        Self {
            search_window: 50_000,
            unique_genes: true,
            compat: CompatMode::Clean,
        }
    }
}

/// Enhancer-centric gene assignment row (`*_ENHANCER_TO_GENE.txt`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EnhancerGeneRow {
    /// Underlying ranked enhancer call.
    pub enhancer: RankedEnhancer,
    /// Comma-joined unique gene symbols whose transcript body overlaps the enhancer.
    pub overlap_genes: Vec<String>,
    /// Comma-joined unique gene symbols whose TSS lies within `±search_window` (excluding overlapping transcripts).
    pub proximal_genes: Vec<String>,
    /// Gene symbol of the transcript whose TSS is closest to the enhancer midpoint.
    pub closest_gene: Option<String>,
}

/// Gene-centric enhancer summary row (`*_GENE_TO_ENHANCER.txt`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeneEnhancerRow {
    /// Gene symbol (`GENE_NAME`).
    pub gene_name: String,
    /// Representative transcript accession (`REFSEQ_ID`).
    pub transcript_id: String,
    /// Overlapping + proximal enhancer locus strings (`chr:start-end`).
    pub proximal_enhancers: Vec<String>,
    /// Corresponding 1-based enhancer ranks (`ENHANCER_RANKS`).
    pub enhancer_ranks: Vec<usize>,
    /// `true` if at least one associated enhancer is a super-enhancer (`IS_SUPER == 1`).
    pub is_super: bool,
}

/// Maps ranked enhancers to overlapping, proximal (`±search_window`), and closest genes in
/// $O(M \log N)$ time using [`AnnotationIndex`].
#[must_use]
pub fn map_enhancers_to_genes(
    ranked: &[RankedEnhancer],
    annot: &AnnotationIndex,
    cfg: &GeneMapConfig,
) -> (Vec<EnhancerGeneRow>, Vec<GeneEnhancerRow>) {
    let mut enh_rows = Vec::with_capacity(ranked.len());

    // Track insertion order of transcript IDs across enhancers (matching `ROSE_geneMapper.py:167-184`).
    let mut seen_tx_order: Vec<String> = Vec::new();
    let mut seen_tx_set: HashSet<String> = HashSet::new();
    let mut tx_to_gene: HashMap<String, String> = HashMap::new();
    let mut overlap_loci_by_tx: HashMap<String, Vec<String>> = HashMap::new();
    let mut prox_loci_by_tx: HashMap<String, Vec<String>> = HashMap::new();
    let mut ranks_by_tx: HashMap<String, Vec<usize>> = HashMap::new();
    let mut super_by_tx: HashMap<String, bool> = HashMap::new();

    for enh in ranked {
        let reg = &enh.item.region;
        let iv = reg.interval;
        let locus_str = format!("{}:{}-{}", reg.chrom, iv.start(), iv.end());

        let overlap_txs = annot.transcripts_overlapping(&reg.chrom, iv);
        let mut overlap_tx_ids: Vec<&str> = Vec::new();
        let mut overlap_tx_set: HashSet<&str> = HashSet::new();
        for tx in &overlap_txs {
            if overlap_tx_set.insert(tx.transcript_id.as_str()) {
                overlap_tx_ids.push(tx.transcript_id.as_str());
                tx_to_gene
                    .entry(tx.transcript_id.clone())
                    .or_insert_with(|| tx.gene_name.clone());
            }
        }

        let (lo, hi) = match cfg.compat {
            CompatMode::Clean => (
                iv.start().saturating_sub(cfg.search_window),
                iv.end().saturating_sub(1).saturating_add(cfg.search_window),
            ),
            CompatMode::RoseLegacy => (
                iv.start().saturating_sub(cfg.search_window),
                iv.end().saturating_add(cfg.search_window),
            ),
        };
        let prox_txs = annot.transcripts_with_tss_in(&reg.chrom, lo, hi);
        let mut prox_tx_ids: Vec<&str> = Vec::new();
        let mut prox_tx_set: HashSet<&str> = HashSet::new();
        for tx in &prox_txs {
            let tid = tx.transcript_id.as_str();
            if !overlap_tx_set.contains(tid) && prox_tx_set.insert(tid) {
                prox_tx_ids.push(tid);
                tx_to_gene
                    .entry(tx.transcript_id.clone())
                    .or_insert_with(|| tx.gene_name.clone());
            }
        }

        let closest_gene = match cfg.compat {
            CompatMode::Clean => {
                // Consider overlapping + proximal + chromosome-wide closest TSS.
                annot
                    .closest_tss(&reg.chrom, iv.center(), None)
                    .map(|tx| tx.gene_name.clone())
            }
            CompatMode::RoseLegacy => {
                // `ROSE_geneMapper.py:114-147` searches overlapping + proximal + 1 Mb distal TSSs.
                let legacy_center = iv.center();
                let mut best: Option<(u64, &str)> = None;
                let lo_1m = iv.start().saturating_sub(1_000_000);
                let hi_1m = iv.end().saturating_add(1_000_000);
                let distal_txs = annot.transcripts_with_tss_in(&reg.chrom, lo_1m, hi_1m);
                for tx in overlap_txs
                    .iter()
                    .chain(prox_txs.iter())
                    .chain(distal_txs.iter())
                {
                    let d = tx.legacy_rose_tss_coord().abs_diff(legacy_center);
                    if best.is_none_or(|(bd, _)| d < bd) {
                        best = Some((d, tx.gene_name.as_str()));
                    }
                }
                best.map(|(_, g)| g.to_owned())
            }
        };

        let overlap_genes = uniquify_gene_names(&overlap_txs, |tx| &tx.gene_name);
        let proximal_genes = uniquify_gene_names(
            &prox_txs
                .into_iter()
                .filter(|tx| !overlap_tx_set.contains(tx.transcript_id.as_str()))
                .collect::<Vec<_>>(),
            |tx| &tx.gene_name,
        );

        for &tid in &overlap_tx_ids {
            if seen_tx_set.insert(tid.to_owned()) {
                seen_tx_order.push(tid.to_owned());
            }
            overlap_loci_by_tx
                .entry(tid.to_owned())
                .or_default()
                .push(locus_str.clone());
            ranks_by_tx
                .entry(tid.to_owned())
                .or_default()
                .push(enh.enhancer_rank);
            let entry = super_by_tx.entry(tid.to_owned()).or_insert(false);
            *entry = *entry || enh.is_super;
        }

        for &tid in &prox_tx_ids {
            if seen_tx_set.insert(tid.to_owned()) {
                seen_tx_order.push(tid.to_owned());
            }
            prox_loci_by_tx
                .entry(tid.to_owned())
                .or_default()
                .push(locus_str.clone());
            ranks_by_tx
                .entry(tid.to_owned())
                .or_default()
                .push(enh.enhancer_rank);
            let entry = super_by_tx.entry(tid.to_owned()).or_insert(false);
            *entry = *entry || enh.is_super;
        }

        enh_rows.push(EnhancerGeneRow {
            enhancer: enh.clone(),
            overlap_genes,
            proximal_genes,
            closest_gene,
        });
    }

    enh_rows.sort_by_key(|r| r.enhancer.enhancer_rank);

    let mut tx_indices: Vec<usize> = (0..seen_tx_order.len()).collect();
    tx_indices.sort_by_key(|&i| {
        let tid = &seen_tx_order[i];
        let min_r = ranks_by_tx[tid].iter().copied().min().unwrap_or(usize::MAX);
        (min_r, i)
    });

    let mut used_genes: HashSet<String> = HashSet::new();
    let mut gene_rows = Vec::new();
    for idx in tx_indices {
        let tid = &seen_tx_order[idx];
        let gene_name = tx_to_gene.get(tid).cloned().unwrap_or_default();
        if cfg.unique_genes && !used_genes.insert(gene_name.clone()) {
            continue;
        }
        let mut loci = overlap_loci_by_tx.remove(tid).unwrap_or_default();
        if let Some(prox) = prox_loci_by_tx.remove(tid) {
            loci.extend(prox);
        }
        let ranks = ranks_by_tx.get(tid).cloned().unwrap_or_default();
        let is_sup = super_by_tx.get(tid).copied().unwrap_or(false);
        gene_rows.push(GeneEnhancerRow {
            gene_name,
            transcript_id: tid.clone(),
            proximal_enhancers: loci,
            enhancer_ranks: ranks,
            is_super: is_sup,
        });
    }

    (enh_rows, gene_rows)
}

fn uniquify_gene_names<'a, T, F>(items: &[T], mut get_name: F) -> Vec<String>
where
    F: FnMut(&T) -> &'a str,
{
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for item in items {
        let name = get_name(item);
        if !name.is_empty() && seen.insert(name) {
            out.push(name.to_owned());
        }
    }
    out
}
