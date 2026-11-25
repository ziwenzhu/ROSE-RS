---
name: rose-rs
description: >-
  Calls, ranks, compares, and maps super-enhancers and typical enhancers from
  ChIP-seq, CUT&RUN, CRAM, or bigWig tracks using the `rose` (`rose-rs`) CLI.
  Performs 12.5 kb constituent stitching, TSS promoter exclusion (`contained` or
  LILY `subtract`), paired-end fragment quantification, amplicon copy-number
  correction (`--cn-correct` and `--cn-segments`), bootstrap/downsampling cutoff
  stability (`--stability`), multi-sample replicate cohort comparison (`rose multi`),
  enhancer-to-gene assignment, and Young Lab ROSE 0.1 parity (`--compat rose`).
  Use whenever a task involves super-enhancer calling, ROSE, LILY, hockey-stick
  enhancer ranking, enhancer-gene mapping, or differential super-enhancer analysis.
---

# `rose-rs` Agent Skill (`rose` CLI)

## Overview

`rose` (`rose-rs`) is a self-contained, pure-Rust replacement for Young Lab **ROSE** and **LILY**. It requires **no** Python, R, `samtools`, or C libraries at runtime, natively handles paired-end fragments, BAM/CRAM/BigWig tracks, copy-number amplicons, bootstrap stability, and multi-sample cohorts, and writes all tables, UCSC BED tracks, SVG/PNG hockey-stick plots, and `run.json` provenance atomically.

---

## 1. Binary Resolution & Pre-Flight Checks

1. **Locate `rose`:**
   - Check `rose --version` on `$PATH`, then `<repo>/target/release/rose` (`rose.exe` on Windows).
   - If not yet built, run `cargo build --release` in the `ROSE-RS` repository root.
2. **Verify BAM/CRAM indexes before invoking indexed modes:**
   - Standard `--norm-denom filtered` (default) and `--cn-correct` can stream unindexed BAMs sequentially, **but** parallel per-chromosome `rayon` acceleration and `--norm-denom index` / `--scale-factor` window queries require a companion `.bam.bai` (or `.bai` / `.csi`).
3. **Never make implicit network assumptions:**
   - `rose run` and `rose multi` perform zero HTTP calls. Pass `-a <refseq.ucsc|genes.gtf>` explicitly (or pre-install local UCSC `refGene` files via `rose fetch-annotation -g hg38 --from-dir <dir> -o annotations/`).

---

## 2. Protocol Selection & Biological Guardrails

Always match flags to the experimental assay and genome state:

### A. Standard ChIP-seq (`--assay chip`, default)
- **Promoter exclusion (`-t 2500`):** The CLI default `--tss` is `0` (matching legacy ROSE's raw flag default), **but** the standard Whyte et al. super-enhancer protocol excludes `±2,500 bp` around TSSs so promoter H3K27ac peaks are not called as enhancers unless stitched into a broader domain. **Always pass `-t 2500 -a <annotation>`** unless the user explicitly asks to include isolated promoters.
- **Full GENCODE / Ensembl GTFs:** Full GTFs contain ~250k lncRNA and pseudogene TSSs that will over-trigger `--max-tss-genes 2` un-stitching. Pass **`--gene-biotype protein-coding`** when `-a` is an unfiltered GENCODE/Ensembl GTF (keep `--gene-biotype all` for UCSC RefSeq `*.ucsc` tables).
- **Blacklist:** Pass `--blacklist <blacklist.bed>` whenever ENCODE blacklist coordinates are available for the assembly.

### B. Paired-End CUT&RUN (`--assay cutrun`)
- Pass **`--assay cutrun`**, which automatically sets:
  - `--floor 0` (avoids zeroing shallow library signal),
  - `--fragment auto` (`[POS, POS + TLEN)` proper-pair fragments counted once from the leftmost mate),
  - `--subtract-control off` (sparse IgG background belongs at peak calling, not raw per-bp subtraction).
- **Always supply `--blacklist <blacklist.bed>`** for CUT&RUN; unmasked telomeric/centromeric repeats otherwise rank #1.
- **IgG sensitivity check:** If an IgG BAM is available (`-c igg.bam`), pass **`--compare-control`** (keeps primary calls unsubtracted while recording concordance against IgG-subtracted calls in `run.json`).
- **Fragment size tuning:** Adjust `--min-fragment` / `--max-fragment` to match biology (default `0..=1000`; e.g., `--min-fragment 120 --max-fragment 1000` for nucleosomal H3K27ac, or `--max-fragment 120` for sub-nucleosomal transcription factor footprints).

### C. Copy-Number Amplified Genomes / LILY Workflow
Focal amplifications (`MYCN`, `MYC`, `ERBB2`, ecDNA) inflate ChIP and Input equally. Choose one of three routes:
1. **Control BAM binning (`--cn-correct`):** When a whole-cell Input BAM has sufficient depth (`>= 50` median post-filter reads per `100 kb` bin), pass `-c input.bam --cn-correct`. `rose` subtracts background *before* dividing by local copy-number ratio: `corrected = max(0, treat - ctrl) / cn_ratio`, leaving unmappable/centromeric gap bins at `1.0`.
2. **External CNV segments (`--cn-segments <file>`):** Required when no deep Input BAM exists (e.g., **always use `--cn-segments` for CUT&RUN**, because sparse IgG triggers `InsufficientControlDepth` under `--cn-correct`). Accepts 6-column CBS `.seg` (`ID, chrom, start, end, num_mark, log2_ratio`, converted via `2^seg.mean`) or 4-column BED (`chrom, start, end, linear_cn_ratio`, where neutral diploid = `1.0`).
3. **Pre-corrected BigWig tracks (`--signal bigwig`):** For HMCan / deepTools continuous density tracks, pass `-r treat.bw [-c ctrl.bw] --signal bigwig --promoter-mode subtract -t 2500`. Use `--drop-top N` only when extreme unmasked artifact spikes compress the tangent slope, and `--min-se N` only if mimicking LILY's minimum-SE floor.

### D. Multi-Sample & Replicate Cohorts (`rose multi`)
- **Do not** run `rose run` separately on replicates/conditions and intersect heterogeneous stitched intervals ad hoc.
- Write a tab-separated sample sheet (`sample`, `group`, `peaks`, `bam`, optional `control`, optional `scale_factor`) and invoke **`rose multi`**.
- Pass **`--min-peak-samples 2`** when replicates exist so singleton peak-caller noise cannot chain adjacent TAD enhancers into megablocks.
- `rose multi` normalizes across samples via DESeq2 median-of-ratios over shared positive enhancers (`>= 80%` positive samples) unless per-sample `scale_factor` values (e.g., *E. coli* / Drosophila spike-in factors) are supplied.

### E. Stitch-Distance Sensitivity & `--stitch auto`
- **Default `12,500 bp` + automatic `0..=25 kb` sensitivity report:** Every `rose run` and `rose multi` run quantifies elementary peak and inter-peak gap pieces once, making multi-distance evaluation free (zero extra BAM reads), and writes `<prefix>_Stitch_sensitivity.tsv` across `0, 1, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25 kb`. Check `MEDIAN_SE_PEAK_COVERAGE` and `BASE_JACCARD_VS_12_5KB` first to see whether stitch choice actually changes calls!
- **When to use `-s auto` / `--stitch auto`:** Recommended when peaks come from **broad peak callers** (such as MACS3 `--broad` in CUT&RUN, where sub-peaks are already partially merged so `12.5 kb` over-stitches) or on **compact / non-mammalian genomes** (*Drosophila*, *C. elegans*, yeast) and sparse/dense factors. `--stitch auto` estimates the clustering scale strictly from the excess of short inter-peak gaps over the unclustered Poisson gap survival tail (never touching signal -> zero circular tuning), clamps to `[2 kb, 25 kb]`, falls back to `12,500 bp` when `< 100` gaps or `< 5%` short-gap excess exists, and in `rose multi` selects **one shared distance** from the pooled consensus peaks.

### F. Legacy Young Lab ROSE 0.1 Reproduction (`--compat rose`)
- Pass **`--compat rose`** *only* when regression-testing against legacy Python/R ROSE outputs (enables `idxstats` denominator, single-read `SEQ + 200 bp` extension, counting `0x400` duplicates, 4-decimal MMR/density rounding, `chrY` drop, and R Brent `optimize` step-function cutoff).

---

## 3. Canonical Command Templates

### Full Single-Sample ChIP-seq Run (with Stability Bootstrap)
```bash
rose --threads 8 run \
  -i H3K27ac_peaks.narrowPeak \
  -r H3K27ac.bam \
  -c Input.bam \
  -a hg38_refseq.ucsc \
  --blacklist hg38-blacklist.v2.bed \
  -t 2500 \
  --stability 500 \
  --seed 42 \
  -o results/H3K27ac
```

### Paired-End CUT&RUN Run
```bash
rose --threads 8 run \
  --assay cutrun \
  -i H3K27ac_cutrun.broadPeak \
  -r H3K27ac_cutrun.bam \
  -c IgG_cutrun.bam \
  --compare-control \
  -a genes.gtf \
  --gene-biotype protein-coding \
  --blacklist hg38-blacklist.v2.bed \
  -t 2500 \
  -o results/cutrun_H3K27ac
```

### Multi-Sample Differential / Replicate Cohort (`rose multi`)
`samples.tsv`:
```tsv
sample	group	peaks	bam	control	scale_factor
DMSO_1	DMSO	dmso1.narrowPeak	dmso1.bam	input1.bam	
DMSO_2	DMSO	dmso2.narrowPeak	dmso2.bam	input2.bam	
JQ1_1	JQ1	jq1_1.narrowPeak	jq1_1.bam	input1.bam	
JQ1_2	JQ1	jq1_2.narrowPeak	jq1_2.bam	input2.bam	
```
*(Omit the `scale_factor` column or leave cells blank for automatic DESeq2 median-of-ratios scaling; populate numeric spike-in multipliers when calibrated.)*

```bash
rose --threads 8 multi \
  --samples samples.tsv \
  --min-peak-samples 2 \
  -a hg38_refseq.ucsc \
  --blacklist hg38-blacklist.v2.bed \
  -t 2500 \
  -o results/cohort_JQ1_vs_DMSO
```

### Modular Step-by-Step Execution
```bash
rose stitch -i peaks.bed -a hg38_refseq.ucsc -t 2500 -o work/stitched.gff
rose signal -i work/stitched.gff -r chip.bam -c input.bam -o work/region_map.tsv
rose call   -i work/region_map.tsv -n SAMPLE1 --stability 500 -o work/calls
rose genes  -i work/calls/SAMPLE1_SuperEnhancers.table.txt -a hg38_refseq.ucsc -o work/genes
```

---

## 4. Output Files & Column Schemas

`<prefix>` defaults to the file stem of `-i` (or `-n <name>` in `rose call`):

| File | Description & Exact Columns |
|---|---|
| `<prefix>_AllEnhancers.table.txt` & `<prefix>_SuperEnhancers.table.txt` | 4 `#` metadata header lines, followed by TSV columns: `REGION_ID`, `CHROM`, `START`, `STOP`, `NUM_LOCI`, `CONSTITUENT_SIZE`, `<rankby_name>`, `[<control_name>]`, `enhancerRank` (`1` = highest), `isSuper` (`1` or `0`), `[CLASS]` (`super`/`typical` in clean mode), `[RAW_NET_SIGNAL, CORRECTED_SIGNAL, CN_RATIO]` (when CN active), `[BOOTSTRAP_FREQ, BORDERLINE]` (when `--stability > 0`). Coordinates are 0-based half-open in `--compat clean` and GFF-style in `--compat rose`. |
| `<prefix>_Enhancers_withSuper.bed`, `<prefix>_Gateway_Enhancers.bed`, `<prefix>_Gateway_SuperEnhancers.bed` | UCSC 0-based half-open BED tracks with valid integer scores in `0..=1000` (`1000` = rank 1). Safe to re-run in place (written via `AtomicWriter` without duplicate append). |
| `<prefix>_SuperEnhancers_ENHANCER_TO_GENE.txt` *(and `_AllEnhancers_`)* | `REGION_ID`, `CHROM`, `START`, `STOP`, `NUM_LOCI`, `CONSTITUENT_SIZE`, `<signal>`, `enhancerRank`, `isSuper`, `OVERLAP_GENES`, `PROXIMAL_GENES`, `CLOSEST_GENE`. |
| `<prefix>_SuperEnhancers_GENE_TO_ENHANCER.txt` *(and `_AllEnhancers_`)* | `GENE_NAME`, `REFSEQ_ID`, `PROXIMAL_ENHANCERS`, `ENHANCER_RANKS`, `IS_SUPER`. |
| `<prefix>_Plot_points.svg`, `<prefix>_Plot_points.png`, `<prefix>_Plot_points.tsv` | Publication-ready vector SVG + raster PNG hockey-stick curves and underlying `(RANK_ASCENDING, SIGNAL, IS_SUPER, CUTOFF)` coordinates. |
| `<prefix>_Stitch_sensitivity.tsv` | Emitted on every `run` / `multi` (unless `--no-stitch-report`): `STITCH_BP`, `IS_ACTIVE`, `NUM_REGIONS`, `NUM_SUPER`, `CUTOFF_SIGNAL`, `MEDIAN_SE_WIDTH_BP`, `P90_SE_WIDTH_BP`, `MEDIAN_SE_PEAK_COVERAGE`, `MEDIAN_ALL_PEAK_COVERAGE`, `BASE_JACCARD_VS_12_5KB`, `SE_OVERLAP_FRAC_VS_12_5KB`. |
| `<prefix>_Stability_saturation.tsv` | Emitted when `--stability > 0`: `READ_FRACTION` (`0.25..1.00`), `CUTOFF_SIGNAL`, `SUPER_ENHANCER_COUNT`, `JACCARD_VS_FULL`. |
| `multi_signal_matrix.tsv` & `multi_consensus_calls.tsv` | Emitted by `rose multi`: normalized signal matrix across all samples, per-sample `*_isSuper`, per-group `*_replicates` & `*_consensus`, plus `DIFF_STATUS` (`gained`, `lost`, `shared`, `typical`) and `LOG2_FC` (`log2((mean_group2 + 1) / (mean_group1 + 1))`) when 2 groups exist. |
| `run.json` | Machine-readable provenance: `version`, `command`, `arguments`, SHA-256 `input_checksums`, `ranking_read_stats` `(index_mapped, filtered_units, effective_mmr)`, `control_read_stats`, `cutoff_signal`, `num_regions`, `num_super_enhancers`, `control_comparison_note`, `active_stitch_bp`, and optional `auto_stitch_report`. |

---

## 5. Post-Run Verification & Troubleshooting

After running `rose`, always inspect **`run.json`** and **`*_Stitch_sensitivity.tsv`**:
1. Exit code is `0` and `num_regions > 0`.
2. `num_super_enhancers` is biologically plausible (typically `2%`–`8%` of stitched enhancer regions, e.g., `200`–`900` SEs on human H3K27ac). If `num_super_enhancers` is `1`–`3` with an astronomically high `cutoff_signal`, an unmasked artifact or high-copy amplicon is compressing the slope — re-run with `--blacklist`, `--cn-correct` / `--cn-segments`, or `--drop-top 1`.
3. Inspect `*_Stitch_sensitivity.tsv`: if `MEDIAN_SE_PEAK_COVERAGE` drops below `0.20` at `12,500 bp` (typical on MACS3 `--broad` CUT&RUN peaks), over-stitching is absorbing mostly background gap — compare against `5,000 bp` or re-run with `--stitch auto`.
4. When `--stability` is enabled, check `JACCARD_VS_FULL` at `0.75` in `*_Stability_saturation.tsv` (`>= 0.85` indicates sequencing depth saturation) and flag `BORDERLINE == 1` enhancers (`0.10 <= BOOTSTRAP_FREQ <= 0.90`).

### Common Errors & Fixes

| Error Message | Root Cause | Fix |
|---|---|---|
| `ContigMismatch` (`< 50%` regions matched) | Peak file uses non-standard scaffolds or wrong assembly relative to BAM/BigWig (`chr1` <-> `1` is already reconciled automatically) | Check assembly build (`hg19` vs `hg38` vs `mm10`) between peak BED and BAM `@SQ` header |
| `MissingIndex` | `--norm-denom index` or `--scale-factor` requested on an unindexed BAM | Index the BAM (`samtools index`) or use default `--norm-denom filtered` sequential streaming |
| `DuplicateRegionId` | Input peak file repeats the same column-4 (BED) or column-2 (GFF) identifier | Deduplicate or strip column 4 so `rose` auto-assigns unique `PEAK_<k>` IDs |
| `InsufficientControlDepth` | `--cn-correct` used on a shallow control BAM (`< 50` median reads per `100 kb` bin, such as CUT&RUN IgG) | Drop `--cn-correct` and pass external segments via `--cn-segments <file.seg>` instead |
| `--tss > 0 requires -a <annotation>` | Promoter exclusion requested without a transcript table | Pass `-a <refseq.ucsc|genes.gtf>` (or `-t 0` if intentionally ranking without TSS exclusion) |
