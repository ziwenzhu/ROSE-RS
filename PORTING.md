# Porting Notes (`PORTING.md`)

This document maps original Young Lab **ROSE** (`bitbucket.org/young_computation/rose`, Python 2.7 + R) and **LILY** conceptual features to `rose-rs` (`rose-core`, `rose-io`, `rose-cli`), records deliberate improvements vs. preserved `--compat rose` quirks, and summarizes verification results.

---

## 1. Module & API Mapping

| Legacy ROSE Source | Original Function / Class | `rose-rs` Equivalent | Algorithmic Notes |
|---|---|---|---|
| `ROSE_utils.py` | `Locus`, `LocusCollection` | `rose_core::Interval`, `rose_core::Constituent`, `rose_core::StitchedRegion` | Canonical 0-based half-open `[start, end)` `Copy` value type with checked construction (`Interval::try_new`), $O(\log N)$ sorted slices instead of Python bin-hash dicts |
| `ROSE_utils.py` | `importRefseq`, `makeStartDict`, `makeTranscriptCollection` | `rose_io::read_annotation_index`, `rose_core::AnnotationIndex` | Supports both UCSC `refGene` (`*.ucsc`) and GTF/GFF3 (including `exon`-only UCSC GTFs); primary chromosomes preferred over `_alt`/`_random`/`chrUn` in clean mode |
| `ROSE_main.py` | `regionStitching`, `LocusCollection.stitchCollection` | `rose_core::stitch_constituents` | Single-pass coordinate sort + sweep merge (`<= stitch_distance`), blacklist removal, TSS `contained` & `subtract` modes, and `> max_tss_genes` disaggregation by distinct gene symbol |
| `ROSE_bamToGFF.py` & `ROSE_main.py` | `mapToGFF`, `Bam.getReadsLocus`, `mapCollection` | `rose_io::compute_bam_signals`, `rose_core::SignalScratch`, `rose_core::normalize_coverage` | In-process `noodles` BAM/CRAM decoding across `rayon` workers; $O(1)$ clamped overlap (`floor = 0`) or reusable per-thread difference array (`floor > 0`); proper-pair `[POS, POS + TLEN)` fragments |
| *(LILY concept)* | HMCan bigWig scoring & copy-number normalization | `rose_io::BigWigFile`, `rose_core::CopyNumberProfile` | Pure-Rust UCSC BigWig v4 (`0x888FFC26`) R+ tree reader with exact width-weighted integration; fused single-pass 100 kb control BAM binning (`--cn-correct`) and `.seg`/BED support (`--cn-segments`) |
| `ROSE_callSuper.R` | `calculate_cutoff`, `numPts_below_line`, `writeSuperEnhancer_table`, `convert_stitched_to_bed` | `rose_core::call_super_enhancers`, `rose_core::compute_cutoff_on_slice`, `rose_io::write_enhancer_call_tables`, `rose_io::write_enhancer_beds`, `rose_cli::plot` | Exact $O(N)$ chord-tangent scan (plus exact R Brent `optimize` port under `CutoffMethod::RoseCompat`); static-safe `plotters` SVG + PNG with embedded `DejaVuSans.ttf`; atomic BED/TSV writes |
| `ROSE_geneMapper.py` | `mapEnhancerToGene` | `rose_core::map_enhancers_to_genes`, `rose_io::write_gene_tables` | $O(M \log N)$ binary search over sorted transcript bodies and TSS arrays; typed headers; preserves all `'R'`-prefixed IDs and gene symbols |
| *(New in `rose-rs`)* | Multi-sample cohort & cutoff stability | `rose_core::evaluate_multi_cohort`, `rose_core::evaluate_stability` | `--min-peak-samples` consensus peak merging, `>= 80%` positive-region DESeq2 median-of-ratios scaling, seeded SplitMix64 region bootstrap, and `xxh3_64` multi-tier read downsampling |
| *(New in `rose-rs`)* | Zero-I/O stitch-distance sweep & `--stitch auto` | `rose_core::PiecewiseStitchIndex`, `rose_core::estimate_auto_stitch_distance`, `rose_io::write_stitch_sensitivity_tsv` | Strictly disjoint `d = 0` peak + inter-peak gap decomposition (`0..=25 kb` sensitivity table on every run with zero extra BAM I/O) + signal-blind ECDF short-gap excess estimator clamped to `[2 kb, 25 kb]` with `12.5 kb` fallback |

---

## 2. Test & Validation Mapping

| Validation Target (`PLAN.md` §10) | `rose-rs` Test Suite | What Is Verified |
|---|---|---|
| §10.1 Golden `--compat rose` parity (`M0`–`M4`) | `crates/rose-io/tests/golden_and_pipeline_tests.rs::test_m0_m4_golden_rose_compat_end_to_end` | Exact parity against Young Lab ROSE logic + unmodified R 4.6.1 `ROSE_callSuper.R` on `tests/golden/` (`28` regions, `PROMO_PEAK_1` TSS drop, `MULTITSS` disaggregation, `chrY` drop, 4-decimal signals, `762650.7858` Brent cutoff, `4` super-enhancers, `MYC` & `RERE` gene mapping) |
| §10.2 Coordinate, stitching & difference-array invariants | `crates/rose-core/tests/property_tests.rs` & `test_m1_bed_vs_gff_equivalence` | `proptest` suites comparing `SignalScratch` + multi-tier downsampling against a brute-force per-base oracle, connected-component stitching invariants, base-exact `PromoterMode::Subtract`, and 0-based BED == 1-based GFF equivalence |
| §10.3 & §10.4 Paired-end BAM, CIGAR & CUT&RUN | `test_m2_m10_paired_end_fragments_cigar_and_cutrun_filters` & `cli_e2e_tests.rs` | `750 bp` proper pair starting `550 bp` upstream counted once (`200 bp` overlap), `1400 bp` & `0x400` duplicate excluded, `36M4D10S` CIGAR span (`40 bp`), `@PG` dedup check, indexed vs streaming parity, `--compare-control` |
| §10.4 BigWig & LILY outlier options (`M6`) | `test_m6_bigwig_width_weighted_and_outlier_controls` | Binary `.bw` roundtrip verifying width-weighted `900*2 + 100*50 = 6800` vs unweighted `26.0`, `--drop-top 1` artifact protection, and `--min-se` |
| §10.5 Multi-sample cohort (`M7`) | `test_m7_multi_sample_consensus_deseq2_and_self_comparison` & `cli_e2e_tests.rs` | 1-sample `multi == run` identity, noisy multi-sample consensus merging without count inflation, zero-inflated DESeq2 median-of-ratios (`2.0x`), and `0` gained / `0` lost on self-comparison |
| §10.6 Stability report (`M8`) | `test_m8_stability_determinism_and_downsampling` | Bit-exact determinism across runs per `--seed` (`SplitMix64` + `xxh3_64`) across both `floor = 0` and `floor = 1` |
| §10.7 Copy-number correction (`M9`) | `test_m9_copy_number_correction_recovers_true_ranking` | Synthetic `5x` ChIP + Input focal amplification + `0`-read centromeric gap bins recovering true ranking via `--cn-correct` and `.seg` `--cn-segments`, plus `InsufficientControlDepth` safeguard on shallow IgG |
| §10.8 Edge cases & atomic I/O | `test_section_10_8_edge_cases` & `cli_e2e_tests.rs` | Flat `0`-SE signal, all-zero plateau (`1` SE), `chr1` <-> `1` Ensembl reconciliation, `< 50%` `ContigMismatch`, `MissingIndex`, `DuplicateRegionId`, exon-only UCSC GTF fallback, `0..=1000` BED scores, and Gateway BED rerun idempotency |

---

## 3. Deliberate Differences (Default `--compat clean` Mode)

1. **Paired-End Fragment Spans & QC Filtering:** Proper pairs (`0x2`) emit `[POS, POS + TLEN)` once from the canonical leftmost mate (`TLEN` within `0..=1000`), and alignments with `0x4 | 0x100 | 0x200 | 0x400 | 0x800` (unmapped, secondary, QC-fail, PCR duplicate, supplementary) are excluded by default, scaling by million post-filter fragments/reads (`--norm-denom filtered`).
2. **Exact $O(N)$ Tangent Cutoff:** Default `--cutoff exact` minimizes `y[i] - s * i` (`s = (y[n] - y[1]) / (n - 1)`) in one linear pass with rightmost tie-breaking across zero-signal plateaus, avoiding R `optimize()` local-minimum traps on step functions.
3. **Contig Name Reconciliation & Validation:** `chr1` and `1` are reconciled automatically against BAM/CRAM/BigWig headers; if fewer than 50% of input regions match any contig, `rose` aborts with a non-zero `ContigMismatch` error instead of writing silent all-zero tables.
4. **Typed Gene Tables & Atomic Output Writes:** No enhancers or genes starting with `'R'` (`REG_*`, `RERE`, `RUNX1`) are dropped; all output files are written via `AtomicWriter` (temp file + atomic rename), preventing partial files and `append=TRUE` Gateway BED row duplication on reruns; BED column 5 scores are scaled to valid UCSC `0..=1000`.
5. **Primary-Assembly RefSeq Preference:** When the same transcript ID appears on both `_alt`/`_random` patches and primary chromosomes (`chr1..chr22, chrX, chrY`), clean mode retains the primary chromosome coordinate.

---

## 4. Preserved Quirks (Behind `--compat rose`)

Do **not** alter these behaviors in `CompatMode::RoseLegacy` — golden regression tests depend on them:
- **Post-stitching multi-TSS disaggregation (`ROSE_main.py:101, 113`):** Checks `±50 bp` around `txStart` (`+`) / `txEnd` (`-`) and un-stitches strictly when `> 2` (`>= 3`) distinct gene symbols overlap (even though ROSE's inline comment says "2").
- **Depth floor summation (`ROSE_bamToGFF.py:113, 135`):** With `floor = 1`, bases with `depth > 1` (`depth >= 2`) contribute their **full `depth`** (not `depth - 1`).
- **4-decimal intermediate rounding (`ROSE_bamToGFF.py:40, 136`):** `MMR = round(idxstats_mapped / 1e6, 4)` and `density = round((bp_sum / width) / MMR, 4)` before multiplying by `width`.
- **Brent Step-Function Cutoff (`ROSE_callSuper.R:6-31, 159`):** Uses `slope = (max - min) / n`, golden-section + parabolic interpolation (`stats::optimize` on `[1, n]`, `tol = eps^0.25`), `floor(minimum)` index, and strict inequality `signal > cutoff`.
- **First-Row RefSeq Deduplication & `chrY` Drop:** Keeps the first row encountered per RefSeq ID in `*.ucsc` and strips `chrY` regions prior to ranking.

---

## 5. Verification & Performance Benchmarks

### End-to-End Scale Benchmark (`PLAN.md` §3 & M5)

Measured on `target/release/rose.exe` (`--threads 8`, pure-Rust `noodles 0.119` + `flate2` `zlib-rs` backend, thin LTO) across `20,020` constituent peaks (`22` human autosomes spanning `2.97 Gb`), `30,000,014` mapped ChIP alignments + `30,000,014` mapped Input alignments (`60M` total BAM records), full `58,350`-transcript human `hg38_refseq.ucsc` annotation, and full table + BED + SVG/PNG plot generation:

| Workload (`0..=25 kb` stitch sensitivity sweep active on all runs) | `PLAN.md` §3 Target | Measured Wall Time (`8` threads) | Peak Working Set RSS |
|---|---|---|---|
| CUT&RUN (`8M` reads, `15k` peaks, `--assay cutrun -s auto -t 2500`) | — | **`2.20 s`** | **`68.05 MiB`** |
| ChIP-seq single track (`30M` reads, `20,020` peaks, `-t 2500`) | — | **`5.71 s`** | **`74.23 MiB`** |
| ChIP + Input (`30M` + `30M` reads, `20,020` peaks, `-s auto -t 2500`) | `< 15.0 s` | **`9.41 s`** | **`73.94 MiB`** |
| ChIP + Input + `--cn-correct` + `--stability 200` (`5`-tier `xxh3` downsampling) | `< 30.0 s` | **`9.55 s`** | **`79.91 MiB`** |
| `rose multi` 4-sample cohort (`240M` total reads, `-s auto --min-peak-samples 2`) | — | **`35.18 s` (`8.80 s / sample`)** | **`112.87 MiB`** |

*Note on `noodles` (`zlib-rs`) vs optional C `htslib` (`PLAN.md` M5):* Because pure-Rust `noodles` + `zlib-rs` + BGZF seek-vs-stream auto-switching (`<= 512` coalesced windows) already achieves **`7.91 s` / `8.22 s`** (`~44–52 MiB` RSS) with zero C compiler or FFI `unsafe` code, the optional `htslib` C binding is unnecessary by empirical benchmark data, preserving `#![forbid(unsafe_code)]` across the entire workspace.

### Criterion Micro-Benchmarks (`crates/rose-core/benches/signal_and_stitch.rs`)

| Kernel | Median Time |
|---|---|
| `signal_floor0_5k_reads` ($O(1)$ interval overlap, 5,000 reads) | **`4.09 µs`** |
| `signal_floor1_diffarray_5k_reads` (difference-array prefix sum, 5,000 reads on 25 kb SE) | **`26.78 µs`** |
| `cutoff_exact_15k` (sort + exact chord-tangent scan, 15,000 enhancers) | **`96.63 µs`** |
| `stitch_20k_peaks` (20,000 genome-wide peaks -> stitched loci) | **`27.64 ms`** |

### Quality Gates Executed

- `cargo fmt --all -- --check` — `0` formatting diffs
- `cargo clippy --offline --workspace --all-targets -- -D warnings` — `0` warnings (`#![forbid(unsafe_code)]`, `clippy::unwrap_used = "deny"`)
- `cargo test --offline --workspace` — `15/15` unit, `proptest` property, golden parity, and CLI e2e suites passing
- `RUSTDOCFLAGS="-D warnings" cargo doc --offline --workspace --no-deps` — `0` rustdoc warnings
