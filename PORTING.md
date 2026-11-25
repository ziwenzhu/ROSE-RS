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

| Validation target | Test | What is verified |
|---|---|---|
| §10.1 `--compat rose` parity | `crates/rose-cli/tests/golden.rs::compat_rose_region_map_matches_rose_reference` | Stitched regions, names and ranking/control signal of `run --compat rose` equal `tests/golden/expected_rose_region_map.tsv` (relative error < 1e-6). The expected table comes from `tests/golden/rose_reference.py`, an independent re-implementation of ROSE 0.1's stitching and `ROSE_bamToGFF.py` counting rules, not from running the original Python 2 / R / samtools pipeline. |
| §10.3 Signal kernel | `rose-io` `bam::tests`, `rose-core` `signal::tests`; `golden.rs::nested_region_does_not_steal_reads_from_enclosing_region` | Overlap scan with nested regions; legacy strict-interior counting and floor handling |
| §10.3 CRAM | `golden.rs::cram_with_reference_matches_bam`, `cram_without_reference_is_a_clean_error` | CRAM decoded with `--reference` equals BAM; missing reference is an error, not a panic |
| §10.2 Stitching | `rose-core` `stitch::tests` | Distance stitching and naming, single-peak multi-TSS un-stitching, duplicate IDs rejected |
| Cutoff | `rose-core` `cutoff::tests` | Exact scan and R Brent port on simple curves |

Not yet covered by tests in this repository: a run of the original ROSE scripts (M0), bigWig (M6), `rose multi` (M7), stability (M8), copy-number correction (M9) and public CUT&RUN data (M10). Regenerate the fixture with `python3 tests/golden/make_fixture.py tests/golden` (needs `pysam`).

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
- **Depth floor summation (`ROSE_bamToGFF.py:113, 135`):** With `floor = 1`, bases with `depth > 1` (`depth >= 2`) contribute their **full `depth`** (not `depth - 1`). Clean mode uses the same rule.
- **Read and region coordinates (`ROSE_utils.py::readsToLoci`, `ROSE_bamToGFF.py:115`):** each read is the closed locus `[POS, POS + len(SEQ)]` (`len + 1` bases) before the 200 bp extension, and only bases strictly inside the GFF locus (`start < x < end`) are counted, while density is divided by `end - start`. Because this is not additive across stitch pieces, `run --compat rose` quantifies the active stitched regions directly; the sensitivity sweep and stability tiers still use the pieces.
- **4-decimal intermediate rounding (`ROSE_bamToGFF.py:40, 136`):** `MMR = round(idxstats_mapped / 1e6, 4)` and `density = round((bp_sum / width) / MMR, 4)` before multiplying by `width`.
- **Brent Step-Function Cutoff (`ROSE_callSuper.R:6-31, 159`):** Uses `slope = (max - min) / n`, golden-section + parabolic interpolation (`stats::optimize` on `[1, n]`, `tol = eps^0.25`), `floor(minimum)` index, and strict inequality `signal > cutoff`.
- **First-Row RefSeq Deduplication & `chrY` Drop:** Keeps the first row encountered per RefSeq ID in `*.ucsc` and strips `chrY` regions prior to ranking.

---

## 5. Verification & Performance Benchmarks

### End-to-End Scale Benchmark

Measured by the author on `target/release/rose.exe` (`--threads 8`, pure-Rust `noodles 0.119` + `flate2` `zlib-rs` backend, thin LTO) across `20,020` constituent peaks (`22` human autosomes spanning `2.97 Gb`), `30,000,014` mapped ChIP alignments + `30,000,014` mapped Input alignments (`60M` total BAM records), full `58,350`-transcript human `hg38_refseq.ucsc` annotation, and full table + BED + SVG/PNG plot generation:

| Workload (`0..=25 kb` stitch sensitivity sweep active on all runs) | Target | Measured Wall Time (`8` threads) | Peak Working Set RSS |
|---|---|---|---|
| CUT&RUN (`8M` reads, `15k` peaks, `--assay cutrun -s auto -t 2500`) | — | **`2.20 s`** | **`68.05 MiB`** |
| ChIP-seq single track (`30M` reads, `20,020` peaks, `-t 2500`) | — | **`5.71 s`** | **`74.23 MiB`** |
| ChIP + Input (`30M` + `30M` reads, `20,020` peaks, `-s auto -t 2500`) | `< 15.0 s` | **`9.41 s`** | **`73.94 MiB`** |
| ChIP + Input + `--cn-correct` + `--stability 200` (`5`-tier `xxh3` downsampling) | `< 30.0 s` | **`9.55 s`** | **`79.91 MiB`** |
| `rose multi` 4-sample cohort (`240M` total reads, `-s auto --min-peak-samples 2`) | — | **`35.18 s` (`8.80 s / sample`)** | **`112.87 MiB`** |

*Note on `noodles` (`zlib-rs`) vs optional C `htslib`:* Because pure-Rust `noodles` + `zlib-rs` + BGZF seek-vs-stream auto-switching (`<= 512` coalesced windows) already achieves **`7.91 s` / `8.22 s`** (`~44–52 MiB` RSS) with zero C compiler or FFI `unsafe` code, the optional `htslib` C binding is unnecessary by empirical benchmark data, preserving `#![forbid(unsafe_code)]` across the entire workspace.

### Quality Gates Executed

- `cargo fmt --all -- --check` — `0` formatting diffs
- `cargo clippy --offline --workspace --all-targets -- -D warnings` — `0` warnings (`#![forbid(unsafe_code)]`, `clippy::unwrap_used = "deny"`)
- `cargo test --workspace` — unit tests in `rose-core` and `rose-io`, plus the CLI golden tests in `crates/rose-cli/tests/golden.rs`
- `RUSTDOCFLAGS="-D warnings" cargo doc --offline --workspace --no-deps` — `0` rustdoc warnings
