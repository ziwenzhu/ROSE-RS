# Plan: a fast rewrite of ROSE

A plan to reimplement ROSE (Rank Ordering of Super-Enhancers) as a single compiled
tool, fixing the defects found in the review of
[young_computation/rose](https://bitbucket.org/young_computation/rose) and adopting
the useful ideas from [BoevaLab/LILY](https://github.com/BoevaLab/LILY).

Working name: `rose-rs`. This is the original design document; the implementation now
exists in this repository. Performance figures below were the targets; measured results
and the current test coverage are in `README.md` and `PORTING.md`.

## 1. Goals

- **Same answers as ROSE by default.** Given the same inputs, a compatibility mode
  reproduces the original's tables closely enough to be a drop-in replacement.
- **Much faster.** Minutes-to-hours becomes seconds on a laptop.
- **Fails loudly.** No silent zeros, stale outputs, or success exit codes on error.
- **One self-contained binary.** No Python 2, no R, no `samtools` on `PATH`, no need
  to run from the install directory, and native builds across Linux (`musl`), macOS
  and Windows.

Supported assays: **ChIP-seq and CUT&RUN**. The algorithm is the same for both; they
differ in how signal is counted and normalised (see "Assay presets" in section 5).

Beyond a faster ROSE, three additions (section 6): a multi-sample and replicate
mode, a cutoff stability report, and built-in copy-number correction.

Non-goals: CUT&Tag (its Tn5 accessibility bias needs separate treatment and is out
of scope), peak calling, and statistical testing of differential signal (the signal
matrix is exported for tools such as DESeq2 or edgeR).

## 2. Language and toolchain: Rust

| Need | Choice |
|---|---|
| BAM/CRAM access | `noodles` + `flate2` (`zlib-rs` SIMD backend) by default; optional `rust-htslib` (`libdeflate`) feature behind an internal `AlignmentReader` trait |
| bigWig access | `bigtools` |
| Parallelism | `rayon`, one reader per thread |
| Deterministic hashing | `xxhash-rust` (`xxh3` / `xxh64`), explicitly seeded |
| Command line | `clap` |
| Plotting | `plotters` (`default-features = false`: `svg_backend`, `bitmap_backend`, `ab_glyph`, embedded font) |
| Testing and benchmarks | `cargo test`, `proptest`, `criterion`, `hyperfine` |
| Distribution | Static Linux `musl` binary, macOS and Windows binaries; Bioconda recipe |

The hot loop is "decompress BGZF blocks, decode BAM records, accumulate interval
coverage". Starting with pure-Rust `noodles` + `zlib-rs` keeps the default build
100% safe Rust (`#![forbid(unsafe_code)]`), avoids `hts-sys` / `libcurl` / `openssl`
C toolchain breakage on Windows and static `musl`, and still benefits from SIMD
deflate. Because BAM/CRAM access sits behind a narrow `AlignmentReader` trait in
`rose-io`, an optional `htslib` Cargo feature can be benchmarked in M5 and dropped
if `noodles` already saturates I/O.

`plotters` must disable default features so it does not link system `fontconfig` and
`freetype` on Linux, which would break static `musl` builds. Workspace-wide lints in
the root `Cargo.toml` (`[workspace.lints]`) enforce `unsafe_code = "forbid"` on all
safe crates and `clippy::unwrap_used = "deny"`.

## 3. Where the original loses time, and the fix

| Original behaviour | Cost | Replacement |
|---|---|---|
| One `samtools view` subprocess per region, output parsed as text | Process spawn and SAM text formatting for every region | In-process binary record decoding |
| Per-base coverage held in Python dicts | Millions of hash inserts per large region | Exact overlap arithmetic (`floor = 0`), or a reusable per-thread difference array (`floor > 0`) |
| Each BAM mapped twice (stitched regions and constituents) | 2x the work; the constituent table is never used for ranking | One pass computes both; the constituent table becomes opt-in |
| `samtools flagstat` for the mapped-read count | A full scan of every BAM | `.bai`/`.csi` index metadata under `--compat rose` (with scan fallback for `.crai`); single-pass filtered count in clean modes |
| Background jobs polled every 5 minutes | Any run longer than 10 seconds takes at least 5 minutes | Threads joined directly |
| Everything serial | One core | `rayon` across chromosomes and BAMs, one reader per thread |
| Gene mapper uses `list.count` per annotation row | Quadratic in transcript count | Sorted transcript and TSS arrays with binary search |
| Stitching by repeated dict lookups | Hard to reason about, order-dependent names | Canonical 0-based half-open intervals, sort once, single sweep |
| `optimize()` on a step function for the cutoff | Can stop at a local minimum | Exact linear scan with deterministic tie-breaking (see section 5) |

Read access and denominator strategy:
- **Indexed range mode:** When regions cover a small share of the genome and the
  scaling denominator does not require scanning every read (`--norm-denom index`,
  `--compat rose`, or user `--scale-factor`), merge nearby fetch windows per
  chromosome—widened upstream by `max(extend, max_fragment)` in paired-end mode—and
  read only overlapping BGZF blocks.
- **Single-pass streaming mode:** When regions are dense, when reading CRAM (`.crai`
  stores no `idxstats` read counts), when counting exact post-filter fragments/reads
  (`--norm-denom filtered`), or when `--cn-correct` scans the control BAM, sweep each
  chromosome once in parallel across `rayon` workers so region coverage, filtered
  library size, read-downsampling tiers, and 100 kb copy-number bins are all
  collected in a single pass.

Target to validate: a ranking BAM plus control of about 30 million reads each, with
about 20,000 constituent peaks, finishes in under 15 seconds for indexed runs (and
under 30 seconds with a full `--cn-correct` sweep) on 8 cores.

## 4. Command-line design

```
rose run               -i peaks.bed -r chip.bam [-c input.bam] -a refseq.ucsc -o outdir
rose stitch            -i peaks.bed -a refseq.ucsc           # regions only
rose signal            -i regions.bed -r chip.bam            # signal table only
rose call              -i region_map.tsv                     # cutoff, tables, plot
rose genes             -i enhancers.tsv -a refseq.ucsc       # gene assignment
rose multi             --samples sheet.tsv -a refseq.ucsc -o outdir
rose fetch-annotation  --genome hg38 -o annotations/         # explicit download helper
```

`run` chains the first four stages in memory; `multi` applies them to many samples
over one shared region set (section 6). `fetch-annotation` isolates network access
in its own subcommand so `run` and `multi` never make implicit HTTP requests on
offline HPC compute nodes (`--genome` looks up a local cache path only).

Main options:

- `--assay chip|cutrun`: sets the signal defaults listed in section 5; every
  individual option below still overrides it.
- `--stitch 12500|auto`, `--tss 0` as in ROSE (`12500` default for cross-study
  comparability; `auto` opt-in selects the stitch scale purely from the inter-peak
  gap distribution without touching signal; see section 6.4).
- `--promoter-mode contained|subtract`: drop peaks wholly inside the TSS window
  (ROSE) or cut the window out of peaks (LILY).
- `--max-tss-genes 2`: disaggregate stitched regions spanning strictly more distinct
  gene symbols than this back into their constituent peaks (matching ROSE's code,
  whose comment says 2 while implementing `> 2`).
- `--gene-biotype all|protein_coding`: optional transcript filter when GTF
  annotations are supplied, preventing ~250k GENCODE lncRNA/pseudogene entries from
  over-triggering `--max-tss-genes` un-stitching.
- `--signal bam|bigwig`: bigWig input lets a copy-number-corrected density track
  (for example from HMCan) drive the ranking, which is LILY's main idea.
- `--extend 200`, `--floor 1`, `--min-mapq`, `--exclude-flags`, `--fragment
  auto|cigar|seq|tlen`, `--min-fragment 0`, `--max-fragment 1000`.
- `--norm-denom filtered|index`: divide by post-filter fragment/read count
  (`filtered`, default in clean modes so high-duplicate libraries are not
  under-scaled) or by `.bai`/`.csi` mapped alignment records (`index`, used under
  `--compat rose`).
- `--blacklist regions.bed`: constituents overlapping these regions are removed
  before stitching.
- `--scale-factor X` (and `--control-scale-factor`): replaces reads-per-million
  scaling, for spike-in or carry-over calibration.
- `--subtract-control on|off`, and `--compare-control` to run both ways and report
  the difference in calls.
- `--exclude-chrom chrY` (ROSE drops chrY silently; here it is an explicit default
  that can be turned off).
- `--cutoff exact|rose-compat`, and `--drop-top N` / `--min-se N` as explicit,
  off-by-default versions of LILY's outlier handling.
- `--stability N` and `--seed 42`: resample N times and report cutoff and call
  stability deterministically.
- `--cn-correct` (estimate from the control BAM) or `--cn-segments file.seg`
  (use externally computed segments), with `--cn-bin 100000`.
- `--min-peak-samples 1` (`rose multi`): minimum number of samples supporting a
  merged constituent peak before stitching.
- `--compat rose`: turn on every legacy quirk at once for regression comparison.
- `--threads`, `--genome` (resolves from local annotation directory).

Inputs accepted: BED, GFF, narrowPeak and broadPeak for regions; UCSC refGene tables
and GTF for annotation; BAM, CRAM and bigWig for signal. Annotation files are
supplied locally or fetched beforehand via `rose fetch-annotation`, never committed
to the repository (the original carries 62 MB of them).

## 5. Algorithms

**Coordinate normalisation.** Every parser converts incoming intervals immediately
into a canonical 0-based half-open `Interval` (`[start, end)`): BED, narrowPeak,
broadPeak and UCSC `refGene` (`txStart..txEnd`) map directly, while 1-based closed
GFF and GTF intervals convert via `(start - 1)..end`. Original ROSE's 1-bp GFF
`end - start` width and gap quirks are confined to `--compat rose`.

**Stitching and un-stitching.** Remove blacklisted peaks, apply `--promoter-mode`
against TSS windows, sort remaining constituents by `(chrom, start, end, id)`, and
sweep once, merging adjacent or overlapping intervals when `next.start - cur.end <=
stitch`. Regions with identical coordinates but different IDs remain distinct
constituents prior to stitching. After merging, count **distinct gene symbols**
(`name2` in refGene, `gene_name` in GTF—never transcript IDs) whose TSS (`1 bp`
coordinate) falls inside `[stitched.start, stitched.end)`. Any stitched region
spanning strictly more than `--max-tss-genes` distinct genes is disaggregated back
into its TSS-filtered constituent peaks. Stitched regions are named
`<count>_<leftmost id>_lociStitched`, matching what the ROSE README promises.

**TSS handling.** Build TSS points and promoter windows into sorted per-chromosome
arrays. For duplicated transcript IDs, prefer primary chromosomes over `_alt`,
`_random` and `chrUn` contigs (in hg38 the original takes whichever row appears
first). For GTF files that omit `transcript` feature rows (common in UCSC GTF
exports), reconstruct transcript bounds from the `min(start)..max(end)` of `exon`
features grouped by `transcript_id`.

**Signal from BAM and unit contract.**
- **Units across subcommands:** `rose signal` (`region_map.tsv`) writes each
  region's `width` (`end - start`) and **mean density (`rpm/bp`)**, matching
  `*_ENHANCER_REGION_MAP.txt`. `rose call` computes total integrated ranking signal
  as `(treat_density - ctrl_density).max(0.0) * width` (`rpm · bp`), matching
  `*_AllEnhancers.table.txt`.
- **Single-end reads:** Read span is derived from reference-consuming CIGAR
  operations (`M`, `=`, `X`, `D`, `N`) by default (or `SEQ` length under `--compat
  rose`), then extended in the 5'->3' strand direction by `--extend`.
- **Paired-end fragments:** Proper pairs (`0x2`) are counted **once per fragment**
  without per-read hash tables by emitting `[POS, POS + TLEN)` only from the
  leftmost mate (`TLEN > 0`, breaking `POS == PNEXT` ties via `READ1` `0x40`) and
  filtering `TLEN` to `--min-fragment..=--max-fragment`; mates are never extended
  separately. Indexed fetch windows are widened upstream by `max(extend,
  max_fragment)` so long fragments whose leftmost mate starts upstream of the peak
  are never missed.
- **Coverage accumulation:** With `--floor 0`, total coverage in `[region.start,
  region.end)` is the exact sum of clamped interval overlaps. With `--floor f`
  (`f >= 1`), increments `+1` / `-1` are written into a reusable per-thread
  difference array, prefix-summed once, and bases with `depth > f` contribute their full
  depth, as in `ROSE_bamToGFF.py` (implemented this way in both modes).
- **Normalisation:** Divide by `scale_factor` if supplied; otherwise divide by
  million post-filter units (`--norm-denom filtered`, where 1 paired-end fragment =
  1 unit and 1 single-end read = 1 unit) or by `.bai`/`.csi` `idxstats` mapped
  alignments (`--norm-denom index`, default under `--compat rose`, with a streaming
  fallback when reading CRAM `.crai`).

**Assay presets.** `--assay` only changes defaults:

| | `chip` | `cutrun` |
|---|---|---|
| Fragment | Paired-end: real fragment (`--compat rose`: single-read `--extend 200`). Single-end: read extended by `--extend` | Real fragment from the read pair; single-end input is accepted with a warning |
| Depth floor | 1, as in ROSE | 0, because libraries are shallow |
| Control | Subtracted when `-c` is given | Not subtracted unless `--subtract-control on`; sparse IgG belongs at the peak-calling step |
| Scaling | Million filtered fragments/reads (`idxstats` under `--compat rose`) | Million filtered fragments, or `--scale-factor` for spike-in calibration |
| Blacklist | Optional | Optional, with a warning when absent (artifact regions otherwise rank at the top) |
| Duplicates | Excluded (`0x400` counted under `--compat rose`) | Excluded (warns if `@PG` headers show no duplicate-marking tool) |
| Fragment length | `0..=1000` on proper pairs (unfiltered single-end) | `--min-fragment` / `--max-fragment` applied per pair; set them to match the target (section 13) |

**Signal from bigWig.** Sum `value * overlap_width` over `[region.start,
region.end)` for total signal, and divide by `region.width` for mean density, so
variable-width bigWig intervals are weighted accurately (fixing LILY's unweighted
bin sum).

**Chromosome names.** Reconcile `chr1` and `1` once against the BAM/bigWig header,
and abort with a clear diagnostic if fewer than a configurable share (default 50%)
of input regions match any contig. This eliminates the original's silent all-zero
output on Ensembl-named BAMs.

**Cutoff.** Clamp negative background-subtracted signals to `0` and sort ascending
as `y[1..=n]`. For `n >= 2` and `y[n] > y[1]`, define chord slope `s = (y[n] - y[1])
/ (n - 1)` (`(y[n] - y[1]) / n` under `--compat rose`) and find `i* = argmin(y[i] -
s * i)`, breaking ties in favour of the **largest** `i` so zero-signal plateaus break
at lift-off. The cutoff is `y[i*]`, and regions with **`signal > cutoff`** (strict
inequality, matching `ROSE_callSuper.R` so the tangent point itself is typical) are
classified as super-enhancers. `rose-compat` also ports R's `stats::optimize`
(Brent's method on the step count function) for exact legacy regression checks.

**Gene assignment.** Report overlapping, proximal (default 50 kb) and closest genes
per region via binary search over sorted transcript intervals and TSS arrays. Table
headers are written from typed schema definitions, and no rows are skipped based on
the first character of an identifier.

## 6. Added features

### 6.1 Multi-sample and replicate mode

Original ROSE handles one sample at a time, so comparing samples means reconciling
different region sets by hand.

- **Input.** A sample sheet with columns `sample`, `group`, `peaks`, `bam`, and
  optional `control` and `scale_factor`.
- **Shared regions without megablock chaining.** Constituent peaks from all samples
  are pooled, filtered against `--blacklist`, and **collapsed (`bedtools merge`
  semantics) into non-overlapping consensus constituent elements** supported by at
  least `--min-peak-samples` samples (default `1`). TSS handling and 12.5 kb
  stitching then run once on those consensus constituents, preventing sample count
  from inflating `<count>_..._lociStitched` names or chaining noise peaks across
  entire TADs.
- **Signal matrix.** Regions by samples, computed in one parallel pass per BAM.
- **Per-sample calls.** The cutoff is computed separately for each sample on the
  shared regions.
- **Consensus.** Within a group, a region is a consensus super-enhancer when it is
  called in at least `--min-replicates` samples (default: all).
- **Comparison.** Between two groups, regions are labelled gained, lost or shared,
  with the log2 ratio of mean normalised signal (with a small pseudocount reported
  in `run.json`). No p-values are produced.
- **Normalisation between samples.** Per-sample `scale_factor` if given; otherwise a
  DESeq2-style median-of-ratios size factor computed over the subset of regions with
  strictly positive signal in at least 80% of samples (avoiding geometric-mean
  collapse from lineage-specific zero-signal enhancers; falls back to library-size
  RPM with a warning if too few shared positive regions exist).

### 6.2 Cutoff stability report

The tangent cutoff gives a hard yes or no, with nothing to show how fragile a call
is. `--stability N` adds:

- **Region bootstrap.** Resample regions with replacement `N` times using a PRNG
  seeded by `--seed` and recompute the cutoff. Needs no extra BAM reading.
- **Read downsampling.** Each fragment/read is assigned to nested depth tiers (for
  example 50%, 75%, 90%, 100%) by hashing `(seed, QNAME)` with `xxhash-rust`
  (`xxh3_64`), never Rust's randomized `DefaultHasher`. When `--floor 0` is active,
  per-tier overlap sums accumulate directly; when `--floor > 0` (`chip` default),
  depth thresholding is nonlinear, so the region worker maintains one difference
  array per downsampling tier during the single BAM pass.
- **Outputs.** A confidence interval for the cutoff; for each region, its bootstrap
  call frequency and a `borderline` flag (e.g. frequency in `[0.1, 0.9]`); and a
  saturation TSV/plot of super-enhancer count across downsampled read depths.

### 6.3 Built-in copy-number correction

Amplified regions gain signal simply because there is more DNA, which is the
problem LILY solves through HMCan. Two routes:

- **From the control BAM (`--cn-correct`).** Count post-filter control reads/fragments
  in fixed bins (`--cn-bin`, default 100 kb as in HMCan) during the **same single
  parallel pass** over the control BAM used for region background counting. Exclude
  blacklisted bins and unmappable/near-zero bins (`< 10%` of the positive-bin median)
  before computing the autosomal median (handling sex chromosomes separately), smooth
  with a running median across neighbouring mappable bins, and clamp extreme ratios
  to a configurable `[min_cn_ratio, max_cn_ratio]` range (leaving unmappable bins at
  `1.0` with an `uncorrected` flag so centromeric gaps never amplify noise).
- **From external segments (`--cn-segments`).** Use a `.seg` or BED file of
  copy-number segments from whole-genome sequencing or a dedicated caller. Needs no
  control BAM scan.
- **Order of operations.** Copy-number amplification inflates both target ChIP and
  control Input equally, so background subtraction happens **before** dividing by the
  local copy-number ratio: `corrected_signal = (treat_signal - ctrl_signal).max(0.0)
  / cn_ratio`.
- **Safeguards.** Refuse `--cn-correct` when the control BAM is too shallow for
  stable 100 kb bin counts (for example `< 50` median reads per bin—always the case
  for CUT&RUN IgG, where `--cn-segments` is the only route). Both raw signal,
  corrected signal and local `cn_ratio` are written to the output table.

### 6.4 Stitch-distance sensitivity and `--stitch auto`

ROSE's fixed 12.5 kb stitch window was chosen empirically on mouse ESC Mediator
ChIP-seq: broad peak callers (such as MACS3 `--broad` in CUT&RUN) are already
partly merged so 12.5 kb over-stitches, whereas compact genomes or sparse/dense
factors cluster at different scales. Conversely, tuning the stitch distance per
sample to "improve" the ranking curve introduces circularity, replicate jitter, and
incomparable boundaries across samples. `rose-rs` addresses both sides cleanly:

- **Zero-I/O stitch scanning (piecewise peak + gap decomposition).** After blacklist
  and TSS filtering, collapse overlapping/abutting (`d = 0`) constituents into
  strictly disjoint base peak blocks $B_0, \dots, B_{K-1}$ per chromosome and emit
  every inter-peak gap interval $G_i = [B_i.\text{end}, B_{i+1}.\text{start})$ up to
  `max(active_stitch, 25 kb)` (plus individual constituent intervals for any rare
  overlapping inputs needed under `--max-tss-genes` disaggregation). Because base-pair
  coverage (`floor = 0`, `floor > 0` difference-array sums, and width-weighted
  bigWig integrals) is strictly additive across disjoint half-open intervals,
  quantifying these elementary peak and gap pieces **once** during the single BAM or
  bigWig pass makes evaluating *any* stitch distance exact in milliseconds with zero
  further BAM reading.
- **Sensitivity report on every run (`*_Stitch_sensitivity.tsv`).** On every `run`
  and `multi` invocation (unless `--no-stitch-report` is passed), evaluate a grid of
  stitch distances from `0` to `25 kb` (`0, 1, 2.5, 5, 7.5, 10, 12.5, 15, 20, 25 kb`
  plus the active distance) and report: super-enhancer count, tangent cutoff, median
  and 90th-percentile super-enhancer width, the share of each region covered by
  constituent peaks (`peak_coverage_frac`, diagnosing broad-peak over-stitching), and
  base-pair Jaccard + recall against the canonical `12.5 kb` result.
- **`--stitch auto` (opt-in, signal-blind gap excess).** Estimate the background
  Poisson inter-peak gap rate $\hat{\lambda}_{\text{bg}}$ from the log-survival slope
  of unclustered euchromatic gaps ($q_{0.65}\dots q_{0.95}$ of gaps $\le 1\text{ Mb}$,
  excluding centromeric deserts and short-gap clustering) and choose $d^*$ maximizing
  the cumulative short-gap excess $D(d) = \hat{F}_{\text{obs}}(d) - (1 -
  e^{-\hat{\lambda}_{\text{bg}} d})$—the exact scale where clustered short gaps end
  and gap density returns to random expectation. Using the ECDF rather than a binned
  histogram avoids smoothing jitter across replicates; never touching signal avoids
  circular tuning.
- **Guard rails.** Clamp $d^*$ to `[2 kb, 25 kb]`; fall back to `12,500 bp` (logged in
  `stderr` and `run.json`) when fewer than `100` euchromatic gaps exist or when peak
  short-gap excess is `< 5%` (`max D(d) < 0.05`); and in `rose multi`, compute one
  shared $d^*$ from the pooled consensus peaks so all samples share identical region
  boundaries.

## 7. Outputs

- The ROSE tables under their original names and column orders
  (`*_AllEnhancers.table.txt`, `*_SuperEnhancers.table.txt`,
  `*_ENHANCER_REGION_MAP.txt`, `*_Enhancers_withSuper.bed`, gene tables).
- A class column (`super_enhancer`, `enhancer`, `promoter`), as LILY reports.
- Spec-compliant BED files with 0-based half-open coordinates and column-5 integer
  scores scaled/clamped to `0..=1000` (raw unscaled rank/signal preserved under
  `--compat rose`).
- The hockey-stick plot rendered natively (`plotters` with embedded font, SVG and
  PNG), plus the underlying curve points as TSV.
- `run.json` recording tool version, resolved arguments, input checksums, raw and
  post-filter read/fragment counts, normalisation factors and the cutoff, for
  reproducibility.
- With the added features: the signal matrix, consensus and group-comparison tables
  (6.1); bootstrap stability columns and depth-saturation table (6.2); and raw
  signal, corrected signal and copy-number ratio per region (6.3).

All files are written to a sibling temporary file in the target directory and
atomically renamed on success, so a rerun can never read a stale or half-written
file. Any failure returns a non-zero exit code.

## 8. Defects from the reviews this design removes

| Source | Defect | Addressed by |
|---|---|---|
| ROSE | All-zero signal for BAMs without `chr` | Contig-name reconciliation with a hard minimum-match check |
| ROSE | Stale or partial files used on rerun; 12-hour wait on a crashed job | Direct thread joins; atomic temp-file writes; `Result` error propagation |
| ROSE | Gateway BED rows duplicated on rerun | Files are always atomically rewritten |
| ROSE | Gene mapper drops IDs starting with `R`; broken header; duplicated columns | Typed table reader and writer |
| ROSE | Standard and turbo paths disagree (depth floor) | One unified signal kernel; `--floor` is an explicit documented option |
| ROSE | Same-coordinate regions merged silently | Constituent IDs are preserved through pre-stitch sorting |
| ROSE | Unquoted shell commands; must run from the repo directory | Zero shell subprocesses; no working-directory assumptions |
| ROSE | Python 2 only | Self-contained compiled Rust binary |
| LILY | Hardcoded floor of 600 super-enhancers | `--min-se`, off by default |
| LILY | Signal not width-weighted | Width-weighted bigWig integration |
| LILY | Threshold assumes pre-sorted peaks | Inputs are sorted internally |
| LILY | GFF annotation treats every feature as a TSS | Transcript/gene feature filtering with `exon`-grouping fallback |
| LILY | Errors exit with status 0 | Non-zero exit on every failure |
| LILY | Whole bigWig loaded into memory twice | Indexed per-chromosome/per-region bigWig queries |

## 9. Project layout

```
rose-rs/
  Cargo.toml          workspace dependencies and [workspace.lints]
  crates/rose-core/   pure algorithms (#![forbid(unsafe_code)], zero I/O or C deps):
                        intervals, stitch, tss, signal, cutoff, stability, cn, genes
  crates/rose-io/     file formats & atomic writers:
                        bam/cram (noodles + optional htslib), bigwig, bed/gff/gtf, tsv
  crates/rose-cli/    clap subcommands, rayon orchestration, plotters SVG/PNG
  tests/golden/       small inputs with expected outputs from original ROSE
  benches/            criterion micro-benchmarks
```

Keeping `rose-core` completely free of filesystem, BAM and bigWig dependencies makes
unit and property testing instant on every OS and lets Python (`PyO3`) or R
bindings wrap the algorithms directly without pulling in format libraries.

## 10. Validation

1. **Golden comparison.** Run original ROSE in a pinned Python 2.7 / R / samtools
   container on the Young lab example data and two or three public H3K27ac
   datasets. `--compat rose` must match region sets exactly, densities and signals
   within rounding (the original rounds to 4 decimals), and the super-enhancer set
   and gene assignments exactly.
2. **Coordinate and stitching property tests.** `proptest` suites checking 0-based
   BED vs 1-based GFF equivalence, stitching + `--max-tss-genes` disaggregation, and
   interval/difference-array coverage against a brute-force per-base oracle.
3. **Signal and paired-end cross-check.** Single-end coverage compared with
   `samtools depth` / `bedtools coverage`; paired-end fragment coverage (including
   long fragments starting upstream of `[region.start - extend, ...)`) compared with
   `bedtools genomecov -pc`; CRAM `.crai` streaming fallback verified against BAM.
4. **CUT&RUN.** On a public paired-end H3K27ac CUT&RUN dataset: fragment-based
   coverage checked against `bedtools genomecov -pc`, and super-enhancer calls
   compared between replicates, with `--compare-control` on/off, and with/without the
   blacklist.
5. **Multi-sample mode.** Running one sample through `multi` (`--min-peak-samples 1`)
   must reproduce `run` exactly; pooled noisy multi-sample peaks must not inflate
   constituent counts or chain across TADs; zero-inflated enhancer matrices must
   yield finite median-of-ratios size factors; a sample compared with itself must
   give zero gained or lost regions.
6. **Stability report.** Fixed `--seed` values produce byte-identical outputs across
   platforms; downsampled signals (both `--floor 0` and `--floor 1`) match BAMs
   subsampled with `samtools view -s` within binomial sampling noise.
7. **Copy-number correction.** A synthetic test (reads in one locus amplified `5x`
   in **both** target and control BAMs, plus simulated centromeric zero-count bins)
   must recover the uncorrected ranking without amplifying gap noise. On public
   neuroblastoma data (`MYCN`-amplified), compare `--cn-correct`, `--cn-segments`,
   and LILY run on HMCan bigWigs.
8. **Edge cases.** Zero or one super-enhancer, all-zero signal plateaus, fewer
   regions than `--min-se`, empty chromosomes, missing `.bai`/`.crai` index,
   mismatched contig names, duplicate peak IDs, UCSC GTFs without `transcript` rows.
9. **Benchmarks.** `hyperfine` for end-to-end wall time and peak RSS against original
   ROSE; `criterion` for stitching, difference-array coverage, and `noodles`
   (`zlib-rs`) vs optional `rust-htslib`.

## 11. Milestones

| | Deliverable | Done when |
|---|---|---|
| M0 | Reference container and golden outputs | Original ROSE runs reproducibly; expected files are stored |
| M1 | `rose-core` intervals/TSS/stitch + `rose-io` BED/GFF/GTF/refGene parsers + `rose stitch` | Stitched and un-stitched regions match golden output across BED/GFF and TSS modes |
| M2 | Unified `rose signal` BAM/CRAM kernel (`noodles`, single-end + paired-end fragments, flags/dedup/MAPQ/length/blacklist filters, `filtered` & `index` denominators) | Signal table matches golden output in `--compat rose` and matches `bedtools` on paired-end BAMs |
| M3 | `rose call` (exact + `rose-compat` Brent cutoff) and static-safe `plotters` SVG/PNG | Super-enhancer set and tables match golden output; plot renders in static `musl` and Windows builds |
| M4 | `rose genes`, `rose run`, `rose fetch-annotation`, and `run.json` | Full single-sample pipeline matches original ROSE end to end |
| M5 | `rayon` parallelism, indexed vs single-pass read-strategy auto-switching, optional `htslib` benchmark | Speed and RSS targets in section 3 met and documented; `htslib` feature kept or dropped based on data |
| M6 | `bigtools` bigWig signal, promoter `subtract` mode, outlier options (`--drop-top`, `--min-se`) | LILY-style run on HMCan bigWig output passes validation |
| M7 | Multi-sample and replicate mode (`rose multi` with consensus peak merging and robust size factors) | Multi-sample validation in section 10 passes |
| M8 | Cutoff stability report (`--stability`, seeded `xxh3` multi-tier downsampling + region bootstrap) | Stability validation passes; overhead stays small next to the base run |
| M9 | Copy-number correction (`--cn-correct` fused control pass + `--cn-segments`) | Synthetic and neuroblastoma copy-number validations pass |
| M10 | `--assay cutrun` preset tuning, `@PG` dedup-header check, `--compare-control` report, end-to-end CUT&RUN validation | The workflow in section 13 runs end to end and passes CUT&RUN replicate/blacklist validation |

M1 to M4 lock correctness and build the unified single-end/paired-end signal model up
front, so M5, M8 and M10 never have to rewrite the BAM inner loop.

## 12. Open decisions

- **Licence.** ROSE is MIT; LILY's scripts are GPL v2 or later. Reimplementing
  LILY's ideas from their published method description is compatible with an MIT
  licence; porting its code is not. Recommendation: MIT, clean-room implementation.
- **Default behaviour.** Recommendation: clean semantics by default (`noodles`
  CIGAR/fragment spans, `--norm-denom filtered`, exact cutoff with rightmost
  tie-break, leftmost `<count>_<id>_lociStitched` naming, `0..=1000` BED scores) and
  all legacy quirks behind `--compat rose`.
- **Read filtering and scaling defaults.** Original ROSE counts duplicates and
  secondary alignments and scales by `flagstat` mapped alignments. Recommendation:
  exclude `0x4` (unmapped), `0x100` (secondary), `0x200` (QC-fail), `0x400`
  (duplicate) and `0x800` (supplementary) and scale by post-filter count by default;
  restore ROSE's behaviour under `--compat rose`.
- **Paired-end ChIP-seq.** Proposed default is the real paired-end fragment span
  `[POS, POS + TLEN)`, with single-read `--extend 200` retained under `--compat
  rose`.
- **GTF biotype default.** Proposed default is `--gene-biotype all` for backward
  compatibility with RefSeq GTFs, while documenting `--gene-biotype protein_coding`
  (or RefSeq `NM_`/`NR_` equivalents) when full GENCODE GTFs are used.
- **Copy-number smoothing.** Proposed first version is a running median over
  mappable 100 kb bins. Circular binary segmentation or GC-content regression (as in
  HMCan) can be evaluated after M9 comparison with LILY.
- **Differential calls.** The comparison in 6.1 is descriptive (labels and log2 fold
  change); statistical testing stays in downstream tools (DESeq2 / edgeR).

## 13. CUT&RUN workflow

CUT&RUN preset validation is the final milestone (M10), building on the paired-end
fragment and filtering engine implemented in M2. This section records the workflow
it serves and the adjustments it makes.

### Typical CUT&RUN to super-enhancer workflow today

1. **Call constituent enhancers.** Run H3K27ac (or Mediator / BRD4) CUT&RUN data
   through nf-core/cutandrun using MACS3 in broad peak mode, producing `.broadPeak`
   / `.bed` peak files and coordinate-sorted `.bam` files.
2. **Convert BED to GFF.** Format the peaks as the 9-column `.gff` ROSE requires:
   column 1 is the chromosome (`chr#`), columns 4 and 5 are 1-based start and end,
   and columns 2 and 9 hold a unique peak ID.
3. **Run ROSE.** Run `ROSE_main.py` with the genome build (`hg38`, `mm10`, ...), the
   `.gff` peak file, the target `.bam`, and optionally an IgG control `.bam`.

### The same workflow with `rose-rs`

| Step today | With `rose-rs` |
|---|---|
| 1. nf-core/cutandrun with MACS3 broad peaks | Unchanged; peak calling stays outside the tool |
| 2. Convert BED to GFF | Not needed. BED and MACS3 `broadPeak` files are read directly as 0-based half-open intervals; missing peak IDs are generated deterministically and duplicate IDs are rejected with a clear message |
| 3. `ROSE_main.py -g ... -i peaks.gff -r target.bam [-c igg.bam]` | `rose run --assay cutrun -i peaks.broadPeak -r target.bam -a refseq.ucsc --blacklist blacklist.bed -o outdir` |

### Key adjustments for CUT&RUN data

- **IgG control handling.** Unlike ChIP-seq input controls, CUT&RUN IgG controls
  have very sparse reads. Passing one through ROSE's `-c` can over-penalise
  low-background regions and add noise, so many researchers run ROSE without `-c`
  or check results both ways. (ROSE subtracts the control and clamps negatives to
  zero; it does not divide by it, so the risk is over-subtraction, not division by
  zero.) In `rose-rs` the `cutrun` preset does not subtract the control unless
  `--subtract-control on` is given, and `--compare-control` runs both ways in one
  pass and reports how many super-enhancer calls change.
- **Deduplication.** Signal should be computed from deduplicated fragments. The
  `cutrun` preset skips reads flagged as duplicates (`0x400`), scales by post-filter
  fragment count (`--norm-denom filtered`), and warns if the BAM `@PG` headers show
  no duplicate-marking step, since unmarked duplicates cannot be detected from flags.
- **Fragment size filtering.** Fragments are filtered to lengths appropriate for the
  target during counting (`--min-fragment` and `--max-fragment`), so the BAM does
  not need to be pre-filtered on disk: short sub-nucleosomal fragments (e.g.
  `<= 120 bp`) for transcription factors, or nucleosomal fragments (e.g.
  `150..=500 bp`) for histone marks such as H3K27ac. Exact preset suggestions will be
  validated against public CUT&RUN datasets in M10.
