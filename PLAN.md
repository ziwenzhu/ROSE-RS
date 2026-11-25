# Plan: a fast rewrite of ROSE

A plan to reimplement ROSE (Rank Ordering of Super-Enhancers) as a single compiled
tool, fixing the defects found in the review of
[young_computation/rose](https://bitbucket.org/young_computation/rose) and adopting
the useful ideas from [BoevaLab/LILY](https://github.com/BoevaLab/LILY).

Working name: `rose-rs`. Nothing here has been built or measured yet; performance
figures below are targets to validate, not results.

## 1. Goals

- **Same answers as ROSE by default.** Given the same inputs, a compatibility mode
  reproduces the original's tables closely enough to be a drop-in replacement.
- **Much faster.** Minutes-to-hours becomes seconds-to-minutes on a laptop.
- **Fails loudly.** No silent zeros, stale outputs, or success exit codes on error.
- **One static binary.** No Python 2, no R, no `samtools` on `PATH`, no need to run
  from the install directory.

Supported assays: **ChIP-seq and CUT&RUN**. The algorithm is the same for both; they
differ in how signal is counted and normalised (see "Assay presets" in section 5).

Beyond a faster ROSE, four additions (section 6): a multi-sample and replicate
mode, a cutoff stability report, built-in copy-number correction, and a stitch
distance sensitivity report with optional automatic selection.

Non-goals: CUT&Tag (its Tn5 accessibility bias needs separate treatment and is out
of scope), peak calling, and statistical testing of differential signal (the signal
matrix is exported for tools such as DESeq2 or edgeR).

## 2. Language and toolchain: Rust

| Need | Choice |
|---|---|
| BAM/CRAM access | `rust-htslib` (htslib bindings); `noodles` (pure Rust) as a possible replacement |
| bigWig access | `bigtools` |
| Parallelism | `rayon`, one reader per thread |
| Command line | `clap` |
| Plotting | `plotters` |
| Testing and benchmarks | `cargo test`, `proptest`, `criterion` |
| Distribution | Static musl binary; Bioconda recipe |

The hot loop is "decompress BAM blocks, decode records, add to counters". htslib
with libdeflate is the fastest well-tested implementation of the first two steps, so
the plan starts with `rust-htslib` and keeps the rest of the code in safe Rust.

To confirm at kickoff: `rust-htslib`'s libdeflate and static-build options, and
whether `noodles` is fast enough to drop the C dependency entirely.

## 3. Where the original loses time, and the fix

| Original behaviour | Cost | Replacement |
|---|---|---|
| One `samtools view` subprocess per region, output parsed as text | Process spawn and SAM text formatting for every region | In-process indexed reads of binary records |
| Per-base coverage held in Python dicts | Millions of hash inserts per large region | Overlap arithmetic, or a reusable difference array when a depth floor is needed |
| Each BAM mapped twice (stitched regions and constituents) | 2x the work; the constituent table is never used for ranking | One pass computes both; the constituent table becomes opt-in |
| `samtools flagstat` for the mapped-read count | A full scan of every BAM | Mapped counts from the index (as `idxstats` does), no scan |
| Background jobs polled every 5 minutes | Any run longer than 10 seconds takes at least 5 minutes | Threads joined directly |
| Everything serial | One core | `rayon` across chromosomes and BAMs, one reader per thread |
| Gene mapper uses `list.count` per annotation row | Quadratic in transcript count | Sorted TSS arrays with binary search |
| Stitching by repeated dict lookups | Hard to reason about, order-dependent names | Sort, then a single sweep |
| `optimize()` on a step function for the cutoff | Can stop at a local minimum | Exact linear scan (see section 5) |

Read access strategy: merge the fetch windows per chromosome and use indexed
queries when regions cover a small share of the genome; fall back to one sequential
pass per chromosome when they cover most of it.

Target to validate: a ranking BAM plus control of about 30 million reads each, with
about 20,000 constituent peaks, finishes in under a minute on 8 cores.

## 4. Command-line design

```
rose run     -i peaks.bed -r chip.bam [-c input.bam] -a refseq.ucsc -o outdir
rose stitch  -i peaks.bed -a refseq.ucsc           # regions only
rose signal  -i regions.bed -r chip.bam            # signal table only
rose call    -i region_map.tsv                     # cutoff, tables, plot
rose genes   -i enhancers.tsv -a refseq.ucsc       # gene assignment
rose multi   --samples sheet.tsv -a refseq.ucsc -o outdir   # several samples
```

`run` is the first four chained in memory; `multi` applies them to many samples
over one shared region set (section 6). The subcommands exist so each stage can
be tested, cached and rerun alone.

Main options:

- `--assay chip|cutrun`: sets the signal defaults listed in section 5; every
  individual option below still overrides it.
- `--stitch 12500|auto`, `--tss 0`. The default stays 12.5 kb as in ROSE; `auto`
  chooses the distance from the peak gaps (section 6.4).
- `--stitch-scan 0:25000:500` (on by default, `--no-stitch-scan` to skip): report
  how the result changes across stitch distances.
- `--promoter-mode contained|subtract`: drop peaks inside the TSS window (ROSE) or
  cut the window out of peaks (LILY).
- `--max-tss-genes 2`: un-stitch regions spanning more genes than this (the ROSE
  code uses "more than 2" although its comment says 2).
- `--signal bam|bigwig`: bigWig input lets a copy-number-corrected density track
  (for example from HMCan) drive the ranking, which is LILY's main idea.
- `--extend 200`, `--floor 1`, `--min-mapq`, `--exclude-flags`, `--fragment
  auto|cigar|seq|tlen`, `--min-fragment 0`, `--max-fragment 1000`.
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
- `--stability N`: resample N times and report cutoff and call stability.
- `--cn-correct` (estimate from the control BAM) or `--cn-segments file.seg`
  (use externally computed segments), with `--cn-bin 100000`.
- `--compat rose`: turn on every legacy quirk at once for regression comparison.
- `--threads`, `--genome` (optional shortcut to a downloaded annotation).

Inputs accepted: BED, GFF, narrowPeak and broadPeak for regions; UCSC refGene tables and GTF for
annotation; BAM and CRAM for reads. Annotation files are downloaded or supplied by
the user, not committed to the repository (the original carries 62 MB of them).

## 5. Algorithms

**Stitching.** Sort constituents by chromosome and start, sweep once, merge when the
gap is at most `--stitch`. Name each region `<count>_<leftmost id>_lociStitched`,
which is what the ROSE README promises and the code does not deliver. Regions with
identical coordinates but different IDs are kept as separate constituents.

**TSS handling.** Build TSS windows from the annotation into a sorted array per
chromosome. For duplicated transcript IDs, prefer primary chromosomes over alt and
random contigs (in hg38 the original takes whichever row comes first).

**Signal from BAM.** For each region, visit reads overlapping the region widened by
the extension. With no floor, signal is the sum of overlap lengths between extended
reads and the region. With a floor, fill a per-thread difference array, prefix-sum
it, and add positions whose depth exceeds the floor. Normalise by million mapped
reads. Fragment length comes from the CIGAR by default, not the sequence length.

Signal is stored per piece, not per stitched region: one value for every constituent
peak and one for every gap between neighbouring peaks up to the largest distance of
interest. Signal is a sum over positions (with or without a floor), so a stitched
region's signal is the sum of its pieces. Any stitch distance can then be evaluated
from prefix sums without reading the BAM again, and the constituent table comes for
free. Under `--compat rose` regions are scored whole, to reproduce the original's
endpoint handling exactly.

**Assay presets.** `--assay` only changes defaults:

| | `chip` | `cutrun` |
|---|---|---|
| Fragment | Paired-end: real fragment. Single-end: read extended by `--extend` | Real fragment from the read pair; single-end input is accepted with a warning |
| Depth floor | 1, as in ROSE | 0, because libraries are shallow |
| Control | Subtracted when `-c` is given | Not subtracted unless `--subtract-control on`; sparse IgG belongs at the peak-calling step |
| Scaling | Reads per million | Reads per million, or `--scale-factor` for spike-in calibration |
| Blacklist | Optional | Optional, with a warning when absent (artifact regions otherwise rank at the top) |
| Duplicates | Excluded (counted under `--compat rose`) | Excluded |
| Fragment length | No filter | `--min-fragment` / `--max-fragment` applied per pair; set them to match the target (section 13) |

With `--fragment auto`, pairing is detected from the BAM flags. For paired-end data
each proper pair is counted once as one fragment spanning both mates, and pairs
outside `--min-fragment`..`--max-fragment` are skipped; mates are never extended
separately.

**Signal from bigWig.** Sum value times overlap width, so the result is correct for
variable-width intervals (LILY sums values without weighting by width).

**Chromosome names.** Reconcile `chr1` and `1` once, against the BAM header, and
abort with a clear message if fewer than a set share of regions match any contig.
This removes the original's all-zero output for BAMs without `chr`.

**Cutoff.** With scores sorted ascending as `y[1..n]` and slope
`s = (max - min) / n`, the tangent point is exactly `argmin(y[i] - s*i)`: a line of
slope `s` through that point has the fewest points beneath it. This is a single
linear pass and has no local-minimum problem. `rose-compat` instead ports R's
`optimize` (Brent's method) to reproduce the original cutoff exactly.

**Gene assignment.** Overlapping, proximal (default 50 kb) and closest gene per
region, via binary search over sorted transcript and TSS arrays. Headers are written
from real column names, and no rows are skipped based on the first letter of an ID.

## 6. Added features

### 6.1 Multi-sample and replicate mode

Original ROSE handles one sample at a time, so comparing samples means reconciling
different region sets by hand.

- **Input.** A sample sheet with columns `sample`, `group`, `peaks`, `bam`, and
  optional `control` and `scale_factor`.
- **Shared regions.** Constituent peaks from all samples are pooled and stitched
  once, so every sample is scored on identical regions.
- **Signal matrix.** Regions by samples, computed in one parallel pass per BAM.
- **Per-sample calls.** The cutoff is computed separately for each sample on the
  shared regions.
- **Consensus.** Within a group, a region is a consensus super-enhancer when it is
  called in at least `--min-replicates` samples (default: all).
- **Comparison.** Between two groups, regions are labelled gained, lost or shared,
  with the log2 ratio of mean normalised signal. No p-values are produced.
- **Normalisation between samples.** Per-sample `scale_factor` if given; otherwise
  a median ratio over regions present in most samples (the idea behind LILY's
  renormalisation script, computed on width-weighted signal).

### 6.2 Cutoff stability report

The tangent cutoff gives a hard yes or no, with nothing to show how fragile a call
is. `--stability N` adds:

- **Region bootstrap.** Resample regions with replacement N times and recompute the
  cutoff. Needs no extra BAM reading.
- **Read downsampling.** Each fragment is assigned to nested subsets (for example
  50%, 75%, 90% of reads) by hashing its name during the normal counting pass, so
  signals at several depths come from the same single pass.
- **Outputs.** A confidence interval for the cutoff; for each region, the share of
  resamples in which it is called and a `borderline` flag; a saturation curve of
  super-enhancer count against sequencing depth.

### 6.3 Built-in copy-number correction

Amplified regions gain signal simply because there is more DNA, which is the
problem LILY solves through HMCan. Two routes:

- **From the control BAM (`--cn-correct`).** Count control fragments in fixed bins
  (`--cn-bin`, default 100 kb as in HMCan), divide by the genome-wide median to get
  a copy-number ratio, smooth it, and divide each region's signal by the local
  ratio (clamped to a sane range). This needs one full pass over the control BAM,
  so it adds to the run time in section 3.
- **From external segments (`--cn-segments`).** Use a `.seg` or BED file of
  copy-number segments from whole-genome sequencing or a dedicated caller. No extra
  BAM pass.

Safeguards: refuse `--cn-correct` when the control is too shallow to give stable
bin counts (always the case for CUT&RUN IgG, where `--cn-segments` is the only
route); handle sex chromosomes separately from the autosomal median; exclude
blacklisted bins. Both raw and corrected signal are written, with the ratio used.

### 6.4 Stitch distance: sensitivity report and automatic selection

The 12.5 kb distance is an empirical choice from one dataset. It is sensitive to
the peak caller (broad peaks are already partly merged), and to how densely a factor
binds. Choosing it per sample from the data, however, costs comparability between
samples and with published work, and can tune the answer if the choice looks at the
ranking. The design therefore separates showing the dependence from acting on it.

- **Default unchanged.** 12.5 kb unless the user asks otherwise.
- **Sensitivity report (on by default).** Using the per-piece signals from section
  5, re-stitch and re-call at each distance in `--stitch-scan`. For each distance
  report: number of stitched regions, number of super-enhancers, median region
  size, share of each region covered by peaks, and overlap of the super-enhancer
  set with the result at the distance actually used. TSS un-stitching is reapplied
  at each distance. A warning is raised when the calls change sharply near the
  distance used.
- **`--stitch auto` (opt-in).** Chosen from the gaps between neighbouring peaks
  only, never from signal or from the ranking, so it cannot tune the result.
  Randomly placed peaks at the observed density give an expected gap distribution;
  clustered peaks show an excess of short gaps. The distance is the gap length at
  which that excess ends (observed gap frequency falls to the expected level),
  pooled over chromosomes.
- **Guard rails.** The automatic value is clamped to a range (proposed 2 to 25 kb);
  if there is no clear excess of short gaps the tool falls back to 12.5 kb and says
  so; the choice is deterministic and written to `run.json` with the gap curve.
- **Multi-sample mode.** One distance for all samples, chosen from the pooled peaks,
  so the shared region set stays comparable.

Not planned: ROSE2-style optimisation on region coverage, any criterion based on
super-enhancer count or curve shape, and stitching by signal in the gap instead of
by distance (a different algorithm, possible later experiment).

## 7. Outputs

- The ROSE tables under their original names and column orders
  (`*_AllEnhancers.table.txt`, `*_SuperEnhancers.table.txt`,
  `*_ENHANCER_REGION_MAP.txt`, `*_Enhancers_withSuper.bed`, gene tables).
- A class column (super-enhancer, enhancer, promoter), as LILY reports.
- BED files with 0-based starts and integer scores.
- The hockey-stick plot rendered natively (`plotters`, PNG and SVG), plus the
  underlying points as TSV.
- `run.json` recording version, arguments, input checksums, read counts and the
  cutoff, for reproducibility.

- With the added features: the signal matrix and consensus and comparison tables
  (6.1), stability columns and the saturation curve (6.2), raw signal,
  corrected signal and copy-number ratio per region (6.3), and the stitch-distance
  scan table and plot, with the gap curve when `--stitch auto` is used (6.4).

All files are written to a temporary name and renamed on success, so a rerun can
never read a stale or half-written file. Any failure returns a non-zero exit code.

## 8. Defects from the reviews this design removes

| Source | Defect | Addressed by |
|---|---|---|
| ROSE | All-zero signal for BAMs without `chr` | Name reconciliation with a hard check |
| ROSE | Stale or partial files used on rerun; 12-hour wait on a crashed job | No polling; atomic writes; error propagation |
| ROSE | Gateway BED rows duplicated on rerun | Files are always rewritten |
| ROSE | Gene mapper drops IDs starting with `R`; broken header; duplicated columns | Typed table reader and writer |
| ROSE | Standard and turbo paths disagree (depth floor) | One code path; floor is a documented option |
| ROSE | Same-coordinate regions merged silently | IDs are part of constituent identity |
| ROSE | Unquoted shell commands; must run from the repo directory | No shell calls; no dependence on the working directory |
| ROSE | Python 2 only | Compiled binary |
| LILY | Hardcoded floor of 600 super-enhancers | `--min-se`, off by default |
| LILY | Signal not width-weighted | Weighted bigWig sums |
| LILY | Threshold assumes pre-sorted peaks | Inputs are sorted internally |
| LILY | GFF annotation treats every feature as a TSS | Only transcript or gene features are used |
| LILY | Errors exit with status 0 | Non-zero exit on every failure |
| LILY | Whole bigWig loaded into memory twice | Region queries only |

## 9. Project layout

```
rose-rs/
  crates/rose-core/   intervals, stitch, annotation, signal, cutoff, genes, io
  crates/rose-cli/    argument parsing, subcommands, plotting
  tests/golden/       small inputs with expected outputs from original ROSE
  benches/            criterion micro-benchmarks
```

`rose-core` is a library with no I/O assumptions, so Python bindings (PyO3) or an R
wrapper can be added later without touching the algorithms.

## 10. Validation

1. **Golden comparison.** Run original ROSE in a pinned Python 2.7 / R / samtools
   container on the Young lab example data and two or three public H3K27ac
   datasets. `--compat rose` must match region sets exactly, signals within rounding
   (the original rounds to 4 decimals), and the super-enhancer set exactly.
2. **Property tests.** Stitching and overlap code checked against a brute-force
   implementation on random intervals.
3. **Signal cross-check.** Per-region coverage compared with `samtools depth` /
   `bedtools coverage` on a small BAM.
4. **CUT&RUN.** On a public paired-end H3K27ac CUT&RUN dataset: fragment-based
   coverage checked against `bedtools genomecov -pc`, and super-enhancer calls
   compared between replicates and with and without the blacklist. There is no
   original-ROSE reference for this mode, so agreement between replicates is the
   acceptance measure.
5. **Multi-sample mode.** Running one sample through `multi` must reproduce `run`
   exactly; replicate pairs should share most calls; a sample compared with itself
   must give no gained or lost regions.
6. **Stability report.** Fixed seeds give identical output; downsampled signals
   match those from a BAM subsampled with `samtools view -s`, within sampling noise.
7. **Copy-number correction.** A synthetic test (reads in one region duplicated in
   both target and control) must return the uncorrected ranking after correction.
   On public neuroblastoma data with a known amplification, results are compared
   with LILY run on HMCan output, and the two routes are compared with each other.
8. **Stitch distance.** Per-piece sums must equal whole-region signal at every
   scanned distance. The scan row at 12.5 kb must equal a normal run. On the
   original mouse ES cell Mediator data, `--stitch auto` should land near 10 to
   15 kb; if it does not, the criterion is reworked before release. The automatic
   value should agree between replicates, and should be smaller for broad peaks
   than for narrow peaks called on the same data.
9. **Edge cases.** Zero or one super-enhancer, fewer regions than `--min-se`, empty
   chromosomes, missing index, mismatched contig names, duplicate IDs.
10. **Benchmarks.** `hyperfine` for end-to-end time against original ROSE and the
   bamliquidator path; `criterion` for the signal loop; peak memory tracked.

## 11. Milestones

| | Deliverable | Done when |
|---|---|---|
| M0 | Reference container and golden outputs | Original ROSE runs reproducibly; expected files are stored |
| M1 | Parsers, interval core, `rose stitch` | Stitched regions match golden output, including TSS modes |
| M2 | `rose signal` from BAM, single-threaded | Signal table matches golden output in compat mode |
| M3 | `rose call` and plot | Super-enhancer set and tables match; plot renders |
| M4 | `rose genes` and `rose run` | Full pipeline matches end to end |
| M5 | Parallelism and read-strategy tuning | Speed target in section 3 met or revised with measurements |
| M6 | bigWig signal, promoter subtract mode, outlier options | LILY-style run on HMCan output works |
| M7 | Multi-sample and replicate mode (`rose multi`) | Multi-sample validation in section 10 passes |
| M8 | Cutoff stability report | Stability validation passes; added time stays small next to the base run |
| M9 | Copy-number correction, both routes | Copy-number validation passes |
| M10 | Stitch distance: per-piece signal, sensitivity report, `--stitch auto` | Stitch-distance validation in section 10 passes |
| M11 | CUT&RUN preset: fragment counting and length filter, deduplication, control handling, blacklist, scale factor | The workflow in section 13 runs end to end and the CUT&RUN validation in section 10 passes |

M1 to M4 give a correct tool; speed work starts only after correctness is locked by
the golden tests.

## 12. Open decisions

- **Licence.** ROSE is MIT; LILY's scripts are GPL v2 or later. Reimplementing
  LILY's ideas from their description is compatible with an MIT licence; porting its
  code is not. Recommendation: MIT, with LILY-inspired features written from the
  method description.
- **Default behaviour.** Recommendation: clean semantics by default (CIGAR-based
  fragments, exact cutoff, leftmost naming) and the legacy quirks behind
  `--compat rose`. The alternative is legacy by default, which favours reproducing
  old papers over new analyses.
- **Read filtering defaults.** The original counts duplicates and secondary
  alignments. Recommendation: keep that under `--compat rose`, exclude them
  otherwise.
- **Paired-end ChIP-seq.** Proposed default is the real fragment, which is more
  accurate than a fixed 200 bp extension but differs from original ROSE; the fixed
  extension stays available under `--compat rose`.
- **Copy-number smoothing.** Proposed first version is a running median over bins.
  Proper segmentation, and GC-content correction as HMCan does, would be more
  accurate but are substantially more work; decide after the M9 comparison with
  LILY.
- **Automatic stitch distance.** The gap-based criterion is untested. Its exact
  form (how the expected gap distribution is built, the smoothing, the clamp range)
  is to be fixed during M10, and `auto` stays opt-in unless the validation shows it
  is stable.
- **Differential calls.** The comparison in 6.1 is descriptive (labels and fold
  change). Adding a statistical test is deliberately left to downstream tools.
- **Control handling.** ROSE subtracts raw RPM. A scaled subtraction could be added
  later; not planned for the first release.

## 13. CUT&RUN workflow

CUT&RUN support is the last milestone (M11). This section records the workflow it
has to serve and the adjustments it makes.

### Typical CUT&RUN to super-enhancer workflow today

1. **Call constituent enhancers.** Run H3K27ac (or Mediator / BRD4) CUT&RUN data
   through nf-core/cutandrun using MACS3 in broad peak mode, producing `.bed` peak
   files and coordinate-sorted `.bam` files.
2. **Convert BED to GFF.** Format the peaks as the 9-column `.gff` ROSE requires:
   column 1 is the chromosome (`chr#`), columns 4 and 5 are start and end, and
   columns 2 and 9 hold a unique peak ID.
3. **Run ROSE.** Run `ROSE_main.py` with the genome build (hg38, mm10, ...), the
   `.gff` peak file, the target `.bam`, and optionally an IgG control `.bam`.

### The same workflow with `rose-rs`

| Step today | With `rose-rs` |
|---|---|
| 1. nf-core/cutandrun with MACS3 broad peaks | Unchanged; peak calling stays outside the tool |
| 2. Convert BED to GFF | Not needed. BED and MACS3 `broadPeak` files are read directly; missing or duplicate peak IDs are generated or rejected with a clear message |
| 3. `ROSE_main.py -g ... -i peaks.gff -r target.bam [-c igg.bam]` | `rose run --assay cutrun -i peaks.broadPeak -r target.bam -a annotation -o outdir` |

### Key adjustments for CUT&RUN data

- **IgG control handling.** Unlike ChIP-seq input controls, CUT&RUN IgG controls
  have very sparse reads. Passing one through ROSE's `-c` can over-penalise
  low-background regions and add noise, so many researchers run ROSE without `-c`
  or check results both ways. (ROSE subtracts the control and clamps negatives to
  zero; it does not divide by it, so the risk is over-subtraction, not a
  division by zero.) In `rose-rs` the `cutrun` preset does not subtract the control
  unless `--subtract-control on` is given, and `--compare-control` runs both ways
  and reports how many super-enhancer calls change.
- **Deduplication.** Signal should be computed from deduplicated reads. The `cutrun`
  preset skips reads flagged as duplicates, and warns if the BAM header shows no
  duplicate-marking step, since unmarked duplicates cannot be detected from flags.
- **Fragment size filtering.** Fragments should be filtered to lengths appropriate
  for the target before signal is computed: short fragments for transcription
  factors and cofactors, nucleosome-sized fragments for histone marks.
  `--min-fragment` and `--max-fragment` apply this while counting, so the BAM does
  not need to be pre-filtered. The bounds for each target class are set by the
  user; suggested values are to be settled during M11 against real data.

