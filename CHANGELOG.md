# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/). From 1.0.0 on, the `rose` command
line, its output file names and columns, and `run.json` are the stable public
interface. The `rose-core` and `rose-io` library APIs may still change in minor
releases.

## [Unreleased]

## [1.0.0] - 2026-11-25

First public release.

### Added

- `rose run`, which runs the whole stitch, signal, call and genes pipeline in one
  process and writes SVG and PNG hockey-stick plots and `run.json` provenance.
- `rose stitch`, `rose signal`, `rose call` and `rose genes` for each stage on
  its own, plus `rose multi` for replicate and multi-sample cohorts and
  `rose fetch-annotation` for a local annotation cache.
- Signal from BAM, CRAM (with `--reference genome.fa`) and bigWig, with
  paired-end fragment counting, MAPQ and duplicate filters, and a
  `--assay cutrun` preset.
- Copy-number correction (`--cn-correct`, `--cn-segments`), bootstrap and
  downsampling cutoff stability (`--stability`), and a 0 to 25 kb stitching
  sensitivity report on every run, with `--stitch auto`.
- `--compat rose`, which reproduces ROSE 0.1 stitching, read and region
  counting, 4-decimal rounding and the R Brent cutoff.
- Golden-fixture tests in `tests/golden/`, and CI on Linux, macOS and Windows.
- Prebuilt binaries for Linux (x86_64 and aarch64, static musl), macOS (Apple
  silicon and Intel) and Windows (x86_64), attached to each GitHub release.

### Known limitations

- `--compat rose` is checked against an independent re-implementation of the
  ROSE counting rules on synthetic data, not yet against the original ROSE
  scripts on real data.
- Mapped-read totals come from the `.bai` index rather than
  `samtools flagstat`.

[Unreleased]: https://github.com/ziwenzhu/ROSE-RS/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/ziwenzhu/ROSE-RS/releases/tag/v1.0.0
