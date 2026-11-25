//! Subcommand execution pipelines (`run`, `stitch`, `signal`, `call`, `genes`, `multi`,
//! `fetch-annotation`).

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rose_core::{
    call_super_enhancers, estimate_auto_stitch_distance, evaluate_multi_cohort, evaluate_stability,
    map_enhancers_to_genes, merge_consensus_peaks, reconcile_chrom_names, stitch_constituents,
    AnnotationIndex, AssayPreset, AutoStitchConfig, CallConfig, CompatMode, ControlSubtraction,
    CopyNumberConfig, CopyNumberProfile, CutoffMethod, ElementaryTrackValues, FragmentMode,
    GeneMapConfig, Interval, NormDenom, PiecewiseStitchIndex, RegionOrigin, RegionSignal,
    SampleSpec, StitchConfig, StitchParam, StitchedRegion, CANONICAL_STITCH_BP,
    DEFAULT_DOWNSAMPLE_FRACTIONS, DEFAULT_SENSITIVITY_DISTANCES,
};
use rose_io::{
    compute_bam_signals, inspect_alignment_header, read_annotation_index, read_blacklist_bed,
    read_cn_segments, read_constituents, read_ranked_enhancers_tsv, read_region_map_tsv,
    read_sample_sheet, sha256_file, write_enhancer_beds, write_enhancer_call_tables,
    write_gene_tables, write_multi_outputs, write_points_and_manifest, write_region_map_tsv,
    write_stitch_sensitivity_tsv, write_stitched_gff, BamSignalConfig, BigWigFile, RunManifest,
    DEFAULT_EXCLUDE_FLAGS,
};

use crate::args::{
    CallCmdArgs, Commands, FetchAnnotationArgs, GenesCmdArgs, MultiCmdArgs, RunArgs,
    SharedCallOpts, SharedSignalOpts, SharedStitchOpts, SignalCmdArgs, SignalSourceArg,
    StitchCmdArgs,
};
use crate::plot::render_hockey_stick_plots;

/// Dispatches the parsed CLI subcommand.
///
/// # Errors
///
/// Returns an `anyhow::Error` on invalid inputs, mismatched contigs, or I/O failures.
pub fn dispatch(cmd: Commands, raw_args: Vec<String>) -> Result<()> {
    match cmd {
        Commands::Run(args) => exec_run(args, raw_args),
        Commands::Stitch(args) => exec_stitch(args),
        Commands::Signal(args) => exec_signal(args),
        Commands::Call(args) => exec_call(args, raw_args),
        Commands::Genes(args) => exec_genes(args),
        Commands::Multi(args) => exec_multi(args, raw_args),
        Commands::FetchAnnotation(args) => exec_fetch_annotation(args),
    }
}

fn resolve_annotation_path(
    explicit: Option<&Path>,
    genome: Option<&str>,
) -> Result<Option<PathBuf>> {
    if let Some(p) = explicit {
        return Ok(Some(p.to_path_buf()));
    }
    let Some(g) = genome else {
        return Ok(None);
    };
    let lower = g.to_ascii_lowercase();
    for dir in ["annotations", "annotation"] {
        let cand = PathBuf::from(dir).join(format!("{lower}_refseq.ucsc"));
        if cand.exists() {
            return Ok(Some(cand));
        }
    }
    bail!(
        "annotation for genome '{g}' not found in ./annotations/; pass -a <file> or run `rose fetch-annotation -g {lower}` first"
    )
}

fn build_stitch_config(opts: &SharedStitchOpts, compat: CompatMode) -> StitchConfig {
    let exclude_chroms: Vec<String> = opts
        .exclude_chrom
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect();
    let initial_bp = match opts.stitch {
        StitchParam::Fixed(bp) => bp,
        StitchParam::Auto => CANONICAL_STITCH_BP,
    };
    StitchConfig {
        stitch_distance: initial_bp,
        tss_window: opts.tss,
        promoter_mode: opts.promoter_mode.into(),
        max_tss_genes: opts.max_tss_genes,
        exclude_chroms,
        compat,
    }
}

fn detect_is_bigwig(path: &Path, explicit: Option<SignalSourceArg>) -> bool {
    if let Some(s) = explicit {
        return s == SignalSourceArg::Bigwig;
    }
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| matches!(e.to_ascii_lowercase().as_str(), "bw" | "bigwig"))
}

fn build_bam_config(
    sig: &SharedSignalOpts,
    call: &SharedCallOpts,
    scale_override: Option<f64>,
    collect_cn_bin: Option<u64>,
    collect_downsample: bool,
) -> BamSignalConfig {
    let compat: CompatMode = call.compat.into();
    let assay: AssayPreset = sig.assay.into();

    let floor = sig.floor.unwrap_or(match assay {
        AssayPreset::Chip => 1,
        AssayPreset::CutRun => 0,
    });
    let exclude_flags = sig.exclude_flags.unwrap_or(match compat {
        CompatMode::RoseLegacy => 0x0004,
        CompatMode::Clean => DEFAULT_EXCLUDE_FLAGS,
    });
    let norm_denom = sig.norm_denom.map(NormDenom::from).unwrap_or(match compat {
        CompatMode::RoseLegacy => NormDenom::Index,
        CompatMode::Clean => NormDenom::Filtered,
    });
    let fragment_mode = match compat {
        CompatMode::RoseLegacy => FragmentMode::Seq,
        CompatMode::Clean => sig.fragment.into(),
    };
    let downsample_fractions = if collect_downsample && call.stability > 0 {
        DEFAULT_DOWNSAMPLE_FRACTIONS.to_vec()
    } else {
        Vec::new()
    };

    BamSignalConfig {
        extend: sig.extend,
        floor,
        min_mapq: sig.min_mapq,
        exclude_flags,
        fragment_mode,
        min_fragment: sig.min_fragment,
        max_fragment: sig.max_fragment,
        norm_denom,
        scale_factor: scale_override,
        compat,
        downsample_fractions,
        seed: call.seed,
        collect_cn_bin_size: collect_cn_bin,
        reference_fasta: sig.reference.clone(),
    }
}

fn build_call_config(assay: AssayPreset, has_control: bool, opts: &SharedCallOpts) -> CallConfig {
    let compat: CompatMode = opts.compat.into();
    let method = opts.cutoff.map(CutoffMethod::from).unwrap_or(match compat {
        CompatMode::RoseLegacy => CutoffMethod::RoseCompat,
        CompatMode::Clean => CutoffMethod::Exact,
    });
    let subtract_control = if !has_control {
        ControlSubtraction::Off
    } else if let Some(toggle) = opts.subtract_control {
        toggle.into()
    } else {
        match assay {
            AssayPreset::Chip => ControlSubtraction::On,
            AssayPreset::CutRun => ControlSubtraction::Off,
        }
    };
    CallConfig {
        method,
        subtract_control,
        drop_top: opts.drop_top,
        min_super_enhancers: opts.min_se,
        compat,
    }
}

struct QuantifiedTrack {
    signals: Vec<(f64, f64)>,
    elementary: ElementaryTrackValues,
    cn_bins: BTreeMap<String, Vec<u64>>,
    contig_lens: BTreeMap<String, u64>,
    read_stats: Option<(Option<u64>, Option<u64>, f64)>,
}

struct QuantifySpec<'a> {
    sig_opts: &'a SharedSignalOpts,
    call_opts: &'a SharedCallOpts,
    scale_factor: Option<f64>,
    collect_cn_bin: Option<u64>,
    collect_downsample: bool,
    emit_cutrun_warnings: bool,
}

fn quantify_track(
    track_path: &Path,
    regions: &[StitchedRegion],
    spec: &QuantifySpec<'_>,
) -> Result<QuantifiedTrack> {
    let region_chroms: Vec<String> = regions.iter().map(|r| r.chrom.clone()).collect();
    if detect_is_bigwig(track_path, spec.sig_opts.signal) {
        if spec.collect_cn_bin.is_some() {
            bail!("--cn-correct requires a BAM/CRAM control track; use --cn-segments with bigWig");
        }
        let bw = BigWigFile::open(track_path)
            .with_context(|| format!("opening bigWig '{}'", track_path.display()))?;
        let contigs = bw.contig_names();
        let chrom_map = reconcile_chrom_names(&region_chroms, &contigs, 0.50)
            .with_context(|| format!("reconciling contigs against '{}'", track_path.display()))?;
        let mut sigs = bw.compute_region_signals(regions, &chrom_map)?;
        if let Some(sf) = spec.scale_factor {
            for (d, s) in &mut sigs {
                *d *= sf;
                *s *= sf;
            }
        }
        let total_signals = sigs.iter().map(|&(_, s)| s).collect();
        return Ok(QuantifiedTrack {
            signals: sigs,
            elementary: ElementaryTrackValues::BigWig { total_signals },
            cn_bins: BTreeMap::new(),
            contig_lens: BTreeMap::new(),
            read_stats: None,
        });
    }

    let hdr = inspect_alignment_header(track_path)
        .with_context(|| format!("reading alignment header '{}'", track_path.display()))?;
    let contig_names: Vec<String> = hdr.contigs.iter().map(|(c, _)| c.clone()).collect();
    let contig_lens: BTreeMap<String, u64> = hdr.contigs.into_iter().collect();
    let chrom_map = reconcile_chrom_names(&region_chroms, &contig_names, 0.50)
        .with_context(|| format!("reconciling contigs against '{}'", track_path.display()))?;

    let bam_cfg = build_bam_config(
        spec.sig_opts,
        spec.call_opts,
        spec.scale_factor,
        spec.collect_cn_bin,
        spec.collect_downsample,
    );
    let out = compute_bam_signals(track_path, regions, &chrom_map, &bam_cfg)
        .with_context(|| format!("quantifying BAM/CRAM '{}'", track_path.display()))?;

    if spec.emit_cutrun_warnings && AssayPreset::from(spec.sig_opts.assay) == AssayPreset::CutRun {
        if !out.has_duplicate_marker_pg {
            eprintln!(
                "rose warning [cutrun]: '{}' @PG headers show no duplicate-marking step; unmarked PCR duplicates cannot be skipped from flags alone",
                track_path.display()
            );
        }
        if out.saw_single_end_reads {
            eprintln!(
                "rose warning [cutrun]: single-end alignments detected in '{}'; paired-end fragments are strongly recommended for CUT&RUN",
                track_path.display()
            );
        }
    }

    let elementary = ElementaryTrackValues::Bam {
        bp_sums: out.raw_bp_sums,
        eff_mmr: out.effective_mmr,
        scale_factor: bam_cfg.scale_factor,
        compat: bam_cfg.compat,
        tier_fractions: bam_cfg.downsample_fractions,
        tier_bp_sums: out.raw_tier_bp_sums,
    };

    Ok(QuantifiedTrack {
        signals: out.region_signals,
        elementary,
        cn_bins: out.cn_bin_counts,
        contig_lens,
        read_stats: Some((
            out.index_mapped_count,
            out.filtered_unit_count,
            out.effective_mmr,
        )),
    })
}

fn file_stem_prefix(path: &Path) -> String {
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("ROSE")
        .to_owned()
}

fn file_basename(path: &Path) -> String {
    path.file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("track")
        .to_owned()
}

fn load_blacklist_opt(path: Option<&Path>) -> Result<BTreeMap<String, Vec<Interval>>> {
    match path {
        Some(p) => {
            read_blacklist_bed(p).with_context(|| format!("reading blacklist '{}'", p.display()))
        }
        None => Ok(BTreeMap::new()),
    }
}

fn exec_stitch(args: StitchCmdArgs) -> Result<()> {
    let compat: CompatMode = args.compat.into();
    let constituents = read_constituents(&args.input, compat)
        .with_context(|| format!("reading peaks '{}'", args.input.display()))?;
    let blacklist = load_blacklist_opt(args.stitch_opts.blacklist.as_deref())?;
    let annot_opt: Option<AnnotationIndex> = match &args.annotation {
        Some(ap) => Some(
            read_annotation_index(ap, args.stitch_opts.gene_biotype.into(), compat)
                .with_context(|| format!("reading annotation '{}'", ap.display()))?,
        ),
        None => {
            if args.stitch_opts.tss > 0 {
                bail!("--tss > 0 requires -a <annotation>");
            }
            None
        }
    };
    let mut cfg = build_stitch_config(&args.stitch_opts, compat);
    if args.stitch_opts.stitch == StitchParam::Auto {
        let rep = estimate_auto_stitch_distance(
            &constituents,
            annot_opt.as_ref(),
            &blacklist,
            &cfg,
            &AutoStitchConfig::default(),
        )?;
        eprintln!("rose [stitch=auto]: {}", rep.note);
        cfg.stitch_distance = rep.selected_stitch_bp;
    }
    let (regions, _stats) =
        stitch_constituents(&constituents, annot_opt.as_ref(), &blacklist, &cfg)?;
    let out_path = if args.out.extension().is_some() {
        args.out
    } else {
        fs::create_dir_all(&args.out)?;
        let prefix = file_stem_prefix(&args.input);
        args.out.join(format!("{prefix}_STITCHED.gff"))
    };
    write_stitched_gff(&out_path, &regions, compat)?;
    Ok(())
}

fn exec_signal(args: SignalCmdArgs) -> Result<()> {
    let compat: CompatMode = args.compat.into();
    let constituents = read_constituents(&args.input, compat)?;
    let regions: Vec<StitchedRegion> = constituents
        .into_iter()
        .map(|c| StitchedRegion {
            constituent_size: c.interval.width(),
            constituent_ids: vec![c.id.clone()],
            id: c.id,
            chrom: c.chrom,
            interval: c.interval,
            strand: c.strand,
            num_loci: 1,
            origin: RegionOrigin::Enhancer,
        })
        .collect();

    let dummy_call = SharedCallOpts {
        subtract_control: None,
        compare_control: false,
        cutoff: None,
        drop_top: 0,
        min_se: 0,
        stability: 0,
        seed: 42,
        cn_correct: false,
        cn_segments: None,
        cn_bin: 100_000,
        compat: args.compat,
    };

    let treat_q = quantify_track(
        &args.rankby,
        &regions,
        &QuantifySpec {
            sig_opts: &args.signal_opts,
            call_opts: &dummy_call,
            scale_factor: args.signal_opts.scale_factor,
            collect_cn_bin: None,
            collect_downsample: false,
            emit_cutrun_warnings: true,
        },
    )?;
    let ctrl_q = match &args.control {
        Some(cp) => Some(quantify_track(
            cp,
            &regions,
            &QuantifySpec {
                sig_opts: &args.signal_opts,
                call_opts: &dummy_call,
                scale_factor: args.signal_opts.control_scale_factor,
                collect_cn_bin: None,
                collect_downsample: false,
                emit_cutrun_warnings: false,
            },
        )?),
        None => None,
    };

    let items: Vec<RegionSignal> = regions
        .into_iter()
        .enumerate()
        .map(|(i, region)| RegionSignal {
            region,
            treat_density: treat_q.signals[i].0,
            treat_signal: treat_q.signals[i].1,
            ctrl_density: ctrl_q.as_ref().map(|c| c.signals[i].0),
            ctrl_signal: ctrl_q.as_ref().map(|c| c.signals[i].1),
            cn_ratio: None,
        })
        .collect();

    let treat_name = file_basename(&args.rankby);
    let ctrl_name = args.control.as_deref().map(file_basename);
    write_region_map_tsv(&args.out, &items, &treat_name, ctrl_name.as_deref(), compat)?;
    Ok(())
}

fn exec_call(args: CallCmdArgs, raw_args: Vec<String>) -> Result<()> {
    let compat: CompatMode = args.call_opts.compat.into();
    let (mut items, treat_name, ctrl_name) = read_region_map_tsv(&args.input)?;
    if let Some(seg_path) = &args.call_opts.cn_segments {
        let segs = read_cn_segments(seg_path)?;
        let prof = CopyNumberProfile::from_segments(segs, 0.25, 20.0);
        for it in &mut items {
            it.cn_ratio = Some(prof.ratio_for_region(&it.region));
        }
    }

    let prefix = args
        .name
        .clone()
        .unwrap_or_else(|| file_stem_prefix(&args.input));
    let call_cfg = build_call_config(AssayPreset::Chip, ctrl_name.is_some(), &args.call_opts);
    let call_res = call_super_enhancers(&items, &call_cfg);

    let stab_opt = (args.call_opts.stability > 0).then(|| {
        let ordered_sigs: Vec<f64> = call_res.ranked.iter().map(|r| r.rank_signal).collect();
        evaluate_stability(
            &ordered_sigs,
            args.call_opts.stability,
            args.call_opts.seed,
            call_cfg.method,
            call_cfg.drop_top,
            compat,
            None,
        )
    });

    fs::create_dir_all(&args.out)?;
    let all_tbl = args.out.join(format!("{prefix}_AllEnhancers.table.txt"));
    let sup_tbl = args.out.join(format!("{prefix}_SuperEnhancers.table.txt"));
    let has_cn = args.call_opts.cn_segments.is_some() || items.iter().any(|x| x.cn_ratio.is_some());
    write_enhancer_call_tables(
        &all_tbl,
        &sup_tbl,
        &prefix,
        &args.input.display().to_string(),
        &treat_name,
        ctrl_name.as_deref(),
        &call_res,
        stab_opt.as_ref(),
        has_cn,
        compat,
    )?;
    write_enhancer_beds(&args.out, &prefix, &call_res)?;
    render_hockey_stick_plots(
        &args.out,
        &prefix,
        &treat_name,
        ctrl_name.as_deref(),
        &call_res,
    )?;

    let pts_tsv = args.out.join(format!("{prefix}_Plot_points.tsv"));
    let sat_path = args.out.join(format!("{prefix}_Stability_saturation.tsv"));
    let manifest = RunManifest {
        version: env!("CARGO_PKG_VERSION").to_owned(),
        command: "call".to_owned(),
        arguments: raw_args,
        input_checksums: vec![(
            "region_map".to_owned(),
            args.input.display().to_string(),
            sha256_file(&args.input)?,
        )],
        ranking_read_stats: None,
        control_read_stats: None,
        cutoff_signal: Some(call_res.summary.cutoff_signal),
        num_regions: call_res.summary.num_regions,
        num_super_enhancers: Some(call_res.summary.num_super_enhancers),
        control_comparison_note: None,
        active_stitch_bp: None,
        auto_stitch_report: None,
    };
    write_points_and_manifest(
        &pts_tsv,
        &call_res,
        stab_opt.as_ref().map(|s| (sat_path.as_path(), s)),
        &args.out.join("run.json"),
        &manifest,
    )?;
    Ok(())
}

fn exec_genes(args: GenesCmdArgs) -> Result<()> {
    let compat: CompatMode = args.compat.into();
    let ranked = read_ranked_enhancers_tsv(&args.input)?;
    let annot = read_annotation_index(&args.annotation, args.gene_biotype.into(), compat)?;
    let g_cfg = GeneMapConfig {
        search_window: args.window,
        unique_genes: true,
        compat,
    };
    let (enh_rows, gene_rows) = map_enhancers_to_genes(&ranked, &annot, &g_cfg);
    fs::create_dir_all(&args.out)?;
    let prefix = file_stem_prefix(&args.input);
    let suffix = if args.window == 50_000 {
        String::new()
    } else {
        format!("_{}KB", args.window / 1_000)
    };
    write_gene_tables(
        args.out
            .join(format!("{prefix}_ENHANCER_TO_GENE{suffix}.txt")),
        args.out
            .join(format!("{prefix}_GENE_TO_ENHANCER{suffix}.txt")),
        &enh_rows,
        &gene_rows,
        compat,
    )?;
    Ok(())
}

fn exec_run(args: RunArgs, raw_args: Vec<String>) -> Result<()> {
    let compat: CompatMode = args.call_opts.compat.into();
    let assay: AssayPreset = args.signal_opts.assay.into();

    if assay == AssayPreset::CutRun && args.stitch_opts.blacklist.is_none() {
        eprintln!(
            "rose warning [cutrun]: no --blacklist supplied; high-signal telomeric/centromeric artifacts may rank at the top"
        );
    }

    let annot_path_opt =
        resolve_annotation_path(args.annotation.as_deref(), args.genome.as_deref())?;
    let annot_opt = match &annot_path_opt {
        Some(ap) => Some(
            read_annotation_index(ap, args.stitch_opts.gene_biotype.into(), compat)
                .with_context(|| format!("reading annotation '{}'", ap.display()))?,
        ),
        None => {
            if args.stitch_opts.tss > 0 {
                bail!("--tss > 0 requires -a <annotation> or -g <genome>");
            }
            None
        }
    };

    let blacklist = load_blacklist_opt(args.stitch_opts.blacklist.as_deref())?;
    let constituents = read_constituents(&args.input, compat)
        .with_context(|| format!("reading constituent peaks '{}'", args.input.display()))?;
    let stitch_cfg = build_stitch_config(&args.stitch_opts, compat);
    let sweep_horizon = if args.stitch_opts.no_stitch_report {
        0
    } else {
        25_000
    };
    let (piecewise_idx, active_stitch_bp) = PiecewiseStitchIndex::build(
        &constituents,
        annot_opt.as_ref(),
        &blacklist,
        args.stitch_opts.stitch,
        &stitch_cfg,
        sweep_horizon,
    )?;
    if let Some(rep) = piecewise_idx.auto_report() {
        eprintln!("rose [stitch=auto]: {}", rep.note);
    }

    let collect_cn_bin = if args.call_opts.cn_correct {
        if args.control.is_none() {
            bail!("--cn-correct requires a control BAM (-c <control.bam>)");
        }
        Some(args.call_opts.cn_bin)
    } else {
        None
    };

    // Under `--compat rose`, ROSE counts only bases strictly inside each stitched locus, which is
    // not additive across peak/gap pieces. Quantify the active regions directly in the same pass
    // (appended after the pieces) so the primary tables match ROSE exactly; the sensitivity sweep
    // and stability tiers keep using the pieces.
    let direct_regions: Vec<StitchedRegion> = if compat == CompatMode::RoseLegacy {
        let mut active_cfg = stitch_cfg.clone();
        active_cfg.stitch_distance = active_stitch_bp;
        stitch_constituents(&constituents, annot_opt.as_ref(), &blacklist, &active_cfg)?.0
    } else {
        Vec::new()
    };
    let n_elem = piecewise_idx.elementary_regions().len();
    let mut quant_regs: Vec<StitchedRegion> = piecewise_idx.elementary_regions().to_vec();
    quant_regs.extend(direct_regions.iter().cloned());
    let elem_regs = quant_regs.as_slice();
    let treat_q = quantify_track(
        &args.rankby,
        elem_regs,
        &QuantifySpec {
            sig_opts: &args.signal_opts,
            call_opts: &args.call_opts,
            scale_factor: args.signal_opts.scale_factor,
            collect_cn_bin: None,
            collect_downsample: args.call_opts.stability > 0,
            emit_cutrun_warnings: true,
        },
    )?;

    let ctrl_q = match &args.control {
        Some(cp) => Some(quantify_track(
            cp,
            elem_regs,
            &QuantifySpec {
                sig_opts: &args.signal_opts,
                call_opts: &args.call_opts,
                scale_factor: args.signal_opts.control_scale_factor,
                collect_cn_bin,
                collect_downsample: args.call_opts.stability > 0,
                emit_cutrun_warnings: false,
            },
        )?),
        None => None,
    };

    let cn_profile_opt: Option<CopyNumberProfile> = if args.call_opts.cn_correct {
        let Some(cq) = ctrl_q.as_ref() else {
            bail!("--cn-correct requires a control BAM (-c <control.bam>)");
        };
        let cn_cfg = CopyNumberConfig {
            bin_size: args.call_opts.cn_bin,
            ..CopyNumberConfig::default()
        };
        Some(CopyNumberProfile::from_bin_counts(
            &cq.cn_bins,
            &cq.contig_lens,
            &blacklist,
            &cn_cfg,
        )?)
    } else if let Some(seg_path) = &args.call_opts.cn_segments {
        let segs = read_cn_segments(seg_path)?;
        Some(CopyNumberProfile::from_segments(segs, 0.25, 20.0))
    } else {
        None
    };

    let mat = piecewise_idx.materialize_at_distance(
        active_stitch_bp,
        annot_opt.as_ref(),
        &treat_q.elementary,
        ctrl_q.as_ref().map(|c| &c.elementary),
        cn_profile_opt.as_ref(),
    )?;
    let regions = mat.regions;
    let mut items = mat.signals;
    if !direct_regions.is_empty() {
        let direct_idx: BTreeMap<&str, usize> = direct_regions
            .iter()
            .enumerate()
            .map(|(i, r)| (r.id.as_str(), n_elem + i))
            .collect();
        for it in &mut items {
            let Some(&qi) = direct_idx.get(it.region.id.as_str()) else {
                bail!(
                    "internal error: region '{}' missing from direct quantification",
                    it.region.id
                );
            };
            (it.treat_density, it.treat_signal) = treat_q.signals[qi];
            if let Some(cq) = &ctrl_q {
                it.ctrl_density = Some(cq.signals[qi].0);
                it.ctrl_signal = Some(cq.signals[qi].1);
            }
        }
    }
    let treat_downsampled = mat.downsampled_totals;
    let ctrl_downsampled = if args.call_opts.stability > 0 {
        if let Some(cq) = &ctrl_q {
            Some(
                piecewise_idx
                    .materialize_at_distance(
                        active_stitch_bp,
                        annot_opt.as_ref(),
                        &cq.elementary,
                        None,
                        None,
                    )?
                    .downsampled_totals,
            )
        } else {
            None
        }
    } else {
        None
    };

    let call_cfg = build_call_config(assay, ctrl_q.is_some(), &args.call_opts);
    let call_res = call_super_enhancers(&items, &call_cfg);

    let sensitivity_opt = if args.stitch_opts.no_stitch_report {
        None
    } else {
        Some(piecewise_idx.evaluate_sensitivity_sweep(
            &DEFAULT_SENSITIVITY_DISTANCES,
            active_stitch_bp,
            annot_opt.as_ref(),
            &treat_q.elementary,
            ctrl_q.as_ref().map(|c| &c.elementary),
            cn_profile_opt.as_ref(),
            &call_cfg,
        )?)
    };

    let compare_note = if args.call_opts.compare_control && ctrl_q.is_some() {
        let mut alt_cfg = call_cfg.clone();
        alt_cfg.subtract_control = match call_cfg.subtract_control {
            ControlSubtraction::On => ControlSubtraction::Off,
            ControlSubtraction::Off => ControlSubtraction::On,
        };
        let alt_res = call_super_enhancers(&items, &alt_cfg);
        let cur_ids: std::collections::HashSet<&str> = call_res
            .ranked
            .iter()
            .filter(|r| r.is_super)
            .map(|r| r.item.region.id.as_str())
            .collect();
        let alt_ids: std::collections::HashSet<&str> = alt_res
            .ranked
            .iter()
            .filter(|r| r.is_super)
            .map(|r| r.item.region.id.as_str())
            .collect();
        let shared = cur_ids.intersection(&alt_ids).count();
        let diff = cur_ids.symmetric_difference(&alt_ids).count();
        Some(format!(
            "primary={} SEs, alternate={} SEs, shared={shared}, changed={diff}",
            cur_ids.len(),
            alt_ids.len()
        ))
    } else {
        None
    };

    let stab_opt = (args.call_opts.stability > 0).then(|| {
        // Align downsampled tiers to `call_res.ranked` order (`enhancer_rank` ascending).
        let id_to_orig: BTreeMap<&str, usize> = regions
            .iter()
            .enumerate()
            .map(|(idx, r)| (r.id.as_str(), idx))
            .collect();
        let ordered_sigs: Vec<f64> = call_res.ranked.iter().map(|r| r.rank_signal).collect();
        let tiers: Vec<(f64, Vec<f64>)> = treat_downsampled
            .iter()
            .enumerate()
            .map(|(t_idx, (frac, t_sigs))| {
                let c_sigs_opt = ctrl_downsampled
                    .as_ref()
                    .and_then(|cd| cd.get(t_idx))
                    .map(|(_, v)| v.as_slice());
                let aligned: Vec<f64> = call_res
                    .ranked
                    .iter()
                    .map(|rk| {
                        let orig_i = id_to_orig[rk.item.region.id.as_str()];
                        let ts = t_sigs[orig_i];
                        let cs = match call_cfg.subtract_control {
                            ControlSubtraction::On => c_sigs_opt.map_or(0.0, |cv| cv[orig_i]),
                            ControlSubtraction::Off => 0.0,
                        };
                        let cn = rk.item.cn_ratio.unwrap_or(1.0).max(1e-6);
                        (ts - cs).max(0.0) / cn
                    })
                    .collect();
                (*frac, aligned)
            })
            .collect();

        evaluate_stability(
            &ordered_sigs,
            args.call_opts.stability,
            args.call_opts.seed,
            call_cfg.method,
            call_cfg.drop_top,
            compat,
            (!tiers.is_empty()).then_some(tiers.as_slice()),
        )
    });

    fs::create_dir_all(&args.out)?;
    let gff_dir = args.out.join("gff");
    fs::create_dir_all(&gff_dir)?;

    let prefix = file_stem_prefix(&args.input);
    let stitched_tag = if args.stitch_opts.tss > 0 {
        format!(
            "{prefix}_{}KB_STITCHED_TSS_DISTAL",
            active_stitch_bp / 1_000
        )
    } else {
        format!("{prefix}_{}KB_STITCHED", active_stitch_bp / 1_000)
    };
    if let Some(sens_rows) = &sensitivity_opt {
        write_stitch_sensitivity_tsv(
            args.out.join(format!("{prefix}_Stitch_sensitivity.tsv")),
            sens_rows,
        )?;
    }

    write_stitched_gff(
        gff_dir.join(format!("{stitched_tag}.gff")),
        &regions,
        compat,
    )?;

    let treat_name = file_basename(&args.rankby);
    let ctrl_name = args.control.as_deref().map(file_basename);
    let reg_map_path = args
        .out
        .join(format!("{stitched_tag}_ENHANCER_REGION_MAP.txt"));
    write_region_map_tsv(
        &reg_map_path,
        &items,
        &treat_name,
        ctrl_name.as_deref(),
        compat,
    )?;

    let all_tbl = args.out.join(format!("{prefix}_AllEnhancers.table.txt"));
    let sup_tbl = args.out.join(format!("{prefix}_SuperEnhancers.table.txt"));
    write_enhancer_call_tables(
        &all_tbl,
        &sup_tbl,
        &prefix,
        &reg_map_path.display().to_string(),
        &treat_name,
        ctrl_name.as_deref(),
        &call_res,
        stab_opt.as_ref(),
        cn_profile_opt.is_some(),
        compat,
    )?;
    write_enhancer_beds(&args.out, &prefix, &call_res)?;
    render_hockey_stick_plots(
        &args.out,
        &prefix,
        &treat_name,
        ctrl_name.as_deref(),
        &call_res,
    )?;

    if let Some(annot) = &annot_opt {
        let g_cfg = GeneMapConfig {
            search_window: args.gene_window,
            unique_genes: true,
            compat,
        };
        let super_only: Vec<_> = call_res
            .ranked
            .iter()
            .filter(|r| r.is_super)
            .cloned()
            .collect();
        let (enh_rows, gene_rows) = map_enhancers_to_genes(&super_only, annot, &g_cfg);
        write_gene_tables(
            args.out
                .join(format!("{prefix}_SuperEnhancers_ENHANCER_TO_GENE.txt")),
            args.out
                .join(format!("{prefix}_SuperEnhancers_GENE_TO_ENHANCER.txt")),
            &enh_rows,
            &gene_rows,
            compat,
        )?;
        let (all_enh_rows, all_gene_rows) = map_enhancers_to_genes(&call_res.ranked, annot, &g_cfg);
        write_gene_tables(
            args.out
                .join(format!("{prefix}_AllEnhancers_ENHANCER_TO_GENE.txt")),
            args.out
                .join(format!("{prefix}_AllEnhancers_GENE_TO_ENHANCER.txt")),
            &all_enh_rows,
            &all_gene_rows,
            compat,
        )?;
    }

    let mut checksums = vec![
        (
            "peaks".to_owned(),
            args.input.display().to_string(),
            sha256_file(&args.input)?,
        ),
        (
            "rankby".to_owned(),
            args.rankby.display().to_string(),
            sha256_file(&args.rankby)?,
        ),
    ];
    if let Some(cp) = &args.control {
        checksums.push((
            "control".to_owned(),
            cp.display().to_string(),
            sha256_file(cp)?,
        ));
    }
    if let Some(ap) = &annot_path_opt {
        checksums.push((
            "annotation".to_owned(),
            ap.display().to_string(),
            sha256_file(ap)?,
        ));
    }

    let pts_tsv = args.out.join(format!("{prefix}_Plot_points.tsv"));
    let sat_path = args.out.join(format!("{prefix}_Stability_saturation.tsv"));
    let manifest = RunManifest {
        version: env!("CARGO_PKG_VERSION").to_owned(),
        command: "run".to_owned(),
        arguments: raw_args,
        input_checksums: checksums,
        ranking_read_stats: treat_q.read_stats,
        control_read_stats: ctrl_q.and_then(|c| c.read_stats),
        cutoff_signal: Some(call_res.summary.cutoff_signal),
        num_regions: call_res.summary.num_regions,
        num_super_enhancers: Some(call_res.summary.num_super_enhancers),
        control_comparison_note: compare_note,
        active_stitch_bp: Some(active_stitch_bp),
        auto_stitch_report: piecewise_idx.auto_report().cloned(),
    };
    write_points_and_manifest(
        &pts_tsv,
        &call_res,
        stab_opt.as_ref().map(|s| (sat_path.as_path(), s)),
        &args.out.join("run.json"),
        &manifest,
    )?;
    Ok(())
}

fn exec_multi(args: MultiCmdArgs, raw_args: Vec<String>) -> Result<()> {
    let compat: CompatMode = args.call_opts.compat.into();
    let assay: AssayPreset = args.signal_opts.assay.into();
    let sheet = read_sample_sheet(&args.samples)?;
    if sheet.is_empty() {
        bail!(
            "sample sheet '{}' contains zero samples",
            args.samples.display()
        );
    }

    let blacklist = load_blacklist_opt(args.stitch_opts.blacklist.as_deref())?;
    let annot_opt = match &args.annotation {
        Some(ap) => Some(read_annotation_index(
            ap,
            args.stitch_opts.gene_biotype.into(),
            compat,
        )?),
        None => {
            if args.stitch_opts.tss > 0 {
                bail!("--tss > 0 requires -a <annotation>");
            }
            None
        }
    };

    let per_sample_peaks: Vec<_> = sheet
        .iter()
        .map(|row| read_constituents(&row.peaks, compat))
        .collect::<std::result::Result<Vec<_>, _>>()?;

    // When `sheet.len() == 1` and `--min-peak-samples 1`, preserve original constituent IDs if
    // non-overlapping so 1-sample `rose multi` matches `rose run` identically!
    let pooled_constituents = if sheet.len() == 1 && args.min_peak_samples == 1 {
        per_sample_peaks.into_iter().next().unwrap_or_default()
    } else {
        merge_consensus_peaks(&per_sample_peaks, args.min_peak_samples)?
    };

    let stitch_cfg = build_stitch_config(&args.stitch_opts, compat);
    let sweep_horizon = if args.stitch_opts.no_stitch_report {
        0
    } else {
        25_000
    };
    let (piecewise_idx, active_stitch_bp) = PiecewiseStitchIndex::build(
        &pooled_constituents,
        annot_opt.as_ref(),
        &blacklist,
        args.stitch_opts.stitch,
        &stitch_cfg,
        sweep_horizon,
    )?;
    if let Some(rep) = piecewise_idx.auto_report() {
        eprintln!("rose multi [stitch=auto]: {}", rep.note);
    }

    let elem_regs = piecewise_idx.elementary_regions();
    let any_ctrl = sheet.iter().any(|r| r.control.is_some());
    let call_cfg = build_call_config(assay, any_ctrl, &args.call_opts);

    let per_sample_tracks: Result<Vec<(QuantifiedTrack, Option<QuantifiedTrack>)>> = sheet
        .par_iter()
        .map(|row| {
            let tq = quantify_track(
                &row.bam,
                elem_regs,
                &QuantifySpec {
                    sig_opts: &args.signal_opts,
                    call_opts: &args.call_opts,
                    scale_factor: row.scale_factor.or(args.signal_opts.scale_factor),
                    collect_cn_bin: None,
                    collect_downsample: false,
                    emit_cutrun_warnings: false,
                },
            )?;
            let cq = match &row.control {
                Some(cp) => Some(quantify_track(
                    cp,
                    elem_regs,
                    &QuantifySpec {
                        sig_opts: &args.signal_opts,
                        call_opts: &args.call_opts,
                        scale_factor: args.signal_opts.control_scale_factor,
                        collect_cn_bin: None,
                        collect_downsample: false,
                        emit_cutrun_warnings: false,
                    },
                )?),
                None => None,
            };
            Ok((tq, cq))
        })
        .collect();

    let sample_tracks = per_sample_tracks?;
    let mut shared_regions = Vec::new();
    let mut sig_vecs = Vec::with_capacity(sample_tracks.len());
    for (idx, (tq, cq)) in sample_tracks.iter().enumerate() {
        let mat = piecewise_idx.materialize_at_distance(
            active_stitch_bp,
            annot_opt.as_ref(),
            &tq.elementary,
            cq.as_ref().map(|c| &c.elementary),
            None,
        )?;
        if idx == 0 {
            shared_regions = mat.regions;
        }
        sig_vecs.push(mat.signals);
    }

    let specs: Vec<SampleSpec> = sheet
        .iter()
        .map(|r| SampleSpec {
            sample_id: r.sample.clone(),
            group: r.group.clone(),
            scale_factor: r.scale_factor.or(args.signal_opts.scale_factor),
        })
        .collect();

    let cohort_res = evaluate_multi_cohort(
        &shared_regions,
        &specs,
        &sig_vecs,
        &call_cfg,
        args.min_replicates,
    )?;

    fs::create_dir_all(&args.out)?;
    write_multi_outputs(&args.out, &cohort_res)?;

    // Also emit per-sample ROSE tables, stitch sensitivity, & hockey-stick plots under `<out>/<sample>/`
    for (s_idx, row) in sheet.iter().enumerate() {
        let s_dir = args.out.join(&row.sample);
        fs::create_dir_all(&s_dir)?;
        let call_res = &cohort_res.per_sample_calls[s_idx];
        let treat_name = file_basename(&row.bam);
        let ctrl_name = row.control.as_deref().map(file_basename);
        write_enhancer_call_tables(
            s_dir.join(format!("{}_AllEnhancers.table.txt", row.sample)),
            s_dir.join(format!("{}_SuperEnhancers.table.txt", row.sample)),
            &row.sample,
            &args.samples.display().to_string(),
            &treat_name,
            ctrl_name.as_deref(),
            call_res,
            None,
            false,
            compat,
        )?;
        render_hockey_stick_plots(
            &s_dir,
            &row.sample,
            &treat_name,
            ctrl_name.as_deref(),
            call_res,
        )?;
        if !args.stitch_opts.no_stitch_report {
            let (tq, cq) = &sample_tracks[s_idx];
            let sens = piecewise_idx.evaluate_sensitivity_sweep(
                &DEFAULT_SENSITIVITY_DISTANCES,
                active_stitch_bp,
                annot_opt.as_ref(),
                &tq.elementary,
                cq.as_ref().map(|c| &c.elementary),
                None,
                &call_cfg,
            )?;
            write_stitch_sensitivity_tsv(
                s_dir.join(format!("{}_Stitch_sensitivity.tsv", row.sample)),
                &sens,
            )?;
        }
    }

    let manifest = RunManifest {
        version: env!("CARGO_PKG_VERSION").to_owned(),
        command: "multi".to_owned(),
        arguments: raw_args,
        input_checksums: vec![(
            "samples".to_owned(),
            args.samples.display().to_string(),
            sha256_file(&args.samples)?,
        )],
        ranking_read_stats: None,
        control_read_stats: None,
        cutoff_signal: None,
        num_regions: shared_regions.len(),
        num_super_enhancers: None,
        control_comparison_note: None,
        active_stitch_bp: Some(active_stitch_bp),
        auto_stitch_report: piecewise_idx.auto_report().cloned(),
    };
    let json_bytes = serde_json::to_vec_pretty(&manifest)?;
    rose_io::write_atomic_bytes(args.out.join("run.json"), &json_bytes)?;
    Ok(())
}

fn exec_fetch_annotation(args: FetchAnnotationArgs) -> Result<()> {
    let lower = args.genome.to_ascii_lowercase();
    if !matches!(
        lower.as_str(),
        "hg18" | "hg19" | "hg38" | "mm8" | "mm9" | "mm10"
    ) {
        bail!(
            "unsupported genome build '{}'; expected one of: hg18, hg19, hg38, mm8, mm9, mm10",
            args.genome
        );
    }
    fs::create_dir_all(&args.out)?;
    let target = args.out.join(format!("{lower}_refseq.ucsc"));

    if let Some(src_dir) = &args.from_dir {
        let src = src_dir.join(format!("{lower}_refseq.ucsc"));
        if !src.exists() {
            bail!("source annotation file '{}' does not exist", src.display());
        }
        let bytes = fs::read(&src)?;
        rose_io::write_atomic_bytes(&target, &bytes)?;
        println!("Installed {} -> {}", src.display(), target.display());
        return Ok(());
    }

    bail!(
        "offline build has no bundled HTTP client by design; download https://hgdownload.soe.ucsc.edu/goldenPath/{lower}/database/refGene.txt.gz, gunzip to {lower}_refseq.ucsc, and pass `--from-dir <dir>` or `-a {}`",
        target.display()
    )
}
