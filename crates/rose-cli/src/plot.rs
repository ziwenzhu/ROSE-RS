//! Native SVG and PNG hockey-stick super-enhancer plot renderer using `plotters` (`ab_glyph`).

use std::path::Path;
use std::sync::OnceLock;

use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters::style::register_font;
use rose_core::CallResult;
use rose_io::write_atomic_bytes;

static EMBEDDED_FONT: &[u8] = include_bytes!("../assets/DejaVuSans.ttf");
static FONT_INIT: OnceLock<()> = OnceLock::new();

fn ensure_font_registered() {
    FONT_INIT.get_or_init(|| {
        let _ = register_font("sans-serif", FontStyle::Normal, EMBEDDED_FONT);
        let _ = register_font("sans-serif", FontStyle::Bold, EMBEDDED_FONT);
    });
}

/// Renders the ROSE hockey-stick curve atomically to both `<prefix>_Plot_points.svg` and
/// `<prefix>_Plot_points.png`.
///
/// # Errors
///
/// Returns an error if drawing or atomic file output fails.
pub fn render_hockey_stick_plots(
    out_dir: &Path,
    prefix_name: &str,
    treat_label: &str,
    ctrl_label: Option<&str>,
    call: &CallResult,
) -> Result<()> {
    ensure_font_registered();

    let svg_path = out_dir.join(format!("{prefix_name}_Plot_points.svg"));
    let png_path = out_dir.join(format!("{prefix_name}_Plot_points.png"));

    // 1. Render SVG into an in-memory String, then write atomically.
    let mut svg_buf = String::new();
    {
        let backend = SVGBackend::with_string(&mut svg_buf, (720, 600));
        draw_chart(backend, treat_label, ctrl_label, call)?;
    }
    write_atomic_bytes(&svg_path, svg_buf.as_bytes())
        .with_context(|| format!("writing SVG plot '{}'", svg_path.display()))?;

    // 2. Render PNG via temporary sibling buffer on disk, then atomic rename via staging path.
    let tmp_png = out_dir.join(format!(
        ".rose_plot_tmp_{}_{prefix_name}.png",
        std::process::id()
    ));
    {
        let backend = BitMapBackend::new(&tmp_png, (720, 600));
        draw_chart(backend, treat_label, ctrl_label, call)?;
    }
    let png_bytes = std::fs::read(&tmp_png)
        .with_context(|| format!("reading staged PNG '{}'", tmp_png.display()))?;
    let _ = std::fs::remove_file(&tmp_png);
    write_atomic_bytes(&png_path, &png_bytes)
        .with_context(|| format!("writing PNG plot '{}'", png_path.display()))?;

    Ok(())
}

fn draw_chart<DB: DrawingBackend>(
    backend: DB,
    treat_label: &str,
    ctrl_label: Option<&str>,
    call: &CallResult,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    let root = backend.into_drawing_area();
    root.fill(&WHITE).map_err(|e| anyhow::anyhow!("{e:?}"))?;

    let n = call.sorted_signals_asc.len();
    let x_max = (n.max(1) as f64) * 1.03;
    let y_max_raw = call
        .sorted_signals_asc
        .last()
        .copied()
        .unwrap_or(1.0)
        .max(1.0);
    let y_max = y_max_raw * 1.08;

    let y_desc = match ctrl_label {
        Some(c) => format!("{treat_label} Signal - {c}"),
        None => format!("{treat_label} Signal"),
    };
    let caption = format!(
        "{treat_label} Enhancers (Cutoff: {:.1}, Super-Enhancers: {})",
        call.summary.cutoff_signal, call.summary.num_super_enhancers
    );

    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", 18).into_font())
        .margin(24)
        .x_label_area_size(44)
        .y_label_area_size(68)
        .build_cartesian_2d(0.0..x_max, 0.0..y_max)
        .map_err(|e| anyhow::anyhow!("{e:?}"))?;

    chart
        .configure_mesh()
        .x_desc("Enhancers Ranked by Signal (Ascending)")
        .y_desc(y_desc)
        .label_style(("sans-serif", 13).into_font())
        .light_line_style(RGBColor(240, 240, 240))
        .draw()
        .map_err(|e| anyhow::anyhow!("{e:?}"))?;

    if n > 0 {
        let cutoff = call.summary.cutoff_signal;
        let x_cut = call.summary.tangent_index_1based as f64;

        // Horizontal and vertical dashed-style reference lines at `(x_cut, cutoff)`
        chart
            .draw_series(LineSeries::new(
                vec![(0.0, cutoff), (x_max, cutoff)],
                RGBColor(150, 150, 150).stroke_width(1),
            ))
            .map_err(|e| anyhow::anyhow!("{e:?}"))?;

        chart
            .draw_series(LineSeries::new(
                vec![(x_cut, 0.0), (x_cut, y_max)],
                RGBColor(150, 150, 150).stroke_width(1),
            ))
            .map_err(|e| anyhow::anyhow!("{e:?}"))?;

        // Sliding diagonal tangent line `y = slope * x + b` clipped to `[0, y_max]`
        if call.summary.slope > 0.0 {
            let slope = call.summary.slope;
            let intercept = cutoff - slope * x_cut;
            let x0 = ((0.0 - intercept) / slope).clamp(0.0, x_max);
            let x1 = ((y_max - intercept) / slope).clamp(0.0, x_max);
            if x1 > x0 {
                chart
                    .draw_series(LineSeries::new(
                        vec![
                            (x0, (slope * x0 + intercept).max(0.0)),
                            (x1, (slope * x1 + intercept).min(y_max)),
                        ],
                        RGBColor(180, 60, 60).stroke_width(1),
                    ))
                    .map_err(|e| anyhow::anyhow!("{e:?}"))?;
            }
        }

        // Typical enhancers curve (dark slate) + super-enhancers curve (crimson red)
        let typical_pts: Vec<(f64, f64)> = call
            .sorted_signals_asc
            .iter()
            .enumerate()
            .filter(|(_, &y)| y <= cutoff)
            .map(|(i, &y)| ((i + 1) as f64, y))
            .collect();
        let super_pts: Vec<(f64, f64)> = call
            .sorted_signals_asc
            .iter()
            .enumerate()
            .filter(|(_, &y)| y > cutoff)
            .map(|(i, &y)| ((i + 1) as f64, y))
            .collect();

        if !typical_pts.is_empty() {
            chart
                .draw_series(LineSeries::new(
                    typical_pts.iter().copied(),
                    RGBColor(75, 85, 99).stroke_width(2),
                ))
                .map_err(|e| anyhow::anyhow!("{e:?}"))?;
        }
        if !super_pts.is_empty() {
            chart
                .draw_series(LineSeries::new(
                    super_pts.iter().copied(),
                    RGBColor(220, 38, 38).stroke_width(3),
                ))
                .map_err(|e| anyhow::anyhow!("{e:?}"))?;

            // Subsample up to 250 scatter circles on super-enhancers so dense cohorts stay crisp
            let step = (super_pts.len() / 250).max(1);
            chart
                .draw_series(
                    super_pts
                        .iter()
                        .step_by(step)
                        .map(|&(x, y)| Circle::new((x, y), 3, RGBColor(220, 38, 38).filled())),
                )
                .map_err(|e| anyhow::anyhow!("{e:?}"))?;
        }
    }

    root.present().map_err(|e| anyhow::anyhow!("{e:?}"))?;
    Ok(())
}
