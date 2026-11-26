//! Native SVG and PNG hockey-stick super-enhancer plot renderer using `plotters` (`ab_glyph`).

use std::path::Path;
use std::sync::OnceLock;

use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters::style::register_font;
use plotters::style::text_anchor::{HPos, Pos, VPos};
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

// Palette shared with the README figures: dark navy text and typical enhancers, coral-red
// super-enhancers.
const NAVY: RGBColor = RGBColor(0x17, 0x2B, 0x46);
const CORAL: RGBColor = RGBColor(0xE2, 0x65, 0x66);
const GUIDE: RGBColor = RGBColor(0x8A, 0x94, 0xA3);
const GRID: RGBColor = RGBColor(0xE6, 0xE9, 0xEE);

/// Strips common alignment / track extensions so plot text reads `H3K27ac`, not `H3K27ac.bam`.
fn display_label(name: &str) -> &str {
    const EXTS: [&str; 6] = [".bam", ".cram", ".bigwig", ".bigWig", ".bw", ".sam"];
    EXTS.iter()
        .find_map(|ext| name.strip_suffix(ext))
        .filter(|s| !s.is_empty())
        .unwrap_or(name)
}

/// Formats an axis value as a whole number with thousands separators (`12,500`).
fn fmt_thousands(v: f64) -> String {
    let n = v.round() as i64;
    let digits = n.unsigned_abs().to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3 + 1);
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    if n < 0 {
        out.insert(0, '-');
    }
    out
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
    let err = |e| anyhow::anyhow!("{e:?}");
    let root = backend.into_drawing_area();
    root.fill(&WHITE).map_err(err)?;

    let n = call.sorted_signals_asc.len();
    let x_max = (n.max(1) as f64) * 1.03;
    let y_max_raw = call
        .sorted_signals_asc
        .last()
        .copied()
        .unwrap_or(1.0)
        .max(1.0);
    let y_max = y_max_raw * 1.08;

    let treat = display_label(treat_label);
    let y_desc = match ctrl_label {
        Some(c) => format!("{treat} signal \u{2212} {}", display_label(c)),
        None => format!("{treat} signal"),
    };
    let n_super = call.summary.num_super_enhancers;
    let caption = format!(
        "{treat} enhancers: {n_super} super-enhancer{}",
        if n_super == 1 { "" } else { "s" }
    );

    let text = |size: u32| ("sans-serif", size).into_font().color(&NAVY);

    let mut chart = ChartBuilder::on(&root)
        .caption(caption, text(18))
        .margin(24)
        .x_label_area_size(44)
        .y_label_area_size(92)
        .build_cartesian_2d(0.0..x_max, 0.0..y_max)
        .map_err(err)?;

    if n > 0 {
        // Light coral band over the super-enhancer ranks, drawn first so the grid and curve sit on top.
        let x_cut = call.summary.tangent_index_1based as f64;
        chart
            .draw_series(std::iter::once(Rectangle::new(
                [(x_cut, 0.0), (x_max, y_max)],
                CORAL.mix(0.10).filled(),
            )))
            .map_err(err)?;
    }

    chart
        .configure_mesh()
        .x_desc("Enhancers ranked by signal")
        .y_desc(y_desc)
        .label_style(text(13))
        .axis_desc_style(text(14))
        .axis_style(NAVY)
        .bold_line_style(GRID)
        .light_line_style(TRANSPARENT)
        .x_label_formatter(&|v| fmt_thousands(*v))
        .y_label_formatter(&|v| fmt_thousands(*v))
        .draw()
        .map_err(err)?;

    if n > 0 {
        let cutoff = call.summary.cutoff_signal;
        let x_cut = call.summary.tangent_index_1based as f64;
        let guide = GUIDE.stroke_width(1);

        // Dashed reference lines through the tangent point `(x_cut, cutoff)`.
        chart
            .draw_series(DashedLineSeries::new(
                vec![(0.0, cutoff), (x_max, cutoff)],
                6,
                4,
                guide,
            ))
            .map_err(err)?;
        chart
            .draw_series(DashedLineSeries::new(
                vec![(x_cut, 0.0), (x_cut, y_max)],
                6,
                4,
                guide,
            ))
            .map_err(err)?;

        // Tangent line `y = slope * x + b` clipped to `[0, y_max]`.
        if call.summary.slope > 0.0 {
            let slope = call.summary.slope;
            let intercept = cutoff - slope * x_cut;
            let x0 = ((0.0 - intercept) / slope).clamp(0.0, x_max);
            let x1 = ((y_max - intercept) / slope).clamp(0.0, x_max);
            if x1 > x0 {
                chart
                    .draw_series(DashedLineSeries::new(
                        vec![
                            (x0, (slope * x0 + intercept).max(0.0)),
                            (x1, (slope * x1 + intercept).min(y_max)),
                        ],
                        3,
                        3,
                        CORAL.stroke_width(1),
                    ))
                    .map_err(err)?;
            }
        }

        // Labels for the guides and the shaded zone, anchored to data points and offset in pixels.
        let label =
            |txt: String, at: (f64, f64), dx: i32, dy: i32, color: RGBColor, h: HPos, v: VPos| {
                EmptyElement::at(at)
                    + Text::new(
                        txt,
                        (dx, dy),
                        ("sans-serif", 12)
                            .into_font()
                            .color(&color)
                            .pos(Pos::new(h, v)),
                    )
            };
        let band_px = chart.backend_coord(&(x_max, 0.0)).0 - chart.backend_coord(&(x_cut, 0.0)).0;
        chart
            .draw_series([
                label(
                    format!("Cutoff {}", fmt_thousands(cutoff)),
                    (0.0, cutoff),
                    6,
                    -4,
                    NAVY,
                    HPos::Left,
                    VPos::Bottom,
                ),
                label(
                    format!("Rank {}", fmt_thousands(x_cut)),
                    (x_cut, y_max),
                    -6,
                    6,
                    NAVY,
                    HPos::Right,
                    VPos::Top,
                ),
                // Inside the band when it is wide enough, otherwise stacked under the rank label.
                if band_px >= 110 {
                    label(
                        "Super-enhancers".to_owned(),
                        (x_cut, y_max),
                        6,
                        6,
                        CORAL,
                        HPos::Left,
                        VPos::Top,
                    )
                } else {
                    label(
                        "Super-enhancers \u{2192}".to_owned(),
                        (x_cut, y_max),
                        -6,
                        22,
                        CORAL,
                        HPos::Right,
                        VPos::Top,
                    )
                },
            ])
            .map_err(err)?;

        // Typical enhancers curve (navy) + super-enhancers curve (coral).
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
                    NAVY.stroke_width(2),
                ))
                .map_err(err)?;
        }
        if !super_pts.is_empty() {
            chart
                .draw_series(LineSeries::new(
                    super_pts.iter().copied(),
                    CORAL.stroke_width(3),
                ))
                .map_err(err)?;

            // Subsample up to 250 scatter circles on super-enhancers so dense cohorts stay crisp
            let step = (super_pts.len() / 250).max(1);
            chart
                .draw_series(
                    super_pts
                        .iter()
                        .step_by(step)
                        .map(|&(x, y)| Circle::new((x, y), 3, CORAL.filled())),
                )
                .map_err(err)?;
        }
    }

    root.present().map_err(err)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{display_label, fmt_thousands};

    #[test]
    fn labels_drop_track_extensions() {
        assert_eq!(display_label("H3K27ac.bam"), "H3K27ac");
        assert_eq!(display_label("input.cram"), "input");
        assert_eq!(display_label("sample.bw"), "sample");
        assert_eq!(display_label("plain"), "plain");
        assert_eq!(display_label(".bam"), ".bam");
    }

    #[test]
    fn thousands_separators() {
        assert_eq!(fmt_thousands(0.0), "0");
        assert_eq!(fmt_thousands(999.4), "999");
        assert_eq!(fmt_thousands(2232.3), "2,232");
        assert_eq!(fmt_thousands(1_234_567.0), "1,234,567");
        assert_eq!(fmt_thousands(-12_500.0), "-12,500");
    }
}
