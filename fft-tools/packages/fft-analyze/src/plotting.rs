use anyhow::{bail, Context, Result};
use charming::{
    component::{AngleAxis, Axis, Legend, PolarCoordinate, RadiusAxis, Title},
    datatype::CompositeValue,
    element::{AxisType, Color, CoordinateSystem, Tooltip, Trigger, TriggerOn},
    series::{Bar, Line},
    Chart, HtmlRenderer,
};
use ndarray::Array1;
use std::path::Path;

/// Generates and saves a rose plot for angular mean.
pub fn generate_and_save_rose_chart(
    path: &Path,
    title_text: &str,
    power_vs_angle: &Array1<f64>,
    angles: &Array1<f64>,
    header: &str,
) -> Result<()> {
    if angles.len() != 36 {
        bail!(
            "Rose chart generation requires 36 5-degree bins from 0-180 degrees, but received {} bins.",
            angles.len()
        );
    }

    // Step 1: Duplicate 0-180 data to 0-360
    let full_power: Vec<f64> = power_vs_angle
        .iter()
        .chain(power_vs_angle.iter())
        .cloned()
        .collect();

    // Step 2: Combine 72 5-degree bins into 36 10-degree bins
    let mut p_10deg = Vec::with_capacity(36);
    if !full_power.is_empty() {
        // Handle the wrap-around case for the 0-degree bin
        p_10deg.push((full_power.last().unwrap() + full_power[0]) / 2.0);
        // Handle the rest
        for i in (1..full_power.len().saturating_sub(1)).step_by(2) {
            p_10deg.push((full_power[i] + full_power[i + 1]) / 2.0);
        }
    }

    let subtitle = header
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| line.trim_start_matches('#').trim())
        .collect::<Vec<&str>>()
        .join("\n");
    let full_title = format!("{}\n{}", title_text, subtitle);

    let angle_labels: Vec<_> = (0..36)
        .map(|i| {
            let angle = i * 10;
            match angle {
                0 => "N".to_string(),
                90 => "E".to_string(),
                180 => "S".to_string(),
                270 => "W".to_string(),
                _ => format!("{}", angle),
            }
        })
        .collect();

    let chart = Chart::new()
        .title(Title::new().text(full_title).left("center"))
        .legend(Legend::new().show(false))
        .polar(PolarCoordinate::new().center(vec![
            CompositeValue::from("50%"),
            CompositeValue::from("50%"),
        ]))
        .angle_axis(
            AngleAxis::new()
                .start_angle(95)
                .type_(AxisType::Category)
                .data(angle_labels),
        )
        .radius_axis(RadiusAxis::new())
        .series(
            Bar::new()
                .coordinate_system(CoordinateSystem::Polar)
                .data(p_10deg)
                .item_style(charming::element::ItemStyle::new().color(Color::from("#2b6cb0"))),
        )
        .tooltip(
            Tooltip::new()
                .trigger(Trigger::Axis)
                .trigger_on(TriggerOn::Mousemove),
        );

    let mut renderer = HtmlRenderer::new(title_text, 1024, 768);
    renderer
        .save(&chart, path)
        .context("Failed to save rose chart to file")?;

    Ok(())
}

/// Generates and saves a line chart plot of the analysis results using charming.
/// NOTE: Due to limitations in charming 0.6.0, a true log-log plot is not possible.
/// This function creates a plot with a categorical X-axis and a logarithmic Y-axis.
pub fn generate_and_save_line_chart(
    path: &Path,
    title_text: &str,
    labels: &Array1<f64>, // Use labels for the X-axis
    values: &Array1<f64>,
    header: &str,
) -> Result<()> {
    let subtitle = header
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| line.trim_start_matches('#').trim())
        .collect::<Vec<&str>>()
        .join("\n");
    let full_title = format!("{}\n{}", title_text, subtitle);

    // Since we are using a category axis, we need to provide the labels as strings.
    // We also need to filter the values and labels together to keep them in sync.
    let (mut filtered_labels, mut data): (Vec<String>, Vec<f64>) = labels // Added mut
        .iter()
        .zip(values.iter())
        .filter(|(_, &v)| v.is_finite() && v > 0.0)
        .map(|(&lbl, &val)| (format!("{:.2}", lbl), val)) // Format labels to string
        .unzip();

    // Reverse the order of labels and data to ensure rising logical order on X-axis
    filtered_labels.reverse();
    data.reverse();

    let chart = Chart::new()
        .title(Title::new().text(full_title).left("center"))
        .legend(Legend::new().show(false))
        .x_axis(
            Axis::new()
                .type_(AxisType::Category) // Use Category axis
                .name("Wavelength (m)")
                .data(filtered_labels), // Pass labels as data
        )
        .y_axis(
            Axis::new()
                .type_(AxisType::Log)
                .name("Power")
                .split_number(10), // Added split_number for more grid lines
        )
        .series(
            Line::new()
                .data(data)
                .item_style(charming::element::ItemStyle::new().color(Color::from("#c53030"))),
        )
        .tooltip(
            Tooltip::new()
                .trigger(Trigger::Axis)
                .trigger_on(TriggerOn::Mousemove),
        );

    let mut renderer = HtmlRenderer::new(title_text, 1024, 768);
    renderer
        .save(&chart, path)
        .context("Failed to save line chart to file")?;

    Ok(())
}

/// Generates and saves a line chart plot for residuals vs. wavelength.
pub fn generate_and_save_residuals_chart(
    path: &Path,
    title_text: &str,
    labels: &Array1<f64>, // Wavelengths
    residuals: &Array1<f64>,
    header: &str,
) -> Result<()> {
    let subtitle = header
        .lines()
        .filter(|line| line.starts_with('#'))
        .map(|line| line.trim_start_matches('#').trim())
        .collect::<Vec<&str>>()
        .join("\n");
    let full_title = format!("{}\n{}", title_text, subtitle);

    // Filter residuals for finite values. Wavelengths are filtered in main.rs
    let (mut filtered_labels, mut filtered_residuals): (Vec<String>, Vec<f64>) = labels
        .iter()
        .zip(residuals.iter())
        .filter(|(_, &r)| r.is_finite())
        .map(|(&lbl, &res)| (format!("{:.2}", lbl), res))
        .unzip();

    // Reverse the order of labels and data to ensure rising logical order on X-axis
    filtered_labels.reverse();
    filtered_residuals.reverse();

    let chart = Chart::new()
        .title(Title::new().text(full_title).left("center"))
        .legend(Legend::new().show(false))
        .x_axis(
            Axis::new()
                .type_(AxisType::Category)
                .name("Wavelength (m)")
                .data(filtered_labels),
        )
        .y_axis(
            Axis::new()
                .type_(AxisType::Value) // Linear axis for residuals
                .name("Relative Residual (%)")
                .split_number(10),
        )
        .series(
            Line::new()
                .data(filtered_residuals)
                .item_style(charming::element::ItemStyle::new().color(Color::from("#c53030"))),
        )
        .tooltip(
            Tooltip::new()
                .trigger(Trigger::Axis)
                .trigger_on(TriggerOn::Mousemove),
        );

    let mut renderer = HtmlRenderer::new(title_text, 1024, 768);
    renderer
        .save(&chart, path)
        .context("Failed to save residuals chart to file")?;

    Ok(())
}
