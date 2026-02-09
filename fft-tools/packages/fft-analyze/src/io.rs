use anyhow::Result;
use ndarray::{Array1, Axis};
use std::io::Write;
use std::path::Path;

use crate::cli::AnalysisMode;
use crate::processing::{find_nearest_index, read_polar_metadata};

/// Retrieves axis labels and constructs a header string for CSV/plot, applying bounds.
///
/// # Arguments
/// * `polar_path` - Path to a polar spectrum TIFF file (used to load metadata).
/// * `mode` - The analysis mode.
/// * `wavelength_bounds` - Optional tuple (min, max) for wavelength filtering.
/// * `angle_bounds` - Optional tuple (min, max) for angle filtering.
///
/// # Returns
/// A `Result` containing a tuple of `(Array1<f64>, String)`:
///   - `Array1<f64>`: The subsetted axis labels (wavelengths or angles).
///   - `String`: The header for CSV/plot, including bounds information as comments.
pub fn get_axis_info(
    polar_path: &Path,
    mode: &AnalysisMode,
    wavelength_bounds: Option<(f64, f64)>,
    angle_bounds: Option<(f64, f64)>,
    detrend_order: Option<usize>,
) -> Result<(Array1<f64>, String)> {
    // Load metadata from the polar spectrum file.
    let polar_metadata = read_polar_metadata(polar_path)?;
    let all_wavelengths_arr = polar_metadata.wavelengths;
    let all_angles_arr = polar_metadata.angles;

    // Generate bounds information strings for the header.
    let mut w_bounds_info = "".to_string();
    if let Some((min_w, max_w)) = wavelength_bounds {
        let min_idx = find_nearest_index(&all_wavelengths_arr, min_w);
        let max_idx = find_nearest_index(&all_wavelengths_arr, max_w);
        w_bounds_info = format!(
            "# Wavelength bounds: {:.2}m to {:.2}m",
            all_wavelengths_arr[min_idx], all_wavelengths_arr[max_idx]
        );
    }

    let mut a_bounds_info = "".to_string();
    if let Some((min_a, max_a)) = angle_bounds {
        let min_idx = find_nearest_index(&all_angles_arr, min_a);
        let max_idx = find_nearest_index(&all_angles_arr, max_a);
        a_bounds_info = format!(
            "# Angle bounds: {:.1}° to {:.1}°",
            all_angles_arr[min_idx], all_angles_arr[max_idx]
        );
    }

    let mut detrend_info = "".to_string();
    if let Some(order) = detrend_order {
        detrend_info = format!("# Detrending of order: {}", order);
    }

    // Combine bounds information into a single string for the header.
    let mut comments = vec![];
    if !w_bounds_info.is_empty() {
        comments.push(w_bounds_info);
    }
    if !a_bounds_info.is_empty() {
        comments.push(a_bounds_info);
    }
    if !detrend_info.is_empty() {
        comments.push(detrend_info);
    }
    let comments_str = comments.join("\n");

    match mode {
        AnalysisMode::RadialMean => {
            // Determine wavelength indices to include.
            let mut w_indices: Vec<usize> = (0..all_wavelengths_arr.len()).collect();
            if let Some((min_w, max_w)) = wavelength_bounds {
                let min_idx = find_nearest_index(&all_wavelengths_arr, min_w);
                let max_idx = find_nearest_index(&all_wavelengths_arr, max_w);
                w_indices = (min_idx..=max_idx).collect();
            }
            let subset_wavelengths = all_wavelengths_arr.select(Axis(0), &w_indices);

            let mut header_cols = vec!["wavelength", "average_psd"];
            if detrend_order.is_some() {
                header_cols.push("trend_power");
                header_cols.push("relative_residual_pct");
            }

            let header = format!("{}\n{}", comments_str, header_cols.join(","));
            Ok((subset_wavelengths, header))
        }
        AnalysisMode::AngularMean => {
            // Determine angle indices to include, handling circular wrap-around.
            let mut a_indices: Vec<usize> = (0..all_angles_arr.len()).collect();
            if let Some((min_a, max_a)) = angle_bounds {
                let min_idx = find_nearest_index(&all_angles_arr, min_a);
                let max_idx = find_nearest_index(&all_angles_arr, max_a);
                if min_idx <= max_idx {
                    a_indices = (min_idx..=max_idx).collect();
                } else {
                    a_indices = (min_idx..all_angles_arr.len()).chain(0..=max_idx).collect();
                }
            }
            let subset_angles = all_angles_arr.select(Axis(0), &a_indices);
            let header = format!("{}\nangle,average_psd", comments_str);
            Ok((subset_angles, header))
        }
    }
}

/// Saves the analysis results to a CSV file.
///
/// # Arguments
/// * `path` - The path to save the CSV file.
/// * `header` - The header string, which may include comment lines starting with '#'.
/// * `labels` - An `Array1<f64>` containing the axis labels (wavelengths or angles).
/// * `values` - An `Array1<f64>` containing the corresponding average PSD values.
pub fn save_to_csv(
    path: &Path,
    header: &str,
    labels: &Array1<f64>,
    values: &Array1<f64>,
    trend: Option<&Array1<f64>>,
    residuals: Option<&Array1<f64>>,
) -> Result<()> {
    let mut file = std::fs::File::create(path)?;

    // Write comment lines (starting with '#') directly to the file before CSV writer.
    let header_lines: Vec<&str> = header.lines().collect();
    for line in &header_lines {
        if line.starts_with('#') {
            writeln!(file, "{}", line)?;
        }
    }

    // Initialize CSV writer and write the actual header row.
    let mut wtr = csv::Writer::from_writer(file);
    if let Some(last_line) = header_lines.last() {
        if !last_line.starts_with('#') {
            wtr.write_record(last_line.split(','))?;
        }
    }

    // Write data rows.
    for i in 0..labels.len() {
        let mut record = vec![labels[i].to_string(), values[i].to_string()];
        if let Some(trend_vals) = trend {
            record.push(trend_vals[i].to_string());
        }
        if let Some(residual_vals) = residuals {
            record.push(residual_vals[i].to_string());
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    Ok(())
}
