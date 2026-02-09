use anyhow::{Context, Result};
use gdal::Dataset;
use ndarray::{Array1, Array2, Axis};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use crate::cli::AnalysisMode;

/// Represents the metadata extracted from a polar FFT block.
pub struct PolarMetadata {
    pub wavelengths: Array1<f64>,
    pub angles: Array1<f64>,
}

/// Reads and parses the metadata JSON file associated with a polar FFT block.
pub fn read_polar_metadata(polar_path: &Path) -> Result<PolarMetadata> {
    let meta_path = polar_path
        .with_file_name(
            polar_path
                .file_name()
                .context("Invalid polar_path")?
                .to_str()
                .context("Invalid UTF-8 in polar_path")?
                .replace("fft_polar_block_", "fft_polar_metadata_block_"),
        )
        .with_extension("json");

    let meta_file = fs::File::open(&meta_path)
        .context(format!("Failed to open metadata file: {:?}", meta_path))?;
    let metadata: BTreeMap<String, serde_json::Value> = serde_json::from_reader(meta_file)?;

    let wavelengths: Vec<f64> = serde_json::from_value(
        metadata
            .get("wavelengths")
            .context("Missing 'wavelengths' in metadata")?
            .clone(),
    )?;
    let angles: Vec<f64> = serde_json::from_value(
        metadata
            .get("angles")
            .context("Missing 'angles' in metadata")?
            .clone(),
    )?;

    Ok(PolarMetadata {
        wavelengths: Array1::from(wavelengths),
        angles: Array1::from(angles),
    })
}

/// Finds the index of the value in `values` that is nearest to `target`.
///
/// # Arguments
/// * `values` - An `Array1<f64>` of values to search within.
/// * `target` - The target `f64` value.
///
/// # Returns
/// The `usize` index of the nearest value.
pub fn find_nearest_index(values: &Array1<f64>, target: f64) -> usize {
    values
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            (*a - target)
                .abs()
                .partial_cmp(&(*b - target).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(index, _)| index)
        .unwrap_or(0)
}

/// Parses a comma-separated string of two numbers into an optional tuple of f64.
///
/// # Arguments
/// * `bounds_str` - An `Option<String>` containing the bounds string (e.g., "10,20").
///
/// # Returns
/// A `Result` containing an `Option<(f64, f64)>` representing the parsed min and max bounds.
pub fn parse_bounds(bounds_str: Option<String>) -> Result<Option<(f64, f64)>> {
    match bounds_str {
        Some(s) => {
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() != 2 {
                anyhow::bail!("Bounds must be two comma-separated numbers.");
            }
            let min = parts[0]
                .parse::<f64>()
                .context("Failed to parse min bound.")?;
            let max = parts[1]
                .parse::<f64>()
                .context("Failed to parse max bound.")?;
            Ok(Some((min, max)))
        }
        None => Ok(None),
    }
}

/// Processes a single polar spectrum file, applies bounds, and calculates a NaN-aware mean profile.
///
/// # Arguments
/// * `polar_path` - Path to the polar spectrum TIFF file.
/// * `mode` - The analysis mode (RadialMean or AngularMean).
/// * `wavelength_bounds` - Optional tuple (min, max) for wavelength filtering.
/// * `angle_bounds` - Optional tuple (min, max) for angle filtering.
///
/// # Returns
/// A `Result` containing an `Array1<f64>` representing the mean power profile for the block.
pub fn process_file(
    polar_path: &Path,
    mode: AnalysisMode,
    wavelength_bounds: Option<(f64, f64)>,
    angle_bounds: Option<(f64, f64)>,
) -> Result<Array1<f64>> {
    let dataset = Dataset::open(polar_path)?;
    let band = dataset.rasterband(1)?;
    let (cols, rows) = band.size();
    let nodata_value = band.no_data_value();

    // Read log-transformed power spectrum.
    let power_log = band
        .read_as::<f64>((0, 0), (cols, rows), (cols, rows), None)?
        .data()
        .to_vec();

    // Convert to linear scale, explicitly handling nodata values.
    let power_spectrum = Array2::from_shape_vec((rows, cols), power_log)?.mapv(|v| {
        let is_nodata = match nodata_value {
            Some(ndv) if ndv.is_nan() => v.is_nan(),
            Some(ndv) => v == ndv,
            None => v.is_nan(), // If no nodata is set, still treat NaNs as nodata
        };

        if is_nodata {
            f64::NAN
        } else {
            10.0_f64.powf(v)
        }
    });

    let polar_metadata = read_polar_metadata(polar_path)?;
    let wavelengths = polar_metadata.wavelengths;
    let angles_arr = polar_metadata.angles;

    // --- Subsetting Logic ---
    // Determine which wavelength indices to include based on bounds.
    let mut w_indices: Vec<usize> = (0..wavelengths.len()).collect();
    if let Some((min_w, max_w)) = wavelength_bounds {
        let mut i1 = find_nearest_index(&wavelengths, min_w);
        let mut i2 = find_nearest_index(&wavelengths, max_w);
        if i1 > i2 {
            std::mem::swap(&mut i1, &mut i2);
        }
        w_indices = (i1..=i2).collect();
    }
    // Determine which angle indices to include based on bounds, handling circular wrap-around.
    let mut a_indices: Vec<usize> = (0..angles_arr.len()).collect();
    if let Some((min_a, max_a)) = angle_bounds {
        let min_idx = find_nearest_index(&angles_arr, min_a);
        let max_idx = find_nearest_index(&angles_arr, max_a);
        if min_idx <= max_idx {
            a_indices = (min_idx..=max_idx).collect();
        } else {
            // Circular selection: from min_idx to end, and from start to max_idx.
            a_indices = (min_idx..angles_arr.len()).chain(0..=max_idx).collect();
        }
    }

    // Select the subset of the power spectrum based on the determined indices.
    let subset_power = power_spectrum
        .select(Axis(0), &w_indices)
        .select(Axis(1), &a_indices);

    // Perform NaN-aware mean calculation based on the analysis mode.
    let mean_profile = match mode {
        AnalysisMode::RadialMean => {
            let mut means = Array1::<f64>::zeros(subset_power.nrows());
            for (i, row) in subset_power.outer_iter().enumerate() {
                let valid_values: Vec<f64> = row.iter().filter(|v| !v.is_nan()).cloned().collect();
                if !valid_values.is_empty() {
                    means[i] = valid_values.iter().sum::<f64>() / valid_values.len() as f64;
                } else {
                    means[i] = f64::NAN;
                }
            }
            means
        }
        AnalysisMode::AngularMean => {
            let mut means = Array1::<f64>::zeros(subset_power.ncols());
            for (i, col) in subset_power.axis_iter(Axis(1)).enumerate() {
                let valid_values: Vec<f64> = col.iter().filter(|v| !v.is_nan()).cloned().collect();
                if !valid_values.is_empty() {
                    means[i] = valid_values.iter().sum::<f64>() / valid_values.len() as f64;
                } else {
                    means[i] = f64::NAN;
                }
            }
            means
        }
    };

    Ok(mean_profile)
}
