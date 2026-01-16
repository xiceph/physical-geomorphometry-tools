use gdal::{Dataset, DriverManager};
use gdal::raster::Buffer;
use ndarray::Array1;
use rayon::prelude::*;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub fn normalize_raster(
    input: PathBuf,
    output: PathBuf,
    csv: Option<PathBuf>,
    initial_k: f64,
    tolerance: f64,
) -> gdal::errors::Result<()> {
    let dataset = Dataset::open(&input)?;
    let band = dataset.rasterband(1)?;
    let (width, height) = band.size();

    // Use the raster's NoData value if it exists, otherwise fallback to a default.
    let no_data_value = band.no_data_value().unwrap_or(-9999.0);

    let buffer: Buffer<f64> = band.read_as((0, 0), (width, height), (width, height), None)?;
    let raster_data_slice = buffer.data();

    // Separate valid data from NoData values, keeping original indices for reconstruction.
    let mut valid_values = Vec::with_capacity(raster_data_slice.len());
    let mut valid_indices = Vec::with_capacity(raster_data_slice.len());

    for (i, &val) in raster_data_slice.iter().enumerate() {
        if (val - no_data_value).abs() >= 1e-9 {
            valid_values.push(val);
            valid_indices.push(i);
        }
    }

    if valid_values.is_empty() {
        panic!("No valid data found in raster.");
    }

    let data_for_optimization = Array1::from(valid_values);

    let (optimal_k, final_kurt, k_list) =
        optimize_kurtosis(data_for_optimization.view(), initial_k, tolerance);

    println!(
        "Optimal k = {:.5}, Final kurtosis = {:.5}",
        optimal_k, final_kurt
    );

    // Apply the arctan transformation to the valid data using the optimal k.
    let transformed_valid_data: Vec<f64> = data_for_optimization
        .into_par_iter()
        .map(|val| (val * optimal_k).atan())
        .collect();

    // Reconstruct the full raster data by placing transformed values at their original indices.
    let mut transformed_full_raster: Vec<f64> = vec![no_data_value; raster_data_slice.len()];

    // To perform a safe "scatter" write, we pair indices with their new values,
    // sort by index, and then perform a fast sequential write. This avoids
    // complex parallel write synchronization and is highly efficient.
    let mut indexed_data: Vec<(usize, f64)> = valid_indices
        .into_par_iter()
        .zip(transformed_valid_data.into_par_iter())
        .collect();

    indexed_data.par_sort_unstable_by_key(|(idx, _)| *idx);

    indexed_data.into_iter().for_each(|(idx, val)| transformed_full_raster[idx] = val);

    let driver = DriverManager::get_driver_by_name("GTiff")?;
    let mut out_ds = driver.create_with_band_type::<f64, _>(&output, width, height, 1)?;
    out_ds.set_projection(&dataset.projection())?;
    out_ds.set_geo_transform(&dataset.geo_transform()?)?;

    let mut out_band = out_ds.rasterband(1)?;
    out_band.set_no_data_value(Some(no_data_value))?;
    let mut out_buffer = Buffer::new((width, height), transformed_full_raster);
    out_band.write((0, 0), (width, height), &mut out_buffer)?;

    if let Some(csv_path) = csv {
        let mut f = File::create(csv_path).unwrap();
        writeln!(f, "k,Kurtosis").unwrap();
        for (k, kurt) in k_list {
            writeln!(f, "{:.6},{:.6}", k, kurt).unwrap();
        }
    }

    Ok(())
}

/// Finds the optimal `k` value by iteratively searching for the minimum positive kurtosis.
fn optimize_kurtosis(
    data: ndarray::ArrayView1<f64>,
    initial_k: f64,
    tol: f64,
) -> (f64, f64, Vec<(f64, f64)>) {
    let mut k = initial_k;
    let mut step = 10000.0; 
    let max_k = 1e7;
    let mut history = Vec::new();

    // Store the kurtosis of the original, untransformed data (k=0) as a baseline.
    let kurt_at_zero = calculate_kurtosis_untransformed(data);
    history.push((0.0, kurt_at_zero));

    // If initial_k is 0, start the search at the first step to avoid a pointless calculation.
    let current_search_k = if k == 0.0 { step } else { k };
    let mut current_kurt = calculate_kurtosis_transformed(data, current_search_k);
    history.push((current_search_k, current_kurt));

    k = current_search_k;

    // Outer loop reduces the step size to refine the search for `k`.
    while step > tol {
        let mut temp_k = k;
        let mut temp_kurt = current_kurt;

        // Inner loop performs a linear search for a better `k` at the current step size.
        loop {
            let next_k = temp_k + step;
            if next_k > max_k {
                break;
            }

            let next_kurt = calculate_kurtosis_transformed(data, next_k);
            history.push((next_k, next_kurt));

            // Stop searching if kurtosis becomes non-positive or starts increasing,
            // as we've likely found or passed the minimum for this step.
            if next_kurt <= 0.0 || next_kurt > temp_kurt {
                break;
            }

            // `next_k` is better, so update `temp_k` and `temp_kurt`
            temp_k = next_k;
            temp_kurt = next_kurt;
        }

        // Update the best `k` found and reduce the step size for the next refinement pass.
        k = temp_k;
        current_kurt = temp_kurt;
        step /= 10.0;
    }

    (k, current_kurt, history)
}

/// Calculates kurtosis for the original, untransformed data.
fn calculate_kurtosis_untransformed(data: ndarray::ArrayView1<f64>) -> f64 {
    let n = data.len() as f64;
    if n < 4.0 {
        return f64::NAN;
    }

    let mean = data.mean().unwrap_or(0.0);

    let (m2_sum, m4_sum) = data.into_par_iter().fold(|| (0.0, 0.0), |(m2_acc, m4_acc), &x| {
        let diff = x - mean;
        (m2_acc + diff * diff, m4_acc + diff.powi(4))
    }).reduce(|| (0.0, 0.0), |(m2_a, m4_a), (m2_b, m4_b)| (m2_a + m2_b, m4_a + m4_b));

    let m2 = m2_sum / n;
    let m4 = m4_sum / n;

    if m2 == 0.0 {
        return f64::NAN;
    }

    (m4 / (m2 * m2)) - 3.0
}

/// Calculates kurtosis for data after applying the `(k * x).atan()` transformation.
/// This is calculated on-the-fly to avoid allocating a new array for each `k` value.
fn calculate_kurtosis_transformed(data: ndarray::ArrayView1<f64>, k: f64) -> f64 {
    let n = data.len() as f64;
    if n < 4.0 {
        return f64::NAN;
    }

    // First pass: calculate mean of transformed data
    let transformed_mean = data.into_par_iter()
        .map(|&x| (k * x).atan())
        .sum::<f64>() / n;

    // Second pass: calculate m2 and m4 of transformed data
    let (m2_sum, m4_sum) = data.into_par_iter()
        .fold(|| (0.0, 0.0), |(m2_acc, m4_acc), &x| {
            let transformed_val = (k * x).atan();
            let diff = transformed_val - transformed_mean;
            (m2_acc + diff * diff, m4_acc + diff.powi(4))
        })
        .reduce(|| (0.0, 0.0), |(m2_a, m4_a), (m2_b, m4_b)| (m2_a + m2_b, m4_a + m4_b));

    let m2 = m2_sum / n;
    let m4 = m4_sum / n;

    if m2 == 0.0 {
        return f64::NAN;
    }

    (m4 / (m2 * m2)) - 3.0
}
