use assert_cmd::Command;
use gdal::raster::Buffer;
use gdal::DriverManager;
use std::path::Path;
use serde_json::json;
use approx::assert_relative_eq;

/// Creates a mock Cartesian PSD GeoTIFF and metadata for testing.
/// 
/// This helper simulates the output of `fft-process`, allowing for verification 
/// of the polar transformation logic in `fft-polar`.
fn create_mock_psd(dir: &Path, row: usize, col: usize, size: usize, pixel_size: f64) {
    let psd_path = dir.join(format!("fft_psd_block_{}_{}.tif", row, col));
    let meta_path = dir.join(format!("fft_metadata_block_{}_{}.json", row, col));

    let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
    let ds = driver.create_with_band_type::<f64, _>(&psd_path, size, size, 1).unwrap();

    // Fill with a uniform linear power spectral density (P = 1.0).
    // In `fft-process`, the values are stored as log10(P + 1e-12).
    let p_linear = 1.0;
    let p_log = (p_linear + 1e-12_f64).log10();
    let data = vec![p_log; size * size];

    let mut band = ds.rasterband(1).unwrap();
    let mut buffer = Buffer::new((size, size), data);
    band.write((0, 0), (size, size), &mut buffer).unwrap();

    let metadata = json!({
        "block_position": [row, col],
        "original_size": [size, size],
        "padded_size": [size, size],
        "processing_params": {},
        "normalization": {
            "original_block_size": [size, size],
            "pixel_size_x": pixel_size,
            "pixel_size_y": pixel_size,
            "psd_normalization": "physical_2D"
        },
        "statistics": {
            "f_nyquist": 1.0 / (2.0 * pixel_size)
        }
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();
}

/// Verifies conservation of power during the Cartesian-to-Polar transformation.
/// 
/// Scientifically, the total topographic variance (power) must be preserved 
/// when re-binning the spectrum from a rectangular grid to polar coordinates.
/// This test integrates the power density across the polar bins and compares 
/// it to the analytical sum of the Cartesian grid.
#[test]
fn test_polar_energy_conservation_integration() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_dir = temp_dir.path().join("input");
    let output_dir = temp_dir.path().join("output");
    std::fs::create_dir_all(&input_dir).unwrap();

    let size = 64;
    let pixel_size = 1.0;
    create_mock_psd(&input_dir, 0, 0, size, pixel_size);

    // Execute the `fft-polar` CLI tool.
    let mut cmd = Command::new(assert_cmd::cargo_bin!("fft-polar"));
    cmd.arg("--input").arg(&input_dir)
       .arg("--output").arg(&output_dir)
       .arg("--n-angles").arg("36")
       .arg("--n-wavenumbers").arg("64")
       .assert()
       .success();

    // Load the generated polar PSD and metadata.
    // fft-polar will create output_dir/fft_polar_block_0_0.tif
    let polar_tif = output_dir.join("fft_polar_block_0_0.tif");
    let ds = gdal::Dataset::open(polar_tif).unwrap();
    let band = ds.rasterband(1).unwrap();
    let (cols, rows) = band.size();
    let data_log = band.read_as::<f64>((0, 0), (cols, rows), (cols, rows), None).unwrap().data().to_vec();
    
    // Convert power values back to the linear domain.
    let polar_psd: Vec<f64> = data_log.iter().map(|&p| if p.is_nan() { 0.0 } else { 10.0_f64.powf(p) }).collect();

    // Metadata filename follows replace("fft_psd", "fft_polar_metadata") logic.
    let meta_path = output_dir.join("fft_polar_metadata_block_0_0.json");
    let meta_content = std::fs::read_to_string(meta_path).unwrap();
    let meta_json: serde_json::Value = serde_json::from_str(&meta_content).unwrap();
    
    let n_angles = meta_json["n_angles"].as_u64().unwrap() as usize;
    let n_wavenumbers = meta_json["n_wavelengths"].as_u64().unwrap() as usize;

    // Calculate total power in polar space by integrating PSD over bin areas.
    // Bin Area = 0.5 * (k_outer^2 - k_inner^2) * Δθ.
    
    let k_min = 2.0 / (size as f64 * pixel_size);
    let k_max = 1.0 / (2.0 * pixel_size);
    
    // Replicate the logarithmic wavenumber spacing used in the tool.
    let log_k_min = k_min.log10();
    let log_k_max = k_max.log10();
    let mut k_edges = Vec::new();
    for i in 0..=n_wavenumbers {
        let t = i as f64 / n_wavenumbers as f64;
        k_edges.push(10.0_f64.powf(log_k_min + t * (log_k_max - log_k_min)));
    }
    
    let d_theta = std::f64::consts::PI / n_angles as f64;
    
    let mut total_polar_power = 0.0;
    for kw in 0..n_wavenumbers {
        let k_inner = k_edges[kw];
        let k_outer = k_edges[kw+1];
        for a in 0..n_angles {
            let p = polar_psd[kw * n_angles + a];
            let bin_area = 0.5 * (k_outer.powi(2) - k_inner.powi(2)) * d_theta;
            total_polar_power += p * bin_area;
        }
    }

    // Calculate the expected total power from the original Cartesian grid.
    // Only frequencies within the [k_min, k_max] range are included.
    let dkx = 1.0 / (size as f64 * pixel_size);
    let mut expected_power = 0.0;
    for y in 0..size {
        for x in 0..size {
            let fx = (x as f64 - size as f64 / 2.0) * dkx;
            let fy = (y as f64 - size as f64 / 2.0) * dkx;
            let k = (fx.powi(2) + fy.powi(2)).sqrt();
            if k >= k_min && k <= k_max { 
                expected_power += 1.0 * dkx * dkx;
            }
        }
    }

    // Verify that total power is conserved within a small tolerance.
    // Minor deviations are expected due to discrete binning and interpolation.
    assert_relative_eq!(total_polar_power, expected_power, epsilon = 1e-3);
}
