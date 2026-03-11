use assert_cmd::Command;
use gdal::raster::Buffer;
use gdal::DriverManager;
use std::path::Path;
use serde_json::json;
use approx::assert_relative_eq;

/// Creates a mock polar PSD GeoTIFF and associated metadata for testing.
/// 
/// This helper simulates the output of `fft-polar`, allowing for controlled 
/// testing of the statistical aggregation logic in `fft-analyze`.
fn create_mock_polar(dir: &Path, block_id: &str, data: Vec<f64>, rows: usize, cols: usize) {
    let polar_path = dir.join(format!("fft_polar_block_{}.tif", block_id));
    let meta_path = dir.join(format!("fft_polar_metadata_block_{}.json", block_id));

    let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
    let ds = driver.create_with_band_type::<f64, _>(&polar_path, cols, rows, 1).unwrap();
    
    let mut band = ds.rasterband(1).unwrap();
    let mut buffer = Buffer::new((cols, rows), data);
    band.write((0, 0), (cols, rows), &mut buffer).unwrap();
    band.set_no_data_value(Some(f64::NAN)).unwrap();

    let wavelengths: Vec<f64> = (0..rows).map(|i| (i + 1) as f64 * 10.0).collect();
    let angles: Vec<f64> = (0..cols).map(|i| (i as f64 * 10.0)).collect();

    let metadata = json!({
        "wavelengths": wavelengths,
        "angles": angles,
        "n_wavelengths": rows,
        "n_angles": cols
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();
}

/// Verifies that `fft-analyze` correctly performs NaN-aware averaging across blocks.
/// 
/// In real-world DEM analysis, blocks near the data boundaries or containing no-data
/// values might result in NaNs in the polar representation. Scientific analysis 
/// requires these to be handled gracefully (ignored in the mean) rather than 
/// propagating the NaN to the final result.
#[test]
fn test_analyze_averaging_with_nans() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_dir = temp_dir.path().join("input");
    let output_file = temp_dir.path().join("summary.csv");
    std::fs::create_dir(&input_dir).unwrap();

    // Define two mock blocks with 2 wavelengths and 1 angle each.
    // Block 1 contains a NaN value.
    // Block 2 contains valid data.
    // Note: Input data to `fft-analyze` is log-transformed.
    // Block 1: [log10(1.0)=0.0, NaN]
    // Block 2: [log10(2.0)≈0.301, log10(3.0)≈0.477]
    create_mock_polar(&input_dir, "0_0", vec![0.0, f64::NAN], 2, 1);
    create_mock_polar(&input_dir, "0_1", vec![0.30102999566, 0.47712125472], 2, 1);

    // Execute the `fft-analyze` CLI tool in radial-mean mode.
    let mut cmd = Command::new(assert_cmd::cargo_bin!("fft-analyze"));
    cmd.arg("--input").arg(&input_dir)
       .arg("--output").arg(&output_file)
       .arg("--mode").arg("radial-mean")
       .assert()
       .success();

    // Read the resulting CSV summary.
    let mut rdr = csv::Reader::from_path(output_file).unwrap();
    let records: Vec<_> = rdr.records().collect();
    
    // The output should contain 2 records (one for each wavelength bin).
    // The mean should be calculated in the linear domain, then optionally logged back.
    
    // Bin 1: Mean of 1.0 and 2.0 = 1.5.
    let val1: f64 = records[0].as_ref().unwrap()[1].parse().unwrap();
    assert_relative_eq!(val1, 1.5, epsilon = 1e-4);
    
    // Bin 2: Mean of NaN and 3.0 = 3.0 (NaN is ignored).
    let val2: f64 = records[1].as_ref().unwrap()[1].parse().unwrap();
    assert_relative_eq!(val2, 3.0, epsilon = 1e-4);
}
