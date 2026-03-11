use assert_cmd::Command;
use gdal::raster::Buffer;
use gdal::DriverManager;
use gdal::spatial_ref::SpatialRef;
use std::path::Path;

/// Creates a synthetic DEM (GeoTIFF) containing a superposition of sine and cosine waves.
/// 
/// This helper generates a deterministic terrain surface for testing spatial decomposition
/// and FFT accuracy. It includes a standard GeoTransform and WKT projection.
fn create_synthetic_dem(path: &Path, width: usize, height: usize) {
    let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
    let mut ds = driver.create_with_band_type::<f64, _>(path, width, height, 1).unwrap();
    
    // Set GeoTransform: origin x, pixel width, rotation, origin y, rotation, pixel height.
    let geo_transform = [0.0, 1.0, 0.0, 0.0, 0.0, -1.0];
    ds.set_geo_transform(&geo_transform).unwrap();

    let srs = SpatialRef::from_epsg(4326).unwrap();
    ds.set_spatial_ref(&srs).unwrap();
    
    // Generate wave pattern: z(x,y) = sin(x/10) + cos(y/10).
    let mut data = vec![0.0; width * height];
    for y in 0..height {
        for x in 0..width {
            data[y * width + x] = (x as f64 / 10.0).sin() + (y as f64 / 10.0).cos();
        }
    }
    
    let mut band = ds.rasterband(1).unwrap();
    let mut buffer = Buffer::new((width, height), data);
    band.write((0, 0), (width, height), &mut buffer).unwrap();
    band.set_no_data_value(Some(-9999.0)).unwrap();
}

/// Verifies the full execution lifecycle of the `fft-process` tool.
/// 
/// This integration test ensures that:
/// 1. The sliding window decomposition correctly calculates block coordinates.
/// 2. All required output artifacts (PSD, metadata) are generated for each block.
/// 3. Scientific validity (e.g., energy conservation) is maintained across the transform.
#[test]
fn test_fft_process_execution() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_path = temp_dir.path().join("input.tif");
    let output_dir = temp_dir.path().join("output");

    // Create a 64x64 synthetic DEM.
    create_synthetic_dem(&input_path, 64, 64);

    // Execute `fft-process` with a 32x32 window and 0 overlap.
    // This should result in exactly 4 non-overlapping blocks.
    let mut cmd = Command::new(assert_cmd::cargo_bin!("fft-process"));
    cmd.arg("--input").arg(&input_path)
       .arg("--output").arg(&output_dir)
       .arg("--window-size").arg("32")
       .arg("--overlap").arg("0")
       .arg("--detrend").arg("1")
       .assert()
       .success();

    // Verify the existence of output files for each block.
    let expected_files = [
        "fft_psd_block_0_0.tif",
        "fft_metadata_block_0_0.json",
        "fft_psd_block_0_32.tif",
        "fft_psd_block_32_0.tif",
        "fft_psd_block_32_32.tif",
    ];

    for file in expected_files {
        assert!(output_dir.join(file).exists(), "Expected artifact {} missing from output directory.", file);
    }
    
    // Perform a scientific validity check on the generated metadata.
    let meta_path = output_dir.join("fft_metadata_block_0_0.json");
    let meta_content = std::fs::read_to_string(meta_path).unwrap();
    let meta_json: serde_json::Value = serde_json::from_str(&meta_content).unwrap();
    
    // 1. Verify Normalization Metadata
    let norm = meta_json.get("normalization").expect("Normalization object missing from metadata.");
    assert_eq!(norm["psd_normalization"], "physical_2D");
    assert_eq!(norm["pixel_size_x"], 1.0);
    assert_eq!(norm["pixel_size_y"], 1.0);
    assert_eq!(norm["original_block_size"], serde_json::json!([32, 32]));

    // 2. Parseval's check: Ensure energy conservation between spatial and frequency domains.
    if let Some(stats) = meta_json.get("statistics") {
        if let Some(error_val) = stats.get("parseval_error") {
            let parseval_error = error_val.as_f64().unwrap();
            // Tolerance is slightly higher than core check to account for multi-tool I/O float precision.
            assert!(parseval_error < 1e-5, "Parseval energy conservation check failed (error: {:.2e}).", parseval_error);
        } else {
             panic!("'parseval_error' field missing from block statistics.");
        }
    } else {
        panic!("'statistics' object missing from block metadata JSON.");
    }
}
