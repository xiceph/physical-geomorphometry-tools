use assert_cmd::Command;
use gdal::raster::Buffer;
use gdal::DriverManager;
use gdal::spatial_ref::SpatialRef;
use std::path::Path;

fn create_synthetic_dem(path: &Path, width: usize, height: usize) {
    let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
    let mut ds = driver.create_with_band_type::<f64, _>(path, width, height, 1).unwrap();
    
    // Set GeoTransform and Projection (WGS84)
    // origin x, pixel width, rotation, origin y, rotation, pixel height
    let geo_transform = [0.0, 1.0, 0.0, 0.0, 0.0, -1.0];
    ds.set_geo_transform(&geo_transform).unwrap();

    let srs = SpatialRef::from_epsg(4326).unwrap();
    ds.set_spatial_ref(&srs).unwrap();
    
    // Simple sine wave
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

#[test]
fn test_fft_process_execution() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_path = temp_dir.path().join("input.tif");
    // Ensure output directory name is unique or clean.
    // fft-process appends .1, .2 if dir exists, but here we provide a path that shouldn't exist yet?
    // temp_dir exists. We want a subdirectory.
    let output_dir = temp_dir.path().join("output");

    create_synthetic_dem(&input_path, 64, 64);

    let mut cmd = Command::cargo_bin("fft-process").unwrap();
    cmd.arg("--input").arg(&input_path)
       .arg("--output").arg(&output_dir)
       .arg("--window-size").arg("32")
       .arg("--overlap").arg("0")
       .arg("--detrend").arg("1")
       .assert()
       .success();

    // Check outputs
    // 64x64 DEM, 32x32 window, 0 overlap -> 4 blocks: (0,0), (0,32), (32,0), (32,32)
    // Wait, the logic loop:
    // for r in 0..n_rows { for c in 0..n_cols { ... } }
    // n_cols = (width - window) / step + 1
    // (64 - 32) / 32 + 1 = 1 + 1 = 2.
    // So 0 and 32. Correct.

    let expected_files = [
        "fft_psd_block_0_0.tif",
        "fft_metadata_block_0_0.json",
        "fft_psd_block_0_32.tif",
        "fft_psd_block_32_0.tif",
        "fft_psd_block_32_32.tif",
    ];

    for file in expected_files {
        assert!(output_dir.join(file).exists(), "File {} missing", file);
    }
    
    // Check metadata for scientific validity
    let meta_path = output_dir.join("fft_metadata_block_0_0.json");
    let meta_content = std::fs::read_to_string(meta_path).unwrap();
    let meta_json: serde_json::Value = serde_json::from_str(&meta_content).unwrap();
    
    // Check Parseval Error
    // The key might be in "statistics" object
    if let Some(stats) = meta_json.get("statistics") {
        if let Some(error_val) = stats.get("parseval_error") {
            let parseval_error = error_val.as_f64().unwrap();
             // 1e-6 is the tolerance in the code, but we use slightly looser here just in case of float diffs
            assert!(parseval_error < 1e-5, "Parseval error too high: {}", parseval_error);
        } else {
             panic!("parseval_error missing in statistics");
        }
    } else {
        panic!("statistics missing in metadata");
    }
}
