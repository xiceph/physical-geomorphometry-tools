use assert_cmd::Command;
use ndarray::Array2;
use std::fs::File;
use std::io::Write;
use serde_json::json;
use approx::assert_relative_eq;

#[test]
fn test_round_trip_accuracy() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_dir = temp_dir.path().join("input");
    let output_file = temp_dir.path().join("reconstructed.tif");
    std::fs::create_dir(&input_dir).unwrap();

    let size = 32;
    let mut data = Array2::<f64>::zeros((size, size));
    for r in 0..size {
        for c in 0..size {
            data[[r, c]] = (r as f64 + c as f64).sin();
        }
    }

    // 1. Manually compute FFT and save like fft-process would
    // We use fft-core directly to generate the mock data
    let fft_result = fft_core::compute_fft(&data, 1.0, 0, 0).unwrap();
    
    // Save complex bin
    let bin_path = input_dir.join("fft_complex_block_0_0.bin");
    let mut f = File::create(bin_path).unwrap();
    for c in fft_result.spectrum.iter() {
        f.write_all(&c.re.to_le_bytes()).unwrap();
        f.write_all(&c.im.to_le_bytes()).unwrap();
    }
    
    // Save metadata json
    let meta_path = input_dir.join("fft_metadata_block_0_0.json");
    let metadata = json!({
        "block_position": [0, 0],
        "original_size": [size, size],
        "padded_size": [size, size],
        "geo_transform": [0.0, 1.0, 0.0, 0.0, 0.0, -1.0],
        "wkt": null,
        "processing_params": {
            "detrend_order": null
        },
        "statistics": {}
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();

    // 2. Run fft-inverse
    let mut cmd = Command::cargo_bin("fft-inverse").unwrap();
    cmd.arg("--input").arg(&input_dir)
       .arg("--output").arg(&output_file)
       .assert()
       .success();

    // 3. Load result and compare
    let ds = gdal::Dataset::open(output_file).unwrap();
    let band = ds.rasterband(1).unwrap();
    let (cols, rows) = band.size();
    let recon_data = band.read_as::<f64>((0, 0), (cols, rows), (cols, rows), None).unwrap().data().to_vec();
    
    for r in 0..size {
        for c in 0..size {
            assert_relative_eq!(recon_data[r * size + c], data[[r, c]], epsilon = 1e-10);
        }
    }
}
