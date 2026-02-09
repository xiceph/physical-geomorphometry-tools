use assert_cmd::Command;
use gdal::raster::Buffer;
use gdal::DriverManager;
use std::path::Path;
use serde_json::json;
use approx::assert_relative_eq;
use std::fs::File;
use std::io::Write;

fn create_mock_results(dir: &Path, scale: f64, size: usize) {
    let psd_path = dir.join("fft_psd_block_0_0.tif");
    let complex_path = dir.join("fft_complex_block_0_0.bin");
    let meta_path = dir.join("fft_metadata_block_0_0.json");

    let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
    let mut ds = driver.create_with_band_type::<f64, _>(&psd_path, size, size, 1).unwrap();
    
    // Power = P. Log-power saved is log10(P + 1e-12).
    // Let's use P = 1.0 * scale^2.
    let p_linear = 1.0 * scale * scale;
    let p_log = (p_linear + 1e-12_f64).log10();
    let data = vec![p_log; size * size];
    let mut band = ds.rasterband(1).unwrap();
    let mut buffer = Buffer::new((size, size), data);
    band.write((0, 0), (size, size), &mut buffer).unwrap();

    // Complex = scale + 0i
    let mut f = File::create(complex_path).unwrap();
    for _ in 0..(size * size) {
        f.write_all(&scale.to_le_bytes()).unwrap();
        f.write_all(&0.0_f64.to_le_bytes()).unwrap();
    }

    let metadata = json!({
        "block_position": [0, 0],
        "original_size": [size, size],
        "padded_size": [size, size],
        "statistics": {
            "f_nyquist": 0.5
        }
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();
}

#[test]
fn test_compare_identity() {
    let temp_dir = tempfile::tempdir().unwrap();
    let dir_a = temp_dir.path().join("a");
    let dir_b = temp_dir.path().join("b");
    let output_dir = temp_dir.path().join("output");
    std::fs::create_dir(&dir_a).unwrap();
    std::fs::create_dir(&dir_b).unwrap();

    create_mock_results(&dir_a, 1.0, 32);
    create_mock_results(&dir_b, 1.0, 32);

    let mut cmd = Command::cargo_bin("fft-compare").unwrap();
    cmd.arg("--input-a").arg(&dir_a)
       .arg("--input-b").arg(&dir_b)
       .arg("--output").arg(&output_dir)
       .assert()
       .success();

    // Check summary.csv
    let csv_path = output_dir.join("comparison_summary.csv");
    let mut rdr = csv::Reader::from_path(csv_path).unwrap();
    for result in rdr.records() {
        let record = result.unwrap();
        let mean_ratio: f64 = record[3].parse().unwrap();
        let coherence: f64 = record[4].parse().unwrap();
        
        assert_relative_eq!(mean_ratio, 1.0, epsilon = 1e-5);
        assert_relative_eq!(coherence, 1.0, epsilon = 1e-5);
    }
}

#[test]
fn test_compare_scaled() {
    let temp_dir = tempfile::tempdir().unwrap();
    let dir_a = temp_dir.path().join("a");
    let dir_b = temp_dir.path().join("b");
    let output_dir = temp_dir.path().join("output");
    std::fs::create_dir(&dir_a).unwrap();
    std::fs::create_dir(&dir_b).unwrap();

    create_mock_results(&dir_a, 1.0, 32);
    create_mock_results(&dir_b, 2.0, 32); // Scale complex by 2, power by 4

    let mut cmd = Command::cargo_bin("fft-compare").unwrap();
    cmd.arg("--input-a").arg(&dir_a)
       .arg("--input-b").arg(&dir_b)
       .arg("--output").arg(&output_dir)
       .assert()
       .success();

    let csv_path = output_dir.join("comparison_summary.csv");
    let mut rdr = csv::Reader::from_path(csv_path).unwrap();
    for result in rdr.records() {
        let record = result.unwrap();
        let mean_ratio: f64 = record[3].parse().unwrap();
        let coherence: f64 = record[4].parse().unwrap();
        
        assert_relative_eq!(mean_ratio, 4.0, epsilon = 1e-5);
        assert_relative_eq!(coherence, 1.0, epsilon = 1e-5);
    }
}
