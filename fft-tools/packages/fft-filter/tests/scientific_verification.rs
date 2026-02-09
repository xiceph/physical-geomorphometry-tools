use assert_cmd::Command;
use ndarray::Array2;
use num_complex::Complex;
use std::fs::File;
use std::io::{Write, Read};
use serde_json::json;
use approx::assert_relative_eq;

#[test]
fn test_filter_functionality() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_dir = temp_dir.path().join("input");
    let output_dir = temp_dir.path().join("output");
    std::fs::create_dir(&input_dir).unwrap();

    let size = 32;
    // let pixel_size = 1.0;
    
    // Create a spectrum with two peaks: 
    // 1. Low frequency (long wavelength): f = 1/16 = 0.0625 -> wavelength = 16
    // 2. High frequency (short wavelength): f = 1/4 = 0.25 -> wavelength = 4
    
    let mut spectrum = Array2::<Complex<f64>>::zeros((size, size));
    
    // d_f = 1 / (size * pixel_size) = 1/32 = 0.03125.
    // f1 = 0.0625 = 2 * d_f. Index offset = 2.
    // f2 = 0.25 = 8 * d_f. Index offset = 8.
    
    // We assume the spectrum is saved AFTER fftshift, so DC is at size/2.
    spectrum[[size/2, size/2 + 2]] = Complex::new(10.0, 0.0);
    spectrum[[size/2, size/2 - 2]] = Complex::new(10.0, 0.0);
    
    spectrum[[size/2, size/2 + 8]] = Complex::new(5.0, 0.0);
    spectrum[[size/2, size/2 - 8]] = Complex::new(5.0, 0.0);

    // Save complex bin
    let bin_path = input_dir.join("fft_complex_block_0_0.bin");
    let mut f = File::create(bin_path).unwrap();
    for c in spectrum.iter() {
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
        "statistics": {
            "f_nyquist": 0.5
        },
        "processing_params": {}
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();

    // Run fft-filter: Low-pass (remove wavelengths shorter than 10m)
    // f_lp = 1/10 = 0.1.
    // Our f1 (0.0625) < 0.1 -> Keep.
    // Our f2 (0.25) > 0.1 -> Remove.
    let mut cmd = Command::cargo_bin("fft-filter").unwrap();
    cmd.arg("--input").arg(&input_dir)
       .arg("--output").arg(&output_dir)
       .arg("--min-wavelength").arg("10")
       .arg("--taper-width").arg("0.1") // Narrow taper to be sure
       .assert()
       .success();

    // Load filtered complex and verify
    let filtered_bin_path = output_dir.join("fft_complex_block_0_0.bin");
    let mut f = File::open(filtered_bin_path).unwrap();
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer).unwrap();
    
    let filtered_data: Vec<f64> = buffer.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
    let filtered_spectrum: Vec<Complex<f64>> = filtered_data.chunks_exact(2).map(|c| Complex::new(c[0], c[1])).collect();
    let filtered_array = Array2::from_shape_vec((size, size), filtered_spectrum).unwrap();
    
    // Check f1 peak (should be kept)
    assert_relative_eq!(filtered_array[[size/2, size/2 + 2]].re, 10.0, epsilon = 1e-5);
    // Check f2 peak (should be removed)
    assert_relative_eq!(filtered_array[[size/2, size/2 + 8]].re, 0.0, epsilon = 1e-5);
}
