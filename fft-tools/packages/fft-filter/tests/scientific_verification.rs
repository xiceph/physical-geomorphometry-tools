use assert_cmd::Command;
use ndarray::Array2;
use num_complex::Complex;
use std::fs::File;
use std::io::{Write, Read};
use serde_json::json;
use approx::assert_relative_eq;

/// Verifies that `fft-filter` correctly removes spectral components based on wavelength.
/// 
/// This test simulates a multi-scale signal consisting of two distinct landform 
/// wavelengths. It then applies a low-pass filter and verifies that only the 
/// long-wavelength (low-frequency) components are preserved in the complex spectrum.
#[test]
fn test_filter_functionality() {
    let temp_dir = tempfile::tempdir().unwrap();
    let input_dir = temp_dir.path().join("input");
    let output_dir = temp_dir.path().join("output");
    std::fs::create_dir(&input_dir).unwrap();

    let size = 32;
    let pixel_size = 1.0;
    
    // Create a spectrum with two distinct peaks: 
    // 1. Long wavelength (low frequency): λ = 16m -> f = 1/16 = 0.0625 Hz.
    // 2. Short wavelength (high frequency): λ = 4m -> f = 1/4 = 0.25 Hz.
    
    let mut spectrum = Array2::<Complex<f64>>::zeros((size, size));
    
    // Calculate frequency resolution: Δf = 1 / (N * Δx) = 1 / (32 * 1.0) = 0.03125 Hz.
    // Frequency bin for f1: 0.0625 / 0.03125 = 2. Index offset from DC = 2.
    // Frequency bin for f2: 0.25 / 0.03125 = 8. Index offset from DC = 8.
    
    // We assume the spectrum is saved AFTER fftshift (DC at center index).
    // f1 peak: preserved by filter (λ=16 > λ_min=10).
    spectrum[[size/2, size/2 + 2]] = Complex::new(10.0, 0.0);
    spectrum[[size/2, size/2 - 2]] = Complex::new(10.0, 0.0);
    
    // f2 peak: removed by filter (λ=4 < λ_min=10).
    spectrum[[size/2, size/2 + 8]] = Complex::new(5.0, 0.0);
    spectrum[[size/2, size/2 - 8]] = Complex::new(5.0, 0.0);

    // Save the complex spectrum to a binary file.
    let bin_path = input_dir.join("fft_complex_block_0_0.bin");
    let mut f = File::create(bin_path).unwrap();
    for c in spectrum.iter() {
        f.write_all(&c.re.to_le_bytes()).unwrap();
        f.write_all(&c.im.to_le_bytes()).unwrap();
    }
    
    // Save associated metadata.
    let meta_path = input_dir.join("fft_metadata_block_0_0.json");
    let metadata = json!({
        "block_position": [0, 0],
        "original_size": [size, size],
        "padded_size": [size, size],
        "geo_transform": [0.0, pixel_size, 0.0, 0.0, 0.0, -pixel_size],
        "statistics": {
            "f_nyquist": 0.5 / pixel_size
        },
        "processing_params": {}
    });
    std::fs::write(meta_path, serde_json::to_string_pretty(&metadata).unwrap()).unwrap();

    // Run the `fft-filter` CLI tool.
    // We apply a low-pass filter by keeping only wavelengths >= 10m.
    // Transition band is tapered to mitigate Gibbs phenomenon artifacts.
    let mut cmd = Command::new(assert_cmd::cargo_bin!("fft-filter"));
    cmd.arg("--input").arg(&input_dir)
       .arg("--output").arg(&output_dir)
       .arg("--min-wavelength").arg("10")
       .arg("--taper-width").arg("0.1") 
       .assert()
       .success();

    // Read the filtered binary spectrum and verify coefficients.
    let filtered_bin_path = output_dir.join("fft_complex_block_0_0.bin");
    let mut f = File::open(filtered_bin_path).unwrap();
    let mut buffer = Vec::new();
    f.read_to_end(&mut buffer).unwrap();
    
    let filtered_data: Vec<f64> = buffer.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().unwrap())).collect();
    let filtered_spectrum: Vec<Complex<f64>> = filtered_data.chunks_exact(2).map(|c| Complex::new(c[0], c[1])).collect();
    let filtered_array = Array2::from_shape_vec((size, size), filtered_spectrum).unwrap();
    
    // Validate that the long-wavelength component (f1) was preserved.
    assert_relative_eq!(filtered_array[[size/2, size/2 + 2]].re, 10.0, epsilon = 1e-5);
    
    // Validate that the short-wavelength component (f2) was effectively removed.
    assert_relative_eq!(filtered_array[[size/2, size/2 + 8]].re, 0.0, epsilon = 1e-5);
}
