use fft_core::{compute_fft, apply_hann_window, detrend, reapply_trend};
use ndarray::Array2;
use approx::assert_relative_eq;

#[test]
fn test_parsevals_theorem_fft() {
    // Create a simple signal: sine wave
    let size = 32;
    let pixel_size = 2.0; // Non-unit pixel size
    let mut data = Array2::<f64>::zeros((size, size));
    for r in 0..size {
        for c in 0..size {
            data[[r, c]] = (2.0 * std::f64::consts::PI * r as f64 / size as f64).sin() 
                         + (4.0 * std::f64::consts::PI * c as f64 / size as f64).cos();
        }
    }

    // Spatial variance (mean of squares for zero-mean signal)
    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();

    // Compute FFT
    let result = compute_fft(&data, pixel_size, 0, 0, (size, size)).unwrap();

    // Physically-correct Parseval: Σ PSD * dkx * dky = Variance
    let (padded_rows, padded_cols) = result.power_spectrum.dim();
    let dkx = 1.0 / (padded_cols as f64 * pixel_size);
    let dky = 1.0 / (padded_rows as f64 * pixel_size);
    let psd_integral = result.power_spectrum.sum() * dkx * dky;

    assert_relative_eq!(psd_integral, spatial_variance, epsilon = 1e-10);
}

#[test]
fn test_physical_psd_normalization_padding_independence() {
    let size = 32;
    let pixel_size = 1.0;
    let mut data = Array2::<f64>::zeros((size, size));
    data[[size/2, size/2]] = 1.0; // Impulse

    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();

    // Case 1: No additional padding
    let result1 = compute_fft(&data, pixel_size, 0, 0, (size, size)).unwrap();
    let dkx1 = 1.0 / (result1.power_spectrum.dim().1 as f64 * pixel_size);
    let dky1 = 1.0 / (result1.power_spectrum.dim().0 as f64 * pixel_size);
    let psd_integral1 = result1.power_spectrum.sum() * dkx1 * dky1;

    // Case 2: Padded data (manually padded with zeros)
    let padded_size = 64;
    let mut padded_data = Array2::<f64>::zeros((padded_size, padded_size));
    padded_data.slice_mut(ndarray::s![0..size, 0..size]).assign(&data);

    let result2 = compute_fft(&padded_data, pixel_size, 0, 0, (size, size)).unwrap();
    let dkx2 = 1.0 / (result2.power_spectrum.dim().1 as f64 * pixel_size);
    let dky2 = 1.0 / (result2.power_spectrum.dim().0 as f64 * pixel_size);
    let psd_integral2 = result2.power_spectrum.sum() * dkx2 * dky2;

    // Both should yield the same variance
    assert_relative_eq!(psd_integral1, spatial_variance, epsilon = 1e-10);
    assert_relative_eq!(psd_integral2, spatial_variance, epsilon = 1e-10);
    
    // Crucially, the PSD amplitude should be independent of padding size.
    assert_relative_eq!(result1.power_spectrum.mean().unwrap(), result2.power_spectrum.mean().unwrap(), epsilon = 1e-10);
    // Wait, let's re-check:
    // result1: PSD1 = |Z1|^2 * area / N_orig
    // result2: PSD2 = |Z2|^2 * area / N_orig
    // Z1_k = Σ_{n=0..N1-1} z_n exp(...)
    // Z2_k = Σ_{n=0..N2-1} z_n exp(...) where z_n=0 for n>=N1.
    // At the same frequencies, Z1_k == Z2_k.
    // So PSD1 == PSD2 at the same frequencies.
    // This is EXACTLY what we want: padding doesn't change the PSD value, it just interpolates.
}

#[test]
fn test_energy_conservation_impulse() {
    let size = 16;
    let pixel_size = 1.0;
    let mut data = Array2::<f64>::zeros((size, size));
    data[[size/2, size/2]] = 10.0;

    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();
    let result = compute_fft(&data, pixel_size, 0, 0, (size, size)).unwrap();
    
    let dkx = 1.0 / (result.power_spectrum.dim().1 as f64 * pixel_size);
    let dky = 1.0 / (result.power_spectrum.dim().0 as f64 * pixel_size);
    let psd_integral = result.power_spectrum.sum() * dkx * dky;

    assert_relative_eq!(psd_integral, spatial_variance, epsilon = 1e-10);
}

#[test]
fn test_hann_window_reduces_energy() {
    let size = 32;
    let mut data = Array2::<f64>::ones((size, size));
    let initial_energy = data.mapv(|x| x.powi(2)).sum();

    apply_hann_window(&mut data, 4).unwrap();

    let final_energy = data.mapv(|x| x.powi(2)).sum();
    assert!(final_energy < initial_energy, "Hann window should reduce total energy");
}

#[test]
fn test_detrend_reapply_linear() {
    let size = 10;
    let mut data = Array2::<f64>::zeros((size, size));
    // Add linear trend: z = 2x + 3y + 5
    for r in 0..size {
        for c in 0..size {
            let x = (2.0 * c as f64 / (size - 1) as f64) - 1.0;
            let y = (2.0 * r as f64 / (size - 1) as f64) - 1.0;
            data[[r, c]] = 2.0 * x + 3.0 * y + 5.0;
        }
    }
    
    let original = data.clone();

    // Detrend (Order 1)
    let coeffs = detrend(&mut data, 1).unwrap();
    
    // Residuals should be near zero
    let max_residual = data.mapv(|x| x.abs()).fold(0.0f64, |a, b| a.max(*b));
    assert_relative_eq!(max_residual, 0.0, epsilon = 1e-10);

    // Reapply
    reapply_trend(&mut data, &coeffs, 1).unwrap();

    // Should match original
    for r in 0..size {
        for c in 0..size {
            assert_relative_eq!(data[[r, c]], original[[r, c]], epsilon = 1e-10);
        }
    }
}
