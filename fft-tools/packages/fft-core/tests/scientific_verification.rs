use fft_core::{compute_fft, apply_hann_window, detrend, reapply_trend};
use ndarray::Array2;
use approx::assert_relative_eq;

/// Verifies that the computed PSD satisfies the physical form of Parseval's Theorem.
///
/// For a physically-normalized PSD, the integral of the PSD over all frequencies
/// (Σ PSD * dkx * dky) must equal the spatial variance of the signal.
#[test]
fn test_parsevals_theorem_fft() {
    // Create a simple synthetic signal: a combination of sine and cosine waves.
    let size = 32;
    let pixel_size_x = 2.0; 
    let pixel_size_y = 3.0; // Test anisotropic scaling
    let mut data = Array2::<f64>::zeros((size, size));
    for r in 0..size {
        for c in 0..size {
            data[[r, c]] = (2.0 * std::f64::consts::PI * r as f64 / size as f64).sin() 
                         + (4.0 * std::f64::consts::PI * c as f64 / size as f64).cos();
        }
    }

    // Calculate spatial variance. Since the signal is constructed to have zero mean,
    // variance is simply the mean of the squares.
    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();

    // Compute the FFT with physical normalization.
    let result = compute_fft(&data, pixel_size_x, pixel_size_y, 0, 0, (size, size)).unwrap();

    // Physical Parseval's check: Σ PSD * dkx * dky = Variance.
    // frequency spacing dkx = 1 / (N_padded * dx).
    let (padded_rows, padded_cols) = result.power_spectrum.dim();
    let dkx = 1.0 / (padded_cols as f64 * pixel_size_x);
    let dky = 1.0 / (padded_rows as f64 * pixel_size_y);
    let psd_integral = result.power_spectrum.sum() * dkx * dky;

    assert_relative_eq!(psd_integral, spatial_variance, epsilon = 1e-10);
}

/// Verifies that the PSD amplitude (spectral density) is independent of the zero-padding width.
///
/// Increasing the FFT grid size through zero-padding should interpolate the spectrum
/// but must not change the absolute spectral density values at existing frequencies.
#[test]
fn test_physical_psd_normalization_padding_independence() {
    let size = 32;
    let pixel_size_x = 1.0;
    let pixel_size_y = 1.0;
    let mut data = Array2::<f64>::zeros((size, size));
    data[[size/2, size/2]] = 1.0; // Use an impulse signal for wide-band spectral coverage.

    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();

    // Case 1: Compute FFT with no additional padding.
    let result1 = compute_fft(&data, pixel_size_x, pixel_size_y, 0, 0, (size, size)).unwrap();
    let dkx1 = 1.0 / (result1.power_spectrum.dim().1 as f64 * pixel_size_x);
    let dky1 = 1.0 / (result1.power_spectrum.dim().0 as f64 * pixel_size_y);
    let psd_integral1 = result1.power_spectrum.sum() * dkx1 * dky1;

    // Case 2: Manually pad the same data with zeros to a larger grid.
    let padded_size = 64;
    let mut padded_data = Array2::<f64>::zeros((padded_size, padded_size));
    padded_data.slice_mut(ndarray::s![0..size, 0..size]).assign(&data);

    // Use the SAME original_size (32, 32) even though input data is now (64, 64).
    let result2 = compute_fft(&padded_data, pixel_size_x, pixel_size_y, 0, 0, (size, size)).unwrap();
    let dkx2 = 1.0 / (result2.power_spectrum.dim().1 as f64 * pixel_size_x);
    let dky2 = 1.0 / (result2.power_spectrum.dim().0 as f64 * pixel_size_y);
    let psd_integral2 = result2.power_spectrum.sum() * dkx2 * dky2;

    // Both cases must yield the same physical variance.
    assert_relative_eq!(psd_integral1, spatial_variance, epsilon = 1e-10);
    assert_relative_eq!(psd_integral2, spatial_variance, epsilon = 1e-10);
    
    // Crucially, the absolute PSD values (density) should be identical.
    // For an impulse signal, Z_k values are identical regardless of padding,
    // so the PSD mean should remain constant.
    assert_relative_eq!(result1.power_spectrum.mean().unwrap(), result2.power_spectrum.mean().unwrap(), epsilon = 1e-10);
}

/// Verifies energy conservation for an impulse signal.
#[test]
fn test_energy_conservation_impulse() {
    let size = 16;
    let pixel_size_x = 1.0;
    let pixel_size_y = 1.0;
    let mut data = Array2::<f64>::zeros((size, size));
    data[[size/2, size/2]] = 10.0;

    let spatial_variance = data.mapv(|x| x.powi(2)).mean().unwrap();
    let result = compute_fft(&data, pixel_size_x, pixel_size_y, 0, 0, (size, size)).unwrap();
    
    let dkx = 1.0 / (result.power_spectrum.dim().1 as f64 * pixel_size_x);
    let dky = 1.0 / (result.power_spectrum.dim().0 as f64 * pixel_size_y);
    let psd_integral = result.power_spectrum.sum() * dkx * dky;

    assert_relative_eq!(psd_integral, spatial_variance, epsilon = 1e-10);
}

/// Verifies that applying a Hann window reduces the total signal energy as expected.
#[test]
fn test_hann_window_reduces_energy() {
    let size = 32;
    let mut data = Array2::<f64>::ones((size, size));
    let initial_energy = data.mapv(|x| x.powi(2)).sum();

    apply_hann_window(&mut data, 4).unwrap();

    let final_energy = data.mapv(|x| x.powi(2)).sum();
    assert!(final_energy < initial_energy, "Hann windowing must reduce total energy due to edge attenuation.");
}

/// Verifies the round-trip consistency of polynomial detrending and reapplication.
#[test]
fn test_detrend_reapply_linear() {
    let size = 10;
    let mut data = Array2::<f64>::zeros((size, size));
    // Add a known linear trend: z = 2x + 3y + 5.
    for r in 0..size {
        for c in 0..size {
            // Internal detrending uses normalized coordinates from -1 to 1.
            let x = (2.0 * c as f64 / (size - 1) as f64) - 1.0;
            let y = (2.0 * r as f64 / (size - 1) as f64) - 1.0;
            data[[r, c]] = 2.0 * x + 3.0 * y + 5.0;
        }
    }
    
    let original = data.clone();

    // Perform 1st order (planar) detrending.
    let coeffs = detrend(&mut data, 1).unwrap();
    
    // After perfect detrending of a planar surface, residuals must be near zero.
    let max_residual = data.mapv(|x| x.abs()).fold(0.0f64, |a, b| a.max(*b));
    assert_relative_eq!(max_residual, 0.0, epsilon = 1e-10);

    // Reapply the original trend.
    reapply_trend(&mut data, &coeffs, 1).unwrap();

    // The reconstructed surface must match the original input.
    for r in 0..size {
        for c in 0..size {
            assert_relative_eq!(data[[r, c]], original[[r, c]], epsilon = 1e-10);
        }
    }
}
