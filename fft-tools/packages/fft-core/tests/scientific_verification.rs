use fft_core::{compute_fft, apply_hann_window, detrend, reapply_trend};
use ndarray::Array2;
use approx::assert_relative_eq;

#[test]
fn test_parsevals_theorem_fft() {
    // Create a simple signal: sine wave
    let size = 32;
    let mut data = Array2::<f64>::zeros((size, size));
    for r in 0..size {
        for c in 0..size {
            data[[r, c]] = (2.0 * std::f64::consts::PI * r as f64 / size as f64).sin() 
                         + (4.0 * std::f64::consts::PI * c as f64 / size as f64).cos();
        }
    }

    let input_mean_sqr = data.mapv(|x| x.powi(2)).mean().unwrap();

    // Compute FFT
    let result = compute_fft(&data, 1.0, 0, 0).unwrap();

    // Sum of PSD should equal Mean Square of Input (Parseval's Theorem for this normalization)
    // PSD = |X|^2 / N^2
    // Sum(PSD) = Sum(|X|^2) / N^2
    // Parseval: Sum(|X|^2) = N * Sum(|x|^2)
    // Sum(PSD) = (N * Sum(|x|^2)) / N^2 = Sum(|x|^2) / N = Mean(|x|^2)
    let psd_sum = result.power_spectrum.sum();

    assert_relative_eq!(psd_sum, input_mean_sqr, epsilon = 1e-10);
}

#[test]
fn test_energy_conservation_impulse() {
    // Test with a delta function (impulse)
    // Parseval's holds for any signal.
    let size = 16;
    let mut data = Array2::<f64>::zeros((size, size));
    data[[size/2, size/2]] = 10.0;

    let input_mean_sqr = data.mapv(|x| x.powi(2)).mean().unwrap();
    let result = compute_fft(&data, 1.0, 0, 0).unwrap();
    let psd_sum = result.power_spectrum.sum();

    assert_relative_eq!(psd_sum, input_mean_sqr, epsilon = 1e-10);
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
            // Normalized coords used in detrend are -1 to 1
            let x = (2.0 * c as f64 / (size - 1) as f64) - 1.0;
            let y = (2.0 * r as f64 / (size - 1) as f64) - 1.0;
            data[[r, c]] = 2.0 * x + 3.0 * y + 5.0;
        }
    }
    
    let original = data.clone();

    // Detrend (Order 1)
    let coeffs = detrend(&mut data, 1).unwrap();
    
    // Residuals should be near zero (machine epsilon)
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
