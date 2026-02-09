use anyhow::{Context, Result};
use clap::Parser;
use console::Term;
use fft_core::{prepare_output_dir, save_gdal_raster, text};
use gdal::Dataset;
use ndarray::{Array1, Array2};
use rayon::prelude::*;
use serde::Deserialize;
use serde_json::json;
use std::f64::consts::PI;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

/// Helper function to round a f64 value to a specified number of decimal places.
fn round_f64(value: f64, decimal_places: u32) -> f64 {
    let factor = 10.0f64.powi(decimal_places as i32);
    (value * factor).round() / factor
}

/// Command-line arguments for the fft-polar tool.
#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about,
    long_about = "Converts Cartesian FFT power spectrums to a polar representation."
)]
struct Args {
    /// Directory containing the fft_psd_block_*.tif files.
    #[arg(long)]
    input: PathBuf,

    /// Directory to save the polar spectrum files.
    #[arg(long)]
    output: PathBuf,

    /// The number of angular bins for the polar spectrum (degrees).
    #[arg(long, default_value_t = 36)]
    n_angles: usize,

    /// The number of wavenumber (radial) bins for the polar spectrum.
    #[arg(long, default_value_t = 64)]
    n_wavenumbers: usize,

    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    jobs: usize,
}

/// Structure to deserialize block metadata from JSON files.
#[derive(Deserialize, Debug)]
struct BlockMetadata {
    original_size: (usize, usize),
    #[serde(default)]
    statistics: Stats,
}

/// Statistics loaded from block metadata, including Nyquist frequency.
#[derive(Deserialize, Debug, Default)]
struct Stats {
    f_nyquist: f64,
}

/// Holds the calculated polar spectrum data.
struct PolarSpectrumData {
    power: Array2<f64>,
    wavelengths: Array1<f64>,
    angles: Array1<f64>,
}

/// Main entry point for the fft-polar tool.
///
/// Parses command-line arguments, prepares output directory,
/// and processes each PSD file in parallel.
fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for converting Cartesian PSDs to a polar representation.\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Polar Transformer"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    println!("{} Configuration:", text::bold("Transformation"));
    println!("  {:<20} {}", "Input Directory:", args.input.display());
    println!("  {:<20} {}", "Output Directory:", args.output.display());
    println!("  {:<20} {}", "Angular Bins:", args.n_angles);
    println!("  {:<20} {}", "Wavenumber Bins:", args.n_wavenumbers);
    println!(
        "  {:<20} {}",
        "Parallel Jobs:",
        if args.jobs == 0 {
            "all available cores".to_string()
        } else {
            args.jobs.to_string()
        }
    );
    println!("{}\n", dline);

    if args.jobs > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.jobs)
            .build_global()?;
    }

    // Prepare the output directory, handling existing directories by appending an index.
    let output_dir = prepare_output_dir(args.output.clone())?;

    // Collect all PSD files from the input directory.
    let entries: Vec<_> = fs::read_dir(&args.input)?
        .filter_map(Result::ok)
        .filter(|e| {
            e.path().extension().is_some_and(|ext| ext == "tif")
                && e.file_name()
                    .to_string_lossy()
                    .starts_with("fft_psd_block_")
        })
        .collect();

    let n_files = entries.len();
    if n_files == 0 {
        println!("{} No PSD files found to process.", text::warning("!"));
        return Ok(());
    }

    let progress_counter = Arc::new(Mutex::new(0));

    print!("Transforming to polar coordinates...");
    std::io::stdout().flush().unwrap();

    // Process each PSD file in parallel using Rayon.
    let args_arc = std::sync::Arc::new(args);
    entries.par_iter().for_each(|entry| {
        let args_clone = args_arc.clone();
        if let Err(e) = process_file(&entry.path(), &output_dir, &args_clone) {
            eprintln!("\n{} Failed to process {}: {}", text::error("Error"), entry.path().display(), e);
        }

        let mut count = progress_counter.lock().unwrap();
        *count += 1;
        let term = Term::stdout();
        let _ = term.clear_line();
        print!("\rTransforming to polar coordinates... {:.0}%", (*count as f32 / n_files as f32) * 100.0);
        let _ = std::io::stdout().flush();
    });

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All files transformed successfully.", text::check_icon());
    println!("{}", line);
    println!("{}", text::success("Polar transformation completed."));
    println!("");

    Ok(())
}

type PsdData = (Array2<f64>, Array1<f64>, Array1<f64>, f64, (usize, usize));

/// Processes a single PSD file: loads data, calculates polar spectrum,
/// and saves results.
///
/// # Arguments
/// * `psd_path` - Path to the input PSD TIFF file.
/// * `output_dir` - Path to the directory where output polar spectrums will be saved.
/// * `args` - Command-line arguments for configuration.
///
/// # Returns
/// A `Result` indicating success or failure.
fn process_file(psd_path: &Path, output_dir: &Path, args: &Args) -> Result<()> {
    let (psd, freqs_x, freqs_y, pixel_size, original_size) = load_psd(psd_path)?;

    // Calculate the polar representation of the PSD.
    let polar_spectrum_data = calculate_polar_spectrum(
        &psd,
        &freqs_x,
        &freqs_y,
        pixel_size,
        original_size.0,
        args.n_angles,
        args.n_wavenumbers,
    )?;

    // Save the polar spectrum to a GeoTIFF. Log transform for better dynamic range storage.
    let file_stem = psd_path.file_stem().unwrap().to_str().unwrap();
    let polar_path = output_dir.join(format!("{}.tif", file_stem.replace("fft_psd", "fft_polar")));
    save_gdal_raster(
        &polar_spectrum_data.power.mapv(|p| {
            if p.is_nan() || p <= 0.0 {
                f64::NAN
            } else {
                p.log10()
            }
        }),
        &polar_path,
        None,
        None,
        Some(f64::NAN),
    )?;
    // Save metadata about the polar axes (wavelengths and angles) to a JSON file.
    let meta_path = output_dir.join(format!(
        "{}.json",
        file_stem.replace("fft_psd", "fft_polar_metadata")
    ));
    let rounded_wavelengths: Vec<f64> = polar_spectrum_data
        .wavelengths
        .iter()
        .map(|&v| round_f64(v, 4))
        .collect();
    let rounded_angles: Vec<f64> = polar_spectrum_data
        .angles
        .iter()
        .map(|&v| round_f64(v, 1))
        .collect();

    let metadata = json!({
        "wavelengths": rounded_wavelengths,
        "angles": rounded_angles,
        "n_wavelengths": rounded_wavelengths.len(),
        "n_angles": rounded_angles.len(),
    });
    let json_string = serde_json::to_string_pretty(&metadata)?;
    std::fs::write(meta_path, json_string)?;

    Ok(())
}

/// Loads a Cartesian PSD spectrum from a TIFF file and its associated metadata.
///
/// # Arguments
/// * `psd_path` - Path to the input PSD TIFF file.
///
/// # Returns
/// A `Result` containing a tuple:
///   - `Array2<f64>`: The PSD data (linear scale).
///   - `Array1<f64>`: X-frequencies.
///   - `Array1<f64>`: Y-frequencies.
///   - `f64`: Pixel size used during FFT.
///   - `(usize, usize)`: The original (unpadded) size of the block.
fn load_psd(psd_path: &Path) -> Result<PsdData> {
    let dataset = Dataset::open(psd_path)?;
    let band = dataset.rasterband(1)?;
    let (cols, rows) = band.size();

    // Read log-transformed PSD and convert back to linear scale.
    let psd_log = band
        .read_as::<f64>((0, 0), (cols, rows), (cols, rows), None)?
        .data()
        .to_vec();
    let psd = Array2::from_shape_vec((rows, cols), psd_log)?.mapv(|p| 10.0_f64.powf(p) - 1e-12);

    // Construct path to the corresponding metadata JSON file.
    let meta_path = psd_path
        .with_file_name(
            psd_path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .replace("fft_psd_block_", "fft_metadata_block_"),
        )
        .with_extension("json");

    // Load metadata to retrieve Nyquist frequency and pixel size.
    let meta_file = fs::File::open(&meta_path)
        .with_context(|| format!("Failed to open metadata file: {:?}", meta_path))?;
    let metadata: BlockMetadata = serde_json::from_reader(meta_file)?;
    let f_nyquist = metadata.statistics.f_nyquist;
    let pixel_size = 1.0 / (2.0 * f_nyquist);
    let original_size = metadata.original_size;

    // Calculate and shift frequencies for X and Y axes.
    let mut freqs_x = fft_core::fftfreq(cols, pixel_size);
    let mut freqs_y = fft_core::fftfreq(rows, pixel_size);
    fft_core::fftshift_1d(&mut freqs_x);
    fft_core::fftshift_1d(&mut freqs_y);

    Ok((
        psd,
        Array1::from(freqs_x),
        Array1::from(freqs_y),
        pixel_size,
        original_size,
    ))
}

/// Calculates the polar representation of a 2D Cartesian PSD using Jacobian-weighted area correction.
///
/// This function bins the Cartesian PSD values into concentric wavenumber rings
/// and angular sectors, normalizing by the physical area of each polar bin
/// to ensure conservation of power.
///
/// # Arguments
/// * `psd` - The 2D Cartesian Power Spectral Density data.
/// * `freqs_x` - Array of frequencies along the X-axis.
/// * `freqs_y` - Array of frequencies along the Y-axis.
/// * `pixel_size` - The physical size of a pixel in the original spatial data.
/// * `original_window_size` - The original (unpadded) window size in pixels.
/// * `n_angles` - The number of angular bins.
/// * `n_wavenumbers` - The number of radial wavenumber (k) bins.
///
/// # Returns
/// A `Result` containing `PolarSpectrumData` with the polar PSD, wavelengths, and angles.
fn calculate_polar_spectrum(
    psd: &Array2<f64>,
    freqs_x: &Array1<f64>,
    freqs_y: &Array1<f64>,
    pixel_size: f64,
    original_window_size: usize,
    n_angles: usize,
    n_wavenumbers: usize,
) -> Result<PolarSpectrumData> {
    let (rows, cols) = psd.dim();
    let dkx = (freqs_x[1] - freqs_x[0]).abs();
    let dky = (freqs_y[1] - freqs_y[0]).abs();
    let cartesian_cell_area = dkx * dky;

    // 1. Define bins in k-space (wavenumber) and theta-space (angle)
    // k_max corresponds to the Nyquist frequency (1 / (2 * pixel_size)).
    let k_max = 1.0 / (2.0 * pixel_size);
    // k_min corresponds to the largest wavelength to be analyzed (window_size / 2).
    let k_min = 2.0 / (original_window_size as f64 * pixel_size);

    let log_k_min = k_min.log10();
    let log_k_max = k_max.log10();

    let k_edges = Array1::logspace(10.0, log_k_min, log_k_max, n_wavenumbers + 1);
    let theta_edges = Array1::linspace(0.0, PI, n_angles + 1);

    // 2. Pre-calculate the physical area of each polar bin
    let mut polar_bin_areas = Array2::<f64>::zeros((n_wavenumbers, n_angles));
    for k_idx in 0..n_wavenumbers {
        let k_inner = k_edges[k_idx];
        let k_outer = k_edges[k_idx + 1];
        for a_idx in 0..n_angles {
            let theta_start = theta_edges[a_idx];
            let theta_end = theta_edges[a_idx + 1];
            let d_theta = theta_end - theta_start;
            // Area of an annulus sector: 0.5 * (r_outer^2 - r_inner^2) * d_theta
            polar_bin_areas[[k_idx, a_idx]] = 0.5 * (k_outer.powi(2) - k_inner.powi(2)) * d_theta;
        }
    }

    // 3. Accumulate power into polar bins
    let mut power_sum = Array2::<f64>::zeros((n_wavenumbers, n_angles));
    for y_idx in 0..rows {
        for x_idx in 0..cols {
            let fx = freqs_x[x_idx];
            let fy = freqs_y[y_idx];
            let power = psd[[y_idx, x_idx]];

            let k = (fx.powi(2) + fy.powi(2)).sqrt();
            if k < k_min {
                continue;
            }

            // atan2 returns values in [-PI, PI]. We want [0, PI]
            let mut theta = fy.atan2(fx);
            if theta < 0.0 {
                theta += PI;
            }
            if theta >= PI {
                theta = PI - 1e-9;
            }

            // Find bins using searchsorted.
            let k_idx = match k_edges
                .as_slice()
                .unwrap()
                .binary_search_by(|v| v.partial_cmp(&k).unwrap())
            {
                Ok(i) => i.saturating_sub(1),
                Err(i) => i.saturating_sub(1),
            };
            let a_idx = match theta_edges
                .as_slice()
                .unwrap()
                .binary_search_by(|v| v.partial_cmp(&theta).unwrap())
            {
                Ok(i) => i.saturating_sub(1),
                Err(i) => i.saturating_sub(1),
            };

            if k_idx < n_wavenumbers && a_idx < n_angles {
                power_sum[[k_idx, a_idx]] += power * cartesian_cell_area;
            }
        }
    }

    // 4. Normalize by polar bin area
    // This is the "Jacobian" part of the transformation. Because polar bins 
    // increase in size as k increases (Area ~ k * dk * dtheta), we must divide 
    // by the area to convert from "Total Power in Bin" back to "Power Density".
    // This ensures that the resulting PSD is physically comparable across all scales.
    let mut polar_power = Array2::<f64>::from_elem((n_wavenumbers, n_angles), f64::NAN);
    for k_idx in 0..n_wavenumbers {
        for a_idx in 0..n_angles {
            let area = polar_bin_areas[[k_idx, a_idx]];
            if area > 1e-12 {
                polar_power[[k_idx, a_idx]] = power_sum[[k_idx, a_idx]] / area;
            }
        }
    }

    // 5. Convert k-bin centers to wavelengths and theta-bin centers to degrees for display
    let k_centers: Array1<f64> = k_edges
        .windows(2)
        .into_iter()
        .map(|w| (w[0] + w[1]) / 2.0)
        .collect();
    let final_wavelengths = k_centers.mapv(|k| 1.0 / k);

    let theta_centers: Array1<f64> = theta_edges
        .windows(2)
        .into_iter()
        .map(|w| (w[0] + w[1]) / 2.0)
        .collect();
    let final_angles = theta_centers.mapv(|t| t.to_degrees());

    Ok(PolarSpectrumData {
        power: polar_power,
        wavelengths: final_wavelengths,
        angles: final_angles,
    })
}
