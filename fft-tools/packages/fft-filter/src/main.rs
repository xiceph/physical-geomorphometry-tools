use anyhow::{Context, Result};
use clap::Parser;
use console::Term;
use fft_core::{prepare_output_dir, text, BlockMetadata};
use ndarray::{Array1, Array2, Zip};
use num_complex::Complex;
use rayon::prelude::*;
use serde_json::json;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use walkdir::WalkDir;

/// Command-line arguments for the fft-filter tool.
#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about,
    long_about = "Applies frequency-domain filters to complex FFT data."
)]
struct Args {
    /// Directory containing the fft_complex_block_*.bin and fft_metadata_block_*.json files.
    #[arg(long)]
    input: PathBuf,

    /// Directory to save the filtered complex data and metadata.
    #[arg(long)]
    output: PathBuf,

    /// The minimum spatial wavelength to preserve (meters). Wavelengths SHORTER than this will be attenuated (implements a low-pass filter).
    #[arg(long)]
    min_wavelength: Option<f64>,

    /// The maximum spatial wavelength to preserve (meters). Wavelengths LONGER than this will be attenuated (implements a high-pass filter).
    #[arg(long)]
    max_wavelength: Option<f64>,

    /// The width of the filter's cosine-tapered transition band, as a proportion of the cutoff frequency.
    #[arg(long, default_value_t = 0.5)]
    taper_width: f64,

    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    jobs: usize,
}

/// Main entry point for the fft-filter tool.
///
/// Parses command-line arguments, finds all metadata files in the input directory,
/// and processes each block in parallel to apply the specified frequency-domain filter.
fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for applying frequency-domain filters to complex FFT data.\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Filter"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    println!("{} Configuration:", text::bold("Filtering"));
    println!("  {:<20} {}", "Input Directory:", args.input.display());
    println!("  {:<20} {}", "Output Directory:", args.output.display());
    if let Some(min_w) = args.min_wavelength {
        println!("  {:<20} {} m", "Min Wavelength:", min_w);
    }
    if let Some(max_w) = args.max_wavelength {
        println!("  {:<20} {} m", "Max Wavelength:", max_w);
    }
    println!("  {:<20} {}", "Taper Width:", args.taper_width);
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

    let output_dir = prepare_output_dir(args.output.clone())?;

    let metadata_paths: Vec<PathBuf> = WalkDir::new(&args.input)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| {
            e.file_name()
                .to_string_lossy()
                .starts_with("fft_metadata_block_")
                && e.file_name().to_string_lossy().ends_with(".json")
        })
        .map(|e| e.path().to_path_buf())
        .collect();

    let n_blocks = metadata_paths.len();
    if n_blocks == 0 {
        println!("{} No blocks found to filter.", text::warning("!"));
        return Ok(());
    }

    let progress_counter = Arc::new(Mutex::new(0));

    print!("Filtering complex spectra...");
    io::stdout().flush().unwrap();

    let args_arc = std::sync::Arc::new(args);
    metadata_paths.par_iter().for_each(|path| {
        if let Err(e) = process_block(path, &output_dir, &args_arc) {
            eprintln!("\n{} Failed to process block for {}: {}", text::error("Error"), path.display(), e);
        }

        let mut count = progress_counter.lock().unwrap();
        *count += 1;
        let term = Term::stdout();
        let _ = term.clear_line();
        print!("\rFiltering complex spectra... {:.0}%", (*count as f32 / n_blocks as f32) * 100.0);
        let _ = io::stdout().flush();
    });

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All blocks filtered successfully.", text::check_icon());
    println!("{}", line);
    println!("{}", text::success("Filtering completed."));
    println!("");

    Ok(())
}

/// Processes a single FFT block: loads data, creates a filter, applies it, and saves the result.
///
/// # Arguments
/// * `metadata_path` - Path to the metadata file for the block to process.
/// * `output_dir` - The directory where the filtered output files will be saved.
/// * `args` - A reference to the command-line arguments containing filter parameters.
fn process_block(metadata_path: &Path, output_dir: &Path, args: &Args) -> Result<()> {
    // 1. Load metadata and complex spectrum data
    let metadata_str = fs::read_to_string(metadata_path)
        .with_context(|| format!("Failed to read metadata file: {:?}", metadata_path))?;
    let mut metadata: BlockMetadata = serde_json::from_str(&metadata_str)
        .with_context(|| format!("Failed to parse metadata: {:?}", metadata_path))?;

    let complex_path = metadata_path.with_file_name(
        metadata_path
            .file_name()
            .unwrap()
            .to_str()
            .unwrap()
            .replace("fft_metadata_block_", "fft_complex_block_")
            .replace(".json", ".bin"),
    );

    let mut file = fs::File::open(&complex_path)
        .with_context(|| format!("Failed to open complex data file: {:?}", complex_path))?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let (padded_rows, padded_cols) = metadata.padded_size;
    let complex_data_f64: Vec<f64> = buffer
        .chunks_exact(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect();

    let mut spectrum: Array2<Complex<f64>> =
        Array2::from_shape_fn((padded_rows, padded_cols), |(r, c)| {
            let idx = (r * padded_cols + c) * 2;
            Complex::new(complex_data_f64[idx], complex_data_f64[idx + 1])
        });

    // 2. Create the filter mask
    let pixel_size = 1.0 / (2.0 * metadata.statistics["f_nyquist"].as_f64().unwrap_or(0.5));
    let filter_mask = create_filter_mask(
        padded_rows,
        padded_cols,
        pixel_size,
        args.min_wavelength,
        args.max_wavelength,
        args.taper_width,
    )?;

    // 3. Apply the filter
    Zip::from(&mut spectrum)
        .and(&filter_mask)
        .for_each(|s, &m| *s *= m);

    // 4. Save the new filtered data and updated metadata
    let output_filename_base = metadata_path.file_name().unwrap().to_str().unwrap();

    // Save filtered complex spectrum
    let filtered_complex_path = output_dir.join(
        output_filename_base
            .replace("fft_metadata_block_", "fft_complex_block_")
            .replace(".json", ".bin"),
    );
    let mut file = fs::File::create(filtered_complex_path)?;
    for c in spectrum.iter() {
        file.write_all(&c.re.to_le_bytes())?;
        file.write_all(&c.im.to_le_bytes())?;
    }

    // Update and save metadata
    let filter_info = json!({
        "min_wavelength_meters": args.min_wavelength,
        "max_wavelength_meters": args.max_wavelength,
        "taper_width_proportion": args.taper_width,
    });
    // This assumes `processing_params` is a serde_json::Value::Object
    if let Some(obj) = metadata.processing_params.as_object_mut() {
        obj.insert("filter".to_string(), filter_info);
    }

    let filtered_metadata_path = output_dir.join(output_filename_base);
    let new_metadata_str = serde_json::to_string_pretty(&metadata)?;
    fs::write(filtered_metadata_path, new_metadata_str)?;

    Ok(())
}

/// Creates a 2D frequency-domain filter mask.
///
/// This function generates a 2D array (mask) that can be multiplied with a 2D spectrum.
/// It supports high-pass, low-pass, and band-pass filtering with a smooth cosine taper
/// at the cutoff frequencies to prevent ringing artifacts.
///
/// # Arguments
/// * `rows` - The number of rows in the spectrum (and the mask to be created).
/// * `cols` - The number of columns in the spectrum.
/// * `pixel_size` - The physical size of a pixel in the original spatial data (in meters).
/// * `min_wavelength` - If Some, specifies the cutoff for a high-pass filter in meters.
/// * `max_wavelength` - If Some, specifies the cutoff for a low-pass filter in meters.
/// * `taper_width` - The width of the transition band as a proportion of the cutoff frequency.
///
/// # Returns
/// A `Result` containing the 2D filter mask as an `Array2<f64>`.
fn create_filter_mask(
    rows: usize,
    cols: usize,
    pixel_size: f64,
    min_wavelength: Option<f64>,
    max_wavelength: Option<f64>,
    taper_width: f64,
) -> Result<Array2<f64>> {
    let mut freqs_x = fft_core::fftfreq(cols, pixel_size);
    fft_core::fftshift_1d(&mut freqs_x);
    let freqs_x = Array1::from_vec(freqs_x);

    let mut freqs_y = fft_core::fftfreq(rows, pixel_size);
    fft_core::fftshift_1d(&mut freqs_y);
    let freqs_y = Array1::from_vec(freqs_y);

    let mut mask = Array2::<f64>::ones((rows, cols));

    // Convert wavelength (m) cutoffs to frequency (cycles/m) cutoffs.
    // min-wavelength sets the LOW-PASS filter cutoff (passes wavelengths LONGER than this).
    let f_low_pass = min_wavelength.map(|l| 1.0 / l);
    // max-wavelength sets the HIGH-PASS filter cutoff (passes wavelengths SHORTER than this).
    let f_high_pass = max_wavelength.map(|l| 1.0 / l);

    for r in 0..rows {
        for c in 0..cols {
            let f_rad = (freqs_y[r].powi(2) + freqs_x[c].powi(2)).sqrt();

            let mut gain_lp = 1.0;
            let mut gain_hp = 1.0;

            // Low-pass filter logic (from --min-wavelength):
            // Attenuates frequencies HIGHER than the cutoff.
            if let Some(f_lp) = f_low_pass {
                let f_start_taper = f_lp * (1.0 - taper_width / 2.0);
                let f_end_taper = f_lp * (1.0 + taper_width / 2.0);

                if f_rad > f_start_taper {
                    gain_lp = if f_rad >= f_end_taper {
                        0.0
                    } else {
                        // Cosine taper from 1 down to 0
                        0.5 * (1.0
                            + (std::f64::consts::PI * (f_rad - f_start_taper)
                                / (f_end_taper - f_start_taper))
                                .cos())
                    };
                }
            }

            // High-pass filter logic (from --max-wavelength):
            // Attenuates frequencies LOWER than the cutoff.
            if let Some(f_hp) = f_high_pass {
                let f_start_taper = f_hp * (1.0 - taper_width / 2.0);
                let f_end_taper = f_hp * (1.0 + taper_width / 2.0);

                if f_rad < f_end_taper {
                    gain_hp = if f_rad <= f_start_taper {
                        0.0
                    } else {
                        // Cosine taper from 0 up to 1
                        0.5 * (1.0
                            - (std::f64::consts::PI * (f_rad - f_start_taper)
                                / (f_end_taper - f_start_taper))
                                .cos())
                    };
                }
            }

            // Combine gains for band-pass filtering.
            mask[[r, c]] = gain_lp * gain_hp;
        }
    }

    Ok(mask)
}
