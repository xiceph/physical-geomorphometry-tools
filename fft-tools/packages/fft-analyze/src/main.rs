use anyhow::Result;
use clap::Parser;
use console::Term;
use fft_core::{polyfit_1d, text};
use ndarray::{stack, Array1, Axis};
use rayon::prelude::*;
use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

mod cli;
mod io;
mod plotting;
mod processing;

use cli::{AnalysisMode, Args};
use io::{get_axis_info, save_to_csv};
use plotting::{
    generate_and_save_line_chart, generate_and_save_residuals_chart, generate_and_save_rose_chart,
};
use processing::{parse_bounds, process_file};

/// Finds all fft_polar_block_*.tif files in the given directory.
fn find_polar_spectrum_files(input_dir: &PathBuf) -> Result<Vec<PathBuf>> {
    let entries: Vec<_> = fs::read_dir(input_dir)?
        .filter_map(Result::ok)
        .filter(|e| {
            e.path().extension().is_some_and(|ext| ext == "tif")
                && e.file_name()
                    .to_string_lossy()
                    .starts_with("fft_polar_block_")
        })
        .map(|e| e.path())
        .collect();

    if entries.is_empty() {
        anyhow::bail!(
            "No polar spectrum files found to process in {}",
            input_dir.display()
        );
    }
    Ok(entries)
}

/// Aggregates results from all blocks by performing a NaN-aware mean.
fn aggregate_results(results: Vec<Array1<f64>>) -> Result<Array1<f64>> {
    let stacked_results = stack(
        Axis(0),
        &results.iter().map(|a| a.view()).collect::<Vec<_>>(),
    )?;
    let mut final_average_psd = Array1::<f64>::zeros(stacked_results.ncols());
    for (j, col) in stacked_results.axis_iter(Axis(1)).enumerate() {
        let valid_values: Vec<f64> = col.iter().filter(|v| !v.is_nan()).cloned().collect();
        if !valid_values.is_empty() {
            final_average_psd[j] = valid_values.iter().sum::<f64>() / valid_values.len() as f64;
        } else {
            final_average_psd[j] = f64::NAN;
        }
    }
    Ok(final_average_psd)
}

/// Main entry point for the fft-analyze tool.
///
/// Parses command-line arguments, processes polar spectrum files, aggregates results,
/// and saves them to a CSV file and an optional SVG plot.
fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for statistical analysis and visualization of polar FFT spectra.\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Analyzer"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    // Validate that --detrend is only used with --mode radial-mean
    if args.detrend.is_some() && !matches!(args.mode, AnalysisMode::RadialMean) {
        anyhow::bail!("The --detrend option can only be used with --mode radial-mean.");
    }

    if args.jobs > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.jobs)
            .build_global()?;
    }

    let entries = find_polar_spectrum_files(&args.input)?;
    let n_files = entries.len();

    // Parse wavelength and angle bounds from command-line arguments.
    let wavelength_bounds = parse_bounds(args.wavelength_bounds.clone())?;
    let angle_bounds = parse_bounds(args.angle_bounds.clone())?;

    println!("{} Configuration:", text::bold("Analysis"));
    println!("  {:<20} {}", "Input Directory:", args.input.display());
    println!("  {:<20} {}", "Output File:", args.output.display());
    println!(
        "  {:<20} {}",
        "Analysis Mode:",
        match args.mode {
            AnalysisMode::RadialMean => "Radial Mean",
            AnalysisMode::AngularMean => "Angular Mean",
        }
    );
    if let Some(bounds) = wavelength_bounds {
        println!("  {:<20} {:?}", "Wavelength Bounds:", bounds);
    }
    if let Some(bounds) = angle_bounds {
        println!("  {:<20} {:?}", "Angle Bounds:", bounds);
    }
    if let Some(order) = args.detrend {
        println!("  {:<20} {}", "Detrend Order:", order);
    }
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

    let progress_counter = Arc::new(Mutex::new(0));

    print!("Analyzing polar spectra...");
    std::io::stdout().flush().unwrap();

    // Process each polar spectrum file in parallel.
    let results: Vec<Array1<f64>> = entries
        .par_iter()
        .map(|entry| {
            let profile = process_file(entry, args.mode.clone(), wavelength_bounds, angle_bounds)?;

            if args.save_partials {
                // Get axis info and header for this specific file
                let (axis_labels, header) = get_axis_info(
                    entry,
                    &args.mode,
                    wavelength_bounds,
                    angle_bounds,
                    None, // Detrending is not applied to partials
                )?;

                // Derive suffix from input filename (e.g., "fft_polar_block_0.0" -> "block_0.0")
                let input_stem = entry
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .replace("fft_polar_", "");

                // Construct partial CSV path
                let output_stem = args
                    .output
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("output");
                let output_dir = args.output.parent().unwrap_or(std::path::Path::new("."));
                let partial_csv_name = format!("{}_{}.csv", output_stem, input_stem);
                let partial_csv_path = output_dir.join(partial_csv_name);

                save_to_csv(
                    &partial_csv_path,
                    &header,
                    &axis_labels,
                    &profile,
                    None,
                    None,
                )?;

                // Construct and save partial plot if requested
                if let Some(ref plot_path) = args.plot {
                    let plot_stem = plot_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("plot");
                    let plot_dir = plot_path.parent().unwrap_or(std::path::Path::new("."));
                    let partial_plot_name = format!("{}_{}.html", plot_stem, input_stem); // Changed extension to .html
                    let partial_plot_path = plot_dir.join(partial_plot_name);

                    match args.mode {
                        AnalysisMode::RadialMean => {
                            generate_and_save_line_chart(
                                &partial_plot_path,
                                &format!("Radial PSD - {}", input_stem),
                                &axis_labels,
                                &profile,
                                &header,
                            )?;
                        }
                        AnalysisMode::AngularMean => {
                            generate_and_save_rose_chart(
                                &partial_plot_path,
                                &format!("Angular Power - {}", input_stem),
                                &profile,
                                &axis_labels,
                                &header,
                            )?;
                        }
                    };
                }
            }

            let mut count = progress_counter.lock().unwrap();
            *count += 1;
            let term = Term::stdout();
            let _ = term.clear_line();
            print!("\rAnalyzing polar spectra... {:.0}%", (*count as f32 / n_files as f32) * 100.0);
            let _ = std::io::stdout().flush();

            Ok(profile)
        })
        .collect::<Result<_>>()?;

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All blocks analyzed and aggregated.", text::check_icon());

    if args.save_partials {
        println!("{} Partial CSV and plot files saved for each block.", text::check_icon());
    }

    let final_average_psd = aggregate_results(results)?;

    // Get the corresponding axis labels (wavelengths or angles) and header information
    // from the first processed file, applying any specified bounds.
    let (axis_labels, header) = get_axis_info(
        &entries[0],
        &args.mode,
        wavelength_bounds,
        angle_bounds,
        args.detrend,
    )?;

    // --- Detrending Logic (if applicable) ---
    let mut trend_power: Option<Array1<f64>> = None;
    let mut residual_power: Option<Array1<f64>> = None;

    if let Some(order) = args.detrend {
        if matches!(args.mode, AnalysisMode::RadialMean) {
            // 1. Filter data: Keep only finite, positive values for both labels and psd
            let (filtered_labels, filtered_psd): (Vec<_>, Vec<_>) = axis_labels
                .iter()
                .zip(final_average_psd.iter())
                .filter(|(&lbl, &psd)| lbl.is_finite() && lbl > 0.0 && psd.is_finite() && psd > 0.0)
                .map(|(&lbl, &psd)| (lbl, psd))
                .unzip();

            if filtered_labels.len() > order {
                // 2. Log-log transform
                let log_labels =
                    Array1::from_vec(filtered_labels.iter().map(|&v| v.log10()).collect());
                let log_psd = Array1::from_vec(filtered_psd.iter().map(|&v| v.log10()).collect());

                // 3. Fit polynomial
                let coeffs = polyfit_1d(&log_labels, &log_psd, order)?;

                // 4. Calculate trend and residuals
                let mut trend = Array1::zeros(axis_labels.len());
                let mut residuals = Array1::zeros(axis_labels.len());
                trend.fill(f64::NAN);
                residuals.fill(f64::NAN);

                let mut trend_idx = 0;
                for (i, &label) in axis_labels.iter().enumerate() {
                    if trend_idx < filtered_labels.len()
                        && (label - filtered_labels[trend_idx]).abs() < 1e-9
                    {
                        if final_average_psd[i].is_finite() && final_average_psd[i] > 0.0 {
                            let log_label = label.log10();
                            let mut log_trend_val = 0.0;
                            for j in 0..=order {
                                log_trend_val += coeffs[j] * log_label.powi((order - j) as i32);
                            }
                            let trend_val = 10.0_f64.powf(log_trend_val);
                            trend[i] = trend_val;
                            residuals[i] = (final_average_psd[i] / trend_val - 1.0) * 100.0;
                        }
                        trend_idx += 1;
                    }
                }
                trend_power = Some(trend);
                residual_power = Some(residuals);
            } else {
                println!("  {} Not enough data points ({}) to fit a polynomial of order {}.", text::warning("!"), filtered_labels.len(), order);
            }
        }
    }

    // Save the aggregated results to a CSV file.
    save_to_csv(
        &args.output,
        &header,
        &axis_labels,
        &final_average_psd,
        trend_power.as_ref(),
        residual_power.as_ref(),
    )?;

    if let Some(plot_path) = args.plot {
        match args.mode {
            AnalysisMode::RadialMean => {
                // Plot PSD data
                generate_and_save_line_chart(
                    &plot_path,
                    "Radial Power Spectral Density",
                    &axis_labels,
                    &final_average_psd,
                    &header,
                )?;
                println!("{} Radial PSD plot saved to {}", text::check_icon(), plot_path.display());

                // If residuals are calculated, generate a separate plot for them
                if let Some(res_power) = residual_power {
                    let residuals_plot_path = plot_path.with_file_name(format!(
                        "{}_residuals.html",
                        plot_path.file_stem().unwrap().to_str().unwrap()
                    ));
                    generate_and_save_residuals_chart(
                        &residuals_plot_path,
                        "Residual Power",
                        &axis_labels,
                        &res_power,
                        &header,
                    )?;
                    println!("{} Residuals plot saved to {}", text::check_icon(), residuals_plot_path.display());
                }
            }
            AnalysisMode::AngularMean => {
                generate_and_save_rose_chart(
                    &plot_path,
                    "Angular Power Distribution",
                    &final_average_psd,
                    &axis_labels,
                    &header,
                )?;
                println!("{} Angular power plot saved to {}", text::check_icon(), plot_path.display());
            }
        };
    }

    println!("{}", line);
    println!("{}", text::success("Analysis completed successfully."));
    println!("");

    Ok(())
}
