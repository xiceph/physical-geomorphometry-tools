use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use console::Term;
use fft_core::{
    apply_hann_window, apply_zero_padding, compute_fft, handle_nodata, prepare_output_dir,
    save_fft_results, save_intermediate_block, text, BlockProcessor, BlockStatistics, ProcessConfig,
};
use gdal::Dataset;
use rayon::prelude::*;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

#[derive(ValueEnum, Clone, Debug, PartialEq)]
enum TaperType {
    Inner,
    Outer,
}

/// Command-line arguments for the fft-process tool.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "A command-line tool to perform FFT analysis on Digital Elevation Models (DEMs)."
)]
struct Args {
    /// Path to the input Digital Elevation Model (DEM) file.
    #[arg(long)]
    input: PathBuf,
    /// Path to the output directory where results will be saved.
    #[arg(long)]
    output: PathBuf,
    /// Size of the square window (block) to process, in pixels. If omitted, defaults to the smaller dimension of the input raster.
    #[arg(long)]
    window_size: Option<usize>,
    /// Overlap between adjacent windows, in pixels.
    #[arg(long)]
    overlap: usize,
    /// Optional: Apply detrending. Specify order (1 or 2). If no order, defaults to 1.
    #[arg(long, num_args(0..=1))]
    detrend: Option<Option<usize>>,
    /// Type of taper to apply.
    #[arg(long, value_enum, default_value_t = TaperType::Outer)]
    taper_type: TaperType,
    /// Optional: Width of the cosine taper for padding.
    #[arg(long, num_args(0..=1))]
    taper: Option<Option<usize>>,
    /// Optional: Apply gradient-based tapered padding. Specify minimum padding.
    #[arg(long, num_args(0..=1), default_value_if("min_pad", "", Some("0")))]
    min_pad: Option<usize>,
    /// Optional: Force a specific total padded dimension.
    #[arg(long)]
    force_padding_size: Option<usize>,
    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    jobs: usize,
    /// If set, intermediate preprocessed block data will be saved.
    #[arg(long)]
    save_intermediate: bool,
}

/// Main entry point for the fft-process tool.
///
/// Parses command-line arguments, configures the processing pipeline,
/// and processes DEM blocks in parallel.
fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for performing FFT analysis on Digital Elevation Models (DEMs).\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Processing Engine"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    // Parse and validate detrending configuration.
    let detrend_config = match args.detrend {
        Some(Some(order)) if order == 1 || order == 2 => Some(order),
        Some(Some(order)) => anyhow::bail!("Detrend order must be 1 or 2, but received {}", order),
        Some(None) => Some(1), // Default detrend order
        None => None,
    };

    let ds = Dataset::open(&args.input).context("Failed to open input DEM")?;
    let (dem_width, dem_height) = (ds.raster_size().0, ds.raster_size().1);

    // Determine the actual window size, defaulting to the smaller DEM dimension if not provided,
    // and clamping to the smaller DEM dimension if provided but larger.
    let final_window_size = match args.window_size {
        Some(ws) => ws.min(dem_width).min(dem_height),
        None => dem_width.min(dem_height),
    };

    let taper_width_config = match args.taper {
        Some(Some(width)) => Some(width),
        Some(None) => Some(final_window_size / 10), // Default to 1/10 of window_size
        None => None,
    };

    // Create a ProcessConfig struct from parsed arguments.
    let config = ProcessConfig {
        input: args.input.clone(),   // Clone to move into config, args.input is used again later
        output: args.output.clone(), // Clone to move into config
        window_size: final_window_size,
        overlap: args.overlap,
        detrend: detrend_config,
        hann_window: if args.taper_type == TaperType::Inner {
            taper_width_config
        } else {
            None
        },
        min_padding: args.min_pad,
        taper_width: if args.taper_type == TaperType::Outer {
            taper_width_config
        } else {
            None
        },
        force_padding_size: args.force_padding_size,
        jobs: args.jobs,
        save_intermediate: args.save_intermediate,
    };

    println!("{} Configuration:", text::bold("Processing"));
    println!("  {:<20} {}", "Input DEM:", config.input.display());
    println!("  {:<20} {}", "Output Directory:", config.output.display());
    println!("  {:<20} {}", "Window Size:", config.window_size);
    println!("  {:<20} {}", "Overlap:", config.overlap);
    if let Some(order) = config.detrend {
        println!("  {:<20} {}", "Detrend Order:", order);
    }
    println!("  {:<20} {:?}", "Taper Type:", args.taper_type);
    if let Some(width) = taper_width_config {
        println!("  {:<20} {}", "Taper Width:", width);
    }
    if let Some(min_pad) = config.min_padding {
        println!("  {:<20} {}", "Min Padding:", min_pad);
    }
    if let Some(force_size) = config.force_padding_size {
        println!("  {:<20} {}", "Forced Padding Size:", force_size);
    }
    println!(
        "  {:<20} {}",
        "Parallel Jobs:",
        if config.jobs == 0 {
            "all available cores".to_string()
        } else {
            config.jobs.to_string()
        }
    );
    println!("{}\n", dline);

    let mut final_config = config.clone();
    final_config.output = prepare_output_dir(config.output)?;

    // Initialize the BlockProcessor, which handles DEM loading and block iteration.
    let processor = BlockProcessor::new(final_config)?;
    let pixel_size = processor.pixel_size();
    let config_arc = processor.config(); // Get an Arc to the config for parallel processing.
    let geo_transform_arc = Arc::new(*processor.geo_transform());
    let wkt_arc = Arc::new(processor.wkt().cloned());

    // Get all block coordinates to process
    let blocks = processor.get_all_blocks();
    let n_blocks = blocks.len();
    let processor_arc = Arc::new(processor);

    let progress_counter = Arc::new(Mutex::new(0));
    let validation_results = Arc::new(Mutex::new(Vec::with_capacity(n_blocks)));

    print!("Processing blocks...");
    io::stdout().flush().unwrap();

    // Process each block in parallel.
    blocks.par_iter().for_each(|&(row_start, col_start)| {
        let result = (|| -> Result<()> {
            let mut stats = BlockStatistics::default();

            // Load Smart Block: Handles real data fetching, detrending, and TPS extrapolation for gaps.
            let (mut block_data, trend_coeffs_option) =
                processor_arc.load_smart_padded_block(row_start, col_start)?;

            stats.original_mean = block_data.mean().unwrap_or(0.0); // Approximation as it's already detrended
            stats.original_std = block_data.std(0.0);

            // Handle nodata values: check for NaNs and abort if any are found.
            // Note: load_smart_padded_block copies valid NaNs if they exist in source.
            let nan_count = handle_nodata(&mut block_data)?;
            stats.nodata_percentage = (nan_count as f64 / block_data.len() as f64) * 100.0;

            stats.residual_mean = block_data.mean().unwrap_or(0.0);
            stats.detrended_std = block_data.std(0.0);

            // Determine dimensions
            let original_size = (config_arc.window_size, config_arc.window_size);

            // Apply tapering if Inner (Outer is handled by load_smart_padded_block)
            let tapered_data = match args.taper_type {
                TaperType::Inner => {
                    if let Some(taper_width) = config_arc.hann_window {
                        apply_hann_window(&mut block_data, taper_width)?;
                    }
                    block_data.clone()
                }
                TaperType::Outer => block_data.clone(),
            };

            // Apply zero-padding
            // If TaperType::Outer, tapered_data is already (N+2T)x(N+2T).
            // apply_zero_padding adds *additional* padding if min_pad > 0 or force_padding_size > dim.
            let (final_data, pad_rows, pad_cols) = apply_zero_padding(
                &tapered_data,
                config_arc.min_padding.unwrap_or(0),
                config_arc.force_padding_size,
            )?;
            stats.variance = final_data.var(0.0);

            // Correctly calculate the geotransform for the intermediate block.
            let outer_taper_offset = if args.taper_type == TaperType::Outer {
                config_arc.taper_width.unwrap_or(0)
            } else {
                0
            };

            // Total pixel offset from the DEM origin to the top-left of the final_data array.
            let total_col_offset =
                col_start as isize - outer_taper_offset as isize - pad_cols as isize;
            let total_row_offset =
                row_start as isize - outer_taper_offset as isize - pad_rows as isize;

            let mut block_geo_transform = *geo_transform_arc;
            block_geo_transform[0] =
                geo_transform_arc[0] + total_col_offset as f64 * geo_transform_arc[1];
            block_geo_transform[3] =
                geo_transform_arc[3] + total_row_offset as f64 * geo_transform_arc[5];

            // Save intermediate preprocessed data if enabled.
            if config_arc.save_intermediate {
                save_intermediate_block(
                    &final_data,
                    row_start,
                    col_start,
                    &config_arc,
                    Some(block_geo_transform),
                    wkt_arc.as_deref(),
                )?;
            }

            // Compute the 2D FFT.
            let fft_result = compute_fft(&final_data, pixel_size, row_start, col_start)?;

            // --- Statistics and Validation ---
            stats.total_power = fft_result.power_spectrum.sum();

            let spatial_mean_square = final_data.mapv(|v| v.powi(2)).mean().unwrap_or(0.0);

            if spatial_mean_square > 1e-9 {
                stats.parseval_error =
                    (stats.total_power - spatial_mean_square).abs() / spatial_mean_square;
            }

            let (rows, cols) = fft_result.power_spectrum.dim();
            let dc_power = fft_result.power_spectrum[[rows / 2, cols / 2]];
            if stats.total_power > 1e-9 {
                stats.dc_power_percentage = (dc_power / stats.total_power) * 100.0;
            }

            // Save the FFT results (PSD, complex spectrum, metadata).
            save_fft_results(
                &fft_result,
                pixel_size,
                &config_arc,
                original_size,
                trend_coeffs_option,
                Some(*geo_transform_arc),
                (*wkt_arc).clone(),
                Some(stats.clone()),
            )?;

            // Collect stats for final report
            validation_results.lock().unwrap().push((row_start, col_start, stats));

            Ok(())
        })();

        if let Err(e) = result {
            eprintln!("\n{} Failed to process block at ({}, {}): {}", text::error("Error"), row_start, col_start, e);
        }

        let mut count = progress_counter.lock().unwrap();
        *count += 1;
        let term = Term::stdout();
        let _ = term.clear_line();
        print!("\rProcessing blocks... {:.0}%", (*count as f32 / n_blocks as f32) * 100.0);
        let _ = io::stdout().flush();
    });

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All blocks processed.", text::check_icon());

    // --- Final Validation Summary ---
    println!("\n{} Scientific Validation Summary:", text::bold("Final"));
    println!("{}", line);

    let stats_vec = validation_results.lock().unwrap();
    
    const PARSEVAL_TOLERANCE: f64 = 1e-6;
    const RESIDUAL_MEAN_TOLERANCE: f64 = 1e-9;
    const DC_POWER_TOLERANCE: f64 = 1.0; // 1%

    let mut parseval_fails = 0;
    let mut residual_fails = 0;
    let mut dc_high = 0;

    for (r, c, s) in stats_vec.iter() {
        if s.parseval_error >= PARSEVAL_TOLERANCE {
            parseval_fails += 1;
            println!("  {} Block {}-{} failed Parseval check (error: {:.2e})", text::warning("!"), r, c, s.parseval_error);
        }
        if s.residual_mean.abs() >= RESIDUAL_MEAN_TOLERANCE {
            residual_fails += 1;
        }
        if s.dc_power_percentage >= DC_POWER_TOLERANCE {
            dc_high += 1;
        }
    }

    if parseval_fails == 0 {
        println!("  {} All blocks passed Parseval energy conservation validation.", text::check_icon());
    } else {
        println!("  {} {} blocks failed Parseval validation.", text::error("✗"), parseval_fails);
    }

    if residual_fails == 0 {
        println!("  {} All blocks successfully detrended (residual mean < 1e-9).", text::check_icon());
    } else {
        println!("  {} {} blocks have high residual mean after detrending.", text::warning("!"), residual_fails);
    }

    if dc_high == 0 {
        println!("  {} DC component power is low across all blocks.", text::check_icon());
    } else {
        println!("  {} {} blocks have significant DC power (> 1%).", text::warning("!"), dc_high);
    }

    println!("{}", line);
    println!("{}", text::success("FFT processing completed."));
    println!("");

    Ok(())
}
