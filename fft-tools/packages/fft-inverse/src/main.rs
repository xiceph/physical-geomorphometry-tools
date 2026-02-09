use anyhow::{Context, Result};
use clap::Parser;
use console::Term;
use fft_core::{reapply_trend, save_gdal_raster, text, BlockMetadata};
use ndarray::{s, Array1, Array2, Zip};
use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::collections::HashMap;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use walkdir::WalkDir;

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Performs an inverse FFT on processed blocks to reconstruct a Digital Elevation Model (DEM)."
)]
struct Args {
    /// Input directory containing FFT result blocks.
    #[arg(long)]
    input: PathBuf,

    /// Path for the output reconstructed raster file (.tif).
    #[arg(long)]
    output: PathBuf,

    /// If set, padding is removed to restore the original dimensions.
    #[arg(long)]
    remove_padding: bool,

    /// If set, the original trend will NOT be reapplied to the inverse FFT result.
    #[arg(long)]
    no_reapply_trend: bool,

    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    jobs: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for performing inverse FFT and reconstructing DEMs from spectral blocks.\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Inverse Engine"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    println!("{} Configuration:", text::bold("Reconstruction"));
    println!("  {:<20} {}", "Input Directory:", args.input.display());
    println!("  {:<20} {}", "Output File:", args.output.display());
    println!("  {:<20} {}", "Remove Padding:", args.remove_padding);
    println!("  {:<20} {}", "Re-apply Trend:", !args.no_reapply_trend);
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

    // Find all metadata files to determine the scope of the reconstruction.
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
        anyhow::bail!("No metadata files found in the input directory.");
    }

    // Read the geo_transform from the first metadata file to apply to the final stitched image.
    let first_metadata_path = metadata_paths
        .first()
        .context("No metadata files found, cannot determine geo_transform")?;
    let first_metadata_str = fs::read_to_string(first_metadata_path).with_context(|| {
        format!(
            "Failed to read first metadata file: {:?}",
            first_metadata_path
        )
    })?;
    let first_metadata: BlockMetadata =
        serde_json::from_str(&first_metadata_str).with_context(|| {
            format!(
                "Failed to parse first metadata file: {:?}",
                first_metadata_path
            )
        })?;
    let output_geo_transform = first_metadata.geo_transform;
    let output_wkt = first_metadata.wkt.as_deref();

    let progress_counter = Arc::new(Mutex::new(0));

    print!("Processing inverse FFT for blocks...");
    io::stdout().flush().unwrap();

    // Process all blocks in parallel and collect the results.
    let processed_blocks: Vec<Result<(Array2<f64>, usize, usize)>> = metadata_paths
        .par_iter()
        .map(|path| {
            let res = process_block(path, args.remove_padding, args.no_reapply_trend);
            
            let mut count = progress_counter.lock().unwrap();
            *count += 1;
            let term = Term::stdout();
            let _ = term.clear_line();
            print!("\rProcessing inverse FFT for blocks... {:.0}%", (*count as f32 / n_blocks as f32) * 100.0);
            let _ = io::stdout().flush();
            
            res
        })
        .collect();

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All blocks transformed back to spatial domain.", text::check_icon());

    // Handle any errors that occurred during block processing.
    let (successful_blocks, failed_blocks): (Vec<_>, Vec<_>) =
        processed_blocks.into_iter().partition(Result::is_ok);
    if !failed_blocks.is_empty() {
        for error in failed_blocks {
            eprintln!("{} Error processing block: {:?}", text::error("Error"), error.unwrap_err());
        }
        anyhow::bail!("Failed to process one or more blocks.");
    }

    let blocks: Vec<(Array2<f64>, usize, usize)> =
        successful_blocks.into_iter().map(Result::unwrap).collect();

    if blocks.is_empty() {
        anyhow::bail!("No blocks were successfully processed.");
    }

    // Stitch the blocks together.
    print!("Stitching blocks into final image...");
    io::stdout().flush().unwrap();
    let final_image = stitch_blocks(blocks)?;
    
    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} Blocks stitched seamlessly.", text::check_icon());

    let final_geo_transform = if !args.remove_padding {
        if let Some(mut gt) = output_geo_transform {
            let (original_rows, original_cols) = first_metadata.original_size;
            let (padded_rows, padded_cols) = first_metadata.padded_size;

            let pad_rows_start = (padded_rows - original_rows) / 2;
            let pad_cols_start = (padded_cols - original_cols) / 2;

            gt[0] -= pad_cols_start as f64 * gt[1]; // Adjust x-origin
            gt[3] -= pad_rows_start as f64 * gt[5]; // Adjust y-origin

            Some(gt)
        } else {
            output_geo_transform
        }
    } else {
        output_geo_transform
    };

    // Save the final image.
    save_gdal_raster(
        &final_image,
        &args.output,
        final_geo_transform.as_ref(),
        output_wkt,
        None,
    )?;
    println!("{} Reconstructed spatial data saved to {}", text::check_icon(), args.output.display());

    println!("{}", line);
    println!("{}", text::success("Inverse FFT reconstruction completed."));
    println!("");

    Ok(())
}

fn process_block(
    metadata_path: &Path,
    remove_padding: bool,
    no_reapply_trend: bool,
) -> Result<(Array2<f64>, usize, usize)> {
    let metadata_str = fs::read_to_string(metadata_path)
        .with_context(|| format!("Failed to read metadata file: {:?}", metadata_path))?;
    let metadata: BlockMetadata = serde_json::from_str(&metadata_str)
        .with_context(|| format!("Failed to parse metadata file: {:?}", metadata_path))?;

    let (padded_rows, padded_cols) = metadata.padded_size;

    // Construct the complex data file path from the metadata path.
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

    let complex_data_f64: Vec<f64> = buffer
        .chunks_exact(8)
        .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
        .collect();

    let mut complex_data: Array2<Complex<f64>> =
        Array2::from_shape_fn((padded_rows, padded_cols), |(r, c)| {
            let idx = (r * padded_cols + c) * 2;
            Complex::new(complex_data_f64[idx], complex_data_f64[idx + 1])
        });

    // The data was shifted before saving, so we need to shift it back before IFFT.
    fft_core::fftshift_2d(&mut complex_data);

    // Inverse FFT
    let mut planner = FftPlanner::new();
    let ifft_rows = planner.plan_fft_inverse(padded_cols);
    let ifft_cols = planner.plan_fft_inverse(padded_rows);

    // IFFT rows.
    complex_data
        .axis_iter_mut(ndarray::Axis(0))
        .for_each(|mut row| {
            ifft_rows.process(row.as_slice_mut().unwrap());
        });

    // IFFT columns.
    let mut contig_transposed = complex_data.t().as_standard_layout().to_owned();
    contig_transposed
        .axis_iter_mut(ndarray::Axis(0))
        .for_each(|mut row| {
            ifft_cols.process(row.as_slice_mut().unwrap());
        });
    let ifft_result = contig_transposed.t().to_owned();

    // Manually normalize the IFFT result.
    let n_elements = (padded_rows * padded_cols) as f64;
    let real_part = ifft_result.mapv(|c| c.re / n_elements);

    if remove_padding {
        let (original_rows, original_cols) = metadata.original_size;
        let r_start = (padded_rows - original_rows) / 2;
        let c_start = (padded_cols - original_cols) / 2;
        let mut unpadded = real_part
            .slice(s![
                r_start..r_start + original_rows,
                c_start..c_start + original_cols
            ])
            .to_owned();

        // Reapply trend if coefficients are present and not explicitly disabled.
        if !no_reapply_trend {
            if let Some(coeffs_vec) = metadata.trend_coeffs {
                let detrend_order_value = metadata.processing_params["detrend_order"]
                    .as_u64()
                    .context("detrend_order not found or invalid in metadata")?;
                let detrend_order = detrend_order_value as usize;

                let trend_coeffs = Array1::from_vec(coeffs_vec);
                reapply_trend(&mut unpadded, &trend_coeffs, detrend_order)?;
            }
        }

        Ok((
            unpadded,
            metadata.block_position.0,
            metadata.block_position.1,
        ))
    } else {
        Ok((
            real_part,
            metadata.block_position.0,
            metadata.block_position.1,
        ))
    }
}

/// Generates a 2D Hann window for blending.
///
/// The Hann window is a cosine-based window that smoothly tapers to zero
/// at the edges, which provides a smoother transition between blocks
/// than a simple linear (pyramid) weighting. The 2D Hann window is
/// created by multiplying two 1D Hann windows.
///
/// # Arguments
/// * `rows` - The number of rows for the weight mask.
/// * `cols` - The number of columns for the weight mask.
///
/// # Returns
/// A 2D array of `f64` weights.
fn generate_hann_weights(rows: usize, cols: usize) -> Array2<f64> {
    let hann_rows: Array1<f64> = Array1::from_shape_fn(rows, |i| {
        0.5 * (1.0 - (2.0 * std::f64::consts::PI * (i as f64 + 0.5) / rows as f64).cos())
    });
    let hann_cols: Array1<f64> = Array1::from_shape_fn(cols, |i| {
        0.5 * (1.0 - (2.0 * std::f64::consts::PI * (i as f64 + 0.5) / cols as f64).cos())
    });

    // Create 2D Hann window by multiplying the 1D windows.
    let mut hann_2d = Array2::zeros((rows, cols));
    for r in 0..rows {
        for c in 0..cols {
            hann_2d[[r, c]] = hann_rows[r] * hann_cols[c];
        }
    }
    hann_2d
}

/// Stitches processed blocks back into a single raster using weighted averaging.
///
/// This function uses a pyramid weighting scheme to smoothly blend overlapping
/// areas, preventing the blocky artifacts that result from simple averaging.
///
/// # Arguments
/// * `blocks` - A vector of tuples, where each contains a processed data block,
///   its starting row, and its starting column.
///
/// # Returns
/// A `Result` containing the final, blended `Array2<f64>` raster.
fn stitch_blocks(blocks: Vec<(Array2<f64>, usize, usize)>) -> Result<Array2<f64>> {
    if blocks.is_empty() {
        return Ok(Array2::zeros((0, 0)));
    }
    if blocks.len() == 1 {
        return Ok(blocks[0].0.clone());
    }

    // Determine the full size of the final canvas.
    let mut max_r = 0;
    let mut max_c = 0;
    for (block, r_start, c_start) in &blocks {
        let (rows, cols) = block.dim();
        max_r = max_r.max(r_start + rows);
        max_c = max_c.max(c_start + cols);
    }

    let mut canvas = Array2::<f64>::zeros((max_r, max_c));
    let mut total_weights = Array2::<f64>::zeros((max_r, max_c));

    // Cache for weight masks to avoid re-generating for same-sized blocks.
    let mut weight_masks: HashMap<(usize, usize), Array2<f64>> = HashMap::new();

    for (block, r_start, c_start) in blocks {
        let (rows, cols) = block.dim();

        // Get or create the weight mask for the current block's dimensions.
        let weights = weight_masks
            .entry((rows, cols))
            .or_insert_with(|| generate_hann_weights(rows, cols));

        // Get views into the canvas for the current block's area.
        let mut canvas_view =
            canvas.slice_mut(s![r_start..r_start + rows, c_start..c_start + cols]);
        let mut weights_view =
            total_weights.slice_mut(s![r_start..r_start + rows, c_start..c_start + cols]);

        // Use Zip to apply weights and accumulate values.
        Zip::from(&mut canvas_view)
            .and(&mut weights_view)
            .and(&block)
            .and(weights)
            .for_each(|canvas_pixel, total_weight_pixel, &block_pixel, weight| {
                *canvas_pixel += block_pixel * (*weight);
                *total_weight_pixel += *weight;
            });
    }

    // Normalize the canvas by the total weights.
    Zip::from(&mut canvas)
        .and(&total_weights)
        .for_each(|canvas_pixel, &total_weight| {
            if total_weight > 0.0 {
                *canvas_pixel /= total_weight;
            }
        });

    Ok(canvas)
}
