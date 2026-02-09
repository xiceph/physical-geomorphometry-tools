//! # FFT Core Library
//!
//! This library provides the core data structures and processing functions for the
//! FFT DEM Analysis Toolkit. It is designed to be used by various binary crates
//! that implement specific analysis pipelines.
//!
//! The main components are:
//! - `ProcessConfig`: A struct to configure the processing pipeline.
//! - `BlockProcessor`: An iterator that reads a DEM and yields processed blocks.
//! - Various public functions for individual processing steps like `detrend`,
//!   `apply_hann_window`, `compute_fft`, etc.

pub mod tps;
pub mod text;

use anyhow::{Context, Result};
use gdal::raster::Buffer;
use gdal::Dataset;
use ndarray::{s, Array1, Array2, Axis};
use ndarray_linalg::Solve;
use num_complex::Complex;
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

// --- Filesystem Utilities ---

/// Checks if an output directory exists. If it does, it finds a new, unused
/// directory name by appending an index (e.g., `output`, `output.1`, `output.2`).
/// It then creates the directory and returns its path.
///
/// # Arguments
/// * `path` - The desired path for the output directory.
///
/// # Returns
/// A `Result` containing the final, unique path to the output directory.
pub fn prepare_output_dir(path: PathBuf) -> Result<PathBuf> {
    if !path.exists() {
        std::fs::create_dir_all(&path)
            .with_context(|| format!("Failed to create directory: {:?}", path))?;
        return Ok(path);
    }

    let mut index = 1;
    loop {
        let new_name = match path.extension() {
            Some(ext) => {
                let mut new_name = path
                    .file_stem()
                    .context("Failed to get file stem")?
                    .to_str()
                    .context("Failed to convert stem to string")?
                    .to_string();
                new_name.push_str(&format!(".{}", index));
                new_name.push('.');
                new_name.push_str(
                    ext.to_str()
                        .context("Failed to convert extension to string")?,
                );
                new_name
            }
            None => format!(
                "{}.{}",
                path.file_name()
                    .context("Failed to get file name")?
                    .to_str()
                    .context("Failed to convert filename to string")?,
                index
            ),
        };

        let new_path = path.with_file_name(new_name);

        if !new_path.exists() {
            std::fs::create_dir_all(&new_path)
                .with_context(|| format!("Failed to create directory: {:?}", new_path))?;
            println!(
                "{}: Output directory {:?} already exists. Using {:?} instead.\n",
                text::warning("Warning"),
                path,
                new_path
            );
            return Ok(new_path);
        }
        index += 1;
    }
}

// --- Public-Facing Structs ---

/// Configuration for the FFT processing pipeline.
#[derive(Debug, Clone)]
pub struct ProcessConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub window_size: usize,
    pub overlap: usize,
    pub detrend: Option<usize>,
    pub hann_window: Option<usize>,
    // If provided, forces specific total dimensions.
    pub force_padding_size: Option<usize>,
    // If None, calculates nearest "5-Smooth Number" where
    // padding >= min_padding.
    pub min_padding: Option<usize>,
    // Width of the cosine decay.
    pub taper_width: Option<usize>,
    pub jobs: usize,
    pub save_intermediate: bool,
}

/// Holds the result of an FFT analysis for a single block.
#[derive(Debug, Clone)]
pub struct FFTResult {
    pub spectrum: Array2<Complex<f64>>,
    pub power_spectrum: Array2<f64>,
    pub frequencies_x: Array1<f64>,
    pub frequencies_y: Array1<f64>,
    pub row_start: usize,
    pub col_start: usize,
}

/// Holds the result of a polar transformation of the power spectrum.
#[derive(Debug, Clone)]
pub struct PolarSpectrum {
    pub wavelengths: Array1<f64>,
    pub angles: Array1<f64>,
    pub power: Array2<f64>,
}

/// A simple struct to hold the data for a single block.
#[derive(Debug)]
pub struct Block {
    pub data: Array2<f64>,
}

/// Detailed statistics calculated for each processed block.
#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct BlockStatistics {
    pub original_mean: f64,
    pub original_std: f64,
    pub detrended_std: f64,
    pub nodata_percentage: f64,
    pub dc_power_percentage: f64,
    pub total_power: f64,
    pub variance: f64,
    pub parseval_error: f64,
    pub residual_mean: f64,
}

/// Struct for serializing block metadata to a JSON file.
#[derive(Debug, Serialize, Deserialize)]
pub struct BlockMetadata {
    pub block_position: (usize, usize),
    pub original_size: (usize, usize),
    pub padded_size: (usize, usize),
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trend_coeffs: Option<Vec<f64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub geo_transform: Option<[f64; 6]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub wkt: Option<String>,
    pub processing_params: serde_json::Value,
    pub statistics: serde_json::Value,
}

// --- Block Processor Iterator ---

/// An iterator that processes a DEM in blocks.
pub struct BlockProcessor {
    config: Arc<ProcessConfig>,
    input_path: PathBuf,
    blocks_to_process: Vec<(usize, usize)>,
    pixel_size: f64,
    nodata_value: Option<f64>,
    current_block: usize,
    geo_transform: [f64; 6],
    wkt: Option<String>,
    dataset_size: (usize, usize),
}

impl BlockProcessor {
    /// Creates a new BlockProcessor.
    ///
    /// Initializes the processor by checking the input DEM metadata, setting up 
    /// the thread pool, and determining the blocks to be processed.
    ///
    /// # Arguments
    /// * `config` - The `ProcessConfig` containing all processing parameters.
    ///
    /// # Returns
    /// A `Result` containing a new `BlockProcessor` instance.
    pub fn new(config: ProcessConfig) -> Result<Self> {
        if config.jobs > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(config.jobs)
                .build_global()?;
        }

        let ds = Dataset::open(&config.input).context("Failed to open input DEM")?;
        let band = ds.rasterband(1)?;
        let (width, height) = (ds.raster_size().0, ds.raster_size().1);
        let geo_transform = ds.geo_transform()?;
        let wkt = ds.spatial_ref()?.to_wkt().ok();
        let pixel_size = (geo_transform[1].abs() + geo_transform[5].abs()) / 2.0;
        let nodata_value = band.no_data_value();

        println!("{} Input DEM metadata loaded ({}x{}).", text::check_icon(), width, height);

        let step = config.window_size - config.overlap;
        let n_cols = (width - config.window_size) / step + 1;
        let n_rows = (height - config.window_size) / step + 1;
        let mut blocks_to_process = Vec::with_capacity(n_rows * n_cols);
        for r in 0..n_rows {
            for c in 0..n_cols {
                blocks_to_process.push((r * step, c * step));
            }
        }
        println!(
            "Processing {} blocks ({}x{})",
            blocks_to_process.len(),
            n_cols,
            n_rows
        );

        Ok(BlockProcessor {
            input_path: config.input.clone(),
            config: Arc::new(config),
            blocks_to_process,
            pixel_size,
            nodata_value,
            current_block: 0,
            geo_transform,
            wkt,
            dataset_size: (width, height),
        })
    }

    /// Returns the pixel size of the DEM.
    pub fn pixel_size(&self) -> f64 {
        self.pixel_size
    }

    /// Returns a clone of the Arc-wrapped configuration.
    pub fn config(&self) -> Arc<ProcessConfig> {
        self.config.clone()
    }

    /// Returns a reference to the geotransform array.
    pub fn geo_transform(&self) -> &[f64; 6] {
        &self.geo_transform
    }

    /// Returns a reference to the WKT string.
    pub fn wkt(&self) -> Option<&String> {
        self.wkt.as_ref()
    }

    /// Returns a list of all blocks to process (row_start, col_start).
    pub fn get_all_blocks(&self) -> Vec<(usize, usize)> {
        self.blocks_to_process.clone()
    }

    /// Loads a block with "smart" padding: real data where available, TPS extrapolation where not.
    ///
    /// # Arguments
    /// * `row_start` - The starting row of the inner window.
    /// * `col_start` - The starting column of the inner window.
    ///
    /// # Returns
    /// A `Result` containing the padded data and optional trend coefficients.
    pub fn load_smart_padded_block(
        &self,
        row_start: usize,
        col_start: usize,
    ) -> Result<(Array2<f64>, Option<Array1<f64>>)> {
        let config = &self.config;
        let window_size = config.window_size;
        let taper_width = config.taper_width.unwrap_or(0);
        let (cols, rows) = self.dataset_size;

        // 1. Initialize Padded Array
        let padded_rows = window_size + 2 * taper_width;
        let padded_cols = window_size + 2 * taper_width;
        let mut padded_data = Array2::<f64>::zeros((padded_rows, padded_cols));
        let mut accumulated_weights = Array2::<f64>::zeros((padded_rows, padded_cols));

        // 2. Identify Valid Region (Global Coords)
        let pad_r_start = row_start as isize - taper_width as isize;
        let pad_c_start = col_start as isize - taper_width as isize;

        // Open dataset handle locally for this thread.
        let ds = Dataset::open(&self.input_path)?;
        let band = ds.rasterband(1)?;

        // Calculate intersection of padded window and dataset extent
        let global_r_min = pad_r_start.max(0);
        let global_r_max = (pad_r_start + padded_rows as isize).min(rows as isize);
        let global_c_min = pad_c_start.max(0);
        let global_c_max = (pad_c_start + padded_cols as isize).min(cols as isize);

        if global_r_max > global_r_min && global_c_max > global_c_min {
            let read_rows = (global_r_max - global_r_min) as usize;
            let read_cols = (global_c_max - global_c_min) as usize;

            let buffer: Buffer<f64> = band.read_as(
                (global_c_min as isize, global_r_min as isize),
                (read_cols, read_rows),
                (read_cols, read_rows),
                None,
            )?;

            let read_data = Array2::from_shape_vec((read_rows, read_cols), buffer.data().to_vec())?;

            // Map read_data into padded_data
            let local_r_start = (global_r_min - pad_r_start) as usize;
            let local_c_start = (global_c_min - pad_c_start) as usize;

            let mut data_slice = padded_data.slice_mut(s![
                local_r_start..local_r_start + read_rows,
                local_c_start..local_c_start + read_cols
            ]);
            data_slice.assign(&read_data);

            let mut weights_slice = accumulated_weights.slice_mut(s![
                local_r_start..local_r_start + read_rows,
                local_c_start..local_c_start + read_cols
            ]);
            weights_slice.fill(1.0);

            // Handle nodata
            if let Some(ndv) = self.nodata_value {
                data_slice.mapv_inplace(|x| if x == ndv { f64::NAN } else { x });
                // Reset weights where NaN
                for ((r, c), &val) in data_slice.indexed_iter() {
                    if val.is_nan() {
                        weights_slice[[r, c]] = 0.0;
                    }
                }
            }
        }

        // 3. Detrending
        let inner_slice = padded_data.slice(s![
            taper_width..taper_width + window_size,
            taper_width..taper_width + window_size
        ]);
        let mut detrended_inner = inner_slice.to_owned();
        let mut trend_coeffs = None;

        if let Some(order) = config.detrend {
            let coeffs = calculate_trend_coeffs(&detrended_inner, order)?;

            // Remove trend from detrended_inner (for TPS source)
            remove_trend(
                &mut detrended_inner,
                &coeffs,
                order,
                window_size,
                window_size,
                0,
                0,
            )?;

            // Remove trend from the FULL padded_data (only where valid)
            let offset_r = -(taper_width as isize);
            let offset_c = -(taper_width as isize);
            remove_trend(
                &mut padded_data,
                &coeffs,
                order,
                window_size,
                window_size,
                offset_r,
                offset_c,
            )?;

            // Reset invalid pixels to 0.0
            for r in 0..padded_rows {
                for c in 0..padded_cols {
                    if accumulated_weights[[r, c]] == 0.0 {
                        padded_data[[r, c]] = 0.0;
                    }
                }
            }
            trend_coeffs = Some(coeffs);
        }

        // 4. Targeted Extrapolation (Missing Sides)
        let is_top_missing = pad_r_start < 0;
        let is_bottom_missing = (pad_r_start + padded_rows as isize) > rows as isize;
        let is_left_missing = pad_c_start < 0;
        let is_right_missing = (pad_c_start + padded_cols as isize) > cols as isize;

        let mut tps_sides = Vec::new();
        if is_top_missing {
            tps_sides.push(tps::Side::Top);
        }
        if is_bottom_missing {
            tps_sides.push(tps::Side::Bottom);
        }
        if is_left_missing {
            tps_sides.push(tps::Side::Left);
        }
        if is_right_missing {
            tps_sides.push(tps::Side::Right);
        }

        // Identify the bounds of ALL valid data in padded_data (including valid padding).
        // If a side is NOT missing, the valid data extends to the edge of padded_data (or 2T+W).
        // If a side IS missing, valid data starts/ends at the boundary of the inner window (T / T+W).
        let valid_r_start = if is_top_missing { taper_width } else { 0 };
        let valid_r_end = if is_bottom_missing {
            taper_width + window_size
        } else {
            padded_rows
        };
        let valid_c_start = if is_left_missing { taper_width } else { 0 };
        let valid_c_end = if is_right_missing {
            taper_width + window_size
        } else {
            padded_cols
        };

        // Create a view of the maximally valid data block to serve as the source for TPS.
        // This allows TPS to see adjacent valid padding (e.g., Right Padding when processing Top Edge),
        // enabling correct "Edge-like" extrapolation for one-sided corners.
        // We clone it immediately to release the borrow on `padded_data` so we can mutate `padded_data` later.
        let detrended_source = padded_data
            .slice(s![valid_r_start..valid_r_end, valid_c_start..valid_c_end])
            .to_owned();
        let (source_rows, source_cols) = detrended_source.dim();

        for side in tps_sides {
            // Expand into corners ONLY if the adjacent side is also missing.
            // If the adjacent side is valid, our `detrended_source` includes it, so we extrapolate
            // *over* it (as part of the edge) rather than expanding *into* it.
            let expand_corners = match side {
                tps::Side::Top | tps::Side::Bottom => (is_left_missing, is_right_missing),
                tps::Side::Left | tps::Side::Right => (is_top_missing, is_bottom_missing),
            };

            let blocks = tps::process_side_blocks(
                &detrended_source,
                taper_width,
                side,
                source_rows,
                source_cols,
                expand_corners,
                true,
            )?;

            for (block_surf, r_start_rel, c_start_rel, side_enum) in blocks {
                let (b_rows, b_cols) = block_surf.dim();
                let is_horizontal = matches!(side_enum, tps::Side::Top | tps::Side::Bottom);

                // Map coordinates from `detrended_source` frame to `padded_data` frame.
                // detrended_source starts at (valid_r_start, valid_c_start) in padded_data.
                //
                // process_side_blocks returns coordinates relative to the source.
                // However, if expansion occurred (e.g. Left missing), the relative coord might be 0
                // but physically represent a position -T relative to the source start.
                // Wait, process_side_blocks logic:
                // If expanded Left (Side::Top), c_start is 0. Source starts at T (valid_c_start).
                // We want result at 0 (Global).
                // Global = valid_c_start + relative_c - T (if expanded)?
                //
                // Let's check `tps.rs` logic again.
                // Side::Top: `c_start_in_padded = 0` if expanded start.
                // `c_start_in_padded = taper_width` if not expanded start.
                // This "padded" in tps.rs refers to a hypothetical frame where the source is centered?
                // No, tps.rs assumes the source is the "Inner" part of a (N+2T) frame.
                // So it returns coordinates assuming the source starts at (T, T).

                // Here, our "source" `detrended_source` is NOT necessarily (N)x(N). It's variable.
                // But `process_side_blocks` treats it as the "Inner Data".
                // So it returns coordinates relative to a frame where `detrended_source` is at offset (T, T).
                // i.e., Output (0,0) corresponds to Source (-T, -T).

                // So: Global_Pos = (Source_Start_Global - T) + Output_Pos_Relative_to_TPS_Frame.
                // valid_r_start is Source_Start_Global_Row.
                // r_start_rel is Output_Pos_Row.
                //
                // Let's verify:
                // Case 1: Top Edge, Left Missing.
                // valid_c_start = T.
                // expand_corners.0 = true.
                // tps.rs: c_start_in_padded = 0.
                // Global Col = (T - T) + 0 = 0. Correct.
                //
                // Case 2: Top Edge, Left Valid (One-sided corner).
                // valid_c_start = 0.
                // expand_corners.0 = false.
                // tps.rs: c_start_in_padded = taper_width = T.
                // Global Col = (0 - T) + T = 0. Correct.

                let r_global =
                    (valid_r_start as isize - taper_width as isize) + r_start_rel as isize;
                let c_global =
                    (valid_c_start as isize - taper_width as isize) + c_start_rel as isize;

                // Bounds check and Add
                let r_start = r_global.max(0) as usize;
                let c_start = c_global.max(0) as usize;
                // If negative global, we must slice the block.
                let block_r_offset = (r_start as isize - r_global) as usize;
                let block_c_offset = (c_start as isize - c_global) as usize;

                let r_end = (r_global + b_rows as isize).min(padded_rows as isize) as usize;
                let c_end = (c_global + b_cols as isize).min(padded_cols as isize) as usize;

                if r_start >= r_end || c_start >= c_end {
                    continue;
                }

                let len_r = r_end - r_start;
                let len_c = c_end - c_start;

                let block_slice = block_surf.slice(s![
                    block_r_offset..block_r_offset + len_r,
                    block_c_offset..block_c_offset + len_c
                ]);

                let mut weights_2d = Array2::ones((b_rows, b_cols));
                if is_horizontal {
                    let hann = hann_window_1d(b_cols);
                    for (c, &val) in hann.iter().enumerate() {
                        weights_2d.column_mut(c).fill(val);
                    }
                } else {
                    let hann = hann_window_1d(b_rows);
                    for (r, &val) in hann.iter().enumerate() {
                        weights_2d.row_mut(r).fill(val);
                    }
                }
                let weights_slice = weights_2d.slice(s![
                    block_r_offset..block_r_offset + len_r,
                    block_c_offset..block_c_offset + len_c
                ]);

                let mut val_view =
                    padded_data.slice_mut(s![r_start..r_start + len_r, c_start..c_start + len_c]);
                let mut w_view = accumulated_weights
                    .slice_mut(s![r_start..r_start + len_r, c_start..c_start + len_c]);

                for br in 0..len_r {
                    for bc in 0..len_c {
                        let r_global = r_start + br;
                        let c_global = c_start + bc;

                        let mut weight = weights_slice[[br, bc]];

                        let taper_width_f = taper_width as f64;

                        // Top-left corner
                        if is_top_missing
                            && is_left_missing
                            && r_global < taper_width
                            && c_global < taper_width
                        {
                            if r_global == taper_width - 1 && c_global == taper_width - 1 {
                                weight = 1.0;
                            } else {
                                let weight_top = (taper_width_f - 1.0 - r_global as f64).powi(2);
                                let weight_left = (taper_width_f - 1.0 - c_global as f64).powi(2);
                                weight = match side_enum {
                                    tps::Side::Top => weight_top,
                                    tps::Side::Left => weight_left,
                                    _ => 0.0,
                                };
                            }
                        }
                        // Top-right corner
                        else if is_top_missing
                            && is_right_missing
                            && r_global < taper_width
                            && c_global >= padded_cols - taper_width
                        {
                            if r_global == taper_width - 1 && c_global == padded_cols - taper_width
                            {
                                weight = 1.0;
                            } else {
                                let c_rel = c_global - (padded_cols - taper_width);
                                let weight_top = (taper_width_f - 1.0 - r_global as f64).powi(2);
                                let weight_right = (c_rel as f64).powi(2);
                                weight = match side_enum {
                                    tps::Side::Top => weight_top,
                                    tps::Side::Right => weight_right,
                                    _ => 0.0,
                                };
                            }
                        }
                        // Bottom-left corner
                        else if is_bottom_missing
                            && is_left_missing
                            && r_global >= padded_rows - taper_width
                            && c_global < taper_width
                        {
                            if r_global == padded_rows - taper_width && c_global == taper_width - 1
                            {
                                weight = 1.0;
                            } else {
                                let r_rel = r_global - (padded_rows - taper_width);
                                let weight_bottom = (r_rel as f64).powi(2);
                                let weight_left = (taper_width_f - 1.0 - c_global as f64).powi(2);
                                weight = match side_enum {
                                    tps::Side::Left => weight_left,
                                    tps::Side::Bottom => weight_bottom,
                                    _ => 0.0,
                                };
                            }
                        }
                        // Bottom-right corner
                        else if is_bottom_missing
                            && is_right_missing
                            && r_global >= padded_rows - taper_width
                            && c_global >= padded_cols - taper_width
                        {
                            if r_global == padded_rows - taper_width
                                && c_global == padded_cols - taper_width
                            {
                                weight = 1.0;
                            } else {
                                let r_rel = r_global - (padded_rows - taper_width);
                                let c_rel = c_global - (padded_cols - taper_width);
                                let weight_bottom = (r_rel as f64).powi(2);
                                let weight_right = (c_rel as f64).powi(2);
                                weight = match side_enum {
                                    tps::Side::Bottom => weight_bottom,
                                    tps::Side::Right => weight_right,
                                    _ => 0.0,
                                };
                            }
                        }

                        if w_view[[br, bc]] == 0.0 {
                            val_view[[br, bc]] = 0.0;
                        }
                        val_view[[br, bc]] += block_slice[[br, bc]] * weight;
                        w_view[[br, bc]] += weight;
                    }
                }
            }
        }

        // Normalize blended regions
        for r in 0..padded_rows {
            for c in 0..padded_cols {
                let w = accumulated_weights[[r, c]];
                if w > 1e-9 {
                    padded_data[[r, c]] /= w;
                }
            }
        }

        // 5. Final Taper with Zero Mean Enforcement
        // This function calculates the weighted mean of the data (using the taper weights),
        // subtracts it to remove DC bias, and then applies the taper.
        // This ensures the edges decay exactly to 0.0 (no artifact) and the weighted mean is 0.0.
        tps::apply_radial_taper_with_zero_mean(
            &mut padded_data,
            taper_width,
            window_size,
            window_size,
        );

        Ok((padded_data, trend_coeffs))
    }
}

/// Helper function to create a 1D Hann window for blending.
fn hann_window_1d(len: usize) -> Vec<f64> {
    if len <= 1 {
        return vec![1.0; len];
    }
    (0..len)
        .map(|i| 0.5 * (1.0 - (2.0 * std::f64::consts::PI * i as f64 / (len - 1) as f64).cos()))
        .collect()
}

impl Iterator for BlockProcessor {
    type Item = Result<(Block, usize, usize)>; // Returns Block and its coordinates

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_block >= self.blocks_to_process.len() {
            return None;
        }

        let (row_start, col_start) = self.blocks_to_process[self.current_block];
        self.current_block += 1;

        let block_result = load_block_data(
            &self.input_path,
            row_start,
            col_start,
            self.config.window_size,
            self.nodata_value,
        );

        Some(block_result.map(|block| (block, row_start, col_start)))
    }
}

// --- Public Processing Functions ---

/// Extracts a block of data from the DEM file.
///
/// # Arguments
/// * `input_path` - Path to the input DEM file.
/// * `row_start` - The starting row index of the block.
/// * `col_start` - The starting column index of the block.
/// * `window_size` - The size of the square window (block).
/// * `nodata_value` - An optional nodata value to convert to `f64::NAN`.
///
/// # Returns
/// A `Result` containing the extracted `Block` data.
pub fn load_block_data(
    input_path: &Path,
    row_start: usize,
    col_start: usize,
    window_size: usize,
    nodata_value: Option<f64>,
) -> Result<Block> {
    let ds = Dataset::open(input_path)?;
    let band = ds.rasterband(1)?;

    let buffer: Buffer<f64> = band.read_as(
        (col_start as isize, row_start as isize),
        (window_size, window_size),
        (window_size, window_size),
        None,
    )?;

    let mut data = Array2::from_shape_vec((window_size, window_size), buffer.data().to_vec())?;

    // Convert nodata values to NaN for consistent handling.
    if let Some(ndv) = nodata_value {
        data.mapv_inplace(|x| if x == ndv { f64::NAN } else { x });
    }

    Ok(Block { data })
}

/// Checks for NaN (nodata) values in a block.
///
/// # Arguments
/// * `data` - A mutable reference to the `Array2<f64>` block data.
///
/// # Returns
/// Returns `Ok(0)` if no NaN values are found. If any `NaN` values are present,
/// an `Err` is returned with a message that includes the count of `NaN` pixels.
pub fn handle_nodata(data: &mut Array2<f64>) -> Result<usize> {
    // TODO: More advanced NaN filling methods may be implemented here in the future,
    // such as interpolation or inpainting. For now, we simply check for NaNs
    // and expect the user to handle them upstream.
    let nan_count = data.iter().filter(|&&v| v.is_nan()).count();
    if nan_count > 0 {
        anyhow::bail!(
            "Block contains {} nodata (NaN) pixels. Input data must not contain nodata values.",
            nan_count
        );
    }
    Ok(0)
}

/// Calculates the polynomial trend coefficients for a given block of data.
///
/// # Arguments
/// * `data` - The `Array2<f64>` block data.
/// * `order` - The order of the polynomial trend (1 for linear, 2 for quadratic).
///
/// # Returns
/// A `Result` containing the coefficients as an `Array1<f64>`.
pub fn calculate_trend_coeffs(data: &Array2<f64>, order: usize) -> Result<Array1<f64>> {
    let (rows, cols) = data.dim();
    let n_points = rows * cols;
    let z: Vec<f64> = data.iter().cloned().collect();

    let n_coeffs = match order {
        1 => 3, // 1 (constant) + x + y
        2 => 6, // 1 (constant) + x + y + x^2 + y^2 + xy
        _ => anyhow::bail!("Detrend order must be 1 or 2"),
    };
    let mut a = Array2::<f64>::zeros((n_points, n_coeffs));
    for r in 0..rows {
        for c in 0..cols {
            let i = r * cols + c;
            // Normalize x and y coordinates to range from -1 to 1.
            // This normalization is critical for numerical stability when solving 
            // the linear system, especially for higher-order polynomials, as it 
            // prevents the Vandermonde matrix entries from exploding or vanishing.
            let x = (2.0 * c as f64 / (cols - 1) as f64) - 1.0;
            let y = (2.0 * r as f64 / (rows - 1) as f64) - 1.0;

            a[[i, 0]] = 1.0; // Constant term
            if order >= 1 {
                a[[i, 1]] = x;
                a[[i, 2]] = y;
            } // Linear terms
            if order >= 2 {
                a[[i, 3]] = x * x;
                a[[i, 4]] = y * y;
                a[[i, 5]] = x * y;
            } // Quadratic terms
        }
    }

    // Solve the linear system A * coeffs = Z to find the trend coefficients.
    let at = a.t();
    let ata = at.dot(&a);
    let atz = at.dot(&Array1::from_vec(z));
    solve_linear_system(ata, atz)
}

/// Removes a polynomial trend from the data, supporting an offset coordinate system.
///
/// # Arguments
/// * `data` - A mutable reference to the `Array2<f64>` data to detrend.
/// * `coeffs` - The trend coefficients calculated from the reference window.
/// * `order` - The order of the polynomial trend.
/// * `ref_rows` - The number of rows in the reference window.
/// * `ref_cols` - The number of columns in the reference window.
/// * `offset_rows` - The row offset of the `data` relative to the reference window.
/// * `offset_cols` - The column offset of the `data` relative to the reference window.
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn remove_trend(
    data: &mut Array2<f64>,
    coeffs: &Array1<f64>,
    order: usize,
    ref_rows: usize,
    ref_cols: usize,
    offset_rows: isize,
    offset_cols: isize,
) -> Result<()> {
    let (rows, cols) = data.dim();

    for r in 0..rows {
        for c in 0..cols {
            let r_ref = r as isize + offset_rows;
            let c_ref = c as isize + offset_cols;

            let x = (2.0 * c_ref as f64 / (ref_cols - 1) as f64) - 1.0;
            let y = (2.0 * r_ref as f64 / (ref_rows - 1) as f64) - 1.0;

            let mut trend_val = coeffs[0];
            if order >= 1 {
                trend_val += coeffs[1] * x + coeffs[2] * y;
            }
            if order >= 2 {
                trend_val += coeffs[3] * x * x + coeffs[4] * y * y + coeffs[5] * x * y;
            }
            data[[r, c]] -= trend_val;
        }
    }

    Ok(())
}

/// Removes a polynomial trend from the data.
///
/// # Arguments
/// * `data` - A mutable reference to the `Array2<f64>` block data.
/// * `order` - The order of the polynomial trend to remove (1 for linear, 2 for quadratic).
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn detrend(data: &mut Array2<f64>, order: usize) -> Result<Array1<f64>> {
    let coeffs = calculate_trend_coeffs(data, order)?;
    let (rows, cols) = data.dim();
    remove_trend(data, &coeffs, order, rows, cols, 0, 0)?;
    Ok(coeffs)
}

/// Reapplies a polynomial trend to the data.
///
/// # Arguments
/// * `data` - A mutable reference to the `Array2<f64>` block data.
/// * `coeffs` - A reference to the `Array1<f64>` containing the polynomial coefficients.
/// * `order` - The order of the polynomial trend to reapply (1 for linear, 2 for quadratic).
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn reapply_trend(data: &mut Array2<f64>, coeffs: &Array1<f64>, order: usize) -> Result<()> {
    let (rows, cols) = data.dim();

    let n_coeffs = coeffs.len();
    let expected_n_coeffs = match order {
        1 => 3, // 1 (constant) + x + y
        2 => 6, // 1 (constant) + x + y + x^2 + y^2 + xy
        _ => anyhow::bail!("Reapply trend order must be 1 or 2"),
    };

    if n_coeffs != expected_n_coeffs {
        anyhow::bail!(
            "Mismatched number of trend coefficients. Expected {}, got {}",
            expected_n_coeffs,
            n_coeffs
        );
    }

    for r in 0..rows {
        for c in 0..cols {
            // Normalize x and y coordinates to range from -1 to 1.
            let x = (2.0 * c as f64 / (cols - 1) as f64) - 1.0;
            let y = (2.0 * r as f64 / (rows - 1) as f64) - 1.0;

            let mut trend_val = coeffs[0];
            if order >= 1 {
                trend_val += coeffs[1] * x + coeffs[2] * y;
            }
            if order >= 2 {
                trend_val += coeffs[3] * x * x + coeffs[4] * y * y + coeffs[5] * x * y;
            }
            data[[r, c]] += trend_val;
        }
    }

    Ok(())
}

/// Solves a system of linear equations Ax = b using ndarray-linalg.
///
/// # Arguments
/// * `a` - The coefficient matrix `A`.
/// * `b` - The right-hand side vector `b`.
///
/// # Returns
/// A `Result` containing the solution vector `x`.
pub fn solve_linear_system(a: Array2<f64>, b: Array1<f64>) -> Result<Array1<f64>> {
    a.solve_into(b)
        .context("Failed to solve linear system for detrending")
}

/// Applies a 2D cosine taper to the edges of the data.
///
/// This function applies a cosine taper, similar to a Hann window function,
/// to the border regions of the input data, leaving the central part
/// unaffected. This is commonly referred to as a Tukey window or a Hann taper.
/// It helps to mitigate spectral leakage by smoothly transitioning the signal
/// to zero at the boundaries of the tapered region.
///
/// # Arguments
/// * `data` - A mutable reference to the `Array2<f64>` block data.
/// * `taper_width` - The width of the cosine taper applied at each edge.
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn apply_hann_window(data: &mut Array2<f64>, taper_width: usize) -> Result<()> {
    let (rows, cols) = data.dim();
    if taper_width == 0 {
        return Ok(());
    }

    if 2 * taper_width >= rows || 2 * taper_width >= cols {
        anyhow::bail!(
            "Taper width ({}) is too large for window dimensions ({}x{}).",
            taper_width,
            rows,
            cols
        );
    }

    let mut window = Array2::<f64>::ones((rows, cols));
    let flat_top_start_r = taper_width as f64;
    let flat_top_end_r = (rows - 1 - taper_width) as f64;
    let flat_top_start_c = taper_width as f64;
    let flat_top_end_c = (cols - 1 - taper_width) as f64;

    for r_idx in 0..rows {
        for c_idx in 0..cols {
            let r = r_idx as f64;
            let c = c_idx as f64;
            // Calculate distance from the flat-top region.
            let dx = (flat_top_start_c - c).max(0.0).max(c - flat_top_end_c);
            let dy = (flat_top_start_r - r).max(0.0).max(r - flat_top_end_r);
            let dist = if dx <= 0.0 && dy <= 0.0 {
                0.0
            } else if dx <= 0.0 {
                dy
            } else if dy <= 0.0 {
                dx
            } else {
                (dx * dx + dy * dy).sqrt()
            };

            // Apply Hann window function if outside the flat-top.
            if dist > 0.0 {
                let t = (dist / taper_width as f64).min(1.0);
                window[[r_idx, c_idx]] = 0.5 * (1.0 + (std::f64::consts::PI * t).cos());
            }
        }
    }
    *data = data.clone() * window;
    Ok(())
}

/// Checks if a number is a 5-smooth number (only prime factors are 2, 3, 5).
fn is_smooth_number(mut n: usize) -> bool {
    if n == 0 {
        return false;
    }
    for d in [2, 3, 5] {
        while n % d == 0 {
            n /= d;
        }
    }
    n == 1
}

/// Applies zero-padding to the data to achieve a desired final size.
///
/// # Arguments
/// * `data` - The input `Array2<f64>` data.
/// * `min_padding` - The minimum padding to add on each side.
/// * `forced_size` - Optional forced total size for the padded array.
///
/// # Returns
/// A `Result` containing a tuple `(padded_data, row_offset, col_offset)`, where:
/// - `padded_data`: The new, zero-padded `Array2<f64>`.
/// - `row_offset`: The starting row index of the original data within `padded_data`.
/// - `col_offset`: The starting column index of the original data within `padded_data`.
pub fn apply_zero_padding(
    data: &Array2<f64>,
    min_padding: usize,
    forced_size: Option<usize>,
) -> Result<(Array2<f64>, usize, usize)> {
    let (rows, cols) = data.dim();

    // Determine the final padded size.
    let total_size = match forced_size {
        Some(size) => {
            if size < rows || size < cols {
                anyhow::bail!(
                    "Forced size ({}) must be larger than data dimensions ({}x{}).",
                    size,
                    rows,
                    cols
                );
            }
            size
        }
        None => {
            let min_dim = rows.max(cols);
            // The total padding must be at least min_padding.
            // We search for the nearest "5-smooth" number (whose prime factors are 
            // only 2, 3, or 5). FFT algorithms are significantly faster when the 
            // input dimensions are composite numbers with small prime factors.
            let min_total_size = min_dim + 2 * min_padding;
            (min_total_size..)
                .find(|&s| is_smooth_number(s))
                .unwrap_or(min_total_size)
        }
    };

    let mut padded_data = Array2::<f64>::zeros((total_size, total_size));
    let pad_rows_start = (total_size - rows) / 2;
    let pad_cols_start = (total_size - cols) / 2;

    let mut view = padded_data.slice_mut(s![
        pad_rows_start..pad_rows_start + rows,
        pad_cols_start..pad_cols_start + cols
    ]);
    view.assign(data);

    Ok((padded_data, pad_rows_start, pad_cols_start))
}

/// Computes the 2D Fast Fourier Transform of the data and calculates Power Spectral Density (PSD).
///
/// # Arguments
/// * `data` - The input `Array2<f64>` data.
/// * `pixel_size` - The physical size of a pixel in the spatial domain.
/// * `row_start` - The starting row of the block (for metadata).
/// * `col_start` - The starting column of the block (for metadata).
///
/// # Returns
/// A `Result` containing an `FFTResult` struct with the complex spectrum, PSD, and frequencies.
pub fn compute_fft(
    data: &Array2<f64>,
    pixel_size: f64,
    row_start: usize,
    col_start: usize,
) -> Result<FFTResult> {
    let (rows, cols) = data.dim();
    let mut complex_data: Array2<Complex<f64>> = data.mapv(|v| Complex::new(v, 0.0));

    // Perform 2D FFT.
    let mut planner = FftPlanner::new();
    let fft_rows = planner.plan_fft_forward(cols);
    let fft_cols = planner.plan_fft_forward(rows);

    // FFT rows.
    complex_data.axis_iter_mut(Axis(0)).for_each(|mut row| {
        fft_rows.process(row.as_slice_mut().unwrap());
    });

    // FFT columns (after transposing for efficient contiguous access).
    let mut contig_transposed = complex_data.t().as_standard_layout().to_owned();

    contig_transposed
        .axis_iter_mut(Axis(0))
        .for_each(|mut row| {
            fft_cols.process(row.as_slice_mut().unwrap());
        });
    complex_data = contig_transposed.t().to_owned();

    // Calculate Power Spectral Density (PSD).
    // We normalize by the square of the number of elements (N^2). 
    // This choice ensures that the sum of all PSD elements equals the mean of the 
    // squared spatial signal (Mean Power), satisfying Parseval's theorem for this 
    // specific implementation of the Discrete Fourier Transform.
    let n_elements = (rows * cols) as f64;
    let mut power_spectral_density = complex_data.mapv(|c| c.norm_sqr() / n_elements.powi(2));

    // Shift zero-frequency component to the center.
    // This is essential for isotropic analysis (polar transformation) and makes 
    // the resulting power maps intuitive to interpret visually.
    fftshift_2d(&mut complex_data);
    fftshift_2d(&mut power_spectral_density);

    // Calculate and shift frequency bins.
    let mut freqs_x = fftfreq(cols, pixel_size);
    let mut freqs_y = fftfreq(rows, pixel_size);
    fftshift_1d(&mut freqs_x);
    fftshift_1d(&mut freqs_y);

    Ok(FFTResult {
        spectrum: complex_data,
        power_spectrum: power_spectral_density,
        frequencies_x: Array1::from_vec(freqs_x),
        frequencies_y: Array1::from_vec(freqs_y),
        row_start,
        col_start,
    })
}

/// Calculates the frequency bins for an FFT.
///
/// # Arguments
/// * `n` - The number of samples.
/// * `d` - The sample spacing (e.g., pixel size).
///
/// # Returns
/// A `Vec<f64>` containing the frequency bins.
pub fn fftfreq(n: usize, d: f64) -> Vec<f64> {
    let val = 1.0 / (n as f64 * d);
    let n_half = n.div_ceil(2);
    let mut results: Vec<f64> = (0..n_half).map(|i| i as f64).collect();
    results.extend((-(n as isize / 2)..0).map(|i| i as f64));
    results.iter_mut().for_each(|x| *x *= val);
    results
}

/// Shifts the zero-frequency component to the center of a 1D array.
///
/// # Arguments
/// * `array` - A mutable slice of any type `T` that implements `Copy`.
pub fn fftshift_1d<T: Copy>(array: &mut [T]) {
    let n = array.len();
    let n2 = n.div_ceil(2);
    array.rotate_left(n2);
}

/// Shifts the zero-frequency component to the center of a 2D array.
///
/// # Arguments
/// * `array` - A mutable `Array2<T>` of any type `T`.
pub fn fftshift_2d<T>(array: &mut Array2<T>) {
    let (rows, cols) = array.dim();
    let (rows2, cols2) = (rows.div_ceil(2), cols.div_ceil(2));

    // Split the array horizontally
    let (top, bottom) = array.view_mut().split_at(Axis(0), rows2);

    // Split each horizontal part vertically
    let (mut top_left, mut top_right) = top.split_at(Axis(1), cols2);
    let (mut bottom_left, mut bottom_right) = bottom.split_at(Axis(1), cols2);

    // Swap quadrants
    for r in 0..rows2 {
        let max_c = if r < bottom_right.dim().0 {
            bottom_right.dim().1
        } else {
            0
        };
        for c in 0..max_c {
            if r < top_left.dim().0 && c < top_left.dim().1 {
                std::mem::swap(&mut top_left[[r, c]], &mut bottom_right[[r, c]]);
            }
        }
    }

    for r in 0..rows2 {
        let max_c = if r < bottom_left.dim().0 {
            bottom_left.dim().1
        } else {
            0
        };
        for c in 0..max_c {
            if r < top_right.dim().0 && c < top_right.dim().1 {
                std::mem::swap(&mut top_right[[r, c]], &mut bottom_left[[r, c]]);
            }
        }
    }
}

// --- File I/O Functions ---

/// Saves the results of an FFT computation for a single block.
///
/// # Arguments
/// * `result` - The `FFTResult` containing the computed spectrum, PSD, and frequencies.
/// * `pixel_size` - The physical size of a pixel.
/// * `config` - The `ProcessConfig` containing output directory and other parameters.
///
/// # Returns
/// A `Result` indicating success or failure.
#[allow(clippy::too_many_arguments)]
pub fn save_fft_results(
    result: &FFTResult,
    pixel_size: f64,
    config: &ProcessConfig,
    original_size: (usize, usize),
    trend_coeffs: Option<Array1<f64>>,
    geo_transform: Option<[f64; 6]>,
    wkt: Option<String>,
    stats: Option<BlockStatistics>,
) -> Result<()> {
    let output_dir = &config.output;

    // Save power spectral density (log-transformed for dynamic range).
    let psd_path = output_dir.join(format!(
        "fft_psd_block_{}_{}.tif",
        result.row_start, result.col_start
    ));
    let log_psd = result.power_spectrum.mapv(|p| (p + 1e-12).log10());
    save_gdal_raster(&log_psd, &psd_path, None, None, None)?;

    // Save complex spectrum to a binary file.
    let complex_path = output_dir.join(format!(
        "fft_complex_block_{}_{}.bin",
        result.row_start, result.col_start
    ));
    let mut file = File::create(complex_path)?;
    for c in result.spectrum.iter() {
        file.write_all(&c.re.to_le_bytes())?;
        file.write_all(&c.im.to_le_bytes())?;
    }

    // Save metadata to a JSON file.
    let meta_path = output_dir.join(format!(
        "fft_metadata_block_{}_{}.json",
        result.row_start, result.col_start
    ));
    let psd_max = result
        .power_spectrum
        .iter()
        .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
    let psd_mean = result.power_spectrum.mean().unwrap_or(0.0);
    let f_nyquist = 1.0 / (2.0 * pixel_size);

    let mut statistics_map = serde_json::Map::new();
    statistics_map.insert("psd_max".to_string(), json!(psd_max));
    statistics_map.insert("psd_mean".to_string(), json!(psd_mean));
    statistics_map.insert("f_nyquist".to_string(), json!(f_nyquist));

    if let Some(block_stats) = stats {
        if let Ok(stats_val) = serde_json::to_value(block_stats) {
            if let Some(stats_obj) = stats_val.as_object() {
                statistics_map.extend(stats_obj.clone());
            }
        }
    }

    let (padded_rows, padded_cols) = result.spectrum.dim();
    let metadata = BlockMetadata {
        block_position: (result.row_start, result.col_start),
        original_size,
        padded_size: (padded_rows, padded_cols),
        trend_coeffs: trend_coeffs.map(|c| c.to_vec()),
        geo_transform,
        wkt,
        processing_params: json!({
            "detrend_order": config.detrend,
            "min_padding": config.min_padding,
            "taper_width": config.taper_width,
            "force_padding_size": config.force_padding_size,
        }),
        statistics: json!(statistics_map),
    };
    let json_string = serde_json::to_string_pretty(&metadata)?;
    std::fs::write(meta_path, json_string)?;

    Ok(())
}

/// Saves the intermediate, preprocessed block data to a GeoTIFF.
///
/// # Arguments
/// * `data` - The `Array2<f64>` data to save.
/// * `row_start` - The starting row of the block (for filename).
/// * `col_start` - The starting column of the block (for filename).
/// * `config` - The `ProcessConfig` containing the output directory.
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn save_intermediate_block(
    data: &Array2<f64>,
    row_start: usize,
    col_start: usize,
    config: &ProcessConfig,
    block_geo_transform: Option<[f64; 6]>,
    wkt: Option<&str>,
) -> Result<()> {
    let filename = format!("intermediate_block_{}_{}.tif", row_start, col_start);
    let path = config.output.join(filename);
    save_gdal_raster(data, &path, block_geo_transform.as_ref(), wkt, None)
}

/// A generic function to save an `Array2<f64>` to a GeoTIFF file.
///
/// # Arguments
/// * `data` - The `Array2<f64>` data to save.
/// * `path` - The path to the output GeoTIFF file.
///
/// # Returns
/// A `Result` indicating success or failure.
pub fn save_gdal_raster(
    data: &Array2<f64>,
    path: &Path,
    geo_transform: Option<&[f64; 6]>,
    wkt: Option<&str>,
    nodata_value: Option<f64>,
) -> Result<()> {
    let (rows, cols) = data.dim();
    let driver = gdal::DriverManager::get_driver_by_name("GTiff")?;
    let mut ds = driver.create_with_band_type::<f64, _>(path, cols, rows, 1)?;

    if let Some(gt) = geo_transform {
        ds.set_geo_transform(gt)?;
    }

    if let Some(wkt_str) = wkt {
        let srs = gdal::spatial_ref::SpatialRef::from_wkt(wkt_str)?;
        ds.set_spatial_ref(&srs)?;
    }

    let mut band = ds.rasterband(1)?;

    if let Some(ndv) = nodata_value {
        band.set_no_data_value(Some(ndv))?;
    }

    // Collect into a Vec to ensure data is contiguous.
    let data_vec: Vec<f64> = data.iter().cloned().collect();
    let mut buffer = Buffer::new((cols, rows), data_vec);
    band.write((0, 0), (cols, rows), &mut buffer)?;
    Ok(())
}

/// Performs a 1D polynomial fit (polyfit) on the given data.
///
/// This function fits a polynomial of a specified order to a set of 1D data points (x, y)
/// using a least-squares approach. It constructs a Vandermonde matrix for the `x` data
/// and solves the linear system `A * coeffs = y` to find the best-fit coefficients.
///
/// # Arguments
/// * `x` - An `Array1<f64>` containing the x-coordinates of the data points.
/// * `y` - An `Array1<f64>` containing the y-coordinates of the data points.
/// * `order` - The `usize` order of the polynomial to fit.
///
/// # Returns
/// A `Result` containing an `Array1<f64>` with the polynomial coefficients,
/// from the highest power to the constant term.
pub fn polyfit_1d(x: &Array1<f64>, y: &Array1<f64>, order: usize) -> Result<Array1<f64>> {
    let n_points = x.len();
    if n_points != y.len() {
        anyhow::bail!("Input arrays x and y must have the same length.");
    }
    if n_points == 0 {
        anyhow::bail!("Input arrays cannot be empty.");
    }

    let n_coeffs = order + 1;
    let mut a = Array2::<f64>::zeros((n_points, n_coeffs));

    // Create the Vandermonde matrix.
    for i in 0..n_points {
        for j in 0..n_coeffs {
            a[[i, j]] = x[i].powi((n_coeffs - 1 - j) as i32);
        }
    }

    // Solve the linear system A * coeffs = y
    let at = a.t();
    let ata = at.dot(&a);
    let aty = at.dot(y);

    let coeffs = ata
        .solve_into(aty)
        .context("Failed to solve linear system for polyfit.")?;

    Ok(coeffs)
}
