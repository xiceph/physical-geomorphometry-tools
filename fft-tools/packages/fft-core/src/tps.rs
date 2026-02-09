//! Thin-Plate Spline (TPS) extrapolation for outer tapering.
//!
//! This module implements a sophisticated method for extrapolating data into a padded
//! region using Thin-Plate Splines. It is designed to be both performant and accurate,
//! providing a smooth continuation of the data that tapers to zero at the boundaries.
//!
//! The core of the algorithm is a parallelized, block-based approach:
//! 1.  Each edge of the data is divided into overlapping blocks.
//! 2.  For each block, a TPS model is fitted in parallel using a subsampled set of control points
//!     from the original data.
//! 3.  The TPS model is used to extrapolate a smooth surface into the padding area.
//! 4.  Residuals (original data - smooth surface) are calculated and added back with a
//!     taper to preserve texture.
//! 5.  Blocks at the ends of edges extend their extrapolation into corner regions.
//! 6.  All extrapolated blocks (sides and corners) are blended together using a weighted
//!     average (with Hann window weights) to create the final continuous surface.
//! 7.  A final radial taper is applied to the entire border to ensure it smoothly goes to zero.

use anyhow::{anyhow, Result};
use dashmap::DashMap;
use nalgebra::{DMatrix, DVector};
use ndarray::{s, Array2};
use once_cell::sync::Lazy;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::f64::consts::PI;

/// Global thread pool for parallel computations.
/// This is initialized once and reused across all TPS padding operations to avoid the
/// overhead of creating new thread pools repeatedly.
static GLOBAL_THREAD_POOL: Lazy<ThreadPool> = Lazy::new(|| {
    rayon::ThreadPoolBuilder::new()
        .num_threads(std::thread::available_parallelism().map_or(1, |x| x.get()))
        .build()
        .expect("Failed to create Rayon thread pool")
});

/// Cache for TPS LU decompositions to avoid re-computing for identical block geometries.
/// Key: (taper_width, slice_rows, slice_cols, side)
#[derive(Hash, Eq, PartialEq, Debug, Clone, Copy)]
struct TPSKey {
    taper_width: usize,
    rows: usize,
    cols: usize,
    side: Side,
}

static TPS_CACHE: Lazy<DashMap<TPSKey, nalgebra::LU<f64, nalgebra::Dyn, nalgebra::Dyn>>> =
    Lazy::new(DashMap::new);

/// The TPS basis function: U(r) = r^2 * log(r).
/// This function defines the shape of the influence of each control point.
/// It returns 0 if r is 0 to handle the case of a point being evaluated at its own location.
fn tps_basis(r: f64) -> f64 {
    if r == 0.0 {
        0.0
    } else {
        r.powi(2) * r.ln()
    }
}

/// Enum to specify which side of the DEM is being processed.
#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Side {
    Top,
    Bottom,
    Left,
    Right,
}

/// Processes all blocks for a single side of the data in parallel.
///
/// This function divides the specified `side` into overlapping blocks, fits a TPS model
/// for each block, and generates the extrapolated surface.
///
/// # Arguments
/// * `data` - The source data (typically the detrended inner window).
/// * `taper_width` - The depth of the padding and the control point band.
/// * `side` - Which of the four sides to process.
/// * `original_rows`, `original_cols` - Dimensions of the source data.
/// * `expand_corners` - A tuple `(expand_start, expand_end)` indicating whether to extrapolate
///   into the adjacent corner regions. `start` refers to top/left, `end` to bottom/right.
///
/// # Returns
/// A `Result` containing a vector of tuples. Each tuple contains:
/// - `Array2<f64>`: The extrapolated surface for a single block (including corner extensions).
/// - `usize`, `usize`: The `(r, c)` coordinates in the final padded array where this block should be placed.
/// - `Side`: The side this block belongs to.
#[allow(clippy::type_complexity)]
pub fn process_side_blocks(
    data: &Array2<f64>,
    taper_width: usize,
    side: Side,
    original_rows: usize,
    original_cols: usize,
    expand_corners: (bool, bool),
    apply_residuals: bool,
) -> Result<Vec<(Array2<f64>, usize, usize, Side)>> {
    // Blocks are 3T long with 2T stride, creating a T overlap on each side.
    // This provides a stable 1.5 aspect ratio for TPS extrapolation, robustly handling
    // corner extensions.
    let block_len_along_edge = 3 * taper_width;
    let block_stride_along_edge = 2 * taper_width;

    let (max_coord_along_edge, is_horizontal) = match side {
        Side::Top | Side::Bottom => (original_cols, true),
        Side::Left | Side::Right => (original_rows, false),
    };

    // Calculate the number of blocks needed to cover the entire edge.
    let num_blocks = if max_coord_along_edge <= block_len_along_edge {
        1
    } else {
        // Ceiling division to ensure the last part of the edge is covered.
        ((max_coord_along_edge - block_len_along_edge) as f64 / block_stride_along_edge as f64)
            .ceil() as usize
            + 1
    };

    // Process all blocks in parallel using the global thread pool.
    let results = GLOBAL_THREAD_POOL.install(|| {
        (0..num_blocks)
            .into_par_iter()
            .map(|block_idx| -> Result<(Array2<f64>, usize, usize, Side)> {
                // Determine the segment of the original data edge this block will process.
                let start_along = block_idx * block_stride_along_edge;
                let end_along = (start_along + block_len_along_edge).min(max_coord_along_edge);

                // Slice the band of data that will serve as control points for this block's TPS fit.
                // Offsets are the global coordinates of the top-left corner of this slice.
                let (control_slice, r_offset, c_offset) = match side {
                    Side::Top => (
                        data.slice(s![0..taper_width, start_along..end_along]),
                        0,
                        start_along,
                    ),
                    Side::Bottom => (
                        data.slice(s![
                            original_rows - taper_width..original_rows,
                            start_along..end_along
                        ]),
                        original_rows - taper_width,
                        start_along,
                    ),
                    Side::Left => (
                        data.slice(s![start_along..end_along, 0..taper_width]),
                        start_along,
                        0,
                    ),
                    Side::Right => (
                        data.slice(s![
                            start_along..end_along,
                            original_cols - taper_width..original_cols
                        ]),
                        start_along,
                        original_cols - taper_width,
                    ),
                };

                // --- Subsample control points for performance ---
                let mut control_points = Vec::new();
                let mut control_values = Vec::new();
                let mut control_lambdas = Vec::new();

                for ((r_local, c_local), &val) in control_slice.indexed_iter() {
                    // Distance from the inner edge (boundary with the rest of the data).
                    let dist_from_inner_edge = match side {
                        Side::Top => (taper_width - 1) - r_local,
                        Side::Bottom => r_local,
                        Side::Left => (taper_width - 1) - c_local,
                        Side::Right => c_local,
                    };

                    // Dynamic sparsity: For large tapers, we don't need every single pixel
                    // even on the boundary to define a smooth trend.
                    // A base stride reduces the dense matrix size N significantly (N^3 speedup).
                    // This is an engineering trade-off: higher sparsity allows for much larger 
                    // taper widths without hitting the O(N^3) bottleneck of the LU decomposition.
                    let base_stride = (taper_width / 20).max(1);

                    // Increase sampling factor with distance to reduce control point density.
                    // We define the thresholds relative to base_stride to ensure we always
                    // capture points near the edge (e.g. at dist=base_stride) to define the gradient.
                    let sampling_factor = if dist_from_inner_edge < 4 * base_stride {
                        base_stride
                    } else if dist_from_inner_edge < 16 * base_stride {
                        2 * base_stride
                    } else {
                        4 * base_stride
                    };

                    let along_edge_coord = if is_horizontal { c_local } else { r_local };

                    // Add the point if it's on the crucial inner boundary or meets the subsampling criteria.
                    if dist_from_inner_edge == 0
                        || (dist_from_inner_edge % sampling_factor == 0
                            && along_edge_coord % sampling_factor == 0)
                    {
                        // Convert local slice coordinates to global data coordinates.
                        control_points
                            .push(((c_local + c_offset) as f64, (r_local + r_offset) as f64));
                        control_values.push(val);

                        // Tighter constraint on the boundary AND the first gradient row
                        // to strongly lock the trend near the edge.
                        if dist_from_inner_edge <= base_stride {
                            control_lambdas.push(1e-9);
                        } else {
                            control_lambdas.push(1e-6);
                        }
                    }
                }

                if control_points.is_empty() {
                    return Ok((Array2::zeros((0, 0)), 0, 0, side)); // Should not happen
                }

                // Construct cache key for this block's geometry
                let (cs_rows, cs_cols) = control_slice.dim();
                let cache_key = TPSKey {
                    taper_width,
                    rows: cs_rows,
                    cols: cs_cols,
                    side,
                };

                // Fit the TPS model to the subsampled control points.
                let weights = fit_tps(
                    &control_points,
                    &control_values,
                    &control_lambdas,
                    Some(cache_key),
                )?;

                let (extrap_rows_base, extrap_cols_base) = if is_horizontal {
                    (taper_width, end_along - start_along)
                } else {
                    (end_along - start_along, taper_width)
                };

                let is_first_block = block_idx == 0;
                let is_last_block = num_blocks > 0 && block_idx == num_blocks - 1;

                let mut extrap_rows = extrap_rows_base;
                let mut extrap_cols = extrap_cols_base;

                // Determine the starting (r,c) position for this block in the final padded array.
                let (mut r_start_in_padded, mut c_start_in_padded) = match side {
                    Side::Top => (0, taper_width + start_along),
                    Side::Bottom => (original_rows + taper_width, taper_width + start_along),
                    Side::Left => (taper_width + start_along, 0),
                    Side::Right => (taper_width + start_along, original_cols + taper_width),
                };

                // --- Handle Corners by Expanding Evaluation Area ---
                let mut r_min = 0; // Local coordinate range for evaluation points
                let mut r_max = extrap_rows as isize;
                let mut c_min = 0;
                let mut c_max = extrap_cols as isize;

                // If this is a first/last block AND we are allowed to expand into that corner,
                // expand its evaluation area into the adjacent corner
                // and adjust its placement coordinates in the final padded array.
                match side {
                    Side::Top | Side::Bottom => {
                        if is_first_block && expand_corners.0 {
                            c_min = -(taper_width as isize);
                            c_start_in_padded = 0;
                            extrap_cols += taper_width;
                        }
                        if is_last_block && expand_corners.1 {
                            c_max = (end_along - start_along) as isize + taper_width as isize;
                            extrap_cols += taper_width;
                        }
                    }
                    Side::Left | Side::Right => {
                        if is_first_block && expand_corners.0 {
                            r_min = -(taper_width as isize);
                            r_start_in_padded = 0;
                            extrap_rows += taper_width;
                        }
                        if is_last_block && expand_corners.1 {
                            r_max = (end_along - start_along) as isize + taper_width as isize;
                            extrap_rows += taper_width;
                        }
                    }
                }

                // Generate the grid of points where the TPS surface will be evaluated.
                let mut eval_points =
                    Vec::with_capacity((r_max - r_min) as usize * (c_max - c_min) as usize);
                for r_loc in r_min..r_max {
                    for c_loc in c_min..c_max {
                        // Convert local block evaluation coordinates to global coordinates.
                        let (x_abs, y_abs) = match side {
                            Side::Top => (
                                (start_along as isize + c_loc) as f64,
                                (r_loc - taper_width as isize) as f64,
                            ),
                            Side::Bottom => (
                                (start_along as isize + c_loc) as f64,
                                (original_rows as isize + r_loc) as f64,
                            ),
                            Side::Left => (
                                (c_loc - taper_width as isize) as f64,
                                (start_along as isize + r_loc) as f64,
                            ),
                            Side::Right => (
                                (original_cols as isize + c_loc) as f64,
                                (start_along as isize + r_loc) as f64,
                            ),
                        };
                        eval_points.push((x_abs, y_abs));
                    }
                }

                // Evaluate the TPS to get the smooth extrapolated surface.
                let extrapolated_vals = evaluate_tps(&eval_points, &control_points, &weights);
                let mut extrap_surf =
                    Array2::from_shape_vec((extrap_rows, extrap_cols), extrapolated_vals)?;

                if apply_residuals {
                    // --- Calculate and Add Back Tapered Residuals ---
                    // Evaluate the TPS surface on the original (full) control point grid.
                    let mut all_control_slice_points = Vec::with_capacity(control_slice.len());
                    for ((r_local, c_local), _) in control_slice.indexed_iter() {
                        all_control_slice_points
                            .push(((c_local + c_offset) as f64, (r_local + r_offset) as f64));
                    }
                    let tps_on_band_vec =
                        evaluate_tps(&all_control_slice_points, &control_points, &weights);
                    let tps_on_band = Array2::from_shape_vec(control_slice.dim(), tps_on_band_vec)?;

                    // The residual is the high-frequency detail not captured by the smooth TPS surface.
                    let residuals = &control_slice - &tps_on_band;

                    // Mirror and taper these residuals into the padding area to preserve texture.
                    // We only mirror into corners if we expanded into them.
                    let do_mirror_start = is_first_block && expand_corners.0;
                    let do_mirror_end = is_last_block && expand_corners.1;

                    let mirrored_residuals = mirror_and_taper_residuals(
                        &residuals,
                        taper_width,
                        &side,
                        do_mirror_start,
                        do_mirror_end,
                    )?;

                    // Add residuals to the entire extrapolated surface (including corners).
                    // The mirrored_residuals array is now sized to match extrap_surf.
                    extrap_surf += &mirrored_residuals;
                }

                Ok((extrap_surf, r_start_in_padded, c_start_in_padded, side))
            })
            .collect::<Result<Vec<_>>>()
    })?;

    // Filter out any empty results, just in case.
    Ok(results
        .into_iter()
        .filter(|(a, _, _, _)| !a.is_empty())
        .collect())
}

/// Helper to calculate the radial taper weight for a given pixel.
/// Returns a value between 0.0 and 1.0.
fn get_radial_taper_weight(
    r: usize,
    c: usize,
    taper_width: usize,
    original_rows: usize,
    original_cols: usize,
) -> f64 {
    // Convert padded coordinates to coordinates relative to the original data's top-left corner.
    let r_orig = r as isize - taper_width as isize;
    let c_orig = c as isize - taper_width as isize;

    // Check if the point is in the padding area.
    if r_orig < 0
        || r_orig >= original_rows as isize
        || c_orig < 0
        || c_orig >= original_cols as isize
    {
        // Calculate the 1-based distance from the nearest edge of the original data.
        let dist_y = if r_orig < 0 {
            -r_orig
        } else if r_orig >= original_rows as isize {
            r_orig - original_rows as isize + 1
        } else {
            0
        };
        let dist_x = if c_orig < 0 {
            -c_orig
        } else if c_orig >= original_cols as isize {
            c_orig - original_cols as isize + 1
        } else {
            0
        };

        // Calculate the Euclidean distance.
        let dist_radial = ((dist_x.pow(2) + dist_y.pow(2)) as f64).sqrt();

        // Apply cosine taper based on the radial distance.
        if dist_radial > taper_width as f64 {
            0.0
        } else {
            let taper_ratio = dist_radial / taper_width as f64;
            0.5 * (1.0 + (PI * taper_ratio).cos())
        }
    } else {
        // Inside original data
        1.0
    }
}

/// Applies a radial cosine taper to the padded border.
///
/// This function iterates through the entire padded array, identifies points within the
/// border region, calculates their radial distance from the original data, and applies a
/// cosine taper.
pub fn apply_radial_taper(
    padded_data: &mut Array2<f64>,
    taper_width: usize,
    original_rows: usize,
    original_cols: usize,
) {
    let (padded_rows, padded_cols) = padded_data.dim();
    for r in 0..padded_rows {
        for c in 0..padded_cols {
            let w = get_radial_taper_weight(r, c, taper_width, original_rows, original_cols);
            padded_data[[r, c]] *= w;
        }
    }
}

/// Applies a radial cosine taper AND enforces a zero weighted mean.
///
/// This calculates the mean of the data *weighted by the final taper mask*,
/// subtracts that mean from the data, and *then* applies the taper.
///
/// Mathematical result:
/// 1. The weighted sum of the result is exactly 0.0 (DC bias removal).
/// 2. The edges smoothly decay to exactly 0.0 (No step artifact).
///
/// This is critical because any non-zero mean in the tapered padding would 
/// introduce a significant DC component and low-frequency "leakage" into the 
/// FFT, potentially obscuring the real topographic signal.
pub fn apply_radial_taper_with_zero_mean(
    padded_data: &mut Array2<f64>,
    taper_width: usize,
    original_rows: usize,
    original_cols: usize,
) {
    let (padded_rows, padded_cols) = padded_data.dim();
    let mut sum_wx = 0.0;
    let mut sum_w = 0.0;

    // Pass 1: Calculate Weighted Mean
    for r in 0..padded_rows {
        for c in 0..padded_cols {
            let w = get_radial_taper_weight(r, c, taper_width, original_rows, original_cols);
            sum_wx += padded_data[[r, c]] * w;
            sum_w += w;
        }
    }

    // Avoid division by zero (though w sum should be large for any useful block)
    let weighted_mean = if sum_w > 1e-9 { sum_wx / sum_w } else { 0.0 };

    // Pass 2: Subtract and Apply
    for r in 0..padded_rows {
        for c in 0..padded_cols {
            let w = get_radial_taper_weight(r, c, taper_width, original_rows, original_cols);
            // Subtract mean first, then taper.
            // Edges: w=0 => result=0.
            // Center: w=1 => result = val - mean.
            padded_data[[r, c]] = (padded_data[[r, c]] - weighted_mean) * w;
        }
    }
}

/// Mirrors residuals into the padding area and applies a cosine taper.
///
/// The residuals represent the high-frequency texture of the data not captured by the
/// smooth TPS fit. Mirroring them preserves this texture in the padded region.
///
/// # Arguments
/// * `residuals` - The `original_data - tps_surface` for a given block.
/// * `taper_width` - The width of the padding.
/// * `side` - The side being processed.
/// * `expand_start` - Whether to mirror into the start corner.
/// * `expand_end` - Whether to mirror into the end corner.
///
/// # Returns
/// A `Result` containing the array of tapered, mirrored residuals.
fn mirror_and_taper_residuals(
    residuals: &Array2<f64>,
    taper_width: usize,
    side: &Side,
    expand_start: bool,
    expand_end: bool,
) -> Result<Array2<f64>> {
    let (res_rows, res_cols) = residuals.dim();

    // Determine output dimensions
    let (out_rows, out_cols) = match side {
        Side::Top | Side::Bottom => {
            let mut w = res_cols;
            if expand_start {
                w += taper_width;
            }
            if expand_end {
                w += taper_width;
            }
            (taper_width, w)
        }
        Side::Left | Side::Right => {
            let mut h = res_rows;
            if expand_start {
                h += taper_width;
            }
            if expand_end {
                h += taper_width;
            }
            (h, taper_width)
        }
    };

    let mut padded_residuals = Array2::zeros((out_rows, out_cols));

    for r_pad in 0..out_rows {
        for c_pad in 0..out_cols {
            // Coordinate mapping logic:
            // We need to map (r_pad, c_pad) in the Padded Residual frame
            // to (r_source, c_source) in the Input Residual frame.
            // AND determine the distance for tapering.

            // 1. Identify "Along" coordinate offset caused by expansion.
            let (r_eff, c_eff) = match side {
                Side::Top | Side::Bottom => {
                    let c_shifted = if expand_start {
                        c_pad as isize - taper_width as isize
                    } else {
                        c_pad as isize
                    };
                    (r_pad as isize, c_shifted)
                }
                Side::Left | Side::Right => {
                    let r_shifted = if expand_start {
                        r_pad as isize - taper_width as isize
                    } else {
                        r_pad as isize
                    };
                    (r_shifted, c_pad as isize)
                }
            };

            // 2. Identify if we are in the main body or a corner extension.
            // Main body: 0 <= along_coord < res_len
            // Start corner: along_coord < 0
            // End corner: along_coord >= res_len

            let (along_idx, _across_idx, res_len) = match side {
                Side::Top | Side::Bottom => (c_eff, r_eff, res_cols as isize),
                Side::Left | Side::Right => (r_eff, c_eff, res_rows as isize),
            };

            let is_start_corner = along_idx < 0;
            let is_end_corner = along_idx >= res_len;

            // If in a corner but expansion wasn't requested, we theoretically shouldn't be here
            // due to loop bounds, but strictly speaking it's 0.0.
            if (is_start_corner && !expand_start) || (is_end_corner && !expand_end) {
                continue;
            }

            // 3. Mirroring Logic
            // Across: mirror as usual.
            // Along: if in corner, mirror the nearest edge of the residual block.
            // For corner extension, we can mirror the `along` coordinate too, effectively creating a
            // "corner mirror" (diagonal flip).

            let (mirror_across, dist_across) = match side {
                Side::Top => (taper_width - 1 - r_pad, (taper_width - 1 - r_pad) as f64),
                Side::Bottom => (res_rows - 1 - r_pad, r_pad as f64),
                Side::Left => (taper_width - 1 - c_pad, (taper_width - 1 - c_pad) as f64),
                Side::Right => (res_cols - 1 - c_pad, c_pad as f64),
            };

            // Calculate 'mirror_along' and 'dist_along' for corners
            let (mirror_along, dist_along) = if is_start_corner {
                // Mirror the start: -1 -> 0, -2 -> 1 ... i.e. -along - 1
                (-along_idx - 1, (-along_idx) as f64)
            } else if is_end_corner {
                // Mirror the end: len -> len - 1, len + 1 -> len - 2 ...
                // diff = along - len. target = len - 1 - diff = len - 1 - (along - len) = 2*len - 1 - along
                (
                    2 * res_len - 1 - along_idx,
                    (along_idx - res_len + 1) as f64,
                )
            } else {
                (along_idx, 0.0)
            };

            // Clamp coordinates to be safe (though logic should hold)
            let m_r = if matches!(side, Side::Top | Side::Bottom) {
                mirror_across
            } else {
                mirror_along as usize
            };
            let m_c = if matches!(side, Side::Top | Side::Bottom) {
                mirror_along as usize
            } else {
                mirror_across
            };

            // Bounds check for safety (esp with float/isize conversions)
            if m_r >= res_rows || m_c >= res_cols {
                continue;
            }

            let mirrored_val = residuals[[m_r, m_c]];

            // Taper: max distance (L-infinity) or Euclidean?
            // L-infinity (box) is better for filling square corners without circular artifacts.
            let dist = dist_across.max(dist_along);

            let taper_ratio = dist / taper_width as f64;
            let taper = if taper_ratio > 1.0 {
                0.0
            } else {
                0.5 * (1.0 + (std::f64::consts::PI * taper_ratio).cos())
            };

            padded_residuals[[r_pad, c_pad]] = mirrored_val * taper;
        }
    }
    Ok(padded_residuals)
}

/// Fits a Thin-Plate Spline (TPS) to a set of control points.
///
/// This function solves the linear system of equations to find the weights (`w_i`)
/// and affine coefficients (`a_1`, `a_x`, `a_y`) for the TPS model.
/// The system is represented in matrix form as:
/// | K  P | | w | = | v |
/// | P.T 0 | | a |   | 0 |
/// where K is the matrix of basis function evaluations, P is the polynomial term matrix,
/// w are the weights, a are the affine coefficients, and v are the control point values.
///
/// # Arguments
/// * `points` - A vector of `(x, y)` coordinates for the control points.
/// * `values` - A vector of `z` values at the control points.
/// * `lambdas` - A vector of regularization parameters for each point.
///   Smaller values enforce tighter interpolation.
/// * `cache_key` - Optional key to identify the block geometry for caching the LU decomposition.
///
/// # Returns
/// A `Result` containing the vector of TPS weights and coefficients `[w..., a_1, a_x, a_y]`.
fn fit_tps(
    points: &[(f64, f64)],
    values: &[f64],
    lambdas: &[f64],
    cache_key: Option<TPSKey>,
) -> Result<DVector<f64>> {
    let n = points.len();
    if n == 0 {
        return Err(anyhow!("Cannot fit TPS with zero points."));
    }
    if lambdas.len() != n {
        return Err(anyhow!(
            "Mismatch between number of points and regularization parameters."
        ));
    }
    let m_size = n + 3; // n weights + 3 affine coefficients

    // --- Build the y vector ---
    let mut y_vector = DVector::<f64>::zeros(m_size);
    // The first n values are the control point values.
    for i in 0..n {
        y_vector[i] = values[i];
    }
    // The last 3 values are zero, corresponding to the polynomial constraints.

    // Check cache if key is provided
    if let Some(key) = cache_key {
        if let Some(lu) = TPS_CACHE.get(&key) {
            let weights = lu.solve(&y_vector).ok_or_else(|| {
                anyhow!("Failed to solve TPS linear system (cached). The matrix may be singular.")
            })?;
            return Ok(weights);
        }
    }

    // --- Build the L matrix ---
    let mut l_matrix = DMatrix::<f64>::zeros(m_size, m_size);

    // Build K matrix (top-left, n x n)
    for r in 0..n {
        for c in 0..n {
            let p1 = points[r];
            let p2 = points[c];
            let dist = ((p1.0 - p2.0).powi(2) + (p1.1 - p2.1).powi(2)).sqrt();
            l_matrix[(r, c)] = tps_basis(dist);
        }
    }
    // Add variable regularization to the diagonal of K.
    for i in 0..n {
        l_matrix[(i, i)] += lambdas[i];
    }

    // Build P matrix (top-right, n x 3)
    for i in 0..n {
        l_matrix[(i, n)] = 1.0;
        l_matrix[(i, n + 1)] = points[i].0; // x
        l_matrix[(i, n + 2)] = points[i].1; // y
    }

    // Build P.T matrix (bottom-left, 3 x n)
    for i in 0..n {
        l_matrix[(n, i)] = 1.0;
        l_matrix[(n + 1, i)] = points[i].0;
        l_matrix[(n + 2, i)] = points[i].1;
    }
    // The bottom-right 3x3 block remains zero.

    // --- Solve the linear system Lw = y ---
    let lu = l_matrix.lu();

    // Store in cache if key provided
    if let Some(key) = cache_key {
        TPS_CACHE.insert(key, lu.clone());
    }

    let weights = lu
        .solve(&y_vector)
        .ok_or_else(|| anyhow!("Failed to solve TPS linear system. The matrix may be singular."))?;

    Ok(weights)
}

/// Evaluates a fitted Thin-Plate Spline at new locations.
///
/// Uses the fitted weights and coefficients to calculate the TPS surface value
/// at a new set of evaluation points.
///
/// # Arguments
/// * `eval_points` - A vector of `(x, y)` coordinates to evaluate the spline at.
/// * `control_points` - The original control points used for fitting.
/// * `weights` - The fitted weights from `fit_tps`, a vector of `[w..., a_1, a_x, a_y]`.
///
/// # Returns
/// A `Vec<f64>` of the evaluated spline values at `eval_points`.
fn evaluate_tps(
    eval_points: &[(f64, f64)],
    control_points: &[(f64, f64)],
    weights: &DVector<f64>,
) -> Vec<f64> {
    let n = control_points.len();
    // Extract the affine coefficients from the end of the weights vector.
    let c_a1 = weights[n];
    let c_ax = weights[n + 1];
    let c_ay = weights[n + 2];

    eval_points
        .iter()
        .map(|&(ex, ey)| {
            // Start with the affine part of the model.
            let mut val = c_a1 + c_ax * ex + c_ay * ey;
            // Add the contribution from each control point's radial basis function.
            for i in 0..n {
                let p_control = control_points[i];
                let dist = ((ex - p_control.0).powi(2) + (ey - p_control.1).powi(2)).sqrt();
                val += weights[i] * tps_basis(dist);
            }
            val
        })
        .collect()
}
