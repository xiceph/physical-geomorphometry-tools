# FFT Inverse

## Overview

`fft-inverse` performs the Inverse Fast Fourier Transform (IFFT) to reconstruct a spatial Digital Elevation Model (DEM) from frequency-domain data. It takes the (potentially filtered) complex blocks from `fft-process` or `fft-filter` and stitches them back together into a single raster image.

## Features

- **Inverse FFT:** Converts complex frequency spectra back to spatial elevation data.
- **Seamless Stitching:** Uses a pyramid-weighting (blending) scheme to merge overlapping blocks, eliminating blocking artifacts.
- **Trend Re-application:** Optionally adds the original detrended surface (planar or quadratic) back to the reconstructed signal.
- **De-padding:** Removes the zero-padding added during the forward process to restore the original block dimensions.

## Usage

### Command Line Arguments

```bash
fft-inverse [OPTIONS] --input <DIR> --output <FILE>
```

#### Required Arguments
- `--input <DIR>`: Directory containing `fft_complex_block_*.bin` and metadata files.
- `--output <FILE>`: Path to the output reconstructed GeoTIFF.

#### Optional Arguments
- `--remove-padding`: Restore the original block dimensions by removing padding. **Strongly Recommended** for final output. Default is `false`.
- `--no-reapply-trend`: Do NOT add the original trend back to the data. Use this if you want to analyze the roughness/residuals only. Default is `false`.
- `--jobs <NUM>`: Number of parallel jobs to run. Defaults to 0 (all available cores).

### Example

Reconstruct a filtered DEM, removing padding and restoring the trend:

```bash
fft-inverse --input ./results/filtered --output reconstructed_dem.tif --remove-padding
```

## Output

- **`reconstructed_dem.tif`**: A single GeoTIFF file covering the extent of the processed blocks.

## Analyzing the Output

### Signal Reconstruction
The output of `fft-inverse` is a realization of the frequency-domain data in physical space.
- **Interpreting Filtered Terrain:** If you have applied a filter (via `fft-filter`), the reconstructed DEM represents the "component" of the landscape within that scale range. For example, a high-pass filtered reconstruction highlights small-scale roughness and texture, while removing the underlying topography.
- **Energy Conservation:** Because the forward and inverse transforms are mathematically linked, the total variance in the reconstructed image (if no filtering was done) should match the variance of the input blocks (minus small windowing losses).

### Blending Artifacts
- **Stitching:** `fft-inverse` uses a weighted blending approach. If you notice "grid-like" patterns in the output, it may be due to insufficient overlap during the `fft-process` step or extreme filtering that has created inconsistent values at block edges. Increasing overlap or using smoother filters (larger taper width) can mitigate this.
- **Trend Discontinuities:** If `--no-reapply-trend` is used, the blocks are stitched as "roughness maps". If the original terrain had sharp trend changes between blocks, the resulting residuals might show subtle jumps at block boundaries. Re-applying the trend usually masks these.
