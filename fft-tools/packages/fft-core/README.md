# FFT Core Library

## Overview

`fft-core` is the foundational library for the FFT DEM Analysis Toolkit. It contains the shared data structures, algorithms, and I/O utilities used by the binary tools (`fft-process`, `fft-polar`, etc.).

This crate is not intended to be run directly but is a dependency for all other packages in the workspace.

## Key Components

### Data Structures
- **`ProcessConfig`**: A centralized configuration struct that holds parameters for window size, overlap, tapering, and file paths.
- **`BlockMetadata`**: A serializable struct that stores essential information about each processed block (position, size, statistics, geotransform).
- **`FFTResult`**: A container for the outputs of a single block's FFT operation (power spectrum, complex spectrum).

### Algorithms
- **`BlockProcessor`**: The engine responsible for iterating over a large raster dataset in chunks, handling spatial referencing and data loading.
- **Tapering & Padding**: Functions for applying Hann windows (`apply_hann_window`) and zero-padding (`apply_zero_padding`).
- **Detrending**: Least-squares fitting and removal of 1st and 2nd-order polynomial surfaces.
- **FFT Utilities**: Wrappers around `RustFFT`, plus helper functions like `fftshift` (1D and 2D) and frequency frequency generation (`fftfreq`).

### I/O
- **GDAL Integration**: Utilities for reading DEMs and writing GeoTIFFs (`save_gdal_raster`).
- **File Management**: Helpers for preparing output directories and saving intermediate binary files.
