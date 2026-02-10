# FFT DEM Analysis Toolkit

## Overview

The **FFT DEM Analysis Toolkit** is a suite of interoperable command-line tools designed for performing advanced **Fast Fourier Transform (FFT)** analysis on **Digital Elevation Models (DEMs)**.

Unlike monolithic software, this toolkit follows the Unix philosophy: it provides small, single-purpose tools that do one thing well. These tools can be chained together to form flexible, reproducible, and scalable analysis pipelines. This approach is particularly suited for scientific researchers and geomorphologists analyzing terrain patterns, land surface structures, and roughness at various scales.

The toolkit addresses common challenges in geospatial FFT analysis, such as:
- **Large Dataset Handling:** Efficient block-based processing to handle DEMs larger than available memory.
- **Edge Effects:** Configurable windowing (tapering) and padding strategies to minimize spectral leakage.
- **Physical Correctness:** Area-weighted polar transformations to ensure conservation of power when converting 2D spectra to 1D radial profiles.

## Tools

The project is organized as a Cargo workspace with the following tools:

- [**fft-process**](./packages/fft-process) – The entry point. Decomposes a DEM into overlapping blocks, applies detrending and windowing, and computes the 2D FFT (Power Spectral Density).
- [**fft-polar**](./packages/fft-polar) – Transforms the 2D Cartesian PSDs produced by `fft-process` into a polar representation (Angle vs. Wavenumber) using a Jacobian-weighted method.
- [**fft-analyze**](./packages/fft-analyze) – Performs statistical analysis on the polar spectra (e.g., radial mean, angular mean) and generates plots and summaries.
- [**fft-compare**](./packages/fft-compare) – Analyzes spectral differences and coherence between two datasets (e.g., reference vs. generalized DEM).
- [**fft-filter**](./packages/fft-filter) – Applies frequency-domain filters (Low-pass, High-pass, Band-pass) to the complex FFT data.
- [**fft-inverse**](./packages/fft-inverse) – Reconstructs a DEM from the processed (and potentially filtered) FFT blocks, seamlessly stitching them back together.
- [**fft-core**](./packages/fft-core) – The shared library containing the core algorithms and data structures.

## Installation

### Prerequisites
- **Rust:** You will need a Rust toolchain installed (stable channel). [Install Rust](https://www.rust-lang.org/tools/install).
- **GDAL:** The tools depend on the GDAL library for reading and writing geospatial raster files. Ensure GDAL development headers are installed on your system.

### Building from Source
To build all tools in the workspace with optimizations enabled:

```bash
cargo build --release
```

The compiled binaries will be available in `target/release/`.

## Example Workflow

A typical analysis pipeline might look like this:

1.  **Process the DEM:** Compute the PSD for 512x512 blocks with 50% overlap.
    ```bash
    ./target/release/fft-process --input input_dem.tif --output ./results/fft --window-size 512 --overlap 256
    ```

2.  **Polar Transformation:** Convert the Cartesian PSDs to polar coordinates.
    ```bash
    ./target/release/fft-polar --input ./results/fft --output ./results/polar
    ```

3.  **Analyze & Plot:** Generate a radial power spectrum summary and plot.
    ```bash
    ./target/release/fft-analyze --input ./results/polar --output summary.csv --plot spectrum_plot.html --mode radial-mean
    ```

4.  **Filter & Reconstruct (Optional):** Remove high-frequency noise and reconstruct the DEM.
    ```bash
    # Apply a low-pass filter (keep wavelengths > 50m)
    ./target/release/fft-filter --input ./results/fft --output ./results/filtered --min-wavelength 50

    # Inverse FFT to get the filtered DEM
    ./target/release/fft-inverse --input ./results/filtered --output filtered_dem.tif
    ```

## License

This project is open-source. Please refer to the [LICENSE](LICENSE) file for details.
