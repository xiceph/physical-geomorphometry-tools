# FFT DEM Analysis Toolkit

## Overview

The **FFT DEM Analysis Toolkit** is a suite of interoperable command-line tools designed for performing advanced **Fast Fourier Transform (FFT)** analysis on **Digital Elevation Models (DEMs)**.

Unlike monolithic software, this toolkit follows the Unix philosophy: it provides small, single-purpose tools that do one thing well. These tools can be chained together to form flexible, reproducible, and scalable analysis pipelines. This approach is particularly suited for scientific researchers and geomorphologists analyzing terrain patterns, land surface structures, and roughness at various scales.

The toolkit addresses common challenges in geospatial FFT analysis, such as:
- **Large Dataset Handling:** Efficient block-based processing to handle DEMs larger than available memory.
- **Edge Effects:** Configurable windowing (tapering) and padding strategies to minimize spectral leakage.
- **Physical Correctness:** 
    - **2D PSD Normalization:** Implements physically-correct normalization ([m^4]) ensuring that the integral of the PSD equals the spatial variance of the terrain, and that PSD values remain independent of the FFT padding width.
    - **Polar Transformation:** Area-weighted polar transformations to ensure conservation of power when converting 2D spectra to 1D radial profiles.

## Quick Start

**No background in Fourier analysis is required to run a standard analysis.**
The typical workflow reduces to three commands and five arguments. All other
parameters have sensible defaults and can be ignored until you need them.

### Step 1 — Compute the Power Spectral Density
```bash
fft-process --input my_dem.tif --output ./results/fft --window-size 512
```
*Processes your DEM in overlapping blocks and computes the 2D power spectrum.
The `--window-size` should be smaller than the shortest side of your input (e.g., 512 or 1024) to ensure multiple regions are sampled. Default overlap (50%) and edge tapering are applied automatically.*

### Step 2 — Convert to Polar Coordinates
```bash
fft-polar --input ./results/fft --output ./results/polar
```
*Reprojects the 2D spectra into wavenumber vs. angle space.
Default binning (36 angular bins, 64 wavenumber bins) is applied automatically.*

### Step 3 — Summarize and Plot
```bash
fft-analyze --input ./results/polar --output summary.csv \
            --mode radial-mean --plot spectrum.html
```
*Produces a CSV radial power spectrum and an interactive HTML plot.*

**Choosing `--mode`** — this is the one conceptual choice required:
| Mode | Use when you want to know… |
|---|---|
| `radial-mean` | How roughness or energy varies with spatial scale (wavelength). Appropriate for most geomorphological studies. |
| `angular-mean` | Whether terrain has a preferred orientation (e.g., fault lineaments, dune crests, glacial striae). |

If in doubt, start with `radial-mean`.

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

## Example Workflow (Advanced)

The Quick Start above covers most use cases. The extended workflow below
demonstrates optional filtering and reconstruction, and exposes additional
parameters for users who need methodological control.

1. **Process the DEM** — 512×512 blocks, 50% overlap (explicit):
```bash
   fft-process --input input_dem.tif --output ./results/fft \
               --window-size 512 --overlap 256
```

2. **Polar Transformation:**
```bash
   fft-polar --input ./results/fft --output ./results/polar
```

3. **Analyze & Plot** — radial mean power spectrum:
```bash
   fft-analyze --input ./results/polar --output summary.csv \
               --mode radial-mean --plot spectrum_plot.html
```

4. **Filter & Reconstruct** *(optional — e.g., to remove high-frequency noise)*:
```bash
   # Keep only wavelengths longer than 50 m
   fft-filter --input ./results/fft --output ./results/filtered \
              --min-wavelength 50

   # Reconstruct the filtered DEM
   fft-inverse --input ./results/filtered --output filtered_dem.tif
```

## License

This project is open-source. Please refer to the [LICENSE](LICENSE) file for details.
