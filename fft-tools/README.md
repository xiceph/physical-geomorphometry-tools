# FFT DEM Analysis Toolkit

## Overview

The **FFT DEM Analysis Toolkit** is a suite of interoperable command-line tools designed for performing advanced **Fast Fourier Transform (FFT)** analysis on **Digital Elevation Models (DEMs)**.

Unlike monolithic software, this toolkit follows the Unix philosophy: it provides small, single-purpose tools that do one thing well. These tools can be chained together to form flexible, reproducible, and scalable analysis pipelines. This approach is particularly suited for scientific researchers and geomorphologists analyzing terrain patterns, land surface structures, and roughness at various scales.

The toolkit addresses common challenges in geospatial FFT analysis, such as:
- **Large Dataset Handling:** Efficient block-based processing to handle DEMs larger than available memory.
- **Edge Effects:** Configurable windowing (tapering) and padding strategies to minimize spectral leakage.
- **Physical Correctness:** 
    - **2D PSD Normalization:** Implements physically-correct normalization ([m⁴]) ensuring that the integral of the PSD equals the spatial variance of the terrain, and that PSD values remain independent of the FFT padding width.
    - **Polar Transformation:** Area-weighted polar transformations to ensure conservation of power when converting 2D spectra to 1D radial profiles.

## Theoretical Background

The core methodology implemented in this toolkit is described in *A modular toolkit for block-based spatially explicit spectral analysis of digital elevation models* by Feciskanin (2026)[^1]. The paper presents a framework for spatially localized FFT analysis of digital elevation models, including moving-window spectral decomposition, mitigation of spectral leakage through extrapolated tapering, physically consistent power spectral density estimation, anisotropy analysis, spectral filtering, and quantitative comparison of DEM datasets in the frequency domain. For a comprehensive description of the theoretical foundations and methodological design, please refer to the original publication.

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

You can either build the toolkit from source or download precompiled binaries for Windows.

### Option 1: Download Standalone Executables (Windows)

Standalone Windows executables are available for download on the [Releases page](https://github.com/xiceph/physical-geomorphometry-tools/releases). No installation is needed — simply download the desired `.exe` files and run them from the command line.

### Option 2: Build from Source

To build the toolkit from source, follow these steps:

#### Prerequisites

- **Rust**: [Install Rust](https://www.rust-lang.org/tools/install).
- **GDAL**: Geospatial data abstraction library (required for raster I/O).
- **Git**: For cloning the repository.

#### Building on Linux

Ensure you have GDAL development headers installed (e.g., `libgdal-dev` on Ubuntu or `gdal-devel` on Fedora).

```bash
git clone https://github.com/xiceph/physical-geomorphometry-tools.git
cd physical-geomorphometry-tools/fft-tools/
cargo build --release
```

The compiled binaries will be located in the `target/release/` directory.

#### Building on Windows

Building on Windows requires the MSVC toolchain and a static GDAL build. The toolkit uses **Intel MKL** (`intel-mkl-static`) to ensure reliable linking and performance.

1.  **Prerequisites:**
    - **Rust**: [Install Rust](https://www.rust-lang.org/tools/install) (MSVC toolchain).
    - **Visual Studio 2022**: Install the "Desktop development with C++" workload.
    - **vcpkg**: [Install and bootstrap](https://vcpkg.io/en/getting-started.html) vcpkg.

2.  **Install GDAL:**
    ```powershell
    .\vcpkg\vcpkg install gdal:x64-windows-static
    ```

3.  **Setup Environment:**
    Open a **Developer PowerShell for VS 2022** and set the required variables:
    ```powershell
    $env:VCPKG_ROOT = "C:\path\to\vcpkg"
    $env:VCPKGRS_TRIPLET = "x64-windows-static"
    ```

4.  **Build:**
    ```powershell
    cd fft-tools
    cargo build --release
    ```

**Troubleshooting:**
- **Linking errors (`dgetrf_`):** Ensure the `intel-mkl-static` feature is enabled in `Cargo.toml`.
- **VcpkgNotFound:** Verify `$env:VCPKG_ROOT` is set correctly.
- **Compiler errors:** Always use a **Developer Shell** to ensure MSVC tools (like `dumpbin.exe`) are in your PATH.

## Example Workflow (Advanced)

The Quick Start covers standard cases. This extended workflow demonstrates how to exercise precise methodological control over the spectral pipeline. **For a full list of parameters and detailed explanations, refer to the README in each package's directory.**

1. **Process with High-Precision Preprocessing**
   Use 2nd-order detrending and explicit tapered padding (64-pixel taper with a minimum 64-pixel additional zero-padding) to minimize edge artifacts and spectral leakage.
   ```bash
   fft-process --input input_dem.tif --output ./results/fft \
               --window-size 512 --overlap 256 \
               --detrend 2 --taper-type outer --taper 64 --min-pad 64
   ```

2. **High-Resolution Polar Transformation**
   Refine the angular and radial resolution for a more detailed isotropic analysis while ensuring bins remain well-populated.
   ```bash
   fft-polar --input ./results/fft --output ./results/polar \
              --n-angles 45 --n-wavenumbers 80
   ```

3. **Subsetting and Detrending Analysis**
   Focus on specific spatial scales (e.g., 10m to 1000m) and apply log-log detrending to the radial profile to better identify slope breaks.
   ```bash
   fft-analyze --input ./results/polar --output summary.csv \
               --mode radial-mean --plot spectrum.html \
               --wavelength-bounds 10,1000 --detrend 1
   ```

4. **Filter with Precise Transition Bands**
   Apply a band-pass filter with a custom cosine taper (`--taper-width 0.2`) to the transition zones to minimize the Gibbs phenomenon while preserving target frequencies.
   ```bash
   fft-filter --input ./results/fft --output ./results/filtered \
              --min-wavelength 20 --max-wavelength 500 --taper-width 0.2

   # Reconstruct the filtered terrain
   fft-inverse --input ./results/filtered --output filtered_dem.tif
   ```

5. **Quantify Spectral Information Loss**
   Use `fft-compare` to quantify the difference between the original and processed (e.g., filtered or generalized) DEMs, identifying the "effective resolution" where coherence drops.
   ```bash
   fft-compare --input-a ./results/fft --input-b ./results/filtered \
               --output ./results/comparison --plot comparison.html
   ```

## License

This project is open-source. Please refer to the [LICENSE](LICENSE) file for details.



[^1]: Feciskanin, R., 2026. *A modular toolkit for block-based spatially explicit spectral analysis of digital elevation models*. Computers & Geosciences, 215, 106216. https://doi.org/10.1016/j.cageo.2026.106216
