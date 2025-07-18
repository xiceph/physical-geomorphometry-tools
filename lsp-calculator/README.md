# Land Surface Parameters Calculator

## Overview

`Land Surface Parameters Calculator` is a tool for calculating **Land Surface Parameters (LSPs)** from an elevation raster (DEM - Digital Elevation Model). These parameters are essential in geospatial analysis, particularly in terrain modeling, hydrology, and geomorphology.

This tool is part of the larger **physical-geomorphometry** project, which focuses on physically based methods for analyzing landforms and land surface dynamics. It integrates seamlessly with other modules in the project, contributing to the comprehensive calculation of geomorphometric characteristics.

The tool leverages polynomial approximations to derive various local LSPs, such as slope, aspect, curvatures, and changes of curvatures from raster data. Users can select specific parameters to compute or opt to calculate all supported parameters at once.

## Features

- **Land Surface Parameter Calculations**:

  - **First-Order Parameters**:
    - Slope
    - Aspect
    - Sine of Slope
    - Sine of Aspect
    - Cosine of Aspect

  - **Basic Trio of Curvatures**:
    - Normal slope line (profile) curvature
    - Normal contour (tangential) curvature
    - Contour torsion

  - **Subforms of Basic Trio**:
    - Second slope line derivative
    - Slope line torsion
    - Second contour derivative
    - Projected contour curvature
    - Projected slope line curvature
    - Contour change of sine slope

  - **Other Gravity-Specific Curvatures**:
    - Difference curvature
    - Total accumulation curvature
    - Total ring curvature
    - Horizontal excess curvature
    - Vertical excess curvature

  - **Principal Curvatures**:
    - Maximal curvature
    - Minimal curvature

  - **Other Gravity-Invariant Curvatures**:
    - Gaussian curvature
    - Elevation laplacian
    - Unsphericity curvature
    - Mean curvature
    - Casorati curvature

  - **Changes of Curvatures**:
    - Contour change of normal contour curvature
    - Slope line change of normal contour curvature
    - Slope line change of normal slope line curvature

  _For more detailed information on each parameter, please refer to the [LSPs description page](./PARAMETERS.md)._

- **Partial Derivatives Calculations**:
  - Can output partial derivatives (and interpolated elevation) instead of LSPs.

- **Input and Output**:
  - Supports input from a GeoTIFF raster file.
  - Outputs parameter-specific GeoTIFF files based on selected computations.

- **Polynomial Approximation**:
  - Employs a least squares regression approach using 3rd or 4th-degree polynomial surfaces for LSPs estimation.

- **Parallel Processing**:
  - Utilizes multiple CPU cores to speed up computations.

## Usage

### Command-Line Arguments

```bash
lsp_calculator [OPTIONS]
```

### Required Arguments

- `-i, --input-file <FILE>`: Path to the input raster file (DEM).
- `-o, --output-prefix <PREFIX>`: Prefix for the output GeoTIFF files.

### Optional Arguments

- `-d, --degree <degree>`: Specify the polynomial degree (3 or 4). Default is `3`.
- `-j, --jobs <jobs>`: Specify the number of threads to use. If omitted, all available processors are used.
- `-h, --help`: Print help.
- `-V, --version`: Print version.

### Batch Selection of Parameters

- `--all`: Compute all Land Surface Parameters offered by the tool.

Or select a subset of parameters to compute (one or more options):

- `--first`: First-order Land Surface Parameters
- `--second`: Second-order Land Surface Parameters
- `--third`: Third-order Land Surface Parameters
- `--for_segmentation`: Land Surface Parameters for Land Surface Segmentation (Sine of Slope, Sine of Aspect, Cosine of Aspect, Normal slope line curvature, Normal contour curvature, Contour torsion, Contour change of normal contour curvature, Slope line change of normal contour curvature,  Slope line change of normal slope line curvature)

### Individual Selection of Calculated Parameters

You can choose one or more specific parameters to calculate using the following options:

- `--slope`: Slope
- `--aspect`: Aspect
- `--sin_slope`: Sine of Slope
- `--sin_aspect`: Sine of Aspect
- `--cos_aspect`: Cosine of Aspect
- `--kns`: Normal slope line (profile) curvature
- `--knc`: Normal contour (tangential) curvature
- `--tc`: Contour torsion
- `--zss`: Second slope line derivative
- `--ts`: Slope line torsion
- `--zcc`: Second contour derivative
- `--kpc`: Projected contour curvature
- `--kps`: Projected slope line curvature
- `--sin_sc`: Contour change of sine slope
- `--kd`: Difference curvature
- `--ka`: Total accumulation curvature
- `--kr`: Total ring curvature
- `--khe`: Horizontal excess curvature
- `--kve`: Vertical excess curvature
- `--kmax`: Maximal curvature
- `--kmin`: Minimal curvature
- `--k`: Gaussian curvature
- `--el`: Elevation laplacian
- `--ku`: Unsphericity curvature
- `--kmean`: Mean curvature
- `--kc`: Casorati curvature
- `--kncc`: Contour change of normal contour curvature
- `--kncs`: Slope line change of normal contour curvature
- `--knss`: Slope line change of normal slope line curvature

### Output of Partial Derivatives Only

- `-p, --partials`

When this flag is provided, the tool does **not** compute or output any Land Surface Parameters (LSPs), regardless of other parameter selection options.

Instead, it calculates and outputs the following partial derivatives derived from the local polynomial surface fitted to each raster cell (up to the third order), along with the interpolated elevation value _z_ from the fitted polynomial. These outputs are useful for custom calculations of land surface characteristics.

**Outputs include:**

- `z_fit` – Interpolated elevation _z_ from fitted polynomial
- `zx` – First-order derivative $\frac{\partial z}{\partial x}$
- `zy` – First-order derivative $\frac{\partial z}{\partial y}$
- `zxx` – Second-order derivative $\frac{\partial^2 z}{\partial x^2}$
- `zxy` – Second-order mixed derivative $\frac{\partial^2 z}{\partial x \partial y}$
- `zyy` – Second-order derivative $\frac{\partial^2 z}{\partial y^2}$
- `zxxx` – Third-order derivative $\frac{\partial^3 z}{\partial x^3}$
- `zxxy` – Third-order mixed derivative $\frac{\partial^3 z}{\partial x^2 \partial y}$
- `zxyy` – Third-order mixed derivative $\frac{\partial^3 z}{\partial x \partial y^2}$
- `zyyy` – Third-order derivative $\frac{\partial^3 z}{\partial y^3}$

> **Note:** This mode overrides any other parameter selection flags (`--all`, `--slope`, etc.). Only the partial derivatives and interpolated elevation will be computed and saved.

### Example Usage

#### Calculate all available parameters

```bash
lsp_calculator -i dem.tif -o output -a
```

#### Calculate specific parameters (e.g., slope and aspect)

```bash
lsp_calculator -i dem.tif -o output --slope --aspect
```

#### Calculate partial derivatives instead of LSPs

```bash
lsp_calculator -i dem.tif -o output --partials
```

#### Example with full paths (Windows)

```bash
C:\tools\lsp_calculator.exe -i C:\data\dems\sample_dem.tif -o C:\results\sample_output --all
```

#### Example with full paths (Linux)

```bash
/home/user/tools/lsp_calculator -i /home/user/data/sample_dem.tif -o /home/user/results/sample_output --all
```

### Output

The output consists of separate GeoTIFF files for each calculated parameter. Each file is named using the full output path followed by a suffix indicating the parameter. The `.tif` extension is automatically added. The format is:

```
<output-path>_<parameter>.tif
```

For example, using:

```bash
lsp_calculator -i dem.tif -o output --slope --aspect
```

Produces:
- `output_slope.tif`
- `output_aspect.tif`

Or using:

```bash
-o C:\results\sample_output
```

Produces:
- `C:\results\sample_output_slope.tif`
- `C:\results\sample_output_aspect.tif`

Or on Linux:

```bash
-o /home/user/results/sample_output
```

Produces:
- `/home/user/results/sample_output_slope.tif`
- `/home/user/results/sample_output_aspect.tif`


## How it Works

1. **Raster Input**:
   - The tool reads the input raster (DEM) file and its metadata, including geotransform and projection.

2. **Polynomial Approximation**:
   - A 5x5 neighborhood of elevation points is used for each raster cell.
   - The least squares method is employed to fit a polynomial surface to these points.
   - The degree of the polynomial can be configured to either 3 or 4.

3. **Land Surface Parameter Calculation**:
   - Based on the polynomial surface, partial derivatives are calculated, and from these, various land surface parameters (LSPs) are determined.
   - Calculations are performed in parallel using multiple CPU cores to improve performance.

4. **Output**:
   - The results for each calculated parameter are written to separate GeoTIFF files.

## Installation

You can either build the tool from source or download a precompiled binary for Windows.

### Option 1: Download Standalone Executable (Windows)

A standalone Windows executable is available for download on the [Releases page](https://github.com/xiceph/physical-geomorphometry-tools/releases). No installation is needed — simply download the `lsp_calculator.exe` file and run it from the command line.

### Option 2: Build from Source

To build the tool from source, follow these steps:

#### Prerequisites

- **GDAL**: A geospatial data abstraction library required for handling raster files.
- **Rust**: The tool is written in Rust, so you'll need a Rust development environment.

#### Building

Clone the repository, navigate to the tool, and compile it:

```bash
git clone https://github.com/xiceph/physical-geomorphometry-tools.git
cd physical-geomorphometry-tools/lsp-calculator/
cargo build --release
```
The compiled binary will be located in the `target/release` directory.

### Running the Tool

After building, you can run the tool using the command as shown in the usage section.

## Future Plans

We intend to add more local LSPs as we identify them to be of interest.

We plan to develop a **QGIS plugin** for easier integration of the tool into common geospatial workflows. This will allow users to calculate land surface parameters directly within the QGIS environment, enhancing accessibility and usability across platforms.

## License

This project is licensed under the MIT License.
