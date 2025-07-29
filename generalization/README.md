# QEM Generalization – Land Surface Generalization Tool

## Overview

`QEM Generalization` is a tool for generating a **generalized raster** from an elevation raster (DEM - Digital Elevation Model) using the **Quadric Error Metric (QEM)**. This process progressively simplifies the land surface shape while preserving key land surface features. The tool aims to enhance support for digital geomorphological mapping, but it can also be useful for a broader range of geoscientific research.

This tool is part of the larger **physical-geomorphometry** project, which provides physically based methods for analyzing landforms and land surface dynamics. `QEM Generalization` complements other tools in the project by offering efficient raster generalization technique.

## Features

- **Quadric Error Metric-based Generalization**:
  - Uses QEM to simplify elevation rasters.
- **Adjustable Parameters**:
  - Control the number of iterations, resolution reduction factor, and sharpness.
- **Multithreading Support**:
  - Utilize multiple CPU cores for faster computation.
- **Customizable Outputs**:
  - Specify input and output file paths, and fine-tune generalization parameters.

## Theoretical Background

The core algorithm implemented in this project is based on the principles outlined in the article _Advancing raster DEM generalization with a quadric error metric approach_ by Feciskanin & Minár. (2025)[^1]. This paper introduces a unique adaptation of Quadric Error Metrics to operate on raster DEM, which helps preserving essential landforms, even at high levels of generalization. For a comprehensive understanding, please refer to the original publication.

## Usage

### Command-Line Arguments

| Argument              | Short | Description                                   | Default Value                  |
|-----------------------|-------|-----------------------------------------------|--------------------------------|
| `--input-file`        | `-i`  | Input raster file (required).                 | N/A                            |
| `--output-file`       | `-o`  | Output raster file (required).                | N/A                            |
| `--iterations`        | `-n`  | Number of generalization iterations.          | `10`                           |
| `--reduction`         | `-r`  | Resolution reduction factor (≥ 1.0).          | `1.0`                          |
| `--sharpness`         | `-s`  | Sharpness for edge retention (1-9).           | `5`                            |
| `--jobs`              | `-j`  | Number of threads. See Configuration section. | From config or all available   |

### Example Usage

#### Generalize with default settings:
```bash
qem_generalization -i dem.tif -o output.tif
```

#### Specify custom iterations, reduction, and sharpness:
```bash
qem_generalization -i dem.tif -o output.tif -n 50 -r 2.0 -s 7
```

#### Limit the number of threads used:
```bash
qem_generalization -i dem.tif -o output.tif -j 4
```

#### Example with full paths (Windows)

```bash
C:\tools\qem_generalization.exe -i C:\data\dems\sample_dem.tif -o C:\results\generalized_output.tif
```

#### Example with full paths (Linux)

```bash
/home/user/tools/qem_generalization -i /home/user/data/sample_dem.tif -o /home/user/results/generalized_output.tif
```


### Output

The tool generates a single GeoTIFF raster file at the location specified by the `--output-file` (`-o`) argument. The output raster reflects the generalized land surface based on the selected parameters.

## Configuration

The tool can be configured using an `app.config` file located in the same directory as the executable. This file allows you to set parameters that control the tool's behavior.

### Available Parameters

- `max_cores`: Specifies the maximum number of CPU cores the tool is allowed to use. If the `--jobs` argument is provided and its value exceeds `max_cores`, the value from the configuration file will be used instead.

### Example `app.config`

To limit the tool to using a maximum of 4 cores, create an `app.config` file with the following content:
```ini
max_cores = 4
```

## How it Works

1. **Raster Input**:
   - Reads the input raster (DEM) and validates its format.
2. **Generalization**:
   - Applies the Quadric Error Metric (QEM) to simplify the raster:
     - **Iterations**: Controls the depth of simplification.
     - **Reduction Factor**: Defines the extent of resolution reduction.
     - **Sharpness**: Influences the retention of sharp land surface features.
3. **Parallel Processing**:
   - Automatically distributes workload across available CPU cores (configurable).
4. **Output**:
   - Writes the generalized raster to the specified file path.

## Installation

You can either build the tool from source or download a precompiled binary for Windows.

### Option 1: Download Standalone Executable (Windows)

A standalone Windows executable is available for download on the [Releases page](https://github.com/xiceph/physical-geomorphometry-tools/releases). No installation is needed — simply download the `qem_generalization.exe` file and run it from the command line.

### Option 2: Build from Source

To build the tool from source, follow these steps:

### Prerequisites

- **Rust**: Install a Rust development environment to compile the source code.
- **GDAL**: Ensure GDAL is installed to handle raster data.

### Building

Clone the repository and compile the tool:
```bash
git clone https://github.com/xiceph/physical-geomorphometry-tools.git
cd physical-geomorphometry-tools/generalization/
cargo build --release
```

The compiled binary will be located in the `target/release` directory.

### Running the Tool

After building, you can run the tool as demonstrated in the usage examples.

## Future Plans

We plan to develop a **QGIS plugin** for easier integration of the tool into common geospatial workflows. This will allow users to calculate land surface parameters directly within the QGIS environment, enhancing accessibility and usability across platforms.

## License

This project is licensed under the MIT License.

 

[^1]: Feciskanin, R., Minár, J. (2025). Advancing raster DEM generalization with a quadric error metric approach. Computers and Geosciences, 202, 105963. https://doi.org/10.1016/j.cageo.2025.105963.
