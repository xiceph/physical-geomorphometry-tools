# Raster Tools Pipeline

## Overview

`Raster Tools` is a command-line pipeline utility for standardized preprocessing of raster data within the **physical geomorphometry workflow**. It is designed to support reproducible and methodologically consistent raster transformations commonly required before further terrain analysis, segmentation, or visualization.

The tool implements a sequential processing pipeline that operates on a single raster input and applies a predefined set of transformations: cropping, normalization, and rescaling, in a fixed and transparent order. Each step is optional and can be combined with others, allowing users to construct concise preprocessing workflows directly from the command line. 

This tool is part of the larger **physical-geomorphometry** project, which focuses on physically based methods for analyzing landforms and land surface dynamics.



## Features

- **Cropping**: Symmetrically crop a specified number of pixels from all sides of a raster.
- **Normalization**: Apply an arctan transformation to normalize pixel values, optimized by minimizing kurtosis.
- **Rescaling**: Rescale pixel values to a standard range (e.g., 0-255).
- **Pipeline Architecture**: Chain multiple operations in a single command.
- **Efficient**: Built with performance in mind, using parallel processing where possible.


## Usage

The tool is designed as a pipeline. You provide an input file and an output file, and specify the processing steps to apply.

```
raster-tools.exe -i <input.tif> -o <output.tif> [OPTIONS]
```

### Arguments

| Short | Long            | Description                                                              |
| :---- | :-------------- | :----------------------------------------------------------------------- |
| `-i`  | `--input-file`  | **Required.** Specify the input raster file path.                        |
| `-o`  | `--output-file` | **Required.** Specify the output raster file path.                       |
| `-c`  | `--crop`        | Crop N pixels from all sides. Expects an integer value (e.g., `-c 10`).  |
| `-n`  | `--normalize`   | Normalize the raster using the arctan transformation method.             |
| `-r`  | `--rescale`     | Rescale raster values to a standard range.                               |
|       | `--keep-temps`  | Keep intermediate files (`_cropped.tif`, `_normalized.tif`, `_k.csv`).   |
| `-h`  | `--help`        | Print help information.                                                  |
| `-V`  | `--version`     | Print version information.                                               |

### Pipeline Steps

The operations are executed in the following order:

1.  **Crop** (`--crop`)
2.  **Normalize** (`--normalize`)
3.  **Rescale** (`--rescale`)

The output of one step becomes the input for the next.

#### Normalization Details

The normalization process applies an arctangent transformation (`atan(k · value)`), where the scaling parameter *k* is selected iteratively so that the **excess kurtosis of the transformed pixel-value distribution approaches zero**. This heuristic aims to reduce heavy-tailed distributions caused by extreme values while avoiding degenerate over-compression. This approach is inspired by the methodology described in _Transformation (normalization) of slope gradient and surface curvatures, automated for statistical analyses from DEMs_ by Csillik, Evans and Drăguţ. (2015)[^1].

If `--keep-temps` is enabled, a CSV file named `<basename>_k.csv` will be created, logging the `k` values and corresponding kurtosis during the optimization process.

### Intermediate Files

By default, the tool creates temporary files for each step and cleans them up at the end, leaving only the final output.

If you want to inspect the result of each step, use the `--keep-temps` flag. This will save the intermediate files in the same directory as the output file with suffixes like `_cropped.tif` and `_normalized.tif`.

### Examples

#### 1. Crop and Normalize a Raster

This command crops 5 pixels from the border of `input.tif`, normalizes the result, and saves the final output to `processed.tif`.

```sh
./target/release/raster-tools -i C:\data\input.tif -o C:\data\processed.tif --crop 5 --normalize
```

#### 2. Full Pipeline with Intermediate Files

This command runs all three steps (crop, normalize, rescale) and keeps all the intermediate files for inspection.

```sh
./target/release/raster-tools `
  --input-file "input.tif" `
  --output-file "final_output.tif" `
  --crop 10 `
  --normalize `
  --rescale `
  --keep-temps
```

This will produce:
- `input_cropped.tif`
- `input_normalized.tif`
- `input_k.csv`
- `final_output.tif` (a copy of `input_rescaled.tif`)


## Installation

You can either build the tool from source or download a precompiled binary for Windows.

### Option 1: Download Standalone Executable (Windows)

A standalone Windows executable is available for download on the [Releases page](https://github.com/xiceph/physical-geomorphometry-tools/releases). No installation is needed — simply download the `raster_tools.exe` file and run it from the command line.

### Option 2: Build from Source

To build the tool from source, follow these steps:

### Prerequisites

- **Rust**: Install a Rust development environment to compile the source code.
- **GDAL**: Ensure GDAL is installed to handle raster data.

### Building

Clone the repository and compile the tool:
```bash
git clone https://github.com/xiceph/physical-geomorphometry-tools.git
cd physical-geomorphometry-tools/segmentation/raster-tools/
cargo build --release
```

The compiled binary will be located in the `target/release` directory.

### Running the Tool

After building, you can run the tool as demonstrated in the usage examples.



## License

This project is licensed under the MIT License.


[^1]: Csillik, O., Evans, I.S., Drăguţ, L. (2015). Transformation (normalization) of slope gradient and surface curvatures, automated for statistical analyses from DEMs. Geomorphology, Volume 232, Pages 65-77, ISSN 0169-555X,. https://doi.org/10.1016/j.geomorph.2014.12.038. 
