# Physically-based land surface segmentation

## Overview

`Physically-based land surface segmentation` provides a two-step command-line workflow to produce meaningful land surface objects from an elevation raster (DEM - Digital Elevation Model).
The process involves calculating land surface parameters (LSPs) from a DEM and then perform land surface segmentation based on these parameters using Trimble eCognition server. Alternatively, user can perform segmentation in eCognition Developer software on his own desktop machine.

This tools are part of the larger **physical-geomorphometry** project, which focuses on physically based methods for analyzing landforms and land surface dynamics. It builds on the other parts of the project by allowing the creation of meaningful land surface segments based on a generalized DEM raster.

For further details, please refer to the article _Physical geomorphometry for elementary land surface segmentation and digital geomorphological mapping_ by Minár et al. (2024)[^1].

The procedure consists of two steps:
1. **Calculating of the land surface parameters for segmentation**:
   A python script `process_dem.py` takes a DEM raster and generates LSPs optimized for land surface segmentation.
3. **Performing land surface segmentation**:
   A python script `segment_dem.py` takes the LSPs and uses an eCognition Server running in a Docker container to perform land surface segmentation and export the results as a shapefile.

## Prerequisites

Before using this workflow, ensure you have the following software installed and configured:

*   **Python 3**: For running the scripts.
*   **`lsp_calculator`**: A command-line tool for calculating land surface parameters. This tool must be installed and accessible in your system's PATH.
*   **`raster_tool`**: A command-line tool for raster manipulation. This tool must be installed and accessible in your system's PATH.
*   **Docker**: For running the eCognition command-line engine.
*   **Trimble eCognition Server**: You must have a valid eCognition Server license and an instance of the eCognition Command Line Engine (`ecognition-cle`) running in a Docker container.

## Workflow steps

### Step 1: Calculate LSPs for land surface segmentation

**Script:** `process_dem.py`

This script performs the following actions:
-   Calculates a suite of LSPs for land surface segmentation (Sine of Slope, Sine of Aspect, Cosine of Aspect, Normal slope line curvature, Normal contour curvature, Contour torsion, Contour change of normal contour curvature, Slope line change of normal contour curvature,  Slope line change of normal slope line curvature) using `lsp_calculator`.
-   Processes the generated rasters using `raster_tool` to crop boundaries, normalize curvature values using arctan trasnformation, and rescale values to a standard range (0-255).
-   The processed LSP rasters are saved into a `LSP\for_seg` subdirectory at the same location as input DEM.

#### Usage

```
python process_dem.py /path/to/your/input_dem.tif
```

### Step 2: Perform land surface segmentation

**Script:** `segment_dem.py`

This script executes a Docker command that calls the eCognition command-line engine to:
-   Import the LSP rasters using a predefined import template.
-   Run a multiresolution segmentation algorithm based on a predefined ruleset.
-   Export the final segmentation results as a shapefile.

Alternatively, user can perform segmentation in eCognition Developer software on his own desktop machine using a predefined import template `import_dem_seg.xml` and a predefined ruleset `dem_seg.dcp`.

**Important**: This script communicates with a Docker container. All file paths provided as arguments (`--img_dir`, `--results_dir`, `--import_tmpl_file`, `--ruleset`) must be absolute paths *inside the Docker container's filesystem*. You must have your data directories mounted into the container for this to work.

#### Usage

```
python segment_dem.py \
    --img_dir /mnt/img_data/your_project/LSP/for_seg \
    --results_dir /mnt/results/your_project \
    --sp 100 \
    --import_tmpl_file /mnt/import_dem_seg.xml \
    --ruleset /mnt/dem_seg.dcp
```
### Required arguments

| Argument                | Description                                                                                                     |
| :---------------------- | :-------------------------------------------------------------------------------------------------------------- |
| `--img_dir`             | The directory *inside the container* where your LSP rasters are located.                                        |
| `--results_dir`         | The directory *inside the container* where the output shapefile will be saved.                                  |
| `--sp`                  | The Scale Parameter for multiresolution segmentation. A larger value results in larger segments. (Default: 100).|
| `--import_tmpl_file`    | Path *inside the container* to the eCognition import template XML file.                                         | 
| `--ruleset`             | Path *inside the container* to the eCognition ruleset (`.dcp`) file.                                            |

## Platform Compatibility

-   **Linux**: The scripts are developed and tested for a Linux environment.
-   **Windows**: The scripts can be adapted for Windows, provided that:
    -  You have Windows-compatible versions of `lsp_calculator` and `raster_tool` installed and available in your PATH.
    -  You use Windows-style paths when running `process_dem.py`.
    -  Trimble eCognition docker for Windows is installed and configured correctly. Note that the paths passed to `segment_dem.py` must be in the Windows-style format that the container expects.

### Examples

#### 1. Calculate LSPs for segmentation.

This command calculates LSPs for land surface segmentation using DEM, processes the LSPs using raster_tool and saves the output to `LSP\for_seg` subdirectory at the same location as input DEM.

```bash
python3 process_dem.py dem.tif
```
#### 2. Perform land surface segmentation.

This command executes land surface segmentation based on predefined import template and predefined ruleset, with scale parameter 150. It uses Trimble eCognition server via command-line engine in docker container. The final shapefile is saved in directory `/mnt/results/sample_project`.

```bash
python3 segment_dem.py --img_dir /mnt/img_data/sample_project/LSP/for_seg \
--results_dir /mnt/results/sample_project --sp 150 \
--import_tmpl_file /mnt/import_dem_seg.xml --ruleset /mnt/dem_seg.dcp
```
## License

This project is licensed under the MIT License.

[^1]: Minár, J., Drăguţ, L., Evans, I. S., Feciskanin, R., Gallay, M., Jenčo, M., & Popov, A. (2024). Physical geomorphometry for elementary land surface segmentation and digital geomorphological mapping. Earth-Science Reviews, 248, 104631. https://doi.org/10.1016/J.EARSCIREV.2023.104631
