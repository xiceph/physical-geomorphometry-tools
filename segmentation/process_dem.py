
import os
import shutil
import subprocess
import argparse
import glob

def run_command(command):
    """Runs a command and prints its output."""
    print(f"Executing: {' '.join(command)}")
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        print("Stdout:", result.stdout)
        print("Stderr:", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command)}")
        print("Stdout:", e.stdout)
        print("Stderr:", e.stderr)
        raise

def calculate_LSPs(input_raster):
    """
    Calculates Land Surface Parameters (LSPs) for segmentation from input DEM TIF file. Prior generalization of DEM is recommended.
    LSP rasters are then adjusted to achieve better segmentation results. At first, all LSP rasters are clipped by 1px on each side to prevent undefined parameter values. The curvature rasters are transformed to approximate a normal distribution. Finally, the rasters are rescaled to a range of 0-255 for standardization.
    Requires lsp_calculator and raster_tool installed.
    Final LSP rasters are saved to the directory LSP/for_seg in base directory of the input DEM file.
    LSP rasters can be segmented in two ways: 1. in eCognition Developer using the ruleset file dem_seg.dcp (for correct layer names, please use the import template import_dem_seg.xml), 2. in eCognition server using the script segment_dem.py (requires Trimble eCognition server license and eCognition command line engine in docker installed).
    """
    prefix = os.path.basename(input_raster).split('.')[0]
    out_dir = os.path.join(os.path.dirname(input_raster),"LSP")
    final_out_dir = os.path.join(out_dir, "for_seg")

    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(final_out_dir, exist_ok=True)

    # Copy and rename DEM
    dem_dest_path = os.path.join(out_dir, f"{prefix}_dem.tif")
    print(f"Copying {input_raster} to {dem_dest_path}")
    shutil.copy(input_raster, dem_dest_path)

    # Run lsp_calculator
    lsp_out_prefix = os.path.join(out_dir, prefix)
    lsp_command = [
        "lsp_calculator",
        "-i", input_raster,
        "-o", lsp_out_prefix,
        "--for_segmentation"
    ]
    run_command(lsp_command)

    # Process LSPs with raster_tool
    lsp_files = glob.glob(os.path.join(out_dir, "*.tif"))
    for lsp_file in lsp_files:
        basename = os.path.basename(lsp_file)
        out_file = os.path.join(final_out_dir, basename)
        
        command = ["raster_tools", "-i", lsp_file, "-o", out_file, "-c", "1"]
        if "_kn" in basename or "_tc" in basename:
            command.extend(["-n", "-r"])
        else:
            command.append("-r")
        
        run_command(command)

    # Rename files
    rename_map = {
        "sin_aspect.tif": "sinA.tif",
        "cos_aspect.tif": "cosA.tif",
        "sin_slope.tif": "sinS.tif"
    }
    for old_suffix, new_suffix in rename_map.items():
        old_name = f"{prefix}_{old_suffix}"
        new_name = f"{prefix}_{new_suffix}"
        old_path = os.path.join(final_out_dir, old_name)
        new_path = os.path.join(final_out_dir, new_name)
        if os.path.exists(old_path):
            print(f"Renaming {old_path} to {new_path}")
            os.rename(old_path, new_path)

    # return final_out_dir

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a DEM raster for segmentation.")
    parser.add_argument("input_raster", help="Path to the input DEM .tif file.")
    
    args = parser.parse_args()

    # Part 1: Process TIFs
    # The output of this function is the directory to be processed by the next step.
    # However, the docker command in the second part seems to take the parent of the "final" directory.
    # The original script's logic for process_dem_seg uses the directory containing the "final" directory
    # let's get the parent of the final_dir

    #final_dir = calculate_LSPs(args.input_raster)
    calculate_LSPs(args.input_raster)
    #processing_dir = os.path.dirname(final_dir)
    
    print("Processing complete.")

