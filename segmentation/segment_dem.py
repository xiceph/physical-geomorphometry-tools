
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

def process_dem_seg(img_dir, results_dir, sp, import_tmpl_file, ruleset):
    """
    Executes multiresolution segmentation of LSPs on the eCognition server and exports the results into shapefile.
    The size of the segments can be influenced by the scale parameter (the larger the value, the larger the objects).
    Requires Trimble eCognition server license and eCognition command line engine in docker installed.
    Ruleset file dem_seg.dcp, import template import_dem_seg.xml and input LSP rasters must be placed in the correct location for the docker container to access them (directory mounted in docker container).
    Paths to input image directory, results directory, import template file and ruleset must be defined relative to docker container filesystem.

    """
    #dir_name = os.path.basename(input_dir)
    #img_dir = os.path.join(container_path, dir_name)
    #output_dir = os.path.join(base_results_dir, dir_name)

    # The original script does mkdir -p, but then comments it out.
    # We will ensure the directory exists.
    # os.makedirs(output_dir, exist_ok=True)

    docker_command = [
        "docker", "exec", "-it", "ecognition-cle",
        "./DIACmdEngine",
        f"param:SP={sp}",
        f"image-dir={img_dir}",
        "import-connector=import_dem_seg",
        f"import-connector-file={import_tmpl_file}",
        f"ruleset={ruleset}",
        f"output-dir={results_dir}/"
    ]

    run_command(docker_command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Physically-based DEM segmentation.")
    parser.add_argument("--img_dir", help="Input directory for image data (in docker container).")
    parser.add_argument("--results_dir", help="Base directory for the output (in docker container).")
    parser.add_argument("--sp", type=int, default=100, help="Scale parameter for multiresolution segmentation.")
    parser.add_argument("--import_tmpl_file", default="/mnt/import_dem_seg.xml", help="Import template XML file (in docker container) for customized import of image layers.")
    parser.add_argument("--ruleset", default="/mnt/dem_seg.dcp", help="Ruleset for physically-based DEM segmentation and export to shapefile.")

    args = parser.parse_args()


    # Part 2: Process with Docker
    # The process_dem_seg.sh script seems to have hardcoded paths, which we now
    # pass as arguments. The `input_dir` for process_dem_seg should be `processing_dir`.
    # Based on the shell script, the docker command constructs paths inside the container
    # so we need to pass the correct directory names.
    #dir_name = os.path.basename(processing_dir)

    # We need to make sure the directory structure matches what the docker command expects.
    # The docker command uses `image-dir=$IMG_DIR/$dir_name`.
    # Let's check how the paths are constructed.
    # The `process_tifs` created `for_seg/<prefix>/final`.
    # The `process_dem_seg` script seems to expect a directory structure like `/mnt/img_data/<dir_name>`.
    # It appears that the output of the first script needs to be in a location accessible
    # to the docker container, under the path specified by --img_dir.

    # For simplicity, this script will just call the docker command with the
    # directory name, assuming the user has placed the processed files in the correct
    # location for the docker container to access them.

    # The `input_dir` for `process_dem_seg` is the directory name that will be appended to `img_dir_base`.
    # In `process_dem_seg.sh`, `dir_name` is `basename "$dir_path"`.
    # `dir_path` is like `$BASE_DIR/*`.
    # In our case, the equivalent of `$dir_path` is `processing_dir`.

    process_dem_seg(
        img_dir=args.img_dir,
        results_dir=args.results_dir,
        sp=args.sp,
        import_tmpl_file=args.import_tmpl_file,
        ruleset=args.ruleset
    )

    print("Processing complete.")

