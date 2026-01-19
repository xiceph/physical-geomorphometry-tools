
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
    Performs physically-based elementary land surface segmentation using land surface parameters (LSPs).
    The function executes multiresolution segmentation of LSPs on the eCognition server and exports the results into shapefile.
    The size of the segments can be influenced by the scale parameter (the larger the value, the larger the objects).
    Requires Trimble eCognition server license and eCognition command line engine in docker installed.
    Ruleset file dem_seg.dcp, import template import_dem_seg.xml and input LSP rasters must be placed in the correct location for the docker container to access them (directory mounted in docker container).
    Paths to input image directory, results directory, import template file and ruleset must be defined relative to docker container filesystem.

    """

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

    process_dem_seg(
        img_dir=args.img_dir,
        results_dir=args.results_dir,
        sp=args.sp,
        import_tmpl_file=args.import_tmpl_file,
        ruleset=args.ruleset
    )

    print("Processing complete.")

