use clap::{Arg, Command};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

mod crop;
mod normalize;
mod rescale;
mod text;

/// Create a path for intermediate files. If `keep_temps` is true the name will be
/// `<base_stem>_<step>.tif`. If false the name will be a unique tmp name:
/// `<base_stem>.tmp.<step>.<pid>.<timestamp>.tif`.
fn make_intermediate_path(output_dir: &Path, base_stem: &str, step: &str, keep_temps: bool) -> PathBuf {
    if keep_temps {
        output_dir.join(format!("{}_{}.tif", base_stem, step))
    } else {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Time went backwards")
            .as_millis();
        let pid = std::process::id();
        output_dir.join(format!("{}.tmp.{}.{}.{}.tif", base_stem, step, pid, ts))
    }
}

/// CSV path for normalization when keep_temps is true (otherwise None)
fn make_k_csv_path(output_dir: &Path, base_stem: &str, keep_temps: bool) -> Option<PathBuf> {
    if keep_temps {
        Some(output_dir.join(format!("{}_k.csv", base_stem)))
    } else {
        None
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    let mut app = Command::new("Raster Tools Pipeline")
        .version(env!("CARGO_PKG_VERSION"))
        .author("Veronika Hajdúchová <veronika.hajduchova@uniba.sk>");
    app = app.about("A command-line pipeline tool to crop, normalize, and rescale rasters.");

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!("\n\
    {}\n\
    {}\n\
    {}\n\
    Part of a {} project.\n\n\
    Author:\n{}\n\
    {}\n",
    format!("{} {}", text::highlight("Raster Tools Pipeline"), text::highlight(app.get_version().unwrap())),
    line,
    app.get_about().unwrap(),
    text::highlight("physical-geomorphometry"),
    app.get_author().unwrap(),
    dline);

    app = app.arg(
            Arg::new("input_file")
                .short('i')
                .long("input-file")
                .required(true)
                .help("Specify the input file path"),
        )
        .arg(
            Arg::new("output_file")
                .short('o')
                .long("output-file")
                .required(true)
                .help("Specify the output file path"),
        )
        .arg(
            Arg::new("crop")
                .short('c')
                .long("crop")
                .num_args(1)
                .value_parser(clap::value_parser!(usize))
                .help("Crop N pixels from all sides"),
        )
        .arg(
            Arg::new("normalize")
                .short('n')
                .long("normalize")
                .action(clap::ArgAction::SetTrue)
                .help("Normalize the raster"),
        )
        .arg(
            Arg::new("rescale")
                .short('r')
                .long("rescale")
                .action(clap::ArgAction::SetTrue)
                .help("Rescale raster values"),
        )
        .arg(
            Arg::new("keep_temps")
                .long("keep-temps")
                .action(clap::ArgAction::SetTrue)
                .help("Keep intermediate rasters (_cropped, _normalized, CSV)"),
        );

    let matches = app.get_matches();

    let input_path = PathBuf::from(matches.get_one::<String>("input_file").unwrap());
    let output_path = PathBuf::from(matches.get_one::<String>("output_file").unwrap());
    let crop_margin = matches.get_one::<usize>("crop").copied();
    let do_normalize = matches.get_flag("normalize");
    let do_rescale = matches.get_flag("rescale");
    let keep_temps = matches.get_flag("keep_temps");

    // Where to place intermediates: use the directory of the output file (or current dir).
    let output_dir = output_path
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));

    // base stem derived from input filename (so intermediate names relate to input)
    let base_stem = input_path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "raster".to_string());

    println!("Starting pipeline...\n");

    // --- Verify input file exists ---
    if !fs::metadata(&input_path).is_ok() {
        eprintln!("{} Input file {:?} does not exist.", text::error_icon(), input_path);
        std::process::exit(1);
    }
    let mut part_time = Instant::now();

    // --- Raster metadata info ---
    println!("Reading input raster...");
    let dataset = gdal::Dataset::open(&input_path)?;
    let band = dataset.rasterband(1)?;
    let (width, height) = band.size();
    let gt = dataset.geo_transform()?;

    let pixel_size_x = gt[1];
    let pixel_size_y = gt[5]; // usually negative

    let elapsed_time = part_time.elapsed();
    println!("{} Input raster ({} x {}) read in {:.2} seconds.", text::check_icon(), width, height, elapsed_time.as_secs_f64());
    println!("  Pixel size: {} x {}", pixel_size_x, pixel_size_y.abs());
    /*if proj.is_empty() {
        print_success("Projection: None");
    } else {
        print_success(&format!("Projection: {}", proj));
    }*/

    // Track the file flowing through the pipeline
    let mut current_input = input_path.clone();
    let mut did_any_step = false;
    let mut temp_files: Vec<PathBuf> = Vec::new();

    // 1. Crop
    if let Some(margin) = crop_margin {
        let tmp_output = make_intermediate_path(&output_dir, &base_stem, "cropped", keep_temps);
        part_time = Instant::now();

        if keep_temps {
            println!("\nCropping raster ({} px margin) {} {:?}", margin, text::arrow_icon(), tmp_output);
        } else {
            println!("\nCropping raster ({} px margin)...", margin);
        }
        crop::crop_raster(current_input.clone(), tmp_output.clone(), margin)?;
        let elapsed_time = part_time.elapsed();
        println!("{} Cropped raster created in {:.2} seconds.", text::check_icon(), elapsed_time.as_secs_f64());

        current_input = tmp_output;
        if !keep_temps {
            temp_files.push(current_input.clone());
        }
        did_any_step = true;
    }

    // 2. Normalize
    if do_normalize {
        let tmp_output = make_intermediate_path(&output_dir, &base_stem, "normalized", keep_temps);
        let tmp_csv = make_k_csv_path(&output_dir, &base_stem, keep_temps);
        part_time = Instant::now();

        if keep_temps {
            println!("\nNormalizing raster {} {:?}", text::arrow_icon(), tmp_output);
        } else {
            println!("\nNormalizing raster...");
        }
        normalize::normalize_raster(
            current_input.clone(),
            tmp_output.clone(),
            tmp_csv,
            0.1,  // initial_k
            1e-2, // tolerance
        )?;
        let elapsed_time = part_time.elapsed();
        println!("{} Normalized raster created in {:.2} seconds.", text::check_icon(), elapsed_time.as_secs_f64());
        current_input = tmp_output;
        if !keep_temps {
            temp_files.push(current_input.clone());
        }
        did_any_step = true;
    }

    // 3. Rescale
    if do_rescale {
        let tmp_output = make_intermediate_path(&output_dir, &base_stem, "rescaled", keep_temps);
        part_time = Instant::now();

        if keep_temps {
            println!("\nRescaling raster {} {:?}", text::arrow_icon(), tmp_output);
        } else {
            println!("\nRescaling raster...");
        }
        rescale::rescale_raster(current_input.clone(), tmp_output.clone())?;
        let elapsed_time = part_time.elapsed();
        println!("{} Rescaled raster created in {:.2} seconds.", text::check_icon(), elapsed_time.as_secs_f64());
        current_input = tmp_output;
        if !keep_temps {
            temp_files.push(current_input.clone());
        }
        did_any_step = true;
    }

    // === Finalize: ensure final output_path contains the final raster ===
    println!("\nFinalizing output file...");
    if output_path == current_input {
        // Already the desired path (rare with this scheme), nothing to do
        println!("{} Output already at {:?}", text::check_icon(), output_path);
    } else {
        if !did_any_step {
            // No processing requested. If input != output, copy input -> output (do not remove input).
            if input_path != output_path {
                fs::copy(&input_path, &output_path)?;
                println!("{} No processing requested — copied input {} {:?}", text::check_icon(), text::arrow_icon(), output_path);
            } else  {
                println!("{} No processing requested and input == output — nothing to do", text::check_icon());
            }
        } else {
            // Did perform steps. If keep_temps, copy the final intermediate to output_path (keep intermediate).
            // If not keep_temps, attempt an atomic rename (move) of the temp into output_path; fallback to copy+remove.
            if keep_temps {
                fs::copy(&current_input, &output_path)?;
                println!("{} {}", text::check_icon(), &format!(
                    "Kept intermediate files; copied final result → {:?}",
                    output_path
                ));
            } else {
                // The last temp file is the final result, so don't remove it from the list.
                // try rename first
                match fs::rename(&current_input, &output_path) {
                    Ok(_) => {
                        println!("{} Moved temporary final {} {:?}", text::check_icon(), text::arrow_icon(), output_path);
                    }
                    Err(e) => {
                        // fallback: copy then remove the temp
                        println!("Rename failed ({}) — falling back to copy + remove", e);
                        fs::copy(&current_input, &output_path)?;
                        // ignore remove error but try
                        if let Err(rem_e) = fs::remove_file(&current_input) {
                            eprintln!("{}: failed to remove temp file {:?}: {}", text::warning("Warning"), current_input, rem_e);
                        }
                        println!("{} Copied temporary final {} {:?}", text::check_icon(), text::arrow_icon(), output_path);
                    }
                }

                // Clean up all other intermediate temp files
                for temp_file in temp_files.iter().filter(|p| **p != current_input) {
                    if let Err(rem_e) = fs::remove_file(temp_file) {
                        eprintln!("{}: failed to remove temp file {:?}: {}", text::warning("Warning"), temp_file, rem_e);
                    }
                }
            }
        }
    }

    // --- Final raster info ---
    println!("\nReading final raster info...");
    let final_ds = gdal::Dataset::open(&output_path)?;
    let final_band = final_ds.rasterband(1)?;
    let (final_width, final_height) = final_band.size();
    let final_gt = final_ds.geo_transform()?;

    let final_pixel_size_x = final_gt[1];
    let final_pixel_size_y = final_gt[5];

    println!("{} Final size: {} x {} px", text::check_icon(), final_width, final_height);
    println!("{} {}", text::check_icon(), &format!(
        "Final pixel size: {} x {}",
        final_pixel_size_x,
        final_pixel_size_y.abs()
    ));
    /*if final_proj.is_empty() {
        print_success("Projection: None");
    } else {
        print_success(&format!("Projection: {}", final_proj));
    }*/

    let elapsed_time = start_time.elapsed();
    println!("\n{}", dline);
    println!("{}", text::success("Pipeline completed successfully."));
    println!("Final output: {:?}", output_path);
    println!("Total elapsed time: {:.2} seconds.", elapsed_time.as_secs_f64());

    Ok(())
}
