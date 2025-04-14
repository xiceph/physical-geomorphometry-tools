use gdal::raster::{RasterBand};
use gdal::{Dataset, DriverManager};
use gdal::cpl::{CslStringList};
use clap::{Command, Arg};
use std::fs;
use std::time::Instant;
use std::error::Error;
use std::io::{self, Write};
use std::sync::{mpsc, Arc};
use std::thread;
use num_cpus;
use rand::Rng;

mod math;
mod raster;
mod qem;
mod data;
mod text;

fn main() -> Result<(), Box<dyn Error>> {
  let start_time = Instant::now();

  // Define the command line arguments
  let app = Command::new("Raster QEM Generalize")
    .version(env!("CARGO_PKG_VERSION"))
    .author("Richard Feciskanin <richard.feciskanin@uniba.sk>")
    //.about("The program generates a generalized raster using the quadric error metric.")
    .arg(
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
      Arg::new("iterations")
        .short('n')
        .long("iterations")
        .default_value("10")
        .help("Specify the number of iterations"),
    )
    .arg(
      Arg::new("reduction")
        .short('r')
        .long("reduction")
        .default_value("1.0")
        .help("Specify the resolution reduction factor"),
    )
    .arg(
      Arg::new("sharpness")
        .short('s')
        .long("sharpness")
        .default_value("5")
        .help("Specify the sharpness value 1 - 9"),
    )
    .arg(
      Arg::new("jobs")
      .short('j')
      .long("jobs")
      .help("Specify the number of threads to use (if omitted, all available processors are used)"),
    );

  let line = "-".repeat(72);
  let dline = "=".repeat(72);

  println!("\n\
  {}\n\
  {}\n\
  Tool for generalizing an elevation raster using Quadric Error Metric.\n\
  Part of a {} project.\n\n\
  Author:\n{}\n\
  {}\n",
  text::highlight(format!("{} {}", text::bold("Raster QEM Generalize".to_string()), app.get_version().unwrap())),
  line,
  text::highlight("physical-geomorphometry".to_string()),
  app.get_author().unwrap(),
  dline);

  // Parse the command line arguments
  let matches = app.get_matches();

  let input_file = matches.get_one::<String>("input_file").unwrap();
  let output_file = matches.get_one::<String>("output_file").unwrap();
  let sharpness = matches.get_one::<String>("sharpness").unwrap().parse::<u32>()?;

  // Compute percentile range for thresholding based on sharpness
  let min_percentile: f32 = 1.0 - (0.1 * sharpness as f32);
  const MAX_PERCENTILE: f32 = 1.01;
  const MAX_ITERATIONS: u32 = 1999;

  // Parsing and validating 'iterations'
  let iterations = match matches.get_one::<String>("iterations").unwrap().parse::<u32>() {
    Ok(value) if value > 0 && value <= MAX_ITERATIONS => {
      if value > 500 {
        println!(
          "{}: 'iterations' value is unusually high ({}). This might cause long computation times.\n",
          text::warning("Warning".to_string()),
          value
        );
      }
      value
    }
    Ok(value) if value > MAX_ITERATIONS => {
      panic!("'iterations' must not exceed {}. Please provide lower value.", MAX_ITERATIONS);
    }
    Ok(_) => {
      panic!("'iterations' must be greater than 0.");
    }
    Err(_) => {
      panic!("'iterations' must be a valid positive integer.");
    }
  };

  // Parsing and validating 'reduction'
  let reduction = match matches.get_one::<String>("reduction").unwrap().parse::<f64>() {
    Ok(value) if value >= 1.0 => {
      if value == 1.0 && iterations > 200 {
        println!(
          "{}: 'iterations' value is high ({}) and not using resolution reduction. Consider using resolution reduction for faster calculation.\n",
          text::warning("Warning".to_string()),
          iterations
        );
      }
      value
    }
    Ok(_) => {
      panic!("'reduction' must be equal or greater than 1.0.");
    }
    Err(_) => {
      panic!("'reduction' must be a valid positive number.");
    }
  };

  // Determine the number of threads to use
  let num_procs = num_cpus::get() as usize;
  let jobs = if let Some(jobs_str) = matches.get_one::<String>("jobs") {
    match jobs_str.parse::<usize>() {
      Ok(max_jobs) if max_jobs > 0 => std::cmp::min(max_jobs, num_procs), // If valid and > 0, use the smaller value
      Ok(_) => {
        println!("{}: 'jobs' value must be greater than 0. Using the number of processors.\n", text::warning("Warning".to_string()));
        num_procs
      },
      Err(_) => {
        println!("{}: 'jobs' value is not a valid number. Using the number of processors.\n", text::warning("Warning".to_string()));
        num_procs
      }
    }
  } else {
    // If no value is provided, default to using the number of processors
    num_procs
  };

  // Summarize initial settings to the user
  println!("Generalize starts with settings:");
  println!("{}", line);
  println!("Number of iterations: {}", iterations);
  println!("Sharpness level: {}", sharpness);
  if reduction > 1.0 {
    println!("Resolution reduction factor: {}", reduction);
  } else {
    println!("Resolution remains unchanged");
  }
  println!("{}\n", dline);

  let mut part_time = Instant::now();

  // Open the input raster file
  let dataset = Dataset::open(input_file)
  .unwrap_or_else(|_| {
    let output = format!("{}: Failed to open input file: {}", text::error("Error".to_string()), input_file);
    eprintln!("{}\n", text::bold(output));
    std::process::exit(1);
  });

  // Get geotransform and projection of the dataset
  let geotransform = dataset.geo_transform().unwrap();
  let projection = dataset.projection();

  let rasterband: RasterBand = dataset.rasterband(1).unwrap();
  let no_data = rasterband.no_data_value();

  // Get the width and height (in pixels) of the original raster dataset.
  let width = dataset.raster_size().0;
  let height = dataset.raster_size().1;

  // Calculate the reduced width and height for the raster after applying the downscaling factor (reduction).
  let width_red: usize = (width as f64 / reduction as f64).ceil() as usize;
  let height_red: usize = (height as f64 / reduction as f64).ceil() as usize;

  // Create a new geotransform matrix for the reduced raster. This adjusts the pixel size based on the reduction.
  let mut new_geotransform = geotransform;
  new_geotransform[1] = geotransform[1] * (width as f64 / width_red as f64); // Pixel size on X
  new_geotransform[5] = geotransform[5] * (height as f64 / height_red as f64); // Pixel size on Y

  // Define the parameters for the new reduced raster
  let raster_params = raster::RasterParams {
    width: width_red,
    height: height_red,
    origin: [ geotransform[0], geotransform[3] ],
    resolution: [ new_geotransform[1], new_geotransform[5] ],
    nodata: no_data.unwrap_or(-9999.0),
  };

  // Read the elevation data into a grid structure
  let elevation = raster::Grid::<f32>::from_raster_band(raster_params, &rasterband);

  let elapsed_time1 = part_time.elapsed();
  if reduction > 1.0 {
    println!("{} Input raster ({}x{}) read and resampled to {}x{} in {:.2} s.", text::success("✓".to_string()), width, height, width_red, height_red, elapsed_time1.as_secs_f64());
  } else {
    println!("{} Input raster ({}x{}) read in {:.2} s.", text::success("✓".to_string()), width, height, elapsed_time1.as_secs_f64());
  }
  part_time = Instant::now();

  // Function to retrieve neighbors in a circular pattern around a given position in a 2D grid
  // It handles the rows before, after, and the current row, and applies a calculation function to each neighbor.
  fn get_circular_neighbors<T, U>(
      prev_row: Option<&Vec<T>>,                      // Previous row of the grid (optional)
      current_row: Option<&Vec<T>>,                   // Current row of the grid (optional)
      next_row: Option<&Vec<T>>,                      // Next row of the grid (optional)
      row_index: usize,                               // Index of the current row
      col: usize,                                     // Column index for the current element
      calculation_fn: impl Fn(&T, usize, usize) -> U, // Function to calculate each neighbor's value based on its position
  ) -> Vec<Option<U>> { // Return a vector of optional calculated neighbor values

    let mut neighbors: Vec<Option<U>> = Vec::with_capacity(8);
    for _ in 0..8 {
      neighbors.push(None);
    }

    // Check and retrieve neighbors in the previous row (if available)
    if let Some(row) = prev_row {
      if col > 0 {
        if let Some(value) = row.get(col - 1) {
          neighbors[0] = Some(calculation_fn(value, row_index - 1, col - 1));
        }
      }
      if let Some(value) = row.get(col) {
        neighbors[1] = Some(calculation_fn(value, row_index - 1, col));
      }
      if let Some(value) = row.get(col + 1) {
        neighbors[2] = Some(calculation_fn(value, row_index - 1, col + 1));
      }
    }
    // Add neighbors from the current row (if available)
    if let Some(row) = current_row {
      if col > 0 {
        if let Some(value) = row.get(col - 1) {
          neighbors[7] = Some(calculation_fn(value, row_index, col - 1));
        }
      }
      if let Some(value) = row.get(col + 1) {
        neighbors[3] = Some(calculation_fn(value, row_index, col + 1));
      }
    }
    // Add neighbors from the next row (if available)
    if let Some(row) = next_row {
      if col > 0 {
        if let Some(value) = row.get(col - 1) {
          neighbors[6] = Some(calculation_fn(value, row_index + 1, col - 1));
        }
      }
      if let Some(value) = row.get(col) {
        neighbors[5] = Some(calculation_fn(value, row_index + 1, col));
      }
      if let Some(value) = row.get(col + 1) {
        neighbors[4] = Some(calculation_fn(value, row_index + 1, col + 1));
      }
    }
    neighbors
  }

  // Retrieves the neighbors of a given column in a 7x7 grid, excluding the center.
  fn get_neighbors_7x7<U>(
      rows: [Option<&Vec<U>>; 7], // 7 rows of data, each row can be Some(Vec<U>) or None
      col_index: usize,           // current column index to fetch the neighbors for
  ) -> Vec<Option<U>> //Return a vector of `Option<U>` representing valid neighbors or `None` for missing values.
  where
    U: Clone,
  {
    // Initialize a vector to store neighbors
    let mut neighbors: Vec<Option<U>> = Vec::with_capacity(48);

    // Loop over all 7 rows and columns
    for i in 0..7 {
      if let Some(row) = rows[i] {
        for j in 0..7 {
          if !(i == 3 && j == 3) { // Skip the center element (the current row and column)
            let col_offset = col_index as isize + j as isize - 3;
            if col_offset >= 0 && col_offset < row.len() as isize {
              let neighbor = row.get(col_offset as usize).cloned();
              neighbors.push(neighbor);
            } else {
              neighbors.push(None);
            }
          }
        }
      } else {
        for j in 0..7 {
          if !(i == 3 && j == 3) {
            neighbors.push(None);
          }
        }
      }
    }

    neighbors
  }

  // Calculates the Q matrix for the given row by considering neighboring elevations.
  // Iterates through each column, checking for valid data points, and computes a Q matrix for each based on the neighbors.
  fn calc_row_q_matrix(
    current_row_index: usize,                 // The index of the current row being processed
    current_row: Option<&Vec<f32>>,           // Current row
    prev_row: Option<&Vec<f32>>,              // Previos row
    next_row: Option<&Vec<f32>>,              // Next row
    raster_params: raster::RasterParams // Raster parameters
  ) -> Vec<(usize, qem::SymmetricMatrix)> {

    let mut row_data = Vec::with_capacity(raster_params.width);

    // Helper closure to calculate 3D vector from raster value
    let calculation_fn = |value: &f32, row: usize, col: usize| -> math::Vector3 {
      let [x, y] = raster_params.get_xy_coords(row, col);
      math::Vector3 { x, y, z: *value as f64 }
    };

    for col in 0..raster_params.width {

      // If there's a value at the current position
      if let Some(elevation) = current_row.unwrap().get(col) {

        if !raster_params.is_nodata(*elevation as f64) {

          let mut q = qem::SymmetricMatrix::new();

          // Get the 3D coordinates of the current point
          let [x, y] = raster_params.get_xy_coords(current_row_index, col);
          let v = math::Vector3 { x, y, z: *elevation as f64 };

          let neighbors = get_circular_neighbors(prev_row, current_row, next_row, current_row_index, col, calculation_fn);

          // Iterate over neighbors and calculate Q matrix for the plane formed by current and adjacent points
          for i in 0..8 {
            // Check for valid neighboring values
            if let (Some(v1), Some(v2)) = (
              neighbors[i],
              if i < 7 {
                neighbors[i + 1]
              } else {
                neighbors[0]
              }
            ) {
              if !raster_params.is_nodata(v1.z) && !raster_params.is_nodata(v2.z) {
                // Calculate the plane parameters and update the Q matrix
                let [a, b, c, d] = math::calc_plane_params(v, v1, v2);
                let mut qn = qem::SymmetricMatrix::new();
                qn.make_plane(a, b, c, d);
                q.add_self(&qn);
              }
            }
          }
          // Store the computed Q matrix for the current column
          row_data.push((col, q));
        }
      }
    }
    row_data
  }

  // Calculates the generalized Q matrix for the given row using neighboring Q matrices and elevation values.
  // For each position in the row, it computes a new Q matrix based on its neighbors, the current elevation, and a specified error threshold.
  fn calc_row_generalized(
    current_row_index: usize,                             // The index of the current row being processed
    rows: [Option<&Vec<qem::SymmetricMatrix>>; 7],  // 7 rows as references
    current_elevation: Option<&Vec<f32>>,                 // Current row's elevation values, optional
    raster_params: raster::RasterParams,            // Raster parameters
    threshold: f64,                                       // Error threshold for updating the Q matrix
  ) -> Vec<(usize, qem::SymmetricMatrix, f32, f64)> {

    let mut row_data = Vec::with_capacity(raster_params.width);

    let current_row = rows[3];

    for col in 0..raster_params.width {

      // If there's a value at the current position
      if let Some(elevation) = current_elevation.unwrap().get(col) {

        if !raster_params.is_nodata(*elevation as f64) {

          let [x, y] = raster_params.get_xy_coords(current_row_index, col);

          // New Q matrix for grid node
          let mut q_next = qem::SymmetricMatrix::new();
          if let Some(matrix) = current_row.and_then(|row| row.get(col)) {
            q_next = matrix.clone();
          }

          // Get the neighboring Q matrices (from surrounding rows)
          let neighbors = get_neighbors_7x7(rows, col);

          let mut used_count = 0;
          let mut errors = Vec::new();

          // Iterate through the neighbors and calculate errors
          for neighbor_q in &neighbors {
            if let Some(q) = neighbor_q {
              let error = qem::calc_vertex_error(q, x, y, *elevation as f64);
              errors.push(error);

              // If the error is within the threshold, add the neighbor's Q matrix to the current one
              if error > 0.0 && error < threshold {
                q_next.add_self(&neighbor_q.unwrap());
                used_count += 1;
              }
            }
          }

          // If any valid neighbors were used, average the Q matrix
          if used_count > 0 {
            q_next.divide((used_count + 1) as f64);
          }

          // Calculate the optimal Z value (elevation) from the new Q matrix
          let new_z = qem::calc_optimal_z(&q_next, x, y);

          // Calculate the 75th percentile error from the errors array
          let q3_error = data::percentile(&mut errors, 75.0);

          // Store the results for this column
          row_data.push((col, q_next, new_z as f32, q3_error));
        }
      }
    }
    row_data // Return the computed data for the row
  }

  match elevation {
    Ok(grid) => {
      // Create a grid of Q matrices with default values
      let default_value = qem::SymmetricMatrix::new();
      let mut q_matrices_next = raster::Grid::<qem::SymmetricMatrix>::from_params(raster_params, default_value);

      let mut new_elevation = grid.clone();
      let grid_arc = Arc::new(grid);
      let height = grid_arc.params().height;

      let (tx, rx) = mpsc::channel();
      let mut handles = vec![];

      // Process rows in parallel using threads
      for tid in 0..jobs {
        let tx = tx.clone();
        let grid_clone = Arc::clone(&grid_arc);
        let handle = thread::spawn(move || {
          for row in (0..height).filter(|r| r % jobs == tid) {
            let current_row = Some(grid_clone.get_row(row).unwrap());
            let prev_row = if row > 0 { Some(grid_clone.get_row(row - 1).unwrap()) } else { None };
            let next_row = if row < height - 1 { Some(grid_clone.get_row(row + 1).unwrap()) } else { None };

            let row_data = calc_row_q_matrix(row, current_row, prev_row, next_row, raster_params);
            tx.send((row, row_data)).unwrap();
          }
        });
        handles.push(handle);
      }
      // Wait for all threads to finish
      for handle in handles {
        handle.join().unwrap();
      }

      // Update grid with calculated data
      for _ in 0..height {
        let (row, row_data) = rx.recv().unwrap();
        for (col, q) in row_data {
          q_matrices_next.set_item(row, col, q);
        }
      }

      let elapsed_time = part_time.elapsed();
      println!("{} Q matrices initialized in {:.2} s.", text::success("✓".to_string()), elapsed_time.as_secs_f64());
      part_time = Instant::now();

      let mut current_error = 0.0;

      // Main loop for performing generalized raster processing using parallelized threading and error thresholding
      // The loop iterates over the specified number of iterations, updating Q matrices and elevations for each row in the raster.
      for iteration in 0..=iterations {

        // Set the error threshold
        let threshold = if iteration == 0 {
          999_999.9
        } else {
          current_error
        };

        let q_matrices_arc = Arc::new(q_matrices_next.clone());

        let (tx, rx) = mpsc::channel();
        let mut handles = vec![];

        // Process rows in parallel using threads
        for tid in 0..jobs {
          let tx = tx.clone();
          let q_matrices_ref = Arc::clone(&q_matrices_arc);
          let elevation_ref = Arc::clone(&grid_arc);
          let handle = thread::spawn(move || {
            for row in (0..height).filter(|r| r % jobs == tid) {
              let current_row = Some(q_matrices_ref.get_row(row).unwrap());
              let prev_row = if row > 0 { Some(q_matrices_ref.get_row(row - 1).unwrap()) } else { None };
              let prev2_row = if row > 1 { Some(q_matrices_ref.get_row(row - 2).unwrap()) } else { None };
              let prev3_row = if row > 2 { Some(q_matrices_ref.get_row(row - 3).unwrap()) } else { None };
              let next_row = if row < height - 1 { Some(q_matrices_ref.get_row(row + 1).unwrap()) } else { None };
              let next2_row = if row < height - 2 { Some(q_matrices_ref.get_row(row + 2).unwrap()) } else { None };
              let next3_row = if row < height - 3 { Some(q_matrices_ref.get_row(row + 3).unwrap()) } else { None };
              let current_elevation = Some(elevation_ref.get_row(row).unwrap());

              let rows = [prev3_row, prev2_row, prev_row, current_row, next_row, next2_row, next3_row];

              let row_data = calc_row_generalized(row, rows, current_elevation, raster_params, threshold);
              tx.send((row, row_data)).unwrap();
            }
          });
          handles.push(handle);
        }

        // Wait for all threads to finish
        for handle in handles {
          handle.join().unwrap();
        }

        let mut errors = Vec::new();

        // Update grid with calculated data
        for _ in 0..height {
          let (row, row_data) = rx.recv().unwrap();
          for (col, q, elevation, error) in row_data {
            // Update the Q matrix and elevation values for each column in the row
            q_matrices_next.set_item(row, col, q);
            new_elevation.set_item(row, col, elevation);

            // Collect the error values to compute error metrics
            errors.push(error);
          }
        }

        // Calculate the current percentile based on iteration and threshold
        let mut percentile = min_percentile + ((iteration + 1) as f32 / iterations as f32) * (MAX_PERCENTILE - min_percentile);
        if percentile > 1.0 {
          percentile = 1.0;
        }

        // Calculate the current error as the percentile of collected errors
        current_error = data::percentile(&mut errors, (percentile * 100.0) as f64);

        print!("\r{}- Generalization progress: {} %", "\x1B[K", (iteration * 100) / iterations);
        io::stdout().flush().unwrap();
      }

      let elapsed_time = part_time.elapsed();
      println!("\r\x1B[K{} Generalized in {:.2} s.", text::success("✓".to_string()), elapsed_time.as_secs_f64());

      // Helper function to generate a random string for unique raster naming
      fn generate_random_string() -> String {
        let letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        let mut rng = rand::thread_rng();
        let letter1 = letters.chars().nth(rng.gen_range(0..letters.len())).unwrap();
        let letter2 = letters.chars().nth(rng.gen_range(0..letters.len())).unwrap();
        format!("{}{}", letter1, letter2)
      }

      // Save the generalized raster to an output file (GeoTIFF format)
      let result = new_elevation.to_gdal_buffer();
      match result {
        Ok(mut buffer) => {
          // Create the output file path, appending a unique string if file already exists
          let mut output_file_safe = output_file.to_string();
          if fs::metadata(&output_file).is_ok() {
            let new_part = format!("_{}", generate_random_string());
            if let Some(index) = output_file.rfind('.') {
              output_file_safe.insert_str(index, &new_part);
            } else {
              output_file_safe.push_str(&new_part);
            }
          }

          // Create the GeoTIFF file and set its properties
          let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
          let options = CslStringList::from_iter(["TILED=YES", "BLOCKXSIZE=16", "BLOCKYSIZE=16"]);
          let mut ds = driver
            .create_with_band_type_with_options::<f32, _>(
              &output_file_safe,
              width_red,
              height_red,
              1,
              &options,
            ).expect("Failed to create GeoTIFF file.");

          ds.set_projection(&projection)?;
          ds.set_geo_transform(&new_geotransform).expect("Failed to set geotransform");

           // Write the buffer data to the raster band
          let mut band1 = ds.rasterband(1).unwrap();
          band1.set_no_data_value(Some(raster_params.nodata))?;
          band1.write((0, 0), (width_red, height_red), &mut buffer).expect("Failed to write data");
          band1.compute_raster_min_max(true)?;
        }
        Err(err) => {
          println!("Failed to create output raster: {}", err);
        }
      }
    }
    Err(err) => {
      println!("Error: {}", err);
    }
  }

  let elapsed_time = start_time.elapsed();
  println!("{}", line);
  println!("{}", text::success("Calculations completed successfully.".to_string()));
  println!("Total elapsed time: {:.2} s.", elapsed_time.as_secs_f64());
  println!("");

  Ok(())
}
