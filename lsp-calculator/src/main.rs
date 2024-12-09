use gdal::raster::{RasterBand, RasterCreationOption};
use gdal::{Dataset, DriverManager};
use clap::{Command, Arg, ArgAction};
use std::time::Instant;
use std::error::Error;
use std::io::{self, Write};
use std::sync::{mpsc, Arc, Mutex};
use std::thread;
use std::collections::HashSet;
use num_cpus;
use nalgebra::DMatrix;

mod poly;
mod lsp;
mod raster;
mod text;

// The function generates text for individual parameter selection.
fn generate_individual_params_help(parameters: &Vec<(&str, &str)>) -> String {
  parameters
  .iter()
  .map(|(name, title)| {
    // Maximum length of parameter name plus spaces
    let param_display_width: usize = 28;
    // Calculate spaces needed to align
    let spaces_needed = param_display_width.saturating_sub(name.len());
    let spaces = " ".repeat(spaces_needed);
    format!("  {}{}{}\n", text::bold(format!("--{}", name.to_string())), spaces, title)
  })
  .collect::<String>()
}

fn main() -> Result<(), Box<dyn Error>> {
  let start_time = Instant::now();
  
  // Define the command line arguments
  let parameters = vec![
    // First order parameters
    ("slope", "Slope"),
    ("aspect", "Aspect"),
    ("sin_slope", "Sine of Slope"),
    ("sin_aspect", "Sine of Aspect"),
    ("cos_aspect", "Cosine of Aspect"),
    // Basic trio of curvatures
    ("kns", "Normal slope line (profile) curvature"),
    ("knc", "Normal contour (tangential) curvature"),
    ("tc", "Contour torsion"),
    // Subforms of Basic trio
    ("zss", "Second slope line derivative"),
    ("ts", "Slope line torsion"),
    ("zcc", "Second contour derivative"),
    ("kpc", "Projected contour curvature"),
    ("kps", "Projected slope line curvature"),
    ("sin_sc", "Contour change of sin slope"),
    // Other gravity specific curvatures
    ("kd", "Difference curvature"),
    ("ka", "Total accumulation curvature"),
    ("kr", "Total ring curvature"),
    ("khe", "Horizontal excess curvature"),
    ("kve", "Vertical excess curvature"),
    // Principal curvatures
    ("k_max", "Maximal curvature"),
    ("k_min", "Minimal curvature"),
    // Other gravity-invariant curvatures
    ("k", "Gaussian curvature"),
    ("el", "Elevation laplacian"),
    ("ku", "Unsphericity curvature"),
    ("k_mean", "Mean curvature"),
    ("kc", "Casorati curvature"),
    // Changes of curvatures
    ("knss", "Slope line change of normal slope line curvature"),
    ("kncc", "Contour change of normal contour curvature"),
    ("kncs", "Slope line change of normal contour curvature"),
  ];
  //let total_parameters = parameters.len();

  let group_first = vec!["slope", "aspect", "sin_slope", "sin_aspect", "cos_aspect"];
  let group_second = vec!["kns", "knc", "tc", "zss", "ts", "zcc", "kpc", "kps", "sin_sc", "kd", "ka", "kr", "khe", "kve", "k_max", "k_min", "k", "el", "ku", "k_mean", "kc"];
  let group_third = vec!["knss", "kncc", "kncs"];
  let group_segmentation = vec!["sin_slope", "sin_aspect", "cos_aspect", "kns", "knc", "tc", "knss", "kncc", "kncs"];
  
  let mut app = Command::new("Land Surface Parameters Calculator")
  .version(env!("CARGO_PKG_VERSION"))
  .author("Richard Feciskanin <richard.feciskanin@uniba.sk>\nVeronika Hajdúchová <hajduchova18@uniba.sk>")
  //.about("The program calculates Land Surface Parameters for an elevation raster.")
  .after_help(&format!("{}\n  \
    {}{}Compute all Land Surface Parameters offered by the tool (listed below)\n  \
    Or select a subset of parameters to compute (one or more options):\n  \
    {}{}First-order Land Surface Parameters\n  \
    {}{}Second-order Land Surface Parameters\n  \
    {}{}Second-order Land Surface Parameters\n  \
    {}{}Land Surface Parameters for Land Surface Segmentation\n\n\
    {}\n{}",
    text::bold(format!("{} of calculated parameters:", text::underline("Batch selection".to_string()))),
    text::bold("--all".to_string()), " ".repeat(25),
    text::bold("--first".to_string()), " ".repeat(23),
    text::bold("--second".to_string()), " ".repeat(22),
    text::bold("--third".to_string()), " ".repeat(23),
    text::bold("--for_segmentation".to_string()), " ".repeat(12),
    text::bold(format!("{} of calculated parameters (choose one or more specific parameters):", text::underline("Individual selection".to_string()))),
    generate_individual_params_help(&parameters)
  ))
  .arg(
    Arg::new("input_file")
    .short('i')
    .long("input-file")
    .value_name("file")
    .required(true)
    .help("Specify the input file path"),
  )
  .arg(
    Arg::new("output_prefix")
    .short('o')
    .long("output-prefix")
    .value_name("prefix")
    .required(true)
    .help("Specify the output file(s) prefix"),
  )
  .arg(
    Arg::new("degree")
    .short('d')
    .long("degree")
    .default_value("3")
    .help("Specify the polynomial degree (3 or 4)"),
  )
  .arg(
    Arg::new("jobs")
    .short('j')
    .long("jobs")
    .help("Specify the number of threads to use (if omitted, all available processors are used)"),
  )
  .arg(
    Arg::new("all")
    .short('a')
    .long("all")
    .help("Compute every Land Surface Parameter offered by the tool")
    .action(ArgAction::SetTrue)
    .hide(true),
  )
  .arg(
    Arg::new("first")
    .long("first")
    .help("Calculate the First-order Land Surface Parameters")
    .action(ArgAction::SetTrue)
    .hide(true),
  )
  .arg(
    Arg::new("second")
    .long("second")
    .help("Calculate the Second-order Land Surface Parameters")
    .action(ArgAction::SetTrue)
    .hide(true),
  )
  .arg(
    Arg::new("third")
    .long("third")
    .help("Calculate the Third-order Land Surface Parameters")
    .action(ArgAction::SetTrue)
    .hide(true),
  )
  .arg(
    Arg::new("for_segmentation")
    .long("for_segmentation")
    .help("Calculate Land Surface Parameters for Land Surface Segmentation")
    .action(ArgAction::SetTrue)
    .hide(true),
  );
  
  let mut list = vec![];
  for (name, title) in &parameters {
    app = app.arg(
      Arg::new(name)
      .long(name)
      .help(&format!("Calculate parameter {}", &title))
      .action(ArgAction::SetTrue)
      .hide(true)
    );
    list.push(*name);
  }

  let line = "-".repeat(72);
  let dline = "=".repeat(72);
  
  println!("\n\
  {}\n\
  {}\n\
  Tool for calculating Land Surface Parameters from an elevation raster.\n\
  Part of a {} project.\n\n\
  Authors:\n{}\n\
  {}\n",
  text::highlight(format!("{} {}", text::bold("Land Surface Parameters Calculator".to_string()), app.get_version().unwrap())),
  line,
  text::highlight("physical-geomorphometry".to_string()),
  app.get_author().unwrap(),
  dline);

  // Parse command-line arguments based on the options defined above
  let matches = app.get_matches();
  
  let input_file = matches.get_one::<String>("input_file").unwrap();
  let output_prefix = matches.get_one::<String>("output_prefix").unwrap();
  
  let degree = matches.get_one::<String>("degree").unwrap().parse::<usize>()?;
  if degree != 3 && degree != 4 {
    panic!("Unsupported polynomial degree. Please use 3 or 4.");
  }

  let num_procs = num_cpus::get() as usize;
  let jobs = if let Some(jobs_str) = matches.get_one::<String>("jobs") {
    // Attempt to convert the string to usize
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

  let all = matches.get_flag("all");

  let mut selected_params_set: HashSet<String> = HashSet::new();

  if matches.get_flag("first"){
    selected_params_set.extend(group_first.into_iter().map(|s| s.to_string()));
  }
  if matches.get_flag("second"){
    selected_params_set.extend(group_second.into_iter().map(|s| s.to_string()));
  }
  if matches.get_flag("third"){
    selected_params_set.extend(group_third.into_iter().map(|s| s.to_string()));
  }
  if matches.get_flag("for_segmentation"){
    selected_params_set.extend(group_segmentation.into_iter().map(|s| s.to_string()));
  }

  for (name, _) in &parameters {
    if matches.get_flag(name) || all {
      selected_params_set.insert(name.to_string());
    }
  }
  let total_selected = selected_params_set.len();

  println!("Starting calculation using a degree {} polynomial approximation.\n", degree);
  println!("The following parameters will be calculated [{}]:", total_selected);
  println!("{}", line);
  
  for (name, title) in &parameters {
    if selected_params_set.contains(&name.to_string()) {
      println!("{}\n  {}", title, text::light(format!("└─▶ {}_{}.tif", output_prefix, name)));
    }
  }
  println!("{}\n", dline);
  
  let selected_params: Vec<String> = selected_params_set.into_iter().collect();

  let mut part_time = Instant::now();
  
  // Open the input raster file
  let dataset = Dataset::open(input_file)
  .unwrap_or_else(|_| {
    let output = format!("{}: Failed to open input file: {}", text::error("Error".to_string()), input_file);
    eprintln!("{}\n", text::bold(output));
    std::process::exit(1);
  });
  
  let geotransform = dataset.geo_transform().unwrap();
  let projection = dataset.projection();
  
  let rasterband: RasterBand = dataset.rasterband(1).unwrap();
  let no_data = rasterband.no_data_value();
  
  let width = dataset.raster_size().0;
  let height = dataset.raster_size().1;
  
  let raster_params = raster::RasterParams {
    width: width,
    height: height,
    origin: [ geotransform[0], geotransform[3] ],
    resolution: [ geotransform[1], geotransform[5] ],
    nodata: no_data.unwrap_or(-9999.0),
  };
  
  let elevation = raster::Grid::<f32>::from_raster_band(raster_params, &rasterband);
  
  let elapsed_time1 = part_time.elapsed();
  println!("{} Input raster ({} x {}) read in {:.2} seconds.", text::success("✓".to_string()), width, height, elapsed_time1.as_secs_f64());
  part_time = Instant::now();
  
  // Retrieves the 5x5 neighborhood around the specified row and column index,
  // returning the indices of used nodes and Z (elevation) coordinates needed for the polynomial approximation.
  fn get_neighbors_5x5_indices_and_z(
    rows: [Option<&Vec<f32>>; 5],         // 5 rows as an array of references
    col_index: usize,                     // Index of the current column
    raster_params: raster::RasterParams,  // Raster metadata
  ) -> (Vec<usize>, Vec<f64>) {
    let mut indices: Vec<usize> = Vec::with_capacity(25);
    let mut z: Vec<f64> = Vec::with_capacity(25);
    
    // Loop through the 5x5 grid
    for i in 0..5 {
      if let Some(row) = rows[i] {
        for j in 0..5 {
          let col_offset = col_index as isize + j as isize - 2;
          // Check if the column index is within bounds
          if col_offset >= 0 && col_offset < row.len() as isize {
            let z_value = row.get(col_offset as usize).unwrap();
            // Skip NoData values
            if !raster_params.is_nodata((*z_value).into()) {
              // Calculate the grid index based on i and j
              let idx = i * 5 + j;
              // Add the index and corresponding z value
              indices.push(idx);
              z.push((*z_value).into());
            }
          }
        }
      }
    }
    
    // Return the valid indices and corresponding z values
    (indices, z)
  }
  
  // Calculate the Land Surface Parameters (LSPs) for each row using its neighboring rows
  fn calc_row_lsp(
    rows: [Option<&Vec<f32>>; 5],         // Array of 5 optional references to row vectors
    raster_params: raster::RasterParams,  // Raster metadata
    selected_params: &Vec<String>,        // A vector of strings representing the selected LSPs
    precomputed_bf: &DMatrix<f64>,        // Precomputed basis functions matrix
    precomputed_m: &DMatrix<f64>          // Precomputed matrix M
  ) -> Vec<(usize, Option<Vec<f32>>)> {
    
    let mut row_data = Vec::with_capacity(raster_params.width);
    
    let current_row = rows[2];
    
    for col in 0..raster_params.width {
      // Check if there is a valid elevation value at the current position
      if let Some(elevation) = current_row.unwrap().get(col) {
        
        if !raster_params.is_nodata(*elevation as f64) {
          
          // Get the z values and corresponding indices of the 5x5 neighboring grid around the current pixel
          let (indices, z) = get_neighbors_5x5_indices_and_z(rows, col, raster_params);
          
          if indices.len() >= precomputed_bf.ncols() + 6 {
            // Extract only the relevant rows from the precomputed basis function matrix (based on valid indices)
            let mut filtered_bf = DMatrix::zeros(indices.len(), precomputed_bf.ncols());
            for (i, &idx) in indices.iter().enumerate() {
              filtered_bf.set_row(i, &precomputed_bf.row(idx));
            }
            
            // Perform polynomial least squares fitting to calculate the derivatives
            let derivatives = if indices.len() == 25 {
              // Use the precomputed matrix M if we have all 25 points
              poly::poly_lsq_der_with_precomputed(&z, &filtered_bf, &precomputed_m)
            } else {
              // Dynamically compute the matrix M for the filtered basis functions
              let filtered_m = poly::compute_matrix_m(&filtered_bf);

              // Perform polynomial least squares fitting with the newly computed matrix M
              poly::poly_lsq_der_with_precomputed(&z, &filtered_bf, &filtered_m)
            };
            
            match derivatives {
              Some(derivatives) => {
                let mut param_values = Vec::with_capacity(selected_params.len());
                // Calculate values of selected LSPs 
                for name in selected_params {
                  let value = lsp::calculate(name, &derivatives);
                  param_values.push(value);
                }
                
                row_data.push((col, Some(param_values)));
              }
              None => {
                row_data.push((col, None));
              }
            }
          } 
        }
      }
    }
    row_data
  }
  
  match elevation {
    Ok(grid) => {
      
      print!("Calculating Land Surface Parameters...");
      io::stdout().flush().unwrap();
      
      let d_x = raster_params.resolution[0];
      let d_y = raster_params.resolution[1];
      let basis_functions = Arc::new(poly::precompute_basis_functions(degree, d_x, d_y));
      let m = Arc::new(poly::compute_matrix_m(&basis_functions));
      
      let grid_arc = Arc::new(grid);
      let height = grid_arc.params().height;
      
      let selected_params_arc = Arc::new(selected_params.clone());

      let count = Arc::new(Mutex::new(0));
      
      let (tx, rx) = mpsc::channel();
      let mut handles = vec![];
      
      // Spawn multiple threads for parallel row processing.
      // Each thread processes rows that match its thread ID using a modulo operation.
      // The neighborhood of 5 rows around the current row is gathered, and the selected parameters are calculated.
      // Results are sent back to the main thread through an MPSC channel.
      for tid in 0..jobs {
        
        let tx = tx.clone();
        let grid_clone = Arc::clone(&grid_arc);
        let selected_params_clone = Arc::clone(&selected_params_arc);
        let basis_functions_clone = Arc::clone(&basis_functions);
        let m_clone = Arc::clone(&m);
        let count_clone = Arc::clone(&count);
        
        let handle = thread::spawn(move || {
          for row in (0..height).filter(|r| r % jobs == tid) {
            let current_row = Some(grid_clone.get_row(row).unwrap());
            let prev_row = if row > 0 { Some(grid_clone.get_row(row - 1).unwrap()) } else { None };
            let prev2_row = if row > 1 { Some(grid_clone.get_row(row - 2).unwrap()) } else { None };
            let next_row = if row < height - 1 { Some(grid_clone.get_row(row + 1).unwrap()) } else { None };
            let next2_row = if row < height - 2 { Some(grid_clone.get_row(row + 2).unwrap()) } else { None };
            
            let rows = [prev2_row, prev_row, current_row, next_row, next2_row];
            
            let row_data = calc_row_lsp(rows, raster_params, &*selected_params_clone, &*basis_functions_clone, &*m_clone);
            tx.send((row, row_data)).unwrap();

            let mut num = count_clone.lock().unwrap();
            *num += 1;
            print!("\r{}Calculating Land Surface Parameters... {:.0}%", "\x1B[K", *num as f32/height as f32 * 100.0);
            io::stdout().flush().unwrap();
          }
        });
        handles.push(handle);
      }
      // Wait for all threads to finish
      for handle in handles {
        handle.join().unwrap();
      }
      
      let selected_count = selected_params.len();
      let mut param_grids = Vec::with_capacity(selected_count);
      
      for _ in &selected_params {
        let grid = raster::Grid::<f32>::from_params(raster_params, raster_params.nodata as f32);
        param_grids.push(grid);
      }
      
      // Receive processed row data from the threads and update the parameter grids accordingly.
      // For each row, the calculated Land Surface Parameters are applied to the corresponding grid.
      for _ in 0..height {
        
        let (row, row_data) = rx.recv().unwrap();
        for (col, lsp) in row_data {
          for (index, _) in selected_params.iter().enumerate() {
            match lsp {
              Some(ref lsp) => {
                if lsp.len() != selected_count {
                  panic!("Return parameters count mismatch!");
                }
                param_grids[index].set_item(row, col, lsp[index]);
              }
              None => {
                //param_grids[index].set_item(row, col, None);
              }
            }
          }
        }
      }
      
      let elapsed_time = part_time.elapsed();
      println!("\r\x1B[K{} Land Surface Parameters calculated in {:.2} seconds.", text::success("✓".to_string()), elapsed_time.as_secs_f64());
      part_time = Instant::now();
      
      // For each selected parameter, create an output GeoTIFF file.
      // The projection and geotransform from the input dataset are applied to each output.
      // The calculated values are written into the raster band of the output files.
      for (index, name) in selected_params.iter().enumerate() {
        
        let grid = &param_grids[index];
        let result = (*grid).to_gdal_buffer();
        match result {
          Ok(buffer) => {
            // output
            let driver = DriverManager::get_driver_by_name("GTiff").unwrap();
            let options = [
              RasterCreationOption {
                key: "TILED",
                value: "YES"
              },
              RasterCreationOption {
                key: "BLOCKXSIZE",
                value: "16"
              },
              RasterCreationOption {
                key: "BLOCKYSIZE",
                value: "16"
              }
            ];
            
            let file_name = format!("{}_{}.tif", output_prefix, name);
            
            let mut ds = driver
            .create_with_band_type_with_options::<f32, _>(
              &file_name,
              width as isize,
              height as isize,
              1,
              &options,
            ).expect("Failed to create GeoTIFF file.");
            
            ds.set_projection(&projection)?;
            ds.set_geo_transform(&geotransform).expect("Failed to set geotransform");
            
            let mut band1 = ds.rasterband(1).unwrap();
            band1.set_no_data_value(no_data)?;
            band1.write((0, 0), (width as usize, height as usize), &buffer).expect("Failed to write data");
            band1.compute_raster_min_max(true)?;
          }
          Err(err) => {
            println!("Failed to create output raster: {}", err);
          }
        }
      }
      let elapsed_time = part_time.elapsed();
      println!("{} Output files written in {:.2} seconds.", text::success("✓".to_string()), elapsed_time.as_secs_f64());
    }
    Err(err) => {
      println!("Chyba: {}", err);
    }
  }
  
  let elapsed_time = start_time.elapsed();
  println!("{}", line);
  println!("{}", text::success("Calculations completed successfully.".to_string()));
  println!("Total elapsed time: {:.2} seconds.", elapsed_time.as_secs_f64());
  println!("");
  
  Ok(())
}
