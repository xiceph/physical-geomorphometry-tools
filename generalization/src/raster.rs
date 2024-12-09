use gdal::raster::{RasterBand, ResampleAlg, Buffer};
use crate::qem::SymmetricMatrix;

// Trait to check if a value represents nodata
pub trait NodataComparable {
  fn is_nodata(&self, nodata: f64) -> bool;
}

// Implement NodataComparable for f64, f32, and SymmetricMatrix
impl NodataComparable for f64 {
  fn is_nodata(&self, nodata: f64) -> bool {
    *self == nodata
  }
}
impl NodataComparable for f32 {
  fn is_nodata(&self, nodata: f64) -> bool {
    *self as f64 == nodata
  }
}
impl NodataComparable for SymmetricMatrix {
  fn is_nodata(&self, _nodata: f64) -> bool {
    false
  }
}

// Raster parameters with width, height, origin, resolution, and nodata value
#[derive(Debug, Clone, Copy)]
pub struct RasterParams {
  pub width: usize,
  pub height: usize,
  pub origin: [f64; 2],
  pub resolution: [f64; 2],
  pub nodata: f64,
}

impl RasterParams {
  // Get XY coordinates for a given row and column
  pub fn get_xy_coords(&self, row: usize, col: usize) -> [f64; 2] {
    if row >= self.height || col >= self.width {
      println!("Varovanie: Zadané súradnice mimo rozsahu rastra.");
    }

    let center: [f64; 2] = [
      (self.origin[0] + self.width as f64 * 0.5 * self.resolution[0]).floor(),
      (self.origin[1] + self.height as f64 * 0.5 * self.resolution[1]).floor(),
    ];

    let x = self.origin[0] + col as f64 * self.resolution[0] - center[0];
    let y = self.origin[1] + row as f64 * self.resolution[1] - center[1];

    [x, y]
  }
  // Check if a value is nodata
  pub fn is_nodata(&self, value: f64) -> bool {
    value == self.nodata
  }
}

// Two-dimensional grid with raster parameters and data
pub struct Grid<T> {
  params: RasterParams,
  data: Vec<Vec<T>>,
}

impl<T: Clone> Clone for Grid<T> {
  fn clone(&self) -> Self {
    Grid {
      params: self.params.clone(),
      data: self.data.clone(),
    }
  }
}

impl<T: PartialEq + Clone + NodataComparable> Grid<T> {

  // Create a new 2D grid with RasterParams and an initial value
  pub fn from_params(params: RasterParams, initial_value: T) -> Self {
    let mut data = Vec::with_capacity(params.height);
    for _ in 0..params.height {
      let row = vec![initial_value.clone(); params.width];
      data.push(row);
    }
    Grid { params, data }
  }

  // Set value at a specific row and column
  pub fn set_item(&mut self, row_index: usize, col_index: usize, value: T) {
    if let Some(row) = self.data.get_mut(row_index) {
      if let Some(cell) = row.get_mut(col_index) {
        *cell = value;
      }
    }
  }

  // Get a reference to a row in the grid
  pub fn get_row(&self, row_index: usize) -> Option<&Vec<T>> {
    self.data.get(row_index)
  }

  // Get grid parameters
  pub fn params(&self) -> &RasterParams {
    &self.params
  }

}

impl Grid<f32> {

  // Load data from a GDAL RasterBand
  pub fn from_raster_band(params: RasterParams, rasterband: &RasterBand) -> Result<Self, String> {

    let (width, height) = rasterband.size();
    let offset = (0, 0);
    let buffer_size = (params.width, params.height);

    let data = rasterband
      .read_as::<f32>(offset, (width, height), buffer_size, Some(ResampleAlg::Lanczos))
      .map_err(|e| format!("Failed to read data: {}", e))
      .expect("Failed to read data");

    // Initialize Grid with data read from RasterBand
    let mut grid_data = Vec::with_capacity(params.height);
    for row in 0..params.height {
      let mut row_data = Vec::with_capacity(params.width);
      for col in 0..params.width {
        let index = (row * params.width + col) as usize;
        row_data.push(data.data()[index]);
      }
      grid_data.push(row_data);
    }

    Ok(Grid {
      params,
      data: grid_data,
    })
  }

  // Convert the grid to a GDAL buffer
  pub fn to_gdal_buffer(&self) -> Result<Buffer<f32>, String> {

    let size = (self.params.width, self.params.height);

    let pixels = size.0 * size.1;
    let mut geotiff_data = Vec::with_capacity(pixels);

    for row in &self.data {
      for value in row {
        geotiff_data.push(*value);
      }
    }

    Ok(Buffer::new(size, geotiff_data))
  }
}
