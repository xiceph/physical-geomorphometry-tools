use gdal::raster::{RasterBand, ResampleAlg, Buffer};

pub trait NodataComparable {
  fn is_nodata(&self, nodata: f64) -> bool;
}

impl NodataComparable for f32 {
  fn is_nodata(&self, nodata: f64) -> bool {
    *self as f64 == nodata
  }
}

#[derive(Debug, Clone, Copy)]
pub struct RasterParams {
  pub width: usize,
  pub height: usize,
  pub origin: [f64; 2],
  pub resolution: [f64; 2],
  pub nodata: f64,
}

impl RasterParams {
  // Method to check if a value is a NoData value
  pub fn is_nodata(&self, value: f64) -> bool {
    value == self.nodata
  }
}

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

  // Create a new two-dimensional grid from RasterParams with a specified initial value
  pub fn from_params(params: RasterParams, initial_value: T) -> Self {
    let mut data = Vec::with_capacity(params.height);
    for _ in 0..params.height {
      let row = vec![initial_value.clone(); params.width];
      data.push(row);
    }
    Grid { params, data }
  }

  // Set the value at a specific row and column
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

  // Retrieve the grid parameters
  pub fn params(&self) -> &RasterParams {
    &self.params
  }

}

impl Grid<f32> {

  // Load data from Gdal RasterBand
  pub fn from_raster_band(params: RasterParams, rasterband: &RasterBand) -> Result<Self, String> {

    let (width, height) = rasterband.size();
    let offset = (0, 0);
    let buffer_size = (params.width, params.height);

    let data = rasterband
    .read_as::<f32>(offset, (width, height), buffer_size, Some(ResampleAlg::NearestNeighbour))
    .map_err(|e| format!("Failed to read data: {}", e))
    .expect("Failed to read data");

    // Initialize Grid with data read from rasterband
    let mut grid_data = Vec::with_capacity(params.height);
    for row in 0..params.height {
      let mut row_data = Vec::with_capacity(params.width);
      for col in 0..params.width {
        let index = (row * params.width + col) as usize;
        row_data.push(data.data[index]);
      }
      grid_data.push(row_data);
    }

    Ok(Grid {
      params,
        data: grid_data,
    })
  }

  // Convert grid data to GDAL buffer
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

