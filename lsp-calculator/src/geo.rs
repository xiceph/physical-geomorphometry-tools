use std::f64::consts::PI;

const RAD: f64 = PI / 180.0;  // Constant for converting degrees to radians

// Function to calculate the distance corresponding to one degree of longitude at a given latitude.
// The formula uses the cosine of latitude in radians.
pub fn get_degree_lon( lat: f64 ) -> f64 {
  111319.4908 * (lat * RAD).cos()
}

// Function to calculate the distance corresponding to one degree of latitude at a given latitude.
// A small adjustment is made to account for Earth's slight oblateness (flattening at the poles).
// The formula uses a cosine series to introduce minor corrections based on latitude in radians.
pub fn get_degree_lat( lat: f64 ) -> f64 {
  111132.09 - 0.55982 * (2.0 * lat * RAD).cos() + 0.00117 * (4.0 * lat * RAD).cos()
}
