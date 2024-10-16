use std::f64::consts::PI;

const RAD: f64 = PI / 180.0;  // Constant for converting degrees to radians
const E: f64 = 0.08181919;    // Earth's eccentricity (WGS-84 model)

// Function to calculate the distance corresponding to one degree of longitude at a given latitude.
// The formula uses the cosine of latitude in radians.
pub fn get_degree_lon( lat: f64 ) -> f64 {

  let reduce_lat = ((lat * RAD).tan() * (1.0 - E * E).sqrt()).atan();

  111319.4908 * (reduce_lat).cos()
}

// Function to calculate the distance corresponding to one degree of latitude at a given latitude.
// A small adjustment is made to account for Earth's slight oblateness (flattening at the poles).
// The formula uses a cosine series to introduce minor corrections based on latitude in radians.
pub fn get_degree_lat( lat: f64 ) -> f64 {

  111132.09 - 0.55982 * (2.0 * lat * RAD).cos() + 0.00117 * (4.0 * lat * RAD).cos()
}
