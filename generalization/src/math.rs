#[derive(Debug, Clone, Copy)]
pub struct Vector3 {
  pub x: f64,
  pub y: f64,
  pub z: f64,
}

// Helper function to create a Vector3 instance
fn vector3(x: f64, y: f64, z: f64) -> Vector3 {
  Vector3 { x, y, z }
}

// Calculate the parameters of a plane (a, b, c, d) passing through three points in 3D space
pub fn calc_plane_params(point_a: Vector3, point_b: Vector3, point_c: Vector3) -> [f64; 4] {

  // Compute vectors
  let vector_ab = vector3(
    point_b.x - point_a.x,
    point_b.y - point_a.y,
    point_b.z - point_a.z,
  );
  let vector_ac = vector3(
    point_c.x - point_a.x,
    point_c.y - point_a.y,
    point_c.z - point_a.z,
  );

  // Compute the cross product
  let normal_vector = vector3(
    vector_ab.y * vector_ac.z - vector_ab.z * vector_ac.y,
    vector_ab.z * vector_ac.x - vector_ab.x * vector_ac.z,
    vector_ab.x * vector_ac.y - vector_ab.y * vector_ac.x,
  );

  // Calculate the squared length of the normal vector
  let length_squared = normal_vector.x * normal_vector.x
    + normal_vector.y * normal_vector.y
    + normal_vector.z * normal_vector.z;

  // Normalize the normal vector to ensure its length is 1
  let normalized_normal_vector = if length_squared != 0.0 {
    let length = length_squared.sqrt();
    vector3(
      normal_vector.x / length,
      normal_vector.y / length,
      normal_vector.z / length,
    )
  } else {
     // Handle the case where the points are collinear (normal vector length is zero)
    panic!("Division by zero encountered when calculating normalized normal vector");
  };

  // Extract components of the normalized normal vector
  let a = normalized_normal_vector.x;
  let b = normalized_normal_vector.y;
  let c = normalized_normal_vector.z;

  // Compute the plane constant `d` using point A
  let d = -a * point_a.x - b * point_a.y - c * point_a.z;

  [a, b, c, d] // Return the plane equation coefficients
}

