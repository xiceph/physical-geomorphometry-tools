#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymmetricMatrix {
  pub m: [f64; 10], // Array to store matrix coefficients
}

impl SymmetricMatrix {
  // Initialize a new SymmetricMatrix with all coefficients set to 0
  pub fn new() -> Self {
    SymmetricMatrix { m: [0.0; 10] }
  }

  // Set values for the 10 independent elements of a 4x4 symmetric matrix
  pub fn set(&mut self, m11: f64, m12: f64, m13: f64, m14: f64, m22: f64, m23: f64, m24: f64, m33: f64, m34: f64, m44: f64) -> &mut Self {
    self.m[0] = m11;
    self.m[1] = m12;
    self.m[2] = m13;
    self.m[3] = m14;
    self.m[4] = m22;
    self.m[5] = m23;
    self.m[6] = m24;
    self.m[7] = m33;
    self.m[8] = m34;
    self.m[9] = m44;
    self
  }

  // Construct a symmetric matrix representing a plane ax + by + cz + d = 0
  pub fn make_plane(&mut self, a: f64, b: f64, c: f64, d: f64) -> &mut Self {
    self.set(
      a * a, a * b, a * c, a * d,
      b * b, b * c, b * d,
      c * c, c * d,
      d * d
    )
  }

  // Add another SymmetricMatrix to the current one (element-wise)
  pub fn add_self(&mut self, n: &SymmetricMatrix) {
    for i in 0..10 {
      self.m[i] += n.m[i];
    }
  }

  // Divide all elements of the matrix by a scalar value
  pub fn divide(&mut self, c: f64) {
    for item in &mut self.m {
      *item /= c;
    }
  }
}

// Calculate the error at a vertex (x, y, z) using the matrix Q
pub fn calc_vertex_error(q: &SymmetricMatrix, x: f64, y: f64, z: f64) -> f64 {
  let m = &q.m;
  m[0] * x * x + 2.0 * m[1] * x * y + 2.0 * m[2] * x * z + 2.0 * m[3] * x +
  m[4] * y * y + 2.0 * m[5] * y * z + 2.0 * m[6] * y +
  m[7] * z * z + 2.0 * m[8] * z + m[9]
}

// Calculate the optimal z-coordinate for a grid node (x, y) using the matrix Q
pub fn calc_optimal_z(q: &SymmetricMatrix, x: f64, y: f64) -> f64 {
  let m = &q.m;
  -((m[2]) * x + (m[5]) * y + (m[8])) / (m[7])
}


