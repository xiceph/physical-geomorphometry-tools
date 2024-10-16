extern crate nalgebra as na;
use na::{DMatrix, DVector};
use std::f64::EPSILON;

// Function to precompute basis functions matrix with resolutions for x and y
pub fn precompute_basis_functions(degree: usize, d_x: f64, d_y: f64) -> DMatrix<f64> {
  let grid_coords_x = [-2.0 * d_x, -1.0 * d_x, 0.0, 1.0 * d_x, 2.0 * d_x];
  let grid_coords_y = [-2.0 * d_y, -1.0 * d_y, 0.0, 1.0 * d_y, 2.0 * d_y];
  let num_points = grid_coords_x.len() * grid_coords_y.len();  // 25 points (5x5 grid)
  let num_basis_functions = (degree + 1) * (degree + 2) / 2;   // Number of basis functions
  
  // Initialize the basis functions matrix with zeros (size: 25 x num_basis_functions)
  let mut basis_functions = DMatrix::zeros(num_points, num_basis_functions);
  
  // Loop through grid coordinates for y and x
  for (y_idx, &y) in grid_coords_y.iter().enumerate() {
    for (x_idx, &x) in grid_coords_x.iter().enumerate() {
      let idx = y_idx * 5 + x_idx;  // Index for the point in the grid

      // Compute basis functions based on polynomial degree
      let mut func_idx = 0;
      for sum in 0..=degree {
        for py in 0..=sum {
          let px = sum - py;
          // Compute x^px * y^py for each basis function
          basis_functions[(idx, func_idx)] = x.powi(px as i32) * y.powi(py as i32);
          func_idx += 1;
        }
      }
    }
  }
  
  basis_functions
}


// Function to compute matrix M
pub fn compute_matrix_m(basis_functions: &DMatrix<f64>) -> DMatrix<f64> {
  let num_basis_functions = basis_functions.ncols();  // Number of basis functions
  
  // Matrix M will have dimensions (num_basis_functions x num_basis_functions)
  let mut m = DMatrix::zeros(num_basis_functions, num_basis_functions);
  
  // Fill matrix M by computing dot products of basis functions for each point (symmetric calculation)
  for i in 0..num_basis_functions {
    let v1 = basis_functions.column(i);
    for j in i..num_basis_functions {
      let v2 = basis_functions.column(j);
      let dot_product = v1.dot(&v2);
      m[(i, j)] = dot_product;
      m[(j, i)] = dot_product; // Matrix M is symmetric
    }
  }
  
  m
}

// Main function to compute derivatives using precomputed values
pub fn poly_lsq_der_with_precomputed(
  z: &[f64],                     // One-dimensional array of z values
  precomputed_bf: &DMatrix<f64>, // Precomputed basis functions (DMatrix)
  precomputed_m: &DMatrix<f64>,  // Precomputed matrix M (DMatrix)
) -> Option<Vec<f64>> {
  // Derive the number of derivatives directly from the matrix precomputed_bf
  let num_derivatives = precomputed_bf.ncols();
  
  // Check that the number of z points matches the number of rows in precomputed_bf
  let num_points = z.len();
  if num_points != precomputed_bf.nrows() {
    panic!("Array lengths do not match.");
  }
  
  // Ensure that the dimensions of precomputed_bf and precomputed_m are consistent
  if precomputed_bf.ncols() != precomputed_m.nrows() || precomputed_m.ncols() != precomputed_bf.ncols() {
    panic!("Dimension mismatch between precomputed basis functions and matrix M.");
  }
  
  // Definition of derivatives for a 3rd degree
  let deriv_x = vec![0, 1, 0, 2, 1, 0, 3, 2, 1, 0];
  let deriv_y = vec![0, 0, 1, 0, 1, 2, 0, 1, 2, 3];
  let mut derivatives = vec![0.0; 10];  // We want to compute only 10 derivatives
  
  // Vector v
  let mut v = DVector::zeros(num_derivatives);
  for i in 0..num_derivatives {
    v[i] = precomputed_bf.column(i).dot(&DVector::from_row_slice(z));
  }
  
  // Solve the system Mx = v
  let c = if let Some(cholesky) = precomputed_m.clone().cholesky() {
    cholesky.solve(&v)
  } else if is_matrix_well_conditioned(precomputed_m, EPSILON) {
    match precomputed_m.clone().lu().solve(&v) {
      Some(solution) => solution,
      None => return None,
    }
  } else {
    return None;
  };
  
  // Precomputed factorials for efficient derivative computation
  const FACTORIALS: [usize; 4] = [1, 1, 2, 6];
  
  // Compute derivatives for the first 10 values
  for k in 0..10 {
    let dx = deriv_x[k];
    let dy = deriv_y[k];
    derivatives[k] = c[k] * FACTORIALS[dx] as f64 * FACTORIALS[dy] as f64;
  }
  
  Some(derivatives)
}

// Function to check if the matrix is well-conditioned
fn is_matrix_well_conditioned(m: &DMatrix<f64>, tolerance: f64) -> bool {
  let det = m.determinant();
  det.abs() >= tolerance
}

