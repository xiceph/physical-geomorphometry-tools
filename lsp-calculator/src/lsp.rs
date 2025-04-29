use std::f64::consts::PI;

const DEG: f64 = 180.0 / PI;  // Constant for converting radians to degrees

// Function to calculate land surface parameters based on derivatives
pub fn calculate(
  params: &Vec<String>,         // A vector of parameter names to calculate
  derivatives: &Vec<f64>,     // A vector containing the derivative values
) -> Vec<f32> {
  
  // Initialize a vector to store the results
  let mut param_values = Vec::with_capacity(params.len());
  
  if derivatives.len() == 10 {
    // Destructure the derivatives for easy access
    let [_z, zx, zy, zxx, zxy, zyy, zxxx, zxxy, zxyy, zyyy] = &derivatives[..] else {
      panic!("Bad derivatives data");  
    };
    
    //TODO check needed calculations

    let zx2 = zx * zx;
    let zy2 = zy * zy;
    let zxy2 = zxy * zxy;
    let zx_zy = zx * zy;

    let a = zxy2 - zxx * zyy;
    let b = 2.0 * zx * (zxxy * zy + a) - zxxx * zy2 - zxyy * zx2;
    let c = 2.0 * zy * (zxyy * zx + a) - zxxy * zy2 - zyyy * zx2;
    let d = zxy * (zxx + zyy);
    let e = -2.0 * (zx * (zxxy * zy + zxy2 + zxx.powi(2)) + d * zy) - zxxx * zx2 - zxyy * zy2;
    let f = -2.0 * (zy * (zxyy * zx + zxy2 + zyy.powi(2)) + d * zx) - zyyy * zy2 - zxxy * zx2;
    let k = 2.0 * zxy * zx_zy - zxx * zy2 - zyy * zx2;
    let l = - zx2 * zxx - 2.0 * zxy * zx_zy - zy2 * zyy;
    let m = zx2 + zy2;
    let n = zxx * zx + zxy * zy;
    let o = zyy * zy + zxy * zx;
    let p = 1.0 + m;

    let norm = m.sqrt();
    let slope = norm.atan();
    let sin_slope = slope.sin();
    let cos_slope = slope.cos();
    let cos_slope2 = cos_slope * cos_slope;

    let angle = (-zx).atan2(-*zy);
    let sin_aspect = angle.sin();
    let cos_aspect = angle.cos();

    let kns_numerator = -1.0 * (zxx * zx2 + 2.0 * zxy * zx_zy + zyy * zy2);
    let kns = kns_numerator / (m * p.powf(1.5));

    let knc_numerator = -1.0 * (zxx * zy2 - 2.0 * zxy * zx_zy + zyy * zx2);
    let knc = knc_numerator / (m * p.sqrt());

    let tc_numerator = zx_zy * (zxx - zyy) - zxy * (zx2 - zy2);
    let tc = tc_numerator / (m * p);

    let numerator1 = (1.0 + zy2) * zxx - 2.0 * zxy * zx_zy + (1.0 + zx2) * zyy;
    let term1 = numerator1 / (2.0 * p.powf(1.5));
    let term2 = (zxx * zyy - zxy2) / p.powi(2);
    let sqrt = (term1.powi(2) - term2).sqrt();
    let k_max = - term1 + sqrt;
    let k_min = - term1 - sqrt;
    let k_min2 = k_min * k_min;

    let bmp = b * m * p;
    let npmn = 2.0 * n * p + m * n;
    let cmpk = c * m * p - k * (2.0 * o * p + m * o);
    let p3 = p.powi(3);
    let m5 = m.powi(5);

    let kncc_numerator = zy * (bmp - k * npmn) - zx * cmpk;
    let kncc_denominator = p3.sqrt() * m5.sqrt();
    let kncc = kncc_numerator / kncc_denominator;

    let kncs_numerator = zx * (k * npmn - bmp) - zy * cmpk;
    let kncs = kncs_numerator / kncc_denominator;

    let knss_numerator = p * (zx * (2.0 * l * n - e * m ) - zy * (f * m - 2.0 * l * o)) + 3.0 * l * m * (zy * o + zx * n);
    let knss_denominator = (p3 * p * p).sqrt() * m5.sqrt();
    let knss = knss_numerator / knss_denominator;

    for param in params {
      let value = match param.as_str() {
        "slope" => {
          // Calculate the Slope in degrees
          slope * DEG
        },
        "sin_slope" => {
          // Use sine of slope
          sin_slope
        }
        "aspect" => {
          // Calculate the Aspect
          let angle_deg = angle * DEG;
          if angle_deg < 0.0 {
            angle_deg + 360.0
          } else {
            angle_deg
          }
        },
        "sin_aspect" => {
          // Use the Sine of Aspect
          sin_aspect
        },
        "cos_aspect" => {
          // Use the Cosine of Aspect
          cos_aspect
        },
        "kns" => {
          // Use the Normal slope line (profile) curvature
          kns
        },
        "zss" => {
          // Calculate Second slope line derivative
          kns / (cos_slope2 * cos_slope)
        },
        "ts" => {
          // Calculate Slope line torsion
          kns / cos_slope
        },
        "knc" => {
          // Use the Normal contour (tangential) curvature
          knc
        },
        "zcc" => {
          // Calculate Second contour derivative
          knc / cos_slope
        },
        "kpc" => {
          // Calculate Projected contour curvature
          knc / sin_slope
        },
        "tc" => {
          // Use the Contour geodesic torsion
          tc
        },
        "kps" => {
          // Calculate Projected slope line curvature
          tc / (sin_slope * cos_slope)
        },
        "sin_sc" => {
          // Calculate Contour change of sin slope
          tc * cos_slope
        },
        "el" => {
          // Calculate Elevation laplacian
          (knc * cos_slope2 + kns) / (cos_slope2 * cos_slope)
        },
        "kmax" => {
          // Use Maximal curvature
          k_max
        },
        "kmin" => {
          // Use Minimal curvature
          k_min
        },
        "kmean" => {
          // Calculate Mean curvature
          (kns + knc) / 2.0
        },
        "kc" => {
          // Calculate Casorati curvature
          ((k_max.powi(2) + k_min2) / 2.0).sqrt()
        },
        "kr" => {
          // Calculate Total ring curvature
          knc * kns - k_min * (knc + kns) + k_min2
        },
        "kd" => {
          // Calculate Difference curvature
          (kns - knc) / 2.0
        },
        "ka" => {
          // Calculate Total accumulation curvature
          kns * knc
        },
        "khe" => {
          // Calculate Horizontal excess curvature
          knc - k_min
        },
        "kve" => {
          // Calculate Vertical excess curvature
          kns - k_min
        },
        "k" => {
          // Calculate Gaussian curvature
          k_max * k_min
        },
        "ku" => {
          // Calculate Unsphericity curvature
          (k_max - k_min) / 2.0
        },
        "kncc" => {
          // Use Contour change of normal contour curvature
          kncc
        },
        "kncs" => {
          // Use Slope line change of normal contour curvature
          kncs
        },
        "knss" => {
          // Use Slope line change of normal slope line curvature
          knss
        },
        _ => panic!("Parameter calculation for {} is not implemented!", param),
      };

      param_values.push(value as f32);
    }

  } else {
    panic!("Bad derivatives length");  
  }
  
  param_values
}

