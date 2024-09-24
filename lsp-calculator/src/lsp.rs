use std::f64::consts::PI;

const DEG: f64 = 180.0 / PI;  // Constant for converting radians to degrees

// Function to calculate land surface parameters based on derivatives
pub fn calculate(
  param: &str,                // Name of the parameter to calculate
  derivatives: &Vec<f64>,     // A vector containing the derivative values
) -> f32 {
  
  let value: f64;  // Variable to store the computed value
  
  if derivatives.len() == 10 {
    // Destructure the derivatives for easy access
    let [_z, zx, zy, zxx, zxy, zyy, zxxx, zxxy, zxyy, zyyy] = &derivatives[..] else {
      panic!("Bad derivatives data");  
    };
    
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

    let angle = (zy).atan2(*zx);
    let sin_aspect = angle.sin();
    let cos_aspect = angle.cos();

    let kns_numerator = -1.0 * (zxx * zx2 + 2.0 * zxy * zx_zy + zyy * zy2);
    let kns = kns_numerator / (m * p.powf(1.5));

    let knc_numerator = -1.0 * (zxx * zy2 - 2.0 * zxy * zx_zy + zyy * zx2);
    let knc = knc_numerator / (m * p.sqrt());

    let tgc_numerator = zx_zy * (zxx - zyy) - zxy * (zx2 - zy2);
    let tgc = tgc_numerator / (m * p);

    let numerator1 = (1.0 + zy2) * zxx - 2.0 * zxy * zx_zy + (1.0 + zx2) * zyy;
    let term1 = numerator1 / (2.0 * p.powf(1.5));
    let term2 = zxx * zyy - zxy2 / p.powi(2);
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



    // Match to determine the parameter to calculate
    match param {
      "slope" => {
        // Calculate the Slope in degrees
        value = slope * DEG;
      },
      "sin_slope" => {
        // Use sine of slope
        value = sin_slope;
      }
      "aspect" => {
        // Calculate the Aspect
        let angle_deg = angle * DEG;
        value = if angle_deg < 0.0 {
          angle_deg + 360.0
        } else {
          angle_deg
        };
      },
      "sin_aspect" => {
        // Use the Sine of Aspect
        value = sin_aspect;
      },
      "cos_aspect" => {
        // Use the Cosine of Aspect
        value = cos_aspect;
      },
      "kns" => {
        // Use the Normal slope line (profile) curvature
        value = kns;
      },
      "zss" => {
        // Calculate Second slope line derivative
        value = kns / (cos_slope2 * cos_slope);
      },
      "ts" => {
        // Calculate Slope line torsion
        value = kns / cos_slope;
      },
      "knc" => {
        // Use the Normal contour (tangential) curvature
        value = knc;
      },
      "zcc" => {
        // Calculate Second contour derivative
        value = knc / cos_slope;
      },
      "kpc" => {
        // Calculate Projected contour curvature
        value = knc / sin_slope;
      },
      "tgc" => {
        // Use the Contour geodesic torsion
        value = tgc;
      },
      "kps" => {
        // Calculate Projected slope line curvature
        value = tgc / (sin_slope * cos_slope);
      },
      "sin_sc" => {
        // Calculate Contour change of sin slope
        value = tgc * cos_slope;
      },
      "el" => {
        // Calculate Elevation laplacian
        value = (knc * cos_slope2 + kns) / (cos_slope2 * cos_slope);
      },
      "k_max" => {
        // Use Maximal curvature
        value = k_max;
      },
      "k_min" => {
        // Use Minimal curvature
        value = k_min;
      },
      "k_mean" => {
        // Calculate Mean curvature
        value = (kns + knc) / 2.0;
      },
      "kc" => {
        // Calculate Casorati curvature
        value = ((k_max.powi(2) + k_min2) / 2.0).sqrt();
      },
      "kr" => {
        // Calculate Total ring curvature
        value = knc * kns - k_min * (knc + kns) + k_min2;
      },
      "kd" => {
        // Calculate Difference curvature
        value = (kns - knc) / 2.0;
      },
      "ka" => {
        // Calculate Total accumulation curvature
        value = kns * knc;
      },
      "khe" => {
        // Calculate Horizontal excess curvature
        value = knc - k_min;
      },
      "kve" => {
        // Calculate Vertical excess curvature
        value = kns - k_min;
      },
      "k" => {
        // Calculate Gaussian curvature
        value = k_max * k_min;
      },
      "ku" => {
        // Calculate Unsphericity curvature
        value = (k_max - k_min) / 2.0;
      },
      "kncc" => {
        // Use Contour change of normal contour curvature
        value = kncc;
      },
      "kncs" => {
        // Use Slope line change of normal contour curvature
        value = kncs;
      },
      "knss" => {
        // Use Slope line change of normal slope line curvature
        value = knss;
      },
      _ => panic!("Parameter calculation for {} is not implemented!", param),  
    }
  } else {
    panic!("Bad derivatives length");  
  }
  
  value as f32 
}

