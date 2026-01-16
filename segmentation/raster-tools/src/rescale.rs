use gdal::{Dataset, DriverManager};
use gdal::raster::Buffer;
use std::path::PathBuf;

pub fn rescale_raster(input: PathBuf, output: PathBuf) -> gdal::errors::Result<()> {
    let dataset = Dataset::open(&input)?;
    let band = dataset.rasterband(1)?;
    let (width, height) = band.size();

    // Use existing NoData or fallback to -9999
    let no_data_value = band.no_data_value().unwrap_or(-9999.0);

    // Read raster into buffer
    let buffer: Buffer<f32> = band.read_as(
        (0, 0),
        (width, height),
        (width, height),
        None,
    )?;

    // Collect valid values
    let valid_values: Vec<f32> = buffer
        .data()
        .iter()
        .cloned()
        .filter(|v| !v.is_nan() && (*v as f64 - no_data_value).abs() > 1e-9)
        .collect();

    if valid_values.is_empty() {
        panic!("✗ No valid data found for rescaling!");
    }

    let min_val = valid_values.iter().cloned().fold(f32::MAX, f32::min);
    let max_val = valid_values.iter().cloned().fold(f32::MIN, f32::max);

    println!("ⓘ Rescale range: min = {}, max = {}", min_val, max_val);

    // Handle flat raster (min == max)
    let mut buffer_f32: Buffer<f32> = Buffer::new((width, height), vec![0f32; width * height]);
    if (max_val - min_val).abs() < f32::EPSILON {
        println!("‼ Flat raster detected (min = max). Filling with 128.0.");
        buffer_f32.data_mut().iter_mut().for_each(|v| *v = 128.0);
    } else {
        for (i, &val) in buffer.data().iter().enumerate() {
            if val.is_nan() || (val as f64 - no_data_value).abs() < 1e-9 {
                buffer_f32.data_mut()[i] = no_data_value as f32; // keep NoData
            } else {
                buffer_f32.data_mut()[i] = (((val - min_val) / (max_val - min_val)) * 255.0).clamp(0.0, 255.0);
            }
        }
    }

    // Save output raster
    let driver = DriverManager::get_driver_by_name("GTiff")?;
    let mut out_ds = driver.create_with_band_type::<f32, _>(&output, width, height, 1)?;
    out_ds.set_projection(&dataset.projection())?;
    out_ds.set_geo_transform(&dataset.geo_transform()?)?;

    let mut out_band = out_ds.rasterband(1)?;
    out_band.write((0, 0), (width, height), &mut buffer_f32)?;
    out_band.set_no_data_value(Some(no_data_value))?;

    Ok(())
}
