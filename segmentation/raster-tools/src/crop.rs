use gdal::{Dataset, DriverManager};
use gdal::raster::{Buffer, ResampleAlg};
use std::path::PathBuf;

pub fn crop_raster(input: PathBuf, output: PathBuf, margin: usize) -> gdal::errors::Result<()> {
    let dataset = Dataset::open(&input)?;
    let band = dataset.rasterband(1)?;
    let (width, height) = band.size();
    let nodata = band.no_data_value().unwrap_or(-9999.0);

    let x_off = margin as isize;
    let y_off = margin as isize;
    let new_width = width - 2 * margin;
    let new_height = height - 2 * margin;

    let mut buffer: Buffer<f32> = band.read_as(
        (x_off, y_off),
        (new_width, new_height),
        (new_width, new_height),
        Some(ResampleAlg::NearestNeighbour),
    )?;

    let driver = DriverManager::get_driver_by_name("GTiff")?;
    let mut out_ds = driver.create_with_band_type::<f32, _>(&output, new_width, new_height, 1)?;
    let mut gt = dataset.geo_transform()?;
    gt[0] += x_off as f64 * gt[1]; // shift origin X
    gt[3] += y_off as f64 * gt[5]; // shift origin Y 
    out_ds.set_geo_transform(&gt)?;

    // Copy projection
    out_ds.set_projection(&dataset.projection())?;

    // Write data
    let mut out_band = out_ds.rasterband(1)?;
    out_band.write((0, 0), (new_width, new_height), &mut buffer)?;

    out_band.set_no_data_value(Some(nodata))?;
    // Explicit flush (optional, ensures disk write now)
    out_ds.flush_cache()?;

    Ok(())
}
