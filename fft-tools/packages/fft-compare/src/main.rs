use anyhow::{Context, Result};
use charming::{
    component::{Axis as ChartAxis, Legend, Title},
    element::{AxisType, Tooltip, Trigger, TriggerOn},
    series::Line,
    Chart, HtmlRenderer,
};
use clap::Parser;
use console::Term;
use fft_core::text;
use gdal::Dataset;
use ndarray::{Array1, Array2};
use num_complex::Complex;
use rayon::prelude::*;
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Compares two sets of FFT analysis results to quantify information loss and spectral differences."
)]
struct Args {
    /// Reference directory (Input A) containing fft-process outputs.
    #[arg(long)]
    input_a: PathBuf,

    /// Comparison directory (Input B) containing fft-process outputs.
    #[arg(long)]
    input_b: PathBuf,

    /// Output directory for comparison results.
    #[arg(long)]
    output: PathBuf,

    /// Target power retention percentage (e.g., 50.0).
    #[arg(long, default_value = "50.0")]
    retention_threshold: f64,

    /// Target coherence threshold (e.g., 0.5).
    #[arg(long, default_value = "0.5")]
    coherence_threshold: f64,

    /// Save partial results (CSV) for each block.
    #[arg(long, default_value_t = false)]
    save_partials: bool,

    /// Path to save the HTML plot.
    #[arg(long)]
    plot: Option<PathBuf>,

    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    jobs: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct BlockCoord {
    row: usize,
    col: usize,
}

struct BlockData {
    psd: Array2<f64>,
    complex: Option<Array2<Complex<f64>>>,
}

#[derive(Clone)]
struct BlockProfiles {
    ratio: Array1<f64>,
}

struct GlobalAccumulator {
    block_profiles: Vec<BlockProfiles>,
    cross_psd_total: Array2<Complex<f64>>,
    psd_a_total: Array2<f64>,
    psd_b_total: Array2<f64>,
}

fn find_block_coords(dir: &Path) -> Result<Vec<BlockCoord>> {
    let mut coords = Vec::new();
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("fft_metadata_block_") && name.ends_with(".json") {
            let parts: Vec<&str> = name
                .trim_start_matches("fft_metadata_block_")
                .trim_end_matches(".json")
                .split('_')
                .collect();
            if parts.len() == 2 {
                coords.push(BlockCoord {
                    row: parts[0].parse()?,
                    col: parts[1].parse()?,
                });
            }
        }
    }
    coords.sort();
    Ok(coords)
}

fn find_first_metadata(dir: &Path) -> Result<PathBuf> {
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_file()
            && path
                .file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("fft_metadata_block_")
        {
            return Ok(path);
        }
    }
    anyhow::bail!("No metadata files found.")
}

fn load_psd(path: &Path) -> Result<Array2<f64>> {
    let dataset = Dataset::open(path)?;
    let band = dataset.rasterband(1)?;
    let (cols, rows) = band.size();
    let power_log = band
        .read_as::<f64>((0, 0), (cols, rows), (cols, rows), None)?
        .data()
        .to_vec();
    Ok(Array2::from_shape_vec((rows, cols), power_log)?.mapv(|p| 10.0_f64.powf(p) - 1e-12))
}

fn load_complex(path: &Path, rows: usize, cols: usize) -> Result<Array2<Complex<f64>>> {
    let mut file = fs::File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let mut data = Vec::with_capacity(rows * cols);
    for chunk in buffer.chunks_exact(16) {
        let re = f64::from_le_bytes(chunk[0..8].try_into()?);
        let im = f64::from_le_bytes(chunk[8..16].try_into()?);
        data.push(Complex::new(re, im));
    }
    Array2::from_shape_vec((rows, cols), data).context("Reshape failed")
}

fn load_block(dir: &Path, coord: &BlockCoord, load_complex_data: bool) -> Result<BlockData> {
    let psd = load_psd(&dir.join(format!("fft_psd_block_{}_{}.tif", coord.row, coord.col)))?;
    let mut complex = None;
    if load_complex_data {
        let (rows, cols) = psd.dim();
        complex = Some(load_complex(
            &dir.join(format!("fft_complex_block_{}_{}.bin", coord.row, coord.col)),
            rows,
            cols,
        )?);
    }
    Ok(BlockData { psd, complex })
}

fn load_metadata_info(meta_path: &Path) -> Result<(Array1<f64>, Array1<f64>, f64, usize)> {
    let metadata: BTreeMap<String, Value> = serde_json::from_reader(fs::File::open(meta_path)?)?;
    let stats = &metadata["statistics"];
    let pixel_size = 1.0 / (2.0 * stats["f_nyquist"].as_f64().unwrap());
    let (p_rows, p_cols) = (
        metadata["padded_size"][0].as_u64().unwrap() as usize,
        metadata["padded_size"][1].as_u64().unwrap() as usize,
    );
    let window_size = metadata["original_size"][0].as_u64().unwrap() as usize;
    let mut fx = fft_core::fftfreq(p_cols, pixel_size);
    let mut fy = fft_core::fftfreq(p_rows, pixel_size);
    fft_core::fftshift_1d(&mut fx);
    fft_core::fftshift_1d(&mut fy);
    Ok((Array1::from(fx), Array1::from(fy), pixel_size, window_size))
}

fn compute_radial_mean(
    psd: &Array2<f64>,
    fx: &Array1<f64>,
    fy: &Array1<f64>,
    n_bins: usize,
    k_min_limit: f64,
) -> (Array1<f64>, Array1<f64>) {
    let (rows, cols) = psd.dim();
    let k_max = (fx.iter().cloned().fold(f64::NEG_INFINITY, f64::max).powi(2)
        + fy.iter().cloned().fold(f64::NEG_INFINITY, f64::max).powi(2))
    .sqrt();
    let k_min = (fx[1] - fx[0]).abs().min((fy[1] - fy[0]).abs());
    let k_bins = Array1::logspace(10.0, k_min.log10(), k_max.log10(), n_bins + 1);
    let mut sums = Array1::<f64>::zeros(n_bins);
    let mut counts = Array1::<f64>::zeros(n_bins);
    for r in 0..rows {
        for c in 0..cols {
            let k = (fx[c].powi(2) + fy[r].powi(2)).sqrt();
            if k < k_min || k > k_max {
                continue;
            }
            let idx = match k_bins
                .as_slice()
                .unwrap()
                .binary_search_by(|v| v.partial_cmp(&k).unwrap())
            {
                Ok(i) => i.min(n_bins - 1),
                Err(i) => i.saturating_sub(1).min(n_bins - 1),
            };
            sums[idx] += psd[[r, c]];
            counts[idx] += 1.0;
        }
    }
    let (mut rad, mut waves) = (Array1::<f64>::zeros(n_bins), Array1::<f64>::zeros(n_bins));
    for i in 0..n_bins {
        let k_mid = (k_bins[i] + k_bins[i + 1]) / 2.0;
        rad[i] = if counts[i] > 0.0 && k_mid >= k_min_limit {
            sums[i] / counts[i]
        } else {
            f64::NAN
        };
        waves[i] = 1.0 / k_mid;
    }
    (waves, rad)
}

fn main() -> Result<()> {
    let args = Args::parse();

    let line = "-".repeat(72);
    let dline = "=".repeat(72);

    println!(
        "\n{}\n{}\nTool for comparing two sets of FFT results to quantify differences.\nPart of the {} toolkit.\n\nAuthors:\n{}\n{}\n",
        format!(
            "{} {}",
            text::highlight("FFT Comparator"),
            env!("CARGO_PKG_VERSION")
        ),
        line,
        text::highlight("fft-tools"),
        env!("CARGO_PKG_AUTHORS"),
        dline
    );

    if args.jobs > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.jobs)
            .build_global()?;
    }

    println!("{} Configuration:", text::bold("Comparison"));
    println!("  {:<20} {}", "Reference (A):", args.input_a.display());
    println!("  {:<20} {}", "Comparison (B):", args.input_b.display());
    println!("  {:<20} {}", "Output Directory:", args.output.display());
    println!("  {:<20} {}%", "Retention Target:", args.retention_threshold);
    println!("  {:<20} {}", "Coherence Target:", args.coherence_threshold);
    println!(
        "  {:<20} {}",
        "Parallel Jobs:",
        if args.jobs == 0 {
            "all available cores".to_string()
        } else {
            args.jobs.to_string()
        }
    );
    println!("{}\n", dline);

    let meta_path = find_first_metadata(&args.input_a)?;
    let (fx, fy, pixel_size, window_size) = load_metadata_info(&meta_path)?;
    let k_min_limit = 2.0 / (window_size as f64 * pixel_size);
    let coords = find_block_coords(&args.input_a)?;
    let n_blocks = coords.len();

    if !args.output.exists() {
        fs::create_dir_all(&args.output)?;
    }

    let global_acc = Mutex::new(GlobalAccumulator {
        block_profiles: Vec::with_capacity(n_blocks),
        cross_psd_total: Array2::zeros((fy.len(), fx.len())),
        psd_a_total: Array2::zeros((fy.len(), fx.len())),
        psd_b_total: Array2::zeros((fy.len(), fx.len())),
    });

    let progress_counter = Arc::new(Mutex::new(0));

    print!("Comparing blocks...");
    io::stdout().flush().unwrap();

    coords.par_iter().for_each(|coord| {
        let res: Result<()> = (|| {
            // Always load complex data for coherence
            let b_a = load_block(&args.input_a, coord, true)?;
            let b_b = load_block(&args.input_b, coord, true)?;
            let (_, r_a) = compute_radial_mean(&b_a.psd, &fx, &fy, 64, k_min_limit);
            let (waves, r_b) = compute_radial_mean(&b_b.psd, &fx, &fy, 64, k_min_limit);
            let ratio = &r_b / &r_a;

            let mut coh_2d = Array2::<f64>::zeros((fy.len(), fx.len()));
            let ca = b_a.complex.as_ref().unwrap();
            let cb = b_b.complex.as_ref().unwrap();
            for r in 0..fy.len() {
                for c in 0..fx.len() {
                    let den = b_a.psd[[r, c]] * b_b.psd[[r, c]];
                    let cross = ca[[r, c]] * cb[[r, c]].conj();
                    coh_2d[[r, c]] = if den > 1e-18 {
                        (cross.norm_sqr() / den).min(1.0)
                    } else {
                        0.0
                    };
                }
            }
            let (_, r_coh) = compute_radial_mean(&coh_2d, &fx, &fy, 64, k_min_limit);

            if args.save_partials {
                let mut wtr = csv::Writer::from_path(
                    args.output
                        .join(format!("comparison_block_{}_{}.csv", coord.row, coord.col)),
                )?;
                wtr.write_record(["wavelength", "psd_a", "psd_b", "ratio", "coherence"])?;
                for i in 0..waves.len() {
                    if r_a[i].is_finite() {
                        wtr.write_record([
                            waves[i].to_string(),
                            r_a[i].to_string(),
                            r_b[i].to_string(),
                            ratio[i].to_string(),
                            r_coh[i].to_string(),
                        ])?;
                    }
                }
            }

            let mut acc = global_acc.lock().unwrap();
            acc.psd_a_total += &b_a.psd;
            acc.psd_b_total += &b_b.psd;
            acc.block_profiles.push(BlockProfiles { ratio });
            let cross_total = &mut acc.cross_psd_total;
            for r in 0..fy.len() {
                for c in 0..fx.len() {
                    cross_total[[r, c]] += ca[[r, c]] * cb[[r, c]].conj();
                }
            }
            Ok(())
        })();
        if let Err(e) = res {
            eprintln!("\n{} Error at {:?}: {}", text::error("Error"), coord, e);
        }

        let mut count = progress_counter.lock().unwrap();
        *count += 1;
        let term = Term::stdout();
        let _ = term.clear_line();
        print!("\rComparing blocks... {:.0}%", (*count as f32 / n_blocks as f32) * 100.0);
        let _ = io::stdout().flush();
    });

    let term = Term::stdout();
    let _ = term.clear_line();
    println!("\r{} All blocks compared and aggregated.", text::check_icon());

    let acc = global_acc.into_inner().unwrap();
    let n_blocks_processed = acc.block_profiles.len();
    if n_blocks_processed == 0 {
        anyhow::bail!("No blocks processed.");
    }

    let global_psd_a = acc.psd_a_total / n_blocks_processed as f64;
    let global_psd_b = acc.psd_b_total / n_blocks_processed as f64;
    let (waves, radial_a_global) = compute_radial_mean(&global_psd_a, &fx, &fy, 64, k_min_limit);
    let (_, radial_b_global) = compute_radial_mean(&global_psd_b, &fx, &fy, 64, k_min_limit);

    let avg_cross = acc.cross_psd_total / n_blocks_processed as f64;
    let mut coh_2d_global = Array2::<f64>::zeros((fy.len(), fx.len()));
    for r in 0..fy.len() {
        for c in 0..fx.len() {
            let den = global_psd_a[[r, c]] * global_psd_b[[r, c]];
            coh_2d_global[[r, c]] = if den > 1e-18 {
                (avg_cross[[r, c]].norm_sqr() / den).min(1.0)
            } else {
                0.0
            };
        }
    }
    let (_, global_coh_profile) = compute_radial_mean(&coh_2d_global, &fx, &fy, 64, k_min_limit);

    let mut header = vec!["wavelength", "mean_psd_a", "mean_psd_b", "mean_ratio"];
    if n_blocks_processed > 1 {
        header.push("std_ratio");
    }
    header.push("coherence");

    let mut wtr = csv::Writer::from_path(args.output.join("comparison_summary.csv"))?;
    wtr.write_record(&header)?;

    let mut mean_ratio_for_retention = vec![];
    for i in 0..waves.len() {
        if !radial_a_global[i].is_finite() {
            continue;
        }
        let block_ratios: Vec<f64> = acc
            .block_profiles
            .iter()
            .map(|p| p.ratio[i])
            .filter(|v| v.is_finite())
            .collect();
        if block_ratios.is_empty() {
            continue;
        }

        let mean_r = block_ratios.iter().sum::<f64>() / block_ratios.len() as f64;
        mean_ratio_for_retention.push((waves[i], mean_r));

        let mut row = vec![
            waves[i].to_string(),
            radial_a_global[i].to_string(),
            radial_b_global[i].to_string(),
            mean_r.to_string(),
        ];
        if n_blocks_processed > 1 {
            let variance = block_ratios
                .iter()
                .map(|v| (v - mean_r).powi(2))
                .sum::<f64>()
                / block_ratios.len() as f64;
            row.push(variance.sqrt().to_string());
        }
        row.push(global_coh_profile[i].to_string());
        wtr.write_record(&row)?;
    }

    if let Some(plot_path) = args.plot {
        let (mut l, mut da, mut db) = (vec![], vec![], vec![]);
        for i in 0..waves.len() {
            if radial_a_global[i].is_finite()
                && radial_a_global[i] > 0.0
                && radial_b_global[i] > 0.0
            {
                l.push(format!("{:.2}", waves[i]));
                da.push(radial_a_global[i]);
                db.push(radial_b_global[i]);
            }
        }
        l.reverse();
        da.reverse();
        db.reverse();
        let chart = Chart::new()
            .title(Title::new().text("Spectral Comparison").left("center"))
            .legend(Legend::new().top("bottom"))
            .x_axis(
                ChartAxis::new()
                    .type_(AxisType::Category)
                    .name("Wavelength (m)")
                    .data(l),
            )
            .y_axis(ChartAxis::new().type_(AxisType::Log).name("Power"))
            .series(Line::new().name("Mean PSD A").data(da))
            .series(Line::new().name("Mean PSD B").data(db))
            .tooltip(
                Tooltip::new()
                    .trigger(Trigger::Axis)
                    .trigger_on(TriggerOn::Mousemove),
            );
        HtmlRenderer::new("Spectral Comparison", 1024, 768).save(&chart, &plot_path)?;
        println!(
            "{} Comparison plot saved to {}",
            text::check_icon(),
            plot_path.display()
        );
    }

    let mut thresholds_found = false;

    if let Some(w) = mean_ratio_for_retention
        .iter()
        .find(|s| s.1.is_finite() && s.1 < args.retention_threshold / 100.0)
        .map(|s| s.0)
    {
        if !thresholds_found {
            println!("\n{} Spectral Thresholds:", text::bold("Summary"));
            println!("{}", line);
            thresholds_found = true;
        }
        println!(
            "  {:<30} {}",
            format!("{}% Power Retention:", args.retention_threshold),
            text::highlight(format!("{:.2} m", w))
        );
    }

    if let Some(w) = waves
        .iter()
        .zip(global_coh_profile.iter())
        .find(|(_, &c)| c.is_finite() && c < args.coherence_threshold)
        .map(|(w, _)| *w)
    {
        if !thresholds_found {
            println!("\n{} Spectral Thresholds:", text::bold("Summary"));
            println!("{}", line);
        }
        println!(
            "  {:<30} {}",
            format!("Coherence ({:.2}):", args.coherence_threshold),
            text::highlight(format!("{:.2} m", w))
        );
    }

    println!("{}", line);
    println!("{}", text::success("Comparison completed successfully."));
    println!("");

    Ok(())
}
