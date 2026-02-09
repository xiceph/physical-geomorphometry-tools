use clap::{Parser, ValueEnum};
use std::path::PathBuf;

/// Command-line arguments for the fft-analyze tool.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = "Analyzes polar FFT spectrums.")]
pub struct Args {
    /// Directory containing the fft_polar_block_*.tif files.
    #[arg(long)]
    pub input: PathBuf,

    /// Path to save the output CSV file.
    #[arg(long)]
    pub output: PathBuf,

    /// Optional path to save the output plot file.
    #[arg(long)]
    pub plot: Option<PathBuf>,

    /// The analysis mode to perform.
    #[arg(long, value_enum)]
    pub mode: AnalysisMode,

    /// Optional wavelength bounds (min,max) for filtering.
    #[arg(long)]
    pub wavelength_bounds: Option<String>,

    /// Optional angle bounds (min,max) for filtering.
    #[arg(long)]
    pub angle_bounds: Option<String>,

    /// Save partial results for each block.
    #[arg(long)]
    pub save_partials: bool,

    /// Polynomial order for detrending in radial-mean mode.
    #[arg(long)]
    pub detrend: Option<usize>,

    /// Number of parallel jobs to run. Defaults to 0 (Rayon chooses).
    #[arg(long, default_value_t = 0)]
    pub jobs: usize,
}

/// Defines the available analysis modes.
#[derive(Debug, Clone, ValueEnum)]
pub enum AnalysisMode {
    /// Calculate mean power along the angle axis to get a 1D radial profile.
    RadialMean,
    /// Calculate mean power along the wavelength axis to get a 1D angular profile.
    AngularMean,
}
