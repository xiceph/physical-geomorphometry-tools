# FFT Compare

## Overview

`fft-compare` is a tool designed to analyze the spectral differences between two datasets processed by `fft-process`. It is particularly useful for quantifying information loss or signal attenuation by comparing a "reference" dataset with a "modified" (e.g., generalized, interpolated, or filtered) dataset.

## Features

- **PSD Ratio Analysis:** Compute the ratio of Power Spectral Density between two datasets to identify scales where energy is lost or gained.
- **Spectral Coherence:** Calculate the magnitude squared coherence to measure the linear correlation and phase consistency between datasets.
- **Retention Analysis:** Determine the wavelength threshold at which a certain percentage (e.g., 50%) of power is retained.
- **Resolution Assessment:** Automatically identify the wavelength where coherence drops below a threshold (default 0.5), defining the "effective resolution".
- **Visualization:** Generate interactive spectral plots (HTML) showing Mean PSD, Ratios, and Coherence.

## Usage

### Command Line Arguments

```bash
fft-compare [OPTIONS] --input-a <DIR> --input-b <DIR> --output <DIR>
```

#### Required Arguments
- `--input-a <DIR>`: Directory containing the reference `fft-process` results.
- `--input-b <DIR>`: Directory containing the comparison `fft-process` results.
- `--output <DIR>`: Directory where CSV summaries will be saved.

#### Optional Arguments
- `--retention-threshold <PERCENT>`: Power percentage (default `50.0`) to find the retention wavelength.
- `--coherence-threshold <VALUE>`: Coherence value (default `0.5`) to find the effective resolution wavelength.
- `--plot <FILE>`: Path to save an interactive HTML plot (e.g., `comparison.html`).
- `--save-partials`: Save individual comparison CSVs for every block. Default is `false`.
- `--jobs <NUM>`: Number of parallel jobs to run. Defaults to 0 (all available cores).

### Example

Compare a reference dataset with a generalized one:

```bash
fft-compare --input-a ./results/original --input-b ./results/modified \
    --output ./results/comparison --plot ./results/comparison/spectral_comparison.html \
    --retention-threshold 50.0 --coherence-threshold 0.5
```

## Output

1.  **`comparison_summary.csv`**: A global summary containing wavelengths, mean PSDs, mean ratios, and coherence values.
2.  **`comparison_block_ROW_COL.csv`** (Optional): Partial results for individual blocks if `--save-partials` is used.
3.  **Spectral Comparison Plot**: An interactive HTML file visualizing the spectral differences.

## Scientific Relevance

In topographic analysis, quantifying the difference between two DEMs is often done using simple RMSE (Root Mean Square Error). However, RMSE is a spatial-domain metric that provides no information about *which scales* are being affected. 

`fft-compare` moves this analysis into the frequency domain, allowing for a much more nuanced assessment:

### 1. Information Retention (PSD Ratio)
The ratio of the Power Spectral Density ($PSD_B / PSD_A$) tells us how much topographic variance is preserved at each wavelength. A ratio of 1.0 means perfect power preservation. A drop below 0.5 (the "50% retention" mark) is often used to define the scale at which a generalization or smoothing process has removed half of the original topographic signal.

### 2. Spectral Coherence
While the PSD ratio measures power, **Spectral Coherence** ($C_{AB}$) measures the linear relationship between the two signals. It ranges from 0 to 1.
*   **$C_{AB} \approx 1$:** The signals are perfectly correlated (same features, same phases).
*   **$C_{AB} \approx 0$:** The signals are uncorrelated.
A coherence threshold of **0.5** is a physically significant benchmark. If coherence drops below 0.5, the features at that scale in the modified dataset no longer reliably represent the reference features, even if the total power (PSD) remains high.

## Analyzing the Output

When comparing a reference DEM with a modified one (e.g., after filtering or interpolation):

| Scenario | Interpretation |
| :--- | :--- |
| **High PSD Ratio & High Coherence** | Data is well-preserved. The modification has not significantly altered features at this scale. |
| **Low PSD Ratio & Low Coherence** | Signal Loss. The modification has smoothed out or removed features at this scale (e.g., low-pass filtering). |
| **High PSD Ratio & Low Coherence** | **Signal Replacement / Noise.** This is a critical edge case. It indicates that the modified dataset has power at this scale, but the features don't match the reference. This often happens due to **aliasing**, **interpolation artifacts**, or the introduction of synthetic noise. |
| **Coherence Threshold Wavelength** | This defines the **Effective Resolution** of the dataset. Wavelengths shorter than this should not be used for quantitative geomorphic analysis (e.g., slope or curvature calculation). |