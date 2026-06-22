# DoubleHeartChecksheet

## Overview
Complete analytical pipeline for the study *"Dissecting the Psychometabolic-Cardiac Axis"*.

This repository contains all R scripts required to reproduce the analyses, from data download to final figures.

## How to Run
1. Set working directory to the repository root.
2. Run `00_setup.R` to install required packages.
3. Run `E01_import_data.R` to download NHANES 2009–2018 data.
4. Run `E02_analysis.R` for survival analyses.
5. Run `03_mr_analysis.R` for Mendelian Randomization (requires OpenGWAS token).
6. Run `E04_clinical_tools.R` for the PMRS, nomogram, DCA, and NRI.
7. Run `E05_figures.R` to generate all figures.
8. Run `E06_results.R` to extract key numerical results.

## Data Sources
All data are downloaded automatically by the scripts:
- **NHANES**: Accessed via the `nhanesA` package from the CDC public API.
- **GWAS summary statistics**: Accessed via `TwoSampleMR` from the IEU OpenGWAS platform.

No pre-downloaded data files are required.

## Requirements
- R ≥ 4.2.2
- R packages: nhanesA, survival, rms, timeROC, TwoSampleMR, ieugwasr, nricens, survminer, ggplot2, dplyr, broom, readr, tidyr

## License
MIT
