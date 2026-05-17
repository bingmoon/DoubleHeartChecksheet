# Double Heart (DH) Clinical Screening Checksheet

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%3E%3D%204.2.2-brightgreen)]()

## 📖 Overview
This repository contains the complete analytical code for the manuscript:
**"Elucidating the Psychometabolic-Cardiac Axis via Multivariable Mendelian Randomization: Development and Validation of a Bedside Nursing Checksheet from the NHANES Cohort"**

The project provides a comprehensive tri-pillared programmatic pipeline:
1. **Real-world Evidence:** Mining and survival analysis of the NHANES (2015-2018) cardiovascular cohort.
2. **Genetic Causal Inference:** Bidirectional and Multivariable Mendelian Randomization (MVMR) to elucidate the psychometabolic-cardiac axis.
3. **Clinical Translation:** Development of a high-precision nomogram and a simplified 17-point bedside checksheet with time-dependent ROC and DCA validations.

## 📂 Repository Structure
* `/Scripts/`: Contains the core R scripts divided into three sequential parts.
* `/Data/`: Directory intended for raw data placement (e.g., CDC NDI linked mortality `.dat` files).
* `/Output/`: Default directory where all generated plots (PDFs), tables (CSVs), and clinical HTML forms are exported.

## 💻 Prerequisites & Environment
All analyses were conducted in the **R programming environment (v4.2.2)**. Ensure the following R packages are installed before running the scripts:

```R
# Core Data & Survival Analysis
install.packages(c("nhanesA", "dplyr", "tidyr", "readr", "survival", "survminer"))

# Mendelian Randomization (requires devtools)
devtools::install_github("MRCIEU/TwoSampleMR")
install.packages("ieugwasr")

# Clinical Modeling & Visualization
install.packages(c("rms", "dcurves", "timeROC", "ggplot2", "grid", "gridExtra"))
