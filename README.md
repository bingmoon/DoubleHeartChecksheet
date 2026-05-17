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

    # Core Data & Survival Analysis
    install.packages(c("nhanesA", "dplyr", "tidyr", "readr", "survival", "survminer"))

    # Mendelian Randomization (requires devtools)
    devtools::install_github("MRCIEU/TwoSampleMR")
    install.packages("ieugwasr")

    # Clinical Modeling & Visualization
    install.packages(c("rms", "dcurves", "timeROC", "ggplot2", "grid", "gridExtra"))

*Note: A valid API token from the [MRC IEU OpenGWAS project](https://gwas.mrcieu.ac.uk/) may be required for the `TwoSampleMR` extraction steps.*

## 🚀 Usage Guide (Pipeline)

### Part 1: NHANES Observational Analysis (`Part1_NHANES_Analysis.R`)
* Connects to the NHANES API to fetch demographics, MCQ, DPQ, and BMX data.
* Parses CDC NDI fixed-width mortality data.
* Executes univariable and multivariable Cox proportional hazards regressions.
* Generates Lollipop charts, Forest plots, and Kaplan-Meier survival curves.

### Part 2: Mendelian Randomization (`Part2_MR_and_MVMR.R`)
* Executes high-throughput bidirectional MR screening.
* Performs sensitivity validations (Leave-one-out, Funnel plots, Pleiotropy tests).
* Runs Multivariable MR (MVMR) to test the metabolic mediation effect of BMI.

### Part 3: Clinical Modeling & Tool Generation (`Part3_Clinical_Modeling.R`)
* Builds the integrated prediction model using `rms`.
* Plots the predictive nomogram, calibration curves, and Decision Curve Analysis (DCA).
* Evaluates the simplified Double Heart (DH) checksheet using time-dependent ROC.
* Auto-generates the clinical supplementary forms in HTML format (ready for PDF conversion).

## ⚠️ Data Availability Statement
* The **NHANES** questionnaire and examination data are publicly available at the [CDC NHANES website](https://wwwn.cdc.gov/nchs/nhanes/). 
* The **NDI linked mortality files** (`.dat`) must be downloaded independently according to the NCHS data use agreement and placed in the `/Data/` directory.
* **GWAS summary statistics** are automatically fetched via the `ieugwasr` API during script execution.

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact
For any programmatic issues, data requests, or methodological inquiries, please contact the corresponding author via the email provided in the manuscript.
