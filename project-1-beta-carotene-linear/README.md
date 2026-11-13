# Plasma Beta-Carotene: Regression Analysis in R

University group project (2025) | Linear and multiple regression | R

## Overview

This project investigates associations between plasma beta-carotene concentrations and demographic, lifestyle, and dietary variables in a cross-sectional dataset of 315 participants. The analysis compares models on the original and log-transformed response scales, assesses regression assumptions, and studies the effect of influential observations on model selection.

## Methods

- Simple linear regression and log-response modelling
- Confidence and prediction intervals; coefficient interpretation and hypothesis tests
- Categorical predictors and multiple linear regression
- Correlation analysis and variance inflation factors (VIF/GVIF)
- Studentized residuals, leverage, Cook's distance, and DFBETAS
- Model comparison and BIC-based stepwise selection

## Main findings

The log-response models were more suitable for the observed skewness and regression diagnostics. BMI, dietary fiber, calories, and vitamin-use categories were included in the final reduced-data model reported in the coursework. After excluding one influential observation and combining the vitamin-use categories into Yes/No, the selected model had approximately **R² = 0.200** and **adjusted R² = 0.190**. These are associations in an observational study, not causal effects.

**Interpretation note:** Exponentiating fitted values from a regression for `log(betaplasma)` gives estimates on the original scale corresponding to conditional medians under a log-normal model, not automatically conditional arithmetic means.

## Files

| Path | Description |
| --- | --- |
| `analysis.R` | R analysis and visualizations |
| `report.pdf` | Group project report |
| `data/carotene.xlsx` | Course-provided, adapted dataset **only if redistribution is authorized** |

The original assignment brief can be summarized here rather than copied into the repository.

## Running the analysis

1. Install R and the required packages: `readxl`, `tidyverse`, `car`, `rstatix`, `GGally`, `patchwork`, and `plyr` (the code uses `plyr::count()`).
2. Make sure you have an authorized copy of the adapted dataset at `data/carotene.xlsx`.
3. Open R or RStudio with the project folder as the working directory and run `source("analysis.R")`.

The code expects the **course-adapted** spreadsheet: its column names and measurement units may differ from those in other public versions of the study data. Before publishing, ensure `analysis.R` reads `data/carotene.xlsx` using a relative path. Re-run the script and check the output after any cleanup or model-selection changes.

## Data source

Nierenberg, D. W., Stukel, T. A., Baron, J. A., Dain, B. J., & Greenberg, E. R. (1989). *Determinants of plasma levels of beta-carotene and retinol*. American Journal of Epidemiology, 130, 511–521.

Related documentation: [R `CAMAN::betaplasma` dataset](https://search.r-project.org/CRAN/refmans/CAMAN/html/betaplasma.html). The associated original dataset's redistribution permissions should be checked separately before including a copy of `carotene.xlsx`.

## Contributions

This was a **group project** with Kseniia Sharnina and Konstantinos Vitos. My contributions included technical writing and revision of the introduction and analysis sections, interpretation of results, and preparation/revision of plots and tables based on the group's code and analyses. Programming and the statistical analysis were collaborative; the original report documents the division of work in more detail.

## Limitations

The study is observational, the explanatory power of the selected model is moderate, and stepwise selection does not guarantee a globally optimal predictor subset. Results should be interpreted with these limitations in mind.
