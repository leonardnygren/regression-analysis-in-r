# Logistic Regression Analysis of Plasma β-Carotene Levels

**University coursework | R | May 2025**

## Overview

This project investigates the association between plasma β-carotene concentrations and demographic, lifestyle, and dietary variables using logistic regression.

The analysis uses a dataset containing 315 observations. Plasma β-carotene concentrations are classified into low and high categories using a threshold of 225.5 ng/mL.

## Methods

The project explores binary logistic regression, odds ratios and confidence intervals, multivariable modelling with interaction terms, and model selection using AIC and BIC.

Model diagnostics include leverage, standardized deviance residuals, and Cook's distance. Classification performance is evaluated using confusion matrices, sensitivity, specificity, ROC curves, and AUC.

## Key Findings

Vitamin use and BMI showed statistically significant associations with low plasma β-carotene concentrations in the initial models.

The multivariable analysis compared several candidate models. A reduced model combining regular and irregular vitamin users into one category achieved an in-sample AUC of approximately 0.788.

The results illustrate the trade-offs between model complexity, statistical fit, and classification performance.

## Tools

**Language:** R

**Libraries:** tidyverse, readxl, broom, car, caret, pROC

## Limitations

The classification performance was evaluated using the same data employed for fitting and model selection. Therefore, reported classification metrics should not be interpreted as independently validated predictive performance.

The dataset is observational, and associations should not be interpreted as causal effects.

## Report and Contributions

This project was completed jointly by Leonard Nygren Löhndorf and Christian Bjur as part of a university course in linear and logistic regression.

Leonard was responsible for report Sections 1 and 3, while Christian was responsible for Sections 2 and 4. Both authors proofread the report.

Some code was adapted from course lectures, and AI assistance was used for tables and parts of the plotting code, as disclosed in the original report.

The complete methodology and results are provided in the accompanying PDF report.