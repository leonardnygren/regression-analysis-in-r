# Regression Analysis in R

Two university group projects (2025) applying linear and logistic regression to a dataset of 315 participants. Both investigate associations between plasma beta-carotene concentrations and demographic, lifestyle, and dietary variables, with different modelling objectives.

## Projects

### [Project 1 — Linear Regression](project-1-beta-carotene-linear/)

Models plasma beta-carotene concentration as a continuous response using simple and multiple linear regression. The analysis compares untransformed and log-transformed responses, examines regression assumptions and influential observations, and explores model selection.

**Methods:** Confidence and prediction intervals, hypothesis testing, categorical predictors, VIF/GVIF, studentized residuals, leverage, Cook's distance, DFBETAS, and BIC-based model selection.

[Code, report, findings, and contributions](project-1-beta-carotene-linear/README.md)

### [Project 2 — Logistic Regression](project-2-beta-carotene-logistic/)

Models the probability of low plasma beta-carotene concentration (below 225.5 ng/mL) as a binary response. The analysis compares logistic regression models, investigates interactions and variable selection, and evaluates classification performance.

**Methods:** Odds ratios and confidence intervals, likelihood-ratio and Wald tests, AIC/BIC, deviance residuals, leverage, Cook's distance, confusion matrices, sensitivity and specificity, ROC curves, AUC, and classification thresholds.

[Code, report, findings, and contributions](project-2-beta-carotene-logistic/README.md)

## Tools

**R**, including `tidyverse`, `ggplot2`, `readxl`, `car`, `broom`, `caret`, and `pROC` (packages vary by project).

## Notes

- Both projects use the same course-adapted dataset and were completed collaboratively, with different co-authors. Authorship and individual contributions are documented in their respective project reports and READMEs.
- Each project folder contains its analysis and report; consult its README for data requirements and instructions for reproducing the analysis.
- The data are observational: reported associations do not establish causation.
- Classification results in Project 2 are evaluated in-sample and should not be interpreted as independently validated predictive performance.
- Any reuse or redistribution restrictions for the course-provided data should be checked before publishing a copy.
