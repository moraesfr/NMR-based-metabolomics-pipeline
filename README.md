# Metabolomics Machine Learning Analysis Pipeline

This repository contains an R script for advanced statistical and machine learning analysis of metabolomics data. The pipeline covers univariate testing, multivariate analysis (PCA, PLS-DA), model building (Random Forest, SVM), cross-validation, and visualization.

## Overview

- **Input:** A CSV file with metabolite measurements (samples as rows, metabolites as columns, last column as group/class).
- **Processing:**
  - Welch's t-test for univariate analysis.
  - Pareto scaling for normalization.
  - Principal Component Analysis (PCA) and loadings/importance (WSSL).
    -  Weighted Sum of Squared Loadings (WSSL) are squared loading by the proportion of variance explained by the corresponding PC. This accounts for the fact that PC1 usually explains more variance than PC5, giving more weight to early components.
  - Partial Least Squares Discriminant Analysis (PLS-DA) and VIP scores.
  - Leave-One-Out Cross-Validation (LOOCV) for Random Forest and SVM classifiers.
  - Calculation of model metrics (AUC, MCC, Precision, Recall).
  - Variable importance ranking for both RF and SVM.
  - Publication-ready visualization of results.
- **Output:** PNG image files for plots, and CSV files with results and variable rankings.

## Requirements

- R (version ≥ 4.0 recommended)
- R packages:
  - `MLmetrics`
  - `tidyverse`
  - `caret`
  - `pROC`
  - `e1071`
  - `randomForest`
  - `ggplot2`
  - `patchwork`
  - `ropls` (Bioconductor)

Install required packages in R with:
```r
install.packages(c("MLmetrics", "tidyverse", "caret", "pROC", "e1071", "randomForest", "ggplot2", "patchwork"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ropls")
```

## Usage

1. **Set Up Working Directory**

   Edit the script's `setwd()` line to your data folder:
   ```r
   setwd("YOUR/DATA/PATH")
   ```

2. **Prepare Input File**

   Ensure your CSV file (e.g., `G1_G3_12_DSS.csv`) is in the working directory. The file should have:
   - Rows: samples
   - Columns: metabolites (last column = group/class label)
   - Header row with column names

3. **Run the Script**

   Execute the script in R or RStudio.

4. **Outputs**

   - `welch_test_results.csv`: Welch t-test results for all variables.
   - `boxplots/`: Folder with boxplots of each metabolite.
   - `PCA_pairs/`: PCA score and WSSL plots for PC pairs.
   - `PLSDA_pairs/`: PLS-DA score plots and VIP importance plots for component pairs.
   - `RandomForest_LOOCV_Metrics_and_Importance.png`: RF LOOCV performance metrics and top variables.
   - `SVM_LOOCV_Metrics_and_Importance.png`: SVM LOOCV performance metrics and top variables.

## Script Structure

- **Load Packages:** Ensures required libraries are loaded.
- **Read Data:** Loads CSV and separates features and group information.
- **Univariate Analysis:** Welch's t-test for group differences.
- **Boxplots:** Creates and saves boxplots for each metabolite.
- **Pareto Scaling:** Normalizes data.
- **PCA:** Dimensionality reduction, variance explained, and variable importance (WSSL).
- **PLS-DA:** Supervised multivariate analysis with VIP scores.
- **Random Forest & SVM with LOOCV:** Model evaluation (AUC, MCC, Precision, Recall) and variable importance.
- **Visualization:** Outputs are saved as CSV and publication-ready PNG images.

## Notes

- The script uses leave-one-out cross-validation for robust model evaluation.
- If your group variable is not binary, you may need to adapt the script for multiclass analysis.
- Ensure that the group labels are factors and that the positive class is the second level for metric calculation.

## Example

```r
# Set working directory (edit this line)
setwd("C:/Users/YourUsername/your_data_folder")

# Run the entire script in R or RStudio.
```

## License

MIT License (or specify your own).

## Contact

For questions or suggestions, open an issue or contact moraesfr@gmail.com.
