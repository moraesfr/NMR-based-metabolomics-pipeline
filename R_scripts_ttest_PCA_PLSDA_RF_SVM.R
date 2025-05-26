# Load required packages
library(MLmetrics)
library(tidyverse)
library(caret)
library(pROC)
library(e1071)
library(randomForest)
library(ggplot2)
library(patchwork)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("ropls")
library(ropls)

# Set the working directory
setwd("YOUR_DIRECTORY_HERE")

# Read the data
# Data should be organized as the first column to be the metabolite names,
# the second column and so on the metabolite concentration data,
# the last column should be the class label: control, disease, for example.

data <- read.csv("YOUR_DATA.csv", header = TRUE, row.names = 1)
group <- as.factor(data[, ncol(data)])
features <- data[, -ncol(data)]

# Univariate Welch Two-Sample t-tests
pvals <- apply(features, 2, function(x) t.test(x ~ group)$p.value)
signif_vars <- names(pvals)[pvals < 0.05]
write.csv(data.frame(Variable=names(pvals), P.Value=pvals), "welch_test_results.csv")

# Boxplots
# Create a folder to store the plots (if not exists)
dir.create("boxplots", showWarnings = FALSE)

# Loop over variables
for (var in colnames(features)) {
  df_plot <- data.frame(Value = features[[var]], Group = group)

  p <- ggplot(df_plot, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot() +
    theme_classic() +  # white background
    ggtitle(paste("Boxplot of", var)) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(filename = paste0("boxplots/", var, "_boxplot.png"), plot = p, width = 6, height = 4)
}

# Pareto scaling function
pareto_scale <- function(df) {
  scaled <- sweep(df, 2, colMeans(df), "-")
  scaled <- sweep(scaled, 2, sqrt(apply(df, 2, sd)), "/")
  return(scaled)
}

# Apply scaling
features <- pareto_scale(features)

# Perform PCA
pca_res <- prcomp(features, scale. = TRUE)

# PCA scores and loadings
scores <- as.data.frame(pca_res$x)
scores$Group <- group  # Add group info for coloring
loadings <- pca_res$rotation[, 1:5]

# Calculate variance explained
var_explained <- (pca_res$sdev[1:5])^2 / sum(pca_res$sdev^2)

# Compute WSSL
WSSL <- rowSums((loadings^2) %*% diag(var_explained[1:5]))
wssl_df <- data.frame(Variable = names(WSSL), WSSL = WSSL)
wssl_df <- wssl_df[order(-wssl_df$WSSL), ]

# Loop over PC pairs (PC1-PC2, PC1-PC3, ..., PC2-PC3, etc.)
pc_indices <- combn(1:5, 2)

for (i in 1:ncol(pc_indices)) {
  pc_x <- pc_indices[1, i]
  pc_y <- pc_indices[2, i]

  # Score plot
  p1 <- ggplot(scores, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), color = "Group")) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    theme_classic() +
    ggtitle(paste0("PCA Score Plot (PC", pc_x, " vs PC", pc_y, ")")) +
    xlab(paste0("PC", pc_x, " (", round(var_explained[pc_x]*100, 1), "%)")) +
    ylab(paste0("PC", pc_y, " (", round(var_explained[pc_y]*100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5))

  # WSSL plot (top 10 variables)
  top10 <- wssl_df[1:10, ]
  p2 <- ggplot(top10, aes(x = reorder(Variable, WSSL), y = WSSL)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_classic() +
    ggtitle("Top 10 Variables by WSSL") +
    xlab("Variable") + ylab("WSSL") +
    theme(plot.title = element_text(hjust = 0.5))

  # Combine and save
  combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  filename <- paste0("PCA_pairs/PCA_PC", pc_x, "_PC", pc_y, "_WSSL.png")
  ggsave(filename, combined_plot, width = 10, height = 5, dpi = 300)
}

# PLS-DA
# Run PLS-DA with 5 components
plsda <- opls(features, group, predI = 5, orthoI = 0)

# Get score matrix and VIP
scores <- as.data.frame(plsda@scoreMN)
scores$Group <- group
vip <- plsda@vipVn
vip_df <- data.frame(Variable = names(vip), VIP = vip)

# Create output folder
dir.create("PLSDA_pairs", showWarnings = FALSE)

# Generate all pairs of the first 5 PLS components
combs <- combn(1:5, 2)

# Loop over combinations
for (i in 1:ncol(combs)) {
  comp_x <- combs[1, i]
  comp_y <- combs[2, i]
  x_col <- paste0("p", comp_x)
  y_col <- paste0("p", comp_y)

  # Score plot
  p1 <- ggplot(scores, aes_string(x = x_col, y = y_col, color = "Group")) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    theme_classic() +
    ggtitle(paste("PLS-DA Score Plot:", x_col, "vs", y_col)) +
    xlab(x_col) + ylab(y_col) +
    theme(plot.title = element_text(hjust = 0.5))

  # Top 10 VIPs
  top10 <- vip_df[order(-vip_df$VIP), ][1:10, ]

  p2 <- ggplot(top10, aes(x = reorder(Variable, VIP), y = VIP)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_classic() +
    ggtitle("Top 10 VIP Scores") +
    xlab("Variable") + ylab("VIP") +
    theme(plot.title = element_text(hjust = 0.5))

  # Combine and save
  combined_plot <- p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))
  fname <- paste0("PLSDA_pairs/", x_col, "_vs_", y_col, ".png")
  ggsave(fname, plot = combined_plot, width = 10, height = 5, dpi = 300)
}


# --- Adjusted Random Forest and SVM with Leave-One-Out Cross-Validation ---
# Custom MCC function
mcc_custom <- function(actual, predicted) {
  cm <- table(actual, predicted)
  if (nrow(cm) < 2 || ncol(cm) < 2) return(NA)
  TP <- cm[2, 2]
  TN <- cm[1, 1]
  FP <- cm[1, 2]
  FN <- cm[2, 1]
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (denominator == 0) return(0)
  return(numerator / denominator)
}

# Leave-One-Out CV for Random Forest
rf_probs <- numeric(nrow(features))
rf_preds <- character(nrow(features))
for (i in 1:nrow(features)) {
  train_idx <- setdiff(1:nrow(features), i)
  model <- randomForest(x = features[train_idx, ], y = group[train_idx], importance = TRUE)
  rf_probs[i] <- predict(model, features[i, , drop = FALSE], type = "prob")[, 2]
  rf_preds[i] <- as.character(predict(model, features[i, , drop = FALSE]))
}

actual <- factor(group)
predicted <- factor(rf_preds, levels = levels(actual))
auc_rf <- auc(actual, rf_probs)
precision_rf <- sum(predicted == "1" & actual == "1") / sum(predicted == "1")
recall_rf <- sum(predicted == "1" & actual == "1") / sum(actual == "1")
mcc_rf <- mcc_custom(actual, predicted)

metrics_df <- data.frame(
  Metric = c("AUC", "MCC", "Precision", "Recall"),
  Value = c(auc_rf, mcc_rf, precision_rf, recall_rf)
)

p1 <- ggplot(metrics_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  ylim(0, 1) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.5) +
  theme_classic() +
  ggtitle("Random Forest - LOOCV Metrics")

# Fit full model for importance
rf_model <- randomForest(x = features, y = group, importance = TRUE)
imp_df <- data.frame(
  Variable = rownames(importance(rf_model)),
  Importance = importance(rf_model)[, "MeanDecreaseGini"]
)
top10_imp <- imp_df[order(-imp_df$Importance), ][1:10, ]

p2 <- ggplot(top10_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  ggtitle("Top 10 Important Variables") +
  xlab("Variable") + ylab("Mean Decrease Gini")

ggsave("RandomForest_LOOCV_Metrics_and_Importance.png", p1 + p2, width = 10, height = 5, dpi = 300)

# --- Support Vector Machine with LOOCV and Variable Importance ---
# Initialize vectors to store predictions and probabilities
svm_probs <- numeric(nrow(features))
svm_preds <- character(nrow(features))

# LOOCV loop
for (i in 1:nrow(features)) {
  train_idx <- setdiff(1:nrow(features), i)

  # Train SVM with probability enabled
  model <- svm(x = features[train_idx, ], y = group[train_idx],
               probability = TRUE, kernel = "radial")

  # Predict the left-out sample
  pred <- predict(model, features[i, , drop = FALSE], probability = TRUE)
  prob <- attr(pred, "probabilities")

  # Store predicted probability for the positive class
  svm_probs[i] <- prob[1, levels(group)[2]]
  svm_preds[i] <- as.character(pred)
}

# Convert predictions to factor
actual <- factor(group)
predicted <- factor(svm_preds, levels = levels(actual))

# Compute metrics
tp <- sum(predicted == levels(actual)[2] & actual == levels(actual)[2])
fp <- sum(predicted == levels(actual)[2] & actual != levels(actual)[2])
fn <- sum(predicted != levels(actual)[2] & actual == levels(actual)[2])

precision_svm <- if ((tp + fp) > 0) tp / (tp + fp) else NA
recall_svm <- if ((tp + fn) > 0) tp / (tp + fn) else NA
auc_svm <- auc(actual, svm_probs)
mcc_svm <- mcc_custom(actual, predicted)

# Create metrics data frame
svm_metrics_df <- data.frame(
  Metric = c("AUC", "MCC", "Precision", "Recall"),
  Value = c(auc_svm, mcc_svm, precision_svm, recall_svm)
)

# Metrics bar plot
p3 <- ggplot(svm_metrics_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  ylim(0, 1) +
  geom_text(aes(label = round(Value, 3)), vjust = -0.5) +
  theme_classic() +
  ggtitle("SVM - LOOCV Metrics")

# --- Variable Importance using caret::train on full dataset ---
# Wrap features and group in a data frame for caret
data_for_caret <- data.frame(features, Class = group)

# Train SVM model with caret
svm_caret_model <- train(Class ~ ., data = data_for_caret,
                         method = "svmRadial",
                         trControl = trainControl(method = "none"),
                         preProcess = c("center", "scale"))

# Compute variable importance
svm_imp <- varImp(svm_caret_model)$importance

# If importance has multiple columns (e.g., one per class), compute row means
if (ncol(svm_imp) > 1) {
  svm_imp$Overall <- rowMeans(svm_imp)
} else {
  svm_imp$Overall <- svm_imp[,1]
}

# Add variable names
svm_imp$Variable <- rownames(svm_imp)

# Sort and select top 10
top10_svm <- svm_imp[order(-svm_imp$Overall), ][1:10, ]

# Importance plot
p4 <- ggplot(top10_svm, aes(x = reorder(Variable, Overall), y = Overall)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_classic() +
  xlab("Variable") +
  ylab("Importance") +
  ggtitle("Top 10 Important Variables (SVM)")

# Combine plots side by side
combined_plot <- p3 + p4

# Save plot
ggsave("SVM_LOOCV_Metrics_and_Importance.png", combined_plot, width = 10, height = 5, dpi = 300)
rm(list=ls())
