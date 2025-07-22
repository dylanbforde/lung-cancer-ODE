# =============================================================================
# TCGA-LUAD Survival Analysis Pipeline
# Based on your original code with improvements
# =============================================================================

# Load required libraries
library(TCGAbiolinks)
library(survival)
library(survminer)
library(dplyr)
library(glmnet)


# =============================================================================
# 1. DATA LOADING
# =============================================================================

# Set your file path here - UPDATED PATH
expr_path <- "C:/Users/mahen/Documents/TCGA.LUAD.sampleMap_HiSeqV2/HiSeqV2"

# Check if file exists and is readable
if (!file.exists(expr_path)) {
  stop("File not found at: ", expr_path)
}

cat("File found at:", expr_path, "\n")

# Try to read first few lines to verify file format
tryCatch({
  test_read <- read.delim(expr_path, header = TRUE, nrows = 2)
  cat("File format verified - Columns:", ncol(test_read), "Rows in preview:", nrow(test_read), "\n")
}, error = function(e) {
  cat("Error reading file:", e$message, "\n")
  stop("Please check if the file is a valid tab-delimited format")
})

# Download clinical data from GDC
cat("Downloading clinical data...\n")
luad_clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

# Read UCSC Xena FPKM expression data
cat("Loading expression data...\n")
expr_data <- read.delim(expr_path, header = TRUE, row.names = 1, check.names = FALSE)

# Keep only tumor samples ending with "-01"
tumor_samples <- grep("-01$", colnames(expr_data), value = TRUE)
expr_data <- expr_data[, tumor_samples]

# Shorten sample IDs to match clinical barcodes (first 12 characters)
colnames(expr_data) <- substr(colnames(expr_data), 1, 12)

cat("Expression data loaded:", nrow(expr_data), "genes,", ncol(expr_data), "samples\n")

# =============================================================================
# 2. MATCH CLINICAL AND EXPRESSION DATA
# =============================================================================

# Match clinical data
sample_ids <- colnames(expr_data)
clinical_data <- luad_clinical[luad_clinical$bcr_patient_barcode %in% sample_ids, ]

# Ensure same order
expr_data <- expr_data[, sample_ids %in% clinical_data$bcr_patient_barcode]
clinical_data <- clinical_data[match(colnames(expr_data), clinical_data$bcr_patient_barcode), ]

cat("Matched samples:", ncol(expr_data), "\n")

# =============================================================================
# 3. BUILD SURVIVAL OBJECT
# =============================================================================

# Create survival time
clinical_data$time <- ifelse(is.na(clinical_data$days_to_death),
                             clinical_data$days_to_last_follow_up,
                             clinical_data$days_to_death)

# Create survival status: 1 = dead, 0 = alive
clinical_data$status <- ifelse(clinical_data$vital_status == "Dead", 1, 0)

# Remove samples with missing survival info
valid_idx <- which(!is.na(clinical_data$time) & !is.na(clinical_data$status))
expr_data <- expr_data[, valid_idx]
clinical_data <- clinical_data[valid_idx, ]

cat("Final sample size:", ncol(expr_data), "\n")
cat("Events (deaths):", sum(clinical_data$status), "\n")
cat("Event rate:", round(mean(clinical_data$status) * 100, 1), "%\n")

# =============================================================================
# 4. COX REGRESSION ANALYSIS
# =============================================================================

cat("Starting Cox regression analysis...\n")

# Transpose expression matrix: rows = samples, columns = genes
expr_t <- t(expr_data)

# Combine expression with survival info
combine_data <- cbind(clinical_data[, c("time", "status")], expr_t)

# Remove samples with missing survival info
combine_data <- na.omit(combine_data)

# Filter out low-expression genes (expressed in <10% of samples)
expr_matrix <- combine_data[, -(1:2)]
keep_genes <- apply(expr_matrix, 2, function(x) mean(x > 0) > 0.1)
combine_data_filtered <- cbind(combine_data[, 1:2], expr_matrix[, keep_genes])
colnames(combine_data_filtered) <- make.names(colnames(combine_data_filtered), unique = TRUE)

cat("Genes after filtering:", ncol(combine_data_filtered) - 2, "\n")

# Initialize result container
cox_results <- data.frame(
  Gene = character(),
  HR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all genes
total_genes <- ncol(combine_data_filtered) - 2
cat("Analyzing", total_genes, "genes...\n")

for (i in 1:total_genes) {
  if (i %% 1000 == 0) cat("Processed", i, "genes\n")

  gene <- colnames(combine_data_filtered)[i + 2]
  gene <- trimws(gene)

  fit <- tryCatch({
    formula <- as.formula(paste0("Surv(time, status) ~ `", gene, "`"))
    coxph(formula, data = combine_data_filtered)
  }, error = function(e) {
    return(NULL)
  })

  if (!is.null(fit)) {
    summary_fit <- summary(fit)
    cox_results <- rbind(cox_results, data.frame(
      Gene = gene,
      HR = round(summary_fit$coefficients[1, 2], 3),
      CI_lower = round(summary_fit$conf.int[1, 3], 3),
      CI_upper = round(summary_fit$conf.int[1, 4], 3),
      p.value = summary_fit$coefficients[1, 5],
      stringsAsFactors = FALSE
    ))
  }
}

# Adjust p-values for multiple testing
cox_results$p.adj <- p.adjust(cox_results$p.value, method = "BH")

# Sort by p-value
cox_results <- cox_results %>% arrange(p.value)

# Top 10 genes
top10_genes <- head(cox_results, 10)
cat("\nTop 10 prognostic genes:\n")
print(top10_genes)

# Save results
write.csv(cox_results, "cox_results_all_genes.csv", row.names = FALSE)
write.csv(top10_genes, "cox_results_top10.csv", row.names = FALSE)

# =============================================================================
# 5. KAPLAN-MEIER PLOTS
# =============================================================================

cat("\nCreating Kaplan-Meier plots...\n")

# Loop over top 10 genes - FIXED VERSION
for (i in 1:nrow(top10_genes)) {
  gene <- top10_genes$Gene[i]

  tryCatch({
    expr_gene <- combine_data_filtered[[gene]]

    surv_df <- data.frame(
      time = combine_data_filtered$time,
      status = combine_data_filtered$status,
      gene_expr = expr_gene
    )

    surv_df <- na.omit(surv_df)

    # Skip if constant expression
    if (length(unique(surv_df$gene_expr)) < 2) {
      cat("Skipped gene (constant expression):", gene, "\n")
      next
    }

    # Create groups (median split)
    surv_df$group <- ifelse(surv_df$gene_expr > median(surv_df$gene_expr), "High", "Low")

    # Skip if single group
    if (length(unique(surv_df$group)) < 2) {
      cat("Skipped gene (single group after split):", gene, "\n")
      next
    }

    # Fit survival model
    fit <- survfit(Surv(time, status) ~ group, data = surv_df)

    # Create plot
    cat("Plotting KM curve for:", gene, "\n")
    plot_km <- ggsurvplot(
      fit,
      data = surv_df,
      title = paste("Kaplan-Meier Curve -", gene),
      subtitle = paste("HR =", top10_genes$HR[i], ", p =",
                       format(top10_genes$p.value[i], scientific = TRUE)),
      pval = TRUE,
      risk.table = TRUE,
      xlab = "Time (days)",
      ylab = "Survival Probability",
      legend.labs = c("High Expression", "Low Expression"),
      palette = c("#E64B35", "#4DBBD5"),
      ggtheme = theme_minimal()
    )

    print(plot_km)

    # Save plot
    ggsave(
      filename = paste0("KM_", gene, ".pdf"),
      plot = plot_km$plot,
      width = 8,
      height = 6
    )

  }, error = function(e) {
    cat("Error plotting gene", gene, ":", e$message, "\n")
  })
}

# =============================================================================
# 6. LASSO-COX REGRESSION
# =============================================================================

cat("\nPerforming Lasso-Cox regression...\n")

set.seed(2025)

# Filter genes with p < 0.01
signif_genes <- cox_results$Gene[cox_results$p.value < 0.01]
cat("Using", length(signif_genes), "significant genes for Lasso\n")

# Remove samples with time <= 0
combine_data_filtered <- combine_data_filtered[combine_data_filtered$time > 0, ]

# Construct matrix
x <- as.matrix(combine_data_filtered[, signif_genes])
y <- Surv(combine_data_filtered$time, combine_data_filtered$status)

# Lasso-Cox with cross-validation
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
best_lambda <- cv_fit$lambda.min

cat("Best lambda:", best_lambda, "\n")

# Plot cross-validation
plot(cv_fit)
title("Cross-Validation for Lasso-Cox")

# Fit final model
lasso_model <- glmnet(x, y, family = "cox", alpha = 1, lambda = best_lambda)
coef_df <- as.data.frame(as.matrix(coef(lasso_model)))
selected_genes <- rownames(coef_df)[coef_df[, 1] != 0]

cat("Genes selected by Lasso-Cox model:", length(selected_genes), "\n")
print(selected_genes)

# =============================================================================
# 7. RISK SCORE CALCULATION
# =============================================================================

cat("\nCalculating risk scores...\n")

# Calculate risk score
risk_score <- predict(lasso_model, newx = x, type = "link")

# Create risk groups
combine_risk <- data.frame(
  time = combine_data_filtered$time,
  status = combine_data_filtered$status,
  risk_score = as.numeric(risk_score)
)
combine_risk$risk_group <- ifelse(combine_risk$risk_score > median(combine_risk$risk_score),
                                  "High Risk", "Low Risk")

# Survival analysis for risk groups
fit <- survfit(Surv(time, status) ~ risk_group, data = combine_risk)

# Plot KM curve for risk groups
risk_km_plot <- ggsurvplot(
  fit,
  data = combine_risk,
  pval = TRUE,
  risk.table = TRUE,
  title = "Kaplan-Meier Curve by Lasso-based Risk Score",
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "Risk Group",
  palette = c("#4DBBD5", "#E64B35"),
  ggtheme = theme_minimal()
)

print(risk_km_plot)

# Save risk stratification plot
ggsave("risk_stratification_km.pdf", plot = risk_km_plot$plot, width = 8, height = 6)

# =============================================================================
# 8. COEFFICIENT ANALYSIS
# =============================================================================

cat("\nAnalyzing Lasso coefficients...\n")

# Check the Lasso coefficient for each gene
coef_matrix <- as.matrix(coef(lasso_model))
coef_table <- data.frame(
  Gene = rownames(coef_matrix),
  Coefficient = coef_matrix[, 1]
)
coef_table <- coef_table[coef_table$Coefficient != 0, ]  # Filter non-zero
coef_table <- coef_table[order(abs(coef_table$Coefficient), decreasing = TRUE), ]

cat("Top 10 genes by coefficient magnitude:\n")
print(head(coef_table, 10))

# Save coefficient table
write.csv(coef_table, "lasso_coefficients.csv", row.names = FALSE)

# Create risk score dataframe
risk_score_df <- data.frame(
  sample = rownames(x),
  risk_score = as.numeric(risk_score),
  risk_group = combine_risk$risk_group
)

# Save risk scores
write.csv(risk_score_df, "lasso_risk_scores.csv", row.names = FALSE)

# =============================================================================
# 9. MODEL PERFORMANCE EVALUATION
# =============================================================================

cat("\nEvaluating model performance...\n")

# Alternative C-index calculation function
calculate_cindex <- function(risk_scores, time, status) {
  # Create all pairs of observations
  n <- length(risk_scores)
  concordant <- 0
  discordant <- 0
  tied <- 0

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Only consider pairs where we can determine order
      if (status[i] == 1 || status[j] == 1) {
        # If both have events, use the one with shorter time
        if (status[i] == 1 && status[j] == 1) {
          if (time[i] < time[j]) {
            if (risk_scores[i] > risk_scores[j]) concordant <- concordant + 1
            else if (risk_scores[i] < risk_scores[j]) discordant <- discordant + 1
            else tied <- tied + 1
          } else if (time[i] > time[j]) {
            if (risk_scores[i] < risk_scores[j]) concordant <- concordant + 1
            else if (risk_scores[i] > risk_scores[j]) discordant <- discordant + 1
            else tied <- tied + 1
          }
        }
        # If only one has event, compare with censored
        else if (status[i] == 1 && status[j] == 0 && time[i] < time[j]) {
          if (risk_scores[i] > risk_scores[j]) concordant <- concordant + 1
          else if (risk_scores[i] < risk_scores[j]) discordant <- discordant + 1
          else tied <- tied + 1
        } else if (status[i] == 0 && status[j] == 1 && time[j] < time[i]) {
          if (risk_scores[j] > risk_scores[i]) concordant <- concordant + 1
          else if (risk_scores[j] < risk_scores[i]) discordant <- discordant + 1
          else tied <- tied + 1
        }
      }
    }
  }

  # Calculate C-index
  c_index <- (concordant + 0.5 * tied) / (concordant + discordant + tied)
  return(c_index)
}

# Try to use survcomp if available, otherwise use our function
if (requireNamespace("survcomp", quietly = TRUE)) {
  c_index <- survcomp::concordance.index(
    x = risk_score,
    surv.time = combine_risk$time,
    surv.event = combine_risk$status,
    method = "noether"
  )$c.index
} else {
  c_index <- calculate_cindex(as.numeric(risk_score), combine_risk$time, combine_risk$status)
}

cat("C-index:", round(c_index, 3), "\n")

# Bootstrap validation
set.seed(2025)
bootstrap_cindex <- replicate(100, {
  boot_idx <- sample(nrow(combine_risk), replace = TRUE)
  boot_data <- combine_risk[boot_idx, ]

  if (requireNamespace("survcomp", quietly = TRUE)) {
    survcomp::concordance.index(boot_data$risk_score, boot_data$time,
                                boot_data$status, method = "noether")$c.index
  } else {
    calculate_cindex(boot_data$risk_score, boot_data$time, boot_data$status)
  }
})

cat("Bootstrap C-index:", round(mean(bootstrap_cindex), 3), "±",
    round(sd(bootstrap_cindex), 3), "\n")

