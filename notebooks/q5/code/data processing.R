# ===================== STEP 0: LOAD LIBRARIES =====================

# Load essential libraries
library(readr)
library(dplyr)
library(TCGAbiolinks)
library(maftools)

# ===================== STEP 1: LOAD EXPRESSION DATA FROM UCSC XENA =====================

# Define path to UCSC Xena expression matrix
expr_path <- "../../dat/TCGA.LUAD.sampleMap_HiSeqV2"

# Check file existence before reading
if (!file.exists(expr_path)) {
  stop("❌ Expression file not found at: ", expr_path)
}

cat("📂 Reading expression data from:\n", expr_path, "\n")

# Load expression matrix (genes as rows, samples as columns)
expr_data <- read.delim(expr_path, header = TRUE, row.names = 1, check.names = FALSE)

# Filter to tumor samples only (sample barcode ends in '-01')
tumor_samples <- grep("-01$", colnames(expr_data), value = TRUE)
expr_data <- expr_data[, tumor_samples]

# Format TCGA sample barcodes (first 12 characters)
colnames(expr_data) <- substr(colnames(expr_data), 1, 12)

# Convert matrix to data frame format, with gene symbols as a column
expr_df <- data.frame(sample = rownames(expr_data), expr_data)

# Save cleaned expression data
write_csv(expr_df, "LUAD_expression_data.csv")
cat("✅ Expression data saved to LUAD_expression_data.csv with", nrow(expr_df), "genes and", ncol(expr_df) - 1, "samples.\n\n")

# ===================== STEP 2: DOWNLOAD CLINICAL DATA FROM GDC =====================

# Query TCGA-LUAD clinical supplement (Biotab format)
query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# Download and prepare data
GDCdownload(query)
clinical_data <- GDCprepare(query)

# Extract patient-level clinical table
patient_df <- clinical_data$clinical_patient_luad

# ===================== STEP 3: CLEAN CLINICAL DATA =====================

# Clean and harmonize clinical variables
real_clinical <- patient_df %>%
  transmute(
    tcga_id = bcr_patient_barcode,
    age = as.numeric(age_at_initial_pathologic_diagnosis),
    gender = gender,
    ecog_ps = as.numeric(ecog_score),  # ✅ Now using real ECOG scores!
    stage = coalesce(ajcc_pathologic_tumor_stage, ajcc_clinical_tumor_stage),
    smoking_status = case_when(
      tobacco_smoking_history_indicator == "1" ~ "Never",
      tobacco_smoking_history_indicator == "2" ~ "Current",
      tobacco_smoking_history_indicator == "3" ~ "Former",
      tobacco_smoking_history_indicator == "4" ~ "Occasional",
      TRUE ~ NA_character_
    ),
    vital_status = vital_status,
    days_to_death = as.numeric(death_days_to),
    days_to_last_followup = as.numeric(last_contact_days_to)
  )


# ===================== STEP 4: COMPUTE PD-L1 (CD274) EXPRESSION =====================

# Load expression matrix again (if needed)
expr_data <- read.delim(expr_path, header = TRUE, row.names = 1, check.names = FALSE)

# Keep only tumor samples
tumor_samples <- grep("-01$", colnames(expr_data), value = TRUE)
expr_data <- expr_data[, tumor_samples]

# Format sample barcodes
colnames(expr_data) <- substr(colnames(expr_data), 1, 12)

# Locate CD274 (PD-L1) row in matrix
cd274_row <- which(rownames(expr_data) == "CD274")

# Check if CD274 exists in the dataset
if (length(cd274_row) == 0) {
  stop("❌ CD274 not found in expression matrix. Check gene names.")
}

# Extract and transpose PD-L1 expression values
cd274_expr <- expr_data[cd274_row, , drop = FALSE] %>% t()

# Convert to data frame format
cd274_df <- data.frame(
  sample = rownames(cd274_expr),
  expression = as.numeric(cd274_expr[,1]),
  stringsAsFactors = FALSE
)

# Quantile binning of PD-L1 expression into High/Medium/Low
cd274_df <- cd274_df %>%
  mutate(
    pdl1_expression = case_when(
      expression >= quantile(expression, 0.67, na.rm = TRUE) ~ "High",
      expression <= quantile(expression, 0.33, na.rm = TRUE) ~ "Low",
      TRUE ~ "Medium"
    )
  )

# Save PD-L1 expression levels
write_csv(cd274_df, "cd274_pdl1_expression.csv")
cat("✅ PD-L1 expression levels saved as cd274_pdl1_expression.csv with", nrow(cd274_df), "samples.\n")

# ===================== STEP 5: SAVE BIOMARKER FILES =====================

# Save EGFR mutation and ALK fusion data (assumes created earlier)
write_csv(egfr_status, "egfr_mutation_status.csv")
write_csv(alk_status,  "alk_fusion_status.csv")

# Save PD-L1 expression (redundant save to ensure availability)
write_csv(cd274_df, "cd274_pdl1_expression.csv")
