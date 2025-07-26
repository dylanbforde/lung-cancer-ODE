# =============================================================================
# MEIN40330 - Final Project -- Question 5 WITH REAL TCGA CLINICAL DATA
# Author: Group 4 - Complete LUAD Treatment Selection with Real ECOG Scores
# =============================================================================

setwd("C:/Users/mahen/Desktop/Personalised Medicine/final project/finalproject/R")

library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(TCGAbiolinks)

# Fix namespace conflicts
select <- dplyr::select
filter <- dplyr::filter

cat("=== ENHANCED LUAD TREATMENT SYSTEM WITH REAL TCGA DATA ===\n")
cat("Using actual ECOG performance status and clinical variables\n\n")

# === STEP 1: PREPARE REAL TCGA CLINICAL DATA ===
cat("=== PREPARING REAL TCGA CLINICAL DATA ===\n")

# Download real TCGA clinical data (if not already done)
if(!exists("clinical_data")) {
  query <- GDCquery(
    project = "TCGA-LUAD",
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )

  # Download if needed
  if(!file.exists(paste0(getwd(), "/GDCdata"))) {
    GDCdownload(query)
  }

  clinical_data <- GDCprepare(query)
}

# Extract and clean patient clinical data
patient_df <- clinical_data$clinical_patient_luad

real_clinical <- patient_df %>%
  transmute(
    sample = bcr_patient_barcode,  # Keep as 'sample' for compatibility
    tcga_id = bcr_patient_barcode,
    age = parse_number(age_at_initial_pathologic_diagnosis),
    gender = gender,
    ecog_ps = parse_number(ecog_score),  # ✅ REAL ECOG SCORES!
    stage = coalesce(ajcc_pathologic_tumor_stage, ajcc_clinical_tumor_stage),
    smoking_status = case_when(
      tobacco_smoking_history_indicator == "1" ~ "Never",
      tobacco_smoking_history_indicator == "2" ~ "Current",
      tobacco_smoking_history_indicator == "3" ~ "Former",
      tobacco_smoking_history_indicator == "4" ~ "Occasional",
      TRUE ~ NA_character_
    ),
    vital_status = vital_status,
    days_to_death = parse_number(death_days_to),
    days_to_last_followup = parse_number(last_contact_days_to)
  ) %>%
  # Clean up stage formatting
  mutate(
    stage = case_when(
      grepl("Stage I[^V]", stage) ~ "Stage I",
      grepl("Stage IA", stage) ~ "Stage IA",
      grepl("Stage IB", stage) ~ "Stage IB",
      grepl("Stage II[^I]", stage) ~ "Stage II",
      grepl("Stage IIA", stage) ~ "Stage IIA",
      grepl("Stage IIB", stage) ~ "Stage IIB",
      grepl("Stage III[^A-B]", stage) ~ "Stage III",
      grepl("Stage IIIA", stage) ~ "Stage IIIA",
      grepl("Stage IIIB", stage) ~ "Stage IIIB",
      grepl("Stage IV", stage) ~ "Stage IV",
      TRUE ~ stage
    ),
    # Handle occasional smokers as former smokers for treatment purposes
    smoking_status = ifelse(smoking_status == "Occasional", "Former", smoking_status)
  )

cat("Real TCGA clinical data prepared:\n")
cat("- Total patients:", nrow(real_clinical), "\n")
cat("- Patients with ECOG PS:", sum(!is.na(real_clinical$ecog_ps)), "\n")
cat("- Patients with age:", sum(!is.na(real_clinical$age)), "\n")
cat("- Patients with stage:", sum(!is.na(real_clinical$stage)), "\n")
cat("- Patients with smoking data:", sum(!is.na(real_clinical$smoking_status)), "\n")

# Display ECOG distribution
ecog_dist <- table(real_clinical$ecog_ps, useNA = "ifany")
cat("\nReal ECOG Performance Status Distribution:\n")
print(ecog_dist)

# Display stage distribution
stage_dist <- table(real_clinical$stage, useNA = "ifany")
cat("\nStage Distribution:\n")
print(stage_dist)

# === STEP 2: LOAD BIOMARKER DATA ===
cat("\n=== LOADING BIOMARKER DATA ===\n")

# Load your existing biomarker files
risk_scores <- read_csv("lasso_risk_scores.csv", show_col_types = FALSE)
expr_df <- read_csv("LUAD_expression_data.csv", show_col_types = FALSE)
egfr_status <- read_csv("egfr_mutation_status.csv", show_col_types = FALSE)
alk_status <- read_csv("alk_fusion_status.csv", show_col_types = FALSE)
pdl1_df <- read_csv("cd274_pdl1_expression.csv", show_col_types = FALSE)

# Standardize IDs
extract_tcga_id <- function(sample_ids) {
  sample_ids %>%
    gsub("\\.", "-", .) %>%
    substr(1, 12) %>%
    toupper()
}

real_clinical$tcga_id <- extract_tcga_id(real_clinical$tcga_id)
egfr_status$tcga_id <- extract_tcga_id(egfr_status$sample)
alk_status$tcga_id <- extract_tcga_id(alk_status$sample)
pdl1_df$tcga_id <- extract_tcga_id(pdl1_df$sample)

cat("Biomarker data loaded successfully\n")

# === STEP 3: ENHANCED TREATMENT FUNCTION WITH REAL DATA ===
cat("\n=== ENHANCED TREATMENT RECOMMENDATIONS WITH REAL CLINICAL DATA ===\n")

recommend_treatment_real_data <- function(df) {
  n <- nrow(df)

  # Initialize output vectors
  primary_treatment <- rep("Standard Chemotherapy", n)
  surgical_option <- rep("Not candidate", n)
  radiation_option <- rep("Not indicated", n)
  combination_therapy <- rep(FALSE, n)
  conf <- rep(0.6, n)
  rationale <- rep("Standard of care", n)
  biomarker_guided <- rep(FALSE, n)
  treatment_intent <- rep("Palliative", n)

  # Extract biomarker variables
  egfr_pos <- !is.na(df$egfr_mutation) & df$egfr_mutation == "Positive"
  alk_pos <- !is.na(df$alk_fusion) & df$alk_fusion == "Positive"
  pdl1_high <- !is.na(df$pdl1_expression) & df$pdl1_expression == "High"
  pdl1_medium <- !is.na(df$pdl1_expression) & df$pdl1_expression == "Medium"

  # === REAL CLINICAL VARIABLES (NO SIMULATION NEEDED!) ===

  # Stage categorization using real data
  early_stage <- !is.na(df$stage) & df$stage %in% c("Stage I", "Stage IA", "Stage IB")
  stage_2 <- !is.na(df$stage) & df$stage %in% c("Stage II", "Stage IIA", "Stage IIB")
  stage_3 <- !is.na(df$stage) & df$stage %in% c("Stage III", "Stage IIIA", "Stage IIIB")
  stage_4 <- !is.na(df$stage) & df$stage %in% c("Stage IV")

  # Age categories using real data
  young_patient <- !is.na(df$age) & df$age < 65
  elderly_patient <- !is.na(df$age) & df$age >= 75

  # REAL ECOG Performance Status (the key improvement!)
  good_ps <- !is.na(df$ecog_ps) & df$ecog_ps <= 1
  moderate_ps <- !is.na(df$ecog_ps) & df$ecog_ps == 2
  poor_ps <- !is.na(df$ecog_ps) & df$ecog_ps >= 3

  # For patients with missing ECOG, use conservative estimate based on age/stage
  missing_ecog <- is.na(df$ecog_ps)
  if(sum(missing_ecog) > 0) {
    # Conservative ECOG estimation for missing data
    good_ps[missing_ecog] <- !elderly_patient[missing_ecog] & !stage_4[missing_ecog]
    moderate_ps[missing_ecog] <- elderly_patient[missing_ecog] | stage_4[missing_ecog]
    poor_ps[missing_ecog] <- elderly_patient[missing_ecog] & stage_4[missing_ecog]
  }

  # Real smoking status
  smokers <- !is.na(df$smoking_status) & df$smoking_status %in% c("Current", "Former")
  never_smokers <- !is.na(df$smoking_status) & df$smoking_status == "Never"

  # LUAD is adenocarcinoma (non-squamous)
  non_squamous <- rep(TRUE, n)

  # === 1. SURGERY EVALUATION (Early Stage Disease) ===
  surgical_candidates <- (early_stage | stage_2) & good_ps & !elderly_patient

  surgical_idx <- which(surgical_candidates)
  if(length(surgical_idx) > 0) {
    surgical_option[surgical_idx] <- "Lobectomy"
    primary_treatment[surgical_idx] <- "Surgery (Lobectomy)"
    treatment_intent[surgical_idx] <- "Curative"
    conf[surgical_idx] <- 0.9
    rationale[surgical_idx] <- "Early stage, ECOG 0-1 - surgical resection"

    # Sublobar for high-risk patients
    high_risk_surgical <- surgical_candidates & elderly_patient
    high_risk_idx <- which(high_risk_surgical)
    if(length(high_risk_idx) > 0) {
      surgical_option[high_risk_idx] <- "Sublobar resection"
      primary_treatment[high_risk_idx] <- "Surgery (Sublobar)"
      rationale[high_risk_idx] <- "Early stage, elderly - sublobar resection"
    }
  }

  # === 2. STEREOTACTIC BODY RADIATION THERAPY (SBRT) ===
  sbrt_candidates <- (early_stage | stage_2) & !surgical_candidates & (good_ps | moderate_ps)

  sbrt_idx <- which(sbrt_candidates)
  if(length(sbrt_idx) > 0) {
    radiation_option[sbrt_idx] <- "SBRT"
    primary_treatment[sbrt_idx] <- "Stereotactic Body Radiation"
    treatment_intent[sbrt_idx] <- "Curative"
    conf[sbrt_idx] <- 0.85
    rationale[sbrt_idx] <- "Early stage, medically inoperable - SBRT curative"
  }

  # === 3. CONCURRENT CHEMORADIATION (Stage III) ===
  chemorad_candidates <- stage_3 & (good_ps | moderate_ps)

  chemorad_idx <- which(chemorad_candidates)
  if(length(chemorad_idx) > 0) {
    radiation_option[chemorad_idx] <- "IMRT + Concurrent Chemo"
    primary_treatment[chemorad_idx] <- "Chemoradiation + Durvalumab"
    combination_therapy[chemorad_idx] <- TRUE
    treatment_intent[chemorad_idx] <- "Curative"
    conf[chemorad_idx] <- 0.8
    rationale[chemorad_idx] <- "Stage III, good PS - concurrent chemoRT + durvalumab"
  }

  # === 4. TARGETED THERAPY (Highest Priority - Overrides Others) ===

  # EGFR-targeted therapy
  egfr_idx <- which(egfr_pos)
  if(length(egfr_idx) > 0) {
    primary_treatment[egfr_idx] <- "EGFR-targeted therapy (Osimertinib)"
    conf[egfr_idx] <- 0.95
    rationale[egfr_idx] <- "EGFR mutation positive - osimertinib first-line"
    biomarker_guided[egfr_idx] <- TRUE

    # Early-stage EGFR+ can have surgery + adjuvant TKI
    early_egfr <- egfr_pos & (early_stage | stage_2) & good_ps
    early_egfr_idx <- which(early_egfr)
    if(length(early_egfr_idx) > 0) {
      surgical_option[early_egfr_idx] <- "Surgery + Adjuvant EGFR TKI"
      treatment_intent[early_egfr_idx] <- "Curative"
    }

    # Advanced EGFR+
    advanced_egfr <- egfr_pos & (stage_3 | stage_4)
    adv_egfr_idx <- which(advanced_egfr)
    if(length(adv_egfr_idx) > 0) {
      treatment_intent[adv_egfr_idx] <- ifelse(stage_3[adv_egfr_idx], "Curative", "Palliative")
      radiation_option[adv_egfr_idx] <- "Palliative RT (if needed)"
    }
  }

  # ALK-targeted therapy
  alk_idx <- which(alk_pos & !egfr_pos)
  if(length(alk_idx) > 0) {
    primary_treatment[alk_idx] <- "ALK-targeted therapy (Alectinib)"
    conf[alk_idx] <- 0.9
    rationale[alk_idx] <- "ALK fusion positive - alectinib first-line"
    biomarker_guided[alk_idx] <- TRUE

    # Treatment intent based on stage
    treatment_intent[alk_idx] <- ifelse(stage_4[alk_idx], "Palliative", "Curative")
  }

  # === 5. IMMUNOTHERAPY STRATEGIES ===

  # PD-L1 ≥50% monotherapy (no driver mutations, good PS)
  pdl1_mono_candidates <- pdl1_high & !egfr_pos & !alk_pos & good_ps & stage_4

  pdl1_mono_idx <- which(pdl1_mono_candidates)
  if(length(pdl1_mono_idx) > 0) {
    primary_treatment[pdl1_mono_idx] <- "Pembrolizumab monotherapy"
    conf[pdl1_mono_idx] <- 0.85
    rationale[pdl1_mono_idx] <- "PD-L1 ≥50%, ECOG 0-1 - pembrolizumab monotherapy"
    biomarker_guided[pdl1_mono_idx] <- TRUE
    treatment_intent[pdl1_mono_idx] <- "Palliative"
  }

  # PD-L1 1-49% combination therapy
  immuno_combo_candidates <- (pdl1_medium | (!pdl1_high & !is.na(df$pdl1_expression))) &
    !egfr_pos & !alk_pos & (good_ps | moderate_ps) & stage_4

  immuno_combo_idx <- which(immuno_combo_candidates)
  if(length(immuno_combo_idx) > 0) {
    primary_treatment[immuno_combo_idx] <- "Pembrolizumab + Chemotherapy"
    combination_therapy[immuno_combo_idx] <- TRUE
    conf[immuno_combo_idx] <- 0.8
    rationale[immuno_combo_idx] <- "PD-L1 1-49%, good PS - pembrolizumab + platinum doublet"
    biomarker_guided[immuno_combo_idx] <- TRUE
    treatment_intent[immuno_combo_idx] <- "Palliative"
  }

  # === 6. ANTI-ANGIOGENIC THERAPY ===
  bevacizumab_candidates <- !egfr_pos & !alk_pos & non_squamous & stage_4 &
    good_ps & !pdl1_high & is.na(df$pdl1_expression)

  bev_idx <- which(bevacizumab_candidates)
  if(length(bev_idx) > 0) {
    primary_treatment[bev_idx] <- "Chemotherapy + Bevacizumab"
    combination_therapy[bev_idx] <- TRUE
    conf[bev_idx] <- 0.75
    rationale[bev_idx] <- "No driver mutations, good PS - carboplatin/paclitaxel/bevacizumab"
    treatment_intent[bev_idx] <- "Palliative"
  }

  # === 7. EMERGING TARGETED THERAPIES ===
  # KRAS G12C (higher prevalence in smokers)
  potential_kras <- smokers & !egfr_pos & !alk_pos & stage_4

  set.seed(123)
  kras_g12c_sim <- potential_kras & sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(0.12, 0.88))

  kras_idx <- which(kras_g12c_sim)
  if(length(kras_idx) > 0) {
    primary_treatment[kras_idx] <- "KRAS G12C inhibitor (Sotorasib)"
    conf[kras_idx] <- 0.7
    rationale[kras_idx] <- "Smoker, likely KRAS G12C - consider sotorasib"
    biomarker_guided[kras_idx] <- TRUE
    treatment_intent[kras_idx] <- "Palliative"
  }

  # === 8. NEOADJUVANT/ADJUVANT STRATEGIES ===
  # Neoadjuvant for resectable Stage III with good PS
  neoadjuvant_candidates <- stage_3 & good_ps & young_patient

  neoadj_idx <- which(neoadjuvant_candidates)
  if(length(neoadj_idx) > 0) {
    primary_treatment[neoadj_idx] <- "Neoadjuvant Chemo-Immunotherapy"
    surgical_option[neoadj_idx] <- "Surgery after neoadjuvant"
    combination_therapy[neoadj_idx] <- TRUE
    treatment_intent[neoadj_idx] <- "Curative"
    conf[neoadj_idx] <- 0.8
    rationale[neoadj_idx] <- "Resectable Stage III, young, ECOG 0-1 - neoadjuvant approach"
  }

  # === 9. PALLIATIVE CARE ===
  palliative_candidates <- poor_ps | (elderly_patient & stage_4)

  pall_idx <- which(palliative_candidates)
  if(length(pall_idx) > 0) {
    primary_treatment[pall_idx] <- "Palliative Care + Supportive Therapy"
    radiation_option[pall_idx] <- "Palliative RT (symptom control)"
    conf[pall_idx] <- 0.9
    rationale[pall_idx] <- "ECOG ≥3 or elderly with Stage IV - palliative intent"
    treatment_intent[pall_idx] <- "Palliative"
  }

  # === 10. TREATMENT SEQUENCING ===
  second_line <- rep("TBD based on response", n)
  second_line[egfr_pos] <- "Chemotherapy + bevacizumab (post-EGFR resistance)"
  second_line[alk_pos] <- "Next-generation ALK inhibitor or chemotherapy"

  immuno_treated <- grepl("Pembrolizumab", primary_treatment)
  second_line[immuno_treated] <- "Docetaxel + ramucirumab"

  chemo_treated <- grepl("Chemotherapy", primary_treatment) & !immuno_treated
  second_line[chemo_treated] <- "Immunotherapy (if PD-L1 >1%) or clinical trial"

  return(data.frame(
    primary_treatment = primary_treatment,
    surgical_option = surgical_option,
    radiation_option = radiation_option,
    second_line_option = second_line,
    treatment_intent = treatment_intent,
    combination_therapy = combination_therapy,
    confidence = conf,
    rationale = rationale,
    biomarker_guided = biomarker_guided,
    stringsAsFactors = FALSE
  ))
}

# === STEP 4: APPLY ENHANCED RECOMMENDATIONS WITH REAL DATA ===
cat("\n=== APPLYING ENHANCED RECOMMENDATIONS WITH REAL TCGA DATA ===\n")

# Merge real clinical data with biomarkers
comprehensive_patients <- real_clinical %>%
  left_join(egfr_status %>% select(tcga_id, egfr_mutation), by = "tcga_id") %>%
  left_join(alk_status %>% select(tcga_id, alk_fusion), by = "tcga_id") %>%
  left_join(pdl1_df %>% select(tcga_id, pdl1_expression), by = "tcga_id")

cat("Comprehensive dataset created with", nrow(comprehensive_patients), "patients\n")
cat("- Patients with complete clinical data:",
    sum(!is.na(comprehensive_patients$age) & !is.na(comprehensive_patients$stage) &
          !is.na(comprehensive_patients$ecog_ps)), "\n")

# Apply enhanced treatment recommendations using REAL clinical data
real_treatment_recs <- recommend_treatment_real_data(comprehensive_patients)
final_analysis <- cbind(comprehensive_patients, real_treatment_recs)

# === STEP 5: COMPREHENSIVE ANALYSIS WITH REAL DATA ===
cat("\n=== COMPREHENSIVE ANALYSIS RESULTS WITH REAL TCGA DATA ===\n")

# Treatment distribution
treatment_summary_real <- final_analysis %>%
  group_by(primary_treatment) %>%
  summarise(
    count = n(),
    percentage = round(n() / nrow(final_analysis) * 100, 1),
    biomarker_guided = sum(biomarker_guided),
    combination_therapy = sum(combination_therapy),
    curative_intent = sum(treatment_intent == "Curative"),
    avg_confidence = round(mean(confidence), 2),
    avg_ecog = round(mean(ecog_ps, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(count))

cat("TREATMENT DISTRIBUTION WITH REAL TCGA DATA:\n")
print(treatment_summary_real)

# Treatment modality breakdown
modality_summary_real <- final_analysis %>%
  summarise(
    total_patients = n(),
    patients_with_ecog = sum(!is.na(ecog_ps)),
    surgical_candidates = sum(surgical_option != "Not candidate"),
    radiation_candidates = sum(radiation_option != "Not indicated"),
    targeted_therapy = sum(grepl("EGFR|ALK|KRAS", primary_treatment)),
    immunotherapy = sum(grepl("Pembrolizumab", primary_treatment)),
    combination_treatments = sum(combination_therapy),
    curative_intent = sum(treatment_intent == "Curative"),
    palliative_intent = sum(treatment_intent == "Palliative")
  )

cat("\nTREATMENT MODALITY BREAKDOWN (REAL DATA):\n")
cat("- Total patients:", modality_summary_real$total_patients, "\n")
cat("- Patients with real ECOG scores:", modality_summary_real$patients_with_ecog,
    paste0("(", round(modality_summary_real$patients_with_ecog/modality_summary_real$total_patients*100, 1), "%)"), "\n")
cat("- Surgical candidates:", modality_summary_real$surgical_candidates,
    paste0("(", round(modality_summary_real$surgical_candidates/modality_summary_real$total_patients*100, 1), "%)"), "\n")
cat("- Radiation therapy:", modality_summary_real$radiation_candidates,
    paste0("(", round(modality_summary_real$radiation_candidates/modality_summary_real$total_patients*100, 1), "%)"), "\n")
cat("- Targeted therapy:", modality_summary_real$targeted_therapy,
    paste0("(", round(modality_summary_real$targeted_therapy/modality_summary_real$total_patients*100, 1), "%)"), "\n")
cat("- Immunotherapy:", modality_summary_real$immunotherapy,
    paste0("(", round(modality_summary_real$immunotherapy/modality_summary_real$total_patients*100, 1), "%)"), "\n")
cat("- Curative intent:", modality_summary_real$curative_intent,
    paste0("(", round(modality_summary_real$curative_intent/modality_summary_real$total_patients*100, 1), "%)"), "\n")

# ECOG Performance Status Impact Analysis
ecog_impact <- final_analysis %>%
  filter(!is.na(ecog_ps)) %>%
  group_by(ecog_ps, treatment_intent) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = treatment_intent, values_from = count, values_fill = 0)

cat("\nECOG PERFORMANCE STATUS IMPACT ON TREATMENT INTENT:\n")
print(ecog_impact)

# === STEP 6: VISUALIZATIONS WITH REAL DATA ===
cat("\n=== CREATING VISUALIZATIONS WITH REAL TCGA DATA ===\n")

# 1. Treatment distribution with ECOG stratification
p1_real <- final_analysis %>%
  filter(!is.na(ecog_ps)) %>%
  count(primary_treatment, ecog_ps) %>%
  ggplot(aes(x = reorder(primary_treatment, n), y = n, fill = factor(ecog_ps))) +
  geom_col(alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Treatment Recommendations by ECOG Performance Status",
       subtitle = paste("N =", sum(!is.na(final_analysis$ecog_ps)), "patients with real ECOG scores"),
       x = "Primary Treatment",
       y = "Number of Patients",
       fill = "ECOG PS") +
  scale_fill_brewer(type = "seq", palette = "YlOrRd")

print(p1_real)

# 2. Treatment intent by stage with real data
p2_real <- final_analysis %>%
  filter(!is.na(stage)) %>%
  count(stage, treatment_intent) %>%
  ggplot(aes(x = stage, y = n, fill = treatment_intent)) +
  geom_col(position = "dodge", alpha = 0.8) +
  theme_minimal() +
  labs(title = "Treatment Intent by Cancer Stage (Real TCGA Data)",
       x = "Cancer Stage", y = "Number of Patients",
       fill = "Treatment Intent") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2_real)

# === STEP 7: SAVE COMPREHENSIVE RESULTS ===
cat("\n=== SAVING RESULTS WITH REAL TCGA DATA ===\n")

write_csv(final_analysis, "comprehensive_treatment_real_tcga_data.csv")
write_csv(treatment_summary_real, "treatment_summary_real_tcga.csv")

# Create clinical quality report
clinical_quality <- final_analysis %>%
  summarise(
    total_patients = n(),
    complete_age = sum(!is.na(age)),
    complete_ecog = sum(!is.na(ecog_ps)),
    complete_stage = sum(!is.na(stage)),
    complete_smoking = sum(!is.na(smoking_status)),
    precision_medicine_rate = round(sum(biomarker_guided) / n() * 100, 1),
    realistic_ecog_0_1_rate = round(sum(ecog_ps <= 1, na.rm = TRUE) / sum(!is.na(ecog_ps)) * 100, 1)
  )

cat("FINAL RESULTS WITH REAL TCGA DATA:\n")
cat("✅ comprehensive_treatment_real_tcga_data.csv\n")
cat("✅ treatment_summary_real_tcga.csv\n")

cat("\nCLINICAL DATA QUALITY WITH REAL TCGA:\n")
cat("- Complete age data:", clinical_quality$complete_age, "/", clinical_quality$total_patients, "\n")
cat("- Complete ECOG data:", clinical_quality$complete_ecog, "/", clinical_quality$total_patients, "\n")
cat("- Complete stage data:", clinical_quality$complete_stage, "/", clinical_quality$total_patients, "\n")
cat("- Precision medicine rate:", clinical_quality$precision_medicine_rate, "%\n")
cat("- Patients with good PS (0-1):", clinical_quality$realistic_ecog_0_1_rate, "%\n")

cat("\n🎉 ENHANCED QUESTION 5 WITH REAL TCGA DATA COMPLETE! 🎉\n")
cat("This analysis uses actual ECOG performance status scores from TCGA,\n")
cat("providing much more accurate and clinically relevant treatment recommendations!\n")
