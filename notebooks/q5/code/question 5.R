# =============================================================================
# IMPROVED LUAD TREATMENT RECOMMENDATION SYSTEM
# Modular Architecture with Better Code Organization
# =============================================================================

# =============================================================================
# 1. SETUP AND CONFIGURATION
# =============================================================================

# Project configuration
project_config <- list(
  base_dir = here::here(),
  data_dir = file.path(here::here(), "data"),
  results_dir = file.path(here::here(), "results"),
  plots_dir = file.path(here::here(), "results", "plots"),
  processed_dir = file.path(here::here(), "processed"),  # Changed path

  # Analysis parameters
  biomarker_thresholds = list(
    pdl1_high = 50,
    pdl1_medium = 1,
    egfr_confidence = 0.95,
    alk_confidence = 0.90
  ),

  # Clinical thresholds
  clinical_thresholds = list(
    elderly_age = 75,
    young_age = 65,
    good_ecog = 1,
    moderate_ecog = 2,
    poor_ecog = 3
  )
)

# Create directory structure
create_project_structure <- function(config) {
  dirs_to_create <- c(
    config$results_dir,
    config$plots_dir,
    config$processed_dir,  # Add this
    file.path(config$results_dir, "reports"),
    file.path(config$results_dir, "exports")
  )

  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Load required packages with better error handling
load_required_packages <- function() {
  required_packages <- c(
    "dplyr", "readr", "ggplot2", "tidyr", "stringr",
    "TCGAbiolinks", "survival", "survminer", "here"
  )

  # Optional packages (will work without these but with reduced functionality)
  optional_packages <- c("DT", "plotly", "knitr", "rmarkdown", "jsonlite", "R6")

  # Install missing required packages
  missing_required <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_required) > 0) {
    cat("Installing missing required packages:", paste(missing_required, collapse = ", "), "\n")
    if ("BiocManager" %in% missing_required || "TCGAbiolinks" %in% missing_required) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(missing_required[missing_required %in% c("TCGAbiolinks", "survival", "survminer")])
      cran_packages <- missing_required[!missing_required %in% c("TCGAbiolinks", "survival", "survminer")]
      if (length(cran_packages) > 0) {
        install.packages(cran_packages)
      }
    } else {
      install.packages(missing_required)
    }
  }

  # Install missing optional packages
  missing_optional <- optional_packages[!sapply(optional_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_optional) > 0) {
    cat("Installing missing optional packages:", paste(missing_optional, collapse = ", "), "\n")
    tryCatch({
      install.packages(missing_optional)
    }, error = function(e) {
      cat("Could not install optional packages:", e$message, "\n")
      cat("System will work with reduced functionality\n")
    })
  }

  # Load packages
  for (pkg in required_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      stop(paste("Required package", pkg, "could not be loaded"))
    }
  }

  # Load optional packages
  for (pkg in optional_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }

  # Fix namespace conflicts
  select <- dplyr::select
  filter <- dplyr::filter

  cat("All required packages loaded successfully\n")
}

# =============================================================================
# 2. DATA LOADING MODULE
# =============================================================================

DataLoader <- R6::R6Class("DataLoader",
                          public = list(
                            config = NULL,

                            initialize = function(config) {
                              self$config <- config
                            },

                            # Load and clean TCGA clinical data
                            load_tcga_clinical = function(force_reload = FALSE) {
                              # Ensure processed directory exists
                              if (!dir.exists(self$config$processed_dir)) {
                                dir.create(self$config$processed_dir, recursive = TRUE)
                                cat("Created processed directory:", self$config$processed_dir, "\n")
                              }

                              clinical_file <- file.path(self$config$processed_dir, "tcga_clinical_clean.rds")

                              if (file.exists(clinical_file) && !force_reload) {
                                cat("Loading cached clinical data...\n")
                                return(readRDS(clinical_file))
                              }

                              cat("Processing TCGA clinical data...\n")

                              # Load raw clinical data
                              tryCatch({
                                query <- GDCquery(
                                  project = "TCGA-LUAD",
                                  data.category = "Clinical",
                                  data.type = "Clinical Supplement",
                                  data.format = "BCR Biotab"
                                )

                                if (!file.exists(file.path(getwd(), "GDCdata"))) {
                                  cat("Downloading clinical data...\n")
                                  GDCdownload(query)
                                }

                                clinical_data <- GDCprepare(query)
                                patient_df <- clinical_data$clinical_patient_luad

                              }, error = function(e) {
                                cat("Error loading TCGA data:", e$message, "\n")
                                cat("Attempting to use existing local data...\n")

                                # Try to load from existing GDCdata folder
                                gdc_path <- file.path(getwd(), "GDCdata")
                                if (dir.exists(gdc_path)) {
                                  # Look for existing clinical data files
                                  clinical_files <- list.files(gdc_path, pattern = "clinical.*luad", recursive = TRUE, full.names = TRUE)
                                  if (length(clinical_files) > 0) {
                                    # Use the first available clinical file
                                    patient_df <- read.delim(clinical_files[1], stringsAsFactors = FALSE)
                                  } else {
                                    stop("No clinical data files found. Please ensure TCGA data is available.")
                                  }
                                } else {
                                  stop("Unable to load clinical data and no local files found.")
                                }
                              })

                              # Clean and standardize clinical data with better error handling
                              clean_clinical <- patient_df %>%
                                transmute(
                                  patient_id = bcr_patient_barcode,
                                  tcga_id = bcr_patient_barcode,
                                  age = suppressWarnings(as.numeric(as.character(age_at_initial_pathologic_diagnosis))),
                                  gender = case_when(
                                    tolower(as.character(gender)) == "male" ~ "Male",
                                    tolower(as.character(gender)) == "female" ~ "Female",
                                    TRUE ~ NA_character_
                                  ),
                                  ecog_ps = suppressWarnings(as.numeric(as.character(ecog_score))),
                                  stage = self$standardize_stage(coalesce(
                                    as.character(ajcc_pathologic_tumor_stage),
                                    as.character(ajcc_clinical_tumor_stage)
                                  )),
                                  smoking_status = self$standardize_smoking(as.character(tobacco_smoking_history_indicator)),
                                  vital_status = as.character(vital_status),
                                  days_to_death = suppressWarnings(as.numeric(as.character(death_days_to))),
                                  days_to_last_followup = suppressWarnings(as.numeric(as.character(last_contact_days_to)))
                                ) %>%
                                # Add derived variables
                                mutate(
                                  age_group = case_when(
                                    is.na(age) ~ "Unknown",
                                    age < self$config$clinical_thresholds$young_age ~ "Young",
                                    age >= self$config$clinical_thresholds$elderly_age ~ "Elderly",
                                    TRUE ~ "Middle-aged"
                                  ),
                                  ecog_group = case_when(
                                    is.na(ecog_ps) ~ "Unknown",
                                    ecog_ps <= self$config$clinical_thresholds$good_ecog ~ "Good",
                                    ecog_ps == self$config$clinical_thresholds$moderate_ecog ~ "Moderate",
                                    ecog_ps >= self$config$clinical_thresholds$poor_ecog ~ "Poor",
                                    TRUE ~ "Unknown"
                                  ),
                                  survival_months = pmax(
                                    coalesce(days_to_death, 0),
                                    coalesce(days_to_last_followup, 0),
                                    na.rm = TRUE
                                  ) / 30.44,
                                  death_event = as.numeric(vital_status == "Dead")
                                ) %>%
                                filter(!is.na(patient_id) & patient_id != "")

                              # Save cleaned data
                              tryCatch({
                                saveRDS(clean_clinical, clinical_file)
                                cat("Clinical data processed and cached to:", clinical_file, "\n")
                              }, error = function(e) {
                                cat("Warning: Could not cache clinical data:", e$message, "\n")
                              })

                              cat("Clinical data processing complete:", nrow(clean_clinical), "patients\n")
                              return(clean_clinical)
                            },

                            # Load biomarker data
                            load_biomarkers = function() {
                              # Check for biomarker files in multiple possible locations
                              possible_locations <- c(
                                self$config$processed_dir,
                                file.path(getwd()),  # Current directory
                                file.path(getwd(), "R"),  # R subdirectory
                                file.path(dirname(getwd()), "luad_real_only_data", "processed")  # Original location
                              )

                              biomarker_files <- list(
                                egfr = "egfr_mutation_status.csv",
                                alk = "alk_fusion_status.csv",
                                pdl1 = "cd274_pdl1_expression.csv",
                                risk_scores = "lasso_risk_scores.csv"
                              )

                              biomarkers <- list()

                              for (marker in names(biomarker_files)) {
                                file_found <- FALSE

                                # Search for the file in possible locations
                                for (location in possible_locations) {
                                  file_path <- file.path(location, biomarker_files[[marker]])

                                  if (file.exists(file_path)) {
                                    tryCatch({
                                      biomarkers[[marker]] <- read_csv(file_path, show_col_types = FALSE) %>%
                                        mutate(patient_id = self$standardize_tcga_id(sample)) %>%
                                        distinct(patient_id, .keep_all = TRUE)
                                      cat("Loaded", marker, "data:", nrow(biomarkers[[marker]]), "patients from", file_path, "\n")
                                      file_found <- TRUE
                                      break
                                    }, error = function(e) {
                                      cat("Error loading", marker, "from", file_path, ":", e$message, "\n")
                                    })
                                  }
                                }

                                if (!file_found) {
                                  cat("Warning: Biomarker file not found:", biomarker_files[[marker]], "\n")
                                  cat("Searched in:", paste(possible_locations, collapse = ", "), "\n")
                                  # Create empty placeholder
                                  biomarkers[[marker]] <- data.frame(
                                    patient_id = character(0),
                                    stringsAsFactors = FALSE
                                  )
                                }
                              }

                              return(biomarkers)
                            },

                            # Standardize TCGA IDs
                            standardize_tcga_id = function(ids) {
                              ids %>%
                                gsub("\\.", "-", .) %>%
                                substr(1, 12) %>%
                                toupper()
                            },

                            # Standardize cancer stages
                            standardize_stage = function(stages) {
                              case_when(
                                grepl("Stage I[^V]|^I[^V]", stages, ignore.case = TRUE) ~ "Stage I",
                                grepl("Stage IA|^IA", stages, ignore.case = TRUE) ~ "Stage IA",
                                grepl("Stage IB|^IB", stages, ignore.case = TRUE) ~ "Stage IB",
                                grepl("Stage II[^I]|^II[^I]", stages, ignore.case = TRUE) ~ "Stage II",
                                grepl("Stage IIA|^IIA", stages, ignore.case = TRUE) ~ "Stage IIA",
                                grepl("Stage IIB|^IIB", stages, ignore.case = TRUE) ~ "Stage IIB",
                                grepl("Stage III[^AB]|^III[^AB]", stages, ignore.case = TRUE) ~ "Stage III",
                                grepl("Stage IIIA|^IIIA", stages, ignore.case = TRUE) ~ "Stage IIIA",
                                grepl("Stage IIIB|^IIIB", stages, ignore.case = TRUE) ~ "Stage IIIB",
                                grepl("Stage IV|^IV", stages, ignore.case = TRUE) ~ "Stage IV",
                                TRUE ~ as.character(stages)
                              )
                            },

                            # Standardize smoking status
                            standardize_smoking = function(smoking_codes) {
                              case_when(
                                smoking_codes == "1" ~ "Never",
                                smoking_codes == "2" ~ "Current",
                                smoking_codes == "3" ~ "Former",
                                smoking_codes == "4" ~ "Former", # Occasional -> Former
                                TRUE ~ NA_character_
                              )
                            }
                          )
)

# =============================================================================
# 3. TREATMENT RECOMMENDATION ENGINE
# =============================================================================

TreatmentEngine <- R6::R6Class("TreatmentEngine",
                               public = list(
                                 config = NULL,

                                 initialize = function(config) {
                                   self$config <- config
                                 },

                                 # Main treatment recommendation function
                                 recommend_treatments = function(patient_data) {
                                   cat("Generating treatment recommendations for", nrow(patient_data), "patients...\n")

                                   # Initialize recommendation data frame
                                   recommendations <- patient_data %>%
                                     mutate(
                                       primary_treatment = "Standard Chemotherapy",
                                       surgical_option = "Not candidate",
                                       radiation_option = "Not indicated",
                                       second_line_treatment = "TBD",
                                       treatment_intent = "Palliative",
                                       combination_therapy = FALSE,
                                       biomarker_guided = FALSE,
                                       confidence_score = 0.6,
                                       treatment_rationale = "Standard of care",
                                       evidence_level = "2A"
                                     )

                                   # Apply treatment algorithms
                                   recommendations <- recommendations %>%
                                     self$apply_surgical_algorithms(.) %>%
                                     self$apply_radiation_algorithms(.) %>%
                                     self$apply_targeted_therapy(.) %>%
                                     self$apply_immunotherapy(.) %>%
                                     self$apply_chemotherapy(.) %>%
                                     self$apply_palliative_care(.) %>%
                                     self$assign_second_line_treatments(.)

                                   return(recommendations)
                                 },

                                 # Surgical treatment algorithms
                                 apply_surgical_algorithms = function(data) {
                                   data %>%
                                     mutate(
                                       surgical_candidate = case_when(
                                         # Early stage with good performance status
                                         stage %in% c("Stage I", "Stage IA", "Stage IB") & ecog_ps <= 1 ~ "Primary",
                                         stage %in% c("Stage II", "Stage IIA", "Stage IIB") & ecog_ps <= 1 ~ "Primary",

                                         # Resectable Stage III
                                         stage %in% c("Stage IIIA") & ecog_ps <= 1 & age < 75 ~ "Post-neoadjuvant",

                                         TRUE ~ "Not candidate"
                                       ),

                                       # Update treatment for surgical candidates
                                       primary_treatment = case_when(
                                         surgical_candidate == "Primary" & age < 75 ~ "Surgery (Lobectomy)",
                                         surgical_candidate == "Primary" & age >= 75 ~ "Surgery (Sublobar)",
                                         surgical_candidate == "Post-neoadjuvant" ~ "Neoadjuvant Chemotherapy → Surgery",
                                         TRUE ~ primary_treatment
                                       ),

                                       surgical_option = case_when(
                                         surgical_candidate == "Primary" & age < 75 ~ "Lobectomy/VATS",
                                         surgical_candidate == "Primary" & age >= 75 ~ "Sublobar resection",
                                         surgical_candidate == "Post-neoadjuvant" ~ "Surgery after neoadjuvant",
                                         TRUE ~ surgical_option
                                       ),

                                       treatment_intent = case_when(
                                         surgical_candidate != "Not candidate" ~ "Curative",
                                         TRUE ~ treatment_intent
                                       ),

                                       confidence_score = case_when(
                                         surgical_candidate == "Primary" ~ 0.9,
                                         surgical_candidate == "Post-neoadjuvant" ~ 0.8,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         surgical_candidate == "Primary" ~ paste0("Early stage (", stage, "), ECOG ≤1 - surgical resection"),
                                         surgical_candidate == "Post-neoadjuvant" ~ "Resectable Stage IIIA - neoadjuvant approach",
                                         TRUE ~ treatment_rationale
                                       ),

                                       evidence_level = case_when(
                                         surgical_candidate != "Not candidate" ~ "1",
                                         TRUE ~ evidence_level
                                       )
                                     )
                                 },

                                 # Radiation therapy algorithms
                                 apply_radiation_algorithms = function(data) {
                                   data %>%
                                     mutate(
                                       radiation_candidate = case_when(
                                         # SBRT for medically inoperable early stage
                                         stage %in% c("Stage I", "Stage IA", "Stage IB") & surgical_candidate == "Not candidate" & ecog_ps <= 2 ~ "SBRT",

                                         # Concurrent chemoradiation for unresectable Stage III
                                         stage %in% c("Stage III", "Stage IIIA", "Stage IIIB") & ecog_ps <= 2 ~ "Concurrent_ChemoRT",

                                         # Palliative radiation
                                         stage == "Stage IV" & ecog_ps <= 3 ~ "Palliative",

                                         TRUE ~ "Not indicated"
                                       ),

                                       # Update treatments for radiation candidates
                                       primary_treatment = case_when(
                                         radiation_candidate == "SBRT" ~ "Stereotactic Body Radiation (SBRT)",
                                         radiation_candidate == "Concurrent_ChemoRT" ~ "Concurrent Chemoradiation → Durvalumab",
                                         TRUE ~ primary_treatment
                                       ),

                                       radiation_option = case_when(
                                         radiation_candidate == "SBRT" ~ "SBRT (curative)",
                                         radiation_candidate == "Concurrent_ChemoRT" ~ "Concurrent chemoradiation",
                                         radiation_candidate == "Palliative" ~ "Palliative radiation",
                                         TRUE ~ radiation_option
                                       ),

                                       treatment_intent = case_when(
                                         radiation_candidate %in% c("SBRT", "Concurrent_ChemoRT") ~ "Curative",
                                         TRUE ~ treatment_intent
                                       ),

                                       confidence_score = case_when(
                                         radiation_candidate == "SBRT" ~ 0.85,
                                         radiation_candidate == "Concurrent_ChemoRT" ~ 0.8,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         radiation_candidate == "SBRT" ~ "Medically inoperable early stage - SBRT",
                                         radiation_candidate == "Concurrent_ChemoRT" ~ "Unresectable Stage III - definitive chemoRT",
                                         TRUE ~ treatment_rationale
                                       )
                                     )
                                 },

                                 # Targeted therapy algorithms
                                 apply_targeted_therapy = function(data) {
                                   data %>%
                                     mutate(
                                       # EGFR targeted therapy (highest priority)
                                       primary_treatment = case_when(
                                         !is.na(egfr_mutation) & egfr_mutation == "Positive" ~ "Osimertinib (EGFR TKI)",
                                         TRUE ~ primary_treatment
                                       ),

                                       # ALK targeted therapy (second priority)
                                       primary_treatment = case_when(
                                         !is.na(alk_fusion) & alk_fusion == "Positive" & !grepl("EGFR", primary_treatment) ~ "Alectinib (ALK inhibitor)",
                                         TRUE ~ primary_treatment
                                       ),

                                       # Update other fields for targeted therapy
                                       biomarker_guided = case_when(
                                         grepl("Osimertinib|Alectinib", primary_treatment) ~ TRUE,
                                         TRUE ~ biomarker_guided
                                       ),

                                       confidence_score = case_when(
                                         grepl("Osimertinib", primary_treatment) ~ 0.95,
                                         grepl("Alectinib", primary_treatment) ~ 0.90,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         grepl("Osimertinib", primary_treatment) ~ "EGFR mutation positive - osimertinib first-line",
                                         grepl("Alectinib", primary_treatment) ~ "ALK fusion positive - alectinib first-line",
                                         TRUE ~ treatment_rationale
                                       ),

                                       evidence_level = case_when(
                                         grepl("Osimertinib|Alectinib", primary_treatment) ~ "1",
                                         TRUE ~ evidence_level
                                       )
                                     )
                                 },

                                 # Immunotherapy algorithms
                                 apply_immunotherapy = function(data) {
                                   data %>%
                                     mutate(
                                       # PD-L1 high (≥50%) - monotherapy
                                       primary_treatment = case_when(
                                         !grepl("Osimertinib|Alectinib|Surgery|SBRT|Chemoradiation", primary_treatment) &
                                           !is.na(pdl1_expression) & pdl1_expression == "High" &
                                           stage == "Stage IV" & ecog_ps <= 1 ~ "Pembrolizumab monotherapy",
                                         TRUE ~ primary_treatment
                                       ),

                                       # PD-L1 moderate (1-49%) - combination
                                       primary_treatment = case_when(
                                         !grepl("Osimertinib|Alectinib|Surgery|SBRT|Chemoradiation|Pembrolizumab monotherapy", primary_treatment) &
                                           !is.na(pdl1_expression) & pdl1_expression == "Medium" &
                                           stage == "Stage IV" & ecog_ps <= 2 ~ "Pembrolizumab + Chemotherapy",
                                         TRUE ~ primary_treatment
                                       ),

                                       # Update combination therapy flag
                                       combination_therapy = case_when(
                                         grepl("Pembrolizumab \\+ Chemotherapy", primary_treatment) ~ TRUE,
                                         TRUE ~ combination_therapy
                                       ),

                                       biomarker_guided = case_when(
                                         grepl("Pembrolizumab", primary_treatment) ~ TRUE,
                                         TRUE ~ biomarker_guided
                                       ),

                                       confidence_score = case_when(
                                         grepl("Pembrolizumab monotherapy", primary_treatment) ~ 0.85,
                                         grepl("Pembrolizumab \\+ Chemotherapy", primary_treatment) ~ 0.8,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         grepl("Pembrolizumab monotherapy", primary_treatment) ~ "PD-L1 ≥50%, good ECOG - pembrolizumab monotherapy",
                                         grepl("Pembrolizumab \\+ Chemotherapy", primary_treatment) ~ "PD-L1 1-49% - pembrolizumab + platinum doublet",
                                         TRUE ~ treatment_rationale
                                       )
                                     )
                                 },

                                 # Chemotherapy algorithms
                                 apply_chemotherapy = function(data) {
                                   data %>%
                                     mutate(
                                       # Standard chemotherapy for patients without other specific treatments
                                       primary_treatment = case_when(
                                         primary_treatment == "Standard Chemotherapy" & stage == "Stage IV" & ecog_ps <= 2 ~ "Carboplatin + Paclitaxel",
                                         primary_treatment == "Standard Chemotherapy" & stage == "Stage IV" & ecog_ps >= 3 ~ "Single-agent chemotherapy",
                                         TRUE ~ primary_treatment
                                       ),

                                       confidence_score = case_when(
                                         grepl("Carboplatin", primary_treatment) ~ 0.7,
                                         grepl("Single-agent", primary_treatment) ~ 0.6,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         grepl("Carboplatin", primary_treatment) ~ "Good ECOG, no driver mutations - platinum doublet",
                                         grepl("Single-agent", primary_treatment) ~ "Poor ECOG - single-agent palliative chemotherapy",
                                         TRUE ~ treatment_rationale
                                       )
                                     )
                                 },

                                 # Palliative care algorithms
                                 apply_palliative_care = function(data) {
                                   data %>%
                                     mutate(
                                       primary_treatment = case_when(
                                         ecog_ps >= 4 | (age >= 80 & stage == "Stage IV") ~ "Best Supportive Care",
                                         TRUE ~ primary_treatment
                                       ),

                                       treatment_intent = case_when(
                                         grepl("Best Supportive Care", primary_treatment) ~ "Supportive",
                                         stage == "Stage IV" ~ "Palliative",
                                         TRUE ~ treatment_intent
                                       ),

                                       confidence_score = case_when(
                                         grepl("Best Supportive Care", primary_treatment) ~ 0.9,
                                         TRUE ~ confidence_score
                                       ),

                                       treatment_rationale = case_when(
                                         grepl("Best Supportive Care", primary_treatment) ~ "ECOG ≥4 or elderly Stage IV - supportive care",
                                         TRUE ~ treatment_rationale
                                       )
                                     )
                                 },

                                 # Assign second-line treatments
                                 assign_second_line_treatments = function(data) {
                                   data %>%
                                     mutate(
                                       second_line_treatment = case_when(
                                         grepl("Osimertinib", primary_treatment) ~ "Chemotherapy + Bevacizumab (post-EGFR resistance)",
                                         grepl("Alectinib", primary_treatment) ~ "Next-generation ALK inhibitor or chemotherapy",
                                         grepl("Pembrolizumab", primary_treatment) ~ "Docetaxel + Ramucirumab",
                                         grepl("Chemotherapy", primary_treatment) & !grepl("Pembrolizumab", primary_treatment) ~ "Immunotherapy (if eligible)",
                                         grepl("Surgery", primary_treatment) ~ "Adjuvant chemotherapy (if high-risk)",
                                         TRUE ~ "Clinical trial or supportive care"
                                       )
                                     )
                                 }
                               )
)

# =============================================================================
# 4. ANALYSIS AND VALIDATION MODULE
# =============================================================================

AnalysisEngine <- R6::R6Class("AnalysisEngine",
                              public = list(
                                config = NULL,

                                initialize = function(config) {
                                  self$config <- config
                                },

                                # Comprehensive analysis of treatment recommendations
                                analyze_recommendations = function(recommendations) {
                                  cat("Performing comprehensive analysis...\n")

                                  analysis_results <- list(
                                    summary = self$generate_summary_statistics(recommendations),
                                    treatment_distribution = self$analyze_treatment_distribution(recommendations),
                                    biomarker_impact = self$analyze_biomarker_impact(recommendations),
                                    clinical_factors = self$analyze_clinical_factors(recommendations),
                                    stage_analysis = self$analyze_stage_patterns(recommendations),
                                    quality_metrics = self$calculate_quality_metrics(recommendations)
                                  )

                                  return(analysis_results)
                                },

                                # Generate summary statistics
                                generate_summary_statistics = function(data) {
                                  list(
                                    total_patients = nrow(data),
                                    patients_with_ecog = sum(!is.na(data$ecog_ps)),
                                    patients_with_biomarkers = sum(!is.na(data$egfr_mutation) | !is.na(data$alk_fusion) | !is.na(data$pdl1_expression)),
                                    biomarker_guided_rate = round(mean(data$biomarker_guided) * 100, 1),
                                    curative_intent_rate = round(mean(data$treatment_intent == "Curative") * 100, 1),
                                    combination_therapy_rate = round(mean(data$combination_therapy) * 100, 1),
                                    average_confidence = round(mean(data$confidence_score), 2)
                                  )
                                },

                                # Analyze treatment distribution
                                analyze_treatment_distribution = function(data) {
                                  data %>%
                                    count(primary_treatment, sort = TRUE) %>%
                                    mutate(
                                      percentage = round(n / sum(n) * 100, 1),
                                      cumulative_percentage = round(cumsum(n) / sum(n) * 100, 1)
                                    )
                                },

                                # Analyze biomarker impact
                                analyze_biomarker_impact = function(data) {
                                  biomarker_analysis <- list()

                                  # EGFR analysis
                                  if (sum(!is.na(data$egfr_mutation)) > 0) {
                                    biomarker_analysis$egfr <- data %>%
                                      filter(!is.na(egfr_mutation)) %>%
                                      group_by(egfr_mutation) %>%
                                      summarise(
                                        n_patients = n(),
                                        targeted_therapy_rate = round(mean(grepl("Osimertinib", primary_treatment)) * 100, 1),
                                        avg_confidence = round(mean(confidence_score), 2),
                                        .groups = "drop"
                                      )
                                  }

                                  # ALK analysis
                                  if (sum(!is.na(data$alk_fusion)) > 0) {
                                    biomarker_analysis$alk <- data %>%
                                      filter(!is.na(alk_fusion)) %>%
                                      group_by(alk_fusion) %>%
                                      summarise(
                                        n_patients = n(),
                                        targeted_therapy_rate = round(mean(grepl("Alectinib", primary_treatment)) * 100, 1),
                                        avg_confidence = round(mean(confidence_score), 2),
                                        .groups = "drop"
                                      )
                                  }

                                  # PD-L1 analysis
                                  if (sum(!is.na(data$pdl1_expression)) > 0) {
                                    biomarker_analysis$pdl1 <- data %>%
                                      filter(!is.na(pdl1_expression)) %>%
                                      group_by(pdl1_expression) %>%
                                      summarise(
                                        n_patients = n(),
                                        immunotherapy_rate = round(mean(grepl("Pembrolizumab", primary_treatment)) * 100, 1),
                                        avg_confidence = round(mean(confidence_score), 2),
                                        .groups = "drop"
                                      )
                                  }

                                  return(biomarker_analysis)
                                },

                                # Analyze clinical factors
                                analyze_clinical_factors = function(data) {
                                  list(
                                    ecog_impact = data %>%
                                      filter(!is.na(ecog_ps)) %>%
                                      group_by(ecog_ps) %>%
                                      summarise(
                                        n_patients = n(),
                                        curative_rate = round(mean(treatment_intent == "Curative") * 100, 1),
                                        surgical_rate = round(mean(grepl("Surgery", primary_treatment)) * 100, 1),
                                        .groups = "drop"
                                      ),

                                    age_impact = data %>%
                                      filter(!is.na(age)) %>%
                                      group_by(age_group) %>%
                                      summarise(
                                        n_patients = n(),
                                        surgical_rate = round(mean(grepl("Surgery", primary_treatment)) * 100, 1),
                                        targeted_rate = round(mean(grepl("Osimertinib|Alectinib", primary_treatment)) * 100, 1),
                                        .groups = "drop"
                                      )
                                  )
                                },

                                # Analyze stage-specific patterns
                                analyze_stage_patterns = function(data) {
                                  data %>%
                                    filter(!is.na(stage)) %>%
                                    group_by(stage) %>%
                                    summarise(
                                      n_patients = n(),
                                      surgical_rate = round(mean(grepl("Surgery", primary_treatment)) * 100, 1),
                                      radiation_rate = round(mean(grepl("Radiation|SBRT", primary_treatment)) * 100, 1),
                                      systemic_rate = round(mean(grepl("Osimertinib|Alectinib|Pembrolizumab|Chemotherapy", primary_treatment)) * 100, 1),
                                      curative_intent_rate = round(mean(treatment_intent == "Curative") * 100, 1),
                                      .groups = "drop"
                                    )
                                },

                                # Calculate quality metrics
                                calculate_quality_metrics = function(data) {
                                  list(
                                    data_completeness = list(
                                      age = round(mean(!is.na(data$age)) * 100, 1),
                                      ecog = round(mean(!is.na(data$ecog_ps)) * 100, 1),
                                      stage = round(mean(!is.na(data$stage)) * 100, 1),
                                      biomarkers = round(mean(!is.na(data$egfr_mutation) | !is.na(data$alk_fusion) | !is.na(data$pdl1_expression)) * 100, 1)
                                    ),

                                    clinical_appropriateness = list(
                                      evidence_level_1_rate = round(mean(data$evidence_level == "1") * 100, 1),
                                      high_confidence_rate = round(mean(data$confidence_score >= 0.8) * 100, 1),
                                      guideline_concordance = round(mean(data$biomarker_guided | data$evidence_level == "1") * 100, 1)
                                    )
                                  )
                                }
                              )
)

# =============================================================================
# 5. VISUALIZATION MODULE
# =============================================================================

VisualizationEngine <- R6::R6Class("VisualizationEngine",
                                   public = list(
                                     config = NULL,

                                     initialize = function(config) {
                                       self$config <- config
                                     },

                                     # Create all visualizations
                                     create_all_plots = function(data, analysis_results) {
                                       plots <- list(
                                         treatment_distribution = self$plot_treatment_distribution(data),
                                         stage_treatment_heatmap = self$plot_stage_treatment_heatmap(data),
                                         ecog_impact = self$plot_ecog_impact(data),
                                         biomarker_outcomes = self$plot_biomarker_outcomes(data),
                                         confidence_distribution = self$plot_confidence_distribution(data)
                                       )

                                       # Save all plots
                                       self$save_plots(plots)

                                       return(plots)
                                     },

                                     # Treatment distribution plot
                                     plot_treatment_distribution = function(data) {
                                       data %>%
                                         count(primary_treatment, sort = TRUE) %>%
                                         mutate(percentage = n / sum(n) * 100) %>%
                                         ggplot(aes(x = reorder(primary_treatment, n), y = n)) +
                                         geom_col(fill = "steelblue", alpha = 0.8) +
                                         geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")),
                                                   hjust = -0.1, size = 3) +
                                         coord_flip() +
                                         theme_minimal() +
                                         labs(
                                           title = "Distribution of Primary Treatment Recommendations",
                                           subtitle = paste("N =", nrow(data), "patients"),
                                           x = "Primary Treatment",
                                           y = "Number of Patients"
                                         ) +
                                         theme(
                                           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                           plot.subtitle = element_text(hjust = 0.5, size = 12),
                                           axis.text.y = element_text(size = 10)
                                         )
                                     },

                                     # Stage-treatment heatmap
                                     plot_stage_treatment_heatmap = function(data) {
                                       stage_treatment_data <- data %>%
                                         filter(!is.na(stage)) %>%
                                         count(stage, primary_treatment) %>%
                                         group_by(stage) %>%
                                         mutate(
                                           total_stage = sum(n),
                                           percentage = n / total_stage * 100
                                         ) %>%
                                         ungroup()

                                       ggplot(stage_treatment_data, aes(x = stage, y = primary_treatment, fill = percentage)) +
                                         geom_tile(color = "white", size = 0.5) +
                                         geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
                                                   size = 3, color = "white") +
                                         scale_fill_gradient(low = "lightblue", high = "darkblue", name = "% of Stage") +
                                         theme_minimal() +
                                         labs(
                                           title = "Treatment Distribution by Cancer Stage",
                                           x = "Cancer Stage",
                                           y = "Primary Treatment"
                                         ) +
                                         theme(
                                           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                           axis.text.x = element_text(angle = 45, hjust = 1),
                                           axis.text.y = element_text(size = 9)
                                         )
                                     },

                                     # ECOG performance status impact
                                     plot_ecog_impact = function(data) {
                                       if (sum(!is.na(data$ecog_ps)) == 0) {
                                         return(ggplot() + theme_void() + labs(title = "ECOG data not available"))
                                       }

                                       ecog_data <- data %>%
                                         filter(!is.na(ecog_ps)) %>%
                                         mutate(
                                           treatment_category = case_when(
                                             grepl("Surgery", primary_treatment) ~ "Surgery",
                                             grepl("Osimertinib|Alectinib", primary_treatment) ~ "Targeted Therapy",
                                             grepl("Pembrolizumab", primary_treatment) ~ "Immunotherapy",
                                             grepl("Chemotherapy|Carboplatin", primary_treatment) ~ "Chemotherapy",
                                             grepl("Radiation|SBRT", primary_treatment) ~ "Radiation",
                                             TRUE ~ "Supportive Care"
                                           )
                                         ) %>%
                                         count(ecog_ps, treatment_category) %>%
                                         group_by(ecog_ps) %>%
                                         mutate(percentage = n / sum(n) * 100)

                                       ggplot(ecog_data, aes(x = factor(ecog_ps), y = percentage, fill = treatment_category)) +
                                         geom_col(position = "stack", alpha = 0.8) +
                                         scale_fill_brewer(type = "qual", palette = "Set2", name = "Treatment Category") +
                                         theme_minimal() +
                                         labs(
                                           title = "Treatment Selection by ECOG Performance Status",
                                           subtitle = paste("N =", sum(!is.na(data$ecog_ps)), "patients with ECOG scores"),
                                           x = "ECOG Performance Status",
                                           y = "Percentage of Patients"
                                         ) +
                                         theme(
                                           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                           plot.subtitle = element_text(hjust = 0.5, size = 12)
                                         )
                                     },

                                     # Biomarker outcomes visualization
                                     plot_biomarker_outcomes = function(data) {
                                       biomarker_data <- data %>%
                                         mutate(
                                           biomarker_status = case_when(
                                             !is.na(egfr_mutation) & egfr_mutation == "Positive" ~ "EGFR+",
                                             !is.na(alk_fusion) & alk_fusion == "Positive" ~ "ALK+",
                                             !is.na(pdl1_expression) & pdl1_expression == "High" ~ "PD-L1 High",
                                             !is.na(pdl1_expression) & pdl1_expression == "Medium" ~ "PD-L1 Medium",
                                             TRUE ~ "No actionable biomarker"
                                           ),
                                           targeted_treatment = case_when(
                                             grepl("Osimertinib", primary_treatment) ~ "EGFR TKI",
                                             grepl("Alectinib", primary_treatment) ~ "ALK inhibitor",
                                             grepl("Pembrolizumab", primary_treatment) ~ "Immunotherapy",
                                             TRUE ~ "Other"
                                           )
                                         ) %>%
                                         count(biomarker_status, targeted_treatment) %>%
                                         group_by(biomarker_status) %>%
                                         mutate(percentage = n / sum(n) * 100)

                                       ggplot(biomarker_data, aes(x = biomarker_status, y = percentage, fill = targeted_treatment)) +
                                         geom_col(position = "stack", alpha = 0.8) +
                                         scale_fill_brewer(type = "qual", palette = "Set1", name = "Treatment Type") +
                                         theme_minimal() +
                                         labs(
                                           title = "Biomarker-Treatment Matching",
                                           subtitle = "Precision medicine approach validation",
                                           x = "Biomarker Status",
                                           y = "Percentage of Patients"
                                         ) +
                                         theme(
                                           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                           plot.subtitle = element_text(hjust = 0.5, size = 12),
                                           axis.text.x = element_text(angle = 45, hjust = 1)
                                         )
                                     },

                                     # Confidence score distribution
                                     plot_confidence_distribution = function(data) {
                                       ggplot(data, aes(x = confidence_score)) +
                                         geom_histogram(bins = 20, fill = "lightcoral", alpha = 0.7, color = "white") +
                                         geom_vline(xintercept = mean(data$confidence_score),
                                                    linetype = "dashed", color = "red", size = 1) +
                                         annotate("text", x = mean(data$confidence_score) + 0.05,
                                                  y = Inf, vjust = 2,
                                                  label = paste("Mean =", round(mean(data$confidence_score), 2)),
                                                  color = "red") +
                                         theme_minimal() +
                                         labs(
                                           title = "Distribution of Treatment Confidence Scores",
                                           subtitle = "Higher scores indicate greater evidence-based confidence",
                                           x = "Confidence Score",
                                           y = "Number of Patients"
                                         ) +
                                         theme(
                                           plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                                           plot.subtitle = element_text(hjust = 0.5, size = 12)
                                         )
                                     },

                                     # Save all plots
                                     save_plots = function(plots) {
                                       for (plot_name in names(plots)) {
                                         filename <- file.path(self$config$plots_dir, paste0(plot_name, ".pdf"))
                                         ggsave(filename, plots[[plot_name]], width = 12, height = 8, dpi = 300)
                                         cat("Saved plot:", filename, "\n")
                                       }
                                     }
                                   )
)

# =============================================================================
# 6. REPORTING MODULE
# =============================================================================

ReportGenerator <- R6::R6Class("ReportGenerator",
                               public = list(
                                 config = NULL,

                                 initialize = function(config) {
                                   self$config <- config
                                 },

                                 # Generate comprehensive report
                                 generate_comprehensive_report = function(data, analysis_results, plots) {
                                   report_file <- file.path(self$config$results_dir, "reports", "comprehensive_treatment_report.html")

                                   # Create report content
                                   report_content <- self$create_report_content(data, analysis_results)

                                   # Write HTML report
                                   self$write_html_report(report_content, report_file)

                                   # Generate executive summary
                                   self$generate_executive_summary(analysis_results)

                                   # Generate patient treatment cards
                                   self$generate_patient_cards(data)

                                   cat("Comprehensive report generated:", report_file, "\n")
                                 },

                                 # Create report content
                                 create_report_content = function(data, analysis) {
                                   list(
                                     title = "LUAD Treatment Recommendation System - Comprehensive Report",
                                     date = Sys.Date(),
                                     summary = analysis$summary,
                                     treatment_distribution = analysis$treatment_distribution,
                                     biomarker_impact = analysis$biomarker_impact,
                                     clinical_factors = analysis$clinical_factors,
                                     stage_analysis = analysis$stage_analysis,
                                     quality_metrics = analysis$quality_metrics,
                                     recommendations = self$generate_clinical_recommendations(analysis)
                                   )
                                 },

                                 # Write HTML report
                                 write_html_report = function(content, filename) {
                                   # HTML template
                                   html_content <- paste0(
                                     '<!DOCTYPE html>
        <html>
        <head>
          <title>', content$title, '</title>
          <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .header { background-color: #f8f9fa; padding: 20px; border-radius: 5px; }
            .section { margin: 20px 0; }
            .metric { background-color: #e9ecef; padding: 15px; margin: 10px 0; border-radius: 5px; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .highlight { background-color: #fff3cd; padding: 10px; border-radius: 5px; }
          </style>
        </head>
        <body>
          <div class="header">
            <h1>', content$title, '</h1>
            <p>Generated on: ', content$date, '</p>
          </div>

          <div class="section">
            <h2>Executive Summary</h2>
            <div class="metric">
              <strong>Total Patients Analyzed:</strong> ', content$summary$total_patients, '<br>
              <strong>Biomarker-Guided Treatments:</strong> ', content$summary$biomarker_guided_rate, '%<br>
              <strong>Curative Intent Rate:</strong> ', content$summary$curative_intent_rate, '%<br>
              <strong>Average Confidence Score:</strong> ', content$summary$average_confidence, '
            </div>
          </div>

          <div class="section">
            <h2>Treatment Distribution</h2>
            ', self$create_table_html(content$treatment_distribution), '
          </div>

          <div class="section">
            <h2>Clinical Recommendations</h2>
            <div class="highlight">
              ', paste(content$recommendations, collapse = "<br>"), '
            </div>
          </div>

        </body>
        </html>'
                                   )

                                   writeLines(html_content, filename)
                                 },

                                 # Generate executive summary
                                 generate_executive_summary = function(analysis) {
                                   summary_file <- file.path(self$config$results_dir, "reports", "executive_summary.txt")

                                   summary_content <- paste0(
                                     "LUAD TREATMENT RECOMMENDATION SYSTEM - EXECUTIVE SUMMARY\n",
                                     "=======================================================\n\n",
                                     "Analysis Date: ", Sys.Date(), "\n",
                                     "Total Patients: ", analysis$summary$total_patients, "\n\n",
                                     "KEY FINDINGS:\n",
                                     "- Biomarker-guided treatments: ", analysis$summary$biomarker_guided_rate, "%\n",
                                     "- Curative intent achievable: ", analysis$summary$curative_intent_rate, "%\n",
                                     "- Average treatment confidence: ", analysis$summary$average_confidence, "\n\n",
                                     "CLINICAL IMPACT:\n",
                                     "- Precision medicine approach implemented\n",
                                     "- Evidence-based treatment selection\n",
                                     "- Stage-appropriate therapeutic strategies\n",
                                     "- Performance status integration\n\n",
                                     "RECOMMENDATIONS:\n",
                                     "1. Implement comprehensive biomarker testing\n",
                                     "2. Consider multidisciplinary tumor board reviews\n",
                                     "3. Monitor treatment response and adapt accordingly\n",
                                     "4. Evaluate for clinical trial opportunities\n"
                                   )

                                   writeLines(summary_content, summary_file)
                                   cat("Executive summary saved:", summary_file, "\n")
                                 },

                                 # Generate patient treatment cards
                                 generate_patient_cards = function(data, n_cards = 20) {
                                   cards_file <- file.path(self$config$results_dir, "reports", "patient_treatment_cards.txt")

                                   # Select diverse patients
                                   selected_patients <- data %>%
                                     group_by(primary_treatment) %>%
                                     slice_head(n = 2) %>%
                                     ungroup() %>%
                                     head(n_cards)

                                   cards_content <- c()

                                   for (i in 1:nrow(selected_patients)) {
                                     patient <- selected_patients[i, ]

                                     card <- paste0(
                                       "\n", paste(rep("=", 70), collapse = ""), "\n",
                                       "PATIENT TREATMENT CARD #", i, "\n",
                                       paste(rep("=", 70), collapse = ""), "\n",
                                       "Patient ID: ", patient$patient_id, "\n",
                                       "Age: ", ifelse(is.na(patient$age), "Unknown", paste(patient$age, "years")),
                                       " | Gender: ", ifelse(is.na(patient$gender), "Unknown", patient$gender), "\n",
                                       "Stage: ", ifelse(is.na(patient$stage), "Unknown", patient$stage), "\n",
                                       "ECOG PS: ", ifelse(is.na(patient$ecog_ps), "Unknown", patient$ecog_ps), "\n\n",
                                       "BIOMARKERS:\n",
                                       "- EGFR: ", ifelse(is.na(patient$egfr_mutation), "Unknown", patient$egfr_mutation), "\n",
                                       "- ALK: ", ifelse(is.na(patient$alk_fusion), "Unknown", patient$alk_fusion), "\n",
                                       "- PD-L1: ", ifelse(is.na(patient$pdl1_expression), "Unknown", patient$pdl1_expression), "\n\n",
                                       "TREATMENT RECOMMENDATIONS:\n",
                                       "Primary: ", patient$primary_treatment, "\n",
                                       "Second-line: ", patient$second_line_treatment, "\n",
                                       "Evidence Level: ", patient$evidence_level, "\n",
                                       "Confidence: ", round(patient$confidence_score, 2), "\n\n",
                                       "Rationale: ", patient$treatment_rationale, "\n",
                                       paste(rep("=", 70), collapse = ""), "\n"
                                     )

                                     cards_content <- c(cards_content, card)
                                   }

                                   writeLines(cards_content, cards_file)
                                   cat("Patient treatment cards saved:", cards_file, "\n")
                                 },

                                 # Helper function to create HTML tables
                                 create_table_html = function(data) {
                                   if (is.null(data) || nrow(data) == 0) return("")

                                   table_html <- "<table>"
                                   table_html <- paste0(table_html, "<tr>")

                                   # Headers
                                   for (col in colnames(data)) {
                                     table_html <- paste0(table_html, "<th>", col, "</th>")
                                   }
                                   table_html <- paste0(table_html, "</tr>")

                                   # Rows
                                   for (i in 1:nrow(data)) {
                                     table_html <- paste0(table_html, "<tr>")
                                     for (col in colnames(data)) {
                                       table_html <- paste0(table_html, "<td>", data[i, col], "</td>")
                                     }
                                     table_html <- paste0(table_html, "</tr>")
                                   }

                                   table_html <- paste0(table_html, "</table>")
                                   return(table_html)
                                 },

                                 # Generate clinical recommendations
                                 generate_clinical_recommendations = function(analysis) {
                                   recommendations <- c(
                                     "1. Implement comprehensive molecular profiling for all LUAD patients",
                                     "2. Establish multidisciplinary tumor boards for complex cases",
                                     "3. Consider clinical trial enrollment for novel therapies",
                                     "4. Regular monitoring and adaptation of treatment plans",
                                     "5. Ensure adequate supportive care integration"
                                   )

                                   # Add specific recommendations based on analysis
                                   if (analysis$summary$biomarker_guided_rate < 50) {
                                     recommendations <- c(recommendations,
                                                          "6. PRIORITY: Increase biomarker testing rates - currently suboptimal")
                                   }

                                   if (analysis$summary$curative_intent_rate < 30) {
                                     recommendations <- c(recommendations,
                                                          "7. Focus on early detection and surgical evaluation")
                                   }

                                   return(recommendations)
                                 }
                               )
)

# =============================================================================
# 7. MAIN ORCHESTRATOR CLASS
# =============================================================================

LUADTreatmentSystem <- R6::R6Class("LUADTreatmentSystem",
                                   public = list(
                                     config = NULL,
                                     data_loader = NULL,
                                     treatment_engine = NULL,
                                     analysis_engine = NULL,
                                     visualization_engine = NULL,
                                     report_generator = NULL,

                                     initialize = function(config = project_config) {
                                       self$config <- config
                                       self$data_loader <- DataLoader$new(config)
                                       self$treatment_engine <- TreatmentEngine$new(config)
                                       self$analysis_engine <- AnalysisEngine$new(config)
                                       self$visualization_engine <- VisualizationEngine$new(config)
                                       self$report_generator <- ReportGenerator$new(config)

                                       # Create project structure
                                       create_project_structure(config)

                                       cat("LUAD Treatment System initialized successfully\n")
                                     },

                                     # Run complete analysis pipeline
                                     run_complete_analysis = function(force_reload = FALSE) {
                                       cat("\n=== STARTING COMPREHENSIVE LUAD TREATMENT ANALYSIS ===\n\n")

                                       # Step 1: Load and prepare data
                                       cat("Step 1: Loading and preparing data...\n")
                                       clinical_data <- self$data_loader$load_tcga_clinical(force_reload)
                                       biomarkers <- self$data_loader$load_biomarkers()

                                       # Merge all data
                                       comprehensive_data <- clinical_data %>%
                                         left_join(biomarkers$egfr %>% select(patient_id, egfr_mutation), by = "patient_id") %>%
                                         left_join(biomarkers$alk %>% select(patient_id, alk_fusion), by = "patient_id") %>%
                                         left_join(biomarkers$pdl1 %>% select(patient_id, pdl1_expression), by = "patient_id") %>%
                                         left_join(biomarkers$risk_scores %>% select(patient_id, risk_score, risk_group), by = "patient_id")

                                       cat("Data integration complete:", nrow(comprehensive_data), "patients\n\n")

                                       # Step 2: Generate treatment recommendations
                                       cat("Step 2: Generating treatment recommendations...\n")
                                       recommendations <- self$treatment_engine$recommend_treatments(comprehensive_data)
                                       cat("Treatment recommendations complete\n\n")

                                       # Step 3: Analyze results
                                       cat("Step 3: Analyzing results...\n")
                                       analysis_results <- self$analysis_engine$analyze_recommendations(recommendations)
                                       cat("Analysis complete\n\n")

                                       # Step 4: Create visualizations
                                       cat("Step 4: Creating visualizations...\n")
                                       plots <- self$visualization_engine$create_all_plots(recommendations, analysis_results)
                                       cat("Visualizations complete\n\n")

                                       # Step 5: Generate reports
                                       cat("Step 5: Generating comprehensive reports...\n")
                                       self$report_generator$generate_comprehensive_report(recommendations, analysis_results, plots)
                                       cat("Reports complete\n\n")

                                       # Step 6: Save results
                                       cat("Step 6: Saving results...\n")
                                       self$save_results(recommendations, analysis_results)
                                       cat("Results saved\n\n")

                                       # Print summary
                                       self$print_analysis_summary(analysis_results)

                                       cat("=== ANALYSIS COMPLETE ===\n")
                                       cat("Check the 'results' folder for all outputs\n\n")

                                       return(list(
                                         data = recommendations,
                                         analysis = analysis_results,
                                         plots = plots
                                       ))
                                     },

                                     # Save all results
                                     save_results = function(recommendations, analysis_results) {
                                       # Main results
                                       write_csv(recommendations, file.path(self$config$results_dir, "comprehensive_treatment_recommendations.csv"))

                                       # Analysis results
                                       write_csv(analysis_results$treatment_distribution,
                                                 file.path(self$config$results_dir, "treatment_distribution.csv"))

                                       if (!is.null(analysis_results$biomarker_impact$egfr)) {
                                         write_csv(analysis_results$biomarker_impact$egfr,
                                                   file.path(self$config$results_dir, "egfr_impact_analysis.csv"))
                                       }

                                       if (!is.null(analysis_results$biomarker_impact$pdl1)) {
                                         write_csv(analysis_results$biomarker_impact$pdl1,
                                                   file.path(self$config$results_dir, "pdl1_impact_analysis.csv"))
                                       }

                                       write_csv(analysis_results$stage_analysis,
                                                 file.path(self$config$results_dir, "stage_analysis.csv"))

                                       # Save summary as JSON for programmatic access
                                       jsonlite::write_json(analysis_results$summary,
                                                            file.path(self$config$results_dir, "analysis_summary.json"),
                                                            pretty = TRUE)
                                     },

                                     # Print analysis summary
                                     print_analysis_summary = function(analysis) {
                                       cat("ANALYSIS SUMMARY\n")
                                       cat("================\n")
                                       cat("Total patients analyzed:", analysis$summary$total_patients, "\n")
                                       cat("Patients with ECOG scores:", analysis$summary$patients_with_ecog, "\n")
                                       cat("Patients with biomarkers:", analysis$summary$patients_with_biomarkers, "\n")
                                       cat("Biomarker-guided treatment rate:", analysis$summary$biomarker_guided_rate, "%\n")
                                       cat("Curative intent rate:", analysis$summary$curative_intent_rate, "%\n")
                                       cat("Average confidence score:", analysis$summary$average_confidence, "\n")
                                       cat("================\n\n")

                                       cat("TOP 5 TREATMENTS:\n")
                                       top_treatments <- head(analysis$treatment_distribution, 5)
                                       for (i in 1:nrow(top_treatments)) {
                                         cat(i, ".", top_treatments$primary_treatment[i], ":",
                                             top_treatments$n[i], "patients (", top_treatments$percentage[i], "%)\n")
                                       }
                                     }
                                   )
)

# =============================================================================
# 8. MAIN EXECUTION FUNCTIONS
# =============================================================================

# Initialize and run the system
initialize_luad_system <- function() {
  # Load required packages
  load_required_packages()

  # Initialize system
  system <- LUADTreatmentSystem$new()

  return(system)
}

# Quick start function
run_luad_analysis <- function(force_reload = FALSE) {
  system <- initialize_luad_system()
  results <- system$run_complete_analysis(force_reload)
  return(results)
}

# =============================================================================
# 9. USAGE INSTRUCTIONS
# =============================================================================
results <- run_luad_analysis()
