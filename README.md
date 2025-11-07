## Set-up Steps

1. Install UV: <https://github.com/astral-sh/uv>
2. Install git: <https://git-scm.com/downloads>
3. Navigate to where you want the folder
4. Run `git clone https://github.com/dylanbforde/lung-cancer-ODE.git` in your terminal.
5. Run `uv sync` in project directory
6. Next run `source .venv/bin/activate`

If you are using VSCode, you will need to select .venv as kernel for notebooks.

---

## ODE References

#### Single-cell RNA
<https://github.com/theislab/scvelo>
<https://github.com/theislab/scvelo_notebooks>
<https://github.com/basilkhuder/Seurat-to-RNA-Velocity>

#### Math-biochemical systems
<https://github.com/pysb/pysb>
<https://github.com/alubbock/thunor>
<https://github.com/RuleWorld/bionetgen>

#### DepMap Analysis
<https://github.com/johnbachman/depmap_analysis>

#### Omics Integration
<https://github.com/fraenkel-lab/OmicsIntegrator2>
<https://github.com/bowang-lab/IntegrAO>

#### Regulatory networks
<https://github.com/Murali-group/Beeline>
<https://github.com/asor-1/Gene-Regulatory-Network-Dynamics-Analysis>

Data For Beeline Downloaded Here:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907>

---
## LUAD Survival Prediction Model (Question 1)

This folder contains the code, data, and results for Question 1 of the Final Project:  
**Can we build a generalizable model to predict patient survival across LUAD datasets?**

#### Folder Structure
| File / Folder | Description |
|---------------|-------------|
| Final Project.qmd | Main analysis script (R Quarto)| 
| updated SA.R | Survival analysis script (R), added hazard ratios| 
| Question1.ipynb | Python-based version of the pipeline| 
| plots/ | KM curves and model plots (Lasso CV, top genes)| 
| output files/ | Processed CSV outputs (risk scores, top genes)| 

#### How the Model Was Built
1. **Data Source**:  
   - Expression data from TCGA-LUAD (FPKM, downloaded from UCSC Xena)  
   - Clinical data (overall survival time + status) via `TCGAbiolinks`

2. **Preprocessing**:  
   - Kept only tumor samples (`-01`)
   - Removed low-expression genes and samples with missing survival data

3. **Analysis Pipeline**:
   - Univariate Cox regression to filter top survival-associated genes  
   - Lasso-Cox regression (via `glmnet`) to build the final model  
   - Risk score calculated for each patient  
   - Patients grouped into High/Low risk using median cutoff  
   - KM curves and C-index used to assess model performance


# LUAD Treatment Recommendation System - Question 5

## Overview

This project implements an intelligent treatment recommendation system for Lung Adenocarcinoma (LUAD) patients using TCGA (The Cancer Genome Atlas) data. The system provides evidence-based treatment recommendations by analyzing clinical characteristics, biomarker profiles, and treatment outcomes.

## Project Structure

```
notebooks/q5/
├── code/
│   └── question 5.R           # Main R script with complete system implementation
├── results/
│   ├── plots/                 # Generated visualization plots
│   ├── reports/              # Analysis reports and summaries
│   │   ├── executive_summary.txt
│   │   └── comprehensive_treatment_report.html
│   └── comprehensive_treatment_recommendations.csv
└── processed/                 # Processed datasets
```

## System Features

### Core Capabilities
- **Biomarker-Guided Treatment Selection**: Analyzes EGFR mutations, ALK fusions, and PD-L1 expression
- **Clinical Factor Integration**: Considers age, ECOG performance status, smoking history, and cancer stage
- **Treatment Confidence Scoring**: Assigns confidence scores based on evidence quality
- **Comprehensive Reporting**: Generates detailed HTML reports and executive summaries
- **Visualization Suite**: Creates treatment distribution plots and confidence score analyses

### Key Components

1. **Data Processing Module**: 
   - TCGA data acquisition and cleaning
   - Biomarker standardization
   - Clinical variable normalization

2. **Treatment Recommendation Engine**:
   - Rule-based treatment selection algorithms
   - Evidence level assessment (Level 1-3)
   - Multidisciplinary treatment planning

3. **Quality Metrics Calculator**:
   - Data completeness assessment
   - Clinical appropriateness evaluation
   - Guideline concordance tracking

4. **Reporting System**:
   - HTML report generation
   - Executive summary creation
   - Treatment distribution analysis

## Installation & Setup

### Required R Packages
```r
# Core packages
install.packages(c("dplyr", "readr", "ggplot2", "tidyr", "stringr", "here"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("TCGAbiolinks", "survival", "survminer"))

# Optional packages (enhanced functionality)
install.packages(c("DT", "plotly", "knitr", "rmarkdown", "jsonlite", "R6"))
```

### Directory Setup
The system automatically creates the required directory structure:
- `results/` - Analysis outputs
- `results/plots/` - Generated visualizations
- `results/reports/` - Summary reports
- `processed/` - Processed datasets

## Usage

### Quick Start
```r
# Load the script
source("notebooks/q5/code/question 5.R")

# Run complete analysis
results <- run_luad_analysis()

# Force data reload (optional)
results <- run_luad_analysis(force_reload = TRUE)
```

### System Initialization
```r
# Initialize system manually
system <- initialize_luad_system()

# Run individual components
clinical_data <- system$load_clinical_data()
recommendations <- system$generate_recommendations(clinical_data)
```

## Key Findings (Latest Analysis)

### Patient Demographics
- **Total Patients Analyzed**: 524
- **Age Distribution**: Young (≤65): 45.2%, Middle-aged (66-74): 32.1%, Elderly (≥75): 22.7%
- **Gender Distribution**: Female: 51.3%, Male: 48.7%

### Treatment Patterns
- **Biomarker-Guided Treatments**: 21.6%
- **Curative Intent Achievable**: 38.7%
- **Average Treatment Confidence**: 0.76

### Treatment Distribution
| Treatment | Patients | Percentage |
|-----------|----------|------------|
| Standard Chemotherapy | 253 | 48.3% |
| Surgery (Lobectomy) | 123 | 23.5% |
| Osimertinib (EGFR TKI) | 73 | 13.9% |
| Alectinib (ALK inhibitor) | 30 | 5.7% |
| Surgery (Sublobar) | 21 | 4.0% |
| SBRT | 13 | 2.5% |
| Pembrolizumab combinations | 9 | 1.7% |

## Clinical Recommendations

Based on the analysis, the system generates evidence-based recommendations:

1. **Implement comprehensive molecular profiling** for all LUAD patients
2. **Establish multidisciplinary tumor boards** for complex cases
3. **Consider clinical trial enrollment** for novel therapies
4. **Regular monitoring and adaptation** of treatment plans
5. **Ensure adequate supportive care integration**
6. **PRIORITY**: Increase biomarker testing rates (currently suboptimal)

## Configuration Options

### Biomarker Thresholds
```r
biomarker_thresholds = list(
  pdl1_high = 50,      # PD-L1 high expression threshold
  pdl1_medium = 1,     # PD-L1 medium expression threshold
  egfr_confidence = 0.95,  # EGFR mutation confidence
  alk_confidence = 0.90    # ALK fusion confidence
)
```

### Clinical Thresholds
```r
clinical_thresholds = list(
  elderly_age = 75,     # Elderly patient threshold
  young_age = 65,       # Young patient threshold
  good_ecog = 1,        # Good performance status
  moderate_ecog = 2,    # Moderate performance status
  poor_ecog = 3         # Poor performance status
)
```

## Output Files

### Generated Reports
1. **executive_summary.txt**: Concise analysis overview
2. **comprehensive_treatment_report.html**: Detailed HTML report
3. **comprehensive_treatment_recommendations.csv**: Patient-level recommendations

### Visualization Outputs
- Treatment distribution plots
- Confidence score distributions
- Biomarker impact analyses
- Stage-specific treatment patterns

## Technical Architecture

The system uses an object-oriented R6 class structure:

```
LUADTreatmentSystem
├── DataProcessor (data cleaning and standardization)
├── TreatmentRecommendationEngine (algorithm implementation)
├── AnalyticsReporter (metrics and summaries)
└── VisualizationGenerator (plot creation)
```

## Data Sources

- **Primary Data**: TCGA-LUAD dataset
- **Clinical Variables**: Age, gender, smoking history, ECOG performance status
- **Molecular Data**: EGFR mutations, ALK fusions, PD-L1 expression
- **Outcome Data**: Overall survival, treatment response

## Performance Metrics

### Data Quality
- **Age Completeness**: >95%
- **ECOG Completeness**: Variable (system handles missing data)
- **Stage Completeness**: >90%
- **Biomarker Completeness**: Variable by marker type

### Clinical Appropriateness
- **Evidence Level 1 Rate**: Tracked per treatment
- **High Confidence Rate**: 76% average
- **Guideline Concordance**: Monitored continuously

## Future Enhancements

1. **Real-time Clinical Decision Support**: Integration with EHR systems
2. **Machine Learning Models**: Predictive outcome modeling
3. **Multi-cancer Support**: Extension beyond LUAD
4. **Interactive Dashboard**: Web-based interface
5. **Clinical Trial Matching**: Automated enrollment suggestions

## Dependencies

- R version ≥ 4.0.0
- Internet connection (for TCGA data download)
- Sufficient memory for genomic data processing (≥8GB recommended)

## License

This project is designed for research and educational purposes. Please ensure compliance with TCGA data usage policies.

## Support

For questions or issues:
1. Check the generated reports in `results/reports/`
2. Review the comprehensive analysis output
3. Consult the inline code documentation

---

**Last Updated**: November 2025  
**System Version**: 1.0  
**Analysis Date**: August 2025 (Latest run)
