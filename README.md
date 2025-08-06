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
