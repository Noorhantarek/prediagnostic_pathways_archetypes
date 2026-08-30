# Pre-diagnostic Pathway Archetypes for Breast and Lung Cancer

Patients follow very different routes to a cancer diagnosis — some direct, others long and circuitous through repeated GP visits, referrals, and investigations. This project asks whether those different routes fall into distinguishable, clinically meaningful patterns ("archetypes"), and whether patients from more deprived areas are disproportionately likely to experience the more complex ones. It applies multichannel sequence analysis and unsupervised clustering to two years of pre-diagnostic primary-care activity for breast and lung cancer patients using CPRD data.

This repository contains annotated R code from an MSc dissertation exploring whether clinically meaningful **pre-diagnostic pathway archetypes** can be identified from longitudinal primary-care data using multichannel sequence analysis and unsupervised clustering.

The public repository focuses on the **analytical methodology** rather than reproducing the restricted-data analysis. No patient-level data, patient identifiers, source extracts, internal directory paths, project IDs, restricted medical-dictionary extracts, or disclosure-sensitive results are included.

## Analytical workflow

```text
Cohort construction
        ↓
Pathway exploration
        ↓
Clinical feature engineering
        ↓
Three-channel state construction
        ↓
Multichannel sequence analysis
        ↓
Optimal Matching distance estimation
        ↓
Hierarchical clustering and validation
        ↓
Pathway archetype interpretation
        ↓
Deprivation / equity analysis
```

## Repository structure

```text
R/
├── 01_cohort_construction.R
├── 02_pathway_exploration.R
├── 03_clinical_feature_engineering.R
├── 04_sequence_state_construction.R
├── 05_multichannel_sequence_analysis.R
├── 06_clustering_archetypes.R
└── 07_deprivation_analysis.R
```

### 01 — Cohort construction
Identifies eligible cancer patients, defines the index diagnosis, applies inclusion/exclusion criteria, and establishes the pre-diagnostic observation window.

In this analysis, out of 7 million uniqe patients from 1,127 compiled CPRD aurum raw files; 156,280 patients, 26.7 million observation records were identifed and included in the analysis. 

### 02 — Pathway exploration
Examines the observation, consultation, and referral channels before sequence construction, including temporal activity patterns and aggregate code-frequency summaries.

### 03 — Clinical feature engineering
Maps high-dimensional coded observations into clinically interpretable analytical features such as cancer-relevant red-flag symptoms, imaging/investigation activity, emergency presentation, and referral urgency. Exact production codelists used in the restricted environment are not distributed publicly.

### 04 — Sequence-state construction
Divides the 24-month pre-diagnostic period into eight non-uniform time slots and creates three patient-level state channels:

- **Observation:** no activity / other clinical activity / imaging / red flag / emergency
- **Consultation:** no consultation / low / medium / high intensity
- **Referral:** no referral / routine / urgent / emergency

Missing activity is explicitly represented as a `NONE` state rather than treated as missing data.

### 05 — Multichannel sequence analysis
Creates TraMineR sequence objects for the three channels and estimates multichannel Optimal Matching distances using transition-rate-derived substitution costs and equal channel weights.

### 06 — Clustering and pathway archetypes
Applies Ward.D2 hierarchical clustering to the multichannel distance matrices. Candidate solutions are evaluated using silhouette width, cluster size, multiple internal validity indices, sequence plots, temporal composition, and clinical interpretability.

### 07 — Deprivation analysis
Examines whether pathway archetype membership varies across deprivation quintiles using descriptive distributions, chi-squared tests, standardised residuals, and ordered trend tests.

The analysis identified 4 distinct pre-diagnostic pathway archetypes for Breast Cancer, and 3 for Lung cancer cohort, ranging from streamlined, low-activity routes to markedly longer and more complex ones. Full results, including the equity analysis by deprivation quintile, are reported in the dissertation itself.

## Software

The workflow uses R packages including:

- `data.table`
- `TraMineR`
- `cluster`
- `WeightedCluster`
- `ggplot2`

## Data governance

The original research was conducted using restricted-access health data. This public code has therefore been deliberately generalised:

- data paths and filenames are placeholders;
- no row-level data are supplied;
- no patient identifiers are shown;
- no restricted source-database codelists or dictionary extracts are distributed;
- source-derived cohort counts and statistical results are not hard-coded;
- code that would save restricted outputs is disabled/commented in the public version.

The scripts are intended to demonstrate the analytical workflow and will not run as-is without appropriately structured, authorised source data.

## What I'd do differently
1. Restrict to a coding-stable period. Archetype membership tracked calendar year almost perfectly — median diagnosis year ran from 2001 to 2013 across breast archetypes. Restricting to 2015 onwards would test whether the archetypes are pathway features or artefacts of changing recording practice.
2. Adjust for confounders. Archetypes differed by twelve years in median age. Multinomial regression of archetype membership on deprivation, adjusted for age, sex and calendar period, would establish whether the deprivation association is independent.
3. Cluster stability. Internal validity indices assess coherence within the data analysed; they say nothing about whether the same solution would emerge from a different sample. Bootstrap resampling would test this.

## Author

Nourhan Ibrahim  
MSc dissertation, 2026
