# Cohort Construction for Longitudinal Primary-Care Data

This repository contains an **annotated R implementation of a cohort-construction workflow** developed for research using longitudinal primary-care electronic health record (EHR) data.

The public code focuses on the analytical methodology: mapping diagnosis code lists, identifying index diagnoses, applying eligibility criteria, excluding prior cancer, and extracting a 24-month pre-diagnostic observation window.

## What this repository demonstrates

- Working with large longitudinal EHR extracts in R
- Mapping SNOMED CT concepts to source-database medical codes
- Memory-conscious, file-by-file processing with `data.table`
- Defining an index diagnosis from longitudinal records
- Building an auditable inclusion/exclusion pipeline
- Identifying prior diagnoses relative to a patient-specific index date
- Defining and extracting a pre-diagnostic observation window
- Applying reproducible quality-control checks

## Cohort-construction workflow

```text
Diagnosis code lists
        |
        v
Map clinical concepts to source medical codes
        |
        v
Search longitudinal observation records
        |
        v
Identify earliest qualifying diagnosis
        |
        v
Define patient-specific index date
        |
        v
Link patient eligibility information
        |
        v
Apply inclusion / exclusion criteria
        |
        +-- acceptable patient record
        +-- eligible age range
        +-- valid index date
        +-- registered at index date
        +-- >= 24 months prior registration
        +-- no prior cancer diagnosis
        +-- cancer-specific eligibility criteria
        |
        v
Final analytical cohort
        |
        v
Define index - 24 months to index - 1 day
        |
        v
Extract pre-diagnostic healthcare events
```

## Repository structure

```text
.
├── README.md
├── .gitignore
└── R/
    └── cohort_construction_public.R
```

## Data availability and confidentiality

The original research was conducted using **restricted-access patient-level health data**.

**No patient-level data, patient identifiers, restricted data extracts, derived patient-level outputs, internal directory paths, project identifiers, or credentials are included in this repository.**

The code has been deliberately generalised for public demonstration. Source paths and filenames are placeholders, and some variable names have been made more descriptive. The script therefore **will not run as-is** without appropriately structured source data and approved code lists.

Researchers using restricted-access data should follow the governance, disclosure-control and code-sharing requirements of their own data provider and institution.

## Software

The workflow uses:

- R
- `data.table`
- `lubridate`

## Main script

[`R/cohort_construction_public.R`](R/cohort_construction_public.R) contains the annotated public implementation.

The script is organised into:

1. Configuration
2. Diagnosis code-list mapping
3. Candidate diagnosis identification
4. Index-date definition
5. Patient-data linkage
6. Cohort eligibility criteria
7. Prior-cancer exclusion
8. Cancer-specific eligibility
9. Pre-diagnostic window definition
10. Event extraction
11. Quality-control checks

## Reproducibility note

The public repository demonstrates the **analytical logic rather than reproducing the restricted analysis**. Counts and patient-level examples from the source data are intentionally omitted.

## Author

Nourhan Ibrahim

MSc research project, 2026.
