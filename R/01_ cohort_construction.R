# ==============================================================================
# Cohort construction for longitudinal primary-care data
# ==============================================================================
# Purpose
# -------
# Public, annotated implementation of the cohort-construction workflow used in
# a research project examining pre-diagnostic pathways in longitudinal EHR data.
#
# IMPORTANT:
# - No patient-level data are included.
# - File paths, project identifiers and source-specific filenames are generalised.
# - Variable names are made more descriptive where practical.
# - This script is illustrative and will not run without appropriately structured
#   source data and code lists.
# - Check the data provider's governance requirements before adapting this code
#   to any restricted-access dataset.
# ==============================================================================

library(data.table)
library(lubridate)

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

# Keep paths outside the analytical logic so the workflow is portable.
# Never hard-code restricted environment paths in a public repository.
data_dir     <- "path/to/restricted/source/data"
codelist_dir <- "path/to/codelists"
output_dir   <- "path/to/approved/output"

# Study boundaries should be set from the approved study protocol.
study_start     <- as.Date("YYYY-MM-DD")
study_end       <- as.Date("YYYY-MM-DD")
extraction_date <- as.Date("YYYY-MM-DD")

minimum_age    <- 18L
maximum_age    <- 100L
lookback_days  <- 730L   # approximately 24 months


# ------------------------------------------------------------------------------
# 2. Load and map diagnosis code lists
# ------------------------------------------------------------------------------

# Read public/approved SNOMED CT code lists as character values. Keeping clinical
# codes as character strings prevents long identifiers being converted to
# scientific notation.

lung_codes <- fread(
  file.path(codelist_dir, "lung_cancer_snomed.csv"),
  colClasses = "character"
)

breast_codes <- fread(
  file.path(codelist_dir, "breast_cancer_snomed.csv"),
  colClasses = "character"
)

# The source data use an internal medical-code identifier. Map SNOMED CT concepts
# to that identifier using the approved medical dictionary.
medical_dictionary <- fread(
  file.path(data_dir, "medical_dictionary.txt"),
  sep = "\t",
  header = TRUE,
  colClasses = "character",
  encoding = "UTF-8"
)

lung_medcodes <- medical_dictionary[
  SnomedCTConceptId %in% lung_codes$snomed_id
]

breast_medcodes <- medical_dictionary[
  SnomedCTConceptId %in% breast_codes$snomed_id
]

# Remove duplicate internal medical codes before searching the longitudinal data.
lung_medcodes   <- unique(lung_medcodes,   by = "MedCodeId")
breast_medcodes <- unique(breast_medcodes, by = "MedCodeId")

lung_medcodes[,   cancer_type := "lung cancer"]
breast_medcodes[, cancer_type := "breast cancer"]

diagnosis_codes <- rbindlist(
  list(
    lung_medcodes[,   .(medical_code = MedCodeId, cancer_type)],
    breast_medcodes[, .(medical_code = MedCodeId, cancer_type)]
  ),
  use.names = TRUE
)


# ------------------------------------------------------------------------------
# 3. Identify candidate diagnosis records
# ------------------------------------------------------------------------------

# Large EHR extracts may be split across many files. Processing files
# sequentially reduces peak memory use compared with loading the full dataset.

observation_files <- list.files(
  data_dir,
  pattern = "observation.*\\.txt$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

diagnosis_records <- vector("list", length(observation_files))

for (i in seq_along(observation_files)) {

  message(sprintf(
    "Processing observation file %d of %d",
    i, length(observation_files)
  ))

  chunk <- fread(
    observation_files[i],
    sep = "\t",
    header = TRUE,
    select = c("patid", "obsdate", "medcodeid"),
    na.strings = c("", "NA"),
    colClasses = list(character = c("patid", "medcodeid")),
    nThread = max(1L, parallel::detectCores() - 1L)
  )

  # Restrict each chunk to diagnosis codes of interest.
  matched <- chunk[medcodeid %in% diagnosis_codes$medical_code]

  if (nrow(matched) > 0L) {
    diagnosis_records[[i]] <- matched
  }

  rm(chunk, matched)
  gc()
}

cancer_observations <- rbindlist(
  diagnosis_records,
  use.names = TRUE,
  fill = TRUE
)

# Attach cancer-type labels.
cancer_observations <- merge(
  cancer_observations,
  diagnosis_codes,
  by.x = "medcodeid",
  by.y = "medical_code",
  all.x = TRUE
)

# Parse dates immediately after loading/combining them.
cancer_observations[, obsdate := as.Date(obsdate, format = "%d/%m/%Y")]


# ------------------------------------------------------------------------------
# 4. Define the index diagnosis
# ------------------------------------------------------------------------------

# Order diagnosis records chronologically within patient and retain the earliest
# qualifying diagnosis as the index event.
setorder(cancer_observations, patid, obsdate)

index_dates <- cancer_observations[
  ,
  .SD[1L],
  by = patid
][
  ,
  .(
    patid,
    cancer_type,
    index_date = obsdate,
    index_medical_code = medcodeid
  )
]

# Restrict index diagnoses to the approved study period.
index_dates <- index_dates[
  !is.na(index_date) &
    index_date >= study_start &
    index_date <= study_end
]


# ------------------------------------------------------------------------------
# 5. Link patient-level eligibility information
# ------------------------------------------------------------------------------

patient <- fread(
  file.path(data_dir, "patient.txt"),
  sep = "\t",
  header = TRUE,
  colClasses = list(character = c("patid", "pracid"))
)

date_columns <- c("regstartdate", "regenddate", "cprd_ddate")

patient[
  ,
  (date_columns) := lapply(
    .SD,
    as.Date,
    format = "%d/%m/%Y"
  ),
  .SDcols = date_columns
]

# A missing registration end date can represent ongoing registration at the
# extraction date. Confirm this interpretation for the dataset being used.
patient[
  ,
  regenddate_clean := fifelse(
    is.na(regenddate),
    extraction_date,
    regenddate
  )
]

cohort <- merge(
  index_dates,
  patient[
    ,
    .(
      patid,
      pracid,
      gender,
      yob,
      regstartdate,
      regenddate_clean,
      acceptable,
      cprd_ddate
    )
  ],
  by = "patid",
  all.x = TRUE
)


# ------------------------------------------------------------------------------
# 6. Apply cohort eligibility criteria
# ------------------------------------------------------------------------------

# Helper function for an auditable filtering workflow. It reports how many rows
# remain after each criterion without hard-coding any restricted-data results.
log_step <- function(step_name, dt, previous_n) {

  current_n  <- nrow(dt)
  n_excluded <- previous_n - current_n

  message(sprintf(
    "%-45s | remaining: %s | excluded: %s",
    step_name,
    format(current_n, big.mark = ","),
    format(n_excluded, big.mark = ",")
  ))

  current_n
}

n <- nrow(cohort)

# Criterion 1: retain records meeting the source database's acceptable-patient
# quality criterion.
cohort <- cohort[acceptable == 1]
n <- log_step("Acceptable patient record", cohort, n)

# Criteria 2-3: adult patients within the pre-specified age range.
# Only year of birth is available here, so age is approximate to calendar year.
cohort[, age_at_index := year(index_date) - yob]

cohort <- cohort[
  !is.na(age_at_index) &
    age_at_index >= minimum_age &
    age_at_index <= maximum_age
]
n <- log_step("Age within eligible range", cohort, n)

# Criterion 4: non-missing index date.
cohort <- cohort[!is.na(index_date)]
n <- log_step("Non-missing index date", cohort, n)

# Criterion 5: index diagnosis must occur while the patient is registered.
cohort[
  ,
  within_registration :=
    index_date >= regstartdate &
    index_date <= regenddate_clean
]

cohort <- cohort[within_registration == TRUE]
n <- log_step("Index within registration period", cohort, n)

# Criterion 6: at least 24 months of observable history before index.
# Using a day-based threshold avoids relying on an average month length.
cohort[, history_days := as.integer(index_date - regstartdate)]

cohort <- cohort[
  !is.na(history_days) &
    history_days >= lookback_days
]
n <- log_step("At least 24 months prior registration", cohort, n)


# ------------------------------------------------------------------------------
# 7. Exclude patients with a prior cancer diagnosis
# ------------------------------------------------------------------------------

# Load an approved all-cancer code list. The exact exclusions should follow the
# study protocol; for example, some studies exclude non-melanoma skin cancer
# from the definition of prior malignancy.
prior_cancer_codes <- fread(
  file.path(codelist_dir, "prior_cancer_codes.csv"),
  colClasses = "character"
)

cohort_patients <- unique(cohort$patid)
index_lookup <- cohort[, .(patid, index_date)]

prior_cancer_matches <- vector("list", length(observation_files))

for (i in seq_along(observation_files)) {

  chunk <- fread(
    observation_files[i],
    sep = "\t",
    header = TRUE,
    select = c("patid", "obsdate", "medcodeid"),
    na.strings = c("", "NA"),
    colClasses = list(character = c("patid", "medcodeid")),
    nThread = max(1L, parallel::detectCores() - 1L)
  )

  # Restrict early: first to the candidate cohort, then to prior-cancer codes.
  chunk <- chunk[
    patid %in% cohort_patients &
      medcodeid %in% prior_cancer_codes$medcodeid
  ]

  if (nrow(chunk) > 0L) {
    prior_cancer_matches[[i]] <- chunk
  }

  rm(chunk)
  gc()
}

prior_cancer_observations <- rbindlist(
  prior_cancer_matches,
  use.names = TRUE,
  fill = TRUE
)

prior_cancer_observations[
  ,
  obsdate := as.Date(obsdate, format = "%d/%m/%Y")
]

prior_cancer_observations <- merge(
  prior_cancer_observations,
  index_lookup,
  by = "patid",
  all.x = TRUE
)

prior_cancer_patients <- unique(
  prior_cancer_observations[
    !is.na(obsdate) &
      !is.na(index_date) &
      obsdate < index_date,
    patid
  ]
)

cohort <- cohort[!patid %in% prior_cancer_patients]
n <- log_step("No prior cancer diagnosis", cohort, n)


# ------------------------------------------------------------------------------
# 8. Apply cancer-specific eligibility rule
# ------------------------------------------------------------------------------

# Example of a pre-specified cancer-specific criterion used in the study.
# Confirm coding conventions and protocol requirements before reuse.
cohort <- cohort[
  !(cancer_type == "breast cancer" & gender == 1)
]
n <- log_step("Cancer-specific sex criterion", cohort, n)


# ------------------------------------------------------------------------------
# 9. Define the pre-diagnostic observation window
# ------------------------------------------------------------------------------

# The analytical window spans approximately 24 months before diagnosis and
# excludes the index date itself.
cohort[
  ,
  `:=`(
    window_start = index_date - lookback_days,
    window_end   = index_date - 1L
  )
]

window_lookup <- cohort[
  ,
  .(patid, index_date, window_start, window_end)
]


# ------------------------------------------------------------------------------
# 10. Extract events occurring within the pre-diagnostic window
# ------------------------------------------------------------------------------

pathway_observations <- vector("list", length(observation_files))

for (i in seq_along(observation_files)) {

  chunk <- fread(
    observation_files[i],
    sep = "\t",
    header = TRUE,
    select = c(
      "patid", "consid", "obsid",
      "obsdate", "medcodeid", "obstypeid"
    ),
    na.strings = c("", "NA"),
    colClasses = list(
      character = c("patid", "consid", "obsid", "medcodeid")
    ),
    nThread = max(1L, parallel::detectCores() - 1L)
  )

  # Restrict to cohort members before merging to minimise memory use.
  chunk <- chunk[patid %in% cohort$patid]

  if (nrow(chunk) > 0L) {

    chunk[, obsdate := as.Date(obsdate, format = "%d/%m/%Y")]

    chunk <- merge(
      chunk,
      window_lookup,
      by = "patid",
      all.x = TRUE
    )

    chunk <- chunk[
      !is.na(obsdate) &
        obsdate >= window_start &
        obsdate <= window_end
    ]

    if (nrow(chunk) > 0L) {
      pathway_observations[[i]] <- chunk
    }
  }

  rm(chunk)
  gc()
}

pathway_observations <- rbindlist(
  pathway_observations,
  use.names = TRUE,
  fill = TRUE
)


# ------------------------------------------------------------------------------
# 11. Quality-control checks
# ------------------------------------------------------------------------------

# Keep QC checks aggregate in public code. Do not print individual patient
# records or restricted-data results to a public repository.

stopifnot(uniqueN(cohort$patid) == nrow(cohort))
stopifnot(all(cohort$index_date >= study_start))
stopifnot(all(cohort$index_date <= study_end))
stopifnot(all(cohort$window_end < cohort$index_date))
stopifnot(all(cohort$history_days >= lookback_days))

message("Cohort-construction workflow completed successfully.")


# ------------------------------------------------------------------------------
# 12. Outputs
# ------------------------------------------------------------------------------

# Do not commit restricted-data outputs to GitHub.
# Save only to an approved secure location when running in the authorised
# research environment.

# saveRDS(cohort, file.path(output_dir, "analytical_cohort.rds"))
# saveRDS(pathway_observations,
#         file.path(output_dir, "prediagnostic_observations.rds"))
