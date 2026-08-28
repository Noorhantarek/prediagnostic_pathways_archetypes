# ==============================================================================
# 03_clinical_feature_engineering.R
# Deriving clinically interpretable pathway features from coded observations/ consultations and referrals
# ==============================================================================
# Public portfolio version.
# The exact source-database codelists and restricted medical-dictionary extracts
# used in the research are not included. The pattern lists below are illustrative
# representations of the study logic and must not be treated as validated codelists.
# ===============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

derived_dir <- file.path('data', 'derived')

pathway_obs <- fread(
  file.path(derived_dir, 'pathway_observations_with_terms.csv'),
  colClasses = list(character = c('patid', 'medcodeid', 'obsid'))
)
pathway_consult <- fread(
  file.path(derived_dir, 'pathway_consultations.csv'),
  colClasses = list(character = c('patid', 'consid'))
)
pathway_referral <- fread(
  file.path(derived_dir, 'pathway_referrals.csv'),
  colClasses = list(character = c('patid', 'obsid'))
)

# Expected fields used below:
# observations:  patid, obsid, medcodeid, Term, cancer_type, days_from_index
# consultations: patid, consid, consdate, index_date, cancer_type
# referrals:     patid, obsid, refurgencyid, cancer_type

# ------------------------------------------------------------------------------
# 2. Guideline-informed clinical concepts
# ------------------------------------------------------------------------------
# The research grouped raw coded observations into a small number of clinically
# interpretable states. Exact production codelists were manually reviewed in the
# secure environment and are not distributed here.

lung_red_flag_patterns <- paste(c(
  'haemoptysis', 'cough', 'dyspnoea', 'shortness of breath',
  'chest pain', 'weight loss', 'loss of appetite', 'fatigue',
  'pneumonia', 'hoarse', 'wheez', 'pleural effusion'
), collapse = '|')

breast_red_flag_patterns <- paste(c(
  'breast lump', 'breast mass', 'nipple discharge', 'nipple retraction',
  'nipple change', 'breast pain', 'skin dimpling', 'axillary lump',
  'breast asymmetry', 'breast abnormal'
), collapse = '|')

imaging_patterns <- paste(c(
  'chest x-ray', 'radiograph', 'computed tomogra', 'CT scan',
  'mammogra', 'MRI', 'ultrasound', 'biopsy', 'bronchoscop', 'PET scan'
), collapse = '|')

screening_patterns <- paste(c(
  'breast screen', 'screening.*breast', 'screening mammog',
  'breast screening programme'
), collapse = '|')

emergency_patterns <- paste(c(
  'A&E attendance', 'accident and emergency', 'emergency hospital admission',
  'seen in accident and emergency'
), collapse = '|')

# Negation / context terms were used during manual review to reduce false-positive
# concept matches. Exact production exclusions may differ.
context_exclusions <- paste(c(
  'family history', 'denies', 'declined', 'refused', 'history of', 'leaflet'
), collapse = '|')

# ------------------------------------------------------------------------------
# 3. Create observation states
# ------------------------------------------------------------------------------

pathway_obs[, term_clean := tolower(trimws(Term))]

pathway_obs[, obs_state_raw := fcase(
  grepl(emergency_patterns, term_clean, ignore.case = TRUE), 'EMERG',

  cancer_type == 'lung cancer' &
    grepl(lung_red_flag_patterns, term_clean, ignore.case = TRUE) &
    !grepl(context_exclusions, term_clean, ignore.case = TRUE), 'RED',

  cancer_type == 'breast cancer' &
    grepl(breast_red_flag_patterns, term_clean, ignore.case = TRUE) &
    !grepl(context_exclusions, term_clean, ignore.case = TRUE), 'RED',

  grepl(imaging_patterns, term_clean, ignore.case = TRUE), 'IMG',

  cancer_type == 'breast cancer' &
    grepl(screening_patterns, term_clean, ignore.case = TRUE), 'SCR',

  default = 'OTH'
)]

# Screening was retained during clinical review but later collapsed into OTH for
# the final sequence alphabet because it was not used as a separate pathway state, according to my research objectives

# ------------------------------------------------------------------------------
# 4. Consultation-intensity feature
# ------------------------------------------------------------------------------
# Consultation activity required a different type of feature engineering from
# the observation channel. Rather than classifying coded clinical concepts, the
# aim was to represent the intensity and level of engagment with primary-care over time.
#
# The analytical timeline contains intervals of unequal duration. Raw numbers of
# consultations therefore cannot be compared directly across slots: a six-month
# interval would naturally contain more consultations than a one-month interval.
# Distinct consultations were consequently standardised to an approximate
# 30-day rate before assigning LOW / MED / HIGH activity states.
#
# The thresholds were informed by the empirical distribution of consultation
# rates in the analytical cohort. Median, 75th-percentile and 90th-percentile
# rates were inspected across temporal slots before selecting fixed cut-offs.
# Fixed rather than slot-specific thresholds were retained so that LOW, MED and
# HIGH have the same interpretation throughout the 24-month pathway.
#
# The numerical summaries used to inform this decision are calculated at run
# time and are not hard-coded in this public version.

pathway_consult[, consdate   := as.Date(consdate)]
pathway_consult[, index_date := as.Date(index_date)]
pathway_consult[, days_from_index := as.numeric(index_date - consdate)]

# Temporal intervals used in the study. These become progressively more granular
# as diagnosis approaches, allowing the analysis to capture late escalation.
slot_def <- data.table(
  slot      = 1:8,
  day_start = c(548, 366, 274, 183,  92, 61, 31,  1),
  day_end   = c(730, 547, 365, 273, 182, 91, 60, 30)
)
slot_def[, n_days := day_end - day_start + 1]

assign_slot <- function(days) {
  fcase(
    days >= 548 & days <= 730, 1L,
    days >= 366 & days <= 547, 2L,
    days >= 274 & days <= 365, 3L,
    days >= 183 & days <= 273, 4L,
    days >=  92 & days <= 182, 5L,
    days >=  61 & days <=  91, 6L,
    days >=  31 & days <=  60, 7L,
    days >=   1 & days <=  30, 8L,
    default = NA_integer_
  )
}

pathway_consult[, slot := assign_slot(days_from_index)]

# Count distinct consultation encounters rather than rows. Using unique
# consultation identifiers avoids treating multiple records attached to the
# same encounter as separate consultations.
consult_counts <- pathway_consult[
  !is.na(slot),
  .(n_consult = uniqueN(consid)),
  by = .(patid, cancer_type, slot)
]

consult_counts <- merge(
  consult_counts,
  slot_def[, .(slot, n_days)],
  by = 'slot',
  all.x = TRUE
)

# Standardise activity to an approximate 30-day rate so that consultation
# intensity is comparable across intervals of different lengths.
consult_counts[, consultation_rate := n_consult / (n_days / 30)]

# Inspect the empirical rate distribution before defining states. These summaries
# informed the final fixed thresholds but are not themselves patient-level data.
consult_rate_summary <- consult_counts[, .(
  median = round(median(consultation_rate, na.rm = TRUE), 2),
  p75    = round(quantile(consultation_rate, 0.75, na.rm = TRUE), 2),
  p90    = round(quantile(consultation_rate, 0.90, na.rm = TRUE), 2)
), by = slot][order(slot)]

print(consult_rate_summary)

# Final study-specific intensity states:
#   LOW  = fewer than 2 consultations per 30 days
#   MED  = 2 to fewer than 4 consultations per 30 days
#   HIGH = 4 or more consultations per 30 days
#
# These cut-offs provide a stable interpretation across all sequence positions.
consult_counts[, cons_state := fcase(
  consultation_rate < 2, 'LOW',
  consultation_rate < 4, 'MED',
  default = 'HIGH'
)]

# ------------------------------------------------------------------------------
# 5. Referral urgency feature
# ------------------------------------------------------------------------------
# Source urgency coding was supplemented by linked observation text where a clear
# fast-track / suspected-cancer referral indicator was present. This is a
# study-specific operationalisation, not a redefinition of the source database.

referral_terms <- unique(pathway_obs[, .(obsid, Term)])
referral_urgency <- merge(pathway_referral, referral_terms, by = 'obsid', all.x = TRUE)

referral_urgency[, ref_state := fcase(
  refurgencyid == 5, 'EMERGENCY',
  refurgencyid == 1, 'URGENT',
  grepl('fast track|2 week rule|suspected.*cancer', Term, ignore.case = TRUE), 'URGENT',
  refurgencyid %in% c(2, 3, 4), 'ROUTINE',
  default = 'ROUTINE'
)]

# ------------------------------------------------------------------------------
# 6. Aggregate QC
# ------------------------------------------------------------------------------

print(pathway_obs[, .N, by = .(cancer_type, obs_state_raw)][order(cancer_type, -N)])

consult_qc <- consult_counts[, .N, by = .(cancer_type, slot, cons_state)]
consult_qc[, pct := round(N / sum(N) * 100, 1), by = .(cancer_type, slot)]
print(dcast(
  consult_qc,
  cancer_type + slot ~ cons_state,
  value.var = 'pct',
  fill = 0
))

print(referral_urgency[, .N, by = .(cancer_type, ref_state)][order(cancer_type, -N)])

# Save only inside an approved secure environment. These outputs are excluded by
# .gitignore and should never be committed to GitHub.
# fwrite(pathway_obs, file.path(derived_dir, 'observations_featured.csv'))
# fwrite(consult_counts, file.path(derived_dir, 'consultations_featured.csv'))
# fwrite(referral_urgency, file.path(derived_dir, 'referrals_featured.csv'))

 
 



