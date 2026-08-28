# ==============================================================================
# 03_clinical_feature_engineering.R
# Deriving clinically interpretable pathway features from coded observations
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
pathway_referral <- fread(
  file.path(derived_dir, 'pathway_referrals.csv'),
  colClasses = list(character = c('patid', 'obsid'))
)

# Expected observation fields used below:
# patid, obsid, medcodeid, Term, cancer_type, days_from_index

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
# the final sequence alphabet because it was not used as a separate pathway state. (Based on my research objectives)

# ------------------------------------------------------------------------------
# 4. Referral urgency feature
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
# 5. Aggregate QC
# ------------------------------------------------------------------------------

print(pathway_obs[, .N, by = .(cancer_type, obs_state_raw)][order(cancer_type, -N)])
print(referral_urgency[, .N, by = .(cancer_type, ref_state)][order(cancer_type, -N)])

# Save only inside an approved secure environment. These outputs are excluded by
# .gitignore and should never be committed to GitHub.
# fwrite(pathway_obs, file.path(derived_dir, 'observations_featured.csv'))
# fwrite(referral_urgency, file.path(derived_dir, 'referrals_featured.csv'))
