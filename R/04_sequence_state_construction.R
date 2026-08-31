# ==============================================================================
# 04_sequence_state_construction.R
# Construct patient-by-time-slot states for three pathway channels
# ==============================================================================
# Public portfolio version.
# ===============================================================================

library(data.table)

# ------------------------------------------------------------------------------
# 1. Define non-uniform temporal slots
# ------------------------------------------------------------------------------
# The 24-month pre-diagnostic period is represented using broader intervals early
# in the pathway and finer monthly intervals closer to diagnosis/index date.

slot_def <- data.table(
  slot = 1:8,
  day_start = c(548, 366, 274, 183, 92, 61, 31, 1),
  day_end   = c(730, 547, 365, 273, 182, 91, 60, 30),
  label = c('m24-19', 'm18-13', 'm12-10', 'm9-7',
            'm6-4', 'm3', 'm2', 'm1')
)
slot_def[, n_days := day_end - day_start + 1L]

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

# ------------------------------------------------------------------------------
# 2. Load featured pathway tables
# ------------------------------------------------------------------------------

derived_dir <- file.path('data', 'derived')

pathway_obs <- fread(file.path(derived_dir, 'observations_featured.csv'),
                     colClasses = list(character = 'patid'))
pathway_consult <- fread(file.path(derived_dir, 'pathway_consultations.csv'),
                         colClasses = list(character = c('patid', 'consid')))
referral_urgency <- fread(file.path(derived_dir, 'referrals_featured.csv'),
                          colClasses = list(character = 'patid'))
cohort <- fread(file.path(derived_dir, 'analytical_cohort.csv'),
                colClasses = list(character = 'patid'))

# ------------------------------------------------------------------------------
# 3. Channel 1: observation state
# ------------------------------------------------------------------------------
# If multiple observation states occur in one patient-slot, retain the most
# clinically salient state according to an explicit priority rule.

pathway_obs[, slot := assign_slot(days_from_index)]
pathway_obs[obs_state_raw == 'SCR', obs_state_raw := 'OTH']

obs_priority <- c(EMERG = 1, RED = 2, IMG = 3, OTH = 4)

obs_state <- pathway_obs[!is.na(slot), .(
  obs_state = obs_state_raw[which.min(obs_priority[obs_state_raw])]
), by = .(patid, cancer_type, slot)]

# ------------------------------------------------------------------------------
# 4. Channel 2: consultation intensity
# ------------------------------------------------------------------------------
# Consultation counts are normalised to a 30-day rate because the sequence slots
# have unequal lengths. Distinct consultation IDs are counted to avoid duplicate
# records inflating activity.

pathway_consult[, slot := assign_slot(days_from_index)]

consult_counts <- pathway_consult[!is.na(slot), .(
  n_consult = uniqueN(consid)
), by = .(patid, cancer_type, slot)]

consult_counts <- merge(
  consult_counts,
  slot_def[, .(slot, n_days)],
  by = 'slot', all.x = TRUE
)

consult_counts[, rate_30d := n_consult / (n_days / 30)]

# Fixed thresholds were used in the study to create comparable intensity states
# across all time slots.
# Thresholds for LOW / MEDIUM / HIGH consultation intensity were selected after
# inspecting the empirical distribution of consultation rates across the
# analytical cohort.

consult_counts[, cons_state := fcase(
  rate_30d < 2, 'LOW',
  rate_30d < 4, 'MED',
  default = 'HIGH'
)]

# ------------------------------------------------------------------------------
# 5. Channel 3: referral state
# ------------------------------------------------------------------------------

referral_urgency[, slot := assign_slot(days_from_index)]
ref_priority <- c(EMERGENCY = 1, URGENT = 2, ROUTINE = 3)

ref_state <- referral_urgency[!is.na(slot), .(
  ref_state = ref_state[which.min(ref_priority[ref_state])]
), by = .(patid, cancer_type, slot)]

# ------------------------------------------------------------------------------
# 6. Aggregate prevalence QC
# ------------------------------------------------------------------------------

check_prevalence <- function(dt, state_col, label) {
  p <- dt[, .N, by = c('cancer_type', 'slot', state_col)]
  p[, pct := 100 * N / sum(N), by = .(cancer_type, slot)]
  message('\n=== ', label, ' ===')
  print(dcast(p, cancer_type + slot ~ get(state_col), value.var = 'pct', fill = 0))
}

check_prevalence(obs_state, 'obs_state', 'OBSERVATION')
check_prevalence(consult_counts, 'cons_state', 'CONSULTATION')
check_prevalence(ref_state, 'ref_state', 'REFERRAL')

# ------------------------------------------------------------------------------
# 7. Build complete patient x slot grids
# ------------------------------------------------------------------------------
# Explicitly creating all eight slots ensures that absence of activity is encoded
# as a meaningful NONE state rather than as missing data.

build_wide <- function(cancer, cohort, obs_state, consult_counts, ref_state) {
  patients <- cohort[cancer_type == cancer, unique(patid)]
  grid <- CJ(patid = patients, slot = 1:8)

  grid <- merge(grid,
                obs_state[cancer_type == cancer, .(patid, slot, obs_state)],
                by = c('patid', 'slot'), all.x = TRUE)
  grid <- merge(grid,
                consult_counts[cancer_type == cancer, .(patid, slot, cons_state)],
                by = c('patid', 'slot'), all.x = TRUE)
  grid <- merge(grid,
                ref_state[cancer_type == cancer, .(patid, slot, ref_state)],
                by = c('patid', 'slot'), all.x = TRUE)

  grid[is.na(obs_state),  obs_state  := 'NONE']
  grid[is.na(cons_state), cons_state := 'NONE']
  grid[is.na(ref_state),  ref_state  := 'NONE']

  list(
    obs  = dcast(grid, patid ~ slot, value.var = 'obs_state'),
    cons = dcast(grid, patid ~ slot, value.var = 'cons_state'),
    ref  = dcast(grid, patid ~ slot, value.var = 'ref_state')
  )
}

wide_breast <- build_wide('breast cancer', cohort, obs_state, consult_counts, ref_state)
wide_lung   <- build_wide('lung cancer', cohort, obs_state, consult_counts, ref_state)

# Save only in the approved secure environment.
# saveRDS(wide_breast, file.path(derived_dir, 'wide_breast_sequences.rds'))
# saveRDS(wide_lung,   file.path(derived_dir, 'wide_lung_sequences.rds'))
