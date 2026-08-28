# ==============================================================================
# 02_pathway_exploration.R
# Exploratory assessment of three pre-diagnostic pathway channels (Observations, consultations and referrals)
# ==============================================================================
# Public portfolio version.
# No patient-level data, restricted file paths, project IDs, or derived results
# from the source dataset are included in this repository.
# ===============================================================================

library(data.table)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load derived pathway tables
# ------------------------------------------------------------------------------
# These files are expected to have been created in the secure environment by the
# cohort-construction / pathway-extraction workflow code.

derived_dir <- file.path('data', 'derived')

pathway_obs <- fread(
  file.path(derived_dir, 'pathway_observations.csv'),
  colClasses = list(character = c('patid', 'medcodeid', 'obsid', 'consid'))
)
pathway_consult <- fread(
  file.path(derived_dir, 'pathway_consultations.csv'),
  colClasses = list(character = c('patid', 'consid'))
)
pathway_referral <- fread(
  file.path(derived_dir, 'pathway_referrals.csv'),
  colClasses = list(character = c('patid', 'obsid'))
)

# Dates should be parsed immediately after import. #Crucial Step
pathway_obs[, obsdate := as.Date(obsdate)]
pathway_consult[, consdate := as.Date(consdate)]
pathway_referral[, obsdate := as.Date(obsdate)]

# ------------------------------------------------------------------------------
# 2. Time relative to the index diagnosis
# ------------------------------------------------------------------------------
# Positive values denote days before diagnosis.

pathway_obs[, days_from_index := as.numeric(index_date - obsdate)]
pathway_consult[, days_from_index := as.numeric(index_date - consdate)]
pathway_referral[, days_from_index := as.numeric(index_date - obsdate)]

# Keep the pre-diagnostic window used by the study. (2 years prior to index/diagnosis)
pathway_obs <- pathway_obs[days_from_index >= 1 & days_from_index <= 730]
pathway_consult <- pathway_consult[days_from_index >= 1 & days_from_index <= 730]
pathway_referral <- pathway_referral[days_from_index >= 1 & days_from_index <= 730]

# ------------------------------------------------------------------------------
# 3. Channel-level descriptive checks
# ------------------------------------------------------------------------------
# Aggregate QC only. 

summarise_channel <- function(dt, date_col, label) {
  by_patient <- dt[, .N, by = patid]
  out <- data.table(
    channel = label,
    records = nrow(dt),
    patients = uniqueN(dt$patid),
    median_records_per_patient = median(by_patient$N),
    q1_records_per_patient = quantile(by_patient$N, 0.25),
    q3_records_per_patient = quantile(by_patient$N, 0.75)
  )
  out
}

channel_summary <- rbindlist(list(
  summarise_channel(pathway_obs, 'obsdate', 'observation'),
  summarise_channel(pathway_consult, 'consdate', 'consultation'),
  summarise_channel(pathway_referral, 'obsdate', 'referral')
))
print(channel_summary)

# ------------------------------------------------------------------------------
# 4. Temporal distribution by month before diagnosis
# ------------------------------------------------------------------------------

monthly_distribution <- function(dt, value_name) {
  dt[
    , .(n_events = .N, n_patients = uniqueN(patid)),
    by = .(cancer_type, month_pre = ceiling(days_from_index / 30.44))
  ][order(cancer_type, month_pre)]
}

monthly_obs <- monthly_distribution(pathway_obs)
monthly_consult <- monthly_distribution(pathway_consult)
monthly_referral <- monthly_distribution(pathway_referral)

# Convert event counts to within-cancer percentages so temporal shapes are
# comparable across cancer groups.
monthly_obs[, pct := 100 * n_events / sum(n_events), by = cancer_type]
monthly_consult[, pct := 100 * n_events / sum(n_events), by = cancer_type]
monthly_referral[, pct := 100 * n_events / sum(n_events), by = cancer_type]

plot_monthly <- function(dt, title, y_label) {
  ggplot(dt, aes(x = month_pre, y = pct, colour = cancer_type, group = cancer_type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    scale_x_reverse(breaks = seq(1, 24, 2)) +
    labs(
      title = title,
      subtitle = 'Month 1 = the month immediately preceding the index diagnosis',
      x = 'Months before diagnosis',
      y = y_label,
      colour = 'Cancer type'
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
}

print(plot_monthly(monthly_obs,
                   'Primary-care observations before diagnosis',
                   '% of observations'))
print(plot_monthly(monthly_consult,
                   'Consultations before diagnosis',
                   '% of consultations'))
print(plot_monthly(monthly_referral,
                   'Referrals before diagnosis',
                   '% of referrals'))

# ------------------------------------------------------------------------------
# 5. Broad pre-diagnostic periods for term exploration
# ------------------------------------------------------------------------------
# These periods were used as an exploratory aid before the final sequence slots
# were defined. They are not the final sequence-state intervals.

pathway_obs[, exploratory_period := fcase(
  days_from_index >= 1   & days_from_index <= 90,  'last_3_months',
  days_from_index >= 91  & days_from_index <= 180, 'months_4_6',
  days_from_index >= 181 & days_from_index <= 365, 'months_7_12',
  days_from_index >= 366 & days_from_index <= 730, 'months_13_24',
  default = NA_character_
)]

term_frequency <- pathway_obs[
  !is.na(exploratory_period) & !is.na(cancer_type),
  .(n_patients = uniqueN(patid), n_obs = .N),
  by = .(cancer_type, exploratory_period, medcodeid)
][order(cancer_type, exploratory_period, -n_obs)]

# If an approved medical dictionary can be used inside the secure environment,
# term labels can be joined there for clinical review.

# End of public pathway-exploration workflow.
