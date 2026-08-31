# ==============================================================================
# 05_multichannel_sequence_analysis.R
# Multichannel sequence construction and Optimal Matching distances
# ==============================================================================
# Public portfolio version.
# ===============================================================================

library(data.table)
library(TraMineR)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

derived_dir <- file.path('data', 'derived')
wide_breast <- readRDS(file.path(derived_dir, 'wide_breast_sequences.rds'))
wide_lung   <- readRDS(file.path(derived_dir, 'wide_lung_sequences.rds'))

slot_labels <- c('m24-19', 'm18-13', 'm12-10', 'm9-7', 'm6-4', 'm3', 'm2', 'm1')

alph_obs  <- c('NONE', 'OTH', 'IMG', 'EMERG', 'RED')
alph_cons <- c('NONE', 'LOW', 'MED', 'HIGH')
alph_ref  <- c('NONE', 'ROUTINE', 'URGENT', 'EMERGENCY')

# ------------------------------------------------------------------------------
# 2. Sampling for pairwise distance estimation
# ------------------------------------------------------------------------------
# Pairwise sequence distances scale quadratically with sample size. The research
# therefore used a reproducible random sample for the multichannel distance and
# clustering stage.

set.seed(42)

sample_wide <- function(wide, n = 10000L) {
  idx <- sort(sample(nrow(wide$obs), min(n, nrow(wide$obs))))
  list(
    obs  = wide$obs[idx],
    cons = wide$cons[idx],
    ref  = wide$ref[idx]
  )
}

samp_breast <- sample_wide(wide_breast)
samp_lung   <- sample_wide(wide_lung)

# Aggregate comparison of the full and sampled observation-state distributions.
compare_state_prevalence <- function(full, sample) {
  full_tab <- prop.table(table(unlist(full$obs[, 2:9])))
  sample_tab <- prop.table(table(unlist(sample$obs[, 2:9])))
  rbind(full = full_tab, sample = sample_tab)
}

print(compare_state_prevalence(wide_breast, samp_breast))
print(compare_state_prevalence(wide_lung, samp_lung))

# ------------------------------------------------------------------------------
# 3. Create TraMineR sequence objects
# ------------------------------------------------------------------------------

make_sequences <- function(sample_data) {
  list(
    obs = seqdef(
      sample_data$obs, var = 2:9,
      alphabet = alph_obs, states = alph_obs,
      labels = c('No activity', 'Other clinical', 'Imaging',
                 'Emergency attendance', 'Red flag symptom'),
      cnames = slot_labels, xtstep = 1
    ),
    cons = seqdef(
      sample_data$cons, var = 2:9,
      alphabet = alph_cons, states = alph_cons,
      labels = c('No consultation', 'Low', 'Medium', 'High'),
      cnames = slot_labels, xtstep = 1
    ),
    ref = seqdef(
      sample_data$ref, var = 2:9,
      alphabet = alph_ref,
      states = c('NONE', 'ROUT', 'URG', 'EMERG'),
      labels = c('No referral', 'Routine', 'Urgent/2WW', 'Emergency'),
      cnames = slot_labels, xtstep = 1
    )
  )
}

seq_breast <- make_sequences(samp_breast)
seq_lung   <- make_sequences(samp_lung)

# ------------------------------------------------------------------------------
# 4. Sequence diagnostics
# ------------------------------------------------------------------------------
# Distribution and frequent-sequence plots provide an interpretable check of the
# state construction before computing distances.

seqdplot(seq_breast$obs,  border = NA, main = 'Breast: observations')
seqdplot(seq_breast$cons, border = NA, main = 'Breast: consultations')
seqdplot(seq_breast$ref,  border = NA, main = 'Breast: referrals')

seqdplot(seq_lung$obs,  border = NA, main = 'Lung: observations')
seqdplot(seq_lung$cons, border = NA, main = 'Lung: consultations')
seqdplot(seq_lung$ref,  border = NA, main = 'Lung: referrals')

# ------------------------------------------------------------------------------
# 5. Multichannel Optimal Matching
# ------------------------------------------------------------------------------
# Transition-rate-derived substitution costs (TRATE) are estimated separately
# for each channel. Channels receive equal weights in the multichannel distance.

compute_multichannel_om <- function(sequences) {
  substitution_costs <- list(
    seqsubm(sequences$obs,  method = 'TRATE'),
    seqsubm(sequences$cons, method = 'TRATE'),
    seqsubm(sequences$ref,  method = 'TRATE')
  )

  seqdistmc(
    channels = list(sequences$obs, sequences$cons, sequences$ref),
    method = 'OM',
    sm = substitution_costs,
    indel = 1,
    cweight = c(1, 1, 1)
  )
}

# These objects can be computationally expensive to create.
dist_breast <- compute_multichannel_om(seq_breast)
dist_lung   <- compute_multichannel_om(seq_lung)

# ------------------------------------------------------------------------------
# 6. Aggregate distance QC
# ------------------------------------------------------------------------------

summarise_distances <- function(d) {
  values <- as.vector(as.matrix(d))
  summary(values)
}

print(summarise_distances(dist_breast))
print(summarise_distances(dist_lung))


# saveRDS(dist_breast, file.path(derived_dir, 'distance_breast.rds'))
# saveRDS(dist_lung,   file.path(derived_dir, 'distance_lung.rds'))
