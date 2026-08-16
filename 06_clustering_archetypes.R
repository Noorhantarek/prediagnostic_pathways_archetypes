# ==============================================================================
# 06_clustering_archetypes.R
# Hierarchical clustering, internal validation, and archetype profiling
# ==============================================================================
# Public portfolio version.
# Selected cluster numbers and restricted-data results are intentionally not
# hard-coded in this public script.
# ===============================================================================

library(data.table)
library(cluster)
library(WeightedCluster)
library(TraMineR)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

derived_dir <- file.path('data', 'derived')

dist_breast <- readRDS(file.path(derived_dir, 'distance_breast.rds'))
dist_lung   <- readRDS(file.path(derived_dir, 'distance_lung.rds'))
seq_breast  <- readRDS(file.path(derived_dir, 'sequence_objects_breast.rds'))
seq_lung    <- readRDS(file.path(derived_dir, 'sequence_objects_lung.rds'))

# ------------------------------------------------------------------------------
# 2. Ward hierarchical clustering
# ------------------------------------------------------------------------------

fit_hierarchical <- function(distance_matrix) {
  d <- as.dist(distance_matrix)
  hc <- hclust(d, method = 'ward.D2')
  list(d = d, hc = hc)
}

breast_hc <- fit_hierarchical(dist_breast)
lung_hc   <- fit_hierarchical(dist_lung)

plot(breast_hc$hc, labels = FALSE, hang = -1,
     main = 'Breast: Ward.D2 hierarchical clustering', xlab = '', sub = '')
plot(lung_hc$hc, labels = FALSE, hang = -1,
     main = 'Lung: Ward.D2 hierarchical clustering', xlab = '', sub = '')

# ------------------------------------------------------------------------------
# 3. Candidate-k validation
# ------------------------------------------------------------------------------
# Cluster selection in the research combined statistical internal-validation
# criteria with clinical interpretability of the resulting pathway profiles.

evaluate_k <- function(hc, d, k_values = 2:10) {
  silhouette_results <- rbindlist(lapply(k_values, function(k) {
    cl <- cutree(hc, k = k)
    data.table(
      k = k,
      avg_silhouette = mean(silhouette(cl, d)[, 3]),
      smallest_cluster = min(table(cl)),
      largest_cluster = max(table(cl))
    )
  }))

  quality <- as.clustrange(hc, diss = d, ncluster = max(k_values))
  quality_dt <- as.data.table(quality$stats, keep.rownames = 'cluster_solution')

  list(silhouette = silhouette_results, quality = quality, quality_table = quality_dt)
}

breast_validation <- evaluate_k(breast_hc$hc, breast_hc$d)
lung_validation   <- evaluate_k(lung_hc$hc, lung_hc$d)

print(breast_validation$silhouette)
print(lung_validation$silhouette)

# WeightedCluster provides multiple internal criteria, including ASW, PBC,
# Hubert's Gamma and Calinski-Harabasz, used alongside clinical interpretation.
plot(breast_validation$quality, stat = c('ASW', 'PBC', 'HG', 'CH'), norm = 'zscore')
plot(lung_validation$quality,   stat = c('ASW', 'PBC', 'HG', 'CH'), norm = 'zscore')

# ------------------------------------------------------------------------------
# 4. Inspect candidate solutions
# ------------------------------------------------------------------------------

plot_candidate_solution <- function(sequences, hc, k) {
  cl <- cutree(hc, k = k)

  seqdplot(sequences$obs,  group = cl, border = NA, withlegend = 'right')
  seqdplot(sequences$cons, group = cl, border = NA, withlegend = 'right')
  seqdplot(sequences$ref,  group = cl, border = NA, withlegend = 'right')

  invisible(cl)
}

# Example usage inside the secure environment:
# plot_candidate_solution(seq_breast, breast_hc$hc, k = <candidate_k>)
# plot_candidate_solution(seq_lung,   lung_hc$hc,   k = <candidate_k>)

# ------------------------------------------------------------------------------
# 5. Numeric cluster fingerprints
# ------------------------------------------------------------------------------
# A compact state-composition table complements sequence plots when assessing
# whether candidate clusters represent distinct and clinically coherent pathways.

profile_clusters <- function(sample_wide, hc, k) {
  cl <- cutree(hc, k = k)

  to_long <- function(w, channel) {
    melt(
      data.table(cluster = cl, w[, 2:9]),
      id.vars = 'cluster',
      variable.name = 'slot',
      value.name = 'state'
    )[, channel := channel]
  }

  long <- rbindlist(list(
    to_long(sample_wide$obs, 'obs'),
    to_long(sample_wide$cons, 'cons'),
    to_long(sample_wide$ref, 'ref')
  ))

  out <- long[, .N, by = .(cluster, channel, state)]
  out[, pct := 100 * N / sum(N), by = .(cluster, channel)]
  out
}

# ------------------------------------------------------------------------------
# 6. Fit the final selected solution
# ------------------------------------------------------------------------------
# The selected k is deliberately supplied as an argument rather than hard-coded
# here because it is a result of the restricted-data analysis.

fit_selected_solution <- function(hc, patient_ids, k) {
  stopifnot(length(patient_ids) == length(hc$order))
  data.table(
    patid = as.character(patient_ids),
    cluster = cutree(hc, k = k)
  )
}

# Example:
# breast_clusters <- fit_selected_solution(breast_hc$hc, breast_patient_ids, k = selected_k)
# lung_clusters   <- fit_selected_solution(lung_hc$hc,   lung_patient_ids,   k = selected_k)

# Archetype names should be assigned only after reviewing the temporal and
# clinical characteristics of each final cluster.
