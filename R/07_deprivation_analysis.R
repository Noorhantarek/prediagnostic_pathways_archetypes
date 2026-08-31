# ==============================================================================
# 07_deprivation_analysis.R
# Equity analysis: pathway archetypes and patient-level deprivation (IMD)
# ==============================================================================
# Public portfolio version.
# No source-data counts, p-values, or final archetype labels are hard-coded.
# ===============================================================================

library(data.table)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------------------------

derived_dir <- file.path('data', 'derived')

cohort <- fread(
  file.path(derived_dir, 'analytical_cohort_with_imd.csv'),
  colClasses = list(character = 'patid')
)
cluster_assignments <- fread(
  file.path(derived_dir, 'cluster_assignments.csv'),
  colClasses = list(character = 'patid')
)

# Expected fields:
# cohort: patid, cancer_type, imd_quintile, age_at_index, gender, index_date
# cluster_assignments: patid, cancer_type, cluster

profile <- merge(
  cluster_assignments,
  cohort,
  by = c('patid', 'cancer_type'),
  all.x = TRUE
)

# ------------------------------------------------------------------------------
# 2. Deprivation quintiles (Patient-level IMD) across archetypes clusters
# ------------------------------------------------------------------------------

imd_distribution <- profile[
  !is.na(imd_quintile),
  .N,
  by = .(cancer_type, imd_quintile, cluster)
]
imd_distribution[, pct := 100 * N / sum(N), by = .(cancer_type, imd_quintile)]

print(imd_distribution[order(cancer_type, imd_quintile, cluster)])

# ------------------------------------------------------------------------------
# 3. Chi-squared association and standardised residuals
# ------------------------------------------------------------------------------

analyse_imd_association <- function(dt) {
  tab <- dcast(
    dt[!is.na(imd_quintile), .N, by = .(cluster, imd_quintile)],
    cluster ~ imd_quintile,
    value.var = 'N',
    fill = 0
  )

  mat <- as.matrix(tab[, -1])
  rownames(mat) <- tab$cluster

  test <- chisq.test(mat)
  residuals <- as.data.table(as.table(test$stdres))
  setnames(residuals, c('cluster', 'imd_quintile', 'std_residual'))

  list(test = test, residuals = residuals)
}

results_by_cancer <- lapply(
  split(profile, profile$cancer_type),
  analyse_imd_association
)

# ------------------------------------------------------------------------------
# 4. Ordinal trend of IMD quintiles across clusters
# ------------------------------------------------------------------------------
# For each IMD quintile, test whether its share changes monotonically across archetypes
# clusters. Interpretation should consider effect size and pattern,
# not the p-value alone.

trend_test_by_cluster <- function(dt) {
  clusters <- sort(unique(dt$cluster))

  rbindlist(lapply(clusters, function(k) {
    x <- dt[
      !is.na(imd_quintile),
      .(in_cluster = sum(cluster == k), total = .N),
      by = imd_quintile
    ][order(imd_quintile)]

    test <- prop.trend.test(x$in_cluster, x$total)

    data.table(
      cluster = k,
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  }))
}

trend_results <- rbindlist(
  lapply(split(profile, profile$cancer_type), trend_test_by_cluster),
  idcol = 'cancer_type'
)
print(trend_results)

# ------------------------------------------------------------------------------
# 5. Visualise quintiles distribution across clusters
# ------------------------------------------------------------------------------

plot_imd_trends <- function(dt, cancer_label) {
  ggplot(
    dt[cancer_type == cancer_label],
    aes(x = imd_quintile, y = pct, colour = factor(cluster), group = cluster)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c('Q1\nleast deprived', 'Q2', 'Q3', 'Q4', 'Q5\nmost deprived')
    ) +
    labs(
      title = paste0(cancer_label, ': Deprivation quintiles distribution across archetypes'),
      x = NULL,
      y = '% of patients within quintile',
      colour = 'Cluster'
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = 'bottom')
}

print(plot_imd_trends(imd_distribution, 'breast cancer'))
print(plot_imd_trends(imd_distribution, 'lung cancer'))

# ------------------------------------------------------------------------------
# 6. Residual heatmap
# ------------------------------------------------------------------------------

plot_residuals <- function(residual_dt, title) {
  ggplot(residual_dt,
         aes(x = imd_quintile, y = factor(cluster), fill = std_residual)) +
    geom_tile() +
    geom_text(aes(label = round(std_residual, 1)), size = 3) +
    scale_fill_gradient2(midpoint = 0) +
    labs(title = title,
         subtitle = 'Positive = over-represented; negative = under-represented',
         x = 'IMD quintile', y = 'Cluster', fill = 'Std. residual') +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank())
}

# Example:
# print(plot_residuals(results_by_cancer[['breast cancer']]$residuals,
#                      'Breast cancer: standardised residuals'))
