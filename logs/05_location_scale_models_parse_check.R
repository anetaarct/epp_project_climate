library(dplyr)
library(readr)
library(tibble)
library(knitr)

model_dir <- "models/main_environmental_models"

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}

data_checks_file <- first_existing(c(
  file.path(model_dir, "main_environmental_model_data_checks.csv")
))

fixed_effects_file <- first_existing(c(
  file.path(model_dir, "fixed_effects_main_environmental_models.csv")
))

summary_files <- list(
  pca_mu = first_existing(c(file.path(model_dir, "summary_m_main_pca_mu.txt"))),
  pca_mu_phi = first_existing(c(file.path(model_dir, "summary_m_main_pca_mu_phi.txt"))),
  short_mu = first_existing(c(file.path(model_dir, "summary_m_main_short_mu.txt"))),
  short_mu_phi = first_existing(c(file.path(model_dir, "summary_m_main_short_mu_phi.txt"))),
  long_mu = first_existing(c(file.path(model_dir, "summary_m_main_long_mu.txt"))),
  long_mu_phi = first_existing(c(file.path(model_dir, "summary_m_main_long_mu_phi.txt")))
)

# m_main_pca_mu <- brms::brm(
#   n_epp_broods | trials(n_broods_sampled) ~
#     PC1_thermal_geographic +
#     PC2_wet_productive +
#     nest_box +
#     marker_type +
#     (1 | population_id) +
#     (1 | species_nonphylo) +
#     (1 | gr(species_phylo, cov = phylo_mat_fit)),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1401,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_main_pca_mu_phi <- brms::brm(
#   brms::bf(
#     n_epp_broods | trials(n_broods_sampled) ~
#       PC1_thermal_geographic +
#       PC2_wet_productive +
#       nest_box +
#       marker_type +
#       (1 | population_id) +
#       (1 | species_nonphylo) +
#       (1 | gr(species_phylo, cov = phylo_mat_fit)),
# 
#     phi ~
#       PC1_thermal_geographic +
#       PC2_wet_productive +
#       nest_box +
#       marker_type
#   ),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu_phi,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1402,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_main_short_mu <- brms::brm(
#   n_epp_broods | trials(n_broods_sampled) ~
#     tmp_short_anom_5yr_C_first_month_z +
#     pre_short_anom_5yr_mm_first_month_z +
#     abs_latitude_z +
#     distance_to_coast_z +
#     nest_box +
#     marker_type +
#     (1 | population_id) +
#     (1 | species_nonphylo) +
#     (1 | gr(species_phylo, cov = phylo_mat_fit)),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1403,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_main_short_mu_phi <- brms::brm(
#   brms::bf(
#     n_epp_broods | trials(n_broods_sampled) ~
#       tmp_short_anom_5yr_C_first_month_z +
#       pre_short_anom_5yr_mm_first_month_z +
#       abs_latitude_z +
#       distance_to_coast_z +
#       nest_box +
#       marker_type +
#       (1 | population_id) +
#       (1 | species_nonphylo) +
#       (1 | gr(species_phylo, cov = phylo_mat_fit)),
# 
#     phi ~
#       tmp_short_anom_5yr_C_first_month_z +
#       pre_short_anom_5yr_mm_first_month_z +
#       nest_box +
#       marker_type
#   ),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu_phi,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1404,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_main_long_mu <- brms::brm(
#   n_epp_broods | trials(n_broods_sampled) ~
#     tmp_long_anom_1961_1989_C_first_month_z +
#     pre_long_anom_1961_1989_mm_first_month_z +
#     abs_latitude_z +
#     distance_to_coast_z +
#     nest_box +
#     marker_type +
#     (1 | population_id) +
#     (1 | species_nonphylo) +
#     (1 | gr(species_phylo, cov = phylo_mat_fit)),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1405,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_main_long_mu_phi <- brms::brm(
#   brms::bf(
#     n_epp_broods | trials(n_broods_sampled) ~
#       tmp_long_anom_1961_1989_C_first_month_z +
#       pre_long_anom_1961_1989_mm_first_month_z +
#       abs_latitude_z +
#       distance_to_coast_z +
#       nest_box +
#       marker_type +
#       (1 | population_id) +
#       (1 | species_nonphylo) +
#       (1 | gr(species_phylo, cov = phylo_mat_fit)),
# 
#     phi ~
#       tmp_long_anom_1961_1989_C_first_month_z +
#       pre_long_anom_1961_1989_mm_first_month_z +
#       nest_box +
#       marker_type
#   ),
#   data = dat_fit,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log"),
#   prior = priors_mu_phi,
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(adapt_delta = 0.99, max_treedepth = 15),
#   seed = 1406,
#   save_pars = brms::save_pars(all = TRUE)
# )

if (!is.na(data_checks_file)) {
  read_csv(data_checks_file, show_col_types = FALSE) %>%
    kable()
} else {
  cat("Main environmental model data checks are not available yet. Run the main environmental model script first.")
}

if (!is.na(fixed_effects_file)) {
  read_csv(fixed_effects_file, show_col_types = FALSE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    kable()
} else {
  cat("Fixed-effect summaries for the main environmental models are not available yet.")
}

summary_labels <- c(
  pca_mu = "Environmental PCA mean-only model",
  pca_mu_phi = "Environmental PCA mean-scale model",
  short_mu = "Short-anomaly mean-only model",
  short_mu_phi = "Short-anomaly mean-scale model",
  long_mu = "Long-anomaly mean-only model",
  long_mu_phi = "Long-anomaly mean-scale model"
)

for (summary_name in names(summary_files)) {
  summary_file <- summary_files[[summary_name]]

  if (!is.na(summary_file)) {
    cat("### ", summary_labels[[summary_name]], "\n\n", sep = "")
    cat("<details><summary>Show full model summary</summary>\n\n```\n")
    cat(readLines(summary_file), sep = "\n")
    cat("\n```\n</details>\n")
  }
}

if (all(is.na(unlist(summary_files)))) {
  cat("Model summaries are not available yet. Run the main environmental model script first.")
}
