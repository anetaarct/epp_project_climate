library(dplyr)
library(readr)
library(tibble)
library(knitr)

model_dir <- "models/sensitivity_models"

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}

data_checks_file <- first_existing(c(
  file.path(model_dir, "sensitivity_model_data_checks.csv")
))


summary_files <- list(
  first_month_mu = first_existing(c(file.path(model_dir, "summary_m_sensitivity_first_month_mu.txt"))),
  first_month_mu_phi = first_existing(c(file.path(model_dir, "summary_m_sensitivity_first_month_mu_phi.txt")))
)

# m_sensitivity_first_month_mu <- brms::brm(
#   n_epp_broods | trials(n_broods_sampled) ~
#     tmp_cru_C_first_month_z +
#     pre_cru_mm_first_month_log_z +
#     dtr_cru_C_first_month_z +
#     ndvi_first_month_median_30km_z +
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
#   seed = 1501,
#   save_pars = brms::save_pars(all = TRUE)
# )

# m_sensitivity_first_month_mu_phi <- brms::brm(
#   brms::bf(
#     n_epp_broods | trials(n_broods_sampled) ~
#       tmp_cru_C_first_month_z +
#       pre_cru_mm_first_month_log_z +
#       dtr_cru_C_first_month_z +
#       ndvi_first_month_median_30km_z +
#       abs_latitude_z +
#       distance_to_coast_z +
#       nest_box +
#       marker_type +
#       (1 | population_id) +
#       (1 | species_nonphylo) +
#       (1 | gr(species_phylo, cov = phylo_mat_fit)),
# 
#     phi ~
#       tmp_cru_C_first_month_z +
#       pre_cru_mm_first_month_log_z +
#       dtr_cru_C_first_month_z +
#       ndvi_first_month_median_30km_z +
#       abs_latitude_z +
#       distance_to_coast_z +
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
#   seed = 1502,
#   save_pars = brms::save_pars(all = TRUE)
# )

if (!is.na(data_checks_file)) {
  read_csv(data_checks_file, show_col_types = FALSE) %>%
    kable()
} else {
  cat("Sensitivity model data checks are not available yet. Run the sensitivity model script first.")
}

if (file.exists("figures/sensitivity_environmental_effects.png")) {
  knitr::include_graphics("figures/sensitivity_environmental_effects.png")
} else {
  cat("Sensitivity environmental effect panel is not available yet. Render the Figures chapter to create it.")
}

if (file.exists("figures/sensitivity_fixed_effects_forest_plot.png")) {
  knitr::include_graphics("figures/sensitivity_fixed_effects_forest_plot.png")
} else {
  cat("Sensitivity fixed-effect plot is not available yet. Render the Figures chapter to create it.")
}

summary_labels <- c(
  first_month_mu = "Mean-only sensitivity model",
  first_month_mu_phi = "Mean-scale sensitivity model"
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
  cat("Model summaries are not available yet. Run the sensitivity model script first.")
}
