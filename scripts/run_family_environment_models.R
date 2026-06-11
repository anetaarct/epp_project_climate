#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
  library(dplyr)
})

dir.create("models_family", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

message("Loading phylogeny-ready model data...")
dat_model_phylo <- readRDS("dat_model_phylo.rds")
phylo_mat <- readRDS("phylo_mat.rds")

z <- function(x) {
  as.numeric(scale(x))
}

complete_model_data <- function(data, vars) {
  data <- data[
    !is.na(data$n_epp_broods) &
      !is.na(data$n_broods_sampled) &
      data$n_broods_sampled > 0 &
      data$n_epp_broods <= data$n_broods_sampled,
  ]

  data <- data[complete.cases(data[, vars]), ]
  droplevels(data)
}

family_epp_scale <- dat_model_phylo %>%
  filter(
    !is.na(family),
    !is.na(species_phylo),
    !is.na(n_epp_broods),
    !is.na(n_broods_sampled),
    n_broods_sampled > 0
  ) %>%
  mutate(epp_rate = n_epp_broods / n_broods_sampled) %>%
  group_by(family, species_phylo) %>%
  summarise(
    species_mean_epp = weighted.mean(
      epp_rate,
      n_broods_sampled,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  group_by(family) %>%
  summarise(
    family_n_species = n_distinct(species_phylo),
    family_mean_epp = mean(species_mean_epp, na.rm = TRUE),
    family_scale = sd(species_mean_epp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(family_n_species >= 3) %>%
  mutate(
    family_mean_epp_z = z(family_mean_epp),
    family_scale_z = z(family_scale)
  )

dat_family_env <- dat_model_phylo %>%
  left_join(
    family_epp_scale,
    by = "family"
  )

dat_family_env$nest_box <- as.factor(dat_family_env$nest_box)
dat_family_env$marker_type <- as.factor(dat_family_env$marker_type)
dat_family_env$population_id <- as.factor(dat_family_env$population_id)
dat_family_env$species_nonphylo <- as.factor(dat_family_env$species_nonphylo)
dat_family_env$species_phylo <- as.factor(dat_family_env$species_phylo)

dat_family_env$abs_latitude_z <- z(dat_family_env$abs_latitude)
dat_family_env$distance_to_coast_z <- z(dat_family_env$distance_to_coast_km)

dat_family_env$tmp_cru_C_first_month_z <- z(dat_family_env$tmp_cru_C_first_month)
dat_family_env$pre_cru_mm_first_month_log_z <- z(log1p(dat_family_env$pre_cru_mm_first_month))
dat_family_env$dtr_cru_C_first_month_z <- z(dat_family_env$dtr_cru_C_first_month)
dat_family_env$ndvi_first_month_median_30km_z <- z(dat_family_env$ndvi_first_month_median_30km)

dat_family_env$tmp_short_anom_5yr_C_first_month_z <- z(dat_family_env$tmp_short_anom_5yr_C_first_month)
dat_family_env$pre_short_anom_5yr_mm_first_month_z <- z(dat_family_env$pre_short_anom_5yr_mm_first_month)

dat_family_env$tmp_long_anom_1961_1989_C_first_month_z <- z(dat_family_env$tmp_long_anom_1961_1989_C_first_month)
dat_family_env$pre_long_anom_1961_1989_mm_first_month_z <- z(dat_family_env$pre_long_anom_1961_1989_mm_first_month)

common_vars <- c(
  "n_epp_broods",
  "n_broods_sampled",
  "family_mean_epp_z",
  "family_scale_z",
  "abs_latitude_z",
  "distance_to_coast_z",
  "nest_box",
  "marker_type",
  "population_id",
  "species_nonphylo",
  "species_phylo"
)

absolute_vars <- c(
  common_vars,
  "tmp_cru_C_first_month_z",
  "pre_cru_mm_first_month_log_z",
  "dtr_cru_C_first_month_z",
  "ndvi_first_month_median_30km_z"
)

short_vars <- c(
  common_vars,
  "tmp_short_anom_5yr_C_first_month_z",
  "pre_short_anom_5yr_mm_first_month_z",
  "dtr_cru_C_first_month_z",
  "ndvi_first_month_median_30km_z"
)

long_vars <- c(
  common_vars,
  "tmp_long_anom_1961_1989_C_first_month_z",
  "pre_long_anom_1961_1989_mm_first_month_z",
  "dtr_cru_C_first_month_z",
  "ndvi_first_month_median_30km_z"
)

dat_family_absolute <- complete_model_data(dat_family_env, absolute_vars)
dat_family_short <- complete_model_data(dat_family_env, short_vars)
dat_family_long <- complete_model_data(dat_family_env, long_vars)

family_model_data_checks <- data.frame(
  model = c(
    "absolute climate",
    "short anomalies",
    "long anomalies"
  ),
  records = c(
    nrow(dat_family_absolute),
    nrow(dat_family_short),
    nrow(dat_family_long)
  ),
  families = c(
    n_distinct(dat_family_absolute$family),
    n_distinct(dat_family_short$family),
    n_distinct(dat_family_long$family)
  ),
  species = c(
    n_distinct(dat_family_absolute$species_phylo),
    n_distinct(dat_family_short$species_phylo),
    n_distinct(dat_family_long$species_phylo)
  )
)

write.csv(
  family_model_data_checks,
  "models_family/family_environment_model_data_checks.csv",
  row.names = FALSE
)

message("Family environment model data:")
print(family_model_data_checks)

priors_family_env <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd")
)

message("Fitting family sensitivity model: absolute climate...")
m_family_absolute <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    (
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_scale_z +
    (
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_mean_epp_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_family_absolute,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_family_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 271,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_family_absolute, "models_family/m_family_absolute.rds")
capture.output(
  summary(m_family_absolute),
  file = "models_family/summary_m_family_absolute.txt"
)

message("Fitting family sensitivity model: short anomalies...")
m_family_short <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    (
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_scale_z +
    (
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_mean_epp_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_family_short,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_family_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 272,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_family_short, "models_family/m_family_short.rds")
capture.output(
  summary(m_family_short),
  file = "models_family/summary_m_family_short.txt"
)

message("Fitting family sensitivity model: long anomalies...")
m_family_long <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    (
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_scale_z +
    (
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
    ) * family_mean_epp_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_family_long,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_family_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 273,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_family_long, "models_family/m_family_long.rds")
capture.output(
  summary(m_family_long),
  file = "models_family/summary_m_family_long.txt"
)

message("Computing LOO for family sensitivity models...")
loo_family_absolute <- loo(m_family_absolute)
loo_family_short <- loo(m_family_short)
loo_family_long <- loo(m_family_long)

loo_comparison_family_environment <- loo_compare(
  loo_family_absolute,
  loo_family_short,
  loo_family_long
)

saveRDS(loo_family_absolute, "models_family/loo_family_absolute.rds")
saveRDS(loo_family_short, "models_family/loo_family_short.rds")
saveRDS(loo_family_long, "models_family/loo_family_long.rds")
saveRDS(
  loo_comparison_family_environment,
  "models_family/loo_comparison_family_environment.rds"
)
capture.output(
  loo_comparison_family_environment,
  file = "models_family/loo_comparison_family_environment.txt"
)

capture.output(sessionInfo(), file = "logs/session_family_environment_models.txt")

message("Done. Family environment model outputs saved in models_family/.")
