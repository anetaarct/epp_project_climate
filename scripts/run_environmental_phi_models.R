#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

dir.create("models_phi", recursive = TRUE, showWarnings = FALSE)
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

dat_env <- dat_model_phylo

dat_env$nest_box <- as.factor(dat_env$nest_box)
dat_env$marker_type <- as.factor(dat_env$marker_type)
dat_env$population_id <- as.factor(dat_env$population_id)
dat_env$species_nonphylo <- as.factor(dat_env$species_nonphylo)
dat_env$species_phylo <- as.factor(dat_env$species_phylo)

dat_env$abs_latitude_z <- z(dat_env$abs_latitude)
dat_env$distance_to_coast_z <- z(dat_env$distance_to_coast_km)

dat_env$tmp_cru_C_first_month_z <- z(dat_env$tmp_cru_C_first_month)
dat_env$pre_cru_mm_first_month_log_z <- z(log1p(dat_env$pre_cru_mm_first_month))
dat_env$dtr_cru_C_first_month_z <- z(dat_env$dtr_cru_C_first_month)
dat_env$ndvi_first_month_median_30km_z <- z(dat_env$ndvi_first_month_median_30km)

dat_env$tmp_short_anom_5yr_C_first_month_z <- z(dat_env$tmp_short_anom_5yr_C_first_month)
dat_env$pre_short_anom_5yr_mm_first_month_z <- z(dat_env$pre_short_anom_5yr_mm_first_month)

dat_env$tmp_long_anom_1961_1989_C_first_month_z <- z(dat_env$tmp_long_anom_1961_1989_C_first_month)
dat_env$pre_long_anom_1961_1989_mm_first_month_z <- z(dat_env$pre_long_anom_1961_1989_mm_first_month)

common_vars <- c(
  "n_epp_broods",
  "n_broods_sampled",
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

dat_env_absolute <- complete_model_data(dat_env, absolute_vars)
dat_env_short <- complete_model_data(dat_env, short_vars)
dat_env_long <- complete_model_data(dat_env, long_vars)

phi_model_data_checks <- data.frame(
  model = c(
    "absolute climate phi",
    "short anomalies phi",
    "long anomalies phi"
  ),
  records = c(
    nrow(dat_env_absolute),
    nrow(dat_env_short),
    nrow(dat_env_long)
  ),
  species = c(
    length(unique(dat_env_absolute$species_phylo)),
    length(unique(dat_env_short$species_phylo)),
    length(unique(dat_env_long$species_phylo))
  ),
  populations = c(
    length(unique(dat_env_absolute$population_id)),
    length(unique(dat_env_short$population_id)),
    length(unique(dat_env_long$population_id))
  )
)

write.csv(
  phi_model_data_checks,
  "models_phi/environmental_phi_model_data_checks.csv",
  row.names = FALSE
)

message("Environmental phi model data:")
print(phi_model_data_checks)

priors_env_phi <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 2), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "b", dpar = "phi")
)

message("Fitting phi model: absolute climate and NDVI...")
m_env_absolute_phi <- brm(
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat)),
    phi ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
  ),
  data = dat_env_absolute,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env_phi,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 241,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_absolute_phi, "models_phi/m_env_absolute_phi.rds")
capture.output(
  summary(m_env_absolute_phi),
  file = "models_phi/summary_m_env_absolute_phi.txt"
)

message("Fitting phi model: short-term anomalies and NDVI...")
m_env_short_phi <- brm(
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat)),
    phi ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
  ),
  data = dat_env_short,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env_phi,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 251,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_short_phi, "models_phi/m_env_short_phi.rds")
capture.output(
  summary(m_env_short_phi),
  file = "models_phi/summary_m_env_short_phi.txt"
)

message("Fitting phi model: long-term anomalies and NDVI...")
m_env_long_phi <- brm(
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat)),
    phi ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z
  ),
  data = dat_env_long,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env_phi,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 261,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_long_phi, "models_phi/m_env_long_phi.rds")
capture.output(
  summary(m_env_long_phi),
  file = "models_phi/summary_m_env_long_phi.txt"
)

message("Computing LOO for environmental phi models...")
loo_env_absolute_phi <- loo(m_env_absolute_phi)
loo_env_short_phi <- loo(m_env_short_phi)
loo_env_long_phi <- loo(m_env_long_phi)

loo_comparison_environmental_phi <- loo_compare(
  loo_env_absolute_phi,
  loo_env_short_phi,
  loo_env_long_phi
)

saveRDS(loo_env_absolute_phi, "models_phi/loo_env_absolute_phi.rds")
saveRDS(loo_env_short_phi, "models_phi/loo_env_short_phi.rds")
saveRDS(loo_env_long_phi, "models_phi/loo_env_long_phi.rds")
saveRDS(
  loo_comparison_environmental_phi,
  "models_phi/loo_comparison_environmental_phi.rds"
)
capture.output(
  loo_comparison_environmental_phi,
  file = "models_phi/loo_comparison_environmental_phi.txt"
)

capture.output(sessionInfo(), file = "logs/session_environmental_phi_models.txt")

message("Done. Environmental phi model outputs saved in models_phi/.")
