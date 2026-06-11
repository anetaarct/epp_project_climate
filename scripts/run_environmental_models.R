#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
  library(performance)
})

dir.create("models", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

message("Loading phylogeny-ready model data...")
dat_model_phylo <- readRDS("data/processed/dat_model_phylo.rds")
phylo_mat <- readRDS("data/processed/phylo_mat.rds")

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

model_data_checks <- data.frame(
  model = c(
    "absolute climate",
    "short anomalies",
    "long anomalies"
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
  model_data_checks,
  "models/environmental_model_data_checks.csv",
  row.names = FALSE
)

message("Environmental model data:")
print(model_data_checks)

priors_env <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd")
)

message("Fitting environmental model: absolute climate and NDVI...")
m_env_absolute <- brm(
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
  data = dat_env_absolute,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 201,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_absolute, "models/m_env_absolute.rds")
capture.output(summary(m_env_absolute), file = "models/summary_m_env_absolute.txt")

message("Fitting environmental model: short-term anomalies and NDVI...")
m_env_short <- brm(
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
  data = dat_env_short,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 211,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_short, "models/m_env_short.rds")
capture.output(summary(m_env_short), file = "models/summary_m_env_short.txt")

message("Fitting environmental model: long-term anomalies and NDVI...")
m_env_long <- brm(
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
  data = dat_env_long,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_env,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 221,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_env_long, "models/m_env_long.rds")
capture.output(summary(m_env_long), file = "models/summary_m_env_long.txt")

message("Computing LOO for environmental models...")
loo_env_absolute <- loo(m_env_absolute)
loo_env_short <- loo(m_env_short)
loo_env_long <- loo(m_env_long)

loo_comparison_environmental <- loo_compare(
  loo_env_absolute,
  loo_env_short,
  loo_env_long
)

saveRDS(loo_env_absolute, "models/loo_env_absolute.rds")
saveRDS(loo_env_short, "models/loo_env_short.rds")
saveRDS(loo_env_long, "models/loo_env_long.rds")
saveRDS(
  loo_comparison_environmental,
  "models/loo_comparison_environmental.rds"
)
capture.output(
  loo_comparison_environmental,
  file = "models/loo_comparison_environmental.txt"
)

capture.output(sessionInfo(), file = "logs/session_environmental_models.txt")

message("Done. Environmental model outputs saved in models/.")
