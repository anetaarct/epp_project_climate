#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(posterior)
})

dir.create("models_scale", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

message("Loading model data...")
dat_model_phylo <- readRDS("dat_model_phylo.rds")

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

dat_scale <- dat_model_phylo

dat_scale$species_nonphylo <- as.factor(dat_scale$species_nonphylo)
dat_scale$population_id <- as.factor(dat_scale$population_id)

dat_scale$year_of_study_start <- as.numeric(dat_scale$year_of_study_start)
dat_scale$year_of_study_end <- as.numeric(dat_scale$year_of_study_end)
dat_scale$study_year <- ifelse(
  !is.na(dat_scale$year_of_study_end),
  rowMeans(
    cbind(dat_scale$year_of_study_start, dat_scale$year_of_study_end),
    na.rm = TRUE
  ),
  dat_scale$year_of_study_start
)
dat_scale$study_year_factor <- as.factor(round(dat_scale$study_year))

scale_model_vars <- c(
  "n_epp_broods",
  "n_broods_sampled",
  "species_nonphylo",
  "population_id",
  "study_year_factor"
)

dat_scale_model <- complete_model_data(
  dat_scale,
  scale_model_vars
)

scale_data_support <- data.frame(
  quantity = c(
    "Records",
    "Species",
    "Populations",
    "Study years"
  ),
  value = c(
    nrow(dat_scale_model),
    dplyr::n_distinct(dat_scale_model$species_nonphylo),
    dplyr::n_distinct(dat_scale_model$population_id),
    dplyr::n_distinct(dat_scale_model$study_year_factor)
  )
)

write.csv(
  scale_data_support,
  "models_scale/scale_model_data_support.csv",
  row.names = FALSE
)

message("Scale model data:")
print(scale_data_support)

priors_scale <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 2), class = "Intercept", dpar = "phi"),
  prior(exponential(1), class = "sd", dpar = "phi")
)

message("Fitting beta-binomial scale model...")
m_scale_species_population_year <- brm(
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      1 +
      (1 | mu_species | species_nonphylo) +
      (1 | mu_population | population_id),
    phi ~
      1 +
      (1 | phi_species | species_nonphylo) +
      (1 | phi_population | population_id) +
      (1 | phi_year | study_year_factor)
  ),
  data = dat_scale_model,
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_scale,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 301,
  save_pars = save_pars(all = TRUE)
)

saveRDS(
  m_scale_species_population_year,
  "models_scale/m_scale_species_population_year.rds"
)

capture.output(
  summary(m_scale_species_population_year),
  file = "models_scale/summary_m_scale_species_population_year.txt"
)

posterior_summary <- as.data.frame(
  posterior::summarise_draws(
    posterior::as_draws_df(m_scale_species_population_year),
    mean,
    sd,
    ~quantile(.x, probs = 0.025),
    ~quantile(.x, probs = 0.975)
  )
)

names(posterior_summary) <- c(
  "parameter",
  "estimate",
  "est_error",
  "lower_95",
  "upper_95"
)

scale_parameter_components <- posterior_summary %>%
  filter(
    grepl("^sd_.*__phi_Intercept$", parameter)
  ) %>%
  mutate(
    scale_component = case_when(
      grepl("species_nonphylo", parameter) ~ "Species",
      grepl("population_id", parameter) ~ "Population",
      grepl("study_year_factor", parameter) ~ "Study year",
      TRUE ~ parameter
    )
  ) %>%
  select(
    scale_component,
    estimate,
    est_error,
    lower_95,
    upper_95
  )

write.csv(
  scale_parameter_components,
  "models_scale/scale_parameter_components.csv",
  row.names = FALSE
)

capture.output(sessionInfo(), file = "logs/session_scale_parameter_model.txt")

message("Done. Scale model outputs saved in models_scale/.")
