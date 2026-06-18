suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(posterior)
  library(loo)
})

out_dir <- "models/publication_environmental_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

chains <- 4
cores <- 4
iter <- 10000
warmup <- 5000

z <- function(x) as.numeric(scale(x))

message("Loading model-ready PCA and phylogenetic data...")

dat <- readRDS("data/processed/dat_model_phylo_pca.rds")
phylo_mat <- readRDS("data/processed/phylo_mat.rds")

dat_fit <- dat %>%
  mutate(
    population_id = factor(population_id),
    species_phylo = factor(species_phylo),
    species_nonphylo = factor(species_nonphylo)
  ) %>%
  filter(
    !is.na(n_epp_broods),
    !is.na(n_broods_sampled),
    n_broods_sampled > 0,
    n_epp_broods >= 0,
    n_epp_broods <= n_broods_sampled,
    !is.na(population_id),
    !is.na(species_nonphylo),
    !is.na(species_phylo),
    !is.na(PC1_thermal_geographic),
    !is.na(PC2_wet_productive),
    !is.na(tmp_short_anom_5yr_C_first_month),
    !is.na(pre_short_anom_5yr_mm_first_month),
    !is.na(tmp_long_anom_1961_1989_C_first_month),
    !is.na(pre_long_anom_1961_1989_mm_first_month),
    !is.na(abs_latitude),
    !is.na(distance_to_coast_km),
    as.character(species_phylo) %in% rownames(phylo_mat)
  ) %>%
  mutate(
    tmp_short_anom_5yr_C_first_month_z = z(tmp_short_anom_5yr_C_first_month),
    pre_short_anom_5yr_mm_first_month_z = z(pre_short_anom_5yr_mm_first_month),
    tmp_long_anom_1961_1989_C_first_month_z = z(tmp_long_anom_1961_1989_C_first_month),
    pre_long_anom_1961_1989_mm_first_month_z = z(pre_long_anom_1961_1989_mm_first_month),
    abs_latitude_z = z(abs_latitude),
    distance_to_coast_z = z(distance_to_coast_km)
  ) %>%
  droplevels()

phylo_mat_fit <- phylo_mat[
  levels(dat_fit$species_phylo),
  levels(dat_fit$species_phylo)
]

saveRDS(dat_fit, file.path(out_dir, "dat_publication_environmental_models.rds"))
saveRDS(phylo_mat_fit, file.path(out_dir, "phylo_mat_publication_environmental_models.rds"))

data_checks <- tibble(
  n_records = nrow(dat_fit),
  n_populations = n_distinct(dat_fit$population_id),
  n_species_phylo = n_distinct(dat_fit$species_phylo),
  n_species_nonphylo = n_distinct(dat_fit$species_nonphylo),
  total_broods = sum(dat_fit$n_broods_sampled),
  total_epp_broods = sum(dat_fit$n_epp_broods),
  raw_epp_rate = total_epp_broods / total_broods
)

write_csv(data_checks, file.path(out_dir, "publication_environmental_model_data_checks.csv"))
print(data_checks)

priors_mu <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "phi")
)

fit_model <- function(model_name, formula, seed) {
  message("Fitting model: ", model_name)

  fit <- brm(
    formula = formula,
    data = dat_fit,
    data2 = list(phylo_mat_fit = phylo_mat_fit),
    family = beta_binomial(link = "logit", link_phi = "log"),
    prior = priors_mu,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = seed,
    save_pars = save_pars(all = TRUE)
  )

  saveRDS(fit, file.path(out_dir, paste0(model_name, ".rds")))

  capture.output(
    summary(fit),
    file = file.path(out_dir, paste0("summary_", model_name, ".txt"))
  )

  fit
}

m_pub_pca_mu <- fit_model(
  "m_pub_pca_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    PC1_thermal_geographic +
    PC2_wet_productive +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  seed = 1801
)

m_pub_short_mu <- fit_model(
  "m_pub_short_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    tmp_short_anom_5yr_C_first_month_z +
    pre_short_anom_5yr_mm_first_month_z +
    abs_latitude_z +
    distance_to_coast_z +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  seed = 1802
)

m_pub_long_mu <- fit_model(
  "m_pub_long_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    tmp_long_anom_1961_1989_C_first_month_z +
    pre_long_anom_1961_1989_mm_first_month_z +
    abs_latitude_z +
    distance_to_coast_z +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  seed = 1803
)

diagnostic_summary <- function(fit, model_name) {
  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  sampler_diagnostics <- brms::nuts_params(fit)

  tibble(
    model = model_name,
    max_rhat = max(draw_diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(draw_diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_diagnostics$ess_tail, na.rm = TRUE),
    divergent_transitions = sampler_diagnostics %>%
      filter(Parameter == "divergent__") %>%
      summarise(n = sum(Value == 1)) %>%
      pull(n),
    treedepth_hits = sampler_diagnostics %>%
      filter(Parameter == "treedepth__") %>%
      summarise(n = sum(Value >= 15)) %>%
      pull(n)
  )
}

model_list <- list(
  m_pub_pca_mu = m_pub_pca_mu,
  m_pub_short_mu = m_pub_short_mu,
  m_pub_long_mu = m_pub_long_mu
)

diagnostics <- bind_rows(lapply(names(model_list), function(model_name) {
  diagnostic_summary(model_list[[model_name]], model_name)
}))

write_csv(diagnostics, file.path(out_dir, "diagnostics_publication_environmental_models.csv"))
print(diagnostics)

loo_list <- lapply(names(model_list), function(model_name) {
  message("Computing PSIS-LOO: ", model_name)
  loo_value <- loo(model_list[[model_name]])
  saveRDS(loo_value, file.path(out_dir, paste0("loo_", model_name, ".rds")))
  loo_value
})

names(loo_list) <- names(model_list)

loo_diagnostics <- bind_rows(lapply(names(loo_list), function(model_name) {
  tibble(
    model = model_name,
    max_pareto_k = max(loo_list[[model_name]]$diagnostics$pareto_k, na.rm = TRUE),
    n_pareto_k_gt_0.7 = sum(loo_list[[model_name]]$diagnostics$pareto_k > 0.7, na.rm = TRUE),
    n_pareto_k_gt_1 = sum(loo_list[[model_name]]$diagnostics$pareto_k > 1, na.rm = TRUE)
  )
}))

write_csv(loo_diagnostics, file.path(out_dir, "loo_diagnostics_publication_environmental_models.csv"))
print(loo_diagnostics)

loo_comparison <- loo_compare(loo_list) %>%
  as.data.frame() %>%
  rownames_to_column("model")

write_csv(loo_comparison, file.path(out_dir, "loo_compare_publication_environmental_models.csv"))

capture.output(
  loo_compare(loo_list),
  file = file.path(out_dir, "loo_compare_publication_environmental_models.txt")
)

fixed_effect_terms <- c(
  b_PC1_thermal_geographic = "PC1 thermal-geographic",
  b_PC2_wet_productive = "PC2 wet-productivity",
  b_tmp_short_anom_5yr_C_first_month_z = "Short temperature anomaly",
  b_pre_short_anom_5yr_mm_first_month_z = "Short precipitation anomaly",
  b_tmp_long_anom_1961_1989_C_first_month_z = "Long temperature anomaly",
  b_pre_long_anom_1961_1989_mm_first_month_z = "Long precipitation anomaly",
  b_abs_latitude_z = "Absolute latitude",
  b_distance_to_coast_z = "Distance to coast"
)

extract_fixef <- function(fit, model_name) {
  draws <- posterior::as_draws_df(fit)
  available_terms <- intersect(names(fixed_effect_terms), names(draws))

  bind_rows(lapply(available_terms, function(term) {
    x <- draws[[term]]
    tibble(
      model = model_name,
      term = term,
      predictor = fixed_effect_terms[[term]],
      estimate = median(x),
      lower = quantile(x, 0.025),
      upper = quantile(x, 0.975),
      prob_positive = mean(x > 0),
      prob_negative = mean(x < 0)
    )
  }))
}

fixed_effects <- bind_rows(lapply(names(model_list), function(model_name) {
  extract_fixef(model_list[[model_name]], model_name)
}))

write_csv(fixed_effects, file.path(out_dir, "fixed_effects_publication_environmental_models.csv"))

capture.output(
  sessionInfo(),
  file = "logs/session_publication_environmental_models.txt"
)

message("Done. Publication environmental models saved in ", out_dir, ".")
