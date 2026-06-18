suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(posterior)
  library(loo)
})

out_dir <- "models/main_environmental_models"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

chains <- 4
cores <- 4
iter <- 10000
warmup <- 5000

z <- function(x) as.numeric(scale(x))

message("Loading model-ready phylogenetic data...")

dat <- readRDS("data/processed/dat_model_phylo.rds")
phylo_mat <- readRDS("data/processed/phylo_mat.rds")

dat_fit <- dat %>%
  mutate(
    population_id = factor(population_id),
    species_phylo = factor(species_phylo),
    species_nonphylo = factor(species_nonphylo),
    nest_box = factor(nest_box),
    marker_type = factor(marker_type)
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
    !is.na(nest_box),
    !is.na(marker_type),
    !is.na(tmp_cru_C_first_month),
    !is.na(pre_cru_mm_first_month),
    !is.na(dtr_cru_C_first_month),
    !is.na(ndvi_first_month_median_30km),
    !is.na(abs_latitude),
    !is.na(distance_to_coast_km),
    !is.na(tmp_short_anom_5yr_C_first_month),
    !is.na(pre_short_anom_5yr_mm_first_month),
    !is.na(tmp_long_anom_1961_1989_C_first_month),
    !is.na(pre_long_anom_1961_1989_mm_first_month),
    as.character(species_phylo) %in% rownames(phylo_mat)
  ) %>%
  mutate(
    tmp_cru_C_first_month_z = z(tmp_cru_C_first_month),
    pre_cru_mm_first_month_log_z = z(log1p(pre_cru_mm_first_month)),
    dtr_cru_C_first_month_z = z(dtr_cru_C_first_month),
    ndvi_first_month_median_30km_z = z(ndvi_first_month_median_30km),
    abs_latitude_z = z(abs_latitude),
    distance_to_coast_z = z(distance_to_coast_km),
    tmp_short_anom_5yr_C_first_month_z = z(tmp_short_anom_5yr_C_first_month),
    pre_short_anom_5yr_mm_first_month_z = z(pre_short_anom_5yr_mm_first_month),
    tmp_long_anom_1961_1989_C_first_month_z = z(tmp_long_anom_1961_1989_C_first_month),
    pre_long_anom_1961_1989_mm_first_month_z = z(pre_long_anom_1961_1989_mm_first_month)
  ) %>%
  droplevels()

pca_vars <- dat_fit %>%
  transmute(
    tmp_z = tmp_cru_C_first_month_z,
    pre_log_z = pre_cru_mm_first_month_log_z,
    dtr_z = dtr_cru_C_first_month_z,
    ndvi_z = ndvi_first_month_median_30km_z,
    abs_latitude_z = abs_latitude_z,
    distance_to_coast_z = distance_to_coast_z
  )

pca <- prcomp(pca_vars, center = FALSE, scale. = FALSE)

pc1_direction <- ifelse(pca$rotation["tmp_z", "PC1"] < 0, -1, 1)
pc2_direction <- ifelse(sum(pca$rotation[c("pre_log_z", "ndvi_z"), "PC2"]) < 0, -1, 1)

pca$x[, "PC1"] <- pc1_direction * pca$x[, "PC1"]
pca$rotation[, "PC1"] <- pc1_direction * pca$rotation[, "PC1"]

pca$x[, "PC2"] <- pc2_direction * pca$x[, "PC2"]
pca$rotation[, "PC2"] <- pc2_direction * pca$rotation[, "PC2"]

dat_fit <- dat_fit %>%
  bind_cols(
    as_tibble(pca$x[, 1:2]) %>%
      transmute(
        PC1_thermal_geographic = PC1,
        PC2_wet_productive = PC2
      )
  ) %>%
  droplevels()

phylo_mat_fit <- phylo_mat[
  levels(dat_fit$species_phylo),
  levels(dat_fit$species_phylo)
]

saveRDS(dat_fit, file.path(out_dir, "dat_main_environmental_models.rds"))
saveRDS(phylo_mat_fit, file.path(out_dir, "phylo_mat_main_environmental_models.rds"))

data_checks <- tibble(
  n_records = nrow(dat_fit),
  n_populations = n_distinct(dat_fit$population_id),
  n_species_phylo = n_distinct(dat_fit$species_phylo),
  n_species_nonphylo = n_distinct(dat_fit$species_nonphylo),
  total_broods = sum(dat_fit$n_broods_sampled),
  total_epp_broods = sum(dat_fit$n_epp_broods),
  raw_epp_rate = total_epp_broods / total_broods
)

write_csv(data_checks, file.path(out_dir, "main_environmental_model_data_checks.csv"))
print(data_checks)

priors_mu <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "phi")
)

priors_mu_phi <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(normal(0, 1), class = "b", dpar = "phi")
)

fit_model <- function(model_name, formula, priors, seed) {
  message("Fitting model: ", model_name)

  fit <- brm(
    formula = formula,
    data = dat_fit,
    data2 = list(phylo_mat_fit = phylo_mat_fit),
    family = beta_binomial(link = "logit", link_phi = "log"),
    prior = priors,
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

m_main_pca_mu <- fit_model(
  "m_main_pca_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    PC1_thermal_geographic +
    PC2_wet_productive +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1401
)

m_main_pca_mu_phi <- fit_model(
  "m_main_pca_mu_phi",
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      PC1_thermal_geographic +
      PC2_wet_productive +
      nest_box +
      marker_type +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat_fit)),
    phi ~
      PC1_thermal_geographic +
      PC2_wet_productive +
      nest_box +
      marker_type
  ),
  priors_mu_phi,
  seed = 1402
)

m_main_short_mu <- fit_model(
  "m_main_short_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    tmp_short_anom_5yr_C_first_month_z +
    pre_short_anom_5yr_mm_first_month_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1403
)

m_main_short_mu_phi <- fit_model(
  "m_main_short_mu_phi",
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat_fit)),
    phi ~
      tmp_short_anom_5yr_C_first_month_z +
      pre_short_anom_5yr_mm_first_month_z +
      nest_box +
      marker_type
  ),
  priors_mu_phi,
  seed = 1404
)

m_main_long_mu <- fit_model(
  "m_main_long_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    tmp_long_anom_1961_1989_C_first_month_z +
    pre_long_anom_1961_1989_mm_first_month_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1405
)

m_main_long_mu_phi <- fit_model(
  "m_main_long_mu_phi",
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat_fit)),
    phi ~
      tmp_long_anom_1961_1989_C_first_month_z +
      pre_long_anom_1961_1989_mm_first_month_z +
      nest_box +
      marker_type
  ),
  priors_mu_phi,
  seed = 1406
)

diagnostic_summary <- function(fit, model_name) {
  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  sampler_diagnostics <- brms::nuts_params(fit)

  divergent_transitions <- sampler_diagnostics %>%
    filter(Parameter == "divergent__") %>%
    summarise(n = sum(Value == 1)) %>%
    pull(n)

  treedepth_hits <- sampler_diagnostics %>%
    filter(Parameter == "treedepth__") %>%
    summarise(n = sum(Value >= 15)) %>%
    pull(n)

  tibble(
    model = model_name,
    max_rhat = max(draw_diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(draw_diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_diagnostics$ess_tail, na.rm = TRUE),
    divergent_transitions = divergent_transitions,
    treedepth_hits = treedepth_hits
  )
}

model_list <- list(
  m_main_pca_mu = m_main_pca_mu,
  m_main_pca_mu_phi = m_main_pca_mu_phi,
  m_main_short_mu = m_main_short_mu,
  m_main_short_mu_phi = m_main_short_mu_phi,
  m_main_long_mu = m_main_long_mu,
  m_main_long_mu_phi = m_main_long_mu_phi
)

diagnostics <- bind_rows(lapply(names(model_list), function(model_name) {
  diagnostic_summary(model_list[[model_name]], model_name)
}))

write_csv(diagnostics, file.path(out_dir, "diagnostics_main_environmental_models.csv"))
print(diagnostics)

loo_list <- lapply(names(model_list), function(model_name) {
  message("Computing LOO: ", model_name)
  loo_value <- loo(model_list[[model_name]], moment_match = TRUE)
  saveRDS(loo_value, file.path(out_dir, paste0("loo_", model_name, ".rds")))
  loo_value
})

names(loo_list) <- names(model_list)

loo_pairs <- list(
  "Environmental PCA" = c("m_main_pca_mu", "m_main_pca_mu_phi"),
  "Short anomaly" = c("m_main_short_mu", "m_main_short_mu_phi"),
  "Long anomaly" = c("m_main_long_mu", "m_main_long_mu_phi")
)

loo_comparison <- bind_rows(lapply(names(loo_pairs), function(hypothesis_group) {
  pair_names <- loo_pairs[[hypothesis_group]]
  pair_comparison <- loo_compare(loo_list[pair_names])

  as.data.frame(pair_comparison) %>%
    rownames_to_column("model") %>%
    mutate(
      hypothesis_group = hypothesis_group,
      model_version = ifelse(grepl("_mu_phi$", model), "mean-scale", "mean-only"),
      .before = model
    )
}))

write_csv(
  loo_comparison,
  file.path(out_dir, "loo_compare_main_environmental_models.csv")
)

capture.output(
  lapply(names(loo_pairs), function(hypothesis_group) {
    cat("\n", hypothesis_group, "\n", sep = "")
    print(loo_compare(loo_list[loo_pairs[[hypothesis_group]]]))
  }),
  file = file.path(out_dir, "loo_compare_main_environmental_models.txt")
)

fixed_effect_terms <- c(
  b_PC1_thermal_geographic = "PC1 thermal-geographic",
  b_PC2_wet_productive = "PC2 wet-productivity",
  b_tmp_short_anom_5yr_C_first_month_z = "Short temperature anomaly",
  b_pre_short_anom_5yr_mm_first_month_z = "Short precipitation anomaly",
  b_tmp_long_anom_1961_1989_C_first_month_z = "Long temperature anomaly",
  b_pre_long_anom_1961_1989_mm_first_month_z = "Long precipitation anomaly",
  b_abs_latitude_z = "Absolute latitude",
  b_distance_to_coast_z = "Distance to coast",
  b_nest_boxyes = "Nest box",
  b_marker_typemicrosatellite = "Marker: microsatellite",
  b_marker_typeSNP = "Marker: SNP",
  b_marker_typeSNPs = "Marker: SNPs",
  b_marker_typeother = "Marker: other",
  b_phi_PC1_thermal_geographic = "phi: PC1 thermal-geographic",
  b_phi_PC2_wet_productive = "phi: PC2 wet-productivity",
  b_phi_tmp_short_anom_5yr_C_first_month_z = "phi: short temperature anomaly",
  b_phi_pre_short_anom_5yr_mm_first_month_z = "phi: short precipitation anomaly",
  b_phi_tmp_long_anom_1961_1989_C_first_month_z = "phi: long temperature anomaly",
  b_phi_pre_long_anom_1961_1989_mm_first_month_z = "phi: long precipitation anomaly",
  b_phi_nest_boxyes = "phi: nest box",
  b_phi_marker_typemicrosatellite = "phi: marker: microsatellite",
  b_phi_marker_typeSNP = "phi: marker: SNP",
  b_phi_marker_typeSNPs = "phi: marker: SNPs",
  b_phi_marker_typeother = "phi: marker: other"
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
      component = ifelse(grepl("^b_phi_", term), "phi", "mu"),
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

write_csv(fixed_effects, file.path(out_dir, "fixed_effects_main_environmental_models.csv"))

capture.output(
  sessionInfo(),
  file = "logs/session_main_environmental_models.txt"
)

message("Done. Main environmental models saved in ", out_dir, ".")
