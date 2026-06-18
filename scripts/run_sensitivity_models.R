suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(posterior)
  library(loo)
})

out_dir <- "models/sensitivity_models"
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
    as.character(species_phylo) %in% rownames(phylo_mat)
  ) %>%
  mutate(
    tmp_cru_C_first_month_z = z(tmp_cru_C_first_month),
    pre_cru_mm_first_month_log_z = z(log1p(pre_cru_mm_first_month)),
    dtr_cru_C_first_month_z = z(dtr_cru_C_first_month),
    ndvi_first_month_median_30km_z = z(ndvi_first_month_median_30km),
    abs_latitude_z = z(abs_latitude),
    distance_to_coast_z = z(distance_to_coast_km)
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
    as_tibble(pca$x[, 1:3]) %>%
      transmute(
        PC1_thermal_geographic = PC1,
        PC2_wet_productive = PC2,
        PC3_environmental_residual = PC3
      )
  ) %>%
  droplevels()

phylo_mat_fit <- phylo_mat[
  levels(dat_fit$species_phylo),
  levels(dat_fit$species_phylo)
]

saveRDS(dat_fit, file.path(out_dir, "dat_sensitivity_models.rds"))
saveRDS(phylo_mat_fit, file.path(out_dir, "phylo_mat_sensitivity_models.rds"))

data_checks <- tibble(
  n_records = nrow(dat_fit),
  n_populations = n_distinct(dat_fit$population_id),
  n_species_phylo = n_distinct(dat_fit$species_phylo),
  n_species_nonphylo = n_distinct(dat_fit$species_nonphylo),
  total_broods = sum(dat_fit$n_broods_sampled),
  total_epp_broods = sum(dat_fit$n_epp_broods),
  raw_epp_rate = total_epp_broods / total_broods
)

write_csv(data_checks, file.path(out_dir, "sensitivity_model_data_checks.csv"))
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

m_sensitivity_first_month_mu <- fit_model(
  "m_sensitivity_first_month_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    tmp_cru_C_first_month_z +
    pre_cru_mm_first_month_log_z +
    dtr_cru_C_first_month_z +
    ndvi_first_month_median_30km_z +
    abs_latitude_z +
    distance_to_coast_z +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1501
)

m_sensitivity_first_month_mu_phi <- fit_model(
  "m_sensitivity_first_month_mu_phi",
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z +
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat_fit)),
    phi ~
      tmp_cru_C_first_month_z +
      pre_cru_mm_first_month_log_z +
      dtr_cru_C_first_month_z +
      ndvi_first_month_median_30km_z +
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type
  ),
  priors_mu_phi,
  seed = 1502
)

m_sensitivity_pca12_mu <- fit_model(
  "m_sensitivity_pca12_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    PC1_thermal_geographic +
    PC2_wet_productive +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1503
)

m_sensitivity_pca123_mu <- fit_model(
  "m_sensitivity_pca123_mu",
  n_epp_broods | trials(n_broods_sampled) ~
    PC1_thermal_geographic +
    PC2_wet_productive +
    PC3_environmental_residual +
    nest_box +
    marker_type +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1504
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
  m_sensitivity_first_month_mu = m_sensitivity_first_month_mu,
  m_sensitivity_first_month_mu_phi = m_sensitivity_first_month_mu_phi,
  m_sensitivity_pca12_mu = m_sensitivity_pca12_mu,
  m_sensitivity_pca123_mu = m_sensitivity_pca123_mu
)

diagnostics <- bind_rows(lapply(names(model_list), function(model_name) {
  diagnostic_summary(model_list[[model_name]], model_name)
}))

write_csv(diagnostics, file.path(out_dir, "diagnostics_sensitivity_models.csv"))
print(diagnostics)

loo_list <- lapply(names(model_list), function(model_name) {
  message("Computing LOO: ", model_name)
  loo_value <- loo(model_list[[model_name]], moment_match = TRUE)
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

write_csv(
  loo_diagnostics,
  file.path(out_dir, "loo_diagnostics_sensitivity_models.csv")
)

print(loo_diagnostics)

loo_pairs <- list(
  "Direct first-month predictors" = c(
    "m_sensitivity_first_month_mu",
    "m_sensitivity_first_month_mu_phi"
  ),
  "PCA dimensionality" = c(
    "m_sensitivity_pca12_mu",
    "m_sensitivity_pca123_mu"
  )
)

loo_comparison <- bind_rows(lapply(names(loo_pairs), function(hypothesis_group) {
  pair_names <- loo_pairs[[hypothesis_group]]
  pair_comparison <- loo_compare(loo_list[pair_names])

  as.data.frame(pair_comparison) %>%
    rownames_to_column("model") %>%
    mutate(
      hypothesis_group = hypothesis_group,
      model_version = case_when(
        grepl("_mu_phi$", model) ~ "mean-scale",
        grepl("pca123", model) ~ "PC1+PC2+PC3",
        grepl("pca12", model) ~ "PC1+PC2",
        TRUE ~ "mean-only"
      ),
      .before = model
    )
}))

write_csv(
  loo_comparison,
  file.path(out_dir, "loo_compare_sensitivity_models.csv")
)

capture.output(
  lapply(names(loo_pairs), function(hypothesis_group) {
    cat("\n", hypothesis_group, "\n", sep = "")
    print(loo_compare(loo_list[loo_pairs[[hypothesis_group]]]))
  }),
  file = file.path(out_dir, "loo_compare_sensitivity_models.txt")
)

fixed_effect_terms <- c(
  b_tmp_cru_C_first_month_z = "First-month mean temperature",
  b_pre_cru_mm_first_month_log_z = "First-month precipitation (log)",
  b_dtr_cru_C_first_month_z = "First-month DTR",
  b_ndvi_first_month_median_30km_z = "First-month NDVI",
  b_abs_latitude_z = "Absolute latitude",
  b_distance_to_coast_z = "Distance to coast",
  b_nest_boxyes = "Nest box",
  b_marker_typemicrosatellite = "Marker: microsatellite",
  b_marker_typeSNP = "Marker: SNP",
  b_marker_typeSNPs = "Marker: SNPs",
  b_marker_typeother = "Marker: other",
  b_PC1_thermal_geographic = "PC1 thermal-geographic",
  b_PC2_wet_productive = "PC2 wet-productivity",
  b_PC3_environmental_residual = "PC3 environmental residual",
  b_phi_tmp_cru_C_first_month_z = "phi: first-month mean temperature",
  b_phi_pre_cru_mm_first_month_log_z = "phi: first-month precipitation (log)",
  b_phi_dtr_cru_C_first_month_z = "phi: first-month DTR",
  b_phi_ndvi_first_month_median_30km_z = "phi: first-month NDVI",
  b_phi_abs_latitude_z = "phi: absolute latitude",
  b_phi_distance_to_coast_z = "phi: distance to coast",
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

write_csv(fixed_effects, file.path(out_dir, "fixed_effects_sensitivity_models.csv"))

capture.output(
  sessionInfo(),
  file = "logs/session_sensitivity_models.txt"
)

message("Done. Sensitivity models saved in ", out_dir, ".")
