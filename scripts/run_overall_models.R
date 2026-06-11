#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(brms)
  library(loo)
})

dir.create("models", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

message("Loading phylogeny-ready model data...")
dat_model_phylo <- readRDS("data/processed/dat_model_phylo.rds")
phylo_mat <- readRDS("data/processed/phylo_mat.rds")

dat_basic_phylo <- dat_model_phylo

dat_basic_phylo$nest_box <- as.factor(dat_basic_phylo$nest_box)
dat_basic_phylo$marker_type <- as.factor(dat_basic_phylo$marker_type)
dat_basic_phylo$id <- as.factor(dat_basic_phylo$id)
dat_basic_phylo$population_id <- as.factor(dat_basic_phylo$population_id)
dat_basic_phylo$species_nonphylo <- as.factor(dat_basic_phylo$species_nonphylo)
dat_basic_phylo$species_phylo <- as.factor(dat_basic_phylo$species_phylo)

dat_basic_phylo <- dat_basic_phylo[
  !is.na(dat_basic_phylo$n_epp_broods) &
    !is.na(dat_basic_phylo$n_broods_sampled) &
    dat_basic_phylo$n_broods_sampled > 0 &
    dat_basic_phylo$n_epp_broods <= dat_basic_phylo$n_broods_sampled &
    !is.na(dat_basic_phylo$id) &
    !is.na(dat_basic_phylo$population_id) &
    !is.na(dat_basic_phylo$species_nonphylo) &
    !is.na(dat_basic_phylo$species_phylo),
]

dat_basic_phylo <- droplevels(dat_basic_phylo)

stopifnot(all(dat_basic_phylo$species_phylo %in% rownames(phylo_mat)))

message("Overall model data:")
message("  records: ", nrow(dat_basic_phylo))
message("  species: ", length(unique(dat_basic_phylo$species_phylo)))
message("  populations: ", length(unique(dat_basic_phylo$population_id)))
message("  studies: ", length(unique(dat_basic_phylo$id)))

priors_basic <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

message("Fitting overall model: population identity only...")
m_basic_pop <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    1 +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_basic_phylo,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_basic,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_basic_pop, "models/m_basic_pop.rds")
capture.output(summary(m_basic_pop), file = "models/summary_m_basic_pop.txt")

message("Fitting overall model: population identity plus study identity...")
m_basic_pop_study <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    1 +
    (1 | id) +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_basic_phylo,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_basic,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 124,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_basic_pop_study, "models/m_basic_pop_study.rds")
capture.output(summary(m_basic_pop_study), file = "models/summary_m_basic_pop_study.txt")

message("Computing LOO for overall models...")
loo_basic_pop <- loo(m_basic_pop)
loo_basic_pop_study <- loo(m_basic_pop_study)

loo_comparison_basic <- loo_compare(
  loo_basic_pop,
  loo_basic_pop_study
)

saveRDS(loo_basic_pop, "models/loo_basic_pop.rds")
saveRDS(loo_basic_pop_study, "models/loo_basic_pop_study.rds")
saveRDS(loo_comparison_basic, "models/loo_comparison_basic.rds")
capture.output(loo_comparison_basic, file = "models/loo_comparison_basic.txt")

capture.output(sessionInfo(), file = "logs/session_overall_models.txt")

message("Done. Overall model outputs saved in models/.")
