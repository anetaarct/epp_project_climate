suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(tibble)
  library(posterior)
})

out_dir <- "models/overall_population"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

chains <- 4
cores <- 4
iter <- 10000
warmup <- 5000

message("Loading model-ready phylogenetic data...")

dat <- readRDS("data/processed/dat_model_phylo.rds")
phylo_mat <- readRDS("data/processed/phylo_mat.rds")

dat_model <- dat %>%
  mutate(
    population_id = factor(population_id),
    species_phylo = factor(species_phylo),
    species_nonphylo = factor(species_nonphylo),
    id_reference = if ("id_reference" %in% names(.)) {
      as.character(id_reference)
    } else {
      as.character(id)
    }
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
    as.character(species_phylo) %in% rownames(phylo_mat)
  ) %>%
  droplevels()

phylo_mat_fit <- phylo_mat[
  levels(dat_model$species_phylo),
  levels(dat_model$species_phylo)
]

saveRDS(dat_model, file.path(out_dir, "dat_overall_population.rds"))
saveRDS(phylo_mat_fit, file.path(out_dir, "phylo_mat_overall_population.rds"))

data_checks <- tibble(
  n_records = nrow(dat_model),
  n_populations = n_distinct(dat_model$population_id),
  n_references = n_distinct(dat_model$id_reference),
  n_species_phylo = n_distinct(dat_model$species_phylo),
  n_species_nonphylo = n_distinct(dat_model$species_nonphylo),
  total_broods = sum(dat_model$n_broods_sampled),
  total_epp_broods = sum(dat_model$n_epp_broods),
  raw_epp_rate = total_epp_broods / total_broods
)

write_csv(data_checks, file.path(out_dir, "overall_population_data_checks.csv"))
print(data_checks)

overlap_dat <- dat_model %>%
  mutate(
    id_reference = as.character(id_reference),
    population_id = as.character(population_id)
  ) %>%
  filter(!is.na(id_reference), !is.na(population_id))

population_overlap <- overlap_dat %>%
  group_by(population_id) %>%
  summarise(
    n_records = n(),
    n_references = n_distinct(id_reference),
    id_references = paste(sort(unique(id_reference)), collapse = "; "),
    .groups = "drop"
  )

reference_overlap <- overlap_dat %>%
  group_by(id_reference) %>%
  summarise(
    n_records = n(),
    n_populations = n_distinct(population_id),
    .groups = "drop"
  )

overlap_summary <- tibble(
  metric = c(
    "Records",
    "References",
    "Populations",
    "Populations linked to one reference",
    "Populations linked to more than one reference",
    "References linked to one population",
    "References linked to more than one population"
  ),
  value = c(
    nrow(overlap_dat),
    n_distinct(overlap_dat$id_reference),
    n_distinct(overlap_dat$population_id),
    sum(population_overlap$n_references == 1),
    sum(population_overlap$n_references > 1),
    sum(reference_overlap$n_populations == 1),
    sum(reference_overlap$n_populations > 1)
  )
)

write_csv(overlap_summary, file.path(out_dir, "reference_population_overlap_summary.csv"))
write_csv(population_overlap, file.path(out_dir, "population_reference_overlap.csv"))
write_csv(reference_overlap, file.path(out_dir, "reference_population_overlap.csv"))

priors_overall <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "phi")
)

message("Fitting overall population model...")

m_overall_population <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    1 +
    (1 | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  data = dat_model,
  data2 = list(phylo_mat_fit = phylo_mat_fit),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_overall,
  chains = chains,
  cores = cores,
  iter = iter,
  warmup = warmup,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 1301,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_overall_population, file.path(out_dir, "m_overall_population.rds"))

capture.output(
  summary(m_overall_population),
  file = file.path(out_dir, "summary_m_overall_population.txt")
)

draw_diagnostics <- posterior::summarise_draws(
  posterior::as_draws_df(m_overall_population),
  "rhat",
  "ess_bulk",
  "ess_tail"
)

sampler_diagnostics <- brms::nuts_params(m_overall_population)

overall_diagnostics <- tibble(
  model = "m_overall_population",
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

write_csv(
  overall_diagnostics,
  file.path(out_dir, "overall_population_model_diagnostics.csv")
)

print(overall_diagnostics)

capture.output(
  sessionInfo(),
  file = "logs/session_overall_model.txt"
)

message("Done. Overall model saved in ", out_dir, ".")
