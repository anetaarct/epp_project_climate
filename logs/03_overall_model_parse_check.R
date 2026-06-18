library(dplyr)
library(readr)
library(tibble)
library(knitr)

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}

overall_summary_file <- first_existing(c(
  "models/overall_population/summary_m_overall_population.txt",
  "models/overall_model/summary_m_overall_population.txt",
  "models/summary_m_overall_population.txt"
))

overall_diagnostics_file <- first_existing(c(
  "models/overall_population/overall_population_model_diagnostics.csv",
  "models/overall_population/diagnostics_m_overall_population.csv",
  "models/overall_model/overall_model_diagnostics.csv"
))

overall_model_file <- first_existing(c(
  "models/overall_population/m_overall_population.rds",
  "models/overall_model/m_overall_population.rds",
  "models/m_overall_population.rds"
))

# m_overall_population <- brms::brm(
#   n_epp_broods | trials(n_broods_sampled) ~
#     1 +
#     (1 | population_id) +
#     (1 | species_nonphylo) +
#     (1 | gr(species_phylo, cov = phylo_mat_fit)),
#   data = dat_model,
#   data2 = list(phylo_mat_fit = phylo_mat_fit),
#   family = brms::beta_binomial(link = "logit", link_phi = "log")
# )

if (!is.na(overall_summary_file)) {
  cat("<details><summary>Show overall model summary</summary>\n\n```\n")
  cat(readLines(overall_summary_file), sep = "\n")
  cat("\n```\n</details>\n")
} else {
  cat("Overall model summary is not available yet.")
}

if (file.exists("data/processed/dat_model_phylo.rds")) {
  overlap_dat <- readRDS("data/processed/dat_model_phylo.rds")
} else if (file.exists("data/processed/dat_model.rds")) {
  overlap_dat <- readRDS("data/processed/dat_model.rds")
} else {
  overlap_dat <- read_csv("epp_data_June2026.csv", show_col_types = FALSE) %>%
    mutate(
      lat_round = round(as.numeric(lat), 2),
      long_round = round(as.numeric(long), 2),
      population_id = paste(scientific_name, lat_round, long_round, sep = "_")
    )
}

overlap_dat <- overlap_dat %>%
  mutate(
    id_reference = if ("id_reference" %in% names(.)) {
      as.character(id_reference)
    } else {
      as.character(id)
    }
  ) %>%
  mutate(
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

population_overlap %>%
  count(n_references, name = "n_populations") %>%
  mutate(percent_populations = 100 * n_populations / sum(n_populations)) %>%
  knitr::kable(digits = 2)

if (!is.na(overall_diagnostics_file)) {
  read_csv(overall_diagnostics_file, show_col_types = FALSE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    knitr::kable()
} else if (!is.na(overall_model_file) &&
           requireNamespace("brms", quietly = TRUE) &&
           requireNamespace("posterior", quietly = TRUE)) {
  m_overall <- readRDS(overall_model_file)

  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(m_overall),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  sampler_diagnostics <- brms::nuts_params(m_overall)

  diagnostic_table <- tibble(
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

  knitr::kable(diagnostic_table, digits = 3)
} else {
  cat("Overall model diagnostics are not available yet.")
}
