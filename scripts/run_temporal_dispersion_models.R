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

dat <- dat_model_phylo

needs_taxonomy <- !all(c("order", "family", "family_common") %in% names(dat)) ||
  ("family" %in% names(dat) && all(is.na(dat$family)))

if (needs_taxonomy && file.exists("epp_data_June2026.csv")) {
  message("Adding family taxonomy from epp_data_June2026.csv...")
  epp_taxonomy <- read.csv(
    "epp_data_June2026.csv",
    stringsAsFactors = FALSE
  )

  taxonomy_cols <- c(
    "record_id",
    "epp_row_id",
    "order",
    "family",
    "family_common"
  )

  epp_taxonomy <- epp_taxonomy[
    ,
    taxonomy_cols[taxonomy_cols %in% names(epp_taxonomy)],
    drop = FALSE
  ]

  if ("record_id" %in% names(dat) && "record_id" %in% names(epp_taxonomy)) {
    dat$.row_order <- seq_len(nrow(dat))
    dat <- dat[
      ,
      !names(dat) %in% c("order", "family", "family_common"),
      drop = FALSE
    ]
    epp_taxonomy <- epp_taxonomy[!duplicated(epp_taxonomy$record_id), ]
    dat <- merge(
      dat,
      epp_taxonomy[
        ,
        c("record_id", "order", "family", "family_common"),
        drop = FALSE
      ],
      by = "record_id",
      all.x = TRUE,
      sort = FALSE
    )
    dat <- dat[order(dat$.row_order), ]
    dat$.row_order <- NULL
  } else if ("epp_row_id" %in% names(dat) && "epp_row_id" %in% names(epp_taxonomy)) {
    dat$.row_order <- seq_len(nrow(dat))
    dat <- dat[
      ,
      !names(dat) %in% c("order", "family", "family_common"),
      drop = FALSE
    ]
    epp_taxonomy <- epp_taxonomy[!duplicated(epp_taxonomy$epp_row_id), ]
    dat <- merge(
      dat,
      epp_taxonomy[
        ,
        c("epp_row_id", "order", "family", "family_common"),
        drop = FALSE
      ],
      by = "epp_row_id",
      all.x = TRUE,
      sort = FALSE
    )
    dat <- dat[order(dat$.row_order), ]
    dat$.row_order <- NULL
  } else {
    warning(
      "Could not add family taxonomy because no shared record_id or epp_row_id was found."
    )
  }
}

dat$nest_box <- as.factor(dat$nest_box)
dat$marker_type <- as.factor(dat$marker_type)
dat$id <- as.factor(dat$id)
dat$population_id <- as.factor(dat$population_id)
dat$species_nonphylo <- as.factor(dat$species_nonphylo)
dat$species_phylo <- as.factor(dat$species_phylo)

dat$year_of_study_start <- as.numeric(dat$year_of_study_start)
dat$year_of_study_end <- as.numeric(dat$year_of_study_end)
dat$study_year <- ifelse(
  !is.na(dat$year_of_study_end),
  rowMeans(
    cbind(dat$year_of_study_start, dat$year_of_study_end),
    na.rm = TRUE
  ),
  dat$year_of_study_start
)

dat$abs_latitude_z <- z(dat$abs_latitude)
dat$distance_to_coast_z <- z(dat$distance_to_coast_km)
dat$ndvi_first_month_median_30km_z <- z(dat$ndvi_first_month_median_30km)
dat$era5_tmp_sd_C_first_month_z <- z(dat$era5_tmp_sd_C_first_month)
dat$era5_pre_sd_mm_day_first_month_log_z <- z(log1p(dat$era5_pre_sd_mm_day_first_month))
dat$era5_pre_wet_days_gt1mm_first_month_z <- z(dat$era5_pre_wet_days_gt1mm_first_month)

base_vars <- c(
  "n_epp_broods",
  "n_broods_sampled",
  "population_id",
  "species_nonphylo",
  "species_phylo"
)

taxonomy_vars <- c(
  "order",
  "family",
  "family_common"
)

taxonomy_vars <- taxonomy_vars[
  taxonomy_vars %in% names(dat)
]

temporal_vars <- c(
  base_vars,
  taxonomy_vars,
  "study_year"
)

dispersion_vars <- c(
  base_vars,
  "abs_latitude_z",
  "distance_to_coast_z",
  "nest_box",
  "marker_type",
  "ndvi_first_month_median_30km_z",
  "era5_tmp_sd_C_first_month_z",
  "era5_pre_sd_mm_day_first_month_log_z",
  "era5_pre_wet_days_gt1mm_first_month_z"
)

dat_temporal <- complete_model_data(dat, temporal_vars)

population_year_counts <- aggregate(
  study_year ~ population_id,
  data = dat_temporal,
  FUN = function(x) length(unique(x))
)

names(population_year_counts) <- c("population_id", "n_years")

repeated_populations <- population_year_counts$population_id[
  population_year_counts$n_years > 1
]

dat_temporal_repeated <- dat_temporal[
  dat_temporal$population_id %in% repeated_populations,
]

family_candidates <- c(
  "family",
  "Family",
  "family_name",
  "Family_name",
  "BIRDBASE_family"
)

family_col <- family_candidates[
  family_candidates %in% names(dat_temporal_repeated)
][1]

repeated_family_count <- if (!is.na(family_col)) {
  length(unique(na.omit(dat_temporal_repeated[[family_col]])))
} else {
  NA_integer_
}

repeated_order_count <- if ("order" %in% names(dat_temporal_repeated)) {
  length(unique(na.omit(dat_temporal_repeated$order)))
} else {
  NA_integer_
}

repeated_data_summary <- data.frame(
  quantity = c(
    "populations with >1 study year",
    "species represented by repeated populations",
    "families represented by repeated populations",
    "orders represented by repeated populations",
    "records in repeated populations"
  ),
  value = c(
    length(unique(dat_temporal_repeated$population_id)),
    length(unique(dat_temporal_repeated$species_phylo)),
    repeated_family_count,
    repeated_order_count,
    nrow(dat_temporal_repeated)
  )
)

population_year_distribution <- as.data.frame(
  table(population_year_counts$n_years)
)

names(population_year_distribution) <- c(
  "n_distinct_study_years",
  "n_populations"
)

write.csv(
  repeated_data_summary,
  "models/repeated_population_summary.csv",
  row.names = FALSE
)

write.csv(
  population_year_distribution,
  "models/population_year_distribution.csv",
  row.names = FALSE
)

message("Repeated-population summary:")
print(repeated_data_summary)

if (is.na(repeated_family_count)) {
  message(
    "Family-level counts were not calculated because no family column is present in the current dataset."
  )
}

dat_temporal <- dat_temporal[
  dat_temporal$population_id %in% repeated_populations,
]

dat_temporal <- droplevels(dat_temporal)
dat_temporal$population_mean_year <- ave(
  dat_temporal$study_year,
  dat_temporal$population_id,
  FUN = mean
)
dat_temporal$year_within_population <- dat_temporal$study_year -
  dat_temporal$population_mean_year
dat_temporal$year_within_population_z <- z(dat_temporal$year_within_population)

dat_dispersion <- complete_model_data(dat, dispersion_vars)

stopifnot(all(dat_temporal$species_phylo %in% rownames(phylo_mat)))
stopifnot(all(dat_dispersion$species_phylo %in% rownames(phylo_mat)))

model_data_checks <- data.frame(
  model = c("within-population temporal", "environmental dispersion"),
  records = c(nrow(dat_temporal), nrow(dat_dispersion)),
  species = c(
    length(unique(dat_temporal$species_phylo)),
    length(unique(dat_dispersion$species_phylo))
  ),
  populations = c(
    length(unique(dat_temporal$population_id)),
    length(unique(dat_dispersion$population_id))
  )
)

write.csv(
  model_data_checks,
  "models/temporal_dispersion_model_data_checks.csv",
  row.names = FALSE
)

message("Temporal and dispersion model data:")
print(model_data_checks)

priors_temporal <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd")
)

priors_dispersion <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 1), class = "b", dpar = "phi"),
  prior(normal(0, 2), class = "Intercept", dpar = "phi")
)

message("Fitting temporal model: within-population year effect...")
m_temporal_population_year <- brm(
  n_epp_broods | trials(n_broods_sampled) ~
    year_within_population_z +
    (1 + year_within_population_z | population_id) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat)),
  data = dat_temporal,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_temporal,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 301,
  save_pars = save_pars(all = TRUE)
)

saveRDS(
  m_temporal_population_year,
  "models/m_temporal_population_year.rds"
)
capture.output(
  summary(m_temporal_population_year),
  file = "models/summary_m_temporal_population_year.txt"
)

message("Fitting dispersion model: phi as a function of environment...")
m_phi_environment <- brm(
  bf(
    n_epp_broods | trials(n_broods_sampled) ~
      abs_latitude_z +
      distance_to_coast_z +
      nest_box +
      marker_type +
      ndvi_first_month_median_30km_z +
      era5_tmp_sd_C_first_month_z +
      era5_pre_sd_mm_day_first_month_log_z +
      era5_pre_wet_days_gt1mm_first_month_z +
      (1 | population_id) +
      (1 | species_nonphylo) +
      (1 | gr(species_phylo, cov = phylo_mat)),
    phi ~
      ndvi_first_month_median_30km_z +
      era5_tmp_sd_C_first_month_z +
      era5_pre_sd_mm_day_first_month_log_z +
      era5_pre_wet_days_gt1mm_first_month_z
  ),
  data = dat_dispersion,
  data2 = list(phylo_mat = phylo_mat),
  family = beta_binomial(link = "logit", link_phi = "log"),
  prior = priors_dispersion,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 302,
  save_pars = save_pars(all = TRUE)
)

saveRDS(m_phi_environment, "models/m_phi_environment.rds")
capture.output(
  summary(m_phi_environment),
  file = "models/summary_m_phi_environment.txt"
)

message("Computing LOO for temporal and dispersion models...")
loo_temporal_population_year <- loo(m_temporal_population_year)
loo_phi_environment <- loo(m_phi_environment)

saveRDS(
  loo_temporal_population_year,
  "models/loo_temporal_population_year.rds"
)
saveRDS(
  loo_phi_environment,
  "models/loo_phi_environment.rds"
)

capture.output(sessionInfo(), file = "logs/session_temporal_dispersion_models.txt")

message("Done. Temporal and dispersion model outputs saved in models/.")
