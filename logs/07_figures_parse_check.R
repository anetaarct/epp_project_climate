library(knitr)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(tidyr)
library(purrr)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}

sensitivity_model_dir <- "models/sensitivity_models"

dat_location_file <- first_existing(c(
  "data/processed/dat_model_phylo.rds",
  "data/processed/dat_model.rds"
))

if (!is.na(dat_location_file)) {
  dat_locations <- readRDS(dat_location_file) %>%
    mutate(
      observed_epp = n_epp_broods / n_broods_sampled
    ) %>%
    filter(
      !is.na(lat),
      !is.na(long),
      !is.na(observed_epp),
      n_broods_sampled > 0
    )

  if (requireNamespace("maps", quietly = TRUE)) {
    world_map <- ggplot2::map_data("world")

    fig_sampling_locations <- ggplot() +
      geom_polygon(
        data = world_map,
        aes(x = long, y = lat, group = group),
        fill = "grey92",
        colour = "white",
        linewidth = 0.2
      ) +
      geom_point(
        data = dat_locations,
        aes(x = long, y = lat, colour = observed_epp, size = n_broods_sampled),
        alpha = 0.72
      ) +
      coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
      scale_colour_gradient(
        low = "#2f6f73",
        high = "#b35c35",
        labels = scales::percent_format(accuracy = 1)
      ) +
      scale_size_continuous(range = c(1, 4), trans = "sqrt") +
      labs(
        x = NULL,
        y = NULL,
        colour = "Observed EPP",
        size = "Sampled broods"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  } else {
    fig_sampling_locations <- ggplot(
      dat_locations,
      aes(x = long, y = lat, colour = observed_epp, size = n_broods_sampled)
    ) +
      geom_point(alpha = 0.72) +
      coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
      scale_colour_gradient(
        low = "#2f6f73",
        high = "#b35c35",
        labels = scales::percent_format(accuracy = 1)
      ) +
      scale_size_continuous(range = c(1, 4), trans = "sqrt") +
      labs(
        x = "Longitude",
        y = "Latitude",
        colour = "Observed EPP",
        size = "Sampled broods"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
  }

  ggsave(
    "figures/sampling_locations_map.png",
    fig_sampling_locations,
    width = 9,
    height = 5.2,
    dpi = 300
  )

  ggsave(
    "figures/sampling_locations_map.pdf",
    fig_sampling_locations,
    width = 9,
    height = 5.2
  )
} else {
  cat("Model-ready data are not available yet.")
}

if (file.exists("figures/sampling_locations_map.png")) {
  include_graphics("figures/sampling_locations_map.png")
} else {
  cat("Sampling-location map is not available yet.")
}

tree_file <- "data/processed/tree_pruned.rds"
dat_phylo_file <- "data/processed/dat_model_phylo.rds"

if (file.exists(tree_file) &&
    file.exists(dat_phylo_file) &&
    requireNamespace("ape", quietly = TRUE)) {
  tree_pruned <- readRDS(tree_file)
  dat_phylo <- readRDS(dat_phylo_file)

  high_epp_species <- dat_phylo %>%
    filter(
      !is.na(species_phylo),
      !is.na(n_epp_broods),
      !is.na(n_broods_sampled),
      n_broods_sampled > 0
    ) %>%
    group_by(species_phylo) %>%
    summarise(
      n_records = n(),
      total_broods = sum(n_broods_sampled),
      total_epp_broods = sum(n_epp_broods),
      species_epp_rate = total_epp_broods / total_broods,
      .groups = "drop"
    ) %>%
    filter(species_phylo %in% tree_pruned$tip.label) %>%
    arrange(desc(species_epp_rate), desc(total_broods)) %>%
    slice_head(n = 10)

  write_csv(
    high_epp_species,
    "figures/high_epp_species_for_phylogeny.csv"
  )

  top_species <- high_epp_species$species_phylo
  tip_colours <- ifelse(tree_pruned$tip.label %in% top_species, "#b35c35", "grey55")

  make_phylogeny_plot <- function() {
    ape::plot.phylo(
      tree_pruned,
      type = "fan",
      show.tip.label = FALSE,
      edge.color = "grey70",
      edge.width = 0.7,
      no.margin = TRUE
    )

    ape::tiplabels(
      pch = 21,
      bg = tip_colours,
      col = "white",
      cex = ifelse(tree_pruned$tip.label %in% top_species, 1.1, 0.45),
      lwd = 0.25
    )

    legend(
      "topleft",
      legend = c("Top species-level EPP", "Other species"),
      pt.bg = c("#b35c35", "grey55"),
      pch = 21,
      pt.cex = c(1.1, 0.7),
      bty = "n",
      cex = 0.85
    )

    legend(
      "bottomleft",
      legend = gsub("_", " ", top_species),
      pch = 21,
      pt.bg = "#b35c35",
      bty = "n",
      cex = 0.55,
      title = "Highest EPP species"
    )
  }

  png("figures/high_epp_species_phylogeny.png", width = 2400, height = 2400, res = 300)
  make_phylogeny_plot()
  dev.off()

  pdf("figures/high_epp_species_phylogeny.pdf", width = 8, height = 8)
  make_phylogeny_plot()
  dev.off()
} else {
  cat("Phylogenetic tree, model-ready data, or the ape package is not available yet.")
}

if (file.exists("figures/high_epp_species_phylogeny.png")) {
  include_graphics("figures/high_epp_species_phylogeny.png")
} else {
  cat("High-EPP phylogeny figure is not available yet.")
}

dat_model_file <- "data/processed/dat_model.rds"

if (file.exists(dat_model_file)) {
  dat_model <- readRDS(dat_model_file)

  environment_window_pairs <- tribble(
    ~variable, ~first_month, ~full_season,
    "Temperature", "tmp_cru_C_first_month", "tmp_cru_C_full_season",
    "Precipitation", "pre_cru_mm_first_month", "pre_cru_mm_full_season",
    "DTR", "dtr_cru_C_first_month", "dtr_cru_C_full_season",
    "NDVI", "ndvi_first_month_median_30km", "ndvi_full_season_median_30km"
  )

  environment_window_plot_data <- purrr::pmap_dfr(
    environment_window_pairs,
    function(variable, first_month, full_season) {
      dat_model %>%
        transmute(
          variable = variable,
          first_month = .data[[first_month]],
          full_season = .data[[full_season]]
        ) %>%
        drop_na()
    }
  )

  fig_first_full_season <- ggplot(
    environment_window_plot_data,
    aes(x = first_month, y = full_season)
  ) +
    geom_point(alpha = 0.35, size = 1.4, colour = "#2f6f73") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, colour = "#b35c35") +
    facet_wrap(~ variable, scales = "free", ncol = 2) +
    labs(
      x = "First-month value",
      y = "Full-season value"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/first_month_full_season_correlations.png",
    fig_first_full_season,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggsave(
    "figures/first_month_full_season_correlations.pdf",
    fig_first_full_season,
    width = 8,
    height = 6
  )
} else {
  cat("data/processed/dat_model.rds is not available yet. Render Data preparation first.")
}

if (file.exists("figures/first_month_full_season_correlations.png")) {
  include_graphics("figures/first_month_full_season_correlations.png")
} else {
  cat("Environmental-window correlation figure is not available yet.")
}

pca_model_set_dir <- "models/pca_environment_model_set"
pca_legacy_dir <- "models/pca_environment"

variance_file <- first_existing(c(
  file.path(pca_model_set_dir, "pca_variance_explained_environment_geography.csv"),
  file.path(pca_legacy_dir, "pca_variance_explained_environment_geography.csv")
))

loadings_file <- first_existing(c(
  file.path(pca_model_set_dir, "pca_loadings_environment_geography.csv"),
  file.path(pca_legacy_dir, "pca_loadings_environment_geography.csv")
))

if (!is.na(variance_file)) {
  variance_explained <- read_csv(variance_file, show_col_types = FALSE)

  fig_pca_variance <- variance_explained %>%
    mutate(PC = factor(PC, levels = PC)) %>%
    ggplot(aes(x = PC, y = variance_explained)) +
    geom_col(width = 0.7, fill = "#2f6f73") +
    geom_point(aes(y = cumulative_variance), colour = "#b35c35", size = 2) +
    geom_line(aes(y = cumulative_variance, group = 1), colour = "#b35c35") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Variance explained") +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/pca_variance_explained.png",
    fig_pca_variance,
    width = 6.5,
    height = 4.5,
    dpi = 300
  )

  ggsave(
    "figures/pca_variance_explained.pdf",
    fig_pca_variance,
    width = 6.5,
    height = 4.5
  )
} else {
  cat("PCA variance file is not available yet. Run the PCA model script first.")
}

if (!is.na(loadings_file)) {
  loadings <- read_csv(loadings_file, show_col_types = FALSE)

  if (!"PC1_thermal_geographic" %in% names(loadings) && "PC1" %in% names(loadings)) {
    loadings <- loadings %>%
      mutate(PC1_thermal_geographic = PC1)
  }

  if (!"PC2_wet_productive" %in% names(loadings) && "PC2" %in% names(loadings)) {
    loadings <- loadings %>%
      mutate(PC2_wet_productive = -PC2)
  }

  fig_pca_loadings <- loadings %>%
    select(variable, PC1_thermal_geographic, PC2_wet_productive) %>%
    pivot_longer(
      c(PC1_thermal_geographic, PC2_wet_productive),
      names_to = "axis",
      values_to = "loading"
    ) %>%
    mutate(
      axis = recode(axis,
        PC1_thermal_geographic = "PC1 thermal-geographic",
        PC2_wet_productive = "PC2 wet-productivity"
      )
    ) %>%
    ggplot(aes(x = loading, y = reorder(variable, loading), fill = loading)) +
    geom_col(width = 0.7) +
    facet_wrap(~ axis, scales = "free_x") +
    scale_fill_gradient2(
      low = "#6b8fb5",
      mid = "#f7f4ef",
      high = "#b35c35",
      midpoint = 0
    ) +
    labs(x = "Loading", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")

  ggsave(
    "figures/pca_loadings.png",
    fig_pca_loadings,
    width = 8,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    "figures/pca_loadings.pdf",
    fig_pca_loadings,
    width = 8,
    height = 4.8
  )
} else {
  cat("PCA loading file is not available yet. Run the PCA model script first.")
}

if (exists("fig_pca_variance") && exists("fig_pca_loadings")) {
  panel_a <- fig_pca_variance +
    ggplot2::labs(title = "A. Variance explained") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "#1f4d3a"))

  panel_b <- fig_pca_loadings +
    ggplot2::labs(title = "B. PCA loadings") +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "#1f4d3a"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    fig_environmental_pca <- panel_a + panel_b +
      patchwork::plot_layout(widths = c(0.8, 1.45))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    fig_environmental_pca <- gridExtra::arrangeGrob(
      panel_a,
      panel_b,
      ncol = 2,
      widths = c(0.8, 1.45)
    )
  } else {
    fig_environmental_pca <- grid::grid.grabExpr({
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(
        layout = grid::grid.layout(1, 2, widths = grid::unit(c(0.8, 1.45), "null"))
      ))
      print(panel_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(panel_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    })
  }

  ggsave(
    "figures/pca_environmental_pca.png",
    fig_environmental_pca,
    width = 10.5,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    "figures/pca_environmental_pca.pdf",
    fig_environmental_pca,
    width = 10.5,
    height = 4.8
  )
}

if (file.exists("figures/pca_environmental_pca.png")) {
  include_graphics("figures/pca_environmental_pca.png")
} else {
  cat("Environmental PCA figure is not available yet.")
}

sensitivity_mu_model_file <- file.path(
  "models/sensitivity_models",
  "m_sensitivity_first_month_mu.rds"
)

sensitivity_data_file <- first_existing(c(
  file.path("models/sensitivity_models", "dat_sensitivity_models.rds"),
  file.path("models/sensitivity_models", "dat_fit_sensitivity_models.rds")
))

if (file.exists(sensitivity_mu_model_file) &&
    !is.na(sensitivity_data_file) &&
    requireNamespace("brms", quietly = TRUE)) {
  sensitivity_mu_model <- readRDS(sensitivity_mu_model_file)
  dat_sensitivity_effects <- readRDS(sensitivity_data_file) %>%
    mutate(
      observed_epp = n_epp_broods / n_broods_sampled,
      nest_box = factor(nest_box),
      marker_type = factor(marker_type)
    )

  sensitivity_effect_variables <- tribble(
    ~variable, ~label,
    "tmp_cru_C_first_month_z", "Temperature",
    "pre_cru_mm_first_month_log_z", "Precipitation",
    "dtr_cru_C_first_month_z", "DTR",
    "ndvi_first_month_median_30km_z", "NDVI"
  )

  make_sensitivity_effect <- function(variable, label) {
    x_range <- range(dat_sensitivity_effects[[variable]], na.rm = TRUE)
    x_grid <- seq(x_range[1], x_range[2], length.out = 80)

    newdata <- dat_sensitivity_effects[rep(1, length(x_grid)), , drop = FALSE]
    newdata$tmp_cru_C_first_month_z <- 0
    newdata$pre_cru_mm_first_month_log_z <- 0
    newdata$dtr_cru_C_first_month_z <- 0
    newdata$ndvi_first_month_median_30km_z <- 0
    newdata$abs_latitude_z <- 0
    newdata$distance_to_coast_z <- 0
    newdata$nest_box <- factor(levels(dat_sensitivity_effects$nest_box)[1], levels = levels(dat_sensitivity_effects$nest_box))
    newdata$marker_type <- factor(levels(dat_sensitivity_effects$marker_type)[1], levels = levels(dat_sensitivity_effects$marker_type))
    newdata$n_broods_sampled <- 1
    newdata[[variable]] <- x_grid

    fitted_values <- fitted(
      sensitivity_mu_model,
      newdata = newdata,
      re_formula = NA,
      scale = "response"
    )

    tibble(
      variable = label,
      x = x_grid,
      estimate = fitted_values[, "Estimate"],
      lower = fitted_values[, "Q2.5"],
      upper = fitted_values[, "Q97.5"]
    )
  }

  sensitivity_effect_plot_data <- pmap_dfr(
    sensitivity_effect_variables,
    make_sensitivity_effect
  )

  sensitivity_effect_raw_data <- pmap_dfr(
    sensitivity_effect_variables,
    function(variable, label) {
      dat_sensitivity_effects %>%
        transmute(
          variable = label,
          x = .data[[variable]],
          observed_epp = observed_epp
        ) %>%
        drop_na()
    }
  )

  sensitivity_effect_levels <- sensitivity_effect_variables$label

  fig_sensitivity_environmental_effects <- ggplot() +
    geom_point(
      data = sensitivity_effect_raw_data,
      aes(x = x, y = observed_epp),
      alpha = 0.12,
      size = 0.8,
      colour = "grey35"
    ) +
    geom_ribbon(
      data = sensitivity_effect_plot_data,
      aes(x = x, ymin = lower, ymax = upper),
      fill = "#9fb9bd",
      alpha = 0.35
    ) +
    geom_line(
      data = sensitivity_effect_plot_data,
      aes(x = x, y = estimate),
      colour = "#b35c35",
      linewidth = 0.8
    ) +
    facet_wrap(
      ~ factor(variable, levels = sensitivity_effect_levels),
      scales = "free_x",
      ncol = 2
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      x = "Standardised predictor value",
      y = "Predicted EPP"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(colour = "grey25")
    )

  ggsave(
    "figures/sensitivity_environmental_effects.png",
    fig_sensitivity_environmental_effects,
    width = 8.5,
    height = 5.8,
    dpi = 300
  )

  ggsave(
    "figures/sensitivity_environmental_effects.pdf",
    fig_sensitivity_environmental_effects,
    width = 8.5,
    height = 5.8
  )
} else {
  cat("Sensitivity mean-only model or data are not available yet.")
}

if (file.exists("figures/sensitivity_environmental_effects.png")) {
  include_graphics("figures/sensitivity_environmental_effects.png")
} else {
  cat("Sensitivity environmental effect panel is not available yet.")
}

sensitivity_mu_model_file <- file.path(
  sensitivity_model_dir,
  "m_sensitivity_first_month_mu.rds"
)

sensitivity_fixed_effect_terms <- c(
  b_tmp_cru_C_first_month_z = "First-month mean temperature",
  b_pre_cru_mm_first_month_log_z = "First-month precipitation (log)",
  b_dtr_cru_C_first_month_z = "First-month DTR",
  b_ndvi_first_month_median_30km_z = "First-month NDVI",
  b_abs_latitude_z = "Absolute latitude",
  b_distance_to_coast_z = "Distance to coast",
  b_nest_boxyes = "Nest box",
  b_marker_typemicrosatellite = "Marker: microsatellite",
  b_marker_typeSNPs = "Marker: SNPs",
  b_marker_typeother = "Marker: other"
)

if (file.exists(sensitivity_mu_model_file) &&
    requireNamespace("posterior", quietly = TRUE) &&
    requireNamespace("bayesplot", quietly = TRUE)) {
  sensitivity_mu_model <- readRDS(sensitivity_mu_model_file)
  sensitivity_draws <- posterior::as_draws_df(sensitivity_mu_model)
  available_terms <- intersect(names(sensitivity_fixed_effect_terms), names(sensitivity_draws))

  sensitivity_fixed_effects <- map_dfr(available_terms, function(term) {
    x <- sensitivity_draws[[term]]
    tibble(
      model = "Mean-only sensitivity",
      term = term,
      predictor = sensitivity_fixed_effect_terms[[term]],
      component = "mu",
      estimate = median(x),
      lower = quantile(x, 0.025),
      upper = quantile(x, 0.975),
      prob_positive = mean(x > 0),
      prob_negative = mean(x < 0)
    )
  })

  write_csv(
    sensitivity_fixed_effects,
    file.path(sensitivity_model_dir, "fixed_effects_sensitivity_models.csv")
  )

  term_order <- names(sensitivity_fixed_effect_terms)[
    names(sensitivity_fixed_effect_terms) %in% available_terms
  ]

  sensitivity_draws_matrix <- as.matrix(sensitivity_draws[, term_order, drop = FALSE])

  fig_sensitivity_fixed_effects <- bayesplot::mcmc_intervals(
    sensitivity_draws_matrix,
    pars = term_order,
    prob = 0.8,
    prob_outer = 0.95,
    point_est = "median",
    outer_size = 0.7,
    inner_size = 1.8,
    point_size = 2.4
  ) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
    scale_y_discrete(
      limits = rev(term_order),
      labels = sensitivity_fixed_effect_terms
    ) +
    labs(
      title = "First-month mean-only",
      x = NULL,
      y = NULL
    ) +
    bayesplot::theme_default(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, colour = "grey25", face = "plain"),
      axis.text.y = element_text(colour = "grey25"),
      axis.text.x = element_text(colour = "grey25"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      panel.grid.major.x = element_line(colour = "grey90")
    )

  ggsave(
    "figures/sensitivity_fixed_effects_forest_plot.png",
    fig_sensitivity_fixed_effects,
    width = 9,
    height = 3.2,
    dpi = 300
  )

  ggsave(
    "figures/sensitivity_fixed_effects_forest_plot.pdf",
    fig_sensitivity_fixed_effects,
    width = 9,
    height = 3.2
  )
} else {
  cat("Sensitivity fixed-effect figure will be generated after the model object and bayesplot are available.")
}
if (file.exists("figures/sensitivity_fixed_effects_forest_plot.png")) {
  include_graphics("figures/sensitivity_fixed_effects_forest_plot.png")
} else {
  cat("Sensitivity fixed-effect figure is not available yet.")
}
