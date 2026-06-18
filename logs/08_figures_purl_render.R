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

publication_model_dir <- "models/publication_environmental_models"

pca_model_file <- first_existing(c(
  file.path(publication_model_dir, "m_pub_pca_mu.rds")
))

pca_publication_data_file <- first_existing(c(
  file.path(publication_model_dir, "dat_publication_environmental_models.rds")
))

if (!is.na(pca_model_file) &&
    !is.na(pca_publication_data_file) &&
    requireNamespace("brms", quietly = TRUE)) {
  fit_pca_publication <- readRDS(pca_model_file)
  dat_pca_publication <- readRDS(pca_publication_data_file) %>%
    filter(
      !is.na(PC1_thermal_geographic),
      !is.na(PC2_wet_productive),
      !is.na(n_epp_broods),
      !is.na(n_broods_sampled),
      n_broods_sampled > 0
    ) %>%
    mutate(observed_epp = n_epp_broods / n_broods_sampled)

  make_pca_effect_grid <- function(axis_name, axis_label) {
    axis_seq <- seq(
      min(dat_pca_publication[[axis_name]], na.rm = TRUE),
      max(dat_pca_publication[[axis_name]], na.rm = TRUE),
      length.out = 100
    )

    newdat <- tibble(
      PC1_thermal_geographic = 0,
      PC2_wet_productive = 0,
      n_broods_sampled = 1,
      population_id = levels(factor(dat_pca_publication$population_id))[1],
      species_nonphylo = levels(factor(dat_pca_publication$species_nonphylo))[1],
      species_phylo = levels(factor(dat_pca_publication$species_phylo))[1]
    ) %>%
      slice(rep(1, length(axis_seq)))

    newdat[[axis_name]] <- axis_seq

    pred <- as_tibble(
      fitted(fit_pca_publication, newdata = newdat, re_formula = NA, scale = "response")
    ) %>%
      bind_cols(newdat) %>%
      transmute(
        axis = axis_label,
        axis_value = .data[[axis_name]],
        estimate = pmin(1, Estimate),
        lower = pmax(0, Q2.5),
        upper = pmin(1, Q97.5)
      )

    raw <- dat_pca_publication %>%
      transmute(
        axis = axis_label,
        axis_value = .data[[axis_name]],
        observed_epp = observed_epp,
        n_broods_sampled = n_broods_sampled
      )

    list(pred = pred, raw = raw)
  }

  pc1_effect <- make_pca_effect_grid(
    "PC1_thermal_geographic",
    "A. PC1 thermal-geographic gradient"
  )
  pc2_effect <- make_pca_effect_grid(
    "PC2_wet_productive",
    "B. PC2 wet-productivity gradient"
  )

  pca_effect_pred <- bind_rows(pc1_effect$pred, pc2_effect$pred)
  pca_effect_raw <- bind_rows(pc1_effect$raw, pc2_effect$raw)

  fig_main_pca_effects <- ggplot() +
    geom_point(
      data = pca_effect_raw,
      aes(x = axis_value, y = observed_epp, size = n_broods_sampled),
      alpha = 0.16,
      colour = "#2f6f73"
    ) +
    geom_ribbon(
      data = pca_effect_pred,
      aes(x = axis_value, ymin = lower, ymax = upper),
      fill = "#9fb9bd",
      alpha = 0.45
    ) +
    geom_line(
      data = pca_effect_pred,
      aes(x = axis_value, y = estimate),
      colour = "#b35c35",
      linewidth = 0.9
    ) +
    facet_wrap(~ axis, scales = "free_x", ncol = 2) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_size_continuous(range = c(0.5, 2.5), trans = "sqrt", guide = "none") +
    labs(
      x = NULL,
      y = "Predicted EPP probability"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", colour = "#1f4d3a"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    "figures/main_model_pca_effects.png",
    fig_main_pca_effects,
    width = 9,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    "figures/main_model_pca_effects.pdf",
    fig_main_pca_effects,
    width = 9,
    height = 4.8
  )
} else {
  cat("Main PCA publication model or data are not available yet.")
}

if (file.exists("figures/main_model_pca_effects.png")) {
  include_graphics("figures/main_model_pca_effects.png")
} else {
  cat("Main model PCA effect figure is not available yet.")
}

publication_model_dir <- "models/publication_environmental_models"

publication_model_files <- c(
  "Environmental PCA" = file.path(publication_model_dir, "m_pub_pca_mu.rds"),
  "Short anomaly" = file.path(publication_model_dir, "m_pub_short_mu.rds"),
  "Long anomaly" = file.path(publication_model_dir, "m_pub_long_mu.rds")
)

publication_funnel_files <- c(
  "Environmental PCA" = "figures/publication_funnel_pca_mu.png",
  "Short anomaly" = "figures/publication_funnel_short_mu.png",
  "Long anomaly" = "figures/publication_funnel_long_mu.png"
)

publication_data_file <- first_existing(c(
  file.path(publication_model_dir, "dat_publication_environmental_models.rds")
))

if (!is.na(publication_data_file) &&
    all(file.exists(publication_model_files)) &&
    requireNamespace("brms", quietly = TRUE)) {
  dat_publication <- readRDS(publication_data_file)

  make_publication_funnel_data <- function(model_file, model_label) {
    fit <- readRDS(model_file)

    epred <- brms::posterior_epred(fit, re_formula = NA)
    expected_counts <- colMeans(epred)
    expected <- ifelse(
      expected_counts > 1,
      expected_counts / dat_publication$n_broods_sampled,
      expected_counts
    )

    dat_publication %>%
      transmute(
        model = model_label,
        n_broods_sampled = n_broods_sampled,
        observed = n_epp_broods / n_broods_sampled,
        expected = expected,
        lower = pmax(0, expected - 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled)),
        upper = pmin(1, expected + 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled))
      )
  }

  publication_funnel_plot_data <- imap_dfr(
    publication_model_files,
    make_publication_funnel_data
  )

  make_publication_funnel_plot <- function(plot_data, model_label) {
    ggplot(plot_data, aes(x = n_broods_sampled, y = observed)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9fb9bd", alpha = 0.35) +
      geom_line(aes(y = expected), colour = "#b35c35", linewidth = 0.7) +
      geom_point(alpha = 0.45, size = 1.2, colour = "#2f6f73") +
      scale_x_log10() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        x = "Number of sampled broods (log scale)",
        y = "Observed EPP",
        title = model_label
      ) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", colour = "#1f4d3a"))
  }

  for (model_label in names(publication_funnel_files)) {
    fig_publication_funnel <- publication_funnel_plot_data %>%
      filter(model == model_label) %>%
      make_publication_funnel_plot(model_label)

    ggsave(
      publication_funnel_files[[model_label]],
      fig_publication_funnel,
      width = 6.5,
      height = 4.8,
      dpi = 300
    )

    ggsave(
      sub("\\.png$", ".pdf", publication_funnel_files[[model_label]]),
      fig_publication_funnel,
      width = 6.5,
      height = 4.8
    )
  }
} else {
  cat("Publication model objects or model data are not available yet.")
}

publication_funnel_files <- c(
  "figures/publication_funnel_pca_mu.png",
  "figures/publication_funnel_short_mu.png",
  "figures/publication_funnel_long_mu.png"
)

existing_publication_funnels <- publication_funnel_files[file.exists(publication_funnel_files)]

if (length(existing_publication_funnels) > 0) {
  include_graphics(existing_publication_funnels)
} else {
  cat("Publication model funnel plots are not available yet.")
}

publication_fixed_effects_file <- file.path(
  "models/publication_environmental_models",
  "fixed_effects_publication_environmental_models.csv"
)

if (file.exists(publication_fixed_effects_file)) {
  publication_fixed_effects <- read_csv(publication_fixed_effects_file, show_col_types = FALSE) %>%
    mutate(
      model = recode(
        model,
        m_pub_pca_mu = "Environmental PCA",
        m_pub_short_mu = "Short anomaly",
        m_pub_long_mu = "Long anomaly"
      ),
      predictor = factor(predictor, levels = rev(unique(predictor)))
    )

  fig_publication_fixed_effects <- ggplot(
    publication_fixed_effects,
    aes(x = estimate, y = predictor)
  ) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.75, colour = "grey45") +
    geom_point(shape = 21, size = 2.5, fill = "white", colour = "grey25", stroke = 0.35) +
    facet_wrap(~ model, scales = "free_y", ncol = 1) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(colour = "grey25"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    "figures/publication_fixed_effects_summary.png",
    fig_publication_fixed_effects,
    width = 8.5,
    height = 4.8,
    dpi = 300
  )

  ggsave(
    "figures/publication_fixed_effects_summary.pdf",
    fig_publication_fixed_effects,
    width = 8.5,
    height = 4.8
  )
} else {
  cat("Publication fixed-effect summaries are not available yet.")
}

if (file.exists("figures/publication_fixed_effects_summary.png")) {
  include_graphics("figures/publication_fixed_effects_summary.png")
} else {
  cat("Publication fixed-effect summary figure is not available yet.")
}

main_model_dir <- "models/main_environmental_models"

main_model_files <- c(
  "PCA mean-only" = file.path(main_model_dir, "m_main_pca_mu.rds"),
  "PCA mean-scale" = file.path(main_model_dir, "m_main_pca_mu_phi.rds"),
  "Short anomaly mean-only" = file.path(main_model_dir, "m_main_short_mu.rds"),
  "Short anomaly mean-scale" = file.path(main_model_dir, "m_main_short_mu_phi.rds"),
  "Long anomaly mean-only" = file.path(main_model_dir, "m_main_long_mu.rds"),
  "Long anomaly mean-scale" = file.path(main_model_dir, "m_main_long_mu_phi.rds")
)

main_model_files <- main_model_files[file.exists(main_model_files)]

main_data_file <- first_existing(c(
  file.path(main_model_dir, "dat_main_environmental_models.rds")
))

if (length(main_model_files) > 0 &&
    !is.na(main_data_file) &&
    requireNamespace("brms", quietly = TRUE)) {
  dat_main <- readRDS(main_data_file)

  make_model_diagnostics <- function(model_file, model_label) {
    fit <- readRDS(model_file)

    yrep <- brms::posterior_predict(fit, ndraws = 200)
    yrep_mean <- colMeans(yrep)
    yrep_lower <- apply(yrep, 2, quantile, probs = 0.025)
    yrep_upper <- apply(yrep, 2, quantile, probs = 0.975)

    tibble(
      model = model_label,
      n_broods_sampled = dat_main$n_broods_sampled,
      observed = dat_main$n_epp_broods / dat_main$n_broods_sampled,
      predicted = yrep_mean / dat_main$n_broods_sampled,
      lower = yrep_lower / dat_main$n_broods_sampled,
      upper = yrep_upper / dat_main$n_broods_sampled
    )
  }

  main_diagnostics_plot_data <- imap_dfr(
    main_model_files,
    make_model_diagnostics
  )

  fig_main_diagnostics <- ggplot(
    main_diagnostics_plot_data,
    aes(x = predicted, y = observed)
  ) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.14, colour = "#6b8fb5") +
    geom_point(alpha = 0.45, size = 1.2, colour = "#2f6f73") +
    facet_wrap(~ model, ncol = 2) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Posterior predictive mean EPP",
      y = "Observed EPP"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/main_environmental_model_diagnostics.png",
    fig_main_diagnostics,
    width = 9,
    height = 8,
    dpi = 300
  )

  ggsave(
    "figures/main_environmental_model_diagnostics.pdf",
    fig_main_diagnostics,
    width = 9,
    height = 8
  )

  make_funnel_data <- function(model_file, model_label) {
    fit <- readRDS(model_file)

    epred <- brms::posterior_epred(fit, re_formula = NA)
    expected_counts <- colMeans(epred)
    expected <- ifelse(
      expected_counts > 1,
      expected_counts / dat_main$n_broods_sampled,
      expected_counts
    )

    dat_main %>%
      transmute(
        model = model_label,
        n_broods_sampled = n_broods_sampled,
        observed = n_epp_broods / n_broods_sampled,
        expected = expected,
        lower = pmax(0, expected - 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled)),
        upper = pmin(1, expected + 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled))
      )
  }

  main_funnel_plot_data <- imap_dfr(
    main_model_files,
    make_funnel_data
  )

  fig_main_funnels <- ggplot(
    main_funnel_plot_data,
    aes(x = n_broods_sampled, y = observed)
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9fb9bd", alpha = 0.35) +
    geom_line(aes(y = expected), colour = "#b35c35", linewidth = 0.7) +
    geom_point(alpha = 0.45, size = 1.2, colour = "#2f6f73") +
    facet_wrap(~ model, ncol = 2) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      x = "Number of sampled broods (log scale)",
      y = "Observed EPP"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/main_environmental_funnel_plots.png",
    fig_main_funnels,
    width = 9,
    height = 8,
    dpi = 300
  )

  ggsave(
    "figures/main_environmental_funnel_plots.pdf",
    fig_main_funnels,
    width = 9,
    height = 8
  )
} else {
  cat("Main environmental model objects or model data are not available yet.")
}

if (file.exists("figures/main_environmental_model_diagnostics.png")) {
  include_graphics("figures/main_environmental_model_diagnostics.png")
} else {
  cat("Main environmental diagnostic plots are not available yet.")
}

if (file.exists("figures/main_environmental_funnel_plots.png")) {
  include_graphics("figures/main_environmental_funnel_plots.png")
} else {
  cat("Main environmental funnel plots are not available yet.")
}

sensitivity_model_dir <- "models/sensitivity_models"

sensitivity_model_files <- c(
  "Mean-only sensitivity" = file.path(sensitivity_model_dir, "m_sensitivity_first_month_mu.rds"),
  "Mean-scale sensitivity" = file.path(sensitivity_model_dir, "m_sensitivity_first_month_mu_phi.rds")
)

sensitivity_model_files <- sensitivity_model_files[file.exists(sensitivity_model_files)]

sensitivity_data_file <- first_existing(c(
  file.path(sensitivity_model_dir, "dat_sensitivity_models.rds"),
  file.path(sensitivity_model_dir, "dat_fit_sensitivity_models.rds"),
  "data/processed/dat_model_with_environment_geography_pca.rds",
  "data/processed/dat_model_phylo.rds"
))

if (length(sensitivity_model_files) > 0 &&
    !is.na(sensitivity_data_file) &&
    requireNamespace("brms", quietly = TRUE)) {
  dat_sensitivity <- readRDS(sensitivity_data_file)

  make_model_diagnostics <- function(model_file, model_label) {
    fit <- readRDS(model_file)
    model_data <- dat_sensitivity

    yrep <- brms::posterior_predict(fit, ndraws = 200)
    yrep_mean <- colMeans(yrep)
    yrep_lower <- apply(yrep, 2, quantile, probs = 0.025)
    yrep_upper <- apply(yrep, 2, quantile, probs = 0.975)

    tibble(
      model = model_label,
      n_broods_sampled = model_data$n_broods_sampled,
      observed = model_data$n_epp_broods / model_data$n_broods_sampled,
      predicted = yrep_mean / model_data$n_broods_sampled,
      lower = yrep_lower / model_data$n_broods_sampled,
      upper = yrep_upper / model_data$n_broods_sampled
    )
  }

  sensitivity_diagnostics_plot_data <- imap_dfr(
    sensitivity_model_files,
    make_model_diagnostics
  )

  fig_sensitivity_diagnostics <- ggplot(
    sensitivity_diagnostics_plot_data,
    aes(x = predicted, y = observed)
  ) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = "#6b8fb5") +
    geom_point(alpha = 0.45, size = 1.2, colour = "#2f6f73") +
    facet_wrap(~ model, ncol = 1) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = "Posterior predictive mean EPP",
      y = "Observed EPP"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/sensitivity_model_diagnostics.png",
    fig_sensitivity_diagnostics,
    width = 8,
    height = 5,
    dpi = 300
  )

  ggsave(
    "figures/sensitivity_model_diagnostics.pdf",
    fig_sensitivity_diagnostics,
    width = 8,
    height = 5
  )

  make_funnel_data <- function(model_file, model_label) {
    fit <- readRDS(model_file)
    model_data <- dat_sensitivity

    epred <- brms::posterior_epred(fit, re_formula = NA)
    expected_counts <- colMeans(epred)
    expected <- ifelse(
      expected_counts > 1,
      expected_counts / model_data$n_broods_sampled,
      expected_counts
    )

    model_data %>%
      transmute(
        model = model_label,
        n_broods_sampled = n_broods_sampled,
        observed = n_epp_broods / n_broods_sampled,
        expected = expected,
        lower = pmax(0, expected - 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled)),
        upper = pmin(1, expected + 1.96 * sqrt(expected * (1 - expected) / n_broods_sampled))
      )
  }

  sensitivity_funnel_plot_data <- imap_dfr(
    sensitivity_model_files,
    make_funnel_data
  )

  fig_sensitivity_funnels <- ggplot(
    sensitivity_funnel_plot_data,
    aes(x = n_broods_sampled, y = observed)
  ) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#9fb9bd", alpha = 0.35) +
    geom_line(aes(y = expected), colour = "#b35c35", linewidth = 0.7) +
    geom_point(alpha = 0.45, size = 1.2, colour = "#2f6f73") +
    facet_wrap(~ model, ncol = 1) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      x = "Number of sampled broods (log scale)",
      y = "Observed EPP"
    ) +
    theme_minimal(base_size = 11)

  ggsave(
    "figures/sensitivity_funnel_plots.png",
    fig_sensitivity_funnels,
    width = 8,
    height = 5,
    dpi = 300
  )

  ggsave(
    "figures/sensitivity_funnel_plots.pdf",
    fig_sensitivity_funnels,
    width = 8,
    height = 5
  )
} else {
  cat("Sensitivity model objects or model data are not available yet.")
}

if (file.exists("figures/sensitivity_model_diagnostics.png")) {
  include_graphics("figures/sensitivity_model_diagnostics.png")
} else {
  cat("Sensitivity model diagnostic plots are not available yet.")
}

if (file.exists("figures/sensitivity_funnel_plots.png")) {
  include_graphics("figures/sensitivity_funnel_plots.png")
} else {
  cat("Sensitivity funnel plots are not available yet.")
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
  b_phi_tmp_cru_C_first_month_z = "phi: mean temperature",
  b_phi_pre_cru_mm_first_month_log_z = "phi: precipitation (log)",
  b_phi_dtr_cru_C_first_month_z = "phi: DTR",
  b_phi_ndvi_first_month_median_30km_z = "phi: NDVI",
  b_phi_abs_latitude_z = "phi: absolute latitude",
  b_phi_distance_to_coast_z = "phi: distance to coast",
  b_phi_nest_boxyes = "phi: nest box",
  b_phi_marker_typemicrosatellite = "phi: marker: microsatellite",
  b_phi_marker_typeSNP = "phi: marker: SNP",
  b_phi_marker_typeSNPs = "phi: marker: SNPs",
  b_phi_marker_typeother = "phi: marker: other"
)

extract_sensitivity_fixef <- function(model_file, model_label) {
  fit <- readRDS(model_file)
  draws <- posterior::as_draws_df(fit)
  available_terms <- intersect(names(fixed_effect_terms), names(draws))

  map_dfr(available_terms, function(term) {
    x <- draws[[term]]
    tibble(
      model = model_label,
      term = term,
      predictor = fixed_effect_terms[[term]],
      component = ifelse(grepl("^b_phi_", term), "phi", "mu"),
      estimate = median(x),
      lower = quantile(x, 0.025),
      upper = quantile(x, 0.975),
      prob_positive = mean(x > 0),
      prob_negative = mean(x < 0)
    )
  })
}

if (length(sensitivity_model_files) > 0 &&
    requireNamespace("posterior", quietly = TRUE)) {
  sensitivity_fixed_effects <- imap_dfr(
    sensitivity_model_files,
    extract_sensitivity_fixef
  ) %>%
    mutate(
      direction = case_when(
        lower > 0 ~ "Positive",
        upper < 0 ~ "Negative",
        TRUE ~ "Uncertain"
      )
    )

  write_csv(
    sensitivity_fixed_effects,
    file.path(sensitivity_model_dir, "fixed_effects_sensitivity_models.csv")
  )

  predictor_order <- c(
    "First-month mean temperature",
    "First-month precipitation (log)",
    "First-month DTR",
    "First-month NDVI",
    "Absolute latitude",
    "Distance to coast",
    "Nest box",
    "Marker: microsatellite",
    "Marker: SNPs",
    "Marker: other"
  )

  sensitivity_fixed_effects_mean_only <- sensitivity_fixed_effects %>%
    filter(model == "Mean-only sensitivity", component == "mu") %>%
    mutate(
      model = "First-month mean-only",
      predictor = factor(predictor, levels = rev(predictor_order))
    )

  fig_sensitivity_fixed_effects <- ggplot(
    sensitivity_fixed_effects_mean_only,
    aes(x = estimate, y = predictor)
  ) +
    geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey45") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 0.75, colour = "grey45") +
    geom_point(shape = 21, size = 2.5, fill = "white", colour = "grey25", stroke = 0.35) +
    facet_wrap(~ model, ncol = 1) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(colour = "grey25"),
      panel.grid.minor = element_blank()
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
  cat("Sensitivity fixed-effect figure will be generated after the model objects are available.")
}

if (file.exists("figures/sensitivity_fixed_effects_forest_plot.png")) {
  include_graphics("figures/sensitivity_fixed_effects_forest_plot.png")
} else {
  cat("Sensitivity fixed-effect figure is not available yet.")
}
