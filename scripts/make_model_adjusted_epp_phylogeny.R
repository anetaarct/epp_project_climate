suppressPackageStartupMessages({
  library(ape)
  library(brms)
})

model_file <- "models/overall_population/m_overall_population.rds"
data_file <- "models/overall_population/dat_overall_population.rds"
tree_file <- "data/processed/tree_pruned.rds"
phylo_file <- "models/overall_population/phylo_mat_overall_population.rds"
output_png <- "figures/model_adjusted_epp_phylogeny.png"
output_pdf <- "figures/model_adjusted_epp_phylogeny.pdf"
output_csv <- "figures/model_adjusted_epp_by_species.csv"

required_files <- c(model_file, data_file, tree_file, phylo_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Cannot create the phylogenetic EPP figure. Missing: ",
    paste(missing_files, collapse = ", ")
  )
}

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

message("Loading the overall beta-binomial model and phylogeny...")
model <- readRDS(model_file)
dat <- readRDS(data_file)
tree <- readRDS(tree_file)
phylo_mat_fit <- readRDS(phylo_file)

required_columns <- c(
  "species_phylo", "species_nonphylo", "population_id", "family",
  "n_epp_broods", "n_broods_sampled"
)
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0L) {
  stop("Model data lack required columns: ", paste(missing_columns, collapse = ", "))
}

dat$species_phylo <- droplevels(factor(dat$species_phylo))
dat$species_nonphylo <- droplevels(factor(dat$species_nonphylo))
dat$population_id <- droplevels(factor(dat$population_id))

model_species <- levels(dat$species_phylo)
if (!setequal(model_species, tree$tip.label)) {
  stop("The overall-model species and pruned-tree tips are not identical.")
}
if (!all(model_species %in% rownames(phylo_mat_fit))) {
  stop("The fitted phylogenetic covariance matrix lacks model species.")
}

family_values <- lapply(split(as.character(dat$family), dat$species_phylo), function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x) != 1L) {
    stop("Each species must map to exactly one non-missing family.")
  }
  x
})
family_lookup <- unlist(family_values, use.names = TRUE)

# One prediction row per species. Setting trials to one makes posterior_epred()
# return a probability. The population random intercept is deliberately omitted,
# while both species terms are retained. The result is therefore a partially
# pooled species-level EPP frequency, standardized to a zero population effect.
tree <- ape::ladderize(tree, right = FALSE)
species_rows <- match(tree$tip.label, as.character(dat$species_phylo))
newdata <- dat[species_rows, , drop = FALSE]
newdata$species_phylo <- factor(
  tree$tip.label,
  levels = levels(dat$species_phylo)
)
newdata$species_nonphylo <- factor(
  tree$tip.label,
  levels = levels(dat$species_nonphylo)
)
newdata$n_epp_broods <- 0L
newdata$n_broods_sampled <- 1L

set.seed(20260819)
message("Calculating posterior species-level EPP frequencies...")
epred <- brms::posterior_epred(
  model,
  newdata = newdata,
  re_formula = ~
    (1 | species_nonphylo) +
    (1 | gr(species_phylo, cov = phylo_mat_fit)),
  ndraws = 4000
)
epred <- as.matrix(epred)

if (ncol(epred) != nrow(newdata)) {
  stop("Unexpected posterior_epred dimensions for the species prediction grid.")
}

raw_by_species <- aggregate(
  cbind(n_epp_broods, n_broods_sampled) ~ species_phylo,
  data = dat,
  FUN = sum
)
n_records <- table(dat$species_phylo)

species_results <- data.frame(
  species_phylo = tree$tip.label,
  family = unname(family_lookup[tree$tip.label]),
  n_records = as.integer(n_records[tree$tip.label]),
  adjusted_epp_mean = colMeans(epred),
  adjusted_epp_median = apply(epred, 2L, stats::median),
  adjusted_epp_lower_95 = apply(epred, 2L, stats::quantile, probs = 0.025),
  adjusted_epp_upper_95 = apply(epred, 2L, stats::quantile, probs = 0.975),
  stringsAsFactors = FALSE
)

raw_match <- match(species_results$species_phylo, as.character(raw_by_species$species_phylo))
species_results$total_epp_broods <- raw_by_species$n_epp_broods[raw_match]
species_results$total_broods <- raw_by_species$n_broods_sampled[raw_match]
species_results$raw_epp_rate <-
  species_results$total_epp_broods / species_results$total_broods

utils::write.csv(
  species_results,
  output_csv,
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
)

draw_arc <- function(angle_min, angle_max, radius, colour, line_width = 8) {
  angle_seq <- seq(
    angle_min,
    angle_max,
    length.out = max(4L, ceiling((angle_max - angle_min) * 180 / pi))
  )
  graphics::lines(
    radius * cos(angle_seq),
    radius * sin(angle_seq),
    col = colour,
    lwd = line_width,
    lend = 1,
    xpd = NA
  )
}

make_tree_plot <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA)

  tree_radius <- max(ape::node.depth.edgelength(tree)[seq_len(ape::Ntip(tree))])
  plot_limit <- tree_radius * 1.53

  ape::plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    edge.color = "#AEB8BC",
    edge.width = 0.62,
    no.margin = TRUE,
    x.lim = c(-plot_limit, plot_limit),
    y.lim = c(-plot_limit, plot_limit)
  )

  last_plot <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  n_tips <- ape::Ntip(tree)
  tip_x <- last_plot$xx[seq_len(n_tips)]
  tip_y <- last_plot$yy[seq_len(n_tips)]
  tip_radius <- sqrt(tip_x^2 + tip_y^2)
  radius <- max(tip_radius)
  tip_angle <- atan2(tip_y, tip_x)
  tip_angle_2pi <- ifelse(tip_angle < 0, tip_angle + 2 * pi, tip_angle)

  # The outer band contains every family. Alternating neutral tones preserve the
  # continuous EPP colour scale as the only quantitative colour encoding.
  family_centres <- vapply(
    split(tip_angle, species_results$family),
    function(x) atan2(mean(sin(x)), mean(cos(x))),
    numeric(1)
  )
  family_centre_2pi <- ifelse(family_centres < 0, family_centres + 2 * pi, family_centres)
  family_order <- names(sort(family_centre_2pi))
  family_band_colours <- stats::setNames(
    rep(c("#D9DEE0", "#AAB7BA"), length.out = length(family_order)),
    family_order
  )

  angle_order <- order(tip_angle_2pi)
  ordered_angles <- tip_angle_2pi[angle_order]
  previous_angles <- c(ordered_angles[length(ordered_angles)] - 2 * pi,
                       ordered_angles[-length(ordered_angles)])
  next_angles <- c(ordered_angles[-1], ordered_angles[1] + 2 * pi)
  lower_bounds <- (previous_angles + ordered_angles) / 2
  upper_bounds <- (ordered_angles + next_angles) / 2
  ring_radius <- radius * 1.085

  for (i in seq_along(angle_order)) {
    tip_index <- angle_order[i]
    family_i <- species_results$family[tip_index]
    draw_arc(
      lower_bounds[i],
      upper_bounds[i],
      ring_radius,
      family_band_colours[[family_i]],
      line_width = 18
    )
  }

  palette <- viridisLite::viridis(256, option = "C")
  colour_index <- pmax(
    1L,
    pmin(256L, floor(species_results$adjusted_epp_mean * 255) + 1L)
  )
  tip_colours <- palette[colour_index]

  graphics::points(
    tip_x,
    tip_y,
    pch = 21,
    bg = tip_colours,
    col = "white",
    cex = 1.25,
    lwd = 0.45
  )

  family_counts <- sort(table(species_results$family), decreasing = TRUE)
  labelled_families <- names(family_counts[family_counts >= 5L])

  for (family_i in labelled_families) {
    angle_i <- family_centres[[family_i]]
    angle_deg <- angle_i * 180 / pi
    rotation <- angle_deg - 90
    while (rotation > 90) rotation <- rotation - 180
    while (rotation < -90) rotation <- rotation + 180

    graphics::text(
      ring_radius * cos(angle_i),
      ring_radius * sin(angle_i),
      labels = family_i,
      srt = rotation,
      adj = c(0.5, 0.5),
      cex = 0.57,
      col = "white",
      font = 2,
      xpd = NA
    )
    graphics::text(
      ring_radius * cos(angle_i),
      ring_radius * sin(angle_i),
      labels = family_i,
      srt = rotation,
      adj = c(0.5, 0.5),
      cex = 0.51,
      col = "#29383B",
      font = 2,
      xpd = NA
    )
  }

  # Fixed 0-100% scale: colours have the same interpretation across figures.
  legend_y_bottom <- -radius * 1.43
  legend_y_top <- -radius * 1.385
  legend_x <- seq(-radius * 0.43, radius * 0.43, length.out = 257)
  for (i in seq_len(256L)) {
    graphics::rect(
      legend_x[i], legend_y_bottom,
      legend_x[i + 1L], legend_y_top,
      col = palette[i], border = NA, xpd = NA
    )
  }
  tick_values <- seq(0, 1, by = 0.25)
  tick_x <- -radius * 0.43 + tick_values * radius * 0.86
  graphics::segments(
    tick_x, legend_y_bottom,
    tick_x, legend_y_bottom - radius * 0.012,
    col = "#29383B", lwd = 0.7, xpd = NA
  )
  graphics::text(
    tick_x,
    legend_y_bottom - radius * 0.035,
    labels = paste0(tick_values * 100, "%"),
    cex = 0.66,
    col = "#29383B",
    xpd = NA
  )
  graphics::text(
    0,
    legend_y_top + radius * 0.035,
    labels = "Model-adjusted EPP frequency",
    cex = 0.76,
    col = "#29383B",
    font = 2,
    xpd = NA
  )
}

message("Writing figure files...")
grDevices::png(output_png, width = 3600, height = 3600, res = 400)
make_tree_plot()
grDevices::dev.off()

grDevices::pdf(output_pdf, width = 9, height = 9, useDingbats = FALSE)
make_tree_plot()
grDevices::dev.off()

message(
  "Done: ", nrow(species_results), " species; ",
  length(unique(species_results$family)), " families; ",
  sum(table(species_results$family) >= 5L), " labelled families."
)
