#!/usr/bin/env Rscript

# Reconstruct and audit distance_to_coast_km.
#
# Usage from the repository root:
#   Rscript scripts/reconstruct_distance_to_coast.R
#
# Optional arguments:
#   1. input CSV (default: epp_data_June2026.csv)
#   2. output directory (default: data/processed)
#
# This script never overwrites the analysis-ready dataset. It reconstructs the
# geographic covariate independently and stops if it fails to match the stored
# values within 0.001 km (1 m).

required_packages <- c("sf", "s2", "rnaturalearthdata")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this reconstruction."
  )
}

required_naturalearthdata_version <- "1.0.0"
installed_naturalearthdata_version <- as.character(
  utils::packageVersion("rnaturalearthdata")
)
if (!identical(
  installed_naturalearthdata_version,
  required_naturalearthdata_version
)) {
  stop(
    "This audit is pinned to rnaturalearthdata ",
    required_naturalearthdata_version,
    "; installed version is ", installed_naturalearthdata_version, "."
  )
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1]] else "epp_data_June2026.csv"
output_dir <- if (length(args) >= 2L) args[[2]] else "data/processed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- utils::read.csv(input_file, check.names = FALSE)
required_columns <- c("location_id", "lat", "long", "distance_to_coast_km")
missing_columns <- setdiff(required_columns, names(d))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

valid_rows <- stats::complete.cases(d[required_columns])
locations <- unique(d[valid_rows, required_columns])
locations <- locations[order(locations$location_id), ]

coordinate_keys <- paste(locations$lat, locations$long, sep = "|")
stored_values_per_coordinate <- tapply(
  locations$distance_to_coast_km,
  coordinate_keys,
  function(x) length(unique(x))
)
if (any(stored_values_per_coordinate != 1L)) {
  stop("At least one coordinate has multiple stored distance-to-coast values.")
}
if (anyDuplicated(locations[c("lat", "long")])) {
  stop("Location identifiers do not map one-to-one to unique coordinates.")
}

points <- sf::st_as_sf(
  locations,
  coords = c("long", "lat"),
  crs = 4326,
  remove = FALSE
)
coastline <- rnaturalearthdata::coastline50

old_s2 <- sf::sf_use_s2(TRUE)
on.exit(sf::sf_use_s2(old_s2), add = TRUE)

distance_matrix_m <- sf::st_distance(points, coastline)
reconstructed_km <- apply(
  distance_matrix_m,
  1,
  function(x) min(as.numeric(x), na.rm = TRUE)
) / 1000

audit <- data.frame(
  location_id = locations$location_id,
  lat = locations$lat,
  long = locations$long,
  stored_distance_to_coast_km = locations$distance_to_coast_km,
  reconstructed_distance_to_coast_km = reconstructed_km,
  difference_km = reconstructed_km - locations$distance_to_coast_km,
  absolute_difference_km = abs(
    reconstructed_km - locations$distance_to_coast_km
  )
)

tolerance_km <- 0.001
summary_table <- data.frame(
  source = "Natural Earth 1:50m coastline",
  object = "rnaturalearthdata::coastline50",
  rnaturalearthdata_version = installed_naturalearthdata_version,
  sf_version = as.character(utils::packageVersion("sf")),
  s2_version = as.character(utils::packageVersion("s2")),
  coordinate_crs = "EPSG:4326",
  distance_engine = "sf::st_distance with sf_use_s2(TRUE)",
  n_dataset_rows = nrow(d),
  n_rows_with_complete_coordinates_and_distance = sum(valid_rows),
  n_unique_locations_compared = nrow(audit),
  mean_absolute_difference_km = mean(audit$absolute_difference_km),
  maximum_absolute_difference_km = max(audit$absolute_difference_km),
  n_difference_gt_0_001_km = sum(audit$absolute_difference_km > 0.001),
  n_difference_gt_0_01_km = sum(audit$absolute_difference_km > 0.01),
  hard_stop_tolerance_km = tolerance_km,
  validation_passed = max(audit$absolute_difference_km) <= tolerance_km
)

utils::write.csv(
  audit,
  file.path(output_dir, "distance_to_coast_reconstruction_audit.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  summary_table,
  file.path(output_dir, "distance_to_coast_reconstruction_summary.csv"),
  row.names = FALSE,
  na = ""
)

if (!summary_table$validation_passed) {
  stop(
    "Distance-to-coast reconstruction failed: maximum absolute difference = ",
    format(summary_table$maximum_absolute_difference_km, scientific = TRUE),
    " km, exceeding the ", tolerance_km, " km tolerance."
  )
}

cat("Distance-to-coast reconstruction passed.\n")
cat("Unique locations compared:", nrow(audit), "\n")
cat(
  "Maximum absolute difference (km):",
  format(summary_table$maximum_absolute_difference_km, scientific = TRUE),
  "\n"
)
cat("Outputs written to:", normalizePath(output_dir), "\n")
