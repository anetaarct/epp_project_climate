# NDVI extraction workflow

This folder contains the files used to extract NDVI values for the EPP mixed-paternity dataset using Google Earth Engine.

## Input file

`EPP_NDVI_GEE_input_all_records.csv`

This file contains one row per study record and was used as the input table asset in Google Earth Engine. It includes:

- `record_id`: unique record identifier used to merge exported NDVI values back to the main dataset
- `id`: study identifier
- `scientific_name`: species name
- `year`: year of study
- `start_month`: first month of the breeding season
- `end_month`: last month of the breeding season
- `can_first_month`: whether first-month NDVI can be calculated
- `can_full_season`: whether full-season NDVI can be calculated
- `lon`, `lat`: rounded geographic coordinates used for extraction

## Earth Engine script

`EPP_NDVI_NOAA_AVHRR_VIIRS_median_30x30km_first_full_records.js`

The script extracts NDVI from daily NOAA CDR datasets:

- `NOAA/CDR/AVHRR/NDVI/V5` for 1981-2013
- `NOAA/CDR/VIIRS/NDVI/V1` for 2014 onward

For each record, the script creates a 30 × 30 km square around the study coordinates. Pixels flagged as cloud, cloud shadow, water, sunglint, night, or invalid channel data are masked.

NDVI was calculated in two windows:

- `ndvi_first_month_median_30km`: median NDVI during the first breeding month, calculated per pixel and then averaged across the 30 × 30 km square
- `ndvi_full_season_median_30km`: median NDVI across the full breeding season, calculated per pixel and then averaged across the 30 × 30 km square

## Output file

`EPP_NOAA_NDVI_median_30x30km_first_full.csv`

This is the exported Google Earth Engine output. It contains the original record identifiers and the extracted NDVI variables.

The output can be joined back to the main EPP dataset using:

```r
left_join(main_data, ndvi_data, by = "record_id")
