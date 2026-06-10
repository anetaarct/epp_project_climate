# Daily weather variability extraction workflow

This folder contains the files used to extract daily weather variability metrics for the EPP dataset using Google Earth Engine.

## Input file

`epp_locations_for_gee_daily_weather.csv`

This file contains one row per EPP record and was uploaded as a table asset in Google Earth Engine. It includes:

* `epp_row_id`: unique record identifier used to merge exported weather variables back to the main dataset
* `year`: study year
* `first_month`: first month of the breeding season
* `latitude`, `longitude`: study coordinates used for environmental extraction

## Earth Engine script

`gee_era5land_daily_weather_variability_first_month.js`

The script extracts daily weather data from:

* `ECMWF/ERA5_LAND/DAILY_AGGR`

For each record, the script:

1. Creates a point geometry from the study coordinates.
2. Identifies the first breeding month based on the study year and breeding phenology.
3. Extracts daily temperature and precipitation values for the entire first-month breeding window.
4. Calculates summary statistics describing within-month weather conditions and variability.

Temperature values are converted from Kelvin to degrees Celsius, and precipitation values are converted from meters to millimetres.

## Extracted variables

### Temperature

* `tmp_mean_C_first_month`: mean daily temperature during the first breeding month
* `tmp_sd_C_first_month`: standard deviation of daily temperature during the first breeding month
* `tmp_range_C_first_month`: range of daily temperature values during the first breeding month

### Precipitation

* `pre_sum_mm_first_month`: total precipitation during the first breeding month
* `pre_mean_mm_day_first_month`: mean daily precipitation during the first breeding month
* `pre_sd_mm_day_first_month`: standard deviation of daily precipitation during the first breeding month
* `pre_max_mm_day_first_month`: maximum daily precipitation during the first breeding month
* `pre_wet_days_gt1mm_first_month`: number of days with precipitation greater than 1 mm during the first breeding month

### Metadata

* `weather_start_date`: first day of the extraction window
* `weather_end_date`: last day of the extraction window
* `weather_dataset`: source dataset used for extraction

## Spatial extraction

Environmental values were extracted directly at the study coordinates using the native ERA5-Land spatial resolution (approximately 11 km).

## Output file

`epp_era5_land_daily_variability_first_month.csv`

This file contains the original record identifiers together with the extracted weather metrics.

The output can be merged back to the main EPP dataset using:

```r
left_join(main_data,
          weather_data,
          by = "epp_row_id")
```

## Variables used in the final models

The daily weather variability model uses the following variables from the extraction:

* `tmp_sd_C_first_month`
* `pre_sd_mm_day_first_month`
* `pre_wet_days_gt1mm_first_month`

These variables quantify temperature variability, precipitation variability, and wet-day frequency during the first breeding month.

