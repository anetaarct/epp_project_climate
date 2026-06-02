// =====================================================
// EPP climate extraction from ERA5 DAILY global
// For comparison with ERA5-Land and gap filling
// 30-km buffer, two windows
// =====================================================

var points = ee.FeatureCollection(
  'projects/epp-project-498119/assets/epp_climate_input'
);

var era5 = ee.ImageCollection('ECMWF/ERA5/DAILY')
  .select([
    'mean_2m_air_temperature',
    'minimum_2m_air_temperature',
    'maximum_2m_air_temperature',
    'total_precipitation'
  ]);

var extractClimate = function(feature) {
  var year = ee.Number(feature.get('year'));
  var startMonth = ee.Number(feature.get('breeding_start_month'));
  var previousMonth = ee.Number(feature.get('previous_month'));
  var previousYear = ee.Number(feature.get('previous_month_year'));
  var window = ee.String(feature.get('climate_window'));

  var startDate = ee.Date.fromYMD(year, startMonth, 1);
  var endDate = startDate.advance(1, 'month');

  var prevDate = ee.Date.fromYMD(previousYear, previousMonth, 1);
  var prevEndDate = startDate.advance(1, 'month');

  var dateStart = ee.Date(
    ee.Algorithms.If(window.equals('start_only'), startDate, prevDate)
  );

  var dateEnd = ee.Date(
    ee.Algorithms.If(window.equals('start_only'), endDate, prevEndDate)
  );

  var point = ee.Geometry.Point([
    ee.Number(feature.get('long')),
    ee.Number(feature.get('lat'))
  ]);

  var buffer = point.buffer(30000);

  var daily = era5.filterDate(dateStart, dateEnd);

  var tempMean = daily.select('mean_2m_air_temperature').mean()
    .subtract(273.15)
    .rename('era5_temp_mean_C');

  var tempSd = daily.select('mean_2m_air_temperature')
    .reduce(ee.Reducer.stdDev())
    .rename('era5_temp_sd_daily_K');

  var tminMean = daily.select('minimum_2m_air_temperature').mean()
    .subtract(273.15)
    .rename('era5_tmin_mean_C');

  var tmaxMean = daily.select('maximum_2m_air_temperature').mean()
    .subtract(273.15)
    .rename('era5_tmax_mean_C');

  var dtrMean = daily.map(function(img) {
      return img.select('maximum_2m_air_temperature')
        .subtract(img.select('minimum_2m_air_temperature'))
        .rename('dtr_daily');
    })
    .mean()
    .rename('era5_dtr_mean_C');

  var precipTotal = daily.select('total_precipitation').sum()
    .multiply(1000)
    .rename('era5_precip_total_mm');

  var precipSd = daily.select('total_precipitation')
    .map(function(img) {
      return img.multiply(1000);
    })
    .reduce(ee.Reducer.stdDev())
    .rename('era5_precip_sd_daily_mm');

  var rainDays = daily.select('total_precipitation')
    .map(function(img) {
      return img.multiply(1000).gt(1).rename('rain_day');
    })
    .sum()
    .rename('era5_rain_days_gt1mm');

  var climateImage = tempMean
    .addBands(tempSd)
    .addBands(tminMean)
    .addBands(tmaxMean)
    .addBands(dtrMean)
    .addBands(precipTotal)
    .addBands(precipSd)
    .addBands(rainDays);

  var values = climateImage.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: buffer,
    scale: 30000,
    maxPixels: 1e9
  });

  return feature.set(values)
    .set('source_gapfill', 'ERA5_DAILY_global')
    .set('buffer_km', 30)
    .set('date_start', dateStart.format('YYYY-MM-dd'))
    .set('date_end', dateEnd.format('YYYY-MM-dd'));
};

var climateERA5 = points.map(extractClimate);

print('Input records:', points.size());
print('ERA5 global preview:', climateERA5.limit(10));

Export.table.toDrive({
  collection: climateERA5,
  description: 'epp_era5_global_climate_30km',
  fileFormat: 'CSV'
});
