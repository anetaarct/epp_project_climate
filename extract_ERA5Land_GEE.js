// =====================================================
// EPP climate extraction from ERA5-Land Daily Aggregated
// 30-km buffer, three breeding windows
// =====================================================

var points = ee.FeatureCollection(
  'projects/epp-project-498119/assets/epp_climate_input_3windows'
);

var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
  .select([
    'temperature_2m',
    'temperature_2m_min',
    'temperature_2m_max',
    'total_precipitation_sum'
  ]);

var extractClimate = function(feature) {
  var year = ee.Number(feature.get('year'));
  var startMonth = ee.Number(feature.get('breeding_start_month'));
  var endMonth = ee.Number(feature.get('breeding_end_month'));
  var previousMonth = ee.Number(feature.get('previous_month'));
  var previousYear = ee.Number(feature.get('previous_month_year'));
  var window = ee.String(feature.get('climate_window'));

  var startDate = ee.Date.fromYMD(year, startMonth, 1);
  var endDate = startDate.advance(1, 'month');

  var prevDate = ee.Date.fromYMD(previousYear, previousMonth, 1);
  var prevEndDate = startDate.advance(1, 'month');

  var seasonEndDate = ee.Date.fromYMD(year, endMonth, 1)
    .advance(1, 'month');

  var dateStart = ee.Date(
    ee.Algorithms.If(
      window.equals('start_only'),
      startDate,
      ee.Algorithms.If(
        window.equals('previous_plus_start'),
        prevDate,
        startDate
      )
    )
  );

  var dateEnd = ee.Date(
    ee.Algorithms.If(
      window.equals('start_only'),
      endDate,
      ee.Algorithms.If(
        window.equals('previous_plus_start'),
        prevEndDate,
        seasonEndDate
      )
    )
  );

  var point = ee.Geometry.Point([
    ee.Number(feature.get('long')),
    ee.Number(feature.get('lat'))
  ]);

  var buffer = point.buffer(30000);

  var daily = era5.filterDate(dateStart, dateEnd);

  var tempMean = daily.select('temperature_2m').mean()
    .subtract(273.15)
    .rename('temp_mean_C');

  var tempSd = daily.select('temperature_2m').reduce(ee.Reducer.stdDev())
    .rename('temp_sd_daily_K');

  var tminMean = daily.select('temperature_2m_min').mean()
    .subtract(273.15)
    .rename('tmin_mean_C');

  var tmaxMean = daily.select('temperature_2m_max').mean()
    .subtract(273.15)
    .rename('tmax_mean_C');

  var dtrMean = daily.map(function(img) {
      return img.select('temperature_2m_max')
        .subtract(img.select('temperature_2m_min'))
        .rename('dtr_daily');
    })
    .mean()
    .rename('dtr_mean_C');

  var precipTotal = daily.select('total_precipitation_sum').sum()
    .multiply(1000)
    .rename('precip_total_mm');

  var precipSd = daily.select('total_precipitation_sum')
    .map(function(img) {
      return img.multiply(1000);
    })
    .reduce(ee.Reducer.stdDev())
    .rename('precip_sd_daily_mm');

  var rainDays = daily.select('total_precipitation_sum')
    .map(function(img) {
      return img.multiply(1000).gt(1).rename('rain_day');
    })
    .sum()
    .rename('rain_days_gt1mm');

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
    scale: 10000,
    maxPixels: 1e9
  });

  return feature.set(values)
    .set('buffer_km', 30)
    .set('date_start', dateStart.format('YYYY-MM-dd'))
    .set('date_end', dateEnd.format('YYYY-MM-dd'));
};

var climateExtracted = points.map(extractClimate);

print('Input records:', points.size());
print('Extracted climate:', climateExtracted.limit(10));

Export.table.toDrive({
  collection: climateExtracted,
  description: 'epp_era5land_climate_30km_3windows',
  fileFormat: 'CSV'
});
