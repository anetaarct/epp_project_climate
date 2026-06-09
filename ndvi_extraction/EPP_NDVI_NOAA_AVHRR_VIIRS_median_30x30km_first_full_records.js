var records = ee.FeatureCollection(
  'projects/epp-project-498119/assets/EPP_NDVI_GEE_input_all_records'
);

// 30 x 30 km square
var halfSideMeters = 15000;

var avhrr = ee.ImageCollection('NOAA/CDR/AVHRR/NDVI/V5');
var viirs = ee.ImageCollection('NOAA/CDR/VIIRS/NDVI/V1');

function prepNdvi(img) {
  var qa = img.select('QA');

  var cloud = qa.bitwiseAnd(1 << 1).eq(0);
  var shadow = qa.bitwiseAnd(1 << 2).eq(0);
  var water = qa.bitwiseAnd(1 << 3).eq(0);
  var sunglint = qa.bitwiseAnd(1 << 4).eq(0);
  var night = qa.bitwiseAnd(1 << 6).eq(0);
  var channelsValid = qa.bitwiseAnd(1 << 7).neq(0);

  var mask = cloud
    .and(shadow)
    .and(water)
    .and(sunglint)
    .and(night)
    .and(channelsValid);

  return img.select('NDVI')
    .multiply(0.0001)
    .rename('ndvi')
    .updateMask(mask)
    .copyProperties(img, ['system:time_start']);
}

var ndvi = avhrr.map(prepNdvi).merge(viirs.map(prepNdvi));

function medianNdviForPeriod(startDate, endDate) {
  return ndvi.filterDate(startDate, endDate).median().rename('ndvi');
}

function seasonEndDate(year, startMonth, endMonth) {
  var endYear = ee.Number(
    ee.Algorithms.If(ee.Number(endMonth).lt(startMonth), ee.Number(year).add(1), year)
  );
  return ee.Date.fromYMD(endYear, endMonth, 1).advance(1, 'month');
}

function getNdviMean(img, geom) {
  return img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom,
    scale: 5566,
    bestEffort: true,
    maxPixels: 1e8
  }).get('ndvi');
}

function addFirstNdvi(feature) {
  var year = ee.Number.parse(feature.get('year'));
  var startMonth = ee.Number.parse(feature.get('start_month'));
  var lon = ee.Number.parse(feature.get('lon'));
  var lat = ee.Number.parse(feature.get('lat'));

  var geom = ee.Geometry.Point([lon, lat])
    .buffer(halfSideMeters)
    .bounds();

  var firstStart = ee.Date.fromYMD(year, startMonth, 1);
  var firstEnd = firstStart.advance(1, 'month');

  var firstImg = medianNdviForPeriod(firstStart, firstEnd);
  var ndviFirst = getNdviMean(firstImg, geom);

  return feature.set('ndvi_first_month_median_30km', ndviFirst);
}

function addFullNdvi(feature) {
  var canFull = ee.Number.parse(feature.get('can_full_season'));

  return ee.Feature(
    ee.Algorithms.If(
      canFull.eq(1),
      addFullNdviReal(feature),
      feature.set('ndvi_full_season_median_30km', null)
    )
  );
}

function addFullNdviReal(feature) {
  var year = ee.Number.parse(feature.get('year'));
  var startMonth = ee.Number.parse(feature.get('start_month'));
  var endMonth = ee.Number.parse(feature.get('end_month'));
  var lon = ee.Number.parse(feature.get('lon'));
  var lat = ee.Number.parse(feature.get('lat'));

  var geom = ee.Geometry.Point([lon, lat])
    .buffer(halfSideMeters)
    .bounds();

  var fullStart = ee.Date.fromYMD(year, startMonth, 1);
  var fullEnd = seasonEndDate(year, startMonth, endMonth);

  var fullImg = medianNdviForPeriod(fullStart, fullEnd);
  var ndviFull = getNdviMean(fullImg, geom);

  return feature.set('ndvi_full_season_median_30km', ndviFull);
}

var recordsFirst = records
  .filter(ee.Filter.eq('can_first_month', 1));

var result = recordsFirst
  .map(addFirstNdvi)
  .map(addFullNdvi);

print('Count input', records.size());
print('Count first ready', recordsFirst.size());
print('Result count', result.size());
print('Preview result', result.limit(5));

Export.table.toDrive({
  collection: result,
  description: 'EPP_NOAA_NDVI_median_30x30km_first_full',
  fileFormat: 'CSV'
});
