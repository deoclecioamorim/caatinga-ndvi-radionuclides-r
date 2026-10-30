/**
 * NDVI v2.0.0 | Bimonthly Landsat NDVI composites, Areas A and B, 2000-2022 | 8 October 2026
 * Paste this complete file into the Google Earth Engine Code Editor.
 * Default: 276 S1/OLI-reference GeoTIFF tasks and nine CSV tasks, 2000-2022.
 * The GeoTIFFs of run 20261008_v2_bimonthly_updated_areas (scenario S1, band NDVI_harmOLI)
 * are the inputs of the R scripts (data/raw/ndvi_area_A and data/raw/ndvi_area_B).
 * Tables include S1/S2, common-support comparisons, sensor overlaps and inventory.
 * All products use the same configured polygons, grid and bimonthly processing.
 * Masked pixels become declared NoData (-9999), never measured NDVI equal to zero.
 * No gap filling or minimum-coverage exclusion; all calendar positions are retained.
 * Study area polygons: set CONFIG.areaAssets to your own Earth Engine assets; the polygons
 * are available from the corresponding author upon reasonable request.
 * Roy et al. (2016), Table 3: doi:10.1016/j.rse.2015.12.024.
 */

// USER SETTINGS: use separate switches to export full tables and smaller raster batches.
var CONFIG = {
  driveFolder: 'data_earth_SEMINE_bimonthly_updated_areas',
  runId: '20261008_v2_bimonthly_updated_areas',
  exportTables: true, exportRasters: true,
  exportAreas: ['A', 'B'],
  firstYear: 2000, lastYear: 2022,
  rasterScenario: 'S1', rasterBand: 'NDVI_harmOLI',
  areaAssets: {
    A: 'projects/<your-project>/assets/area_A', // polygon of Area A (not distributed)
    B: 'projects/<your-project>/assets/area_B'  // polygon of Area B (not distributed)
  },

  // STUDY SETTINGS: processing retained from the executed bimonthly pipeline.
  version: 'ndvi-v2.0.0',
  scriptRevision: 'bimonthly-complete-updated-areas-20261008',
  scriptFile: 'gee_ndvi_v2_bimonthly_complete.js',
  sourceRunId: '20261007_v2_primary', // Code lineage, not the source of new polygon results.
  sourceScriptSha256: 'eda170ff61ddb6742c1f1ffdf2f1a04865b935bd58529562f5c0f3fb4990057e',
  parentScriptSha256: '046f7d8180f7f9ecc49e95e64f11ffced5a7846ce1f18284c56e6933670aaf4f',
  processingSourceScriptSha256: 'a7ac2241791b9a496533b5db1d204b7281ffd09d34bec171873ab994d9aa82e2',
  canonicalSha256: '4b2d450c2adb8d5c3b00ef6fc8539421cc9c2ef160b3b488a821a4b1c808ebb4',
  start: '2000-01-01', endExclusive: '2023-01-01',
  periodMonths: 2, frequencyPerYear: 6,
  crs: 'EPSG:32724', crsTransform: [30, 0, 0, 0, -30, 10000000],
  tileScale: 4, maxPixels: 1e8, denominatorEpsilon: 1e-6,
  medianMaxRaw: 1000, noData: -9999,
  coverageThresholds: [0, 0.3, 0.5, 0.8],
  overlapMinCentroidPixels: 100,
  enabledScenarios: ['S1', 'S2']
};
if (CONFIG.periodMonths !== 2 || CONFIG.frequencyPerYear !== 6 ||
    CONFIG.start !== '2000-01-01' || CONFIG.endExclusive !== '2023-01-01') {
  throw new Error('This edition requires the fixed 2000-2022 bimonthly calendar.');
}
if (CONFIG.firstYear % 1 !== 0 || CONFIG.lastYear % 1 !== 0 ||
    CONFIG.firstYear < 2000 || CONFIG.lastYear > 2022 || CONFIG.firstYear > CONFIG.lastYear) {
  throw new Error('Choose whole years with 2000 <= firstYear <= lastYear <= 2022.');
}
if (!CONFIG.exportAreas.length || CONFIG.exportAreas.some(function(key, i) {
  return ['A','B'].indexOf(key) === -1 || CONFIG.exportAreas.indexOf(key) !== i;
})) {
  throw new Error('exportAreas must contain A, B or both, without duplicates.');
}
if (typeof CONFIG.exportTables !== 'boolean' || typeof CONFIG.exportRasters !== 'boolean') {
  throw new Error('exportTables and exportRasters must be true or false.');
}
if (CONFIG.noData !== -9999 || CONFIG.enabledScenarios.join(',') !== 'S1,S2') {
  throw new Error('Keep NoData = -9999 and the S1/S2 table scenarios for this edition.');
}

// Roy et al. (2016), Table 3: SURFACE NDVI, two separate OLS regressions.
// CONUS/pre-Collection-2/ETM+-OLI calibration is a transfer hypothesis here.
var HARM = {ETM_TO_OLI: {a: 0.0235, b: 0.9723}, OLI_TO_ETM: {a: 0.0029, b: 0.9589}};
var SENSORS = {
  L5: {id: 'LANDSAT/LT05/C02/T1_L2', red: 'SR_B3', nir: 'SR_B4', oli: false, radBits: 524},
  L7: {id: 'LANDSAT/LE07/C02/T1_L2', red: 'SR_B3', nir: 'SR_B4', oli: false, radBits: 524},
  L8: {id: 'LANDSAT/LC08/C02/T1_L2', red: 'SR_B4', nir: 'SR_B5', oli: true, radBits: 2072},
  L9: {id: 'LANDSAT/LC09/C02/T1_L2', red: 'SR_B4', nir: 'SR_B5', oli: true, radBits: 2072}
};
var SCENARIOS = [
  {id: 'S1', sensors: ['L5','L7','L8','L9'], l7End: null, strictSR: true, familyTransfer: true},
  // Historical sensitivity cutoff; March-April 2013 remains a partial transition bin.
  // 11 April is not a claim about the operational start of Landsat 8.
  {id: 'S2', sensors: ['L5','L7','L8','L9'], l7End: '2013-04-11', strictSR: true, familyTransfer: true}
];
var BANDS = ['NDVI_raw','NDVI_harmOLI','NDVI_harmETM'];
var PROJ = ee.Projection(CONFIG.crs, CONFIG.crsTransform);
var AREAS = CONFIG.exportAreas.map(function(a) {
  var asset = CONFIG.areaAssets[a];
  var fc = ee.FeatureCollection(asset);
  return {name: 'Area_' + a, asset: asset, fc: fc, geom: fc.geometry()};
});
var ROI = AREAS[0].geom;
AREAS.slice(1).forEach(function(a) {ROI = ROI.union(a.geom, ee.ErrorMargin(0.1));});
// GRID AND SPATIAL SUPPORT
function grid(image) { return ee.Image(image).reproject({crs: PROJ}); }
function zero(name) { return grid(ee.Image.constant(0).rename(name)); }
function empty(bands) { return grid(ee.Image.constant(bands.map(function() { return 0; })).rename(bands).toFloat()).updateMask(zero('mask')); }
function args(area, reducer) {
  return {geometry: area.geom, reducer: reducer, crs: CONFIG.crs,
    crsTransform: CONFIG.crsTransform, maxPixels: CONFIG.maxPixels,
    tileScale: CONFIG.tileScale, bestEffort: false};
}
function safeNumber(value) { return ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(value, null), 0, value)); }
function safeRatio(a, b) { return ee.Algorithms.If(ee.Number(b).gt(0), ee.Number(a).divide(ee.Number(b).max(1e-12)), null); }
function sum(image, area) { return safeNumber(ee.Image(image).rename('v').reduceRegion(args(area, ee.Reducer.sum())).get('v')); }
var PIXEL_AREA = grid(ee.Image.pixelArea());
function areaTotal(area) { return sum(PIXEL_AREA, area); }
function validMask(image) { return grid(ee.Image(image).mask().reduce(ee.Reducer.min()).unmask(0, false).gt(0)); }
function supportArea(mask, area) { return sum(PIXEL_AREA.multiply(ee.Image(mask).unmask(0, false)), area); }
function areaMean(image, area, mask) {
  var m = mask || validMask(image);
  var w = PIXEL_AREA.updateMask(m);
  return safeRatio(sum(ee.Image(image).multiply(w), area), sum(w, area));
}
// SENSOR PREPROCESSING
function qaBits(key) { return SENSORS[key].oli ? 63 : 59; }
function preprocess(key, scenario) {
  var sensor = SENSORS[key];
  return function(image) {
    var base = image.select('QA_PIXEL').bitwiseAnd(qaBits(key)).eq(0)
      .and(image.select('QA_RADSAT').bitwiseAnd(sensor.radBits).eq(0));
    var red = image.select(sensor.red).multiply(0.0000275).add(-0.2);
    var nir = image.select(sensor.nir).multiply(0.0000275).add(-0.2);
    var denom = nir.add(red);
    var strict = red.gt(0).and(red.lte(1)).and(nir.gt(0)).and(nir.lte(1));
    var goodDenom = denom.abs().gt(CONFIG.denominatorEpsilon);
    var raw = nir.subtract(red).divide(denom).updateMask(base.and(goodDenom));
    if (scenario.strictSR) raw = raw.updateMask(strict);
    raw = raw.rename('NDVI_raw');
    var domain = raw.gte(0).and(raw.lte(1));
    var applyETM = !sensor.oli && (key === 'L7' || scenario.familyTransfer);
    var applyOLI = sensor.oli && (key === 'L8' || scenario.familyTransfer);
    // Outside calibration domain, retain raw (do not extrapolate, clip or drop it).
    var oli = applyETM ? raw.where(domain, raw.multiply(HARM.ETM_TO_OLI.b).add(HARM.ETM_TO_OLI.a)) : raw;
    var etm = applyOLI ? raw.where(domain, raw.multiply(HARM.OLI_TO_ETM.b).add(HARM.OLI_TO_ETM.a)) : raw;
    var qa = raw.multiply(0).add(1).rename('accepted');
    var preRange = red.multiply(0).add(1).updateMask(base.and(goodDenom)).rename('qa_denom');
    var reject = preRange.updateMask(strict.not()).rename('range_reject');
    var outside = qa.updateMask(domain.not()).rename('outside_domain');
    var acquisitionDay=raw.multiply(0).add(ee.Date(ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'))
      .difference(ee.Date(CONFIG.start),'day')).rename('acq_day');
    // All bands use nearest-neighbor reprojection before mosaicking/compositing.
    return grid(ee.Image.cat([raw, oli.rename(BANDS[1]), etm.rename(BANDS[2]), qa, preRange, reject, outside,acquisitionDay]).toFloat())
      .set({day: ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
        scene_id: image.get('system:index'), 'system:time_start': image.get('system:time_start')});
  };
}
// SENSOR-DAY OBSERVATIONS AND BIMONTHLY COMPOSITES
function rawCollection(key, area, start, end, scenario) {
  var c = ee.ImageCollection(SENSORS[key].id).filterBounds(area.geom).filterDate(start,end);
  if (key === 'L7' && scenario.l7End) c = c.filterDate(CONFIG.start, scenario.l7End);
  return c.sort('system:index');
}
function dailyCollection(key, area, start, end, scenario) {
  var c = rawCollection(key, area, start, end, scenario).map(preprocess(key,scenario));
  var days = ee.List(c.aggregate_array('day')).distinct().sort();
  return ee.ImageCollection.fromImages(days.map(function(day) {
    // Each sensor/UTC date has at most one observation per pixel. Last valid
    // system:index wins in a row overlap; no cross-date spatial gap filling.
    return c.filter(ee.Filter.eq('day',day)).sort('scene_id').mosaic().setDefaultProjection(PROJ)
      .set({day: day, sensor: key, 'system:time_start': ee.Date(day).millis()});
  }));
}
var AUX = ['accepted','qa_denom','range_reject','outside_domain','acq_day'];
function count(col, band, name) {
  return grid(col.merge(ee.ImageCollection.fromImages([empty([band])])).select(band).count().unmask(0,false)).rename(name);
}
function composite(area, start, scenario) {
  start = ee.Date(start); var end = start.advance(2,'month');
  var merged = ee.ImageCollection.fromImages([empty(BANDS.concat(AUX))]);
  var density = [];
  scenario.sensors.forEach(function(key) {
    var c = dailyCollection(key,area,start,end,scenario);
    merged = merged.merge(c);
    density.push(count(c,'NDVI_raw','n_obs_'+key));
  });
  // At most 4 sensors × 62 UTC dates = 248 inputs in a two-month bin.
  // Explicit maxRaw guarantees direct median computation for sensor-day inputs.
  var median = merged.select(BANDS).reduce(ee.Reducer.median(256,1e-6,CONFIG.medianMaxRaw))
    .rename(BANDS).setDefaultProjection(PROJ);
  var counts = ee.Image.cat(density);
  var total = counts.reduce(ee.Reducer.sum()).rename('n_obs');
  return median.addBands(counts).addBands(total)
    .addBands(merged.select('acq_day').mean().rename('mean_observation_day_since_2000'))
    .addBands(count(merged,'qa_denom','n_obs_qa_denom'))
    .addBands(count(merged,'range_reject','n_obs_range_reject'))
    .addBands(count(merged,'outside_domain','n_obs_outside_domain'));
}
// CALENDAR AND AREA SERIES
function periodProperties(start) {
  start=ee.Date(start); var end=start.advance(2,'month');
  // Same display convention as received code: first day of second month.
  var representative=start.advance(1,'month');
  return {period_start:start.format('YYYY-MM-dd'), period_end_exclusive:end.format('YYYY-MM-dd'),
    date_representative:representative.format('YYYY-MM-dd'),
    time_representative_ms:representative.millis(), year:start.get('year'),
    period:start.get('month').subtract(1).divide(2).floor().add(1),
    frequency_per_year:6, version:CONFIG.version, run_id:CONFIG.runId};
}
function periodStarts(firstYear, lastYear) {
  return ee.List.sequence(0,(lastYear-firstYear+1)*6-1).map(function(i) {
    return ee.Date.fromYMD(firstYear,1,1).advance(ee.Number(i).multiply(2),'month').millis();
  });
}
function periodFeature(area, start, scenario) {
  var c=composite(area,start,scenario), mask=validMask(c.select('NDVI_raw'));
  // Bundle the sufficient statistics into one reduction per period, plus one SD.
  var terms=[PIXEL_AREA.rename('area_all'),PIXEL_AREA.multiply(mask).rename('area_valid')];
  ['n_obs','n_obs_qa_denom','n_obs_range_reject','n_obs_outside_domain'].forEach(function(b) {
    terms.push(c.select(b).multiply(PIXEL_AREA).rename(b+'_sum'));
  });
  scenario.sensors.forEach(function(k) {terms.push(c.select('n_obs_'+k).multiply(PIXEL_AREA).rename('n_obs_'+k+'_sum'));});
  BANDS.forEach(function(b){terms.push(c.select(b).multiply(PIXEL_AREA).rename(b+'_sum'));});
  terms.push(c.select('mean_observation_day_since_2000').multiply(PIXEL_AREA).rename('mean_day_sum'));
  var sums=ee.Image.cat(terms).reduceRegion(args(area,ee.Reducer.sum()));
  var total=safeNumber(sums.get('area_all')),valid=safeNumber(sums.get('area_valid')),fraction=valid.divide(total);
  var sd=c.select(BANDS).reduceRegion(args(area,ee.Reducer.stdDev()));
  var props=periodProperties(start);
  props.area=area.name; props.scenario=scenario.id; props.area_asset=area.asset;
  props.valid_area_m2=valid; props.grid_area_m2=total; props.valid_fraction=fraction;
  props.strict_sr=scenario.strictSR; props.family_transfer=scenario.familyTransfer;
  props.observation_unit='sensor_UTC_day';
  props.l7_end_exclusive=scenario.l7End || CONFIG.endExclusive;
  props.n_obs_mean_all=safeNumber(sums.get('n_obs_sum')).divide(total);
  props.n_obs_mean_valid=safeRatio(safeNumber(sums.get('n_obs_sum')),valid);
  props.mean_observation_day_since_2000=safeRatio(safeNumber(sums.get('mean_day_sum')),valid);
  ['n_obs_qa_denom','n_obs_range_reject','n_obs_outside_domain'].forEach(function(b) {
    props[b+'_mean_all']=safeNumber(sums.get(b+'_sum')).divide(total);
  });
  var n=ee.Number(props.n_obs_mean_all), oli=ee.Number(0), scenes=ee.Number(0);
  ['L5','L7','L8','L9'].forEach(function(key) {
    var enabled=scenario.sensors.indexOf(key)!==-1;
    var sceneCount=enabled ? rawCollection(key,area,ee.Date(start),ee.Date(start).advance(2,'month'),scenario).size() : ee.Number(0);
    var obsMean=enabled ? safeNumber(sums.get('n_obs_'+key+'_sum')).divide(total) : ee.Number(0);
    props['n_scenes_'+key]=sceneCount; props['n_obs_'+key+'_mean_all']=obsMean;
    props['valid_observation_share_'+key]=safeRatio(obsMean,n);
    if(key==='L8'||key==='L9') oli=oli.add(obsMean);
    scenes=scenes.add(sceneCount);
  });
  props.n_scenes_total=scenes; props.oli_valid_observation_share=safeRatio(oli,n);
  // These shares describe accepted sensor-day observations, not which observation
  // supplied the median or a weighting of the final median.
  props.no_valid_data=valid.eq(0);
  CONFIG.coverageThresholds.forEach(function(t) {props['eligible_c'+Math.round(t*100)]=valid.gt(0).and(fraction.gte(t));});
  BANDS.forEach(function(b) {
    props[b+'_mean']=safeRatio(safeNumber(sums.get(b+'_sum')),valid);
    props[b+'_spatial_sd']=sd.get(b);
  });
  return ee.Feature(null,props);
}
function table(scenario, firstYear, lastYear) {
  var starts=periodStarts(firstYear,lastYear);
  return ee.FeatureCollection(AREAS.map(function(area) {
    return ee.FeatureCollection(starts.map(function(start) {return periodFeature(area,start,scenario);}));
  })).flatten();
}
// COMMON-SUPPORT COMPARISONS
function pairFeature(area,start,leftScenario,rightScenario,label,leftSensor,rightSensor) {
  var left=composite(area,start,leftScenario), right=composite(area,start,rightScenario);
  var lm=validMask(left.select(BANDS[0])), rm=validMask(right.select(BANDS[0]));
  var common=lm.and(rm),terms=[PIXEL_AREA.rename('area_all'),PIXEL_AREA.multiply(lm).rename('area_left'),
    PIXEL_AREA.multiply(rm).rename('area_right'),PIXEL_AREA.multiply(common).rename('area_common')];
  BANDS.forEach(function(b) {
    var li=left.select(b).updateMask(common),ri=right.select(b).updateMask(common),delta=ri.subtract(li);
    [li.rename(b+'_left_common'),ri.rename(b+'_right_common'),delta.rename(b+'_delta'),
      delta.abs().rename(b+'_abs'),delta.pow(2).rename(b+'_sq'),
      left.select(b).rename(b+'_left_native'),right.select(b).rename(b+'_right_native')].forEach(function(img) {
        terms.push(img.multiply(PIXEL_AREA));
    });
  });
  if(leftSensor && rightSensor) {
    terms.push(left.select('n_obs').updateMask(common).multiply(PIXEL_AREA).rename('left_n_obs_common'));
    terms.push(right.select('n_obs').updateMask(common).multiply(PIXEL_AREA).rename('right_n_obs_common'));
    terms.push(left.select('mean_observation_day_since_2000').updateMask(common).multiply(PIXEL_AREA).rename('left_day_common'));
    terms.push(right.select('mean_observation_day_since_2000').updateMask(common).multiply(PIXEL_AREA).rename('right_day_common'));
  }
  var sums=ee.Image.cat(terms).reduceRegion(args(area,ee.Reducer.sum()));
  var all=safeNumber(sums.get('area_all')),la=safeNumber(sums.get('area_left')),ra=safeNumber(sums.get('area_right'));
  var ca=safeNumber(sums.get('area_common')),props=periodProperties(start);
  props.area=area.name;props.comparison=label;props.left=leftScenario.id;props.right=rightScenario.id;
  props.left_valid_fraction=la.divide(all);props.right_valid_fraction=ra.divide(all);
  props.common_valid_fraction=ca.divide(all); props.common_area_m2=ca;
  props.common_fraction_of_left=safeRatio(ca,la);props.common_fraction_of_right=safeRatio(ca,ra);
  props.n_common_centroid=safeNumber(zero('n').add(1).updateMask(common).reduceRegion(args(area,ee.Reducer.count())).get('n'));
  props.enough_common_centroid=ee.Number(props.n_common_centroid).gte(CONFIG.overlapMinCentroidPixels);
  BANDS.forEach(function(b) {
    props[b+'_left_common']=safeRatio(safeNumber(sums.get(b+'_left_common')),ca);
    props[b+'_right_common']=safeRatio(safeNumber(sums.get(b+'_right_common')),ca);
    props[b+'_delta_mean']=safeRatio(safeNumber(sums.get(b+'_delta')),ca);
    props[b+'_mae']=safeRatio(safeNumber(sums.get(b+'_abs')),ca);
    props[b+'_mse']=safeRatio(safeNumber(sums.get(b+'_sq')),ca);
    props[b+'_rmsd']=ee.Algorithms.If(ca.gt(0),safeNumber(props[b+'_mse']).sqrt(),null);
    props[b+'_delta_native']=ee.Algorithms.If(la.gt(0).and(ra.gt(0)),
      safeNumber(safeRatio(safeNumber(sums.get(b+'_right_native')),ra))
        .subtract(safeNumber(safeRatio(safeNumber(sums.get(b+'_left_native')),la))),null);
  });
  // For L7/L8 these compare L8 raw with transformed L7 and transformed L8 with L7 raw.
  if(leftSensor && rightSensor) {
    props.left_sensor=leftSensor;props.right_sensor=rightSensor;
    props['left_n_obs_common']=safeRatio(safeNumber(sums.get('left_n_obs_common')),ca);
    props['right_n_obs_common']=safeRatio(safeNumber(sums.get('right_n_obs_common')),ca);
    props.mean_observation_day_difference_common=ee.Algorithms.If(ca.gt(0),
      safeNumber(safeRatio(safeNumber(sums.get('right_day_common')),ca))
        .subtract(safeNumber(safeRatio(safeNumber(sums.get('left_day_common')),ca))),null);
    props.epoch=ee.Algorithms.If(ee.Date(start).millis().gte(ee.Date('2017-01-01').millis()),'2017plus','pre2017');
  }
  return ee.Feature(null,props);
}
function comparisons(right,firstYear,lastYear) {
  return ee.FeatureCollection(AREAS.map(function(area) {
    return ee.FeatureCollection(periodStarts(firstYear,lastYear).map(function(s) {
      return pairFeature(area,s,SCENARIOS[0],right,'S1_vs_'+right.id);
    }));
  })).flatten();
}
function overlap(pair,firstYear,lastYear) {
  var left={id:pair[0],sensors:[pair[0]],strictSR:true,familyTransfer:true,l7End:null};
  var right={id:pair[1],sensors:[pair[1]],strictSR:true,familyTransfer:true,l7End:null};
  var starts=periodStarts(firstYear,lastYear);
  if(pair[0]==='L7') starts=starts.filter(ee.Filter.gte('item',ee.Date('2013-05-01').millis()));
  if(pair[1]==='L9') starts=starts.filter(ee.Filter.gte('item',ee.Date('2021-11-01').millis()));
  return ee.FeatureCollection(AREAS.map(function(area) {
    return ee.FeatureCollection(starts.map(function(s) {
      return pairFeature(area,s,left,right,pair.join('_vs_'),pair[0],pair[1]);
    }));
  })).flatten();
}
// SCENE INVENTORY
function sceneInventory(firstYear,lastYear) {
  var out=ee.FeatureCollection([]);
  AREAS.forEach(function(area) {Object.keys(SENSORS).forEach(function(key) {
    var c=rawCollection(key,area,ee.Date.fromYMD(firstYear,1,1),ee.Date.fromYMD(lastYear+1,1,1),SCENARIOS[0]);
    // toList needs a positive limit even when this sensor has no scenes.
    var list=c.toList(c.size().max(1)).map(function(img) {
      img=ee.Image(img);
      return ee.Feature(null,{area:area.name,sensor:key,collection:SENSORS[key].id,
        scene_id:img.get('system:index'), product_id:img.get('LANDSAT_PRODUCT_ID'),
        acquisition_ms:img.get('system:time_start'), day:ee.Date(img.get('system:time_start')).format('YYYY-MM-dd'),
        wrs_path:img.get('WRS_PATH'),wrs_row:img.get('WRS_ROW'),
        sun_elevation:img.get('SUN_ELEVATION'),processing_level:img.get('PROCESSING_LEVEL'),
        sr_algorithm:img.get('ALGORITHM_SOURCE_SURFACE_REFLECTANCE'),
        generated_ms:img.get('DATE_PRODUCT_GENERATED'),cloud_cover_metadata:img.get('CLOUD_COVER'),version:CONFIG.version});
    }); out=out.merge(ee.FeatureCollection(list));
  });});return out;
}
// Clip first: unmasking later must not expose valid data outside the polygon.
function rasterForExport(image, area) {
  return ee.Image(image).select(CONFIG.rasterBand).clip(area.geom)
    .unmask({value: CONFIG.noData, sameFootprint: false}).toFloat();
}
function twoDigits(value) {return ('0' + value).slice(-2);}
var RASTER_SCENARIO = SCENARIOS.filter(function(s) {return s.id === CONFIG.rasterScenario;})[0];
if (!RASTER_SCENARIO || BANDS.indexOf(CONFIG.rasterBand) === -1) {
  throw new Error('Choose S1/S2 and an existing NDVI band.');
}
var PREFIX = 'Cabrobo_' + CONFIG.runId;
// Full-run names match the existing table schema; partial batches get distinct suffixes.
var BATCH_SUFFIX = CONFIG.exportAreas.join(',') === 'A,B' && CONFIG.firstYear === 2000 && CONFIG.lastYear === 2022
  ? '' : '_' + CONFIG.exportAreas.join('_') + '_' + CONFIG.firstYear + '_' + CONFIG.lastYear;
function exportTable(collection, suffix) {
  var name = PREFIX + '_' + suffix + BATCH_SUFFIX;
  Export.table.toDrive({collection: collection, description: name,
    fileNamePrefix: name, folder: CONFIG.driveFolder, fileFormat: 'CSV'});
}
function exportAnalysisTables() {
  SCENARIOS.forEach(function(s) {
    exportTable(table(s, CONFIG.firstYear, CONFIG.lastYear), s.id + '_series');
    if (s.id !== 'S1') exportTable(comparisons(s, CONFIG.firstYear, CONFIG.lastYear), 'S1_vs_' + s.id + '_common');
  });
  [['L5','L7',2000,2011],['L7','L8',2013,2022],['L8','L9',2021,2022]].forEach(function(pair) {
    var first = Math.max(CONFIG.firstYear, pair[2]), last = Math.min(CONFIG.lastYear, pair[3]);
    if (first <= last) exportTable(overlap(pair.slice(0,2), first, last), pair[0] + '_vs_' + pair[1] + '_overlap');
  });
  exportTable(sceneInventory(CONFIG.firstYear, CONFIG.lastYear), 'scene_inventory');
  exportTable(ee.FeatureCollection([ee.Feature(null, {
    version: CONFIG.version, run_id: CONFIG.runId,
    config_json: JSON.stringify(CONFIG), scenarios_json: JSON.stringify(SCENARIOS),
    harmonization_json: JSON.stringify(HARM), sensor_json: JSON.stringify(SENSORS)
  })]), 'configuration');
}
function exportRasterCalendar() {
  var calendar = [];
  AREAS.forEach(function(area) {
    // The same projected rectangle is reused for all periods of this area.
    var region = area.geom.bounds(1, ee.Projection(CONFIG.crs));
    for (var year = CONFIG.firstYear; year <= CONFIG.lastYear; year++) {
      for (var month = 1; month <= 11; month += 2) {
        var start = ee.Date.fromYMD(year, month, 1);
        var basename = PREFIX + '_' + area.name + '_' + CONFIG.rasterScenario + '_' +
          CONFIG.rasterBand + '_' + year + '_' + twoDigits(month) + '_' + twoDigits(month + 1);
        if (CONFIG.exportRasters) Export.image.toDrive({
          image: rasterForExport(composite(area, start, RASTER_SCENARIO), area),
          description: basename, fileNamePrefix: basename, folder: CONFIG.driveFolder,
          region: region, crs: CONFIG.crs, crsTransform: CONFIG.crsTransform,
          maxPixels: CONFIG.maxPixels, fileFormat: 'GeoTIFF', skipEmptyTiles: false,
          formatOptions: {cloudOptimized: true, noData: CONFIG.noData}
        });
        if (CONFIG.exportTables) {
          var properties = periodProperties(start);
          properties.area = area.name; properties.area_asset = area.asset;
          properties.scenario = CONFIG.rasterScenario; properties.ndvi_band = CONFIG.rasterBand;
          properties.period_months = CONFIG.periodMonths;
          properties.period_index = (year - 2000) * 6 + (month + 1) / 2;
          properties.file_name = basename + '.tif'; properties.no_data = CONFIG.noData;
          properties.crs = CONFIG.crs; properties.crs_transform = JSON.stringify(CONFIG.crsTransform);
          properties.observation_unit = 'sensor_UTC_day';
          properties.series_csv = PREFIX + '_' + CONFIG.rasterScenario + '_series' + BATCH_SUFFIX + '.csv';
          properties.source_run_id = CONFIG.sourceRunId;
          calendar.push(ee.Feature(null, properties));
        }
      }
    }
  });
  if (CONFIG.exportTables) exportTable(ee.FeatureCollection(calendar),
    CONFIG.rasterScenario + '_' + CONFIG.rasterBand + '_calendar');
}
var NDVI = {config:CONFIG, scenarios:SCENARIOS, areas:AREAS, table:table,
  comparisons:comparisons, overlap:overlap, sceneInventory:sceneInventory,
  composite:composite, preprocess:preprocess, periodStarts:periodStarts};

// Run creates tasks; submit the requested exports in the Tasks tab.
if (typeof NDVI_HEADLESS === 'undefined' || !NDVI_HEADLESS) {
  print('Export configuration', CONFIG);
  print('Default full run: 276 GeoTIFF tasks and 9 CSV tasks; no tasks are submitted automatically.');
  print('All calendar positions are retained, including periods or pixels with no valid NDVI.');
  Map.centerObject(ROI, 13);
  AREAS.forEach(function(a) {Map.addLayer(a.fc, {}, a.name);});
  if (CONFIG.exportTables) exportAnalysisTables();
  if (CONFIG.exportTables || CONFIG.exportRasters) exportRasterCalendar();
}
