//Google Earth Engine Script to export monthly NDWI-trend statistics point-wise
//(output point attributes: somenumber_statistic_monthnameabbrevation, statistic::one of the following 5 variables:
//NDWI_tau:Mann Kendall's tau, count:count of available months per time series, p_val:p-value, 
//Sslope:Sen's Slope, Sint: Intercept) 

///---------------------user input -----------------------///
// set path variable and choose spatial peatland shape package
//see https://github.com/lgiese/MoistureTrends_EUpeat/blob/main/scripts/MoistureTrendsEUpeatlands_datapreparation_p0.Rmd how to prepare mulitpoint dataset
var aoi_all_path = "projects/gee_path_asset/to/multipoint_csv";

// generate folder with specific name in your google drive for each point shape package
// set folder name to folder in your Google drive
var folder = ee.String('some_folder_name').getInfo()

//click 'Run' and wait until all tasks show up in Tasks (upper right corner)
//press F12 to open browser console and follow: https://github.com/gee-hydro/gee_monkey
///-----------------------------------------------//

//All functions for Mann Kendall Trend Analysis of NDWI based on Landsat satellite data


var lsNDWI = function(aoi){
  //cloud and saturation mask adapted from
  //https://gis.stackexchange.com/questions/425159/how-to-make-a-cloud-free-composite-for-landsat-8-collection-2-surface-reflectanc
  // implemented by Marvin Ludwig 
  
  function maskLcl(image) {
  // Bit 0 - Fill
  // Bit 1 - Dilated Cloud
  // Bit 2 - Unused
  // Bit 3 - Cloud
  // Bit 4 - Cloud Shadow
  // Bit 5 - Snow
  //var snowBit = 1 << 6;
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('111111', 2)).eq(0)
  //  .and(image.select('QA_PIXEL').bitwiseAnd(snowBit).eq(0));
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // Apply the scaling factors to the appropriate bands.
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // Replace the original bands with the scaled ones and apply the masks.
  return image.addBands(opticalBands, null, true)
      .updateMask(qaMask)
      .updateMask(saturationMask)
}



// Define functions for NDWI
var L457NDWI = function(img){
    var ndwi = img.expression('(NIR-SWIR)/(NIR+SWIR)', {
              'NIR': img.select('SR_B4'),
              'SWIR': img.select('SR_B5')
              }).rename('NDWI');
    return img.addBands(ndwi)
  }
  
  
  
var L8NDWI = function(img){

    var ndwi = img.expression('(NIR-SWIR)/(NIR+SWIR)', {
              'NIR': img.select('SR_B5'),
              'SWIR': img.select('SR_B6')
              }).rename('NDWI');
    return img.addBands(ndwi)
  }


var l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2")
  .filterBounds(aoi)
  //.filterDate("1982-09-01", "1984-03-31")
  .filter(ee.Filter.lt('CLOUD_COVER', 80))
  .select("SR_B4", "SR_B5", "QA_PIXEL", "QA_RADSAT")
  .map(maskLcl)
  .map(L457NDWI)
  .select("NDWI")


var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
  .filterBounds(aoi)
  //.filterDate("1984-04-01", "2012-04-30")
  .filter(ee.Filter.lt('CLOUD_COVER', 80))
  .select("SR_B4", "SR_B5", "QA_PIXEL", "QA_RADSAT")
  .map(maskLcl)
  .map(L457NDWI)
  .select("NDWI")
  
var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
  .filterBounds(aoi)
  //.filterDate("2012-05-01", "2013-03-17")
  .filter(ee.Filter.lt('CLOUD_COVER', 80))
  .select("SR_B4", "SR_B5", "QA_PIXEL", "QA_RADSAT")
  .map(maskLcl)
  .map(L457NDWI)
  .select("NDWI")

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(aoi)
    //.filterDate('2019-04-01', '2019-04-30')
    .filter(ee.Filter.lt('CLOUD_COVER', 80))
    .select(['SR_B5', 'SR_B6', 'QA_PIXEL', 'QA_RADSAT'])
    .map(maskLcl)
    .map(L8NDWI)
    .select("NDWI")


var L = l4.merge(l5).merge(l7).merge(l8);

L = L.sort("DATE_ACQUIRED");

  
return(L);


};

//The following code has been produced by Laura Giese
//function for building monthly medians

var YMed = function(monthColl, from, to, month){


// iterate over year list
  
  var year_list = ee.List.sequence(from, to);
  var coll_mon_median = year_list.map(function(year){
  // filter specific year
  var yearly_med = monthColl.select("NDWI").filter(ee.Filter.calendarRange(year, year, 'year')).median();
  
    
    var start = ee.Date.fromYMD(year, month, 1);
    
  
    //reduce collection to one image using median of the rasters
    var img = yearly_med.set({
      'month': month,
      "year": year,
      'system:time_start': start});
    
    var bandnum = img.bandNames().length();
    
    return ee.Algorithms.If(ee.Number(bandnum).eq(0), 0, img);
  
  })

  var coll_mon_median_clean = coll_mon_median.removeAll([0])
  var coll_mm =ee.ImageCollection(coll_mon_median_clean)
return coll_mm
  

}



//calculate Kendalls statistics and counts. 


var funTA = function(coll, aoi){
// by Laura Giese (after https://developers.google.com/earth-engine/tutorials/community/nonparametric-trends)
var afterFilter = ee.Filter.lessThan({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
  primary: coll,
  secondary: coll,
  condition: afterFilter
}));

var sign = function(i, j) { // i and j are images
  return ee.Image(j).neq(i) // Zero case
      .multiply(ee.Image(j).subtract(i).clamp(-1, 1)).int();
};

var kendall = ee.ImageCollection(joined.map(function(current) {
  var afterCollection = ee.ImageCollection.fromImages(current.get('after'));
  return afterCollection.map(function(image) {
    // The unmask is to prevent accumulation of masked pixels that
    // result from the undefined case of when either current or image
    // is masked.  It won't affect the sum, since it's unmasked to zero.
    return ee.Image(sign(current, image)).unmask(0);
  });
  // Set parallelScale to avoid User memory limit exceeded.
}).flatten()).reduce('sum', 2);


var kendall_aoi = kendall.clip(aoi);

//Alternative: compute SensSlope: https://developers.google.com/earth-engine/apidocs/ee-reducer-sensslope
//https://courses.spatialthoughts.com/gee-water-resources-management.html#trend-analysis
//https://rdrr.io/github/USGS-R/smwrStats/man/senSlope.html

//divide by 10000 since we devided by 10000 before and set year image meta information
var processImage = function(image) {
  var image_ndwi = image.divide(10000)
  var year = image.get('year');

  var yearImage = ee.Image.constant(ee.Number(year)).toShort();
  return ee.Image.cat(yearImage, image_ndwi).rename(['year', 'NDWI']).set('year', year);
}

var processedCol = coll.map(processImage)
//print(processedCol, 'procColl')

// Calculate time series slope using sensSlope() reducer function.
var sens = processedCol.reduce(ee.Reducer.sensSlope());

// The resulting image has 2 bands: slope and intercept
// We select the 'slope' and 'offset' (intercept)
var sensSlope = sens.select('slope')
var sensInt = sens.select('offset')

//compute Kendalls' tau
var count_t = coll.count().double()//.clip(aoi);
var count_t_named = count_t.rename("count");
//since definitions how to calculate kendall test statistics and sens solpe exactly slightly vary 
//we use this approach to provide reproducable calculation of tau. taus are exactly the same values as reducer.KendallsCorr
var tau = kendall_aoi.divide(count_t.subtract(1).multiply(count_t).divide(2)).float().rename('tau');


///calculate p-values based on input data - variance
// Values that are in a group (ties).  Set all else to zero.
var groups = coll.map(function(i) {
  var matches = coll.map(function(j) {
    return i.eq(j); // i and j are images.
  }).sum();
  return i.multiply(matches.gt(1));
});

///// Compute tie group sizes in a sequence.  The first group is discarded.
var group = function(array) {
  var length = array.arrayLength(0);
  // Array of indices.  These are 1-indexed.
  var indices = ee.Image([1])
      .arrayRepeat(0, length)
      .arrayAccum(0, ee.Reducer.sum())
      .toArray(1);
  var sorted = array.arraySort();
  var left = sorted.arraySlice(0, 1);
  var right = sorted.arraySlice(0, 0, -1);
  // Indices of the end of runs.
  var mask = left.neq(right)
  // Always keep the last index, the end of the sequence.
      .arrayCat(ee.Image(ee.Array([[1]])), 0);
  var runIndices = indices.arrayMask(mask);
  // Subtract the indices to get run lengths.
  var groupSizes = runIndices.arraySlice(0, 1)
      .subtract(runIndices.arraySlice(0, 0, -1));
  return groupSizes;
};

// See equation 2.6 in Sen (1968).
var factors = function(image) {
  return image.expression('b() * (b() - 1) * (b() * 2 + 5)');
};

var groupSizes = group(groups.toArray());
var groupFactors = factors(groupSizes);
var groupFactorSum = groupFactors.arrayReduce('sum', [0])
      .arrayGet([0, 0]);

var count = joined.count();

var kendallVariance = factors(count)
    .subtract(groupFactorSum)
    .divide(18)
    .float();


///// Compute Z-statistics.
var zero = kendall.multiply(kendall.eq(0));
var pos = kendall.multiply(kendall.gt(0)).subtract(1);
var neg = kendall.multiply(kendall.lt(0)).add(1);

var z = zero
    .add(pos.divide(kendallVariance.sqrt()))
    .add(neg.divide(kendallVariance.sqrt()));
//var z_aoi = z.clip(aoi);
//Map.addLayer(z_aoi, {min: -2, max: 2}, 'z');

// https://en.wikipedia.org/wiki/Error_function#Cumulative_distribution_function
function eeCdf(z) {
  return ee.Image(0.5)
      .multiply(ee.Image(1).add(ee.Image(z).divide(ee.Image(2).sqrt()).erf()));
}

function invCdf(p) {
  return ee.Image(2).sqrt()
      .multiply(ee.Image(p).multiply(2).subtract(1).erfInv());
}

// Compute P-values.
var p = ee.Image(1).subtract(eeCdf(z.abs())).float();
var p_aoi = p.rename('p_val'); //.clip(aoi).
//Map.addLayer(p_aoi, {min: 0, max: 0.5}, 'p');

var monthMK_add = tau.addBands(count_t_named).addBands(p_aoi).addBands(sensSlope).addBands(sensInt).float(); //


//print(monthMK_add)
return monthMK_add.float()
}



// wrapper function for trend analysis and handling na's 

var fImg = function(aoi, from_year, to_year){
// by Laura Giese


var L_start = lsNDWI(aoi)
//multiply by 10000
var L = L_start.map(function(img){
  var img_10000 = img.multiply(10000).copyProperties(img, ["system:time_start"])
  return(img_10000)
});
//create list of month to iterte over
var month_list = ee.List.sequence(1,12)
var month_list_i = ee.List.sequence(-1,11)
//print(month_list.get(0), 'month_list_el')
var month_name_list = ee.List(['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'])
  var from = from_year
  var to = to_year
//var month=2

var monthcoll_ndwi = month_list.map(function(month){
  //select data for one month
  var monthCollection = L.select("NDWI").filter(ee.Filter.calendarRange(month,month,"month"));
  var monthColl = ee.ImageCollection(monthCollection);
  //print(monthColl, 'collection_band')
  
  //check if enough bands exist for calculating median
  var monthcoll_img = monthColl.toBands();
  //print('monthcoll_img',monthcoll_img)
  var bandnum_base = monthcoll_img.bandNames().length();
  var list_excl = ee.List([0, 1])

  
  //call function to calculate monthly medians
  var yearMK = ee.Algorithms.If(list_excl.contains(ee.Number(bandnum_base)), 0, YMed(monthColl, from, to, month));
  var yearColl = ee.ImageCollection(yearMK);
  //print('yearColl', yearColl)

  //check if enough bands exist to continue
  var yearcoll_img_if = ee.Algorithms.If(list_excl.contains(ee.Number(bandnum_base)), 0, yearColl.toBands());
  var yearcoll_img = ee.Image(yearcoll_img_if)
  var bandnum = ee.Algorithms.If(list_excl.contains(ee.Number(bandnum_base)), 0, yearcoll_img.bandNames().length());
  //var check_cond = ee.Algorithms.If(list_excl.contains(bandnum), bandnum, 'multiband trend');
  //print('check_cond',check_cond)
  
  //get index for monthname
  var month_index = month_list_i.get(month);
  var mon_nam = ee.String(month_name_list.get(month_index));
  //built name strings
  var count_str = ee.String('count_').cat(mon_nam);
  var ndwi_str = ee.String('NDWI_tau_').cat(mon_nam);
  var p_str = ee.String('p_val_').cat(mon_nam); 
  var slope_str = ee.String('Sslope_').cat(mon_nam);
  var ic_str = ee.String('Sint_').cat(mon_nam); 
  var name_list = ee.List([ndwi_str,count_str, p_str, slope_str, ic_str]); //
  
  // calculate p-value, kendalls tau, count variable (after: https://developers.google.com/earth-engine/tutorials/community/nonparametric-trends)
  var monthMK_add = ee.Algorithms.If(list_excl.contains(ee.Number(bandnum)), 0,funTA(yearColl, aoi));

  var monthMK_add_image = ee.Image(monthMK_add); //.float();
  var monthMK_named = ee.Algorithms.If(list_excl.contains(ee.Number(bandnum)), 0, monthMK_add_image.rename(name_list));
  //print(monthMK_named,'monthMK_named')

  return ee.Algorithms.If(list_excl.contains(ee.Number(bandnum)), 0, monthMK_named);
  //var monthcoll_ndwi = ee.Algorithms.If(ee.Number(bandnum).eq(0), 0, monthMK_named);
});
//print('monthcount', monthcoll_ndwi)


//print(monthcount, "Monthcount")
var monthcount_cl = monthcoll_ndwi.removeAll([0]);
var monthcount = ee.ImageCollection(monthcount_cl).toBands();//.clipToCollection(aoi);
return monthcount
  
}

///wrapper function to export results in loop (using apply)  


var funExport = function(aoi_all, number, ID_numbers, from_year, to_year, folder, aoi_all_path){
// by Laura Giese

  //iteration variable (multipoint feature number)      
  var iter = number

//call feature by system index
var aoi_num = ID_numbers.get(iter); 
  //print(aoi_num, 'aoi_num')
var aoi = aoi_all.filter(ee.Filter.eq("system:index", aoi_num)); 
  //print(aoi,'aoi')


//multipoint

var multiPoints = aoi.geometry();
//print(multiPoints);

// Apply the coordinates method to the MultiPoint object.
var multiPointCoordinates = multiPoints.coordinates();


var iter_len_m = multiPointCoordinates.length()
var i_len_m = ee.Number(iter_len_m).subtract(ee.Number(1)).getInfo()
//shorten list for tests
//var i_len = ee.Number(3).getInfo()
//print('i_len', i_len)

var iter_list_multi = ee.List.sequence(0, i_len_m)

// add empty band for all band names which are not included in image with trend statistics to keep NA information in table and
// to prevent loosing columns
var trend_iter = iter_list_multi.map(function(iter_multi){


var multiPointCoordinates5 = multiPointCoordinates.get(iter_multi)
var aoi_geom = ee.Geometry.Point(multiPointCoordinates5)
var aoi_buff = aoi_geom.buffer(30)

//read all bandnames already included in trend image and assign unified bandnames (later: column names)
var coll_trend_img = fImg(aoi_geom, from_year, to_year, aoi_num)
//bandname
var bN = coll_trend_img.bandNames()
var new_bN = bN.map(function(bn){
  var bn_string = ee.String(bn)
  var split_bn = bn_string.split('_');
  var slit_len=split_bn.length()
  var bn_new_str = ee.Algorithms.If(slit_len.eq(3), ee.String(split_bn.get(1)).cat('_').cat(split_bn.get(2)), ee.String(split_bn.get(1)).cat('_').cat(split_bn.get(2)).cat('_').cat(split_bn.get(3)))
  return(bn_new_str)
})
var coll_trend_img_renamed = coll_trend_img.rename(new_bN).unmask()
var exp_names_all = ee.List(['NDWI_tau_apr','NDWI_tau_aug','NDWI_tau_feb','NDWI_tau_jul','NDWI_tau_jun',
'NDWI_tau_mar','NDWI_tau_may','NDWI_tau_nov','NDWI_tau_oct','NDWI_tau_sep','Sint_apr','Sint_aug',
'Sint_feb','Sint_jul','Sint_jun','Sint_mar','Sint_may','Sint_nov','Sint_oct','Sint_sep','Sslope_apr',
'Sslope_aug','Sslope_feb','Sslope_jul','Sslope_jun','Sslope_mar','Sslope_may','Sslope_nov','Sslope_oct',
'Sslope_sep','count_apr','count_aug','count_feb','count_jul','count_jun','count_mar','count_may','count_nov',
'count_oct','count_sep','p_val_apr','p_val_aug','p_val_feb','p_val_jul','p_val_jun','p_val_mar','p_val_may',
'p_val_nov','p_val_oct','p_val_sep','NDWI_tau_jan','Sint_jan','Sslope_jan','count_jan','p_val_jan','NDWI_tau_dec',
'Sint_dec','Sslope_dec','count_dec','p_val_dec'])
var missing_bands = exp_names_all.map(function(l){
  return ee.Algorithms.If(new_bN.contains(ee.String(l)), 0, l)
});

var missing_bands_clean = missing_bands.removeAll([0])
var missing_bands_len = missing_bands_clean.length()
var const_bands_string = ee.String('0x00')
var empty_band_vals = ee.List.repeat(const_bands_string, missing_bands_len)
var empty_img = ee.Image();
var empty_multiB_img = missing_bands_clean.map(function(addB){
  var eBn_List = ee.List(addB)
  var em_img = ee.Image(0).selfMask().toFloat().unmask()//.rename(eBn_List)
  return em_img
})
var empty_img_named = ee.ImageCollection.fromImages(empty_multiB_img).toBands().rename(missing_bands_clean)

var coll_trend_img_renamed_allBands = coll_trend_img_renamed.addBands(empty_img_named)

//print(coll_trend_img_renamed)
var newft = coll_trend_img_renamed_allBands.sample({
                                  region: aoi_geom, 
                                  //collection: aoi_geom, 
                                  //properties: ['FID'], 
                                  geometries: true,
                                  projection: 'EPSG:4326',
                                  scale: 30,
                                  dropNulls: false
                                  })

var newft_single = ee.Feature(newft.first())                                
return(newft_single)
})      


var trend_iter_coll=ee.FeatureCollection(trend_iter)


var aoi_all_path_str = ee.String(aoi_all_path)
var pkg_n_split = aoi_all_path_str.split('/');

  var split_len=pkg_n_split.length()
var split_num = ee.Number(split_len).subtract(ee.Number(1)).getInfo()

var id_str = ee.String(pkg_n_split.get(split_num))
var iter_str= iter.toString()
var gap = ee.String('_')
var mp_str = ee.String('mp')
print(id_str)

var filename = ee.String(id_str).cat(mp_str).cat(gap).cat(iter_str).getInfo()

//export ndwi trend statistics for all month as csv table
Export.table.toDrive({
  collection: trend_iter_coll,
  description: filename,
  folder: folder,
  fileFormat: 'CSV'
});


};

//----------------------------------------------------
// set start & end year for trend analysis
var from_year = 1984
var to_year = 2022

//Create Feature Collection of All Input multipoint features
var aoi_all = ee.FeatureCollection(aoi_all_path)

//convert to point list to apply function on each shape
var aoi_L = aoi_all.toList(60000)
print(aoi_L)

//get peatland index in export list
var ID_numbers= aoi_L.map(function(aoi){
  var feat_aoi = ee.Feature(aoi)
  var sys_id = feat_aoi.id()
  return sys_id
}); 

//get length of export list
var iter_len = ID_numbers.length();
var i_len = ee.Number(iter_len).getInfo();//.subtract(ee.Number(1)).getInfo()

print('i_len', i_len);

//run export task for all multipoints
//var i = 0
var tasks_exported = Array.apply(null, {length: i_len}).map(Number.call, Number) 
            .map(function(number){
              return funExport(aoi_all, number, ID_numbers, from_year, to_year, folder, aoi_all_path)});

