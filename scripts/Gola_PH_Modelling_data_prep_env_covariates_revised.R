#### data preparation for environmental covariates script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in Feb 2025

#### completely revised in May 2025 to make the script shorter and improve readability


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(units)
library(terra)
library(raster)
library(stringr)
library(exactextractr) # for very fast implementation of exact extraction of values from raster 
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in data and tidy everything up #### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set path
path <- 'C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/raster data'

# raster data - capital letters
DTM <- rast(x = paste0(path, '/Gola90mDTM.tif'))
SLOPE <- rast(x = paste0(path, '/ASTER30SLOPE.tif'))
VEG_aug_2013 <- rast(x = paste0(path, '/MOD13Q1.A2013225.h16v08.061.2021237170725.hdf')) # from 13.08.2013 0:00 to 28.08.2013 23:59
VEG_feb_2013 <- rast(x = paste0(path, '/MOD13Q1.A2013033.h16v08.061.2021226153239.hdf')) # from 02.02.2013 0:00 to 17.02.2013 23:59
VEG_aug_2017 <- rast(x = paste0(path, '/MOD13Q1.A2017225.h16v08.061.2021280152407.hdf')) # from 13.08.2017 0:00 to 28.08.2017 23:59
VEG_feb_2017 <- rast(x = paste0(path, '/MOD13Q1.A2017033.h16v08.061.2021266145520.hdf')) # from 02.02.2017 0:00 to 17.08.2017 23:59
VEG_feb_2021 <- rast(x = paste0(path, '/MOD13Q1.A2021033.h16v08.061.2021049223006.hdf')) # from 02.02.2021 0:00 to 17.02.2021 23:59
VEG_aug_2021 <- rast(x = paste0(path, '/MOD13Q1.A2021225.h16v08.061.2021243200737.hdf')) # from 13.08.2021 0:00 to 28.08.2021 23:59
JRC_ANN_CHANGE <- terra::rast(x = paste0(path, '/JRC_TMF_AnnualChanges2021_GGL.tif'))
JRC_TRANSITION <- terra::rast(raster::raster(x = paste0(path, '/JRC_TMF_Transition2021_GGL.tif'))) 
# JRC_DEGRAD <- terra::rast(raster::raster(x = paste0(path, '/JRC_TMF_DegradationYear2021_GGL.tif'))) # more information on the JRC data at https://forobs.jrc.ec.europa.eu/TMF/resources/tutorial/gee


# immediately throw away the stuff which isn't needed
JRC_ANN_CHANGE_2013 <- JRC_ANN_CHANGE[[c('Dec2013')]] # here select the year that is needed, but then also change it in the if else statement, I think data is only available until Dez2021
JRC_ANN_CHANGE_2021 <- JRC_ANN_CHANGE[[c('Dec2021')]] # here select the year that is needed, but then also change it in the if else statement, I think data is only available until Dez2021
JRC_ANN_CHANGE_2017 <- JRC_ANN_CHANGE[[c('Dec2017')]] # as intermediate year 
rm(JRC_ANN_CHANGE)
VEG_aug_2013 <-  VEG_aug_2013[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_feb_2013 <-  VEG_feb_2013[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_aug_2017 <-  VEG_aug_2017[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_feb_2017 <-  VEG_feb_2017[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_aug_2021 <-  VEG_aug_2021[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_feb_2021 <-  VEG_feb_2021[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]

# vector data 
rivers <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_rivers.shp')
study_area <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_PH_study_area.shp')
# settlements <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Settlements_2021_Points.shp')
reserves <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Greater_Gola_Landscape_2024_polygons.shp')
leakage_belt <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_Leakage_Belt.shp') # this file does not really help due to its mismatch with the reserves layer:/
roads <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Roads_InProgressGola.shp')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Set WGS 84 / UTM zone 29N CRS for all spatial objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set crs to WGS 84 / UTM zone 29N, 32629 for all data
spat_rasters <- ls()[sapply(mget(ls()), inherits, "SpatRaster")] # get all spat rasters 
target_crs <- crs('EPSG:32629') # set target crs 
for(i in 1:length(spat_rasters)){ # loop over all spatRasters and project to new crs 
    raster <- get(spat_rasters[i])  # get raster from 
    if(crs(raster) != target_crs){
      if(names(raster)[1] == "TransitionClass" | names(raster)[1] == "DegradationYear" | names(raster)[1] == "Dec2021"| names(raster)[1] == "Dec2017"| names(raster)[1] == "Dec2013"){ # change Dez2021 here if needed
        raster <- terra::project(raster, target_crs, , method = 'near')  # project to new crs using terra which is rather slow and appropriate handling for classes
        assign(spat_rasters[i], raster, envir = .GlobalEnv)  # save to global environment again
        message(spat_rasters[i], " has been projected to target crs.")
      } else{
        raster <- terra::project(raster, target_crs)  # project to new crs using terra 
        assign(spat_rasters[i], raster, envir = .GlobalEnv)  # save to global environment again
        message(spat_rasters[i], " has been projected to target crs.")
      }
    } else{
      message(spat_rasters[i], " already in target crs, thus, not projected to save time.")
    }
}

# set crs to WGS 84 / UTM zone 29N, 32629 for all data
spat_vectors <- ls()[sapply(mget(ls()), inherits, "sf")]
target_crs <- 32629 # set target crs 
for(i in 1:length(spat_vectors)){ # loop over all spatVectors and project to new crs 
  vector <- get(spat_vectors[i])  # get vectors 
  if(st_crs(vector)$epsg != target_crs){
    vector <- st_transform(vector, target_crs)  # project to new crs using terra which is rather slow
    assign(spat_vectors[i], vector, envir = .GlobalEnv)  # save to global environment again
    rm(vector)
  } else {
   message(spat_vectors[i], " already in ", target_crs, ", thus, not projected to save time.")
  }
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create a grid and buffered study area##########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create grid, for ~2km2 choose 1519 as cellsize
grid_sf <- st_make_grid(study_area, cellsize = 1519.671, square = F) %>% # creates a st object, for squares: square = T, F is polygons - here make sure that we create the appropriate grid cell size!
  st_as_sf() %>% # store as sf object 
  st_filter(study_area, .predicate = st_intersects) %>% # filter all polygons that lie within or touch the study area
  mutate(CellID = row_number())

# create buffered study area 
study_area_buffer <- study_area %>% st_buffer(2000)  # buffers by 2000m 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Crop all rasters by buffered study area ##########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

spat_rasters # all spatialRasters from global environment
for(i in 1:length(spat_rasters)){
  raster <- get(spat_rasters[i])
  raster <- terra::crop(raster, study_area_buffer)
  assign(spat_rasters[i], raster, envir = .GlobalEnv)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Calculate fraction/mean per grid cell of raster data classes from different sources ##########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract area-weighted means of vegetation indices for all related layers at once

veg_layers <- ls()[grep('VEG', ls())] # get all layers from global environment with VEG 

for(veg in veg_layers){
  raster <- get(veg) # get different rasters
  parts <- strsplit(veg, "_")[[1]] # get the different months and years for saving the data later again
  year <- parts[2]
  month <- parts[3]
  
  extracted <- exact_extract(raster, grid_sf, fun = 'weighted_mean', weights = 'area', append_cols = 'CellID') %>% 
    rename(!!paste0('NDVI_', year, "_", month) := `weighted_mean."250m 16 days NDVI"`, 
           !!paste0('EVI_', year, "_", month) :=  `weighted_mean."250m 16 days EVI"`)
  assign(paste0(veg, '_extracted'), extracted, envir = .GlobalEnv) # save them back to global environment 
}


# extract elevation data from DTM and slope from SLOPE
elev_mean <- exactextractr::exact_extract(DTM, grid_sf, fun = 'weighted_mean', weights = 'area')
slope_mean <- exactextractr::exact_extract(SLOPE, grid_sf, fun = 'weighted_mean', weights = 'area')

# extract data from JRC annual changes layers
# more information on JRC data here: https://forobs.jrc.ec.europa.eu/static/tmf/TMF_DataUsersGuide.pdf

# prepare data extraction from JRC annual changes dataset
ann_changes <- ls()[grep('JRC_ANN_CHANGE', ls())]
change_names <- data.frame(ID = 1:6, 
                           label = c('Undisturbed_tropical_moist_forest', 
                                     'Degraded_tropical_moist_forest', 
                                     'Deforested_land', 
                                     'Tropical_moist_forest_regrowth', 
                                     'Permanent_and_seasonal_water', 
                                     'Other_land_cover'))

# call loop to extract data from JRC annual data raster layers
for(ann in ann_changes){
  raster <- get(ann)
  year <- strsplit(ann, "_")[[1]][4] # get the years the data is available for 
  
  raster <- as.factor(raster) # save as raster to avoid issues with projection
  
  extracted <- exactextractr::exact_extract(raster, grid_sf, fun = 'weighted_frac', weights = 'area', append_cols = 'CellID')
  names(extracted) <- c('CellID', paste('JRC_ann_changes', change_names$label, year, sep = '_')) 
  
  assign(paste0('JRC_ann_changes_', year, '_stats'), extracted, envir = .GlobalEnv) # save extracted data back to global environment 
}


# get data from JRC Transition
JRC_TRANSITION <- as.factor(JRC_TRANSITION) # improve data structure of JRC_TRANSITION
class_names <- data.frame(Value = levels(JRC_TRANSITION)[[1]]$ID, 
                          label = c('Undisturbed_tropical_moist_forest', ## class 10
                                    'Degraded_forest_short_duration_disturbance_before_2014', # class 21
                                    'Degraded_forest_short_duration_disturbance_after_2014', # class 22
                                    'Degraded_forest_long_duration_disturbance_before_2014', # class 23
                                    'Degraded_forest_long_duration_disturbance_after_2014', # class 24
                                    'Degraded_forest_2_3_degradation_periods_before_2014', # class 25
                                    'Degraded_forest_2_3_degradation_periods_after_2014', # class 26
                                    'Regrowth_desturbed_before_2004', # class 31
                                    'Regrowth_desturbed_between_2004_2013', # class 32
                                    'Regrowth_desturbed_between_2014_2020', # class 33
                                    'Deforestation_started_before_2013', # class 41
                                    'Deforestation_started_between_2013_2020', # class 42
                                    'Deforestation_started_2021', 'Deforestation_started_2022', 'Deforestation_started_2023', 'Degradation_started_2023', # class 51 to 54  
                                    'Permanent_water','Seasonal_water', 'Deforestation_to_permanent_water', 'Deforestation_to_seasonal_water', # class 71 to 74
                                    'Old_plantation', 'Plantation_regrowing_disturbed_before_2014', 'Plantation_regrowing_disturbed_between_2014_2020', 'Conversion_to_plantation_before_2014', 'Conversion_to_plantation_between_2015_2021', 'Recent_conversion_to_plantation_started_2021',
                                    'Other_lc_without_afforestation','Afforestation_young', 'Afforestation_old', 'Water_recently_converted_into_forest_regrowth')) 
levels(JRC_TRANSITION) <- class_names

# extract data from JRC transition map -sub types
JRC_transition_stats <- exactextractr::exact_extract(JRC_TRANSITION, grid_sf, fun = 'weighted_frac', weights = 'area', append_cols = 'CellID')
names(JRC_transition_stats) <- c('CellID', paste0('JRC_transition_', class_names$label)) # set appropriate class names, taken from here: https://forobs.jrc.ec.europa.eu/static/tmf/TMF_DataUsersGuide.pdf#[{%22num%22%3A16%2C%22gen%22%3A0}%2C{%22name%22%3A%22XYZ%22}%2C69%2C736%2C0]

# select only those variables which might be of interest
JRC_transition_stats <- JRC_transition_stats %>% select(CellID, JRC_transition_Undisturbed_tropical_moist_forest, JRC_transition_Degraded_forest_short_duration_disturbance_before_2014, 
                                                        JRC_transition_Degraded_forest_short_duration_disturbance_after_2014, JRC_transition_Degraded_forest_long_duration_disturbance_before_2014, 
                                                        JRC_transition_Degraded_forest_long_duration_disturbance_after_2014, JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014, 
                                                        JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014, JRC_transition_Regrowth_desturbed_between_2004_2013, JRC_transition_Regrowth_desturbed_between_2014_2020)




# join all data to grid_sf
grid_sf <- grid_sf %>% 
  left_join(VEG_feb_2013_extracted, join_by(CellID)) %>% 
  left_join(VEG_aug_2013_extracted, join_by(CellID)) %>% 
  left_join(VEG_feb_2017_extracted, join_by(CellID)) %>% 
  left_join(VEG_aug_2017_extracted, join_by(CellID)) %>% 
  left_join(VEG_feb_2021_extracted, join_by(CellID)) %>% 
  left_join(VEG_aug_2021_extracted, join_by(CellID)) %>%
  left_join(JRC_ann_changes_2013_stats, join_by(CellID)) %>% 
  left_join(JRC_ann_changes_2017_stats, join_by(CellID)) %>% 
  left_join(JRC_ann_changes_2021_stats, join_by(CellID)) %>% 
  left_join(JRC_transition_stats, join_by(CellID)) %>% 
  mutate(mean_elev = elev_mean, 
         mean_slope = slope_mean)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Calculate metrics using vector layer ##########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# river density of medium and large rivers per grid cell 
rivers_cell <- st_intersection(rivers, grid_sf) # Intersect rivers with grid to get river segments within each grid cell, drops all the cells where no rivers are

# define river density as length of rivers within each grid cell 
river_length_med_large <- rivers_cell %>% 
  filter(River_Size %in% c('Large', 'Medium')) %>% # filter out small rivers 
  mutate(length_m = st_length(geometry, which = 'Euclidean')) %>%  # compute length of each river segment per grid cell
  group_by(CellID) %>%  # group by grid cell
  summarise(river_length_med_large = sum(length_m, na.rm = TRUE)) %>%
  st_drop_geometry()

# transfer river length to grid_sf and calculate density
grid_sf <- grid_sf %>% left_join(river_length_med_large, join_by(CellID)) %>% 
  mutate(area = set_units(st_area(x), km^2),
         river_density_med_large = river_length_med_large/area) %>% 
  replace(is.na(.), 0) %>% select(-river_length_med_large) # NA results from grid cells without rivers



# distance to large river, and roads
dist_large_river <- st_distance(grid_sf, rivers %>% filter(River_Size == 'Large'), which = 'Euclidean') %>% # calc dist to large, main river streams
  apply(1, min) 
dist_road <- st_distance(grid_sf, roads %>% filter(Category == 'Road'), which = 'Euclidean') %>% # calc dist to roads
  apply(1, min) 


grid_sf <- grid_sf %>% mutate(Distance_large_river = as.numeric(set_units(dist_large_river, m)), 
                              Distance_road = as.numeric(set_units(dist_road, m))) 


# location in forest reserve, national park, community forest or outside 
reserves <- reserves %>% rename(Name = NAME, Reserve_Type = DESIG) %>% 
  filter(Name %in% c('Gola North', 'Gola South', 'Gola Central', 'Tiwai Island Sanctuary', 'Kambui South', 'Kambui Hills and Extensions')) %>% 
  select(Name, Reserve_Type) %>% 
  st_zm() %>% # remove z dimension
  st_cast('POLYGON')# clean shp file up 

# tidy up the leakage belt shp and make it fit the reserves sf
leakage_belt <- leakage_belt %>% select(geometry) %>% mutate(Name = 'Gola', Reserve_Type = 'Leakage Belt')

# bring leakage belt and 
reserves_lk_belt <- bind_rows(reserves, leakage_belt)

# prepare area of reserve types per grid cell 
location <- st_intersection(reserves_lk_belt %>% st_transform(32629), grid_sf) %>% select(CellID, Reserve_Type) %>% 
  mutate(area = set_units(st_area(.), km^2))%>% 
  group_by(CellID, Reserve_Type) %>% 
  summarise(area = sum(area)) %>% 
  group_by(CellID) %>%
  slice_max(area, n = 1)#  %>% mutate(CellID = as.numeric(CellID))

# join back to grid and assign Reserve Type for outside if the area which is outside any reserve or leakage belt is larger than inside
loc <- grid_sf %>% mutate(area = set_units(st_area(.), km^2)) %>%
  select(area, CellID) %>%
  left_join(location %>% select(Reserve_Type, CellID, area_l = area) %>% st_drop_geometry(), join_by(CellID)) %>% 
  mutate(area = as.numeric(area), area_l = as.numeric(area_l),
         area_out = area- area_l, 
         area_out = as.numeric(if_else(is.na(area_out), area, area_out)), 
         Reserve_Type = if_else(is.na(Reserve_Type), 'Outside', Reserve_Type), 
         Reserve_Type = if_else(as.numeric(area_out) > (area/2), 'Outside', Reserve_Type)) 

# join data back to grid
grid_sf <- grid_sf %>% left_join(loc %>% select(CellID, Reserve_Type) %>% st_drop_geometry(), join_by(CellID))

# check that everything went okay 
tmap_mode(mode = 'view')
tm_shape(reserves_lk_belt) +
  tm_polygons(fill='Reserve_Type') +
tm_shape(grid_sf) +
  tm_polygons(fill = 'Reserve_Type', fill_alpha = 0.4)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Save full set of environmental covariates as shp #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# export all environmental covariates as csv and grid as shp file
st_write(grid_sf %>% select(CellID, x), "data/PH_grid.shp", append = F)
write.csv(grid_sf %>% st_drop_geometry(), "data/PH_prepared_env_covariates.csv")



