#### data preparation for environmental covariates script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in Feb 2025

#### completely revised in March 2025 to include new raster data, make script run as background process and improve readability


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
LANDCOVER_2021 <- rast(x = paste0(path, '/Gola_lc_all_rf50_classification_200122_1341.tif'))
LANDCOVER_2023 <- rast(x = paste0(path, '/classification_2023.tif')) # land use map from 2023 which was created from Brit, Felicity thinks, that it overpredicts cocoa
DTM <- rast(x = paste0(path, '/Gola90mDTM.tif'))
SLOPE <- rast(x = paste0(path, '/ASTER30SLOPE.tif'))
VEG_2013 <- rast(x = paste0(path, '/MOD13Q1.A2013001.h16v08.061.2021224063452.hdf'))
VEG_2020 <- rast(x = paste0(path, '/MOD13Q1.A2020353.h16v08.061.2021012024041.hdf'))
JRC_TRANSITION <- terra::rast(raster::raster(x = paste0(path, '/JRC_TMF_Transition2021_GGL.tif'))) # there were issues with terra::rast, thus using raster::raster for all JRC tif data  
JRC_DEGRAD <- terra::rast(raster::raster(x = paste0(path, '/JRC_TMF_DegradationYear2021_GGL.tif'))) # more information on the JRC data at https://forobs.jrc.ec.europa.eu/TMF/resources/tutorial/gee
JRC_ANN_CHANGE <- terra::rast(x = paste0(path, '/JRC_TMF_AnnualChanges2021_GGL.tif'))

# immediately throw away the stuff which isn't needed
JRC_ANN_CHANGE_2013 <- JRC_ANN_CHANGE[[c('Dec2013')]] # here select the year that is needed, but then also change it in the if else statement, I think data is only available until Dez2021
JRC_ANN_CHANGE_2021 <- JRC_ANN_CHANGE[[c('Dec2021')]] # here select the year that is needed, but then also change it in the if else statement, I think data is only available until Dez2021
rm(JRC_ANN_CHANGE)
VEG_2013 <-  VEG_2013[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]
VEG_2020 <-  VEG_2020[[c('"250m 16 days NDVI"', '"250m 16 days EVI"')]]

# vector data 
rivers <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_rivers.shp')
study_area <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_PH_study_area.shp')
settlements <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Settlements_2021_Points.shp')
reserves <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Greater_Gola_Landscape_2024_polygons.shp')
# roads <- st_read()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Set WGS 84 / UTM zone 29N CRS for all spatial objects ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set crs to WGS 84 / UTM zone 29N, 32629 for all data
spat_rasters <- ls()[sapply(mget(ls()), inherits, "SpatRaster")] # get all spat rasters 
target_crs <- crs('EPSG:32629') # set target crs 
for(i in 1:length(spat_rasters)){ # loop over all spatRasters and project to new crs 
    raster <- get(spat_rasters[i])  # get raster from 
    if(crs(raster) != target_crs){
      if(names(raster)[1] == "TransitionClass" | names(raster)[1] == "DegradationYear" | names(raster)[1] == "Dec2021"| names(raster)[1] == "Dec2013"){ # change Dez2021 here if needed
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

# save all projected spatial data as Rdata to not have to do all this again and again
save(DTM, JRC_ANN_CHANGE_2013, JRC_ANN_CHANGE_2021, JRC_DEGRAD, JRC_TRANSITION, LANDCOVER_2021, LANDCOVER_2023, 
     LANDCOVER_2021, rivers, settlements, SLOPE, study_area, VEG_2013, VEG_2020, file = 'data/PH_projected_spatial_data.Rdata')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create a grid and buffered study area##########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create grid, for 2km2 choose 1519 as cellsize, 500 is rather small, nearly breakes my laptop
grid_sf <- st_make_grid(study_area, cellsize = 1519.671, square = F) %>% # creates a st object, for squares: square = T, F is polygons - here make sure that we create the appropriate grid cell size!
  st_as_sf() %>% # store as sf object 
  st_filter(study_area, .predicate = st_intersects) %>% # filter all polygons that lie within or touch the study area
  mutate(CellID = row_number())

# create points in the center of each grid cell
grid_point_sf <- st_make_grid(study_area, cellsize = 1519.671, square = F, what = 'centers') %>% # make sure that cellsize is the same as above!
  st_as_sf() %>% # store as sf object 
  st_filter(grid_sf, .predicate = st_intersects) %>% # filter all polygons that lie within or touch the study area
  mutate(CellID = row_number())

# create buffered study area 
study_area_buffer <- study_area %>% st_buffer(2000)  # buffers by 2000m 

# plot 
plot(st_geometry(grid_sf), main = 'Study area')
plot(st_geometry(study_area), add = T)
plot(st_geometry(grid_point_sf), add = T)

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

# plot land cover data 
plot(LANDCOVER_2021)
plot(LANDCOVER_2023)
plot(study_area, add = T)

# improve LANDCOVER_2023 data structure
# LANDCOVER_2023 <- as.factor(LANDCOVER_2023)
# levels(LANDCOVER_2023) <- data.frame(Value=1:7, label = c('bareground', 'regrowth1', 'farm_bush', 'forest', 'water', 'swamp', 'cocoa')) # watch out, these are not the actual labels but spaceholder
# cats(LANDCOVER_2023)

# extract data from LANDCOVER_2021
# landcover_21_stats <- exactextractr::exact_extract(x = LANDCOVER_2021, grid_sf, weights = 'area', fun = 'weighted_frac')
# names(landcover_21_stats) <- paste0(str_replace(cats(LANDCOVER_2021)[[1]]$label, ' ', '_'), '_2021')
# landcover_21_stats <- landcover_21_stats %>% mutate(CellID = 1:nrow(.))

# extract data from LANDCOVER_2023
# landcover_23_stats <- exactextractr::exact_extract(x = LANDCOVER_2023, grid_sf, weights = 'area', fun = 'weighted_frac')
# names(landcover_23_stats) <- paste0(str_replace(cats(LANDCOVER_2023)[[1]]$label, ' ', '_'), '_2023')
# landcover_23_stats <- landcover_23_stats %>% mutate(CellID = 1:nrow(.))

# extract data from MODIS for vegatation in 2013
veg_13 <- exact_extract(VEG_2013, grid_sf, fun = 'weighted_mean', weights = 'area', append_cols = 'CellID') %>% 
  rename(NDVI_2013 = `weighted_mean."250m 16 days NDVI"`, 
         EVI_2013 = `weighted_mean."250m 16 days EVI"`)

# extract data from MODIS for vegatation in 2020
veg_20 <- exact_extract(VEG_2020, grid_sf, fun = 'weighted_mean', weights = 'area', append_cols = 'CellID') %>% 
  rename(NDVI_2020 = `weighted_mean."250m 16 days NDVI"`, 
         EVI_2020 = `weighted_mean."250m 16 days EVI"`)

# extract data from SLOPE and DTM
elev_mean <- exactextractr::exact_extract(DTM, grid_sf, fun = 'weighted_mean', weights = 'area')
slope_mean <- exactextractr::exact_extract(SLOPE, grid_sf, fun = 'weighted_mean', weights = 'area')

# more information on JRC data here: https://forobs.jrc.ec.europa.eu/static/tmf/TMF_DataUsersGuide.pdf

# improve data structure of JRC_TRANSITION
JRC_TRANSITION <- as.factor(JRC_TRANSITION)
class_names <- data.frame(Value = levels(JRC_TRANSITION)[[1]]$ID, 
                          label = c('Undisturbed_tropical_moist_forest', ## class 10
                                    'Degraded_forest_short_duration_disturbance_before_2014', # class 21
                                    'Degraded_forest_short_duration_disturbance_after_2014', # class 22
                                    'Degraded_forest_long_duration_disturbance_before_2014', # class 23
                                    'Degraded_forest_long_duration_disturbance_after_2014', # class 24
                                    'Degraded_forest_2_3_degradation_periods_before_2014', # class 25
                                    'Degraded_forest_2_3_degradation_periods_after_2014', # class 26
                                    'Regrowth_desturbed_before_2004', 'Regrowth_desturbed_between_2004_2013', 'Regrowth_desturbed_between_2014_2020', # class 31, 32, 33, respectively
                                    'Deforestation_started_before_2013','Deforestation_started_between_2013_2020', 
                                    'Deforestation_started_2021', 'Deforestation_started_2022', 'Deforestation_started_2023', 'Degradation_started_2023',  
                                    'Permanent_water','Seasonal_water', 'Deforestation_to_permanent_water', 'Deforestation_to_seasonal_water', 
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

# extract data from JRC_DEGRAD, omitted for now since I daut this is needed!
# unique(JRC_DEGRAD) # this raster has only one value per cell , possibly it gives the year of degradation for the first time?
# hist(as.data.frame(JRC_DEGRAD)$DegradationYear)
# plot(ifel(JRC_DEGRAD != 0, 1, 0))

# improve data structure from JRC_ANN_CHANGE - extracted data for Dec2021 which then has values from 1 to 6, legend in https://forobs.jrc.ec.europa.eu/static/tmf/TMF_DataUsersGuide.pdf
JRC_ANN_CHANGE_2013 <- as.factor(JRC_ANN_CHANGE_2013)
JRC_ANN_CHANGE_2021 <- as.factor(JRC_ANN_CHANGE_2021)
change_names <- data.frame(ID = 1:nrow(levels(JRC_ANN_CHANGE_2013)[[1]]), 
                           label = c('Undisturbed_tropical_moist_forest', 
                                     'Degraded_tropical_moist_forest', 
                                     'Deforested_land', 
                                     'Tropical_moist_forest_regrowth', 
                                     'Permanent_and_seasonal_water', 
                                     'Other_land_cover'))
levels(JRC_ANN_CHANGE_2013)[[1]] <- change_names %>% mutate(label = paste(change_names$label, 'Dec2013', sep = '_'))
levels(JRC_ANN_CHANGE_2021)[[1]] <- change_names %>% mutate(label = paste(change_names$label, 'Dec2021', sep = '_'))
# extract data from JRC_ANN_CHANGE, which are basically rather broad land cover classes!
plot(JRC_ANN_CHANGE_2021); plot(JRC_ANN_CHANGE_2013)
JRC_ann_changes_2013_stats <- exactextractr::exact_extract(JRC_ANN_CHANGE_2013, grid_sf, fun = 'weighted_frac', weights = 'area', append_cols = 'CellID')
JRC_ann_changes_2021_stats <- exactextractr::exact_extract(JRC_ANN_CHANGE_2021, grid_sf, fun = 'weighted_frac', weights = 'area', append_cols = 'CellID')
names(JRC_ann_changes_2013_stats) <- c('CellID', paste('JRC_ann_changes', change_names$label, 'Dec2013', sep = '_'))
names(JRC_ann_changes_2021_stats) <- c('CellID', paste('JRC_ann_changes', change_names$label, 'Dec2021', sep = '_'))


# join data to grid_sf
grid_sf <- grid_sf %>% 
  # left_join(landcover_21_stats, join_by(CellID)) %>% 
  # left_join(landcover_23_stats, join_by(CellID)) %>% 
  left_join(veg_13, join_by(CellID)) %>% 
  left_join(veg_20, join_by(CellID)) %>% 
  left_join(JRC_ann_changes_2013_stats, join_by(CellID)) %>% 
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



# distance to settlement, to large river
dist_settle <- st_distance(grid_sf, settlements, which = 'Euclidean') %>% # gives a distance matrix
  apply(1, min) 
dist_large_river <- st_distance(grid_sf, rivers %>% filter(River_Size == 'Large'), which = 'Euclidean') %>% # calc dist to large, main river streams
  apply(1, min) 
grid_sf <- grid_sf %>% mutate(Distance_settlement = set_units(dist_settle, m), 
                              Distance_large_river = set_units(dist_large_river, m)) 



###############################################################################
##### CONTINUE HERE TO CREATE A LEAKAGE BELT ##################################
###############################################################################

# location in forest reserve, national park, community forest or outside 
reserves <- reserves %>% rename(Name = NAME, Reserve_Type = DESIG) %>% 
  filter(Name %in% c('Gola North', 'Gola South', 'Gola Central', 'Bunumbu', 'Tiwai Island Sanctuary', 'Kambui South', 'Kambui Hills and Extensions')) %>% 
  select(Name, Reserve_Type) # clean shp file up 
reserves_buffered <- vect(st_buffer(reserves %>% filter(Name %in% c('Gola North', 'Gola South', 'Gola Central')), dist = 4000))

str(reserves_buffered)
reserves_buffered <- terra::union(reserves_buffered)
reserves_buffered

plot(reserves_buffered)
  st_union()  %>% st_cast('POLYGON')
  
  st_difference(reserves) %>% st_cast('POLYGON')
test[[1]]
st_geometry_type(test)
plot(st_geometry(test))
plot(st_geometry(test[[7]]))

location <- st_intersection(reserves, grid_sf) %>% select(CellID, Reserve_Type) %>% 
  mutate(area = set_units(st_area(.), km^2))%>% 
  group_by(CellID, Reserve_Type) %>% 
  summarise(sum_area = sum(area)) %>% 
  group_by(CellID) %>%
  slice_max(sum_area, n = 1) # this removes the reserve types which are not dominant / have the most area in a grid cell

grid_sf <- grid_sf %>% left_join(location %>% st_drop_geometry() %>% select(Reserve_Type), join_by(CellID)) %>% 
  mutate(if_else(is.na(Reserve_Type), 'Outside', Reserve_Type))
  
library(tmap)
tm_shape(reserves) +
  tm_polygons()

################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Save full set of environmental covariates as shp #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# export all environmental covariates as csv and grid as shp file
st_write(grid_sf %>% select(CellID, x), "data/PH_grid.shp", append = F)
write.csv(grid_sf %>% st_drop_geometry(), "data/PH_prepared_env_covariates.csv")













#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### x. Calculate distance from every camera to nearest river ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in all data sets and remove unneeded data - these are only example site locations 
#deploy_cam <- read_excel(path = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/combined_PH_data2025_draft7.xlsx', sheet = 'locations_cameras') %>% 
#  select(ProjectIDName, PlotIDName, CameraIDName, UTM_X_meters, UTM_Y_meters) %>% 
#  distinct() %>%
#  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), crs = 32629) # transform into sf object


# calculate euclidean distance - choose euclidean or circle dist, ATTENTION ON CRS - please doublecheck again!
#dist_river<- st_distance(deploy_cam, rivers, which = 'Euclidean') %>% # gives a distance matrix
#  apply(1, min) # find the min dist from each point [point,] to each river [,river]

# save distance to nearest point in deploy_cam df 
#deploy_cam <- deploy_cam %>% mutate(Dist_river_m = dist_river)











# make an extensive plot of results 
#forest <- ggplot(grid_sf) +
#  geom_sf(aes(fill = forest), color = "black") +
#  scale_fill_viridis_c(option = "plasma") +
#  theme_minimal() + 
#  labs(title = 'Percentage of Forest Cover per Grid Cell in Gola, SL') 
#swamp <- ggplot(grid_sf) +
#  geom_sf(aes(fill = swamp), color = "black") +
#  scale_fill_viridis_c(option = "plasma") +
#  theme_minimal() + 
#  labs(title = 'Percentage of Swamp Cover per Grid Cell in Gola, SL')
#water <- ggplot(grid_sf) +
#  geom_sf(aes(fill = water), color = "black") +
#  scale_fill_viridis_c(option = "plasma") +
#  theme_minimal() + 
#  labs(title = 'Percentage of Water Cover per Grid Cell in Gola, SL')
#river <- ggplot(grid_sf) +
#  geom_sf(aes(fill = as.numeric(river_density)), color = "black") +
#  scale_fill_viridis_c(option = "plasma", name = 'River Density') +
#  theme_minimal() + 
#  labs(title = 'River Density in km/km2 in Gola, SL') 
#large_river <- ggplot(grid_sf) +
#  geom_sf(aes(fill = as.numeric(Distance_large_river)), color = "black") +
#  scale_fill_viridis_c(option = "plasma", name = 'Distance to Large River') +
#  theme_minimal() + 
#  labs(title = 'Distance to Large Rivers in Gola, SL') 
#camps <- ggplot(grid_sf) +
#  geom_sf(aes(fill = as.numeric(Distance_campsite)), color = "black") +
#  scale_fill_viridis_c(option = "plasma", name = 'Distance to Campsites') +
#  theme_minimal() + 
#  labs(title = 'Distance to Campsites in Gola, SL')
#settlements <- ggplot(grid_sf) +
#  geom_sf(aes(fill = as.numeric(Distance_settlement)), color = "black") +
#  scale_fill_viridis_c(option = "plasma", name = 'Distance to Settlements') +
#  theme_minimal() + 
#  labs(title = 'Distance to Settlements in Gola, SL')
#elevation <- ggplot(grid_sf) +
#  geom_sf(aes(fill = as.numeric(mean_elev_m)), color = "black") +
#  scale_fill_viridis_c(option = "plasma", name = 'Mean Elevation') +
#  theme_minimal() + 
#  labs(title = 'Mean Elevation per Cell in Gola, SL')
#library(ggpubr)
#ggarrange(forest, swamp, water, river)
#ggarrange(large_river, camps, settlements)

#plot(grid_sf %>% select(-CellID, -area, -river_length, -river_length_med_large, -bareground, -farm_bush, -regrowth1, -Distance_campsite, -oil_palm))


