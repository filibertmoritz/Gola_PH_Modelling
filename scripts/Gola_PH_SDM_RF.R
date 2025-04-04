#### pygmy hippo random forest sdm to explore important variables for integrated occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025 
#### Thanks for hints and comments from Felicity Edwards and Steffen Oppel

#### Background: We have data from various PH surveys which can be mainly distinguished in camera trap surveys, transect surveys and opportunistic presence-only data 
#### Issue: The data is rather scattered and we are unsure, if its enough to fit a more sophisticated, integrated model 
#### Solution: pool all data and explore, a) if we can fit a more basic random forest model (to afterwards, possibly go on to the integrated model) and b) explore which predictors could be relevant for the integrated model 
#### Target: Built a RF model with all pooled data 
#### Steps: a) load elsewhere prepared data sets and fit them into one df, b) load elsewhere prepared environmental covariates, c) create/design an effort variable to use as an offset in the RF model
####        d) fit RF model, predict/plot distribution of PH and e) explore the effect of the different variables 


##### THINGS TO IMPROVE 
# read through the random forest model literature and Steffen's scripts to check that I do not make any silly mistakes 
# add plot that explore the effect of all the variables 

##### BACKGROUND READING 
# https://christophm.github.io/interpretable-ml-book/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(sf)
library(terra)
library(tidyverse)
library(lubridate)
library(hms)
library(data.table)
library(ranger)
library(scales)
library(tmap)
library(units)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

getwd()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. load data which has been prepared elsewhere #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# presence and deploy/survey data 
pres_opp <- fread(file = 'data/PH_prepared_pres_opp_data.csv', stringsAsFactors = T) %>% select(-V1)
pres_cam <- fread(file = 'data/PH_prepared_pres_cam_data.csv', stringsAsFactors = T) %>% select(-V1)
pres_transects <- fread(file = 'data/PH_prepared_pres_transect_data.csv', stringsAsFactors = T) %>% select(-V1)
deploy_cam <- fread(file = 'data/PH_prepared_deploy_cam_data.csv', stringsAsFactors = T) %>% select(-V1)
locs_transects_paths <- st_read('data/PH_prepared_transect_paths.shp')
locs_transects_data <- fread(file = 'data/PH_prepared_transect_paths.csv', stringsAsFactors = T) %>% select(-V1)
locs_transects <- locs_transects_paths %>% left_join(locs_transects_data, join_by(uniqueID)) # combine shp and csv file to one object 
rm(locs_transects_data, locs_transects_paths) # remove unneeded objects 

# grid and environmental covariates 
grid_sf <- st_read('data/PH_grid.shp') %>% st_as_sf() # laod shp and transform to sf 
envCovs <- fread('data/PH_prepared_env_covariates.csv') %>% 
  select(-V1) %>% as_tibble()
envCovs_sf <- grid_sf %>% left_join(envCovs, join_by(CellID)) %>% st_as_sf()# transfer geometry from grid to envCovs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Merge all presence data into one big data frame, beforehand filter for independence #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare opportunistic presences as sf
head(pres_opp)
str(pres_opp)
pres_opp <- pres_opp %>% 
  rename(UTM_X_meters = UTM_X_m, UTM_Y_meters = UTM_Y_m) %>% 
  mutate(Obs_Method = 'Opportunistic') %>%
  select(Project, Country, Obs_DateTime, UTM_X_meters, UTM_Y_meters, Sign, Obs_Method) 

# prepare transect presences as sf
pres_transects <- pres_transects %>% 
  mutate(Country = 'SierraLeone', 
         Obs_Method = 'Transect Survey')

# prepare camera trap presences as sf
head(pres_cam)
pres_cam <- pres_cam %>% mutate(Obs_Method = 'Camera Trap') %>% 
  select(Project, Country, Obs_DateTime, UTM_X_meters , UTM_Y_meters, SiteID, Obs_Method)

# merge data together 
pres_sf <- bind_rows(pres_cam, pres_opp,
                     pres_transects) %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Create effort variable #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save camera deployment data as sf 
deploy_cam_sf <- deploy_cam  %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)
rm(deploy_cam)


# effort variable for camera trap data 
deploy_cam_sf <- deploy_cam_sf %>% mutate(Deployment_Time = difftime(Collection, Deployment, units = 'days')) # calculate deployment time 

effort_cam <- deploy_cam_sf %>% # sum up deployment time and number of cameras per cell
  st_join(grid_sf) %>% group_by(CellID) %>% 
  summarise(CameraEffort_Time = sum(Deployment_Time), 
            CameraEffort_N = n(), 
            .groups = 'keep') %>% 
  st_drop_geometry()

# effort variable for transect data 
effort_transects <- locs_transects %>% 
  st_intersection(grid_sf) %>%  # drops all the cells where no rivers are
  mutate(TransectEffort_Length = st_length(geometry, which = 'Euclidean')) %>%  # compute length of each river segment per grid cell
  group_by(CellID) %>%  # group by grid cell
  summarise(TransectEffort_Length = sum(TransectEffort_Length, na.rm = TRUE)) %>%
  st_drop_geometry()

# bring effort data together 
effort_sf <- grid_sf %>% 
  left_join(effort_cam, join_by(CellID)) %>% 
  left_join(effort_transects, join_by(CellID))
effort_sf[is.na(effort_sf)] <- 0 # replace all NAs with 0 

# calculate one effort index with rescaled variables to between 0 and 1
effort_sf <- effort_sf %>%
 mutate(CameraEffort = (rescale(as.numeric(CameraEffort_Time), to = c(0,1)) + 
                          rescale(as.numeric(CameraEffort_N), to = c(0,1)) / 2), # one Camera Trap effort variable as mean from both N and Time
         TransectEffort = rescale(as.numeric(TransectEffort_Length), to = c(0,1)), # rescaled Transect length per cell
         Effort = (CameraEffort + (TransectEffort))/2) # mean of both efforts


# calculate an effort variable for presence_only data 

# distance weighted influence 
# weighted euclidean distance

### FILL SOMETHING SMART IN - possibly the effort from other survey methods as benchmark where people went?

# plot 
tmap_mode("view")  # interactive mode
tm_shape(effort_sf) + 
  tm_polygons(fill = "Effort", lwd = .2, 
              fill.scale = tm_scale_continuous(values = "viridis", midpoint = 0)) +
tm_shape(locs_transects) + 
  tm_lines(col = 'red', lwd = 1.5) +
tm_shape(deploy_cam_sf) +
  tm_dots(fill = 'orange', size = .3) +
tm_shape(pres_sf) +
  tm_symbols(fill = 'brown', shape = 3)

# make a plot which is easier to export 
# plot raw observations
effort_plot<-ggplot(effort_sf) +
  geom_sf(aes(fill = Effort), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Survey Effort in Gola, Sierra Leone",
       fill = "Effort", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Create presence-absence information per grid cell #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

occu <- pres_sf %>% select(geometry) %>%
  st_join(effort_sf) %>%  # effort data transferred as well to possibly infer on 0's if there was effort, but no presence in a grid cell
  st_drop_geometry() %>% 
  group_by(CellID) %>% 
  mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
  distinct() 
occu_sf <- grid_sf %>% left_join(occu, join_by(CellID)) # transfer to full grid 

occu_sf[is.na(occu_sf)] <- 0 # this replaces all NAs - places where no species where seen and no effort occured - to 0's

plot(occu_sf)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Bring all data together #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data overview
effort_sf %>% class()
envCovs_sf %>% class()
occu_sf %>% class()

# create one massive data frame 
data <- envCovs_sf %>% 
  left_join(effort_sf %>% st_drop_geometry(), join_by(CellID)) %>% 
  left_join(occu_sf %>% select(CellID, Occu) %>% st_drop_geometry(), join_by(CellID))
plot(data)

# create a training data set where either a presence was observed or surveys were performed!
traindata <- data %>% filter(Occu == 1 | Effort > 0) %>% st_drop_geometry()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Look at correlation of features/predictors #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data overview
library(GGally)
data %>% 
  st_drop_geometry() %>%
  select(river_density, river_density_med_large, Distance_campsite, Distance_settlement, 
         Distance_large_river, bareground, regrowth1, farm_bush, forest, water, swamp, 
         cocoa, oil_palm, mean_elev, NDVI_2013, NDVI_2020, EVI_2013, EVI_2020, Occu) %>% 
  ggpairs() 


# correlation
library(corrplot)
data %>% 
  st_drop_geometry() %>% 
  select(# bareground_2021, regrowth1_2021, farm_bush_2021,forest_2021, water_2021, swamp_2021, cocoa_2021, oil_palm_2021, # land cover data 2021 from ?
        # bareground_2023, regrowth1_2023, farm_bush_2023, regrowth1_2023, forest_2023, water_2023, swamp_2023, cocoa_2023, # land cover data 2023 from Britt, not sure if classes are assigned correctly!
        NDVI_2013, NDVI_2020, EVI_2013, EVI_2020, # vegetation indices for two times from MODIS
        mean_elev, mean_slope, 
        JRC_ann_changes_2021Undisturbed_tropical_moist_forest, JRC_ann_changes_2021Degraded_tropical_moist_forest, JRC_ann_changes_2021Deforested_land, JRC_ann_changes_2021Tropical_moist_forest_regrowth, JRC_ann_changes_2021Undisturbed_tropical_moist_forest, 
        JRC_transition_Undisturbed_tropical_moist_forest, JRC_transition_Degraded_forest_short_duration_disturbance_before_2014, 
        JRC_transition_Degraded_forest_short_duration_disturbance_after_2014, JRC_transition_Degraded_forest_long_duration_disturbance_before_2014, 
        JRC_transition_Degraded_forest_long_duration_disturbance_after_2014, JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014, 
        JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014, JRC_transition_Regrowth_desturbed_between_2004_2013, JRC_transition_Regrowth_desturbed_between_2014_2020, 
        Distance_settlement, Distance_large_river, river_density_med_large, Occu) %>%  
  cor() %>% 
  corrplot(tl.pos = "l")

# pca
pc <- prcomp(data %>% 
               st_drop_geometry() %>%
               select(NDVI_2013, NDVI_2020, EVI_2013, EVI_2020, # vegetation indices for two times from MODIS
                      mean_elev, mean_slope, 
                      JRC_ann_changes_2021Undisturbed_tropical_moist_forest, JRC_ann_changes_2021Degraded_tropical_moist_forest, JRC_ann_changes_2021Deforested_land, JRC_ann_changes_2021Tropical_moist_forest_regrowth, JRC_ann_changes_2021Undisturbed_tropical_moist_forest, 
                      JRC_transition_Undisturbed_tropical_moist_forest, JRC_transition_Degraded_forest_short_duration_disturbance_before_2014, 
                      JRC_transition_Degraded_forest_short_duration_disturbance_after_2014, JRC_transition_Degraded_forest_long_duration_disturbance_before_2014, 
                      JRC_transition_Degraded_forest_long_duration_disturbance_after_2014, JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014, 
                      JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014, JRC_transition_Regrowth_desturbed_between_2004_2013, JRC_transition_Regrowth_desturbed_between_2014_2020, 
                      Distance_settlement, Distance_large_river, river_density_med_large, Occu),
             center = TRUE,
             scale. = TRUE)

# biplot 
library('ggfortify')
autoplot(pc, data = data %>% 
           st_drop_geometry() %>%
           select(NDVI_2013, NDVI_2020, EVI_2013, EVI_2020, # vegetation indices for two times from MODIS
                  mean_elev, mean_slope, 
                  JRC_ann_changes_2021Undisturbed_tropical_moist_forest, JRC_ann_changes_2021Degraded_tropical_moist_forest, JRC_ann_changes_2021Deforested_land, JRC_ann_changes_2021Tropical_moist_forest_regrowth, JRC_ann_changes_2021Undisturbed_tropical_moist_forest, 
                  JRC_transition_Undisturbed_tropical_moist_forest, JRC_transition_Degraded_forest_short_duration_disturbance_before_2014, 
                  JRC_transition_Degraded_forest_short_duration_disturbance_after_2014, JRC_transition_Degraded_forest_long_duration_disturbance_before_2014, 
                  JRC_transition_Degraded_forest_long_duration_disturbance_after_2014, JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014, 
                  JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014, JRC_transition_Regrowth_desturbed_between_2004_2013, JRC_transition_Regrowth_desturbed_between_2014_2020, 
                  Distance_settlement, Distance_large_river, river_density_med_large, Occu) %>% mutate(Occu = as.factor(Occu)), 
         loadings = TRUE, loadings.label = TRUE, colour = 'Occu') + 
  theme_bw()

# another biplot
library(factoextra)
fviz_pca_biplot(pc, label = "var", habillage = data$Occu, addEllipses = TRUE)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_data_exploration_pca.jpg', plot = last_plot())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Train random forest model #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Steffen recommended two options to include effort as an 'offset' into the random forest model: 
#   1)  Incorporate effort in the response: divide the number of hippo sightings by the effort (in days or so), 
#       and model this ratio. The resulting output predictions could then be multiplied by the max effort (days) 
#       to get the ‘true’ occupancy. This is crude, but should work hopefully.
#   2)  Use effort as a normal predictor variable, with the response being just 0/1 or the number of hippos. This may 
#       however result in garbage output, because you would then essentially be modelling the occurrence of a sighting and 
#       NOT the occurrence of an animal (which may be high even when there is no effort and no sighting).
#       where you artificially set ‘effort’ to the maximum recorded anywhere, which would in theory give you the probability 
#       of occurrence if there had been exhaustive survey effort. However, this gets complicated very quickly, because if the 
#       effort variable is very strong and dominant and explains most of the variation, then all the other forest structure 
#       variables will basically be ignored by the random forest (or have only very mild influence) which would then result 
#       in a map that predicts that hippos are everywhere when there is maximum effort everywhere.

# However, here I will also try all this WITHOUT ANY EFFORT variables as option 3)

# Things to consider: 
#     a)  Craig and Huettmann 2009 used binary logistic regression in TreeNet and balanced option under class weights, which automatically reweighted each class to 
#         account unequal sample size between presence and absence points
#     b)  validate predictions by a confusion matrix which report precentages of presences and absences correctly classified (cutoff at 0.5)
#     c)  Use effort as weights? 
#     d)  Valavi et al. 2023: class overlap (presence/absence is not clearly distinguishable by one variable - common for changes in occupancy across ecological gradients and rare species)
#         - solutions in paper, also example code


# 1) Incorporate effort in the response

# check for distribution of data 
traindata %>% select(Occu, Effort) %>% summary() # there are lots of 0 in both Effort and Occu!

# create a separate response variable which a) incorporates the Effort but b) doesn't become Inf from dividing by 0
traindata <- traindata %>% 
  mutate(Effort_manipulated = if_else(Occu == 1 & Effort == 0, 0.0001, Effort), # if effort is 0 but an opportunistic presence occured, this ensures, that no Inf occures! ITS A BIT CRUDE
         CameraEffort_Time = as.numeric(CameraEffort_Time), 
         TransectEffort_Length = as.numeric(TransectEffort_Length),
         Occu_Effort = Occu*Effort_manipulated) # divide Occu by Effort

rf1 <- ranger(Occu_Effort ~ NDVI_2013+NDVI_2020+EVI_2013+ EVI_2020+ # vegetation indices for two times from MODIS
              mean_elev+mean_slope+ # topography
              JRC_ann_changes_2021Undisturbed_tropical_moist_forest+JRC_ann_changes_2021Degraded_tropical_moist_forest+JRC_ann_changes_2021Deforested_land+JRC_ann_changes_2021Tropical_moist_forest_regrowth+ JRC_ann_changes_2021Undisturbed_tropical_moist_forest+ 
              JRC_transition_Undisturbed_tropical_moist_forest+ JRC_transition_Degraded_forest_short_duration_disturbance_before_2014+ 
              JRC_transition_Degraded_forest_short_duration_disturbance_after_2014+JRC_transition_Degraded_forest_long_duration_disturbance_before_2014+ 
              JRC_transition_Degraded_forest_long_duration_disturbance_after_2014+JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014+ 
              JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014+JRC_transition_Regrowth_desturbed_between_2004_2013+ JRC_transition_Regrowth_desturbed_between_2014_2020+ 
              Distance_settlement+ Distance_large_river+river_density_med_large, # rivers
              data = traindata, num.trees = 1000, 
              #mtry = T, 
              oob.error = T,
              importance = 'permutation', # choose one of: "impurity" (Based on Gini impurity reduction), "impurity_corrected" (bias corrected), "permutation" (Based on permuting feature values and measuring performance drop)
              sample.fraction = 0.5) # perform down-sampling, as suggested in Valavi et al. 2022
summary(rf1)
rf1


# this is the shallow-tuned model from Valavi et al 2022
#rf2 <- ranger(Occu ~ NDVI_2013+NDVI_2020+EVI_2013+ EVI_2020+ # vegetation indices for two times from MODIS
#                mean_elev+mean_slope+ # topography
#                JRC_ann_changes_2021Undisturbed_tropical_moist_forest+JRC_ann_changes_2021Degraded_tropical_moist_forest+JRC_ann_changes_2021Deforested_land+JRC_ann_changes_2021Tropical_moist_forest_regrowth+ JRC_ann_changes_2021Undisturbed_tropical_moist_forest+ 
#                JRC_transition_Undisturbed_tropical_moist_forest+ JRC_transition_Degraded_forest_short_duration_disturbance_before_2014+ 
#                JRC_transition_Degraded_forest_short_duration_disturbance_after_2014+JRC_transition_Degraded_forest_long_duration_disturbance_before_2014+ 
#                JRC_transition_Degraded_forest_long_duration_disturbance_after_2014+JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014+ 
#                JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014+JRC_transition_Regrowth_desturbed_between_2004_2013+ JRC_transition_Regrowth_desturbed_between_2014_2020+ 
#                Distance_settlement+ Distance_large_river+river_density_med_large + # rivers
#                Effort,
#              data = traindata,
#              num.trees = 2000, 
#              splitrule = 'hellinger',
#              max.depth = 4, 
#              probability = T,
#              replace = T,
#              oob.error = T,
#              importance = 'permutation', # choose one of: "impurity" (Based on Gini impurity reduction), "impurity_corrected" (bias corrected), "permutation" (Based on permuting feature values and measuring performance drop)
#              #sample.fraction = 0.5 # perform down-sampling, as suggested in Valavi et al. 2022
#              ) 
#rf2

# 2) use effort as response variable 
rf2 <- ranger(Occu ~ NDVI_2013+NDVI_2020+EVI_2013+ EVI_2020+ # vegetation indices for two times from MODIS
                mean_elev+mean_slope+ # topography
                JRC_ann_changes_2021Undisturbed_tropical_moist_forest+JRC_ann_changes_2021Degraded_tropical_moist_forest+JRC_ann_changes_2021Deforested_land+JRC_ann_changes_2021Tropical_moist_forest_regrowth+ JRC_ann_changes_2021Undisturbed_tropical_moist_forest+ 
                JRC_transition_Undisturbed_tropical_moist_forest+ JRC_transition_Degraded_forest_short_duration_disturbance_before_2014+ 
                JRC_transition_Degraded_forest_short_duration_disturbance_after_2014+JRC_transition_Degraded_forest_long_duration_disturbance_before_2014+ 
                JRC_transition_Degraded_forest_long_duration_disturbance_after_2014+JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014+ 
                JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014+JRC_transition_Regrowth_desturbed_between_2004_2013+ JRC_transition_Regrowth_desturbed_between_2014_2020+ 
                Distance_settlement+ Distance_large_river+river_density_med_large+ # rivers
                Effort, 
              data = traindata, num.trees = 1000, 
              #mtry = T, 
              importance = 'permutation',
              sample.fraction = 0.5) # perform down-sampling, as suggested in Valavi et al. 2022
rf2

# 3) model occu without considering effort 
rf3 <- ranger(Occu ~ NDVI_2013+NDVI_2020+EVI_2013+ EVI_2020+ # vegetation indices for two times from MODIS
                mean_elev+mean_slope+ # topography
                JRC_ann_changes_2021Undisturbed_tropical_moist_forest+JRC_ann_changes_2021Degraded_tropical_moist_forest+JRC_ann_changes_2021Deforested_land+JRC_ann_changes_2021Tropical_moist_forest_regrowth+ JRC_ann_changes_2021Undisturbed_tropical_moist_forest+ 
                JRC_transition_Undisturbed_tropical_moist_forest+ JRC_transition_Degraded_forest_short_duration_disturbance_before_2014+ 
                JRC_transition_Degraded_forest_short_duration_disturbance_after_2014+JRC_transition_Degraded_forest_long_duration_disturbance_before_2014+ 
                JRC_transition_Degraded_forest_long_duration_disturbance_after_2014+JRC_transition_Degraded_forest_2_3_degradation_periods_before_2014+ 
                JRC_transition_Degraded_forest_2_3_degradation_periods_after_2014+JRC_transition_Regrowth_desturbed_between_2004_2013+ JRC_transition_Regrowth_desturbed_between_2014_2020+ 
                Distance_settlement+ Distance_large_river+river_density_med_large, 
              data = traindata, num.trees = 1000, 
              #mtry = T, 
              importance = 'permutation',
              sample.fraction = 0.5) # perform down-sampling, as suggested in Valavi et al. 2022
rf3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Get Impression of Importance of the environmental covariates #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reading on partial dependence plots
# https://www.rpubs.com/vishal1310/QuickIntroductiontoPartialDependencePlots

# alternatively, use vivid package with nice tutorials: https://alaninglis.github.io/vivid/articles/rangerVivid.html
# however, a few computations are enormously computation intensive!

# get variable importance plots
library(vip)
vip(rf1, all_permutations = T, include_type = T) + labs(title = 'Variable Importance in RF1') + theme_bw()
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_vip_rf1.jpg', plot = last_plot())
vip(rf2, all_permutations = T, include_type = T) + labs(title = 'Variable Importance in RF2') + theme_bw()
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_vip_rf2.jpg', plot = last_plot())
vip(rf3, all_permutations = T, include_type = T) + labs(title = 'Variable Importance in RF3') + theme_bw()
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_vip_rf3.jpg', plot = last_plot())

# get most imortant variables 
variables <- names(sort(importance(rf1), decreasing = TRUE))

# plot partial predictive plots, first using vivid
library(vivid)
gc() # its rather computation intensive
viv_plot1 <- vivid::pdpPairs(data = traindata, 
                fit = rf1, 
                response = 'Occu_Effort', 
                nmax = 100, # samples rows from data to display, NULL means all rows/observations, 50 is good choice for computation time 
                gridSize = 13, # huge influence on computation time, set to 20 for higher resolution in pixel-plot
                nIce = 100, # Number of ice curves to be plotted, defaults to 30.
                vars = variables[1:10]) # choose important variables here 
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_feature_influence_rf1_new.jpg', plot = viv_plot1, width = 20, height = 16)

gc()
variables2 <- names(sort(importance(rf2), decreasing = TRUE))
viv_plot2 <- vivid::pdpPairs(data = traindata, 
                fit = rf2, 
                response = 'Occu', 
                nmax = 100, # samples rows from data to display, NULL means all rows/observations, 50 is good choice for computation time 
                gridSize = 10, # huge influence on computation time, set to 20 for higher resolution in pixel-plot
                nIce = 100, # Number of ice curves to be plotted, defaults to 30.
                vars = variables2[1:10]) # choose important variables here 
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_feature_influence_rf2_new.jpg', plot = viv_plot2, width = 20, height = 16)

gc()
variables3 <- names(sort(importance(rf3), decreasing = TRUE))
viv_plot3 <- vivid::pdpPairs(data = traindata, 
                fit = rf3, 
                response = 'Occu', 
                nmax = 100, # samples rows from data to display, NULL means all rows/observations, 50 is good choice for computation time 
                gridSize = 13, # huge influence on computation time, set to 20 for higher resolution in pixel-plot
                nIce = 100, # Number of ice curves to be plotted, defaults to 30.
                vars = variables3[1:10]) # choose important variables here 
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_feature_influence_rf3_new.jpg', plot = viv_plot3, width = 20, height = 16)
?vivid::pdpVars()

# plot partial predictive plots for only one variable manually for quicker computation
variables <- names(sort(importance(rf1), decreasing = TRUE))
variables
pdp::partial(rf1, pred.var = 'cocoa', train = traindata, prob = T) %>% 
  ggplot(aes(x = cocoa, y = yhat)) +
  geom_line(size = 1, color = "blue") +
  theme_bw() +
  labs(title = "Partial Dependence of cocoa",
       x = "cocoa",
       y = "Predicted Probability")

pdp::plotPartial(pdp::partial(rf1, pred.var = variables, train = traindata, prob = T))

# make heatmap of predictors - I'm not sure if this is actually very helpful, but it breakes my computer!
viviHeatmap(mat = vivi(data = traindata, fit = rf1, response = 'Occu_Effort')) # yhat represents the Variable Interaction Value Importance - quantifies how much a pair of predictor variables interact in influencing the model's predictions


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Predict with all 3 random forest models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# predict and store data in data frame 
data$pred1 <- predict(rf1, data = data %>% st_drop_geometry())$predictions / max(traindata$Effort_manipulated)
data$pred2 <- predict(rf2, data = data %>% st_drop_geometry() %>% mutate(Effort = max(Effort)))$predictions #[,1]
data$pred3 <- predict(rf3, data = data %>% st_drop_geometry())$predictions

# confusion matrix - doesn't work yet
# library(caret)
# caret::confusionMatrix(data$pred1, reference = data$Occu)

# plot
plot(data %>% select(pred1, pred2, pred3, Occu))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Plot predictions  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot option 1
opt1 <- ggplot(data) +
  geom_sf(aes(fill = pred1), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Distributon in Gola, Sierra Leone",
       subtitle = 'Predictions from a random forest SDM with survey effort included in response variable', 
       fill = "Prediction", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()


# plot option 2
opt2 <- ggplot(data) +
  geom_sf(aes(fill = pred2), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Distributon in Gola, Sierra Leone",
       subtitle = 'Predictions from a random forest SDM with survey effort included as a predictor variable', 
       fill = "Prediction", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()

# plot option 3
opt3 <- ggplot(data) +
  geom_sf(aes(fill = pred3), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Distributon in Gola, Sierra Leone",
       subtitle = 'Predictions from a random forest SDM without considering survey effort', 
       fill = "Prediction", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()

# plot raw observations
raw<-ggplot(data) +
  geom_sf(aes(fill = Occu), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Occupancy in Gola, Sierra Leone",
       subtitle = 'Raw observations aggregated per grid cell', 
       fill = "Occurence", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()

# create an arranged plot which compares the different methods 
library(ggpubr)
ggarrange(opt1, opt2, opt3, raw)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_predictions_big_grid.jpg', plot = last_plot(), width = 15*2, height = 10*2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 10. Plot influence of predictor against prediction  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ggplot(data %>% filter(Occu == 1 | Effort > 0)) +
  geom_point(mapping = aes(y = Distance_settlement, x = pred1, colour = factor(Occu))) +
  labs(title = 'Distance to Settlements', 
       subtitle = 'Predictions from a random forst SDM with survey effort included in response variable', 
       x = 'Predicted Occupancy', y = 'Distance to Settlements', color = 'True Occurence') +
  theme_bw()
ggplot(data %>% filter(Occu == 1 | Effort > 0)) +
  geom_point(mapping = aes(y = Distance_large_river, x = pred1, colour = factor(Occu))) +
  labs(title = 'Distance to Large Rivers', 
       subtitle = 'Predictions from a random forst SDM with survey effort included in response variable', 
       x = 'Predicted Occupancy', y = 'Distance to Large Rivers', color = 'True Occurence') +
  theme_bw()
ggplot(data %>% filter(Occu == 1 | Effort > 0)) +
  geom_point(mapping = aes(y = water, x = pred1, colour = factor(Occu))) +
  labs(title = 'Water', 
       subtitle = 'Predictions from a random forst SDM with survey effort included in response variable', 
       x = 'Predicted Occupancy', y = 'Percent Water per Grid Cell', color = 'True Occurence') +
  theme_bw()
ggplot(data %>% filter(Occu == 1 | Effort > 0)) +
  geom_point(mapping = aes(y = forest, x = pred1, colour = factor(Occu))) +
  labs(title = 'Forest', 
       subtitle = 'Predictions from a random forst SDM with survey effort included in response variable', 
       x = 'Predicted Occupancy', y = 'Percent Forest per Grid Cell', color = 'True Occurence') +
  theme_bw()
ggplot(data %>% filter(Occu == 1 | Effort > 0)) +
  geom_point(mapping = aes(y = mean_elev, x = pred1, colour = factor(Occu))) +
  labs(title = 'Elevation', 
       subtitle = 'Predictions from a random forst SDM with survey effort included in response variable', 
       x = 'Predicted Occupancy', y = 'Mean Elevation per Grid Cell', color = 'True Occurence') +
  theme_bw()

