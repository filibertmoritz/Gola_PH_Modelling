#### pygmy hippo random forest approach to select important variables for integrated occupancy model 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in April 2025 


#### Background:  We have data from various PH surveys which can be mainly distinguished in camera trap surveys, transect surveys and opportunistic presence-only data.
####              We want to built an integrated occupancy model with these 2 or 3 data sources to a) map occupancy, b) explore the effect of the important variables on occupancy and c) potentially explore the temporal change in occupancy.
#### Issue:       We don't know which variables could drive occupancy and thus should be included as predictors into the integrated occupancy model. 
#### Solution:    Use a random forest model on all pooled data (from all 3 sources) to select important variables using the AUC approach        
#### Steps:       a) load elsewhere prepared data sets and bring them into one big df, 
####              b) load elsewhere prepared environmental covariates and remove highly correlated variables  
####              c) create/design an effort variable to include it into the random forest model 
####              d) fit RF model and get variable importance 

#### Literature:  1: Bradter et al. 2022 - https://www.cambridge.org/core/journals/environmental-data-science/article/variable-ranking-and-selection-with-random-forest-for-unbalanced-data/D00D9D74FA395B4FAC8886A84CC2FCCA, I mostly follow the recommendations from this paper
####              2: Genuer et al. 2010 - https://www.sciencedirect.com/science/article/abs/pii/S0167865510000954, general but old paper on variable selection using rf
####              3: Degenhardt et al. 2019 - https://academic.oup.com/bib/article/20/2/492/4554516?login=false, variable selection in omics using rf, not super relevant
####              4: Gregorutti et al. 2016 - https://link.springer.com/article/10.1007/s11222-016-9646-1, correlation and variable selection in rf
####              5: Hanberry 2024 - https://www.sciencedirect.com/science/article/pii/S1574954123004351, about correlation in SDMs and variable importance 


#### Hints:       1: The current version doesn't consider oppotunistic presence only data 
####              2: cforest from the party package might handle multicollinearity better 

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
library(party)
library(varImp) # probably only one of these is needed!
library(moreparty) # not sure if needed!
library(permimp)
library(scales)
library(tmap)
library(stringr)
# library(units)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. load data which has been prepared elsewhere #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# presence data 
# pres_opp <- fread(file = 'data/PH_prepared_pres_opp_data.csv', stringsAsFactors = T) %>% select(-V1)
pres_cam <- fread(file = 'data/PH_prepared_pres_cam_data.csv', stringsAsFactors = T) %>% select(-V1)
pres_transects <- fread(file = 'data/PH_prepared_pres_transect_data.csv', stringsAsFactors = T) %>% select(-V1)

# deploy and visit data 
deploy_cam_sf <- fread(file = 'data/PH_prepared_deploy_cam_data.csv', stringsAsFactors = T) %>% select(-V1) %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)
locs_transects_paths <- st_read('data/PH_prepared_transect_paths.shp')
locs_transects_data <- fread(file = 'data/PH_prepared_transect_paths.csv', stringsAsFactors = T) %>% select(-V1)
locs_transects <- locs_transects_paths %>% left_join(locs_transects_data, join_by(uniqueID)) # combine shp and csv file to one object 
rm(locs_transects_data, locs_transects_paths) # remove unneeded objects 

# grid and environmental covariates 
grid_sf <- st_read('data/PH_grid.shp') %>% st_as_sf() # laod shp and transform to sf 
envCovs <- fread('data/PH_prepared_env_covariates.csv', stringsAsFactors = T) %>% 
  select(-V1) %>% as_tibble()
envCovs_sf <- grid_sf %>% left_join(envCovs, join_by(CellID)) %>% st_as_sf()# transfer geometry from grid to envCovs
rm(envCovs)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Merge all presence data into one big data frame #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare opportunistic presences as sf
# head(pres_opp)
# str(pres_opp)
# pres_opp <- pres_opp %>% 
#   rename(UTM_X_meters = UTM_X_m, UTM_Y_meters = UTM_Y_m) %>% 
#   mutate(Obs_Method = 'Opportunistic') %>%
#   select(Project, Country, Obs_DateTime, UTM_X_meters, UTM_Y_meters, Sign, Obs_Method) 

# prepare transect presences as sf
pres_transects <- pres_transects %>% 
  mutate(Country = 'SierraLeone', 
         Obs_Method = 'Transect Survey')

# prepare camera trap presences as sf
pres_cam <- pres_cam %>% mutate(Obs_Method = 'Camera Trap') %>% 
  select(Project, Country, Obs_DateTime, UTM_X_meters , UTM_Y_meters, SiteID, Obs_Method)

# merge data together and save as sf
pres_sf <- bind_rows(pres_cam, 
                     # pres_opp,
                     pres_transects) %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create effort variable #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# effort variable for camera trap data 
deploy_cam_sf <- deploy_cam_sf %>% 
  mutate(Deployment_Time = as.numeric(difftime(Collection, Deployment, units = 'days'))) # calculate deployment time in days

effort_cam <- deploy_cam_sf %>% # sum up deployment time and number of cameras per cell
  st_join(grid_sf) %>% group_by(CellID) %>% 
  summarise(CameraEffort_Time = sum(Deployment_Time), 
            CameraEffort_N = n(), 
            .groups = 'keep') %>% 
  st_drop_geometry()

# effort variable for transect data 
effort_transects <- locs_transects %>% 
  st_intersection(grid_sf) %>%  # drops all the cells where no trancests are
  mutate(TransectEffort_Length = st_length(geometry, which = 'Euclidean')) %>%  # compute length of each river segment per grid cell, in m (which is the unit of the projected CRS)
  group_by(CellID) %>%  # group by grid cell
  summarise(TransectEffort_Length = as.numeric(sum(TransectEffort_Length, na.rm = TRUE)))%>%
  st_drop_geometry()

# effort variable for opportunistic presence-only data 

# distance weighted influence 
# weighted euclidean distance
# possibly the effort from other survey methods as benchmark where people went?

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

# remove unneeded effort objects 
rm(effort_cam, effort_transects)

# overview effort plot 
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


effort_sf %>% filter(Effort > 0) %>% nrow() # there are 313 cells where surveys were performed


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Create presence-absence information per grid cell #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

occu <- pres_sf %>% select(geometry) %>%
  st_join(effort_sf) %>%  # effort data transferred as well to possibly infer on 0's if there was effort, but no presence in a grid cell (only relevant if pres_opp is included too)
  st_drop_geometry() %>% 
  group_by(CellID) %>% 
  mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
  distinct() 
occu_sf <- grid_sf %>% left_join(occu, join_by(CellID)) # transfer to full grid 
occu_sf[is.na(occu_sf)] <- 0 # this replaces all NAs - places where no species where seen and no effort occured - to 0's


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Deal with highly correlated predictors #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare predictors
str(envCovs_sf)
pred <- envCovs_sf %>% 
  st_drop_geometry() %>% 
  select(-CellID, -area, -NDVI_2013, -EVI_2013, -matches("2021|2013"))  # remove unneeded variables and non-numeric variables 
colnames(pred) <- str_remove_all(colnames(pred), 'JRC_ann_changes_')
colnames(pred) <- str_remove_all(colnames(pred), 'JRC_transition_')

str(pred)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Bring all data together #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data overview
effort_sf %>% class()
envCovs_sf %>% class()
occu_sf %>% class()

# create one massive data frame 
data <- envCovs_sf %>% 
  left_join(effort_sf %>% st_drop_geometry(), join_by(CellID)) %>% 
  left_join(occu_sf %>% select(CellID, Occu) %>% st_drop_geometry(), join_by(CellID))

# create a training data set where either a presence was observed or surveys were performed!
traindata <- data %>% 
  filter(Occu == 1 | Effort > 0) %>% # this filters for all cases where Effort > 0 and also includes all grid cells with presences (to correctly handle presence-only data)
  st_drop_geometry()

# create a separate response variable which a) incorporates the Effort but b) doesn't become Inf from dividing by 0
traindata <- traindata %>% 
  mutate(# Effort_manipulated = if_else(Occu == 1 & Effort == 0, 0.0001, Effort), # if effort is 0 but an opportunistic presence occured, this ensures, that no Inf occures! ITS A BIT CRUDE
         CameraEffort_Time = as.numeric(CameraEffort_Time), 
         TransectEffort_Length = as.numeric(TransectEffort_Length),
         Occu_Effort = Occu*Effort) # divide Occu by Effort (if presence-only data considered - unse Effort_manipulated instead)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Train random forest model and compute variable importance #####
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




# get all predictors from data set
pred <- traindata %>% 
  select(-CellID, -area, -NDVI_2013, -EVI_2013, -matches("2021|2013|Effort"), -Occu) %>% # remove a few predictors which aren't needed
  colnames()

# built formulas
f_without_effort <- as.formula(paste('as.factor(Occu)', "~", paste(pred, collapse = " + "))) # classification RF without consideration of survey effort 
f_effort_pred <- as.formula(paste('as.factor(Occu)', "~", paste(pred, collapse = " + "), '+ Effort')) # classification RF with consideration of survey effort as predictor variable
f_occu_effort_response <- as.formula(paste('Occu_Effort', "~", paste(pred, collapse = " + "))) # regression random forest with consideration by modelling Occu-Effort ratio as response variable

# preparation to built fr in loop
formulas <- list(f_without_effort, f_effort_pred, f_occu_effort_response) # save all formulas in a list
names(formulas) <- c('f_without_effort', 'f_effort_pred', 'f_occu_effort_response')
models <- list() # initial list for model objects 
i <- 1 # as manual counter

# build random forests in a loop 
for(f in formulas){
  model_name <- paste0('cr', names(formulas)[i])
  models[[model_name]] <- cforest(formula =  f,
    data = traindata, 
    controls = cforest_control(ntree = 2000, 
                               mtry = length(attr(terms(f), "term.labels"))/2)) # calculates number of predictors and divides by 2
  i <- i+1 # count 1 further
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Calculate AUC Variable Importance for each RF #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# there are mainly 3 packages to calculate variable importance from cforest objects:
## 1. permimp - https://cran.r-project.org/web/packages/permimp/permimp.pdf
# permimp(crf1)
## 2. party - 
# party::varImpAUC(crf1, conditional = T)# or without conditional
## 3. moreparty - 
# # moreparty::fastvarImpAUC(crf1, conditional = T, parallel = T)# or without conditional



# be aware that this takes some time!

# calculate AUC Variable importance 10 times for each model 
results_vip <- data.frame(model = factor(), iter = numeric(), variable = character(), impoertance = numeric())

for(m_names in names(models)){
  m <- models[[m_names]]
  for(n in 1:10){ # this could also be done in nperm within varImpAUC, however it may not be 
    if(grepl("as\\.factor", as.character(m@data@formula$response[2]))){
      v <- moreparty::fastvarImpAUC(m, conditional = 0, parallel = T)} # v <- party::varimpAUC(m, conditional = T) # alternatively moreparty::varImpAUC(m, conditional = 0, parallel = T) or varImp::varImpAUC
    else{
      v <- party::varimp(m, conditional = T)} # use normal varImp since AUC is not possible in regression classifications 
    r <- data.frame(model = m_names, iter = n, variable = names(v), importance = v)
    results_vip <- rbind(results_vip, r)
  }
}

# calculate mean AUC varImp per model and variable
results_vip <- results_vip %>% 
  group_by(model, variable) %>% 
  mutate(mean_importance = mean(importance))

# make plot out of curiosity

# only the means 
ggplot(results_vip, aes(x = reorder(variable, mean_importance), y = mean_importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  facet_wrap(~ model, scales = "free") +
  labs(
    x = "Variable",
    y = "Mean AUC-based Importance",
    title = "Variable Importance by Model") +
  theme_minimal(base_size = 13)

# means and single varImps together
ggplot(results_vip) +
  # Grey jittered points: individual importance scores
  geom_point(aes(x = reorder(variable, mean_importance), y = importance),
             color = "grey20", alpha = 0.5, size = 1.5,
             position = position_jitter(width = 0.2)) +
  # Red vertical lines from 0 to mean importance
  geom_segment(aes(x = reorder(variable, mean_importance), xend = reorder(variable, mean_importance),
                   y = 0, yend = mean_importance),
               color = "red", linewidth = 1) +
  # Red mean points
  geom_point(aes(x = reorder(variable, mean_importance), y = mean_importance),
             color = "red", size = 1.5) +
  coord_flip() +
  facet_wrap(~ model, scales = "free") +
  labs(x = "Variable",
    y = "AUC-based Importance",
    title = "Variable Importance per Model",
    subtitle = "Grey dots: individual runs (n = 10); Red line + dot: mean importance") +
  theme_minimal(base_size = 13)

ggplot(results_vip) +
  geom_point(aes(x = reorder(variable, mean_importance), y = importance),
             color = "grey20", alpha = 0.5, size = 1.5,
             position = position_jitter(width = 0.2)) +
  geom_segment(aes(x = reorder(variable, mean_importance), xend = reorder(variable, mean_importance),
                   y = 0, yend = mean_importance),
               color = "red", linewidth = 1) +
  geom_point(aes(x = reorder(variable, mean_importance), y = mean_importance),
             color = "red", size = 1.5) +
  coord_flip() +
  facet_grid(. ~ model, scales = "free", switch = "y") +  # key change
  labs(x = "Variable",
       y = "AUC-based Importance",
       title = "Variable Importance per Model",
       subtitle = "Grey dots: individual runs (n = 10); Red line + dot: mean importance") +
  theme_minimal(base_size = 13) +
  theme(strip.placement = "outside")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 10. Forward variable selection #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# consecutively select best variables to add them to the model and compare via AUC 

# follow the structure which has been applied in Bradter et al. 2022 - https://www.cambridge.org/core/journals/environmental-data-science/article/variable-ranking-and-selection-with-random-forest-for-unbalanced-data/D00D9D74FA395B4FAC8886A84CC2FCCA

# this is not implemented here! 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 11. Recursive Feature Selection #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this makes things much easier 

# moreparty::FeatureSelection(Y = traindata$Occu, 
#                             X = traindata %>% select(-CellID, -area, -NDVI_2013, -EVI_2013, -matches("2021|2013|Effort"), -Occu), 
#                             method = 'RFE', ntree = 1000)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 12. Select all 10 best variables per RF #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select 10 best predictors from all models 
best_pred <- results_vip %>% group_by(model, variable) %>% 
  summarise(mean_importance = mean(importance)) %>%
  slice_max(order_by = mean_importance, n = 6) %>% pull(variable) %>% unique()
best_pred <- best_pred[best_pred != 'Effort'] # exclude effort
# best_pred <- best_pred[best_pred != 'Reserve_Type'] # exclude Reserve_Type which is categorial

# check for correlation between best predictors

traindata[, best_pred] %>% st_drop_geometry() %>% cor() %>% corrplot(type = 'lower')

names(traindata[, best_pred])








################################################################################
#########STOP HERE #############################################################
#################################################################################
################################################################################


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











#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 6.1. Correlation analysis using #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#################################################################################
##### FROM HERE CORRELATION ISSUE HAS TO BE REVISED ############################
#################################################################################

# calculate pairwise correlation and remove variables above certain threshold (between 0.5 to 0.7 - Dormann 2017)

# calculate correlation matrix 
cor.mat <- cor(as.matrix(pred[, sapply(pred, is.numeric)])) # uses only numeric predictors

# create corplot
library(corrplot)
corrplot(cor.mat,  type = 'lower')

# set threshold for high correlation
threshold <- 0.6 # conservative value

# get upper triangle of correlation matrix (to avoid duplicates)
upper_tri <- cor.mat2
upper_tri[lower.tri(upper_tri, diag = TRUE)] <- NA

# Find indices of highly correlated pairs
high_corr_idx <- which(abs(upper_tri) > threshold, arr.ind = TRUE) # abs handles both - positive and negative correlation

# bring data together in a df 
high_corr_pairs <- data.frame(
  var1 = rownames(upper_tri)[high_corr_idx[, 1]],
  var2 = colnames(upper_tri)[high_corr_idx[, 2]],
  correlation = cor.mat2[high_corr_idx])

# check result
print(high_corr_pairs)

# remove redundant predictors 
pred_new <- pred %>% select(-EVI_2020, -Undisturbed_tropical_moist_forest, -Degraded_tropical_moist_forest_Dec2017, -Tropical_moist_forest_regrowth_Dec2017, -Distance_thoroughfares, -Deforested_land_Dec2017, -Degraded_forest_long_duration_disturbance_after_2014, -Other_land_cover_Dec2017, -Degraded_forest_2_3_degradation_periods_after_2014)
cor.mat2 <- cor(as.matrix(pred_new[, sapply(pred_new, is.numeric)])) 
corrplot(cor.mat2)

# get names of new predictors 
new_predictors <- colnames(pred_new)

# ALTERNATIVE approach for selecting variables using the correlation - I don't really understand why certain variables are removed and others not! 

# find highly correlated variables 
library(caret)
red_var <- findCorrelation(cor.mat, cutoff = 0.5)
names(pred)[-red_var]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 6.2. Variance inflation factor #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# From the manual for the usdm package "Collinearity causes instability in parameter estimation in regression-type models. The VIF is based on the square of the multiple correlation coefficient resulting from regressing a predictor variable against all other predictor variables. If a variable has a strong linear relationship with at least one of the other variables, the correlation coefficient will be close to 1, and VIF for that variable would be large."
# Two options for VIFs: vifcor and vifstep. The latter is the option used by Geue & Thomassen (2020) Functions "exclude highly collinear variables in a stepwise procedure".
# "vifcor" first finds the pair of variables with the maximum linear correlation (greater than "th" threshold), and excludes the one with the greater VIF. The procedure is repeated until no variable has a high correlation coefficient (greater than threshold) with other variables.
# "vifstep" calculates the VIF for all variables, excludes the one with highest VIF (greater than threshold), and repeats the procedure until no variables with VIF greater than "th" remains.
# In Species Distribution Modelling the norm is to retain predictors with VIFs of < 10. In other modelling areas, then ViFs of 3 or 5 are more normal.

# Check the original paper for the package: Naimi et al. (2014) https://doi.org/10.1111/j.1600-0587.2013.00205.x
# Also Geue & Thomassen (2020), where it is used in a SDM https://doi.org/10.1002/ece3.6232

# load necessary package
library(usdm)

options(scipen = 999) # swith scientific notation off 
set.seed(123) # Not certain if there is random process? Results seem to vary?

# Predictors can only be numeric
VIF.STEP <- vifstep(pred[, sapply(pred, is.numeric)], th = 4)
VIF.COR <- vifcor(pred[, sapply(pred, is.numeric)], th = 0.6)

# have a look at outputs 
summary(VIF.STEP)
VIF.STEP@results
VIF.COR@results

# Returns a VIF object, examine different outputs.
VIF.STEP.MATRIX <- data.frame(VIF.STEP@corMatrix)
VIF.COR.MATRIX <- data.frame(VIF.COR@corMatrix)

# Check which variables are still included.
colnames(pred)
row.names(VIF.STEP.MATRIX)
row.names(VIF.COR.MATRIX)

VIF.COR@excluded

#################################################################################
##### UNTIL HERE CORRELATION ISSUE HAS TO BE REVISED ############################
#################################################################################
