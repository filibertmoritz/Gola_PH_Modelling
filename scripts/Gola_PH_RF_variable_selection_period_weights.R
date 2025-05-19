#### pygmy hippo random forest approach to select important variables for integrated occupancy model 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in April 2025 

# completely revised in May 2025 to remove unnecessary code chunks and improve readability 


#### Background:  We have data from two PH surveys which can be distinguished in camera trap surveys and transect surveys.
####              We want to built an integrated occupancy model with these 2 data sources to a) map occupancy, b) explore the effect of the important variables on occupancy and c) explore the temporal change in occupancy.
#### Issue:       We don't know which variables could drive occupancy and thus should be included as predictors into the integrated occupancy model. 
#### Solution:    Use a conditional random forest model on pooled data to select important variables using the variable importance measures, preferably an AUC based approach        
#### Steps:       a) load elsewhere prepared data sets and bring them into one big df, 
####              b) load elsewhere prepared environmental covariates 
####              c) create/design an effort variable to include it into the random forest model 
####              d) fit RF model and get variable importance 

#### Literature:  1: Bradter et al. 2022 - https://www.cambridge.org/core/journals/environmental-data-science/article/variable-ranking-and-selection-with-random-forest-for-unbalanced-data/D00D9D74FA395B4FAC8886A84CC2FCCA, I mostly follow the recommendations from this paper
####              2: Genuer et al. 2010 - https://www.sciencedirect.com/science/article/abs/pii/S0167865510000954, general but old paper on variable selection using rf
####              3: Degenhardt et al. 2019 - https://academic.oup.com/bib/article/20/2/492/4554516?login=false, variable selection in omics using rf, not super relevant
####              4: Gregorutti et al. 2016 - https://link.springer.com/article/10.1007/s11222-016-9646-1, correlation and variable selection in rf
####              5: Hanberry 2024 - https://www.sciencedirect.com/science/article/pii/S1574954123004351, about correlation in SDMs and variable importance 
####              6: Buston and Elith 2011 - https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2011.01803.x

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(sf)
library(terra)
library(tidyverse)
library(lubridate)
library(hms)
library(data.table)
# library(ranger)
library(party)
library(varImp) # probably only one of these is needed!
library(moreparty) # not sure if needed!
# library(permimp)
library(scales)
library(tmap)
library(stringr)
library(units)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. load data which has been prepared elsewhere #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# presence data 
pres_cam <- fread(file = 'data/PH_prepared_pres_cam_data.csv', stringsAsFactors = T) %>% select(-V1)
pres_transects <- fread(file = 'data/PH_prepared_pres_transect_data.csv', stringsAsFactors = T) %>% select(-V1)

# camera deployment data and transect visits data 
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
##### 3. Create yearly envCov data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create period column
envCovs_sf_period <- envCovs_sf %>% 
  mutate(Period = as.factor('2011-2017')) %>%
  rbind(envCovs_sf %>% mutate(Period = as.factor('2018-2025'))) 

# remove data from 2017 
envCovs_sf_period <- envCovs_sf_period %>%
  select(!matches('2017')) 

# create yearly siteCovs
yearly_envCovs <- sub('_2021', '', names(envCovs_sf_period)[grepl('2021', names(envCovs_sf_period))])

# merge columns depending on period for each yearly/period specific environmental covariate
for(y in yearly_envCovs){
  envCovs_sf_period <- envCovs_sf_period %>% 
    mutate(!!y := case_when(Period == '2011-2017' ~ .[[paste0(y, "_2013")]], 
                                   Period == '2018-2025' ~ .[[paste0(y, "_2021")]], 
                                   TRUE ~ NA_real_))
}

# remove unneeded columns - NDVI, EVI, JRC annual changes 
envCovs_sf_period <- envCovs_sf_period %>%
  select(-matches("^NDVI_.*_(2013|2021)$"), -matches("^EVI_.*_(2013|2021)$"), 
         -matches("^JRC_ann_changes_.*_(2013|2021)$"), 
         -JRC_transition_Undisturbed_tropical_moist_forest, -area) # remove this column because it is probably very similar to JRC_ann_changes_Undisturbed_tropical_moist_forest

# change column names to remove prefix (JRC....)
names(envCovs_sf_period) <- sub('JRC_ann_changes_|JRC_transition_', '', names(envCovs_sf_period))
str(envCovs_sf_period)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Merge all presence data into one big data frame #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare transect presences as sf
pres_transects <- pres_transects %>% 
  mutate(Country = 'SierraLeone', 
         Obs_Method = 'Transect Survey')

# prepare camera trap presences as sf
pres_cam <- pres_cam %>% mutate(Obs_Method = 'Camera Trap') %>% 
  select(Project, Country, Obs_DateTime, UTM_X_meters , UTM_Y_meters, SiteID, Obs_Method)

# merge data together and save as sf
pres_sf <- bind_rows(pres_cam, 
                     pres_transects) %>% 
  mutate(Period = as.factor(if_else(year(Obs_DateTime) <= 2017, '2011-2017', '2018-2025'))) %>%  # create column with time 2 distinct time periods 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create effort variable #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate effort variable for camera trap data 
deploy_cam_sf <- deploy_cam_sf %>% 
  mutate(Deployment_Time = as.numeric(difftime(Collection, Deployment, units = 'days')), # calculate deployment time in days
         Period = as.factor(if_else(year(Collection) <= 2017, '2011-2017', '2018-2025'))) # add period variable 

effort_cam <- deploy_cam_sf %>% # sum up deployment time and number of cameras per cell and period
  st_join(grid_sf) %>% group_by(CellID, Period) %>% 
  summarise(CameraEffort_Time = sum(Deployment_Time), 
            CameraEffort_N = n(), 
            .groups = 'keep') %>% 
  st_drop_geometry()

# calculate effort variable for transect data 
effort_transects <- locs_transects %>% 
  mutate(Period = as.factor(if_else(year(DateTime_End) <= 2017, '2011-2017', '2018-2025'))) %>% # add period variable
  st_intersection(grid_sf) %>%  # drops all the cells where no transects are
  mutate(TransectEffort_Length = set_units(st_length(geometry, which = 'Euclidean'), 'm')) %>%  # compute length of each river segment per grid cell, in m (which is the unit of the projected CRS)
  group_by(CellID, Period) %>%  # group by grid cell and period
  summarise(TransectEffort_Length = as.numeric(sum(TransectEffort_Length, na.rm = TRUE)))%>%
  st_drop_geometry()

# bring effort data together, beforehand create a template grid with a CellID for each Time period
effort_sf <- envCovs_sf_period %>% select(CellID, Period) %>%
  left_join(effort_cam, join_by(CellID, Period)) %>% 
  left_join(effort_transects, join_by(CellID, Period))
effort_sf[is.na(effort_sf)] <- 0 # replace all NAs with 0 

# summarise data to calculate ratio between transect and camera effort effectiveness to presences // THIS WORKS WITH RAW DATA which is not filtered by independence!
effort_comparison <- effort_sf %>% 
  st_drop_geometry() %>%
  group_by(Period) %>% 
  summarise(overall_transect_effort = sum(TransectEffort_Length), # in m 
            overall_camera_time_effort = sum(CameraEffort_Time), # in days
            overall_camera_N_effort = sum(CameraEffort_N)) %>% # Number of Cameras per grid cells
  mutate(pres_trans = c(nrow(pres_transects %>% filter(year(Obs_DateTime)<= 2017)), # this calculates the number of transect presences in the first time period 
                        nrow(pres_transects %>% filter(year(Obs_DateTime)> 2017))), # this calculates the number of transect presences in the second time period 
         pres_cam = c(nrow(pres_cam %>% filter(year(Obs_DateTime)<= 2017)), # this calculates the number of camera presences in the first time period 
                      nrow(pres_cam %>% filter(year(Obs_DateTime)> 2017)))) # this calculates the number of camera presences in the second time period 

ratio <- effort_comparison %>% mutate(ratio_transects = pres_trans/overall_transect_effort, # per m surveyed transect we have these many PH signs 
                 ratio_camera_time = pres_cam/overall_camera_time_effort, # per day camera deployment we have these many PH sightings 
                 ratio_camera_N = pres_cam/overall_camera_N_effort) # per camera we have these many PH sightings

# calculate factor to inflate the camera effort to the scale of meters for each time period 
inflation <- ratio %>% mutate(cam_time_inflation = ratio_camera_time/ratio_transects) %>% pull(cam_time_inflation)
# these factors mean, that one had to survey this many meters to equal out one camera trapping day for the first and the second time period

# adjust effort across survey methods 
effort_sf <- effort_sf %>% 
  mutate(CameraEffort_Time_inflated = if_else(Period == '2011-2017', CameraEffort_Time*inflation[1], CameraEffort_Time*inflation[2]), # use inflation factor for different time periods
         Effort_unscaled = CameraEffort_Time_inflated+TransectEffort_Length, 
         Effort_scaled = rescale(Effort_unscaled, to = c(0,1))) # this scales across periods, not within periods!
  
# this is the old effort calculation

## calculate one effort index with rescaled variables to between 0 and 1
# effort_sf <- effort_sf %>%
#  mutate(CameraEffort = (rescale(as.numeric(CameraEffort_Time), to = c(0,1)) + 
#                           rescale(as.numeric(CameraEffort_N), to = c(0,1)) / 2), # one Camera Trap effort variable as mean from both N and Time
#          TransectEffort = rescale(as.numeric(TransectEffort_Length), to = c(0,1)), # rescaled Transect length per cell
#          Effort = (CameraEffort + (TransectEffort))/2) # mean of both efforts


# remove unneeded effort objects 
rm(effort_cam, effort_transects, ratio, inflation)

# overview effort plot 
tmap_mode("view")  # interactive mode
tm_shape(effort_sf %>% filter(Period == '2018-2025')) + 
  tm_polygons(fill = "Effort_scaled", lwd = .2, 
              fill.scale = tm_scale_continuous(values = "viridis", midpoint = 0)) +
tm_shape(locs_transects) + 
  tm_lines(col = 'red', lwd = 1.5) +
tm_shape(deploy_cam_sf) +
  tm_dots(fill = 'orange', size = .3) +
tm_shape(pres_sf) +
  tm_symbols(fill = 'brown', shape = 3)


effort_sf %>% filter(Effort_scaled > 0) %>% nrow() # there are 420 cells where surveys were performed across periods


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Create presence-absence information per grid cell and period #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

occu_list <- list() # create input list
for(p in unique(pres_sf$Period)){
  occu_list[[p]] <- pres_sf %>% filter(Period == p) %>% # filter presence df for correct period
    select(geometry) %>%
    st_join(effort_sf %>% filter(Period == p)) %>%  # effort data transferred as well to possibly infer on 0's if there was effort, but no presence in a grid cell (only relevant if pres_opp is included too)
    st_drop_geometry() %>% 
    group_by(CellID) %>% 
    mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
    distinct()
}
occu <- bind_rows(occu_list) # bind all df from list together
rm(occu_list) # remove list which isn't needed anymore

occu_sf <- envCovs_sf_period %>% select(CellID, Period) %>% # get grid for each time period
  left_join(occu, join_by(CellID, Period)) # transfer presences per grid cell and period to full grid 
occu_sf[is.na(occu_sf)] <- 0 # this replaces all NAs - places where no species where seen and no effort occured - to 0's


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Bring all data together #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data overview
effort_sf %>% str()
envCovs_sf %>% class()
occu_sf %>% class()

# create one massive data frame with all data
data <- envCovs_sf_period %>% 
  left_join(effort_sf %>% st_drop_geometry(), join_by(CellID, Period)) %>% 
  left_join(occu_sf %>% select(CellID, Occu, Period) %>% st_drop_geometry(), join_by(CellID, Period))

# create a training data set where either a presence was observed or surveys were performed!
traindata <- data %>% 
  filter(Occu == 1 | Effort_scaled > 0) %>% # this filters for all cases where Effort > 0 and also includes all grid cells with presences (to correctly handle presence-only data)
  st_drop_geometry()


# create a separate response variable which a) incorporates the Effort but b) doesn't become Inf from dividing by 0
traindata <- traindata %>% 
  mutate(CameraEffort_Time = as.numeric(CameraEffort_Time), 
         TransectEffort_Length = as.numeric(TransectEffort_Length),
         Occu_Effort = Occu*Effort_scaled) # divide Occu by Effort (if presence-only data considered - unse Effort_manipulated instead)

# the relevant variables are now: 
#     1. Occu_Effort - the occupancy multiplied by effort (Occu-Effort ratio), 
#     2. TransectEffort_Length (Effort for all transects in m), 
#     3. CameraEffort_Time - camera trapping effort in days

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Train random forest model and compute variable importance #####
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
traindata <- traindata  %>% st_drop_geometry() # remove geometry
pred <- names(traindata)[!grepl("Occu|Effort|Period|CellID", names(traindata))] # get all predictor names which are not connected to effort or occu

# built formulas and set parameters
f_without_effort <- as.formula(paste('as.factor(Occu)', "~", paste(pred, collapse = " + "))) # classification RF without consideration of survey effort 
f_effort_pred <- as.formula(paste('as.factor(Occu)', "~", paste(pred, collapse = " + "), '+ TransectEffort_Length + CameraEffort_Time')) # classification RF with consideration of survey effort as predictor variable
f_occu_effort_response <- as.formula(paste('Occu_Effort', "~", paste(pred, collapse = " + "))) # regression random forest with consideration by modelling Occu-Effort ratio as response variable

# sort data and create weights matrix
traindata <- traindata %>% arrange(Period, CellID) %>% select(CellID, Period, everything())
ntree <- 2000
weightsmatrix <- matrix(NA, nrow = nrow(traindata), ncol = ntree) # initialize matrix
period <- levels(traindata$Period) # get different periods
treeseq <- rep(period, ntree/2)
for(i in 1:ntree){weightsmatrix[,i] <- if_else(traindata$Period == treeseq[i], 0, 1)} # create weightmatrix with 0 for first period and 1 for other period, next column the other way round

weightsmatrix[,1]; weightsmatrix[,2] # check that everything went okay - should be complimentary


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
    weights = weightsmatrix,
    controls = cforest_control(ntree = ntree, 
                               mtry = length(attr(terms(f), "term.labels"))/2)) # calculates number of predictors and divides by 2
  i <- i+1 # count 1 further
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Extract and calculate performance measures for all 3 randomforests #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(pROC) # needed for ROC and AUC calculation

# confusion matrix
conf_matrix <- table(traindata$Occu, predict(models$crf_without_effort, OOB = TRUE))

# calculate oob error
oob_error <- 1 - (sum(diag(conf_matrix)) / sum(conf_matrix)) # diag extracts the correctly classified counts from the diagonal of the confusion matrix

# calculate ROC (true-positive-rate as y-axis and false-positive-rate as x-axis) and AUC
roc(traindata$Occu, as.numeric(predict(models$crf_without_effort, OOB = TRUE)))


# predict using the OOB = T and extract predicted probabilities from list
prob <- predict(models$crf_effort_pred, OOB = TRUE, type = 'prob')
pred_prob <- numeric(length(prob)) # create input vector
for (i in 1:length(prob)) {
  pred_prob[i] <- prob[[i]][2] # 1 is 0 (absence) and 2 is 1 (presence)
}

# calculate ROC and AUC from predicted probabilities of the classification tree
roc(traindata$Occu, pred_prob)
auc(roc(traindata$Occu, pred_prob))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Predict Occupancy across Gola as validation #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pred_occu <- data.frame(CellID = numeric(), Occupancy = numeric(), Period = factor(), Model = factor())

for(m in 1:length(models)){
  raw_pred <- predict(models[[m]], newdata = data, type = 'prob')
  
  pred_vect <- numeric(length(prob)) # create input vector
  for (i in 1:length(raw_pred)) {
    pred_vect[i] <- if_else(grepl('crf_occu_effort_response', names(models)[m]), raw_pred[[i]][1],raw_pred[[i]][2])} # 1 is 0 (absence) and 2 is 1 (presence)
  
  p <- data.frame(Occupancy = pred_vect, 
                  Period = factor(data$Period), 
                  CellID = data$CellID, 
                  Model = factor(names(models)[m]))
  pred_occu <- rbind(pred_occu, p)
}

raw_pred[[1]][1]

# make a better df for plotting
pred_occu_plot <- envCovs_sf_period %>% # add geometry
  select(CellID, Period) %>% 
  left_join(pred_occu, join_by(CellID, Period)) %>% 
  group_by(Model) %>% 
  mutate(Occupancy_scaled = rescale(Occupancy,to = c(0,1)))  %>% # scale to 0 and 1 within each model
  mutate(Model = case_when(Model == 'crf_without_effort'~ 'a) Model without consideration of effort.', 
                           Model == 'crf_occu_effort_response' ~ 'b) Model with occupancy-effort ratio as response variable.', 
                           Model == 'crf_effort_pred' ~ 'c) Model with effort as predictor variable.'))
  
# plot predictions - different predictions for different periods originate from the usage of different (period-specific) covariated per period
ggplot(pred_occu_plot) +
  geom_sf(aes(fill = Occupancy_scaled), alpha = 1) +
  facet_grid(Period~Model) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Distributon in Gola, Sierra Leone",
       subtitle = 'Predictions from multiple different conditional random forest SDM`s with period as random effect', 
       fill = "Rescaled \nPredicted \nOccupancy \n", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()


################################################################################
##### CONTINUE HERE ############################################################
################################################################################



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Calculate AUC Variable Importance for each RF #####
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
  mutate(mean_importance = mean(importance)) %>% 
  group_by(model) %>% 
  mutate(variation = abs(min(importance)), 
         # rel_mean_importance = rescale(mean_importance, c(0,1))
         )


# plot variable importance and save
plot <- results_vip %>% 
  mutate(Model = case_when(model == 'crf_without_effort'~ 'a) Model without consideration of effort.', 
                           model == 'crf_occu_effort_response' ~ 'b) Model with occupancy-effort ratio as response variable.', 
                           model == 'crf_effort_pred' ~ 'c) Model with effort as predictor variable.'), 
         Variable = case_when(variable == 'Distance_large_river' ~ 'Distance to large river', 
                              variable == 'mean_elev' ~ 'Elevation', 
                              variable == 'river_density_med_large'~ 'Density of medium and large rivers', 
                              variable == 'Distance_road'~ 'Distance to roads', 
                              variable == 'Reserve_Type'~ 'Reserve Type', 
                              .default = variable))  %>% 
  ggplot() +
  geom_point(aes(x = reorder(Variable, mean_importance), y = importance, color = "Individual replicates"), alpha = 0.5, size = 1.5,
             position = position_jitter(width = 0.2)) +
  geom_segment(aes(x = reorder(Variable, mean_importance), xend = reorder(Variable, mean_importance),
                   y = 0, yend = mean_importance, color = "Mean importance"), linewidth = .9) +
  geom_point(aes(x = reorder(Variable, mean_importance), y = mean_importance, color = "Mean importance"), size = 1.8) +
  geom_hline(aes(yintercept = variation, color = "Random variation"), linetype = 2) +
  coord_flip() +
  facet_grid(. ~ Model, scales = "free", switch = "y") + 
  scale_color_manual(name = "Legend",
                     values = c("Mean importance" = "firebrick", 
                                "Mean importance" = "firebrick", 
                                "Individual replicates" = "gray30",
                                "Random variation" = "grey50")) +
  labs(x = "Variable",
       y = "Relative Variable Importance",
       title = "Variable Importance") +
  theme_bw(base_size = 13) +
  theme(legend.background = element_rect(color = "black"), 
        legend.position = c(0.9, 0.15))
ggsave(plot = plot, filename = 'output/plots/PH_crf_variable_selection.jpg', width = 16, height = 9)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Select all 5 best variables per RF #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select 10 best predictors from all models 
best_pred <- results_vip %>% group_by(model, variable) %>% 
  summarise(mean_importance = mean(importance)) %>%
  slice_max(order_by = mean_importance, n = 6) %>% pull(variable) %>% unique()
best_pred <- best_pred[best_pred != 'Effort'] # exclude effort
best_pred



# check for correlation between best predictors
library(corrplot)

traindata[, best_pred] %>% st_drop_geometry()  %>% select(-Reserve_Type) %>% cor() %>% corrplot(type = 'lower')



