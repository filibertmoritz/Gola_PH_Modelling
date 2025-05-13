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


## this checks out what steffen suggested 
#ntree <- 1500
#weightsmatrix <- matrix(NA, nrow = dim(pres_sf)[1], ncol = ntree) # this creates a matrix with ntree columns and number of sites for rows - for every tree there is a column with all sites 

## fill in the different periods (2011-2017 is 0 and 2018-2025 is 1)
#names <- levels(pres_sf$Period) # use period as grouping factor 
#treeseq <- rep(1:5, 300)
#for(i in 1:ntree){
#  weightsmatrix[,i] <- if_else(pres_sf$Period == names)
#}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create effort variable #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# effort variable for camera trap data 
deploy_cam_sf <- deploy_cam_sf %>% 
  mutate(Deployment_Time = as.numeric(difftime(Collection, Deployment, units = 'days')), # calculate deployment time in days
         Period = as.factor(if_else(year(Collection) <= 2017, '2011-2017', '2018-2025'))) 

effort_cam <- deploy_cam_sf %>% # sum up deployment time and number of cameras per cell per period
  st_join(grid_sf) %>% group_by(CellID, Period) %>% 
  summarise(CameraEffort_Time = sum(Deployment_Time), 
            CameraEffort_N = n(), 
            .groups = 'keep') %>% 
  st_drop_geometry()

# effort variable for transect data 
effort_transects <- locs_transects %>% 
  mutate(Period = as.factor(if_else(year(DateTime_End) <= 2017, '2011-2017', '2018-2025'))) %>%
  st_intersection(grid_sf) %>%  # drops all the cells where no trancests are
  mutate(TransectEffort_Length = set_units(st_length(geometry, which = 'Euclidean'), 'm')) %>%  # compute length of each river segment per grid cell, in m (which is the unit of the projected CRS)
  group_by(CellID, Period) %>%  # group by grid cell
  summarise(TransectEffort_Length = as.numeric(sum(TransectEffort_Length, na.rm = TRUE)))%>%
  st_drop_geometry()

# bring effort data together, beforehand create a template grid with a CellID for each Time period
effort_sf <- grid_sf %>% 
  mutate(Period = as.factor('2011-2017')) %>% rbind(grid_sf %>% mutate(Period = as.factor('2018-2025'))) %>%
  left_join(effort_cam, join_by(CellID, Period)) %>% 
  left_join(effort_transects, join_by(CellID, Period))
effort_sf[is.na(effort_sf)] <- 0 # replace all NAs with 0 

# summarise data to calculate ratio between transect and camera effort effectiveness to presences // THIS WORKS WITH RAW DATA which is not filtered by independence!
effort_comparison <- effort_sf %>% 
  st_drop_geometry() %>%
  group_by(Period) %>% 
  summarise(overall_transect_effort = sum(TransectEffort_Length), # in m 
            overall_camera_time_effort = sum(CameraEffort_Time), # in days
            overall_camera_N_effort = sum(CameraEffort_N)) %>% 
  mutate(pres_trans = c(nrow(pres_transects %>% filter(year(Obs_DateTime)<= 2017)), # this calculates the number of transect presences in the first time period 
                        nrow(pres_transects %>% filter(year(Obs_DateTime)> 2017))), # this calculates the number of transect presences in the second time period 
         pres_cam = c(nrow(pres_cam %>% filter(year(Obs_DateTime)<= 2017)), # this calculates the number of camera presences in the first time period 
                      nrow(pres_cam %>% filter(year(Obs_DateTime)> 2017)))) # this calculates the number of camera presences in the second time period 

ratio <- effort_comparison %>% mutate(ratio_transects = pres_trans/overall_transect_effort, # per m surveyed transect we have these many PH signs 
                 ratio_camera_time = pres_cam/overall_camera_time_effort, # per day camera deployment we have these many PH sightings 
                 ratio_camera_N = pres_cam/overall_camera_N_effort) # per camera we have these many PH sightings

################################################################################
##### FILL SOMETHING APPROPRIATE IN ############################################
#################################################################################

# calculate factor to inflate the camera effort to the scale of meters for each time period 
inflation <- ratio %>% mutate(cam_time_inflation = ratio_camera_time/ratio_transects) %>% pull(cam_time_inflation)
# these factors mean, that one had to survey this many meters to equal out one camera trapping day

# calculate effort across survey methods 
effort_sf <- effort_sf %>% 
  mutate(CameraEffort_Time_inflated = if_else(Period == '2011-2017', CameraEffort_Time*inflation[1], CameraEffort_Time*inflation[2]), 
         Effort_unscaled = CameraEffort_Time_inflated+TransectEffort_Length, 
         Effort_scaled = rescale(Effort_unscaled, to = c(0,1))) # this scales across periods, not within periods!
  


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


effort_sf %>% filter(Effort_scaled > 0) %>% nrow() # there are 420 cells where surveys were performed in both periods


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Create presence-absence information per grid cell and period #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#################################################################################
###### SOMEHOWE THIS DOES NOT REALLY WORk #######################################
##################################################################################

# even if there is no effort in this time period, the loop fills in a 1 for a presence if there is another presence in the other period

occu_list <- list()
for(p in unique(pres_sf$Period)){
  occu_list[[p]] <- pres_sf %>% filter(Period == p) %>%
    select(geometry) %>%
    st_join(effort_sf) %>%  # effort data transferred as well to possibly infer on 0's if there was effort, but no presence in a grid cell (only relevant if pres_opp is included too)
    st_drop_geometry() %>% 
    group_by(CellID) %>% 
    mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
    distinct()
}

occu <- bind_rows(occu_list)

################################################################################
##### unTIL HERE ################################
################################################################################

occu <- pres_sf %>% select(geometry) %>%
  st_join(effort_sf) %>%  # effort data transferred as well to possibly infer on 0's if there was effort, but no presence in a grid cell (only relevant if pres_opp is included too)
  st_drop_geometry() %>% 
  group_by(CellID) %>% 
  mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
  distinct() 
occu_sf <- grid_sf %>% left_join(occu, join_by(CellID, Period)) # transfer to full grid 
occu_sf[is.na(occu_sf)] <- 0 # this replaces all NAs - places where no species where seen and no effort occured - to 0's


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Bring all data together #####
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


names(traindata)

# get all predictors from data set
traindata <- traindata  %>% 
  st_drop_geometry() %>% 
  select(-CellID, -area, -JRC_transition_Undisturbed_tropical_moist_forest,-matches("2013|2021"))  # remove unneeded variables 
colnames(traindata) <- str_remove_all(colnames(traindata), 'JRC_ann_changes_|JRC_transition_')
pred <- names(traindata)[!grepl("Occu|Effort", names(traindata))] # get all predictor names which are not connected to effort or occu

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



