#### Pygmy hippo single season integrated occupancy model in spOccupancy  
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025

# data has roughgly been prepared elsewhere

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(data.table)
library(spOccupancy)
library(tmap)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

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
##### 3. Prepare presence camera trap data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# create visits from camera trap deployment data 
deploy_cam_visit <- deploy_cam %>%
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629) %>% # create an sf
  st_join(grid_sf) %>%
  rowwise() %>%
  mutate(Visit_start = list(seq(Deployment, Collection, by = "20 days"))) %>%  # split in 20 days intervals, however, there are mostly days remaining as tail
  unnest(Visit_start) %>% 
  group_by(SiteID, Deployment, Collection) %>% # combination of SiteID, Deployment and Collection works as deploymentID
  mutate(Visit_end = lead(Visit_start, default = last(Collection)), 
         Visit_start = as.POSIXct(paste0(Visit_start, ' 12:00:00'), format = '%Y-%m-%d %H:%M:%OS'), # make meaningful time up where cameras could have been deployed and collected to avoid issues with matching pictures
         Visit_end = as.POSIXct(paste0(Visit_end, ' 12:00:00'), format = '%Y-%m-%d %H:%M:%OS'), 
         Visit_length = difftime(Visit_end, Visit_start, units = 'days')) %>%  
  filter(!Visit_length == 0) %>% # remove cases where 0 days were left as overshoot, 17 cases
  arrange(CellID, SiteID, Deployment) %>%
  group_by(CellID) %>% # CellID instead of SiteID
  mutate(Cell_visit = row_number()) # this creates a consecutive number for each visit for each site across deployments

# transfer presences to deployment data, this assumes, that the observations lie within the deployment period
#### THIS IS A CRITICAL STEP AND I WANT THIS TO BE DOUBLECHECKED BY SOMEONE ELSE!
deploy_cam_visit_occu <- deploy_cam_visit %>% st_drop_geometry() %>%
  ungroup() %>%
  left_join(pres_cam %>% select(Project, SiteID, Obs_DateTime), join_by(Project, SiteID), multiple = 'all', relationship = 'many-to-many') %>% 
  filter(Visit_start <= Obs_DateTime & Obs_DateTime <= Visit_end) %>% 
  group_by(SiteID, Cell_visit) %>% 
  summarise(Presences = n_distinct(Obs_DateTime), .groups = 'drop') %>% # this produces a count of pictures, but if distinct (not just multiple pictures in a row) detections of pygmy hippos are needed, use n_distinct() here 
  right_join(deploy_cam_visit %>% st_drop_geometry(), join_by(SiteID, Cell_visit)) %>% 
  mutate(Presences = if_else(is.na(Presences), 0, Presences), 
         Occu = if_else(Presences >= 1, 1, 0)) %>%
  arrange(CellID, SiteID, Deployment, Cell_visit) 


# there is an issue with one camera trapping site which is outside of the study area
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% filter(!is.na(CellID )) # this filters out one camera trap which lies outside the study area
deploy_cam_visit <- deploy_cam_visit %>% filter(!is.na(CellID))

# quick plot
tm_shape(grid_sf) +
  tm_polygons()+
tm_shape(deploy_cam_visit) +
  tm_dots(fill = 'red', size = 0.5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Prepare transect survey data   #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this option is not the best, but there could be other ways by frst sampling points on the line and then splitting the lines up at these points
# st_line_sample and lwgeom::st_split()

library(stplanr)

# split transects into subtransects which are interpreted as visits 
locs_transects <- locs_transects %>% 
  st_intersection(grid_sf) %>%
  st_cast(to = 'MULTILINESTRING') %>% # first transform everything to MULTILINESTRING
  st_cast(to ='LINESTRING') %>% # transfer everything pack to linestring
  group_by(uniqueID, CellID) %>% 
  # mutate(transect_length = st_length(geometry, which = 'Euclidean')) %>% 
  line_segment(segment_length = 700) %>% # segment length in the crs unit, which is meters here 
  group_by(uniqueID, CellID) %>% 
  mutate(subtransect = as.factor(row_number())) %>% 
  group_by(uniqueID,CellID,subtransect) %>% 
  mutate(transect_length = st_length(geometry, which = 'Euclidean'))

# plot result  
tmap_mode("view")  # interactive mode
tm_shape(grid_sf) +
  tm_polygons(fill = tm_const()) + 
tm_shape(locs_transects) + 
  tm_lines(col = 'subtransect', lwd = 1.5, id = 'uniqueID') 

# link the presences to the subtransects 
pres_transects_sf <- pres_transects %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629) %>% 
  st_transform(crs = 32629)

unique_ID <- unique(locs_transects$uniqueID) # get uniqueIDs for each transect
joined_pres <- st_join(pres_transects_sf %>% filter(uniqueID == unique_ID[1]), # create df to store data with first iteration
                       locs_transects %>% filter(uniqueID == unique_ID[1]) %>% select(subtransect, CellID), 
                       st_nearest_feature)
for(i in 2:length(unique_ID)) { # loop over all other IDs 
  df <- st_join(pres_transects_sf %>% filter(uniqueID == unique_ID[i]), 
                locs_transects %>% filter(uniqueID == unique_ID[i]) %>% select(subtransect, CellID), 
                st_nearest_feature)
  joined_pres <- bind_rows(joined_pres, df) # bind results together 
}
pres_transects_sf <- joined_pres
rm(joined_pres)

# check, that everything went okay 
tm_shape(grid_sf) +
  tm_polygons() +
tm_shape(locs_transects) + 
  tm_lines(col='subtransect', lwd = 2, palette=c("blue", "darkgreen", "red", "lightblue", 'orange', 'purple')) +
tm_shape(pres_transects_sf) +
  tm_dots(col='subtransect', palette=c("blue", "darkgreen", "red", "lightblue", 'orange', 'purple'), size=0.5)

# add presence absence information to locs transect df
occu_transects <- pres_transects_sf  %>% 
  left_join(locs_transects %>% st_drop_geometry(), join_by(uniqueID, subtransect, CellID)) %>% 
  group_by(CellID, uniqueID, subtransect) %>% 
  summarise(Presences = n(), .groups = 'drop') # this produces a count of observations per subtransect in each gridCell, but if distinct (not just multiple pictures in a row) detections of pygmy hippos are needed, use n_distinct() here 
  # st_cast(to = 'MULTIPOINT')

locs_transects <- locs_transects %>% 
  left_join(occu_transects %>% st_drop_geometry(), join_by(CellID, uniqueID, subtransect)) %>% 
  mutate(Presences = if_else(is.na(Presences), 0, Presences),
         Occu = if_else(Presences >= 1, 1, 0)) %>% 
  group_by(CellID) %>% 
  mutate(Cell_visit = row_number())


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Bring data together for spOccupancy ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# transfer transect data into grid
# transect_grid <- grid_sf %>% 
#   left_join(locs_transects %>% st_drop_geometry(), join_by(CellID)) %>%
#   group_by(CellID) %>% 
#   mutate(Cell_visit = row_number()) %>% 
#   st_drop_geometry() %>% 
#   drop_na(Occu)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1 Det Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate a few det covariates, improve season by assigning the season where the majority was surveyed
locs_transects <- locs_transects %>% ungroup() %>%
  mutate(Julian_Date_Start_Transect = lubridate::yday(DateTime_Start), # improve by taking a mean date?
         Year_Transect = as.numeric(year(DateTime_Start)), # survey year coded as numeric or factor
         Season_Transect= as.factor(if_else(month(DateTime_Start) %in% 5:10, 'Wet', 'Dry')), # wet season from May (05) to October (10), https://doi.org/10.51847/8Wz28ID8Mn
         Transect_Length = as.numeric(transect_length), 
         Project_Transect = as.factor(Project))  # %>% select(-transect_length) 
locs_transects_scaled <- locs_transects %>% mutate(across(c(Year_Transect, Julian_Date_Start_Transect, Transect_Length), ~ scale(.) %>% as.vector())) # scale numeric variables

deploy_cam_visit_occu <- deploy_cam_visit_occu %>% ungroup()%>% 
  mutate(Project_Camera = as.factor(Project), 
         Julian_Date_Start_Camera = yday(Visit_start), 
         Year_Camera = as.numeric(year(Visit_start)), # survey year coded as numeric or factor?
         Season_Camera = as.factor(if_else(month(Visit_start) %in% 5:10, 'Wet', 'Dry')),
         Trapping_Days = as.numeric(Visit_length)) # %>% select(-transect_length) # %>% select(-Visit_length) 
deploy_cam_visit_occu_scaled <- deploy_cam_visit_occu %>% mutate(across(c(Year_Camera, Julian_Date_Start_Camera, Trapping_Days), ~ scale(.) %>% as.vector())) # scale numeric variables

# extract det.covs data
names(locs_transects_scaled)
det.variables.transect <- c('Project_Transect', 'Julian_Date_Start_Transect', 'Transect_Length', 'Season_Transect', 'Year_Transect') # fill in all variables that are interesting
det.covs.transect <- list()
for(det.var in det.variables.transect){
  dat <- locs_transects_scaled %>% st_drop_geometry() %>%
    ungroup() %>% 
    select(CellID, Cell_visit, !!sym(det.var)) %>% 
    pivot_wider(names_from = Cell_visit, names_prefix = 'Visit_',
                values_from = !!sym(det.var))
  det.covs.transect[[det.var]] <- dat[,2:ncol(dat)] # dat # for troubleshooting use solely dat to keep CellID
  rm(dat)
}

names(deploy_cam_visit_occu_scaled)
det.variables.camera <- c('Project_Camera', 'Julian_Date_Start_Camera', 'Trapping_Days', 'Season_Camera', 'Year_Camera') # fill in all variables which could be interesting
det.covs.camera <- list()
for(det.var in det.variables.camera){
  dat <- deploy_cam_visit_occu_scaled %>% ungroup() %>% 
    select(CellID, Cell_visit, !!sym(det.var)) %>% 
    pivot_wider(names_from = Cell_visit, names_prefix = 'Visit_',
                values_from = !!sym(det.var))
  det.covs.camera[[det.var]] <- dat[,2:ncol(dat)] # dat # for troubleshooting use solely dat to keep CellID
  rm(dat)
}

# bind det.covs for both camera trapping and transect data together in one list 
det.covs <- list(det.covs.transect, det.covs.camera)
det.covs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.2 Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# these are environmental covariates, in spOccupancy called occ.covs 

# for transect data 
sites.transect <- unique(locs_transects$CellID) 
occ.covs.transect <- envCovs %>%
  mutate(across(-CellID, ~ scale(.) %>% as.vector())) %>% # scale all numeric values 
  filter(CellID %in% sites.transect) %>% 
  arrange(CellID)

# for camera data 
sites.camera <- unique(deploy_cam_visit_occu$CellID) # there is an NA, and I assume that because of the one location with lies outside the study area
occ.covs.camera <- envCovs %>%
  mutate(across(-CellID, ~ scale(.) %>% as.vector())) %>% # scale all numeric values 
filter(CellID %in% sites.camera) %>% 
  arrange(CellID)

occ.covs <- bind_rows(occ.covs.transect, occ.covs.camera)
occ.covs



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3 Presence-absence data as y ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract occu data as y 
y.transect <- locs_transects %>% 
  ungroup() %>% 
  st_drop_geometry() %>%
  select(CellID, Cell_visit, Occu) %>% 
  pivot_wider(names_from = Cell_visit, names_prefix = 'Visit_',
              values_from = Occu) %>% 
  select(-CellID) # keep CellID for troubleshooting

y.camera <- deploy_cam_visit_occu %>% ungroup() %>% 
  select(CellID, Cell_visit, Occu) %>% 
  pivot_wider(names_from = Cell_visit, names_prefix = 'Visit_',
              values_from = Occu) %>% 
  select(-CellID) # keep CellID for troubleshooting

y <- list(y.transect, y.camera)
y

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.4 Create site data to ensure proper linkage between data sets ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get sites of each data source
sites.transect.ind <- 1:length(sites.transect)
sites.camera.ind <- (max(sites.transect.ind)+1):(max(sites.transect.ind)+length(sites.camera))
sites <- list(sites.transect.ind, sites.camera.ind)
sites

# calc n.sites
n.sites.camera <- length(sites.camera)
n.sites.transect <- length(sites.transect)
n.sites <- n.sites.camera + n.sites.transect

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.5. Bring all data together ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check data 
sites
occ.covs
y
det.covs

# create data list which holds all data needed to run the model
data.list <- list(occ.covs, det.covs, y, sites)
names(data.list) <- c('occ.covs', 'det.covs', 'y', 'sites') # name data.list corresponding to the requirements if spOccupancy::intPGOcc

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Set parameters for integrated occupancy model using spOccupancy::intPGOcc ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set inits, alpha - det.covs, z for 
inits.list <- list(alpha = list(transect = rep(0, length(det.covs.transect)+1), # alpha gives initial values for det with each a list per data source, which hosts a vec with length of n() of predictors for this data source
                                camera = rep(0, 18)), # choose 0 as initial value
                   beta = rep(0, 1), # for ecological state model, start occupancy covariates at 0, length is the number of occu covs (n.occ.covs)
                   #sigma.sq.psi.inits, for random effects in the occurence model
                   sigma.sq.p.inits = list(transect = 0, # for random effects in the det model, wrong list name in the documentation
                                     camera = 0),  
                   z = rep(1, n.sites)) # z is for latent variable (here occupancy), start with 1 for all sites occupied

# set priors 
priors.list <- list(beta.normal = list(mean = 0, var = 2.72), # priors for beta, the ecological state model (occu) given in a list where two vectors are given, first for mean and second for variance, if they are all the same, only one value per tag
                    alpha.normal = list(mean = list(0, 0), 
                                        var = list(2.72, 2.72)),
                    sigma.sq.p.ig = list(shape = c(0.1), # random effects for det, vector in the lists have to have length of number of random effects 
                                       scale = c(0.1)))
n.samples <- 20000

# call simple first model 
gc()
m1 <- intPGOcc(occ.formula = ~ river_density_med_large + Distance_large_river + mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + JRC_transition_Undisturbed_tropical_moist_forest, #occ.cov, 
               det.formula = list(transect = ~ Julian_Date_Start_Transect + I(Julian_Date_Start_Transect^2)  + Transect_Length + Project_Transect + Season_Transect + (1 | Year_Transect), 
                                  camera = ~ Julian_Date_Start_Camera + I(Julian_Date_Start_Camera^2) + Trapping_Days + Project_Camera + Season_Camera + (1 | Year_Camera)), 
               data = data.list,
               inits = inits.list,
               n.samples = n.samples, 
               priors = priors.list, 
               n.omp.threads = 5, # use 5 cores
               verbose = TRUE, 
               n.report = 2500, 
               n.burn = 10000, 
               n.thin = 1, # no thinning
               n.chains = 3)


# access model 
summary(m1) # Rhat looks pretty good, ESS always (mostly far) above 1000
# plot(m1, param = 'alpha')


# produce posterior predictive checks, here using chi-square and freeman tukey - the results aren't very similar
ppc_m1_chi <- ppcOcc(m1, fit.stat = 'chi-squared', group = 1) # group 1 - groups values by site, group 2 - groups values per replicate
ppc_m1_ft <- ppcOcc(m1, fit.stat = 'freeman-tukey', group = 1)
summary(ppc_m1_chi) # bayesian p-values are 0.135 and 0.69
summary(ppc_m1_ft) # here bayesian p-values look much better, both ~0.27-0.28


# produce a model check plot, code taken from https://doserlab.com/files/spoccupancy-web/articles/modelfitting
# chi square stats
ppc_m1_chi_df <- data.frame(fit = ppc_m1_chi$fit.y[[2]], 
                            fit.rep = ppc_m1_chi$fit.y.rep[[2]], 
                            color = 'lightskyblue1')
ppc_m1_chi_df $color[ppc_m1_chi_df$fit.rep > ppc_m1_chi_df$fit] <- 'lightsalmon'
plot(ppc_m1_chi_df$fit, ppc_m1_chi_df$fit.rep, bg = ppc_m1_chi_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_m1_chi_df$fit, ppc_m1_chi_df$fit, col = 'black')
# freeman tukey stats 
ppc_m1_ft_df <- data.frame(fit = ppc_m1_ft$fit.y[[2]], 
                            fit.rep = ppc_m1_ft$fit.y.rep[[2]], 
                            color = 'lightskyblue1')
ppc_m1_ft_df$color[ppc_m1_ft_df$fit.rep > ppc_m1_ft_df$fit] <- 'lightsalmon'
plot(ppc_m1_ft_df$fit, ppc_m1_ft_df$fit.rep, bg = ppc_m1_ft_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc_m1_ft_df$fit, ppc_m1_ft_df$fit, col = 'black')

# lets have a look if there are any very influential data points 
diff_fit_chi <- ppc_m1_chi$fit.y.rep.group.quants[[1]][3, ] - ppc_m1_chi$fit.y.group.quants[[1]][3, ] # change list to 1 or 2, depending on transect or camera, respectively 
plot(diff_fit_chi, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')

diff_fit_ft <- ppc_m1_ft$fit.y.rep.group.quants[[2]][3, ] - ppc_m1_ft$fit.y.group.quants[[2]][3, ] # change list to 1 or 2, depending on transect or camera, respectively 
plot(diff_fit_ft, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Visualise effect sizes ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot effect sizes using MCMCvis
library(MCMCvis)
#jpeg(filename = 'output/plots/Effect_sizes_occu.jpg', height = 800, width = 800)
MCMCplot(m1$beta.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Occupancy Effect Sizes")
#dev.off()

#jpeg(filename = 'output/plots/Effect_sizes_det.jpg', height = 1200, width = 1000)
MCMCplot(m1$alpha.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Detection Effect Sizes")
#dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Plot influence of variables on Occupancy and Detection  ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

new_data <- data.frame(
  `(Intercept)` = rep(1, 1000),
  # Scaled variables
  river_density_med_large = scale(seq(min(occ.covs$river_density_med_large), max(occ.covs$river_density_med_large), length.out = 1000)),
  Distance_large_river = scale(seq(min(occ.covs$Distance_large_river), max(occ.covs$Distance_large_river), length.out = 1000)),
  mean_elev = scale(seq(min(occ.covs$mean_elev), max(occ.covs$mean_elev), length.out = 1000)),
  JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 = scale(seq(min(occ.covs$JRC_transition_Degraded_forest_short_duration_disturbance_after_2014), 
                                                                                   max(occ.covs$JRC_transition_Degraded_forest_short_duration_disturbance_after_2014), length.out = 1000)),
  JRC_transition_Undisturbed_tropical_moist_forest = scale(seq(min(occ.covs$JRC_transition_Undisturbed_tropical_moist_forest), 
                                                               max(occ.covs$JRC_transition_Undisturbed_tropical_moist_forest), length.out = 1000)),
  
  # Unscaled variables for later plotting
  river_density_med_large_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(river_density_med_large)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(river_density_med_large)), 
    length.out = 1000),
  Distance_large_river_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(Distance_large_river)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(Distance_large_river)), 
    length.out = 1000),
  mean_elev_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(mean_elev)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(mean_elev)), 
    length.out = 1000),
  JRC_transition_Degraded_forest_short_duration_disturbance_after_2014_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(JRC_transition_Degraded_forest_short_duration_disturbance_after_2014)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(JRC_transition_Degraded_forest_short_duration_disturbance_after_2014)), 
    length.out = 1000),
  JRC_transition_Undisturbed_tropical_moist_forest_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(JRC_transition_Undisturbed_tropical_moist_forest)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(JRC_transition_Undisturbed_tropical_moist_forest)), 
    length.out = 1000))

new_data <- data.frame(
  `(Intercept)` = rep(1, 1000),
  # Scaled variables
  river_density_med_large = mean(scale(seq(min(occ.covs$river_density_med_large), max(occ.covs$river_density_med_large), length.out = 1000))),
  Distance_large_river = mean(scale(seq(min(occ.covs$Distance_large_river), max(occ.covs$Distance_large_river), length.out = 1000))),
  mean_elev = scale(seq(min(occ.covs$mean_elev), max(occ.covs$mean_elev), length.out = 1000)),
  JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 = mean(scale(seq(min(occ.covs$JRC_transition_Degraded_forest_short_duration_disturbance_after_2014), 
                                                                                   max(occ.covs$JRC_transition_Degraded_forest_short_duration_disturbance_after_2014), length.out = 1000))),
  JRC_transition_Undisturbed_tropical_moist_forest = mean(scale(seq(min(occ.covs$JRC_transition_Undisturbed_tropical_moist_forest), 
                                                               max(occ.covs$JRC_transition_Undisturbed_tropical_moist_forest), length.out = 1000))), 
  mean_elev_unscaled = seq(
    min(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(mean_elev)), 
    max(envCovs_sf %>% filter(CellID %in% c(sites.camera, sites.transect)) %>% pull(mean_elev)), 
    length.out = 1000))

X.0 <- as.matrix(new_data[,1:6])

# predict on new data frame 
pred_m1_effect_occu <- predict(m1, type = 'occupancy', ignore.RE = T, X.0 = X.0)

# summarize prediction 
pred_m1_effect_occu <- apply(pred_m1_effect_occu$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))
# create a df for plotting 
plot_m1_effect_occu <- as.data.frame(t(pred_m1_effect_occu)) %>% 
  rename(pred_mean = `50%`, pred_CI_lower = `2.5%`, pred_CI_upper = `97.5%`) 
plot_m1_effect_occu <- cbind(plot_m1_effect_occu, new_data)


# produce quick plot
ggplot(plot_m1_effect_occu) +
  geom_ribbon(aes(x = mean_elev_unscaled, ymin = pred_CI_lower, ymax = pred_CI_upper), fill = "skyblue", alpha = 0.5) +
  geom_line(aes(x = mean_elev_unscaled, y = pred_mean), color = "blue", size = 1) +
  ggtitle("Effect of Elevation on Predicted Occupancy") +
  xlab("Mean Elevation in meter above NN") + 
  ylab("Predicted Occupancy") +
  theme_bw()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Predict Pygmy hippo Occurence in Gola ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create model matrix
X.0 <- model.matrix(~ river_density_med_large + Distance_large_river + 
                      mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + 
                      JRC_transition_Undisturbed_tropical_moist_forest, 
                    data = envCovs_sf %>% st_drop_geometry() %>% mutate(across(everything(),~ scale(.) %>% as.vector())))

# predict 
pred_m1 <- predict(m1, type = 'occupancy', ignore.RE = T, X.0 = X.0)

# summarise prediction 
psi.hat.quants <- apply(pred_m1$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))

# create a df for plotting 
plot_m1 <- as.data.frame(t(psi.hat.quants)) %>% 
  mutate(CellID = row_number()) %>% 
  rename(pred_mean = `50%`, pred_CI_lower = `2.5%`, pred_CI_upper = `97.5%`)

# create sf
plot_m1_sf <- envCovs_sf %>% left_join(plot_m1, join_by(CellID))

# plot prediction on map
tm_shape(plot_m1_sf) +
  tm_polygons(fill = 'pred_mean')


# fancier plot
library(ggspatial)
ggplot(data = plot_m1_sf) +
  annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_sf(aes(fill = pred_mean), color = NA, alpha = 0.5) +  # Use pred_mean for fill color
  scale_fill_viridis_c(option = "viridis", name = "Predicted Mean Occupancy") + # scale_fill_viridis_c(option = "magma"), or replace "magma" with "inferno", "plasma", "cividis", 
  theme_minimal() +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), # add frame
        axis.line = element_blank()) +
  labs(title = "Predicted Pygmy Hippo Occupancy Probability", 
       subtitle = "Based on an Integrated Occupancy Model fitted in spOccupancy", 
       x = 'Longitude', y = 'Latutude')
#ggsave(filename = 'output/plots/PH_hotspot_map_all_data.jpg', plot = hotspot_map, height = 6, width = 12)

