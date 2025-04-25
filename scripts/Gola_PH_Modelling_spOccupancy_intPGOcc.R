#### Pygmy hippo single season integrated occupancy model in spOccupancy  
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025

# data has roughgly been prepared elsewhere


# check again how the variables are coded and with scaled and unscaled variables, scale only numeric 
# check priors and inits and set them to appropriate values 
# 


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

# make a plot to check if everything worked 
# tmap_mode("view")  # interactive mode
# tm_shape(grid_sf) +
#   tm_polygons(fill = tm_const()) + 
# tm_shape(locs_transects) + 
#   tm_lines(col = 'subtransect', lwd = 1.5, id = 'uniqueID') 

# create a sf object 
pres_transects_sf <- pres_transects %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629) %>% 
  st_transform(crs = 32629)

# link the presences to the subtransects, this is a bit ugly but the best way I could find to do this 
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
rm(joined_pres, df)

# make a plot to check, that everything went okay - there should be subtransects 
# tm_shape(grid_sf) +
#   tm_polygons() +
# tm_shape(locs_transects) + 
#   tm_lines(col='subtransect', lwd = 2, palette=c("blue", "darkgreen", "red", "lightblue", 'orange', 'purple')) +
# tm_shape(pres_transects_sf) +
#   tm_dots(col='subtransect', palette=c("blue", "darkgreen", "red", "lightblue", 'orange', 'purple'), size=0.5)

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

deploy_cam_visit_occu %>% head()
envCovs_sf %>% head()
locs_transects %>% head()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1 Det Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate a few det covariates, improve season by assigning the season where the majority was surveyed
locs_transects <- locs_transects %>% ungroup() %>%
  mutate(Date_Transect_unscaled = lubridate::yday(DateTime_Start), # Julian start date, improve by taking a mean date?
         Year_Transect_fact = as.factor(year(DateTime_Start)), # survey year coded as factor
         Year_Transect_num_unscaled = as.numeric(year(DateTime_Start)-2011), # survey year coded as numeric
         Season_Transect = as.factor(if_else(month(DateTime_Start) %in% 5:10, 'Wet', 'Dry')), # wet season from May (05) to October (10), https://doi.org/10.51847/8Wz28ID8Mn
         Transect_Length_unscaled = as.numeric(transect_length), 
         # Period_Transect = as.factor(if_else(year(DateTime_End) <= 2017, '2011-2017', '2018-2024')),
         Project_Transect_fact = as.factor(Project), 
         Project_Transect_num_unscaled = as.numeric(if_else(Project =='Pygmy Hippo ARTP_REDD 2013-2014', 0, 1))) # %>% select(-transect_length) 
locs_transects <- locs_transects %>% mutate(Year_Transect_num = as.numeric(scale(Year_Transect_num_unscaled)), # scale numeric variables
                                            Date_Transect = as.numeric(scale(Date_Transect_unscaled)), 
                                            Transect_Length = as.numeric(scale(Transect_Length_unscaled)), 
                                            Project_Transect_num = as.numeric(scale(Project_Transect_num_unscaled)))

deploy_cam_visit_occu <- deploy_cam_visit_occu %>% ungroup()%>% 
  mutate(Project_Camera = as.factor(Project), 
         Project_lumped_Camera = as.factor(case_when(Project %in% c('BasalZoo_2024','Basel Zoo PygmyHippo 2018-2020') ~ 'BaselZoo', # fill a more meaningful categorisation in
                                                     Project %in% c('REDD_SP1', 'REDD_SP2', 'REDD_SP3', 'REDD_SP4', 'Pygmy Hippo REDD CT 2019-2021','artp_p1', 'artp_p2', 'REDD_ARTP_PygmyHippo 2013-2014') ~ 'REDD_ARTP', 
                                                     Project %in% c('Darwin19_22', 'Darwin_morroRiver', 'darwin_13_17') ~ 'Darwin', 
                                                     Project %in% c('IWT_CF') ~ 'IWT')),
         Project_Focus_fact = as.ordered(if_else(Project %in% c('BasalZoo_2024', 'Basel Zoo PygmyHippo 2018-2020', 
                                                                'Darwin_morroRiver', 'IWT_CT', 'Pygmy Hippo REDD CT 2019-2021', 
                                                                'REDD_ARTP_PygmyHippo 2013-2014'), 'Pygmy_hippo', 'Other')),
         Project_Focus_num_unscaled = as.numeric(if_else(Project_Focus_fact == 'Pygmy_hippo', 1, 0)), # Pygmy hippo is coded as 1, other as 0
         Date_Camera_unscaled = yday(Visit_start), 
         Year_Camera_fact = as.factor(year(Visit_start)), # survey year coded as factor
         Year_Camera_num_unscaled = as.numeric(year(Visit_start)-2011), # survey year coded as numeric
         Season_Camera = as.factor(if_else(month(Visit_start) %in% 5:10, 'Wet', 'Dry')),
         # Period_Camera = as.factor(if_else(year(Visit_end) <= 2017, '2011-2017', '2018-2024')),
         Trapping_Days_unscaled = as.numeric(Visit_length)) # %>% select(-transect_length) # %>% select(-Visit_length) 
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% mutate(Year_Camera_num = as.numeric(scale(Year_Camera_num_unscaled)), # scale numeric variables
                                            Date_Camera = as.numeric(scale(Date_Camera_unscaled)), 
                                            Trapping_Days = as.numeric(scale(Trapping_Days_unscaled)), 
                                            Project_Focus_num = as.numeric(scale(Project_Focus_num_unscaled)))

# extract det.covs data
names(locs_transects) 
det.variables.transect <- c('Project_Transect_fact', 'Project_Transect_num', 'Project_Transect_num_unscaled', 'Date_Transect', 'Date_Transect_unscaled', 'Transect_Length',  
                            'Transect_Length_unscaled', 'Season_Transect', 'Year_Transect_num', 'Year_Transect_num_unscaled', 'Year_Transect_fact') # fill in all variables that are interesting
det.covs.transect <- list()
for(det.var in det.variables.transect){
  dat <- locs_transects %>% st_drop_geometry() %>%
    ungroup() %>% 
    select(CellID, Cell_visit, !!sym(det.var)) %>% 
    pivot_wider(names_from = Cell_visit, names_prefix = 'Visit_',
                values_from = !!sym(det.var))
  det.covs.transect[[det.var]] <- dat[,2:ncol(dat)] # dat # for troubleshooting use solely dat to keep CellID
  rm(dat)
}

names(deploy_cam_visit_occu)
det.variables.camera <- c('Project_Camera', 'Project_lumped_Camera', 'Project_Focus_num', 'Project_Focus_num_unscaled', 'Project_Focus_fact', 'Date_Camera', 'Date_Camera_unscaled', 'Trapping_Days', 
                          'Trapping_Days_unscaled', 'Season_Camera', 'Year_Camera_num', 'Year_Camera_num_unscaled', 'Year_Camera_fact') # fill in all variables which could be interesting
det.covs.camera <- list()
for(det.var in det.variables.camera){
  dat <- deploy_cam_visit_occu %>% ungroup() %>% 
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
occ.covs.transect.unscaled <- envCovs %>%
  filter(CellID %in% sites.transect) %>% 
  arrange(CellID)
occ.covs.transect <- occ.covs.transect.unscaled %>% 
  mutate(across(!c(CellID, Reserve_Type), ~ scale(.) %>% as.vector()), 
         Reserve_Type=as.factor(Reserve_Type)) # saved as factor, all other as numeric
str(occ.covs.transect)  

# for camera data 
sites.camera <- unique(deploy_cam_visit_occu$CellID) 
occ.covs.camera.unscaled <- envCovs %>%
  filter(CellID %in% sites.camera) %>% 
  arrange(CellID)
occ.covs.camera <- occ.covs.camera.unscaled %>% 
  mutate(across(!c(CellID, Reserve_Type), ~ scale(.) %>% as.vector()), # scale all numeric variables
         Reserve_Type=as.factor(Reserve_Type)) # saved as factor, all other as numeric
  
# bind everything together 
occ.covs <- bind_rows(occ.covs.transect, occ.covs.camera)
occ.covs.unscaled <- bind_rows(occ.covs.camera.unscaled, occ.covs.transect.unscaled)
occ.covs

rm(occ.covs.camera.unscaled, occ.covs.transect.unscaled) # remove unneeded objects 



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
inits.list <- list(alpha = list(transect = 0, # alpha gives initial values for det with each a list per data source, which hosts a vec with length of n() of predictors for this data source
                                camera = 0), # choose 0 as initial value
                   beta = 0, # for ecological state model, start occupancy covariates at 0, length is the number of occu covs (n.occ.covs)
                   # sigma.sq.psi.inits = list, # for random effects in the occurence model
                   sigma.sq.p.inits = list(transect = 0, # for random effects in the det model, wrong list name in the documentation
                                     camera = 0),  
                   z = rep(1, n.sites)) # z is for latent variable (here occupancy), start with 1 for all sites occupied

# set priors 
priors.list <- list(beta.normal = list(mean = 0, var = 2.72), # priors for beta, the ecological state model (occu) given in a list where two vectors are given, first for mean and second for variance, if they are all the same, only one value per tag
                    alpha.normal = list(mean = list(0, 0), 
                                        var = list(2.72, 2.72)),
                    sigma.sq.p.ig = list(shape = c(0.1), # random effects for det all follow inverse Gamma distribution, vector in the lists have to have length of number of random effects 
                                       scale = c(0.1)))
n.samples <- 20000

# call simple global model 

# different models: 
## year as.numeric and as.factor in both det in random effect 
## camera project lumped in det formula, camera project non lumped with all levels as random effects additional to year
## add reserve type into occ model 
## year is not doable in occ model 
## in det try both Date and Date^2
## add nested random effect for Season and year as factor in det models 

names(occ.covs) # variables for occ submodel
names(det.covs.transect) # variables for det transect submodel
names(det.covs.camera) # variables for det camera submodel 


gc()
m1 <- intPGOcc(occ.formula = ~  river_density_med_large + Distance_large_river + mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + JRC_transition_Undisturbed_tropical_moist_forest + EVI_2020 + Distance_settlement, #occ.cov, 
               det.formula = list(transect = ~ Date_Transect + I(Date_Transect^2)  + Transect_Length + Season_Transect + Project_Transect_fact + (1 | Year_Transect_num), 
                                  camera = ~ Date_Camera + I(Date_Camera^2) + Trapping_Days + Season_Camera + Project_Camera + (1 | Year_Camera_num)), 
               data = data.list, inits = inits.list, n.samples = n.samples, priors = priors.list, 
               n.omp.threads = 5, # use 5 cores
               verbose = TRUE, # means that messages about processing and computation are printed to console
               n.report = 2500, n.burn = 10000, 
               n.thin = 1, n.chains = 3)

#m2 <- intPGOcc(occ.formula = ~ river_density_med_large + Distance_large_river + mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + JRC_transition_Undisturbed_tropical_moist_forest , #occ.cov, 
#               det.formula = list(transect = ~ Date_Transect + I(Date_Transect^2)  + Transect_Length + Project_Transect + Season_Transect + (1 | Year_Transect), 
#                                  camera = ~ Date_Camera + I(Date_Camera^2) + Trapping_Days + Season_Camera + (1 | Year_Camera)), 
#               data = data.list,
#               inits = inits.list,
#               n.samples = n.samples, 
#               priors = priors.list, 
#               n.omp.threads = 5, # use 5 cores
#               verbose = TRUE, # means that messages about processing and computation are printed to console
#               n.report = 2500, 
#               n.burn = 10000, 
#               n.thin = 1, # no thinning
#               n.chains = 3)


# access model 
summary(m1) # Rhat looks pretty good, ESS always (mostly far) above 1000
# plot(m1, param = 'alpha')


# produce posterior predictive checks, here using chi-square and freeman tukey - the results aren't very similar
gc()
ppc_m1 <- ppcOcc(m1, fit.stat = 'freeman-tukey', group = 1) # group 1 - groups values by site, group 2 - groups values per replicate, there is also chi squared available
summary(ppc_m1) # get bayes p-value


# produce a model check plot, code taken from https://doserlab.com/files/spoccupancy-web/articles/modelfitting
ppc_m1_df <- data.frame(fit = ppc_m1$fit.y[[2]], 
                            fit.rep = ppc_m1$fit.y.rep[[2]], 
                            color = 'lightskyblue1')
ppc_m1_df$color[ppc_m1_df$fit.rep > ppc_m1_df$fit] <- 'lightsalmon'
plot(ppc_m1_df$fit, ppc_m1_df$fit.rep, bg = ppc_m1_df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True', main = 'Camera Data bay p-val 0.2071')
lines(ppc_m1_df$fit, ppc_m1_df$fit, col = 'black')

# lets have a look if there are any very influential data points 
diff_fit <- ppc_m1$fit.y.rep.group.quants[[2]][3, ] - ppc_m1$fit.y.group.quants[[2]][3, ] # change list to 1 or 2, depending on transect or camera, respectively 
plot(diff_fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy', main = 'Camera Data')


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
##### 8. Plot marginal effects plots influence for Occupancy   ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# predict for each variable while holding the other variables constant (here mean of scaled numeric variables, which would be 0)
n <- 1000 # length of df determines how many predictions should be made, influence on computation time 
const_val <- data.frame(`(Intercept)` = rep(1, n),
                       river_density_med_large = rep(mean(occ.covs$river_density_med_large), length.out = n), 
                       Distance_large_river = rep(mean(occ.covs$Distance_large_river), length.out = n), 
                       mean_elev = rep(mean(occ.covs$mean_elev), length.out = n), 
                       JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 = rep(mean(occ.covs$JRC_transition_Degraded_forest_short_duration_disturbance_after_2014), length.out = n), 
                       JRC_transition_Undisturbed_tropical_moist_forest = rep(mean(occ.covs$JRC_transition_Undisturbed_tropical_moist_forest), length.out = n), 
                       EVI_2020 = rep(mean(occ.covs$EVI_2020), length.out = n), 
                       Distance_settlement = rep(mean(occ.covs$Distance_settlement), length.out = n))
effect_occu_all <- data.frame() # initialise df
occ.var <- colnames(m1$X)[2:dim(m1$X)[2]] # get all variables used in occu fomula, except of intercept
for(var in occ.var){
  X.0 <- const_val %>% # create matrix to hand over to predict function
    mutate(!!var := seq(min(occ.covs[,var]), max(occ.covs[,var]), length.out = n),
           unscaled = seq(min(occ.covs.unscaled[,var]), max(occ.covs.unscaled[,var]), length.out = n)) %>% 
    as.matrix()
  
  pred_occu_effect <- predict(m1, type = 'occupancy', ignore.RE = T, X.0 = X.0[,!grepl('unscaled', colnames(X.0))]) # predict with matrix
  pred_occu_effect <- apply(pred_occu_effect$psi.0.samples, 2, quantile, c(0.025, 0.5, 0.975))   # summarise prediction 
  
  occu_effect <- as.data.frame(t(pred_occu_effect)) %>%   # wrangle data together 
    rename(pred_mean = `50%`, pred_CI_lower = `2.5%`, pred_CI_upper = `97.5%`) 
  occu_effect <- bind_cols(occu_effect,  unscaled = X.0[,'unscaled']) %>% 
    mutate(Variable = var)
  
  effect_occu_all <- bind_rows(effect_occu_all, occu_effect)   # save in massive df
}

# produce quick plot
ggplot(effect_occu_all) +
  geom_ribbon(aes(x = unscaled, ymin = pred_CI_lower, ymax = pred_CI_upper), fill = "skyblue", alpha = 0.5) +
  geom_line(aes(x = unscaled, y = pred_mean), color = "blue", size = 1) +
  facet_wrap(~Variable, scale = 'free_x')+
  ggtitle("Effect of Model Variables on Predicted Occupancy") +
  xlab("Variable") + 
  ylab("Predicted Occupancy") +
  theme_bw()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Predict Pygmy hippo Occurence in Gola ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create model matrix
X.0 <- model.matrix(~ river_density_med_large + Distance_large_river + 
                      mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + 
                      JRC_transition_Undisturbed_tropical_moist_forest + EVI_2020 + Distance_settlement, 
                    data = envCovs_sf %>% select(-Reserve_Type) %>% st_drop_geometry() %>% mutate(across(everything(),~ scale(.) %>% as.vector())))

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
  #annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_sf(aes(fill = pred_mean), # color = NA, 
          alpha = 1) +  # Use pred_mean for fill color
  scale_fill_viridis_c(option = "plasma", name = "Mean Predicted Occupancy", limits = c(0,1), direction = -1) + # scale_fill_viridis_c(option = "magma"), or replace "magma" with "inferno", "plasma", "cividis", 
  theme_bw() +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), # add frame
        axis.line = element_blank()) +
  labs(title = "Predicted Pygmy Hippo Occupancy Probability", 
       subtitle = "Based on an Integrated Occupancy Model fitted in spOccupancy", 
       x = 'Longitude', y = 'Latutude')
#ggsave(filename = 'output/plots/PH_hotspot_map_all_data.jpg', plot = hotspot_map, height = 6, width = 12)

