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
         Season_Transect= as.factor(if_else(month(DateTime_Start) %in% 5:10, 'Wet', 'Dry')), # wet season from May (05) to October (10), https://doi.org/10.51847/8Wz28ID8Mn
         Transect_Length = as.numeric(transect_length), 
         Project_Transect = as.factor(Project)) # %>% select(-transect_length) 

deploy_cam_visit_occu <- deploy_cam_visit_occu %>% ungroup()%>% 
  mutate(Project_Camera = as.factor(Project), 
         Julian_Date_Start_Camera = yday(Visit_start), 
         Season_Camera = as.factor(if_else(month(Visit_start) %in% 5:10, 'Wet', 'Dry')),
         Trapping_Days = as.numeric(Visit_length)) # %>% select(-Visit_length) 


# extract det.covs data
names(locs_transects)
det.variables.transect <- c('Project_Transect', 'Julian_Date_Start_Transect', 'Transect_Length', 'Season_Transect') # fill in all variables that are interesting
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
det.variables.camera <- c('Project_Camera', 'Julian_Date_Start_Camera', 'Trapping_Days', 'Season_Camera') # fill in all variables which could be interesting
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
occ.covs.transect <- envCovs %>% 
  filter(CellID %in% sites.transect) %>% 
  arrange(CellID)

# for camera data 
sites.camera <- unique(deploy_cam_visit_occu$CellID) # there is an NA, and I assume that because of the one location with lies outside the study area
occ.covs.camera <- envCovs %>% 
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

# alpha corresponds to 


# set inits, alpha - det.covs, z for 
inits.list <- list(alpha = list(transect = rep(0, length(det.covs.transect)+1), # alpha gives initial values for det with each a list per data source, which hosts a vec with length of n() of predictors for this data source
                                camera = rep(0, 17)), # choose 0 as initial value
                   beta = rep(0, 1), # for ecological state model, start occupancy covariates at 0, length is the number of occu covs (n.occ.covs)
                   #sigma.sq.psi, for random effects in the occurence model
                   #sigma.sq.p, # for random effects in the det model 
                   z = rep(1, n.sites)) # z is for latent variable (here occupancy), start with 1 for all sites occupied

# set priors 
priors.list <- list(beta.normal = list(mean = 0, var = 2.72), # priors for beta, the ecological state model (occu) given in a list where two vectors are given, first for mean and second for variance, if they are all the same, only one value per tag
                    alpha.normal = list(mean = list(0, 0), 
                                        var = list(2.72, 2.72)))
n.samples <- 70000

# call model 
out <- intPGOcc(occ.formula = ~ 1, #occ.cov, 
                det.formula = list(transect = ~ Julian_Date_Start_Transect + Transect_Length + Project_Transect + Season_Transect, 
                                   camera = ~ Julian_Date_Start_Camera + Trapping_Days + Project_Camera + Season_Camera), 
                data = data.list,
                inits = inits.list,
                n.samples = n.samples, 
                priors = priors.list, 
                n.omp.threads = 1, 
                verbose = TRUE, 
                n.report = 1000, 
                n.burn = 20000, 
                n.thin = 1, # no thinning
                n.chains = 3)




summary(out)









plot(locs_transects_visit_occu)
class(locs_transects_visit_occu)

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
head(pres_cam)
pres_cam <- pres_cam %>% mutate(Obs_Method = 'Camera Trap') %>% 
  select(Project, Country, Obs_DateTime, UTM_X_meters , UTM_Y_meters, SiteID, Obs_Method)

# merge data together 
pres_sf <- bind_rows(pres_cam, pres_opp,
                     pres_transects) %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629)


head(hb.dat)


######
## data for spOccupancy has to be in a 3-dim array with dim corresponding to species (1), site (row), and replicate
## site = row, replicated visits are columns, every value is either presence or absence 

head(pres_transects)
