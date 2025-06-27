#### Pygmy hippo static multi-season integrated occupancy model in spOccupancy
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in April 2025, revised in May


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
library(scales) # only for plot
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. load data which has been prepared elsewhere #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# presence and deploy/survey data 
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
##### 2. Determine an appropriate Date for splitting the data into 2 time periods #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# camera data 
deploy_cam %>% 
  ggplot() +
  geom_segment(data = deploy_cam %>% filter(year(Deployment ) %in% c(2017, 2018)), 
               aes(x = Deployment, xend = Collection, y = SiteID, yend = SiteID), linewidth = 2, color = "blue") +
  geom_vline(xintercept = as.Date('2017-12-01'), linetype = 'dotted', color = 'red', linewidth = 1.5)


# transect data 
locs_transects %>% 
  ggplot() +
  geom_point(aes(x = DateTime_Start, y = uniqueID), size = 1.5, color = 'forestgreen') +
  geom_point(aes(x = DateTime_End, y = uniqueID), size = 1.5, color = 'darkred') +
  geom_vline(xintercept = as.POSIXct('2017-12-01'), linetype = 'dotted', color = 'red', linewidth = 1.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Prepare presence camera trap data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create an index for the period - this would have to be improved for more automation and more survey periods
deploy_cam <- deploy_cam %>% mutate(Period = if_else(year(Deployment) <= 2017, 1,2))

# create list for results 
deploy_cam_visits_periods <- list()

for(p in 1:max(deploy_cam$Period)){
  visits <- deploy_cam %>% filter(Period == p) %>%
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
    group_by(CellID) %>% # CellID instead of SiteID - this is the step, where we decide to define a visit across all projects an
    mutate(Cell_visit = row_number())
  
  deploy_cam_visits_periods[[p]] <- visits
}

# bind df of different periods together and remove unneeded objects 
deploy_cam_visit <- bind_rows(deploy_cam_visits_periods)
rm(visits, deploy_cam_visits_periods)


# transfer presences to deployment data, this assumes, that the observations lie within the deployment period
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

library(stplanr) # package needed for line_segment function 

# create survey periods for transect surveys 
locs_transects <- locs_transects %>% 
  mutate(Period = if_else(year(DateTime_Start) <= 2017, 1, 2))

# split transects into subtransects which are interpreted as visits 
locs_transects <- locs_transects %>% 
  st_intersection(grid_sf) %>%
  st_cast(to = 'MULTILINESTRING') %>% # first transform everything to MULTILINESTRING
  st_cast(to ='LINESTRING') %>% # transfer everything back to linestring
  group_by(Period, uniqueID, CellID) %>% 
  # mutate(transect_length = st_length(geometry, which = 'Euclidean')) %>% 
  line_segment(segment_length = 700) %>% # segment length in the crs unit, which is meters here 
  # group_by(Period, uniqueID, CellID) %>% 
  mutate(subtransect = as.factor(row_number())) %>% 
  group_by(Period, uniqueID, CellID, subtransect) %>% 
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
  mutate(Period = if_else(year(Obs_DateTime) <= 2017, 1, 2)) %>%
  group_by(Period, CellID, uniqueID, subtransect) %>% 
  summarise(Presences = n(), .groups = 'drop') # this produces a count of observations per subtransect in each gridCell, but if distinct (not just multiple signs in a subtransect) detections of pygmy hippos are needed, use n_distinct() here 
  # st_cast(to = 'MULTIPOINT')

locs_transects <- locs_transects %>% 
  left_join(occu_transects %>% st_drop_geometry(), join_by(Period, CellID, uniqueID, subtransect)) %>% 
  mutate(Presences = if_else(is.na(Presences), 0, Presences),
         Occu = if_else(Presences >= 1, 1, 0)) %>% 
  group_by(Period, CellID) %>% 
  mutate(Cell_visit = row_number())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Bring data together for spOccupancy ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# all nessecary data is organised in these data frames
deploy_cam_visit_occu %>% head()
envCovs_sf %>% head()
locs_transects %>% head()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1 Observation level - Det Covs and Occu data  ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# here, I prepare observation level data for a) detection covariates and b) Occupancy data itself since the data structures are identical

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1.1 Transect Surveys Data ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate a few det covariates, improve season by assigning the season where the majority was surveyed
locs_transects <- locs_transects %>% ungroup() %>%
  mutate(Date_Transect_unscaled = lubridate::yday(DateTime_Start), # Julian start date, improve by taking a mean date?
         Year_Transect_fact = as.factor(year(DateTime_Start)), # survey year coded as factor
         Year_Transect_num_unscaled = as.numeric(year(DateTime_Start)-2011), # survey year coded as numeric
         Season_Transect_fact = as.factor(if_else(month(DateTime_Start) %in% 5:10, 'Wet', 'Dry')), # wet season from May (05) to October (10), dry season all other months  https://doi.org/10.51847/8Wz28ID8Mn
         Season_Transect_num_unscaled = as.numeric(if_else(month(DateTime_Start) %in% 5:10, 1, 0)), # wet season coded as 1, dry as 0
         Transect_Length_unscaled = as.numeric(transect_length), 
         Period_Transect_fact = as.factor(if_else(year(DateTime_Start) <= 2017, '2011-2017', '2018-2025')),
         Period_Transect_num_unscaled = as.numeric(if_else(year(DateTime_Start) <= 2017, 1, 2)), # 1 is period from 2011 to 2017
         Project_Transect_fact = as.factor(Project), 
         Project_Transect_num_unscaled = as.numeric(if_else(Project =='Pygmy Hippo ARTP_REDD 2013-2014', 0, 1))) # %>% select(-transect_length) 
locs_transects <- locs_transects %>% mutate(Year_Transect_num = as.numeric(scale(Year_Transect_num_unscaled)), # scale numeric variables
                                            Date_Transect = as.numeric(scale(Date_Transect_unscaled)), 
                                            Transect_Length = as.numeric(scale(Transect_Length_unscaled)), 
                                            Project_Transect_num = as.numeric(scale(Project_Transect_num_unscaled)),
                                            Period_Transect_num = as.numeric(scale(Project_Transect_num_unscaled)),
                                            Season_Transect_num = as.numeric(scale(Season_Transect_num_unscaled)))


# extract det.covs data - observation.covariates
periods <- 1:max(locs_transects$Period)
visits <- 1:max(locs_transects$Cell_visit)
cells <- sort(unique(locs_transects$CellID))
obs.covs.selection <- c('Date_Transect', 'Transect_Length', 'Project_Transect_num', 'Project_Transect_fact', 'Season_Transect_fact', 'Season_Transect_num','Year_Transect_num', 'Year_Transect_fact', 'Period_Transect_fact', 'Period_Transect_num','Occu') 
obs.covs.transects <- replicate(length(obs.covs.selection),
                          array(NA, dim = c(length(cells), max(locs_transects$Period), max(locs_transects$Cell_visit))), 
                          simplify = F)
names(obs.covs.transects) <- obs.covs.selection

for(obs.det in obs.covs.selection){
  for(p in periods){
    for(v in visits){
      for(c in 1:length(cells)){
        value <- locs_transects %>%
          filter(CellID == cells[c], Period == p, Cell_visit == v) %>%
          pull(obs.det)
        if(length(value) == 0){
          obs.covs.transects[[obs.det]][c,p,v] <- NA
        } else {
          obs.covs.transects[[obs.det]][c,p,v] <- value
          }
      }
    }
  }
  dimnames(obs.covs.transects[[obs.det]]) <- list(cells, paste0("Period_", periods), paste0("Visit_", visits)) # name the array
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1.2 Camera Trapping Data ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate covariates 
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% ungroup()%>% 
  mutate(Project_Camera = as.factor(Project), 
         Project_Focus_fact = as.ordered(if_else(Project %in% c('Pygmy Hippo REDD CT 2019-2021', 'REDD_ARTP_PygmyHippo 2013-2014', 
                                                                'Darwin_morroRiver', 'Basel Zoo PygmyHippo 2018-2020', 'BasalZoo_2024', 'CYCV camera trap data 2020'), 
                                                 'Pygmy_hippo', 'Other')),
         Project_Focus_num_unscaled = as.numeric(if_else(Project_Focus_fact == 'Pygmy_hippo', 1, 0)), # Pygmy hippo is coded as 1, other as 0
         Date_Camera_unscaled = yday(Visit_start), 
         Year_Camera_fact = as.factor(year(Visit_start)), # survey year coded as factor
         Year_Camera_num_unscaled = as.numeric(year(Visit_start)-2011), # survey year coded as numeric
         Season_Camera_fact = as.factor(if_else(month(Visit_start) %in% 5:10, 'Wet', 'Dry')),
         Season_Camera_num_unscaled = as.numeric(if_else(month(Visit_start) %in% 5:10, 1, 0)), # wet season coded as 1, dry as 0
         Period_Camera_fact = as.factor(if_else(year(Visit_end) <= 2017, '2011-2017', '2018-2025')),
         Period_Camera_num_unscaled = as.numeric(if_else(year(Visit_end)<= 2017, 1, 2)), # 1 is period 1 from 2011 to 2017
         Trapping_Days_unscaled = as.numeric(Visit_length)) 
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% mutate(Year_Camera_num = as.numeric(scale(Year_Camera_num_unscaled)), # scale numeric variables
                                                          Date_Camera = as.numeric(scale(Date_Camera_unscaled)), 
                                                          Trapping_Days = as.numeric(scale(Trapping_Days_unscaled)), 
                                                          Project_Focus_num = as.numeric(scale(Project_Focus_num_unscaled)), 
                                                          Period_Camera_num = as.numeric(scale(Period_Camera_num_unscaled)),
                                                          Season_Camera_num = as.numeric(scale(Season_Camera_num_unscaled)))


# extract det.covs on observation level - this takes some time due to 53 visits!
periods <- 1:max(deploy_cam_visit_occu$Period)
visits <- 1:max(deploy_cam_visit_occu$Cell_visit)
cells <- sort(unique(deploy_cam_visit_occu$CellID))
obs.covs.selection <- c('Date_Camera', 'Trapping_Days', 'Project_Camera' ,'Project_Focus_fact', 'Project_Focus_num', 'Season_Camera_fact', 'Season_Camera_num', 'Period_Camera_fact','Period_Camera_num', 'Year_Camera_num', 'Year_Camera_fact', 'Occu')
obs.covs.camera <- replicate(length(obs.covs.selection),
                          array(NA, dim = c(length(cells), max(periods), max(visits))), 
                          simplify = F)
names(obs.covs.camera) <- obs.covs.selection

for(obs.det in obs.covs.selection){ # this is very slow and it would be great to find a quicker solution!
  for(p in periods){
    for(v in visits){
      for(c in 1:length(cells)){
        value <- deploy_cam_visit_occu %>%
          filter(CellID == cells[c], Period == p, Cell_visit == v) %>%
          pull(obs.det)
        if(length(value) == 0){
          obs.covs.camera[[obs.det]][c,p,v] <- NA
        } else {
          obs.covs.camera[[obs.det]][c,p,v] <- value
        }
      }
    }
  }
  dimnames(obs.covs.camera[[obs.det]]) <- list(cells, paste0("Period_", periods), paste0("Visit_", visits)) # name the array
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1.3 Bring det.covs from both data sources together ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


det.covs <- list(Transect = obs.covs.transects[names(obs.covs.transects) != "Occu"], 
                 Camera = obs.covs.camera[names(obs.covs.camera) != "Occu"])



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.2 Presence-absence data as y ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract occupancy data from observation level lists and safe in new list
y <- list(obs.covs.transects["Occu"]$Occu, obs.covs.camera["Occu"]$Occu)
names(y) <- c('Transect', 'Camera')
y %>% head()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3 Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3.1 Site-Level Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get sites (cells) for each data source
sites.transect <- sort(unique(locs_transects$CellID))
sites.camera <- sort(unique(deploy_cam_visit_occu$CellID))

# df with all envCovs for all sites (across transect and camera data)
envCovs_all_unscaled <- bind_rows(envCovs %>% filter(CellID %in% sites.transect) %>% arrange(CellID) %>% mutate(Type = 'Transect'), 
          envCovs %>% filter(CellID %in% sites.camera) %>% arrange(CellID) %>% mutate(Type = 'Camera'))

# scale numeric variables  
envCovs_all <- envCovs_all_unscaled %>%
  mutate(across(-c(CellID, Reserve_Type, Type), ~ (. - mean(.)) / sd(.)), # manual scaling (otherwise use: scale(.) %>% as.vector())
         Reserve_Type=as.factor(Reserve_Type)) 

# bring all (scaled and unscaled)
site.covs <- envCovs_all_unscaled %>% 
  select(mean_elev, mean_slope, Distance_large_river, Distance_road, river_density_med_large) 
site.covs <- bind_cols(site.covs %>% rename_with(~paste0(., '_unscaled')), 
                       envCovs_all %>% select(all_of(c(names(site.covs), 'Reserve_Type')))) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3.2 Season-Level/Yearly-Site Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# keep season level covariates unscaled here and scale later on covariate level, not within each year of the covariate seperately

# prepare loop
periods <- 1:max(deploy_cam_visit_occu$Period) # number of periods 
cells <- c(sites.transect, sites.camera) # site/cell indices 
year_preds <- names(envCovs_sf)[grep("JRC_ann_changes|NDVI|EVI", names(envCovs_sf))] # select all variables from annual changes data set and EVI/NDVI
year_preds <- year_preds[grep('2013|2021', year_preds)] # filter out intermediate year 2017
# year_preds <- sub(pattern = 'JRC_ann_changes_', replacement = '',year_preds)
year.site.covs <- unique(gsub('_2013|_2021','', year_preds)) # fill in the years data should be extracted for (which refer to the periods)
year.site.covs.list <- list()


# call loop
for(yr in year.site.covs){
  m <- matrix(data = NA, nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods))) 
  for(p in periods){
      m[,p] <- envCovs_all_unscaled[[ year_preds[grep(yr, year_preds)][p] ]] # take the unscaled df with envCovs
  }
  m_scaled <- (as.vector(m) - mean(as.vector(m)))/sd(as.vector(m)) # scale all variables across years manually (otherwise use scale(as.vector(m)))
  m_scaled <- matrix(m_scaled, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))  
  
  year.site.covs.list[[yr]] <- m_scaled
  year.site.covs.list[[paste0(yr, '_unscaled')]] <- m
}

# rename variables 
names(year.site.covs.list) <- sub(pattern = 'JRC_ann_changes_', replacement = '',names(year.site.covs.list))

# manually add year or year-period as yearly-site covariate - two versions, numeric and scaled and as factor
year.site.covs.list[['Period_num']] <- matrix(data = as.numeric(scale(rep(periods, each = length(cells)))), nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods)))
year.site.covs.list[['Period_fact']] <- matrix(data = rep(c('2011-2017','2018-2025'), each = length(cells)), nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods)))


# create occ.covs list which compromises all data 
occ.covs <- c(as.list(site.covs), year.site.covs.list) # add all matrices from a list to the occ.covs list by using c()
occ.covs %>% str()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.4 Create site data to ensure proper linkage between data sets ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a df from occ.covs which includes information on data source/type and row.numbers
sites.idx <- envCovs_all %>% mutate(row.idx = row_number())

sites <- list(Transect = sites.idx$row.idx[sites.idx$Type == 'Transect'], 
              Camera = sites.idx$row.idx[sites.idx$Type == 'Camera'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.4 Create season data to ensure proper linkage between primary survey periods/years ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seasons <- list(Transect = sort(unique(locs_transects$Period)), 
                Camera = sort(unique(deploy_cam_visit_occu$Period)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.5. Bring all data together ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# check data 
sites
seasons
occ.covs %>% str()
y
det.covs

# create data list which holds all data needed to run the model
data.list <- list(occ.covs, det.covs,y,  sites,  seasons)
names(data.list) <- c('occ.covs', 'det.covs', 'y', 'sites', 'seasons') # name data.list corresponding to the requirements if spOccupancy::intPGOcc
str(data.list)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Set parameters and call integrated occupancy model using spOccupancy::tIntPGOcc ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set inits, alpha - det.covs, z for 
inits.list = list(z = matrix(data = 1, ncol = max(seasons$Transect), nrow = max(sites$Camera)), # z is for latent variable (here occupancy), start with 1 for all sites occupied
                  beta = 0, # for ecological state model, start occupancy covariates at 0, length is the number of occu covs (n.occ.covs)
                  alpha = list(Transect = 0, # for each data source a list with initial values for det model
                               Camera = 0), 
                  sigma.sq.psi = NULL, # for random effects in occurrence model as list
                  sigma.sq.p = c(0), # for random effects detection model as c
                  sigma.sq.t = NULL, # only relevant if ar1 = T
                  rho = NULL) # only relevant if ar1 = T

priors.list <- list(beta.normal = list(mean = 0, var = 2.72), # priors for beta, the ecological state model (occu) given in a list where two vectors are given, first for mean and second for variance, if they are all the same, only one value per tag
                    alpha.normal = list(mean = list(0, 0), 
                                        var = list(2.72, 2.72)),
                    sigma.sq.p.ig = list(shape = c(0.1), # random effects for det all follow inverse Gamma distribution, vector in the lists have to have length of number of random effects 
                                         scale = c(0.1)), 
                    sigma.sq.psi.ig = NULL, 
                    sigma.sq.t.ig = NULL, 
                    rho.unif = NULL)

# set model parameters 
n.report <- 2
batch.length <- 1500
n.batch <- 10
n.thin <- 3
n.chains <- 4

# call global model 
gm <- tIntPGOcc(occ.formula = ~  river_density_med_large + Distance_large_river + Permanent_and_seasonal_water 
                + Undisturbed_tropical_moist_forest 
                + mean_elev + mean_slope
                + Period_fact, 
                det.formula = list(Transect = ~Transect_Length + I(Date_Transect^2)  + as.factor(Season_Transect_fact) + as.factor(Period_Transect_fact), 
                                   Camera = ~Trapping_Days + I(Date_Camera^2) + as.factor(Project_Camera) + as.factor(Season_Camera_fact)  + (1| Year_Camera_num)), 
                data = data.list,
                n.report = n.report,
                n.batch = n.batch,
                batch.length = batch.length,
                n.thin = n.thin, 
                n.chains = n.chains)

# call model summary 
summary(gm)

# check MCMC chains 
plot(gm, param = 'alpha', density = F) # looks good
plot(gm, param = 'beta', density = F)

# check effect sizes 
library(MCMCvis)
MCMCplot(gm$beta.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Occupancy Effect Sizes")
MCMCplot(gm$alpha.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Detection Effect Sizes")

# Goodness-of-fit test, freeman-tukey and grouped by sites, described in https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml
ppc_gm1 <- ppcOcc(gm, fit.stat = 'freeman-tukey', group = 1) # group 1 - groups values by site, group 2 - groups values per replicate, there is also chi squared available
summary(ppc_gm1) # get bayes. p-value

# visualise Goodness-of-fit 

# produce a model check plot, code example taken from https://doserlab.com/files/spoccupancy-web/articles/modelfitting
ppc_result <- data.frame(fit = numeric(),fit.rep = numeric(),season = character(),dataset = character(),color = character())
for(d in 1:length(gm$det.formula)){
  for(s in 1:length(gm$seasons[[d]])){
    ppc_frame <- data.frame(fit = ppc_gm1$fit.y[[d]][,s], fit.rep = ppc_gm1$fit.y.rep[[d]][,s], 
                            season = s, dataset = names(gm$det.formula)[d], color = 'lightskyblue1')
    ppc_frame$color[ppc_frame$fit.rep > ppc_frame$fit] <- 'lightsalmon'
    ppc_result <- rbind(ppc_result, ppc_frame)
  }
}

# plot true vs fitted values 
ppc_result %>% mutate(season = if_else(season == 1, '2011-2017', '2018-2025')) %>% 
  ggplot() +
  geom_point(mapping = aes(x = fit, y = fit.rep, colour = color), size = 0.7, alpha = 0.25, show.legend = F) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1.2) +
  facet_grid(season~dataset, scales = 'free') +
  labs(x = 'True', y = 'Fit', title = 'True vs. Fitted Values for all Data Sources and Time Periods') +
  theme_bw()


# plot influential data points 
fit_result <- data.frame(fit = numeric(), numeric(),season = character(),dataset = character())
for(d in 1:length(gm$det.formula)){
  for(s in 1:length(gm$seasons[[d]])){
    fit_frame <- data.frame(fit = ppc_gm$fit.y.rep.group.quants[[d]][3,,][,s] - ppc_gm$fit.y.group.quants[[d]][3,,][,s], 
                            season = s, dataset = names(gm$det.formula)[d], SiteID = sites[[d]])
    fit_result <- rbind(fit_result, fit_frame)
  }
}

fit_result %>% mutate(season = if_else(season == 1, '2011-2017', '2017-2025')) %>% 
  ggplot() +
  geom_point(mapping = aes(y = fit, x = SiteID)) +
  facet_grid(season~dataset, scale = 'free') + 
  labs(y = 'Fit', title = 'Replicate - True Discrepancy') + 
  theme_bw()

# remove saved objects to avoid memory issues 
rm(gm, ppc_gm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Model selection using wAIC ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set model parameters
n.report <- 2
batch.length <- 1500
n.batch <- 10
n.thin <- 3
n.chains <- 4
# occ.formula <- ~  river_density_med_large + Distance_large_river + Permanent_and_seasonal_water + Undisturbed_tropical_moist_forest + mean_elev + mean_slope + Period_fact
# det.formula <- list(Transect = ~Transect_Length + I(Date_Transect^2)  + as.factor(Season_Transect_fact) + as.factor(Period_Transect_fact), 
#                     Camera = ~Trapping_Days + I(Date_Camera^2) + as.factor(Project_Camera) + as.factor(Season_Camera_fact) + (1| Year_Camera_num))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7.1 Detection submodels ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# only implement two steps here to decide which det models to choose - separately fit models with Date or Date^2 and (1| Year_Transect_num) or as.factor(Period_Transect_fact)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7.1.1 Transect data detection submodels   ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# for date to decide on date or date^2

date.formula <- list(~Transect_Length+as.factor(Season_Transect_fact)+as.factor(Period_Transect_fact)+Date_Transect, # build formulas 
                     ~Transect_Length+as.factor(Season_Transect_fact)+as.factor(Period_Transect_fact)+I(Date_Transect^2))
model.comparison <- data.frame(occ = character(),det.transect = character(),det.camera = character(), # crate input df
                               wAIC = numeric(),wAIC.transect=numeric(),wAIC.camera=numeric(), step = character()) 
for(d in 1:length(date.formula)){
  # fit model
  model <- tIntPGOcc(occ.formula = ~1, det.formula = list(Transect = date.formula[[d]], Camera = ~1), 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wAIC
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'det.transect.date')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

# get best Date model and use this single best model (lowest wAIC) to build formula of next models 
best.date.transect <- model.comparison %>% filter(step == 'det.transect.date') %>% slice_min(wAIC.transect) %>% pull(det.transect) # select best submodel based on wAIC for this particular data source



# decide between period or (1|Year)

# year.formula <- list(as.formula(paste(best.date.transect, '+ as.factor(Period_Transect_fact)')), 
#                      as.formula(paste(best.date.transect, '+ (1| Year_Transect_num)')))

# for(f in 1:length(year.formula)){
#   # fit model
#   model <- tIntPGOcc(occ.formula = ~1, det.formula = list(Transect = year.formula[[f]], Camera = ~1), 
#                      data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
#   # save all values to model comparison table 
#   output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
#                        det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
#                        det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
#                        wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
#                        wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
#                        wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wA
#                         max.Rhat = max(unlist(model$rhat)),
#                         min.ESS = min(unlist(model$ESS)),
#                        step = 'det.transect.year')
#   model.comparison <- rbind(model.comparison, c(output))
#   
#   rm(model) # remove model object to clear memory
# }

# get best year/period model and use this single best model (lowest wAIC) to build formula of next models 
# best.year.transect <- as.formula(model.comparison %>% filter(step == 'det.transect.year') %>% slice_min(wAIC.transect) %>% pull(det.transect)) # select best submodel based on wAIC for this particular data source


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7.1.2 Camera data detection submodels   ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# decide for date or date^2
best.date.transect <- as.formula(best.date.transect)
date.formula <- list(~Trapping_Days+as.factor(Season_Camera_fact)+Date_Camera, # build formulas 
                     ~Trapping_Days+as.factor(Season_Camera_fact)+I(Date_Camera^2))


for(d in 1:length(date.formula)){
  # fit model
  model <- tIntPGOcc(occ.formula = ~1, det.formula = list(Transect = best.date.transect, Camera = date.formula[[d]]), 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wAIC
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'det.camera.date')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

# get best Date model and use this single best model (lowest wAIC) to build formula of next models 
best.date.camera <- model.comparison %>% filter(step == 'det.camera.date') %>% slice_min(wAIC.camera) %>% pull(det.camera)


# decide if project or project focus is better

project.formula <- list(as.formula(paste(best.date.camera, '+ Project_Camera')), 
                        as.formula(paste(best.date.camera, '+ as.factor(Project_Focus_fact)')))

for(d in 1:length(project.formula)){
  # fit model
  model <- tIntPGOcc(occ.formula = ~1, det.formula = list(Transect = best.date.transect, Camera = project.formula[[d]]), 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wA
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'det.camera.project')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

# get best Date model and use this single best model (lowest wAIC) to build formula of next models 
best.period.camera <- model.comparison %>% filter(step == 'det.camera.project') %>% slice_min(wAIC.camera) %>% pull(det.camera)


# decide if period or (1|Year) is better

year.formula <- list(as.formula(paste(best.period.camera, '+ as.factor(Period_Camera_fact)')), 
                     as.formula(paste(best.period.camera, '+ (1| Year_Camera_num)')))

for(d in 1:length(date.formula)){
  # fit model
  model <- tIntPGOcc(occ.formula = ~1, det.formula = list(Transect = best.date.transect, Camera = year.formula[[d]]), 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wA
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'det.camera.year')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

# get best year formula 
best.year.camera <- model.comparison %>% filter(step == 'det.camera.year') %>% slice_min(wAIC.camera) %>% pull(det.camera)



# bring best formulas together 
det.formula <- list(Transect = as.formula(best.date.transect),
                          Camera = as.formula(best.year.camera))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7.2 Occupancy submodel ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# tackle the model selection in 3 consecutive steps - similar variables are grouped together and tested together in all combinations to get the best combination 
# groups are river (river density, Distance large river, permanent and seasonal water), topography (slope and elevation), and vegetation (NDVI and EVI)

# water 
water.formula <- list(~ Undisturbed_tropical_moist_forest + Period_fact + river_density_med_large, 
                    ~ Undisturbed_tropical_moist_forest + Period_fact + river_density_med_large+Distance_large_river, 
                    ~ Undisturbed_tropical_moist_forest + Period_fact + Distance_large_river+Permanent_and_seasonal_water,
                    ~ Undisturbed_tropical_moist_forest + Period_fact + river_density_med_large+Distance_large_river+Permanent_and_seasonal_water)

for(w in 1:length(water.formula)){
  model <- tIntPGOcc(occ.formula = water.formula[[w]], det.formula = det.formula, 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wA
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'occ.water')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

best.water <- model.comparison %>% filter(step == 'occ.water') %>% slice_min(wAIC) %>% pull(occ)


# topography 
topography.formula <- list(as.formula(paste(best.water, '+ mean_elev')), 
                          as.formula(paste(best.water, '+ mean_slope')),
                          as.formula(paste(best.water, '+ mean_elev + mean_slope')))

for(w in 1:length(topography.formula)){
  model <- tIntPGOcc(occ.formula = topography.formula[[w]], det.formula = det.formula, 
                     data = data.list, n.report = n.report, n.batch = n.batch, batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)
  # save all values to model comparison table
  output <- data.frame(occ=paste('~',as.character(model$occ.formula)[2]), # occ
                       det.transect=paste('~', as.character(model$det.formula[[1]])[2]), # det.transect
                       det.camera=paste('~', as.character(model$det.formula[[2]])[2]), # det.camera
                       wAIC=sum(waicOcc(model)[3]), # summed wAIC across both data sources
                       wAIC.transect=waicOcc(model)[3][1,], # 1st column is transect data, wAIC.transect
                       wAIC.camera=waicOcc(model)[3][2,], # 2nd column is camera data, wA
                       max.Rhat = max(unlist(model$rhat)),
                       min.ESS = min(unlist(model$ESS)),
                       step = 'occ.topography')
  model.comparison <- rbind(model.comparison, c(output))
  
  rm(model) # remove model object to clear memory
}

best.topography <- as.formula(model.comparison %>% filter(step == 'occ.topography') %>% slice_min(wAIC) %>% pull(occ))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Fit best model again ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n.report <- 2
batch.length <- 1500
n.batch <- 10
n.thin <- 3
n.chains <- 4


best.model <- tIntPGOcc(occ.formula = best.topography, det.formula = det.formula, 
                         data = data.list, n.report = n.report, n.batch = n.batch, 
                         batch.length = batch.length, n.thin = n.thin, n.chains = n.chains)

# call model
summary(best.model)


# visualise effect sizes using MCMCvis
library(MCMCvis)
MCMCplot(best.model$beta.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Occupancy Effect Sizes")
MCMCplot(best.model$alpha.samples, ref_ovl = TRUE, ci = c(50, 95),  main = "Detection Effect Sizes")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Predict Occupancy throughout the study area  ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# since the predictors were scaled before fitting the model, the envionmental covariates used to predict and extrapolate throughout Gola 
# have to be scaled with the same factors as the scaling was done before

# approach: take mean and sd from old data model was trained with - then take new environmental covariates and manually scale them with these factors ((envCovs - train.data.mean) / train.data.sd)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9.1 Prepare data for prediction ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9.1.1 SiteCovs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get site covariates and normalize them with the same parameters as scaling was done before model fitting 
op <- all.vars(best.model$occ.formula)[!all.vars(best.model$occ.formula) %in% names(year.site.covs.list)]
envCovs_pred <- envCovs_sf[,c(op, 'CellID')] %>% st_drop_geometry()
for(o in op){
  sd.sites <- sd(as.numeric(occ.covs[[paste0(o, '_unscaled')]])) # use the data the model was fitted with for calculation of mean and sd for normalization
  mean.sites <- mean(occ.covs[[paste0(o, '_unscaled')]])
  
  envCovs_pred[[o]] <- (envCovs_sf[[o]] -mean.sites) / sd.sites # apply above calc normalization factors to whole study area envCovs 
}

envCovs_pred <- envCovs_pred %>% 
  mutate(`(Intercept)` = 1)

# create an array
occ.preds <- c('(Intercept)', all.vars(best.model$occ.formula)[!all.vars(best.model$occ.formula) %in% names(year.site.covs.list)]) # get all occu predictors + intercept for prediction which vary at site not year-site level
X.0 <- array(NA, dim = c(nrow(envCovs_pred), length(best.model$seasons[[1]]), length(occ.preds)))
dimnames(X.0) <- list(CellID = as.character(envCovs_pred$CellID), Period = paste0("Period_", 1:length(best.model$seasons[[1]])), Occ.covs = occ.preds)

# loop over all site level covariates
for (o in occ.preds) {
  for (s in 1:length(best.model$seasons[[1]])) {
    X.0[, s, o] <- envCovs_pred[[o]]
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9.1.2 yearly SiteCovs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all.vars(best.model$occ.formula)[all.vars(best.model$occ.formula) %in% names(year.site.covs.list) ] # these are all the variables at yearly site level
year.occ.preds <- all.vars(best.model$occ.formula)[all.vars(best.model$occ.formula) %in% names(year.site.covs.list) & all.vars(best.model$occ.formula) != "Period_fact"] # these are all yearly site covs but without Period
year.occ.preds.array <- array(data = NA, dim = c(nrow(envCovs), length(best.model$seasons[[1]]), length(year.occ.preds)))
dimnames(year.occ.preds.array) <- list(CellID = as.character(envCovs_pred$CellID), Period = paste0("Period_", 1:length(best.model$seasons[[1]])), Occ.covs = year.occ.preds)

for(year.occ in year.occ.preds){
  
  # pool both years since normalization was done across years before model fitting
  year.envCovs <- occ.covs[[paste0(year.occ, '_unscaled')]]
  
  sd.sites <- sd(year.envCovs)
  mean.sites <- mean(year.envCovs)
  
  year.occ.preds.array[,,year.occ] <- (c(envCovs_sf[[paste0('JRC_ann_changes_', year.occ, '_2013')]], 
                                        envCovs_sf[[paste0('JRC_ann_changes_', year.occ, '_2021')]]) -mean.sites) / sd.sites # manual scaling
}

# manually add period_fact
period.fact <- array(data = rep((c(0, 1)), each = dim(year.occ.preds.array)[1]), 
      dim = c(nrow(envCovs), length(best.model$seasons[[1]]), 1), 
      dimnames = list(CellID = as.character(envCovs_pred$CellID), Period = paste0("Period_", 1:length(best.model$seasons[[1]])), Occ.covs = 'Period_fact'))

# bring yeary site covs and site covs together 
library(abind)
X.0.pred <- abind(X.0, year.occ.preds.array, period.fact, along = 3)
X.0.pred <- X.0.pred[,,c('(Intercept)', all.vars(best.model$occ.formula))] # this ensures, that the variables are in the correct order!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9.2 Predict ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# predict 
t.cols <- 1:2 # this indicates for which time periods we are predicting
pred_bm <- predict(best.model, type = 'occupancy', ignore.RE = F, X.0 = X.0.pred, t.cols = t.cols) # check, if 4 for t.cols is correct!

# prediction is a list with two components: psi.0.samples (the occurrence probability predictions) and z.0.samples (the latent occurrence predictions), each as 3D arrays with dimensions corresponding to MCMC sample, site, primary period
str(pred_bm)

plot_bm <- data.frame(CellID = numeric(), pred_mean = numeric(), Period = character())
for(p in 1:length(best.model$seasons[[1]])){
  prediction <- data.frame(CellID = 1:nrow(X.0), 
                           pred_mean = apply(pred_bm$psi.0.samples[, , p], 2, mean), 
                           pred_sd = apply(pred_bm$psi.0.samples[, , p], 2, sd),
                           Period = c('2011-2017', '2018-2025')[p])
  plot_bm <- bind_rows(plot_bm, prediction)
}

# bring into long format 
# plot_bm <- plot_bm %>% pivot_longer(cols = -c(CellID, Period), 
#                          names_to = 'pred.type', 
#                          values_to = 'pred.value')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9.3 Format predicted data and plot ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create sf
plot_bm_sf <- envCovs_sf %>% select(CellID) %>% 
  left_join(plot_bm, join_by(CellID))


# get information on data the model was trained with
preds <- setdiff(all.vars(best.model$occ.formula), "Period_fact") # all predictors except of Period_fact are important here 
pred.limits <- data.frame(variable = character(), Period = character(), train_min = numeric(), train_max = numeric())


for(p in 1:length(best.model$seasons[[1]])){
  for(variable in preds){
    # if else statement to handle yearly site covs and site covs differently
    if (is.matrix(occ.covs[[variable]])) { # yearly site covs should be a matrix 
      train_min <- min(occ.covs[[variable]][, p])
      train_max <- max(occ.covs[[variable]][, p])
    } else {
      train_min <- min(occ.covs[[variable]])
      train_max <- max(occ.covs[[variable]])}
    # save min and max values in a df
    period.limits <- data.frame(variable = variable, Period = c('2011-2017', '2018-2025')[p], 
                                train_min = train_min,train_max = train_max)
    # aggregate results
    pred.limits <- bind_rows(pred.limits, period.limits)
  }
}

# get values predictions were made for
pred.values <- data.frame(Period = character(), CellID = numeric())
for(p in 1:length(best.model$seasons[[1]])){
    df <- data.frame(Period = rep(c('2011-2017', '2018-2025')[p], times = dim(X.0.pred)[1]),
                     X.0.pred[,p,]) %>% mutate(CellID = row_number())
    df$X.Intercept. <- NULL
    df$Period_fact <- NULL
    pred.values <- bind_rows(pred.values, df)
    # rm(df)
}

# bring pred values into long format 
pred.values = pred.values %>%
  pivot_longer(cols = -c(CellID, Period),
    names_to = "variable",
    values_to = "value") 

# put the respective min and max values in there as well, specific to period and variable
pred.values <- pred.values %>% 
  left_join(pred.limits, join_by(variable, Period), relationship = 'many-to-many')


# flag all values that lie within the borders of what the model was trained for, flag extrapolations
flag <- pred.values %>% ungroup() %>% 
  mutate(extrapol = if_else(value >= train_min & value <= train_max, 1, 0.5)) %>% 
  group_by(CellID, Period) %>% 
  reframe(flag = sum(extrapol), 
            alpha = as.numeric(if_else(sum(extrapol) == length(preds), 1, 0.95)))

# plot prediction
library(ggspatial)

pred_map <- plot_bm_sf %>% 
  left_join(flag, join_by(CellID, Period)) %>% mutate(outline = if_else(alpha == 1, 'in', 'outside')) %>%
  ggplot() +
  geom_sf(aes(fill = pred_mean, alpha = alpha, col = outline)) +  # alpha is numeric (e.g., 1 or 0.6)
  #geom_sf(data = locs_transects %>% mutate(Period = if_else(Period == 1, '2011-2017', '2018-2025')), aes(linetype = "Transect lines"), color = "forestgreen", linewidth = 1) +
  #geom_sf(data = deploy_cam_visit %>% mutate(Period = if_else(Period == 1, '2011-2017', '2018-2025')) %>% group_by(SiteID, Period) %>% summarise(deploy_length = as.numeric(sum(Visit_length))), 
  #        aes(size = deploy_length), col = 'black') +
  scale_fill_viridis_c(option = "plasma", name = "Predicted \nOccupancy", limits = c(0, 1), direction = -1) +
  #scale_linetype_manual(name = "Transect surveys", values = c("Transect lines" = "solid")) +
  #scale_size_continuous(name = "Camera deployment \nlength in days", range = c(0.2, 1.5)) +
  scale_alpha_continuous(range = c(0.5, 1), guide = "none") +  # optional: hide alpha legend
  scale_color_manual(name = "Extrapolation",
    values = c("in" = "gray30", "outside" = "grey60")) +
  facet_grid(~Period) +
  theme_bw() +
  theme(legend.position = c(0.95, 0.25), # 
    panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), 
    axis.line = element_blank()) +
  labs(title = "Predicted Pygmy Hippo Occupancy Probability", 
    subtitle = "Based on a Static Multi-Year Integrated Occupancy Model fitted in spOccupancy", 
    x = 'Longitude', y = 'Latitude')


ggsave(filename = 'output/plots/PH_tIntPGOcc_predicted_map.jpg', plot = pred_map, height = 7, width = 16)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 10. Plot marginal effects plots for Occupancy   ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Attention! The covariates for prediction in the X.0 array should be organized in the same order as they were specified in the corresponding formula argument of tIntPGOcc.

# get all occu predictors + intercept for prediction which vary at site not year-site level
occ.preds <- c('(Intercept)', all.vars(best.model$occ.formula))
n <- 2000 # number of predictions per variable 
head(X.0.pred) # this is the old array for prediction in gola 
X.0.mrgnl.effects <- array(data = NA, dim = c(n, length(best.model$seasons[[1]]), length(occ.preds)))
dimnames(X.0.mrgnl.effects) <- list(NULL, Period = paste0("Period_", 1:length(best.model$seasons[[1]])), Occ.covs = occ.preds)


# use old array to built a new one with constant values 
for(o in occ.preds){
  for (s in 1:length(best.model$seasons[[1]])){
      X.0.mrgnl.effects[,s,o] <- rep(mean(X.0.pred[,s, o], na.rm = T), times = n) # this adds constant values for all predictors
  }
}

# check that everything went okay 
head(X.0.mrgnl.effects) 
t.cols <- 1:2 # this indicates for which time periods we are predicting
plot.mrgnl.effects <- data.frame(variable = character(), value = numeric(), pred.mean = numeric(), Period = character(), pred.upper = numeric(), pred.lower = numeric()) 

# actual prediction 
for(o in occ.preds){
  varying.array <- array(data = NA, dim = c(n, length(best.model$seasons[[1]])))
  dimnames(varying.array) <- list(NULL, paste0("Period_", 1:length(best.model$seasons[[1]])))
  for(s in  1:length(best.model$seasons[[1]])){
    varying.array[,s] <- seq(from = min(X.0.pred[,s,o], na.rm = T), to = max(X.0.pred[,s,o], na.rm = T), length.out = n)
  }
  pred.array <- X.0.mrgnl.effects
  pred.array[,,o] <- varying.array
  pred.mrgnl.effects <- predict(best.model, t.cols = t.cols, X.0 = pred.array, type = 'occupancy', ignore.RE = F)
  
  for(p in 1:length(best.model$seasons[[1]])){
    prediction.mrgnl.effects <- data.frame(variable = o, 
                             value = pred.array[,p,o], 
                             pred.mean = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, mean), 
                             pred.lower = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, quantile, probs = 0.025),
                             pred.upper = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, quantile, probs = 0.975),
                             Period = c('2011-2017', '2018-2025')[p])
    plot.mrgnl.effects <- bind_rows(plot.mrgnl.effects, prediction.mrgnl.effects)
  }
  # plot.mrgnl.effects <- bind_rows(plot.mrgnl.effects, prediction.mrgnl.effects)
}

# recode variable names with units
plot.mrgnl.effects <- plot.mrgnl.effects %>%
  mutate(Variable = recode(variable,
                           "Distance_large_river" = "Distance to large river in m",
                           'river_density_med_large' = 'Density of medium and large rivers',
                           "mean_elev" = "Elevation in m asl.",
                           "mean_slope" = "Slope",
                           "Permanent_and_seasonal_water" = "Fraction of permanent and seasonal water per grid cell",
                           "Undisturbed_tropical_moist_forest" = "Fraction of undisturbed tropical moist forest per grid cell")) 

# add unscaled values to marginal effects df
preds <- setdiff(all.vars(best.model$occ.formula), "Period_fact") # all predictors except of Period_fact are important here 
plot.mrgnl.effects$value_unscaled <- NA
for(pred in preds){
  for(p in 1:length(best.model$seasons[[1]])){
    if(is.matrix(occ.covs[[pred]])){ # alll yearly site covs should be a matrix
      value_unscaled <- seq(from = min(occ.covs[[paste0(pred, '_unscaled')]][,p]), to = max(occ.covs[[paste0(pred, '_unscaled')]][,p]), length.out = n)
    } else{
      value_unscaled <- seq(from = min(occ.covs[[paste0(pred, '_unscaled')]]), to = max(occ.covs[[paste0(pred, '_unscaled')]]), length.out = n)
    }
    plot.mrgnl.effects$value_unscaled[plot.mrgnl.effects$variable == pred & plot.mrgnl.effects$Period == c('2011-2017', '2018-2025')[p]] <- value_unscaled
  }
}

# produce plot 
marginal.effects.vis <- plot.mrgnl.effects %>% filter(!variable %in% c('(Intercept)', 'Period_fact')) %>% 
  ggplot() +
  geom_ribbon(aes(x = value_unscaled, ymin = pred.lower, ymax = pred.upper, fill = Period), alpha = 0.3, linewidth = 1.2)+
  geom_line(mapping = aes(x = value_unscaled, y = pred.mean, color = Period), linewidth = 1.1) +
  scale_color_manual(values = c("2011-2017" = "#E69F00", "2018-2025" = "#56B4E9")) +
  scale_fill_manual(values = c("2011-2017" = "#E69F00", "2018-2025" = "#56B4E9")) +
  facet_wrap(~Variable, scale = 'free_x') + 
  theme_bw(base_size = 13) +
  labs(x = 'Environmental covariates', y = 'Predicted occupancy', title = 'Effect of environmental covariates on occupancy') +
  theme(#legend.background = element_rect(color = "gray30"), 
        #legend.position = c(0.85, 0.25)
        )
ggsave(filename = 'output/plots/PH_marginal_effects.jpg', plot = marginal.effects.vis, height = 8, width = 10)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 11. Export all essential data ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# save whole workspace 
save.image(file = "output/PH_workspace.Rdata")


# save all important objects 
save(best.model, ppc_gm, model.comparison, data.list, plot.mrgnl.effects, plot_bm_sf, 
     flag, envCovs_sf, ppc_result, fit_result, pred.limits, plot_bm, marginal.effects.vis, pred_map,
     file = 'PH_tIntPGOcc_output.RData')


