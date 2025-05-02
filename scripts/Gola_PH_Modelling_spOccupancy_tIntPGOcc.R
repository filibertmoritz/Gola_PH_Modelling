#### Pygmy hippo static multi-season integrated occupancy model in spOccupancy  
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in April 2025

# data has roughgly been prepared elsewhere


# check again how the variables are coded and with scaled and unscaled variables, scale only numeric 
# check priors and inits and set them to appropriate values 



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
# pres_opp <- fread(file = 'data/PH_prepared_pres_opp_data.csv', stringsAsFactors = T) %>% select(-V1)
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
deploy_cam %>% # unfortunately, it is not possible to create time periods where projects aren't split up
  group_by(Project) %>% 
  summarise(min_deploy = min(Deployment), max_collect = max(Collection)) %>% 
  ggplot() +
  geom_segment(aes(x = min_deploy, xend = max_collect, y = Project, yend = Project), linewidth = 2, color = "blue") 

# try to not split up deployment periods
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
  
# calculate the effort for each data set and period 

# calculate the number of presences for each data set and period 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Prepare presence camera trap data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create an index for the period - this woueld have to be improved for more automation and more survey periods
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
         # Period_Transect = as.factor(if_else(year(DateTime_End) <= 2017, '2011-2017', '2018-2024')),
         Project_Transect_fact = as.factor(Project), 
         Project_Transect_num_unscaled = as.numeric(if_else(Project =='Pygmy Hippo ARTP_REDD 2013-2014', 0, 1))) # %>% select(-transect_length) 
locs_transects <- locs_transects %>% mutate(Year_Transect_num = as.numeric(scale(Year_Transect_num_unscaled)), # scale numeric variables
                                            Date_Transect = as.numeric(scale(Date_Transect_unscaled)), 
                                            Transect_Length = as.numeric(scale(Transect_Length_unscaled)), 
                                            Project_Transect_num = as.numeric(scale(Project_Transect_num_unscaled)), 
                                            Season_Transect_num = as.numeric(scale(Season_Transect_num_unscaled)))


# extract det.covs data - observation.covariates
periods <- 1:max(locs_transects$Period)
visits <- 1:max(locs_transects$Cell_visit)
cells <- sort(unique(locs_transects$CellID))
obs.covs.selection <- c('Date_Transect', 'Transect_Length', 'Project_Transect_num', 'Project_Transect_fact', 'Season_Transect_fact', 'Season_Transect_num','Year_Transect_num', 'Year_Transect_fact', 'Occu') 
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
 ####### another solution ######
#names(locs_transects)
#obs.det.selection <- c('Date_Transect', 'Date_Transect_unscaled', 'Transect_Length', 'Transect_Length_unscaled')
#obs.det.covs <- list()
#cell.order <- sort(unique(locs_transects$CellID))
#n.periods <- max(locs_transects$Period)
#all.data <- data.frame(CellID = rep(cell.order, times = n.periods), 
#                       Period = rep(1:n.periods, each = length(cell.order)))

#for(obs.det in obs.det.selection){
#  # for each observation covariate create an array with NA's
#  a <- array(NA, dim = c(length(cell.order), n.periods, max(locs_transects$Cell_visit)))
#  dimnames(a) <- list(cell.order, paste0("Period_", 1:n.periods), paste0("Visit_", 1:max(locs_transects$Cell_visit)))

#  for(v in 1:max(locs_transects$Cell_visit)){
#    visits <- locs_transects %>% filter(Cell_visit == v) %>% 
#      st_drop_geometry() %>% 
#      ungroup() %>% 
#      select(CellID,!!sym(obs.det)) %>% 
#      right_join(all.data, join_by(CellID)) %>%
#      arrange(Period, CellID)

#    dat <- visits %>% 
#      pivot_wider(names_from = Period, names_prefix = 'Period_', 
#                  values_from = !!sym(obs.det)) %>% 
#      arrange(CellID)

#    #dat <- dat %>% 
#    #  right_join(tibble(CellID = cell.order), by = "CellID") %>%
#    #  arrange(CellID)

#    a[,,v] <- as.matrix(dat[,2:ncol(dat)])
#  }
#  obs.det.covs[[obs.det]] <- a
#}

# check that everything went right 
#obs.det.covs$Date_Transect_unscaled[,,1]


#####


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.1.2 Camera Trapping Data ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calculate covariates 
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% ungroup()%>% 
  mutate(Project_Camera = as.factor(Project), 
         # Project_lumped_Camera = as.factor(case_when(Project %in% c('BasalZoo_2024','Basel Zoo PygmyHippo 2018-2020') ~ 'BaselZoo', # fill a more meaningful categorisation in
         #                                             Project %in% c('REDD_SP1', 'REDD_SP2', 'REDD_SP3', 'REDD_SP4', 'Pygmy Hippo REDD CT 2019-2021','artp_p1', 'artp_p2', 'REDD_ARTP_PygmyHippo 2013-2014') ~ 'REDD_ARTP', 
         #                                             Project %in% c('Darwin19_22', 'Darwin_morroRiver', 'darwin_13_17') ~ 'Darwin', 
         #                                             Project %in% c('IWT_CF') ~ 'IWT')),
         Project_Focus_fact = as.ordered(if_else(Project %in% c('Pygmy Hippo REDD CT 2019-2021', 'REDD_ARTP_PygmyHippo 2013-2014', 
                                                                'Darwin_morroRiver', 'Basel Zoo PygmyHippo 2018-2020', 'BasalZoo_2024'), 
                                                 'Pygmy_hippo', 'Other')),
         Project_Focus_num_unscaled = as.numeric(if_else(Project_Focus_fact == 'Pygmy_hippo', 1, 0)), # Pygmy hippo is coded as 1, other as 0
         Date_Camera_unscaled = yday(Visit_start), 
         Year_Camera_fact = as.factor(year(Visit_start)), # survey year coded as factor
         Year_Camera_num_unscaled = as.numeric(year(Visit_start)-2011), # survey year coded as numeric
         Season_Camera_fact = as.factor(if_else(month(Visit_start) %in% 5:10, 'Wet', 'Dry')),
         Season_Camera_num_unscaled = as.numeric(if_else(month(Visit_start) %in% 5:10, 1, 0)), # wet season coded as 1, dry as 0
         # Period_Camera = as.factor(if_else(year(Visit_end) <= 2017, '2011-2017', '2018-2024')),
         Trapping_Days_unscaled = as.numeric(Visit_length)) # %>% select(-transect_length) # %>% select(-Visit_length) 
deploy_cam_visit_occu <- deploy_cam_visit_occu %>% mutate(Year_Camera_num = as.numeric(scale(Year_Camera_num_unscaled)), # scale numeric variables
                                                          Date_Camera = as.numeric(scale(Date_Camera_unscaled)), 
                                                          Trapping_Days = as.numeric(scale(Trapping_Days_unscaled)), 
                                                          Project_Focus_num = as.numeric(scale(Project_Focus_num_unscaled)), 
                                                          Season_Camera_num = as.numeric(scale(Season_Camera_num_unscaled)))


# extract det.covs on observation level - this takes some time due to 53 visits!
periods <- 1:max(deploy_cam_visit_occu$Period)
visits <- 1:max(deploy_cam_visit_occu$Cell_visit)
cells <- sort(unique(deploy_cam_visit_occu$CellID))
obs.covs.selection <- c('Date_Camera', 'Trapping_Days', 'Project_Focus_fact', 'Project_Focus_num', 'Season_Camera_fact', 'Season_Camera_num', 'Year_Camera_num', 'Year_Camera_fact', 'Occu')
obs.covs.camera <- replicate(length(obs.covs.selection),
                          array(NA, dim = c(length(cells), max(periods), max(visits))), 
                          simplify = F)
names(obs.covs.camera) <- obs.covs.selection

for(obs.det in obs.covs.selection){
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
y



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3 Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3.1 Site-Level Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get sites (Cells) for each data source
sites.transect <- sort(unique(locs_transects$CellID))
sites.camera <- sort(unique(deploy_cam_visit_occu$CellID))

# df with all envCovs for all sites (across transect and camera data)
envCovs_all_unscaled <- bind_rows(envCovs %>% filter(CellID %in% sites.transect) %>% arrange(CellID) %>% mutate(Type = 'Transect'), 
          envCovs %>% filter(CellID %in% sites.camera) %>% arrange(CellID) %>% mutate(Type = 'Camera'))

# scale numeric variables  
envCovs_all <- envCovs_all_unscaled %>%
  mutate(across(!c(CellID, Reserve_Type, Type), ~ scale(.) %>% as.vector()), 
         Reserve_Type=as.factor(Reserve_Type))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### 5.3.2 Season-Level/Yearly-Site Occ Covs ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# keep season level covariates unscaled here and scale later on covariate level, not within each year of the covariate seperately

# prepare loop
periods <- 1:max(deploy_cam_visit_occu$Period) # number of periods 
cells <- c(sites.transect, sites.camera) # site/cell indices 
# names(envCovs_sf)[grep('JRC_ann_changes', names(envCovs_sf))] # this is if only the years that are mn
year_preds <- names(envCovs_sf)[grep("JRC_ann_changes.*(2013|2021)", names(envCovs_sf))]
year.site.covs <- unique(gsub('_Dec2013|_Dec2021','', year_preds)) # fill in the years data should be extracted for (which refer to the periods)
year.site.covs.list <- list()


# call loop
for(y in year.site.covs){
  m <- matrix(data = NA, nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods))) 
  for(p in periods){
      m[,p] <- envCovs_all_unscaled[[ year_preds[grep(y, year_preds)][p] ]] # take the unscaled df with envCovs
  }
  m_scaled <- scale(as.vector(m)) # scale all variables 
  m <- matrix(m_scaled, nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))  
  year.site.covs.list[[y]] <- m
}

# manually add year or year-period as yearly-site covariate - two versions, numeric and scaled and as factor
year.site.covs.list[['Period_num']] <- matrix(data = as.numeric(scale(rep(periods, each = length(cells)))), nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods)))
year.site.covs.list[['Period_fact']] <- matrix(data = rep(c('2011-2017','2018-2025'), each = length(cells)), nrow = length(cells), ncol = length(periods), dimnames = list(cells, paste0('Period_',periods)))


# create occ.covs list which compromises all data 
occ.covs <- c(as.list(envCovs_all), year.site.covs.list) # add all matrices from a list to the occ.covs list by using c()


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
inits.list = list(z = 1, # z is for latent variable (here occupancy), start with 1 for all sites occupied
                  beta = 0, # for ecological state model, start occupancy covariates at 0, length is the number of occu covs (n.occ.covs)
                  alpha = list(Transect = 0, # for each data source a list with initial values for det model
                               Camera = 0), 
                  sigma.sq.psi = NULL, # for random effects in occurrence model as list
                  sigma.sq.p = NULL, # for random effects detection model as list
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


# call global model 

names(occ.covs)
names(det.covs$Transect)

m1 <- tIntPGOcc(occ.formula = ~ river_density_med_large + Distance_large_river + 
                  mean_elev + JRC_transition_Degraded_forest_short_duration_disturbance_after_2014 + 
                  JRC_ann_changes_Undisturbed_tropical_moist_forest  + Period_num, 
          det.formula = list(Transect = ~Transect_Length + Date_Transect + I(Date_Transect^2) + Project_Transect_fact + (1| Season_Transect_num) + (1| Year_Transect_num), 
                             Camera = ~Trapping_Days + Date_Camera + I(Date_Camera^2) + Project_Focus_fact + (1| Season_Camera_num) + (1| Year_Camera_num)), 
          data = data.list, 
          n.batch = 5, 
          batch.length = 2000,
          n.chains = 3)

# call model summary 
summary(m1)

# Goodness-of-fit test, freeman-tukey and grouped by sites, described in https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml
ppc_m1 <- ppcOcc(m1, fit.stat = 'freeman-tukey', group = 1) # group 1 - groups values by site, group 2 - groups values per replicate, there is also chi squared available
summary(ppc_m1) # get bayes p-value

# visualise Goodness-of-fit 

# produce a model check plot, code taken from https://doserlab.com/files/spoccupancy-web/articles/modelfitting
ppc_result <- data.frame(fit = numeric(),fit.rep = numeric(),season = character(),dataset = character(),color = character())
for(d in 1:length(m1$det.formula)){
  for(s in 1:length(m1$seasons[[d]])){
    ppc_frame <- data.frame(fit = ppc_m1$fit.y[[d]][,s], fit.rep = ppc_m1$fit.y.rep[[d]][,s], 
                            season = s, dataset = names(m1$det.formula)[d], color = 'lightskyblue1')
    ppc_frame$color[ppc_frame$fit.rep > ppc_frame$fit] <- 'lightsalmon'
    ppc_result <- rbind(ppc_result, ppc_frame)
  }
}

# plot true vs fitted values 
ppc_result %>% mutate(season = if_else(season == 1, '2011-2017', '2017-2025')) %>% 
  ggplot() +
  geom_point(mapping = aes(x = fit, y = fit.rep, colour = color), size = 0.8) +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1.2) +
  facet_grid(season~dataset) +
  labs(x = 'True', y = 'Fit', title = 'True vs. Fitted Values for all Data Sources and Time Periods') +
  theme_bw()


# plot influential data points 
fit_result <- data.frame(fit = numeric(), numeric(),season = character(),dataset = character())
for(d in 1:length(m1$det.formula)){
  for(s in 1:length(m1$seasons[[d]])){
    fit_frame <- data.frame(fit = ppc_m1$fit.y.rep.group.quants[[d]][3,,][,s] - ppc_m1$fit.y.group.quants[[d]][3,,][,s], 
                            season = s, dataset = names(m1$det.formula)[d], SiteID = sites[[d]])
    fit_result <- rbind(fit_result, fit_frame)
  }
}

fit_result %>% mutate(season = if_else(season == 1, '2011-2017', '2017-2025')) %>% 
  ggplot() +
  geom_point(mapping = aes(y = fit, x = SiteID)) +
  facet_grid(season~dataset, scale = 'free') + 
  labs(y = 'Fit', title = 'Replicate - True Discrepancy') + 
  theme_bw()

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
##### 8. Predict Occupancy throughout study area  ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8.1 Prepare data for prediction ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# prepare site covariates 

# get all occu predictors + intercept for prediction which vary at site not year-site level
occ.preds <- c('(Intercept)', all.vars(m1$call$occ.formula)[!all.vars(m1$call$occ.formula) %in% names(year.site.covs.list)])

# prepare all occu covariates for prediction
envCovs_pred <- envCovs_sf %>%
  st_drop_geometry() %>%
  select(occ.preds[occ.preds != "(Intercept)"]) %>%
  mutate(across(everything(), ~ scale(.) %>% as.vector())) %>%  
  mutate(`(Intercept)` = 1) # extrac column with 1's for Intercept 

# create an array
X.0 <- array(NA, dim = c(nrow(envCovs_pred), length(m1$seasons[[1]]), length(occ.preds)))
dimnames(X.0) <- list(CellID = as.character(1:nrow(envCovs_pred)), Period = paste0("Period_", 1:length(m1$seasons[[1]])), Occ.covs = occ.preds)

# loop over all site level covariates
for (o in occ.preds) {
  for (s in 1:length(m1$seasons[[1]])) {
    X.0[, s, o] <- envCovs_pred[[o]]
  }
}



# prepare yearly site covariates 

# these are all the variables at site year level
all.vars(m1$call$occ.formula)[all.vars(m1$call$occ.formula) %in% names(year.site.covs.list) ]

occ.preds.year <- array(data = c(as.numeric(scale(c(envCovs$JRC_ann_changes_Undisturbed_tropical_moist_forest_Dec2013, envCovs$JRC_ann_changes_Undisturbed_tropical_moist_forest_Dec2013))), 
                                 # as.factor(rep(c('2011-2017', '2018-2025'), each = nrow(envCovs))), 
                                 as.numeric(scale(rep(c(1, 2), each = nrow(envCovs))))), 
      dim = c(nrow(envCovs), 2, 2), 
      dimnames = list(paste0(1:nrow(envCovs)), paste0("Period_", 1:length(m1$seasons[[1]])),c("JRC_ann_changes_Undisturbed_tropical_moist_forest", 'Period_num')))

# bring all prepared occ.covs for prediction together into one array 
library(abind) # this is a package which makes binding arrays much easier
X.0 <- abind(X.0, occ.preds.year, along = 3) 

# predict 
t.cols <- 1:2 # this indicates for which time periods we are predicting
pred_m1 <- predict(m1, type = 'occupancy', ignore.RE = T, X.0 = X.0, t.cols = t.cols) # check, if 4 for t.cols is correct!

# prediction is a list with two components: psi.0.samples (the occurrence probability predictions) and z.0.samples (the latent occurrence predictions), each as 3D arrays with dimensions corresponding to MCMC sample, site, primary period
str(pred_m1)

plot_m1 <- data.frame(CellID = numeric(), pred_mean = numeric(), Period = character())
for(p in 1:length(m1$seasons[[1]])){
  prediction <- data.frame(CellID = 1:nrow(X.0), 
                           pred_mean = apply(pred_m1$psi.0.samples[, , p], 2, mean), 
                           Period = c('2011-2017', '2018-2021')[p])
  plot_m1 <- bind_rows(plot_m1, prediction)
}

# create sf
plot_m1_sf <- envCovs_sf %>% left_join(plot_m1, join_by(CellID))


# fancier plot
library(ggspatial)
plot_m1_sf %>% 
  ggplot() +
  #annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_sf(aes(fill = pred_mean), # color = NA, 
          alpha = 1) +  # Use pred_mean for fill color
  scale_fill_viridis_c(option = "plasma", name = "Mean Predicted Occupancy", limits = c(0,1), direction = -1) + # scale_fill_viridis_c(option = "magma"), or replace "magma" with "inferno", "plasma", "cividis", 
  facet_grid(~Period) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), # add frame
        axis.line = element_blank()) +
  labs(title = "Predicted Pygmy Hippo Occupancy Probability", 
       subtitle = "Based on an Static Multi Year Integrated Occupancy Model fitted in spOccupancy", 
       x = 'Longitude', y = 'Latutude')
#ggsave(filename = 'output/plots/PH_hotspot_map_all_data.jpg', plot = hotspot_map, height = 6, width = 12)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Plot marginal effects plots for Occupancy   ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Attention! The covariates for prediction in the X.0 array should be organized in the same order as they were specified in the corresponding formula argument of tIntPGOcc.

# get all occu predictors + intercept for prediction which vary at site not year-site level
occ.preds <- c('(Intercept)', all.vars(m1$call$occ.formula))
n <- 500 # number of predictions per variable 
head(X.0) # this is the old array for prediction in gola 
X.0.mrgnl.effects <- array(data = NA, dim = c(n, length(m1$seasons[[1]]), length(occ.preds)))
dimnames(X.0.mrgnl.effects) <- list(NULL, Period = paste0("Period_", 1:length(m1$seasons[[1]])), Occ.covs = occ.preds)

head(X.0.mrgnl.effects)
head(X.0)

# use old array to built a new one with constant values 
for(o in occ.preds){
  for (s in 1:length(m1$seasons[[1]])){
      X.0.mrgnl.effects[,s,o] <- rep(mean(X.0[,s, o], na.rm = T), times = n) # this adds constant values for all predictors
  }
}

# check that everything went okay 
head(X.0.mrgnl.effects) # somehow its a bit weird, yearly site covariates seem to be constant across the periods - doublecheck this 
t.cols <- 1:2 # this indicates for which time periods we are predicting
plot.mrgnl.effects <- data.frame(variable = character(), value = numeric(), pred.mean = numeric(), Period = character(), pred.upper = numeric(), pred.lower = numeric()) 

# actual prediction 
for(o in occ.preds){
  varying.array <- array(data = NA, dim = c(n, length(m1$seasons[[1]])))
  dimnames(varying.array) <- list(NULL, paste0("Period_", 1:length(m1$seasons[[1]])))
  for(s in  1:length(m1$seasons[[1]])){
    varying.array[,s] <- seq(from = min(X.0[,s,o], na.rm = T), to = max(X.0[,s,o], na.rm = T), length.out = n)
  }
  pred.array <- X.0.mrgnl.effects
  pred.array[,,o] <- varying.array
  pred.mrgnl.effects <- predict(m1, t.cols = t.cols, X.0 = pred.array, type = 'occupancy', ignore.RE = F)
  
  for(p in 1:length(m1$seasons[[1]])){
    prediction.mrgnl.effects <- data.frame(variable = o, 
                             value = pred.array[,p,o], 
                             pred.mean = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, mean), 
                             pred.lower = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, quantile, probs = 0.025),
                             pred.upper = apply(pred.mrgnl.effects$psi.0.samples[, ,p], 2, quantile, probs = 0.975),
                             Period = c('2011-2017', '2018-2021')[p])
    plot.mrgnl.effects <- bind_rows(plot.mrgnl.effects, prediction.mrgnl.effects)
  }
  # plot.mrgnl.effects <- bind_rows(plot.mrgnl.effects, prediction.mrgnl.effects)
}


# produce a plot 
plot.mrgnl.effects %>% filter(!variable %in% c('(Intercept)', 'Period_num')) %>% 
  ggplot() +
  geom_ribbon(aes(x = value, ymin = pred.lower, ymax = pred.upper, fill = Period), alpha = 0.3, linewidth = 1.2)+
  geom_line(mapping = aes(x = value, y = pred.mean, color = Period)) +
  facet_wrap(~variable, scale = 'free_x') + 
  theme_bw()


