#### data preparation for camera trap data script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in Feb 2025


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations ######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install needed packages 
#install.packages('readxl')
#install.packages('Microsoft365R')

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in data and tidy everything up #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in all data sets 
excel_sheets('data/combined_PH_data2025_draft7.xlsx')
pres_cam <- read_excel(path = 'data/combined_PH_data2025_draft9.xlsx', sheet = 'combined_camera_data')
deploy_cam <- read_excel(path = 'data/combined_PH_data2025_draft9.xlsx', sheet = 'locations_cameras')
deploy_cam_lookup_IWTDarwin <- read.csv(file = 'data/rowlookup_darwinIWT_upt.csv') # this is a table with all sites and dates where the cameras have been deployed for IWT and Darwin projects

# format data types properly and get rid of unneeded columns - camera deployment data 
str(deploy_cam)
deploy_cam <- deploy_cam %>% mutate(Project = factor(ProjectIDName), 
                                    Country = factor(country), 
                                    River = factor(Riverlocation), 
                                    PlotIDName = factor(PlotIDName),
                                    Placename = factor(placename), 
                                    CameraIDName = factor(CameraIDName), 
                                    Deployment = as.Date(`DeploymentDate_dd/mm/yyyy`, format = '%d-%m-%Y'), 
                                    Collection = as.Date(`CollectionDate_dd/mm/yyyy`, format = '%d-%m-%Y'), 
                                    Camera_Functioning = factor(if_else(Camera_Functioning == "Camera Functioning", 'functioning', 'failure')), 
                                    Habitat = factor(overallhabitat)) %>% 
  select(Project, Country, Placename, River, PlotIDName, CameraIDName, Deployment, Collection, Camera_Functioning, UTM_X_meters, UTM_Y_meters, Habitat)
deploy_cam$CameraIDName[deploy_cam$CameraIDName == "na"] <- NA
deploy_cam$Habitat[deploy_cam$Habitat == "na"] <- NA

# there are a few inconsistencies with the data 
deploy_cam <- deploy_cam %>% mutate(Project = if_else(Project == 'REDD PygmyHippo 2013-2014', 'REDD_ARTP_PygmyHippo 2013-2014', Project), # assuming they are the same projects 
                                    CameraIDName = if_else(Project == 'Pygmy Hippo REDD CT 2019-2021', str_replace_all(CameraIDName, pattern = ' ', replacement = ''), CameraIDName))


# format data types properly and get rid of unneeded columns - camera presence data 
str(pres_cam) # this looks rather scattered, thus do it separately for each data set 
pres_cam <- pres_cam %>% mutate(Project = factor(project), 
                                Country = factor(country), 
                                Placename = factor(placename), 
                                DeploymentID = factor(deploymentID), 
                                CameraIDName = factor(CameraIDName), 
                                Obs_Date = as.Date(date,format = '%d/%m/%Y'), 
                                Obs_Time = format(as.POSIXct(as.numeric(time) * 86400, origin = Obs_Date, tz = 'UTC'), "%H:%M:%OS"), # introduces NAs since there are times labeled as 'unknown'
                                Obs_DateTime = case_when(
                                  is.na(Obs_Time) ~ as.POSIXct(paste0(Obs_Date, " 04:00:00"), format = "%Y-%m-%d %H:%M:%S"), # this makes up a observation time for all cases where no time was available
                                  TRUE ~ as.POSIXct(paste0(Obs_Date, " ", Obs_Time), format = "%Y-%m-%d %H:%M:%S")),
                                UTM_X_meters = x_coord, UTM_Y_meters = y_coord, 
                                habitat = factor(habitat), 
                                Season = factor(str_to_title(season))) %>% 
  select(Project, Country, Placename, DeploymentID, CameraIDName, Obs_DateTime, UTM_X_meters, UTM_Y_meters, Season, count, notes) 
pres_cam$Placename[pres_cam$Placename == 'na'] <- NA
pres_cam$DeploymentID[pres_cam$DeploymentID == 'na'] <- NA 
pres_cam$CameraIDName[pres_cam$CameraIDName == 'unknown'] <- NA 

# format data properly - lookup deployment table for IWT and Darwin morrow river project 
str(deploy_cam_lookup_IWTDarwin)
deploy_cam_lookup_IWTDarwin[deploy_cam_lookup_IWTDarwin ==''] <- NA # label empty rows with NA 
deploy_cam_lookup_IWTDarwin <- deploy_cam_lookup_IWTDarwin %>% drop_na() # remove tail of the table with empty values 
deploy_cam_lookup_IWTDarwin <- deploy_cam_lookup_IWTDarwin %>% 
  mutate(date = as.Date(date, format = '%d/%m/%Y'), 
         Placename = as.factor(placename), 
         DeploymentID = as.factor(deploymentID)) 
deploy_cam_lookup_IWTDarwin <- deploy_cam_lookup_IWTDarwin %>% # calculate deployment period per site and deployment 
  group_by(Placename, DeploymentID) %>% 
  summarise(Deployment = min(date), 
            Collection = max(date), .groups = 'drop')
str(deploy_cam_lookup_IWTDarwin)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Find matching columns between tables to combine them #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Basel Zoo data
pres_cam %>% filter(Project == 'Basel Zoo PygmyHippo 2018-2020') %>% print(n = 24)
deploy_cam %>% filter(Project == 'Basel Zoo PygmyHippo 2018-2020') # its possible to join the data via CameraIDName or coordinates

# ARTP data - project names are different in both tables 
pres_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% select(-notes)
deploy_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') # only join via coordinates possible 

# Pygmy Hippo REDD CT 2019-2021
pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% select(-notes) %>% print(n = 26) # the placename could be an ident representation of the site, but no equivalent in cam_deploy data 
deploy_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% print(n = 60) # connection via CameraIDName, but unfortunately coordinates may not work
deploy_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% left_join(pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), join_by(UTM_X_meters, UTM_Y_meters)) 
# DWCN23 has been in the field for only one day?

pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021', CameraIDName %in% c('PH09', 'DWCN08')) %>% select(-notes) %>% print(n = 26) # the placename could be an ident representation of the site, but no equivalent in cam_deploy data 
deploy_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021', CameraIDName %in% c('PH09', 'DWCN08')) %>% print(n = 60) # connection via CameraIDName, but unfortunately coordinates may not work

ggplot() +
  geom_point(data = deploy_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 2) +
  geom_point(data = pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 1, color = 2)

unique_sites <- pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% distinct(CameraIDName, UTM_X_meters, UTM_Y_meters)
semi_join(unique_sites, deploy_cam, join_by(UTM_X_meters, UTM_Y_meters))  # look up which presences are included in the deployment table - PH09, DWCN23, DWCN18 are not included 
unique_sites
# deploy_cam %>% filter(CameraIDName %in% c('PH09', 'DWCN23', 'DWCN18'))

pres_cam %>% filter(Project == 'Darwin_morroRiver') %>% select(-notes) # no coordinates available 
deploy_cam %>% filter(Project == 'Darwin_morroRiver') # ident placename between both tables (placename and CamerIDName are the same in cam_deploy), deployment data for second deployment periods missing

pres_cam %>% filter(Project == 'IWT_CF') %>% select(-notes) #
deploy_cam %>% filter(Project == 'IWT_CF') %>% print(n = 31)# ident Placenames (which is also the same as PlotIDName), lots of observations are outside the deployment period


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. prepare each data set #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

unique(pres_cam$Project)
unique(deploy_cam$Project)

# create SiteID
pres_cam$SiteID <- NA
deploy_cam$SiteID <- NA

# BasalZoo2024 project, Plotname is not equal between tables, but Placename corresponds to 

pres_cam <- pres_cam %>% mutate(Obs_DateTime = if_else(Project == 'BasalZoo_2024' & Placename == 'BZ404206', update(Obs_DateTime, year = 2024), Obs_DateTime)) # there is one mistake in the year, agreed for with Felicity
deploy_cam <- deploy_cam %>% mutate(Collection = if_else(Project == 'BasalZoo_2024' & PlotIDName == 'BZ403205', update(Collection, year = 2024), Collection)) # there is again a year mistake, agreed for with Felicity

pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'BasalZoo_2024', Placename, SiteID), 
                    CameraIDName = case_when(Placename == 'BZ404206' ~ 'DWCN31', # just transfer two cases where Camera ID name is missing manually 
                                             Placename == 'BZ403205' ~ 'IWT11', 
                                             TRUE ~ CameraIDName)) 
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'BasalZoo_2024', PlotIDName, SiteID), 
                      CameraIDName = if_else(Project == 'BasalZoo_2024', str_replace_all(CameraIDName, ' ', ''), CameraIDName))

# artp1 project, SiteID made from Project name and Placename/PlotIDname/DeploymentID
pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'artp_p1', paste0(Project, '_', Placename), SiteID), 
                                CameraIDName = if_else(Placename == 'P1255' & Project == 'artp_p1', 'C17', CameraIDName))
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'artp_p1', paste0(Project, '_', Placename), SiteID))

# artp_p2, SiteID made from Project name and Placename/PlotIDname/DeploymentID
pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'artp_p2', paste0(Project, '_', Placename), SiteID), 
                                CameraIDName = if_else(Placename == 'P9137' & Project == 'artp_p2', 'C25', CameraIDName)) # transfer CameraID name manually
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'artp_p2', paste0(Project, '_', Placename), SiteID)) 

# redd sp1, I'm not sure if 
pres_cam <- pres_cam %>% mutate(CameraIDName = if_else( Project == 'REDD_SP1' & Placename == 'P5515', 'C10', CameraIDName), 
                    SiteID = if_else(Project == 'REDD_SP1', paste0(Project, '_', Placename), SiteID))
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'REDD_SP1', paste0(Project, '_', Placename), SiteID)) 

# redd sp2, sp3, sp4, darwin_13_17,  Darwin19_22, there are only presences in sp4
deploy_cam <- deploy_cam %>% 
  mutate(CameraIDName = if_else(Project %in% c('REDD_SP2', 'REDD_SP3', 'REDD_SP4', 'darwin_13_17', 'Darwin19_22'), str_replace_all(CameraIDName, ' ', ''), CameraIDName), 
                      SiteID = if_else(Project %in% c('REDD_SP2','REDD_SP3', 'REDD_SP4', 'darwin_13_17',  'Darwin19_22'), paste0(Project, '_', Placename), SiteID)) 
pres_cam <- pres_cam %>% 
  mutate(SiteID = if_else(Project %in% c('REDD_SP2','REDD_SP3', 'REDD_SP4', 'darwin_13_17',  'Darwin19_22'), paste0(Project, '_', Placename), SiteID))

# Basel Zoo project, there are a few issues, solve them
pres_cam <- pres_cam %>% mutate(Placename = if_else(Project == 'Basel Zoo PygmyHippo 2018-2020', 'Moa_Moa River', Placename), # transfer Placename to deployment df 
                    SiteID = if_else(Project == 'Basel Zoo PygmyHippo 2018-2020', paste(Placename, CameraIDName, UTM_X_meters, UTM_Y_meters, sep = '_'), SiteID)) 
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'Basel Zoo PygmyHippo 2018-2020', paste(Placename, CameraIDName, UTM_X_meters, UTM_Y_meters, sep = '_'), SiteID)) 

# REDD ARTP 2023-2014, transfer columns CameraIDName, Placename from deploy_cam to pres_cam - only possible via coordinates (which have one mismatch)
pres_cam$UTM_X_meters[pres_cam$UTM_X_meters == 292289] <- 291289 # replace this singe value , since it must be a typo, agreed with Felicity

pres_artp_replacement <- pres_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% select(-Placename, -CameraIDName) %>%  # transfer columns from deploy_cam to pres_cam
  left_join(deploy_cam %>% select(UTM_X_meters, UTM_Y_meters, CameraIDName, Placename), join_by(UTM_X_meters, UTM_Y_meters))
pres_cam <- pres_cam %>% filter(Project != 'REDD_ARTP_PygmyHippo 2013-2014') %>% rbind(pres_artp_replacement) # save improved columns in pres_cam df
rm(pres_artp_replacement)

pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'REDD_ARTP_PygmyHippo 2013-2014', paste(Placename, CameraIDName, sep = '_'), SiteID)) # create SiteID's 
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'REDD_ARTP_PygmyHippo 2013-2014', paste(Placename, CameraIDName, sep = '_'), SiteID))

# Pygmy Hippo REDD CT 2019-2021, there are a few issues which need troubleshooting!
pres_cam <- pres_cam %>% filter(!Placename %in% c('R377302(1)', 'R377302(2)')) # remove two observations which cause issues 
pres_cam$UTM_X_meters[pres_cam$UTM_X_meters == 240128] <- 240126 # correct very minor mistake in coordinates, agreed for with Felicity
pres_cam$UTM_Y_meters[pres_cam$UTM_Y_meters == 830246] <- 830245 # same here 
pres_cam$UTM_X_meters[pres_cam$UTM_X_meters == 251504] <- 251481 # correct a location error for Placename R514137, CameraIDName DWCN23, agreed for with Felicity
pres_cam$UTM_Y_meters[pres_cam$UTM_Y_meters == 813575] <- 813792 # same here 
# pres_cam$UTM_Y_meters[pres_cam$UTM_Y_meters == 830212] <- 830044 # PH09 has mismatch in presence and deployment data, I just aligned the y value to the deployment data value #### NEEDS GO FROM FELICITY
pres_cam <- pres_cam %>% mutate(Obs_DateTime = if_else(Project == 'Pygmy Hippo REDD CT 2019-2021' & CameraIDName == 'PH09', update(Obs_DateTime, year = 2020), Obs_DateTime)) # agreed for with Felicity

pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'Pygmy Hippo REDD CT 2019-2021', paste(Placename, UTM_X_meters, UTM_Y_meters, CameraIDName, sep = '_'), SiteID)) # unfortunately, I have to take 4 variables in to create an unique SiteID
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'Pygmy Hippo REDD CT 2019-2021', paste(Placename, UTM_X_meters, UTM_Y_meters,CameraIDName, sep = '_'), SiteID)) 

# Darwin_morroRiver, there are a few issues in the data!
removal <- c('D050459', 'D070481', 'IWT025422', 'IWT912319') # sites where cameras failed - exclude
removal_lookup <- c('IWT993381-1dpl') # certain deployments that have to be removed from lookup table
deploy_cam_lookup_IWTDarwin <- deploy_cam_lookup_IWTDarwin %>% filter(!DeploymentID %in% removal_lookup, !Placename %in% removal) # there were a failures and no pictures were taken at all, thus exclude, agreed with Felicity
deploy_cam <- deploy_cam %>% filter(!Placename %in% removal) # These cameras failed, agreed for with Felicity
pres_cam %>% filter(Placename %in% removal, DeploymentID %in% removal_lookup) # check, that no presences were removed - if nothing is filtered out: all good!

# Darwin_morroRiver, create SiteID from Placename in both tables 
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'Darwin_morroRiver', Placename, SiteID)) # define SiteID for this project in both df's just as placename
pres_cam <- pres_cam %>% mutate(SiteID = if_else(Project == 'Darwin_morroRiver', Placename, SiteID)) 

# IWT_CF, just create the SiteID from Placename - attention: the deployments are not transferred yet!
deploy_cam <- deploy_cam %>% mutate(SiteID = if_else(Project == 'IWT_CF', Placename, SiteID))
pres_cam <- pres_cam %>% mutate(SiteID = ifelse(Project == 'IWT_CF', Placename, SiteID))

#### removed a few sites: D050459 (done), D070481 (done), IWT025422, IWT993381-1dpl (done) and IWT912319 (done), 
#### do not remove:  IWT900306, IWT006396, D046462, IWT868306 
#### still answer needed for IWT846293, IWT803255, D116537_1dpl, # once I got a reply, include them!

# compare the two tables deploy_cam_lookupIWTDarwin and deploy_cam, in each table is data which is not in the other one - make them similar!
anti_join(deploy_cam %>% filter(Project %in% c('Darwin_morroRiver', 'IWT_CF')), deploy_cam_lookup_IWTDarwin, join_by(Placename)) # the three sites (IWT846293, IWT803255, IWT025422) are available in deploy cam but not in the lookup table
anti_join(deploy_cam_lookup_IWTDarwin, deploy_cam, join_by(Placename)) # 4 deployments from the lookup table are not in the deployment table: D046462-3dpl, D070481_1dpl, IWT006396-1dpl, IWT868306-1dpl and IWT900306-1dpl

# transfer 3 sites from deploy_cam to lookup table
deploy_cam_lookup_IWTDarwin <- deploy_cam %>% 
  filter(SiteID %in% c('IWT846293', 'IWT803255', 'IWT993381')) %>% 
  select(Placename, Deployment, Collection) %>% 
  mutate(DeploymentID = paste0(Placename, '-1dpl')) %>% 
  rbind(deploy_cam_lookup_IWTDarwin)


##### I THINK THIS IS NOT AN ISSUE ANYMORE; BUT KEEP THE BLOCK ######

# transfer 4 deployments from deploy_cam_lookupIWTDarwin to deploy_cam data frame 
# deploy_cam <- deploy_cam_lookup_IWTDarwin %>% 
#  filter(DeploymentID %in% c('D046462-3dpl', 'IWT006396-1dpl', 'IWT868306-1dpl', 'IWT900306-1dpl')) %>% 
#  mutate(SiteID = Placename, River = 'Moro', Country = 'SierraLeone', Project = factor(if_else(str_starts(Placename, 'IWT'), 'IWT_CF', 'Darwin_morroRiver'))) %>% 
#  bind_rows(deploy_cam) %>% 
#  distinct() %>%
#  select(-DeploymentID)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. ensure, that deploy_cam and the Darwin and IWT lookup table are identical #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this ensures that the part further down runs properly to create absences

head(deploy_cam)
head(deploy_cam_lookup_IWTDarwin)

# transfer camera functioning from deploy_cam to deploy_cam_lookup_IWT
deploy_cam_lookup_IWTDarwin <- deploy_cam_lookup_IWTDarwin %>% mutate(SiteID = Placename) %>% left_join(deploy_cam %>% select(SiteID, Deployment, Collection, Camera_Functioning), join_by(SiteID, Deployment, Collection)) 

# check, that deployment periods are identical between lookup and deploy data

deploy_cam <- deploy_cam %>%  # for D116537_1dpl and IWT993381-1dpl and IWT912319 the collection dates in deploy_cam might not be correct, just take the lookup table version and remove other data from deploy_cam
  filter(!Project %in% c('Darwin_morroRiver', 'IWT_CF')) %>% 
  rbind(deploy_cam_lookup_IWTDarwin %>% select(-Camera_Functioning, -DeploymentID) %>% left_join(deploy_cam %>% select(-Deployment, -Collection), join_by(SiteID, Placename)))
# rm(deploy_cam_lookup_IWTDarwin)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. export/save prepared data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.csv(pres_cam, file = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_pres_cam_data.csv')
write.csv(deploy_cam, file = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_deploy_cam_data.csv')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. create visits from deployment periods  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot an overview of all deployments with their observations of PH 
ggplot() +
  geom_point(data = deploy_cam, mapping = aes(x= UTM_X_meters, y = UTM_Y_meters), size = 2) +
  geom_point(data = pres_cam, mapping = aes(x= UTM_X_meters, y = UTM_Y_meters), color = 'red', size =.7) +
  labs(title = 'Camera Trap Deployments and their Pygmy hippo Observations in Gola')

# create duplicate visits in deployment data 
deploy_cam_expanded <- deploy_cam %>%
  rowwise() %>%
  mutate(Visit_start = list(seq(Deployment, Collection, by = "20 days"))) %>%  # split in 20 days intervals, however, there are mostly days remaining as tail
  unnest(Visit_start) %>% 
  group_by(SiteID, Deployment, Collection) %>% # combination of SiteID, Deployment and Collection works as deploymentID
  mutate(Visit_end = lead(Visit_start, default = last(Collection)), 
         Visit_start = as.POSIXct(paste0(Visit_start, ' 12:00:00'), format = '%Y-%m-%d %H:%M:%OS'), # make meaningful time up where cameras could have been deployed and collected to avoid issues with matching pictures
         Visit_end = as.POSIXct(paste0(Visit_end, ' 12:00:00'), format = '%Y-%m-%d %H:%M:%OS'), 
         Visit_length = difftime(Visit_end, Visit_start, units = 'days')) %>%  
  group_by(SiteID) %>% 
  mutate(Visit = row_number()) %>% # this creates a consecutive number for each visit for each site across deployments
  filter(!Visit_length == 0) # remove cases where 0 days were left as overshoot, 17 cases
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. transfer presences into deployment data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# transfer presences to deployment data, this assumes, that the observations lie within the deployment period
#### THIS IS A CRITICAL STEP AND I WANT THIS TO BE DOUBLECHECKED BY SOMEONE ELSE!
deploy_cam_expanded <- deploy_cam_expanded %>% ungroup() %>%
  left_join(pres_cam %>% select(Project, SiteID, Obs_DateTime), join_by(Project, SiteID), multiple = 'all', relationship = 'many-to-many') %>% 
  filter(Visit_start <= Obs_DateTime & Obs_DateTime <= Visit_end) %>% 
  group_by(SiteID, Visit) %>% 
  summarise(Presences = n_distinct(Obs_DateTime), .groups = 'drop') %>% # this produces a count of pictures, but if distinct (not just multiple pictures in a row) detections of pygmy hippos are needed, use n_distinct() here 
  right_join(deploy_cam_expanded, join_by(SiteID, Visit)) %>% 
  mutate(Presences = if_else(is.na(Presences), 0, Presences)) %>% arrange(Project, SiteID, Visit) 

# this is a check that no presences got lost while #### THIS HAS TO BE WORKED OUT STILL!
pres_cam %>% distinct(SiteID, as.Date(Obs_DateTime)) # this checks the number of distinct presences per day
deploy_cam_expanded %>% pull(Presences) %>% sum()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. create plots of deployment periods with all presences mapped on them ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create plots for each data set as checks 
deploy_pres_plot <- list()
for(i in 1:length(unique(deploy_cam$Project))){
  deploy_pres_plot[[i]] <- ggplot() +
    geom_segment(data = deploy_cam %>% filter(Project == unique(deploy_cam$Project)[i]), # %>% filter(Project %in% c('Pygmy Hippo REDD CT 2019-2021'))
                 aes(x = Deployment, xend = Collection, y = SiteID, yend = SiteID), linewidth = 2, color = "blue") +
    geom_point(data = pres_cam %>% filter(Project == unique(deploy_cam$Project)[i]), # %>% filter(Project %in% c('Pygmy Hippo REDD CT 2019-2021'))
               aes(x = as.Date(Obs_DateTime), y = SiteID), color = "red", size = 3) + 
    facet_wrap(~Project, scales = "free")
  names(deploy_pres_plot)[i] <- unique(deploy_cam$Project)[i]
}

for(i in 1:length(deploy_pres_plot)){
  plot(deploy_pres_plot[[i]])}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. create a plot that shows distribution of camera trapping days ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

days <- deploy_cam %>% rowwise() %>% mutate(Date = list(seq(Deployment, Collection, by = "1 days"))) %>%  # split in 20 days intervals, however, there are mostly days remaining as tail
  unnest(Date) %>% 
  ggplot() +
  geom_histogram(mapping = aes( x = Date), binwidth = 1) + 
  theme_bw() +
  labs(title = "Number of Cameras per Day Across Time", x = "Time", y = "Number of Deployed Cameras per Day")
pres_all <- pres_cam %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Obs_DateTime)) +
  theme_bw()  +
  labs(title = "Number of Presences Across Time", subtitle = "Infinite number of signs allowed per Camera trap & date", x = "Time", y = "Number of Presences")
pres_filt <- pres_cam %>% mutate(Date = as.Date(Obs_DateTime)) %>% group_by(Project, SiteID, Date) %>% slice(1) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Obs_DateTime))+
  theme_bw()  +
  labs(title = "Number of Presences Across Time", subtitle = "Only one presence per camera site and date allowed", x = "Time", y = "Number of Presences")

library(ggpubr)
overview_cam <- ggarrange(days, pres_all, pres_filt, ncol = 2, nrow = 2, common.legend = F)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/Data_overview_cameras.jpg', plot = overview_cam, width = 8, height = 6)


