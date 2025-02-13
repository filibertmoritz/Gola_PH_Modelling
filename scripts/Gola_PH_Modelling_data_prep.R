#### data preparation script for pigmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install needed packages 
#install.packages('readxl')
#install.packages('Microsoft365R')

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
library(spOccupancy)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in data and tidy everything up 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in data 
excel_sheets('data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx')
pres <- read_excel('data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx', sheet = 'PYGMY HIPPO MASTER') # these are all presences we have across different projects
pres_cam <- read.csv('data/combined_PH_data2025_draft2_combined_camera_data.csv') # all presences from the camera trap data
cam_deploy <- read.csv('data/combined_PH_data2025_draft2_locations_cameras.csv') # all camera deployment locations and times of pres_cam

# format data types properly - pres
str(pres)

pres %>%
  mutate(Date = as.Date(VisitDate_dd_mm_yyyy), 
         Time = hms::as_hms(ifelse(is.na(VisitTime_hh_mm_ss), format(VisitDate_dd_mm_yyyy, '%H:%M:%S'), format(VisitTime_hh_mm_ss, '%H:%M:%S')))) %>%
  select(DatasetName, Date, Time, VisitDate_dd_mm_yyyy, VisitTime_hh_mm_ss)

pres[pres == -9999] <- NA # replace all -9999 with NA 

# format data types properly and get rid of unneeded columns - camera deployment data 
str(cam_deploy)
cam_deploy <- cam_deploy %>% mutate(Project = factor(ProjectIDName), 
                      Country = factor(country), 
                      River = factor(Riverlocation), 
                      PlotIDName = factor(PlotIDName),
                      Placename = factor(placename), 
                      CameraIDName = factor(CameraIDName), 
                      Deployment = as.Date(DeploymentDate_dd.mm.yyyy, format = '%d/%m/%Y'), 
                      Collection = as.Date(CollectionDate_dd.mm.yyyy, format = '%d/%m/%Y'), 
                      Camera_Functioning = factor(if_else(Camera_Functioning == "Camera Functioning", 'functioning', 'failure')), 
                      Habitat = factor(overallhabitat)) %>% 
  select(Project, Country, Placename, River, PlotIDName, CameraIDName, Deployment, Collection, Camera_Functioning, UTM_X_meters, UTM_Y_meters, Habitat)

# there are a few inconsistencies with the data 
cam_deploy <- cam_deploy %>% mutate(Project = if_else(Project == 'REDD PygmyHippo 2013-2014', 'REDD_ARTP_PygmyHippo 2013-2014', Project), # assuming they are the same projects 
                                    CameraIDName = if_else(Project == 'Pygmy Hippo REDD CT 2019-2021', str_replace_all(CameraIDName, pattern = ' ', replacement = ''), CameraIDName))

# format data types properly and get rid of unneeded columns - camera presence data 
str(pres_cam) # this looks rather scattered, thus do it separately for each data set 
pres_cam <- pres_cam %>% mutate(Project = factor(project), 
                                Country = factor(country), 
                                Placename = factor(placename), 
                                DeploymentID = factor(deploymentID), 
                                CameraIDName = factor(CameraIDName), 
                                Obs_Date = as.Date(date,format = '%d/%m/%Y'), 
                                UTM_X_meters = x_coord, UTM_Y_meters = y_coord, habitat = factor(habitat)) %>% 
  select(Project, Country, Placename, DeploymentID, CameraIDName, Obs_Date, time, UTM_X_meters, UTM_Y_meters, count, notes)

# solve issues in Basel Zoo PygmyHippo 2018-2020 data set - transfer placename from deployment data to presence data set 
pres_cam <- pres_cam %>% left_join(cam_deploy %>% 
                         filter(Project == 'Basel Zoo PygmyHippo 2018-2020') %>% 
                         select(Project, CameraIDName, UTM_X_meters, UTM_Y_meters, Placename) %>% 
                         rename(placename = Placename), 
                       join_by(Project, CameraIDName, UTM_X_meters, UTM_Y_meters)) %>% 
  mutate(Placename = if_else(Project == 'Basel Zoo PygmyHippo 2018-2020', placename, Placename)) %>% select(-placename)

# solve issues with ARTP data set 
pres_cam <- pres_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% 
  mutate(UTM_X_meters = ifelse(UTM_X_meters == 292289, 291289, UTM_X_meters)) %>% # most likely this is a typo - replaced 2 by one to match the deployments coordoninates 
  select(-notes, -CameraIDName, -Placename) %>% 
  left_join(cam_deploy %>% 
              filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% 
              select(Placename, CameraIDName, UTM_X_meters, UTM_Y_meters), 
            join_by(UTM_X_meters, UTM_Y_meters)) %>% 
  mutate(placename = Placename, cameraIDName = CameraIDName) %>% # create dummy columns which need to be merged with the other data frame pres cam - but only for 'REDD_ARTP_PygmyHippo 2013-2014'
  select(UTM_Y_meters, UTM_X_meters, placename, cameraIDName) %>% 
  distinct() %>% # remove all duplicates which resulted from replicated observations but case many to many relationship in join
  right_join(pres_cam %>%  mutate(UTM_X_meters = ifelse(UTM_X_meters == 292289, 291289, UTM_X_meters)), join_by(UTM_X_meters, UTM_Y_meters)) %>% 
  mutate(Placename = if_else(Project == 'REDD_ARTP_PygmyHippo 2013-2014', placename, Placename), 
         CameraIDName = if_else(Project == 'REDD_ARTP_PygmyHippo 2013-2014', cameraIDName, CameraIDName)) %>% 
  select(Project, Country, Placename, DeploymentID, CameraIDName,Obs_Date, time, count, UTM_X_meters, UTM_Y_meters, notes)

# solve issues in the Pygmy Hippo REDD CT 2019-2021 data set 

# transfer coordinates from deployment data to pres_cam data set for DArwin_morroRiver and IWT_CF data set
pres_cam <- pres_cam %>% left_join(cam_deploy %>% select(Project, Placename, UTM_X_meters, UTM_Y_meters) %>% 
                         filter(Project %in% c('Darwin_morroRiver', 'IWT_CF')) %>%
                         rename(utm_y_meters = UTM_Y_meters, utm_x_meters = UTM_X_meters), 
                       join_by(Project, Placename)) %>%
  mutate(UTM_X_meters = if_else(Project  %in% c('Darwin_morroRiver', 'IWT_CF'), utm_x_meters, UTM_X_meters), 
         UTM_Y_meters = if_else(Project  %in% c('Darwin_morroRiver', 'IWT_CF'), utm_y_meters, UTM_Y_meters)) %>% 
  select(-utm_x_meters, -utm_y_meters) %>% View()



ggplot() +
  geom_point(data = cam_deploy %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014', UTM_X_meters > 290000), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 4) +
  geom_point(data = pres_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014',  UTM_X_meters > 290000), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 2, color = 2)
ggplot() +
  geom_point(data = cam_deploy %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 2) +
  geom_point(data = pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), aes(y = UTM_Y_meters, x = UTM_X_meters), size = 1, color = 2)

# creation of a consistent, unique SiteID, since this is missing so far
cam_deploy <- cam_deploy %>% mutate(SiteID = case_when(Project == 'Basel Zoo PygmyHippo 2018-2020' ~ paste0(Placename, '_', CameraIDName), # Basel Zoo PH Project - River CameraIDName
                                         Project == 'REDD_ARTP_PygmyHippo 2013-2014' ~ paste0(Placename, '_', CameraIDName),  # placename and CameraIDName
                                         # Project == 'Pygmy Hippo REDD CT 2019-2021' ~ paste0(), # until now I did not manage to make sense of this data set
                                         Project == 'Darwin_morroRiver' ~ Placename, 
                                         Project == 'IWT_CF' ~ Placename)) 
pres_cam <- pres_cam %>% mutate(SiteID = case_when(Project == 'Basel Zoo PygmyHippo 2018-2020' ~ paste0(Placename, '_', CameraIDName), # Basel Zoo PH Project - River CameraIDName
                                       Project == 'REDD_ARTP_PygmyHippo 2013-2014' ~ paste0(Placename, '_', CameraIDName),  # placename and CameraIDName
                                       # Project == 'Pygmy Hippo REDD CT 2019-2021' ~ paste0(), # until now I did not manage to make sense of this data set
                                       Project == 'Darwin_morroRiver' ~ Placename, 
                                       Project == 'IWT_CF' ~ Placename))



# take the different data sources and try to identical link columns to each other for a join 

# Basel Zoo data
pres_cam %>% filter(Project == 'Basel Zoo PygmyHippo 2018-2020')
cam_deploy %>% filter(Project == 'Basel Zoo PygmyHippo 2018-2020') # its possible to join the data via CameraIDName or coordinates

# ARTP data - project names are different in both tables 
pres_cam %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% select(-notes)
cam_deploy %>% filter(Project == 'REDD_ARTP_PygmyHippo 2013-2014') # only join via coordinates possible 

# Pygmy Hippo REDD CT 2019-2021
pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% select(-notes) # the placename could be an ident representation of the site, but no equivalent in cam_deploy data 
cam_deploy %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') # connection via CameraIDName, but unfortunately coordinates may not work
cam_deploy %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021') %>% left_join(pres_cam %>% filter(Project == 'Pygmy Hippo REDD CT 2019-2021'), join_by(UTM_X_meters, UTM_Y_meters)) %>% View()
# DWCN23 has been in the field for only one day?

pres_cam %>% filter(Project == 'Darwin_morroRiver') %>% select(-notes) # no coordinates available 
cam_deploy %>% filter(Project == 'Darwin_morroRiver') # ident placename between both tables (placename and CamerIDName are the same in cam_deploy), deployment data for second deployment periods missing

pres_cam %>% filter(Project == 'IWT_CF') %>% select(-notes) #
cam_deploy %>% filter(Project == 'IWT_CF') # ident placenames, but no coordinates in presence data, lots of observations are outside the deployment period

# in darvinmorrow river - what is deploymentID?



# join data together 
presence <- data.frame(row.names = c('Project', 'Country', 'SiteID', 'CameraIDName', 'ObsDate', 'ObsTime', 'Count', 'UTM_X_meters', 'UTM_Y_metsers', 'Notes'))
pres_cam %>% filter(project == 'Basel Zoo PygmyHippo 2018-2020') %>% mutate(SiteID = paste0()) # placename 

deployments <- data.frame(row.names = c('Project', 'Country', 'River', 'SiteID', 'CameraIDName', 'Deployment', 'Collection', 'CameraFunctioning', 'UTM_X_meters', 'UTM_Y_metsers', 'Notes'))
cam_deploy %>% filter(project == 'Basel Zoo PygmyHippo 2018-2020')







unique(pres_cam$placename)

unique(cam_deploy$ProjectIDName)


# data from REDD_ARTP_PygmyHippo 2013-2014 cannot be matched with deployment data because no site IDs exist, probs a spatial join is needed 

# DarwinMorro river data and IWT_CF has no spatial location of deployments 

# visualisize the different deployment periods 
cam_deploy %>% 
  mutate(ID = paste0(PlotIDName, '_', CameraIDName)) %>%
  ggplot(aes(y = ID)) +
  geom_segment(aes(x = Deployment, xend = Collection, yend = ID), linewidth = 1, color = 'blue') 

# visualise the diffferent locations of various monitoring schemes
par(mfrow = c(2,1))
plot(data = cam_deploy, x = cam_deploy$UTM_X_meters, y = cam_deploy$UTM_Y_meters, add = T)
plot(data = pres_cam, x = pres_cam$x_coord, y = pres_cam$y_coord, col = 'red')


################################################################################
# THIS IS JUST RUBBISH STUFF


# Authenticate with SharePoint (First time will ask for browser login)
site - get_sharepoint_site(httpsrspb.sharepoint.comsitesGOLALIBRARY) # this requires Admin rights which is not doable for me

# Get the document library (where the file is stored)
doc_lib - site$get_drive()

# Open file stream from SharePoint
file - doc_lib$get_item(Shared DocumentsPygmy Hippopotamus MASTER 2008-2021 08May2023.xls)

# Read the Excel file directly into R
df - read_excel(file$download())

# View first rows
head(df)