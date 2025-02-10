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
cam_deploy <- cam_deploy %>% mutate(project = factor(ProjectIDName), 
                      country = factor(country), 
                      River = factor(Riverlocation), 
                      PlotIDName = factor(PlotIDName),
                      placename = factor(placename), 
                      CameraIDName = factor(CameraIDName), 
                      Deployment = as.Date(DeploymentDate_dd.mm.yyyy, format = '%d/%m/%Y'), 
                      Collection = as.Date(CollectionDate_dd.mm.yyyy, format = '%d/%m/%Y'), 
                      Camera_Functioning = factor(if_else(Camera_Functioning == "Camera Functioning", 'functioning', 'failure')), 
                      habitat = factor(overallhabitat)) %>% 
  select(project, country, placename, Riverlocation, PlotIDName, CameraIDName, Deployment, Collection, Camera_Functioning, UTM_X_meters, UTM_Y_meters, habitat)

# format data types properly and get rid of uneeded columns - camera presence data 
str(pres_cam) # this looks rather scattered, thus do it separately for each data set 
pres_cam <- pres_cam %>% mutate(project = factor(project), 
                    country = factor(country), 
                    placename = factor(placename), 
                    deploymentID = factor(deploymentID), 
                    CameraIDName = factor(CameraIDName), 
                    Obs_Date = as.Date(date,format = '%d/%m/%Y'), 
                    UTM_X_meters = x_coord, UTM_Y_meters = y_coord, habitat = factor(habitat)) %>% 
  select(project, country, placename, deploymentID, CameraIDName, Obs_Date, time, UTM_X_meters, UTM_Y_meters, count, notes)

# take the different data sources and try to identical link columns to each other for a join 

# Basel Zoo data
pres_cam %>% filter(project == 'Basel Zoo PygmyHippo 2018-2020')
cam_deploy %>% filter(project == 'Basel Zoo PygmyHippo 2018-2020') # its possible to join the data via CameraIDName or coordinates

# ARTP data - project names are different in both tables 
pres_cam %>% filter(project == 'REDD_ARTP_PygmyHippo 2013-2014') %>% select(-notes)
cam_deploy %>% filter(project == 'REDD PygmyHippo 2013-2014') # only join via coordinates possible 






#















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