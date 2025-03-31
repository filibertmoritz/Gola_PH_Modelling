#### data preparation for transects script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in Feb 2025


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(maptiles)
library(stringr)

library(mapview)
library(leaflet)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in table data and tidy everything up #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in all data sets 
excel_sheets('data/combined_PH_data2025_draft7.xlsx')
pres_transects <- read_excel(path = 'data/combined_PH_data2025_draft7.xlsx', sheet = 'combined_transect_data')
locs_transects <- read_excel(path = 'data/combined_PH_data2025_draft7.xlsx', sheet = 'locations_transects') %>% 
  drop_na(DatasetName) # there are lots of NA/empty fields in the table tail

# format data types properly and get rid of unneeded columns - transect location data 
str(locs_transects)
locs_transects <- locs_transects %>% 
  mutate(Project = as.factor(DatasetName), 
         River = as.factor(River_Stream_Name), 
         TransectIDName = as.factor(TransectIDName), 
         uniqueID = as.factor(uniqueID), ### attention, there might be a typo with uniquEID
         Season = str_to_title(Season), 
         start_date = if_else(Project == 'Pygmy Hippo ARTP_REDD 2013-2014', paste(start_date, '00:00:00'), start_date), 
         end_date = if_else(Project == 'Pygmy Hippo ARTP_REDD 2013-2014', paste(end_date, '00:00:00'), end_date), 
         start_date = as.POSIXct(start_date, format = "%d.%m.%Y %H:%M:%OS"), 
         end_date = as.POSIXct(paste0(end_date), format = "%d.%m.%Y %H:%M:%OS")) %>% 
  rename(Start_UTM_X_meters = Start_UTM_X_m, Start_UTM_Y_meters = Start_UTM_Y_m, End_UTM_X_meters = End_UTM_X_m, End_UTM_Y_meters = End_UTM_Y_m) %>% 
  select(Project, TransectIDName, uniqueID, River, Start_UTM_X_meters, Start_UTM_Y_meters, End_UTM_X_meters, End_UTM_Y_meters, start_date, end_date) 
locs_transects$River[locs_transects$River=='Unknown'] <- NA # replace weird na's

# format data types properly and get rid of unneeded columns - transect presence data 
str(pres_transects)
pres_transects <- pres_transects %>% 
  mutate(Project = as.factor(DatasetName), 
         Obs_DateTime = as.POSIXct(paste0(date(VisitDate_dd_mm_yyyy), as_hms(VisitTime_hh_mm_ss)), format = '%Y-%m-%d %H:%M:%OS'), 
         Obs_DateTime = if_else(is.na(Obs_DateTime), as.POSIXct(paste0(date(VisitDate_dd_mm_yyyy), ' 00:00:00'), format = '%Y-%m-%d %H:%M:%OS'), Obs_DateTime), # if there is no time, then just always take 00:00
         Sign = as.factor(Observation_Category_3), 
         River = as.factor(River_Stream_Name),
         uniqueID = factor(uniqueID),
         TransectIDName = as.factor(REDD_SectionIDStreamName)) %>% 
  rename(UTM_X_meters = UTM_X_m, UTM_Y_meters = UTM_Y_m, Distance_River_m = Distance_to_nearest_stream_m, River_Width = Estimated_stream_width_m) %>% 
  select(Project, TransectIDName, uniqueID, River, Obs_DateTime, Sign, UTM_X_meters, UTM_Y_meters)

# create sf for presence points 
pres_transects <- pres_transects %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Read in spatial and tidy everything up #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in spatial transct data 
paths_artp <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/PH_ARTP_REDD_2013_2014_transects.shp') %>% st_as_sf()
paths_redd <- st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/REDD_Pygmy_Hippo_Streams_Transects_UPDATED_2.shp') %>% st_as_sf()

# plot paths
plot(st_geometry(paths_artp), col = 'forestgreen')
plot(st_geometry(paths_redd), add = T, col = 'red')

st_crs(paths_artp) == st_crs(paths_redd)


# get additional data from artp transects shp, especially geometries and survey times 
print(locs_transects %>% filter(Project == 'Pygmy Hippo ARTP_REDD 2013-2014'), n = 35)
head(paths_artp)

paths_artp <- paths_artp %>% # first part only to tidy up the sf
  mutate(Project = 'Pygmy Hippo ARTP_REDD 2013-2014', 
         DateTime_Start = as.POSIXct(paste0(start_date, Start_time), format = "%d.%m.%Y %H:%M:%OS"), 
         DateTime_End = as.POSIXct(paste0(end_date, End_time), format = "%d.%m.%Y %H:%M:%OS"), 
         uniqueID = as.factor(uniqueID)) %>% 
  select(TransectID, uniqueID, DateTime_Start, DateTime_End, geometry) %>% st_as_sf() # only select those columns which contains data to complement the locs_transects data 

locs_transects <- locs_transects %>% # transfer data from shp to overall transect sf
  left_join(paths_artp %>% select(-TransectID), join_by(uniqueID)) %>%
  mutate(DateTime_Start = if_else(Project == 'Pygmy Hippo ARTP_REDD 2013-2014', DateTime_Start, start_date),
         DateTime_End = if_else(Project == 'Pygmy Hippo ARTP_REDD 2013-2014', DateTime_End, end_date)) %>% 
  select(-start_date, -end_date) 

# get additional data from redd transects shp, especially geometries 
print(locs_transects %>% filter(Project == 'Pygmy Hippo REDD 2019-2021'), n = 13)
paths_redd <- paths_redd %>% 
  rename(TransectIDName = StreamIDNa, River = StreamName, Start_UTM_X_m = StreamStar, Start_UTM_Y_m = StreamSt_1, End_UTM_X_m = StreamEnd_,  End_UTM_Y_m = StreamEnd1) %>% 
  select(-StreamID, -Shape__Len, -StreamEn_2, -StreamSt_4, -StreamSt_3, -StreamEn_3, -StreamSt_2, -StreamEn_1)

locs_transects <- locs_transects %>% 
  left_join(paths_redd %>% rename(geometry_redd = geometry) %>% select(TransectIDName, geometry_redd), join_by(TransectIDName)) %>% 
  mutate(geometry = if_else(Project == 'Pygmy Hippo REDD 2019-2021', geometry_redd, geometry)) %>% 
  select(-geometry_redd) %>% 
  st_as_sf(crs = 32629) %>% 
  st_cast('MULTILINESTRING') %>%
  mutate(geometry = st_zm(geometry, drop = T))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Brief data overview ###########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
# check head and crs
head(locs_transects) # locations on all transects
head(pres_transects)

# plot
plot(st_geometry(locs_transects), col = locs_transects$Project)
plot(st_geometry(pres_transects), col = locs_transects$Project, add = T)

tmap_mode("view")  # interactive mode
tm_shape(locs_transects) + 
  tm_lines(col = "Project", lwd = 2) +
tm_shape(pres_transects) + 
  tm_dots(col = "Project", size = 1, alpha = 1) 


#################################################################################
#### I LEFT OUT ALL THE DATA REMOVING STUFF, SINCE I ASSUME THAT FELICITY DID THIS ALREADY ###
################################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Delete weird data as agreed for with Felicity ###########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# remove transects which do not make sense - however, ensure that no presences are connected to these!

# Filiberts suggestions

# Number 15 and 14 as well as 13 and 12 - are they each just duplicates - solve by deleting each one?
# Patch of Number 16, 17,18, 20 - which ones are relevant? - possibly all of them and we could treat them as replicates within one grid cell in different visits (either season or year) and correct for det by length through a transect.
# Number 24 and 25 - are the transects correct in the way I followed the river - other rotes would be possible?
# Number 38: I just made this connection up in case there might be a river which is not included in the river shape used
# Number 53: this transect is very short (also number 47) and not close to a river, how to proceed with this one?

# these are the things we agreed for with Felicity

# The two very short transects, remove these.
# The other transects make sense to me.
# For those transects which link together or slightly cross over we perhaps need to group them under a 
### common river name and have the total length of river that way, and then if there are overlaps then this 
### isn’t counted twice, unless surveying was in different years. 

# one could remove these transects, however, to me this makes slightly sense and I think the I would only delete this one transect 

# removal <- c('Weisei') # TransectIDNames to remove, checked, that we don't loose presences 

# locs_transects <- locs_transects %>% filter(!TransectIDName %in% removal)
# surveys_transects <- surveys_transects %>% filter(!TransectIDName %in% removal)
# paths_artp <-  paths_artp %>% filter(!TransectIDName %in% removal)
# paths_redd <- paths_redd %>% filter(!TransectIDName %in% removal)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Export data  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

st_write(locs_transects %>% select(uniqueID, geometry), 
         dsn = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_transect_paths.shp', 
         append = F)
write.csv(locs_transects %>% st_drop_geometry(), 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_transect_paths.csv')
write.csv(pres_transects %>% st_drop_geometry(), 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_pres_transect_data.csv')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Create a few overview plots ##### 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pres_all <- pres_transects %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Obs_DateTime)) +
  theme_bw()  +
  labs(title = "Number of Presences Across Time", subtitle = "Infinite number of signs allowed per transect & date", x = "Time", y = "Number of Presences")

pres_filt <- pres_transects %>% 
  mutate(Date = as.Date(Obs_DateTime)) %>% 
  group_by(TransectIDName, Date) %>% 
  slice(1) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Date)) +
  theme_bw()  +
  labs(title = "Number of Presences Across Time", subtitle = "Only one presence per transect and date allowed", x = "Time", y = "Number of Presences")

survey <- locs_transects %>%
  ggplot() +
  geom_histogram(mapping = aes(x = DateTime_Start))  +
  theme_bw() +
  labs(title = "Number of Transect Surveys Across Time",x = "Time", y = "Number of Surveys")

# create arranged plot and save
library(ggpubr)
overview <- ggarrange(pres_all, pres_filt, survey, ncol = 2, nrow = 2, common.legend = F)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/Data_overview_transects.jpg', plot = overview, width = 8, height = 6)

