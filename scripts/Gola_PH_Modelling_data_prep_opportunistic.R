#### data preparation for transects script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations #######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install needed packages 
#install.packages('readxl')
#install.packages('Microsoft365R')

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(stringr)
library(mapview)
filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in table data and tidy everything up #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load in data 
excel_sheets('data/combined_PH_data2025_draft5.xlsx')
pres_opp <- read_excel(path = 'data/combined_PH_data2025_draft5.xlsx', sheet = 'opportunistic_data')

# tidy everything a bit up, change data type in columns and remove unneeded columns  
pres_opp <- pres_opp %>% 
  mutate(Project = factor(DatasetName), Country = factor(country), 
         Obs_Date = case_when(grepl('/', VisitDate_dd_mm_yyyy) == F ~ as.Date(as.numeric(VisitDate_dd_mm_yyyy), origin = "1899-12-30"),
                                  grepl('/', VisitDate_dd_mm_yyyy) == T ~ as.Date(VisitDate_dd_mm_yyyy, format = '%d/%m/%Y')), 
         Obs_Time = as_hms(VisitTime_hh_mm_ss), 
         Obs_DateTime = as.POSIXct(if_else(Project == 'Pygmy Hippo REDD 2013-2014', 
                                           paste(Obs_Date, '00:00:00'), 
                                           paste(Obs_Date, Obs_Time)), 
                                   format = "%Y-%m-%d %H:%M:%OS"), 
         Sign = factor(Observation_Category_3), River = factor(River_Stream_Name), Weather = factor(weather)) %>% 
  select(Project, Country, Obs_DateTime, UTM_X_m, UTM_Y_m, Sign, River, Weather, Comments)

# replace na and unknown by NA
pres_opp <- pres_opp %>% mutate(River = str_replace(River, pattern = 'Unknown', replacement = NA_character_), 
                    Weather = str_replace(Weather, pattern = 'na', replacement = NA_character_), 
                    Sign = str_replace(Sign, pattern = 'unknown', replacement = NA_character_)) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. create geometry for locations #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a sf object
pres_opp <- pres_opp %>% 
  st_as_sf(coords = c('UTM_X_m', 'UTM_Y_m'), crs = 32629, remove = F) # transform into sf object with WGS 84 / UTM zone 29N (32629)

#plot data 
mapview(pres_opp, popup = T, zcol = 'Sign', legend = T) 
map <- mapview(pres_opp, popup = T, zcol = 'Sign', legend = T) 

# save prepared data without geometry
write.csv(pres_opp %>% st_drop_geometry(), file = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/data/PH_prepared_pres_opp_data.csv')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. create geometry for locations #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pres_all <- pres_opp %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Obs_DateTime)) +
  theme_bw()  +
  labs(title = "Number of Opportunistic Presences Across Time", subtitle = "Infinite number of signs allowed per date", x = "Time", y = "Number of Presences")
pres_filt <- pres_opp %>% mutate(Date = as.Date(Obs_DateTime)) %>% group_by(Project, Date) %>% slice(1) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = Obs_DateTime))+
  theme_bw()  +
  labs(title = "Number of Opportunistic Presences Across Time", subtitle = "Only one presence per date allowed", x = "Time", y = "Number of Presences")

# arrange subplots 
library(ggpubr)
overview <- ggarrange(pres_all, pres_filt)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/Data_overview_opportunistic.jpg', plot = overview, width = 8, height = 6)
