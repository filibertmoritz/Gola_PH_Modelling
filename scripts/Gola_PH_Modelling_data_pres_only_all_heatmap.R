#### data preparation for transects script for pygmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025

# the aim of this script is to create a hotspot map, where and when observations of Pygmy Hippos occured to inform upcoming monitoring 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations #######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(hms)
library(readxl)
library(tidyverse)
library(lubridate)
library(sf)
library(stringr)
library(mapview)
library(ggmap)
library(ggspatial)
library('prettymapr')

filter <- dplyr::filter
select <- dplyr::select
rename <- dplyr::rename 

# set working directory
# setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. Read in table data and tidy everything up #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load in data - both PH Presence Master Table and data prepared from different data sets
excel_sheets('data/combined_PH_data2025_draft7.xlsx')
pres_opp <- read_excel(path = 'data/combined_PH_data2025_draft7.xlsx', sheet = 'opportunistic_data')

excel_sheets('data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx')
pres_master <- read_excel(path = 'data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx', sheet = 'PYGMY HIPPO MASTER')

# tidy everything a bit up, change data type in columns and remove unneeded columns  
pres_master <- pres_master %>% 
  mutate(Project = factor(DatasetName),
         Obs_Date = as.Date(VisitDate_dd_mm_yyyy), River = factor(River_Stream_Name), 
         Sign = factor(Observation_Category_3), ObservationMethod = factor(ObservationMethod)) %>% 
  select(Project, Obs_Date, UTM_X_m, UTM_Y_m, River, ObservationMethod, Sign)

# create a spatial object 
pres_master_sf <- pres_master %>% 
  st_as_sf(coords = c('UTM_X_m', 'UTM_Y_m'), crs = 32629, remove = F) # transform into sf object with the correct crs

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Create a raster through the whole study area #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create variables that are needed to split up the time period
pres_master_sf <- pres_master_sf %>% 
  mutate(Period = factor(if_else(Obs_Date < as.Date('2015-12-31'), '2008-2015', '2016-2025'))) %>% 
  group_by(Project, Obs_Date) %>% 
  mutate(dist_to_nearest = sapply(1:n(), function(i) {
    dists <- as.numeric(st_distance(geometry[i], geometry))  # pairwise distance
    dists[dists == 0] <- NA  # remove dist to itself
    min(dists, na.rm = TRUE)}), # get min dist
  dist_to_nearest = if_else(dist_to_nearest == 'Inf', NA, dist_to_nearest))

# create a grid to aggregate the presences to these grid cells 
pres_grid <- st_make_grid(pres_master_sf, cellsize = 1500, square = F) %>% st_as_sf() # grid that captures the whole study area

# join data frames and count numbers of points within polygon
pres_per_cell <- st_join(pres_master_sf, pres_grid, join = st_within) %>% group_by(geometry, Period) %>% summarise(n_obs = n(), .groups = 'keep')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Plot all the data #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

hotspot_map <- ggplot() +
  annotation_map_tile(zoom = 10, type = 'cartolight') +  # change map by argument 'type' - osm, cartolight, hotstyle
  geom_point(data = pres_per_cell, aes(x = st_coordinates(geometry)[,1], 
                                     y = st_coordinates(geometry)[,2], 
                                     size = n_obs, color = Period), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +  # Adjust bubble sizes
  scale_color_manual(values = c("2008-2015" = "blue", "2016-2025" = "red")) +
  theme_minimal() +
  labs(title = "Temporal Hotspot Map of Pygmy Hippo Observations in Gola",
       size = "Number of Observations", color = "Time Period", 
       x = 'Longitude', y = 'Latutude' ) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), # add frame
        axis.line = element_blank())
ggsave(filename = 'output/plots/PH_hotspot_map_all_data.jpg', plot = hotspot_map, height = 6, width = 12)


