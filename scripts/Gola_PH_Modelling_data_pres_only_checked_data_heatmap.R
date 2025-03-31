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

# load in data - all 3 different data sets 
excel_sheets('data/combined_PH_data2025_draft7.xlsx')
pres_opp <- read_excel(path = 'data/combined_PH_data2025_draft7.xlsx', sheet = 'opportunistic_data')
pres_transects <- read_excel(path = 'data/combined_PH_data2025_draft5.xlsx', sheet = 'combined_transect_data')
pres_cam <- read_excel(path = 'data/combined_PH_data2025_draft7.xlsx', sheet = 'combined_camera_data')

# prepare presences from transect data 
pres_transects <- pres_transects %>% 
  mutate(Project = as.factor(if_else(DatasetName == 'Pygmy Hippo REDD 2013-2014', 'PH_ARTP_transects', 'PH_REDD_transects')), 
         Obs_Date = date(VisitDate_dd_mm_yyyy), 
         Sign = as.factor(Observation_Category_3), 
         River = as.factor(River_Stream_Name), 
         Obs_Method = factor('Transect Survey'))  %>% 
  rename(UTM_X_meters = UTM_X_m, UTM_Y_meters = UTM_Y_m) %>% 
  select(Project, River, Obs_Date, Sign, Obs_Method, UTM_X_meters, UTM_Y_meters)

# prepare presences from camera data 
pres_cam <- pres_cam %>% mutate(Project = factor(project), 
                                Obs_Date = as.Date(date,format = '%d/%m/%Y'), ,
                                UTM_X_meters = x_coord, UTM_Y_meters = y_coord, 
                                Obs_Method = factor('Camera Trap'), 
                                River = NA, Sign = NA) %>%
  select(Project, River, Obs_Date, Sign, Obs_Method, UTM_X_meters, UTM_Y_meters) 

# prepare opportunistic presences 
pres_opp <- pres_opp %>% 
  mutate(Project = factor(DatasetName), Country = factor(country), 
         Obs_Date = case_when(grepl('/', VisitDate_dd_mm_yyyy) == F ~ as.Date(as.numeric(VisitDate_dd_mm_yyyy), origin = "1899-12-30"),
                              grepl('/', VisitDate_dd_mm_yyyy) == T ~ as.Date(VisitDate_dd_mm_yyyy, format = '%d/%m/%Y')), 
         Sign = factor(Observation_Category_3), 
         River = factor(str_replace(River_Stream_Name, pattern = 'Unknown', replacement = NA_character_)), 
         Obs_Method = factor('Opportunistic')) %>% 
  rename(UTM_X_meters = UTM_X_m, UTM_Y_meters = UTM_Y_m) %>% 
  select(Project, River, Obs_Date, Sign, Obs_Method, UTM_X_meters, UTM_Y_meters)

pres_all_sf <- bind_rows(pres_cam, pres_transects, pres_opp) %>% 
  st_as_sf(coords = c('UTM_X_meters', 'UTM_Y_meters'), remove = F, crs = 32629) 

plot(st_geometry(pres_all_sf))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Create a raster through the whole study area #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create variables that are needed to split up the time period
pres_all_sf <- pres_all_sf %>% 
  mutate(Period = factor(if_else(Obs_Date < as.Date('2015-12-31'), '2008-2015', '2016-2025')), 
         Obs_Year = year(Obs_Date), 
         Obs_Week = week(Obs_Date)) %>% 
  group_by(Project, Obs_Date) %>% 
  mutate(dist_to_nearest = sapply(1:n(), function(i) {
    dists <- as.numeric(st_distance(geometry[i], geometry))  # pairwise distance
    dists[dists == 0] <- NA  # remove dist to itself
    min(dists, na.rm = TRUE)}), # get min dist
  dist_to_nearest = if_else(dist_to_nearest == 'Inf', NA, dist_to_nearest)) # handly Inf, which results from no other point to calc distance to

# apply filter to restrict violations of independence of observations
pres_ind_sf <- pres_all_sf %>% group_by(Project, Obs_Date) %>% 
  arrange(Project, Obs_Date) %>% 
  mutate(first = if_else(row_number()==1, 1, NA), # this flags one observation per project and day
         dist = if_else(dist_to_nearest > 200, 1, NA)) %>% # this flags all observations within the same project and day that are 200m apart from each other
  group_by(Project, Obs_Year, Obs_Week, UTM_X_meters, UTM_Y_meters) %>%
         mutate(camera = if_else(row_number() == 1, 1, NA)) %>% # this flags only one observation within the same project, year and week at the same location - only for camera traps!
  filter(first == 1 | dist == 1, Obs_Method != 'Camera Trap' | camera == 1) %>%# filter out observations that do not have a single flag, last condition only applicable for camera trap data 
  select(-first, -dist, -camera) 


# create a grid to aggregate the presences to these grid cells 
pres_grid <- st_make_grid(pres_all_sf, cellsize = 1500, square = F) %>% st_as_sf() # grid that captures the whole study area

# join data frames and count numbers of points within polygon distinct time periods
pres_per_cell_all <- st_join(pres_all_sf, pres_grid, join = st_within) %>% group_by(geometry, Period) %>% summarise(n_obs = n(), .groups = 'keep') # for all presences
pres_per_cell_ind <- st_join(pres_ind_sf, pres_grid, join = st_within) %>% group_by(geometry, Period) %>% summarise(n_obs = n(), .groups = 'keep') # for independent presences 


# join data frames and count numbers of points within polygon yearly for past 2018 data 
pres_per_cell_all_2018 <- st_join(pres_all_sf %>% filter(Obs_Year > 2017), pres_grid, join = st_within) %>% group_by(geometry, Obs_Year) %>% summarise(n_obs = n(), .groups = 'keep') # for all presences
pres_per_cell_ind_2018 <- st_join(pres_ind_sf %>% filter(Obs_Year > 2017), pres_grid, join = st_within) %>% group_by(geometry, Obs_Year) %>% summarise(n_obs = n(), .groups = 'keep') # for independent presences 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 4. Plot all the data for distinct time periods #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot the presences that are assumened to be independent 
hotspot_map_ind <- ggplot() +
  annotation_map_tile(zoom = 10, type = 'cartolight') +  # change map by argument 'type' - osm, cartolight, hotstyle
  geom_point(data = pres_per_cell_ind, aes(x = st_coordinates(geometry)[,1], 
                                     y = st_coordinates(geometry)[,2], 
                                     size = n_obs, color = Period), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +  # Adjust bubble sizes
  scale_color_manual(values = c("2008-2015" = "blue", "2016-2025" = "red")) +
  theme_minimal() +
  labs(title = "Temporal Hotspot Map of Pygmy Hippo Observations in Gola",
       subtitle = 'Filtered for independent observations',
       size = paste0("Number of Observations\n(Overall: ", pres_per_cell_ind %>% pull(n_obs) %>% sum(), ")"), 
       x = 'Longitude', y = 'Latutude' ) +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1), # add frame
        axis.line = element_blank())
ggsave(filename = 'output/plots/PH_hotspot_map_ind.jpg', plot = hotspot_map_ind, height = 6, width = 12)

# plot just all presences without filtering for indepencence
hotspot_map_all <- ggplot() +
  annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_point(data = pres_per_cell_all, aes(x = st_coordinates(geometry)[,1], 
                                           y = st_coordinates(geometry)[,2], 
                                           size = n_obs, color = Period), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_manual(values = c("2008-2015" = "blue", "2016-2025" = "red")) +
  theme_minimal() +
  labs(title = "Temporal Hotspot Map of Pygmy Hippo Observations in Gola",
       subtitle = 'All presences without filtering for independent observations are shown',
       size = paste0("Number of Observations\n(Overall: ", pres_per_cell_all %>% pull(n_obs) %>% sum(), ")"), 
       color = "Time Period", 
       x = 'Longitude', y = 'Latitude') +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(filename = 'output/plots/PH_hotspot_map_all.jpg', plot = hotspot_map_ind, height = 6, width = 12)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Plot all the data for yearly data  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a plot that shows unfiltered data (possibly not independent) from 2018 onwards with indication of the year (as.numeric)
hotspot_map_all_yearly <- ggplot() +
  annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_point(data = pres_per_cell_all_2018, aes(x = st_coordinates(geometry)[,1], 
                                           y = st_coordinates(geometry)[,2], 
                                           size = n_obs, color = Obs_Year), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "Temporal Hotspot Map of Pygmy Hippo Observations in Gola after 2018",
       subtitle = 'All presences without filtering for independent observations are shown',
       size = paste0("Number of Observations\n(Overall: ", pres_per_cell_all_2018 %>% pull(n_obs) %>% sum(), ")"), 
       color = "Time Period", 
       x = 'Longitude', y = 'Latitude') +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(filename = 'output/plots/PH_hotspot_map_all_yearly.jpg', plot = hotspot_map_all_yearly, height = 6, width = 12)


# create a plot with filtered, independent observations from 2018 onwards with indication of the year (as.numeric)
hotspot_map_ind_yearly <- ggplot() +
  annotation_map_tile(zoom = 10, type = 'cartolight') +
  geom_point(data = pres_per_cell_ind_2018, aes(x = st_coordinates(geometry)[,1], 
                                                y = st_coordinates(geometry)[,2], 
                                                size = n_obs, color = Obs_Year), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_minimal() +
  labs(title = "Temporal Hotspot Map of Pygmy Hippo Observations in Gola after 2018",
       subtitle = 'Filtered for independent observations',
       size = paste0("Number of Observations\n(Overall: ", pres_per_cell_ind_2018 %>% pull(n_obs) %>% sum(), ")"), 
       color = "Time Period", 
       x = 'Longitude', y = 'Latitude') +
  theme(legend.position = "right", 
        panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(filename = 'output/plots/PH_hotspot_map_ind_yearly.jpg', plot = hotspot_map_ind_yearly, height = 6, width = 12)
