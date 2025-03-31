#### pygmy hippo random forest sdm to explore important variables for integrated occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in March 2025 
#### Thanks for hints and comments from Felicity Edwards and Steffen Oppel

#### Background: We have data from various PH surveys which can be mainly distinguished in camera trap surveys, transect surveys and opportunistic presence-only data 
#### Issue: The data is rather scattered and we are unsure, if its enough to fit a more sophisticated, integrated model 
#### Solution: take all pooled data and explore, a) if we can fit a more basic random forest model (to afterwards, possibly go on to the integrated model) and b) explore which predictors could be relevant for the integrated model 
#### Target: Built a RF model with all pooled data - here, I will simply use the Pygmy hippo master table 
#### Steps: a) load pygmy hippo master table , b) load elsewhere prepared environmental covariates, 
####        c) fit RF model, predict/plot distribution of PH and d) explore the effect of the different variables 


##### THINGS TO IMPROVE 
# read through the random forest model literature and Steffen's scripts to check that I do not make any silly mistakes 
# add plot that explore the effect of all the variables 

##### BACKGROUND READING 
# https://christophm.github.io/interpretable-ml-book/

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(sf)
library(terra)
library(tidyverse)
library(lubridate)
library(hms)
library(data.table)
library(ranger)
library(scales)
library(tmap)
library(units)
library(readxl)

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

getwd()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 2. load and prepara data #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# presence data 
excel_sheets('data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx')
pres_master <- read_excel(path = 'data/Pygmy Hippopotamus MASTER 2008-2021 08May2023.xlsx', sheet = 'PYGMY HIPPO MASTER')

# tidy everything a bit up, change data type in columns and remove unneeded columns  
pres_master <- pres_master %>% 
  mutate(Project = factor(DatasetName),
         Obs_Date = as.Date(VisitDate_dd_mm_yyyy), River = factor(River_Stream_Name), 
         Sign = factor(Observation_Category_3), ObservationMethod = factor(ObservationMethod)) %>% 
  select(Project, Obs_Date, UTM_X_m, UTM_Y_m, River, ObservationMethod, Sign)

# create a spatial frame
pres_master_sf <- pres_master %>% 
  st_as_sf(coords = c('UTM_X_m', 'UTM_Y_m'), crs = 32629, remove = F) 


# grid and environmental covariates 
grid_sf <- st_read('data/PH_grid.shp') %>% st_as_sf() # laod shp and transform to sf 
envCovs <- fread('data/PH_prepared_env_covariates.csv') %>% 
  select(-V1) %>% as_tibble()
envCovs_sf <- grid_sf %>% left_join(envCovs, join_by(CellID)) %>% st_as_sf()# transfer geometry from grid to envCovs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 3. Create presence-absence information per grid cell and bring data into one big data frame #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# calc occurence per grid cell
occu <- pres_master_sf %>% select(geometry) %>%
  st_join(grid_sf) %>% 
  st_drop_geometry() %>%
  group_by(CellID) %>% 
  mutate(Occu = if_else(n() > 0, 1, 0)) %>% # this transferres all count data into presence-absence data 
  distinct() %>% ungroup()

# create one massive sf with envCovs, grid cells and occu data 
data_sf <- envCovs_sf %>% left_join(occu, join_by(CellID)) # transfer to full grid 

data_sf$Occu[is.na(data_sf$Occu)] <- 0 # this replaces all NAs to 0

plot(data_sf %>% select(Occu))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 5. Bring all data together #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data overview
data_sf

# create a training data set where either a presence was observed or surveys were performed!
traindata <- data_sf %>% st_drop_geometry()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Look at correlation of features/predictors #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# correlation
library(corrplot)
data_sf %>% ungroup() %>%
  st_drop_geometry() %>% 
  select(river_density, river_density_med_large, Distance_campsite, Distance_settlement, 
         Distance_large_river, bareground, regrowth1, farm_bush, forest, water, swamp, cocoa, 
         oil_palm, mean_elev, NDVI_2013, NDVI_2020, EVI_2013, EVI_2020, Occu) %>%  
  cor() %>% 
  corrplot()

# pca
pc <- prcomp(data_sf %>% 
               st_drop_geometry() %>%
               select(river_density, river_density_med_large, Distance_campsite, Distance_settlement, 
                      Distance_large_river, bareground, regrowth1, farm_bush, forest, water, swamp, 
                      cocoa, oil_palm, mean_elev, NDVI_2013, NDVI_2020, EVI_2013, EVI_2020),
             center = TRUE,
             scale. = TRUE)

# biplot
library(factoextra)
fviz_pca_biplot(pc, label = "var", habillage = data_sf$Occu, addEllipses = TRUE)
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_data_exploration_pca.jpg', plot = last_plot())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 6. Train random forest model #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Things to consider: 
#     a)  Craig and Huettmann 2009 used binary logistic regression in TreeNet and balanced option under class weights, which automatically reweighted each class to 
#         account unequal sample size between presence and absence points
#     b)  validate predictions by a confusion matrix which report precentages of presences and absences correctly classified (cutoff at 0.5)
#     c)  Use effort as weights? 
#     d)  Valavi et al. 2023: class overlap (presence/absence is not clearly distinguishable by one variable - common for changes in occupancy across ecological gradients and rare species)
#         - solutions in paper, also example code


# 1) Incorporate effort in the response
rf1 <- ranger(Occu ~ river_length_med_large+Distance_settlement+Distance_large_river+river_length+forest+swamp+bareground+regrowth1+water+farm_bush+cocoa+oil_palm+mean_elev+NDVI_2020+EVI_2020, 
              data = traindata, num.trees = 1000, 
              #mtry = T, 
              importance = 'permutation', # choose one of: "impurity" (Based on Gini impurity reduction), "impurity_corrected" (bias corrected), "permutation" (Based on permuting feature values and measuring performance drop)
              sample.fraction = 0.5) # perform down-sampling, as suggested in Valavi et al. 2022
summary(rf1)
rf1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 7. Get Impression of Importance of the environmental covariates #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reading on partial dependence plots
# https://www.rpubs.com/vishal1310/QuickIntroductiontoPartialDependencePlots

# alternatively, use vivid package with nice tutorials: https://alaninglis.github.io/vivid/articles/rangerVivid.html
# however, a few computations are enormously computation intensive!

# get variable importance plots
library(vip)
vip(rf1, all_permutations = T, include_type = T) + labs(title = 'Variable Importance in RF1') + theme_bw()
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_vip_rf1.jpg', plot = last_plot())

# get most imortant variables 
variables <- names(sort(importance(rf1), decreasing = TRUE))

# plot partial predictive plots, first using vivid
library(vivid)
gc() # its rather computation intensive
viv_plot1 <- vivid::pdpPairs(data = traindata, 
                             fit = rf1, 
                             response = 'Occu', 
                             nmax = 150, # samples rows from data to display, NULL means all rows/observations, 50 is good choice for computation time 
                             gridSize = 15, # huge influence on computation time, set to 20 for higher resolution in pixel-plot
                             nIce = 100, # Number of ice curves to be plotted, defaults to 30.
                             vars = variables) # choose important variables here 
ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_feature_influence_all_presences.jpg', plot = viv_plot1, width = 20, height = 16)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 8. Predict with all 3 random forest models #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# predict and store data in data frame 
data_sf$pred1 <- predict(rf1, data = data_sf %>% st_drop_geometry())$predictions 

# confusion matrix - doesn't work yet
# library(caret)
# caret::confusionMatrix(data$pred1, reference = data$Occu)

# plot
plot(data_sf %>% select(pred1, Occu))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 9. Plot predictions  #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# plot option 1
opt1 <- ggplot(data_sf) +
  geom_sf(aes(fill = pred1), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Distributon in Gola, Sierra Leone",
       subtitle = 'Predictions from a random forest SDM with survey effort included in response variable', 
       fill = "Prediction", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()

# plot raw observations
raw<-ggplot(data_sf) +
  geom_sf(aes(fill = Occu), alpha = 1) +
  scale_fill_viridis_c(option = 'plasma', direction = -1) +  # Use a color gradient (or use scale_fill_gradient())
  labs(title = "Pygmy Hippo Occupancy in Gola, Sierra Leone",
       subtitle = 'Raw observations aggregated per grid cell', 
       fill = "Occurence", 
       x = 'Longitude', y = 'Latitude') +
  theme_bw()

# create an arranged plot which compares the different methods 
library(ggpubr)
ggarrange(opt1,raw)
# ggsave(filename = 'C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/output/plots/PH_RF_predictions.jpg', plot = last_plot(), width = 15, height = 10)


