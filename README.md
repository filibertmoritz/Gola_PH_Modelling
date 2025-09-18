# Pygmy hippo distribution modelling Gola, Sierra Leone

Pygmy hippos are elusive and rare forest mammals of conservation concern that occur in the Greater Gola Landscape, Sierra Leone. 
Unfortunately, limited knowledge exists on distribution of pygmy hippos across Gola, the drivers of occurence as well as changes in occurence over time. 

Since knowledge on drivers of pygmy hippo occupancy is rather limited, we apply a random forest model to narrow down possible predictors. Afterwards, we use these predictors in a static multi-season integrated occupancy model in spOcuupancy which integrates both camera trapping and riverine transect survey data to model Pygmy hippo distribution.

In more detail we apply the following steps: 

  **1. Monitoring data preparation:** First, prepare the monitoring data - this is organised in two separate scripts for [camera trapping](https://github.com/filibertmoritz/Gola_PH_Modelling/blob/main/scripts/Gola_PH_Modelling_data_prep_camera_traps_revised.R) and [riverine transect](https://github.com/filibertmoritz/Gola_PH_Modelling/blob/main/scripts/Gola_PH_Modelling_data_prep_transects_revised.R) data.
  **2. Environmental covariate preparation:** Second, extract environmental covariates from multiple different sources by using this [script](https://github.com/filibertmoritz/Gola_PH_Modelling/blob/main/scripts/Gola_PH_Modelling_data_prep_env_covariates_revised.R). Therefore, the study area can easily be adjusted by employing another .shp file as study area.
  **3. Find relevant environmental covariates:** Thrid, we fit a conditional random forest model on the pooled data to find relevant predictor variables and you can find this script [here](https://github.com/filibertmoritz/Gola_PH_Modelling/blob/main/scripts/Gola_PH_RF_variable_selection_period_weights.R).
  ***4. Fit occupancy model:** Lastly, we fit a static multi-season occupancy model using solely the most important environmental covariates, please find the script [here](https://github.com/filibertmoritz/Gola_PH_Modelling/blob/main/scripts/Gola_PH_Modelling_spOccupancy_tIntPGOcc_revised2.R).


