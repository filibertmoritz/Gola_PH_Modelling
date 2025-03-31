#### script to get data from .gdb files that are used in ArcGis but not easilly accessable in qGIS #####
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de, in Feb 2025

#### Problem: Data in .gdb file which is not that easily accessable, thus, load data in and store it again in my files in as usable .shp or .tif 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1. load packages and set working directoty ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(sf) # for vector data 
library(terra)
library(dplyr)
select <- dplyr::select
rename <- dplyr::rename
select <- dplyr::select

setwd('C:/Users/filib/Documents/Praktika/RSPB/Gola_PH_Modelling/')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. get elevation data ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# list layers in data base 
st_read(dsn = 'C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', layer = 'ASTER_Contour10m')

# try to load data gives short list with file names - thats the best I can get
# rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'xxx') # produce a warning with all the raster data available inside 

# load data in 
elev <- list()
elev$Gola30mASTER <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'Gola30mASTER')
elev$ASTER30SLOPE <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'ASTER30SLOPE') # gives the slope for each grid cell within a DEM
elev$ASTER30HILLSHADE <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'ASTER30HILLSHADE') # this seems to give hillshades for better visualisation 
elev$Gola90mDTM <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'Gola90mDTM') # this might be the DTM/DEM which contains elevation 
elev$SRTM90m <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'SRTM90m')
elev$ASTER30ASPECT <- rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Elevation.gdb/Elevation.gdb/', 'ASTER30ASPECT')

# plot data for overview
plot(elev[[4]])

# write data to files 
# writeRaster(elev$Gola30mASTER, filename = 'C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/raster data/Gola30mASTER.tif')

for(i in 1:length(elev)){ # save files in a loop
  writeRaster(elev[[i]], filename = paste0('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/raster data/', names(elev)[i], '.tif'), overwrite = T)
}
plot(rast('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/raster data/Gola30mASTER.tif')) # check that everything worked out 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. get data on campsites, here vector layers ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

layers <- st_layers('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Campsites.gdb/Campsites.gdb/') # only for vector data

campsites <- vect('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Campsites.gdb/Campsites.gdb/', layer = layers$name[1]) %>% 
  st_as_sf() %>% 

# write data to files 
st_write(campsites, dsn = 'C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_campsites.shp')

# plot, to be sure that everything worked 
plot(st_read('C:/Users/filib/Documents/Praktika/RSPB/data/spatialdata/Gola_campsites.shp'))





