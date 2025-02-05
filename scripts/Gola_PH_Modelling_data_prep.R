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

# format data types properly
str(pres)

pres %>%
  mutate(Date = as.Date(VisitDate_dd_mm_yyyy), 
         Time = hms::as_hms(ifelse(is.na(VisitTime_hh_mm_ss), format(VisitDate_dd_mm_yyyy, '%H:%M:%S'), format(VisitTime_hh_mm_ss, '%H:%M:%S')))) %>%
  select(DatasetName, Date, Time, VisitDate_dd_mm_yyyy, VisitTime_hh_mm_ss)

# replace all -9999 with NA 
pres[pres == -9999] <- NA 

















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