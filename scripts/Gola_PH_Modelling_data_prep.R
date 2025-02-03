#### data preparation script for pigmy hippo occupancy modelling 
#### script written by Filibert Heim, filibert.heim@rspb.org.uk or filibert.heim@posteo.de


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### 1. Preparatations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install needed packages 
#install.packages('readxl')
#install.packages('Microsoft365R')

# load packages
library(Microsoft365R)
library(readxl)
library(tidyverse)
library(lubridate)
library(spOccupancy)
filter - dplyrfilter
select - dplyrselect
rename - dplyrrename 

# read in data 




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