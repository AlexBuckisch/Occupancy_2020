# Occupancy Modeling for White-tailed deer at Chiricahua National Monument 2017
# Alex Buckisch
# abuckisch@email.arizona.edu
# May 2021

# This script will get reproducible occupancy models for presence/absence data
# of White-tailed deer (Odocoileus virginianus) detected with camera traps at
# Chricahua NM in 2017. By adjusting the "SPECIES" variable other species can be
# selected for modeling.

# All necessary data have already been formatted and uploaded to GitHub for ease
# of use, access, and reproducibility.

# In theory (!), you should be able to Ctrl+A and Ctrl+Enter in R and it should run
# without a hitch. Downloading and extracting data from my GitHub repository and
# creating plots and maps automatically. library(repmis) is important for the
# automatic download and extraction from GitHub, check to see if you need to
# select (Y/n) after loading it.

# Packages ---------------------------------------------------------------------

library(sf)
library(lubridate)
library(unmarked)
library(tidyverse)
library(tidyr)
library(rgdal)
library(units)
library(viridis)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggnewscale)
library(ggsn)
library(raster)
library(broom)
library(RColorBrewer)
library(manipulate)
library(repmis) # possibly need to select (Y/n) after loading for package to work

# Set working directory and clear environment ----------------------------------

rm(list=ls()) # clear work environment
getwd()

setwd("V:/Wildlife/Analysis/Occupancy")

# Loading data from GitHub -----------------------------------------------------

# Some can be loaded through .csv-files and tinyurls
# Complete locations of all cameras deployed
cam_locations_final   <- read.csv(file   = "https://tinyurl.com/RNR620-locations",
                                  header = TRUE)

# Assigning numeric value to sites and storing them for cross-reference
# Numeric values are easier to call up during analysis than character strings
SiteID_names          <- read.csv(file   = "https://tinyurl.com/RNR620-sites",
                                  header = TRUE)

# Others need to be read as RData files because they have date data associated
# with them (.csv would require reformatting). Tinyurl doesn't work for these...
# libary(repmis) is required to work and download properly from GitHub

# Complete deployments across multiple parks
source_data("https://github.com/AlexBuckisch/RNR_620/blob/master/Data/cam_deployments_final.RData?raw=true")

# Complete detections of relevant mammal species 
source_data("https://github.com/AlexBuckisch/RNR_620/blob/master/Data/cam_detections_final.RData?raw=true")  

# Assigning values crucial for analysis ----------------------------------------

# Variables used in code
SPECIES <- "ODVI" # 4-letter species code (ODocoileus VIrginianus)
PARK    <- "CHIR" # 4-letter park code for Chiricahua National Monument
YEAR    <- "2017" # Year of interest

# Variables used for naming plots
Species_name <- "White-tailed deer"
Park_name    <- "Chiricahua National Monument"

# Assigning first and last day of study
day_first  <- as_date(min(cam_deployments_final$DeployDate)) 
day_last   <- as_date(max(cam_deployments_final$RetrievalDate))

# Creating the presence/absence matrix -----------------------------------------

# Populate a matrix of NA's to generate a gross matrix of encounter histories 
# that can be paired down later for occupancy modeling
y_matrix            <- as.data.frame(matrix(data = NA,
                                            nrow = length(unique(cam_deployments_final$SiteID)),
                                            ncol = (day_last - day_first) + 1))
rownames(y_matrix ) <- unique(cam_deployments_final$SiteID)
colnames(y_matrix)  <- paste0(seq(from = 1, to = (day_last - day_first) + 1))

#  Create a list containing the vector of active study days for each site

activelist <- list()

k <- 0

# Iterate over unique SiteIDs to find active days for each site, relative to
# start of study

for (i in unique(cam_deployments_final$SiteID)){
  # i takes value of each SiteID
  
  k          <- k + 1
  activedays <- vector(length = 0)
  
  # get deploy dates for ith site
  deployDates <- cam_deployments_final$DeployDate[cam_deployments_final$SiteID == i]
  deployDates <- unique(deployDates)
  totalDeployDates <- length(deployDates)
  
  # iterate over all deploy dates for site i, extracting deploy + retrieval date
  for (j in 1:totalDeployDates){
    
    # extract deploy date relative to start of study
    deployday <- cam_deployments_final$DeployDate[cam_deployments_final$SiteID == i][j] -
      day_first + 1
    
    # extract retrieval date relative to start of study
    retrievalday <- cam_deployments_final$RetrievalDate[cam_deployments_final$SiteID == i][j] -
      day_first + 1
    
    # if deploy date and retrieval overlap only add deploy date, else add range
    # between both
    if (deployday != retrievalday){
      activedays <- c(activedays, deployday:retrievalday)
    } else {
      activedays <- c(activedays, deployday)
    }
    activedays <- unique(activedays)
    
    activelist[[i]] <- activedays
  }
}

# Change NA's in the y_matrix to 0's for each active day at each site 
# (NA's remain if a site was inactive)
for (i in unique(cam_deployments_final$SiteID)){
  y_matrix[i,activelist[[i]]] <- 0
  
}

# Change any 0's to 1's for the detection of the species of interest (ODVI)
# For each study day
for(i in unique(cam_detections_final$SiteID[cam_detections_final$Species_Code == SPECIES])){
  y_matrix[i, unique(cam_detections_final$Study_day[cam_detections_final$SiteID == i & 
                                                      cam_detections_final$Species_Code == SPECIES])] <- 1
}

# After creation of matrix use StdLocName for rownames which makes it easier for
# sorting and analysis
rownames(y_matrix) <- SiteID_names$StdLocName[match(x     = rownames(y_matrix),
                                                    table = SiteID_names$SiteID)]

# Restrict detections to CHIR, eliminate non-active days, and reorder
# restrict to the park
y_matrix_CHIR <- y_matrix[cam_locations_final$SiteID[cam_locations_final$UnitCode == PARK], ] 

# eliminate non-active days
y_matrix_CHIR <- y_matrix_CHIR[ ,sapply(y_matrix_CHIR, function(x)!all(is.na(x)))]             

# reorder sites by name
y_matrix_CHIR <- y_matrix_CHIR[order(rownames(y_matrix_CHIR)), ]                                

# CHIR 2017 specific matrix
# Select proper time window of deployment
y_CHIR_2017 <- y_matrix_CHIR[ ,c(1:91)]

# Rename columns
colnames(y_CHIR_2017) <- paste0("y.", 1:ncol(y_CHIR_2017))

# Covariate data ---------------------------------------------------------------

# Load site-covariates
SiteCovs <- read.csv(file   = "https://tinyurl.com/RNR620-sitecovs",
                     header = TRUE, fileEncoding = "UTF-8-BOM")

# Load observer-covariates (Temperature, precipitation, and date (time of year))

# Temperature
Temp_2017 <- read.csv(file   = "https://tinyurl.com/RNR620-temp",
                      header = TRUE)

# Precipitation
Precip_2017 <- read.csv(file   = "https://tinyurl.com/RNR620-precip",
                        header = TRUE)

# Short formatting of data frames to make them ready for analysis
Temp_2017           <- Temp_2017[ ,-1]
rownames(Temp_2017) <- rownames(y_CHIR_2017)

Precip_2017           <- Precip_2017[ ,-1]
rownames(Precip_2017) <- rownames(y_CHIR_2017)

# Date matrix needs to be created in R based on cam_detections_final
# All dates are the corresponding julian day of the year

Date_2017   <- matrix(data = NA, nrow = nrow(y_CHIR_2017), ncol = ncol(Temp_2017))
colnames(Date_2017) <- colnames(Temp_2017)
rownames(Date_2017) <- rownames(Temp_2017)

Date_2017[1, ] <- rep(seq(from = 250, length.out = ncol(Date_2017)))
Date_2017      <- Date_2017[rep(1, each = nrow(Date_2017)), ]

# Creating unmarked frame and standardizing covariates -------------------------

ODVI_umf <- unmarkedFrameOccu(
  y        = y_CHIR_2017,
  siteCovs = data.frame(Elevation       = SiteCovs$Elevation,
                        Slope           = SiteCovs$Slope,
                        Water           = SiteCovs$Water,
                        Roads           = SiteCovs$Roads,
                        Trails          = SiteCovs$Trails,
                        Cover_Field     = SiteCovs$Cover_Field,
                        Cover_Subcanopy = SiteCovs$Cover_Subcanopy,
                        Cover_Canopy    = SiteCovs$Cover_Canopy,
                        Cover_Total     = SiteCovs$Cover_Total),
  obsCovs = list(Precip = Precip_2017,
                 Temp   = Temp_2017,
                 Date   = Date_2017))

# Standardize covariates after unmarked frame is assembled
# Since all covariates are numeric, this is easily done
siteCovs(ODVI_umf) <- scale(siteCovs(ODVI_umf))
obsCovs(ODVI_umf)  <- scale(obsCovs(ODVI_umf))

# Look at the complete umarked frame to see if everything is in order
# All covariates should be in the right category and standardized
summary(ODVI_umf)

# Occupancy modeling -----------------------------------------------------------

# Starting with "Null model" add covariates just for occupancy probability first
fm1  <- occu(~1 ~1, data = ODVI_umf)
fm2  <- occu(~1 ~Elevation, data = ODVI_umf)
fm3  <- occu(~1 ~Water, data = ODVI_umf)
fm4  <- occu(~1 ~Slope, data = ODVI_umf)
fm5  <- occu(~1 ~Trails, data = ODVI_umf)
fm6  <- occu(~1 ~Roads, data = ODVI_umf)
fm7  <- occu(~1 ~Cover_Field, data = ODVI_umf)
fm8  <- occu(~1 ~Cover_Subcanopy, data = ODVI_umf)
fm9  <- occu(~1 ~Cover_Canopy, data = ODVI_umf)
fm10 <- occu(~1 ~Cover_Total, data = ODVI_umf)
fm11 <- occu(~1 ~Elevation + Slope, data = ODVI_umf)
fm12 <- occu(~1 ~Elevation + Water, data = ODVI_umf)
fm13 <- occu(~1 ~Elevation + Cover_Total, data = ODVI_umf)
fm14 <- occu(~1 ~Roads + Trails, data = ODVI_umf)
fm15 <- occu(~1 ~Elevation + Roads + Trails, data = ODVI_umf)
fm16 <- occu(~1 ~Water + Cover_Total, ODVI_umf)

models1 <- fitList("psi(.)p(.)" = fm1,
                   "psi(Elevation)p(.)" = fm2,
                   "psi(Water)p(.)" = fm3,
                   "psi(Slope)p(.)" = fm4,
                   "psi(Trails)p(.)" = fm5,
                   "psi(Roads)p(.)" = fm6,
                   "psi(Cover_Field)p(.)" = fm7,
                   "psi(Cover_Subcanopy)p(.)" = fm8,
                   "psi(Cover_Canopy)p(.)" = fm9,
                   "psi(Cover_Total)p(.)" = fm10,
                   "psi(Elevation + Slope)p(.)" = fm11,
                   "psi(Elevation + Water)p(.)" = fm12,
                   "psi(Elevation + Cover_Total)p(.)" = fm13,
                   "psi(Roads + Trails)p(.)" = fm14,
                   "psi(Elevation + Roads + Trails)p(.)" = fm15,
                   "psi(Water + Cover_Totalp(.)" = fm16 )

(MS1 <- modSel(models1))

# Now add observer covariates for effect on detection probability
fm17 <- occu(~Precip ~1, data = ODVI_umf)
fm18 <- occu(~Temp ~1 , data = ODVI_umf)
fm19 <- occu(~Date ~1, data = ODVI_umf)
fm20 <- occu(~Precip + Temp ~1, data = ODVI_umf)
fm21 <- occu(~Precip + Date ~1, data = ODVI_umf)
fm22 <- occu(~Temp + Date ~1, data = ODVI_umf)
fm23 <- occu(~Precip + Temp + Date ~1, data = ODVI_umf)

models2 <- fitList("psi(.)p(.)" = fm1,
                   "psi(.)p(Precipitation)" = fm17,
                   "psi(.)p(Temperature)" = fm18,
                   "psi(.)p(Date)" = fm19,
                   "psi(.)p(Precipitation + Temperature)" = fm20,
                   "psi(.)p(Precipitation + Date)" = fm21,
                   "psi(.)p(Temperature + Date)" = fm22,
                   "psi(.)p(Precipitation + Temperature + Date)" = fm23)

(MS2 <- modSel(models2))

# Discard all models with higher AIC than the null-model (psi(.)p(.)), since they 
# don't offer anything to interpretation.
# Add the best fitting covariates for deteciton probability to models and fit again
# This results in the best models that affect occupancy and detection probability

# Combined models for occupancy and detction probability
fm24 <- occu(~Temp ~Elevation, data = ODVI_umf)
fm25 <- occu(~Temp ~Elevation + Water, data = ODVI_umf)
fm26 <- occu(~Temp ~Elevation + Cover_Total, data = ODVI_umf)
fm27 <- occu(~Temp ~Elevation + Slope, data = ODVI_umf)
fm28 <- occu(~Temp ~Water, data = ODVI_umf)
fm29 <- occu(~Temp ~Roads, data = ODVI_umf)
fm30 <- occu(~Temp ~Elevation + Roads + Trails, data = ODVI_umf)
fm31 <- occu(~Temp ~Roads + Trails, data = ODVI_umf)
fm32 <- occu(~Date ~Elevation, data = ODVI_umf)
fm33 <- occu(~Date ~Elevation + Water, data = ODVI_umf)
fm34 <- occu(~Date ~Elevation + Cover_Total, data = ODVI_umf)
fm35 <- occu(~Date ~Elevation + Slope, data = ODVI_umf)
fm36 <- occu(~Date ~Water, data = ODVI_umf)
fm37 <- occu(~Date ~Roads, data = ODVI_umf)
fm38 <- occu(~Date ~Elevation + Roads + Trails, data = ODVI_umf)
fm39 <- occu(~Date ~Roads + Trails, data = ODVI_umf)
fm40 <- occu(~Precip + Date ~Elevation, data = ODVI_umf)
fm41 <- occu(~Precip + Date ~Elevation + Water, data = ODVI_umf)
fm42 <- occu(~Precip + Date ~Elevation + Cover_Total, data = ODVI_umf)
fm43 <- occu(~Precip + Date ~Elevation + Slope, data = ODVI_umf)
fm44 <- occu(~Precip + Date ~Water, data = ODVI_umf)
fm45 <- occu(~Precip + Date ~Roads, data = ODVI_umf)
fm46 <- occu(~Precip + Date ~Elevation + Roads + Trails, data = ODVI_umf)
fm47 <- occu(~Precip + Date ~Roads + Trails, data = ODVI_umf)
fm48 <- occu(~Temp + Date ~Elevation, data = ODVI_umf)
fm49 <- occu(~Temp + Date ~Elevation + Water, data = ODVI_umf)
fm50 <- occu(~Temp + Date ~Elevation + Cover_Total, data = ODVI_umf)
fm51 <- occu(~Temp + Date ~Elevation + Slope, data = ODVI_umf)
fm52 <- occu(~Temp + Date ~Water, data = ODVI_umf)
fm53 <- occu(~Temp + Date ~Roads, data = ODVI_umf)
fm54 <- occu(~Temp + Date ~Elevation + Roads + Trails, data = ODVI_umf)
fm55 <- occu(~Temp + Date ~Roads + Trails, data = ODVI_umf)


models3 <- fitList("psi(.)p(.)" = fm1,
                   "psi(Elevation)p(.)" = fm2,
                   "psi(Water)p(.)" = fm3,
                   "psi(Slope)p(.)" = fm4,
                   "psi(Trails)p(.)" = fm5,
                   "psi(Roads)p(.)" = fm6,
                   "psi(Cover_Field)p(.)" = fm7,
                   "psi(Cover_Subcanopy)p(.)" = fm8,
                   "psi(Cover_Canopy)p(.)" = fm9,
                   "psi(Cover_Total)p(.)" = fm10,
                   "psi(Elevation + Slope)p(.)" = fm11,
                   "psi(Elevation + Water)p(.)" = fm12,
                   "psi(Elevation + Cover_Total)p(.)" = fm13,
                   "psi(Roads + Trails)p(.)" = fm14,
                   "psi(Elevation + Roads + Trails)p(.)" = fm15,
                   "psi(Water + Cover_Totalp(.)" = fm16,
                   "psi(.)p(Precipitation)" = fm17,
                   "psi(.)p(Temperature)" = fm18,
                   "psi(.)p(Precipitation + Temperature)" = fm20,
                   "psi(Elevation)p(Temperature)" = fm24,
                   "psi(Elevation + Water)p(Temperature)" = fm25,
                   "psi(Elevation + Cover_Total)p(Temperature)" = fm26,
                   "psi(Elevation + Slope)p(Temperature)" = fm27,
                   "psi(Water)p(Temperature)" = fm28,
                   "psi(Roads)p(Temperature)" = fm29,
                   "psi(Elevation + Roads + Trails)p(Temperature)" = fm30,
                   "Psi(Roads + Trails)p(Temperature)" = fm31,
                   "psi(Elevation)p(Date)" = fm32,
                   "psi(Elevation + Water)p(Date)" = fm33,
                   "psi(Elevation + Cover_Total)p(Date)" = fm34,
                   "psi(Elevation + Slope)p(Date)" = fm35,
                   "psi(Water)p(Date)" = fm36,
                   "psi(Roads)p(Date)" = fm37,
                   "psi(Elevation + Roads + Trails)p(Date)" = fm38,
                   "Psi(Roads + Trails)p(Date)" = fm39,
                   "psi(Elevation)p(Precipitation + Date)" = fm40,
                   "psi(Elevation + Water)p(Precipitation + Date)" = fm41,
                   "psi(Elevation + Cover_Total)p(Precipitation + Date)" = fm42,
                   "psi(Elevation + Slope)p(Precipitation + Date)" = fm43,
                   "psi(Water)p(Precipitation + Date)" = fm44,
                   "psi(Roads)p(Precipitation + Date)" = fm45,
                   "psi(Elevation + Roads + Trails)p(Precipitation + Date)" = fm46,
                   "Psi(Roads + Trails)p(Precipitation + Date)" = fm47,
                   "psi(Elevation)p(Temperature + Date)" = fm48,
                   "psi(Elevation + Water)p(Temperature + Date)" = fm49,
                   "psi(Elevation + Cover_Total)p(Temperature + Date)" = fm50,
                   "psi(Elevation + Slope)p(Temperature + Date)" = fm51,
                   "psi(Water)p(Temperature + Date)" = fm52,
                   "psi(Roads)p(Temperature + Date)" = fm53,
                   "psi(Elevation + Roads + Trails)p(Temperature + Date)" = fm54,
                   "Psi(Roads + Trails)p(Temperature + Date)" = fm55)

(MS3 <- modSel(models3))


# FREQUENTIST APPROACH ---------------------------------------------------------

# Using the Frequentist-approach, first fit a 'rich' model that includes every
# covariate. Following, look at the individual p-values and eliminate the
# covariate with the highest p-value. Start with the detection model (i.e.
# observer covariates) and then the occupancy-model (i.e. site covariates). The
# the complete and best model is achieved once all covariates have a p-value of
# less than 0.10 (arbitrary value).

# Null model 
fm1f <- occu(~1 ~1, data = ODVI_umf)

# Rich model
fm2f <- occu(~Precip + Temp + Date ~Elevation + Slope + Water + Roads + Trails +
              Cover_Field + Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm2f) # eliminate Temp (p > 0.6)

# Without Temp in detection model
fm3f <- occu(~Precip + Date ~Elevation + Slope + Water + Roads + Trails +
              Cover_Field + Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm3f) # all detection covariates < 0.10 (detection model is set)

# Eliminate Roads
fm4f <- occu(~Precip + Date ~Elevation + Slope + Water + Trails +
              Cover_Field + Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm4f)

# Eliminate Cover_Field
fm5f <- occu(~Precip + Date ~Elevation + Slope + Water + Trails +
              Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm5f)

# Eliminate Cover_Total
fm6f <- occu(~Precip + Date ~Elevation + Slope + Water + Trails + Cover_Subcanopy +
              Cover_Canopy, data = ODVI_umf)

summary(fm6f)

# Eliminate Cover_Subcanopy
fm7f <- occu(~Precip + Date ~Elevation + Slope + Water + Trails + Cover_Canopy,
            data = ODVI_umf)

summary(fm7f)

# Eliminate Slope
fm8f <- occu(~Precip + Date ~Elevation + Water + Trails + Cover_Canopy,
            data = ODVI_umf)

summary(fm8f)

# Eliminate Trails
fm9f <- occu(~Precip + Date ~Elevation + Water + Cover_Canopy,
            data = ODVI_umf)

summary(fm9f)

# Eliminate Cover_Canopy
fm10f <- occu(~Precip + Date ~Elevation + Water, data = ODVI_umf)

summary(fm10f)

# Eliminate Water
fm11f <- occu(~Precip + Date ~Elevation, data = ODVI_umf)

summary(fm11f)

# Best model is fm11f - only uses Elevation as a site-covariate and Precip, Date
# as a observer covariate

# Predicting occupancy and detection probability -------------------------------

# Occupancy and detection probability for the null model
backTransform(fm1, type = "state")
backTransform(fm1, type = "det")

# Occupancy and detection probability for the top model (fm32)

summary(fm32)
# Elevation is a significant covariate on occupancy probability (p < 0.05)
# Date is a highly significant covariate on detection probability (p < 0.01)

plogis(coef(fm32)) # backtransforming model on probability scale
# psi(Elevation) = 0.31
# p(Date) = 0.43                                        

# When creating data frame for the predicted values, the name of the column needs
# to be the same name as the variable used in the modeling
# i.e. Elevation, Temp, Date, Precip, Roads, etc. otherwise 'predict' won't work

# Predicting occupancy and detection probability for the top model (fm32)
Elev_pred <- data.frame(Elevation = seq(from = min(ODVI_umf@siteCovs$Elevation),
                                        to   = max(ODVI_umf@siteCovs$Elevation),
                                        length.out = 100))

Elev_orig <- data.frame(Elev_orig = seq(from = min(SiteCovs$Elevation),
                                        to   = max(SiteCovs$Elevation),
                                        length.out = 100))

Occu_Elev_pred <- predict(fm32, newdata = Elev_pred, type = "state")

# Predicting detection probability based on date for the top model (fm32)
Date_pred <- data.frame(Date       = seq(from = min(ODVI_umf@obsCovs$Date),
                                         to   = max(ODVI_umf@obsCovs$Date),
                                         length.out = 91)) # days of deployment
Date_orig <- data.frame(Date_orig  = seq(from = 250, to = 340,
                                         length.out = 91)) # julian days of the year (same deployment days)

Det_Date_pred <- predict(fm32, newdata = Date_pred, type = "det")

# Predicting occupancy prob. based on Water from best model with this covariate
# Model fm33
# This model uses two site covariates but only Water will be plotted
# Elevation needs to be assigned a value, usually the mean
Water_pred <- data.frame(Water = seq(from       = min(ODVI_umf@siteCovs$Water),
                                     to         = max(ODVI_umf@siteCovs$Water),
                                     length.out = 100),
                         Elevation   = mean(ODVI_umf@siteCovs$Elevation))

Water_orig <- data.frame(Water_orig  = seq(from       = min(SiteCovs$Water),
                                           to         = max(SiteCovs$Water),
                                           length.out = 100))

Occu_Water_pred <- predict(fm33, newdata = Water_pred, type = "state")

# Predicting occupancy prob. based on total Cover from best model with this covariate
# Model fm34
Cover_pred <- data.frame(Cover_Total = seq(from = min(ODVI_umf@siteCovs$Cover_Total),
                                           to   = max(ODVI_umf@siteCovs$Cover_Total),
                                           length.out = 100),
                         Elevation   = mean(ODVI_umf@siteCovs$Elevation))

Cover_orig <- data.frame(Cover_orig  = seq(from = min(SiteCovs$Cover_Total),
                                           to   = max(SiteCovs$Cover_Total),
                                           length.out = 100))

Occu_Cover_pred <- predict(fm34, newdata = Cover_pred, type = "state")

# Predicting detection prob. based on temperature from best model with these covariates
# Model fm48
Temp_pred <- data.frame(Temp = seq(from       = min(ODVI_umf@obsCovs$Temp), 
                                   to         = max(ODVI_umf@obsCovs$Temp),
                                   length.out = 91),
                        Date = mean(ODVI_umf@obsCovs$Date))

Temp_orig <- data.frame(Temp_orig = seq(from       = min(Temp_2017), 
                                        to         = max(Temp_2017), 
                                        length.out = 91)) # julian day of the year

Det_Temp_pred <- predict(fm48, newdata = Temp_pred, type = "det")

# Plotting occupancy and detection probability against covariates---------------

# Make new data frame consisting of standardized elevation, original elevation,
# and predicted occupancy
gg_data_elev <- data.frame(Elevation_predicted = Elev_pred,
                           Elevation_original  = Elev_orig,
                           Occupancy           = Occu_Elev_pred)

# Make a new data frame consisting of standardized date, original date,
# and predicted detection probability
gg_data_date <- data.frame(Date_predicted = Date_pred,
                           Date_original  = Date_orig,
                           Detection      = Det_Date_pred)

# Make a new data frame consisting of standardized water distance, original water
# distance, and predicted occupancy
gg_data_water <- data.frame(Water_predicted = Water_pred[ ,-2], # delete Elevation
                            Water_original  = Water_orig,
                            Occupancy       = Occu_Water_pred)

# Make a new data frame consisting of standardized total cover, original total
# cover, and predicted occupancy
gg_data_cover <- data.frame(Cover_predicted = Cover_pred[ ,-2], # delete Elevation
                            Cover_original  = Cover_orig,
                            Occupancy       = Occu_Cover_pred)

# Make a new data frame consisting of standardized temperature, original
# temperature, and predicted detection probability
gg_data_temp <- data.frame(Temp_predicted = Temp_pred[ ,-2], # delete Date
                           Temp_original  = Temp_orig,
                           Detection      = Det_Temp_pred)

# Plot predicted occupancy vs original elevation based on predictions with
# standardized elevation
(Elev_gg <- ggplot(data = gg_data_elev, aes(x = Elev_orig, y = Occupancy.Predicted)) +
    geom_line(size = 2, color = "darkblue") +
    geom_ribbon(aes(ymin  = Occupancy.lower,
                    ymax  = Occupancy.upper, 
                    alpha = 0.1), fill = "lightgrey") +
    theme_bw() +
    labs(title    = paste0("Occupancy Probability of ", Species_name, " vs Elevation"),
         subtitle = paste0(Park_name, " ", YEAR),
         x        = "Elevation [m]",
         y        = "Predicted Occupancy Probability") +
    theme(plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_x_continuous(limits = round(c(min(gg_data_elev$Elev_orig),
                                        max(gg_data_elev$Elev_orig))),
                       breaks = c(1550,1750,2000,2250)))

# Plot predicted Occupancy vs original total cover based on predicitions with
# standardized total cover
(Cover_gg <- ggplot(data = gg_data_cover, aes(x = Cover_orig, y = Occupancy.Predicted)) +
    geom_line(size = 2, color = "darkblue") +
    geom_ribbon(aes(ymin  = Occupancy.lower,
                    ymax  = Occupancy.upper,
                    alpha = 0.1),  fill = "lightgrey") +
    theme_bw() +
    labs(title    = paste0("Occupancy Probability of ", Species_name, " vs Total Cover"),
         subtitle = paste0(Park_name, " ", YEAR),
         x        = "Total Vegetation Cover",
         y        = "Predicted Occupancy Probability") +
    theme(plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5),
          legend.position = "none"))

# Plot predicted occupancy probability vs water distance based on prediction
# with standardized water distance
(Wat_gg <- ggplot(data = gg_data_water, aes(x = Water_orig, y = Occupancy.Predicted)) +
    geom_line(size = 2, color = "midnightblue") +
    geom_ribbon(aes(ymin  = Occupancy.lower,
                    ymax  = Occupancy.upper,
                    alpha = 0.1),  fill = "lightgrey") +
    theme_bw() +
    labs(title    = paste0("Occupancy Probability of ", Species_name, " vs Water distance"),
         subtitle = paste0(Park_name, " ", YEAR),
         x        = "Water distance [m]",
         y        = "Predicted Occupancy Probability") +
    theme(plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000)) +
    scale_y_continuous(limits = 0:1))

# Plot predicted detection probability vs date (Time of year) based on predictions
# with standardized dates
(Date_gg <- ggplot(data = gg_data_date, aes(x = Date_orig, y = Detection.Predicted)) +
    geom_line(size = 2, color = "darkgreen") +
    geom_ribbon(aes(ymin  = Detection.lower,
                    ymax  = Detection.upper,
                    alpha = 0.1),  fill = "lightgrey") +
    theme_bw() +
    labs(title    = paste0("Detection Probability of ", Species_name, " vs Time of Year"),
         subtitle = paste0(Park_name, " ", YEAR),
         x        = "Time of Year",
         y        = "Predicted Detection Probability") +
    theme(plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_x_continuous(breaks = c(250, 275, 306, 336),
                       labels = c("September", "October", "November",
                                  "December")))

# Plot predicted detection probability vs temperature based on predictions with
# standardized temperature
(Temp_gg <- ggplot(data = gg_data_temp, aes(x = Temp_orig, y = Detection.Predicted)) +
    geom_line(size = 2, color = "darkgreen") +
    geom_ribbon(aes(ymin  = Detection.lower,
                    ymax  = Detection.upper,
                    alpha = 0.1),  fill = "lightgrey") +
    theme_bw() +
    labs(title    = paste0("Detection Probability of ", Species_name, " vs Temperature"),
         subtitle = paste0(Park_name, " ", YEAR),
         x        = "Temperature [C]",
         y        = "Predicted Detection Probability") +
    theme(plot.title      = element_text(hjust = 0.5),
          plot.subtitle   = element_text(hjust = 0.5),
          legend.position = "none"))

# All plots together
ggarrange(nrow = 2, ncol = 3,
          Elev_gg, Wat_gg, Cover_gg, Date_gg, Temp_gg)

# Distribution mapping ---------------------------------------------------------

# Read in file that contains all XY coordinates from park of interest
# All coordinates have to be in UTM to make multiple layers fit in ggplot

# 30 m fishnet from CHIR
Fishnet30 <- data.frame(read.csv(file   = "https://tinyurl.com/RNR620-fishnet",
                                 header = TRUE))

# Make new data frame with Elevation data from fishnet that can be used for
# predict
Fishnet30_Elev           <- data.frame(Fishnet30$RASTERVALU)
colnames(Fishnet30_Elev) <- "Elevation"


# Standardize elevation values
Fishnet30_zElev <- data.frame(scale(Fishnet30_Elev))

# Use the best model with elevation for predictions
# This generates predictions of occupancy probability for all 76,535 points
# generated by the fishnet. After trying out different resolutions (10, 30, 50 m)
# I decided to use a 30 m fishnet resolution. Better than 50, but reasonably fast
# executed. The 10 m fishnet generated nearly 1 million points, and predicting
# occupancy for those took almost 5 hours. The 30 m will only take a minute to
# simulate.
# For this assignment, the generated .csv-file can be loaded through GitHub,
# no need to simulate.

Map_pred_30 <- data.frame(read.csv(file = "https://tinyurl.com/RNR620-predict")
                          [ ,-1])

# Make gg data frame with predictions and original values
gg_data_map30 <- data.frame(X_coordinate = Fishnet30$POINT_X,
                            Y_coordinate = Fishnet30$POINT_Y,
                            Occupancy    = Map_pred_30$Predicted)

# Download all shapefiles required for mapping from GitHub
download.file(url = "https://tinyurl.com/RNR620-shapefiles",
              destfile = "Shapefiles.zip", mode = "wb")

# Unzip folder
Shapefiles <- unzip(zipfile = "Shapefiles.zip")

# Assign shapefiles to variables for plotting

# Roads
CHIR_Roads <- st_read("./Shapefiles/CHIR_Roads_Park.shp")
# Trails
CHIR_Trails <- st_read("./Shapefiles/CHIR_Trails.shp")
# Boundary
CHIR_Boundary <- st_read("./Shapefiles/CHIR_Boundary.shp")
# DEM
DEM <- st_read("./Shapefiles/CHIR_Fishnet_30m_Elevation.shp")

# DEM
# Formatting DEM requires a copy
DEM2 <- data.frame(DEM)
DEM2 <- DEM2[ ,-4]
colnames(DEM2) <- c("X", "Y", "Elevation")

# Plotting DEM
(gg_DEM <- ggplot(data = DEM2, aes(x = X,
                                   y = Y)) +
    geom_raster(aes(fill = Elevation), alpha = 0.75) +
    scale_fill_gradientn(name = "Elevation", colours = terrain.colors(100)) +
    theme_bw() +
    rremove("grid") +
    coord_equal() +
    xlab("Longitude") +
    ylab("Latitude"))

# Making DEM to raster to get hillshade
DEM_raster <- rasterFromXYZ(DEM2)
plot(DEM_raster)

# Assign correct coordinate reference system (coordinates are in UTM, WGS1984)
crs(DEM_raster) <- "+proj=utm +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# Getting slope and aspect from DEM raster
DEM_terrain <- terrain(DEM_raster, opt = c("slope", "aspect"), unit = "radians")

# Creating hillshade
slope  <- terrain(DEM_raster, opt = "slope")
aspect <- terrain(DEM_raster, opt = "aspect")
hill   <- hillShade(slope, aspect, angle = 60, diretion = 270)

hill_spdf         <- as(hill, "SpatialPixelsDataFrame")
hill_df           <- data.frame(hill_spdf)
colnames(hill_df) <- c("Value", "X_coord", "Y_coord")

# Just Hillshade map
(gg_Hill <- ggplot() + 
    geom_tile(data = hill_df, aes(x = X_coord, y = Y_coord, fill = Value)) + 
    scale_fill_gradient(low = "black", high = "white") +
    new_scale_fill() +
    geom_tile(data = DEM2, aes(x = X, y = Y, fill = Elevation), alpha = 0.4) + 
    scale_fill_gradientn(colours = rev(terrain.colors(10))) +
    theme_bw() + 
    theme(legend.position = "none"))

# Plot complete distribution map (Occupancy + Hillshade + Shapefiles)

OccuMap_30 <- ggplot() +
  # Hillshade layer
  geom_tile(data = hill_df, aes(x    = X_coord,
                                y    = Y_coord,
                                fill = Value)) + 
  scale_fill_gradient(low = "black", high = "white") +
  new_scale_fill() +
  # Occupancy layer
  geom_tile(data = gg_data_map30, aes(x    = X_coordinate,
                                      y    = Y_coordinate,
                                      fill = Occupancy), alpha = 0.75) +
  scale_fill_viridis(option = "viridis", direction = 1) +
  # CHIR roads
  geom_sf(data = CHIR_Roads$geometry, color = "black", size = 1) +
  # CHIR trails
  geom_sf(data = CHIR_Trails$geometry, color = "dimgray", size = 0.7,
          linetype = "dashed") +
  # CHIR boundary
  geom_sf(data = CHIR_Boundary$geometry, color = "black", fill = "NA", size = 1.5) +
  # Plot anatomy        
  theme_bw() + 
  rremove("grid") +
  north(data = CHIR_Boundary, location = "topright", symbol = 3) +
  ggsn::scalebar(data      = CHIR_Boundary, location = "bottomleft", dist = 1,
                 dist_unit = "km", transform = FALSE) +
  labs(title    = "Occupancy probability of White-tailed Deer",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Longitude",
       y        = "Latitude") +
  theme(plot.title       = element_text(hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5),
        legend.position  = "right",
        legend.box       = "vertical",
        legend.key.size  = unit(2, "cm"),
        legend.key.width = unit(0.8, "cm")) +
  guides(fill_new = FALSE, # remove hillshade from legend
         fill     = guide_colorbar(ticks.colour = "black",
                                   frame.colour = "black"))

OccuMap_30 # if there is a warining message ('reached elapsed time limit') and
# the map is missing some layers, run the line again. It should look
# like a real map.
