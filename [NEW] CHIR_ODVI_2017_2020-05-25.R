# SINGLE-SEASON OCCUPANCY MODEL FOR ODVI IN CHIR 2017 --------------------------
# Alex Buckisch
# abuckisch@email.arizona.edu
# June 2020

# This script requires 'Camdata_2020.R' script to be loaded into the R environment
# and run before it will work

# PACKAGES ---------------------------------------------------------------------
library(unmarked)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggsn)
library(viridis)
library(sf)
library(lubridate)
library(unmarked)
library(rgdal)
library(raster)
library(units)
library(broom)
library(RColorBrewer)
library(manipulate)

# LOAD CAMDATA_2020.R ----------------------------------------------------------

# Load Camdata.R so this scripts executes without error
source("C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/R/Models/CAMDATA/Camdata_2020.R")

# SET WORKING DIRECTORY---------------------------------------------------------
setwd("C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/R/Models/CHIR/ODVI")

# SET SPECIES, PARK OF INTEREST, AND SAMPLING PERIOD----------------------------

# Variables used in code
SPECIES <- "ODVI"
PARK    <- "CHIR"
YEAR    <- "2017"

# Variables used for naming plots
Species_name <- "White-tailed deer"
Park_name    <- "Chiricahua National Monument"

# Select sampling period in number of days
SAMPLE <- 7

# CREATE THE PRESENCE/ABSENCE MATRIX -------------------------------------------

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


# For multiple species use %in% operator (substitute above with below

# Change any 0's to 1's for the detection of the species of interest (ODVI and LECA) 
# For each study day
#for(i in unique(cam_detections_final$SiteID[cam_detections_final$Species_Code %in% c(SPECIES, SPECIES2)])){
#  y_matrix[i, unique(cam_detections_final$Study_day[cam_detections_final$SiteID == i & 
#                                                    cam_detections_final$Species_Code %in% c(SPECIES, SPECIES2)])] <- 1
#}



# After creation of matrix use StdLocName for rownames which makes it easier for
# sorting and analysis
rownames(y_matrix) <- SiteID_names$StdLocName[match(x     = rownames(y_matrix),
                                                    table = SiteID_names$SiteID)]


# CREATE PRESENCE/ABSENCE MATRIX FOR SELECTED SAMPLING PERIOD ------------------

# Fill matrix with NAs based on dimensions of y_matrix (divided by SAMPLE)
y_matrix_sample <- data.frame(matrix(data = NA, nrow = nrow(y_matrix),
                                     ncol = floor(ncol(y_matrix) / SAMPLE)))

# Assign new matrix the same rownames
rownames(y_matrix_sample) <- rownames(y_matrix)

# Fill in y_matrix span of days
for(i in 1:ncol(y_matrix_sample)){
  day_end <- i * SAMPLE
  day_start <- day_end - SAMPLE + 1
  
  # only assign NA if all days are NA
  day_sum <- rowSums((unlist(sapply(X = y_matrix[ ,c(day_start:day_end)], FUN = is.na))))
  all_na  <- day_sum == SAMPLE
  #y_matrix_sample[all_na,i] <- NA
  
  # if it's greater than 0 (i.e. animal detected) assign 1 to sampling period
  y_matrix_sample[!all_na,i] <- rowSums(y_matrix[!all_na, c(day_start:day_end)],
                                        na.rm = TRUE) > 0
  y_matrix_sample[ ,i] <- as.numeric(y_matrix_sample[ ,i])
  
}


# NAIVE OCCUPANCY --------------------------------------------------------------

# Count number of detections for each site 
no_detections_1           <- data.frame(rowSums(x = y_matrix, na.rm = TRUE))
colnames(no_detections_1) <- paste0("Number of Detections of ", SPECIES)

# Sampling period 5 days
no_detections_5           <- data.frame(rowSums(x = y_matrix_sample, na.rm = TRUE))
colnames(no_detections_5) <- paste0("Number of Detections of ", SPECIES)

# Sampling period 10 days
no_detections_10           <- data.frame(rowSums(x = y_matrix_sample, na.rm = TRUE))
colnames(no_detections_10) <- paste0("Number of Detections of ", SPECIES)

# CHIR SPECIFIC PRESENCE/ABSENCE MATRIX ----------------------------------------

# Restrict detections to CHIR, eliminate non-active days, and reorder
# restrict to the park

# original (1 day)
y_matrix_CHIR1 <- y_matrix[cam_locations_final$SiteID[cam_locations_final$UnitCode == PARK], ] 

# sampling period specific
y_matrix_CHIR <- y_matrix_sample[cam_locations_final$SiteID[cam_locations_final$UnitCode == PARK], ] 

# eliminate non-active days (1 day)
y_matrix_CHIR1 <- y_matrix_CHIR1[ ,sapply(y_matrix_CHIR1, function(x)!all(is.na(x)))]             

# eliminate non-active days
y_matrix_CHIR <- y_matrix_CHIR[ ,sapply(y_matrix_CHIR, function(x)!all(is.na(x)))]  

# reorder sites by name
y_matrix_CHIR <- y_matrix_CHIR[order(rownames(y_matrix_CHIR)), ]                                

# CHIR 2017 specific matrix
# Select proper time window of deployment
y_CHIR_2017 <- y_matrix_CHIR[ ,c(1:13)] # 7 days
y_CHIR1_2017 <- y_matrix_CHIR1[ ,c(1:91)] # 1 day

# Rename columns
colnames(y_CHIR_2017) <- paste0("Week.", 1:ncol(y_CHIR_2017))

# Number of detections during study periods
# 1 day
sum(x = y_CHIR1_2017, na.rm = TRUE)

# 7 days
sum(x = y_CHIR_2017, na.rm = TRUE)

# SCALE AND FORMAT COVARIATES --------------------------------------------------
# Read in site covariates
SiteCovs <- read.csv(file   = "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/Covariates/All_SiteCovariates/CHIR_Covariates_All.csv",
                     header = TRUE, fileEncoding = "UTF-8-BOM")

# Read in observer-covariates - since these change every time it is best to 
# format them in R individually

# Precipitation and Time need to be formatted first
Precip <- read.csv(file = "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/Covariates/Climate/CHIR_Precip_All.csv")
Temp   <- read.csv(file = "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/Covariates/Climate/CHIR_Temp_All.csv")


# Transpose columns and rows to get them into the right covariate alignment
Precip <- data.frame(t(Precip))
Temp   <- data.frame(t(Temp))

# Rename column names to dates (starting at 01/01/2017 and ending at 31/12/2019)
# Total of 1095 days
colnames(Precip) <- paste0("Date.", 1:ncol(Precip))
colnames(Temp)   <- paste0("Date.", 1:ncol(Temp))

# Autopopulate values for number of plots in CHIR, i.e. 41 rows
Precip <- Precip[rep(1, each = nrow(y_matrix_CHIR)), ]
Temp   <- Temp[rep(1, each = nrow(y_matrix_CHIR)), ]

# Rename rows to StdLocName
rownames(Precip) <- rownames(y_CHIR_2017)
rownames(Temp)   <- rownames(y_CHIR_2017)

# From start date 09/07/2017 to end date 12/06/2017
Precip_2017 <- Precip[ ,c(250:340)]
Temp_2017   <- Temp[ ,c(250:340)] 

# Summarize and average Temp and Precip over individual sampling periods

# Fill matrix with NAs based on dimensions of y_matrix_sample (divided by SAMPLE)
Precip_sample <- data.frame(matrix(data = NA, nrow = nrow(Precip_2017),
                                   ncol = floor(ncol(Precip_2017) / SAMPLE)))

# Assign new matrix the same rownames
rownames(Precip_sample) <- rownames(Precip_2017)

# Fill in y_matrix span of days
for(i in 1:ncol(Precip_sample)){
  day_end <- i * SAMPLE
  day_start <- day_end - SAMPLE + 1
  
  # Error checking
  #cat("day_start = ", day_start, "\n")
  #cat("day_end = ", day_end, "\n")
  
  # Calculate mean for all rows starting at day_start ending at day_end
  sample_mean <- rowMeans(Precip_2017[ ,c(day_start:day_end)])
  Precip_sample[ ,i] <- sample_mean
  
  }

# Fill matrix with NAs based on dimensions of y_matrix_sample (divided by SAMPLE)
Temp_sample <- data.frame(matrix(data = NA, nrow = nrow(Temp_2017),
                                 ncol = floor(ncol(Temp_2017) / SAMPLE)))

# Assign new matrix the same rownames
rownames(Temp_sample) <- rownames(Temp_2017)

# Fill in y_matrix span of days
for(i in 1:ncol(Temp_sample)){
  day_end <- i * SAMPLE
  day_start <- day_end - SAMPLE + 1
  
  # Error checking
  #cat("day_start = ", day_start, "\n")
  #cat("day_end = ", day_end, "\n")
  
  # Calculate mean for all rows starting at day_start ending at day_end
  sample_mean <- rowMeans(Temp_2017[ ,c(day_start:day_end)])
  Temp_sample[ ,i] <- sample_mean
  
}

# Date and Time matrices need to be created in R based on cam_detections_final
# All dates are the corresponding julian day of the year

# Date
Date_2017   <- matrix(data = NA, nrow = nrow(y_CHIR_2017), ncol = ncol(Temp_sample))
colnames(Date_2017) <- paste0("Week.", 1:ncol(Date_2017))
rownames(Date_2017) <- rownames(Temp_2017)

Date_2017[1, ] <- rep(seq(from = 36, length.out = ncol(Date_2017)))
Date_2017      <- Date_2017[rep(1, each = nrow(Date_2017)), ]


# Time (will be added later)
# Time_2017 <- matrix(data = NA, nrow = nrow(y_CHIR_2017), ncol = ncol(Temp_2017))
# colnames(Time_2017) <- colnames(Temp_2017)
# rownames(Time_2017) <- rownames(Temp_2017)

# Make unmarked frame
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
  obsCovs = list(Precip = Precip_sample,
                 Temp   = Temp_sample,
                 Date   = Date_2017))

# Standardize covariates after unmarked frame is assembled
# Since all covariates are numeric, this is easily done
siteCovs(ODVI_umf) <- scale(siteCovs(ODVI_umf))
obsCovs(ODVI_umf)  <- scale(obsCovs(ODVI_umf))

# Look at the complete umarked frame to see if everything is in order
# All covariates should be in the right category and standardized
summary(ODVI_umf)

# OCCUPANCY MODELING -----------------------------------------------------------

# Using the Frequentist-approach, first fit a 'rich' model that includes every
# covariate. Following, look at the individual p-values and eliminate the
# covariate with the highest p-value. Start with the detection model (i.e.
# observer covariates) and then the occupancy-model (i.e. site covariates). The
# the complete and best model is achieved once all covariates have a p-value of
# less than 0.10 (arbitrary value).

# Null model 
fm1 <- occu(~1 ~1, data = ODVI_umf)

# Rich model
fm2 <- occu(~Precip + Temp + Date ~Elevation + Slope + Water + Roads + Trails +
                                   Cover_Field + Cover_Subcanopy + Cover_Canopy +
                                   Cover_Total, data = ODVI_umf)

summary(fm2) # eliminate Temp (p > 0.6)

# Without Temp in detection model
fm3 <- occu(~Precip + Date ~Elevation + Slope + Water + Roads + Trails +
              Cover_Field + Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm3) # all detection covariates < 0.10 (detection model is set)

# Eliminate Roads
fm4 <- occu(~Precip + Date ~Elevation + Slope + Water + Trails +
              Cover_Field + Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm4)

# Eliminate Cover_Field
fm5 <- occu(~Precip + Date ~Elevation + Slope + Water + Trails +
              Cover_Subcanopy + Cover_Canopy +
              Cover_Total, data = ODVI_umf)

summary(fm5)

# Eliminate Cover_Total
fm6 <- occu(~Precip + Date ~Elevation + Slope + Water + Trails + Cover_Subcanopy +
              Cover_Canopy, data = ODVI_umf)

summary(fm6)

# Eliminate Cover_Subcanopy
fm7 <- occu(~Precip + Date ~Elevation + Slope + Water + Trails + Cover_Canopy,
            data = ODVI_umf)

summary(fm7)

# Eliminate Slope
fm8 <- occu(~Precip + Date ~Elevation + Water + Trails + Cover_Canopy,
            data = ODVI_umf)

summary(fm8)

# Eliminate Trails
fm9 <- occu(~Precip + Date ~Elevation + Water + Cover_Canopy,
            data = ODVI_umf)

summary(fm9)

# Eliminate Cover_Canopy
fm10 <- occu(~Precip + Date ~Elevation + Water, data = ODVI_umf)

summary(fm10)

# Eliminate Water
fm11 <- occu(~Precip + Date ~Elevation, data = ODVI_umf)

summary(fm11)

# Best model is fm11 - only uses Elevation as a site-covariate and Precip, Date
# as a observer covariate

# PREDICTING OCCUPANCY AND DETECTION PROBABILITY -------------------------------

# Predicting occupancy and detection probability for the null model
backTransform(fm1, type = "state")
backTransform(fm1, type = "det")

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

# Predicting occupancy prob. based on distance to roads and trails from best model
# with these covariates
# Model fm38
RoaTra_pred <- data.frame(Roads  = seq(from = min(ODVI_umf@siteCovs$Roads),
                                       to    = max(ODVI_umf@siteCovs$Roads),
                                       length.out = 100),
                          Trails = seq(from = min(ODVI_umf@siteCovs$Trails),
                                       to   = max(ODVI_umf@siteCovs$Trails),
                                       length.out = 100),
                          Elevation = mean(ODVI_umf@siteCovs$Elevation))

RoaTra_orig <- data.frame(Roads  = seq(from = min(SiteCovs$Roads),
                                       to   = max(SiteCovs$Roads),
                                       length.out = 100),
                          Trails = seq(from = min(SiteCovs$Trails),
                                       to   = max(SiteCovs$Trails),
                                       length.out = 100))

Occu_RoaTra_pred <- predict(fm38, newdata = RoaTra_pred, type = "state")

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


# PLOTTTING --------------------------------------------------------------------

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

# Make a new data frame consisting of standardized distance to roads and trails,
# original distance, and predicted occupancy
gg_data_RoaTra <- data.frame(Roads_predicted  = RoaTra_pred[ ,1],
                             Trails_predicted = RoaTra_pred[ ,2],
                             Roads_original   = RoaTra_orig[ ,1],
                             Trails_original  = RoaTra_orig[ ,2],
                             Occupancy        = Occu_RoaTra_pred)

# Make a new data frame consisting of standardized temperature, original
# temperature, and predicted detection probability
gg_data_temp <- data.frame(Temp_predicted = Temp_pred[ ,-2], # delete Date
                           Temp_original  = Temp_orig,
                           Detection      = Det_Temp_pred)

# Plot predicted occupancy vs original elevation based on predictions with
# standardized elevation
Elev_gg <- ggplot(data = gg_data_elev, aes(x = Elev_orig, y = Occupancy.Predicted)) +
  geom_line(size = 2, color = "darkblue") +
  geom_ribbon(aes(ymin  = Occupancy.lower,
                  ymax  = Occupancy.upper, 
                  alpha = 0.1), fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Occupancy Probability of White-tailed deer vs Elevation",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Elevation [m]",
       y        = "Predicted Occupancy Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_continuous(limits = round(c(min(gg_data_elev$Elev_orig),
                                      max(gg_data_elev$Elev_orig))),
                     breaks = c(1550,1750,2000,2250))

# Plot predicted Occupancy vs original total cover based on predicitions with
# standardized total cover
Cover_gg <- ggplot(data = gg_data_cover, aes(x = Cover_orig, y = Occupancy.Predicted)) +
  geom_line(size = 2, color = "darkblue") +
  geom_ribbon(aes(ymin  = Occupancy.lower,
                  ymax  = Occupancy.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Occupancy Probability of White-tailed deer vs total Vegetation Cover",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Total Vegetation Cover",
       y        = "Predicted Occupancy Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none")

# Plot predicted detection probability vs date (Time of year) based on predictions
# with standardized dates
Date_gg <- ggplot(data = gg_data_date, aes(x = Date_orig, y = Detection.Predicted)) +
  geom_line(size = 2, color = "darkgreen") +
  geom_ribbon(aes(ymin  = Detection.lower,
                  ymax  = Detection.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Detection probability of White-tailed deer vs Time of Year",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Time of Year",
       y        = "Predicted Detection Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_continuous(breaks = c(250, 275, 306, 336),
                     labels = c("September", "October", "November",
                                "December"))

# Plot predicted occupancy probability vs water distance based on prediction
# with standardized water distance
Wat_gg <- ggplot(data = gg_data_water, aes(x = Water_orig, y = Occupancy.Predicted)) +
  geom_line(size = 2, color = "midnightblue") +
  geom_ribbon(aes(ymin  = Occupancy.lower,
                  ymax  = Occupancy.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Occupancy Probability of White-tailed deer vs Water distance",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Water distance [m]",
       y        = "Predicted Occupancy Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000)) +
  scale_y_continuous(limits = 0:1)

# Plot predicted occupancy probability vs distance to roads and trails based on 
# predictions with standardized distance to roads and trails
RoaTra_gg <- ggplot(data = gg_data_RoaTra, aes(y = Occupancy.Predicted)) +
  geom_line(aes(x = Trails_original), size = 2, color = "midnightblue") +
  geom_line(aes(x = Roads_original),  size = 2, color = "darkgreen") +
  geom_ribbon(aes(x     = Trails_original,
                  ymin  = Occupancy.lower,
                  ymax  = Occupancy.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  geom_ribbon(aes(x     = Roads_original,
                  ymin  = Occupancy.lower,
                  ymax  = Occupancy.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Occupancy Probability of White-tailed deer vs 
                               Distance to Roads and Trails",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Distance to Roads or Trails [m]",
       y        = "Predicted Occupancy Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none")

# Plot predicted detection probability vs temperature based on predicitons with
# standardized temperature
Temp_gg <- ggplot(data = gg_data_temp, aes(x = Temp_orig, y = Detection.Predicted)) +
  geom_line(size = 2, color = "darkgreen") +
  geom_ribbon(aes(ymin  = Detection.lower,
                  ymax  = Detection.upper,
                  alpha = 0.1),  fill = "lightgrey") +
  theme_bw() +
  labs(title    = "Detection probability of White-tailed deer vs Temperature",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Temperature [C]",
       y        = "Predicted Detection Probability") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        legend.position = "none")

# All plots together
ggarrange(nrow = 2, ncol = 3,
          Elev_gg, Wat_gg, Cover_gg, RoaTra_gg, Date_gg, Temp_gg)

# DISTRIBUTION MAPPING ---------------------------------------------------------

# Read in file that contains all XY coordinates from park of interest
# All coordinates have to be in UTM to make multiple layers fit in ggplot

# 30 m fishnet from CHIR
Fishnet30 <- data.frame(read.csv(file   = "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/Covariates/Fishnets/CHIR_Fishnet_30m_2.csv",
                                 header = TRUE))

# Make new data frame with Elevation data from fishnet that can be used for
# predict
Fishnet30_Elev           <- data.frame(Fishnet30$RASTERVALU) # 30 m
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

# Map_pred_30 <- predict(fm32, newdata = Fishnet30_zElev, type = "state")

# Save prediction of best model, so the code doesn't need to be executed again
# write.csv(x = Map_pred_30, file = "Map_pred_30_2.csv")
Map_pred_30 <- data.frame(read.csv(file = "Map_pred_30.csv")[ ,-1])

# Make gg data frame with predictions and original values
gg_data_map30 <- data.frame(X_coordinate = Fishnet30$POINT_X,
                            Y_coordinate = Fishnet30$POINT_Y,
                            Occupancy    = Map_pred_30$Predicted)

# Plot the distribution map
# Get height and width of cell
gg_data_map30$X_coordinate[[2]] - gg_data_map30$X_coordinate[[1]]

# Plotting the boundary
# Reading in the shapefile with the boundary
CHIR_Boundary <- readOGR(dsn = GeoDB, layer = "CHIR_Boundary")
CHIR_Boundary <- tidy(CHIR_Boundary)

ggplot() +
  geom_polygon(data = CHIR_Boundary, aes(x     = long,
                                         y     = lat,
                                         group = group),
               fill = NA, color = "black", size = 1.5) +
  theme_classic()

# Plotting
OccuMap_30 <- ggplot() +
  # Occupancy layer
  geom_tile(data = gg_data_map30, aes(x    = X_coordinate,
                                      y    = Y_coordinate,
                                      fill = Occupancy)) +
  # Boundary layer
  geom_polygon(data = CHIR_Boundary, aes(x     = long,
                                         y     = lat,
                                         group = group),
               fill = NA, color = "black", size = 2) +
  scale_fill_viridis(option = "viridis", direction = 1) +
  theme_bw() +
  rremove("grid") +
  labs(title    = "Occupancy probability of White-tailed Deer",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Easting",
       y        = "Northing") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5))

# Creating hillshade effect

# DEM
# Loading DEM from created 30 m Fishnet that is used for OCcupancy predictions
DEM <- data.frame(readOGR(dsn = GeoDB, layer = "CHIR_Fishnet_30m_Elevation"))
DEM <- (DEM[ ,-c(4:6)])
colnames(DEM) <- c("X", "Y", "Elevation")

# Plotting DEM
(gg_DEM <- ggplot(data = DEM, aes(x = X,
                                  y = Y)) +
    geom_raster(aes(fill = Elevation), alpha = 0.75) +
    scale_fill_gradientn(name = "Elevation", colours = terrain.colors(100)) +
    theme_bw() +
    rremove("grid") +
    coord_equal() +
    xlab("Longitude") +
    ylab("Latitude"))

# Making DEM to raster
DEM_raster <- rasterFromXYZ(DEM)
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

# Plot

# Just Hillshade map
gg_Hill <- ggplot() + 
  geom_tile(data = hill_df, aes(x = X_coord, y = Y_coord, fill = Value)) + 
  scale_fill_gradient(low = "black", high = "white") +
  new_scale_fill() +
  geom_tile(data = DEM, aes(x = X, y = Y, fill = Elevation), alpha=0.4) + 
  scale_fill_gradientn(colours = rev(terrain.colors(10))) +
  theme_bw() + 
  theme(legend.position="none")

# Combine Hillshade and Occupancy predictions
gg_OccuHill <- ggplot() +
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
  # Boundary layer
  geom_polygon(data = CHIR_Boundary, aes(x     = long,
                                         y     = lat,),
               fill = NA, color = "black", size = 1.5) +
  # Plot anatomy        
  theme_bw() + 
  rremove("grid") +
  labs(title    = "Occupancy probability of White-tailed Deer",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Easting",
       y        = "Northing") +
  theme(plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5))

# Plotting roads and trails
# Roads
CHIR_Roads <- readOGR(dsn = GeoDB, layer = "CHIR_Roads_Park")
CHIR_Roads <- tidy(CHIR_Roads)

# See if plot looks as it should
ggplot() +
  geom_path(data = CHIR_Roads, aes(x     = long,
                                   y     = lat,
                                   group = group),
            fill = NA, color = "black", size = 1) +
  theme_classic()

# Trails
CHIR_Trails <- readOGR(dsn = GeoDB, layer = "CHIR_Trails")
CHIR_Trails <- tidy(CHIR_Trails)

# See if plot looks as it should
ggplot() +
  geom_path(data = CHIR_Trails, aes(x     = long,
                                    y     = lat,
                                    group = group),
            fill = NA, color = "grey", size = 1) +
  theme_classic()



# All layers combined
gg_OccuHill_ALL <- ggplot() +
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
  # Roads layer
  geom_path(data = CHIR_Roads, aes(x     = long,
                                   y     = lat,
                                   group = group),
            fill = NA, color = "black", size = 1) +
  # Trails layer
  geom_path(data = CHIR_Trails, aes(x     = long,
                                    y     = lat,
                                    group = group),
            fill = NA, color = "dimgray", size = 0.7, linetype = "dashed") +
  # Boundary layer
  geom_polygon(data = CHIR_Boundary, aes(x     = long,
                                         y     = lat,),
               fill = NA, color = "black", size = 1.5) +
  # Plot anatomy        
  theme_bw() + 
  rremove("grid") +
  north(data = CHIR_Boundary, location = "topright", symbol = 3) +
  ggsn::scalebar(data      = CHIR_Boundary, location = "bottomleft", dist = 1,
                 dist_unit = "km", transform = FALSE) +
  labs(title    = "Occupancy probability of White-tailed Deer",
       subtitle = "Chiricahua National Monument 2017",
       x        = "Easting",
       y        = "Northing") +
  theme(plot.title       = element_text(hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5),
        legend.position  = "right",
        legend.box       = "vertical",
        legend.key.size  = unit(2, "cm"),
        legend.key.width = unit(0.8, "cm")) +
  guides(fill_new = FALSE, # remove hillshade from legend
         fill     = guide_colorbar(ticks.colour = "black",
                                   frame.colour = "black"))

# SAVING PLOTS AND MAPS --------------------------------------------------------

# Set working directory to the Output location
setwd("C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/R/Output/CHIR/ODVI")

# Covariate plots

# ggsave("CHIR_ODVI_2017_Occu_Elev.png", plot  = Elev_gg, width = 12, height = 9, 
#        units = "in", dpi = 600)
# ggsave("CHIR_ODVI_2017_Det_Date.png", plot  = Det_gg, width = 12, height = 9,
#        units = "in", dpi = 600)
# ggsave("CHIR_ODVI_2017_Occu_Wat.png", plot  = Wat_gg, width = 12, height = 9,
#        units = "in", dpi = 600)
# ggsave("CHIR_ODVI_2017_Occu_Cov.png", plot  = Cover_gg, width = 12, height = 9,
#        units = "in", dpi = 600)
# ggsave("CHIR_ODVI_2017_Occu_Roa+Tra.png", plot  = RoaTra_gg, width = 12, height = 9,
#        units = "in", dpi = 600)
# 
# # Distribution maps
# 
# ggsave("CHIR_ODVI_2017_Occu_Map_30.png",
#        plot   = OccuMap_30, width = 12, height = 9, units = "in", dpi = 600)
# ggsave("CHIR_ODVI_2017_OccuMap.png",
#        plot   = gg_OccuHill_ALL, width = 12, height = 9, units = "in", dpi = 600)



