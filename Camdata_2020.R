# CAMDATA ----------------------------------------------------------------------
# Alex Buckisch
# abuckisch@email.arizona.edu
# September 2020

# This script contains code to extract data from NPS ArcGIS .gdb files for 
# occupancy analysis. In this script there are references to the "super-study" 
# such as the "super-study period" and "super-study area". This references the 
# total number of sites and the time span over which cameras were deployed by 
# NPS, useful for organizing the data. In occupancy modeling, the spatial
# scale and temporal scale is usually much narrower.

# PACKAGES ---------------------------------------------------------------------

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

# SET WORKING DIRECTORY---------------------------------------------------------

setwd("C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/R/Models/CAMDATA")


# IMPORTING GEODATABASE FROM ARCGIS PRO ----------------------------------------

GeoDB <- "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/ArcGIS Pro/MyProject/SODN_Data.gdb"

# List all feature classes in a file geodatabase
fc_list <- ogrListLayers(dsn = GeoDB)

# IMPORTING AND FORMATTING LOCATIONS, DEPLOYMENTS, AND DETECTIONS --------------

# Read in the feature class that contains all camera locations (location data)
cam_locations_df <- as.data.frame(readOGR(dsn = GeoDB, layer = "NEW_ECO_WildlifeCameraLocations"))

# Delete everything before 2016 (LCNWR, SBNWR, and CHIR 2014/2016
# side note: TONT 2017 used the same uplands locations for the herpetofauna
# deployment, won't get deleted)
cam_locations_final  <- cam_locations_df

# Read in the feature class that contains all deployment windows (event data)
cam_deployments_df <- as.data.frame(readOGR(dsn = GeoDB,layer= "NEW_ECO_WildlifeCameraEvents"))

# Delete everything before 2016
# (LCNWR, SBNWR, and CHIR 2014/2015/2016) as well as TONT 2017
cam_deployments_sorted <- cam_deployments_df %>% arrange(DeployDate)

# delete all LCNWR, SBNWR, and CHIR 2014
cam_deployments_sorted <- cam_deployments_sorted[-c(1:73), ]             
cam_deployments_sorted <- cam_deployments_sorted %>% arrange(DeployDate)

# delete all CHIR 2016
cam_deployments_sorted <- cam_deployments_sorted[-c(61:84), ]            
cam_deployments_sorted <- cam_deployments_sorted %>% arrange(DeployDate)

# delete all TONT 2017 (Herpetofauna cams)
cam_deployments_sorted <- cam_deployments_sorted[-c(182:205), ]          
cam_deployments_final  <- cam_deployments_sorted %>% arrange(DeployDate)

# Read in the feature class that contains all species detections (species data)
cam_detections <- sf::st_read(dsn = GeoDB, layer = "NEW_WildlifeSpeciesData")

# Delete everything before 2016 (LCNWR, SBNWR, and CHIR 2014/2016)
# Delete all LCNWR, SBNWR, TONT FY17
cam_detections_sorted <- cam_detections[- grep("LCNWR|SBNWR|Herpetofauna",
                                               cam_detections$Filepath), ]
# Delete all CHIR 2014/2016
cam_detections_sorted <- cam_detections_sorted[- (1:3914), ]

# Add a new StdLocName to detections and match it with the corresponding plot
# (based on EventID) 
cam_detections_sorted$StdLocName <- cam_deployments_sorted$StdLocName[match(x     = cam_detections_sorted$Event_ID,
                                                                            table = cam_deployments_sorted$ID)]

# Reorder columns in detections so StdLocName is not at the end for better readability
cam_detections_sorted <- cam_detections_sorted[, c(1,2,7,3,4,5,6)]

# List all Species codes detected (93 total species detected)
codes           <- as.data.frame(unique(x = cam_detections_sorted$Species_Code)) 
colnames(codes) <- "Species_Code"

# Remove all bird species, reptile species, unknown species, rodent species 
# (except squirrels and woodrats), dogs, no species, i.e. humans + setup photos

# Birds: AMRO|BEWR|BHCO|BTSP|CACW|CAKI|CANT|CANW|CASJ|CBTH|CORA|COYE|DEJU|GAQU|
#        GBHE|GHOW|GRRO|LEGO|MEJA|MODO|MONQ|NOFL|NOMO|RCSP|ROWR|RTHA|SPTO|STJA|
#        THCY|TUVU|WESO|WITU
# Mammals: AMHA|CAFA|DIDE|DIME|SPTE|TADO
# Other: CNTI|CRCO|GOAG|NOSP|PICA|SCCL|THCY
# Unknowns: UNBI|UNCA|UNCE|UNCO|UNCR|UNDI|UNDO|UNHU|UNJA|UNKI|UNLE|UNME1|UNNE|
#           UNOD|UNPE|UNRA|UNRO|UNSQ|UNWO (don't delete UNSY)



# 29 occupancy relevant species remain after elimination
cam_detections_final <- cam_detections_sorted[- grep("AMRO|BEWR|BHCO|BTSP|CACW|
                                                     |CAKI|CANT|CANW|CASJ|CBTH|
                                                     |CORA|COYE|DEJU|GAQU|GBHE|
                                                     |GHOW|GRRO|LEGO|MEJA|MODO|
                                                     |MONQ|NOFL|NOMO|RCSP|ROWR|
                                                     |RTHA|SPTO|STJA|THCY|TUVU|
                                                     |WESO|WITU|AMHA|CAFA|DIDE|
                                                     |DIME|SPTE|TADO|CNTI|CRCO|
                                                     |GOAG|NOSP|PICA|SCCL|THCY|
                                                     |UNBI|UNCA|UNCE|UNCO|UNCR|
                                                     |UNDI|UNDO|UNHU|UNJA|UNKI|
                                                     |UNLE|UNME1|UNNE|UNOD|UNPE|
                                                     |UNRA|UNRO|UNSQ|UNWO",
                                                     cam_detections_sorted$Species_Code), ]                   

# Replace outlier COME with COLE (hog-nosed skunk wrong species code)
cam_detections_final$Species_Code <- gsub(pattern = "COME",
                                          replacement = "COLE",
                                          x = cam_detections_final$Species_Code)

# Check if detection matrix now only contains the right species
codes           <- as.data.frame(unique(x = cam_detections_final$Species_Code))
colnames(codes) <- "Species_Code"

# For consistency and simplicity, it is best to rename each StdLocName 
# with a simple unique number, called "SiteID", and store this in a table for 
# cross-referencing
SiteID_names <- data.frame(StdLocName = unique(cam_deployments_final$StdLocName),
                           SiteID = seq(1:length(unique(cam_deployments_final$StdLocName))))

# Add a "SiteID" column to locations, deployments, and detections
cam_locations_final$SiteID   <- SiteID_names$SiteID[match(x = cam_locations_final$StdLocName, table = SiteID_names$StdLocName)]
cam_deployments_final$SiteID <- SiteID_names$SiteID[match(x = cam_deployments_final$StdLocName, table = SiteID_names$StdLocName)]
cam_detections_final$SiteID  <- SiteID_names$SiteID[match(x = cam_detections_final$StdLocName, table = SiteID_names$StdLocName)]

# Add a "UnitCode" column to deployments and detections
cam_deployments_final$UnitCode <- cam_locations_final$UnitCode[match(x = cam_deployments_final$StdLocName, table = cam_locations_final$StdLocName)]
cam_detections_final$UnitCode  <- cam_locations_final$UnitCode[match(x = cam_detections_final$StdLocName, table = cam_locations_final$StdLocName)]

# DATE CONVERSION --------------------------------------------------------------

# Convert the columns with date data into dates that R can recognize
cam_deployments_final$DeployDate    <- as_date(ymd_hms(as.character(cam_deployments_final$DeployDate)))
cam_deployments_final$RetrievalDate <- as_date(ymd_hms(as.character(cam_deployments_final$RetrievalDate)))


# Calculate super-study days of operation (first delployment to last retrieval)
cam_deployments_final$Study_days     <- difftime(time1 = cam_deployments_final$RetrievalDate,
                                                 time2 = cam_deployments_final$DeployDate, units = "days") + 1
cam_deployments_final$DeployJulian   <- yday(cam_deployments_final$DeployDate)
cam_deployments_final$RetrieveJulian <- yday(cam_deployments_final$RetrievalDate)

# Identify first and last day of super-study
day_first  <- min(cam_deployments_final$DeployDate) 
day_last   <- max(cam_deployments_final$RetrievalDate)

# Add seperate date and time columns to detections (easier to analyze ymd and hms than ymd_hms)
cam_detections_final$Date <- as_date(ymd_hms(as.character(cam_detections_final$DateTime_Trig)))              
cam_detections_final$Time <- strftime(as.character(cam_detections_final$DateTime_Trig), format = "%H:%M:%S")

# Add the super-study day and Julian day to detections (useful for analysis)
cam_detections_final$Study_day <- difftime(time1 = cam_detections_final$Date,
                                           time2 = day_first, units = "days") + 1
cam_detections_final$JulianDay <- yday(cam_detections_final$DateTime_Trig)

# Add the super-study day of Deployment date and retrieval date to deployments
cam_deployments_final$Study_day_Deploy   <- difftime(time1 = cam_deployments_final$DeployDate,
                                                     time2 = day_first, units = "days") + 1
cam_deployments_final$Study_day_Retrieve <- difftime(time1 = cam_deployments_final$RetrievalDate,
                                                     time2 = day_first, units = "days") + 1

# FINAL FORMATTING FOR LOCATIONS, DEPLOYMENTS, AND DETECTIONS ------------------

cam_locations_final   <- cam_locations_final[ ,c(1,2,22,3:16)]
cam_deployments_final <- cam_deployments_final[ ,c(30,2,29,4,34,5,35,32,33,31,6:28)]
cam_detections_final  <- cam_detections_final[ ,c(1,2,10,11,12,13,8,9,3,4,5,6,7)]

# SELECTING DATA BASED ON DATE, PARK, ETC. -------------------------------------

# Select all rows of data frame that were in 2018 for total number of photos
cam_detections_2018 <- data.frame(cam_detections_sorted[grep("2018", cam_detections_sorted$DateTime_Trig), ])

# Now remove all species that are not of interest
cam_detections_2018 <- cam_detections_2018[- grep("AMRO|BEWR|BHCO|BTSP|CACW|
                                                     |CAKI|CANT|CANW|CASJ|CBTH|
                                                     |CORA|COYE|DEJU|GAQU|GBHE|
                                                     |GHOW|GRRO|LEGO|MEJA|MODO|
                                                     |MONQ|NOFL|NOMO|RCSP|ROWR|
                                                     |RTHA|SPTO|STJA|THCY|TUVU|
                                                     |WESO|WITU|AMHA|CAFA|
                                                     |SPTE|CNTI|CRCO|
                                                     |GOAG|NOSP|PICA|SCCL|THCY|
                                                     |UNBI|UNCA|UNCE|UNCO|UNCR|
                                                     |UNDI|UNDO|UNHU|UNJA|UNKI|
                                                     |UNLE|UNME1|UNNE|UNOD|UNPE|
                                                     |UNRA|UNRO|UNSQ|UNWO",
                                                  cam_detections_2018$Species_Code), ]

# Add "Unit Code" to detections 2018
cam_detections_2018$UnitCode  <- cam_locations_final$UnitCode[match(x = cam_detections_2018$StdLocName,
                                                                    table = cam_locations_final$StdLocName)]

# List all unique species
data.frame(unique(cam_detections_2018$Species_Code))

# Count how many detections of each species
table(cam_detections_2018$Species_Code)

# Transpose and make into a data frame
Species_detections_2018 <- data.frame(t(table(cam_detections_2018$Species_Code)))

# Delete first column and rename
Species_detections_2018 <- Species_detections_2018[ ,-1]
colnames(Species_detections_2018) <- c("Species_Code", "Detections_in_2018")

# Keep only 19 identified target species
Target_detections_2018 <- Species_detections_2018[- grep("LEAL|OVCA|MEMA|NEAL|
                                                         |SCAR|SPGR|TADO|TATA",
                                                  Species_detections_2018$Species_Code), ]

# Add all detections of target species
sum(Target_detections_2018$Detections_in_2018)

# Save data frame as table
# write.csv(x = Target_detections_2018, file = "C:/Users/Alex/Dropbox/Privat/University of Arizona/THESIS/R/Miscellaneous/Target_Species_Detections_2018.csv")

# Park specific detection data frames for 2018
CAGR_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "CAGR", ]
CHIR_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "CHIR", ]
GICL_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "GICL", ]
MOCC_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "MOCC", ]
MOWE_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "MOWE", ]
MOCA_2018 <- rbind(MOCC_2018, MOWE_2018)
ORPI_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "ORPI", ]
SAGW_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "SAGW", ]
TONT_2018 <- cam_detections_2018[cam_detections_2018$UnitCode == "TONT", ]

# Get species detected at each park and only focus on target species
sort(unique(TONT_2018$Species_Code))

# Count how many detections of each species at specific park unit
nrow(data.frame(cam_detections_2018[cam_detections_2018$UnitCode == "TONT" &
                                      cam_detections_2018$Species_Code == "PETA", ]))


# TESTING SELECT QUERIES -------------------------------------------------------

# Select all detections of ODVI at specific plot
as.data.frame(table(filter(cam_detections_2018[cam_detections_2018$Species_Code == "ODVI" &
                           cam_detections_2018$StdLocName   == "Wildlife_ORPI_V102_005", ])))


# CHIR 2019 June - October
CHIR_2019_June <- data.frame(cam_deployments_final %>%
  #select(UnitCode, StdLocName, DeployDate, RetrievalDate, TotalPics) %>%
  filter(between(DeployDate, as.Date("2019-06-17"), as.Date("2019-06-21"))))

# Total number of pictures taken
sum(as.integer(CHIR_2019_June$TotalPics))

# CHIR 2019 May and June
Det_MayJune <- data.frame(CHIR_2019 %>%
                            filter(between(DateTime_Trig, as.POSIXct("2019-05-01"),as.POSIXct("2019-06-30"))))

# Select just May
df <- Det_MayJune %>% filter(between(DateTime_Trig, as.POSIXct("2019-05-01"),as.POSIXct("2019-05-31")))

unique(df$StdLocName)
rowSums(df)
