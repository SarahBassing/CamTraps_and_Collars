  #'  ============================================
  #'  Occupancy Models (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  Februray 2021
  #'  ============================================
  #'  Script to create unmarked data frames and run single-species, single-season
  #'  occupancy models for deer, elk, cougars, wolves, coyotes, and bobcats for
  #'  summer 2018 (7/1/18 - 9/29/18) and winter 2018-2019 (12/1/18 - 3/1/19),
  #'  respectively. Each single-season occupancy model includes 13 7-day sampling 
  #'  occasions comprising the warmest months in summer and coldest months in 
  #'  winter with the most consistent snow.
  #'  
  #'  Encounter histories are generated with the CameratTrap_DetectionHistories.R
  #'  script. Covariate data included in occupancy models were collected at each
  #'  camera site or extracted from remotely sensed data.
  #'  ============================================

  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(unmarked)
  library(tidyverse)

  #'  Source script that generates detection histories
  #'  Detection histories come trimmed to desired season length
  source("./Scripts/CameraTrap_DetectionHistories.R")
  #'  Alternatively, can read in specific detection histories, e.g., 
  # DH_coug_smr18 <- read.csv("./Data/Detection_Histories/Cougar__detection_history__with_effort__7_days_per_occasion__occasionStart0h__first_day_2018-07-01__2021-02-01.csv")
  # DH_coug_wtr1819 <- read.csv("./Data/Detection_Histories/Cougar__detection_history__with_effort__7_days_per_occasion__occasionStart0h__first_day_2018-12-01__2021-02-01.csv")
  
  #'  Read in covariate data collected during camera deployment & scale
  #'  Canopy cover, land management & owner, habitat type => site-level occ covs
  #'  Dist. to focal pt, height, monitoring = > survey/site-level detection covs
  # stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18_2021-02-15.csv") %>%
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18-20_2021-03-04.csv") %>%
    #'  Get rid of mysterious space after one of the NEs
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    #'  Consolidate Monitoring feature (determined at camera site) into fewer 
    #'  categories that capture whether we expect vehicles to be on the feature
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Trail", "Game trail", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Decommissioned road", "Game trail", as.character(Monitoring))
    ) %>%
    #'  Consolidate Habitat Type into fewer categories (determined at camera site)
    mutate(
      Habitat_Type = ifelse(Habitat_Type == "Agriculture", "Grassland", as.character(Habitat_Type)), # Grassland / Agriculture
      Habitat_Type = ifelse(Habitat_Type == "Riparian", "Mixed conifer", as.character(Habitat_Type)) # Mixed conifer / Riparian
    ) %>%
    #'  Consolidate NLCD and landcover into fewer categories (determined via remote sensing)
    #'  NLCD 21: Developed- open space, 22: Developed- low intensity, 
    #'  42: Evergreen forest, 52: Shrub/scrub, 71: Grassland, 81: Pasture/hay, 
    #'  82: Cultivated crops, 90: Woody wetlands, 95: Emergent herbaceous wetlands
    mutate(
      NLCD_landcov = ifelse(NLCD_landcov == "21", "Developed", as.character(NLCD_landcov)), # Developed- open land / low intensity use
      NLCD_landcov = ifelse(NLCD_landcov == "22", "Developed", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "42", "Forested", as.character(NLCD_landcov)), # Evergreen forest / Woody wetland
      NLCD_landcov = ifelse(NLCD_landcov == "90", "Forested", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "52", "Shrub", as.character(NLCD_landcov)),
      NLCD_landcov = ifelse(NLCD_landcov == "71", "Ag-Grassland", as.character(NLCD_landcov)), # Grassland / Agriculture / Emergent herbaceous wetland
      NLCD_landcov = ifelse(NLCD_landcov == "95", "Ag-Grassland", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "81", "Ag-Grassland", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "82", "Ag-Grassland", as.character(NLCD_landcov)),
      # NLCD_landcov = ifelse(NLCD_landcov == "22", "21", as.character(NLCD_landcov)), # Developed open/low intensity human use
      # NLCD_landcov = ifelse(NLCD_landcov == "90", "42", as.character(NLCD_landcov)), # Evergreen forest/Woody wetland
      # NLCD_landcov = ifelse(NLCD_landcov == "95", "71", as.character(NLCD_landcov)), # Grassland / Agriculture / Emergent herbaceous wetland
      # NLCD_landcov = ifelse(NLCD_landcov == "81", "71", as.character(NLCD_landcov)), 
      # NLCD_landcov = ifelse(NLCD_landcov == "82", "71", as.character(NLCD_landcov)),
    #'  121: barren, 201: emergent wetland, 211: mesic grass, 212: xeric grass,
    #'  221: mesic shrub, 222: xeric shrub, 230: forest, 310: agriculture, 332: residential
      landcov18 = ifelse(landcov18 == "212", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "310", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "332", "211", as.character(landcov18)),  # Grassland / Agriculture / Residential
      landcov18 = ifelse(landcov18 == "222", "221", as.character(landcov18)),  # Shrub
      landcov19 = ifelse(landcov19 == "121", "211", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "201", "211", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "212", "211", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "310", "211", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "332", "211", as.character(landcov19)),  # Grassland / Agriculture / Residential / Other
      landcov19 = ifelse(landcov19 == "222", "221", as.character(landcov19)),  # Shrub
      landcov = ifelse(landcov == "211", "Ag-Grassland", as.character(landcov)),
      landcov = ifelse(landcov == "121", "Ag-Grassland", as.character(landcov)),
      landcov = ifelse(landcov == "201", "Ag-Grassland", as.character(landcov)), 
      landcov = ifelse(landcov == "212", "Ag-Grassland", as.character(landcov)),
      landcov = ifelse(landcov == "310", "Ag-Grassland", as.character(landcov)),
      landcov = ifelse(landcov == "332", "Ag-Grassland", as.character(landcov)),  # Grassland / Emergent wetland / Other: Barren, Agriculture, Residential
      landcov = ifelse(landcov == "221", "Shrub", as.character(landcov)),   # Shrub
      landcov = ifelse(landcov == "222", "Shrub", as.character(landcov)),
      landcov = ifelse(landcov == "230", "Forest", as.character(landcov))   # Forest
      # landcov = ifelse(landcov == "121", "211", as.character(landcov)),
      # landcov = ifelse(landcov == "201", "211", as.character(landcov)), 
      # landcov = ifelse(landcov == "212", "211", as.character(landcov)),
      # landcov = ifelse(landcov == "310", "211", as.character(landcov)),
      # landcov = ifelse(landcov == "332", "211", as.character(landcov)),  # Grassland / Emergent wetland / Other: Barren, Agriculture, Residential
      # landcov = ifelse(landcov == "222", "221", as.character(landcov))   # Shrub
    ) %>%
    #'  Rename, format, and scale as needed
    transmute(
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      # Cell_ID = as.factor(Cell_ID),
      # Camera_ID = as.factor(Camera_ID),
      CameraLocation = as.factor(CameraLocation),
      Distance = scale(Distance_Focal_Point),
      Height = scale(Height_frm_grnd),
      Trail = as.factor(Monitoring),
      Canopy_Cov = scale(Canopy_Cov),           # Fine-scale camera station
      Land_Mgnt = as.factor(Land_Mgnt),
      Land_Owner = as.factor(Land_Owner),
      Habitat_Type = as.factor(Habitat_Type),   # Fine-scale camera station
      NDVI_sp18 = scale(ndvi_sp18),             #  KEEP IN MIND:
      NDVI_sm18 = scale(ndvi_sm18),             # Scaling annual vs all values 
      NDVI_sp19 = scale(ndvi_sp19),             # together changes the values
      NDVI_sm19 = scale(ndvi_sm19),
      NDVI_sp = scale(ndvi_sp),
      NDVI_sm = scale(ndvi_sm),
      dNBR_sp18 = scale(dnbr_sp18),
      dNBR_sm18 = scale(dnbr_sm18),
      dNBR_sp19 = scale(dnbr_sp19),
      dNBR_sm19 = scale(dnbr_sm19),
      dNBR_sp = scale(dnbr_sp),
      dNBR_sm = scale(dnbr_sm),
      Disturb18 = as.factor(disturb18),
      Disturb19 = as.factor(disturb19),
      Disturb = as.factor(disturb),
      Burn18 = as.factor(burnPerim18),
      Burn19 = as.factor(burnPerim19),
      Burn = as.factor(burnPerim),
      Landcov18 = as.factor(landcov18),         # Courser-scale, 30m res
      Landcov19 = as.factor(landcov19),
      Landcov = as.factor(landcov),
      NLCD = as.factor(NLCD_landcov),           # Courser-scale, 30m res
      Elev = scale(elevation),
      Slope = scale(slope), # this is a circular variable so not quite sure how to treat this...
      Aspect = scale(aspect), # this is a circular variable and 90degrees used when slope = 0
      TRI = scale(tri),
      Roughness = scale(roughness),
      Canopy18 = scale(canopy18),               # Courser-scale, 30m res
      Canopy19 = scale(canopy19),
      Canopy = scale(canopy),
      NearestH2o = scale(km2water)
    )
  #WaterDen = scale(water_density)
  
  #'  Identify which sites have missing data (just look at 1st column)
  noCanopy <- which(is.na(stations$Canopy_Cov))
  noDist <- which(is.na(stations$Distance))
  noHgt <- which(is.na(stations$Height))
  missing_dat <- unique(c(noHgt, noDist, noCanopy))
  
  #'  How many detections am I losing as result of this?
  sum(DH_bob_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_bob_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_coug_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_coug_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_coy_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_coy_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_wolf_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_wolf_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_elk_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_elk_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_md_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_md_wtr1819[missing_dat,], na.rm = TRUE)
  sum(DH_wtd_smr18[missing_dat,], na.rm = TRUE)
  sum(DH_wtd_wtr1819[missing_dat,], na.rm = TRUE)
  #'  Bummer dude 
  
  #'  Find median value of these covaraites (mean = 0, sd = 1 since center & scaled)
  median(stations$Distance, na.rm = TRUE)
  median(stations$Height, na.rm = TRUE)
  median(stations$Canopy_Cov, na.rm = TRUE) # data skewed due to many cameras in open/burned areas
  
  #'  Replace missing covariate values with mean value
  #'  Missing 2018 Canopy_Cov obs in forested are so using mean instead of median
  # stations[is.na(stations$Distance),] <- 0
  # stations[is.na(stations$Height),] <- 0
  # stations[is.na(stations$Canopy_Cov),] <- 0
    
  #'  Check for correlation among covaraites
  #'  Watch out for NAs (use="complete.obs")
  cor(stations$Distance, stations$Height, use = "complete.obs")
  cor(stations$Elev, stations$Slope, use = "complete.obs")
  cor(stations$Elev, stations$Aspect, use = "complete.obs")
  cor(stations$Slope, stations$Aspect, use = "complete.obs")
  cor(stations$Elev, stations$TRI, use = "complete.obs")
  cor(stations$TRI, stations$Roughness, use = "complete.obs")  #  EEK!
  cor(stations$Elev, stations$NDVI_sm18, use = "complete.obs")
  cor(stations$Slope, stations$NDVI_sm18, use = "complete.obs")
  cor(stations$Aspect, stations$NDVI_sm18, use = "complete.obs")
  
  #'  Create survey-level covariate matrix
  #'  Requires uniqe column for each sampling occasion and covariate
  nrows <- nrow(stations)
  ncols <- 13
 
  srvy_covs <- list(
    Height = matrix(c(Hgt1 = stations$Height, Hgt2 = stations$Height,
                      Hgt3 = stations$Height, Hgt4 = stations$Height,
                      Hgt5 = stations$Height, Hgt6 = stations$Height,
                      Hgt7 = stations$Height, Hgt8 = stations$Height,
                      Hgt9 = stations$Height, Hgt10 = stations$Height,
                      Hgt11 = stations$Height, Hgt12 = stations$Height,
                      Hgt13 = stations$Height),
                    nrow = nrows, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = stations$Distance, Dist2 = stations$Distance,
                        Dist3 = stations$Distance, Dist4 = stations$Distance,
                        Dist5 = stations$Distance, Dist6 = stations$Distance,
                        Dist7 = stations$Distance, Dist8 = stations$Distance,
                        Dist9 = stations$Distance, Dist10 = stations$Distance,
                        Dist11 = stations$Distance, Dist12 = stations$Distance,
                        Dist13 = stations$Distance),
                      nrow = nrows, ncol = ncols, byrow = FALSE)
    )
  
  #'  Double check it looks OK
  head(srvy_covs[[2]])
  
  #'  NEED TO BRING IN EFFORT DATA AS WELL!
  


  
  #'  Create unmarked dataframes
  ####  BOBCAT UMF  ####
  bob_s1819_UMF <- unmarkedFrameOccu(DH_bob_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Landcov = stations$Landcov,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm,  
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         Tree_cov = stations$Canopy,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  # Mgnt = stations$Land_Mgnt,
  # Habitat = stations$Habitat_Type,
  # dNBR = stations$dNBR_sm,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
  # TRI = stations$TRI,
  
  #'  Remove rows with missing observation covariate data (Height & Distance data)
  nrow(bob_s1819_UMF@y)
  bob_s1819_UMF <- bob_s1819_UMF[-missing_dat]
  nrow(bob_s1819_UMF@y)
  ##### DOES IT MAKE SENSE TO DROP ENTIRE CAMERA SITE DUE TO 1 MISSING SITE COVARIATE THAT OFTEN ISN'T EVEN SIGNIFICANT?

  bob_w1820_UMF <- unmarkedFrameOccu(DH_bob_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Landcov = stations$Landcov,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           Tree_cov = stations$Canopy,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)

  #'  Remove rows with missing observation covariate data (Height & Distance data)
  nrow(bob_w1820_UMF@y)
  bob_w1820_UMF <- bob_w1820_UMF[-missing_dat]
  nrow(bob_w1820_UMF@y)
  
  summary(bob_w1820_UMF)
  
  ####  COUGAR UMF  ####
  coug_s1819_UMF <- unmarkedFrameOccu(DH_coug_smr1819,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Landcov = stations$Landcov,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm,  
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  nrow(coug_s1819_UMF@y)
  coug_s1819_UMF <- coug_s1819_UMF[-missing_dat]
  nrow(coug_s1819_UMF@y)
  
  coug_w1820_UMF <- unmarkedFrameOccu(DH_coug_wtr1820,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Landcov = stations$Landcov,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  nrow(coug_w1820_UMF@y)
  coug_w1820_UMF <- coug_w1820_UMF[-missing_dat]
  nrow(coug_w1820_UMF@y)
  
  # coug_s19_UMF <- unmarkedFrameOccu(DH_coug_smr19,
  #                                   siteCovs = data.frame(Year = stations$Year,
  #                                                         Area = stations$Study_Area,
  #                                                         Trail = stations$Trail,
  #                                                         Canopy_cov = stations$Canopy_Cov,
  #                                                         Landcov = stations$Landcov,
  #                                                         NLCD = stations$NLCD,
  #                                                         NDVI = stations$NDVI_sm,  
  #                                                         Elev = stations$Elev,
  #                                                         Slope = stations$Slope,
  #                                                         Aspect = stations$Aspect,
  #                                                         TRI = stations$TRI,
  #                                                         Tree_cov = stations$Canopy,
  #                                                         NearestH2o = stations$NearestH2o),
  #                                   obsCovs = srvy_covs)
  # nrow(coug_s19_UMF@y)
  # coug_s19_UMF <- coug_s19_UMF[-missing_dat]
  # nrow(coug_s19_UMF@y)
  # coug_w1920_UMF <- unmarkedFrameOccu(DH_coug_wtr1920,
  #                                     siteCovs = data.frame(Year = stations$Year,
  #                                                           Area = stations$Study_Area,
  #                                                           Trail = stations$Trail,
  #                                                           Canopy_cov = stations$Canopy_Cov,
  #                                                           Landcov = stations$Landcov,
  #                                                           NLCD = stations$NLCD,
  #                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
  #                                                           Elev = stations$Elev,
  #                                                           Slope = stations$Slope,
  #                                                           Aspect = stations$Aspect,
  #                                                           TRI = stations$TRI,
  #                                                           Tree_cov = stations$Canopy,
  #                                                           NearestH2o = stations$NearestH2o),
  #                                     obsCovs = srvy_covs)
  # nrow(coug_w1920_UMF@y)
  # coug_w1920_UMF <- coug_w1920_UMF[-missing_dat]
  # nrow(coug_w1920_UMF@y)
  
  ####  COYOTE UMF  ####
  coy_s1819_UMF <- unmarkedFrameOccu(DH_coy_smr1819,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Landcov = stations$Landcov,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm,  
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          Tree_cov = stations$Canopy,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  coy_s1819_UMF <- coy_s1819_UMF[-missing_dat]
  
  coy_w1820_UMF <- unmarkedFrameOccu(DH_coy_wtr1820,
                                      siteCovs = data.frame(Year = stations$Year,
                                                            Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            Canopy_cov = stations$Canopy_Cov,
                                                            Landcov = stations$Landcov,
                                                            NLCD = stations$NLCD,
                                                            NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            Aspect = stations$Aspect,
                                                            Tree_cov = stations$Canopy,
                                                            NearestH2o = stations$NearestH2o),
                                      obsCovs = srvy_covs)
  coy_w1820_UMF <- coy_w1820_UMF[-missing_dat]
  
  ####  WOLF UMF  ####
  wolf_s1819_UMF <- unmarkedFrameOccu(DH_wolf_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Landcov = stations$Landcov,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm,  
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         Tree_cov = stations$Canopy,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  wolf_s1819_UMF <- wolf_s1819_UMF[-missing_dat]
  
  wolf_w1820_UMF <- unmarkedFrameOccu(DH_wolf_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Landcov = stations$Landcov,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           Tree_cov = stations$Canopy,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)
  wolf_w1820_UMF <- wolf_w1820_UMF[-missing_dat]
  
  
  ####  ELK UMF  ####
  #'  Consider removing observations/sites from the OK since so few detections and no collars there
  elk_s1819_UMF <- unmarkedFrameOccu(DH_elk_smr1819,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Landcov = stations$Landcov,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm,  
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          Tree_cov = stations$Canopy,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  elk_s1819_UMF <- elk_s1819_UMF[-missing_dat]
  
  elk_w1820_UMF <- unmarkedFrameOccu(DH_elk_wtr1820,
                                      siteCovs = data.frame(Year = stations$Year,
                                                            Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            Canopy_cov = stations$Canopy_Cov,
                                                            Landcov = stations$Landcov,
                                                            NLCD = stations$NLCD,
                                                            NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            Aspect = stations$Aspect,
                                                            Tree_cov = stations$Canopy,
                                                            NearestH2o = stations$NearestH2o),
                                      obsCovs = srvy_covs)
  elk_w1820_UMF <- elk_w1820_UMF[-missing_dat]
  
  ####  MULE DEER UMF  ####
  md_s1819_UMF <- unmarkedFrameOccu(DH_md_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Landcov = stations$Landcov,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm,  
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         Tree_cov = stations$Canopy,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  md_s1819_UMF <- md_s1819_UMF[-missing_dat]
  
  md_w1820_UMF <- unmarkedFrameOccu(DH_md_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Landcov = stations$Landcov,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           Tree_cov = stations$Canopy,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)
  md_w1820_UMF <- md_w1820_UMF[-missing_dat]
  
  ####  WHITE-TAILED DEER UMF  ####
  wtd_s1819_UMF <- unmarkedFrameOccu(DH_wtd_smr1819,
                                  siteCovs = data.frame(Year = stations$Year,
                                                        Area = stations$Study_Area,
                                                        Trail = stations$Trail,
                                                        Canopy_cov = stations$Canopy_Cov,
                                                        Landcov = stations$Landcov,
                                                        NLCD = stations$NLCD,
                                                        NDVI = stations$NDVI_sm,  
                                                        Elev = stations$Elev,
                                                        Slope = stations$Slope,
                                                        Aspect = stations$Aspect,
                                                        Tree_cov = stations$Canopy,
                                                        NearestH2o = stations$NearestH2o),
                                  obsCovs = srvy_covs)
  wtd_s1819_UMF <- wtd_s1819_UMF[-missing_dat]
  
  wtd_w1820_UMF <- unmarkedFrameOccu(DH_wtd_wtr1820,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Landcov = stations$Landcov,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          Tree_cov = stations$Canopy,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  wtd_w1820_UMF <- wtd_w1820_UMF[-missing_dat]
  
  
  #'  Occupancy models
  #'  =============================
  #'  unmarked formula: ~detection ~occupancy
  #'  ~1 for intercept only
  #'  
  #'  Use chi-sq test to evaluate model fit after model selection (pg. 4 vignette)


  ####  BOBCAT MODELS  ####                   NLCD more signif than Landcov
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (bob_s1819_null <- occu(~1 ~1, bob_s1819_UMF))
  backTransform(bob_s1819_null, 'det')
  backTransform(bob_s1819_null, 'state')
  #'  Null with detection covariates
  (bob_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, bob_s1819_UMF))
  #'  Terrain model
  (bob_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, bob_s1819_UMF))
  (bob_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, bob_s1819_UMF))
  #'  Vegetation model
  (bob_s1819_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, bob_s1819_UMF))
  #'  Habitat model
  (bob_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, bob_s1819_UMF))
  (bob_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, bob_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(bob_s1819_null, bob_s1819_null2, bob_s1819_terrain, bob_s1819_terrain2, bob_s1819_veg, bob_s1819_terrain_veg, bob_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Null
  (bob_w1820_null <- occu(~1 ~1, bob_w1820_UMF))
  backTransform(bob_w1820_null, 'det')
  backTransform(bob_w1820_null, 'state')
  #'  Null with detection covariates
  (bob_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, bob_w1820_UMF))
  #'  Terrain model
  (bob_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, bob_w1820_UMF))
  (bob_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, bob_w1820_UMF))
  #'  Vegetation model
  (bob_w1820_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, bob_w1820_UMF))  
  #'  Habitat model
  (bob_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, bob_w1820_UMF))
  (bob_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, bob_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(bob_w1820_null, bob_w1820_null2, bob_w1820_terrain, bob_w1820_terrain2, bob_w1820_veg, bob_w1820_terrain_veg, bob_w1820_terrain_veg2)
  modSel(mods)
  
  ####  COUGAR MODELS  ####
  #'  Included covariates informed by Dickson & Beier 2002, Kertson et al. 2011, 
  #'  Smereka et al. 2020
                                              #  No real diff btwn NLCD & landcov in summer
  #'  SUMMERS 2018 & 2019
  #'  Null
  (coug_s1819_null <- occu(~1 ~1, coug_s1819_UMF))
  backTransform(coug_s1819_null, 'det')
  backTransform(coug_s1819_null, 'state')
  #'  Null with detection covariates
  (coug_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, coug_s1819_UMF))
  #'  Terrain model
  (coug_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, coug_s1819_UMF))
  (coug_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, coug_s1819_UMF))
  #'  Vegetation model
  (coug_s1819_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, coug_s1819_UMF))
  #'  Habitat model
  (coug_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, coug_s1819_UMF))
  (coug_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, coug_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model

  mods <- fitList(coug_s1819_null, coug_s1819_null2, coug_s1819_terrain, coug_s1819_terrain2, coug_s1819_veg, coug_s1819_terrain_veg, coug_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Null
  (coug_w1820_null <- occu(~1 ~1, coug_w1820_UMF))
  backTransform(coug_w1820_null, 'det')
  backTransform(coug_w1820_null, 'state')
  #'  Null with detection covariates
  (coug_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, coug_w1820_UMF))
  #'  Terrain model
  (coug_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, coug_w1820_UMF))
  (coug_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, coug_w1820_UMF))
  #'  Vegetation model
  (coug_w1820_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, coug_w1820_UMF))  #Aspect makes NDVI very significant
  #'  Habitat model
  (coug_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, coug_w1820_UMF))
  (coug_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, coug_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(coug_w1820_null, coug_w1820_null2, coug_w1820_terrain, coug_w1820_terrain2, coug_w1820_veg, coug_w1820_terrain_veg, coug_w1820_terrain_veg2)
  modSel(mods)

  
  ####  COYOTE MODELS  ####
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (coy_s1819_null <- occu(~1 ~1, coy_s1819_UMF))
  backTransform(coy_s1819_null, 'det')
  backTransform(coy_s1819_null, 'state')
  #'  Null with detection covariates
  (coy_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, coy_s1819_UMF))
  #'  Terrain model
  (coy_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, coy_s1819_UMF))
  (coy_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, coy_s1819_UMF))
  #'  Vegetation model
  (coy_s1819_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, coy_s1819_UMF))
  #'  Habitat model
  (coy_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, coy_s1819_UMF))
  (coy_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, coy_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(coy_s1819_null, coy_s1819_null2, coy_s1819_terrain, coy_s1819_terrain2, coy_s1819_veg, coy_s1819_terrain_veg, coy_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Null
  (coy_w1820_null <- occu(~1 ~1, coy_w1820_UMF))
  backTransform(coy_w1820_null, 'det')
  backTransform(coy_w1820_null, 'state')
  #'  Null with detection covariates
  (coy_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, coy_w1820_UMF))
  #'  Terrain model
  (coy_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, coy_w1820_UMF))
  (coy_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, coy_w1820_UMF))
  #'  Vegetation model
  (coy_w1820_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, coy_w1820_UMF))
  #'  Habitat model
  (coy_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, coy_w1820_UMF))
  (coy_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, coy_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(coy_w1820_null, coy_w1820_null2, coy_w1820_terrain, coy_w1820_terrain2, coy_w1820_veg, coy_w1820_terrain_veg, coy_w1820_terrain_veg2)
  modSel(mods)
  
  
  
  ####  WOLF MODELS  ####
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (wolf_s1819_null <- occu(~1 ~1, wolf_s1819_UMF))
  backTransform(wolf_s1819_null, 'det')
  backTransform(wolf_s1819_null, 'state')
  #'  Null with detection covariates
  (wolf_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, wolf_s1819_UMF))
  #'  Terrain model
  (wolf_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, wolf_s1819_UMF))
  (wolf_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, wolf_s1819_UMF))
  #'  Vegetation model
  (wolf_s1819_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, wolf_s1819_UMF))
  #'  Habitat model
  (wolf_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, wolf_s1819_UMF))
  (wolf_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, wolf_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(wolf_s1819_null, wolf_s1819_null2, wolf_s1819_terrain, wolf_s1819_terrain2, wolf_s1819_veg, wolf_s1819_terrain_veg, wolf_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Null
  (wolf_w1820_null <- occu(~1 ~1, wolf_w1820_UMF))
  backTransform(wolf_w1820_null, 'det')
  backTransform(wolf_w1820_null, 'state')
  #'  Null with detection covariates
  (wolf_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, wolf_w1820_UMF))
  #'  Terrain model
  (wolf_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, wolf_w1820_UMF))
  (wolf_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, wolf_w1820_UMF))
  #'  Vegetation model
  (wolf_w1820_veg <- occu(~Height*Distance + Trail + Year ~NDVI + NLCD, wolf_w1820_UMF))
  #'  Habitat model
  (wolf_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + NDVI + NLCD, wolf_w1820_UMF))
  (wolf_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + NDVI + NLCD, wolf_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(wolf_w1820_null, wolf_w1820_null2, wolf_w1820_terrain, wolf_w1820_terrain2, wolf_w1820_veg, wolf_w1820_terrain_veg, wolf_w1820_terrain_veg2)
  modSel(mods)
  
  
  ####  ELK MODELS ####                          #  DROP OK STUDY AREA FROM THIS SET????
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (elk_s1819_null <- occu(~1 ~1, elk_s1819_UMF))
  backTransform(elk_s1819_null, 'det')
  backTransform(elk_s1819_null, 'state')
  #'  Null with detection covariates
  (elk_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, elk_s1819_UMF))
  #'  Terrain model
  (elk_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, elk_s1819_UMF))
  (elk_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, elk_s1819_UMF))
  #'  Vegetation model
  (elk_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + NLCD, elk_s1819_UMF))
  #'  Habitat model
  (elk_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + NLCD, elk_s1819_UMF))
  (elk_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + NLCD, elk_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(elk_s1819_null, elk_s1819_null2, elk_s1819_terrain, elk_s1819_terrain2, elk_s1819_veg, elk_s1819_terrain_veg, elk_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020                    #NLCD fails to converge but Landcov ok
  #'  Null
  (elk_w1820_null <- occu(~1 ~1, elk_w1820_UMF))
  backTransform(elk_w1820_null, 'det')
  backTransform(elk_w1820_null, 'state')
  #'  Null with detection covariates
  (elk_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, elk_w1820_UMF))
  #'  Terrain model
  (elk_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, elk_w1820_UMF))
  (elk_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, elk_w1820_UMF))
  #'  Vegetation model
  (elk_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + Landcov, elk_w1820_UMF))
  #'  Habitat model
  (elk_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + Landcov, elk_w1820_UMF))
  (elk_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + Landcov, elk_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(elk_w1820_null, elk_w1820_null2, elk_w1820_terrain, elk_w1820_terrain2, elk_w1820_veg, elk_w1820_terrain_veg, elk_w1820_terrain_veg2)
  modSel(mods)
  
  
  ####  MULE DEER MODELS  ####
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (md_s1819_null <- occu(~1 ~1, md_s1819_UMF))
  backTransform(md_s1819_null, 'det')
  backTransform(md_s1819_null, 'state')
  #'  Null with detection covariates
  (md_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, md_s1819_UMF))
  #'  Terrain model
  (md_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, md_s1819_UMF))
  (md_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, md_s1819_UMF))
  #'  Vegetation model
  (md_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + NLCD, md_s1819_UMF))
  #'  Habitat model
  (md_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + NLCD, md_s1819_UMF))
  (md_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + NLCD, md_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(md_s1819_null, md_s1819_null2, md_s1819_terrain, md_s1819_terrain2, md_s1819_veg, md_s1819_terrain_veg, md_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020                    #Landcov signif but NLCD is not
  #'  Null
  (md_w1820_null <- occu(~1 ~1, md_w1820_UMF))
  backTransform(md_w1820_null, 'det')
  backTransform(md_w1820_null, 'state')
  #'  Null with detection covariates
  (md_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, md_w1820_UMF))
  #'  Terrain model
  (md_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, md_w1820_UMF))
  (md_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, md_w1820_UMF))
  #'  Vegetation model
  (md_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + Landcov, md_w1820_UMF))
  #'  Habitat model
  (md_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + Landcov, md_w1820_UMF))
  (md_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + Landcov, md_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(md_w1820_null, md_w1820_null2, md_w1820_terrain, md_w1820_terrain2, md_w1820_veg, md_w1820_terrain_veg, md_w1820_terrain_veg2)
  modSel(mods)
  
  
  ####  WHITE-TAILED DEER MODELS  ####
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (wtd_s1819_null <- occu(~1 ~1, wtd_s1819_UMF))
  backTransform(wtd_s1819_null, 'det')
  backTransform(wtd_s1819_null, 'state')
  #'  Null with detection covariates
  (wtd_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, wtd_s1819_UMF))
  #'  Terrain model
  (wtd_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, wtd_s1819_UMF))
  (wtd_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, wtd_s1819_UMF))
  #'  Vegetation model
  (wtd_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + NLCD, wtd_s1819_UMF))
  #'  Habitat model
  (wtd_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + NLCD, wtd_s1819_UMF))
  (wtd_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + NLCD, wtd_s1819_UMF))
  #'  Anthropogenic model- dist to road/road density, human pop, 
  #'  Combined model
  
  mods <- fitList(wtd_s1819_null, wtd_s1819_null2, wtd_s1819_terrain, wtd_s1819_terrain2, wtd_s1819_veg, wtd_s1819_terrain_veg, wtd_s1819_terrain_veg2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020                    #Landcov signif but NLCD is not
  #'  Null
  (wtd_w1820_null <- occu(~1 ~1, wtd_w1820_UMF))
  backTransform(wtd_w1820_null, 'det')
  backTransform(wtd_w1820_null, 'state')
  #'  Null with detection covariates
  (wtd_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, wtd_w1820_UMF))
  #'  Terrain model
  (wtd_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, wtd_w1820_UMF))
  (wtd_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, wtd_w1820_UMF))
  #'  Vegetation model
  (wtd_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + NDVI + Landcov, wtd_w1820_UMF))
  #'  Habitat model
  (wtd_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + NDVI + Landcov, wtd_w1820_UMF))
  (wtd_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + NDVI + Landcov, wtd_w1820_UMF))
  #'  Anthropogenic model
  #'  Combined model
  
  mods <- fitList(wtd_w1820_null, wtd_w1820_null2, wtd_w1820_terrain, wtd_w1820_terrain2, wtd_w1820_veg, wtd_w1820_terrain_veg, wtd_w1820_terrain_veg2)
  modSel(mods)
  
  
  
  
  