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
    mutate(
      NLCD_landcov = ifelse(NLCD_landcov == "22", "21", as.character(NLCD_landcov)), # Developed open/low intensity human use
      NLCD_landcov = ifelse(NLCD_landcov == "95", "42", as.character(NLCD_landcov)), # Evergreen forest/Emergent herbaceous wetland
      NLCD_landcov = ifelse(NLCD_landcov == "81", "71", as.character(NLCD_landcov)), # Grassland / Agriculture
      NLCD_landcov = ifelse(NLCD_landcov == "82", "71", as.character(NLCD_landcov)),
    #'  121: barren, 201: emergent wetland, 211: mesic grass, 212: xeric grass,
    #'  221: mesic shrub, 222: xeric shrub, 230: forest, 310: agriculture, 332: residential
      landcov18 = ifelse(landcov18 == "212", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "310", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "332", "211", as.character(landcov18)),  # Grassland / Agriculture / Residential
      landcov19 = ifelse(landcov19 == "121", "211", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "201", "211", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "212", "211", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "310", "211", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "332", "211", as.character(landcov19)),  # Grassland / Agriculture / Residential / Other
      landcov = ifelse(landcov == "121", "211", as.character(landcov)), 
      landcov = ifelse(landcov == "201", "211", as.character(landcov)), 
      landcov = ifelse(landcov == "212", "211", as.character(landcov)),
      landcov = ifelse(landcov == "310", "211", as.character(landcov)), 
      landcov = ifelse(landcov == "332", "211", as.character(landcov))  # Grassland / Agriculture / Residential / Other
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
  bob_s18_UMF <- unmarkedFrameOccu(DH_bob_smr18,
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
  nrow(bob_s18_UMF@y)
  bob_s18_UMF <- bob_s18_UMF[-missing_dat]
  nrow(bob_s18_UMF@y)
  ##### DOES IT MAKE SENSE TO DROP ENTIRE CAMERA SITE DUE TO 1 MISSING SITE COVARIATE THAT OFTEN ISN'T EVEN SIGNIFICANT?

  bob_w1819_UMF <- unmarkedFrameOccu(DH_bob_wtr1819,
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
  nrow(bob_w1819_UMF@y)
  bob_w1819_UMF <- bob_w1819_UMF[-missing_dat]
  nrow(bob_w1819_UMF@y)
  
  summary(bob_w1819_UMF)
  
  ####  COUGAR UMF  ####
  coug_s18_UMF <- unmarkedFrameOccu(DH_coug_smr18,
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
  nrow(coug_s18_UMF@y)
  coug_s18_UMF <- coug_s18_UMF[-missing_dat]
  nrow(coug_s18_UMF@y)
  coug_w1819_UMF <- unmarkedFrameOccu(DH_coug_wtr1819,
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
  nrow(coug_w1819_UMF@y)
  coug_w1819_UMF <- coug_w1819_UMF[-missing_dat]
  nrow(coug_w1819_UMF@y)
  coug_s19_UMF <- unmarkedFrameOccu(DH_coug_smr19,
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
  nrow(coug_s19_UMF@y)
  coug_s19_UMF <- coug_s19_UMF[-missing_dat]
  nrow(coug_s19_UMF@y)
  coug_w1920_UMF <- unmarkedFrameOccu(DH_coug_wtr1920,
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
  nrow(coug_w1920_UMF@y)
  coug_w1920_UMF <- coug_w1920_UMF[-missing_dat]
  nrow(coug_w1920_UMF@y)
  
  ####  COYOTE UMF  ####
  coy_s18_UMF <- unmarkedFrameOccu(DH_coy_smr18,
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
  coy_w1819_UMF <- unmarkedFrameOccu(DH_coy_wtr1819,
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
  
  ####  WOLF UMF  ####
  wolf_s18_UMF <- unmarkedFrameOccu(DH_wolf_smr18,
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
  wolf_w1819_UMF <- unmarkedFrameOccu(DH_wolf_wtr1819,
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
  
  ####  ELK UMF  ####
  #'  Consider removing observations/sites from the OK since so few detections and no collars there
  elk_s18_UMF <- unmarkedFrameOccu(DH_elk_smr18,
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
  elk_w1819_UMF <- unmarkedFrameOccu(DH_elk_wtr1819,
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
  
  ####  MULE DEER UMF  ####
  md_s18_UMF <- unmarkedFrameOccu(DH_md_smr18,
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
  md_w1819_UMF <- unmarkedFrameOccu(DH_md_wtr1819,
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
  
  ####  WHITE-TAILED DEER UMF  ####
  wtd_s18_UMF <- unmarkedFrameOccu(DH_wtd_smr18,
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
  wtd_w1819_UMF <- unmarkedFrameOccu(DH_wtd_wtr1819,
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
  
  
  #'  Occupancy models
  #'  =============================
  #'  unmarked formula: ~detection ~occupancy
  #'  ~1 for intercept only
  #'  
  #'  Use chi-sq test to evaluate model fit after model selection (pg. 4 vignette)
  #'  
  #'  Think about: diff sites discarded b/c missing data => different likelihoods
  #'  for each season and species... are these not comparable across seasons for 
  #'  a given species or across species for a given season?

  ####  BOBCAT MODELS  ####
  #'  Null model
  (bob_s18_null <- occu(~1 ~1, bob_s18_UMF))
  backTransform(bob_s18_null, 'det'); plogis(bob_s18_null@estimates@estimates$det@estimates)
  backTransform(bob_s18_null, 'state'); plogis(bob_s18_null@estimates@estimates$state@estimates)

  (bob_w1819_null <- occu(~1 ~1, bob_w1819_UMF))
  backTransform(bob_w1819_null, 'det'); plogis(bob_w1819_null@estimates@estimates$det@estimates)
  backTransform(bob_w1819_null, 'state'); plogis(bob_w1819_null@estimates@estimates$state@estimates)

  
  #'  Detection models
  (bob_s18_hgt <- occu(~Height ~1, bob_s18_UMF))
  (bob_s18_dist <- occu(~Distance ~1, bob_s18_UMF))
  (bob_s18_angle <- occu(~Height*Distance ~1, bob_s18_UMF))
  (bob_s18_trail <- occu(~Trail ~1, bob_s18_UMF))
  (bob_s18_angle.trail <- occu(~Height*Distance + Trail ~1, bob_s18_UMF))
  bob_s18_detlist <- fitList(bob_s18_null, bob_s18_hgt, bob_s18_dist, bob_s18_angle, bob_s18_trail, bob_s18_angle.trail)
  modSel(bob_s18_detlist)

  (bob_w1819_hgt <- occu(~Height ~1, bob_w1819_UMF))
  (bob_w1819_dist <- occu(~Distance ~1, bob_w1819_UMF))
  (bob_w1819_angle <- occu(~Height*Distance ~1, bob_w1819_UMF))
  (bob_w1819_trail <- occu(~Trail ~1, bob_w1819_UMF))
  (bob_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, bob_w1819_UMF))
  bob_w1819_detlist <- fitList(bob_w1819_null, bob_w1819_hgt, bob_w1819_dist, bob_w1819_angle, bob_w1819_trail, bob_w1819_angle.trail) 
  modSel(bob_w1819_detlist)

  #'  Quick notes 2/17: game trails negatively affect detection in summer, no affect in winter
  #'  camera height & distance don't seem to matter
  
  #'  Occupancy & detection models
  (bob_s18_SA <- occu(~Height*Distance + Trail ~Area, bob_s18_UMF))
  (bob_s18_H2o <- occu(~Height*Distance + Trail ~NearestH2o, bob_s18_UMF))
  (bob_s18_Landcov <- occu(~Height*Distance + Trail ~Landcov, bob_s18_UMF))
  (bob_s18_Habitat <- occu(~Height*Distance + Trail ~Habitat, bob_s18_UMF)) # shrub-steppe doesn't converge
  (bob_s18_NLCD <- occu(~Height*Distance + Trail ~NLCD, bob_s18_UMF))       # Fails to converge
  (bob_s18_NDVI <- occu(~Height*Distance + Trail ~NDVI, bob_s18_UMF))
  (bob_s18_Elev <- occu(~Height*Distance + Trail ~Elev, bob_s18_UMF))
  (bob_s18_Elev2 <- occu(~Height*Distance + Trail ~Elev + I(Elev^2), bob_s18_UMF))
  (bob_s18_Slope <- occu(~Height*Distance + Trail ~Slope, bob_s18_UMF))
  (bob_s18_Aspect <- occu(~Height*Distance + Trail ~Aspect, bob_s18_UMF))
  (bob_s18_Aspect2 <- occu(~Height*Distance + Trail ~Aspect + I(Aspect^2), bob_s18_UMF))
  (bob_s18_TRI <- occu(~Height*Distance + Trail ~TRI, bob_s18_UMF))
  (bob_s18_CanopyCov <- occu(~Height*Distance + Trail ~Canopy_cov, bob_s18_UMF))  # MISSING 1 OBSERVATION SO SCREWS THINGS UP- DROP THAT ROW???
  (bob_s18_TreeCov <- occu(~Height*Distance + Trail ~Tree_cov, bob_s18_UMF))
  
  (bob_s18_SA_NLCD <- occu(~Height*Distance + Trail ~Area + NLCD, bob_s18_UMF))        # Fails to converge
  (bob_s18_SA_Habitat <- occu(~Height*Distance + Trail ~Area + Habitat, bob_s18_UMF))  # shrub-steppe doesn't converge
  (bob_s18_SA_Landcov <- occu(~Height*Distance + Trail ~Area + Landcov, bob_s18_UMF))
  (bob_s18_SA_NearestH2o <- occu(~Height*Distance + Trail ~Area + NearestH2o, bob_s18_UMF))  # Fails to converge
  (bob_s18_SA_Elev2 <- occu(~Height*Distance + Trail ~Area + Elev + I(Elev^2), bob_s18_UMF))
  (bob_s18_SA_Aspect <- occu(~Height*Distance + Trail ~Area + Aspect, bob_s18_UMF))
  (bob_s18_SA_TRI <- occu(~Height*Distance + Trail ~Area + TRI, bob_s18_UMF))
  (bob_s18_SA_Slope <- occu(~Height*Distance + Trail ~Area + Slope, bob_s18_UMF))
  (bob_s18_SA_Treecov <- occu(~Height*Distance + Trail ~Area + Tree_cov, bob_s18_UMF))
  (bob_s18_SA_NDVI <- occu(~Height*Distance + Trail ~Area + NDVI, bob_s18_UMF))
  (bob_s18_SA_NLCD_Elev2 <- occu(~Height*Distance + Trail ~Area + NLCD + Elev + I(Elev^2), bob_s18_UMF))

  #'  Quick notes 2/17: 
  #'  Aspect + effect, NDVI almost signif + effect
  #'  OK study area - effect when some additional covs included
  #'  NLCD - effect (and coverges) when OK and Elev2 added
  
  #'  Model selection (exclude models where covaraites clearly did not converge)
  bob_s18_list <- fitList(bob_s18_SA, bob_s18_H2o, bob_s18_Landcov, bob_s18_NDVI, 
                          bob_s18_Elev, bob_s18_Elev2, bob_s18_Slope, bob_s18_Aspect,
                          bob_s18_TRI, bob_s18_TreeCov,
                          bob_s18_SA_Landcov, bob_s18_SA_Elev2, bob_s18_SA_Aspect, 
                          bob_s18_SA_TRI, bob_s18_SA_Slope, bob_s18_SA_Treecov, 
                          bob_s18_SA_NDVI, bob_s18_SA_NLCD_Elev2)  #bob_s18_CanopyCov
  modSel(bob_s18_list)
  
  
  (bob_w1819_SA <- occu(~1 ~Area, bob_w1819_UMF))
  (bob_w1819_H2o <- occu(~1 ~NearestH2o, bob_w1819_UMF))
  (bob_w1819_Landcov <- occu(~1 ~Landcov, bob_w1819_UMF))
  (bob_w1819_Habitat <- occu(~1 ~Habitat, bob_w1819_UMF))  
  (bob_w1819_NLCD <- occu(~1 ~NLCD, bob_w1819_UMF))
  (bob_w1819_NDVI <- occu(~1 ~NDVI, bob_w1819_UMF))
  (bob_w1819_dNBR <- occu(~1 ~dNBR, bob_w1819_UMF))
  (bob_w1819_Elev <- occu(~1 ~Elev, bob_w1819_UMF))
  (bob_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), bob_w1819_UMF))
  (bob_w1819_Slope <- occu(~1 ~Slope, bob_w1819_UMF))
  (bob_w1819_Aspect <- occu(~1 ~Aspect, bob_w1819_UMF))
  (bob_w1819_Aspect2 <- occu(~1 ~Aspect + I(Aspect^2), bob_w1819_UMF))
  (bob_w1819_TRI <- occu(~1 ~TRI, bob_w1819_UMF))
  (bob_w1819_CanopyCov <- occu(~1 ~Canopy_cov, bob_w1819_UMF))
  (bob_w1819_TreeCov <- occu(~1 ~Tree_cov, bob_w1819_UMF))
  
  
  (bob_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, bob_w1819_UMF))
  (bob_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, bob_w1819_UMF))  
  (bob_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, bob_w1819_UMF))
  (bob_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, bob_w1819_UMF))
  (bob_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), bob_w1819_UMF))
  (bob_w1819_SA_Aspect <- occu(~1 ~Area + Aspect, bob_w1819_UMF))
  (bob_w1819_SA_TRI <- occu(~1 ~Area + TRI, bob_w1819_UMF))
  (bob_w1819_SA_Slope <- occu(~1 ~Area + Slope, bob_w1819_UMF))
  (bob_w1819_SA_Tree_cov <- occu(~1 ~Area + Tree_cov, bob_w1819_UMF))
  (bob_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, bob_w1819_UMF))
  (bob_w1819_SA_NLCD_Elev2 <- occu(~1 ~Area + NLCD + Elev + I(Elev^2), bob_w1819_UMF))

  #'  Model selection (exclude models where covaraites clearly did not converge)
  bob_w1819_list <- fitList(bob_w1819_SA, bob_w1819_H2o, bob_w1819_Landcov, 
                          bob_w1819_Habitat, bob_w1819_NDVI, bob_w1819_NLCD, 
                          bob_w1819_dNBR, bob_w1819_Elev, bob_w1819_Elev2, 
                          bob_w1819_Slope, bob_w1819_Aspect, bob_w1819_TRI, 
                          bob_w1819_TreeCov, bob_w1819_SA_NLCD, bob_w1819_SA_Habitat, 
                          bob_w1819_SA_Landcov, bob_w1819_SA_NearestH2o, 
                          bob_w1819_SA_Elev2, bob_w1819_SA_Aspect, 
                          bob_w1819_SA_TRI, bob_w1819_SA_Slope, bob_w1819_SA_NDVI, 
                          bob_w1819_SA_NLCD_Elev2)  #bob_w1819_CanopyCov
  modSel(bob_w1819_list)
  
  ####  COUGAR MODELS  ####
  #'  Null model
  (coug_s18_null <- occu(~1 ~1, coug_s18_UMF))
  backTransform(coug_s18_null, 'det')
  backTransform(coug_s18_null, 'state')
  
  (coug_w1819_null <- occu(~1 ~1, coug_w1819_UMF))
  backTransform(coug_w1819_null, 'det')
  backTransform(coug_w1819_null, 'state')
  
  #'  Detection models
  # (coug_s18_hgt <- occu(~Height ~1, coug_s18_UMF))
  # (coug_s18_dist <- occu(~Distance ~1, coug_s18_UMF))
  # (coug_s18_angle <- occu(~Height*Distance ~1, coug_s18_UMF))
  # (coug_s18_trail <- occu(~Trail ~1, coug_s18_UMF))
  # (coug_s18_angle.trail <- occu(~Height*Distance + Trail ~1, coug_s18_UMF))
  # 
  # (coug_w1819_hgt <- occu(~Height ~1, coug_w1819_UMF))    #  SIGNIF at 0.05
  # (coug_w1819_dist <- occu(~Distance ~1, coug_w1819_UMF)) #  SIGNIF at 0.05
  # (coug_w1819_angle <- occu(~Height*Distance ~1, coug_w1819_UMF)) 
  # (coug_w1819_hgt_dist <- occu(~Height + Distance ~1, coug_w1819_UMF))
  # (coug_w1819_trail <- occu(~Trail ~1, coug_w1819_UMF))
  # (coug_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, coug_w1819_UMF))
  # mods <- fitList(coug_w1819_hgt, coug_w1819_dist, coug_w1819_hgt_dist, coug_w1819_angle) 
  # modsel <- modSel(mods)
  
  #'  Occupancy & detection models
  (coug_s18_SA <- occu(~Height*Distance + Trail ~Area, coug_s18_UMF))
  (coug_s18_H2o <- occu(~Height*Distance + Trail ~NearestH2o, coug_s18_UMF))
  (coug_s18_Landcov <- occu(~Height*Distance + Trail ~Landcov, coug_s18_UMF))  # only partially converges
  (coug_s18_Habitat <- occu(~Height*Distance + Trail ~Habitat, coug_s18_UMF))  # Fails to converge
  (coug_s18_NLCD <- occu(~Height*Distance + Trail ~NLCD, coug_s18_UMF))
  (coug_s18_NDVI <- occu(~Height*Distance + Trail ~NDVI, coug_s18_UMF))
  (coug_s18_Elev <- occu(~Height*Distance + Trail ~Elev, coug_s18_UMF))
  (coug_s18_Elev2 <- occu(~Height*Distance + Trail ~Elev + I(Elev^2), coug_s18_UMF))
  (coug_s18_Slope <- occu(~Height*Distance + Trail ~Slope, coug_s18_UMF))
  (coug_s18_Aspect <- occu(~Height*Distance + Trail ~Aspect, coug_s18_UMF))
  (coug_s18_TRI <- occu(~Height*Distance + Trail ~TRI, coug_s18_UMF))
  (coug_s18_CanopyCov <- occu(~Height*Distance + Trail ~Canopy_cov, coug_s18_UMF))
  (coug_s18_TreeCov <- occu(~Height*Distance + Trail ~Tree_cov, coug_s18_UMF))
  
  (coug_s18_SA_NDVI <- occu(~Trail ~Area + NDVI, coug_s18_UMF))
  (coug_s18_SA_NLCD <- occu(~Trail ~Area + NLCD, coug_s18_UMF))        # partially converges
  (coug_s18_SA_Habitat <- occu(~Trail ~Area + Habitat, coug_s18_UMF))  # fails to converges
  (coug_s18_SA_Landcov <- occu(~Trail ~Area + Landcov, coug_s18_UMF))  # partially converges
  (coug_s18_SA_TRI <- occu(~Trail ~Area + TRI, coug_s18_UMF)) 
  (coug_s18_SA_Slope <- occu(~Trail ~Area + Slope, coug_s18_UMF))
  (coug_s18_SA_Elev <- occu(~Trail ~Area + Elev, coug_s18_UMF))
  (coug_s18_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), coug_s18_UMF))
  (coug_s18_SA_Canopycov <- occu(~Trail ~Area + Canopy_cov, coug_s18_UMF))
  (coug_s18_SA_Treecov <- occu(~Trail ~Area + Tree_cov, coug_s18_UMF))
  (coug_s18_SA_Elev_NDVI <- occu(~Trail ~Area +Elev + NDVI, coug_s18_UMF))
  (coug_s18_Elev_NDVI <- occu(~Trail ~Elev + NDVI, coug_s18_UMF))
  (coug_s18_Elev_NLCD <- occu(~Trail ~Elev + NLCD, coug_s18_UMF))        
  (coug_s18_Elev_Habitat <- occu(~Trail ~Elev + Habitat, coug_s18_UMF))  # partially converges
  (coug_s18_Elev_TRI <- occu(~Trail ~Elev + TRI, coug_s18_UMF)) 
  (coug_s18_Elev_Slope <- occu(~Trail ~Elev + Slope, coug_s18_UMF))
  (coug_s18_Elev_Canopycov <- occu(~Trail ~Elev + Canopy_cov, coug_s18_UMF))
  (coug_s18_Elev_Treecov <- occu(~Trail ~Elev + Tree_cov, coug_s18_UMF))
  (coug_s18_Elev_Slope_Canopycov <- occu(~Trail ~Elev + Slope + Canopy_cov, coug_s18_UMF))
  (coug_s18_Elev_TRI_Canopycov <- occu(~Trail ~Elev + TRI + Canopy_cov, coug_s18_UMF))
  
  #'  NOTES 2/18
  #'  OK study area - effect
  #'  NDVI, tree cov + effect
  #'  Landcov & NLCD partially signif, Elev2
  #'  OK study area - effect with terrain variables, all signif
  
  (coug_w1819_SA <- occu(~Distance ~Area, coug_w1819_UMF))
  (coug_w1819_H2o <- occu(~Distance ~NearestH2o, coug_w1819_UMF))
  (coug_w1819_Landcov <- occu(~Distance ~Landcov, coug_w1819_UMF))  
  (coug_w1819_Habitat <- occu(~Distance ~Habitat, coug_w1819_UMF))  
  (coug_w1819_NLCD <- occu(~Distance ~NLCD, coug_w1819_UMF))
  (coug_w1819_NDVI <- occu(~Distance ~NDVI, coug_w1819_UMF))
  (coug_w1819_dNBR <- occu(~Distance ~dNBR, coug_w1819_UMF))
  (coug_w1819_Elev <- occu(~Distance ~Elev, coug_w1819_UMF))
  (coug_w1819_Elev2 <- occu(~Distance ~Elev + I(Elev^2), coug_w1819_UMF))
  (coug_w1819_Slope <- occu(~Distance ~Slope, coug_w1819_UMF))
  (coug_w1819_Aspect <- occu(~Distance ~Aspect, coug_w1819_UMF))
  (coug_w1819_TRI <- occu(~Distance ~TRI, coug_w1819_UMF))
  (coug_w1819_CanopyCov <- occu(~Distance ~Canopy_cov, coug_w1819_UMF))
  (coug_w1819_TreeCov <- occu(~Distance ~Tree_cov, coug_w1819_UMF))
  
  (coug_w1819_SA_NDVI <- occu(~Distance ~Area + NDVI, coug_w1819_UMF))
  (coug_w1819_SA_NLCD <- occu(~Distance ~Area + NLCD, coug_w1819_UMF))        # partially converges
  (coug_w1819_SA_Habitat <- occu(~Distance ~Area + Habitat, coug_w1819_UMF))  
  (coug_w1819_SA_Landcov <- occu(~Distance ~Area + Landcov, coug_w1819_UMF))  
  (coug_w1819_SA_TRI <- occu(~Distance ~Area + TRI, coug_w1819_UMF)) 
  (coug_w1819_SA_Slope <- occu(~Distance ~Area + Slope, coug_w1819_UMF))
  (coug_w1819_SA_Elev <- occu(~Distance ~Area + Elev, coug_w1819_UMF))
  (coug_w1819_SA_Elev2 <- occu(~Distance ~Area + Elev + I(Elev^2), coug_w1819_UMF))
  (coug_w1819_SA_Canopycov <- occu(~Distance ~Area + Canopy_cov, coug_w1819_UMF))
  (coug_w1819_SA_Treecov <- occu(~Distance ~Area + Tree_cov, coug_w1819_UMF))
  (coug_w1819_SA_Elev_NDVI <- occu(~Distance ~Area +Elev + NDVI, coug_w1819_UMF))
  (coug_w1819_Elev_NDVI <- occu(~Distance ~Elev + NDVI, coug_w1819_UMF))
  (coug_w1819_Elev_NLCD <- occu(~Distance ~Elev + NLCD, coug_w1819_UMF))        
  (coug_w1819_Elev_Habitat <- occu(~Distance ~Elev + Habitat, coug_w1819_UMF)) 
  (coug_w1819_Elev_Landcov <- occu(~Distance ~Elev + Landcov, coug_w1819_UMF))
  (coug_w1819_Elev_TRI <- occu(~Distance ~Elev + TRI, coug_w1819_UMF)) 
  (coug_w1819_Elev_Slope <- occu(~Distance ~Elev + Slope, coug_w1819_UMF))
  (coug_w1819_Elev_Canopycov <- occu(~Distance ~Elev + Canopy_cov, coug_w1819_UMF))
  (coug_w1819_Elev_Treecov <- occu(~Distance ~Elev + Tree_cov, coug_w1819_UMF))
  (coug_w1819_Elev_Slope_Canopycov <- occu(~Distance ~Elev + Slope + Canopy_cov, coug_w1819_UMF))
  (coug_w1819_Elev_TRI_Canopycov <- occu(~Distance ~Elev + TRI + Canopy_cov, coug_w1819_UMF))
  
  #'  Using Distance on detection in winter but Height could also be used- need model selection to make final call
  
  #'  NOTES 2/18:
  #'  NDVI, canopy cov + effect
  #'  Slope + effect with study area
  
  ####  COYOTE MODELS  ####
  #'  Null model
  (coy_s18_null <- occu(~1 ~1, coy_s18_UMF))
  backTransform(coy_s18_null, 'det')
  backTransform(coy_s18_null, 'state')
  
  (coy_w1819_null <- occu(~1 ~1, coy_w1819_UMF))
  backTransform(coy_w1819_null, 'det')
  backTransform(coy_w1819_null, 'state')
  
  #'  Detection models
  (coy_s18_hgt <- occu(~Height ~1, coy_s18_UMF))
  (coy_s18_dist <- occu(~Distance ~1, coy_s18_UMF))
  (coy_s18_angle <- occu(~Height*Distance ~1, coy_s18_UMF))
  (coy_s18_trail <- occu(~Trail ~1, coy_s18_UMF))
  (coy_s18_angle.trail <- occu(~Height*Distance + Trail ~1, coy_s18_UMF))
  
  (coy_w1819_hgt <- occu(~Height ~1, coy_w1819_UMF))
  (coy_w1819_dist <- occu(~Distance ~1, coy_w1819_UMF))
  (coy_w1819_angle <- occu(~Height*Distance ~1, coy_w1819_UMF))
  (coy_w1819_trail <- occu(~Trail ~1, coy_w1819_UMF))
  (coy_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, coy_w1819_UMF))
  
  #'  Trail - effect summer
  
  #'  Occupancy & detection models
  (coy_s18_SA <- occu(~Trail ~Area, coy_s18_UMF))
  (coy_s18_H2o <- occu(~Trail ~NearestH2o, coy_s18_UMF))
  (coy_s18_Landcov <- occu(~Trail ~Landcov, coy_s18_UMF))  # only partially converges
  (coy_s18_Habitat <- occu(~Trail ~Habitat, coy_s18_UMF))  # only partially converges
  (coy_s18_NLCD <- occu(~Trail ~NLCD, coy_s18_UMF))
  (coy_s18_NDVI <- occu(~Trail ~NDVI, coy_s18_UMF))
  (coy_s18_Elev <- occu(~Trail ~Elev, coy_s18_UMF))
  (coy_s18_Elev2 <- occu(~Trail ~Elev + I(Elev^2), coy_s18_UMF))
  (coy_s18_Slope <- occu(~Trail ~Slope, coy_s18_UMF))
  (coy_s18_Aspect <- occu(~Trail ~Aspect, coy_s18_UMF))
  (coy_s18_TRI <- occu(~Trail ~TRI, coy_s18_UMF))
  (coy_s18_CanopyCov <- occu(~Trail ~Canopy_cov, coy_s18_UMF))
  (coy_s18_TreeCov <- occu(~Trail ~Tree_cov, coy_s18_UMF))
  
  (coy_s18_SA_NDVI <- occu(~Trail ~Area + NDVI, coy_s18_UMF))
  (coy_s18_SA_NLCD <- occu(~Trail ~Area + NLCD, coy_s18_UMF))        
  (coy_s18_SA_Habitat <- occu(~Trail ~Area + Habitat, coy_s18_UMF))  # partially converges
  (coy_s18_SA_Landcov <- occu(~Trail ~Area + Landcov, coy_s18_UMF))  # partially converges
  (coy_s18_SA_TRI <- occu(~Trail ~Area + TRI, coy_s18_UMF)) 
  (coy_s18_SA_Slope <- occu(~Trail ~Area + Slope, coy_s18_UMF))
  (coy_s18_SA_Elev <- occu(~Trail ~Area + Elev, coy_s18_UMF))
  (coy_s18_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), coy_s18_UMF))
  (coy_s18_SA_Canopycov <- occu(~Trail ~Area + Canopy_cov, coy_s18_UMF))
  (coy_s18_SA_Treecov <- occu(~Trail ~Area + Tree_cov, coy_s18_UMF))
  (coy_s18_Elev_NDVI <- occu(~Trail ~Elev + NDVI, coy_s18_UMF))
  (coy_s18_Elev_NLCD <- occu(~Trail ~Elev + NLCD, coy_s18_UMF))        
  (coy_s18_Elev_Habitat <- occu(~Trail ~Elev + Habitat, coy_s18_UMF))  # partially converges
  (coy_s18_Elev_TRI <- occu(~Trail ~Elev + TRI, coy_s18_UMF)) 
  (coy_s18_Elev_Slope <- occu(~Trail ~Elev + Slope, coy_s18_UMF))
  (coy_s18_Elev_Canopycov <- occu(~Trail ~Elev + Canopy_cov, coy_s18_UMF))
  (coy_s18_Elev_Treecov <- occu(~Trail ~Elev + Tree_cov, coy_s18_UMF))
  (coy_s18_Elev_Slope_Canopycov <- occu(~Trail ~Elev + Slope + Canopy_cov, coy_s18_UMF))
  (coy_s18_Elev_TRI_Canopycov <- occu(~Trail ~Elev + TRI + Canopy_cov, coy_s18_UMF))

  #'  NOTES 2/18:
  #'  Elev, slope, tri, canopy cov - effect 
  #'  Elev & tri, Elev & slope, Elev & canopy cov - together
  #'  Elev & Slope & Canopy - together, Elev & TRI & Canopy - together
  
  (coy_w1819_SA <- occu(~1 ~Area, coy_w1819_UMF))
  (coy_w1819_H2o <- occu(~1 ~NearestH2o, coy_w1819_UMF))
  (coy_w1819_Landcov <- occu(~1 ~Landcov, coy_w1819_UMF))  
  (coy_w1819_Habitat <- occu(~1 ~Habitat, coy_w1819_UMF))  # Shrub-steppe fails to converge
  (coy_w1819_NLCD <- occu(~1 ~NLCD, coy_w1819_UMF))        
  (coy_w1819_NDVI <- occu(~1 ~NDVI, coy_w1819_UMF))
  (coy_w1819_dNBR <- occu(~1 ~dNBR, coy_w1819_UMF))
  (coy_w1819_Elev <- occu(~1 ~Elev, coy_w1819_UMF))
  (coy_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), coy_w1819_UMF))
  (coy_w1819_Slope <- occu(~1 ~Slope, coy_w1819_UMF))
  (coy_w1819_Aspect <- occu(~1 ~Aspect, coy_w1819_UMF))
  (coy_w1819_TRI <- occu(~1 ~TRI, coy_w1819_UMF))
  (coy_w1819_CanopyCov <- occu(~1 ~Canopy_cov, coy_w1819_UMF))
  (coy_w1819_TreeCov <- occu(~1 ~Tree_cov, coy_w1819_UMF))
  
  (coy_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, coy_w1819_UMF))        
  (coy_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, coy_w1819_UMF))  # Shrub-steppe fails to converge
  (coy_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, coy_w1819_UMF))
  (coy_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, coy_w1819_UMF))
  (coy_w1819_SA_Elev <- occu(~1 ~Area + Elev, coy_w1819_UMF))
  (coy_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), coy_w1819_UMF))
  (coy_w1819_SA_Canopycov <- occu(~1 ~Area + Canopy_cov, coy_w1819_UMF))
  (coy_w1819_SA_Treecov <- occu(~1 ~Area + Tree_cov, coy_w1819_UMF))
  (coy_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, coy_w1819_UMF))
  (coy_w1819_SA_NDVI_Lancov <- occu(~1 ~Area + NDVI + Landcov, coy_w1819_UMF))
  (coy_w1819_SA_Elev_Slope <- occu(~1 ~Area + Elev + Slope, coy_w1819_UMF))
  (coy_w1819_Elev_NLCD <- occu(~1 ~Elev + NLCD, coy_w1819_UMF))        
  (coy_w1819_Elev_Habitat <- occu(~1 ~Elev + Habitat, coy_w1819_UMF))  # Shrub-steppe fails to converge
  (coy_w1819_Elev_Landcov <- occu(~1 ~Elev + Landcov, coy_w1819_UMF))
  (coy_w1819_Elev_Slope <- occu(~1 ~Elev + Slope, coy_w1819_UMF))
  (coy_w1819_Elev_TRI <- occu(~1 ~Elev + TRI, coy_w1819_UMF))
  (coy_w1819_Elev_NearestH2o <- occu(~1 ~Elev + NearestH2o, coy_w1819_UMF))
  (coy_w1819_Elev_Canopycov <- occu(~1 ~Elev + Canopy_cov, coy_w1819_UMF))
  (coy_w1819_Elev_Treecov <- occu(~1 ~Elev + Tree_cov, coy_w1819_UMF))
  (coy_w1819_Elev_NDVI <- occu(~1 ~Elev + NDVI, coy_w1819_UMF))

  #'  NOTES 2/18:
  #'  Grassland & Mixed conifer + effect (but shrub fails to converge)
  #'  Elev, Slope, TRI - effect
  #'  OK study area - effect when NDVI, NLCD or landcov added, some of landcov signif, NDVI - effect
  #'  Tree cov - effect when study area included
  
  
  ####  WOLF MODELS  ####
  #'  Null model
  (wolf_s18_null <- occu(~1 ~1, wolf_s18_UMF))
  backTransform(wolf_s18_null, 'det')
  backTransform(wolf_s18_null, 'state')
  
  (wolf_w1819_null <- occu(~1 ~1, wolf_w1819_UMF))
  backTransform(wolf_w1819_null, 'det')
  backTransform(wolf_w1819_null, 'state')
  
  #'  Detection models
  (wolf_s18_hgt <- occu(~Height ~1, wolf_s18_UMF))
  (wolf_s18_dist <- occu(~Distance ~1, wolf_s18_UMF))
  (wolf_s18_angle <- occu(~Height*Distance ~1, wolf_s18_UMF))
  (wolf_s18_trail <- occu(~Trail ~1, wolf_s18_UMF))
  (wolf_s18_angle.trail <- occu(~Height*Distance + Trail ~1, wolf_s18_UMF))
  
  (wolf_w1819_hgt <- occu(~Height ~1, wolf_w1819_UMF))
  (wolf_w1819_dist <- occu(~Distance ~1, wolf_w1819_UMF))
  (wolf_w1819_angle <- occu(~Height*Distance ~1, wolf_w1819_UMF))
  (wolf_w1819_trail <- occu(~Trail ~1, wolf_w1819_UMF))
  (wolf_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, wolf_w1819_UMF))
  
  #  NOTES 2/18: trail - effect in summer
  
  #'  Occupancy & detection models
  (wolf_s18_SA <- occu(~Trail ~Area, wolf_s18_UMF))
  (wolf_s18_H2o <- occu(~Trail ~NearestH2o, wolf_s18_UMF))
  (wolf_s18_Landcov <- occu(~Trail ~Landcov, wolf_s18_UMF))  # only partially converges
  (wolf_s18_Habitat <- occu(~Trail ~Habitat, wolf_s18_UMF))  # only partially converges
  (wolf_s18_NLCD <- occu(~Trail ~NLCD, wolf_s18_UMF))
  (wolf_s18_NDVI <- occu(~Trail ~NDVI, wolf_s18_UMF))
  (wolf_s18_Elev <- occu(~Trail ~Elev, wolf_s18_UMF))
  (wolf_s18_Elev2 <- occu(~Trail ~Elev + I(Elev^2), wolf_s18_UMF))
  (wolf_s18_Slope <- occu(~Trail ~Slope, wolf_s18_UMF))
  (wolf_s18_Aspect <- occu(~Trail ~Aspect, wolf_s18_UMF))
  (wolf_s18_TRI <- occu(~Trail ~TRI, wolf_s18_UMF))
  (wolf_s18_CanopyCov <- occu(~Trail ~Canopy_cov, wolf_s18_UMF))
  (wolf_s18_TreeCov <- occu(~Trail ~Tree_cov, wolf_s18_UMF))
  
  (wolf_s18_SA_NLCD <- occu(~Trail ~Area + NLCD, wolf_s18_UMF))        
  (wolf_s18_SA_Habitat <- occu(~Trail ~Area + Habitat, wolf_s18_UMF))  # only partially converges
  (wolf_s18_SA_Landcov <- occu(~Trail ~Area + Landcov, wolf_s18_UMF))  # only partially converges
  (wolf_s18_SA_NearestH2o <- occu(~Trail ~Area + NearestH2o, wolf_s18_UMF))
  (wolf_s18_SA_Elev <- occu(~Trail ~Area + Elev, wolf_s18_UMF))
  (wolf_s18_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), wolf_s18_UMF))
  (wolf_s18_SA_Canopycov <- occu(~Trail ~Area + Canopy_cov, wolf_s18_UMF))
  (wolf_s18_SA_Treecov <- occu(~Trail ~Area + Tree_cov, wolf_s18_UMF))
  (wolf_s18_SA_NDVI <- occu(~Trail ~Area + NDVI, wolf_s18_UMF))
  
  #' NOTES 2/18:
  #' Elev almost + effect, Elev + effect with study area included
  
  (wolf_w1819_SA <- occu(~1 ~Area, wolf_w1819_UMF))          # fails to converge
  (wolf_w1819_H2o <- occu(~1 ~NearestH2o, wolf_w1819_UMF))
  (wolf_w1819_Landcov <- occu(~1 ~Landcov, wolf_w1819_UMF))  # only partially converges
  (wolf_w1819_Habitat <- occu(~1 ~Habitat, wolf_w1819_UMF))  # only partially converges
  (wolf_w1819_NLCD <- occu(~1 ~NLCD, wolf_w1819_UMF))        # fails to converge
  (wolf_w1819_NDVI <- occu(~1 ~NDVI, wolf_w1819_UMF))        # fails to converge
  (wolf_w1819_dNBR <- occu(~1 ~dNBR, wolf_w1819_UMF))        # fails to converge
  (wolf_w1819_Elev <- occu(~1 ~Elev, wolf_w1819_UMF))
  (wolf_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), wolf_w1819_UMF))
  (wolf_w1819_Slope <- occu(~1 ~Slope, wolf_w1819_UMF))
  (wolf_w1819_Aspect <- occu(~1 ~Aspect, wolf_w1819_UMF))
  (wolf_w1819_TRI <- occu(~1 ~TRI, wolf_w1819_UMF))
  (wolf_w1819_CanopyCov <- occu(~1 ~Canopy_cov, wolf_w1819_UMF))
  (wolf_w1819_TreeCov <- occu(~1 ~Tree_cov, wolf_w1819_UMF))
  
  # (wolf_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, wolf_w1819_UMF))        
  # (wolf_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, wolf_w1819_UMF))  
  # (wolf_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, wolf_w1819_UMF))  
  # (wolf_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, wolf_w1819_UMF)) 
  # (wolf_w1819_SA_Elev <- occu(~1 ~Area + Elev, wolf_w1819_UMF))
  # (wolf_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), wolf_w1819_UMF))
  # (wolf_w1819_SA_Canopycov <- occu(~1 ~Area + Canopy_cov, wolf_w1819_UMF))
  # (wolf_w1819_SA_Treecov <- occu(~1 ~Area + Tree_cov, wolf_w1819_UMF))
  # (wolf_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, wolf_w1819_UMF))
  # (wolf_w1819_SA_NLCD_Elev2 <- occu(~1 ~Area + Habitat + Elev + I(Elev^2), wolf_w1819_UMF))
  # (wolf_w1819_SA_NLCD_Elev2 <- occu(~1 ~Area + Landcov + Elev + I(Elev^2), wolf_w1819_UMF))
  
  
  ####  ELK MODELS ####
  #'  Null model
  (elk_s18_null <- occu(~1 ~1, elk_s18_UMF))
  backTransform(elk_s18_null, 'det')
  backTransform(elk_s18_null, 'state')
  
  (elk_w1819_null <- occu(~1 ~1, elk_w1819_UMF))
  backTransform(elk_w1819_null, 'det')
  backTransform(elk_w1819_null, 'state')
  
  #'  Detection models
  (elk_s18_hgt <- occu(~Height ~1, elk_s18_UMF))
  (elk_s18_dist <- occu(~Distance ~1, elk_s18_UMF))
  (elk_s18_angle <- occu(~Height*Distance ~1, elk_s18_UMF))
  (elk_s18_trail <- occu(~Trail ~1, elk_s18_UMF))           
  (elk_s18_angle.trail <- occu(~Height*Distance + Trail ~1, elk_s18_UMF))
  
  (elk_w1819_hgt <- occu(~Height ~1, elk_w1819_UMF))
  (elk_w1819_dist <- occu(~Distance ~1, elk_w1819_UMF))
  (elk_w1819_angle <- occu(~Height*Distance ~1, elk_w1819_UMF))
  (elk_w1819_trail <- occu(~Trail ~1, elk_w1819_UMF))       
  (elk_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, elk_w1819_UMF))
  
  #'  Tail - effect summer
  
  #'  Occupancy & detection models
  (elk_s18_SA <- occu(~Trail ~Area, elk_s18_UMF))
  (elk_s18_H2o <- occu(~Trail ~NearestH2o, elk_s18_UMF))
  (elk_s18_Landcov <- occu(~Trail ~Landcov, elk_s18_UMF))
  (elk_s18_Habitat <- occu(~Trail ~Habitat, elk_s18_UMF))  # Shrub-steppe fails to converge
  (elk_s18_NLCD <- occu(~Trail ~NLCD, elk_s18_UMF))
  (elk_s18_NDVI <- occu(~Trail ~NDVI, elk_s18_UMF))
  (elk_s18_Elev <- occu(~Trail ~Elev, elk_s18_UMF))
  (elk_s18_Elev2 <- occu(~Trail ~Elev + I(Elev^2), elk_s18_UMF))
  (elk_s18_Slope <- occu(~Trail ~Slope, elk_s18_UMF))
  (elk_s18_Aspect <- occu(~Trail ~Aspect, elk_s18_UMF))
  (elk_s18_TRI <- occu(~Trail ~TRI, elk_s18_UMF))
  (elk_s18_CanopyCov <- occu(~Trail ~Canopy_cov, elk_s18_UMF))
  (elk_s18_TreeCov <- occu(~Trail ~Tree_cov, elk_s18_UMF))
  
  (elk_s18_SA_NLCD <- occu(~Trail ~Area + NLCD, elk_s18_UMF))        # NLCD FAILS TO CONVERGE
  (elk_s18_SA_Habitat <- occu(~Trail ~Area + Habitat, elk_s18_UMF))  # Shrub-steppe fails to converge
  (elk_s18_SA_Landcov <- occu(~Trail ~Area + Landcov, elk_s18_UMF))
  (elk_s18_SA_NearestH2o <- occu(~Trail ~Area + NearestH2o, elk_s18_UMF))
  (elk_s18_SA_Elev <- occu(~Trail ~Area + Elev, elk_s18_UMF))
  (elk_s18_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), elk_s18_UMF))
  (elk_s18_SA_Canopycov <- occu(~Trail ~Area + Canopy_cov, elk_s18_UMF))
  (elk_s18_SA_Treecov <- occu(~Trail ~Area + Tree_cov, elk_s18_UMF))
  (elk_s18_SA_NDVI <- occu(~Trail ~Area + NDVI, elk_s18_UMF))
  (elk_s18_SA_Habitat_Elev2 <- occu(~Trail ~Area + Habitat + Elev + I(Elev^2), elk_s18_UMF))
  (elk_s18_SA_Landcov_Elev2 <- occu(~Trail ~Area + Landcov + Elev + I(Elev^2), elk_s18_UMF))
  
  #'  NOTES 2/18:
  #'  NDVI, tree cov + effect
  #'  OK study area, Elev, Elev2 - effect
  
  (elk_w1819_SA <- occu(~1 ~Area, elk_w1819_UMF))
  (elk_w1819_H2o <- occu(~1 ~NearestH2o, elk_w1819_UMF))
  (elk_w1819_Landcov <- occu(~1 ~Landcov, elk_w1819_UMF))  # Landcov only partially converged
  (elk_w1819_Habitat <- occu(~1 ~Habitat, elk_w1819_UMF))  # Shrub-steppe fails to converge
  (elk_w1819_NLCD <- occu(~1 ~NLCD, elk_w1819_UMF))        # DOESNT CONVERGE
  (elk_w1819_NDVI <- occu(~1 ~NDVI, elk_w1819_UMF))
  (elk_w1819_dNBR <- occu(~1 ~dNBR, elk_w1819_UMF))
  (elk_w1819_Elev <- occu(~1 ~Elev, elk_w1819_UMF))
  (elk_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), elk_w1819_UMF))
  (elk_w1819_Slope <- occu(~1 ~Slope, elk_w1819_UMF))
  (elk_w1819_Aspect <- occu(~1 ~Aspect, elk_w1819_UMF))
  (elk_w1819_TRI <- occu(~1 ~TRI, elk_w1819_UMF))
  (elk_w1819_CanopyCov <- occu(~1 ~Canopy_cov, elk_w1819_UMF))
  (elk_w1819_TreeCov <- occu(~1 ~Tree_cov, elk_w1819_UMF))
  
  (elk_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, elk_w1819_UMF))        # NLCD FAILS TO CONVERGE
  (elk_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, elk_w1819_UMF))  # Shrub-steppe fails to converge
  (elk_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, elk_w1819_UMF))
  (elk_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, elk_w1819_UMF))
  (elk_w1819_SA_Elev <- occu(~1 ~Area + Elev, elk_w1819_UMF))
  (elk_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), elk_w1819_UMF))
  (elk_w1819_SA_Canopycov <- occu(~1 ~Area + Canopy_cov, elk_w1819_UMF))
  (elk_w1819_SA_Treecov <- occu(~1 ~Area + Tree_cov, elk_w1819_UMF))
  (elk_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, elk_w1819_UMF))
  (elk_w1819_SA_Habitat_Elev2 <- occu(~1 ~Area + Habitat + Elev + I(Elev^2), elk_w1819_UMF))
  (elk_w1819_SA_Landcov_Elev2 <- occu(~1 ~Area + Landcov + Elev + I(Elev^2), elk_w1819_UMF))
  
  #' NOTES 2/18:
  #' OK study area -
  
  
  ####  MULE DEER MODELS  ####
  #'  Null model
  (md_s18_null <- occu(~1 ~1, md_s18_UMF))
  backTransform(md_s18_null, 'det')
  backTransform(md_s18_null, 'state')
  
  (md_w1819_null <- occu(~1 ~1, md_w1819_UMF))
  backTransform(md_w1819_null, 'det')
  backTransform(md_w1819_null, 'state')
  
  #'  Detection models
  (md_s18_hgt <- occu(~Height ~1, md_s18_UMF))
  (md_s18_dist <- occu(~Distance ~1, md_s18_UMF))
  (md_s18_angle <- occu(~Height*Distance ~1, md_s18_UMF))
  (md_s18_trail <- occu(~Trail ~1, md_s18_UMF))            #  at least this converges
  (md_s18_angle.trail <- occu(~Height*Distance + Trail ~1, md_s18_UMF))
  
  (md_w1819_hgt <- occu(~Height ~1, md_w1819_UMF))
  (md_w1819_dist <- occu(~Distance ~1, md_w1819_UMF))
  (md_w1819_angle <- occu(~Height*Distance ~1, md_w1819_UMF))
  (md_w1819_trail <- occu(~Trail ~1, md_w1819_UMF))
  (md_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, md_w1819_UMF))
  
  #'  Occupancy & detection models
  (md_s18_SA <- occu(~1 ~Area, md_s18_UMF))
  (md_s18_H2o <- occu(~1 ~NearestH2o, md_s18_UMF))
  (md_s18_Landcov <- occu(~1 ~Landcov, md_s18_UMF))
  (md_s18_Habitat <- occu(~1 ~Habitat, md_s18_UMF)) 
  (md_s18_NLCD <- occu(~1 ~NLCD, md_s18_UMF))
  (md_s18_NDVI <- occu(~1 ~NDVI, md_s18_UMF))
  (md_s18_Elev <- occu(~1 ~Elev, md_s18_UMF))
  (md_s18_Elev2 <- occu(~1 ~Elev + I(Elev^2), md_s18_UMF))
  (md_s18_Slope <- occu(~1 ~Slope, md_s18_UMF))
  (md_s18_Aspect <- occu(~1 ~Aspect, md_s18_UMF))
  (md_s18_TRI <- occu(~1 ~TRI, md_s18_UMF))
  (md_s18_CanopyCov <- occu(~1 ~Canopy_cov, md_s18_UMF))
  (md_s18_TreeCov <- occu(~1 ~Tree_cov, md_s18_UMF))
  
  (md_s18_SA_NLCD <- occu(~1 ~Area + NLCD, md_s18_UMF))
  (md_s18_SA_Habitat <- occu(~1 ~Area + Habitat, md_s18_UMF))  
  (md_s18_SA_Landcov <- occu(~1 ~Area + Landcov, md_s18_UMF))
  (md_s18_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, md_s18_UMF))
  (md_s18_SA_Elev <- occu(~1 ~Area + Elev, md_s18_UMF))
  (md_s18_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), md_s18_UMF))
  (md_s18_SA_Canopycov <- occu(~1 ~Area + Canopy_cov, md_s18_UMF))
  (md_s18_SA_NDVI <- occu(~1 ~Area + NDVI, md_s18_UMF))
  (md_s18_SA_Habtiat_Elev2 <- occu(~1 ~Area + Habitat + Elev + I(Elev^2), md_s18_UMF))
  (md_s18_SA_Landcov_Elev2 <- occu(~1 ~Area + Landcov + Elev + I(Elev^2), md_s18_UMF))
  
  
  #'  NOTES 2/17:
  #'  OK study area, elev, slope, tri, canopy cov +
  #'  NDVI -
  #'  Some Habtiat, Landcov signif -
  
  
  (md_w1819_SA <- occu(~1 ~Area, md_w1819_UMF))
  (md_w1819_H2o <- occu(~1 ~NearestH2o, md_w1819_UMF))
  (md_w1819_Landcov <- occu(~1 ~Landcov, md_w1819_UMF))
  (md_w1819_Habitat <- occu(~1 ~Habitat, md_w1819_UMF))  
  (md_w1819_NLCD <- occu(~1 ~NLCD, md_w1819_UMF))
  (md_w1819_NDVI <- occu(~1 ~NDVI, md_w1819_UMF))
  (md_w1819_dNBR <- occu(~1 ~dNBR, md_w1819_UMF))
  (md_w1819_Elev <- occu(~1 ~Elev, md_w1819_UMF))
  (md_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), md_w1819_UMF))
  (md_w1819_Slope <- occu(~1 ~Slope, md_w1819_UMF))
  (md_w1819_Aspect <- occu(~1 ~Aspect, md_w1819_UMF))
  (md_w1819_TRI <- occu(~1 ~TRI, md_w1819_UMF))
  (md_w1819_CanopyCov <- occu(~1 ~Canopy_cov, md_w1819_UMF))
  (md_w1819_TreeCov <- occu(~1 ~Tree_cov, md_w1819_UMF))
  
  (md_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, md_w1819_UMF))
  (md_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, md_w1819_UMF))  
  (md_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, md_w1819_UMF))
  (md_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, md_w1819_UMF))
  (md_w1819_SA_Elev <- occu(~1 ~Area + Elev, md_w1819_UMF))
  (md_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), md_w1819_UMF))
  (md_w1819_SA_Slope <- occu(~1 ~Area + Slope, md_w1819_UMF))
  (md_w1819_SA_Aspect <- occu(~1 ~Area + Aspect, md_w1819_UMF))
  (md_w1819_SA_TRI <- occu(~1 ~Area + TRI, md_w1819_UMF))
  (md_w1819_SA_Tree_cov <- occu(~1 ~Area + Tree_cov, md_w1819_UMF))
  (md_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, md_w1819_UMF))
  (md_w1819_SA_NLCD_Elev2 <- occu(~1 ~Area + NLCD + Elev + I(Elev^2), md_w1819_UMF))
  (md_w1819_SA_Landcov_Elev2 <- occu(~1 ~Area + Landcov + Elev + I(Elev^2), md_w1819_UMF))
  (md_w1819_SA_NLCD_Slope <- occu(~1 ~Area + NLCD + Slope, md_w1819_UMF))
  (md_w1819_SA_Landcov_Slope <- occu(~1 ~Area + Landcov + Slope, md_w1819_UMF))
  
  
  #'  NOTE 2/17:
  #'  Ok study area, slope, tri +
  #'  Mixed conifer, NDVI, canopy cov, tree cov -
  #'  Some Landcov, NLCD signif
  #'  Elev - when added with study area
  
  
  ####  WHITE-TAILED DEER MODELS  ####
  #'  Null model
  (wtd_s18_null <- occu(~1 ~1, wtd_s18_UMF))
  backTransform(wtd_s18_null, 'det')
  backTransform(wtd_s18_null, 'state')
  
  (wtd_w1819_null <- occu(~1 ~1, wtd_w1819_UMF))
  backTransform(wtd_w1819_null, 'det')
  backTransform(wtd_w1819_null, 'state')
  
  #'  Detection models
  (wtd_s18_hgt <- occu(~Height ~1, wtd_s18_UMF))
  (wtd_s18_dist <- occu(~Distance ~1, wtd_s18_UMF))
  (wtd_s18_angle <- occu(~Height*Distance ~1, wtd_s18_UMF))
  (wtd_s18_trail <- occu(~Trail ~1, wtd_s18_UMF))           
  (wtd_s18_angle.trail <- occu(~Height*Distance + Trail ~1, wtd_s18_UMF))
  
  (wtd_w1819_hgt <- occu(~Height ~1, wtd_w1819_UMF))
  (wtd_w1819_dist <- occu(~Distance ~1, wtd_w1819_UMF))
  (wtd_w1819_angle <- occu(~Height*Distance ~1, wtd_w1819_UMF))
  (wtd_w1819_trail <- occu(~Trail ~1, wtd_w1819_UMF))       
  (wtd_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, wtd_w1819_UMF))
  
  #'  Notes 2/17:
  #'  trail + effect detection in winter
  
  #'  Occupancy & detection models
  (wtd_s18_SA <- occu(~1 ~Area, wtd_s18_UMF))
  (wtd_s18_H2o <- occu(~1 ~NearestH2o, wtd_s18_UMF))
  (wtd_s18_Landcov <- occu(~1 ~Landcov, wtd_s18_UMF))
  (wtd_s18_Habitat <- occu(~1 ~Habitat, wtd_s18_UMF)) 
  (wtd_s18_NLCD <- occu(~1 ~NLCD, wtd_s18_UMF))
  (wtd_s18_NDVI <- occu(~1 ~NDVI, wtd_s18_UMF))
  (wtd_s18_Elev <- occu(~1 ~Elev, wtd_s18_UMF))
  (wtd_s18_Elev2 <- occu(~1 ~Elev + I(Elev^2), wtd_s18_UMF))
  (wtd_s18_Slope <- occu(~1 ~Slope, wtd_s18_UMF))
  (wtd_s18_Aspect <- occu(~1 ~Aspect, wtd_s18_UMF))
  (wtd_s18_TRI <- occu(~1 ~TRI, wtd_s18_UMF))
  (wtd_s18_CanopyCov <- occu(~1 ~Canopy_cov, wtd_s18_UMF))
  (wtd_s18_TreeCov <- occu(~1 ~Tree_cov, wtd_s18_UMF))
  
  (wtd_s18_SA_NLCD <- occu(~1 ~Area + NLCD, wtd_s18_UMF))
  (wtd_s18_SA_Habitat <- occu(~1 ~Area + Habitat, wtd_s18_UMF))  
  (wtd_s18_SA_Landcov <- occu(~1 ~Area + Landcov, wtd_s18_UMF))
  (wtd_s18_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, wtd_s18_UMF))
  (wtd_s18_SA_Elev <- occu(~1 ~Area + Elev, wtd_s18_UMF))
  (wtd_s18_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), wtd_s18_UMF))
  (wtd_s18_SA_Tree_cov <- occu(~1 ~Area + Tree_cov, wtd_s18_UMF))
  (wtd_s18_SA_NDVI <- occu(~1 ~Area + NDVI, wtd_s18_UMF))
  (wtd_s18_SA_NLCD_Elev2 <- occu(~1 ~Area + NLCD + Elev + I(Elev^2), wtd_s18_UMF))
  
  #' NOTES 2/17:
  #' OK study area, H2o, elev, elev2, slope, TRI - effect occu in summer
  #' some Landcov, NLCD, Mixed conifer + occu in summer
  #' NDVI, canopy cov + occu in summer
  
  (wtd_w1819_SA <- occu(~Trail ~Area, wtd_w1819_UMF))
  (wtd_w1819_H2o <- occu(~Trail ~NearestH2o, wtd_w1819_UMF))
  (wtd_w1819_Landcov <- occu(~Trail ~Landcov, wtd_w1819_UMF))
  (wtd_w1819_Habitat <- occu(~Trail ~Habitat, wtd_w1819_UMF))  # DOES NOT CONVERGE AT ALL?!
  (wtd_w1819_NLCD <- occu(~Trail ~NLCD, wtd_w1819_UMF))
  (wtd_w1819_NDVI <- occu(~Trail ~NDVI, wtd_w1819_UMF))
  (wtd_w1819_dNBR <- occu(~Trail ~dNBR, wtd_w1819_UMF))
  (wtd_w1819_Elev <- occu(~Trail ~Elev, wtd_w1819_UMF))
  (wtd_w1819_Elev2 <- occu(~Trail ~Elev + I(Elev^2), wtd_w1819_UMF))
  (wtd_w1819_Slope <- occu(~Trail ~Slope, wtd_w1819_UMF))
  (wtd_w1819_Aspect <- occu(~Trail ~Aspect, wtd_w1819_UMF))
  (wtd_w1819_TRI <- occu(~Trail ~TRI, wtd_w1819_UMF))
  (wtd_w1819_CanopyCov <- occu(~Trail ~Canopy_cov, wtd_w1819_UMF))
  (wtd_w1819_TreeCov <- occu(~Trail ~Tree_cov, wtd_w1819_UMF))
  
  (wtd_w1819_SA_NLCD <- occu(~Trail ~Area + NLCD, wtd_w1819_UMF))
  (wtd_w1819_SA_Habitat <- occu(~Trail ~Area + Habitat, wtd_w1819_UMF))  # DOES NOT CONVERGE AT ALL?!
  (wtd_w1819_SA_Landcov <- occu(~Trail ~Area + Landcov, wtd_w1819_UMF))
  (wtd_w1819_SA_NearestH2o <- occu(~Trail ~Area + NearestH2o, wtd_w1819_UMF))
  (wtd_w1819_SA_Elev <- occu(~Trail ~Area + Elev, wtd_w1819_UMF))
  (wtd_w1819_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), wtd_w1819_UMF))
  (wtd_w1819_SA_Slope <- occu(~Trail ~Area + Slope, wtd_w1819_UMF))
  (wtd_w1819_SA_Aspect <- occu(~Trail ~Area + Aspect, wtd_w1819_UMF))
  (wtd_w1819_SA_TRI <- occu(~Trail ~Area + TRI, wtd_w1819_UMF))
  (wtd_w1819_SA_Tree_cov <- occu(~Trail ~Area + Tree_cov, wtd_w1819_UMF))
  (wtd_w1819_SA_NDVI <- occu(~Trail ~Area + NDVI, wtd_w1819_UMF))
  (wtd_w1819_SA_NLCD_Elev2 <- occu(~Trail ~Area + NLCD + Elev + I(Elev^2), wtd_w1819_UMF))
  
  #'  REMEMBER TO ADJUST COVS ON DETECTION AFTER RUNNING MODEL SELECTION OF DETECTION MODEL
  
  #'  NOTES 2/17:
  #'  OK study area, elev, slope, TRI - effect in winter
  #'  NDVI, aspect, canopy cov + effect in winter
  #'  Some Landcov, NDVI + effect in winter