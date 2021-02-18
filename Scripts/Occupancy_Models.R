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
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18_2021-02-15.csv") %>%
  # stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/Camera_Station_Covariates_2021-02-05.csv") %>%
    # select("Year", "Study_Area", "Cell_ID", "Camera_ID", "CameraLocation", 
    #        "Latitude", "Longitude", "Distance_Focal_Point", "Height_frm_grnd", 
    #        "Monitoring", "Canopy_Cov", "Land_Mgnt", "Land_Owner", "Habitat_Type") %>%
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
      landcov18 = ifelse(landcov18 == "212", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "310", "211", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "332", "211", as.character(landcov18))  # Grassland / Agriculture / Residential
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
      NDVI_sp18 = scale(ndvi_sp18),
      NDVI_sm18 = scale(ndvi_sm18),
      dNBR_sp18 = scale(dnbr_sp18),
      dNBR_sm18 = scale(dnbr_sm18),
      Disturb18 = as.factor(disturb18),
      Burn18 = as.factor(burnPerim18),
      Landcov18 = as.factor(landcov18),         # Courser-scale, 30m res
      NLCD = as.factor(NLCD_landcov),           # Courser-scale, 30m res
      Elev = scale(elevation),
      Slope = scale(slope), # this is a circular variable so not quite sure how to treat this...
      Aspect = scale(aspect), # this is a circular variable and 90degrees used when slope = 0
      TRI = scale(tri),
      Roughness = scale(roughness),
      Canopy18 = scale(canopy18),               # Courser-scale, 30m res
      NearestH2o = scale(km2water)
    )

    
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
                    nrow = nrows, ncol = ncols, byrow = TRUE),
    Distance = matrix(c(Dist1 = stations$Distance, Dist2 = stations$Distance,
                        Dist3 = stations$Distance, Dist4 = stations$Distance,
                        Dist5 = stations$Distance, Dist6 = stations$Distance,
                        Dist7 = stations$Distance, Dist8 = stations$Distance,
                        Dist9 = stations$Distance, Dist10 = stations$Distance,
                        Dist11 = stations$Distance, Dist12 = stations$Distance,
                        Dist13 = stations$Distance),
                      nrow = nrows, ncol = ncols, byrow = TRUE)
    )
  #'  NEED TO BRING IN EFFORT DATA AS WELL!
  
  
  #'  Create unmarked dataframes
  ####  BOBCAT UMF  ####
  bob_s18_UMF <- unmarkedFrameOccu(DH_bob_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type,
                                                         Landcov = stations$Landcov18,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm18,  #  LOOK INTO CALCULATING CUMULATIVE NDVI  FOR SUMMER MODELS
                                                         dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         TRI = stations$TRI,
                                                         Tree_cov = stations$Canopy18,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  bob_w1819_UMF <- unmarkedFrameOccu(DH_bob_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type,
                                                           Landcov = stations$Landcov18,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                           dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           TRI = stations$TRI,
                                                           Tree_cov = stations$Canopy18,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)
  
  summary(bob_s18_UMF)
  
  ####  COUGAR UMF  ####
  coug_s18_UMF <- unmarkedFrameOccu(DH_coug_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type,
                                                          Landcov = stations$Landcov18,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm18,
                                                          dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy18,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  coug_w1819_UMF <- unmarkedFrameOccu(DH_coug_wtr1819,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type,
                                                          Landcov = stations$Landcov18,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                          dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy18,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  
  ####  COYOTE UMF  ####
  coy_s18_UMF <- unmarkedFrameOccu(DH_coy_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type,
                                                          Landcov = stations$Landcov18,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm18,
                                                          dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy18,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  coy_w1819_UMF <- unmarkedFrameOccu(DH_coy_wtr1819,
                                      siteCovs = data.frame(Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            Canopy_cov = stations$Canopy_Cov,
                                                            Mgnt = stations$Land_Mgnt,
                                                            Habitat = stations$Habitat_Type,
                                                            Landcov = stations$Landcov18,
                                                            NLCD = stations$NLCD,
                                                            NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                            dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            Aspect = stations$Aspect,
                                                            TRI = stations$TRI,
                                                            Tree_cov = stations$Canopy18,
                                                            NearestH2o = stations$NearestH2o),
                                      obsCovs = srvy_covs)
  
  ####  WOLF UMF  ####
  wolf_s18_UMF <- unmarkedFrameOccu(DH_wolf_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type,
                                                         Landcov = stations$Landcov18,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm18,
                                                         dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         TRI = stations$TRI,
                                                         Tree_cov = stations$Canopy18,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  wolf_w1819_UMF <- unmarkedFrameOccu(DH_wolf_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type,
                                                           Landcov = stations$Landcov18,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                           dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           TRI = stations$TRI,
                                                           Tree_cov = stations$Canopy18,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)
  
  ####  ELK UMF  ####
  elk_s18_UMF <- unmarkedFrameOccu(DH_elk_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type,
                                                          Landcov = stations$Landcov18,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm18,
                                                          dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy18,
                                                          NearestH2o = stations$NearestH2o),
                                    obsCovs = srvy_covs)
  elk_w1819_UMF <- unmarkedFrameOccu(DH_elk_wtr1819,
                                      siteCovs = data.frame(Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            Canopy_cov = stations$Canopy_Cov,
                                                            Mgnt = stations$Land_Mgnt,
                                                            Habitat = stations$Habitat_Type,
                                                            Landcov = stations$Landcov18,
                                                            NLCD = stations$NLCD,
                                                            NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                            dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            Aspect = stations$Aspect,
                                                            TRI = stations$TRI,
                                                            Tree_cov = stations$Canopy18,
                                                            NearestH2o = stations$NearestH2o),
                                      obsCovs = srvy_covs)
  
  ####  MULE DEER UMF  ####
  md_s18_UMF <- unmarkedFrameOccu(DH_md_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         Canopy_cov = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type,
                                                         Landcov = stations$Landcov18,
                                                         NLCD = stations$NLCD,
                                                         NDVI = stations$NDVI_sm18,
                                                         dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,
                                                         TRI = stations$TRI,
                                                         Tree_cov = stations$Canopy18,
                                                         NearestH2o = stations$NearestH2o),
                                   obsCovs = srvy_covs)
  md_w1819_UMF <- unmarkedFrameOccu(DH_md_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           Canopy_cov = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type,
                                                           Landcov = stations$Landcov18,
                                                           NLCD = stations$NLCD,
                                                           NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                           dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,
                                                           TRI = stations$TRI,
                                                           Tree_cov = stations$Canopy18,
                                                           NearestH2o = stations$NearestH2o),
                                     obsCovs = srvy_covs)
  
  ####  WHITE-TAILED DEER UMF  ####
  wtd_s18_UMF <- unmarkedFrameOccu(DH_wtd_smr18,
                                  siteCovs = data.frame(Area = stations$Study_Area,
                                                        Trail = stations$Trail,
                                                        Canopy_cov = stations$Canopy_Cov,
                                                        Mgnt = stations$Land_Mgnt,
                                                        Habitat = stations$Habitat_Type,
                                                        Landcov = stations$Landcov18,
                                                        NLCD = stations$NLCD,
                                                        NDVI = stations$NDVI_sm18,
                                                        dNBR = stations$dNBR_sm18,  #  THINK ABOUT USING BURN SEVERITY FROM PREVIOUS YEAR INSTEAD
                                                        Elev = stations$Elev,
                                                        Slope = stations$Slope,
                                                        Aspect = stations$Aspect,
                                                        TRI = stations$TRI,
                                                        Tree_cov = stations$Canopy18,
                                                        NearestH2o = stations$NearestH2o),
                                  obsCovs = srvy_covs)
  wtd_w1819_UMF <- unmarkedFrameOccu(DH_wtd_wtr1819,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          Canopy_cov = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type,
                                                          Landcov = stations$Landcov18,
                                                          NLCD = stations$NLCD,
                                                          NDVI = stations$NDVI_sm18,  #  USE SUMMER NDVI FOR WINTER MODELS
                                                          dNBR = stations$dNBR_sm18,  #  USE SUMMER BURN SEVERITY FOR WINTER MODELS
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,
                                                          TRI = stations$TRI,
                                                          Tree_cov = stations$Canopy18,
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
  bob_s18_list <- fitList(bob_s18_null, bob_s18_hgt, bob_s18_dist, bob_s18_angle, bob_s18_trail, bob_s18_angle.trail)
  modSel(bob_s18_list)
  ####  ERROR: Data are not the same among models due to missing covariate values. Consider removing NAs before analysis.
  
  (bob_w1819_hgt <- occu(~Height ~1, bob_w1819_UMF))
  (bob_w1819_dist <- occu(~Distance ~1, bob_w1819_UMF))
  (bob_w1819_angle <- occu(~Height*Distance ~1, bob_w1819_UMF))
  (bob_w1819_trail <- occu(~Trail ~1, bob_w1819_UMF))
  (bob_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, bob_w1819_UMF))
  bob_w1819_list <- fitList(bob_w1819_null, bob_w1819_hgt, bob_w1819_dist, bob_w1819_angle, bob_w1819_trail, bob_w1819_angle.trail) 
  modSel(bob_w1819_list)
  ####  ERROR: Data are not the same among models due to missing covariate values. Consider removing NAs before analysis.
  
  #'  Quick notes 2/17: game trails negatively affect detection in summer, no affect in winter
  #'  camera height & distance don't seem to matter
  
  #'  Occupancy & detection models
  (bob_s18_SA <- occu(~Trail ~Area, bob_s18_UMF))
  (bob_s18_H2o <- occu(~Trail ~NearestH2o, bob_s18_UMF))
  (bob_s18_Landcov <- occu(~Trail ~Landcov, bob_s18_UMF))
  (bob_s18_Habitat <- occu(~Trail ~Habitat, bob_s18_UMF)) # shrub-steppe doesn't converge
  (bob_s18_NLCD <- occu(~Trail ~NLCD, bob_s18_UMF))
  (bob_s18_NDVI <- occu(~Trail ~NDVI, bob_s18_UMF))
  (bob_s18_Elev <- occu(~Trail ~Elev, bob_s18_UMF))
  (bob_s18_Elev2 <- occu(~Trail ~Elev + I(Elev^2), bob_s18_UMF))
  (bob_s18_Slope <- occu(~Trail ~Slope, bob_s18_UMF))
  (bob_s18_Slope <- occu(~Trail ~Aspect, bob_s18_UMF))
  (bob_s18_TRI <- occu(~Trail ~TRI, bob_s18_UMF))
  (bob_s18_CanopyCov <- occu(~Trail ~Canopy_cov, bob_s18_UMF))
  (bob_s18_TreeCov <- occu(~Trail ~Tree_cov, bob_s18_UMF))
  
  (bob_s18_SA_NLCD <- occu(~Trail ~Area + NLCD, bob_s18_UMF))
  (bob_s18_SA_Habitat <- occu(~Trail ~Area + Habitat, bob_s18_UMF))  # shrub-steppe doesn't converge
  (bob_s18_SA_Landcov <- occu(~Trail ~Area + Landcov, bob_s18_UMF))
  (bob_s18_SA_NearestH2o <- occu(~Trail ~Area + NearestH2o, bob_s18_UMF))
  (bob_s18_SA_Elev2 <- occu(~Trail ~Area + Elev + I(Elev^2), bob_s18_UMF))
  (bob_s18_SA_Tree_cov <- occu(~Trail ~Area + Aspect, bob_s18_UMF))
  (bob_s18_SA_Tree_cov <- occu(~Trail ~Area + Tree_cov, bob_s18_UMF))
  (bob_s18_SA_NDVI <- occu(~Trail ~Area + NDVI, bob_s18_UMF))
  (bob_s18_SA_NLCD_Elev2 <- occu(~Trail ~Area + NLCD + Elev + I(Elev^2), bob_s18_UMF))
  
  #'  Note: using trail type on detection for summer models since only signif cov.
  #'  BUT need to use formal model selection to make this official
  
  #'  Quick notes 2/17: 
  #'  NDVI, Aspect + affect occupancy in summer
  #'  NLCD71 - affect occupancy in summer
  #'  Elevation + affect occupancy (signif), elev^2 - affects occu (non-signif)
  #'  Tree_cov + affect occupancy in summer
  #'  OK - affect occupancy in summer when other covs included
  
  (bob_w1819_SA <- occu(~ ~Area, bob_w1819_UMF))
  (bob_w1819_H2o <- occu(~1 ~NearestH2o, bob_w1819_UMF))
  (bob_w1819_Landcov <- occu(~1 ~Landcov, bob_w1819_UMF))
  (bob_w1819_Habitat <- occu(~1 ~Habitat, bob_w1819_UMF))  # Shrub-steppe converges
  (bob_w1819_NLCD <- occu(~1 ~NLCD, bob_w1819_UMF))
  (bob_w1819_NDVI <- occu(~1 ~NDVI, bob_w1819_UMF))
  (bob_w1819_Elev <- occu(~1 ~Elev, bob_w1819_UMF))
  (bob_w1819_Elev2 <- occu(~1 ~Elev + I(Elev^2), bob_w1819_UMF))
  (bob_w1819_Slope <- occu(~1 ~Slope, bob_w1819_UMF))
  (bob_w1819_TRI <- occu(~1 ~Aspect, bob_w1819_UMF))
  (bob_w1819_CanopyCov <- occu(~1 ~Canopy_cov, bob_w1819_UMF))
  (bob_w1819_TreeCov <- occu(~1 ~Tree_cov, bob_w1819_UMF))
  
  
  (bob_w1819_SA_NLCD <- occu(~1 ~Area + NLCD, bob_w1819_UMF))
  (bob_w1819_SA_Habitat <- occu(~1 ~Area + Habitat, bob_w1819_UMF))  # shrub-steppe converges
  (bob_w1819_SA_Landcov <- occu(~1 ~Area + Landcov, bob_w1819_UMF))
  (bob_w1819_SA_NearestH2o <- occu(~1 ~Area + NearestH2o, bob_w1819_UMF))
  (bob_w1819_SA_Elev2 <- occu(~1 ~Area + Elev + I(Elev^2), bob_w1819_UMF))
  (bob_w1819_SA_Tree_cov <- occu(~1 ~Area + Aspect, bob_w1819_UMF))
  (bob_w1819_SA_Tree_cov <- occu(~1 ~Area + Tree_cov, bob_w1819_UMF))
  (bob_w1819_SA_NDVI <- occu(~1 ~Area + NDVI, bob_w1819_UMF))
  (bob_w1819_SA_NLCD_Elev2 <- occu(~1 ~Area + NLCD + Elev + I(Elev^2), bob_w1819_UMF))
  
  
  #'  Note: using null detection on winter models since no signif cov.
  #'  BUT need to use formal model selection to make this official
  
  #'  Quick notes 2/17: 
  #'  OK study area, aspect - affect occupancy in winter
  #'  Mixed coniger + affect occupancy in winter
  #'  NDVI + affect occupancy in winter
  
  ####  COUGAR MODELS  ####
  #'  Null model
  (coug_s18_null <- occu(~1 ~1, coug_s18_UMF))
  backTransform(coug_s18_null, 'det')
  backTransform(coug_s18_null, 'state')
  
  (coug_w1819_null <- occu(~1 ~1, coug_w1819_UMF))
  backTransform(coug_w1819_null, 'det')
  backTransform(coug_w1819_null, 'state')
  
  #'  Detection models
  (coug_s18_hgt <- occu(~Height ~1, coug_s18_UMF))
  (coug_s18_dist <- occu(~Distance ~1, coug_s18_UMF))
  (coug_s18_angle <- occu(~Height*Distance ~1, coug_s18_UMF))
  (coug_s18_trail <- occu(~Trail ~1, coug_s18_UMF))
  (coug_s18_angle.trail <- occu(~Height*Distance + Trail ~1, coug_s18_UMF))
  
  (coug_w1819_hgt <- occu(~Height ~1, coug_w1819_UMF))    #  SIGNIF at 0.05
  (coug_w1819_dist <- occu(~Distance ~1, coug_w1819_UMF)) #  SIGNIF at 0.05
  (coug_w1819_angle <- occu(~Height*Distance ~1, coug_w1819_UMF)) # Dist SIGNIF at 0.06
  (coug_w1819_trail <- occu(~Trail ~1, coug_w1819_UMF))
  (coug_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, coug_w1819_UMF))
  
  
  #'  Occupancy & detection models
  
  
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
  
  #'  Occupancy & detection models
  
   
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
  
  #'  Occupancy & detection models
  
  
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
  (elk_s18_trail <- occu(~Trail ~1, elk_s18_UMF))           #  SIGNIF at 0.05/0.06, but also fails to converge
  (elk_s18_angle.trail <- occu(~Height*Distance + Trail ~1, elk_s18_UMF))
  
  (elk_w1819_hgt <- occu(~Height ~1, elk_w1819_UMF))
  (elk_w1819_dist <- occu(~Distance ~1, elk_w1819_UMF))
  (elk_w1819_angle <- occu(~Height*Distance ~1, elk_w1819_UMF))
  (elk_w1819_trail <- occu(~Trail ~1, elk_w1819_UMF))       #  SIGNIF at 0.05, but also fails to converge
  (elk_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, elk_w1819_UMF))
  #'  Occupancy & detection models
  
  
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
  (md_s18_SA_Tree_cov <- occu(~1 ~Area + Canopy_cov, md_s18_UMF))
  (md_s18_SA_NDVI <- occu(~1 ~Area + NDVI, md_s18_UMF))
  (md_s18_SA_NLCD_Elev2 <- occu(~1 ~Area + Habitat + Elev + I(Elev^2), md_s18_UMF))
  (md_s18_SA_NLCD_Elev2 <- occu(~1 ~Area + Landcov + Elev + I(Elev^2), md_s18_UMF))
  
  
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
  #'  trail + affect detection in winter
  
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
  #' OK study area, H2o, elev, elev2, slope, TRI - affect occu in summer
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