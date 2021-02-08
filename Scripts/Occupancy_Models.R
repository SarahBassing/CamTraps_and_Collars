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
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/Camera_Station_Covariates_2021-02-05.csv") %>%
    select("Year", "Study_Area", "Cell_ID", "Camera_ID", "CameraLocation", 
           "Latitude", "Longitude", "Distance_Focal_Point", "Height_frm_grnd", 
           "Monitoring", "Canopy_Cov", "Land_Mgnt", "Land_Owner", "Habitat_Type") %>%
    mutate(
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      Cell_ID = as.factor(Cell_ID),
      Camera_ID = as.factor(Camera_ID),
      CameraLocation = as.factor(CameraLocation),
      Distance_Focal_Point = scale(Distance_Focal_Point),
      Height_frm_grnd = scale(Height_frm_grnd),
      Monitoring = as.factor(Monitoring), # consider changing 1 "closed road" to "dirt road"
      Canopy_Cov = scale(Canopy_Cov),
      Land_Mgnt = as.factor(Land_Mgnt),
      Land_Owner = as.factor(Land_Owner),
      Habitat_Type = as.factor(Habitat_Type)
    )
    
  #'  Check for correlation among covaraites
  #'  Watch out for NAs (use="complete.obs")
  cor(stations$Distance_Focal_Point, stations$Height_frm_grnd, use="complete.obs")
  
  #'  Create survey-level covariate matrix
  #'  Requires uniqe column for each sampling occasion and covariate
  nrows <- nrow(stations)
  ncols <- 13
 
  srvy_covs <- list(
    Height = matrix(c(Hgt1 = stations$Height_frm_grnd, Hgt2 = stations$Height_frm_grnd,
                      Hgt3 = stations$Height_frm_grnd, Hgt4 = stations$Height_frm_grnd,
                      Hgt5 = stations$Height_frm_grnd, Hgt6 = stations$Height_frm_grnd,
                      Hgt7 = stations$Height_frm_grnd, Hgt8 = stations$Height_frm_grnd,
                      Hgt9 = stations$Height_frm_grnd, Hgt10 = stations$Height_frm_grnd,
                      Hgt11 = stations$Height_frm_grnd, Hgt12 = stations$Height_frm_grnd,
                      Hgt13 = stations$Height_frm_grnd),
                    nrow = nrows, ncol = ncols, byrow = TRUE),
    Distance = matrix(c(Dist1 = stations$Distance_Focal_Point, Dist2 = stations$Distance_Focal_Point,
                        Dist3 = stations$Distance_Focal_Point, Dist4 = stations$Distance_Focal_Point,
                        Dist5 = stations$Distance_Focal_Point, Dist6 = stations$Distance_Focal_Point,
                        Dist7 = stations$Distance_Focal_Point, Dist8 = stations$Distance_Focal_Point,
                        Dist9 = stations$Distance_Focal_Point, Dist10 = stations$Distance_Focal_Point,
                        Dist11 = stations$Distance_Focal_Point, Dist12 = stations$Distance_Focal_Point,
                        Dist13 = stations$Distance_Focal_Point),
                      nrow = nrows, ncol = ncols, byrow = TRUE)
    )
  
  #'  Create unmarked dataframes

  ####  BOBCAT UMF  ####
  bob_s18_UMF <- unmarkedFrameOccu(DH_bob_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Monitoring,
                                                         Canopy = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type),
                                   obsCovs = srvy_covs)
  bob_w1819_UMF <- unmarkedFrameOccu(DH_bob_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Monitoring,
                                                           Canopy = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type),
                                     obsCovs = srvy_covs)
  
  summary(bob_s18_UMF)
  
  ####  COUGAR UMF  ####
  coug_s18_UMF <- unmarkedFrameOccu(DH_coug_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Monitoring,
                                                          Canopy = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type),
                                    obsCovs = srvy_covs)
  coug_w1819_UMF <- unmarkedFrameOccu(DH_coug_wtr1819,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Monitoring,
                                                          Canopy = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type),
                                    obsCovs = srvy_covs)
  
  ####  COYOTE UMF  ####
  coy_s18_UMF <- unmarkedFrameOccu(DH_coy_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Monitoring,
                                                          Canopy = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type),
                                    obsCovs = srvy_covs)
  coy_w1819_UMF <- unmarkedFrameOccu(DH_coy_wtr1819,
                                      siteCovs = data.frame(Area = stations$Study_Area,
                                                            Trail = stations$Monitoring,
                                                            Canopy = stations$Canopy_Cov,
                                                            Mgnt = stations$Land_Mgnt,
                                                            Habitat = stations$Habitat_Type),
                                      obsCovs = srvy_covs)
  
  ####  WOLF UMF  ####
  wolf_s18_UMF <- unmarkedFrameOccu(DH_wolf_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Monitoring,
                                                         Canopy = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type),
                                   obsCovs = srvy_covs)
  wolf_w1819_UMF <- unmarkedFrameOccu(DH_wolf_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Monitoring,
                                                           Canopy = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type),
                                     obsCovs = srvy_covs)
  
  ####  ELK UMF  ####
  elk_s18_UMF <- unmarkedFrameOccu(DH_elk_smr18,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Monitoring,
                                                          Canopy = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type),
                                    obsCovs = srvy_covs)
  elk_w1819_UMF <- unmarkedFrameOccu(DH_elk_wtr1819,
                                      siteCovs = data.frame(Area = stations$Study_Area,
                                                            Trail = stations$Monitoring,
                                                            Canopy = stations$Canopy_Cov,
                                                            Mgnt = stations$Land_Mgnt,
                                                            Habitat = stations$Habitat_Type),
                                      obsCovs = srvy_covs)
  
  ####  MULE DEER UMF  ####
  md_s18_UMF <- unmarkedFrameOccu(DH_md_smr18,
                                   siteCovs = data.frame(Area = stations$Study_Area,
                                                         Trail = stations$Monitoring,
                                                         Canopy = stations$Canopy_Cov,
                                                         Mgnt = stations$Land_Mgnt,
                                                         Habitat = stations$Habitat_Type),
                                   obsCovs = srvy_covs)
  md_w1819_UMF <- unmarkedFrameOccu(DH_md_wtr1819,
                                     siteCovs = data.frame(Area = stations$Study_Area,
                                                           Trail = stations$Monitoring,
                                                           Canopy = stations$Canopy_Cov,
                                                           Mgnt = stations$Land_Mgnt,
                                                           Habitat = stations$Habitat_Type),
                                     obsCovs = srvy_covs)
  
  ####  WHITE-TAILED DEER UMF  ####
  wtd_s18_UMF <- unmarkedFrameOccu(DH_wtd_smr18,
                                  siteCovs = data.frame(Area = stations$Study_Area,
                                                        Trail = stations$Monitoring,
                                                        Canopy = stations$Canopy_Cov,
                                                        Mgnt = stations$Land_Mgnt,
                                                        Habitat = stations$Habitat_Type),
                                  obsCovs = srvy_covs)
  wtd_w1819_UMF <- unmarkedFrameOccu(DH_wtd_wtr1819,
                                    siteCovs = data.frame(Area = stations$Study_Area,
                                                          Trail = stations$Monitoring,
                                                          Canopy = stations$Canopy_Cov,
                                                          Mgnt = stations$Land_Mgnt,
                                                          Habitat = stations$Habitat_Type),
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
  
  (bob_w1819_hgt <- occu(~Height ~1, bob_w1819_UMF))
  (bob_w1819_dist <- occu(~Distance ~1, bob_w1819_UMF))
  (bob_w1819_angle <- occu(~Height*Distance ~1, bob_w1819_UMF))
  (bob_w1819_trail <- occu(~Trail ~1, bob_w1819_UMF))
  (bob_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, bob_w1819_UMF))
  
  
  #'  Occupancy & detection models
  
  
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
  
  (coug_w1819_hgt <- occu(~Height ~1, coug_w1819_UMF))
  (coug_w1819_dist <- occu(~Distance ~1, coug_w1819_UMF))
  (coug_w1819_angle <- occu(~Height*Distance ~1, coug_w1819_UMF))
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
  (elk_s18_trail <- occu(~Trail ~1, elk_s18_UMF))
  (elk_s18_angle.trail <- occu(~Height*Distance + Trail ~1, elk_s18_UMF))
  
  (elk_w1819_hgt <- occu(~Height ~1, elk_w1819_UMF))
  (elk_w1819_dist <- occu(~Distance ~1, elk_w1819_UMF))
  (elk_w1819_angle <- occu(~Height*Distance ~1, elk_w1819_UMF))
  (elk_w1819_trail <- occu(~Trail ~1, elk_w1819_UMF))
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
  (md_s18_trail <- occu(~Trail ~1, md_s18_UMF))
  (md_s18_angle.trail <- occu(~Height*Distance + Trail ~1, md_s18_UMF))
  
  (md_w1819_hgt <- occu(~Height ~1, md_w1819_UMF))
  (md_w1819_dist <- occu(~Distance ~1, md_w1819_UMF))
  (md_w1819_angle <- occu(~Height*Distance ~1, md_w1819_UMF))
  (md_w1819_trail <- occu(~Trail ~1, md_w1819_UMF))
  (md_w1819_angle.trail <- occu(~Height*Distance + Trail ~1, md_w1819_UMF))
  
  #'  Occupancy & detection models
  
  
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
  
  #'  Occupancy & detection models
  
  
  
  