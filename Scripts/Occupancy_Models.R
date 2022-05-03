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
  library(MuMIn)
  library(condformat)
  library(tidyverse)

  #'  Source script that generates detection histories
  #'  Detection histories come trimmed to desired season length
  source("./Scripts/CameraTrap_DetectionHistories.R")

  #'  Read in covariate data collected during camera deployment & scale
  #'  Canopy cover, land management & owner, habitat type => site-level occ covs
  #'  Dist. to focal pt, height, monitoring = > survey/site-level detection covs
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18-20_2021-08-10.csv") %>% #2021-05-14 update covariate names below if using this older version
  #'  Get rid of mysterious space after one of the NEs (ugh)
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", as.character(Monitoring)),
    ) %>%
    #'  Consolidate Habitat Type into fewer categories (determined at camera site)
    mutate(
      Habitat_Type = ifelse(Habitat_Type == "Agriculture", "Grassland", as.character(Habitat_Type)), # Grassland / Agriculture
      Habitat_Type = ifelse(Habitat_Type == "Riparian", "Mixed conifer", as.character(Habitat_Type)) # Mixed conifer / Riparian
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
      PercForestMix = scale(PercForestMix2),    # Forest + Mesic Shrub mix
      PercXericShrub = scale(PercXericShrub),
      PercXericGrass = scale(PercXericGrass),
      Elev = scale(Elev), 
      Slope = scale(Slope), 
      RoadDensity = scale(RoadDen), 
      HumanModified = scale(HumanMod)
      ) %>%
    arrange(Year) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER 2021-08-10 version

  #'  Adjust reference category for Trail factors
  order_trail <- c("Trail", "Dirt road", "Decommissioned road")
  stations <- stations %>%
    mutate(
      Trail = fct_relevel(Trail, order_trail)
    )
    
  #'  Identify which sites have missing data (just look at 1st column)
  noCanopy <- which(is.na(stations$Canopy_Cov))
  noDist <- which(is.na(stations$Distance))
  noHgt <- which(is.na(stations$Height))
  missing_dat <- unique(c(noHgt, noDist, noCanopy))
  
  #'  Find median value of these covaraites (mean = 0, sd = 1 since center & scaled)
  median(stations$Distance, na.rm = TRUE)
  median(stations$Height, na.rm = TRUE)
  median(stations$Canopy_Cov, na.rm = TRUE) # data skewed due to many cameras in open/burned areas
  
  #'  Replace missing covariate values with mean value
  #'  Missing 2018 Canopy_Cov obs in forested are so using mean instead of median
  stations$Distance[is.na(stations$Distance),] <- 0
  stations$Height[is.na(stations$Height),] <- 0
  stations$Canopy_Cov[is.na(stations$Canopy_Cov),] <- 0
  
  #'  Save
  # write.csv(stations, paste0('./Outputs/Tables/stations18-20_', Sys.Date(), '.csv'))
    
  #'  Save study-area specific covariates (important for ungulate models)
  stations_NE <- filter(stations, Study_Area == "NE")
  stations_OK <- filter(stations, Study_Area == "OK")
  
  #'  Check for correlation among covariates
  #'  Watch out for NAs (use = "complete.obs")
  cov_data <- stations[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity", "HumanModified")]
  (cov_matrix <- cor(cov_data, use = "complete.obs")) # Elev & HM -0.6407; Forest & Grass -0.6348
  cov_data_NE <- stations_NE[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity", "HumanModified")]
  (cov_matrix_NE <- cor(cov_data_NE, use = "complete.obs"))
  cov_data_OK <- stations_OK[ , c("Elev", "Slope", "PercForestMix", "PercXericGrass", "PercXericShrub", "RoadDensity", "HumanModified")]
  (cov_matrix_OK <- cor(cov_data_OK, use = "complete.obs")) # Elev & HM -0.6887
  #'  Elevation and Human Modification are pretty closely correlated, especially 
  #'  in the Okanogan, although they're not quite high enough to keep separate
  #'  based on my arbitrary r = |0.7| cutoff. Will decide to keep or drop HM 
  #'  based on level of correlation in RSF data.
  
  #'  Survey covariates
  #'  Read in & format weather data that was not included in covariates above
  narr <- read.csv("./Data/NARR_weeklymeans_smr18-wtr20.csv") 
  temp <- narr %>%
    dplyr::select(-c(X, MeanDPrecip_mm)) %>%
    group_by(Study_Area, Year, Season) %>%
    spread(Occasion, MeanDTemp_K) %>%
    ungroup() %>%
    arrange(Year, CameraLocation)
  temp <- as.data.frame(temp)
  
  #'  Break up temp data in summer and winter data frames
  temp_smr <- filter(temp, Season == "Summer18" | Season == "Summer19") 
  #'  Create rows for cameras that are missing from seasonal temp data
  #'  Happens when a camera is deployed after specific data range
  missing_smr <- subset(stations, !(CameraLocation %in% temp_smr$CameraLocation)) %>%
    dplyr::select(CameraLocation, Study_Area, Year) %>%
    mutate(
      Season = ifelse(Year == "Year1", "Summer18", "Summer19"),
    )
  missing_occ <- matrix(NA, nrow = nrow(missing_smr), ncol = 13)
  missing_smr <- cbind(missing_smr, missing_occ)
  temp_smr <- rbind(temp_smr, missing_smr) %>%
    arrange(Year, CameraLocation)
  #'  Format for UMFs and standardize each column
  Temp_smr <- temp_smr %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
    
  temp_wtr <- filter(temp, Season == "Winter1819" | Season == "Winter1920") 
  #'  Create rows for cameras that are missing from seasonal temp data
  #'  Happens when a camera was moved part way through the year
  missing_wtr <- subset(stations, !(CameraLocation %in% temp_wtr$CameraLocation)) %>%
    dplyr::select(CameraLocation, Study_Area, Year) %>%
    mutate(
      Season = ifelse(Year == "Year1", "Winter1819", "Winter1920"),
    )
  missing_occ <- matrix(NA, nrow = nrow(missing_wtr), ncol = 13)
  missing_wtr <- cbind(missing_wtr, missing_occ)
  temp_wtr <- rbind(temp_wtr, missing_wtr) %>%
    arrange(Year, CameraLocation)
  #'  Format for UMFs and standardize each column
  Temp_wtr <- temp_wtr %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
  
  
  #'  Create survey-level covariate matrix
  #'  Requires unique column for each sampling occasion and covariate
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
                      nrow = nrows, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr,
    Temp_wtr = Temp_wtr
    )
  # Effort_smr = matrix(Effort_smr1819, nrow = nrows, ncol = ncols, byrow = FALSE),
  # Effort_wtr = matrix(Effort_wtr1820, nrow = nrows, ncol = ncols, byrow = FALSE)
  #'  FYI effort covariate does weird things when scaled so not doing it now
  #'  Very few sampling occasions had low sampling effort (~3%; ignoring occasions
  #'  when camera failed completely) so not including sampling effort as a 
  #'  covariate on detection process- not enough variation to estimate effect
  
  #'  Double check it looks OK
  head(srvy_covs[[2]])


  #'  Save study-area specific survey covariates
  Hgt_NE <- stations$Height[stations$Study_Area == "NE"]
  Hgt_OK <- stations$Height[stations$Study_Area == "OK"]
  Dist_NE <- stations$Distance[stations$Study_Area == "NE"]
  Dist_OK <- stations$Distance[stations$Study_Area == "OK"]
  Temp_smr_NE <- temp_smr[temp_smr$Study_Area == "NE",] %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
  Temp_smr_OK <- temp_smr[temp_smr$Study_Area == "OK",] %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
  Temp_wtr_NE <- temp_wtr[temp_wtr$Study_Area == "NE",] %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
  Temp_wtr_OK <- temp_wtr[temp_wtr$Study_Area == "OK",] %>%
    dplyr::select(-c(CameraLocation, Study_Area, Year, Season)) %>%
    scale(.)
  
  nrows_NE <- nrow(stations_NE)
  nrows_OK <- nrow(stations_OK)
  
  srvy_covs_NE <- list(
    Height = matrix(c(Hgt1 = Hgt_NE, Hgt2 = Hgt_NE, Hgt3 = Hgt_NE, Hgt4 = Hgt_NE,
                      Hgt5 = Hgt_NE, Hgt6 = Hgt_NE, Hgt7 = Hgt_NE, Hgt8 = Hgt_NE,
                      Hgt9 = Hgt_NE, Hgt10 = Hgt_NE, Hgt11 = Hgt_NE, Hgt12 = Hgt_NE,
                      Hgt13 = Hgt_NE),
                    nrow = nrows_NE, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = Dist_NE, Dist2 = Dist_NE, Dist3 = Dist_NE, 
                        Dist4 = Dist_NE, Dist5 = Dist_NE, Dist6 = Dist_NE,
                        Dist7 = Dist_NE, Dist8 = Dist_NE, Dist9 = Dist_NE, 
                        Dist10 = Dist_NE, Dist11 = Dist_NE, Dist12 = Dist_NE,
                        Dist13 = Dist_NE),
                      nrow = nrows_NE, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr_NE,
    Temp_wtr = Temp_wtr_NE
  )
  
  srvy_covs_OK <- list(
    Height = matrix(c(Hgt1 = Hgt_OK, Hgt2 = Hgt_OK, Hgt3 = Hgt_OK, Hgt4 = Hgt_OK,
                      Hgt5 = Hgt_OK, Hgt6 = Hgt_OK, Hgt7 = Hgt_OK, Hgt8 = Hgt_OK,
                      Hgt9 = Hgt_OK, Hgt10 = Hgt_OK, Hgt11 = Hgt_OK, Hgt12 = Hgt_OK,
                      Hgt13 = Hgt_OK),
                    nrow = nrows_OK, ncol = ncols, byrow = FALSE),
    Distance = matrix(c(Dist1 = Dist_OK, Dist2 = Dist_OK, Dist3 = Dist_OK, 
                        Dist4 = Dist_OK, Dist5 = Dist_OK, Dist6 = Dist_OK,
                        Dist7 = Dist_OK, Dist8 = Dist_OK, Dist9 = Dist_OK, 
                        Dist10 = Dist_OK, Dist11 = Dist_OK, Dist12 = Dist_OK,
                        Dist13 = Dist_OK),
                      nrow = nrows_OK, ncol = ncols, byrow = FALSE),
    Temp_smr = Temp_smr_OK,
    Temp_wtr = Temp_wtr_OK
  )
  
  #'  Double check these
  head(srvy_covs_NE[[1]])
  tail(srvy_covs_OK[[1]])

  
  #'  Create unmarked dataframes
  #'  FYI: unmarkedFrameOccu is NOT smart enough to match up DH and covs by CameraLocation
  #'  so BE SURE DH AND STATIONS CAMERA SITES ARE ORDERED IN THE SAME WAY!!!!
  
  ####  BOBCAT UMF  ####
  bob_s1819_UMF <- unmarkedFrameOccu(DH_bob_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         PercForMix = stations$PercForestMix,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         RoadDensity = stations$RoadDensity,
                                                         HumanMod = stations$HumanModified),
                                   obsCovs = srvy_covs)

  nrow(bob_s1819_UMF@y)
  #'  Double check the order of CameraLocations in stations matches the DH
  # bob_station_arranged <- stations
  # bob_covs_arranged <- bob_s1819_UMF@siteCovs
  # bob_DH_arranged <- bob_s1819_UMF@y
  #'  Remove rows with missing observation covariate data (Height & Distance data)
  # bob_s1819_UMF <- bob_s1819_UMF[-missing_dat]
  # nrow(bob_s1819_UMF@y)

  bob_w1820_UMF <- unmarkedFrameOccu(DH_bob_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           PercForMix = stations$PercForestMix,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           RoadDensity = stations$RoadDensity,
                                                           HumanMod = stations$HumanModified),
                                     obsCovs = srvy_covs)

  nrow(bob_w1820_UMF@y)
  #'  Remove rows with missing observation covariate data (Height & Distance data)
  #' bob_w1820_UMF <- bob_w1820_UMF[-missing_dat]
  #' nrow(bob_w1820_UMF@y)
  
  summary(bob_w1820_UMF)
  
  ####  COUGAR UMF  ####
  coug_s1819_UMF <- unmarkedFrameOccu(DH_coug_smr1819,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanMod = stations$HumanModified),
                                    obsCovs = srvy_covs)
  nrow(coug_s1819_UMF@y)
  # coug_s1819_UMF <- coug_s1819_UMF[-missing_dat]
  # nrow(coug_s1819_UMF@y)
  
  coug_w1820_UMF <- unmarkedFrameOccu(DH_coug_wtr1820,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanMod = stations$HumanModified),
                                    obsCovs = srvy_covs)
  nrow(coug_w1820_UMF@y)
  # coug_w1820_UMF <- coug_w1820_UMF[-missing_dat]
  # nrow(coug_w1820_UMF@y)
  
  
  ####  COYOTE UMF  ####
  coy_s1819_UMF <- unmarkedFrameOccu(DH_coy_smr1819,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          PercForMix = stations$PercForestMix,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanMod = stations$HumanModified),
                                    obsCovs = srvy_covs)
  # coy_s1819_UMF <- coy_s1819_UMF[-missing_dat]
  
  coy_w1820_UMF <- unmarkedFrameOccu(DH_coy_wtr1820,
                                      siteCovs = data.frame(Year = stations$Year,
                                                            Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            PercForMix = stations$PercForestMix,
                                                            PercXGrass = stations$PercXericGrass,
                                                            PercXShrub = stations$PercXericShrub,
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            RoadDensity = stations$RoadDensity,
                                                            HumanMod = stations$HumanModified),
                                      obsCovs = srvy_covs)
  # coy_w1820_UMF <- coy_w1820_UMF[-missing_dat]
  
  ####  WOLF UMF  ####
  wolf_s1819_UMF <- unmarkedFrameOccu(DH_wolf_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         PercForMix = stations$PercForestMix,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         RoadDensity = stations$RoadDensity,
                                                         HumanMod = stations$HumanModified),
                                   obsCovs = srvy_covs)
  # wolf_s1819_UMF <- wolf_s1819_UMF[-missing_dat]
  
  wolf_w1820_UMF <- unmarkedFrameOccu(DH_wolf_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           PercForMix = stations$PercForestMix,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           RoadDensity = stations$RoadDensity,
                                                           HumanMod = stations$HumanModified),
                                     obsCovs = srvy_covs)
  # wolf_w1820_UMF <- wolf_w1820_UMF[-missing_dat]
  
  
  ####  ELK UMF  ####
  #'  Data from NE study area only to mirror GPS collar data
  elk_s1819_UMF <- unmarkedFrameOccu(DH_elk_smr1819,
                                    siteCovs = data.frame(Year = stations_NE$Year,
                                                          Trail = stations_NE$Trail,
                                                          PercForMix = stations_NE$PercForestMix,
                                                          PercXGrass = stations_NE$PercXericGrass,
                                                          PercXShrub = stations_NE$PercXericShrub,
                                                          Elev = stations_NE$Elev,
                                                          Slope = stations_NE$Slope,
                                                          RoadDensity = stations_NE$RoadDensity,
                                                          HumanMod = stations_NE$HumanModified),
                                    obsCovs = srvy_covs_NE)
  # elk_s1819_UMF <- elk_s1819_UMF[-missing_dat]
  
  elk_w1820_UMF <- unmarkedFrameOccu(DH_elk_wtr1820,
                                      siteCovs = data.frame(Year = stations_NE$Year,
                                                            Trail = stations_NE$Trail,
                                                            PercForMix = stations_NE$PercForestMix,
                                                            PercXGrass = stations_NE$PercXericGrass,
                                                            PercXShrub = stations_NE$PercXericShrub,
                                                            Elev = stations_NE$Elev,
                                                            Slope = stations_NE$Slope,
                                                            RoadDensity = stations_NE$RoadDensity,
                                                            HumanMod = stations_NE$HumanModified),
                                      obsCovs = srvy_covs_NE)
  # elk_w1820_UMF <- elk_w1820_UMF[-missing_dat]
  
  ####  MULE DEER UMF  ####
  #'  Data from OK study area only to mirror GPS collar data
  md_s1819_UMF <- unmarkedFrameOccu(DH_md_smr1819,
                                   siteCovs = data.frame(Year = stations_OK$Year,
                                                         Trail = stations_OK$Trail,
                                                         PercForMix = stations_OK$PercForestMix,
                                                         PercXGrass = stations_OK$PercXericGrass,
                                                         PercXShrub = stations_OK$PercXericShrub,
                                                         Elev = stations_OK$Elev,
                                                         Slope = stations_OK$Slope,
                                                         RoadDensity = stations_OK$RoadDensity,
                                                         HumanMod = stations_OK$HumanModified),
                                   obsCovs = srvy_covs_OK)
  # md_s1819_UMF <- md_s1819_UMF[-missing_dat]
  
  md_w1820_UMF <- unmarkedFrameOccu(DH_md_wtr1820,
                                     siteCovs = data.frame(Year = stations_OK$Year,
                                                           Trail = stations_OK$Trail,
                                                           PercForMix = stations_OK$PercForestMix,
                                                           PercXGrass = stations_OK$PercXericGrass,
                                                           PercXShrub = stations_OK$PercXericShrub,
                                                           Elev = stations_OK$Elev,
                                                           Slope = stations_OK$Slope,
                                                           RoadDensity = stations_OK$RoadDensity,
                                                           HumanMod = stations_OK$HumanModified),
                                     obsCovs = srvy_covs_OK)
  # md_w1820_UMF <- md_w1820_UMF[-missing_dat]
  
  ####  WHITE-TAILED DEER UMF  ####
  #'  Data from NE study area only to mirror GPS collar data
  wtd_s1819_UMF <- unmarkedFrameOccu(DH_wtd_smr1819,
                                  siteCovs = data.frame(Year = stations_NE$Year,
                                                        Trail = stations_NE$Trail,
                                                        PercForMix = stations_NE$PercForestMix,
                                                        PercXGrass = stations_NE$PercXericGrass,
                                                        PercXShrub = stations_NE$PercXericShrub,
                                                        Elev = stations_NE$Elev,
                                                        Slope = stations_NE$Slope,
                                                        RoadDensity = stations_NE$RoadDensity,
                                                        HumanMod = stations_NE$HumanModified),
                                  obsCovs = srvy_covs_NE)
  # wtd_s1819_UMF <- wtd_s1819_UMF[-missing_dat]
  
  wtd_w1820_UMF <- unmarkedFrameOccu(DH_wtd_wtr1820,
                                    siteCovs = data.frame(Year = stations_NE$Year,
                                                          Trail = stations_NE$Trail,
                                                          PercForMix = stations_NE$PercForestMix,
                                                          PercXGrass = stations_NE$PercXericGrass,
                                                          PercXShrub = stations_NE$PercXericShrub,
                                                          Elev = stations_NE$Elev,
                                                          Slope = stations_NE$Slope,
                                                          RoadDensity = stations_NE$RoadDensity,
                                                          HumanMod = stations_NE$HumanModified),
                                    obsCovs = srvy_covs_NE)
  # wtd_w1820_UMF <- wtd_w1820_UMF[-missing_dat]
  nrow( wtd_w1820_UMF@y)
  
  
  #'  Occupancy models
  #'  =============================
  #'  unmarked formula: ~detection ~occupancy
  #'  ~1 for intercept only
  #'  
  #'  1) Run GLOBAL MODEL for all species and seasons.
  #'  Only removed Study Area covariate if data only collected in one area.
  #'  Only removed PercXShrub if converged poorly (usually due to no/few 
  #'  detections in areas with shrubland).
  #'  2) Calculate VIF (Variance Inflation Factor) for the global model to identify
  #'  if there's any massive correlation once all covariates are included in the
  #'  model. VIF = 1 means no correlation among any predictors; 1 < VIF < 5 means
  #'  moderate correlation among predictors but should not be a problem; VIF > 5
  #'  means severe correlation between predictors- coefficients & p-value likely
  #'  unreliable!
  #'  3) DREDGE global model to identify "best" model for each species & season.
  #'  Goal is to identify which covariates are important in "best" model & compare 
  #'  that to covariates that were significant in the global model. What am I
  #'  missing when I use the global model? Use this to help determine a reasonable
  #'  cutoff when interpreting p-values. Consider anything p < 0.1 to have trend
  #'  worth interpreting given so many non-significant variables in the model 
  #'  can drag down the certainty around more important covariates.
  #'  NOTE: dredging with this many variables TAKES FOREVER so only do once!!!
  #'  =============================

  ####  BOBCAT MODELS  ####                   
  #'  SUMMERS 2018 & 2019
  (bob_s1819_global <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                            ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area, 
                            bob_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(bob_s1819_global, type = "state")
  unmarked::vif(bob_s1819_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' bob_s1819_dd <- dredge(bob_s1819_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_bob_smr_dd <- nrow(bob_s1819_dd)
  #' print(bob_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' bob_s1819_all <- get.models(bob_s1819_dd, subset = delta < 2,)
  #' bob_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (bob_s1819_top <- occu(formula = ~Distance + Height + Trail + Distance:Height  ~Elev + PercForMix + Slope, data = bob_s1819_UMF))
  
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXShrub to be consistent with RSF           
  (bob_w1820_global <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                            ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + Area, 
                            bob_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(bob_w1820_global, type = "state")
  unmarked::vif(bob_w1820_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' bob_w1820_dd <- dredge(bob_w1820_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_bob_wtr_dd <- nrow(bob_w1820_dd)
  #' print(bob_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' bob_w1820_all <- get.models(bob_w1820_dd, subset = delta < 2,)
  #' bob_w1820_all[[1]]
  #' bob_w1820_all[[2]]
  #' #'  Dredge identified top model
  #' (bob_w1820_top <- occu(formula = ~Distance + Trail ~PercForMix + RoadDensity, data = bob_w1820_UMF))
  #' (bob_w1820_top <- occu(formula = ~Distance ~Elev + PercForMix, data = bob_w1820_UMF))

  
  ####  COUGAR MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coug_s1819_global <- occu(~Trail + Temp_smr + Height + Distance + Height*Distance 
                             ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area, 
                             coug_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coug_s1819_global, type = "state")
  unmarked::vif(coug_s1819_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' coug_s1819_dd <- dredge(coug_s1819_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_cg_smr_dd <- nrow(coug_s1819_dd)
  #' print(coug_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' coug_s1819_all <- get.models(coug_s1819_dd, subset = delta < 2,)
  #' coug_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (coug_s1819_top <- occu(formula = ~Height + Trail ~Area + Elev + PercForMix, data = coug_s1819_UMF))
  #' #'  Uni-variate model (on occupancy) for elevation to make sure opposing OccMod & RSF effects aren't due to over-parameterization or confounding variables in OccMod
  #' (coug_s1819_elev <- occu(formula = ~Trail + Temp_smr + Height + Distance + Height*Distance ~Elev, data = coug_s1819_UMF)) # elev has + effect

  #'  WINTERS 2018-2019 & 2019-2020     
  #'  Dropping PercXGrass due to high correlation with PercForMix in RSF      
  (coug_w1820_global <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                             ~Elev + Slope + PercForMix + PercXShrub + RoadDensity + Area, 
                             coug_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coug_w1820_global, type = "state") 
  unmarked::vif(coug_w1820_global, type = "det")
  #' #'  Dredge the global model for all possible combinations  
  #' coug_w1820_dd <- dredge(coug_w1820_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_cg_wtr_dd <- nrow(coug_w1820_dd)
  #' print(coug_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' coug_w1820_all <- get.models(coug_w1820_dd, subset = delta < 2,)
  #' coug_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (coug_w1820_top <- occu(formula = ~Height + Trail ~Elev + HumanMod + PercForMix + Slope, data = coug_w1820_UMF))  
  
  ####  COYOTE MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coy_s1819_global <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                            ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + Area, 
                            coy_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coy_s1819_global, type = "state") 
  unmarked::vif(coy_s1819_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' coy_s1819_dd <- dredge(coy_s1819_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_cy_smr_dd <- nrow(coy_s1819_dd)
  #' print(coy_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' coy_s1819_all <- get.models(coy_s1819_dd, subset = delta < 2,)
  #' coy_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (coy_s1819_top <- occu(formula = ~Height + Trail ~Area + Elev + Slope, data = coy_s1819_UMF))
   
  #'  WINTERS 2018-2019 & 2019-2020    
  #'  Dropping PercXGrass to be consistent with RSF          
  (coy_w1820_global <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                            ~Elev + Slope + PercForMix + PercXShrub + RoadDensity + Area, 
                            coy_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(coy_w1820_global, type = "state") 
  unmarked::vif(coy_w1820_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' coy_w1820_dd <- dredge(coy_w1820_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_cy_wtr_dd <- nrow(coy_w1820_dd)
  #' print(coy_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' coy_w1820_all <- get.models(coy_w1820_dd, subset = delta < 2,)
  #' coy_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (coy_w1820_top <- occu(formula = ~Distance + Height + Temp_wtr + Trail + Distance:Height ~Elev + PercXGrass, data = coy_w1820_UMF))


  
  ####  WOLF MODELS  ####
  #'  SUMMERS 2018 & 2019    
  #'  Removed PercXShrub due to poor convergence                 
  (wolf_s1819_global2 <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                              ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + Area, 
                              wolf_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wolf_s1819_global2, type = "state") 
  unmarked::vif(wolf_s1819_global2, type = "det")
  #' #'  Dredge the global2 model for all possible combinations
  #' wolf_s1819_dd <- dredge(wolf_s1819_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_wf_smr_dd <- nrow(wolf_s1819_dd)
  #' print(wolf_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' wolf_s1819_all <- get.models(wolf_s1819_dd, subset = delta < 2,)
  #' wolf_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (wolf_s1819_top <- occu(formula = ~Trail + Temp_smr ~Area + Elev, data = wolf_s1819_UMF)) 
  
  #'  WINTERS 2018-2019 & 2019-2020 
  #'  Removed PercXShrub due to poor convergence 
  #'  Removed PercXGrass to be consistent with RSF   
  (wolf_w1820_global2 <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                              ~Elev + Slope + PercForMix + RoadDensity + Area, 
                              wolf_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wolf_w1820_global2, type = "state") 
  unmarked::vif(wolf_w1820_global2, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' wolf_w1820_dd <- dredge(wolf_w1820_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_wf_wtr_dd <- nrow(wolf_w1820_dd)
  #' print(wolf_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' wolf_w1820_all <- get.models(wolf_w1820_dd, subset = delta < 2,)
  #' wolf_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (wolf_w1820_top <- occu(formula = ~Trail ~Area, data = wolf_w1820_UMF))
  
  
  ####  ELK MODELS ####                          
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  #'  Removed PercXGrass and PercXShrub models due to poor convergence                             
  (elk_s1819_global2 <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                             ~Elev + Slope + PercForMix + RoadDensity, 
                             elk_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(elk_s1819_global2, type = "state") 
  unmarked::vif(elk_s1819_global2, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' elk_s1819_dd <- dredge(elk_s1819_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_ek_smr_dd <- nrow(elk_s1819_dd)
  #' print(elk_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' elk_s1819_all <- get.models(elk_s1819_dd, subset = delta < 2,)
  #' elk_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (elk_s1819_top <- occu(formula = ~Distance + Height + Trail ~1, data = elk_s1819_UMF))
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  NE study area only so no Area effect 
  #'  Removed PercXGrass and PercXShrub due to poor convergence               
  (elk_w1820_global2 <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                             ~Elev + Slope + PercForMix + RoadDensity, 
                             elk_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(elk_w1820_global2, type = "state") 
  unmarked::vif(elk_w1820_global2, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' elk_w1820_dd <- dredge(elk_w1820_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_ek_wtr_dd <- nrow(elk_w1820_dd)
  #' print(elk_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' elk_w1820_all <- get.models(elk_w1820_dd, subset = delta < 2,)
  #' elk_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (elk_w1820_top <- occu(formula = ~Height ~1, data = elk_w1820_UMF))
  
  
  ####  MULE DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  OK study area only so no Area effect
  #'  Dropped PercXShrub b/c excluded in RSF due to correlation with PercForMix
  (md_s1819_global <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                           ~Elev + Slope + PercForMix + PercXGrass + RoadDensity, #PercXShrub + 
                           md_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(md_s1819_global, type = "state") 
  unmarked::vif(md_s1819_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' md_s1819_dd <- dredge(md_s1819_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_md_smr_dd <- nrow(md_s1819_dd)
  #' print(md_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' md_s1819_all <- get.models(md_s1819_dd, subset = delta < 2,)
  #' md_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (md_s1819_top <- occu(formula = ~Distance + Temp_smr ~PercForMix, data = md_s1819_UMF))

  #'  WINTERS 2018-2019 & 2019-2020, OK study area only so no Area effect                    
  (md_w1820_global <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                           ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity, 
                           md_w1820_UMF))  
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(md_w1820_global, type = "state") 
  unmarked::vif(md_w1820_global, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' md_w1820_dd <- dredge(md_w1820_global, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_md_wtr_dd <- nrow(md_w1820_dd)
  #' print(md_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' md_w1820_all <- get.models(md_w1820_dd, subset = delta < 2,)
  #' md_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (md_w1820_top <- occu(formula = ~Distance + Height + Trail + Distance:Height ~PercXGrass + PercXShrub, data = md_w1820_UMF))
  
  
  ####  WHITE-TAILED DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  #'  Removed PercXGrass & PercXShrub due to poor convergence
  (wtd_s1819_global2 <- occu(~Trail + Temp_smr + Height + Distance + Distance*Height 
                             ~Elev + Slope + PercForMix +  RoadDensity, 
                             wtd_s1819_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wtd_s1819_global2, type = "state") 
  unmarked::vif(wtd_s1819_global2, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' wtd_s1819_dd <- dredge(wtd_s1819_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_wt_smr_dd <- nrow(wtd_s1819_dd)
  #' print(wtd_s1819_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' wtd_s1819_all <- get.models(wtd_s1819_dd, subset = delta < 2,)
  #' wtd_s1819_all[[1]]
  #' #'  Dredge identified top model
  #' (wtd_s1819_top <- occu(formula = ~Height + Trail ~Elev + PercForMix, data = wtd_s1819_UMF))
  
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only so no Area effect     
  #'  Removed PercXGrass & PercXShrub due to poor convergence         
  (wtd_w1820_global2 <- occu(~Trail + Temp_wtr + Height + Distance + Distance*Height 
                             ~Elev + Slope + PercForMix + RoadDensity, 
                             wtd_w1820_UMF))
  #'  Calculate variance inflation factor for each sub-model
  unmarked::vif(wtd_w1820_global2, type = "state") 
  unmarked::vif(wtd_w1820_global2, type = "det")
  #' #'  Dredge the global model for all possible combinations
  #' wtd_w1820_dd <- dredge(wtd_w1820_global2, rank = "AIC")
  #' #'  Count the number of dredged models
  #' tot_wt_wtr_dd <- nrow(wtd_w1820_dd)
  #' print(wtd_w1820_dd[1:5,])
  #' #'  Keep top models (within 2 deltaAIC) & review the top model
  #' wtd_w1820_all <- get.models(wtd_w1820_dd, subset = delta < 2,)
  #' wtd_w1820_all[[1]]
  #' #'  Dredge identified top model
  #' (wtd_w1820_top <- occu(formula = ~Distance + Height + Temp_wtr + Trail + Distance:Height ~Elev, data = wtd_w1820_UMF))
  

  #' #'  Average number of dredged models
  #' ndd <- sum(tot_bob_smr_dd, tot_bob_wtr_dd, tot_cg_smr_dd, tot_cg_wtr_dd, tot_cy_smr_dd, 
  #'   tot_cy_wtr_dd, tot_ek_smr_dd, tot_ek_wtr_dd, tot_md_smr_dd, tot_md_wtr_dd, 
  #'   tot_wf_smr_dd, tot_wf_wtr_dd, tot_wt_smr_dd, tot_wt_wtr_dd)
  #' meandd <- sum(tot_bob_smr_dd, tot_bob_wtr_dd, tot_cg_smr_dd, tot_cg_wtr_dd, tot_cy_smr_dd, 
  #'   tot_cy_wtr_dd, tot_ek_smr_dd, tot_ek_wtr_dd, tot_md_smr_dd, tot_md_wtr_dd, 
  #'   tot_wf_smr_dd, tot_wf_wtr_dd, tot_wt_smr_dd, tot_wt_wtr_dd)/14
  #' mindd <- min(tot_bob_smr_dd, tot_bob_wtr_dd, tot_cg_smr_dd, tot_cg_wtr_dd, tot_cy_smr_dd, 
  #'   tot_cy_wtr_dd, tot_ek_smr_dd, tot_ek_wtr_dd, tot_md_smr_dd, tot_md_wtr_dd, 
  #'   tot_wf_smr_dd, tot_wf_wtr_dd, tot_wt_smr_dd, tot_wt_wtr_dd)
  #' maxdd <- max(tot_bob_smr_dd, tot_bob_wtr_dd, tot_cg_smr_dd, tot_cg_wtr_dd, tot_cy_smr_dd, 
  #'   tot_cy_wtr_dd, tot_ek_smr_dd, tot_ek_wtr_dd, tot_md_smr_dd, tot_md_wtr_dd, 
  #'   tot_wf_smr_dd, tot_wf_wtr_dd, tot_wt_smr_dd, tot_wt_wtr_dd)
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info

  #'  Function to save occupancy results
  occ_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
        # Model = rep(model, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
      # relocate(Model, .before = Species)
    return(out)
  }
  
  #'  Run each model through function
  #'  Full models
  bob_s1819_occ <- occ_out(bob_s1819_global, "Bobcat", "Summer") #, "Global"
  bob_w1820_occ <- occ_out(bob_w1820_global, "Bobcat", "Winter")
  coug_s1819_occ <- occ_out(coug_s1819_global, "Cougar", "Summer")
  coug_w1820_occ <- occ_out(coug_w1820_global, "Cougar", "Winter")
  coy_s1819_occ <- occ_out(coy_s1819_global, "Coyote", "Summer")
  coy_w1820_occ <- occ_out(coy_w1820_global, "Coyote", "Winter")
  wolf_s1819_occ <- occ_out(wolf_s1819_global2, "Wolf", "Summer")
  wolf_w1820_occ <- occ_out(wolf_w1820_global2, "Wolf", "Winter")
  elk_s1819_occ <- occ_out(elk_s1819_global2, "Elk", "Summer")
  elk_w1820_occ <- occ_out(elk_w1820_global2, "Elk", "Winter")
  md_s1819_occ <- occ_out(md_s1819_global, "Mule Deer", "Summer")
  md_w1820_occ <- occ_out(md_w1820_global, "Mule Deer", "Winter")
  wtd_s1819_occ <- occ_out(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_w1820_occ <- occ_out(wtd_w1820_global2, "White-tailed Deer", "Winter")
  #'  Top models identified by dredging
  # bob_s1819_occt <- occ_out(bob_s1819_top, "Bobcat", "Summer", "Top")
  # bob_w1820_occt <- occ_out(bob_w1820_top, "Bobcat", "Winter", "Top")
  # coug_s1819_occt <- occ_out(coug_s1819_top, "Cougar", "Summer", "Top")
  # coug_w1820_occt <- occ_out(coug_w1820_top, "Cougar", "Winter", "Top")
  # coy_s1819_occt <- occ_out(coy_s1819_top, "Coyote", "Summer", "Top")
  # coy_w1820_occt <- occ_out(coy_w1820_top, "Coyote", "Winter", "Top")
  # wolf_s1819_occt <- occ_out(wolf_s1819_top, "Wolf", "Summer", "Top")
  # wolf_w1820_occt <- occ_out(wolf_w1820_top, "Wolf", "Winter", "Top")
  # elk_s1819_occt <- occ_out(elk_s1819_top, "Elk", "Summer", "Top")
  # elk_w1820_occt <- occ_out(elk_w1820_top, "Elk", "Winter", "Top")
  # md_s1819_occt <- occ_out(md_s1819_top, "Mule Deer", "Summer", "Top")
  # md_w1820_occt <- occ_out(md_w1820_top, "Mule Deer", "Winter", "Top")
  # wtd_s1819_occt <- occ_out(wtd_s1819_top, "White-tailed Deer", "Summer", "Top")
  # wtd_w1820_occt <- occ_out(wtd_w1820_top, "White-tailed Deer", "Winter", "Top")
  
  #'  Merge into larger data frames for easy comparison
  #'  Full models
  summer_occ <- rbind(bob_s1819_occ, coug_s1819_occ, coy_s1819_occ, wolf_s1819_occ,
                      elk_s1819_occ, md_s1819_occ, wtd_s1819_occ)
  winter_occ <- rbind(bob_w1820_occ, coug_w1820_occ, coy_w1820_occ, wolf_w1820_occ,
                      elk_w1820_occ, md_w1820_occ, wtd_w1820_occ)
  occ_results <- rbind(summer_occ, winter_occ) %>%
    arrange(Species)
  colnames(occ_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval") #"Model", 
  #' #'  Top models identified by derdging
  #' summer_occ_top <- rbind(bob_s1819_occt, coug_s1819_occt, coy_s1819_occt, wolf_s1819_occt,
  #'                     elk_s1819_occt, md_s1819_occt, wtd_s1819_occt)
  #' winter_occ_top <- rbind(bob_w1820_occt, coug_w1820_occt, coy_w1820_occt, wolf_w1820_occt,
  #'                     elk_w1820_occt, md_w1820_occt, wtd_w1820_occt)
  #' occ_results_top <- rbind(summer_occ_top, winter_occ_top) %>%
  #'   arrange(Species) 
  #' colnames(occ_results_top) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  #' #'  Combine full model and top model results
  #' summer_occ_combo <- rbind(bob_s1819_occ, bob_s1819_occt, coug_s1819_occ, coug_s1819_occt, 
  #'                         coy_s1819_occ, coy_s1819_occt, wolf_s1819_occ, wolf_s1819_occt,
  #'                         elk_s1819_occ, elk_s1819_occt, md_s1819_occ, md_s1819_occt, 
  #'                         wtd_s1819_occ, wtd_s1819_occt)
  #' winter_occ_combo <- rbind(bob_w1820_occ, bob_w1820_occt, coug_w1820_occ, coug_w1820_occt, 
  #'                         coy_w1820_occ, coy_w1820_occt, wolf_w1820_occ, wolf_w1820_occt,
  #'                         elk_w1820_occ, elk_w1820_occt, md_w1820_occ, md_w1820_occt, 
  #'                         wtd_w1820_occ, wtd_w1820_occt)
  #' occ_results_combo <- rbind(summer_occ_combo, winter_occ_combo) %>%
  #'   arrange(Species) 
  #' colnames(occ_results_combo) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are easier to look at
  rounddig <- 2
  results_psi <- occ_results %>%
    mutate(
      Estimate = round(Estimate, rounddig),
      SE = round(SE, rounddig),
      z = round(z, rounddig),
      Pval = round(Pval, rounddig)
    )
  # results_psi_top <- occ_results_top %>%
  #   mutate(
  #     Estimate = round(Estimate, rounddig),
  #     SE = round(SE, rounddig),
  #     z = round(z, rounddig),
  #     Pval = round(Pval, rounddig)
  #   )
  # results_psi_combo <- occ_results_combo %>%
  #   mutate(
  #     Estimate = round(Estimate, rounddig),
  #     SE = round(SE, rounddig),
  #     z = round(z, rounddig),
  #     Pval = round(Pval, rounddig)
  #   )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_psi_wide <- results_psi %>%  #results_psi_combo
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    separate("AreaOK", c("AreaOK (SE)", "AreaOK Pval"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev Pval"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix Pval"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass Pval"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub Pval"), sep = "_") %>%
    separate("RoadDensity", c("RoadDensity (SE)", "RoadDensity Pval"), sep = "_") %>%
    # separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))
  
  #'  Save!
  write.csv(results_psi, paste0("./Outputs/Tables/OccMod_OccProb_Results_matchRSF_", Sys.Date(), ".csv"))  
  write.csv(results_psi_wide, paste0("./Outputs/Tables/OccMod_OccProb_Results_wide_matchRSF_", Sys.Date(), ".csv"))
  
 
  #'  Function to save detection results
  det_out <- function(mod, spp, season, model) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
        # Model = rep(model, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) 
      # relocate(Model, .before = Species)
    return(out)
  }
  
  #'  Run each model through detection function
  bob_s1819_det <- det_out(bob_s1819_global, "Bobcat", "Summer") #, "Global"
  bob_w1820_det <- det_out(bob_w1820_global, "Bobcat", "Winter")
  coug_s1819_det <- det_out(coug_s1819_global, "Cougar", "Summer")
  coug_w1820_det <- det_out(coug_w1820_global, "Cougar", "Winter")
  coy_s1819_det <- det_out(coy_s1819_global, "Coyote", "Summer")
  coy_w1820_det <- det_out(coy_w1820_global, "Coyote", "Winter")
  wolf_s1819_det <- det_out(wolf_s1819_global2, "Wolf", "Summer")
  wolf_w1820_det <- det_out(wolf_w1820_global2, "Wolf", "Winter")
  elk_s1819_det <- det_out(elk_s1819_global2, "Elk", "Summer")
  elk_w1820_det <- det_out(elk_w1820_global2, "Elk", "Winter")
  md_s1819_det <- det_out(md_s1819_global, "Mule Deer", "Summer")
  md_w1820_det <- det_out(md_w1820_global, "Mule Deer", "Winter")
  wtd_s1819_det <- det_out(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_w1820_det <- det_out(wtd_w1820_global2, "White-tailed Deer", "Winter")

  # bob_s1819_dett <- det_out(bob_s1819_top, "Bobcat", "Summer", "Top")
  # bob_w1820_dett <- det_out(bob_w1820_top, "Bobcat", "Winter", "Top")
  # coug_s1819_dett <- det_out(coug_s1819_top, "Cougar", "Summer", "Top")
  # coug_w1820_dett <- det_out(coug_w1820_top, "Cougar", "Winter", "Top")
  # coy_s1819_dett <- det_out(coy_s1819_top, "Coyote", "Summer", "Top")
  # coy_w1820_dett <- det_out(coy_w1820_top, "Coyote", "Winter", "Top")
  # wolf_s1819_dett <- det_out(wolf_s1819_top, "Wolf", "Summer", "Top")
  # wolf_w1820_dett <- det_out(wolf_w1820_top, "Wolf", "Winter", "Top")
  # elk_s1819_dett <- det_out(elk_s1819_top, "Elk", "Summer", "Top")
  # elk_w1820_dett <- det_out(elk_w1820_top, "Elk", "Winter", "Top")
  # md_s1819_dett <- det_out(md_s1819_top, "Mule Deer", "Summer", "Top")
  # md_w1820_dett <- det_out(md_w1820_top, "Mule Deer", "Winter", "Top")
  # wtd_s1819_dett <- det_out(wtd_s1819_top, "White-tailed Deer", "Summer", "Top")
  # wtd_w1820_dett <- det_out(wtd_w1820_top, "White-tailed Deer", "Winter", "Top")
  
  #'  Merge into larger data frames for easy comparison
  summer_det <- rbind(bob_s1819_det, coug_s1819_det, coy_s1819_det, wolf_s1819_det,
                      elk_s1819_det, md_s1819_det, wtd_s1819_det)
  winter_det <- rbind(bob_w1820_det, coug_w1820_det, coy_w1820_det, wolf_w1820_det,
                      elk_w1820_det, md_w1820_det, wtd_w1820_det)
  det_results <- rbind(summer_det, winter_det) %>%
    arrange(Species)
  colnames(det_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval") #"Model", 
  
  # summer_det_top <- rbind(bob_s1819_dett, coug_s1819_dett, coy_s1819_dett, wolf_s1819_dett,
  #                     elk_s1819_dett, md_s1819_dett, wtd_s1819_dett)
  # winter_det_top <- rbind(bob_w1820_dett, coug_w1820_dett, coy_w1820_dett, wolf_w1820_dett,
  #                     elk_w1820_dett, md_w1820_dett, wtd_w1820_dett)
  # det_results_top <- rbind(summer_det_top, winter_det_top) %>%
  #   arrange(Species)
  # colnames(det_results_top) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  # 
  # summer_det_combo <- rbind(bob_s1819_det, bob_s1819_dett, coug_s1819_det, coug_s1819_dett, 
  #                           coy_s1819_det, coy_s1819_dett, wolf_s1819_det, wolf_s1819_dett,
  #                           elk_s1819_det, elk_s1819_dett, md_s1819_det, md_s1819_dett, 
  #                           wtd_s1819_det, wtd_s1819_dett)
  # winter_det_combo <- rbind(bob_w1820_det, bob_w1820_dett, coug_w1820_det, coug_w1820_dett, 
  #                           coy_w1820_det, coy_w1820_dett, wolf_w1820_det, wolf_w1820_dett,
  #                           elk_w1820_det, elk_w1820_dett, md_w1820_det, md_w1820_dett, 
  #                           wtd_w1820_det, wtd_w1820_dett)
  # det_results_combo <- rbind(summer_det_combo, winter_det_combo) %>%
  #   arrange(Species)
  # colnames(det_results_combo) <- c("Model", "Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are a little easier to interpret
  results_det <- det_results %>%
    mutate(
      Estimate = round(Estimate, 2),
      SE = round(SE, 2),
      z = round(z, 2),
      Pval = round(Pval, 2)
    )
  # results_det_top <- det_results_top %>%
  #   mutate(
  #     Estimate = round(Estimate, 3),
  #     SE = round(SE, 3),
  #     z = round(z, 3),
  #     Pval = round(Pval, 3),
  #     Parameter = ifelse(Parameter == "Distance:Height", "Height:Distance", Parameter)
  #   )
  # results_det_combo <- det_results_combo %>%
  #   mutate(
  #     Estimate = round(Estimate, 3),
  #     SE = round(SE, 3),
  #     z = round(z, 3),
  #     Pval = round(Pval, 3),
  #     Parameter = ifelse(Parameter == "Distance:Height", "Height:Distance", Parameter)
  #   )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_det_wide <- results_det %>% #results_det_combo
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    #' #'  Bold significant variables- doesn't work if continue manipulating data frame
    #' condformat(.) %>%
    #' rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    separate("TrailDirt road", c("Road (SE)", "Road Pval"), sep = "_") %>%
    separate("TrailDecommissioned road", c("Decom Road (SE)", "Decom Road Pval"), sep = "_") %>%
    separate("Height", c("Height (SE)", "Height Pval"), sep = "_") %>%
    separate("Temp_smr", c("Summer Temp (SE)", "Smr Temp Pval"), sep = "_") %>%
    separate("Temp_wtr", c("Winter Temp (SE)", "Wtr Temp Pval"), sep = "_") %>%
    separate("Distance", c("Distance (SE)", "Distance Pval"), sep = "_") %>%
    separate("Height:Distance", c("Height*Distance (SE)", "Height*Distance Pval"), sep = "_") %>%
    # separate("YearYear2", c("Year2 (SE)", "Year2 Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Wolf", "Coyote", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))

  #'  Save!
  write.csv(results_det, paste0("./Outputs/Tables/OccMod_DetProb_Results_matchRSF_", Sys.Date(), ".csv"))  
  write.csv(results_det_wide, paste0("./Outputs/Tables/OccMod_DetProb_Results_matchRSF_wide", Sys.Date(), ".csv"))


  #'  Predict probability of occupancy across sites
  mu_occ <- function(mod, species, season) {
    #'  Predict occupancy probability for all camera sties
    occu_mean <- predict(object = mod, type = "state") %>%  # do I provide newdata = sitecovariates if it's the same values as what went into the model?
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE)
    #'  Predict occupancy probability for all camera sties
    det_mean <- predict(object = mod, type = "det") %>%
      #'  Average occupancy probabilities across sites for mean psi
      summarise_at(c("Predicted", "SE"), mean, na.rm = TRUE) 
    predicted <- as.data.frame(rbind(occu_mean, det_mean))
    colnames(predicted) <- c("Mean", "SE")
    Parameter <- c("Occupancy", "Detection")
    Species <- species
    Season <- season
    predicted <- cbind(predicted, Parameter)
    predicted <- cbind(predicted, Species)
    predicted <- cbind(predicted, Season)
    return(predicted)
  }
  #'  Estimate mean probability of occupancy and detection per species and season
  md_predict_smr <- mu_occ(md_s1819_global, "Mule Deer", "Summer")
  md_predict_wtr <- mu_occ(md_w1820_global, "Mule Deer", "Winter")
  elk_predict_smr <- mu_occ(elk_s1819_global2, "Elk", "Summer")
  elk_predict_wtr <- mu_occ(elk_w1820_global2, "Elk", "Winter")
  wtd_predict_smr <- mu_occ(wtd_s1819_global2, "White-tailed Deer", "Summer")
  wtd_predict_wtr <- mu_occ(wtd_w1820_global2, "White-tailed Deer", "Winter")
  coug_predict_smr <- mu_occ(coug_s1819_global, "Cougar", "Summer")
  coug_predict_wtr <- mu_occ(coug_w1820_global, "Cougar", "Winter")
  wolf_predict_smr <- mu_occ(wolf_s1819_global2, "Wolf", "Summer")
  wolf_predict_wtr <- mu_occ(wolf_w1820_global2, "Wolf", "Winter")
  bob_predict_smr <- mu_occ(bob_s1819_global, "Bobcat", "Summer")
  bob_predict_wtr <- mu_occ(bob_w1820_global, "Bobcat", "Winter")
  coy_predict_smr <- mu_occ(coy_s1819_global, "Coyote", "Summer")
  coy_predict_wtr <- mu_occ(coy_w1820_global, "Coyote", "Winter")

  #'  Make a pretty table
  Mean_tbl <- bind_rows(md_predict_smr, md_predict_wtr, elk_predict_smr, elk_predict_wtr, wtd_predict_smr, 
              wtd_predict_wtr, coug_predict_smr, coug_predict_wtr, wolf_predict_smr, 
              wolf_predict_wtr, bob_predict_smr, bob_predict_wtr, coy_predict_smr, 
              coy_predict_wtr) %>%
    relocate(Species, .before = Mean) %>%
    relocate(Season, .after = Species) %>%
    relocate(Parameter, .after = Season) %>%
    arrange(Parameter, Mean, Species)
  
  #'  Save
  write.csv(Mean_tbl, paste0("./Outputs/Tables/OccMod_Mean_Estimates_matchRSF_", Sys.Date(), ".csv"))  

 
  #'  Save workspace
  save.image(file = "./Outputs/OccMod_script_results.RData")
  