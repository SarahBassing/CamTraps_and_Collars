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
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18-20_2021-03-24.csv") %>%
    #'  Get rid of mysterious space after one of the NEs
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    #'  Consolidate Monitoring feature (determined at camera site) into fewer 
    #'  categories that capture whether we expect vehicles to be on the feature
    # mutate(
    #   Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
    #   Monitoring = ifelse(Monitoring == "Trail", "Game trail", as.character(Monitoring)),
    #   Monitoring = ifelse(Monitoring == "Decommissioned road", "Game trail", as.character(Monitoring))
    # ) %>%
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", as.character(Monitoring)),
    ) %>%
    #'  Consolidate Habitat Type into fewer categories (determined at camera site)
    mutate(
      Habitat_Type = ifelse(Habitat_Type == "Agriculture", "Grassland", as.character(Habitat_Type)), # Grassland / Agriculture
      Habitat_Type = ifelse(Habitat_Type == "Riparian", "Mixed conifer", as.character(Habitat_Type)) # Mixed conifer / Riparian
    ) %>%
    #'  Consolidate NLCD and landcover into fewer categories (determined via remote sensing)
    #'  NOTE: Developed includes A LOT of different habitats and seems to be more 
    #'  tied with land ownership than landcover 
    #'  I don't trust NLCD is a good metric of landcover!
    mutate(
      NLCD_landcov = ifelse(NLCD_landcov == "21", "Other", as.character(NLCD_landcov)), # Developed open space seems like a catch all for a mix of habitats
      NLCD_landcov = ifelse(NLCD_landcov == "22", "Forested", as.character(NLCD_landcov)), # Developed low intensity is really forested areas
      NLCD_landcov = ifelse(NLCD_landcov == "42", "Forested", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "90", "Forested", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "95", "Forested", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "52", "Shrub", as.character(NLCD_landcov)),
      NLCD_landcov = ifelse(NLCD_landcov == "71", "Ag-Grassland", as.character(NLCD_landcov)), # Grassland / Agriculture
      NLCD_landcov = ifelse(NLCD_landcov == "81", "Ag-Grassland", as.character(NLCD_landcov)), 
      NLCD_landcov = ifelse(NLCD_landcov == "82", "Ag-Grassland", as.character(NLCD_landcov)),
      landcov18 = ifelse(landcov18 == "201", "ForestMix", as.character(landcov18)),
      landcov18 = ifelse(landcov18 == "230", "ForestMix", as.character(landcov18)),
      landcov18 = ifelse(landcov18 == "211", "MesicGrass", as.character(landcov18)),
      landcov18 = ifelse(landcov18 == "212", "XericGrass", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "310", "Developed", as.character(landcov18)), 
      landcov18 = ifelse(landcov18 == "332", "Developed", as.character(landcov18)),  
      landcov18 = ifelse(landcov18 == "221", "ForestMix", as.character(landcov18)),
      landcov18 = ifelse(landcov18 == "222", "XericShrub", as.character(landcov18)), 
      landcov19 = ifelse(landcov19 == "121", "MesicGrass", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "201", "ForestMix", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "230", "ForestMix", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "211", "MesicGrass", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "212", "XericGrass", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "310", "Developed", as.character(landcov19)), 
      landcov19 = ifelse(landcov19 == "332", "Developed", as.character(landcov19)),  
      landcov19 = ifelse(landcov19 == "221", "ForestMix", as.character(landcov19)),
      landcov19 = ifelse(landcov19 == "222", "XericShrub", as.character(landcov19)),  
      landcov = ifelse(landcov == "211", "MesicGrass", as.character(landcov)),
      landcov = ifelse(landcov == "121", "MesicGrass", as.character(landcov)), # Barren becomes mesic grass
      landcov = ifelse(landcov == "201", "ForestMix", as.character(landcov)),  # Emergent Wetland becomes forest
      landcov = ifelse(landcov == "212", "XericGrass", as.character(landcov)),
      landcov = ifelse(landcov == "310", "Developed", as.character(landcov)), # Agriculture becomes developed
      landcov = ifelse(landcov == "332", "Developed", as.character(landcov)), # Residential becomes developed 
      landcov = ifelse(landcov == "221", "ForestMix", as.character(landcov)), # Mesic shrub becomes forest  
      landcov = ifelse(landcov == "222", "XericShrub", as.character(landcov)),
      landcov = ifelse(landcov == "230", "ForestMix", as.character(landcov))
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
      NDVI_sp18 = scale(ndvi_sp18),             # KEEP IN MIND:
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
      PercForest = scale(PercForest),
      PercForestMix = scale(PercForestMix2),    # Forest + Mesic Shrub mix
      PercMesicShrub = scale(PercMesicShrub),
      PercMesicGrass = scale(PercMesicGrass),
      PercXericShrub = scale(PercXericShrub),
      PercXericGrass = scale(PercXericGrass),
      PercDeveloped = scale(PercDeveloped),
      Elev = scale(elevation),
      Slope = scale(slope), 
      Aspect = scale(aspect), # CIRCULAR! Also, 90-degrees is used when slope = 0
      TRI = scale(tri),
      Roughness = scale(roughness),
      Canopy18 = scale(canopy18),               # Courser-scale, 30m res
      Canopy19 = scale(canopy19),
      Canopy = scale(canopy),           # Skewed but transformations don't help
      Landfire = scale(landfire),       # Skewed but transformations don't help
      WaterDensity = scale(water_density), 
      RoadDensity = scale(road_density),   
      NearestH2o = scale(km2water),        
      NearestRd = scale(km2road),          
      HumanDensity = scale(human_density), 
      HumanModified = scale(modified)      
    )

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
  # write.csv(stations, paste0('./Outputs/stations18-20_', Sys.Date(), '.csv'))
    
  #'  Save study-area specific covaraites (important for ungulate models)
  stations_NE <- filter(stations, Study_Area == "NE")
  stations_OK <- filter(stations, Study_Area == "OK")
  
  #'  Check for correlation among covaraites
  #'  Watch out for NAs (use="complete.obs")
  cor(stations$Distance, stations$Height, use = "complete.obs")
  cor(stations$Elev, stations$Slope, use = "complete.obs")
  cor(stations$Elev, stations$Aspect, use = "complete.obs")
  cor(stations$Slope, stations$Aspect, use = "complete.obs")
  cor(stations$Elev, stations$TRI, use = "complete.obs")
  cor(stations$TRI, stations$Roughness, use = "complete.obs")  # EEK! 0.986
  cor(stations$Elev, stations$NDVI_sm18, use = "complete.obs")
  cor(stations$Slope, stations$NDVI_sm18, use = "complete.obs")
  cor(stations$Aspect, stations$NDVI_sm18, use = "complete.obs")
  cor(stations$TRI, stations$WaterDensity, use = "complete.obs")
  cor(stations$Elev, stations$WaterDensity, use = "complete.obs")
  cor(stations$Elev, stations$RoadDensity, use = "complete.obs")
  cor(stations$Elev, stations$HumanDensity, use = "complete.obs")
  cor(stations$Elev, stations$HumanModified, use = "complete.obs")  # YEP -0.639
  cor(stations$WaterDensity, stations$NearestH2o, use = "complete.obs")
  cor(stations$RoadDensity, stations$NearestRd, use = "complete.obs")
  cor(stations$RoadDensity, stations$HumanDensity, use = "complete.obs")
  cor(stations$WaterDensity, stations$HumanDensity, use = "complete.obs")
  cor(stations$Landfire, stations$Canopy, use = "complete.obs")             # YUP 0.668
  cor(stations$HumanDensity, stations$HumanModified, use = "complete.obs")  # EHH 0.498
  cor(stations$HumanModified, stations$RoadDensity, use = "complete.obs")
  cor(stations$PercForest, stations$NDVI_sm, use = "complete.obs")        # OOF 0.745
  cor(stations$PercForestMix, stations$NDVI_sm, use = "complete.obs")     # GAH 0.815
  cor(stations$PercMesicGrass, stations$NDVI_sm, use = "complete.obs")
  cor(stations$PercXericGrass, stations$NDVI_sm, use = "complete.obs")    # EHH -0.577
  cor(stations$PercMesicShrub, stations$NDVI_sm, use = "complete.obs")
  cor(stations$PercXericShrub, stations$NDVI_sm, use = "complete.obs")    # EHH -0.486
  cor(stations$PercMesicMix, stations$NDVI_sm, use = "complete.obs")
  cor(stations$PercDeveloped, stations$NDVI_sm, use = "complete.obs")
  cor(stations$PercForest, stations$PercMesicMix, use = "complete.obs")   # EHH -0.567
  cor(stations$PercForestMix, stations$Landfire, use = "complete.obs")    # YEAH 0.681
  
  #'  Study area specific correlations
  cor(stations_NE$Elev, stations_NE$HumanModified, use = "complete.obs")  # MEH -0.567
  cor(stations_OK$Elev, stations_OK$HumanModified, use = "complete.obs")  # YEP -0.688
  cor(stations_NE$RoadDensity, stations_NE$HumanModified, use = "complete.obs")
  cor(stations_OK$RoadDensity, stations_OK$HumanModified, use = "complete.obs") # MEH 0.573
  
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
                      nrow = nrows, ncol = ncols, byrow = FALSE)
    )
  # Effort_smr = matrix(Effort_smr1819, nrow = nrows, ncol = ncols, byrow = FALSE),
  # Effort_wtr = matrix(Effort_wtr1820, nrow = nrows, ncol = ncols, byrow = FALSE)
  #'  FYI effort covariate does weird things when scaled so not doing it now
  #'  Very few sampling occasions had low sampling effort (~3%; ignoring occasions
  #'  when camera failed completely) so not including sampling effort as a 
  #'  covariate on detection process- not enough variation to estimate effect
  
  #'  Double check it looks OK
  head(srvy_covs[[2]])
  
  #'  NEED TO BRING IN TEMPERATURE DATA AS WELL
  

  #'  Save study-area specific survey covariates
  Hgt_NE <- stations$Height[stations$Study_Area == "NE"]
  Hgt_OK <- stations$Height[stations$Study_Area == "OK"]
  Dist_NE <- stations$Distance[stations$Study_Area == "NE"]
  Dist_OK <- stations$Distance[stations$Study_Area == "OK"]
  
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
                      nrow = nrows_NE, ncol = ncols, byrow = FALSE)
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
                      nrow = nrows_OK, ncol = ncols, byrow = FALSE)
  )
  
  #'  Double check these
  head(srvy_covs[[1]])
  head(srvy_covs_NE[[1]])
  tail(srvy_covs[[1]])
  tail(srvy_covs_OK[[1]])
  
  
  ####  Naive Occupancy & Covariates  ####
  
  #'  Get a feel for naive occupancy across covaraites
  smbob <- as.numeric(apply(DH_bob_smr1819, 1, max, na.rm=FALSE))
  wtrbob <- as.numeric(apply(DH_bob_wtr1820, 1, max, na.rm=FALSE))
  smcoug <- as.numeric(apply(DH_coug_smr1819, 1, max, na.rm=FALSE))
  wtrcoug <- as.numeric(apply(DH_coug_wtr1820, 1, max, na.rm=FALSE))
  smcoy <- as.numeric(apply(DH_coy_smr1819, 1, max, na.rm=FALSE))
  wtrcoy <- as.numeric(apply(DH_coy_wtr1820, 1, max, na.rm=FALSE))
  smwolf <- as.numeric(apply(DH_wolf_smr1819, 1, max, na.rm=FALSE))
  wtrwolf <- as.numeric(apply(DH_wolf_wtr1820, 1, max, na.rm=FALSE))
  smelk <- as.numeric(apply(DH_elk_smr1819, 1, max, na.rm=FALSE))
  wtrelk <- as.numeric(apply(DH_elk_wtr1820, 1, max, na.rm=FALSE))
  smmd <- as.numeric(apply(DH_md_smr1819, 1, max, na.rm=FALSE))
  wtrmd <- as.numeric(apply(DH_md_wtr1820, 1, max, na.rm=FALSE))
  smwtd <- as.numeric(apply(DH_wtd_smr1819, 1, max, na.rm=FALSE))
  wtrwtd <- as.numeric(apply(DH_wtd_wtr1820, 1, max, na.rm=FALSE))
  
  #'  Raw detections by landcover categories
  landcov <- as.character(as.factor(stations$Landcov))
  landcovNE <- as.character(as.factor(stations_NE$Landcov))
  landcovOK <- as.character(as.factor(stations_OK$Landcov))
  hab_det <- as.data.frame(cbind(landcov, smbob, wtrbob, smcoug, wtrcoug, smcoy, wtrcoy, smwolf, wtrwolf)) %>%
    mutate(
      smbob = as.numeric(smbob),
      wtrbob = as.numeric(wtrbob),
      smcoug = as.numeric(smcoug),
      wtrcoug = as.numeric(wtrcoug),
      smcoy = as.numeric(smcoy),
      wtrcoy = as.numeric(wtrcoy),
      smwolf = as.numeric(smwolf),
      wtrwolf = as.numeric(wtrwolf)
    )
  habNE_det <- as.data.frame(cbind(landcovNE, smelk, wtrelk, smwtd, wtrwtd)) %>%
    mutate(
      smelk = as.numeric(smelk),
      wtrelk = as.numeric(wtrelk),
      smwtd = as.numeric(smwtd),
      wtrwtd = as.numeric(wtrwtd)
    )
  habOK_det <- as.data.frame(cbind(landcovOK, smmd, wtrmd)) %>%
    mutate(
      smmd = as.numeric(smmd),
      wtrmd = as.numeric(wtrmd)
    )
  #'  Count number of raw detections per landcover category
  tbl_pred <- hab_det %>%
    group_by(landcov) %>%
    summarise(across(smbob:wtrwolf, sum, na.rm= TRUE)) %>%
    ungroup()
  tbl_NE <- habNE_det %>%
    group_by(landcovNE) %>%
    summarise(across(smelk:wtrwtd, sum, na.rm= TRUE)) %>%
    ungroup()
  tbl_OK <- habOK_det %>%
    group_by(landcovOK) %>%
    summarise(across(smmd:wtrmd, sum, na.rm= TRUE)) %>%
    ungroup()
  sum_raw_det <- full_join(tbl_pred, tbl_OK, by = c("landcov" = "landcovOK")) %>%
    full_join(tbl_NE, by = c("landcov" = "landcovNE"))
  #'  Visualize another way
  barplot(sum_raw_det$smwtd, names.arg = sum_raw_det$landcov, 
          xlab = "Landcover Type", ylab = "Raw Detections")

  #'  Raw detections and continuous variables (remember- these are center & scaled)
  plot(stations$HumanModified, wtrcoy)
  plot(stations$PercForestMix, wtrcoy)
  hist(stations$PercForestMix[hab_det$wtrcoy == "1",], breaks = 50)
  hist(stations$PercXericShrub[habOK_det$wtrmd == "1",], breaks = 25)
  hist(stations$PercXericShrub[hab_det$wtrcoy == "1",], breaks = 25)

  
  
  #'  Create unmarked dataframes
  ####  BOBCAT UMF  ####
  bob_s1819_UMF <- unmarkedFrameOccu(DH_bob_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         NDVI = stations$NDVI_sm, 
                                                         Landcov = stations$Landcov,
                                                         PercFor = stations$PercForest,
                                                         PercForMix = stations$PercForestMix,
                                                         PercMGrass = stations$PercMesicGrass,
                                                         PercMShrub = stations$PercMesicShrub,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         PercDev = stations$PercDeveloped,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,                             
                                                         Landfire = stations$Landfire,
                                                         NearestH2o = stations$NearestH2o,
                                                         NearestRd = stations$NearestRd,
                                                         WaterDensity = stations$WaterDensity,
                                                         RoadDensity = stations$RoadDensity,
                                                         HumanDensity = stations$HumanDensity,
                                                         HumanMod = stations$HumanModified),
                                   obsCovs = srvy_covs)

  nrow(bob_s1819_UMF@y)
  #'  Remove rows with missing observation covariate data (Height & Distance data)
  # bob_s1819_UMF <- bob_s1819_UMF[-missing_dat]
  # nrow(bob_s1819_UMF@y)

  bob_w1820_UMF <- unmarkedFrameOccu(DH_bob_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Landcov = stations$Landcov,
                                                           PercFor = stations$PercForest,
                                                           PercForMix = stations$PercForestMix,
                                                           PercMGrass = stations$PercMesicGrass,
                                                           PercMShrub = stations$PercMesicShrub,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           PercDev = stations$PercDeveloped,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,                             
                                                           Landfire = stations$Landfire,
                                                           NearestH2o = stations$NearestH2o,
                                                           NearestRd = stations$NearestRd,
                                                           WaterDensity = stations$WaterDensity,
                                                           RoadDensity = stations$RoadDensity,
                                                           HumanDensity = stations$HumanDensity,
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
                                                          NDVI = stations$NDVI_sm,
                                                          Landcov = stations$Landcov,
                                                          PercFor = stations$PercForest,
                                                          PercForMix = stations$PercForestMix,
                                                          PercMGrass = stations$PercMesicGrass,
                                                          PercMShrub = stations$PercMesicShrub,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          PercDev = stations$PercDeveloped,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,                             
                                                          Landfire = stations$Landfire,
                                                          NearestH2o = stations$NearestH2o,
                                                          NearestRd = stations$NearestRd,
                                                          WaterDensity = stations$WaterDensity,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanDensity = stations$HumanDensity,
                                                          HumanMod = stations$HumanModified),
                                    obsCovs = srvy_covs)
  nrow(coug_s1819_UMF@y)
  # coug_s1819_UMF <- coug_s1819_UMF[-missing_dat]
  # nrow(coug_s1819_UMF@y)
  
  coug_w1820_UMF <- unmarkedFrameOccu(DH_coug_wtr1820,
                                    siteCovs = data.frame(Year = stations$Year,
                                                          Area = stations$Study_Area,
                                                          Trail = stations$Trail,
                                                          NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                          Landcov = stations$Landcov,
                                                          PercFor = stations$PercForest,
                                                          PercForMix = stations$PercForestMix,
                                                          PercMGrass = stations$PercMesicGrass,
                                                          PercMShrub = stations$PercMesicShrub,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          PercDev = stations$PercDeveloped,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,                             
                                                          Landfire = stations$Landfire,
                                                          NearestH2o = stations$NearestH2o,
                                                          NearestRd = stations$NearestRd,
                                                          WaterDensity = stations$WaterDensity,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanDensity = stations$HumanDensity,
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
                                                          NDVI = stations$NDVI_sm, 
                                                          Landcov = stations$Landcov,
                                                          PercFor = stations$PercForest,
                                                          PercForMix = stations$PercForestMix,
                                                          PercMGrass = stations$PercMesicGrass,
                                                          PercMShrub = stations$PercMesicShrub,
                                                          PercXGrass = stations$PercXericGrass,
                                                          PercXShrub = stations$PercXericShrub,
                                                          PercDev = stations$PercDeveloped,
                                                          Elev = stations$Elev,
                                                          Slope = stations$Slope,
                                                          Aspect = stations$Aspect,                             
                                                          Landfire = stations$Landfire,
                                                          NearestH2o = stations$NearestH2o,
                                                          NearestRd = stations$NearestRd,
                                                          WaterDensity = stations$WaterDensity,
                                                          RoadDensity = stations$RoadDensity,
                                                          HumanDensity = stations$HumanDensity,
                                                          HumanMod = stations$HumanModified),
                                    obsCovs = srvy_covs)
  # coy_s1819_UMF <- coy_s1819_UMF[-missing_dat]
  
  coy_w1820_UMF <- unmarkedFrameOccu(DH_coy_wtr1820,
                                      siteCovs = data.frame(Year = stations$Year,
                                                            Area = stations$Study_Area,
                                                            Trail = stations$Trail,
                                                            NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                            Landcov = stations$Landcov,
                                                            PercFor = stations$PercForest,
                                                            PercForMix = stations$PercForestMix,
                                                            PercMGrass = stations$PercMesicGrass,
                                                            PercMShrub = stations$PercMesicShrub,
                                                            PercXGrass = stations$PercXericGrass,
                                                            PercXShrub = stations$PercXericShrub,
                                                            PercDev = stations$PercDeveloped,
                                                            Elev = stations$Elev,
                                                            Slope = stations$Slope,
                                                            Aspect = stations$Aspect,                             
                                                            Landfire = stations$Landfire,
                                                            NearestH2o = stations$NearestH2o,
                                                            NearestRd = stations$NearestRd,
                                                            WaterDensity = stations$WaterDensity,
                                                            RoadDensity = stations$RoadDensity,
                                                            HumanDensity = stations$HumanDensity,
                                                            HumanMod = stations$HumanModified),
                                      obsCovs = srvy_covs)
  # coy_w1820_UMF <- coy_w1820_UMF[-missing_dat]
  
  ####  WOLF UMF  ####
  wolf_s1819_UMF <- unmarkedFrameOccu(DH_wolf_smr1819,
                                   siteCovs = data.frame(Year = stations$Year,
                                                         Area = stations$Study_Area,
                                                         Trail = stations$Trail,
                                                         NDVI = stations$NDVI_sm, 
                                                         Landcov = stations$Landcov,
                                                         PercFor = stations$PercForest,
                                                         PercForMix = stations$PercForestMix,
                                                         PercMGrass = stations$PercMesicGrass,
                                                         PercMShrub = stations$PercMesicShrub,
                                                         PercXGrass = stations$PercXericGrass,
                                                         PercXShrub = stations$PercXericShrub,
                                                         PercDev = stations$PercDeveloped,
                                                         Elev = stations$Elev,
                                                         Slope = stations$Slope,
                                                         Aspect = stations$Aspect,                             
                                                         Landfire = stations$Landfire,
                                                         NearestH2o = stations$NearestH2o,
                                                         NearestRd = stations$NearestRd,
                                                         WaterDensity = stations$WaterDensity,
                                                         RoadDensity = stations$RoadDensity,
                                                         HumanDensity = stations$HumanDensity,
                                                         HumanMod = stations$HumanModified),
                                   obsCovs = srvy_covs)
  # wolf_s1819_UMF <- wolf_s1819_UMF[-missing_dat]
  
  wolf_w1820_UMF <- unmarkedFrameOccu(DH_wolf_wtr1820,
                                     siteCovs = data.frame(Year = stations$Year,
                                                           Area = stations$Study_Area,
                                                           Trail = stations$Trail,
                                                           NDVI = stations$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Landcov = stations$Landcov,
                                                           PercFor = stations$PercForest,
                                                           PercForMix = stations$PercForestMix,
                                                           PercMGrass = stations$PercMesicGrass,
                                                           PercMShrub = stations$PercMesicShrub,
                                                           PercXGrass = stations$PercXericGrass,
                                                           PercXShrub = stations$PercXericShrub,
                                                           PercDev = stations$PercDeveloped,
                                                           Elev = stations$Elev,
                                                           Slope = stations$Slope,
                                                           Aspect = stations$Aspect,                             
                                                           Landfire = stations$Landfire,
                                                           NearestH2o = stations$NearestH2o,
                                                           NearestRd = stations$NearestRd,
                                                           WaterDensity = stations$WaterDensity,
                                                           RoadDensity = stations$RoadDensity,
                                                           HumanDensity = stations$HumanDensity,
                                                           HumanMod = stations$HumanModified),
                                     obsCovs = srvy_covs)
  # wolf_w1820_UMF <- wolf_w1820_UMF[-missing_dat]
  
  
  ####  ELK UMF  ####
  #'  Data from NE study area only to mirror GPS collar data
  elk_s1819_UMF <- unmarkedFrameOccu(DH_elk_smr1819,
                                    siteCovs = data.frame(Year = stations_NE$Year,
                                                          Trail = stations_NE$Trail,
                                                          NDVI = stations_NE$NDVI_sm,
                                                          Landcov = stations_NE$Landcov,
                                                          PercFor = stations_NE$PercForest,
                                                          PercForMix = stations_NE$PercForestMix,
                                                          PercMGrass = stations_NE$PercMesicGrass,
                                                          PercMShrub = stations_NE$PercMesicShrub,
                                                          PercXGrass = stations_NE$PercXericGrass,
                                                          PercXShrub = stations_NE$PercXericShrub,
                                                          PercDev = stations_NE$PercDeveloped,
                                                          Elev = stations_NE$Elev,
                                                          Slope = stations_NE$Slope,
                                                          Aspect = stations_NE$Aspect,
                                                          Landfire = stations_NE$Landfire,
                                                          NearestH2o = stations_NE$NearestH2o,
                                                          NearestRd = stations_NE$NearestRd,
                                                          WaterDensity = stations_NE$WaterDensity,
                                                          RoadDensity = stations_NE$RoadDensity,
                                                          HumanDensity = stations_NE$HumanDensity,
                                                          HumanMod = stations_NE$HumanModified),
                                    obsCovs = srvy_covs_NE)
  # elk_s1819_UMF <- elk_s1819_UMF[-missing_dat]
  
  elk_w1820_UMF <- unmarkedFrameOccu(DH_elk_wtr1820,
                                      siteCovs = data.frame(Year = stations_NE$Year,
                                                            Trail = stations_NE$Trail,
                                                            NDVI = stations_NE$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                            Landcov = stations_NE$Landcov,
                                                            PercFor = stations_NE$PercForest,
                                                            PercForMix = stations_NE$PercForestMix,
                                                            PercMGrass = stations_NE$PercMesicGrass,
                                                            PercMShrub = stations_NE$PercMesicShrub,
                                                            PercXGrass = stations_NE$PercXericGrass,
                                                            PercXShrub = stations_NE$PercXericShrub,
                                                            PercDev = stations_NE$PercDeveloped,
                                                            Elev = stations_NE$Elev,
                                                            Slope = stations_NE$Slope,
                                                            Aspect = stations_NE$Aspect,
                                                            Landfire = stations_NE$Landfire,
                                                            NearestH2o = stations_NE$NearestH2o,
                                                            NearestRd = stations_NE$NearestRd,
                                                            WaterDensity = stations_NE$WaterDensity,
                                                            RoadDensity = stations_NE$RoadDensity,
                                                            HumanDensity = stations_NE$HumanDensity,
                                                            HumanMod = stations_NE$HumanModified),
                                      obsCovs = srvy_covs_NE)
  # elk_w1820_UMF <- elk_w1820_UMF[-missing_dat]
  
  ####  MULE DEER UMF  ####
  #'  Data from OK study area only to mirror GPS collar data
  md_s1819_UMF <- unmarkedFrameOccu(DH_md_smr1819,
                                   siteCovs = data.frame(Year = stations_OK$Year,
                                                         Trail = stations_OK$Trail,
                                                         NDVI = stations_OK$NDVI_sm, 
                                                         Landcov = stations_OK$Landcov,
                                                         PercFor = stations_OK$PercForest,
                                                         PercForMix = stations_OK$PercForestMix,
                                                         PercMGrass = stations_OK$PercMesicGrass,
                                                         PercMShrub = stations_OK$PercMesicShrub,
                                                         PercXGrass = stations_OK$PercXericGrass,
                                                         PercXShrub = stations_OK$PercXericShrub,
                                                         PercDev = stations_OK$PercDeveloped,
                                                         Elev = stations_OK$Elev,
                                                         Slope = stations_OK$Slope,
                                                         Aspect = stations_OK$Aspect,
                                                         Landfire = stations_OK$Landfire,
                                                         NearestH2o = stations_OK$NearestH2o,
                                                         NearestRd = stations_OK$NearestRd,
                                                         WaterDensity = stations_OK$WaterDensity,
                                                         RoadDensity = stations_OK$RoadDensity,
                                                         HumanDensity = stations_OK$HumanDensity,
                                                         HumanMod = stations_OK$HumanModified),
                                   obsCovs = srvy_covs_OK)
  # md_s1819_UMF <- md_s1819_UMF[-missing_dat]
  
  md_w1820_UMF <- unmarkedFrameOccu(DH_md_wtr1820,
                                     siteCovs = data.frame(Year = stations_OK$Year,
                                                           Trail = stations_OK$Trail,
                                                           NDVI = stations_OK$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                           Landcov = stations_OK$Landcov,
                                                           PercFor = stations_OK$PercForest,
                                                           PercForMix = stations_OK$PercForestMix,
                                                           PercMGrass = stations_OK$PercMesicGrass,
                                                           PercMShrub = stations_OK$PercMesicShrub,
                                                           PercXGrass = stations_OK$PercXericGrass,
                                                           PercXShrub = stations_OK$PercXericShrub,
                                                           PercDev = stations_OK$PercDeveloped,
                                                           Elev = stations_OK$Elev,
                                                           Slope = stations_OK$Slope,
                                                           Aspect = stations_OK$Aspect,
                                                           Landfire = stations_OK$Landfire,
                                                           NearestH2o = stations_OK$NearestH2o,
                                                           NearestRd = stations_OK$NearestRd,
                                                           WaterDensity = stations_OK$WaterDensity,
                                                           RoadDensity = stations_OK$RoadDensity,
                                                           HumanDensity = stations_OK$HumanDensity,
                                                           HumanMod = stations_OK$HumanModified),
                                     obsCovs = srvy_covs_OK)
  # md_w1820_UMF <- md_w1820_UMF[-missing_dat]
  
  ####  WHITE-TAILED DEER UMF  ####
  #'  Data from NE study area only to mirror GPS collar data
  wtd_s1819_UMF <- unmarkedFrameOccu(DH_wtd_smr1819,
                                  siteCovs = data.frame(Year = stations_NE$Year,
                                                        Trail = stations_NE$Trail,
                                                        NDVI = stations_NE$NDVI_sm,
                                                        Landcov = stations_NE$Landcov,
                                                        PercFor = stations_NE$PercForest,
                                                        PercForMix = stations_NE$PercForestMix,
                                                        PercMGrass = stations_NE$PercMesicGrass,
                                                        PercMShrub = stations_NE$PercMesicShrub,
                                                        PercXGrass = stations_NE$PercXericGrass,
                                                        PercXShrub = stations_NE$PercXericShrub,
                                                        PercDev = stations_NE$PercDeveloped,
                                                        Elev = stations_NE$Elev,
                                                        Slope = stations_NE$Slope,
                                                        Aspect = stations_NE$Aspect,
                                                        Landfire = stations_NE$Landfire,
                                                        NearestH2o = stations_NE$NearestH2o,
                                                        NearestRd = stations_NE$NearestRd,
                                                        WaterDensity = stations_NE$WaterDensity,
                                                        RoadDensity = stations_NE$RoadDensity,
                                                        HumanDensity = stations_NE$HumanDensity,
                                                        HumanMod = stations_NE$HumanModified),
                                  obsCovs = srvy_covs_NE)
  # wtd_s1819_UMF <- wtd_s1819_UMF[-missing_dat]
  
  wtd_w1820_UMF <- unmarkedFrameOccu(DH_wtd_wtr1820,
                                    siteCovs = data.frame(Year = stations_NE$Year,
                                                          Trail = stations_NE$Trail,
                                                          NDVI = stations_NE$NDVI_sm, #  USE SUMMER NDVI FOR WINTER MODELS
                                                          Landcov = stations_NE$Landcov,
                                                          PercFor = stations_NE$PercForest,
                                                          PercForMix = stations_NE$PercForestMix,
                                                          PercMGrass = stations_NE$PercMesicGrass,
                                                          PercMShrub = stations_NE$PercMesicShrub,
                                                          PercXGrass = stations_NE$PercXericGrass,
                                                          PercXShrub = stations_NE$PercXericShrub,
                                                          PercDev = stations_NE$PercDeveloped,
                                                          Elev = stations_NE$Elev,
                                                          Slope = stations_NE$Slope,
                                                          Aspect = stations_NE$Aspect,
                                                          Landfire = stations_NE$Landfire,
                                                          NearestH2o = stations_NE$NearestH2o,
                                                          NearestRd = stations_NE$NearestRd,
                                                          WaterDensity = stations_NE$WaterDensity,
                                                          RoadDensity = stations_NE$RoadDensity,
                                                          HumanDensity = stations_NE$HumanDensity,
                                                          HumanMod = stations_NE$HumanModified),
                                    obsCovs = srvy_covs_NE)
  # wtd_w1820_UMF <- wtd_w1820_UMF[-missing_dat]
  nrow( wtd_w1820_UMF@y)
  
  
  #'  Occupancy models
  #'  =============================
  #'  unmarked formula: ~detection ~occupancy
  #'  ~1 for intercept only
  #'  
  #'  Use chi-sq test to evaluate model fit after model selection (pg. 4 vignette)


  ####  BOBCAT MODELS  ####                   
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019
  #'  Null
  (bob_s1819_null <- occu(~1 ~1, bob_s1819_UMF))
  backTransform(bob_s1819_null, 'det')
  backTransform(bob_s1819_null, 'state')
  #'  Null with detection covariates
  (bob_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, bob_s1819_UMF))
  #'  Terrain model
  (bob_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, bob_s1819_UMF))
  (bob_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, bob_s1819_UMF))
  #'  Vegetation model
  (bob_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, bob_s1819_UMF))
  #'  Habitat model
  (bob_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, bob_s1819_UMF))
  (bob_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, bob_s1819_UMF))
  #'  Anthropogenic model 
  (bob_s1819_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, bob_s1819_UMF))
  #'  Combined model
  (bob_s1819_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, bob_s1819_UMF))
  (bob_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, bob_s1819_UMF))
  
  mods <- fitList(bob_s1819_null, bob_s1819_null2, bob_s1819_terrain, bob_s1819_terrain2, 
                  bob_s1819_veg, bob_s1819_anthro, bob_s1819_terrain_veg, bob_s1819_terrain_veg2, 
                  bob_s1819_full, bob_s1819_full2)
  modSel(mods)
  
  
  #'  WINTERS 2018-2019 & 2019-2020           
  #'  Null
  (bob_w1820_null <- occu(~1 ~1, bob_w1820_UMF))
  backTransform(bob_w1820_null, 'det')
  backTransform(bob_w1820_null, 'state')
  #'  Null with detection covariates
  (bob_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, bob_w1820_UMF))
  #'  Terrain model
  (bob_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, bob_w1820_UMF))
  (bob_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, bob_w1820_UMF))
  #'  Vegetation model
  (bob_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, bob_w1820_UMF))  
  #'  Habitat model
  (bob_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, bob_w1820_UMF))
  (bob_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, bob_w1820_UMF))
  #'  Anthropogenic model 
  (bob_w1820_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, bob_w1820_UMF))
  #'  Combined model
  (bob_w1820_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, bob_w1820_UMF))
  (bob_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, bob_w1820_UMF))
  
  mods <- fitList(bob_w1820_null, bob_w1820_null2, bob_w1820_terrain, bob_w1820_terrain2, 
                  bob_w1820_veg, bob_w1820_anthro, bob_w1820_terrain_veg, bob_w1820_terrain_veg2,
                  bob_w1820_full, bob_w1820_full2)
  modSel(mods)
  
  ####  COUGAR MODELS  ####
  #'  Included covariates informed by Dickson & Beier 2002, Kertson et al. 2011, 
  #'  Smereka et al. 2020
                                              
  #'  SUMMERS 2018 & 2019                     
  #'  Null
  (coug_s1819_null <- occu(~1 ~1, coug_s1819_UMF))
  backTransform(coug_s1819_null, 'det')
  backTransform(coug_s1819_null, 'state')
  #'  Null with detection covariates
  (coug_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, coug_s1819_UMF))
  #'  Terrain model
  (coug_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, coug_s1819_UMF))
  (coug_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, coug_s1819_UMF))
  #'  Vegetation model
  (coug_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, coug_s1819_UMF))
  #'  Habitat model
  (coug_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, coug_s1819_UMF))
  (coug_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, coug_s1819_UMF))
  #'  Anthropogenic model
  (coug_s1819_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, coug_s1819_UMF))
  #'  Combined model
  (coug_s1819_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coug_s1819_UMF))
  (coug_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coug_s1819_UMF))

  mods <- fitList(coug_s1819_null, coug_s1819_null2, coug_s1819_terrain, coug_s1819_terrain2, 
                  coug_s1819_veg, coug_s1819_anthro, coug_s1819_terrain_veg, coug_s1819_terrain_veg2,
                  coug_s1819_full, coug_s1819_full2) 
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Null
  (coug_w1820_null <- occu(~1 ~1, coug_w1820_UMF))
  backTransform(coug_w1820_null, 'det')
  backTransform(coug_w1820_null, 'state')
  #'  Null with detection covariates
  (coug_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, coug_w1820_UMF))
  #'  Terrain model
  (coug_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, coug_w1820_UMF))
  (coug_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, coug_w1820_UMF))
  #'  Vegetation model
  (coug_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, coug_w1820_UMF))
  #'  Habitat model
  (coug_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, coug_w1820_UMF))
  (coug_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, coug_w1820_UMF))
  #'  Anthropogenic model
  (coug_w1820_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, coug_w1820_UMF))
  #'  Combined model
  (coug_w1820_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coug_w1820_UMF))
  (coug_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coug_w1820_UMF))
  
  mods <- fitList(coug_w1820_null, coug_w1820_null2, coug_w1820_terrain, coug_w1820_terrain2, 
                  coug_w1820_veg, coug_w1820_anthro, coug_w1820_terrain_veg, coug_w1820_terrain_veg2,
                  coug_w1820_full, coug_w1820_full2)
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
  (coy_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, coy_s1819_UMF))
  (coy_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, coy_s1819_UMF))
  #'  Vegetation model
  (coy_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, coy_s1819_UMF))
  #'  Habitat model
  (coy_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, coy_s1819_UMF))
  (coy_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, coy_s1819_UMF))
  #'  Anthropogenic model
  (coy_s1819_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, coy_s1819_UMF))
  #'  Combined model
  (coy_s1819_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coy_s1819_UMF))
  (coy_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coy_s1819_UMF))
  
  mods <- fitList(coy_s1819_null, coy_s1819_null2, coy_s1819_terrain, coy_s1819_terrain2, 
                  coy_s1819_veg, coy_s1819_anthro, coy_s1819_terrain_veg, coy_s1819_terrain_veg2, 
                  coy_s1819_full, coy_s1819_full2)
  modSel(mods)
  
  #'  WINTERS 2018-2019 & 2019-2020              
  #'  Null
  (coy_w1820_null <- occu(~1 ~1, coy_w1820_UMF))
  backTransform(coy_w1820_null, 'det')
  backTransform(coy_w1820_null, 'state')
  #'  Null with detection covariates
  (coy_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, coy_w1820_UMF))
  #'  Terrain model
  (coy_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, coy_w1820_UMF))
  (coy_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, coy_w1820_UMF))
  #'  Vegetation model
  (coy_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, coy_w1820_UMF))
  #'  Habitat model
  (coy_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, coy_w1820_UMF))
  (coy_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, coy_w1820_UMF))
  #'  Anthropogenic model
  (coy_w1820_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, coy_w1820_UMF))
  #'  Combined model
  (coy_w1820_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coy_w1820_UMF))
  (coy_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, coy_w1820_UMF))
  
  mods <- fitList(coy_w1820_null, coy_w1820_null2, coy_w1820_terrain, coy_w1820_terrain2, 
                  coy_w1820_veg, coy_w1820_anthro, coy_w1820_terrain_veg, coy_w1820_terrain_veg2,
                  coy_w1820_full, coy_w1820_full2)
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
  (wolf_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, wolf_s1819_UMF))
  (wolf_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, wolf_s1819_UMF))
  #'  Vegetation model------------------ NOTE: Dropped Xeric Shrub from model due to lack of convergence (was messing with intercept too)
  (wolf_s1819_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass, wolf_s1819_UMF))
  #'  Habitat model
  (wolf_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass, wolf_s1819_UMF))
  (wolf_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass, wolf_s1819_UMF))
  #'  Anthropogenic model 
  (wolf_s1819_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, wolf_s1819_UMF))
  #'  Combined model
  (wolf_s1819_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wolf_s1819_UMF))
  (wolf_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wolf_s1819_UMF))
  
  mods <- fitList(wolf_s1819_null, wolf_s1819_null2, wolf_s1819_terrain, wolf_s1819_terrain2, 
                  wolf_s1819_veg, wolf_s1819_anthro, wolf_s1819_terrain_veg, 
                  wolf_s1819_terrain_veg2, wolf_s1819_full, wolf_s1819_full2)
  modSel(mods)
  #  Naive occupancy low (0.139) so may limit estimable effect of most covariates
  
  #'  WINTERS 2018-2019 & 2019-2020            
  #'  Null
  (wolf_w1820_null <- occu(~1 ~1, wolf_w1820_UMF))
  backTransform(wolf_w1820_null, 'det')
  backTransform(wolf_w1820_null, 'state')
  #'  Null with detection covariates
  (wolf_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, wolf_w1820_UMF))
  #'  Terrain model
  (wolf_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope, wolf_w1820_UMF))
  (wolf_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope, wolf_w1820_UMF))
  #'  Vegetation model
  (wolf_w1820_veg <- occu(~Height*Distance + Trail + Year ~Area + PercForMix + PercXGrass + PercXShrub, wolf_w1820_UMF))
  #'  Habitat model
  (wolf_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + PercXShrub, wolf_w1820_UMF))  
  (wolf_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, wolf_w1820_UMF))  
  #'  Anthropogenic model
  (wolf_w1820_anthro <- occu(~Height*Distance + Trail + Year ~Area + RoadDensity + HumanMod, wolf_w1820_UMF))
  #'  Combined model
  (wolf_w1820_full <- occu(~Height*Distance + Trail + Year ~Area + Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wolf_w1820_UMF))  
  (wolf_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Area + Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wolf_w1820_UMF))  #  FAILS to converge
  
  mods <- fitList(wolf_w1820_null, wolf_w1820_null2, wolf_w1820_terrain,  
                  wolf_w1820_veg, wolf_w1820_anthro, wolf_w1820_terrain_veg) 
                  #wolf_w1820_terrain2, #wolf_w1820_terrain_veg2, #wolf_w1820_full, #wolf_w1820_full2 FAILS
  modSel(mods)
  #'  Naive occupancy low (0.122) so may limit estimable effect of most covariates
  #'  Anthro model lowest AIC b/c includes Area covariate- has nothing to do with anthro effects
  
  ####  ELK MODELS ####                          
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019, NE study area only so no Area effect                             
  #'  Null
  (elk_s1819_null <- occu(~1 ~1, elk_s1819_UMF))
  backTransform(elk_s1819_null, 'det')
  backTransform(elk_s1819_null, 'state')
  #'  Null with detection covariates
  (elk_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, elk_s1819_UMF))
  #'  Terrain model
  (elk_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, elk_s1819_UMF))
  (elk_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, elk_s1819_UMF))
  #'  Vegetation model--------- NOTE: removed Xeric Shrub from veg models b/c very little Xeric Shrub in NE study area
  (elk_s1819_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass , elk_s1819_UMF))
  #'  Habitat model
  (elk_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass, elk_s1819_UMF))
  (elk_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass, elk_s1819_UMF))
  #'  Anthropogenic model 
  (elk_s1819_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, elk_s1819_UMF))
  #'  Combined model
  (elk_s1819_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, elk_s1819_UMF))
  (elk_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, elk_s1819_UMF))
  
  mods <- fitList(elk_s1819_null, elk_s1819_null2, elk_s1819_terrain, elk_s1819_terrain2, 
                  elk_s1819_veg, elk_s1819_terrain_veg, elk_s1819_terrain_veg2, 
                  elk_s1819_anthro, elk_s1819_full, elk_s1819_full2)
  modSel(mods)
  # Naive occupancy isn't that low (0.409) so not sure why no psi covariates are significant
  
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only so no Area effect                
  #'  Null                                         
  (elk_w1820_null <- occu(~1 ~1, elk_w1820_UMF))
  backTransform(elk_w1820_null, 'det')
  backTransform(elk_w1820_null, 'state')
  #'  Null with detection covariates
  (elk_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, elk_w1820_UMF))
  #'  Terrain model
  (elk_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, elk_w1820_UMF))
  (elk_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, elk_w1820_UMF))
  #'  Vegetation model--------- NOTE: removed Xeric Shrub from veg models b/c very little Xeric Shrub in NE study area
  (elk_w1820_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass, elk_w1820_UMF))
  #'  Habitat model
  (elk_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass, elk_w1820_UMF))
  (elk_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass, elk_w1820_UMF))
  #'  Anthropogenic model
  (elk_w1820_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, elk_w1820_UMF))
  #'  Combined model
  (elk_w1820_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, elk_w1820_UMF))
  (elk_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, elk_w1820_UMF))
  
  mods <- fitList(elk_w1820_null, elk_w1820_null2, elk_w1820_terrain, elk_w1820_terrain2, 
                  elk_w1820_veg, elk_w1820_terrain_veg, elk_w1820_terrain_veg2, 
                  elk_w1820_anthro, elk_w1820_full, elk_w1820_full2) 
                  #elk_w1820_veg, #elk_w1820_terrain_veg, #elk_w1820_terrain_veg2, #, elk_w1820_full, #elk_w1820_full2 questionable convergence on Grass variable
  modSel(mods)
  #  Naive occupancy pretty low (0.159) which may explain why psi covaraites not significant
  
  
  ####  MULE DEER MODELS  ####
  #'  Included covariates informed by knowledge of study system & 
  
  #'  SUMMERS 2018 & 2019, OK study area only so no Area effect
  #'  Null
  (md_s1819_null <- occu(~1 ~1, md_s1819_UMF))
  backTransform(md_s1819_null, 'det')
  backTransform(md_s1819_null, 'state')
  #'  Null with detection covariates
  (md_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, md_s1819_UMF))
  #'  Terrain model
  (md_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, md_s1819_UMF))
  (md_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, md_s1819_UMF))
  #'  Vegetation model
  (md_s1819_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass + PercXShrub, md_s1819_UMF))
  #'  Habitat model
  (md_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub, md_s1819_UMF))
  (md_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, md_s1819_UMF))
  #'  Anthropogenic model
  (md_s1819_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, md_s1819_UMF))
  #'  Combined model
  (md_s1819_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, md_s1819_UMF))
  (md_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, md_s1819_UMF))
  
  mods <- fitList(md_s1819_null, md_s1819_null2, md_s1819_terrain, md_s1819_terrain2, 
                  md_s1819_veg, md_s1819_terrain_veg, md_s1819_terrain_veg2, 
                  md_s1819_anthro, md_s1819_full, md_s1819_full2)
  modSel(mods)
  #'  Naive occupancy is high (0.84) but can still estimate at least some covariate effects on psi
  #'  Given elev & hum mod are highly correlated in the OK I suspect -effect of 
  #'  hum mod in anthro model is more associated with human-elevation relationship
  #'  than actual avoidance of humans... then again, why wouldn't there be a similarly
  #'  significant -effect of elevation in th terrain model is that were the case?
  
  #'  WINTERS 2018-2019 & 2019-2020, OK study area only so no Area effect                    
  #'  Null
  (md_w1820_null <- occu(~1 ~1, md_w1820_UMF))
  backTransform(md_w1820_null, 'det')
  backTransform(md_w1820_null, 'state')
  #'  Null with detection covariates
  (md_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, md_w1820_UMF))
  #'  Terrain model
  (md_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, md_w1820_UMF))
  (md_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, md_w1820_UMF))
  #'  Vegetation model
  (md_w1820_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass + PercXShrub, md_w1820_UMF))
  #'  Habitat model
  (md_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub, md_w1820_UMF))
  (md_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub, md_w1820_UMF))
  #'  Anthropogenic model
  (md_w1820_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, md_w1820_UMF))
  #'  Combined model
  (md_w1820_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, md_w1820_UMF))
  (md_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + PercXShrub + RoadDensity + HumanMod, md_w1820_UMF))
  
  mods <- fitList(md_w1820_null, md_w1820_null2, md_w1820_terrain, md_w1820_terrain2, 
                  md_w1820_veg, md_w1820_terrain_veg, md_w1820_terrain_veg2, 
                  md_w1820_anthro, md_w1820_full, md_w1820_full2)
  modSel(mods)
  
  
  ####  WHITE-TAILED DEER MODELS  ####
  #'  Included covariates informed by 
  
  #'  SUMMERS 2018 & 2019, NE study area only so no Area effect 
  #'  Null
  (wtd_s1819_null <- occu(~1 ~1, wtd_s1819_UMF))
  backTransform(wtd_s1819_null, 'det')
  backTransform(wtd_s1819_null, 'state')
  #'  Null with detection covariates
  (wtd_s1819_null2 <- occu(~Height*Distance + Trail + Year ~1, wtd_s1819_UMF))
  #'  Terrain model
  (wtd_s1819_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, wtd_s1819_UMF))
  (wtd_s1819_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, wtd_s1819_UMF))
  #'  Vegetation model--------- NOTE: removed Xeric Shrub from veg models b/c very little Xeric Shrub in NE study area
  (wtd_s1819_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass, wtd_s1819_UMF))
  #'  Habitat model
  (wtd_s1819_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass, wtd_s1819_UMF))
  (wtd_s1819_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass, wtd_s1819_UMF))
  #'  Anthropogenic model
  (wtd_s1819_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, wtd_s1819_UMF))
  #'  Combined model
  (wtd_s1819_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wtd_s1819_UMF))
  (wtd_s1819_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wtd_s1819_UMF))
  
  mods <- fitList(wtd_s1819_null, wtd_s1819_null2, wtd_s1819_terrain, wtd_s1819_terrain2, 
                  wtd_s1819_veg, wtd_s1819_terrain_veg, wtd_s1819_terrain_veg2, 
                  wtd_s1819_anthro, wtd_s1819_full, wtd_s1819_full2)
  modSel(mods)
  # Naive occupancy so high (0.92) not able to estimate effects of variables on psi
  
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only so no Area effect              
  #'  Null
  (wtd_w1820_null <- occu(~1 ~1, wtd_w1820_UMF))
  backTransform(wtd_w1820_null, 'det')
  backTransform(wtd_w1820_null, 'state')
  #'  Null with detection covariates
  (wtd_w1820_null2 <- occu(~Height*Distance + Trail + Year ~1, wtd_w1820_UMF))
  #'  Terrain model
  (wtd_w1820_terrain <- occu(~Height*Distance + Trail + Year ~Elev + Slope, wtd_w1820_UMF))
  (wtd_w1820_terrain2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope, wtd_w1820_UMF))
  #'  Vegetation model--------- NOTE: removed Xeric Shrub from veg models b/c very little Xeric Shrub in NE study area
  (wtd_w1820_veg <- occu(~Height*Distance + Trail + Year ~PercForMix + PercXGrass, wtd_w1820_UMF))
  #'  Habitat model
  (wtd_w1820_terrain_veg <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass, wtd_w1820_UMF))
  (wtd_w1820_terrain_veg2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass, wtd_w1820_UMF))
  #'  Anthropogenic model
  (wtd_w1820_anthro <- occu(~Height*Distance + Trail + Year ~RoadDensity + HumanMod, wtd_w1820_UMF))
  #'  Combined model
  (wtd_w1820_full <- occu(~Height*Distance + Trail + Year ~Elev + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wtd_w1820_UMF))
  (wtd_w1820_full2 <- occu(~Height*Distance + Trail + Year ~Elev + I(Elev^2) + Slope + PercForMix + PercXGrass + RoadDensity + HumanMod, wtd_w1820_UMF))
  
  mods <- fitList(wtd_w1820_null, wtd_w1820_null2, wtd_w1820_terrain, wtd_w1820_terrain2, 
                  wtd_w1820_veg, wtd_w1820_terrain_veg, wtd_w1820_terrain_veg2, 
                  wtd_w1820_anthro, wtd_w1820_full, wtd_w1820_full2)
  modSel(mods)
  
  
  
  ####  Goodness of Fit Test  ####
  #'  Checking model fit: Chi^2 test for binary data (1 df)
  #'  Chi^2 GOF is a one-sided test used to determine if a sample matches the 
  #'   population & asks what's the probability that the test statistic is greater  
  #'   than (falls above) the critical value on the tail of the Chi^2 distribution.
  #'  H0: there is no significant difference between the observed and the expected value
  #'  Ha: there is a significant difference between the observed and the expected value
  #'  If the Chi^2 test statistic is greater than the critical value derived from
  #'   the Chi^2 distribution (determined by degrees of freedom), then reject
  #'   the null & conclude there is a significant difference btwn observed & expected. 
  #'  If test statistic is less than the critical value then fail to reject & 
  #'   conclude there is no significant difference btwn observed & expected frequencies.
  #'  Use p-value with alpha-level 0.05 to reject/fail to reject the null hypothesis.
  #'  Failing to reject the null (p-value > 0.05) indicates adequate model fit.
  #'  
  #'  Code from unmarked vignette
  #'  Parametric Bootstrap Statistics:
  #'  t0 = Original statistic computed from data
  #'  t_B = Vector of bootstrap samples
  #'  If Pr(t_B > t0) is > 0.05 then fail to reject null & conclude model fit is adequate
  chisq <- function(fm) {
    umf <- getData(fm)
    y <- getY(umf)
    y[y>1] <- 1
    sr <- fm@sitesRemoved
    if(length(sr)>0)
      y <- y[-sr,,drop=FALSE]
    fv <- fitted(fm, na.rm=TRUE)
    y[is.na(fv)] <- NA
    sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
  }
  
  #'  Chi^2 test for models within 2 delta AIC for each species and season-- RUN WITH nsim = 10000 FOR FINAL ANALYSES
  (pb_bob_s1819 <- parboot(bob_s1819_anthro, statistic = chisq, nsim = 100, parallel = FALSE))         
  (pb_bob_w1820 <- parboot(bob_w1820_anthro, statistic = chisq, nsim = 100, parallel = FALSE))          
  (pb_coug_s1819 <- parboot(coug_s1819_terrain_veg2, statistic = chisq, nsim = 100, parallel = FALSE))         
  (pb_coug_w1820 <- parboot(coug_w1820_terrain_veg2, statistic = chisq, nsim = 100, parallel = FALSE))    #pval= 0.0396 Trash model
  (pb_coy_s1819 <- parboot(coy_s1819_terrain, statistic = chisq, nsim = 100, parallel = FALSE))   
  (pb_coy_w1820 <- parboot(coy_w1820_terrain, statistic = chisq, nsim = 100, parallel = FALSE))   
  (pb_wolf_s1819 <- parboot(wolf_s1819_terrain, statistic = chisq, nsim = 100, parallel = FALSE))      
  (pb_wolf_w1820 <- parboot(wolf_w1820_anthro, statistic = chisq, nsim = 100, parallel = FALSE))       
  (pb_elk_s1819 <- parboot(elk_s1819_null2, statistic = chisq, nsim = 100, parallel = FALSE))         
  (pb_elk_w1820 <- parboot(elk_w1820_null, statistic = chisq, nsim = 100, parallel = FALSE))          
  (pb_md_s1819 <- parboot(md_s1819_anthro, statistic = chisq, nsim = 100, parallel = FALSE))          
  (pb_md_w1820 <- parboot(md_w1820_terrain_veg2, statistic = chisq, nsim = 100, parallel = FALSE))             
  (pb_wtd_s1819 <- parboot(wtd_s1819_null2, statistic = chisq, nsim = 100, parallel = FALSE))          # wtd_s1819_terrain
  (pb_wtd_w1820 <- parboot(wtd_w1820_terrain_veg, statistic = chisq, nsim = 100, parallel = FALSE))        
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info

  #'  Function to save occupancy results
  occ_out <- function(mod, spp, season) {
    out <- summary(mod@estimates)$state %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$state),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter)
    return(out)
  }
  
  #'  Run each model through function
  bob_s1819_occ <- occ_out(bob_s1819_anthro, "Bobcat", "Summer")
  bob_w1820_occ <- occ_out(bob_w1820_anthro, "Bobcat", "Winter")
  coug_s1819_occ <- occ_out(coug_s1819_terrain_veg2, "Cougar", "Summer")
  coug_w1820_occ <- occ_out(coug_w1820_terrain_veg2, "Cougar", "Winter")
  coy_s1819_occ <- occ_out(coy_s1819_terrain, "Coyote", "Summer")
  coy_w1820_occ <- occ_out(coy_w1820_terrain, "Coyote", "Winter")
  wolf_s1819_occ <- occ_out(wolf_s1819_terrain, "Wolf", "Summer")
  wolf_w1820_occ <- occ_out(wolf_w1820_anthro, "Wolf", "Winter")
  elk_s1819_occ <- occ_out(elk_s1819_null2, "Elk", "Summer")
  elk_w1820_occ <- occ_out(elk_w1820_null, "Elk", "Winter")
  md_s1819_occ <- occ_out(md_s1819_anthro, "Mule Deer", "Summer")
  md_w1820_occ <- occ_out(md_w1820_terrain_veg2, "Mule Deer", "Winter")
  wtd_s1819_occ <- occ_out(wtd_s1819_null2, "White-tailed Deer", "Summer")
  wtd_w1820_occ <- occ_out(wtd_w1820_terrain_veg, "White-tailed Deer", "Winter")
  
  #'  Merge into larger data frames for easy comparison
  summer_occ <- rbind(bob_s1819_occ, coug_s1819_occ, coy_s1819_occ, wolf_s1819_occ,
                      elk_s1819_occ, md_s1819_occ, wtd_s1819_occ)
  winter_occ <- rbind(bob_w1820_occ, coug_w1820_occ, coy_w1820_occ, wolf_w1820_occ,
                      elk_w1820_occ, md_w1820_occ, wtd_w1820_occ)
  occ_results <- rbind(summer_occ, winter_occ) %>%
    arrange(Species)
  colnames(occ_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are easier to look at
  occ_results <- occ_results %>%
    mutate(
      Estimate = round(Estimate, 3),
      SE = round(SE, 3),
      z = round(z, 3),
      Pval = round(Pval, 3)
    )
  #'  Save!
  # write.csv(occ_results, paste0("./Outputs/OccMod_OccProb_Results_", Sys.Date(), ".csv"))
  
 
  #'  Function to save detection results
  det_out <- function(mod, spp, season) {
    out <- summary(mod@estimates)$det %>%
      mutate(
        Parameter = row.names(summary(mod@estimates)$det),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.))
      ) %>%
      relocate(Parameter, .before = Estimate) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter)
    return(out)
  }
  
  #'  Run each model through detection function
  bob_s1819_det <- det_out(bob_s1819_anthro, "Bobcat", "Summer")
  bob_w1820_det <- det_out(bob_w1820_anthro, "Bobcat", "Winter")
  coug_s1819_det <- det_out(coug_s1819_veg, "Cougar", "Summer")
  coug_w1820_det <- det_out(coug_w1820_terrain2, "Cougar", "Winter")
  coy_s1819_det <- det_out(coy_s1819_terrain, "Coyote", "Summer")
  coy_w1820_det <- det_out(coy_w1820_terrain, "Coyote", "Winter")
  wolf_s1819_det <- det_out(wolf_s1819_terrain, "Wolf", "Summer")
  wolf_w1820_det <- det_out(wolf_w1820_anthro, "Wolf", "Winter")
  elk_s1819_det <- det_out(elk_s1819_null2, "Elk", "Summer")
  elk_w1820_det <- det_out(elk_w1820_null, "Elk", "Winter")
  md_s1819_det <- det_out(md_s1819_anthro, "Mule Deer", "Summer")
  md_w1820_det <- det_out(md_w1820_terrain_veg2, "Mule Deer", "Winter")
  wtd_s1819_det <- det_out(wtd_s1819_null2, "White-tailed Deer", "Summer")
  wtd_w1820_det <- det_out(wtd_w1820_terrain_veg, "White-tailed Deer", "Winter")
  
  #'  Merge into larger data frames for easy comparison
  summer_det <- rbind(bob_s1819_det, coug_s1819_det, coy_s1819_det, wolf_s1819_det,
                      elk_s1819_det, md_s1819_det, wtd_s1819_det)
  winter_det <- rbind(bob_w1820_det, coug_w1820_det, coy_w1820_det, wolf_w1820_det,
                      elk_w1820_det, md_w1820_det, wtd_w1820_det)
  det_results <- rbind(summer_det, winter_det) %>%
    arrange(Species)
  colnames(det_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are a little easier to interpret
  det_results <- det_results %>%
    mutate(
      Estimate = round(Estimate, 3),
      SE = round(SE, 3),
      z = round(z, 3),
      Pval = round(Pval, 3)
    )
  
  #'  Save!
  # write.csv(det_results, paste0("./Outputs/OccMod_DetProb_Results_", Sys.Date(), ".csv"))







  


  
  
  