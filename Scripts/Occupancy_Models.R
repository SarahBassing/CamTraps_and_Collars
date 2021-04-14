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
    #'  Force km2road = 0 for cameras deployed on roads (not decommissioned roads)
    mutate(
      km2road = ifelse(Monitoring == "Dirt road", 0, km2road)
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
  
  
  plot(stations$Elev, stations$HumanModified) # seems to be a break point around mid-elevation
  
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
  
  #' #'  Get a feel for naive occupancy across covaraites
  #' smbob <- as.numeric(apply(DH_bob_smr1819, 1, max, na.rm=FALSE))
  #' wtrbob <- as.numeric(apply(DH_bob_wtr1820, 1, max, na.rm=FALSE))
  #' smcoug <- as.numeric(apply(DH_coug_smr1819, 1, max, na.rm=FALSE))
  #' wtrcoug <- as.numeric(apply(DH_coug_wtr1820, 1, max, na.rm=FALSE))
  #' smcoy <- as.numeric(apply(DH_coy_smr1819, 1, max, na.rm=FALSE))
  #' wtrcoy <- as.numeric(apply(DH_coy_wtr1820, 1, max, na.rm=FALSE))
  #' smwolf <- as.numeric(apply(DH_wolf_smr1819, 1, max, na.rm=FALSE))
  #' wtrwolf <- as.numeric(apply(DH_wolf_wtr1820, 1, max, na.rm=FALSE))
  #' smelk <- as.numeric(apply(DH_elk_smr1819, 1, max, na.rm=FALSE))
  #' wtrelk <- as.numeric(apply(DH_elk_wtr1820, 1, max, na.rm=FALSE))
  #' smmd <- as.numeric(apply(DH_md_smr1819, 1, max, na.rm=FALSE))
  #' wtrmd <- as.numeric(apply(DH_md_wtr1820, 1, max, na.rm=FALSE))
  #' smwtd <- as.numeric(apply(DH_wtd_smr1819, 1, max, na.rm=FALSE))
  #' wtrwtd <- as.numeric(apply(DH_wtd_wtr1820, 1, max, na.rm=FALSE))
  #' 
  #' #'  Raw detections by landcover categories
  #' landcov <- as.character(as.factor(stations$Landcov))
  #' landcovNE <- as.character(as.factor(stations_NE$Landcov))
  #' landcovOK <- as.character(as.factor(stations_OK$Landcov))
  #' hab_det <- as.data.frame(cbind(landcov, smbob, wtrbob, smcoug, wtrcoug, smcoy, wtrcoy, smwolf, wtrwolf)) %>%
  #'   mutate(
  #'     smbob = as.numeric(smbob),
  #'     wtrbob = as.numeric(wtrbob),
  #'     smcoug = as.numeric(smcoug),
  #'     wtrcoug = as.numeric(wtrcoug),
  #'     smcoy = as.numeric(smcoy),
  #'     wtrcoy = as.numeric(wtrcoy),
  #'     smwolf = as.numeric(smwolf),
  #'     wtrwolf = as.numeric(wtrwolf)
  #'   )
  #' habNE_det <- as.data.frame(cbind(landcovNE, smelk, wtrelk, smwtd, wtrwtd)) %>%
  #'   mutate(
  #'     smelk = as.numeric(smelk),
  #'     wtrelk = as.numeric(wtrelk),
  #'     smwtd = as.numeric(smwtd),
  #'     wtrwtd = as.numeric(wtrwtd)
  #'   )
  #' habOK_det <- as.data.frame(cbind(landcovOK, smmd, wtrmd)) %>%
  #'   mutate(
  #'     smmd = as.numeric(smmd),
  #'     wtrmd = as.numeric(wtrmd)
  #'   )
  #' #'  Count number of raw detections per landcover category
  #' tbl_pred <- hab_det %>%
  #'   group_by(landcov) %>%
  #'   summarise(across(smbob:wtrwolf, sum, na.rm= TRUE)) %>%
  #'   ungroup()
  #' tbl_NE <- habNE_det %>%
  #'   group_by(landcovNE) %>%
  #'   summarise(across(smelk:wtrwtd, sum, na.rm= TRUE)) %>%
  #'   ungroup()
  #' tbl_OK <- habOK_det %>%
  #'   group_by(landcovOK) %>%
  #'   summarise(across(smmd:wtrmd, sum, na.rm= TRUE)) %>%
  #'   ungroup()
  #' sum_raw_det <- full_join(tbl_pred, tbl_OK, by = c("landcov" = "landcovOK")) %>%
  #'   full_join(tbl_NE, by = c("landcov" = "landcovNE"))
  #' #'  Visualize another way
  #' barplot(sum_raw_det$smwtd, names.arg = sum_raw_det$landcov, 
  #'         xlab = "Landcover Type", ylab = "Raw Detections")
  #' 
  #' #'  Raw detections and continuous variables (remember- these are center & scaled)
  #' plot(stations$HumanModified, wtrcoy)
  #' plot(stations$PercForestMix, wtrcoy)
  #' hist(stations$PercForestMix[hab_det$wtrcoy == "1",], breaks = 50)
  #' hist(stations$PercXericShrub[habOK_det$wtrmd == "1",], breaks = 25)
  #' hist(stations$PercXericShrub[hab_det$wtrcoy == "1",], breaks = 25)

  
  
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
  #'  Running global model for all species and seasons.
  #'  Only removed Study Area covariate if data only collected in one area.
  #'  Only removed PercXShrub if converged poorly (usually due to no/few 
  #'  detections in areas with shrubland). 
  #'  =============================

  ####  BOBCAT MODELS  ####                   
  #'  SUMMERS 2018 & 2019
  (bob_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, bob_s1819_UMF))
    #'  WINTERS 2018-2019 & 2019-2020           
  (bob_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, bob_w1820_UMF))

  ####  COUGAR MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coug_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, coug_s1819_UMF))
  #'  WINTERS 2018-2019 & 2019-2020           
  (coug_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, coug_w1820_UMF))
  
  ####  COYOTE MODELS  ####
  #'  SUMMERS 2018 & 2019
  (coy_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, coy_s1819_UMF))
  #'  WINTERS 2018-2019 & 2019-2020              
  (coy_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, coy_w1820_UMF))

  ####  WOLF MODELS  ####
  #'  SUMMERS 2018 & 2019    
  #'  Removed PercXShrub in global2 models due to poor convergence, esp. summer model                 
  (wolf_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, wolf_s1819_UMF))
  (wolf_s1819_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod + Area, wolf_s1819_UMF))
      #'  WINTERS 2018-2019 & 2019-2020            
  (wolf_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area, wolf_w1820_UMF))
  (wolf_w1820_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod + Area, wolf_w1820_UMF))
  
  ####  ELK MODELS ####                          
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  #'  Removed PercXShrub in global2 models due to poor convergence                             
  (elk_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, elk_s1819_UMF))
  (elk_s1819_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod, elk_s1819_UMF))
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only so no Area effect                
  (elk_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, elk_w1820_UMF))
  (elk_w1820_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod, elk_w1820_UMF))
  
  ####  MULE DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  OK study area only so no Area effect
  (md_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, md_s1819_UMF))
  #'  WINTERS 2018-2019 & 2019-2020, OK study area only so no Area effect                    
  (md_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, md_w1820_UMF))

  ####  WHITE-TAILED DEER MODELS  ####
  #'  SUMMERS 2018 & 2019
  #'  NE study area only so no Area effect 
  #'  Removed PercXShrub in global2 models due to poor convergence, esp. on winter model
  (wtd_s1819_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, wtd_s1819_UMF))
  (wtd_s1819_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod, wtd_s1819_UMF))
  ####  PROBABLY NEED TO RUN THIS AS A LOGISTIC REGRESSION SINCE OCCUPANCY IS SO HIGH
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only so no Area effect              
  (wtd_w1820_global <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod, wtd_w1820_UMF))
  (wtd_w1820_global2 <- occu(~Trail + Height + Distance + Height*Distance + Year ~Elev + Slope + PercForMix + PercXGrass + NearestRd + HumanMod, wtd_w1820_UMF))
  

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
  bob_s1819_occ <- occ_out(bob_s1819_global, "Bobcat", "Summer")
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
  
  #'  Merge into larger data frames for easy comparison
  summer_occ <- rbind(bob_s1819_occ, coug_s1819_occ, coy_s1819_occ, wolf_s1819_occ,
                      elk_s1819_occ, md_s1819_occ, wtd_s1819_occ)
  winter_occ <- rbind(bob_w1820_occ, coug_w1820_occ, coy_w1820_occ, wolf_w1820_occ,
                      elk_w1820_occ, md_w1820_occ, wtd_w1820_occ)
  occ_results <- rbind(summer_occ, winter_occ) %>%
    arrange(Species)
  colnames(occ_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are easier to look at
  results_psi <- occ_results %>%
    mutate(
      Estimate = round(Estimate, 3),
      SE = round(SE, 3),
      z = round(z, 3),
      Pval = round(Pval, 3)
    )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_psi_wide <- results_psi %>%
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
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
    separate("NearestRd", c("NearestRd (SE)", "NearestRd Pval"), sep = "_") %>%
    separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Wolf", "Coyote", "Mule Deer", "Elk", "White-tailed Deer")))

  #'  Save!
  # write.csv(results_psi, paste0("./Outputs/OccMod_OccProb_Results_", Sys.Date(), ".csv"))
  # write.csv(results_psi_wide, paste0("./Outputs/OccMod_OccProb_Results_wide_", Sys.Date(), ".csv"))
  
 
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
  bob_s1819_det <- det_out(bob_s1819_global, "Bobcat", "Summer")
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
  
  #'  Merge into larger data frames for easy comparison
  summer_det <- rbind(bob_s1819_det, coug_s1819_det, coy_s1819_det, wolf_s1819_det,
                      elk_s1819_det, md_s1819_det, wtd_s1819_det)
  winter_det <- rbind(bob_w1820_det, coug_w1820_det, coy_w1820_det, wolf_w1820_det,
                      elk_w1820_det, md_w1820_det, wtd_w1820_det)
  det_results <- rbind(summer_det, winter_det) %>%
    arrange(Species)
  colnames(det_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")

  #'  Round so numbers are a little easier to interpret
  results_det <- det_results %>%
    mutate(
      Estimate = round(Estimate, 3),
      SE = round(SE, 3),
      z = round(z, 3),
      Pval = round(Pval, 3)
    )
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  results_det_wide <- results_det %>%
    dplyr::select(-z) %>%
    mutate(
      SE = round(SE, 2),
      SE = paste0("(", SE, ")")
    ) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_Pval, Est_SE, Pval, sep = "_") %>%
    spread(Parameter, Est_SE_Pval) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept Pval"), sep = "_") %>%
    separate("TrailDirt road", c("Road (SE)", "Road Pval"), sep = "_") %>%
    separate("TrailDecommissioned road", c("Decom Road (SE)", "Decom Road Pval"), sep = "_") %>%
    separate("Height", c("Height (SE)", "Height Pval"), sep = "_") %>%
    separate("Distance", c("Distance (SE)", "Distance Pval"), sep = "_") %>%
    separate("Height:Distance", c("Height*Distance (SE)", "Height*Distance Pval"), sep = "_") %>%
    separate("YearYear2", c("Year2 (SE)", "Year2 Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Wolf", "Coyote", "Mule Deer", "Elk", "White-tailed Deer")))

  #'  Save!
  # write.csv(results_det, paste0("./Outputs/OccMod_DetProb_Results_", Sys.Date(), ".csv"))
  # write.csv(results_det_wide, paste0("./Outputs/OccMod_DetProb_Results_wide", Sys.Date(), ".csv"))

  #'  Save workspace
  # save.image(file = "CamvCollar_OccMods.RData")




  


  
  
  