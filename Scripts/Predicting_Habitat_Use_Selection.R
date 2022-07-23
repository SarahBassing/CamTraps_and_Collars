  #'  ============================================================
  #'  Predicting Probability of Habitat Use and Relative Selection
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  September 2021
  #'  ============================================================
  #'  Script to predict probability of use and relative probability of selection
  #'  for each species across the relevant study area(s) by season based on 
  #'  occupancy model and RSF results. Script scales covariate data across both
  #'  study areas based on mean & SD values from each species & season specific 
  #'  RSF and the mean & SD values from the camera sites. Then predicts prob.
  #'  of use and relative selection and calculates pixel by pixel correlation
  #'  between results. Finally, script maps results for visual comparison and
  #'  saves raster files for predicted probabilities.
  #'  ============================================================

  #'  Clear memory
  rm(list=ls())
  
  #'  Load libraries
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(png)
  library(RCurl)
  library(RColorBrewer)
  library(rphylopic)
  library(sf)
  library(raster)
  library(tidyverse)
  
  ####  Helpful spatial data  ####
  #'  Define projections
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Read in study area data and reproject
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK_SA$NAME <- "Okanogan"
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE_SA$NAME <- "Northeast"
  
  projection(OK_SA)
  extent(OK_SA)
  extent(NE_SA)

  
  ####  Load Model Results  ####
  #'  Occupancy model output
  occ_out <- read.csv("./Outputs/Tables/OccMod_OccProb_Results_matchRSF_2022-05-27.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE
    #'  Calculate 90% confidence intervals to mirror alpha = 0.1
    mutate(
      l95 = (Estimate - (1.64 * SE)),  #### REMINDER: this is 90% CI even though column says l95/u95
      u95 = (Estimate + (1.64 * SE))   
    ) %>%
    dplyr::select(-c(X)) #, Model
  
  #'  RSF results output
  rsf_out <- read.csv("./Outputs/Tables/RSF_Results_BuffHR_2022-05-03.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE 2021-09-14
    #'  Calculate 95% confidence intervals to mirror alpha = 0.05
    mutate(
      l95 = (Estimate - (1.96 * SE)),  #### REMINDER: this is 95% CI
      u95 = (Estimate + (1.96 * SE))
    ) %>%
    dplyr::select(-X)


  ####  Center & Scale Covariate Data  ####
  #'  Read in original covariate data from occupancy models
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18-20_2021-08-10.csv") %>%
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    transmute(
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      CameraLocation = as.factor(CameraLocation),
      PercForMix = PercForestMix2, 
      PercXShrub = PercXericShrub,
      PercXGrass = PercXericGrass,
      Elev = as.numeric(Elev), 
      Slope = Slope, 
      RoadDen = RoadDen, 
      HumanMod = HumanMod 
    ) %>%
    arrange(Year)
  
  #'  Read in original covariate data from RSFs
  load("./Outputs/RSF_pts/md_dat_2nd_buffHR_all_2022-04-18.RData")  
  load("./Outputs/RSF_pts/elk_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/wtd_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/coug_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/wolf_dat_2nd_buffHR_all_2022-04-18.RData") 
  load("./Outputs/RSF_pts/bob_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/coy_dat_2nd_buffHR_all_2022-04-18.RData")
  
  #'  Function to find mean and standard deviation for each covariate
  #'  Used when center & scaling covariates for original models and will be used 
  #'  to standardize covariate data when predicting across study areas
  cov_summary <- function(covs) {
    mu.cov <- covs %>%
      summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE)))
    sd.cov <- covs %>%
      summarise(across(where(is.numeric), ~sd(.x, na.rm = TRUE)))
    mu.sd.cov <- rbind(mu.cov, sd.cov)
    parameter <- as.data.frame(c("Mean", "SD"))
    colnames(parameter) <- "Parameter"
    cov_summary <- cbind(parameter, mu.sd.cov)
    return(cov_summary)
  }
  #'  Find mean & SD across all camera locations
  summary_occ_covs <- cov_summary(stations)
  #'  Find mean & SD across all used & available collar data, but separately by season & spp
  summary_md_smr <- cov_summary(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",])
  summary_md_wtr <- cov_summary(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",])
  summary_elk_smr <- cov_summary(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19",])
  summary_elk_wtr <- cov_summary(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920",])
  summary_wtd_smr <- cov_summary(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19",])
  summary_wtd_wtr <- cov_summary(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920",])
  summary_coug_smr <- cov_summary(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",])
  summary_coug_wtr <- cov_summary(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",])
  summary_wolf_smr <- cov_summary(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",])
  summary_wolf_wtr <- cov_summary(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",])
  summary_bob_smr <- cov_summary(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",])
  summary_bob_wtr <- cov_summary(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",])
  summary_coy_smr <- cov_summary(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19",])
  summary_coy_wtr <- cov_summary(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920",])
  
  #'  Read in & scale covariate data for entire study areas (based on 30m & 1km grids)
  #'  Data generated in "Covariate_Extraction.R" script in WPPP_CameraTrapping Project
  # NE_covs_30m <- read.csv("./Outputs/Tables/StudyAreaWide_NE_Covariates_30m_2021-09-15.csv") %>%
  #   dplyr::select(-X) %>%
  #   mutate(Area = 0)
  # OK_covs_30m <- read.csv("./Outputs/Tables/StudyAreaWide_OK_Covariates_30m_2021-09-15.csv") %>%
  #   dplyr::select(-X) %>%
  #   mutate(Area = 1)
  # all_covs_30m <- as.data.frame(rbind(NE_covs_30m, OK_covs_30m)) 
  NE_covs_1km <- read.csv("./Outputs/Tables/StudyAreaWide_NE_Covariates_1km_2021-09-16.csv") %>% 
    dplyr::select(-X) %>%
    mutate(Area = 0)
  OK_covs_1km <- read.csv("./Outputs/Tables/StudyAreaWide_OK_Covariates_1km_2021-09-16.csv") %>% 
  dplyr::select(-X) %>%
    mutate(Area = 1)
  all_covs_1km <- as.data.frame(rbind(NE_covs_1km, OK_covs_1km))
  
  #'  Exclude elevations >2150 for camera covariates since didn't sample above this elevation
  # all_covs_adj_30m <- all_covs_30m %>%
  #   mutate(Elev = ifelse(Elev > 2150, NA, Elev))
  all_covs_adj_1km <- all_covs_1km %>%
    mutate(Elev = ifelse(Elev > 2150, NA, Elev))
  
  #'  Function to scale the individual covariates based on the covariate-specific 
  #'  mean and SD values used to standardize the covariate data in the original
  #'  models. This will differ by species and season for the RSF predictions since
  #'  the used and available locations differed with each model. This will be the
  #'  same for all occupancy model predictions since the same camera sites were
  #'  included in all analyses.
  scaling_covs <- function(covs, mu.sd) {
    scaling_covs <- covs %>%
      transmute(
        obs = obs,
        Elev = (Elev - mu.sd$Elev[1]) / mu.sd$Elev[2],
        Slope = (Slope - mu.sd$Slope[1]) / mu.sd$Slope[2],
        RoadDen = (RoadDen - mu.sd$RoadDen[1]) / mu.sd$RoadDen[2],
        HumanMod = (HumanMod - mu.sd$HumanMod[1]) / mu.sd$HumanMod[2],
        PercForMix = (PercForestMix2 - mu.sd$PercForMix[1]) / mu.sd$PercForMix[2],
        PercXGrass = (PercXericGrass - mu.sd$PercXGrass[1]) / mu.sd$PercXGrass[2],
        PercXShrub = (PercXericShrub - mu.sd$PercXShrub[1]) / mu.sd$PercXShrub[2],
        x = x,
        y = y,
        Area = Area)
  }
  #'  Standardize covariate data based on camera covariate means & SDs
  cam_zcovs <- scaling_covs(all_covs_adj_1km, summary_occ_covs)  # all_covs_adj_30m
  #'  Standardize covariate data based on collar covariate means & SDs for specific species & seasons
  md_smr_zcovs <- scaling_covs(all_covs_1km, summary_md_smr)  # all_covs_30m
  md_wtr_zcovs <- scaling_covs(all_covs_1km, summary_md_wtr)
  elk_smr_zcovs <- scaling_covs(all_covs_1km, summary_elk_smr)
  elk_wtr_zcovs <- scaling_covs(all_covs_1km, summary_elk_wtr)
  wtd_smr_zcovs <- scaling_covs(all_covs_1km, summary_wtd_smr)
  wtd_wtr_zcovs <- scaling_covs(all_covs_1km, summary_wtd_wtr)
  coug_smr_zcovs <- scaling_covs(all_covs_1km, summary_coug_smr)
  coug_wtr_zcovs <- scaling_covs(all_covs_1km, summary_coug_wtr)
  wolf_smr_zcovs <- scaling_covs(all_covs_1km, summary_wolf_smr)
  wolf_wtr_zcovs <- scaling_covs(all_covs_1km, summary_wolf_wtr)
  bob_smr_zcovs <- scaling_covs(all_covs_1km, summary_bob_smr)
  bob_wtr_zcovs <- scaling_covs(all_covs_1km, summary_bob_wtr)
  coy_smr_zcovs <- scaling_covs(all_covs_1km, summary_coy_smr)
  coy_wtr_zcovs <- scaling_covs(all_covs_1km, summary_coy_wtr)
  
  #'  Function to covert all covariate data into SpatialPointsDataFrames
  sp_covs <- function(covs) {
    xy <- covs[,c(10,11)]
    # xy <- covs[,c(9,10)]
    covs_spdf <- SpatialPointsDataFrame(data = covs, coords = xy,
                                        proj4string = CRS(sa_proj))
  }
  #'  Gather all scaled covariate data into a monster list
  zcovs_list <- list(cam_zcovs, md_smr_zcovs, md_wtr_zcovs, 
                elk_smr_zcovs, elk_wtr_zcovs, wtd_smr_zcovs, 
                wtd_wtr_zcovs, coug_smr_zcovs, coug_wtr_zcovs, 
                wolf_smr_zcovs, wolf_wtr_zcovs, bob_smr_zcovs, 
                bob_wtr_zcovs, coy_smr_zcovs, coy_wtr_zcovs)
  #'  Run list of scaled covariates through spatial function
  zcovs <- lapply(zcovs_list, sp_covs) #sp_zcovs
  
  

  ####  Predict probability of use across study areas  ####
  #'  Snag the intercept for each species (needed below)
  occ_int <- occ_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename intercept
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept)
  
  #'  Manipulate Occupancy result tables
  #'  Exclude all non-significant coefficients by forcing to 0
  occ_coefs_signif <- occ_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Problem is many of the intercepts were reduced to 0 which we don't want 
    #'  so we need to swap this intercept column with the original intercepts 
    #'  using the occ_int data frame above 
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = occ_int$alpha,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDensity,
      B.area = AreaOK) %>% #B.hm = HumanMod
    #'  Change NAs to 0 (no effect) for coefficients not included in species-specific models
    mutate(
      B.elev = ifelse(is.na(B.elev), 0, B.elev),
      B.slope = ifelse(is.na(B.slope), 0, B.slope),
      B.for = ifelse(is.na(B.for), 0, B.for),
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.rd = ifelse(is.na(B.rd), 0, B.rd),
      B.area = ifelse(is.na(B.area), 0, B.area)
    )
  
  #'  Function to predict across all grid cells based on occupancy model results
  #'  Should end up with 1 predicted value per grid cell
  predict_occ <- function(cov, coef) {
    predict_odds <- c()
    predict_prob <- c()
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$alpha + coef$B.area*cov$Area[i] + coef$B.elev*cov$Elev[i] + 
                               coef$B.slope*cov$Slope[i]+ coef$B.for*cov$PercForMix[i] + 
                               coef$B.grass*cov$PercXGrass[i] + coef$B.shrub*cov$PercXShrub[i] + 
                               coef$B.rd*cov$RoadDen[i]) #+ coef$B.hm*cov$HumanMod[i]
      predict_prob[i] <- predict_odds[i] / (1 + predict_odds[i])
    }
    predict_prob <- as.data.frame(predict_prob) %>%
      transmute(
        Predicted_Occ = predict_prob
      )
    return(predict_prob)
  }
  #'  Run estimates from occupancy sub-model through function to predict probability of use
  #'  BE SURE TO USE THE RIGHT COVARIATE DATA! sp_zcovs[[1]] are covariates from camera sites
  #'  Includes ONLY significant coefficients per species and season
  md_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Summer",])
  md_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Winter",])
  elk_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Summer",])
  elk_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Winter",])
  wtd_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]][zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Winter",])
  coug_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Summer",])
  coug_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Winter",])
  wolf_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Winter",])
  bob_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Summer",])
  bob_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Winter",])
  coy_smr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Summer",])
  coy_wtr_predict_occ_sgnf <- predict_occ(zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predator data that spans both study areas
  Predicted_occ <- as.data.frame(zcovs[[1]]) %>%  # NOTE WHICH GRID YOU'RE USING!
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(coug_smr_predict_occ_sgnf, coug_wtr_predict_occ_sgnf, 
          wolf_smr_predict_occ_sgnf, wolf_wtr_predict_occ_sgnf, 
          bob_smr_predict_occ_sgnf, bob_wtr_predict_occ_sgnf, 
          coy_smr_predict_occ_sgnf, coy_wtr_predict_occ_sgnf)
  #'  Make sure you have the order right when you change the names!!!
  colnames(Predicted_occ) <- c("obs", "Area", "x", "y",  
                               "COUG_smr_occ", "COUG_wtr_occ", 
                               "WOLF_smr_occ", "WOLF_wtr_occ", "BOB_smr_occ", 
                               "BOB_wtr_occ", "COY_smr_occ", "COY_wtr_occ")
  
  #'  Okanogan-only data (mule deer)
  # OK_rows <- seq(1:nrow(md_smr_predict_occ_sgnf))
  # Area <- rep("Okanogan", length(OK_rows))
  OK_occ <- as.data.frame(cbind(md_smr_predict_occ_sgnf, md_wtr_predict_occ_sgnf)) # OK_rows, Area, 
  colnames(OK_occ) <- c("MD_smr_occ", "MD_wtr_occ") #"obs", "Area", 
  
  #'  Northeast-only data (elk & white-tailed deer)
  # NE_rows <- seq(1:nrow(elk_smr_predict_occ_sgnf))
  # Area <- rep("Northeast", length(NE_rows))
  NE_occ <- as.data.frame(cbind(elk_smr_predict_occ_sgnf, elk_wtr_predict_occ_sgnf, 
                                wtd_smr_predict_occ_sgnf, wtd_wtr_predict_occ_sgnf)) # NE_rows, Area,  
  colnames(NE_occ) <- c("ELK_smr_occ", "ELK_wtr_occ", "WTD_smr_occ", "WTD_wtr_occ") # "obs", "Area", 
  
  #'  Merge ungulate & predator data by study area
  Predicted_occ_OK <- Predicted_occ[Predicted_occ$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_occ = NA,
      ELK_wtr_occ = NA,
      WTD_smr_occ = NA,
      WTD_wtr_occ = NA
    ) %>%
    cbind(OK_occ) 
  Predicted_occ_NE <- Predicted_occ[Predicted_occ$Area == "Northeast",] %>%
    cbind(NE_occ) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_occ = NA,
      MD_wtr_occ = NA
    )
  
  #'  Merge NE and OK predictions together
  Predicted_occ <- as.data.frame(rbind(Predicted_occ_NE, Predicted_occ_OK)) #%>%
    # mutate(Obs = seq(1:nrow(.))) %>%
    # dplyr::select(-obs) %>%
    # relocate(Obs, .before = "Area")
  
  #'  Save
  write.csv(Predicted_occ, paste0("./Outputs/Tables/Predicted_Prob_Occupancy_", Sys.Date(), ".csv"))

  
  ####  Predict relative probability of selesction across study areas  ####
  #'  Snag the intercept for each species (needed below)
  rsf_int <- rsf_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept)
  
  #'  Manipulate RSF result tables
  #'  Exclude all non-significant coefficients by forcing to 0
  rsf_coefs_signif <- rsf_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  No intercepts should have been changed here but just in case- swap current
    #'  intercept column with the intercept from above
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = rsf_int$alpha,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDen) %>% # B.hm = HumanMod
    mutate(
      B.elev = ifelse(is.na(B.elev), 0, B.elev),
      B.slope = ifelse(is.na(B.slope), 0, B.slope),
      B.for = ifelse(is.na(B.for), 0, B.for),
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.rd = ifelse(is.na(B.rd), 0, B.rd)
    )
  
  #'  Function to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  #'  NOTE: I want the predict relative probability of selection from RSF so not 
  #'  using a logit transformation like I normally would with logistic regression.
  #'  Instead, dropping the intercept from the model and just exponentiating the 
  #'  coeffs*covs (Fieberg et al. 2020)
  predict_rsf <- function(cov, coef) {
    predict_odds <- c()
    predict_prob <- c()
    #'  Normalizing constant so values are constrained between 0 & 1
    # c <- exp(coef$alpha)
    #'  Predict across each grid cell
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$B.elev*cov$Elev[i] + coef$B.slope*cov$Slope[i] + 
                               coef$B.for*cov$PercForMix[i] + coef$B.grass*cov$PercXGrass[i] + 
                               coef$B.shrub*cov$PercXShrub[i] + coef$B.rd*cov$RoadDen[i]) #+ coef$B.hm*cov$HumanMod[i])
      #'  If we were back-transforming like a normal logistic regression
      # predict_prob[i] <- predict_odds[i] / (1 + predict_odds[i])
    }
    predict_rsf <- as.data.frame(cbind(predict_prob, predict_odds)) #%>%
    #' transmute(
    #'   # Predicted_RSF = predict_odds,
    #'   # Logit_RSF = predict_prob,
    #'   #'  Multiply predictions by normalizing constant
    #'  normalized_RSF = c*predict_odds
    #'  )
    return(predict_rsf)
  }
  #'  Run estimated coefficients from RSFs through function to predict relative probability of selection
  #'  BE SURE TO USE THE RIGHT COVARIATE DATA! 
  #'  e.g., sp_zcovs[[2]] = covs standardized based on summer mule deer mean & SDs (includes both study areas)
  #'  e.g., sp_zcovs[[15]] = covs standardized based on winter coyote mean & SDs (inclused both study areas)
  #'  Includes ONLY significant coefficients per species and season
  md_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[2]][zcovs[[2]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Summer",])
  md_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[3]][zcovs[[3]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Winter",])
  elk_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[4]][zcovs[[4]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Summer",])
  elk_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[5]][zcovs[[5]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Winter",])
  wtd_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[6]][zcovs[[6]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[7]][zcovs[[7]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Winter",])
  coug_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[8]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Summer",])
  coug_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[9]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Winter",])
  wolf_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[10]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[11]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Winter",])
  bob_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[12]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Summer",])
  bob_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[13]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Winter",])
  coy_smr_predict_rsf_sgnf <- predict_rsf(zcovs[[14]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Summer",])
  coy_wtr_predict_rsf_sgnf <- predict_rsf(zcovs[[15]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predators
  Predicted_rsf <- as.data.frame(all_covs_1km) %>% # NOTE WHICH GRID YOU'RE USING!
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(coug_smr_predict_rsf_sgnf, coug_wtr_predict_rsf_sgnf, # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o non-signif coeffs, w/ or w/o normalizing constant)
          wolf_smr_predict_rsf_sgnf, wolf_wtr_predict_rsf_sgnf, 
          bob_smr_predict_rsf_sgnf, bob_wtr_predict_rsf_sgnf, 
          coy_smr_predict_rsf_sgnf, coy_wtr_predict_rsf_sgnf)
  #'  Make sure you have the order right when you change the names!!!
  colnames(Predicted_rsf) <- c("obs", "Area", "x", "y",  
                               "COUG_smr_rsf", "COUG_wtr_rsf", 
                               "WOLF_smr_rsf", "WOLF_wtr_rsf", "BOB_smr_rsf", 
                               "BOB_wtr_rsf", "COY_smr_rsf", "COY_wtr_rsf")
  
  #'  Okanogan-only data (mule deer)
  # OK_rows <- seq(1:nrow(md_smr_predict_rsf_sgnf))
  # Area <- rep("Okanogan", length(OK_rows))
  OK_rsf <- as.data.frame(cbind(md_smr_predict_rsf_sgnf, md_wtr_predict_rsf_sgnf)) # OK_rows, Area, 
  colnames(OK_rsf) <- c("MD_smr_rsf", "MD_wtr_rsf") #"obs", "Area", 
  
  #'  Northeast-only data (elk & white-tailed deer)
  # NE_rows <- seq(1:nrow(elk_smr_predict_rsf_sgnf))
  # Area <- rep("Northeast", length(NE_rows))
  NE_rsf <- as.data.frame(cbind(elk_smr_predict_rsf_sgnf, elk_wtr_predict_rsf_sgnf, # NE_rows, Area, 
                                wtd_smr_predict_rsf_sgnf, wtd_wtr_predict_rsf_sgnf)) 
  colnames(NE_rsf) <- c("ELK_smr_rsf", "ELK_wtr_rsf", "WTD_smr_rsf", "WTD_wtr_rsf") #"obs", "Area", 
  
  #'  Merge ungulate & predator data by study area
  Predicted_rsf_OK <- Predicted_rsf[Predicted_rsf$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_rsf = NA,
      ELK_wtr_rsf = NA,
      WTD_smr_rsf = NA,
      WTD_wtr_rsf = NA
    ) %>%
    cbind(OK_rsf) 
  Predicted_rsf_NE <- Predicted_rsf[Predicted_rsf$Area == "Northeast",] %>%
    cbind(NE_rsf) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_rsf = NA,
      MD_wtr_rsf = NA
    )
  
  #'  Merge NE and OK predictions togther
  Predicted_rsf <- as.data.frame(rbind(Predicted_rsf_NE, Predicted_rsf_OK)) #%>%
    # mutate(Obs = seq(1:nrow(.))) %>%
    # dplyr::select(-obs) %>%
    # relocate(Obs, .before = "Area")
  
  
  #'  Function to identify any outliers
  outliers <- function(predicted, title) { #, covs_list
    #'  Summarize predicted values
    hist(predicted, breaks = 100, main = title)
    boxplot(predicted, main = title)
    #'  What value represents the 99th percentile in the predicted RSF values
    quant <- quantile(predicted, c(0.99), na.rm = TRUE)
    #'  Print that value and maximum prediction
    print(quant); print(max(predicted, na.rm = TRUE))
    #'  Identify the 1% most extreme values and set to 99th percentile value
    predicted <- as.data.frame(predicted) %>%
      mutate(outlier = ifelse(predicted > quant, "outlier", "not_outlier"),
             adjusted_rsf = ifelse(outlier == "outlier", quant, predicted))
    #'  How many predicted values are above the 99th percentile?
    outlier <- predicted[predicted$outlier == "outlier",]
    outlier <- filter(outlier, !is.na(outlier))
    print(nrow(outlier))
    adjusted_rsf <- predicted$adjusted_rsf
    
    return(adjusted_rsf)
  }
  #'  Identify outlier predictions and possible covariates associated with those
  #'  Be sure to used standardized covariates for evaluation
  Predicted_rsf$MD_smr_rsf2 <- outliers(Predicted_rsf$MD_smr_rsf, "Mule Deer Summer RSF Predictions")
  Predicted_rsf$MD_wtr_rsf2 <- outliers(Predicted_rsf$MD_wtr_rsf, "Mule Deer Winter RSF Predictions")
  Predicted_rsf$ELK_smr_rsf2 <- outliers(Predicted_rsf$ELK_smr_rsf, "Elk Summer RSF Predictions")
  Predicted_rsf$ELK_wtr_rsf2 <- outliers(Predicted_rsf$ELK_wtr_rsf, "Elk Winter RSF Predictions")
  Predicted_rsf$WTD_smr_rsf2 <- outliers(Predicted_rsf$WTD_smr_rsf, "White-tailed Deer Summer RSF Predictions")
  Predicted_rsf$WTD_wtr_rsf2 <- outliers(Predicted_rsf$WTD_wtr_rsf, "White-tailed Deer Winter RSF Predictions")
  Predicted_rsf$COUG_smr_rsf2 <- outliers(Predicted_rsf$COUG_smr_rsf, "Cougar Summer RSF Predictions")
  Predicted_rsf$COUG_wtr_rsf2 <- outliers(Predicted_rsf$COUG_wtr_rsf, "Cougar Winter RSF Predictions")
  Predicted_rsf$WOLF_smr_rsf2 <- outliers(Predicted_rsf$WOLF_smr_rsf, "Wolf Summer RSF Predictions")
  Predicted_rsf$WOLF_wtr_rsf2 <- outliers(Predicted_rsf$WOLF_wtr_rsf, "Wolf Winter RSF Predictions")
  Predicted_rsf$COY_smr_rsf2 <- outliers(Predicted_rsf$COY_smr_rsf, "Coyote Summer RSF Predictions")
  Predicted_rsf$COY_wtr_rsf2 <- outliers(Predicted_rsf$COY_wtr_rsf, "Coyote Winter RSF Predictions")
  Predicted_rsf$BOB_smr_rsf2 <- outliers(Predicted_rsf$BOB_smr_rsf, "Bobcat Summer RSF Predictions")
  Predicted_rsf$BOB_wtr_rsf2 <- outliers(Predicted_rsf$BOB_wtr_rsf, "Bobcat Winter RSF Predictions")
  
  #' #'  Identify any outliers
  #' outliers <- function(predicted, title) {
  #'   print(summary(predicted))
  #'   hist(predicted, breaks = 100, main = title)
  #'   boxplot(predicted, main = title)
  #' }
  #' #'  Identify outlier perdictions
  #' outliers(Predicted_rsf$MD_smr_rsf, "Mule Deer Summer RSF Predictions")
  #' outliers(Predicted_rsf$MD_wtr_rsf, "Mule Deer Winter RSF Prediction")
  #' outliers(Predicted_rsf$ELK_smr_rsf, "Elk Summer RSF Predictions")
  #' outliers(Predicted_rsf$ELK_wtr_rsf, "Elk Winter RSF Prediction")
  #' outliers(Predicted_rsf$WTD_smr_rsf, "White-tailed Deer Summer RSF Predictions")
  #' outliers(Predicted_rsf$WTD_wtr_rsf, "White-tailed Deer Winter RSF Prediction")
  #' outliers(Predicted_rsf$COUG_smr_rsf, "Cougar Summer RSF Predictions")
  #' outliers(Predicted_rsf$COUG_wtr_rsf, "Cougar Winter RSF Prediction")
  #' outliers(Predicted_rsf$WOLF_smr_rsf, "Wolf Summer RSF Predictions")
  #' outliers(Predicted_rsf$WOLF_wtr_rsf, "Wolf Winter RSF Prediction")
  #' outliers(Predicted_rsf$BOB_smr_rsf, "Bobcat Summer RSF Predictions")
  #' outliers(Predicted_rsf$BOB_wtr_rsf, "Bobcat Winter RSF Prediction")
  #' outliers(Predicted_rsf$COY_smr_rsf, "Coyote Summer RSF Predictions")
  #' outliers(Predicted_rsf$COY_wtr_rsf, "Coyote Winter RSF Prediction")
  #' 
  #' #'  Exclude extreme outliers as identified by histograms & boxplots
  #' Predicted_rsf <- Predicted_rsf %>%
  #'   mutate(
  #'     MD_smr_rsf2 = ifelse(MD_smr_rsf > 5, NA, MD_smr_rsf),
  #'     MD_wtr_rsf2 = ifelse(MD_wtr_rsf > 15, NA, MD_wtr_rsf), 
  #'     # ELK_wtr_rsf2 = ifelse(ELK_wtr_rsf > 3, NA, ELK_wtr_rsf),  
  #'     WTD_smr_rsf2 = ifelse(WTD_smr_rsf > 6, NA, WTD_smr_rsf),
  #'     WTD_wtr_rsf2 = ifelse(WTD_wtr_rsf > 6, NA, WTD_wtr_rsf),
  #'     COUG_wtr_rsf2 = ifelse(COUG_wtr_rsf > 10, NA, COUG_wtr_rsf), 
  #'     WOLF_smr_rsf2 = ifelse(WOLF_smr_rsf > 9, NA, WOLF_smr_rsf),
  #'     # WOLF_wtr_rsf2 = ifelse(WOLF_wtr_rsf > 5, NA, WOLF_wtr_rsf),
  #'     BOB_smr_rsf2 = ifelse(BOB_smr_rsf > 6, NA, BOB_smr_rsf),
  #'     BOB_wtr_rsf2 = ifelse(BOB_wtr_rsf > 6, NA, BOB_wtr_rsf),
  #'     COY_wtr_rsf2 = ifelse(COY_wtr_rsf > 5, NA, COY_wtr_rsf)
  #'   )

 
  #'  Merge all predictions together (with unscaled RSF predictions)
  Predicted_Occ_RSF <- Predicted_occ %>%
    full_join(Predicted_rsf, by = c("obs", "Area", "x", "y"))
  write.csv(Predicted_Occ_RSF, paste0("./Outputs/Tables/Predictions_OccMod_v_RSF_buffHR_", Sys.Date(), ".csv"))  
  
  
  
  ####  Calculate Correlations between OccMod & RSF Predictions ####
  #'  Evaluate correlation between predicted space use for each paired set of models
  predict_corr <- function(predict_occu, predict_rsfs) {
    #'  Identify the maximum value in the predicted RSFs
    m <- max(predict_rsfs, na.rm = TRUE)
    #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
    stand_rsf <- as.data.frame(predict_rsfs) %>%
      mutate(scaled_rsf = predict_rsfs/m)
    #'  Combine occupancy and re-scaled RSF values
    predicted <- as.data.frame(cbind(predict_occu, stand_rsf$scaled_rsf)) %>%
      mutate(CellID = seq(1:nrow(.)))
    colnames(predicted) <- c("Occ_predictions", "RSF_rs_predictions", "CellID")
    #'  Calculated correlation between occupancy and rsf predictions
    pred_corr <- cor(predicted$Occ_predictions, predicted$RSF_rs_predictions, use = "complete.obs")
    #' #'  Plot correlation
    #' plot(predicted$Occ_predictions, predicted$RSF_rs_predictions)
    #' abline(a = 0, b = 1, col = "red")
    return(pred_corr)
  }
  #'  Run each predictions through function
  md_smr_corr <- predict_corr(Predicted_occ$MD_smr_occ, Predicted_rsf$MD_smr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  md_wtr_corr <- predict_corr(Predicted_occ$MD_wtr_occ, Predicted_rsf$MD_wtr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  elk_smr_corr <- predict_corr(Predicted_occ$ELK_smr_occ, Predicted_rsf$ELK_smr_rsf2)  
  elk_wtr_corr <- predict_corr(Predicted_occ$ELK_wtr_occ, Predicted_rsf$ELK_wtr_rsf2)  
  wtd_smr_corr <- predict_corr(Predicted_occ$WTD_smr_occ, Predicted_rsf$WTD_smr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  wtd_wtr_corr <- predict_corr(Predicted_occ$WTD_wtr_occ, Predicted_rsf$WTD_wtr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  coug_smr_corr <- predict_corr(Predicted_occ$COUG_smr_occ, Predicted_rsf$COUG_smr_rsf2)
  coug_wtr_corr <- predict_corr(Predicted_occ$COUG_wtr_occ, Predicted_rsf$COUG_wtr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  wolf_smr_corr <- predict_corr(Predicted_occ$WOLF_smr_occ, Predicted_rsf$WOLF_smr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  wolf_wtr_corr <- predict_corr(Predicted_occ$WOLF_wtr_occ, Predicted_rsf$WOLF_wtr_rsf2)  
  bob_smr_corr <- predict_corr(Predicted_occ$BOB_smr_occ, Predicted_rsf$BOB_smr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  bob_wtr_corr <- predict_corr(Predicted_occ$BOB_wtr_occ, Predicted_rsf$BOB_wtr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  coy_smr_corr <- predict_corr(Predicted_occ$COY_smr_occ, Predicted_rsf$COY_smr_rsf2)
  coy_wtr_corr <- predict_corr(Predicted_occ$COY_wtr_occ, Predicted_rsf$COY_wtr_rsf2)  # NOTE WHETHER EXCLUDING OUTLIER PREDICTIONS HERE
  
  #'  Wrangle results into a table
  spp <- rep(c("Mule Deer", "Elk", "White-tailed Deer", "Cougar", "Wolf", "Bobcat", "Coyote"), each = 2)
  season <- rep(c("Summer", "Winter"), 7)
  corr <- c(md_smr_corr, md_wtr_corr, elk_smr_corr, elk_wtr_corr, wtd_smr_corr, 
            wtd_wtr_corr, coug_smr_corr, coug_wtr_corr, wolf_smr_corr, wolf_wtr_corr,
            bob_smr_corr, bob_wtr_corr, coy_smr_corr, coy_wtr_corr)
  corr <- as.data.frame(corr)
  corr_results <- as.data.frame(cbind(spp, season, corr)) %>%
    transmute(
      Species = spp,
      Season = season,
      #Correlation = as.numeric(corr),
      Correlation = round(corr, digits = 2)
    ) %>%
    arrange(Species)
  
  #'  Save correlations
  write.csv(corr_results, paste0("./Outputs/Tables/Correlation_OccMod_RSF_Predictions_BuffHR_1km_", Sys.Date(), ".csv"))  # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o outliers)
  
  
  ####  Re-scale RSF values between 0 & 1 for mapping  ####
  #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
  Predicted_rsf_rescale <- Predicted_rsf %>%
    transmute(
      obs = obs,
      Area = Area,
      x = x,
      y = y,
      COUG_smr_rsf = round(COUG_smr_rsf2/(max(COUG_smr_rsf2, na.rm = T)), digits = 2),
      COUG_wtr_rsf = round(COUG_wtr_rsf2/(max(COUG_wtr_rsf2, na.rm = T)), digits = 2), 
      WOLF_smr_rsf = round(WOLF_smr_rsf2/(max(WOLF_smr_rsf2, na.rm = T)), digits = 2),
      WOLF_wtr_rsf = round(WOLF_wtr_rsf2/(max(WOLF_wtr_rsf2, na.rm = T)), digits = 2),
      BOB_smr_rsf = round(BOB_smr_rsf2/(max(BOB_smr_rsf2, na.rm = T)), digits = 2),
      BOB_wtr_rsf = round(BOB_wtr_rsf2/(max(BOB_wtr_rsf2, na.rm = T)), digits = 2),
      COY_smr_rsf = round(COY_smr_rsf2/(max(COY_smr_rsf2, na.rm = T)), digits = 2),
      COY_wtr_rsf = round(COY_wtr_rsf2/(max(COY_wtr_rsf2, na.rm = T)), digits = 2),
      ELK_smr_rsf = round(ELK_smr_rsf2/(max(ELK_smr_rsf2, na.rm = T)), digits = 2),
      ELK_wtr_rsf = round(ELK_wtr_rsf2/(max(ELK_wtr_rsf2, na.rm = T)), digits = 2), 
      WTD_smr_rsf = round(WTD_smr_rsf2/(max(WTD_smr_rsf2, na.rm = T)), digits = 2), 
      WTD_wtr_rsf = round(WTD_wtr_rsf2/(max(WTD_wtr_rsf2, na.rm = T)), digits = 2),
      MD_smr_rsf = round(MD_smr_rsf2/(max(MD_smr_rsf2, na.rm = T)), digits = 2), 
      MD_wtr_rsf = round(MD_wtr_rsf2/(max(MD_wtr_rsf2, na.rm = T)), digits = 2)  
    )
  
  #'  Save
  write.csv(Predicted_rsf_rescale, paste0("./Outputs/Tables/Predicted_Relative_Selection_rescale_1km_", Sys.Date(), ".csv"))
  
  
  #'  For some reason ggplot is freaking out over plotting actual 0 values and
  #'  represents them as NA in maps. So for plotting purposes only I am changing
  #'  all pixels with value 0 to 0.0001 so they are slightly >0 and will plot.
  Predicted_rsf_rescale <- Predicted_rsf_rescale %>%
    mutate(COUG_smr_rsf = ifelse(COUG_smr_rsf == 0, 0.00001, COUG_smr_rsf),
           WOLF_smr_rsf = ifelse(WOLF_smr_rsf == 0, 0.00001, WOLF_smr_rsf),
           BOB_smr_rsf = ifelse(BOB_smr_rsf == 0, 0.00001, BOB_smr_rsf),
           BOB_wtr_rsf = ifelse(BOB_wtr_rsf == 0, 0.00001, BOB_wtr_rsf),
           ELK_smr_rsf = ifelse(ELK_smr_rsf == 0, 0.00001, ELK_smr_rsf),
           WTD_smr_rsf = ifelse(WTD_smr_rsf == 0, 0.00001, WTD_smr_rsf),
           MD_wtr_rsf = ifelse(MD_wtr_rsf == 0, 0.00001, MD_wtr_rsf))
  
  
  ####  Plot predicted estimates  ####
  #'  Keep in mind the occupancy and RSF results are on very different scales
  #'  so the coloration is going to differ just because of that.
  #'  Is there a way to weight the RSF results so they higher selected areas show
  #'  up better?
  
  ####  MULE DEER  ####
  #'  Summer Occ
  md_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + #low = "azure" #low = "floralwhite"
    #' #'  Define number of ticks per axis
    #' scale_x_discrete(breaks = 3) +
    #' scale_y_discrete(breaks = 4) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer occupancy model") 
  #'  Summer RSF
  md_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill = cut(MD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Summer resource selection function")
  #'  Winter Occ
  md_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter occupancy model") 
  #'  Winter RSF
  md_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill = cut(MD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  md_smr_map <- md_smr_occ_fig + plot_annotation(title = "Predicted Summer Mule Deer Space Use") + md_smr_rsf_fig 
  md_wtr_map <- md_wtr_occ_fig + plot_annotation(title = "Predicted Winter Mule Deer Space Use") + md_wtr_rsf_fig 
  md_predicted_map <- md_smr_occ_fig + plot_annotation(title = "Predicted Mule Deer Space Use") + 
    md_wtr_occ_fig + md_smr_rsf_fig + md_wtr_rsf_fig + plot_layout(ncol = 2)
  
  #'  Visualize
  plot(md_smr_map)
  plot(md_wtr_map)
  plot(md_predicted_map)
  
  #'  How do I thin out x-axis values so they aren't so bunched in panel figure?
  
  
  
  ####  ELK  ####
  #'  Summer Occ
  elk_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  elk_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  elk_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  elk_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  elk_smr_map <- elk_smr_occ_fig + plot_annotation(title = "Predicted Summer Elk Space Use") + elk_smr_rsf_fig 
  elk_wtr_map <- elk_wtr_occ_fig + plot_annotation(title = "Predicted Winter Elk Space Use") + elk_wtr_rsf_fig 
  elk_predicted_map <- elk_smr_occ_fig + plot_annotation(title = "Predicted Elk Space Use") + 
    elk_wtr_occ_fig + elk_smr_rsf_fig + elk_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(elk_smr_map)
  plot(elk_wtr_map)
  plot(elk_predicted_map)
  
  
  ####  WHITE-TAILED DEER  ####
  #'  Summer Occ
  wtd_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wtd_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  wtd_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wtd_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  wtd_smr_map <- wtd_smr_occ_fig + plot_annotation(title = "Predicted Summer White-tailed Deer Space Use") + wtd_smr_rsf_fig 
  wtd_wtr_map <- wtd_wtr_occ_fig + plot_annotation(title = "Predicted Winter White-tailed Deer Space Use") + wtd_wtr_rsf_fig 
  wtd_predicted_map <- wtd_smr_occ_fig + plot_annotation(title = "Predicted White-tailed Deer Space Use") + 
    wtd_wtr_occ_fig + wtd_smr_rsf_fig +  wtd_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(wtd_smr_map)
  plot(wtd_wtr_map)
  plot(wtd_predicted_map)
  
  
  ####  COUGAR  ####
  #'  Summer Occ
  coug_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coug_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  coug_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coug_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  coug_smr_map <- coug_smr_occ_fig + plot_annotation(title = "Predicted Summer Cougar Space Use") + coug_smr_rsf_fig + plot_layout(ncol = 1)
  coug_wtr_map <- coug_wtr_occ_fig + plot_annotation(title = "Predicted Winter Cougar Space Use") + coug_wtr_rsf_fig + plot_layout(ncol = 1)
  coug_predicted_map <- coug_smr_occ_fig + plot_annotation(title = "Predicted Cougar Space Use") + 
    coug_wtr_occ_fig + coug_smr_rsf_fig + coug_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(coug_smr_map)
  plot(coug_wtr_map)
  plot(coug_predicted_map)
  
  
  ####  WOLF  ####
  #'  Summer Occ
  wolf_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wolf_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  wolf_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wolf_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  wolf_smr_map <- wolf_smr_occ_fig + plot_annotation(title = "Predicted Summer Wolf Space Use") + wolf_smr_rsf_fig + plot_layout(ncol = 1)
  wolf_wtr_map <- wolf_wtr_occ_fig + plot_annotation(title = "Predicted Winter Wolf Space Use") + wolf_wtr_rsf_fig + plot_layout(ncol = 1)
  wolf_predicted_map <- wolf_smr_occ_fig + plot_annotation(title = "Predicted Wolf Space Use") + 
    wolf_wtr_occ_fig + wolf_smr_rsf_fig + wolf_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(wolf_smr_map)
  plot(wolf_wtr_map)
  plot(wolf_predicted_map)
  
  
  ####  BOBCAT  ####
  #'  Summer Occ
  bob_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  bob_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  bob_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  bob_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  bob_smr_map <- bob_smr_occ_fig + plot_annotation(title = "Predicted Summer Bobcat Space Use") + bob_smr_rsf_fig + plot_layout(ncol = 1)
  bob_wtr_map <- bob_wtr_occ_fig + plot_annotation(title = "Predicted Winter Bobcat Space Use") + bob_wtr_rsf_fig + plot_layout(ncol = 1)
  bob_predicted_map <- bob_smr_occ_fig + plot_annotation(title = "Predicted Bobcat Space Use") + 
    bob_wtr_occ_fig + bob_smr_rsf_fig + bob_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(bob_smr_map)
  plot(bob_wtr_map)
  plot(bob_predicted_map)
  
  
  ####  COYOTE  ####
  #'  Summer Occ
  coy_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coy_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  coy_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coy_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  coy_smr_map <- coy_smr_occ_fig + plot_annotation(title = "Predicted Summer Coyote Space Use") + coy_smr_rsf_fig + plot_layout(ncol = 1)
  coy_wtr_map <- coy_wtr_occ_fig + plot_annotation(title = "Predicted Winter Coyote Space Use") + coy_wtr_rsf_fig + plot_layout(ncol = 1)
  coy_predicted_map <- coy_smr_occ_fig + plot_annotation(title = "Predicted Coyote Space Use") + 
    coy_wtr_occ_fig + coy_smr_rsf_fig + coy_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(coy_smr_map)
  plot(coy_wtr_map)
  plot(coy_predicted_map)
  
  
  #'  Save figures as PNG images
  ggsave("./Outputs/Figures/Maps/MuleDeer_predict_smr_plot.png", md_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/MuleDeer_predict_wtr_plot.png", md_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/MuleDeer_predicted_plot.png", md_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Elk_predict_smr_plot.png", elk_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/Elk_predict_wtr_plot.png", elk_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/Elk_predicted_plot.png", elk_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/WTDeer_predict_smr_plot.png", wtd_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/WTDeer_predict_wtr_plot.png", wtd_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/WTDeer_predicted_plot.png", wtd_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Cougar_predict_smr_plot.png", coug_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Cougar_predict_wtr_plot.png", coug_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Cougar_predicted_plot.png", coug_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Wolf_predict_smr_plot.png", wolf_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Wolf_predict_wtr_plot.png", wolf_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Wolf_predicted_plot.png", wolf_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Bobcat_predict_smr_plot.png", bob_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Bobcat_predict_wtr_plot.png", bob_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Bobcat_predicted_plot.png", bob_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Coyote_predict_smr_plot.png", coy_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Coyote_predict_wtr_plot.png", coy_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Coyote_predicted_plot.png", coy_predicted_map, width = 18.5, units = "in")
  
  
  ####  Focal Figures for Publication  ####
  ####  COYOTE  ####
  #'  Summer Occ
  coy_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Coyote summer occupancy model") #+
    # theme(text = element_text(size = 8),
    #       plot.title = element_text(size = 10)) +
    # theme(legend.position = "none")
  #'  Summer RSF
  coy_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Coyote summer RSF") #+
    # theme(text = element_text(size = 8),
    #       plot.title = element_text(size = 10)) +
    
  #'  patchwork figures together:
  coy_smr_map <- coy_smr_occ_fig + coy_smr_rsf_fig #+ plot_layout(ncol = 1)
  
  #'  Winter Occ
  coy_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Coyote winter occupancy model") #+
    # theme(text = element_text(size = 8),
    #       plot.title = element_text(size = 10)) 
  #'  Winter RSF
  coy_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Coyote winter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  coy_wtr_map <- coy_wtr_occ_fig / coy_wtr_rsf_fig
  
  ####  WHITE-TAILED DEER  ####
  #'  Summer Occ
  wtd_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    #theme(axis.text.x = element_text(size = 7)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    # theme(text = element_text(size = 14),
    #       plot.title = element_text(size = 14)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("White-tailed deer summer occupancy model")  #White-tailed deer \nsummer occupancy model
  #'  Summer RSF
  wtd_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    #theme(axis.text.x = element_text(size = 7)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    # theme(text = element_text(size = 14),
    #       plot.title = element_text(size = 14)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("White-tailed deer summer RSF") #"White-tailed deer \nsummer RSF"
  #'  patchwork figures together:
  wtd_smr_map <- wtd_smr_occ_fig + wtd_smr_rsf_fig + #guide_area() + 
    plot_layout(guides = 'collect') & theme(legend.box = 'horizontal')
  
  #'  Winter Occ
  wtd_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("White-tailed deer \nwinter occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Winter RSF
  wtd_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("White-tailed deer \nwinter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) 
  #'  patchwork figures together:
  wtd_wtr_map <- wtd_wtr_occ_fig + wtd_wtr_rsf_fig #+ 
  #plot_layout(guides = 'collect') & theme(legend.box = 'horizontal')
  
  ####  MULE DEER  ####
  #'  Summer Occ
  md_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Mule deer \nsummer occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Summer RSF
  md_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill = cut(MD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Mule deer summer RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  md_smr_map <- md_smr_occ_fig + md_smr_rsf_fig #+ 
    #plot_layout(guides = 'collect') & theme(legend.box = 'horizontal')
  
  #'  Winter Occ
  md_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Mule deer winter occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) #+
    #theme(legend.position = "none")
  #'  Winter RSF
  md_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill = cut(MD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Mule deer winter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) #+
    #theme(legend.position = "none")
  #'  patchwork figures together:
  md_wtr_map <- md_wtr_occ_fig + md_wtr_rsf_fig #+ 
  #plot_layout(guides = 'collect') & theme(legend.box = 'horizontal')
  
  ####  COUGAR  ####
  #'  Summer Occ
  coug_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Cougar summer occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Summer RSF
  coug_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Cougar summer RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  coug_smr_map <- coug_smr_occ_fig / coug_smr_rsf_fig
  
  #'  Winter Occ
  coug_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("") + 
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Cougar winter occupancy model") 
  #'  Winter RSF
  coug_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("Longitude") + 
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Cougar winter RSF") 
  #'  patchwork figures together:
  coug_wtr_map <- coug_wtr_occ_fig / coug_wtr_rsf_fig #+ plot_layout(ncol = 1)
  
  ####  WOLF  ####
  #'  Summer Occ
  wolf_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Wolf summer occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Summer RSF
  wolf_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Wolf summer RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  wolf_smr_map <- wolf_smr_occ_fig + wolf_smr_rsf_fig
  
  #'  Winter Occ
  wolf_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("") + 
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Wolf winter occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Winter RSF
  wolf_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("Longitude") + 
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Wolf winter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  wolf_wtr_map <- wolf_wtr_occ_fig / wolf_wtr_rsf_fig
  
  ####  BOBCAT  ####
  #'  Summer Occ
  bob_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Bobcat summer occupancy model") 
  #'  Summer RSF
  bob_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Bobcat summer RSF") 
  #'  patchwork figures together:
  bob_smr_map <- bob_smr_occ_fig + bob_smr_rsf_fig #+ plot_layout(ncol = 1)
  
  #'  Winter Occ
  bob_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(legend.position="none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("") + 
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Bobcat winter occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Winter RSF
  bob_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20", size = 1) +
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Latitude") + xlab("Longitude") + 
    theme(axis.text.x = element_text(size = 7)) +
    theme(legend.position="none") +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Bobcat winter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  bob_wtr_map <- bob_wtr_occ_fig / bob_wtr_rsf_fig
  
  
  ####  ELK  ####
  #'  Summer Occ
  elk_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("Longitude") + xlab("") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Elk summer occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Summer RSF
  elk_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Elk summer RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  elk_smr_map <- elk_smr_occ_fig + elk_smr_rsf_fig
  
  #'  Winter Occ
  elk_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Probability \nof site use')  +
    ggtitle("Elk winter occupancy model") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  Winter RSF
  elk_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F, drop = FALSE) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20", size = 1) +
    #'  Get rid of lines and gray background
    theme_bw(base_size = 10) +
    theme(panel.border = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
    #'  Change legend, axis, & main titles
    ylab("") + xlab("Longitude") + 
    #theme(axis.text.x = element_text(size = 7)) +
    labs(fill = 'Re-scaled relative \nprobability of selection')  +
    ggtitle("Elk winter RSF") #+
    # theme(text = element_text(size = 16),
    #       plot.title = element_text(size = 16))
  #'  patchwork figures together:
  elk_wtr_map <- elk_wtr_occ_fig + elk_wtr_rsf_fig
  
  
  ####  Panels for publication  ####
  p1 <- coy_smr_occ_fig / coy_smr_rsf_fig
  p2 <- wtd_smr_occ_fig / wtd_smr_rsf_fig
  mismatch_maps <- plot_grid(p1, p2, labels = c('a', 'b'), rel_widths = c(2, 1.5)) +
    plot_annotation(title = "Inconsistent predicted space use")
  ggsave("./Outputs/Figures/Maps/Mismatch_Figure5.tiff", mismatch_maps, width = 10, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p3 <- bob_smr_occ_fig / bob_smr_rsf_fig
  p4 <- md_wtr_occ_fig / md_wtr_rsf_fig
  goodmatch_maps <- plot_grid(p3, p4, labels = c('a', 'b'), rel_widths = c(2, 1.25)) +
    plot_annotation(title = "Consistent predicted space use")
  ggsave("./Outputs/Figures/Maps/Match_Figure4.tiff", goodmatch_maps, width = 10, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p5 <- coug_smr_occ_fig / coug_smr_rsf_fig
  p6 <- coug_wtr_occ_fig / coug_wtr_rsf_fig
  appendix_map1 <- plot_grid(p5, p6, labels = 'a', rel_widths = c(1.05, 1)) +
    plot_annotation(title = "Seasonal predicted space use")
  ggsave("./Outputs/Figures/Maps/Appendix_map1.tiff", appendix_map1, width = 12, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p7 <- coy_wtr_occ_fig / coy_wtr_rsf_fig
  appendix_map2 <- plot_grid(p7, labels = 'b')
  ggsave("./Outputs/Figures/Maps/Appendix_map2.tiff", appendix_map2, width = 7, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p8 <- wolf_smr_occ_fig / wolf_smr_rsf_fig
  p9 <- wtd_wtr_occ_fig / wtd_wtr_rsf_fig
  appendix_map3 <- plot_grid(p8, p9, labels = c('c', 'd'), rel_widths = c(2, 1.25))
  ggsave("./Outputs/Figures/Maps/Appendix_map3.tiff", appendix_map3, width = 10, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p10 <- elk_smr_occ_fig / elk_smr_rsf_fig
  p11 <- elk_wtr_occ_fig / elk_wtr_rsf_fig
  appendix_map4 <- plot_grid(p10, p11, labels = 'a', rel_widths = c(1.1, 1)) +
    plot_annotation(title = "Seasonal predicted space use")
  ggsave("./Outputs/Figures/Maps/Appendix_map4.tiff", appendix_map4, width = 8, height = 5, dpi = 800, units = "in", device = 'tiff')
  
  p12 <- bob_wtr_occ_fig / bob_wtr_rsf_fig
  p13 <- wolf_wtr_occ_fig / wolf_wtr_rsf_fig
  appendix_map5 <- plot_grid(p12, p13, labels = c('b', 'c'), rel_widths = c(1, 1.2)) 
  ggsave("./Outputs/Figures/Maps/Appendix_map5.tiff", appendix_map5, width = 12, height = 5, dpi = 800, units = "in", device = 'tiff')
  

  ####  OLD VERSIONS  ####
  #'  ----------------
  
  #'  Example mismatches Option 1
  mismatch_fig <- coy_smr_map / wtd_smr_map + plot_layout(heights = c(1,1)) #+
    #plot_layout(guides = 'collect') & theme(legend.box = 'horizontal')
  mismatch_fig + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 16))
  mismatch_fig[[1]] <- mismatch_fig[[1]] + plot_layout(tag_level = "new")
  mismatch_fig[[2]] <- mismatch_fig[[2]] + plot_layout(tag_level = "new")
  mismatch_fig <- mismatch_fig + plot_annotation(tag_levels = c("a", "1")) + 
    plot_annotation(title = "Inconsistent predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  mismatch_fig
  
  #'  Example mismatches Option 2
  coy_smr_map2 <- coy_smr_occ_fig / coy_smr_rsf_fig 
  wtd_smr_map2 <- wtd_smr_occ_fig / wtd_smr_rsf_fig #+ theme(axis.text.x = element_text(size = 16))
  
  mismatch_fig2 <- (coy_smr_occ_fig / coy_smr_rsf_fig) | (wtd_smr_occ_fig / wtd_smr_rsf_fig) + 
    plot_layout(widths = c(4,1)) + plot_layout(guides = 'collect')
  mismatch_fig2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 16))
  mismatch_fig2[[1]] <- mismatch_fig2[[1]] + plot_layout(tag_level = "new")
  mismatch_fig2[[2]] <- mismatch_fig2[[2]] + plot_layout(tag_level = "new")
  mismatch_fig2 <- mismatch_fig2 + plot_annotation(tag_levels = c("a", "1")) + 
    plot_annotation(title = "Inconsistent predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  mismatch_fig2
  
  #'  Example of consistent maps
  goodmatch_fig <- bob_smr_map / coug_wtr_map #+ plot_layout(heights = c(1,1))
  goodmatch_fig + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 16))
  goodmatch_fig[[1]] <- goodmatch_fig[[1]] + plot_layout(tag_level = "new")
  goodmatch_fig[[2]] <- goodmatch_fig[[2]] + plot_layout(tag_level = "new")
  goodmatch_fig <- goodmatch_fig + plot_annotation(tag_levels = c("a", "1")) +
    plot_annotation(title = "Consistent predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  goodmatch_fig
  
  goodmatch_fig2 <- (bob_smr_occ_fig / bob_smr_rsf_fig) | (md_wtr_occ_fig / md_wtr_rsf_fig) + 
    plot_layout(widths = c(4,1)) + plot_layout(guides = 'collect')
  goodmatch_fig2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 16))
  goodmatch_fig2[[1]] <- goodmatch_fig2[[1]] + plot_layout(tag_level = "new")
  goodmatch_fig2[[2]] <- goodmatch_fig2[[2]] + plot_layout(tag_level = "new")
  goodmatch_fig2 <- goodmatch_fig2 + plot_annotation(tag_levels = c("a", "1")) +
    plot_annotation(title = "Consistent predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  goodmatch_fig2
  
  #' #' Additional plots for supplemental materials
  #' extramaps1_fig <- coug_smr_map / coy_wtr_map / wolf_smr_map
  #' extramaps1_fig + plot_annotation(tag_levels = "A")
  #' extramaps1_fig[[1]] <- extramaps1_fig[[1]] + plot_layout(tag_level = "new")
  #' extramaps1_fig[[2]] <- extramaps1_fig[[2]] + plot_layout(tag_level = "new")
  #' extramaps1_fig[[3]] <- extramaps1_fig[[3]] + plot_layout(tag_level = "new")
  #' extramaps1_fig <- extramaps1_fig + plot_annotation(tag_levels = c("A", "1")) +
  #'   plot_annotation(title = "Seasonal Predicted Space Use")
  #' extramaps1_fig
  #' #'  Keep the tag levels consistent with those above by providing list of custom tags
  #' extramaps2_fig <- md_smr_map / md_wtr_map / wtd_wtr_map #+ plot_layout(guides = 'collect')
  #' extramaps2_fig[[1]] <- extramaps2_fig[[1]] + plot_layout(tag_level = "new")
  #' extramaps2_fig[[2]] <- extramaps2_fig[[2]] + plot_layout(tag_level = "new")
  #' extramaps2_fig[[3]] <- extramaps2_fig[[3]] + plot_layout(tag_level = "new")
  #' extramaps2_fig <- extramaps2_fig + plot_annotation(tag_levels = list(c('D', 'E', 'F'), '1')) #+, 'F'
  #'   # plot_annotation(title = "Mule Deer and White-tailed Deer Predicted Space Use")
  #' extramaps2_fig
  
  #' Additional plots for supplemental materials
  extramaps1_fig <- coug_smr_map | coug_wtr_map 
  extramaps1_fig + plot_layout(heights = c(2, 2)) + plot_annotation(tag_levels = "a")
  extramaps1_fig[[1]] <- extramaps1_fig[[1]] + plot_layout(tag_level = "new")
  extramaps1_fig[[2]] <- extramaps1_fig[[2]] + plot_layout(tag_level = "new")
  extramaps1_fig <- extramaps1_fig + plot_annotation(tag_levels = c("a", "1")) +
    plot_annotation(title = "Seasonal predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  extramaps1_fig
  #'  Keep the tag levels consistent with those above by providing list of custom tags
  extramaps2_fig <- coy_wtr_occ_fig / coy_wtr_rsf_fig #coy_wtr_map 
  extramaps2_fig + plot_annotation(tag_levels = "a")
  # extramaps2_fig[[1]] <- extramaps2_fig[[1]] + plot_layout(tag_level = "new")
  # extramaps2_fig[[2]] <- extramaps2_fig[[2]] + plot_layout(tag_level = "new")
  extramaps2_fig <- extramaps2_fig + plot_annotation(tag_levels = list(c('c')))
  extramaps2_fig
  #'  Keep the tag levels consistent with those above by providing list of custom tags
  extramaps3_fig <- wtd_wtr_map / wolf_smr_map + plot_layout(heights = c(1,2))
  extramaps3_fig[[1]] <- extramaps3_fig[[1]] + plot_layout(tag_level = "new")
  extramaps3_fig[[2]] <- extramaps3_fig[[2]] + plot_layout(tag_level = "new")
  extramaps3_fig <- extramaps3_fig + plot_annotation(tag_levels = list(c('d', 'e'), '1'))
  extramaps3_fig

  extramaps3_figa <- (wolf_smr_occ_fig / wolf_smr_rsf_fig) | (wtd_wtr_occ_fig / wtd_wtr_rsf_fig) + 
    plot_layout(widths = c(4,1)) + plot_layout(guides = 'collect')
  extramaps3_figa + plot_annotation(tag_levels = "a")
  extramaps3_figa[[1]] <- extramaps3_figa[[1]] + plot_layout(tag_level = "new")
  extramaps3_figa[[2]] <- extramaps3_figa[[2]] + plot_layout(tag_level = "new")
  extramaps3_figa <- extramaps3_figa + plot_annotation(tag_levels = list(c('d', 'e'), '1'))
  extramaps3_figa
  
  extramaps4_fig <- (elk_smr_occ_fig / elk_smr_rsf_fig) | (elk_wtr_occ_fig / elk_wtr_rsf_fig) + 
    plot_layout(widths = c(4,1)) + plot_layout(guides = 'collect')
  extramaps4_fig + plot_annotation(tag_levels = "a")
  extramaps4_fig[[1]] <- extramaps4_fig[[1]] + plot_layout(tag_level = "new")
  extramaps4_fig[[2]] <- extramaps4_fig[[2]] + plot_layout(tag_level = "new")
  extramaps4_fig <- extramaps4_fig + plot_annotation(tag_levels = c("a", "1")) +
    plot_annotation(title = "Seasonal predicted space use",
                    theme = theme(plot.title = element_text(size = 18)))
  extramaps4_fig
  
  extramaps5_fig <- bob_wtr_map | wolf_wtr_map + plot_layout(guides = 'collect')
  extramaps5_fig + plot_annotation(tag_levels = "a")
  extramaps5_fig[[1]] <- extramaps5_fig[[1]] + plot_layout(tag_level = "new")
  extramaps5_fig[[2]] <- extramaps5_fig[[2]] + plot_layout(tag_level = "new")
  extramaps5_fig <- extramaps5_fig + plot_annotation(tag_levels = list(c('c', 'd'), 1))
  extramaps5_fig
  
  # extramaps4_5_fig <- ((elk_smr_occ_fig / elk_smr_rsf_fig) | (elk_wtr_occ_fig / elk_wtr_rsf_fig)) /
  #   bob_wtr_map / wolf_wtr_map + plot_layout(widths = c(3,2)) + plot_layout(guides = 'collect') +
  #   plot_layout(heights = c(1, 2, 2))
  # extramaps4_5_fig + plot_annotation(tag_levels = "a")
  # extramaps4_5_fig[[1]] <- extramaps4_5_fig[[1]] + plot_layout(tag_level = "new")
  # extramaps4_5_fig[[2]] <- extramaps4_5_fig[[2]] + plot_layout(tag_level = "new")
  # extramaps4_5_fig[[3]] <- extramaps4_5_fig[[3]] + plot_layout(tag_level = "new")
  # extramaps4_5_fig <- extramaps4_5_fig + plot_annotation(tag_levels = c("a", "1")) +
  #   plot_annotation(title = "Seasonal predicted space use")
  # extramaps4_5_fig
  
  #' #' Lame maps for supplemental materials
  #' rsfonly_figs <- (elk_smr_rsf_fig + elk_wtr_rsf_fig) / bob_wtr_rsf_fig / wolf_wtr_rsf_fig  + 
  #'   plot_layout(guides = 'collect') + plot_layout(heights = c(2, 2, 2))  
  #' rsfonly_figs + plot_annotation(tag_levels = "a")
  #' rsfonly_figs[[1]] <- rsfonly_figs[[1]] + plot_layout(tag_level = "new")
  #' # rsfonly_figs[[2]] <- rsfonly_figs[[2]] + plot_layout(tag_level = "new")
  #' # rsfonly_figs[[3]] <- rsfonly_figs[[3]] + plot_layout(tag_level = "new")
  #' rsfonly_figs <- rsfonly_figs + plot_annotation(tag_levels = c("a", "1")) +
  #'   plot_annotation(title = "Predicted resource selection functions")#, \nIncomparable to Corresponding Occupancy Models
  #' rsfonly_figs
  
  #'  Save 'em
  ggsave("./Outputs/Figures/Maps/Mismatch_figa.tiff", mismatch_fig, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Mismatch_fig4a.tiff", mismatch_fig2, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Match_figa.tiff", goodmatch_fig, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Match_fig5a.tiff", goodmatch_fig2, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra1_fig.tiff", extramaps1_fig, width = 12, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra2_fig.tiff", extramaps2_fig, width = 11, height = 11, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra1_figa.tiff", extramaps1_fig, width = 12, height = 5, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra2_figa.tiff", extramaps2_fig, width = 11, height = 9, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra3_fig.tiff", extramaps3_fig, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra3_figa.tiff", extramaps3_figa, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra4_fig.tiff", extramaps4_fig, width = 12, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures/Maps/Extra5_fig.tiff", extramaps5_fig, width = 15, height = 7, dpi = 800, units = "in", device = 'tiff')
  #ggsave("./Outputs/Figures/Maps/Extra4-5_fig.tiff", extramaps4_5_fig, width = 8, height = 11, dpi = 800, units = "in", device = 'tiff')
  
  
  
  
  