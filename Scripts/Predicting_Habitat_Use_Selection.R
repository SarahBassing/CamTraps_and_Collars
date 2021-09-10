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
  
  projection(WA)
  projection(OK_SA)
  extent(OK_SA)
  extent(NE_SA)

  
  ####  Load Model Results  ####
  #'  Occupancy model output
  occ_out <- read.csv("./Outputs/Tables/OccMod_OccProb_Results_noHM_2021-08-15.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE
    #'  Calculate 90% confidence intervals to mirror alpha = 0.1
    mutate(
      l95 = (Estimate - (1.64 * SE)),  #### REMINDER: this is 90% CI even though column says l95/u95
      u95 = (Estimate + (1.64 * SE))   
    ) %>%
    dplyr::select(-c(X, Model))
  
  #'  RSF results output
  rsf_out <- read.csv("./Outputs/Tables/RSF_Results_noHM_2021-08-15.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE
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
    arrange(Year) #NECESSARY TO MATCH DH's CAMERALOCATION ORDER 2021-08-10 version
  
  #'  Read in original covariate data from RSFs
  load("./Outputs/RSF_pts/md_dat_2nd_all_2021-08-10.RData")  
  load("./Outputs/RSF_pts/elk_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/wtd_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/coug_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/wolf_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/bob_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/coy_dat_2nd_all_2021-08-10.RData")
  
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
  
  #'  Read in & scale covariate data for entire study areas (based on 1km grids)
  #'  Data generated in "Covariate_Extraction.R" script in WPPP_CameraTrapping Project
  NE_covs_1km <- read.csv("./Outputs/Tables/StudyAreaWide_NE_Covariates_2021-08-10.csv") %>%
    dplyr::select(-X) %>%
    mutate(Area = 0)
  OK_covs_1km <- read.csv("./Outputs/Tables/StudyAreaWide_OK_Covariates_2021-08-10.csv") %>%
  dplyr::select(-X) %>%
    mutate(Area = 1)
  all_covs_1km <- as.data.frame(rbind(NE_covs_1km, OK_covs_1km))
  
  #'  Exclude elevations >2150 for camera covariates since didn't sample above this elevation
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
  cam_zcovs_1km <- scaling_covs(all_covs_adj_1km, summary_occ_covs)
  #'  Standardize covariate data based on collar covariate means & SDs for specific species & seasons
  md_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_md_smr)
  md_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_md_wtr)
  elk_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_elk_smr)
  elk_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_elk_wtr)
  wtd_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_wtd_smr)
  wtd_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_wtd_wtr)
  coug_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_coug_smr)
  coug_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_coug_wtr)
  wolf_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_wolf_smr)
  wolf_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_wolf_wtr)
  bob_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_bob_smr)
  bob_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_bob_wtr)
  coy_smr_zcovs_1km <- scaling_covs(all_covs_1km, summary_coy_smr)
  coy_wtr_zcovs_1km <- scaling_covs(all_covs_1km, summary_coy_wtr)
  
  #'  Function to covert all covariate data into SpatialPointsDataFrames
  sp_covs <- function(covs) {
    xy <- covs[,c(10,11)]
    covs_spdf <- SpatialPointsDataFrame(data = covs, coords = xy,
                                        proj4string = CRS(sa_proj))
  }
  #'  Gather all scaled covariate data into a monster list
  zcovs <- list(cam_zcovs_1km, md_smr_zcovs_1km, md_wtr_zcovs_1km, 
                elk_smr_zcovs_1km, elk_wtr_zcovs_1km, wtd_smr_zcovs_1km, 
                wtd_wtr_zcovs_1km, coug_smr_zcovs_1km, coug_wtr_zcovs_1km, 
                wolf_smr_zcovs_1km, wolf_wtr_zcovs_1km, bob_smr_zcovs_1km, 
                bob_wtr_zcovs_1km, coy_smr_zcovs_1km, coy_wtr_zcovs_1km)
  #'  Run list of scaled covariates through spatial function
  sp_zcovs <- lapply(zcovs, sp_covs)
  
  

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
    #'  Use p-values to change non-significant coefficients (alpha-level = 0.1) to 0 so there is no effect
    mutate(Estimate = ifelse(Pval > 0.1, Estimate == 0, Estimate)) %>%
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
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
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
  md_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Summer",])
  md_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 1,], occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Winter",])
  elk_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Summer",])
  elk_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Winter",])
  wtd_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]][sp_zcovs[[1]]$Area == 0,], occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Winter",])
  coug_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Summer",])
  coug_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Winter",])
  wolf_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Winter",])
  bob_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Summer",])
  bob_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Winter",])
  coy_smr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Summer",])
  coy_wtr_predict_occ_sgnf <- predict_occ(sp_zcovs[[1]], occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predator data that spans both study areas
  Predicted_occ <- as.data.frame(all_covs_1km) %>%
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
  OK_rows <- seq(1:nrow(md_smr_predict_occ_sgnf))
  Area <- rep("Okanogan", length(OK_rows))
  OK_occ <- as.data.frame(cbind(OK_rows, Area, md_smr_predict_occ_sgnf, md_wtr_predict_occ_sgnf)) 
  colnames(OK_occ) <- c("obs", "Area", "MD_smr_occ", "MD_wtr_occ")
  
  #'  Northeast-only data (elk & white-tailed deer)
  NE_rows <- seq(1:nrow(elk_smr_predict_occ_sgnf))
  Area <- rep("Northeast", length(NE_rows))
  NE_occ <- as.data.frame(cbind(NE_rows, Area, elk_smr_predict_occ_sgnf, elk_wtr_predict_occ_sgnf, 
                                wtd_smr_predict_occ_sgnf, wtd_wtr_predict_occ_sgnf)) 
  colnames(NE_occ) <- c("obs", "Area", "ELK_smr_occ", "ELK_wtr_occ", "WTD_smr_occ", "WTD_wtr_occ")
  
  #'  Merge ungulate & predator data by study area
  Predicted_occ_OK <- Predicted_occ[Predicted_occ$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_occ = NA,
      ELK_wtr_occ = NA,
      WTD_smr_occ = NA,
      WTD_wtr_occ = NA
    ) %>%
    full_join(OK_occ, by = c("obs", "Area")) 
  Predicted_occ_NE <- Predicted_occ[Predicted_occ$Area == "Northeast",] %>%
    full_join(NE_occ, by = c("obs", "Area")) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_occ = NA,
      MD_wtr_occ = NA
    )
  
  #'  Merge NE and OK predictions together
  Predicted_occ <- as.data.frame(rbind(Predicted_occ_NE, Predicted_occ_OK))
  
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
    #'  Use p-values to change non-significant coefficients (alpha-level = 0.05) to 0 so there is no effect
    mutate(Estimate = ifelse(Pval > 0.05, Estimate == 0, Estimate)) %>%
    #'  For some reason the mutation above changes estimates that are already 0.00 with Pval > 0.05 to equal 1!?!?!
    #'  So changing those back to 0 since the effect is non-significant
    mutate(Estimate = ifelse(Estimate == 1, Estimate == 0, Estimate)) %>%
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
      B.rd = RoadDen) # B.hm = HumanMod
  
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
  md_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[2]][sp_zcovs[[2]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Summer",])
  md_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[3]][sp_zcovs[[3]]$Area == 1,], rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Winter",])
  elk_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[4]][sp_zcovs[[4]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Summer",])
  elk_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[5]][sp_zcovs[[5]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Winter",])
  wtd_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[6]][sp_zcovs[[6]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[7]][sp_zcovs[[7]]$Area == 0,], rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Winter",])
  coug_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[8]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Summer",])
  coug_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[9]], rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Winter",])
  wolf_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[10]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[11]], rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Winter",])
  bob_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[12]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Summer",])
  bob_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[13]], rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Winter",])
  coy_smr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[14]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Summer",])
  coy_wtr_predict_rsf_sgnf <- predict_rsf(sp_zcovs[[15]], rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predators
  Predicted_rsf <- as.data.frame(all_covs_1km) %>%
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
  OK_rows <- seq(1:nrow(md_smr_predict_rsf_sgnf))
  Area <- rep("Okanogan", length(OK_rows))
  OK_rsf <- as.data.frame(cbind(OK_rows, Area, md_smr_predict_rsf_sgnf, md_wtr_predict_rsf_sgnf)) # KEEP TRACK of which version of the predicted results I'm using
  colnames(OK_rsf) <- c("obs", "Area", "MD_smr_rsf", "MD_wtr_rsf")
  
  #'  Northeast-only data (elk & white-tailed deer)
  NE_rows <- seq(1:nrow(elk_smr_predict_rsf_sgnf))
  Area <- rep("Northeast", length(NE_rows))
  NE_rsf <- as.data.frame(cbind(NE_rows, Area, elk_smr_predict_rsf_sgnf, elk_wtr_predict_rsf_sgnf, # KEEP TRACK of which version of the predicted results I'm using
                                wtd_smr_predict_rsf_sgnf, wtd_wtr_predict_rsf_sgnf)) 
  colnames(NE_rsf) <- c("obs", "Area", "ELK_smr_rsf", "ELK_wtr_rsf", "WTD_smr_rsf", "WTD_wtr_rsf")
  
  #'  Merge ungulate & predator data by study area
  Predicted_rsf_OK <- Predicted_rsf[Predicted_rsf$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_rsf = NA,
      ELK_wtr_rsf = NA,
      WTD_smr_rsf = NA,
      WTD_wtr_rsf = NA
    ) %>%
    full_join(OK_rsf, by = c("obs", "Area")) 
  Predicted_rsf_NE <- Predicted_rsf[Predicted_rsf$Area == "Northeast",] %>%
    full_join(NE_rsf, by = c("obs", "Area")) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_rsf = NA,
      MD_wtr_rsf = NA
    )
  
  #'  Merge NE and OK predictions togther
  Predicted_rsf <- as.data.frame(rbind(Predicted_rsf_NE, Predicted_rsf_OK))
  
  #'  Identify any outliers
  outliers <- function(predicted, title) {
    print(summary(predicted))
    hist(predicted, breaks = 100, main = title)
    boxplot(predicted, main = title)
  }
  #'  Identify outlier perdictions
  outliers(Predicted_rsf$MD_smr_rsf, "Mule Deer Summer RSF Predictions")
  outliers(Predicted_rsf$MD_wtr_rsf, "Mule Deer Winter RSF Prediction")
  outliers(Predicted_rsf$ELK_smr_rsf, "Elk Summer RSF Predictions")
  outliers(Predicted_rsf$ELK_wtr_rsf, "Elk Winter RSF Prediction")
  outliers(Predicted_rsf$WTD_smr_rsf, "White-tailed Deer Summer RSF Predictions")
  outliers(Predicted_rsf$WTD_wtr_rsf, "White-tailed Deer Winter RSF Prediction")
  outliers(Predicted_rsf$COUG_smr_rsf, "Cougar Summer RSF Predictions")
  outliers(Predicted_rsf$COUG_wtr_rsf, "Cougar Winter RSF Prediction")
  outliers(Predicted_rsf$WOLF_smr_rsf, "Wolf Summer RSF Predictions")
  outliers(Predicted_rsf$WOLF_wtr_rsf, "Wolf Winter RSF Prediction")
  outliers(Predicted_rsf$BOB_smr_rsf, "Bobcat Summer RSF Predictions")
  outliers(Predicted_rsf$BOB_wtr_rsf, "Bobcat Winter RSF Prediction")
  outliers(Predicted_rsf$COY_smr_rsf, "Coyote Summer RSF Predictions")
  outliers(Predicted_rsf$COY_wtr_rsf, "Coyote Winter RSF Prediction")
  
  #'  Exclude extreme outliers as identified by histograms & boxplots
  Predicted_rsf <- Predicted_rsf %>%
    mutate(
      MD_smr_rsf2 = ifelse(MD_smr_rsf > 15, NA, MD_smr_rsf),
      MD_wtr_rsf2 = ifelse(MD_wtr_rsf > 300, NA, MD_wtr_rsf),
      ELK_wtr_rsf2 = ifelse(ELK_wtr_rsf > 25, NA, ELK_wtr_rsf),
      WTD_smr_rsf2 = ifelse(WTD_smr_rsf > 10, NA, WTD_smr_rsf),
      COUG_wtr_rsf2 = ifelse(COUG_wtr_rsf > 25, NA, COUG_wtr_rsf)
    )

 
  #'  Merge all predictions together (with unscaled RSF predictions)
  Predicted_Occ_RSF <- Predicted_occ %>%
    full_join(Predicted_rsf, by = c("obs", "Area", "x", "y"))
  write.csv(Predicted_Occ_RSF, paste0("./Outputs/Tables/Predictions_OccMod_v_RSF_noHM_", Sys.Date(), ".csv"))  
  
  
  
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
  md_smr_corr <- predict_corr(Predicted_occ$MD_smr_occ, Predicted_rsf$MD_smr_rsf2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
  md_wtr_corr <- predict_corr(Predicted_occ$MD_wtr_occ, Predicted_rsf$MD_wtr_rsf2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
  elk_smr_corr <- predict_corr(Predicted_occ$ELK_smr_occ, Predicted_rsf$ELK_smr_rsf)  
  elk_wtr_corr <- predict_corr(Predicted_occ$ELK_wtr_occ, Predicted_rsf$ELK_wtr_rsf2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
  wtd_smr_corr <- predict_corr(Predicted_occ$WTD_smr_occ, Predicted_rsf$WTD_smr_rsf2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
  wtd_wtr_corr <- predict_corr(Predicted_occ$WTD_wtr_occ, Predicted_rsf$WTD_wtr_rsf)
  coug_smr_corr <- predict_corr(Predicted_occ$COUG_smr_occ, Predicted_rsf$COUG_smr_rsf)
  coug_wtr_corr <- predict_corr(Predicted_occ$COUG_wtr_occ, Predicted_rsf$COUG_wtr_rsf2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
  wolf_smr_corr <- predict_corr(Predicted_occ$WOLF_smr_occ, Predicted_rsf$WOLF_smr_rsf)
  wolf_wtr_corr <- predict_corr(Predicted_occ$WOLF_wtr_occ, Predicted_rsf$WOLF_wtr_rsf)
  bob_smr_corr <- predict_corr(Predicted_occ$BOB_smr_occ, Predicted_rsf$BOB_smr_rsf)
  bob_wtr_corr <- predict_corr(Predicted_occ$BOB_wtr_occ, Predicted_rsf$BOB_wtr_rsf)
  coy_smr_corr <- predict_corr(Predicted_occ$COY_smr_occ, Predicted_rsf$COY_smr_rsf)
  coy_wtr_corr <- predict_corr(Predicted_occ$COY_wtr_occ, Predicted_rsf$COY_wtr_rsf)
  
  #'  Wrangle results into a table
  spp <- rep(c("Mule Deer", "Elk", "White-tailed Deer", "Cougar", "Wolf", "Bobcat", "Coyote"), each = 2)
  season <- rep(c("Summer", "Winter"), 7)
  corr <- c(md_smr_corr, md_wtr_corr, elk_smr_corr, elk_wtr_corr, wtd_smr_corr, 
            wtd_wtr_corr, coug_smr_corr, coug_wtr_corr, wolf_smr_corr, wolf_wtr_corr,
            bob_smr_corr, bob_wtr_corr, coy_smr_corr, coy_wtr_corr)
  corr_results <- as.data.frame(cbind(spp, season, corr)) %>%
    transmute(
      Species = spp,
      Season = season,
      Correlation = as.numeric(corr),
      Correlation = round(Correlation, digits = 2)
    ) %>%
    arrange(Species)
  
  #'  Save correlations
  write.csv(corr_results, paste0("./Outputs/Tables/Correlation_OccMod_RSF_Predictions_noHM_", Sys.Date(), ".csv"))  # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o non-signif coeffs)
  
  
  ####  Re-scale RSF values between 0 & 1 for mapping  ####
  #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
  Predicted_rsf_rescale <- Predicted_rsf %>%
    transmute(
      obs = obs,
      Area = Area,
      x = x,
      y = y,
      COUG_smr_rsf = round(COUG_smr_rsf/(max(COUG_smr_rsf, na.rm = T)), digits = 2),
      COUG_wtr_rsf = round(COUG_wtr_rsf2/(max(COUG_wtr_rsf2, na.rm = T)), digits = 2), # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
      WOLF_smr_rsf = round(WOLF_smr_rsf/(max(WOLF_smr_rsf, na.rm = T)), digits = 2),
      WOLF_wtr_rsf = round(WOLF_wtr_rsf/(max(WOLF_wtr_rsf, na.rm = T)), digits = 2),
      BOB_smr_rsf = round(BOB_smr_rsf/(max(BOB_smr_rsf, na.rm = T)), digits = 2),
      BOB_wtr_rsf = round(BOB_wtr_rsf/(max(BOB_wtr_rsf, na.rm = T)), digits = 2),
      COY_smr_rsf = round(COY_smr_rsf/(max(COY_smr_rsf, na.rm = T)), digits = 2),
      COY_wtr_rsf = round(COY_wtr_rsf/(max(COY_wtr_rsf, na.rm = T)), digits = 2),
      ELK_smr_rsf = round(ELK_smr_rsf/(max(ELK_smr_rsf, na.rm = T)), digits = 2),
      ELK_wtr_rsf = round(ELK_wtr_rsf2/(max(ELK_wtr_rsf2, na.rm = T)), digits = 2), # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
      WTD_smr_rsf = round(WTD_smr_rsf2/(max(WTD_smr_rsf2, na.rm = T)), digits = 2), # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
      WTD_wtr_rsf = round(WTD_wtr_rsf/(max(WTD_wtr_rsf, na.rm = T)), digits = 2),
      MD_smr_rsf = round(MD_smr_rsf2/(max(MD_smr_rsf2, na.rm = T)), digits = 2), # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
      MD_wtr_rsf = round(MD_wtr_rsf2/(max(MD_wtr_rsf2, na.rm = T)), digits = 2)  # NOTE: EXCLUDING OUTLIER PREDICTIONS HERE
    )
  
  #'  Save
  write.csv(Predicted_rsf_rescale, paste0("./Outputs/Tables/Predicted_Relative_Selection_rescale_", Sys.Date(), ".csv"))
  
  
  
  ####  Plot predicted estimates  ####
  #'  Keep in mind the occupancy and RSF results are on very different scales
  #'  so the coloration is going to differ just because of that.
  #'  Is there a way to weight the RSF results so they higher selected areas show
  #'  up better?
  
  ####  MULE DEER  ####
  #'  Summer Occ
  md_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + #low = "azure" #low = "floralwhite"
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  md_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill = cut(MD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) + 
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  md_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  md_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Okanogan",], aes(x = x, y = y, fill =cut(MD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  elk_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  elk_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(ELK_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wtd_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wtd_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale[Predicted_rsf_rescale$Area == "Northeast",], aes(x = x, y = y, fill = cut(WTD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coug_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coug_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wolf_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wolf_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  bob_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  bob_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coy_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coy_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
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
  
  