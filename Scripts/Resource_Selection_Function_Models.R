  #'  ============================================
  #'  Resource Selection Functions (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  ============================================
  #'  Script to build resource selection function models for comparison to 
  #'  occupancy models and HMMs.
  #'  
  #'  Cleaned telemetry and covariate data were prepared with
  #'  Collar_RSF_DataPrep.R script which take FOREVER to run. 
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load packages for selecting available points
  library(tidyverse)
  library(lme4)
  
  #'  Load used and available locations, and covariate data
  # load("./Outputs/RSF_pts/md_dat_all_DATE.RData")
  load("./Outputs/RSF_pts/elk_dat_all_2021-06-17.RData")
  # load("./Outputs/RSF_pts/wtd_dat_all_DATE.RData")
  load("./Outputs/RSF_pts/coug_dat_all_2021-06-17.RData")
  load("./Outputs/RSF_pts/wolf_dat_all_2021-06-17.RData")
  load("./Outputs/RSF_pts/bob_dat_all_2021-06-17.RData")
  load("./Outputs/RSF_pts/coy_dat_all_2021-06-17.RData")
  
  
  #' #'  Covariates at used locations
  #' load("./Outputs/Telemetry_covs/spp_telem_covs_2021-05-10.RData")
  #' load("./Outputs/RSF_available_pts/coy_avail_covs_2021-05-17.RData")
  #' #'  Used locations
  #' # load("./Outputs/Telemetry_tracks/MD_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/MD_wtr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/ELK_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/ELK_wtr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/WTD_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/WTD_wtr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/COUG_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/COUG_wtr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/WOLF_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/WOLF_wtr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/BOB_smr_track.RData")
  #' # load("./Outputs/Telemetry_tracks/BOB_wtr_track.RData")
  #' load("./Outputs/Telemetry_tracks/COY_smr_track.RData")
  #' load("./Outputs/Telemetry_tracks/COY_wtr_track.RData")
  #' #'  Available locations
  #' load("./Outputs/RSF_available_pts/coy_available_2021-05-12.RData")
  #' 
  #' 
  #' #'  Add use classification to the collar locations
  #' used <- function(used_pts, covs) {
  #'   used_locs <- left_join(used_pts, covs, by = c("time", "ID")) %>%
  #'     mutate(
  #'       use = 1,
  #'       Year = ifelse(Season.x == "Summer18" | Season.x == "Winter1819", "Year1", "Year2")
  #'     ) %>%
  #'     dplyr::select(c(AnimalID, Sex, Season.x, x, y, use, Elev, Slope, HumanMod, #NearestRd, 
  #'                     PercForMix, PercXGrass, PercXShrub, Area, Year))
  #'   names(used_locs)[names(used_locs) == "Season.x"] <- "Season"
  #'   return(used_locs)
  #' }
  #' coy_smr_used <- used(COY_smr_track, spp_telem_covs[[13]])
  #' coy_wtr_used <- used(COY_wtr_track, spp_telem_covs[[14]])
  #' # track_list <- list(COY_smr_track, COY_wtr_track)
  #' # coy_used <- lapply(track_list, used)
  #' 
  #' 
  #' #'  Retain only relevant location data
  #' avail <- function(avail_pts, avail_covs) {
  #'   avail_pts$obs <- as.numeric(seq(1:nrow(avail_pts)))
  #'   avail_locs <- full_join(avail_pts, avail_covs, by = c("obs", "ID", "Season")) %>%
  #'     dplyr::select(c(ID, Sex, Season, x, y, use, Elev, Slope, HumanMod, #NearestRd,
  #'                     PercForMix, PercXGrass, PercXShrub, Area, Year))
  #'   names(avail_locs)[names(avail_locs) == "ID"] <- "AnimalID"
  #'   return(avail_locs)
  #' } 
  #' coy_smr_avail <- avail(coy_available[[1]], coy_avail_covs[[1]])
  #' coy_wtr_avail <- avail(coy_available[[2]], coy_avail_covs[[2]])
  #' 
  #' 
  #' #'  Combine into a single dataframe
  #' coy_locs_smr <- rbind(coy_smr_used, coy_smr_avail) 
  #' coy_locs_wtr <- rbind(coy_wtr_used, coy_wtr_avail)
  
  
  
  #'  Center & scale covariates 
  #'  Note: standardizing across all IDs & years, but separately by season & spp
  spp_dataPrep <- function(locs){
    #'  Make categorical variables factors
    locs$ID <- as.factor(locs$ID)
    locs$Used <- as.factor(locs$Used)
    locs$Area <- as.factor(locs$Area)
    # locs$Sex <- as.factor(locs$Sex)
    locs$Year <- as.factor(locs$Year)
    locs$Season <- as.factor(locs$Season)
    #'  Standardize continuous variables
    locs$Elev <- scale(locs$Elev)
    locs$Slope <- scale(locs$Slope)
    locs$RoadDen <- scale(locs$RoadDen)
    locs$HumanMod <- scale(locs$HumanMod)
    # locs$NearestRd <- scale(locs$NearestRd)
    locs$PercForMix <- scale(locs$PercForMix)
    locs$PercXGrass <- scale(locs$PercXGrass)
    locs$PercXShrub <- scale(locs$PercXShrub)
    
    locs <- as.data.frame(locs)
  
    return(locs)
  }
  #'  Run season & species-specific data through prep function
  # mdData_smr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",])
  # mdData_wtr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",])
  elkData_smr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19",])
  elkData_wtr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920",])
  # wtdData_smr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19",])
  # wtdData_wtr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920",])
  cougData_smr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",])
  cougData_wtr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",])
  wolfData_smr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",])
  wolfData_wtr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",])
  bobData_smr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",])
  bobData_wtr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",])
  coyData_smr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19",])
  coyData_wtr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920",])
  
  
  
  #'  Resource Selection Function Models
  #'  ==================================
  #'  Functions to run logistic mixed effects models that include random effect 
  #'  for individual
  #'  Prey RSFs excluded random effect for study area
  #'  Predator RSFs include random effect for study area
  #'  Other habitat covariates excluded based on species and convergence issues
  
  ####  Mule Deer RSF  ####
  md_rsf <- function(dat) {
    #' #'  Null model
    #' mod0 <- glm(Used ~ 1, data = dat, family = binomial(link = "logit"))
    #' print(summary(mod0))
    
    #'  Random effect for individual and year
    mod_reff <- glmer(Used ~ 1 + (1|ID) + (1|Year), data = dat, family = binomial(link = "logit"))
    print(summary(mod_reff))
    
    #' #'  Semi-global model including random effects
    #' mod_semiglobal <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
    #'                     data = dat, family = binomial(link = "logit"))
    #' print(summary(mod_semiglobal))
    
    #'  Global model with random effect for individual and year
    mod_global <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                        data = dat, family = binomial(link = "logit"))
    print(summary(mod_global))
    
    #'  AIC for model selection
    # print(AIC(mod_reff, mod_semiglobal, mod_global))
    print(AIC(mod_reff, mod_global))
    
    #'  List outputs
    # mod_out <- list(mod_reff, mod_semiglobal, mod_global)
    mod_out <- list(mod_reff, mod_global)
    return(mod_out)
    
  }
  #'  Run RSF on seasonal used and available locations for each species
  md_rsf_smr <- md_rsf(mdData_smr)
  md_rsf_wtr <- md_rsf(mdData_wtr)
  
  
  ####  Elk RSF  ####
  #'  Global model with random effect for individual and year
  #'  Summer
  elk_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                        data = elkData_smr, family = binomial(link = "logit"))
  summary(elk_global_smr)
  #'  Winter
  elk_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                      data = elkData_wtr, family = binomial(link = "logit"))
  summary(elk_global_wtr)

  
  ####  White-tailed Deer RSF  ####
  wtd_rsf <- function(dat) {
    
    #'  Random effect for individual and year
    mod_reff <- glmer(Used ~ 1 + (1|ID) + (1|Year), data = dat, family = binomial(link = "logit"))
    print(summary(mod_reff))
    
    #'  Semi-global model including random effects, EXCLUDING shrub & grass
    mod_semiglobal <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + HumanMod + (1|ID) + (1|Year),
                        data = dat, family = binomial(link = "logit"))
    print(summary(mod_semiglobal))
    
    #'  Global model with random effect for individual and year
    mod_global <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                        data = dat, family = binomial(link = "logit"))
    print(summary(mod_global))
    
    #'  AIC for model selection
    print(AIC(mod_reff, mod_semiglobal, mod_global))
    
    #'  List outputs
    # mod_out <- list(mod_reff, mod_semiglobal, mod_global)
    mod_out <- list(mod_reff, mod_global)
    return(mod_out)
    
  }
  #'  Run RSF on seasonal used and available locations for each species
  wtd_rsf_smr <- wtd_rsf(wtdData_smr)
  wtd_rsf_wtr <- wtd_rsf(wtdData_wtr)
  
  
  ####  Cougar RSF  ####
  #'  Summer
  #'  Random effect for individual and year
  # coug_reff_smr <- glmer(Used ~ 1 + (1|ID) + (1|Year) + (1|Area), data = cougData_smr, family = binomial(link = "logit"))
  # summary(coug_reff_smr)
  #' #'  Semi-global model including random effects, EXCLUDING shrub & grass
  #' coug_semiglobal_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                     data = cougData_smr, family = binomial(link = "logit"))
  #' summary(coug_semiglobal_smr)
  #'  Global model with random effect for individual and year
  coug_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
                          data = cougData_smr, family = binomial(link = "logit"))
  summary(coug_global_smr)
  #'  AIC for model selection
  # AIC(coug_reff_smr, coug_semiglobal_smr, coug_global_smr)
  # AIC(coug_reff_smr, coug_global_smr)
  #### WHY ARE THESE FITTED WITH DIFFERENT NUMBER OF OBSERVATIONS?!?!?!
  
  #'  Winter
  #'  Random effect for individual and year
  # coug_reff_wtr <- glmer(Used ~ 1 + (1|ID) + (1|Year) + (1|Area), data = cougData_wtr, family = binomial(link = "logit"))
  # summary(coug_reff_wtr)
  #' #'  Semi-global model including random effects, EXCLUDING shrub & grass
  #' coug_semiglobal_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                     data = cougData_wtr, family = binomial(link = "logit"))
  #' summary(coug_semiglobal_wtr)
  #'  Global model with random effect for individual and year
  coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
                          data = cougData_wtr, family = binomial(link = "logit"))
  summary(coug_global_wtr)
  #'  AIC for model selection
  # AIC(coug_reff_wtr, coug_semiglobal_wtr, coug_global_wtr)
  # AIC(coug_reff_wtr, coug_global_wtr)
  
  
  
  ####  Wolf RSF  ####
  #'  Summer  
  #' #'  Semi-global model including random effects, EXCLUDING shrub
  #' wolf_semiglobal_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
  #'                     data = wolfData_smr, family = binomial(link = "logit"))
  #' summary(wolf_semiglobal_smr)
  #'  Global model with random effect for individual, year, and study area
  wolf_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                       data = wolfData_smr, family = binomial(link = "logit"))
  summary(wolf_global_smr)
  
  #'  Winter  
  #' #'  Semi-global model including random effects, EXCLUDING shrub
  #' wolf_semiglobal_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
  #'                     data = wolfData_wtr, family = binomial(link = "logit"))
  #' summary(wolf_semiglobal_wtr)
  #'  Global model with random effect for individual, year, and study area
  wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = wolfData_wtr, family = binomial(link = "logit"))
  summary(wolf_global_wtr)
  
  
  ####  Bobcat RSF  ####
  #'  Summer  
  #'  Global model with random effect for individual, year, and study area
  bob_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = bobData_smr, family = binomial(link = "logit"))
  summary(bob_global_smr)
  
  #'  Winter  
  #'  Global model with random effect for individual, year, and study area
  bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = bobData_wtr, family = binomial(link = "logit"))
  summary(bob_global_wtr)
  
  
  ####  Coy RSF  ####
  #'  Summer  
  #'  Global model with random effect for individual, year, and study area
  coy_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = coyData_smr, family = binomial(link = "logit"))
  summary(coy_global_smr)
  
  #'  Winter  
  #'  Global model with random effect for individual, year, and study area
  coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = coyData_wtr, family = binomial(link = "logit"))
  summary(coy_global_wtr)
  
  
  
  
  
  #'  Save
  save(coy_rsf_smr, file = "./Outputs/RSF_coy_smr.RData")
  save(coy_rsf_wtr, file = "./Outputs/RSF_coy_wtr.RData")
  load("./Outputs/RSF_coy_smr.RData")
  load("./Outputs/RSF_coy_wrr.RData")
  
  #'  Pull out results
  coy_rsf_smr
  
  
  
  
  
  
  
  
  