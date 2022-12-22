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
  library(car)
  library(lme4)
  
  #'  Load used and available locations, and covariate data
  #' #'  2nd Order Selection (based on single large MCP)
  #' load("./Outputs/RSF_pts/md_dat_2nd_all_2021-09-13.RData")
  #' load("./Outputs/RSF_pts/elk_dat_2nd_all_2021-09-13.RData")
  #' load("./Outputs/RSF_pts/wtd_dat_2nd_all_2021-09-13.RData")
  #' load("./Outputs/RSF_pts/coug_dat_2nd_all_2021-09-13.RData")
  #' load("./Outputs/RSF_pts/wolf_dat_2nd_all_2021-10-29.RData")
  #' load("./Outputs/RSF_pts/bob_dat_2nd_all_2021-09-13.RData")
  #' load("./Outputs/RSF_pts/coy_dat_2nd_all_2021-09-13.RData")
  #'  2nd Order Selection (based on buffered home ranges) 1:20 ratio
  load("./Outputs/RSF_pts/md_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/elk_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/wtd_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/coug_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/wolf_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/bob_dat_2nd_buffHR_all_2022-04-18.RData")
  load("./Outputs/RSF_pts/coy_dat_2nd_buffHR_all_2022-04-18.RData")
  #' #'  2nd Order Selection (based on buffered home ranges) 1:200 ratio
  #' load("./Outputs/RSF_pts/Subsampled_pts/md_dat_2nd_200avail_all_2022-05-11.RData")  
  #' load("./Outputs/RSF_pts/Subsampled_pts/elk_dat_2nd_200avail_all_2022-05-11.RData")
  #' load("./Outputs/RSF_pts/Subsampled_pts/wtd_dat_2nd_200avail_all_2022-05-11.RData")
  #' load("./Outputs/RSF_pts/Subsampled_pts/coug_dat_2nd_200avail_all_2022-05-10.RData")
  #' load("./Outputs/RSF_pts/Subsampled_pts/wolf_dat_2nd_200avail_all_2022-05-10.RData") 
  #' load("./Outputs/RSF_pts/Subsampled_pts/bob_dat_2nd_200avail_all_2022-05-10.RData")
  #' load("./Outputs/RSF_pts/Subsampled_pts/coy_dat_2nd_200avail_all_2022-05-10.RData")
  
  
  #'  Exclude used/available locations above 2100m
  high_elev_locs <- function(locs, spp, season) {
    #'  How many locations total?
    nlocs <- nrow(locs)
    #'  How many locations actually fall above this range?
    n_hi_elev <- nrow(locs[locs$Elev > 2100,])
    #'  What proportion of the total locations are above 2100m?
    prop_hi_elev <- round(n_hi_elev/nlocs, 2)
    
    #'  How many used locationts total?
    nused <- nrow(locs[locs$Used == 1,])
    #'  How many of those locations are used?
    nused_hi_elev <- nrow(locs[locs$Elev > 2100 & locs$Used == 1,])
    #'  What proportion of the total used locations are above 2100m?
    prop_used_hi_elev <- round(nused_hi_elev/nused, 2)
    
    #'  What's the maximum elevation included in the data set (used or available)
    max_elev <- max(locs$Elev)
    max_used_elev <- max(locs$Elev[locs$Used == 1])
    #'  Combine data and rename
    hi_elev <- cbind(spp, season, nlocs, n_hi_elev, prop_hi_elev, nused, nused_hi_elev, 
                     prop_used_hi_elev, max_elev, max_used_elev)
    colnames(hi_elev) <- c("Species", "Season", "total locs", "total high elev locs", 
                           "proportion high elev", "total used locs", "n used high elev locs", "proportion used high elev", 
                           "maximum elev", "maximum used elev")
    
    return(hi_elev)
  }
  md_smr_elev <- high_elev_locs(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",], spp = "Mule Deer", season = "Summer")
  md_wtr_elev <- high_elev_locs(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",], spp = "Mule Deer", season = "Winter")
  elk_smr_elev <- high_elev_locs(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19",], spp = "Elk", season = "Summer")
  elk_wtr_elev <- high_elev_locs(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920",], spp = "Elk", season = "Winter")
  wtd_smr_elev <- high_elev_locs(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19",], spp = "White-tailed Deer", season = "Summer")
  wtd_wtr_elev <- high_elev_locs(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920",], spp = "White-tailed Deer", season = "Winter")
  coug_smr_elev <- high_elev_locs(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",], spp = "Cougar", season = "Summer")
  coug_wtr_elev <- high_elev_locs(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",], spp = "Cougar", season = "Winter")
  wolf_smr_elev <- high_elev_locs(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",], spp = "Wolf", season = "Summer")
  wolf_wtr_elev <- high_elev_locs(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",], spp = "Wolf", season = "Winter")
  bob_smr_elev <- high_elev_locs(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",], spp = "Bobcat", season = "Summer")
  bob_wtr_elev <- high_elev_locs(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",], spp = "Bobcat", season = "Winter")
  coy_smr_elev <- high_elev_locs(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19",], spp = "Coyote", season = "Summer")
  coy_wtr_elev <- high_elev_locs(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920",], spp = "Coyote", season = "Winter")
  
  #'  Create one data frame with all results and save
  high_elev_deets <- as.data.frame(rbind(bob_smr_elev, bob_wtr_elev, coug_smr_elev, coug_wtr_elev,
                           coy_smr_elev, coy_wtr_elev, elk_smr_elev, elk_wtr_elev,
                           md_smr_elev, md_wtr_elev, wtd_smr_elev, wtd_wtr_elev,
                           wolf_smr_elev, wolf_wtr_elev))
  
  # write.csv(high_elev_deets, paste0("./Outputs/Tables/High_Elevation_Stats_BuffHR_", Sys.Date(), ".csv"))
  
  #'  Based on the above results, filter out high elevation sites for select species
  no_high_elev_locs <- function(locs) {
    drop_hi_elev <- locs %>%
      filter(locs$Elev < 2100)
    return(drop_hi_elev)
  }
  mdData_smr_low <- no_high_elev_locs(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",])
  mdData_wtr_low <- no_high_elev_locs(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",])
  cougData_smr_low <- no_high_elev_locs(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",])
  cougData_wtr_low <- no_high_elev_locs(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",])
  wolfData_smr_low <- no_high_elev_locs(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",])
  wolfData_wtr_low <- no_high_elev_locs(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",])
  bobData_smr_low <- no_high_elev_locs(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",])
  bobData_wtr_low <- no_high_elev_locs(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",])
  
  
  
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
    #'  Leave weights as is
    locs$w <- locs$w
    
    locs <- as.data.frame(locs)
  
    return(locs)
  }
  #'  Run season & species-specific data through prep function
  mdData_smr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Summer18" | md_dat_all$Season == "Summer19",])
  mdData_wtr <- spp_dataPrep(md_dat_all[md_dat_all$Season == "Winter1819" | md_dat_all$Season == "Winter1920",])
  elkData_smr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Summer18" | elk_dat_all$Season == "Summer19",])
  elkData_wtr <- spp_dataPrep(elk_dat_all[elk_dat_all$Season == "Winter1819" | elk_dat_all$Season == "Winter1920",])
  wtdData_smr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Summer18" | wtd_dat_all$Season == "Summer19",])
  wtdData_wtr <- spp_dataPrep(wtd_dat_all[wtd_dat_all$Season == "Winter1819" | wtd_dat_all$Season == "Winter1920",])
  cougData_smr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Summer18" | coug_dat_all$Season == "Summer19",])
  cougData_wtr <- spp_dataPrep(coug_dat_all[coug_dat_all$Season == "Winter1819" | coug_dat_all$Season == "Winter1920",])
  wolfData_smr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Summer18" | wolf_dat_all$Season == "Summer19",])
  wolfData_wtr <- spp_dataPrep(wolf_dat_all[wolf_dat_all$Season == "Winter1819" | wolf_dat_all$Season == "Winter1920",])
  bobData_smr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Summer18" | bob_dat_all$Season == "Summer19",])
  bobData_wtr <- spp_dataPrep(bob_dat_all[bob_dat_all$Season == "Winter1819" | bob_dat_all$Season == "Winter1920",])
  coyData_smr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Summer18" | coy_dat_all$Season == "Summer19",])
  coyData_wtr <- spp_dataPrep(coy_dat_all[coy_dat_all$Season == "Winter1819" | coy_dat_all$Season == "Winter1920",])
  
  #'  Format datasets with high elevation locations excluded
  mdData_smr_low <- spp_dataPrep(mdData_smr_low)
  mdData_wtr_low <- spp_dataPrep(mdData_wtr_low)
  cougData_smr_low <- spp_dataPrep(cougData_smr_low)
  cougData_wtr_low <- spp_dataPrep(cougData_wtr_low)
  wolfData_smr_low <- spp_dataPrep(wolfData_smr_low)
  wolfData_wtr_low <- spp_dataPrep(wolfData_wtr_low)
  bobData_smr_low <- spp_dataPrep(bobData_smr_low)
  bobData_wtr_low <- spp_dataPrep(bobData_wtr_low)
  
  #'  Function to create correlation matrix for all covariates at once
  cov_correlation <- function(dat) {
    used <- dat[dat$Used == 1,]
    covs <- used[,c("Elev", "Slope", "PercForMix", "PercXGrass", "PercXShrub", "RoadDen", "HumanMod")]
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (md_smr_corr <- cov_correlation(mdData_smr)) # Elev & HM -0.8111; Shrub & HM 0.8370; Shrub & Elev -0.6981
  (md_wtr_corr <- cov_correlation(mdData_wtr)) # Elev & HM -0.6427; Shrub & HM 0.6281
  (elk_smr_corr <- cov_correlation(elkData_smr)) # Elev & HM -0.6453
  (elk_wtr_corr <- cov_correlation(elkData_wtr)) # Forest & Grass -0.6426
  (wtd_smr_corr <- cov_correlation(wtdData_smr))
  (wtd_wtr_corr <- cov_correlation(wtdData_wtr))
  (coug_smr_corr <- cov_correlation(cougData_smr)) 
  (coug_wtr_corr <- cov_correlation(cougData_wtr)) # Forest & Grass -0.6712
  (wolf_smr_corr <- cov_correlation(wolfData_smr))
  (wolf_wtr_corr <- cov_correlation(wolfData_wtr)) #Forest & Grass -0.6969
  (bob_smr_corr <- cov_correlation(bobData_smr))
  (bob_wtr_corr <- cov_correlation(bobData_wtr)) # Forest & Shrub -0.7068
  (coy_smr_corr <- cov_correlation(coyData_smr)) # Elev & HM -0.7852
  (coy_wtr_corr <- cov_correlation(coyData_wtr)) # Elev & HM -0.7479; Forest & Grass -0.6146
  #'  Human Modified is definitely out

  cov_correlation(mdData_smr_low)
  cov_correlation(mdData_wtr_low)
  cov_correlation(cougData_smr_low)
  cov_correlation(cougData_wtr_low)
  cov_correlation(wolfData_smr_low)
  cov_correlation(wolfData_wtr_low)
  cov_correlation(bobData_smr_low)
  cov_correlation(bobData_wtr_low)
  
  
  #'  Resource Selection Function Models
  #'  ==================================
  #'  Functions to run logistic mixed effects models that include random effect 
  #'  for individual
  #'  Wolf, bobcat, & coyote RSFs exclude random effect for year because very few
  #'  collars were on air in both years
  #'  Mule deer summer RSF excludes random effect for year because too many
  #'  collars were on air in both years
  #'  Do I want study area as a fixed effect? Does that make sense? Including it
  #'  would be consistent with the occupancy models but a collared animal can't
  #'  "select" for one study area over the other since both study areas aren't
  #'  "available" to them. Not including.
  #'  Other habitat covariates excluded based on species and convergence issues
  
  ####  Mule Deer RSF  ####
  #'  SUMMERS 2018 & 2019
  #'  Dropping PercXShrub due to high correlation with PercForMix/Elev
  md_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                         data = mdData_smr, weights = w, family = binomial(link = "logit")) 
  summary(md_global_smr)
  car::vif(md_global_smr)

  #'  WINTERS 2018-2019 & 2019-2020
  md_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID),
                          data = mdData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(md_global_wtr)
  car::vif(md_global_wtr)
  
  ####  Mule Deer Subset RSFs  ####
  #'  Run RSF on summer mule deer data for ONLY locations within the OK study area
  library(sf)
  sa_proj <- st_crs("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  OK_SA <- st_transform(OK_SA, sa_proj)
  #'  Make mule deer input data spatial
  md_pts <- st_as_sf(mdData_smr, coords = c("x", "y"), crs = sa_proj)
  #'  Crop mule deer used/available data to the extent of the OK study area
  md_SA_only <- st_intersection(md_pts, OK_SA)
  #'  Make study area only data non-spatial
  md_SA_only_df <- as.data.frame(md_SA_only)
  #'  Run exact same RSF as above but with only relocation data from within study area
  md_SA_only_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),
                         data = md_SA_only_df, weights = w, family = binomial(link = "logit"))
  summary(md_SA_only_smr)
  car::vif(md_SA_only_smr)
  
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Mule deer summer
  #'  Run exact same RSF as above but with only locations <2100m
  md_lowElev_only_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),
                          data = mdData_smr_low, weights = w, family = binomial(link = "logit"))
  summary(md_lowElev_only_smr)
  car::vif(md_lowElev_only_smr)
  
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Mule deer winter
  md_lowElev_only_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID),
                         data = mdData_wtr_low, weights = w, family = binomial(link = "logit"))  
  summary(md_lowElev_only_wtr)
  car::vif(md_lowElev_only_wtr)
  
  
  ####  Elk RSF  ####
  #'  SUMMERS 2018 & 2019
  #'  Dropped PercXGrass & PercXShrub to be consistent with occupancy model 
  #'  DROPPING RANDOM EFFECT (AND USING GLM) due to singularity issue 
  elk_global_smr <- glm(Used ~ 1 + Elev + Slope + PercForMix + RoadDen,  #+ (1|ID)
                        data = elkData_smr, weights = w, family = binomial(link = "logit")) 
  #' GAH! warning: boundary (singular) fit: see ?isSingular   ---- OK if I drop random effect
  summary(elk_global_smr)
  car::vif(elk_global_smr)
  
  #'  Does something obvious stand out with the used/available locations?
  used_elk <- elkData_smr %>%
    filter(w == 1) %>%
    group_by(ID, Year) %>%
    summarise(n = n()) %>%
    ungroup()
  avail_elk <- elkData_smr %>%
    filter(w == 5000) %>%
    group_by(ID, Year) %>%
    summarise(n = n()) %>%
    ungroup()
  
  #'  WINTERS 2018-2019 & 2019-2020 
  #'  Dropped PercXGrass and PercXShrub due to high correlations with PercForMix 
  #'  and to be consistent with occupancy model
  elk_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID), 
                          data = elkData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(elk_global_wtr)
  car::vif(elk_global_wtr)
  
  ####  White-tailed Deer RSF  ####
  #'  SUMMERS 2018 & 2019
  #'  Dropping PercXGrass & PercXShurb to be consistent with occupancy model
  wtd_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),   
                          data = wtdData_smr, weights = w, family = binomial(link = "logit"))  # + I(Elev^2)
  summary(wtd_global_smr)
  car::vif(wtd_global_smr)

  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass & PercXShurb to be consistent with occupancy model
  wtd_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),  
                          data = wtdData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(wtd_global_wtr)
  car::vif(wtd_global_wtr)
  
  
  ####  Cougar RSF  ####
  #'  SUMMERS 2018 & 2019
  coug_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                          data = cougData_smr, weights = w, family = binomial(link = "logit")) 
  summary(coug_global_smr)
  car::vif(coug_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass due to high correlation with PercForMix
  coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                           data = cougData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(coug_global_wtr)
  car::vif(coug_global_wtr)
  
  ####  Cougar Subset RSFs  ####
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Cougar summer
  #'  Run exact same RSF as above but with only locations <2100m
  coug_lowElev_only_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                                 data = cougData_smr_low, weights = w, family = binomial(link = "logit"))
  summary(coug_lowElev_only_smr)
  car::vif(coug_lowElev_only_smr)
  
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Cougar winter
  coug_lowElev_only_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                                 data = cougData_wtr_low, weights = w, family = binomial(link = "logit"))  
  summary(coug_lowElev_only_wtr)
  car::vif(coug_lowElev_only_wtr)
  
  
  ####  Wolf RSF  ####
  #'  SUMMERS 2018 & 2019
  #'  Dropping PercXShrub to be consistent with occupancy model 
  wolf_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),  
                       data = wolfData_smr, weights = w, family = binomial(link = "logit"))
  summary(wolf_global_smr)
  car::vif(wolf_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXShrub to be consistent with wolf occupancy model
  #'  Dropping PercXGrass due to high correlation with PercForMix  
  wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID), 
                           data = wolfData_wtr, weights = w, family = binomial(link = "logit"))  
  summary(wolf_global_wtr)
  car::vif(wolf_global_wtr)

  ####  Wolf Subset RSFs  ####
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Wolf summer
  #'  Run exact same RSF as above but with only locations <2100m
  wolf_lowElev_only_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID),  
                                 data = wolfData_smr_low, weights = w, family = binomial(link = "logit"))
  summary(wolf_lowElev_only_smr)
  car::vif(wolf_lowElev_only_smr)
  
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Wolf winter
  wolf_lowElev_only_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID),  
                                 data = wolfData_wtr_low, weights = w, family = binomial(link = "logit"))  
  summary(wolf_lowElev_only_wtr)
  car::vif(wolf_lowElev_only_wtr)
  
  
  ####  Bobcat RSF  ####
  #'  SUMMERS 2018 & 2019
  bob_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                           data = bobData_smr, weights = w, family = binomial(link = "logit")) 
  summary(bob_global_smr)
  car::vif(bob_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXShrub due to high correlation with PercForMix
  bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                          data = bobData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(bob_global_wtr)
  car::vif(bob_global_wtr)
  
  ####  Bobcat Subset RSFs  ####
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: bobcat summer
  #'  Run exact same RSF as above but with only locations <2100m
  bob_lowElev_only_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                                 data = bobData_smr_low, weights = w, family = binomial(link = "logit"))
  summary(bob_lowElev_only_smr)
  car::vif(bob_lowElev_only_smr)
  
  #'  HIGH ELEVATION LOCATIONS EXCLUDED: Wolf winter
  bob_lowElev_only_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                                 data = bobData_wtr_low, weights = w, family = binomial(link = "logit"))  
  summary(bob_lowElev_only_wtr)
  car::vif(bob_lowElev_only_wtr)
  
  
  ####  Coy RSF  ####
  #'  SUMMERS 2018 & 2019
  coy_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                           data = coyData_smr, weights = w, family = binomial(link = "logit")) #+ I(Elev^2)
  summary(coy_global_smr)
  car::vif(coy_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXGrass due to high correlation with PercForMix
  coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                          data = coyData_wtr, weights = w, family = binomial(link = "logit")) 
  summary(coy_global_wtr)
  car::vif(coy_global_wtr)
  
  
  #'  Save
  save(md_global_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_BuffHR_", Sys.Date(), ".RData"))  
  save(md_global_wtr, file = paste0("./Outputs/RSF_output/md_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(elk_global_smr, file = paste0("./Outputs/RSF_output/elk_RSF_smr_BuffHR_", Sys.Date(), ".RData")) #' Note, currently excludes the random effect for ID
  save(elk_global_wtr, file = paste0("./Outputs/RSF_output/elk_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(wtd_global_smr, file = paste0("./Outputs/RSF_output/wtd_RSF_smr_BuffHR_", Sys.Date(), ".RData"))
  save(wtd_global_wtr, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(coug_global_smr, file = paste0("./Outputs/RSF_output/coug_RSF_smr_BuffHR_", Sys.Date(), ".RData"))
  save(coug_global_wtr, file = paste0("./Outputs/RSF_output/coug_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(wolf_global_smr, file = paste0("./Outputs/RSF_output/wolf_RSF_smr_BuffHR_", Sys.Date(), ".RData"))
  save(wolf_global_wtr, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(bob_global_smr, file = paste0("./Outputs/RSF_output/bob_RSF_smr_BuffHR_", Sys.Date(), ".RData"))
  save(bob_global_wtr, file = paste0("./Outputs/RSF_output/bob_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  save(coy_global_smr, file = paste0("./Outputs/RSF_output/coy_RSF_smr_BuffHR_", Sys.Date(), ".RData"))
  save(coy_global_wtr, file = paste0("./Outputs/RSF_output/coy_RSF_wtr_BuffHR_", Sys.Date(), ".RData"))
  
  # save(md_SA_only_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_BuffHR_SAonly_", Sys.Date(), ".RData"))
  # 
  # save(md_lowElev_only_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(md_lowElev_only_wtr, file = paste0("./Outputs/RSF_output/md_RSF_wtr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(coug_lowElev_only_smr, file = paste0("./Outputs/RSF_output/coug_RSF_smr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(coug_lowElev_only_wtr, file = paste0("./Outputs/RSF_output/coug_RSF_wtr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(wolf_lowElev_only_smr, file = paste0("./Outputs/RSF_output/wolf_RSF_smr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(wolf_lowElev_only_wtr, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(bob_lowElev_only_smr, file = paste0("./Outputs/RSF_output/bob_RSF_smr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  # save(bob_lowElev_only_wtr, file = paste0("./Outputs/RSF_output/bob_RSF_wtr_BuffHR_lowElev_", Sys.Date(), ".RData"))
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Pull out RSF results
  load("./Outputs/RSF_output/md_RSF_smr_BuffHR_2022-05-03.RData")  #' Make sure I read in the right dataset 
  load("./Outputs/RSF_output/md_RSF_wtr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/elk_RSF_smr_BuffHR_2022-05-03.RData") # Note: excludes random effect
  load("./Outputs/RSF_output/elk_RSF_wtr_BuffHR_2022-05-03.RData") 
  load("./Outputs/RSF_output/wtd_RSF_smr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/wtd_RSF_wtr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/coug_RSF_smr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/coug_RSF_wtr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/wolf_RSF_smr_BuffHR_2022-05-03.RData") 
  load("./Outputs/RSF_output/wolf_RSF_wtr_BuffHR_2022-05-03.RData") 
  load("./Outputs/RSF_output/bob_RSF_smr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/bob_RSF_wtr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/coy_RSF_smr_BuffHR_2022-05-03.RData")
  load("./Outputs/RSF_output/coy_RSF_wtr_BuffHR_2022-05-03.RData")
  
  # load("./Outputs/RSF_output/md_RSF_smr_SAonly_BuffHR_2022-05-03.RData")
  

  #'  Function to save parameter estimates & p-values
  #'  use coef(mod) to look at random effects estimates
  rounddig <- 2
  
  rsf_out <- function(mod, spp, season){
    # betas <- mod@beta
    # se <- sqrt(diag(vcov(mod)))
    betas <- summary(mod)$coef[,1]
    se <- summary(mod)$coef[,2]
    z <- summary(mod)$coef[,3]
    pval <- summary(mod)$coef[,4]
    out <- as.data.frame(cbind(betas, se, pval)) %>%
      transmute(
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Parameter = row.names(.),
        Estimate = round(betas, rounddig),
        SE = round(se, rounddig),
        Z = round(z, rounddig),
        Pval = round(pval, rounddig)) 
    rownames(out) <- NULL
    return(out)
  }
  md_s1819_rsf <- rsf_out(md_global_smr, "Mule Deer", "Summer")
  md_w1820_rsf <- rsf_out(md_global_wtr, "Mule Deer", "Winter")
  elk_s1819_rsf <- rsf_out(elk_global_smr, "Elk", "Summer")
  elk_w1820_rsf <- rsf_out(elk_global_wtr, "Elk", "Winter")
  wtd_s1819_rsf <- rsf_out(wtd_global_smr, "White-tailed Deer", "Summer")
  wtd_w1820_rsf <- rsf_out(wtd_global_wtr, "White-tailed Deer", "Winter")
  coug_s1819_rsf <- rsf_out(coug_global_smr, "Cougar", "Summer")
  coug_w1820_rsf <- rsf_out(coug_global_wtr, "Cougar", "Winter")
  wolf_s1819_rsf <- rsf_out(wolf_global_smr, "Wolf", "Summer")
  wolf_w1820_rsf <- rsf_out(wolf_global_wtr, "Wolf", "Winter")
  bob_s1819_rsf <- rsf_out(bob_global_smr, "Bobcat", "Summer")
  bob_w1820_rsf <- rsf_out(bob_global_wtr, "Bobcat", "Winter")
  coy_s1819_rsf <- rsf_out(coy_global_smr, "Coyote", "Summer")
  coy_w1820_rsf <- rsf_out(coy_global_wtr, "Coyote", "Winter")
  
  md_sSAonly_rsf <- rsf_out(md_SA_only_smr, "Mule Deer", "Summer")
  # write.csv(md_sSAonly_rsf, paste0("./Outputs/Tables/RSF_Results_MuleDeer_SAonly_", Sys.Date(), ".csv"))  

  #'  Merge into larger data frames for easy comparison
  summer_rsf <- rbind(bob_s1819_rsf, coug_s1819_rsf, coy_s1819_rsf, wolf_s1819_rsf,
                      elk_s1819_rsf, md_s1819_rsf, wtd_s1819_rsf) 
  winter_rsf <- rbind(bob_w1820_rsf, coug_w1820_rsf, coy_w1820_rsf, wolf_w1820_rsf,
                      elk_w1820_rsf, md_w1820_rsf, wtd_w1820_rsf) 
  rsf_results <- rbind(summer_rsf, winter_rsf) %>%
    arrange(Species)
  colnames(rsf_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  rsf_results_wide <- rsf_results %>% 
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
    # separate("AreaOK", c("AreaOK (SE)", "AreaOK Pval"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev Pval"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope Pval"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix Pval"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass Pval"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub Pval"), sep = "_") %>%
    separate("RoadDen", c("Road Density (SE)", "Road Density Pval"), sep = "_") %>%
    # separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    arrange(match(Season, c("Summer", "Winter")))
  
  
  #' #'  Save!
  #' write.csv(rsf_results, paste0("./Outputs/Tables/RSF_Results_BuffHR_", Sys.Date(), ".csv"))  
  #' write.csv(rsf_results_wide, paste0("./Outputs/Tables/RSF_Results_wide_BuffHR_", Sys.Date(), ".csv"))
  #' 
  #' save.image("./Outputs/RSF_script_results.RData")
  
  
  
 