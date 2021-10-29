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
  #'  3rd Order Selection
  # load("./Outputs/RSF_pts/md_dat_all_2021-06-30.RData")  # 2021-06-22 uses reprojected rasters
  # load("./Outputs/RSF_pts/elk_dat_all_2021-06-30.RData")
  # load("./Outputs/RSF_pts/wtd_dat_all_2021-06-30.RData")
  # load("./Outputs/RSF_pts/coug_dat_all_2021-06-30.RData")
  # load("./Outputs/RSF_pts/wolf_dat_all_2021-06-30.RData")
  # load("./Outputs/RSF_pts/bob_dat_all_2021-06-30.RData")
  # load("./Outputs/RSF_pts/coy_dat_all_2021-06-30.RData")
  #'  2nd Order Selection
  load("./Outputs/RSF_pts/md_dat_2nd_all_2021-09-13.RData")  # 2021-08-10 uses un-buffered and un-masked MCPs
  load("./Outputs/RSF_pts/elk_dat_2nd_all_2021-09-13.RData")
  load("./Outputs/RSF_pts/wtd_dat_2nd_all_2021-09-13.RData")
  load("./Outputs/RSF_pts/coug_dat_2nd_all_2021-09-13.RData")
  load("./Outputs/RSF_pts/wolf_dat_2nd_all_2021-09-13.RData")
  load("./Outputs/RSF_pts/bob_dat_2nd_all_2021-09-13.RData")
  load("./Outputs/RSF_pts/coy_dat_2nd_all_2021-09-13.RData")
  
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
  #'  Global model with random effect for individual and year (winter only)
  #'  SUMMERS 2018 & 2019
  #'  Dropping PercXShrub due to high correlation with PercForMix
  # md_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
  #                         data = mdData_smr, family = binomial(link = "logit"))
  md_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                         data = mdData_smr, family = binomial(link = "logit"))
  summary(md_global_smr)
  #'  Singularity issues caused when year is included- probably because almost all collars deployed both summers so no variation
  car::vif(md_global_smr)

  #'  WINTERS 2018-2019 & 2019-2020
  md_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year),
                          data = mdData_wtr, family = binomial(link = "logit"))
  summary(md_global_wtr)
  car::vif(md_global_wtr)
  
  
  ####  Elk RSF  ####
  #'  Global model with random effect for individual and year
  #'  SUMMERS 2018 & 2019
  elk_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year),  
                        data = elkData_smr, family = binomial(link = "logit"))
  summary(elk_global_smr)
  car::vif(elk_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020 
  #'  Dropped PercXGrass and PercXShrub due to high correlations with PercForMix 
  #'  causing problems with predictions
  # elk_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
  #                     data = elkData_wtr, family = binomial(link = "logit"))
  elk_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + RoadDen + (1|ID) + (1|Year), 
                          data = elkData_wtr, family = binomial(link = "logit"))
  summary(elk_global_wtr)
  car::vif(elk_global_wtr)
  
  ####  White-tailed Deer RSF  ####
  #'  Global model with random effect for individual and year
  #'  SUMMERS 2018 & 2019
  wtd_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = wtdData_smr, family = binomial(link = "logit"))
  summary(wtd_global_smr)
  car::vif(wtd_global_smr)

  #'  WINTERS 2018-2019 & 2019-2020
  wtd_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = wtdData_wtr, family = binomial(link = "logit"))
  summary(wtd_global_wtr)
  car::vif(wtd_global_wtr)
  
  
  ####  Cougar RSF  ####
  #'  Random effect for individual & year
  #'  SUMMERS 2018 & 2019
  coug_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = cougData_smr, family = binomial(link = "logit"))
  summary(coug_global_smr)
  car::vif(coug_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass due to high correlation with PercForMix
  # coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
  #                         data = cougData_wtr, family = binomial(link = "logit"))
  coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                           data = cougData_wtr, family = binomial(link = "logit"))
  summary(coug_global_wtr)
  car::vif(coug_global_wtr)
  
  
  ####  Wolf RSF  ####
  #'  Random effect for individual; too few collars active both years for year effect
  #'  SUMMERS 2018 & 2019
  wolf_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                       data = wolfData_smr, family = binomial(link = "logit"))  #1/9 collars active both years
  summary(wolf_global_smr)
  car::vif(wolf_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  #'  Dropping PercXGrass due to high correlation with PercForMix  
  # wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
  #                          data = wolfData_wtr, family = binomial(link = "logit")) #no collars active both years
  wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID), 
                           data = wolfData_wtr, family = binomial(link = "logit")) #no collars active both years
  summary(wolf_global_wtr)
  car::vif(wolf_global_wtr)

  
  ####  Bobcat RSF  ####
  #'  Random effect for individual; too few collars active both years for year effect
  #'  SUMMERS 2018 & 2019
  bob_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
                           data = bobData_smr, family = binomial(link = "logit")) #1/10 collars active both years
  summary(bob_global_smr)
  car::vif(bob_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXShrub due to high correlation with PercForMix
  # bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID), 
  #                          data = bobData_wtr, family = binomial(link = "logit")) #1/11 collars active both years
  bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + RoadDen + (1|ID), 
                          data = bobData_wtr, family = binomial(link = "logit")) #1/11 collars active both years
  summary(bob_global_wtr)
  car::vif(bob_global_wtr)
  #  Winter1920 MVBOB71M disperses somewhere on the Colville reservation between 
  #  two study areas and make the OK bobcat MCP super big as a result. Cut his
  #  winter 1920 data completely (exclude from RSF and MCP).
  
  
  ####  Coy RSF  ####
  #'  Random effect for individual and year
  #'  SUMMERS 2018 & 2019
  coy_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                           data = coyData_smr, family = binomial(link = "logit"))  #5/16 collars active both years
  summary(coy_global_smr)
  car::vif(coy_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  #'  Dropping PercXGrass due to high correlation with PercForMix
  # coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + (1|ID) + (1|Year), 
  #                          data = coyData_wtr, family = binomial(link = "logit")) #4/19 collars active both years
  coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXShrub + RoadDen + (1|ID) + (1|Year), 
                          data = coyData_wtr, family = binomial(link = "logit")) #4/19 collars active both years
  summary(coy_global_wtr)
  car::vif(coy_global_wtr)
  
  
  #'  Save
  save(md_global_smr, file = paste0("./Outputs/RSF_output/md_RSF_smr_NoHM_", Sys.Date(), ".RData"))  #'  KEEP TRACK of whether human modified was excluded from models!
  save(md_global_wtr, file = paste0("./Outputs/RSF_output/md_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(elk_global_smr, file = paste0("./Outputs/RSF_output/elk_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(elk_global_wtr, file = paste0("./Outputs/RSF_output/elk_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(wtd_global_smr, file = paste0("./Outputs/RSF_output/wtd_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(wtd_global_wtr, file = paste0("./Outputs/RSF_output/wtd_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(coug_global_smr, file = paste0("./Outputs/RSF_output/coug_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(coug_global_wtr, file = paste0("./Outputs/RSF_output/coug_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(wolf_global_smr, file = paste0("./Outputs/RSF_output/wolf_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(wolf_global_wtr, file = paste0("./Outputs/RSF_output/wolf_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(bob_global_smr, file = paste0("./Outputs/RSF_output/bob_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(bob_global_wtr, file = paste0("./Outputs/RSF_output/bob_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  save(coy_global_smr, file = paste0("./Outputs/RSF_output/coy_RSF_smr_NoHM_", Sys.Date(), ".RData"))
  save(coy_global_wtr, file = paste0("./Outputs/RSF_output/coy_RSF_wtr_NoHM_", Sys.Date(), ".RData"))
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Pull out RSF results
  load("./Outputs/RSF_output/md_RSF_smr_noHM_2021-09-23.RData")  #' Make sure I read in the right dataset 2021-09-13
  load("./Outputs/RSF_output/md_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/elk_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/elk_RSF_wtr_noHM_2021-09-23.RData") 
  load("./Outputs/RSF_output/wtd_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/wtd_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coug_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coug_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/wolf_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/wolf_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/bob_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/bob_RSF_wtr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coy_RSF_smr_noHM_2021-09-23.RData")
  load("./Outputs/RSF_output/coy_RSF_wtr_noHM_2021-09-23.RData")
  

  #'  Function to save parameter estimates & p-values
  #'  use coef(mod) to look at random effects estimates
  rounddig <- 2
  
  rsf_out <- function(mod, spp, season){
    betas <- mod@beta
    se <- sqrt(diag(vcov(mod)))
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
  
  
  #'  Save!
  write.csv(rsf_results, paste0("./Outputs/Tables/RSF_Results_NoHM_", Sys.Date(), ".csv"))  #'  KEEP TRACK of whether Human Modified was excluded from models!
  write.csv(rsf_results_wide, paste0("./Outputs/Tables/RSF_Results_wide_NoHM_", Sys.Date(), ".csv"))
  
  save.image("./Outputs/RSF_script_results.RData")
  
  
  
  #'  SHOULD I BE CONSIDERING A VARIENCE INFLATION FACTOR ON THESE ESTIMATES????
  