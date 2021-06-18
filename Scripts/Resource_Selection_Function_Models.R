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
  
  #' ####  Mule Deer RSF  ####
  #' #'  Global model with random effect for individual and year
  #' #'  SUMMERS 2018 & 2019, OK study area only = no Area effect
  #' md_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                         data = mdData_smr, family = binomial(link = "logit"))
  #' summary(md_global_smr)
  #' 
  #' #'  WINTERS 2018-2019 & 2019-2020, OK study area only = no Area effect
  #' md_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                         data = mdData_wtr, family = binomial(link = "logit"))
  #' summary(md_global_wtr)
  
  
  ####  Elk RSF  ####
  #'  Global model with random effect for individual and year
  #'  SUMMERS 2018 & 2019, NE study area only = no Area effect
  elk_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                        data = elkData_smr, family = binomial(link = "logit"))
  summary(elk_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020, NE study area only = no Area effect
  elk_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
                      data = elkData_wtr, family = binomial(link = "logit"))
  summary(elk_global_wtr)

  
  #' ####  White-tailed Deer RSF  ####
  #' #'  Global model with random effect for individual and year
  #' #'  SUMMERS 2018 & 2019, NE study area only = no Area effect
  #' wtd_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                         data = wtdData_smr, family = binomial(link = "logit"))
  #' summary(wtd_global_smr)
  #' 
  #' #'  WINTERS 2018-2019 & 2019-2020, NE study area only = no Area effect
  #' wtd_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year),
  #'                         data = wtdData_wtr, family = binomial(link = "logit"))
  #' summary(wtd_global_wtr)
  
  
  ####  Cougar RSF  ####
  #'  Random effect for individual, year and study area
  #'  SUMMERS 2018 & 2019
  coug_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
                          data = cougData_smr, family = binomial(link = "logit"))
  summary(coug_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020
  coug_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + RoadDen + HumanMod + (1|ID) + (1|Year) + (1|Area),
                          data = cougData_wtr, family = binomial(link = "logit"))
  summary(coug_global_wtr)
  
  
  ####  Wolf RSF  ####
  #'  Random effect for individual, year and study area
  #'  SUMMERS 2018 & 2019
  wolf_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                       data = wolfData_smr, family = binomial(link = "logit"))
  summary(wolf_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  wolf_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = wolfData_wtr, family = binomial(link = "logit"))
  summary(wolf_global_wtr)
  
  
  ####  Bobcat RSF  ####
  #'  Random effect for individual, year and study area
  #'  SUMMERS 2018 & 2019
  bob_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = bobData_smr, family = binomial(link = "logit"))
  summary(bob_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  bob_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = bobData_wtr, family = binomial(link = "logit"))
  summary(bob_global_wtr)
  
  
  ####  Coy RSF  ####
  #'  Random effect for individual, year and study area
  #'  SUMMERS 2018 & 2019
  coy_global_smr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = coyData_smr, family = binomial(link = "logit"))
  summary(coy_global_smr)
  
  #'  WINTERS 2018-2019 & 2019-2020  
  coy_global_wtr <- glmer(Used ~ 1 + Elev + Slope + PercForMix + PercXGrass + PercXShrub + HumanMod + (1|ID) + (1|Year) + (1|Area),
                           data = coyData_wtr, family = binomial(link = "logit"))
  summary(coy_global_wtr)
  
  
  #'  Save
  # save(md_global_smr, file = "./Outputs/RSF_output/md_RSF_smr.RData")
  # save(md_global_wtr, file = "./Outputs/RSF_output/md_RSF_wtr.RData")
  save(elk_global_smr, file = "./Outputs/RSF_output/elk_RSF_smr.RData")
  save(elk_global_wtr, file = "./Outputs/RSF_output/elk_RSF_wtr.RData")
  # save(wtd_global_smr, file = "./Outputs/RSF_output/wtd_RSF_smr.RData")
  # save(wtd_global_wtr, file = "./Outputs/RSF_output/wtd_RSF_wtr.RData")
  save(coug_global_smr, file = "./Outputs/RSF_output/coug_RSF_smr.RData")
  save(coug_global_wtr, file = "./Outputs/RSF_output/coug_RSF_wtr.RData")
  save(wolf_global_smr, file = "./Outputs/RSF_output/wolf_RSF_smr.RData")
  save(wolf_global_wtr, file = "./Outputs/RSF_output/wolf_RSF_wtr.RData")
  save(bob_global_smr, file = "./Outputs/RSF_output/bob_RSF_smr.RData")
  save(bob_global_wtr, file = "./Outputs/RSF_output/bob_RSF_wtr.RData")
  save(coy_global_smr, file = "./Outputs/RSF_output/coy_RSF_smr.RData")
  save(coy_global_wtr, file = "./Outputs/RSF_output/coy_RSF_wtr.RData")
  
  
  ####  Summary tables  ####
  #'  Save model outputs in table format 
  #'  Functions extract outputs for each sub-model and appends species/season info
  
  #'  Pull out RSF results
  # load("./Outputs/RSF_output/md_RSF_smr.RData")
  # load("./Outputs/RSF_output/md_RSF_wtr.RData")
  # load("./Outputs/RSF_output/elk_RSF_smr.RData")
  # load("./Outputs/RSF_output/elk_RSF_wtr.RData")
  # load("./Outputs/RSF_output/wtd_RSF_smr.RData")
  # load("./Outputs/RSF_output/wtd_RSF_wtr.RData")
  # load("./Outputs/RSF_output/coug_RSF_smr.RData")
  # load("./Outputs/RSF_output/coug_RSF_wtr.RData")
  # load("./Outputs/RSF_output/wolf_RSF_smr.RData")
  # load("./Outputs/RSF_output/wolf_RSF_wtr.RData")
  # load("./Outputs/RSF_output/bob_RSF_smr.RData")
  # load("./Outputs/RSF_output/bob_RSF_wtr.RData")
  # load("./Outputs/RSF_output/coy_RSF_smr.RData")
  # load("./Outputs/RSF_output/coy_RSF_wtr.RData")
  

  #'  Function to save parameter estimates & p-values
  #'  use coef(mod) to look at random effects estimates
  rounddig <- 3
  
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
  # md_s1819_rsf <- rsf_out(md_global_smr, "Mule Deer", "Summer")
  # md_w1820_rsf <- rsf_out(md_global_wtr, "Mule Deer", "Winter")
  elk_s1819_rsf <- rsf_out(elk_global_smr, "Elk", "Summer")
  elk_w1820_rsf <- rsf_out(elk_global_wtr, "Elk", "Winter")
  # wtd_s1819_rsf <- rsf_out(wtd_global_smr, "White-tailed Deer", "Summer")
  # wtd_w1820_rsf <- rsf_out(wtd_global_wtr, "White-tailed Deer", "Winter")
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
                      elk_s1819_rsf) #, md_s1819_rsf, wtd_s1819_rsf
  winter_rsf <- rbind(bob_w1820_rsf, coug_w1820_rsf, coy_w1820_rsf, wolf_w1820_rsf,
                      elk_w1820_rsf) #, md_w1820_rsf, wtd_w1820_rsf
  rsf_results <- rbind(summer_occ, winter_occ) %>%
    arrange(Species)
  colnames(rsf_results) <- c("Species", "Season", "Parameter", "Estimate", "SE", "z", "Pval")
  
  
  #'  Spread this out so the coefficient effects are easier to compare across species
  rsf_results_wide <- coug_s1819_rsf %>% #rsf_results
    dplyr::select(-Z) %>%
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
    separate("HumanMod", c("HumanMod (SE)", "HumanMod Pval"), sep = "_") #%>%
    # arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf", "Mule Deer", "Elk", "White-tailed Deer"))) %>%
    # arrange(match(Season, c("Summer", "Winter")))
  
  
  #'  Save!
  write.csv(rsf_results, paste0("./Outputs/RSF_Results_", Sys.Date(), ".csv"))
  write.csv(rsf_results_wide, paste0("./Outputs/RSF_Results_wide_", Sys.Date(), ".csv"))
  
  
  
  
  
  