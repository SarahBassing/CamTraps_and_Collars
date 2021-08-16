  #'  =====================================
  #'  Hypothesis testing: Cam vs Collars
  #'  Washington Predator-Prey Project
  #'  Sarah B. Bassing
  #'  August 2021
  #'  ======================================
  #'  Using GLMMs to test whether collaring effort or camera deployment influenced
  #'  whether we found agreement between camera trap-based occupancy model
  #'  estimates and satellite collar-based resource selection functions. 
  #'  
  #'  We run two separate tests: one to test for significant effects on the
  #'  number of agreements we found and one to test for significant effects on
  #'  the number of disagreements we found between the two models. These only
  #'  include instances where both models found significant (or trending towards
  #'  significant effects). 
  #'  
  #'  Mean occupancy and detection probabilities and road effect estimates come 
  #'  from the Occupancy_Models.R script and rely on me predicting and 
  #'  back-transforming these estimates correctly.
  #'  ======================================
  
  #'  Clear memory
  rm(list=ls())
  
  #'  Load libraries
  library(lme4)
  library(tidyverse)
  
  #'  Read in occupancy model results
  OccMod_est <- read.csv("./Outputs/Tables/OccMod_Mean_Estimates_noHM_2021-08-15.csv") %>% # make sure you have the right file!
    dplyr::select(-c(X, SE))
  occ_mu <- OccMod_est %>%
    filter(Parameter == "Occupancy") %>%
    transmute(
      Species = Species,
      Season = Season,
      Mean_Occ = Mean
    )
  det_mu <- OccMod_est %>%
    filter(Parameter == "Detection") %>%
    transmute(
      Species = Species,
      Season = Season,
      Mean_Det = Mean
    )
  
  #'  Read in estimated road v trail effect on detection probability 
  detprob <- read.csv("./Outputs/Tables/OccMod_DetProb_Results_noHM_2021-08-15.csv") %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    filter(Parameter  == "(Intercept)" | Parameter == "TrailDirt road" | Parameter == "TrailDecommissioned road") %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Trail_int", Parameter),
           Parameter = ifelse(Parameter == "TrailDirt road", "Road", Parameter),
           Parameter = ifelse(Parameter == "TrailDecommissioned road", "Decommed", Parameter),
           #'  Change estimates to 0 if effect was non-significant
           Estimate = ifelse(Pval > 0.1, Estimate == 0, Estimate)) %>%
    dplyr::select(-Pval) %>%
    pivot_wider(names_from = Parameter, values_from = Estimate)
  
  #'  Species with FEMALE-ONLY collaring efforts (F = 1, M&F = 0)
  Species <- c("Bobcat", "Cougar", "Coyote", "Elk", "Mule Deer", "White-tailed Deer", "Wolf")
  Female_only <- c(0,0,0,1,1,1,0)
  Collared_sex <- as.data.frame(cbind(Species, Female_only))
  
  #'  Build input table of independent variables for hypothesis testing
  dat_x <- full_join(occ_mu, det_mu, by = c("Species", "Season")) %>%
    full_join(detprob, by = c("Species", "Season")) %>%
    arrange(Species, Season) %>%
    left_join(Collared_sex, by = "Species") %>%
    mutate(
      Mean_Occ = scale(Mean_Occ),
      Mean_Det = scale(Mean_Det),
      Trail_int = scale(Trail_int),
      Road = scale(Road),
      Decommed = scale(Decommed),
      Female_only = as.factor(Female_only)
      )
  
  #'  Agreement & disagreement between significant Occupancy & RSF estimates
  #'  Counts are the number of covariates per species & season where we saw 
  #'  agreement (++ or --) or disagreement (+- or -+) between Occ & RSF estimates
  dat_y <- dplyr::select(dat_x, c(Species, Season))
  #'  Be sure to follow correct order of how species & season are listed!  ## MAKE SURE THESE ARE UP TO DATE!!!!! 
  Agree <- c(1,0,1,2,2,1,0,0,1,1,2,1,1,0)  # WITHOUT HM in models
  Disagree <- c(0,0,1,0,0,1,0,0,0,0,0,0,0,0)  # WITHOUT HM in models
  # Agree <- c(1,1,1,3,2,1,0,0,0,1,2,1,1,0)  # WITH HM in models
  # Disagree <- c(0,1,1,0,0,1,0,0,1,0,0,0,1,0)  # WITH HM in models
  #'  Format so each row is 1 instance of an agreement per species & season
  dat_a <- cbind(dat_y, Agree) %>%
    group_by(Species, Season) %>%
    #'  Add rows for a given species & season if Agree > 1 (completes sequence of 1 thru x)
    complete(Agree = full_seq(1:Agree, 1)) %>%
    #'  Add column where all values > 1 become 1 since each row is equivalent to
    #'  as single instance of agreement
    mutate(Agree2 = ifelse(Agree >= 1, 1, 0)) %>%
    #'  full_seq adds a row for species with 0 instances of agreement so must 
    #'  remove those based on the correct species and season
    #'  DOUBLE CHECK THAT THE CORRECT SEASON WAS REMOVED!!!!
    filter(Species != "Elk" | Agree2 < 1) %>%
    filter(Species != "Mule Deer" | Season != "Summer" | Agree2 < 1) %>%
    filter(Species != "Wolf" | Season != "Winter" | Agree2 < 1) %>%
    ungroup() %>%
    transmute(Species = as.factor(Species),
              Season = as.factor(Season),
              Agree = as.numeric(Agree2))
  #'  Format so each row is 1 instance of a disagreement per species & season
  #'  Super easy since never more than 1 instance per species-season combo
  dat_d <- cbind(dat_y, Disagree) 
  
  #'  Combine independent data with reformatted agreement & disagreement data
  dat_agree <- dat_a %>%
    full_join(dat_x, by = c("Species", "Season"))
  dat_disagree <- dat_d %>%
    full_join(dat_x, by = c("Species", "Season"))
  
  
  ####  Hypothesis testing  ####
  #'  Agreement between OccMod and RSF
  #'  Start with GLMMs, includes random effect for species
  mod1 <- glmer(Agree ~ Mean_Det + (1 | Species), family = binomial, data = dat_agree)
  mod2 <- glmer(Agree ~ Mean_Occ + (1 | Species), family = binomial, data = dat_agree)
  mod3 <- glmer(Agree ~ Female_only + (1 | Species), family = binomial, data = dat_agree)
  mod4 <- glmer(Agree ~ Road + (1 | Species), family = binomial, data = dat_agree)
  summary(mod1)
  summary(mod2)
  summary(mod3)
  summary(mod4)
  AIC(mod1, mod2, mod3, mod4)
  #'  No good. 
  
  #'  Try a more basic GLM to estimate probability of an agreement
  mod1 <- glm(Agree ~ Mean_Det, family = binomial, data = dat_agree)
  mod2 <- glm(Agree ~ Mean_Occ , family = binomial, data = dat_agree)
  mod3 <- glm(Agree ~ Female_only, family = binomial, data = dat_agree)
  mod4 <- glm(Agree ~ Road, family = binomial, data = dat_agree)
  summary(mod1)
  summary(mod2)
  summary(mod3)
  summary(mod4)
  AIC(mod1, mod2, mod3, mod4)
  
  #'  Should this be a Poisson regression instead? Estimating the number of agreements?
  mod1 <- glm(Agree ~ Mean_Det, family = poisson, data = dat_agree)
  mod2 <- glm(Agree ~ Mean_Occ , family = poisson, data = dat_agree)
  mod3 <- glm(Agree ~ Female_only, family = poisson, data = dat_agree)
  mod4 <- glm(Agree ~ Road, family = poisson, data = dat_agree)
  summary(mod1)
  summary(mod2)
  summary(mod3)
  summary(mod4)
  AIC(mod1, mod2, mod3, mod4)
  
  
  #'  Disagreement between OccMod and RSF
  #'  Start with GLMMs, includes random effect for species
  mod1 <- glmer(Disagree ~ Mean_Det + (1 | Species), family = binomial, data = dat_disagree)
  mod2 <- glmer(Disagree ~ Mean_Occ + (1 | Species), family = binomial, data = dat_disagree)
  mod3 <- glmer(Disagree ~ Female_only + (1 | Species), family = binomial, data = dat_disagree)
  mod4 <- glmer(Disagree ~ Road + (1 | Species), family = binomial, data = dat_disagree)
  summary(mod1)
  summary(mod2)
  summary(mod3)
  summary(mod4)
  AIC(mod1, mod2, mod3, mod4)
  #'  BAH! no good. Random effect definitely not necessary.
  
  #'  Try a more basic GLM: logistic model
  mod1 <- glm(Disagree ~ Mean_Det, family = binomial, data = dat_disagree)
  mod2 <- glm(Disagree ~ Mean_Occ , family = binomial, data = dat_disagree)
  mod3 <- glm(Disagree ~ Female_only, family = binomial, data = dat_disagree)
  mod4 <- glm(Disagree ~ Road, family = binomial, data = dat_disagree)
  summary(mod1)
  summary(mod2)
  summary(mod3)
  summary(mod4)
  AIC(mod1, mod2, mod3, mod4)
  
  
  
  
  
  
  
  