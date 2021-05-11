  #'  ============================================
  #'  Hidden Markove Movement Models (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2021
  #'  ============================================
  #'  Script to run hidden Markov movement models for deer, elk, cougars, wolves, 
  #'  coyotes, and bobcats for summer 2018/2019 & winter 2018-2019, respectively. 
  #'  Data were collected & generously provided by WPPP collaborators including
  #'  T.Ganz, T.Roussin, L.Satterfield, B.Windell, and others. Code adapted from
  #'  momentuHMM GitHub, J.Merkel Movement Workshop, L.Satterfield, & R.Emmet.
  #'  Time periods and covariates to match up with single-season occupancy models.
  #'  
  #'  Cleaned telemetry and covariate data were prepared for HMMs with the
  #'  Collar_Movement_DataPrep.R script which took FOREVER to run so only due once.
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(momentuHMM)
  library(rgdal)
  library(tidyverse)

  #'  Load crwOut & covaraite data
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_2021-05-03.RData")
  load("./Outputs/Telemetry_covs/spp_telem_covs_2021-05-10.RData")
  # load("./Outputs/Telemetry_covs/coy_telem_covs_smr.RData")
  # load("./Outputs/Telemetry_covs/coy_telem_covs_wtr.RData")
  

  #'  Merge datasets and create momentuHMMData object
  spp_dataPrep <- function(crwOut, telem_covs){
    #'  Merge crawlOut data with extracted covariate data
    crwlMerge <- crawlMerge(crwOut, telem_covs, Time.name = "time")
    #'  Make categorical variables factors
    crwlMerge$crwPredict$Area <- as.factor(crwlMerge$crwPredict$Area)
    crwlMerge$crwPredict$Sex <- as.factor(crwlMerge$crwPredict$Sex)
    crwlMerge$crwPredict$Year <- as.factor(crwlMerge$crwPredict$Year)
    crwlMerge$crwPredict$Season <- as.factor(crwlMerge$crwPredict$Season)
    #'  Standardize continuous variables
    crwlMerge$crwPredict$Elev <- scale(crwlMerge$crwPredict$Elev)
    crwlMerge$crwPredict$Slope <- scale(crwlMerge$crwPredict$Slope)
    crwlMerge$crwPredict$HumanMod <- scale(crwlMerge$crwPredict$HumanMod)
    crwlMerge$crwPredict$NearestRd <- scale(crwlMerge$crwPredict$NearestRd)
    crwlMerge$crwPredict$PercForMix <- scale(crwlMerge$crwPredict$PercForMix)
    crwlMerge$crwPredict$PercXGrass <- scale(crwlMerge$crwPredict$PercXGrass)
    crwlMerge$crwPredict$PercXShrub <- scale(crwlMerge$crwPredict$PercXShrub)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Elev", "Slope", "HumanMod", "NearestRd", 
                                                   "PercForMix", "PercXGrass", "PercXShrub", 
                                                   "Year", "Sex", "Area", "Season"))
    return(Data)
  }
  #'  Run season & species-specific data through prep function
  #'  Warnings are due to missing Sex data for interpolated locations
  mdData_smr <- spp_dataPrep(crwOut_ALL[[1]], spp_telem_covs[[1]])
  mdData_wtr <- spp_dataPrep(crwOut_ALL[[2]], spp_telem_covs[[2]])
  elkData_smr <- spp_dataPrep(crwOut_ALL[[3]], spp_telem_covs[[3]])
  elkData_wtr <- spp_dataPrep(crwOut_ALL[[4]], spp_telem_covs[[4]])
  wtdData_smr <- spp_dataPrep(crwOut_ALL[[5]], spp_telem_covs[[5]])
  wtdData_wtr <- spp_dataPrep(crwOut_ALL[[6]], spp_telem_covs[[6]])
  cougData_smr <- spp_dataPrep(crwOut_ALL[[7]], spp_telem_covs[[7]])
  cougData_wtr <- spp_dataPrep(crwOut_ALL[[8]], spp_telem_covs[[8]])
  wolfData_smr <- spp_dataPrep(crwOut_ALL[[9]], spp_telem_covs[[9]])
  wolfData_wtr <- spp_dataPrep(crwOut_ALL[[10]], spp_telem_covs[[10]])
  bobData_smr <- spp_dataPrep(crwOut_ALL[[11]], spp_telem_covs[[11]])
  bobData_wtr <- spp_dataPrep(crwOut_ALL[[12]], spp_telem_covs[[12]])
  coyData_smr <- spp_dataPrep(crwOut_ALL[[13]], spp_telem_covs[[13]])
  coyData_wtr <- spp_dataPrep(crwOut_ALL[[14]], spp_telem_covs[[14]])
  # tst <- mapply(spp_dataPrep, crwOut_ALL, spp_telem_covs)
  
  
  #'  Visualize data to inform initial parameter specifications
  # plot(mdData_smr)  #250, 500, 250, 500
  # plot(elkData_smr)  #500, 1000, 500, 1000
  # plot(wtdData_smr)  #100, 500, 100, 500
  # plot(cougData_smr)  #500, 1500, 500, 1500
  # plot(wolfData_smr)  #500, 3000, 500, 3000
  # plot(bobData_smr)  #500, 1000, 500, 1000
  # plot(coyData_smr)  #500, 2000, 500, 2000
  
  
  
  ####  Initial model set up  ####
  
  #'  Define initial parameters associated with each distribution & each state
  #'  Species-specific parameters based on viewing plotted data
  Par0_m1_md <- list(step = c(250, 500, 250, 500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_elk <- list(step = c(500, 1000, 500, 1000, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_wtd <- list(step = c(100, 500, 100, 500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_coug <- list(step = c(500, 1500, 500, 1500, 0.01, 0.005), angle = c(0.3, 0.7))  #zero-mass params needed
  Par0_m1_wolf <- list(step = c(500, 3000, 500, 3000), angle = c(0.3, 0.7))  
  Par0_m1_bob <- list(step = c(500, 1000, 500, 1000), angle = c(0.3, 0.7))  
  Par0_m1_coy <- list(step = c(500, 2000, 500, 2000), angle = c(0.3, 0.7))  
  #'  Step arguments: report 2 means then the 2 SD for the two different states
  #'  Gamma distribution: mean & standard deviation of step lengths for each state
  #'  Michelot & Langrock 2019 recommend using same value for mean and SD per state
  #'  Wrapped Cauchy distribution: concentration of turning angles for each state
  #'  Include zero-mass parameters when there are 0s in the data w/gamma, Weibull, etc. distributions
  #'  e.g., zeromass0 <- c(0.1,0.05) # step zero-mass
  # For md:Par0_m1_md <- list(step = c(250, 1000, 250, 500, 0.01, 0.005), angle = c(0.3, 0.7))  # Last 2 parameters are for zero-mass
  # For wolf:Par0_m1 <- list(step = c(1000, 3000, 1000, 3000), angle = c(0.3, 0.7))  #  No zero-mass parameters needed
  # For elk: Par0_m1 <- list(step = c(250, 1000, 250, 500, 0.01, 0.005), angle = c(0.3, 0.7)) 
  
  #'  Label states
  stateNames <- c("encamped", "exploratory")
  
  #' Distributions for observation processes
  dists <- list(step = "gamma", angle = "wrpcauchy")  
  # dist2 <- list(step = "weibull", angle = "wrpcauchy")
  #' Can test out different distributions 
  #' Step length: gamma or Weibull; Turning angle: von Mises or wrapped Cauchy
  #' State dwell time: geometric distribution
  #' Weibull = "weibull"; von Mises = "vm"
  
  #'  Define formula to be applied to transition probabilities
  pred_formula <- ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod + Area
  prey_formula <- ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearestRd + HumanMod
  #'  Covariates affecting probability of transitioning from one state to the other
  #'  Same covariates as on psi in occupancy models
  #'  Predators have study area included since these vary by individual
  #'  Prey lack study area b/c ungulate collars restricted to 1 study area and all female
  
  #'  Include sex on the state-dependent distributions for predators only
  #'  Males tend to have larger home ranges than females (at least felids) and 
  #'  this likely influences their movement behavior
  #'  Only female ungulates collared so not necessary for prey models
  #'  Add zeromass = ~Sex if needed
  pred_DM <- list(step = list(mean = ~Sex, sd = ~Sex), angle = list(concentration = ~1)) 
  # pred_DM <- list(step = list(mean = ~Sex, sd = ~Sex))
  # pred_DM <- list(step = list(mean = ~Sex, sd = ~1))
  # pred_DM <- list(step = matrix(c(1,0,0,0,"Sex",0,0,0,0,1,0,0,0,"Sex",0,0,0,0,1,0,0,0,0,1), 4, 6))
  #'  Matrix format: repeat mean1, mean2, sd1, sd2 minimum of 4 times if using
  #'  intercept-only (intercept = 1); add rows for each additional covariate 
  #'  where covariates are placed in columns corresponding to each parameter
  
  prey_DM <- NULL

  
  
  ####  It's H[a]MM[er] time!  ####
  
  #'  Keep in mind I can fit covariates on the state transition probabilities, 
  #'  meaning the variables that influence whether an animal will transition from
  #'  one state to the other, or on the state-dependent observation distributions,
  #'  meaning variables that influence step length and/or turning angle for each
  #'  of the states. Currently fitting covariates to state transition probabilities.

  #'  Use retryFits argument to specify the number of attempts to minimize the 
  #'  negative log-likelihood based on random perturbations of the parameter 
  #'  estimates at the current minimum- helps ensure convergence
  
  #'  Function to run data through null and global HMM for each species
  HMM_fit <- function(Data, Par0_m1, pformula, dm) {
    
    #' Fit basic model with no covariates
    m1 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m1,
                 estAngleMean = list(angle = FALSE), stateNames = stateNames)

    #'  Compute the most likely state sequence
    states <- viterbi(m1)
    #'  Derive percentage of time spent in each state
    table(states)/nrow(Data)
    
    #'  Get new initial parameter values for global model based on nested m1 model
    Par0_m2 <- getPar0(model = m1, formula = pformula)  
    
    #'  Fit model with sex covariate on transition probability
    m2 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m2$Par,
                 beta0 = Par0_m2$beta, stateNames = stateNames, formula = pformula, 
                 DM = dm)
    
    #'  What proportion of the locations fall within each state?
    states <- viterbi(m2)
    print(table(states)/nrow(Data))
    
    #'  Model selection with AIC
    print(AIC(m1,m2))
    
    #'  Model summary and covariate effects
    print(m2)
    
    global_est <- CIbeta(m2, alpha = 0.95)
    print(global_est[[3]])
    
    return(m2)
    
  }
  
  #'  Run species-specific data through function
  md_HMM_smr <- HMM_fit(mdData_smr, Par0_m1_md, prey_formula, prey_DM) 
  md_HMM_wtr <- HMM_fit(mdData_wtr, Par0_m1_md, prey_formula, prey_DM) 
  elk_HMM_smr <- HMM_fit(elkData_smr, Par0_m1_elk, prey_formula, prey_DM) 
  elk_HMM_wtr <- HMM_fit(elkData_wtr, Par0_m1_elk, prey_formula, prey_DM) 
  # In fitHMM.momentuHMMData(data = Data, nbStates = 2, dist = dist,  :
  #                            ginv of the hessian failed -- Error in svd(X): infinite or missing values in 'x'
  wtd_HMM_smr <- HMM_fit(wtdData_smr, Par0_m1_wtd, prey_formula, prey_DM) 
  wtd_HMM_wtr <- HMM_fit(wtdData_wtr, Par0_m1_wtd, prey_formula, prey_DM) 
  coug_HMM_smr <- HMM_fit(cougData_smr, Par0_m1_coug, pred_formula, prey_DM) #pred_DM
  coug_HMM_wtr <- HMM_fit(cougData_wtr, Par0_m1_coug, pred_formula, prey_DM) #pred_DM
  wolf_HMM_smr <- HMM_fit(wolfData_smr, Par0_m1_wolf, pred_formula, prey_DM) #pred_DM
  wolf_HMM_wtr <- HMM_fit(wolfData_wtr, Par0_m1_wolf, pred_formula, prey_DM) #pred_DM
  bob_HMM_smr <- HMM_fit(bobData_smr, Par0_m1_bob, pred_formula, prey_DM) #pred_DM
  bob_HMM_wtr <- HMM_fit(bobData_wtr, Par0_m1_bob, pred_formula, prey_DM) #pred_DM
  coy_HMM_smr <- HMM_fit(coyData_smr, Par0_m1_coy, pred_formula, prey_DM) #pred_DM
  coy_HMM_wtr <- HMM_fit(coyData_wtr, Par0_m1_coy, pred_formula, prey_DM) #pred_DM
  
  #'  Save model results
  spp_HMM_output <- list(md_HMM_smr, md_HMM_wtr, elk_HMM_smr, elk_HMM_wtr, wtd_HMM_smr, 
                         wtd_HMM_wtr, coug_HMM_smr, coug_HMM_wtr, wolf_HMM_smr, 
                         wolf_HMM_wtr, bob_HMM_smr, bob_HMM_wtr, coy_HMM_smr, coy_HMM_wtr)
  save(spp_HMM_output, file = paste0("./Outputs/spp_HMM_output_", Sys.Date(), ".RData"))
  
  #'  Calculate transition probabilities
  # md_trProbs_smr <- getTrProbs(md_HMM_smr, getCI=TRUE)
  # md_trProbs_wtr <- getTrProbs(md_HMM_wtr, getCI=TRUE)
  # elk_trProbs_smr <- getTrProbs(elk_HMM_smr, getCI=TRUE)
  # elk_trProbs_wtr <- getTrProbs(elk_HMM_wtr, getCI=TRUE)
  # wtd_trProbs_smr <- getTrProbs(wtd_HMM_smr, getCI=TRUE)
  # wtd_trProbs_wtr <- getTrProbs(wtd_HMM_wtr, getCI=TRUE)
  # coug_trProbs_smr <- getTrProbs(coug_HMM_smr, getCI=TRUE)
  # coug_trProbs_wtr <- getTrProbs(coug_HMM_wtr, getCI=TRUE)
  # wolf_trProbs_smr <- getTrProbs(wolf_HMM_smr, getCI=TRUE)
  # wolf_trProbs_wtr <- getTrProbs(wolf_HMM_wtr, getCI=TRUE)
  # bob_trProbs_smr <- getTrProbs(bob_HMM_smr, getCI=TRUE)
  # bob_trProbs_wtr <- getTrProbs(bob_HMM_wtr, getCI=TRUE)
  # coy_trProbs_smr <- getTrProbs(coy_HMM_smr, getCI=TRUE)
  # coy_trProbs_wtr <- getTrProbs(coy_HMM_wtr, getCI=TRUE)
  
  #'  Function to extract model outputs
  rounddig <- 2
  hmm_out <- function(mod, spp, season) {
    #'  Extract estimates, standard error, and 95% Confidence Intervals for effect
    #'  of each covariate on transition probabilities
    est_out <- CIbeta(mod, alpha = 0.95)
    beta1.2 <- formatC(round(est_out[[3]]$est[,1], rounddig), rounddig, format="f")
    beta2.1 <- formatC(round(est_out[[3]]$est[,2], rounddig), rounddig, format="f")
    se1.2 <- formatC(round(est_out[[3]]$se[,1], rounddig), rounddig, format="f")
    se2.1 <- formatC(round(est_out[[3]]$se[,2], rounddig), rounddig, format="f")
    lci1.2 <- formatC(round(est_out[[3]]$lower[,1], rounddig), rounddig, format="f")
    lci2.1 <- formatC(round(est_out[[3]]$lower[,2], rounddig), rounddig, format="f")
    uci1.2 <- formatC(round(est_out[[3]]$upper[,1], rounddig), rounddig, format="f")
    uci2.1 <- formatC(round(est_out[[3]]$upper[,2], rounddig), rounddig, format="f")
    #'  Merge into a data frame and organize
    out1.2 <- as.data.frame(cbind(beta1.2, se1.2, lci1.2, uci1.2)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Transition = rep("Trans.1->2", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta1.2) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out1.2) <- c("Species", "Season", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out2.1 <- as.data.frame(cbind(beta2.1, se2.1, lci2.1, uci2.1)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        Transition = rep("Trans.2->1", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta2.1) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out2.1) <- c("Species", "Season", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out <- as.data.frame(rbind(out1.2, out2.1))
    return(out)
  }
  #'  Run each season and species-specific model through function
  md_s1819_hmm <- hmm_out(md_HMM_smr, "Mule Deer", "Summer")
  md_w1820_hmm <- hmm_out(md_HMM_wtr, "Mule Deer", "Winter")
  elk_s1819_hmm <- hmm_out(elk_HMM_smr, "Elk", "Summer")
  elk_w1820_hmm <- hmm_out(elk_HMM_wtr, "Elk", "Winter")
  wtd_s1819_hmm <- hmm_out(wtd_HMM_smr, "White-tailed Deer", "Summer")
  wtd_w1820_hmm <- hmm_out(wtd_HMM_wtr, "White-tailed Deer", "Winter")
  coug_s1819_hmm <- hmm_out(coug_HMM_smr, "Cougar", "Summer")
  coug_w1820_hmm <- hmm_out(coug_HMM_wtr, "Cougar", "Winter")
  wolf_s1819_hmm <- hmm_out(wolf_HMM_smr, "Wolf", "Summer")
  wolf_w1820_hmm <- hmm_out(wolf_HMM_wtr, "Wolf", "Winter")
  bob_s1819_hmm <- hmm_out(bob_HMM_smr, "Bobcat", "Summer")
  bob_w1820_hmm <- hmm_out(bob_HMM_wtr, "Bobcat", "Winter")
  coy_s1819_hmm <- hmm_out(coy_HMM_smr, "Coyote", "Summer")
  coy_w1820_hmm <- hmm_out(coy_HMM_wtr, "Coyote", "Winter")
  
  #'  Gather prey and predator results to put into a single results table
  results_hmm_prey <- rbind(md_s1819_hmm, md_w1820_hmm, elk_s1819_hmm, elk_w1820_hmm, 
                            wtd_s1819_hmm, wtd_w1820_hmm)
  results_hmm_pred <- rbind(bob_s1819_hmm, bob_w1820_hmm, coug_s1819_hmm, coug_w1820_hmm, 
                            coy_s1819_hmm, coy_w1820_hmm, wolf_s1819_hmm, wolf_w1820_hmm)
  
  #'  Spread results so the coefficient effects are easier to compare between 
  #'  transition probabilities and across species
  #'  Ungulates (no study area covariate included)
  results_hmm_wide_prey <- results_hmm_prey %>% 
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    # separate("AreaOK", c("AreaOK (SE)", "AreaOK 95% CI"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev 95% CI"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope 95% CI"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix 95% CI"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass 95% CI"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub 95% CI"), sep = "_") %>%
    separate("NearestRd", c("NearestRd (SE)", "NearestRd 95% CI"), sep = "_") %>%
    separate("HumanMod", c("HumanMod (SE)", "HumanMod 95% CI"), sep = "_") %>%
    mutate(
      AreaOK = rep("NA", nrow(.)),
      AreaCI = rep("NA", nrow(.))
    ) %>%
    relocate(AreaOK, .before = "Elev (SE)") %>%
    relocate(AreaCI, .before = "Elev (SE)") %>%
    arrange(match(Species, c("Mule Deer", "Elk", "White-tailed Deer")))
  names(results_hmm_wide_prey)[names(results_hmm_wide_prey) == "AreaOK"] <- "AreaOK (SE)"
  names(results_hmm_wide_prey)[names(results_hmm_wide_prey) == "AreaCI"] <- "AreaOK 95% CI"

  #'  Predators (study area covariate included)
  results_hmm_wide_pred <- results_hmm_pred %>% 
    # dplyr::select(-z) %>%
    mutate(
      # SE = round(SE, 2),
      SE = paste0("(", SE, ")"),
    ) %>%
    #'  Bold significant variables- doesn't work if continue manipulating data frame
    # condformat(.) %>%
    # rule_text_bold(c(Estimate, SE, Pval), expression = Pval <= 0.05) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("(Intercept)", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("AreaOK", c("AreaOK (SE)", "AreaOK 95% CI"), sep = "_") %>%
    separate("Elev", c("Elev (SE)", "Elev 95% CI"), sep = "_") %>%
    separate("Slope", c("Slope (SE)", "Slope 95% CI"), sep = "_") %>%
    separate("PercForMix", c("PercForMix (SE)", "PercForMix 95% CI"), sep = "_") %>%
    separate("PercXGrass", c("PercXGrass (SE)", "PercXGrass 95% CI"), sep = "_") %>%
    separate("PercXShrub", c("PercXShrub (SE)", "PercXShrub 95% CI"), sep = "_") %>%
    separate("NearestRd", c("NearestRd (SE)", "NearestRd 95% CI"), sep = "_") %>%
    separate("HumanMod", c("HumanMod (SE)", "HumanMod 95% CI"), sep = "_") %>%
    arrange(match(Species, c("Bobcat", "Cougar", "Coyote", "Wolf"))) 
  # arrange(match(Season, c("Summer", "Winter")))
  
  results_hmm_wide <- rbind(results_hmm_wide_pred, results_hmm_wide_prey)
  
  write.csv(results_hmm_wide, paste0("./Outputs/HMM_Results_wide", Sys.Date(), ".csv"))
  
  #' #' Plot estimates and CIs for Pr(exploratory) at each time step
  #' plot(trProbs$est[1,2,], type="l", ylim=c(0,1), ylab="Pr(exploratory)", xlab="t", col=c("#E69F00", "#56B4E9")[coy_HMM_smr$miSum$Par$states])
  #' arrows(1:dim(trProbs$est)[3],
  #'        trProbs$lower[1,2,],
  #'        1:dim(trProbs$est)[3],
  #'        trProbs$upper[1,2,],
  #'        length=0.025, angle=90, code=3, col=c("#E69F00", "#56B4E9")[coy_HMM_smr$miSum$Par$states], lwd=1.3)
  #' abline(h=0.5,lty=2)
  #' 
  #' # proportion of entire time series spent in each state
  #' coy_HMM_smr$miSum$Par$timeInStates
  #' 
  #' # histograms of distance to water by state
  #' par(mfrow=c(2,1))
  #' hist(coy_HMM_smr$miSum$data$dist2sabie[which(coy_HMM_smr$miSum$Par$states==1)],main=stateNames[1],xlab="distance to water (m)")
  #' hist(coy_HMM_smr$miSum$data$dist2sabie[which(coy_HMM_smr$miSum$Par$states==2)],main=stateNames[2],xlab="distance to water (m)")
  
  #'  CURRENT QUESTIONS:
  #'  1. make sure I don't need to consolidate bursts in any way (I think each new
  #'  burst should be assigned a new starting parameter and movement of each track 
  #'  is estimated from there, but not 100% sure)
  #'  2. do I add a random effect for individual? Seems like that might be difficult
  #'  in momentuHMM and the discrete-random effects seem weird to me unless it's 
  #'  just random effect for sex or year but not sure if that's how the discrete
  #'  rnd eff works
  #'  3. do I move sex covariate to the state-dependent distributions instead of
  #'  on transition probabilities... I think that makes more sense but having 
  #'  trouble with design matrix on state-dep. distributions
  
  
  
  #' #'  Testing with one species
  #' #'  Create momentuHMMData object from crwData object and covariates
  #' #'  Missing values due to sex missing from interpolated locations
  #' coyData_wtr <- prepData(data = coyMerge_wtr,
  #'                         covNames = c("Elev", "Slope", "HumanMod", "NearstRd",
  #'                                      "PercForMix", "PercXGrass", "PercXShrub",
  #'                                      "Year", "Sex", "Area"))
  #' # mdData <- prepData(data = mdMerge, covNames = c("Elev", "Slope", "HumanMod", "Sex", "Season"))
  #' # wolfData <- prepData(data = wolfMerge, covNames = c("Elev", "Slope", "HumanMod", "dist2road")) #"Sex"
  #' # wolfData <- prepData(data = crwOut_WOLF, covNames = "Sex") #covNames = c("Elev", "Slope", "HumanMod")
  #' 
  #' #' Fit basic model with no covariates
  #' m1 <- fitHMM(data = coyData_wtr, nbStates = 2, dist = dist, Par0 = Par0_m1_coy_wtr,
  #'              estAngleMean = list(angle = FALSE), stateNames = stateNames)
  #' #'  Compute the most likely state sequence
  #' states <- viterbi(m1)
  #' #'  Derive percentage of time spent in each state
  #' table(states)/nrow(coyData_wtr)
  #' 
  #' #'  Adding complexity
  #' formula <- ~Elev + Slope + PercForMix + PercXGrass + PercXShrub + NearstRd + HumanMod + Area + Sex  # Remove Area + Sex if using ungulate data!
  #' #'  Consider putting sex (predators only) or season on the state-dependent distributions
  #' #'  I could see season influencing step length at the very least
  #' DM <- list(step = list(mean = ~Season, sd = ~Season), angle = list(concentration = ~1)) # zeromass = ~Season
  #' 
  #' #'  Get new initial parameter values based on nested m1 model
  #' Par0_m2_coy_wtr <- getPar0(model = m1, formula = formula)  #DM = DM
  #' Par0_m2_coy_wtr$beta  # should the covariate betas be 0.00000???
  #' 
  #' #'  Fit model with all covariates on transition probability
  #' m2 <- fitHMM(data = coyData_wtr, nbStates = 2, dist = dist, Par0 = Par0_m2_coy_wtr$Par,
  #'              beta0 = Par0_m2_coy_wtr$beta, stateNames = stateNames, formula = formula, DM = DM) #
  #' states <- viterbi(m2)
  #' table(states)/nrow(coyData_wtr)
  #' 
  #' #'  Model selection with AIC
  #' AIC(m1,m2)
  #' 
  #' Par0_m3_coy_wtr <- getPar0(model = m2, formula = formula)
  #' Par0_m3_coy_wtr$beta
  #' 
  #' #'  Plot model- this will plot every individual track....
  #' # plot(m1, plotCI = TRUE)
  


  
  
  