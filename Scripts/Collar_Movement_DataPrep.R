  #'  ============================================
  #'  Movement Data Prep (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2021
  #'  ============================================
  #'  Script to format cleaned GPS location data for subsequent HMMM analyses 
  #'  for deer, elk, cougars, wolves, coyotes, and bobcats for summer 2018 
  #'  (7/1/18 - 9/29/18) and winter 2018-2019 (12/1/18 - 3/1/19), respectively. 
  #'  Data were collected & generously provided by WPPP collaborators including
  #'  T.Ganz, T.Roussin, L.Satterfield, B.Windell, and others. Code adapted from
  #'  momentuHMM GitHub, J.Merkel Movement Workshop, L.Satterfield, & R.Emmet.
  #'  Time periods and covariates to match up with single-season occupancy models.
  #'  
  #'  Telemetry data initially cleaned with Collar_DataCleaning.R script, 
  #'  Collar_Truncating&Filtering.R script, and by T.Ganz & L.Satterfield.
  #'  Covariate data based remotely sensed data available from various sources
  #'  (noted in occupancy model script).
  #'  ============================================
  
  #'  Clear memory
  rm(list=ls())
  
  #'  Load libraries
  library(momentuHMM)
  library(rgdal)
  library(tidyverse)
  
  #'  Source cleaned telemetry data
  source("./Scripts/Collar_Truncating&Filtering.R")
  
  ####  Data preparation  ####
  #'  Select relevant columns
  #'  Keeping version of datetime that have been floored to beginning of hour
  rawELK <- elk_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt) %>%
    arrange(ID, Floordt)
  colnames(rawELK) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  #'  Only keep first track to practice with
  # rawELK1 <- subset(rawELK, ID == unique(ID)[2])
  rawMD <- md_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawMD) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  rawWTD <- wtd_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawWTD) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  rawCOUG <- coug_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawCOUG) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  rawWOLF <- wolf_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawWOLF) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  rawBOB <- bob_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawBOB) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  rawCOY <- coy_gtg %>%
    dplyr::select(ID, FullID, Sex, Season, Longitude, Latitude, Floordt)%>%
    arrange(ID, Floordt)
  colnames(rawCOY) <- c("ID", "FullID", "Sex", "Season", "Long", "Lat", "time")
  

  
  
  #'  Function to covert times to POSIX & make locations spatial for each species
  prep_raw <- function(raw, plotit = TRUE) {
    raw$time <- as.POSIXct(raw$time, tz = "America/Los_Angeles")
    
    #'  Make locations spatial and project to UTM coordinates with study area projection
    llcoord <- SpatialPoints(raw[,5:6], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    utmcoord <- spTransform(llcoord, CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs "))
    #'  Add UTM locations to data frame
    raw$x <- attr(utmcoord, "coords")[,1]
    raw$y <- attr(utmcoord, "coords")[,2]
    
    return(raw)
  }
  
  #'  Run data from each species through data prep function
  rawMD <- prep_raw(rawMD)
  rawELK <- prep_raw(rawELK)
  rawWTD <- prep_raw(rawWTD)
  rawCOUG <- prep_raw(rawCOUG)
  rawWOLF <- prep_raw(rawWOLF)
  rawBOB <- prep_raw(rawBOB)
  rawCOY <- prep_raw(rawCOY)
  
  #'  Quick peak at example data for each species
  plot_collar <- function(raw) {
    #'  Pull out data from first collar only
    first <- raw %>%
      slice_head(n = 1) %>%
      dplyr::select(ID)
    animal1 <- raw[raw$ID == first$ID,]
    #'  Plot the first individual in data set
    ggplot(animal1, aes(x = Long, y = Lat)) + 
      geom_point() + 
      geom_path() + 
      facet_wrap(FullID~Season)
  }
  #'  Plot seasonal locations for one individual
  plot_collar(rawMD)
  plot_collar(rawELK)
  plot_collar(rawWTD)
  plot_collar(rawCOUG)
  plot_collar(rawWOLF)
  plot_collar(rawBOB)
  plot_collar(rawCOY)
  
  #'  Take a closer look at those seasonal locations, esp. the mule deer
  #'  Are they migrating halfway through a season? Any dispersal?
  #'  Function to plot locations from individual animals
  plot_telem <- function(spdf){
    #'  Split out spatial points df by individual animal ID
    ind_animal <- group_split(spdf, spdf$ID)
    #'  Place holder for unique animal ID
    names <- c()
    #'  Empty list to hold individual maps
    plot_list <- list()
    #'  Loop through all animals one at a time to create maps of their locations
    for(i in 1:length(unique(ind_animal))) {
      names <- c(names, unique(as.character(ind_animal[[i]]$ID)))
      plot <- ggplot(ind_animal[[i]], aes(x = Long, y = Lat, color = time)) + 
        geom_point() + 
        geom_path() + 
        #'  Nifty side-by-side plots
        facet_wrap(FullID~Season)
      plot_list[[i]] <- plot
    }
    return(plot_list)
  }
  #'  Run data from each species through
  #'  Looking for obvious movement from summer to winter range (or vice verse) 
  #'  within a single season's worth of locations, suggesting migration or dispersal
  rawMD_maps <- plot_telem(rawMD)
  rawELK_maps <- plot_telem(rawELK)
  rawWTD_maps <- plot_telem(rawWTD)
  rawCOUG_maps <- plot_telem(rawCOUG)
  rawWOLF_maps <- plot_telem(rawWOLF)
  rawBOB_maps <- plot_telem(rawBOB)
  rawCOY_maps <- plot_telem(rawCOY)
  
  #'  Save plots as PDF to review and look for evidence of migration or dispersal
  # pdf("./Outputs/GPSlocs_byseason_maps.pdf")
  # for (i in 1:length(unique(rawMD_maps))) {
  #   print(rawMD_maps[[i]])
  # }
  # for (i in 1:length(unique(rawELK_maps))) {
  #   print(rawELK_maps[[i]])
  # }
  # for (i in 1:length(unique(rawWTD_maps))) {
  #   print(rawWTD_maps[[i]])
  # }
  # for (i in 1:length(unique(rawCOUG_maps))) {
  #   print(rawCOUG_maps[[i]])
  # }
  # for (i in 1:length(unique(rawWOLF_maps))) {
  #   print(rawWOLF_maps[[i]])
  # }
  # for (i in 1:length(unique(rawBOB_maps))) {
  #   print(rawBOB_maps[[i]])
  # }
  # for (i in 1:length(unique(rawCOY_maps))) {
  #   print(rawCOY_maps[[i]])
  # }
  # dev.off()
    
  
  #'  Identify bursts of sequential locations & where there are prolonged gaps
  #'  Source creat.burst.R function to identify bursts
  #'  Requires location data be ordered by ID and time (do this above)
  #'  Script written by J.Merkle & provided at Movement Workshop
  source("./Scripts/creat.burst.R")
  
  #'  Function to run locations for each species through the creat.burst function
  bursts <- function(rawloc) {
    #'  Tmax = 28.25 hours (in seconds) so that locations are still grouped in a 
    #'  single burst if there's a gap of 24hr or less in the data (up to 6 sequential   
    #'  fixes missed with 4 hours on each end) but a new burst if gap is >24 hrs.
    rawloc$burst <- creat.burst(data = rawloc, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)   
    #'  Exclude super short bursts (where burst length is 3 or less locations) 
    #'  because need at least 3 points to get a turning angle 
    loc_burst <- rawloc[rawloc$burst %in% names(table(rawloc$burst))[table(rawloc$burst) >=3],]
    length(unique(loc_burst$burst)) 
    head(loc_burst)
    #'  Add burst value ids to Unique ID column to create unique IDs for each track
    loc_burst$UniqueID <- with(loc_burst, paste0(ID, "_", burst))
    loc_burst <- loc_burst %>%
      #'  Change ID from unique animal ID to a combo animal ID & track number identifier
      #'  This is important for crawlWrap function to ensure it does not interpolate
      #'  large chunks of missing data
      mutate(
        AnimalID = ID,
        ID = UniqueID
      ) %>%
      dplyr::select(-UniqueID)
    return(loc_burst)
  }
  #'  Run each species through the function that identifies bursts in the data
  MD_track <- bursts(rawMD)
  ELK_track <- bursts(rawELK)
  WTD_track <- bursts(rawWTD)
  COUG_track <- bursts(rawCOUG)
  WOLF_track <- bursts(rawWOLF)
  BOB_track <- bursts(rawBOB)
  COY_track <- bursts(rawCOY)
  
  
  #'  Function to interpolate missing fixes based on regular time intervals (4hr)
  #'  Interpolating up to 6 missed fixes (24 hr gaps) within each track/animal
  #'  Use crawlWrap function from momentuHMM where:
  #'   -theta are starting values; crawlWrap defaults to 0 if none are provided
  #'   -fixPar contain all parameter values to be held fixed, if not specified 
  #'    then none are fixed
  #'   -estimated parameters: sigma & beta intercepts... what are these?!
  crwWrp <- function(track) {
    #'  Vector numbering unique tracks
    tracks <- 1:length(unique(track$ID))
    #'  Create list of crawlWrap arguments for each track
    #'  Remember: ID is a unique identifier for animal ID & track number
    theta <- fixPar <- list()
    for(i in unique(track$ID)) {
      theta[[i]] <- c(0, 0)
      fixPar[[i]] <- c(NA, NA)
    }
    #'  Interpolate missing locations within each track
    crwOut <- crawlWrap(obsData = track[which(track$ID %in% unique(track$ID)[tracks]),], 
                        theta = theta, fixPar = fixPar, attempts = 100, 
                        Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
    
    return(crwOUT) # LOOK INTO RUNNING THIS IN PARALLEL
  }
  #'  Interpolate missing locations for each species
  crwOut_MD <- crwWrp(MD_track)
  crwOut_ELK <- crwWrp(ELK_track)
  crwOut_WTD <- crwWrp(WTD_track)
  crwOut_COUG <- crwWrp(COUG_track)
  crwOut_WOLF <- crwWrp(WOLF_track)
  # crwOut_BOB <- crwWrp(BOB_track)
  # crwOut_COY <- crwWrp(COY_track)
  
  #'  View interpolated data and new data (step length and turning angle)
  md_move <- crwOut_MD[[2]]
  elk_move <- crwOut_ELK[[2]]
  wtd_move <- crwOut_WTD[[2]]
  coug_move <- crwOut_COUG[[2]]
  wolf_move <- crwOut_WOLF[[2]]
  # bob_move <- crwOut_BOB[[2]]
  # coy_move <- crwOut_COY[[2]]
  
  #'  Save individual crwOut datasets
  save(crwOut_MD, file = "./Outputs/crwOut_MD.RData")
  save(crwOut_ELK, file = "./Outputs/crwOut_ELK.RData")
  save(crwOut_WTD, file = "./Outputs/crwOut_WTD.RData")
  save(crwOut_COUG, file = "./Outputs/crwOut_COUG.RData")
  save(crwOut_WOLF, file = "./Outputs/crwOut_WOLF.RData")
  # save(crwOut_BOB, file = "./Outputs/crwOut_BOB.RData")
  # save(crwOut_COY, file = "./Outputs/crwOut_COY.RData")
  
  
  #'  Save workspace so I never need to rerun crawlWrap function again
  # save.image(paste0("Collar_Movement_DataPrep_", Sys.Date(), ".RData"))
  
  
  ####  Run one species through without function  ####
  
  #' # rawELK1$time <- as.POSIXct(rawELK1$time, tz = "America/Los_Angeles")
  #' rawELK$time <- as.POSIXct(rawELK$time, tz = "America/Los_Angeles")
  #' 
  #' #'  Make locations spatial and project to UTM coordinates with study area projection
  #' # llcoord <- SpatialPoints(rawELK1[,5:6], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  #' llcoord <- SpatialPoints(rawELK[,5:6], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  #' utmcoord <- spTransform(llcoord, CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs "))
  #' #'  Add UTM locations to data frame
  #' # rawELK1$x <- attr(utmcoord, "coords")[,1]
  #' # rawELK1$y <- attr(utmcoord, "coords")[,2]
  #' rawELK$x <- attr(utmcoord, "coords")[,1]
  #' rawELK$y <- attr(utmcoord, "coords")[,2]
  #' 
  #' #'  Quick peak
  #' ggplot(rawELK[rawELK$ID == "1397ELK18",], aes(x = Long, y = Lat)) + 
  #'   geom_point() + 
  #'   geom_path() + 
  #'   facet_wrap(FullID~Season)
  #' 
  #' 
  #' #'  Tmax = 28.25 hours (in seconds) so that locations are still grouped in a 
  #' #'  single burst if there's a gap of 24hr or less in the data (up to 6 sequential   
  #' #'  fixes missed with 4 hours on each end) but a new burst if gap is >24 hrs.
  #' # rawELK1$burst <- creat.burst(data = rawELK1, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)
  #' rawELK$burst <- creat.burst(data = rawELK, id = TRUE, id_name = "ID", date_name = "time", Tmax = 87300)   
  #' head(rawELK)
  #' tail(rawELK)
  #' # rawELK <- rawELK1
  #' #'  Summarize length the bursts
  #' res <- rawELK %>% group_by(burst) %>% summarise(Freq = n())
  #' #'  Look at bursts, ordered from shortest to longest
  #' res[order(res$Freq),]
  #' #'  Frequency of counts (5 1-point bursts, 4 2-point bursts, 0 3-point bursts...)
  #' #'  Keep in mind there are multiple bursts per individual animal, some shorter than others
  #' table(res$Freq)
  #' 
  #' #'  Exclude super short bursts (where burst length is 3 or less locations) 
  #' #'  because need at least 3 points to get a turning angle 
  #' ELK_burst <- rawELK[rawELK$burst %in% names(table(rawELK$burst))[table(rawELK$burst) >=3],]
  #' length(unique(ELK_burst$burst)) 
  #' head(ELK_burst)
  #' #'  Add burst value ids to Unique ID column to create unique IDs for each track
  #' ELK_burst$UniqueID <- with(ELK_burst, paste0(ID, "_", burst))
  #' ELK_burst <- ELK_burst %>%
  #'   mutate(
  #'     AnimalID = ID,
  #'     ID = UniqueID
  #'   ) %>%
  #'   dplyr::select(-UniqueID)
  #' #'  Vector for length of unique tracks
  #' tracks <- 1:length(unique(ELK_burst$ID))
  #' 
  #' 
  #' #'  Fit crawl model to interpolate missing fixes
  #' #'  theta are starting values, crawlWrap defaults to 0 if none are provided
  #' #'  fixPar contain all parameter values to be held fixed, if not specified then none are fixed
  #' #'  estimated parameters: sigma & beta intercepts... what are these?!
  #' #'  Vector for length of unique tracks
  #' tracks <- 1:length(unique(ELK_burst$ID))
  #' #'  List of parameters for each individual animal
  #' theta <- fixPar <- list()
  #' for(i in unique(ELK_burst$ID)) {
  #'   theta[[i]] <- c(0, 0)
  #'   fixPar[[i]] <- c(NA, NA)
  #' }
  #' crwOut_ELK <- crawlWrap(ELK_burst, theta = theta, fixPar = fixPar, attempts = 100, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
  #' crwOut_ELK_burst <- crawlWrap(obsData = ELK_burst[which(ELK_burst$ID %in% unique(ELK_burst$ID)[tracks]),], theta = theta, fixPar = fixPar, attempts = 100, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"))
  #' 
  #' tst <- crwOut_ELK[[2]]
  #' tstb <- crwOut_ELK_burst[[2]]
  #' 
  #' # crwOut_ELK1 <- crawlWrap(obsData = rawELK1, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), retryFits = 0)
  #' # crwOut_ELK <- crawlWrap(obsData = ELK_burst[which(ELK_burst$ID %in% unique(ELK_burst$ID)[tracks]),], Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), retryFits = 100)
  #' # crwOut_ELK <- crawlWrap(obsData = ELK_burst, Time.name = "time", timeStep = "4 hours", coord = c("x", "y"), theta = c(0, 0), fixPar = c(NA, NA), initial.state=list(a=c(0,0),P = diag(c(5000 ^ 2, 10 * 3600 ^ 2))), retryFits = 100)
  #' 
  #' head(crwOut_ELK1[[1]])
  #' cW_ELK1 <- crwOut_ELK1[[2]]
  
  
  