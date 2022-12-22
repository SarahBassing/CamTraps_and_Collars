  #'  ============================================
  #'  Resource Selection Functions (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  ============================================
  #'  Script to randomly sample "available" points from the home range of each 
  #'  collared individual and build resource selection function models for each
  #'  species. This will be compared to the occupancy model results to evaluate
  #'  habitat associations derived from different types of data and whether those
  #'  associations vary with sampling method. 
  #'  
  #'  Note, the track data used here exclude collars with poor fix success or were
  #'  deployed less than 20 days within a season of interest. They also exclude 
  #'  locations generated during obvious predator dispersal events or deer 
  #'  migrations.
  #'  
  #'  Cleaned telemetry and covariate data were prepared for with the
  #'  Collar_Movement_DataPrep.R and Collar_Covariate_Extraction.R scripts 
  #'  which take FOREVER to run. 
  #'  ============================================
  
  #'  Clean workspace & load libraries
  rm(list = ls())  

  #'  Load packages for selecting available points
  library(tidyverse)
  library(sf)
  library(sp)
  library(raster)
  library(lme4)
  library(adehabitatHR)
  
  #'  Load telemetry data
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispersal.RData")
  # load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispMig.RData")
  
  #'  Set projection for spatial analyses
  sa_proj <- CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Functions to filter species-specific data by study area
  #'  NE study area collars
  NE_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "NE")
    return(tracks)
  }
  #'  Drop mule deer summer and winter lists since no mulies in NE study area
  noMD <- spp_all_tracks[-c(1:2)]
  #'  Run list through function
  NE_tracks <- lapply(noMD, NE_collars)
  #'  OK study area collars
  OK_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "OK")
    return(tracks)
  }
  #'  Drop elk and white-tailed deer summer and winter lists since no elk or wtd in OK
  noELKWTD <- spp_all_tracks[-c(3:6)]
  #'  Run list through function
  OK_tracks <- lapply(noELKWTD, OK_collars)

  
  #'  Prep lists to generate available locations per individual
  #'  =========================================================
  #'  Function to pull out unique animal IDs
  unq_id <- function(locs) {
    animal <- as.data.frame(locs$AnimalID) %>%
      unique()
    colnames(animal) <- "AnimalID"
    return(animal)
  }
  animal_ID <- lapply(spp_all_tracks, unq_id)
  NE_ID <- lapply(NE_tracks, unq_id)
  OK_ID <- lapply(OK_tracks, unq_id)
  
  #'  Function to split data into list of dataframes by unique animal ID and year
  #'  (Use AnimalID if ignoring year aspect of data)
  split_animal <- function(locs) {
    ind_animal <- group_split(locs, locs$FullID) #FullID so it's by year too (important for landcover data extraction)
    return(ind_animal)
  }
  animal_split <- lapply(spp_all_tracks, split_animal)
  NE_split <- lapply(NE_tracks, split_animal)
  OK_split <- lapply(OK_tracks, split_animal)
  
  
  #'  Generate random "Available" locations for each individual
  #'  =========================================================
  #'  Define availability based on 2nd or 3rd order selection
  
  #' ####  3rd ORDER SELECTION Average 1:20 used:available  ####
  #' #'  ---------------------------------------------------
  #' #'  Function to randomly select "Available" points within each animal's seasonal 
  #' #'  home range (utilization distributions)
  #' avail_pts_3rd <- function(locs, navail, plotit = F) {
  #'   #'  1. Make each animal's locations spatial
  #'   #'  ---------------------------------------------------------
  #'   locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj)
  #'   #'  Plot to make sure step 1 worked
  #'   if(plotit) {
  #'     plot(locs_sp, col = "blue", pch = 19)
  #'   }
  #'   
  #'   #'  2. Create UDs for each animal following Bivariate normal utilization distributions
  #'   #'  ----------------------------------------------------------
  #'   #'  Estimate KDEs for each home range, extend the spatial extent by 1.5 
  #'   #'  when estimating UDs (throws an error about grid being too small to 
  #'   #'  estimate home range otherwise)
  #'   UD <- kernelUD(locs_sp, kern = "bivnorm", extent = 1.5) 
  #'   UD95 <- getverticeshr(UD, 95)
  #'   #'  Plot to make sure step 2 worked
  #'   if(plotit) {
  #'     plot(UD95, border = "darkgreen", col = NA)
  #'   }
  #'   
  #'   #'  3. Randomly select points within each home range
  #'   #'  ------------------------------------------------
  #'   #'  Sampling 20 available points to every 1 used point
  #'   #'  Identify number of used points per individual
  #'   nused <- nrow(locs) 
  #'   #'  Multiply by 20 
  #'   navailable <- nused*navail#20
  #'   #'  Set seed for reproducibility
  #'   # set.seed(2021)
  #'   rndpts <- spsample(UD95, navailable, type = "random")
  #'   #'  Turn them into spatial points
  #'   rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
  #'   #' Plot to make sure step 3 worked
  #'   if(plotit) {
  #'     plot(rndpts_sp, col = "red", pch = 19)
  #'   }
  #'   
  #'   #'  4. Make list of locations non-spatial
  #'   #'  -------------------------------------
  #'   rndpts_df <- as.data.frame(rndpts_sp) 
  #'   ID <- unique(droplevels(locs$AnimalID))
  #'   Season <- unique(locs$Season)
  #'   rndpts_df$ID <- ID
  #'   rndpts_df$Season <- Season
  #'   
  #'   return(rndpts_df)
  #' 
  #' }
  #' #'  Call function
  #' #'  Run list of lists together
  #' # spp_pts <- lapply(animal_split, lapply, avail_pts, T)
  #' #'  Run lists by species and season
  #' md_smr_df <- lapply(animal_split[[1]], avail_pts_3rd, navail = 20, T) #works when migrations are excluded (with _noDispMig.RData)
  #' md_wtr_df <- lapply(animal_split[[2]], avail_pts_3rd, navail = 20, T) #works
  #' elk_smr_df <- lapply(animal_split[[3]], avail_pts_3rd, navail = 20, T) #works
  #' elk_wtr_df <- lapply(animal_split[[4]], avail_pts_3rd, navail = 20, T) #works
  #' wtd_smr_df <- lapply(animal_split[[5]], avail_pts_3rd, navail = 20, T) #works when rando pt removed from 3144WTD18 smr19 (with _noDispMig.RData)
  #' wtd_wtr_df <- lapply(animal_split[[6]], avail_pts_3rd, navail = 20, T) #works
  #' coug_smr_df <- lapply(animal_split[[7]], avail_pts_3rd, navail = 20, T) #works
  #' coug_wtr_df <- lapply(animal_split[[8]], avail_pts_3rd, navail = 20, T) #works
  #' wolf_smr_df <- lapply(animal_split[[9]], avail_pts_3rd, navail = 20, T) #works
  #' wolf_wtr_df <- lapply(animal_split[[10]], avail_pts_3rd, navail = 20, T) #works
  #' bob_smr_df <- lapply(animal_split[[11]], avail_pts_3rd, navail = 20, T) #works
  #' bob_wtr_df <- lapply(animal_split[[12]], avail_pts_3rd, navail = 20, T) #works
  #' coy_smr_df <- lapply(animal_split[[13]], avail_pts_3rd, navail = 20, T) #works
  #' coy_wtr_df <- lapply(animal_split[[14]], avail_pts_3rd, navail = 20, T) #works
  #' 
  #' #'  Convert to dataframes instead of lists
  #' md_smr_df <- do.call(rbind.data.frame, md_smr_df)
  #' md_wtr_df <- do.call(rbind.data.frame, md_wtr_df)
  #' elk_smr_df <- do.call(rbind.data.frame, elk_smr_df)
  #' elk_wtr_df <- do.call(rbind.data.frame, elk_wtr_df)
  #' wtd_smr_df <- do.call(rbind.data.frame, wtd_smr_df)
  #' wtd_wtr_df <- do.call(rbind.data.frame, wtd_wtr_df)
  #' coug_smr_df <- do.call(rbind.data.frame, coug_smr_df)
  #' coug_wtr_df <- do.call(rbind.data.frame, coug_wtr_df)
  #' wolf_smr_df <- do.call(rbind.data.frame, wolf_smr_df)
  #' wolf_wtr_df <- do.call(rbind.data.frame, wolf_wtr_df)
  #' bob_smr_df <- do.call(rbind.data.frame, bob_smr_df)
  #' bob_wtr_df <- do.call(rbind.data.frame, bob_wtr_df)
  #' coy_smr_df <- do.call(rbind.data.frame, coy_smr_df)
  #' coy_wtr_df <- do.call(rbind.data.frame, coy_wtr_df)
  #' 
  #' #'  Gather into one big list
  #' md_available <- list(md_smr_df, md_wtr_df)
  #' elk_available <- list(elk_smr_df, elk_wtr_df)
  #' wtd_available <- list(wtd_smr_df, wtd_wtr_df)
  #' coug_available <- list(coug_smr_df, coug_wtr_df)
  #' wolf_available <- list(wolf_smr_df, wolf_wtr_df)
  #' bob_available <- list(bob_smr_df, bob_wtr_df)
  #' coy_available <- list(coy_smr_df, coy_wtr_df)
  #' 
  #' #'  Save available points based on individual home ranges
  #' save(md_available, file = paste0("./Outputs/RSF_pts/md_available_3rd_", Sys.Date(), ".RData"))
  #' save(elk_available, file = paste0("./Outputs/RSF_pts/elk_available_3rd_", Sys.Date(), ".RData"))
  #' save(wtd_available, file = paste0("./Outputs/RSF_pts/wtd_available_3rd_", Sys.Date(), ".RData"))
  #' save(coug_available, file = paste0("./Outputs/RSF_pts/coug_available_3rd_", Sys.Date(), ".RData"))
  #' save(wolf_available, file = paste0("./Outputs/RSF_pts/wolf_available_3rd_", Sys.Date(), ".RData"))
  #' save(bob_available, file = paste0("./Outputs/RSF_pts/bob_available_3rd_", Sys.Date(), ".RData"))
  #' save(coy_available, file = paste0("./Outputs/RSF_pts/coy_available_3rd_", Sys.Date(), ".RData"))
 
  
  #'  Calculation size of Home Range based on KDEs
  HR_size <- function(locs, plotit = F) {
    #'  1. Make each animal's locations spatial
    #'  ---------------------------------------------------------
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj)
    #'  Estimate KDEs for each home range
    UD <- kernelUD(locs_sp, kern = "bivnorm", extent = 1.5) 
    #'  Get 95% KDE vertices
    UD95 <- getverticeshr(UD, 95)
    #'  Calculate area of 95% KDE home range (original units in meters so area
    #'  expressed in hectares according to adehabitatHR documentation where it
    #'  default with “m” in and “ha” for output)
    UD95_area <- kernel.area(UD, percent = 95)
    #'  Plot to make sure vertices worked
    if(plotit) {
      plot(UD95, border = "darkgreen", col = NA)
    }
    return(UD95_area)
  }
  md_smr_df <- unlist(lapply(animal_split[[1]], HR_size, T)) #works when migrations are excluded (with _noDispMig.RData) 
  md_wtr_df <- unlist(lapply(animal_split[[2]], HR_size, T))
  elk_smr_df <- unlist(lapply(animal_split[[3]], HR_size, T)) 
  elk_wtr_df <- unlist(lapply(animal_split[[4]], HR_size, T))
  wtd_smr_df <- unlist(lapply(animal_split[[5]], HR_size, T)) #works when rando pt removed from 3144WTD18 smr19 (with _noDispMig.RData)
  wtd_wtr_df <- unlist(lapply(animal_split[[6]], HR_size, T))
  coug_smr_df <- unlist(lapply(animal_split[[7]], HR_size, T))
  coug_wtr_df <- unlist(lapply(animal_split[[8]], HR_size, T))
  wolf_smr_df <- unlist(lapply(animal_split[[9]], HR_size, T))
  wolf_wtr_df <- unlist(lapply(animal_split[[10]], HR_size, T))
  bob_smr_df <- unlist(lapply(animal_split[[11]], HR_size, T))
  bob_wtr_df <- unlist(lapply(animal_split[[12]], HR_size, T))
  coy_smr_df <- unlist(lapply(animal_split[[13]], HR_size, T))
  coy_wtr_df <- unlist(lapply(animal_split[[14]], HR_size, T))
  
  #'  Calculate mean area of home range for each species
  mean_HR <- function(HR_area) {
    #'  Mean area in hectares
    mean_HR <- mean(HR_area)
    #'  Divide by 100 to get area in square kilometers
    mean_HR_km <- round(mean_HR/100, 2)
    print(mean_HR_km)
    #'  Approx. distance from HR center to edge if...
    #'  home range was a perfect circle, radius would be...
    print(sqrt(mean_HR_km/pi))
    #'  home range was a perfect square, half of length of one side would be...
    print(sqrt(mean_HR_km)/2)
    #'  Essentially, would we expect >1 camera within a HR based on average 
    #'  distance between cameras and home range size?
    return(mean_HR_km)
  }
  md_smr_muHR <- mean_HR(md_smr_df) #works when migrations are excluded (with _noDispMig.RData)
  md_wtr_muHR <- mean_HR(md_wtr_df)
  elk_smr_muHR <- mean_HR(elk_smr_df)
  elk_wtr_muHR <- mean_HR(elk_wtr_df)
  wtd_smr_muHR <- mean_HR(wtd_smr_df) #works when rando pt removed from 3144WTD18 smr19 (with _noDispMig.RData)
  wtd_wtr_muHR <- mean_HR(wtd_wtr_df)
  coug_smr_muHR <- mean_HR(coug_smr_df)
  coug_wtr_muHR <- mean_HR(coug_wtr_df)
  wolf_smr_muHR <- mean_HR(wolf_smr_df)
  wolf_wtr_muHR <- mean_HR(wolf_wtr_df)
  bob_smr_muHR <- mean_HR(bob_smr_df)
  bob_wtr_muHR <- mean_HR(bob_wtr_df)
  coy_smr_muHR <- mean_HR(coy_smr_df)
  coy_wtr_muHR <- mean_HR(coy_wtr_df)
  
  
  ####  2nd ORDER SELECTION using buffered individual "home range"  ####
  #'  --------------------------------------------------------------
  #'  Load polygons which represent all locations from each individual's home range
  #'  source("./Scripts/Collar_Buffered_HRs.R")
  load("./Outputs/MCPs/MCP_HomeRange_Polygons_allSpp.RData")
  #'  Functions to
  #'    1. Match buffered HR polygon with seasonal location per individual
  #'       HR_poly = annual home range for each animal per study area & includes
  #'       OK list No. 1) md, 4) OK coug, 6) OK wolf, 8) OK bob, 10) OK coy, and 
  #'       NE list No. 2) elk, 3) wtd, 5) NE coug, 7) NE wolf, 9) NE bob, 11) NE coy 
  #'       spp_all_tracks = animal locations for each animal per season and 
  #'       includes 14 lists that are species and season specific, but not study 
  #'       area specific
  #'    2. Generate available locations within the extent of each individuals'
  #'       buffered home range.
  
  #'  Functions to subset lists of HR polygons to match animals in each season
  #'  OKANOGAN DATA 
  subset_poly_OK <- function(locs, poly) {
    
    #'  Filter locations by study area
    locs_OK <- filter(locs, StudyArea == "OK")
    
    #'  Convert list of SpatialPolygonDataFrames into single sf object
    poly_sf_list <- lapply(poly, st_as_sf)
    poly_sf <- do.call(rbind.data.frame, poly_sf_list)
    
    #'  Subset sf object based on individual IDs present in seasonal locations
    shortlist <- subset(poly_sf, FullID %in% locs_OK$FullID)
    
    #'  Convert sf object back to list of SpatialPolygonDataFrames
    sp_poly <- list()
    for(i in 1:nrow(shortlist)){
      sp_poly[[i]] <- as(shortlist[i,], "Spatial")
    }
    
    return(sp_poly)
  }
  #'  Filter HR polygons based on seasonal location data for OKANOGAN study area
  #'  KEEP TRACK OF LIST ORDERS HERE! HR_poly lists noted above
  md_smr_OK_poly <- subset_poly_OK(spp_all_tracks[[1]], poly = HR_poly[[1]])
  md_wtr_OK_poly <- subset_poly_OK(spp_all_tracks[[2]], poly = HR_poly[[1]])
  coug_smr_OK_poly <- subset_poly_OK(spp_all_tracks[[7]], poly = HR_poly[[4]])
  coug_wtr_OK_poly <- subset_poly_OK(spp_all_tracks[[8]], poly = HR_poly[[4]])
  wolf_smr_OK_poly <- subset_poly_OK(spp_all_tracks[[9]], poly = HR_poly[[6]])
  wolf_wtr_OK_poly <- subset_poly_OK(spp_all_tracks[[10]], poly = HR_poly[[6]])
  bob_smr_OK_poly <- subset_poly_OK(spp_all_tracks[[11]], poly = HR_poly[[8]])
  bob_wtr_OK_poly <- subset_poly_OK(spp_all_tracks[[12]], poly = HR_poly[[8]])
  coy_smr_OK_poly <- subset_poly_OK(spp_all_tracks[[13]], poly = HR_poly[[10]])
  coy_wtr_OK_poly <- subset_poly_OK(spp_all_tracks[[14]], poly = HR_poly[[10]])
  
  #'  NORTHEAST DATA
  subset_poly_NE <- function(locs, poly) {
    
    #'  Filter locations to just Northeast study area
    locs_NE <- filter(locs, StudyArea == "NE")
    
    #'  Convert list of SpatialPolygonDataFrames into single sf object
    poly_sf_list <- lapply(poly, st_as_sf)
    poly_sf <- do.call(rbind.data.frame, poly_sf_list)
    
    #'  Subset sf object based on individual IDs present in seasonal locations
    shortlist <- subset(poly_sf, FullID %in% locs_NE$FullID)
    
    #'  Convert sf object back to list of SpatialPolygonDataFrames
    sp_poly <- list()
    for(i in 1:nrow(shortlist)){
      sp_poly[[i]] <- as(shortlist[i,], "Spatial")
    }

    return(sp_poly)
  }
  #'  Filter HR polygons based on seasonal location data for NORTHEAST study area
  #'  KEEP TRACK OF LIST ORDERS HERE! HR_poly lists noted above 
  elk_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[3]], poly = HR_poly[[2]])
  elk_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[4]], poly = HR_poly[[2]])
  wtd_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[5]], poly = HR_poly[[3]])
  wtd_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[6]], poly = HR_poly[[3]])
  coug_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[7]], poly = HR_poly[[5]])
  coug_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[8]], poly = HR_poly[[5]])
  wolf_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[9]], poly = HR_poly[[7]])
  wolf_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[10]], poly = HR_poly[[7]])
  bob_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[11]], poly = HR_poly[[9]])
  bob_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[12]], poly = HR_poly[[9]])
  coy_smr_NE_poly <- subset_poly_NE(spp_all_tracks[[13]], poly = HR_poly[[11]])
  coy_wtr_NE_poly <- subset_poly_NE(spp_all_tracks[[14]], poly = HR_poly[[11]])

  
  #'  Function to randomly select points within each buffered home range
  avail_pts_buff_2nd <- function(locs, navail, buff_poly, plotit = F)  {
    
    #'  Make location data spatial for plotting below
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj) 
    ID <- as.character(droplevels(unique(locs$AnimalID)))
    Season <- unique(locs$Season)

    #'  Identify number of used points per individual
    nused <- nrow(locs)
    #'  Multiply by desired number of available points
    navailable <- nused*navail
    
    #'  Create data frame with AnimalID and Season info based on number of random 
    #'  points per individual
    AnimalID <- as.data.frame(rep(ID, navailable))
    AnimalID$Season <- Season
    colnames(AnimalID) <- c("AnimalID", "Season")
    
    #'  Set seed for reproducibility
    # set.seed(2021)
    #'  Randomly sample "Available" locations within each polygon
    rndpts <- spsample(buff_poly, navailable, type = "random")
    
    #'  Turn them into SpatialPoints
    rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
    #'  Add AnimalID and Season data, making SpatialPointsDataFrame
    rndpts_spdf <- SpatialPointsDataFrame(rndpts_sp, AnimalID)
    
    #' Plot to make sure it worked
    if(plotit) {
      plot(buff_poly, border = "black", col = NA)
      plot(rndpts_spdf, col = "red", pch = 19, cex = 0.70, add = T)
      plot(locs_sp, col = "blue", cex = 0.70, add = T)
    }
    
    #'  Convert list of SpatialPointsDataFrames to list of sf objects
    rndpts_sf_list <- st_as_sf(rndpts_spdf)
    
    return(rndpts_sf_list)
    
  }
  #'  Generate "available" points for each individual per season and study area 
  #'  where extent of availability is based on buffered annual home range of 
  #'  that individual
  md_smr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[1]], buff_poly = md_smr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    #'  Converts list of sf objects to a single data frame with all locations
    do.call(rbind.data.frame, .)
  md_wtr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[2]], buff_poly = md_wtr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  elk_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[1]], buff_poly = elk_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  elk_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[2]], buff_poly = elk_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wtd_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[3]], buff_poly = wtd_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wtd_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[4]], buff_poly = wtd_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coug_smr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[3]], buff_poly = coug_smr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coug_wtr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[4]], buff_poly = coug_wtr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coug_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[5]], buff_poly = coug_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coug_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[6]], buff_poly = coug_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wolf_smr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[5]], buff_poly = wolf_smr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wolf_wtr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[6]], buff_poly = wolf_wtr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wolf_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[7]], buff_poly = wolf_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  wolf_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[8]], buff_poly = wolf_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  bob_smr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[7]], buff_poly = bob_smr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  bob_wtr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[8]], buff_poly = bob_wtr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  bob_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[9]], buff_poly = bob_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  bob_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[10]], buff_poly = bob_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coy_smr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[9]], buff_poly = coy_smr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coy_wtr_2nd_buff_OK_df <- mapply(avail_pts_buff_2nd, locs = OK_split[[10]], buff_poly = coy_wtr_OK_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coy_smr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[11]], buff_poly = coy_smr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  coy_wtr_2nd_buff_NE_df <- mapply(avail_pts_buff_2nd, locs = NE_split[[12]], buff_poly = coy_wtr_NE_poly, navail = 2000, SIMPLIFY = FALSE, T) %>%
    do.call(rbind.data.frame, .)
  
  #'  Merge species-specific data by season
  md_smr_2nd_buff_df <- md_smr_2nd_buff_OK_df
  md_wtr_2nd_buff_df <- md_wtr_2nd_buff_OK_df
  elk_smr_2nd_buff_df <- elk_smr_2nd_buff_NE_df
  elk_wtr_2nd_buff_df <- elk_wtr_2nd_buff_NE_df
  wtd_smr_2nd_buff_df <- wtd_smr_2nd_buff_NE_df
  wtd_wtr_2nd_buff_df <- wtd_wtr_2nd_buff_NE_df
  coug_smr_2nd_buff_df <- rbind(coug_smr_2nd_buff_OK_df, coug_smr_2nd_buff_NE_df)
  coug_wtr_2nd_buff_df <- rbind(coug_wtr_2nd_buff_OK_df, coug_wtr_2nd_buff_NE_df)
  wolf_smr_2nd_buff_df <- rbind(wolf_smr_2nd_buff_OK_df, wolf_smr_2nd_buff_NE_df)
  wolf_wtr_2nd_buff_df <- rbind(wolf_wtr_2nd_buff_OK_df, wolf_wtr_2nd_buff_NE_df)
  bob_smr_2nd_buff_df <- rbind(bob_smr_2nd_buff_OK_df, bob_smr_2nd_buff_NE_df)
  bob_wtr_2nd_buff_df <- rbind(bob_wtr_2nd_buff_OK_df, bob_wtr_2nd_buff_NE_df)
  coy_smr_2nd_buff_df <- rbind(coy_smr_2nd_buff_OK_df, coy_smr_2nd_buff_NE_df)
  coy_wtr_2nd_buff_df <- rbind(coy_wtr_2nd_buff_OK_df, coy_wtr_2nd_buff_NE_df)
  
  #'  Gather into one big list per species
  md_available_2nd_buffHR <- list(md_smr_2nd_buff_df, md_wtr_2nd_buff_df)
  elk_available_2nd_buffHR <- list(elk_smr_2nd_buff_df, elk_wtr_2nd_buff_df)
  wtd_available_2nd_buffHR <- list(wtd_smr_2nd_buff_df, wtd_wtr_2nd_buff_df)
  coug_available_2nd_buffHR <- list(coug_smr_2nd_buff_df, coug_wtr_2nd_buff_df)
  wolf_available_2nd_buffHR <- list(wolf_smr_2nd_buff_df, wolf_wtr_2nd_buff_df)
  bob_available_2nd_buffHR <- list(bob_smr_2nd_buff_df, bob_wtr_2nd_buff_df)
  coy_available_2nd_buffHR <- list(coy_smr_2nd_buff_df, coy_wtr_2nd_buff_df)
  
  #'  Save available points based on individual home ranges
  save(md_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/md_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(elk_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/elk_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(wtd_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/wtd_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(coug_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/coug_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(wolf_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/wolf_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(bob_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/bob_available_2nd_buffHR_", Sys.Date(), ".RData"))
  save(coy_available_2nd_buffHR, file = paste0("./Outputs/RSF_pts/coy_available_2nd_buffHR_", Sys.Date(), ".RData"))
  

  #' ####  2nd ORDER SELECTION using single MCP  ####
  #' #'  ----------------------------------------
  #' #'  Create MCPs based on all locations from all individuals per species per study area
  #' source("./Scripts/Collar_MCPs.R")
  #' 
  #' #'  Function to randomly select "Available" points within MCPs generated from
  #' #'  all individuals of a given species across seasons within a study area.
  #' #'  MPCs created in Collar_MCPs.R script
  #' avail_pts_2nd <- function(locs, mcps, navail, plotit = T) {
  #'   #'  1. Randomly select points within each MCP
  #'   #'  -----------------------------------------
  #'   #'  Identify number of used points per individual
  #'   nused <- nrow(locs) 
  #'   #'  Multiply by desired number of available points 
  #'   navailable <- nused*navail
  #'   #'  Set seed for reproducibility
  #'   # set.seed(2021)
  #'   rndpts <- spsample(mcps, navailable, type = "random")
  #'   #'  Turn them into spatial points
  #'   rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
  #'   #' Plot to make sure step 1 worked
  #'   if(plotit) {
  #'     plot(rndpts_sp, col = "red", pch = 19, cex = 0.70)
  #'   }
  #'   
  #'   #'  2. Make list of locations non-spatial
  #'   #'  -------------------------------------
  #'   rndpts_df <- as.data.frame(rndpts_sp) 
  #'   ID <- unique(droplevels(locs$AnimalID))
  #'   Season <- unique(locs$Season)
  #'   rndpts_df$ID <- ID
  #'   rndpts_df$Season <- Season
  #'   
  #'   return(rndpts_df)
  #' }
  #' #'  Run lists through by species, season, & study area... this is ugly
  #' #'  Keep a close eye on those indices and study area MCPs
  #' #'  Using buffered MCPs with large water bodies masked out   
  #' #'  OKANOGAN COLLARS
  #' md_smr_2nd_OK_df <- lapply(OK_split[[1]], avail_pts_2nd, mcps = md_poly_clip, navail = 200)
  #' md_wtr_2nd_OK_df <- lapply(OK_split[[2]], avail_pts_2nd, mcps = md_poly_clip, navail = 200)
  #' coug_smr_2nd_OK_df <- lapply(OK_split[[3]], avail_pts_2nd, mcps = coug_OK_poly_clip, navail = 200)
  #' coug_wtr_2nd_OK_df <- lapply(OK_split[[4]], avail_pts_2nd, mcps = coug_OK_poly_clip, navail = 200)
  #' wolf_smr_2nd_OK_df <- lapply(OK_split[[5]], avail_pts_2nd, mcps = wolf_OK_poly_clip, navail = 200)
  #' wolf_wtr_2nd_OK_df <- lapply(OK_split[[6]], avail_pts_2nd, mcps = wolf_OK_poly_clip, navail = 200)
  #' bob_smr_2nd_OK_df <- lapply(OK_split[[7]], avail_pts_2nd, mcps = bob_OK_poly_clip, navail = 200)
  #' bob_wtr_2nd_OK_df <- lapply(OK_split[[8]], avail_pts_2nd, mcps = bob_OK_poly_clip, navail = 200)
  #' coy_smr_2nd_OK_df <- lapply(OK_split[[9]], avail_pts_2nd, mcps = coy_OK_poly_clip, navail = 200)
  #' coy_wtr_2nd_OK_df <- lapply(OK_split[[10]], avail_pts_2nd, mcps = coy_OK_poly_clip, navail = 200)
  #' #'  NORTHEAST COLLARS
  #' elk_smr_2nd_NE_df <- lapply(NE_split[[1]], avail_pts_2nd, mcps = elk_poly_clip, navail = 200)
  #' elk_wtr_2nd_NE_df <- lapply(NE_split[[2]], avail_pts_2nd, mcps = elk_poly_clip, navail = 200)
  #' wtd_smr_2nd_NE_df <- lapply(NE_split[[3]], avail_pts_2nd, mcps = wtd_poly_clip, navail = 200)
  #' wtd_wtr_2nd_NE_df <- lapply(NE_split[[4]], avail_pts_2nd, mcps = wtd_poly_clip, navail = 200)
  #' coug_smr_2nd_NE_df <- lapply(NE_split[[5]], avail_pts_2nd, mcps = coug_NE_poly_clip, navail = 200)
  #' coug_wtr_2nd_NE_df <- lapply(NE_split[[6]], avail_pts_2nd, mcps = coug_NE_poly_clip, navail = 200)
  #' wolf_smr_2nd_NE_df <- lapply(NE_split[[7]], avail_pts_2nd, mcps = wolf_NE_poly_clip, navail = 200)
  #' wolf_wtr_2nd_NE_df <- lapply(NE_split[[8]], avail_pts_2nd, mcps = wolf_NE_poly_clip, navail = 200)
  #' bob_smr_2nd_NE_df <- lapply(NE_split[[9]], avail_pts_2nd, mcps = bob_NE_poly_clip, navail = 200)
  #' bob_wtr_2nd_NE_df <- lapply(NE_split[[10]], avail_pts_2nd, mcps = bob_NE_poly_clip, navail = 200)
  #' coy_smr_2nd_NE_df <- lapply(NE_split[[11]], avail_pts_2nd, mcps = coy_NE_poly_clip, navail = 200)
  #' coy_wtr_2nd_NE_df <- lapply(NE_split[[12]], avail_pts_2nd, mcps = coy_NE_poly_clip, navail = 200)
  #' 
  #' #'  Convert to dataframes instead of lists
  #' #'  Okanogan
  #' md_smr_2nd_OK_df <- do.call(rbind.data.frame, md_smr_2nd_OK_df)
  #' md_wtr_2nd_OK_df <- do.call(rbind.data.frame, md_wtr_2nd_OK_df)
  #' coug_smr_2nd_OK_df <- do.call(rbind.data.frame, coug_smr_2nd_OK_df)
  #' coug_wtr_2nd_OK_df <- do.call(rbind.data.frame, coug_wtr_2nd_OK_df)
  #' wolf_smr_2nd_OK_df <- do.call(rbind.data.frame, wolf_smr_2nd_OK_df)
  #' wolf_wtr_2nd_OK_df <- do.call(rbind.data.frame, wolf_wtr_2nd_OK_df)
  #' bob_smr_2nd_OK_df <- do.call(rbind.data.frame, bob_smr_2nd_OK_df)
  #' bob_wtr_2nd_OK_df <- do.call(rbind.data.frame, bob_wtr_2nd_OK_df)
  #' coy_smr_2nd_OK_df <- do.call(rbind.data.frame, coy_smr_2nd_OK_df)
  #' coy_wtr_2nd_OK_df <- do.call(rbind.data.frame, coy_wtr_2nd_OK_df)
  #' #'  Northeast
  #' elk_smr_2nd_NE_df <- do.call(rbind.data.frame, elk_smr_2nd_NE_df)
  #' elk_wtr_2nd_NE_df <- do.call(rbind.data.frame, elk_wtr_2nd_NE_df)
  #' wtd_smr_2nd_NE_df <- do.call(rbind.data.frame, wtd_smr_2nd_NE_df)
  #' wtd_wtr_2nd_NE_df <- do.call(rbind.data.frame, wtd_wtr_2nd_NE_df)
  #' coug_smr_2nd_NE_df <- do.call(rbind.data.frame, coug_smr_2nd_NE_df)
  #' coug_wtr_2nd_NE_df <- do.call(rbind.data.frame, coug_wtr_2nd_NE_df)
  #' wolf_smr_2nd_NE_df <- do.call(rbind.data.frame, wolf_smr_2nd_NE_df)
  #' wolf_wtr_2nd_NE_df <- do.call(rbind.data.frame, wolf_wtr_2nd_NE_df)
  #' bob_smr_2nd_NE_df <- do.call(rbind.data.frame, bob_smr_2nd_NE_df)
  #' bob_wtr_2nd_NE_df <- do.call(rbind.data.frame, bob_wtr_2nd_NE_df)
  #' coy_smr_2nd_NE_df <- do.call(rbind.data.frame, coy_smr_2nd_NE_df)
  #' coy_wtr_2nd_NE_df <- do.call(rbind.data.frame, coy_wtr_2nd_NE_df)
  #' 
  #' #'  Merge species-specific data by season
  #' md_smr_2nd_df <- md_smr_2nd_OK_df
  #' md_wtr_2nd_df <- md_wtr_2nd_OK_df
  #' elk_smr_2nd_df <- elk_smr_2nd_NE_df
  #' elk_wtr_2nd_df <- elk_wtr_2nd_NE_df
  #' wtd_smr_2nd_df <- wtd_smr_2nd_NE_df
  #' wtd_wtr_2nd_df <- wtd_wtr_2nd_NE_df
  #' coug_smr_2nd_df <- rbind(coug_smr_2nd_OK_df, coug_smr_2nd_NE_df)
  #' coug_wtr_2nd_df <- rbind(coug_wtr_2nd_OK_df, coug_wtr_2nd_NE_df)
  #' wolf_smr_2nd_df <- rbind(wolf_smr_2nd_OK_df, wolf_smr_2nd_NE_df)
  #' wolf_wtr_2nd_df <- rbind(wolf_wtr_2nd_OK_df, wolf_wtr_2nd_NE_df)
  #' bob_smr_2nd_df <- rbind(bob_smr_2nd_OK_df, bob_smr_2nd_NE_df)
  #' bob_wtr_2nd_df <- rbind(bob_wtr_2nd_OK_df, bob_wtr_2nd_NE_df)
  #' coy_smr_2nd_df <- rbind(coy_smr_2nd_OK_df, coy_smr_2nd_NE_df)
  #' coy_wtr_2nd_df <- rbind(coy_wtr_2nd_OK_df, coy_wtr_2nd_NE_df)
  #' 
  #' #'  Gather into one big list per species
  #' md_available_2nd <- list(md_smr_2nd_df, md_wtr_2nd_df)
  #' elk_available_2nd <- list(elk_smr_2nd_df, elk_wtr_2nd_df)
  #' wtd_available_2nd <- list(wtd_smr_2nd_df, wtd_wtr_2nd_df)
  #' coug_available_2nd <- list(coug_smr_2nd_df, coug_wtr_2nd_df)
  #' wolf_available_2nd <- list(wolf_smr_2nd_df, wolf_wtr_2nd_df)
  #' bob_available_2nd <- list(bob_smr_2nd_df, bob_wtr_2nd_df)
  #' coy_available_2nd <- list(coy_smr_2nd_df, coy_wtr_2nd_df)
  #' 
  #' #'  Save available points based on individual home ranges
  #' save(md_available_2nd, file = paste0("./Outputs/RSF_pts/md_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(elk_available_2nd, file = paste0("./Outputs/RSF_pts/elk_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(wtd_available_2nd, file = paste0("./Outputs/RSF_pts/wtd_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(coug_available_2nd, file = paste0("./Outputs/RSF_pts/coug_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(wolf_available_2nd, file = paste0("./Outputs/RSF_pts/wolf_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(bob_available_2nd, file = paste0("./Outputs/RSF_pts/bob_available_2nd_200avail_", Sys.Date(), ".RData"))
  #' save(coy_available_2nd, file = paste0("./Outputs/RSF_pts/coy_available_2nd_200avail_", Sys.Date(), ".RData"))


  
  #'  Extract covariate data for each available location
  #'  ==================================================
  #'  This will take awhile!

  #'  Load packages for covariate extraction
  library(sf)
  library(sp)
  library(rgeos)
  library(raster)
  library(parallel)
  library(future.apply)
  library(tidyverse)
  
  #'  Load location data
  #'  3rd Order Selection
  # load("./Outputs/RSF_pts/md_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/elk_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/wtd_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/coug_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/wolf_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/bob_available_2021-06-22.RData")
  # load("./Outputs/RSF_pts/coy_available_2021-06-22.RData")
  #'  2nd Order Selection (buffered HR approach)
  load("./Outputs/RSF_pts/md_available_2nd_buffHR_2022-04-15.RData") 
  load("./Outputs/RSF_pts/elk_available_2nd_buffHR_2022-04-15.RData")
  load("./Outputs/RSF_pts/wtd_available_2nd_buffHR_2022-04-22.RData")
  load("./Outputs/RSF_pts/coug_available_2nd_buffHR_2022-04-15.RData")
  load("./Outputs/RSF_pts/wolf_available_2nd_buffHR_2022-04-15.RData") 
  load("./Outputs/RSF_pts/bob_available_2nd_buffHR_2022-04-15.RData")
  load("./Outputs/RSF_pts/coy_available_2nd_buffHR_2022-04-15.RData")
  #' #'  2nd Order Selection (MCP-based)
  #' load("./Outputs/RSF_pts/md_available_2nd_2022-04-06.RData") 
  #' load("./Outputs/RSF_pts/elk_available_2nd_2022-04-06.RData")
  #' load("./Outputs/RSF_pts/wtd_available_2nd_2022-04-06.RData")
  #' load("./Outputs/RSF_pts/coug_available_2nd_2022-04-06.RData")
  #' load("./Outputs/RSF_pts/wolf_available_2nd_2022-04-06.RData") 
  #' load("./Outputs/RSF_pts/bob_available_2nd_2022-04-06.RData")
  #' load("./Outputs/RSF_pts/coy_available_2nd_2022-04-06.RData")
  
  
  ####  Extract Covariates at Used & Available Locations  ####
  #'  ----------------------------------------------------
  #'  Read in spatial data
  wppp_bound <- st_read("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  #'  Terrain rasters
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m.tif")
  slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect.tif", band = 1)
  #'  Cascadia Biodiveristy Watch rasters & shapefiles
  landcov18 <- raster("./Shapefiles/Cascadia_layers/landcover_2018.tif")
  landcov19 <- raster("./Shapefiles/Cascadia_layers/landcover_2019.tif")
  interp_landcov18 <- raster("./Shapefiles/Cascadia_layers/interpolated_landcover_2018.tif")
  interp_landcov19 <- raster("./Shapefiles/Cascadia_layers/interpolated_landcover_2019.tif")
  formix2prop18 <- raster("./Shapefiles/Cascadia_layers/forestmix2prop_18.tif")
  formix2prop19 <- raster("./Shapefiles/Cascadia_layers/forestmix2prop_19.tif")
  xgrassprop18 <- raster("./Shapefiles/Cascadia_layers/xgrassprop_18.tif")
  xgrassprop19 <- raster("./Shapefiles/Cascadia_layers/xgrassprop_19.tif")
  xshrubprop18 <- raster("./Shapefiles/Cascadia_layers/xshrubprop_18.tif")
  xshrubprop19 <- raster("./Shapefiles/Cascadia_layers/xshrubprop_19.tif")
  rdden <- raster("./Shapefiles/Cascadia_layers/roadsForTaylor/RoadDensity_1km.tif")
  #'  Human density and human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Create raster stack of terrain layers
  terra_stack <- stack(dem, slope)
  
  #'  Create raster stacks of 2018 and  2019 data
  perc_stack18 <- stack(formix2prop18, xgrassprop18, xshrubprop18)
  perc_stack19 <- stack(formix2prop19, xgrassprop19, xshrubprop19)
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  
  #'  Function to make used and available points data a spatial sp object
  #'  Must be a SpatialPointsDF with sp if running in parallel
  #'  No clue why it won't work as an sf object but whatever
  spatial_locs <- function(locs) {
    sp_locs <- SpatialPointsDataFrame(data = locs, coords = locs[,c("x", "y")], proj4string = CRS(sa_proj))
    return(sp_locs)
  }
  #'  Run used & available locations through function
  used_locs <- lapply(spp_all_tracks, spatial_locs)
  #' #'  3rd Order Selection Covariates
  #' md_locs <- lapply(md_available, spatial_locs)
  #' elk_locs <- lapply(elk_available, spatial_locs)
  #' wtd_locs <- lapply(wtd_available, spatial_locs)
  #' coug_locs <- lapply(coug_available, spatial_locs)
  #' wolf_locs <- lapply(wolf_available, spatial_locs)
  #' bob_locs <- lapply(bob_available, spatial_locs)
  #' coy_locs <- lapply(coy_available, spatial_locs)
  #' #'  2nd Order Selection Covariates (MCP)
  #' md_locs <- lapply(md_available_2nd, spatial_locs)
  #' elk_locs <- lapply(elk_available_2nd, spatial_locs)
  #' wtd_locs <- lapply(wtd_available_2nd, spatial_locs)
  #' coug_locs <- lapply(coug_available_2nd, spatial_locs)
  #' wolf_locs <- lapply(wolf_available_2nd, spatial_locs)
  #' bob_locs <- lapply(bob_available_2nd, spatial_locs)
  #' coy_locs <- lapply(coy_available_2nd, spatial_locs)
  
  #'  Function to make available points into a sp object 
  #'  Necessary b/c buffered 2nd order available data are currently structured 
  #'  to be made into an sf object but need these to be an sp object instead
  sp_locs <- function(locs) {
    locs$ID <- locs$AnimalID
    sp_locs <- as(locs, "Spatial")
    return(sp_locs)
  }
  #'  Run 2nd Order Selection Covariates (buffered HR) through function
  md_locs <- lapply(md_available_2nd_buffHR, sp_locs)
  elk_locs <- lapply(elk_available_2nd_buffHR, sp_locs)
  wtd_locs <- lapply(wtd_available_2nd_buffHR, sp_locs)
  coug_locs <- lapply(coug_available_2nd_buffHR, sp_locs)
  wolf_locs <- lapply(wolf_available_2nd_buffHR, sp_locs)
  bob_locs <- lapply(bob_available_2nd_buffHR, sp_locs)
  coy_locs <- lapply(coy_available_2nd_buffHR, sp_locs)
 
  

  #'  COVARIATE EXTRACTION & CALCULATIONS  
  #'  ===========================================
  #'  Takes forever but running in parallel helps 
  
  #'  Monitor time
  start.time <- Sys.time()
  
  #'  Setup script to run in parallel
  #'  Extract covariates for each species at once
  #'  Identify how many cores I want to use
  detectCores(logical = FALSE)
  cl <- parallel::makeCluster(3) 
  #'  Run in parallel on local computer with specified number of cores
  plan(cluster, workers = cl)
  
  #'  Master function to extract and manipulate covaraite data for each species
  cov_extract <- function(locs) {
    #'  Reproject to WGS84 to match unprojected rasters
    locs_reproj <- spTransform(locs, wgs84)
    
    #'  1. Extract data from terrain & anthropogenic rasters
    #'  ----------------------------------------------------
    terrain <- raster::extract(terra_stack, locs_reproj, df = TRUE)
    modified <- raster::extract(HM, locs_reproj, df = TRUE)
    road_den <- raster::extract(rdden, locs, df = TRUE) # note the projection
    #'  Merge into a single data frame of covariates
    join_covs <- terrain %>%
      full_join(road_den, by = "ID") %>%
      full_join(modified, by = "ID") %>%
      transmute(
        obs = ID,
        Elev = round(WPPP_DEM_30m, digits = 2),
        Slope = round(WPPP_slope_aspect, digits = 2),
        RoadDen = round(RoadDensity_1km, digits = 2),
        HumanMod = round(WPPP_gHM, digits = 2)
      ) %>%
      #'  Need to change NAs to 0 for road density (if NA it means there are no
      #'  roads within that 1km pixel and raster pixel was empty)
      mutate(
        RoadDen = ifelse(is.na(RoadDen), 0, RoadDen)
      )
    #'  Make location and animal ID information non-spatial
    #'  Be sure to remove geometry if this is an sf object
    animal <- as.data.frame(locs) #%>%
      #dplyr::select(-geometry)
    #'  Merge animal/time information with covariates
    covs <- as.data.frame(cbind(animal, join_covs))
    
    #'  3. Extract landcover data from within 250m of each location & calculate
    #'     percent habitat type within that buffer
    #'  ------------------------------------------------------------------------
    animal$obs <- as.numeric(seq(1:nrow(animal)))
    #'  Extract percent landcover type using 250m moving window at each camera site
    perc_landcover18 <- raster::extract(perc_stack18, locs, df = TRUE) %>%
      transmute(
        obs = ID,
        PercForestMix2 = round(forestmix2prop_18, 2),
        PercXericGrass = round(xgrassprop_18, 2),
        PercXericShrub = round(xshrubprop_18, 2)
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer18" | Season == "Winter1819")
    perc_landcover19 <- raster::extract(perc_stack19, locs, df = TRUE) %>%
      transmute(
        obs = ID,
        PercForestMix2 = round(forestmix2prop_19, 2),
        PercXericGrass = round(xgrassprop_19, 2),
        PercXericShrub = round(xshrubprop_19, 2)
      ) %>%
      full_join(animal, by = "obs") %>%
      filter(Season == "Summer19" | Season == "Winter1920")
    percHab <- rbind(perc_landcover18, perc_landcover19)
    
    #'  4. Join all covatiates together & clean up for inclusion in HMM
    telem_covs <- covs %>%
      full_join(percHab, by = c("obs", "ID", "Season")) %>%
      transmute(
        ID = ID,
        Season = Season,
        Year = ifelse(Season == "Summer18" | Season == "Winter1819", "Year1", "Year2"),
        Elev = Elev,
        Slope = Slope,
        RoadDen = RoadDen,
        HumanMod = HumanMod,
        # NearestRd = dist2road, 
        PercForMix = PercForestMix2,
        PercXGrass = PercXericGrass,
        PercXShrub = PercXericShrub,
        obs = obs) %>%
      mutate(
        Area = ifelse(grepl("NE", ID), "NE", "OK"),
        Area = ifelse(grepl("MD", ID), "OK", Area),
        Area = ifelse(grepl("EA", ID), "NE", Area),
        Area = ifelse(grepl("ELK", ID), "NE", Area),
        Area = ifelse(grepl("WTD", ID), "NE", Area),
        #'  Indicate whether this location was used = 1 or available = 0
        Used = 0,
        #'  Add weights to used/available locations (used = 1, available = 5000 
        #'  per Fieberg et al. 2021)
        w = 5000
        )
    
    return(telem_covs)
    
  }
  
  #'  Run list of species used & available location data through function in parallel
  #'  This will take AWHILE even in parallel
  #'  Keep track of which order of selection is being extracted here!!!
  used_covs <- future_lapply(used_locs, cov_extract, future.seed = TRUE)
  md_avail_covs <- future_lapply(md_locs, cov_extract, future.seed = TRUE)
  elk_avail_covs <- future_lapply(elk_locs, cov_extract, future.seed = TRUE)
  wtd_avail_covs <- future_lapply(wtd_locs, cov_extract, future.seed = TRUE)
  coug_avail_covs <- future_lapply(coug_locs, cov_extract, future.seed = TRUE)
  wolf_avail_covs <- future_lapply(wolf_locs, cov_extract, future.seed = TRUE)
  bob_avail_covs <- future_lapply(bob_locs, cov_extract, future.seed = TRUE)
  coy_avail_covs <- future_lapply(coy_locs, cov_extract, future.seed = TRUE)
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  Stop running in parallel
  parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  #'  Add study area to wolf data
  #'  No easy way of doing this because ID not associated with WPPP study areas
  #'  Double check lists 9 & 10 are wolf summer & winter data
  wolf_avail_covs[[1]]$Area <- "NE"
  wolf_avail_covs[[1]] <- mutate(wolf_avail_covs[[1]],
                                Area = ifelse(grepl("W61M", ID), "OK", Area),  #double check no "W71F" in here
                                Area = ifelse(grepl("W88M", ID), "OK", Area),
                                Area = ifelse(grepl("W93M", ID), "OK", Area),
                                Area = ifelse(grepl("W94M", ID), "OK", Area))
  wolf_avail_covs[[2]]$Area <- "NE"
  wolf_avail_covs[[2]] <- mutate(wolf_avail_covs[[2]],
                                 Area = ifelse(grepl("W61M", ID), "OK", Area),
                                 Area = ifelse(grepl("W88M", ID), "OK", Area),
                                 Area = ifelse(grepl("W93M", ID), "OK", Area),
                                 Area = ifelse(grepl("W94M", ID), "OK", Area))
  
  
  #'  Merge lists across lists of coordinates & covariate data
  combo_data <- function(pts, covs){
    merged <- mapply(data.frame, pts, covs, SIMPLIFY = FALSE)
    return(merged)
  }
  #'  Run used and available locations and covariates through
  used_dat <- combo_data(used_locs, used_covs)
  #'  2nd order (with MCPs)
  md_avail_dat <- combo_data(md_available_2nd, md_avail_covs)
  elk_avail_dat <- combo_data(elk_available_2nd, elk_avail_covs)
  wtd_avail_dat <- combo_data(wtd_available_2nd, wtd_avail_covs)
  coug_avail_dat <- combo_data(coug_available_2nd, coug_avail_covs)
  wolf_avail_dat <- combo_data(wolf_available_2nd, wolf_avail_covs)
  bob_avail_dat <- combo_data(bob_available_2nd, bob_avail_covs)
  coy_avail_dat <- combo_data(coy_available_2nd, coy_avail_covs)
  #'  2nd order (with buffered home ranges)
  md_avail_dat <- combo_data(md_available_2nd_buffHR, md_avail_covs)
  elk_avail_dat <- combo_data(elk_available_2nd_buffHR, elk_avail_covs)
  wtd_avail_dat <- combo_data(wtd_available_2nd_buffHR, wtd_avail_covs)
  coug_avail_dat <- combo_data(coug_available_2nd_buffHR, coug_avail_covs)
  wolf_avail_dat <- combo_data(wolf_available_2nd_buffHR, wolf_avail_covs)
  bob_avail_dat <- combo_data(bob_available_2nd_buffHR, bob_avail_covs)
  coy_avail_dat <- combo_data(coy_available_2nd_buffHR, coy_avail_covs)
  
  #'  Function to drop unneeded columns from list of used data sets and indicate 
  #'  whether this location was used = 1 or available = 0 and adjust weights for
  #'  used locations
  select_cols <- function(dat) {
    used_skinny <- dat %>%
      dplyr::select(x, y, AnimalID, Season, ID, Season.1, Year, Elev, Slope, RoadDen, 
                    HumanMod, PercForMix, PercXGrass, PercXShrub, obs, Area) %>%
      mutate(
        Used = 1,
        w = 1)
    colnames(used_skinny) <-  c("x", "y", "ID", "Season", "ID.1", "Season.1", 
                                "Year", "Elev", "Slope", "RoadDen", "HumanMod", 
                                "PercForMix", "PercXGrass", "PercXShrub", "obs", 
                                "Area", "Used", "w")
    return(used_skinny)
  }
  #'  Run the list of used locations through function
  used_dat <- lapply(used_dat, select_cols)
  
  #'  Function to drop unneeded columns from list of available data sets when
  #'  generated from the buffered home ranges
  select_avail_cols <- function(dat) {
    avail_skinny <- dat %>%
      #'  Drop unneeded columns
      dplyr::select(geometry, ID, Season, AnimalID, Season.1, Year, Elev, Slope, RoadDen, 
                    HumanMod, PercForMix, PercXGrass, PercXShrub, obs, Area, Used, w) %>%
      #'  Split geometry column (remnant of sf object) into x & y columns
      mutate(geom = gsub(geometry,pattern="(\\))|(\\()|c",replacement = "")) %>%
      tidyr::separate(geom,into=c("x","y"),sep=",") %>%
      #'  Reorder the coordinates to match used datasets
      relocate(x, .before = ID) %>%
      relocate(y, .after = x) %>%
      #'  Remove the unneeded geometry column
      dplyr::select(-c(geometry))
    colnames(avail_skinny) <-  c("x", "y", "ID", "Season", "ID.1", "Season.1",
                                 "Year", "Elev", "Slope", "RoadDen", "HumanMod",
                                 "PercForMix", "PercXGrass", "PercXShrub", "obs",
                                 "Area", "Used", "w")
    return(avail_skinny)
  }
  md_avail_dat <- lapply(md_avail_dat, select_avail_cols)
  elk_avail_dat <- lapply(elk_avail_dat, select_avail_cols)
  wtd_avail_dat <- lapply(wtd_avail_dat, select_avail_cols)
  coug_avail_dat <- lapply(coug_avail_dat, select_avail_cols)
  wolf_avail_dat <- lapply(wolf_avail_dat, select_avail_cols)
  bob_avail_dat <- lapply(bob_avail_dat, select_avail_cols)
  coy_avail_dat <- lapply(coy_avail_dat, select_avail_cols)
  
    
  #'  Save used and available data separately
  #'  Adjust between 3rd, 2nd, and 2nd_buffHR order depending on version of available locs used above
  save(used_dat, file = paste0("./Outputs/RSF_pts/used_dat_", Sys.Date(), ".RData"))
  save(md_avail_dat, file = paste0("./Outputs/RSF_pts/md_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(elk_avail_dat, file = paste0("./Outputs/RSF_pts/elk_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(wtd_avail_dat, file = paste0("./Outputs/RSF_pts/wtd_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(coug_avail_dat, file = paste0("./Outputs/RSF_pts/coug_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(wolf_avail_dat, file = paste0("./Outputs/RSF_pts/wolf_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(bob_avail_dat, file = paste0("./Outputs/RSF_pts/bob_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))
  save(coy_avail_dat, file = paste0("./Outputs/RSF_pts/coy_avail_2nd_buffHR_dat_", Sys.Date(), ".RData"))

  
  #'  Sub-sample data sets to fewer available locations (20, 200, 1000) to evaluate 
  #'  what a sufficient number of locations is for RSF coefficients to stabilize.
  subsample_avial <- function(used, avail, navail) {
    #'  Merge all used and available locations together
    full_dat <- rbind(used, avail) %>%  
      arrange(ID, Season, Year) %>%
      #'  Generate new ID for each grouping of locations based on observations 
      #'  from each unique individual, season, and year
      group_by(ID, Season, Year) %>%
      mutate(NewID = cur_group_id()) %>%
      ungroup()
    #'  Count number of used locations per individual/season/year
    used_locs <- full_dat %>%
      filter(Used == "1") %>%
      group_by(NewID) %>%
      summarise(n = n()) %>%
      ungroup()
    #'  Pull out available locations, included NewID column
    avail_locs <- filter(full_dat, Used == "0")
    #'  Sub-sample number of available locations based on the number of used locations
    #'  per individual/season/year and a desired ratio of used:available locations
    sub_locs <- data.frame()
    for(i in 1:nrow(used_locs)) {
      #'  Grab the number of used locations for each individual
      nused <- used_locs$n[i]
      #'  Subset available locations to that same individual
      avail_pts <- subset(avail_locs, NewID == used_locs$NewID[i])
      #'  Sub-sample available locations based on ratio of used:available locations
      subsamp <- sample_n(avail_pts, nused*navail)
      #'  Append sub-sampled available locations together
      sub_locs <- rbind(sub_locs, subsamp) 
    }
    #'  Drop unneeded columns and merge with used locations
    sub_locs <- dplyr::select(sub_locs, -NewID)
    useavail <- rbind(sub_locs, used) %>%
      dplyr::select(-c("Season.1", "ID.1"))
    
    return(useavail)
  }
  #'  Sub-sample summer available points 10 times for species with smallest (wtd) 
  #'  & largest (wolf) home ranges for varying number of available points
  library(purrr)
  wtd_smr_dat_all_20 <- rerun(10, subsample_avial(used_dat[[5]], wtd_avail_dat[[1]], navail = 20))
  wtd_smr_dat_all_200 <- rerun(10, subsample_avial(used_dat[[5]], wtd_avail_dat[[1]], navail = 200))
  wtd_smr_dat_all_1000 <- rerun(10, subsample_avial(used_dat[[5]], wtd_avail_dat[[1]], navail = 1000))
  wolf_smr_dat_all_20 <- rerun(10, subsample_avial(used_dat[[9]], wolf_avail_dat[[1]], navail = 20))
  wolf_smr_dat_all_200 <- rerun(10, subsample_avial(used_dat[[9]], wolf_avail_dat[[1]], navail = 200))
  wolf_smr_dat_all_1000 <- rerun(10, subsample_avial(used_dat[[9]], wolf_avail_dat[[1]], navail = 1000))
  
  #' #'  Sub-sample data sets to fewer available locations (20, 200, 1000) to evaluate 
  #' #'  what a sufficient number of locations is for RSF coefficients to stabilize.
  #' library(purrr)
  #' 
  #' subsample_avial <- function(avail, navail) {
  #'   sub_locs <- avail %>%
  #'     group_by(ID, Season, Year) %>%
  #'     sample_n(., navail) %>%
  #'     ungroup()
  #'   return(sub_locs)
  #' }
  #' #'  Sub-sample summer available points 10 times for species with smallest (wtd) 
  #' #'  & largest (wolf) home ranges for varying number of available points
  #' wtd_smr_avail_dat_20 <- rerun(10, subsample_avial(wtd_avail_dat[[1]], navail = 20))
  #' wtd_smr_avail_dat_200 <- rerun(10, subsample_avial(wtd_avail_dat[[1]], navail = 200))
  #' wtd_smr_avail_dat_1000 <- rerun(10, subsample_avial(wtd_avail_dat[[1]], navail = 1000))
  #' wolf_smr_avail_dat_20 <- rerun(10, subsample_avial(wolf_avail_dat[[1]], navail = 20))
  #' wolf_smr_avail_dat_200 <- rerun(10, subsample_avial(wolf_avail_dat[[1]], navail = 200))
  #' wolf_smr_avail_dat_1000 <- rerun(10, subsample_avial(wolf_avail_dat[[1]], navail = 1000))
  #' 
  #' #'  Merge 10 sub-samples of available locations with used locations
  #' merge_use_avail <- function(used, avail) {
  #'   full_dat <- rbind(used, avail) %>%
  #'     arrange(ID, Season, Used) %>%
  #'     dplyr::select(-c("Season.1", "ID.1"))
  #'   return(full_dat)
  #' }
  #' #'  Bind same used locations to each resample of available points
  #' wtd_smr_dat_all_20 <- lapply(wtd_smr_avail_dat_20, merge_use_avail, used = used_dat[[5]])
  #' wtd_smr_dat_all_200 <- lapply(wtd_smr_avail_dat_200, merge_use_avail, used = used_dat[[5]])
  #' wtd_smr_dat_all_1000 <- lapply(wtd_smr_avail_dat_1000, merge_use_avail, used = used_dat[[5]])
  #' wolf_smr_dat_all_20 <- lapply(wolf_smr_avail_dat_20, merge_use_avail, used = used_dat[[9]])
  #' wolf_smr_dat_all_200 <- lapply(wolf_smr_avail_dat_200, merge_use_avail, used = used_dat[[9]])
  #' wolf_smr_dat_all_1000 <- lapply(wolf_smr_avail_dat_1000, merge_use_avail, used = used_dat[[9]])
  
  #'  Save
  save(wtd_smr_dat_all_20, file = paste0("./Outputs/RSF_pts/wtd_smr_dat_all_20_", Sys.Date(), ".RData"))
  save(wtd_smr_dat_all_200, file = paste0("./Outputs/RSF_pts/wtd_smr_dat_all_200_", Sys.Date(), ".RData"))
  save(wtd_smr_dat_all_1000, file = paste0("./Outputs/RSF_pts/wtd_smr_dat_all_1000_", Sys.Date(), ".RData"))
  save(wolf_smr_dat_all_20, file = paste0("./Outputs/RSF_pts/wolf_smr_dat_all_20_", Sys.Date(), ".RData"))
  save(wolf_smr_dat_all_200, file = paste0("./Outputs/RSF_pts/wolf_smr_dat_all_200_", Sys.Date(), ".RData"))
  save(wolf_smr_dat_all_1000, file = paste0("./Outputs/RSF_pts/wolf_smr_dat_all_1000_", Sys.Date(), ".RData"))
  
  
  #'  Save the used and full available data set (2000 sampled locations) 
  wtd_smr_dat_all_2K <- rbind(used_dat[[5]], wtd_avail_dat[[1]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  wolf_smr_dat_all_2K <- rbind(used_dat[[9]], wolf_avail_dat[[1]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  
  save(wtd_smr_dat_all_2K, file = paste0("./Outputs/RSF_pts/wtd_smr_dat_2nd_2Kavail_all_", Sys.Date(), ".RData"))
  save(wolf_smr_dat_all_2K, file = paste0("./Outputs/RSF_pts/wolf_smr_dat_2nd_2Kavail_all_", Sys.Date(), ".RData"))
  
  
  #'  Merge used & available data per species with all available locations
  #'  List order: summer used, summer available, winter used, winter available
  md_dat_all <- rbind(used_dat[[1]], md_avail_dat[[1]], used_dat[[2]], md_avail_dat[[2]])  %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  elk_dat_all <- rbind(used_dat[[3]], elk_avail_dat[[1]], used_dat[[4]], elk_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  wtd_dat_all <- rbind(used_dat[[5]], wtd_avail_dat[[1]], used_dat[[6]], wtd_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  coug_dat_all <- rbind(used_dat[[7]], coug_avail_dat[[1]], used_dat[[8]], coug_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  wolf_dat_all <- rbind(used_dat[[9]], wolf_avail_dat[[1]], used_dat[[10]], wolf_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  bob_dat_all <- rbind(used_dat[[11]], bob_avail_dat[[1]], used_dat[[12]], bob_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  coy_dat_all <- rbind(used_dat[[13]], coy_avail_dat[[1]], used_dat[[14]], coy_avail_dat[[2]]) %>%
    arrange(ID, Season, Used) %>%
    dplyr::select(-c("Season.1", "ID.1"))
  
  #'  Save combined data for final RSFs
  #'  Adjust between 3rd, 2nd & 2nd_buffHR order depending on version of available locs used above
  save(md_dat_all, file = paste0("./Outputs/RSF_pts/md_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(elk_dat_all, file = paste0("./Outputs/RSF_pts/elk_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(wtd_dat_all, file = paste0("./Outputs/RSF_pts/wtd_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(coug_dat_all, file = paste0("./Outputs/RSF_pts/coug_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(wolf_dat_all, file = paste0("./Outputs/RSF_pts/wolf_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(bob_dat_all, file = paste0("./Outputs/RSF_pts/bob_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  save(coy_dat_all, file = paste0("./Outputs/RSF_pts/coy_dat_2nd_buffHR_all_", Sys.Date(), ".RData"))
  
  
  