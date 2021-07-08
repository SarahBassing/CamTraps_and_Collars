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
  
  #'  Load packages for selecting available points
  library(tidyverse)
  library(sp)
  library(raster)
  library(lme4)
  library(adehabitatHR)
  
  #'  Load MCPs which represent all locations from all individuals of a species
  #'  Ungulates have 1 MCP each; carnivores have 1 MCP per study area
  # load("./Outputs/RSF_pts/MCP_all_ind.RData")
  source("./Scripts/Collar_MCPs.R")
  
  #'  Load telemetry data
  # load("./Outputs/Telemetry_tracks/spp_all_tracks.RData")
  # load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispersal.RData")
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispMig.RData")
  
  
  #'  Functions to filter species-specific data by study area
  #'  NE study area collars
  NE_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "NE")
    return(tracks)
  }
  #'  Drop mule deer summer and winter lists
  noMD <- spp_all_tracks[-c(1:2)]
  #'  Run list through function
  NE_tracks <- lapply(noMD, NE_collars)
  #'  OK study area collars
  OK_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "OK")
    return(tracks)
  }
  #'  Drop elk and white-tailed deer summer and winter lists
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
  
  #'  Function to calculate mean number of observations per species
  #'  Used to decide how many "available" points to draw for each individual
  #'  Using average across all individuals & species to be consistent across models
  #'  Create empty dataframe to hold mean locations
  avg_locs <- c()
  mean_obs <- function(locs) {
    mean_locs <- (nrow(locs))/(length(unique(locs$FullID))) #FullID? #AnimalID
    avg_locs <- c(avg_locs, mean_locs)
    return(avg_locs)
  }
  mean_locs <- lapply(spp_all_tracks, mean_obs)
  # NE_mean_locs <- lapply(NE_tracks, mean_obs)
  # OK_mean_locs <- lapply(OK_tracks, mean_obs)

  #'  Calculate mean number of used locations for all species
  mean_used <- mean(unlist(mean_locs)); sd_used <- sd(unlist(mean_locs))
  # NE_mean_used <- mean(unlist(NE_mean_locs)); NE_sd_used <- sd(unlist(NE_mean_locs))
  # OK_mean_used <- mean(unlist(OK_mean_locs)); OK_sd_used <- sd(unlist(OK_mean_locs))
  #'  RSF literature suggests 1:20 ratio used:available
  navailable <- mean_used*20
  # NE_navailable <- NE_mean_used*20
  # OK_navailable <- OK_mean_used*20

  
  
  #'  Generate random "Available" locations for each individual
  #'  =========================================================
  #'  Set projection for spatial analyses
  sa_proj <- CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  3rd ORDER SELECTION
  #'  -------------------
  #'  Function to randomly select "Available" points within each animal's seasonal 
  #'  home range (utilization distributions)
  avail_pts_3rd <- function(locs, plotit = F) {
    #'  1. Make each animal's locations spatial
    #'  ---------------------------------------------------------
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj)
    #'  Plot to make sure step 1 worked
    if(plotit) {
      plot(locs_sp, col = "blue", pch = 19)
    }
    
    #'  2. Create UDs for each animal following Bivariate normal utilization distributions
    #'  ----------------------------------------------------------
    #'  Estimate KDEs for each home range, extend the spatial extent by 1.5 
    #'  when estimating UDs (throws an error about grid being too small to 
    #'  estimate home range otherwise)
    UD <- kernelUD(locs_sp, kern = "bivnorm", extent = 1.5) 
    UD95 <- getverticeshr(UD, 95)
    #'  Plot to make sure step 2 worked
    if(plotit) {
      plot(UD95, border = "darkgreen", col = NA)
    }
    
    #'  3. Randomly select points within each home range
    #'  ------------------------------------------------
    set.seed(2021)
    rndpts <- spsample(UD95, navailable, type = "random")
    #'  Turn them into spatial points
    rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
    #' Plot to make sure step 3 worked
    if(plotit) {
      plot(rndpts_sp, col = "red", pch = 19)
    }
    
    #'  4. Make list of locations non-spatial
    #'  -------------------------------------
    rndpts_df <- as.data.frame(rndpts_sp) 
    ID <- unique(droplevels(locs$AnimalID))
    Season <- unique(locs$Season)
    rndpts_df$ID <- ID
    rndpts_df$Season <- Season
    
    return(rndpts_df)

  }
  #'  Call function
  #'  Run list of lists together
  # spp_pts <- lapply(animal_split, lapply, avail_pts, T)
  #'  Run lists by species and season
  md_smr_df <- lapply(animal_split[[1]], avail_pts_3rd, T) #works when migrations are excluded
  md_wtr_df <- lapply(animal_split[[2]], avail_pts_3rd, T) #works
  elk_smr_df <- lapply(animal_split[[3]], avail_pts_3rd, T) #works
  elk_wtr_df <- lapply(animal_split[[4]], avail_pts_3rd, T) #works
  wtd_smr_df <- lapply(animal_split[[5]], avail_pts_3rd, T) #works when rando pt removed from 3144WTD18 smr19
  wtd_wtr_df <- lapply(animal_split[[6]], avail_pts_3rd, T) #works
  coug_smr_df <- lapply(animal_split[[7]], avail_pts_3rd, T) #works
  coug_wtr_df <- lapply(animal_split[[8]], avail_pts_3rd, T) #works
  wolf_smr_df <- lapply(animal_split[[9]], avail_pts_3rd, T) #works
  wolf_wtr_df <- lapply(animal_split[[10]], avail_pts_3rd, T) #works
  bob_smr_df <- lapply(animal_split[[11]], avail_pts_3rd, T) #works
  bob_wtr_df <- lapply(animal_split[[12]], avail_pts_3rd, T) #works
  coy_smr_df <- lapply(animal_split[[13]], avail_pts_3rd, T) #works
  coy_wtr_df <- lapply(animal_split[[14]], avail_pts_3rd, T) #works
  
  #'  Convert to dataframes instead of lists
  md_smr_df <- do.call(rbind.data.frame, md_smr_df)
  md_wtr_df <- do.call(rbind.data.frame, md_wtr_df)
  elk_smr_df <- do.call(rbind.data.frame, elk_smr_df)
  elk_wtr_df <- do.call(rbind.data.frame, elk_wtr_df)
  wtd_smr_df <- do.call(rbind.data.frame, wtd_smr_df)
  wtd_wtr_df <- do.call(rbind.data.frame, wtd_wtr_df)
  coug_smr_df <- do.call(rbind.data.frame, coug_smr_df)
  coug_wtr_df <- do.call(rbind.data.frame, coug_wtr_df)
  wolf_smr_df <- do.call(rbind.data.frame, wolf_smr_df)
  wolf_wtr_df <- do.call(rbind.data.frame, wolf_wtr_df)
  bob_smr_df <- do.call(rbind.data.frame, bob_smr_df)
  bob_wtr_df <- do.call(rbind.data.frame, bob_wtr_df)
  coy_smr_df <- do.call(rbind.data.frame, coy_smr_df)
  coy_wtr_df <- do.call(rbind.data.frame, coy_wtr_df)
  
  #'  Gather into one big list
  md_available <- list(md_smr_df, md_wtr_df)
  elk_available <- list(elk_smr_df, elk_wtr_df)
  wtd_available <- list(wtd_smr_df, wtd_wtr_df)
  coug_available <- list(coug_smr_df, coug_wtr_df)
  wolf_available <- list(wolf_smr_df, wolf_wtr_df)
  bob_available <- list(bob_smr_df, bob_wtr_df)
  coy_available <- list(coy_smr_df, coy_wtr_df)
  
  #'  Save available points based on individual home ranges
  save(md_available, file = paste0("./Outputs/RSF_pts/md_available_", Sys.Date(), ".RData"))
  save(elk_available, file = paste0("./Outputs/RSF_pts/elk_available_", Sys.Date(), ".RData"))
  save(wtd_available, file = paste0("./Outputs/RSF_pts/wtd_available_", Sys.Date(), ".RData"))
  save(coug_available, file = paste0("./Outputs/RSF_pts/coug_available_", Sys.Date(), ".RData"))
  save(wolf_available, file = paste0("./Outputs/RSF_pts/wolf_available_", Sys.Date(), ".RData"))
  save(bob_available, file = paste0("./Outputs/RSF_pts/bob_available_", Sys.Date(), ".RData"))
  save(coy_available, file = paste0("./Outputs/RSF_pts/coy_available_", Sys.Date(), ".RData"))
 
  
  #'  2nd ORDER SELECTION
  #'  -------------------
  #'  Function to randomly select "Available" points within MCPs generated from
  #'  all individuals of a given species across seasons within a study area
  #'  MPCs created in Collar_MCPs.R script
  avail_pts_2nd <- function(locs, mcps, plotit = F) {
    #'  1. Randomly select points within each MCP
    #'  -----------------------------------------
    set.seed(2021)
    rndpts <- spsample(mcps, navailable, type = "random")
    #'  Turn them into spatial points
    rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
    #' Plot to make sure step 1 worked
    if(plotit) {
      plot(rndpts_sp, col = "red", pch = 19)
    }
    
    #'  2. Make list of locations non-spatial
    #'  -------------------------------------
    rndpts_df <- as.data.frame(rndpts_sp) 
    ID <- unique(droplevels(locs$AnimalID))
    Season <- unique(locs$Season)
    rndpts_df$ID <- ID
    rndpts_df$Season <- Season
    
    return(rndpts_df)
  }
  #'  Run lists through by species, season, & study area... this is ugly
  #'  Keep a close eye on those indices and study area MCPs   
  #'  OKANOGAN COLLARS
  md_smr_2nd_OK_df <- lapply(OK_split[[1]], avail_pts_2nd, mcps = md_mcp)
  md_wtr_2nd_OK_df <- lapply(OK_split[[2]], avail_pts_2nd, mcps = md_mcp)
  coug_smr_2nd_OK_df <- lapply(OK_split[[3]], avail_pts_2nd, mcps = coug_OK_mcp)
  coug_wtr_2nd_OK_df <- lapply(OK_split[[4]], avail_pts_2nd, mcps = coug_OK_mcp)
  wolf_smr_2nd_OK_df <- lapply(OK_split[[5]], avail_pts_2nd, mcps = wolf_OK_mcp)
  wolf_wtr_2nd_OK_df <- lapply(OK_split[[6]], avail_pts_2nd, mcps = wolf_OK_mcp)
  bob_smr_2nd_OK_df <- lapply(OK_split[[7]], avail_pts_2nd, mcps = bob_OK_mcp)
  bob_wtr_2nd_OK_df <- lapply(OK_split[[8]], avail_pts_2nd, mcps = bob_OK_mcp)
  coy_smr_2nd_OK_df <- lapply(OK_split[[9]], avail_pts_2nd, mcps = coy_OK_mcp)
  coy_wtr_2nd_OK_df <- lapply(OK_split[[10]], avail_pts_2nd, mcps = coy_OK_mcp)
  #'  NORTHEAST COLLARS
  elk_smr_2nd_NE_df <- lapply(NE_split[[1]], avail_pts_2nd, mcps = elk_mcp)
  elk_wtr_2nd_NE_df <- lapply(NE_split[[2]], avail_pts_2nd, mcps = elk_mcp)
  wtd_smr_2nd_NE_df <- lapply(NE_split[[3]], avail_pts_2nd, mcps = wtd_mcp)
  wtd_wtr_2nd_NE_df <- lapply(NE_split[[4]], avail_pts_2nd, mcps = wtd_mcp)
  coug_smr_2nd_NE_df <- lapply(NE_split[[5]], avail_pts_2nd, mcps = coug_NE_mcp)
  coug_wtr_2nd_NE_df <- lapply(NE_split[[6]], avail_pts_2nd, mcps = coug_NE_mcp)
  wolf_smr_2nd_NE_df <- lapply(NE_split[[7]], avail_pts_2nd, mcps = wolf_NE_mcp)
  wolf_wtr_2nd_NE_df <- lapply(NE_split[[8]], avail_pts_2nd, mcps = wolf_NE_mcp)
  bob_smr_2nd_NE_df <- lapply(NE_split[[9]], avail_pts_2nd, mcps = bob_NE_mcp)
  bob_wtr_2nd_NE_df <- lapply(NE_split[[10]], avail_pts_2nd, mcps = bob_NE_mcp)
  coy_smr_2nd_NE_df <- lapply(NE_split[[11]], avail_pts_2nd, mcps = coy_NE_mcp)
  coy_wtr_2nd_NE_df <- lapply(NE_split[[12]], avail_pts_2nd, mcps = coy_NE_mcp)
  
  #'  Convert to dataframes instead of lists
  #'  Okanogan
  md_smr_2nd_OK_df <- do.call(rbind.data.frame, md_smr_2nd_OK_df)
  md_wtr_2nd_OK_df <- do.call(rbind.data.frame, md_wtr_2nd_OK_df)
  coug_smr_2nd_OK_df <- do.call(rbind.data.frame, coug_smr_2nd_OK_df)
  coug_wtr_2nd_OK_df <- do.call(rbind.data.frame, coug_wtr_2nd_OK_df)
  wolf_smr_2nd_OK_df <- do.call(rbind.data.frame, wolf_smr_2nd_OK_df)
  wolf_wtr_2nd_OK_df <- do.call(rbind.data.frame, wolf_wtr_2nd_OK_df)
  bob_smr_2nd_OK_df <- do.call(rbind.data.frame, bob_smr_2nd_OK_df)
  bob_wtr_2nd_OK_df <- do.call(rbind.data.frame, bob_wtr_2nd_OK_df)
  coy_smr_2nd_OK_df <- do.call(rbind.data.frame, coy_smr_2nd_OK_df)
  coy_wtr_2nd_OK_df <- do.call(rbind.data.frame, coy_wtr_2nd_OK_df)
  #'  Northeast
  elk_smr_2nd_NE_df <- do.call(rbind.data.frame, elk_smr_2nd_NE_df)
  elk_wtr_2nd_NE_df <- do.call(rbind.data.frame, elk_wtr_2nd_NE_df)
  wtd_smr_2nd_NE_df <- do.call(rbind.data.frame, wtd_smr_2nd_NE_df)
  wtd_wtr_2nd_NE_df <- do.call(rbind.data.frame, wtd_wtr_2nd_NE_df)
  coug_smr_2nd_NE_df <- do.call(rbind.data.frame, coug_smr_2nd_NE_df)
  coug_wtr_2nd_NE_df <- do.call(rbind.data.frame, coug_wtr_2nd_NE_df)
  wolf_smr_2nd_NE_df <- do.call(rbind.data.frame, wolf_smr_2nd_NE_df)
  wolf_wtr_2nd_NE_df <- do.call(rbind.data.frame, wolf_wtr_2nd_NE_df)
  bob_smr_2nd_NE_df <- do.call(rbind.data.frame, bob_smr_2nd_NE_df)
  bob_wtr_2nd_NE_df <- do.call(rbind.data.frame, bob_wtr_2nd_NE_df)
  coy_smr_2nd_NE_df <- do.call(rbind.data.frame, coy_smr_2nd_NE_df)
  coy_wtr_2nd_NE_df <- do.call(rbind.data.frame, coy_wtr_2nd_NE_df)
  
  #'  Merge species-specific data by season
  md_smr_2nd_df <- md_smr_2nd_OK_df
  md_wtr_2nd_df <- md_wtr_2nd_OK_df
  elk_smr_2nd_df <- elk_smr_2nd_NE_df
  elk_wtr_2nd_df <- elk_wtr_2nd_NE_df
  wtd_smr_2nd_df <- wtd_smr_2nd_NE_df
  wtd_wtr_2nd_df <- wtd_wtr_2nd_NE_df
  coug_smr_2nd_df <- rbind(coug_smr_2nd_OK_df, coug_smr_2nd_NE_df)
  coug_wtr_2nd_df <- rbind(coug_wtr_2nd_OK_df, coug_wtr_2nd_NE_df)
  wolf_smr_2nd_df <- rbind(wolf_smr_2nd_OK_df, wolf_smr_2nd_NE_df)
  wolf_wtr_2nd_df <- rbind(wolf_wtr_2nd_OK_df, wolf_wtr_2nd_NE_df)
  bob_smr_2nd_df <- rbind(bob_smr_2nd_OK_df, bob_smr_2nd_NE_df)
  bob_wtr_2nd_df <- rbind(bob_wtr_2nd_OK_df, bob_wtr_2nd_NE_df)
  coy_smr_2nd_df <- rbind(coy_smr_2nd_OK_df, coy_smr_2nd_NE_df)
  coy_wtr_2nd_df <- rbind(coy_wtr_2nd_OK_df, coy_wtr_2nd_NE_df)
  
  #'  Gather into one big list per species
  md_available_2nd <- list(md_smr_2nd_df, md_wtr_2nd_df)
  elk_available_2nd <- list(elk_smr_2nd_df, elk_wtr_2nd_df)
  wtd_available_2nd <- list(wtd_smr_2nd_df, wtd_wtr_2nd_df)
  coug_available_2nd <- list(coug_smr_2nd_df, coug_wtr_2nd_df)
  wolf_available_2nd <- list(wolf_smr_2nd_df, wolf_wtr_2nd_df)
  bob_available_2nd <- list(bob_smr_2nd_df, bob_wtr_2nd_df)
  coy_available_2nd <- list(coy_smr_2nd_df, coy_wtr_2nd_df)
  
  #'  Save available points based on individual home ranges
  save(md_available_2nd, file = paste0("./Outputs/RSF_pts/md_available_2nd_", Sys.Date(), ".RData"))
  save(elk_available_2nd, file = paste0("./Outputs/RSF_pts/elk_available_2nd_", Sys.Date(), ".RData"))
  save(wtd_available_2nd, file = paste0("./Outputs/RSF_pts/wtd_available_2nd_", Sys.Date(), ".RData"))
  save(coug_available_2nd, file = paste0("./Outputs/RSF_pts/coug_available_2nd_", Sys.Date(), ".RData"))
  save(wolf_available_2nd, file = paste0("./Outputs/RSF_pts/wolf_available_2nd_", Sys.Date(), ".RData"))
  save(bob_available_2nd, file = paste0("./Outputs/RSF_pts/bob_available_2nd_", Sys.Date(), ".RData"))
  save(coy_available_2nd, file = paste0("./Outputs/RSF_pts/coy_available_2nd_", Sys.Date(), ".RData"))
  
  
  
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
  #'  2nd Order Selection
  load("./Outputs/RSF_pts/md_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/elk_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/wtd_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/coug_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/wolf_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/bob_available_2nd_2021-07-07.RData")
  load("./Outputs/RSF_pts/coy_available_2nd_2021-07-07.RData")
  
  #'  Read in spatial data
  wppp_bound <- st_read("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  #'  Terrain rasters
  # dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m_reproj.tif")
  # Slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect_reproj.tif", band = 1)
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
  # roads <- st_read("./Shapefiles/Cascadia_layers/roadsForTaylor", layer = "roadsForTaylor")
  rdden <- raster("./Shapefiles/roaddensity/road.density_km2_TIF.tif")
  #'  Human density and human modified landscapes
  # HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM_reproj.tif")
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Create raster stack of terrain layers
  terra_stack <- stack(dem, slope)
  
  #'  Create raster stacks of 2018 and  2019 data
  perc_stack18 <- stack(formix2prop18, xgrassprop18, xshrubprop18)
  perc_stack19 <- stack(formix2prop19, xgrassprop19, xshrubprop19)
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Reproject road shapefile to match animal location projection
  # road_reproj <- st_transform(roads, crs = st_crs(sa_proj))
  # projection(road_reproj)

  
  #'  Function to make available points data a spatial sf object
  #'  Actually make it a SpatialPointsDF with sp if running in parallel
  #'  No clue why it won't work as an sf object but whatever
  spatial_locs <- function(locs) {
    sf_locs <- SpatialPointsDataFrame(data = locs, coords = locs[,c("x", "y")], proj4string = CRS(sa_proj))
    # sf_locs <- st_as_sf(locs, coords = c("x", "y"), crs = sa_proj)
    return(sf_locs)
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
  #'  2nd Order Selection Covariates
  md_locs <- lapply(md_available_2nd, spatial_locs)
  elk_locs <- lapply(elk_available_2nd, spatial_locs)
  wtd_locs <- lapply(wtd_available_2nd, spatial_locs)
  coug_locs <- lapply(coug_available_2nd, spatial_locs)
  wolf_locs <- lapply(wolf_available_2nd, spatial_locs)
  bob_locs <- lapply(bob_available_2nd, spatial_locs)
  coy_locs <- lapply(coy_available_2nd, spatial_locs)
  
  
  ####  KEEP TRACK OF WHICH VERSION OF THESE LOCS DATA I USE BELOW!  ####
  

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
        RoadDen = round(road.density_km2_TIF, digits = 2),
        HumanMod = round(WPPP_gHM, digits = 2)
        # Elev = round(WPPP_DEM_30m_reproj, digits = 2),
        # Slope = round(WPPP_slope_aspect_reproj, digits = 2),
        # RoadDen = round(road.density_km2_TIF, digits = 2),
        # HumanMod = round(WPPP_gHM_reproj, digits = 2)
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

    
    #' #'  2. Extract data from roads shapefile & calculate distance to nearest road
    #' #'     for each location
    #' #'  ------------------------------------------------------------------------
    #' dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
    #' dist2road <- as.data.frame(dist2road)
    #' dist2road$ID <- locs$ID
    #' dist2road$obs <- c(1:nrow(locs))
    #' #'  Append to covariate data frame
    #' covs <- covs %>%
    #'   full_join(dist2road, by = c("obs", "ID"))
    
    
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
        Used = 0
        )
    
    return(telem_covs)
    
  }
  
  #'  Run list of species used & available location data through function in parallel
  #'  This will take AWHILE even in parallel
  #'  3rd Order Selection
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
  md_avail_dat <- combo_data(md_available_2nd, md_avail_covs)
  elk_avail_dat <- combo_data(elk_available_2nd, elk_avail_covs)
  wtd_avail_dat <- combo_data(wtd_available_2nd, wtd_avail_covs)
  coug_avail_dat <- combo_data(coug_available_2nd, coug_avail_covs)
  wolf_avail_dat <- combo_data(wolf_available_2nd, wolf_avail_covs)
  bob_avail_dat <- combo_data(bob_available_2nd, bob_avail_covs)
  coy_avail_dat <- combo_data(coy_available_2nd, coy_avail_covs)
  
  #'  Function to drop unneeded columns from list of used data sets and indicate 
  #'  whether this location was used = 1 or available = 0
  select_cols <- function(dat) {
    used_skinny <- dat %>%
      dplyr::select(x, y, AnimalID, Season, ID, Season.1, Year, Elev, Slope, RoadDen, 
                    HumanMod, PercForMix, PercXGrass, PercXShrub, obs, Area) %>%
      mutate(
        Used = 1)
    colnames(used_skinny) <-  c("x", "y", "ID", "Season", "ID.1", "Season.1", 
                                "Year", "Elev", "Slope", "RoadDen", "HumanMod", 
                                "PercForMix", "PercXGrass", "PercXShrub", "obs", 
                                "Area", "Used")
    return(used_skinny)
  }
  #'  Run the list of used locations through function
  used_dat <- lapply(used_dat, select_cols)
    
  #'  Save used and available data separately
  #'  Adjust between 3rd & 2nd order depending on version of available locs used above
  save(used_dat, file = paste0("./Outputs/RSF_pts/used_dat_", Sys.Date(), ".RData"))
  save(md_avail_dat, file = paste0("./Outputs/RSF_pts/md_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(elk_avail_dat, file = paste0("./Outputs/RSF_pts/elk_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(wtd_avail_dat, file = paste0("./Outputs/RSF_pts/wtd_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(coug_avail_dat, file = paste0("./Outputs/RSF_pts/coug_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(wolf_avail_dat, file = paste0("./Outputs/RSF_pts/wolf_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(bob_avail_dat, file = paste0("./Outputs/RSF_pts/bob_avail_2nd_dat_", Sys.Date(), ".RData"))
  save(coy_avail_dat, file = paste0("./Outputs/RSF_pts/coy_avail_2nd_dat_", Sys.Date(), ".RData"))
  
  
  #'  Merge used & available data per species
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
  #'  Adjust between 3rd & 2nd order depending on version of available locs used above
  save(md_dat_all, file = paste0("./Outputs/RSF_pts/md_dat_2nd_all_", Sys.Date(), ".RData"))
  save(elk_dat_all, file = paste0("./Outputs/RSF_pts/elk_dat_2nd_all_", Sys.Date(), ".RData"))
  save(wtd_dat_all, file = paste0("./Outputs/RSF_pts/wtd_dat_2nd_all_", Sys.Date(), ".RData"))
  save(coug_dat_all, file = paste0("./Outputs/RSF_pts/coug_dat_2nd_all_", Sys.Date(), ".RData"))
  save(wolf_dat_all, file = paste0("./Outputs/RSF_pts/wolf_dat_2nd_all_", Sys.Date(), ".RData"))
  save(bob_dat_all, file = paste0("./Outputs/RSF_pts/bob_dat_2nd_all_", Sys.Date(), ".RData"))
  save(coy_dat_all, file = paste0("./Outputs/RSF_pts/coy_dat_2nd_all_", Sys.Date(), ".RData"))

  
  
  # 2021-06-22 uses reprojected rasters
  
  
  
  