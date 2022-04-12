  #'  ========================================
  #'  Minimum Convex Polygons
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  July 3, 2021
  #'  ========================================
  #'  Generate species-specific MCPs for all collared individuals within the 
  #'  study. Ignoring individual home ranges and generating one per species based
  #'  on ALL locations from that species. These MCPs will be used to select
  #'  random available locations for 2nd order RSF analyses. Using telemetry data
  #'  that have been cleaned (spurious locations, dispersals, etc. removed) and  
  #'  truncated around capture and mortality events but no further filtering done.
  #'  These data do include any seasonal migration movement.
  #'  
  #'  Ungulate data cleaned with Collar_DataCleaning.R script, meso-predator data
  #'  cleaned with Collar_DataCleaning_Mesopredator.R script, and large predator
  #'  data cleaned by L.Satterfield.
  #'  ========================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(adehabitatHR)
  library(sf)
  library(sp)
  library(rgdal)
  library(rgeos)
  library(tidyverse)
  
  #'  Define coordinate projection
  # wgs84 <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  sa_proj <- CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Load study area shapefiles
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  OK_wgs84 <- st_transform(OK_SA, sa_proj)
  NE_wgs84 <- st_transform(NE_SA, sa_proj)
  OK_SA_sp <- readOGR(dsn = "./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA_sp <- readOGR(dsn = "./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  OK_SA_sp <- spTransform(OK_SA_sp, sa_proj)
  NE_SA_sp <- spTransform(NE_SA_sp, sa_proj)
    
  #'  Load water bodies shapefile (needed to mask unavailable habitat)
  waterbody <- sf::st_read("./Shapefiles/WA_DeptEcology_HydroWA", layer = "WPPP_waterbody") %>%
    st_transform(crs = sa_proj)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbody[waterbody$AreSqKm > 1,]
  
  #'  Load telemetry data
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispersal_allSeasons.RData") #noDispMig

  #'  Functions to filter species-specific data by study area
  #'  NE study area collars
  NE_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "NE")
    return(tracks)
  }
  #'  Drop mule deer summer and winter lists since no mulies in NE study area
  noMD <- spp_all_tracks_all_seasons[-1]
  #'  Run list through function
  NE_tracks <- lapply(noMD, NE_collars)
  #'  OK study area collars
  OK_collars <- function(locs) {
    tracks <- filter(locs, StudyArea == "OK")
    return(tracks)
  }
  #'  Drop elk and white-tailed deer summer and winter lists since no elk or wtd in OK
  noELKWTD <- spp_all_tracks_all_seasons[-c(2:3)]
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
  animal_ID <- lapply(spp_all_tracks_all_seasons, unq_id) #spp_all_tracks
  NE_ID <- lapply(NE_tracks, unq_id)
  OK_ID <- lapply(OK_tracks, unq_id)
  
  #'  Function to split data into list of dataframes by unique animal ID and year
  #'  (Use AnimalID if ignoring year aspect of data)
  split_animal <- function(locs) {
    ind_animal <- group_split(locs, FullID) #FullID so it's by year too (important for landcover data extraction)
    return(ind_animal)
  }
  animal_split <- lapply(spp_all_tracks_all_seasons, split_animal)
  NE_split <- lapply(NE_tracks, split_animal)
  OK_split <- lapply(OK_tracks, split_animal)
  
  
  avail_pts <- function(locs, ex, plotit = F) { 
    #'  1. Make each animal's locations spatial
    #'  ---------------------------------------------------------
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj) 
    #'  Plot to make sure step 1 worked
    if(plotit) {
      plot(locs_sp, col = "blue", pch = 19, cex = 0.75, main = locs$FullID) 
      # plot(OK_SA_sp, add = T)
    }
    
    #'  Hold out AnimalID info to add back to individual polygons below
    AnimalID <- unique(locs[,"FullID"]) # FullID b/c split data based on this
    AnimalID <- as.data.frame(AnimalID)
    
    #'  2. Create UDs for each animal following Bivariate normal utilization distributions
    #'  ----------------------------------------------------------
    #'  Estimate KDEs for each home range, extend the spatial extent by 1.5
    #'  when estimating UDs (throws an error about grid being too small to
    #'  estimate home range otherwise)
    # MCP95 <- mcp(locs_sp, percent = 95) 
    UD <- kernelUD(locs_sp, kern = "bivnorm", extent = ex)
    UD95 <- getverticeshr(UD, 95)
    UD50 <- getverticeshr(UD, 50)
    UD75 <- getverticeshr(UD, 75)
    UD90 <- getverticeshr(UD, 90)
    
    #'  Calculate area of 95% UD
    UD95_area <- kernel.area(UD, percent = 95)
    #'  Convert area to square-kilometers
    UD95_km2 <- round(UD95_area/100, 2)
    print(UD95_km2)
    #'  Create buffer around home range based on approx. diameter of home range
    #'  (pretending home range is a perfect circle)
    #'  Multiply by 1000 so buffer is measured in meters, not kilometers
    buff <- (sqrt(UD95_km2/pi)*2)*1000
    
    #'  Buffer home range by the area of the 95% KDE
    buffered_poly <- gBuffer(UD95, width = buff)
    
    #'  Intersect and clip water body polygons from buffered polygons so large 
    #'  bodies of water are not available to collared animals
    bigwater_sp <- as(st_geometry(bigwater), "Spatial")
    poly_clip <- rgeos::gDifference(buffered_poly, bigwater_sp)
    
    #'  Plot to make sure buffering and clipping worked
    if(plotit) {
      plot(poly_clip, border = "red", col = NA)
      plot(locs_sp, col = "blue", pch = 19, cex = 0.75, add = T)
      plot(UD50, border = "green", col = NA, add = T)
      plot(UD75, border = "green", col = NA, add = T)
      plot(UD90, border = "green", col = NA, add = T)
      plot(UD95, border = "darkgreen", col = NA, add = T)
      
    }
     
    #'  Append AnimalID to polygon
    poly_clipDF <- SpatialPolygonsDataFrame(poly_clip, AnimalID)
    
    
    #' #'  3. Randomly select points within each home range
    #' #'  ------------------------------------------------
    #' #'  Sampling n available points to every 1 used point
    #' #'  Identify number of used points per individual
    #' nused <- nrow(locs)
    #' #'  Multiply by desired number of available points set by navail
    #' navailable <- nused*navail
    #' #'  Set seed for reproducibility
    #' # set.seed(2021)
    #' rndpts <- spsample(poly_clipDF, navailable, type = "random")
    #' #'  Turn them into spatial points
    #' rndpts_sp <- SpatialPoints(rndpts, proj4string = sa_proj)
    #' #' Plot to make sure step 3 worked
    #' if(plotit) {
    #'   plot(rndpts_sp, col = "red", pch = 19)
    #'   plot(poly_clip, border = "red", col = NA, add = T)
    #'   plot(UD95, border = "darkgreen", col = NA, add = T)
    #' }
    #' 
    #' #'  4. Make list of locations non-spatial
    #' #'  -------------------------------------
    #' rndpts_df <- as.data.frame(rndpts_sp) 
    #' ID <- unique(droplevels(locs$AnimalID))
    #' Season <- unique(locs$Season)
    #' rndpts_df$ID <- ID
    #' rndpts_df$Season <- Season
    #' 
    #' return(rndpts_df)
    return(poly_clipDF)
    
  }
  md_poly <- lapply(animal_split[[1]], avail_pts, ex = 2.5, T) #ex = 2.5 for noDispMig
  elk_poly <- lapply(animal_split[[2]], avail_pts, ex = 1.5, T)
  wtd_poly <- lapply(animal_split[[3]], avail_pts, ex = 1.5, T)
  coug_poly <- lapply(animal_split[[4]], avail_pts, ex = 1.5, T)
  wolf_poly <- lapply(animal_split[[5]], avail_pts, ex = 1.5, T)
  bob_poly <- lapply(animal_split[[6]], avail_pts, ex = 1.5, T)
  coy_poly <- lapply(animal_split[[7]], avail_pts, ex = 1.5, T)
  
  #'  Use this to define available extent for each animal in RSF analyses