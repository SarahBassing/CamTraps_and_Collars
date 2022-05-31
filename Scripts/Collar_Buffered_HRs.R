  #'  ========================================
  #'  Buffered Home Ranges
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2022
  #'  ========================================
  #'  Generate individual home ranges for all collared individuals within the 
  #'  study. Buffer each home range by diameter of home range and mask out areas
  #'  unavailable to animals (waterbodies). These polygons will be used to select
  #'  random available locations for 2nd order RSF analyses. Using telemetry data
  #'  that have been cleaned (spurious locations, dispersals, migrations removed)   
  #'  & truncated around capture and mortality events but no further filtering done.
  #'  These data exclude seasonal migration movements for mule deer.
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
  
  #'  Load telemetry data
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispersal_allSeasons.RData") 
  
  #'  Define coordinate projection
  sa_proj <- CRS("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Load study area shapefiles
  OK_SA_sp <- readOGR(dsn = "./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA_sp <- readOGR(dsn = "./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  OK_SA_sp <- spTransform(OK_SA_sp, sa_proj)
  NE_SA_sp <- spTransform(NE_SA_sp, sa_proj)
    
  #'  Load water bodies shapefile (needed to mask unavailable habitat)
  waterbody <- sf::st_read("./Shapefiles/WA_DeptEcology_HydroWA", layer = "WPPP_waterbody") %>%
    st_transform(crs = sa_proj)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbody[waterbody$AreSqKm > 1,]

  
  #'  Filter data by study area
  #'  =========================
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
  
  #'  Prep lists to estimate annual home range per individual
  #'  =======================================================
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
  
  
  #'  Function to estimate home ranges & create buffered polygons
  #'  ===========================================================
  avail_pts <- function(locs, plotit = F) { #ex, h, if using KDEs
    
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
    
    #'  2. Create MCPs for each animal 
    #'  ------------------------------
    MCP100 <- mcp(locs_sp, percent = 100, unin = "m", unout = "km2")
    #' #'  Notes about estimating KDEs for each home range: 
    #' #'    1) extend the spatial extent by 1.5 - 2.5 when estimating UDs (throws 
    #' #'    an error about grid being too small to estimate home range otherwise) 
    #' #'    2) adjust the smoothing parameter (h) from the reference bandwidth 
    #' #'    ("href") to least squares cross validation ("LSCV") if estimates are
    #' #'    too large (oversmoothed) but LSCV is prone to convergence failure. 
    #' #'    Explore LSCV convergence issues with plotLSCV(UD95)
    #' UD <- kernelUD(locs_sp, kern = "bivnorm", extent = ex, h = h)
    #' UD95 <- getverticeshr(UD, 95)
    #' UD50 <- getverticeshr(UD, 50)
    #' UD75 <- getverticeshr(UD, 75)
    #' UD90 <- getverticeshr(UD, 90)

    #'  Calculate area of 100% MCP (default of “m” in and “ha” for output)
    MCP100_area <- mcp.area(locs_sp, percent = 100, unin = "m", unout = "km2", plotit = FALSE) 
    MCP100_area <- round(MCP100_area$a, 2)
    print(MCP100_area)
    # UD95_area <- kernel.area(UD, percent = 95, unin = "m", unout = "km2")
    # UD95_km2 <- round(UD95_area, 2)
    # print(UD95_km2)
    
    #'  Create buffer around home range based on approx. diameter of home range
    #'  (pretending home range is a perfect circle)
    #'  Multiply by 1000 so buffer is measured in meters, not kilometers
    buff <- (sqrt(MCP100_area/pi)*2)*1000 
    # buff <- (sqrt(UD95_km2/pi)*2)*1000

    #'  Buffer home range polygon by diameter of home range
    buffered_poly <- gBuffer(MCP100, width = buff)
    # buffered_poly <- gBuffer(UD95, width = buff)

    #'  Intersect and clip water body polygons from buffered polygons so large
    #'  bodies of water are not available to collared animals
    bigwater_sp <- as(st_geometry(bigwater), "Spatial")
    poly_clip <- rgeos::gDifference(buffered_poly, bigwater_sp)
    
    #'  Plot to make sure buffering and clipping worked
    if(plotit) {
      plot(poly_clip, border = "red", col = NA)
      plot(OK_SA_sp, add = T)
      plot(NE_SA_sp, add = T)
      plot(locs_sp, col = "blue", pch = 19, cex = 0.75, add = T)
      plot(MCP100, border = "darkgreen", col = NA, add = T)
      # plot(UD50, border = "green", col = NA, add = T)
      # plot(UD75, border = "green", col = NA, add = T)
      # plot(UD90, border = "green", col = NA, add = T)
      # plot(UD95, border = "darkgreen", col = NA, add = T)

    }
    
    #'  Append AnimalID to polygon
    poly_clipDF <- SpatialPolygonsDataFrame(poly_clip, AnimalID)
    
    return(poly_clipDF)
  }
  #'  Estimate annual home range for each individual per study area
  #'  Add ex = 2.5, h = "href" if creating KDEs for mule deer; ex = 1.5 for other spp
  md_OK_poly <- lapply(OK_split[[1]], avail_pts, F) 
  elk_NE_poly <- lapply(NE_split[[1]], avail_pts, F)
  wtd_NE_poly <- lapply(NE_split[[2]], avail_pts, F)
  coug_OK_poly <- lapply(OK_split[[2]], avail_pts, F)
  coug_NE_poly <- lapply(NE_split[[3]], avail_pts, F)
  wolf_OK_poly <- lapply(OK_split[[3]], avail_pts, F)
  wolf_NE_poly <- lapply(NE_split[[4]], avail_pts, F)
  bob_OK_poly <- lapply(OK_split[[4]], avail_pts, F)
  bob_NE_poly <- lapply(NE_split[[5]], avail_pts, F)
  coy_OK_poly <- lapply(OK_split[[5]], avail_pts, F)
  coy_NE_poly <- lapply(NE_split[[6]], avail_pts, F)
  
  #'  List polygons together
  HR_poly <- list(md_OK_poly, elk_NE_poly, wtd_NE_poly, coug_OK_poly, coug_NE_poly, 
                  wolf_OK_poly, wolf_NE_poly, bob_OK_poly, bob_NE_poly, coy_OK_poly, 
                  coy_NE_poly)
  
  #'  Save to use when sampling "available" resources
  save(HR_poly, file = "./Outputs/MCPs/MCP_HomeRange_Polygons_allSpp.RData")
  
  
  
  #'  Use this to define available extent for each animal in RSF analyses
  #'  Next: Collar_RSF_DataPrep.R
  
  
  
  #'  Create spatial polygons for mapping
  load("./Outputs/MCPs/MCP_HomeRange_Polygons_allSpp.RData")
  
  #'  Convert buffered HR polygons to sf objects
  HR_poly_sf <- function(poly) {
    buffHR_sf <- st_as_sf(poly)
    return(buffHR_sf)
  }
  buffHR_sf <- lapply(HR_poly, lapply, HR_poly_sf)
  #'  Split up by species, merge into single df and dissolve inner boundaries of 
  #'  overlappling buffered MCPS
  md_buffHR <- buffHR_sf[[1]] %>% dplyr::bind_rows(.) %>% st_union(.)
  elk_buffHR <- buffHR_sf[[2]] %>% dplyr::bind_rows(.) %>% st_union(.)
  wtd_buffHR <- buffHR_sf[[3]] %>% dplyr::bind_rows(.) %>% st_union(.)
  coug_OK_buffHR <- buffHR_sf[[4]] %>% dplyr::bind_rows(.) %>% st_union(.)
  coug_NE_buffHR <- buffHR_sf[[5]] %>% dplyr::bind_rows(.) %>% st_union(.)
  wolf_OK_buffHR <- buffHR_sf[[6]] %>% dplyr::bind_rows(.) %>% st_union(.)
  wolf_NE_buffHR <- buffHR_sf[[7]] %>% dplyr::bind_rows(.) %>% st_union(.)
  bob_OK_buffHR <- buffHR_sf[[8]] %>% dplyr::bind_rows(.) %>% st_union(.)
  bob_NE_buffHR <- buffHR_sf[[9]] %>% dplyr::bind_rows(.) %>% st_union(.)
  coy_OK_buffHR <- buffHR_sf[[10]] %>% dplyr::bind_rows(.) %>% st_union(.)
  coy_NE_buffHR <- buffHR_sf[[11]] %>% dplyr::bind_rows(.) %>% st_union(.)
  
  #'  List polygons together
  buffHR_poly <- list(md_buffHR, elk_buffHR, wtd_buffHR, coug_OK_buffHR, coug_NE_buffHR, 
                  wolf_OK_buffHR, wolf_NE_buffHR, bob_OK_buffHR, bob_NE_buffHR, coy_OK_buffHR, 
                  coy_NE_buffHR)
  #  Save for plotting
  save(buffHR_poly, file = "./Outputs/MCPs/MCP_BufferedHomeRange_Polygons_allSpp.RData")
  
  
  #'  Generate individual home ranges (MCPs) without buffering
  HR_poly <- function(locs, plotit = F) { 
    
    #'  Make each animal's locations spatial
    locs_sp <- SpatialPointsDataFrame(locs[,c("x", "y")], data = locs[2], proj4string = sa_proj) 
    
    #'  Create MCPs for each animal 
    MCP100 <- mcp(locs_sp, percent = 100, unin = "m", unout = "km2")
    
    #'  Convert to sf object instead of sp and add AnimalID back to df
    MCP100_sf <- st_as_sf(MCP100)
    
    #'  Plot to make sure buffering and clipping worked
    if(plotit) {
      plot(locs_sp, col = "blue", pch = 19, cex = 0.75)
      plot(OK_SA_sp, add = T)
      plot(NE_SA_sp, add = T)
      plot(MCP100, border = "darkgreen", col = NA, add = T)
    }
    
    return(MCP100_sf)
  }
  #'  Estimate annual home range for each individual per study area, save as
  #'  single sf object, and dissolve inner boundaries of overlapping MCPs
  md_OK_HRpoly <- lapply(OK_split[[1]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  elk_NE_HRpoly <- lapply(NE_split[[1]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  wtd_NE_HRpoly <- lapply(NE_split[[2]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  coug_OK_HRpoly <- lapply(OK_split[[2]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  coug_NE_HRpoly <- lapply(NE_split[[3]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  wolf_OK_HRpoly <- lapply(OK_split[[3]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  wolf_NE_HRpoly <- lapply(NE_split[[4]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  bob_OK_HRpoly <- lapply(OK_split[[4]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  bob_NE_HRpoly <- lapply(NE_split[[5]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  coy_OK_HRpoly <- lapply(OK_split[[5]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
  coy_NE_HRpoly <- lapply(NE_split[[6]], HR_poly, F) %>% dplyr::bind_rows(.) %>% st_union(.)
 
  #'  List polygons together
  indHR_poly <- list(md_OK_HRpoly, elk_NE_HRpoly, wtd_NE_HRpoly, coug_OK_HRpoly, coug_NE_HRpoly, 
                  wolf_OK_HRpoly, wolf_NE_HRpoly, bob_OK_HRpoly, bob_NE_HRpoly, coy_OK_HRpoly, 
                  coy_NE_HRpoly)
  
  #'  Save for making figures
  save(indHR_poly, file = "./Outputs/MCPs/MCP_individual_HomeRange_Polygons_allSpp.RData")
  
  
  
  
  
  