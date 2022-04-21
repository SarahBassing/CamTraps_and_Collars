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
  
  #' #'  Clean workspace
  #' rm(list = ls())

  #'  Load libraries
  library(adehabitatHR)
  library(sf)
  library(sp)
  library(rgdal)
  library(rgeos)
  library(tidyverse)
  
  #'  Define coordinate projection
  wgs84 <- st_crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  sa_proj <- st_crs("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Load study area shapefiles
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA")
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA")
  OK_wgs84 <- st_transform(OK_SA, sa_proj)
  NE_wgs84 <- st_transform(NE_SA, sa_proj)
  
  #'  Load water bodies shapefile (needed to mask unavailable habitat)
  waterbody <- sf::st_read("./Shapefiles/WA_DeptEcology_HydroWA", layer = "WPPP_waterbody") %>%
    st_transform(crs = sa_proj)
  #'  Identify large bodies of water (anything larger than 1 sq-km in size)
  bigwater <- waterbody[waterbody$AreSqKm > 1,]
  
  #'  Load cleaned data    ### NEED TO BREAK APART BY STUDY AREA!!!
  load("./Data/Collar_AllSpecies_AllLocations_Clean.RData")
  md_skinny <- clean_data[[1]] %>%
    dplyr::select(Species, Latitude, Longitude) 
  elk_skinny <- clean_data[[2]] %>%
    dplyr::select(Species, Latitude, Longitude) 
  wtd_skinny <- clean_data[[3]] %>%
    dplyr::select(Species, Latitude, Longitude) 
  coug_skinny <- clean_data[[4]] %>%
    mutate(Species = "Cougar") %>%
    dplyr::select(Species, StudyArea, Latitude, Longitude) 
  wolf_skinny <- clean_data[[5]] %>%
    mutate(Species = "WOlf") %>%
    dplyr::select(Species, StudyArea, Latitude, Longitude) 
  meso_skinny <- clean_data[[6]] %>%
    dplyr::select(Species, StudyArea, Latitude, Longitude) 
  bob_skinny <- meso_skinny[meso_skinny$Species == "Bobcat",]
  coy_skinny <- meso_skinny[meso_skinny$Species == "Coyote",]

  #'  Separate predator locations by study area so MCP are study area specific
  coug_NE <- coug_skinny[coug_skinny$StudyArea == "NE",] %>%
    dplyr::select(-StudyArea)
  coug_OK <- coug_skinny[coug_skinny$StudyArea == "OK",] %>%
    dplyr::select(-StudyArea)
  wolf_NE <- wolf_skinny[wolf_skinny$StudyArea == "NE",] %>%
    dplyr::select(-StudyArea)
  wolf_OK <- wolf_skinny[wolf_skinny$StudyArea == "OK",] %>%
    dplyr::select(-StudyArea)
  bob_NE <- bob_skinny[bob_skinny$StudyArea == "NE",] %>%
    dplyr::select(-StudyArea)
  bob_OK <- bob_skinny[bob_skinny$StudyArea == "OK",] %>%
    dplyr::select(-StudyArea)
  coy_NE <- coy_skinny[coy_skinny$StudyArea == "NE",] %>%
    dplyr::select(-StudyArea)
  coy_OK <- coy_skinny[coy_skinny$StudyArea == "OK",] %>%
    dplyr::select(-StudyArea)

  
  #'  Make telemetry locations spatial and transform projection to UTMs
  #'  I don't know why I can't get spTrnasform to work for me but whatever
  sp_object <- function(locs) {
    sf_ob <- st_as_sf(locs, coords = c("Longitude", "Latitude"), crs = wgs84) %>%
      st_transform(crs = sa_proj)
    sp_ob <- as(sf_ob, "Spatial")
    return(sp_ob)
  }
  md_sp <- sp_object(md_skinny)
  elk_sp <- sp_object(elk_skinny)
  wtd_sp <- sp_object(wtd_skinny)
  coug_NE_sp <- sp_object(coug_NE)
  coug_OK_sp <- sp_object(coug_OK)
  wolf_NE_sp <- sp_object(wolf_NE)
  wolf_OK_sp <- sp_object(wolf_OK)
  bob_NE_sp <- sp_object(bob_NE)
  bob_OK_sp <- sp_object(bob_OK)
  coy_NE_sp <- sp_object(coy_NE)
  coy_OK_sp <- sp_object(coy_OK)

  #'  Create massive MCP that includes all locations from a species
  #'  This intentionally ignores individual home ranges, migrations, and year
  md_mcp <- mcp(md_sp, percent = 100) 
  elk_mcp <- mcp(elk_sp, percent = 100)
  wtd_mcp <- mcp(wtd_sp, percent = 100)
  coug_NE_mcp <- mcp(coug_NE_sp, percent = 100)
  coug_OK_mcp <- mcp(coug_OK_sp, percent = 100)
  wolf_NE_mcp <- mcp(wolf_NE_sp, percent = 100)
  wolf_OK_mcp <- mcp(wolf_OK_sp, percent = 100)
  bob_NE_mcp <- mcp(bob_NE_sp, percent = 100)
  bob_OK_mcp <- mcp(bob_OK_sp, percent = 100)
  coy_NE_mcp <- mcp(coy_NE_sp, percent = 100)
  coy_OK_mcp <- mcp(coy_OK_sp, percent = 100)
  
  #'  Expand MCPs by 3km to allow available points to be drawn from just beyond-
  #'  individuals could have moved beyond the used points within 4hr fix interval
  #'  Width is in units of the projection (meters)
  md_poly <- gBuffer(md_mcp, width = 3000)
  elk_poly <- gBuffer(elk_mcp, width = 3000)
  wtd_poly <- gBuffer(wtd_mcp, width = 3000)
  coug_NE_poly <- gBuffer(coug_NE_mcp, width = 3000)
  coug_OK_poly <- gBuffer(coug_OK_mcp, width = 3000)
  wolf_NE_poly <- gBuffer(wolf_NE_mcp, width = 3000)
  wolf_OK_poly <- gBuffer(wolf_OK_mcp, width = 3000)
  bob_NE_poly <- gBuffer(bob_NE_mcp, width = 3000)
  bob_OK_poly <- gBuffer(bob_OK_mcp, width = 3000)
  coy_NE_poly <- gBuffer(coy_NE_mcp, width = 3000)
  coy_OK_poly <- gBuffer(coy_OK_mcp, width = 3000)
  
  #'  Intersect and clip water body polygons from MCP polygons so large bodies
  #'  of water are not available to collared animals
  bigwater_sp <- as(st_geometry(bigwater), "Spatial")
  md_poly_clip <- rgeos::gDifference(md_poly, bigwater_sp)
  elk_poly_clip <- rgeos::gDifference(elk_poly, bigwater_sp)
  wtd_poly_clip <- rgeos::gDifference(wtd_poly, bigwater_sp)
  coug_NE_poly_clip <- rgeos::gDifference(coug_NE_poly, bigwater_sp)
  coug_OK_poly_clip <- rgeos::gDifference(coug_OK_poly, bigwater_sp)
  wolf_NE_poly_clip <- rgeos::gDifference(wolf_NE_poly, bigwater_sp)
  wolf_OK_poly_clip <- rgeos::gDifference(wolf_OK_poly, bigwater_sp)
  bob_NE_poly_clip <- rgeos::gDifference(bob_NE_poly, bigwater_sp)
  bob_OK_poly_clip <- rgeos::gDifference(bob_OK_poly, bigwater_sp)
  coy_NE_poly_clip <- rgeos::gDifference(coy_NE_poly, bigwater_sp)
  coy_OK_poly_clip <- rgeos::gDifference(coy_OK_poly, bigwater_sp)

  
  #' #'  Plot to make sure these make sense
  #' plot(md_mcp)
  #' plot(md_sp, add = TRUE, pch = 16, col = "blue")
  #' plot(md_poly, add = TRUE)
  #' plot(md_poly_clip, add = TRUE)
  #' 
  #' plot(elk_mcp)
  #' plot(elk_sp, add = TRUE, pch = 16, col = "blue")
  #' plot(elk_poly, add = TRUE)
  #' plot(elk_poly_clip, add = TRUE)
  #' 
  #' plot(wtd_mcp)
  #' plot(wtd_sp, add = TRUE, pch = 16, col = "blue")
  #' plot(wtd_poly, add = TRUE)
  #' plot(wtd_poly_clip, add = TRUE)
  #' 
  #' plot(coug_OK_mcp)
  #' plot(coug_OK_sp, add = TRUE, pch = 16, col = "red")
  #' plot(coug_OK_poly, add = TRUE)
  #' plot(coug_OK_poly_clip, add = TRUE)
  #' 
  #' plot(coug_NE_mcp)
  #' plot(coug_NE_sp, add = TRUE, pch = 16, col = "red")
  #' plot(coug_NE_poly, add = TRUE)
  #' plot(coug_NE_poly_clip, add = TRUE)
  #' 
  #' plot(wolf_OK_mcp)
  #' plot(wolf_OK_sp, add = TRUE, pch = 16, col = "red")
  #' plot(wolf_OK_poly, add = TRUE)
  #' plot(wolf_OK_poly_clip, add = TRUE)
  #' 
  #' plot(wolf_NE_mcp)
  #' plot(wolf_NE_sp, add = TRUE, pch = 16, col = "red")
  #' plot(wolf_NE_poly, add = TRUE)
  #' plot(wolf_NE_poly_clip, add = TRUE)
  #' 
  #' plot(bob_OK_mcp)
  #' plot(bob_OK_sp, add = TRUE, pch = 16, col = "red")
  #' plot(bob_OK_poly, add = TRUE)
  #' plot(bob_OK_poly_clip, add = TRUE)
  #' 
  #' plot(bob_NE_mcp)
  #' plot(bob_NE_sp, add = TRUE, pch = 16, col = "red")
  #' plot(bob_NE_poly, add = TRUE)
  #' plot(bob_NE_poly_clip, add = TRUE)
  #' 
  #' plot(coy_OK_mcp)
  #' plot(coy_OK_sp, add = TRUE, pch = 16, col = "red")
  #' plot(coy_OK_poly, add = TRUE)
  #' plot(coy_OK_poly_clip, add = TRUE)
  #' 
  #' plot(coy_NE_mcp)
  #' plot(coy_NE_sp, add = TRUE, pch = 16, col = "red")
  #' plot(coy_NE_poly, add = TRUE)
  #' plot(coy_NE_poly_clip, add = TRUE)

  #' #'  Coerce buffered MCP SpatialPolygons into SpatialPolygonDataFrames and save
  #' md_poly <- as(md_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(md_poly, dsn = "./Outputs/MCPs", layer = "md_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' elk_poly <- as(elk_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(elk_poly, "./Outputs/MCPs", "elk_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' wtd_poly <- as(wtd_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(wtd_poly, "./Outputs/MCPs", "wtd_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' coug_NE_poly <- as(coug_NE_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(coug_NE_poly, "./Outputs/MCPs", "coug_NE_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' coug_OK_poly <- as(coug_OK_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(coug_OK_poly, "./Outputs/MCPs", "coug_OK_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' wolf_NE_poly <- as(wolf_NE_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(wolf_NE_poly, "./Outputs/MCPs", "wolf_NE_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' wolf_OK_poly <- as(wolf_OK_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(wolf_OK_poly, "./Outputs/MCPs", "wolf_OK_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' bob_NE_poly <- as(bob_NE_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(bob_NE_poly, "./Outputs/MCPs", "bob_NE_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' bob_OK_poly <- as(bob_OK_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(bob_OK_poly, "./Outputs/MCPs", "bob_OK_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' coy_NE_poly <- as(coy_NE_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(coy_NE_poly, "./Outputs/MCPs", "coy_NE_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  #' coy_OK_poly <- as(coy_OK_poly, "SpatialPolygonsDataFrame")
  #' writeOGR(coy_OK_poly, "./Outputs/MCPs", "coy_OK_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE)
  
  