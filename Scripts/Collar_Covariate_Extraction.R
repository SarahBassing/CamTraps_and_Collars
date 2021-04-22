  #'  ============================================
  #'  Covariate extraction for telemetry locations
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  April 2021
  #'  ============================================
  #'  Script to extract covariate data at each GPS collar location and interpolated
  #'  locations based on the crawlWrap function. Covariates gathered from 
  #'  multiple sources and described in Covariate_Extract.R script from
  #'  WPPP_CameraTrapping repository. Relevant covariates include:
  #'    -Elevation (30m res)
  #'    -Slope (30m res)
  #'    -Percent Mixed Forest (within 250m of point)
  #'    -Percent Xeric Grass (within 250m of point)
  #'    -Percent Xeric Shrub (within 205m of point)
  #'    -Distance to nearest road
  #'    -Human Modified Landscape (1km res)
  #'    -Study Area
  #'  ============================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(stars)
  library(rgeos)
  library(raster)
  library(tidyverse)
  
  #'  Read in covariate data extracted from other sources
  dist2road <- read.csv("./Output/dist2road18-20.csv") %>%   ### NEED TO GENERATE THIS ON LAB COMPUTER
    mutate(
      km2road = dist2road/1000
    ) %>%
    dplyr::select(-X)
  
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
  roads <- st_read("./Shapefiles/Cascadia_layers/roadsForTaylor", layer = "roadsForTaylor")
  #'  Human density and human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  projection(wppp_bound)
  projection(dem)
  projection(slope)
  projection(landcov18)
  projection(interp_landcov18)
  projection(HM)
  projection(roads)
  
  res(dem)
  res(landcov18)
  res(interp_landcov18)
  res(HM)

  #'  Load animal location data for each species
  load("./Outputs/Telemetry_crwOut/crwOut_MD.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_ELK.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_WTD.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_COUG.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_WOLF.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_BOB.RData")
  load("./Outputs/Telemetry_crwOut/crwOut_COY.RData")
  
  #'  Pull out locations and interpolated data
  md_move <- crwOut_MD[[2]]
  elk_move <- crwOut_ELK[[2]]
  wtd_move <- crwOut_WTD[[2]]
  coug_move <- crwOut_COUG[[2]]
  wolf_move <- crwOut_WOLF[[2]]
  bob_move <- crwOut_BOB[[2]]
  coy_move <- crwOut_COY[[2]]
  
  #'  Make location data spatial----- UPDATE EXACT COLUMNS
  md_locs <- st_as_sf(md_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  elk_locs <- st_as_sf(elk_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  wtd_locs <- st_as_sf(wtd_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  coug_locs <- st_as_sf(coug_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  wolf_locs <- st_as_sf(wolf_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  bob_locs <- st_as_sf(bob_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  coy_locs <- st_as_sf(coy_move, coords = c("mu.x", "mu.y"), crs = sa_proj)
  

  
  #'  Calculate distance from each animal location to nearest road and take the 
  #'  minimum value -- default units of measurement are in METERS
  #'  Make sure features are in either Equidistant or State-specific projection
  #'  to preserve distances between features (DO NOT use WGS84)
  #'  Albers Equal Area projection good for measuring area
  
  #'  Reproject road shapefile to match animal location projection
  road_reproj <- st_transform(roads, crs = st_crs(sa_proj))
  projection(road_reproj)
  
  #'  Calculate distance to nearest road feature
  #'  Test with one location
  dist2road <- min(st_distance(road_reproj, wolf_locs[1,]))
  dist2road <- as.data.frame(dist2road)
  dist2road$AnimalID <- wolf_locs$AnimalID[1]
  
  #'  Function to iterate over all points per species (this will take awhile!)
  RoadDist <- function(locs) {
    dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
    dist2road <- as.data.frame(dist2road)
    dist2road$AnimalID <- locs$AnimalID
    return(dist2road)
  }
  #'  Run each species through function to calculate distance to nearest road
  #'  THIS WILL TAKE FOREVER!!!
  locs <- list(md_locs, elk_locs, wtd_locs, coug_locs, wolf_locs, bob_locs, coy_locs)
  spp_nearestRd <- lapply(locs, RoadDist)
  tst <- RoadDist(wolf_locs)
  
  #' #'  Iterate over all points (this will take awhile!)
  #' dist2road <- sapply(1:nrow(reproj_cams), function(x) min(st_distance(rd_reproj, reproj_cams[x, ])))
  #' dist2road <- as.data.frame(dist2road)
  #' dist2road$CameraLocation <- camera_stations$CameraLocation
  #' 
  #' #'  Save save save
  #' write.csv(dist2road, "./Output/dist2road18-20.csv")
  
  
  #'  Extract covariate data from rasters
  #'  Reproject to match rasters in wgs84
  reproj_md <- st_transform(md_locs, crs = st_crs(wgs84))
  reproj_elk <- st_transform(elk_locs, crs = st_crs(wgs84))
  reproj_wtd <- st_transform(wtd_locs, crs = st_crs(wgs84))
  reproj_coug <- st_transform(coug_locs, crs = st_crs(wgs84))
  reproj_wolf <- st_transform(wolf_locs, crs = st_crs(wgs84))
  reproj_bob <- st_transform(bob_locs, crs = st_crs(wgs84))
  reproj_coy <- st_transform(coy_locs, crs = st_crs(wgs84))
  
  #'  Function to extract covariates from raster data for each species
  rast_extract <- function(locs) {
    #'  Extract covariates that share same projection
    elevation <- raster::extract(dem, locs, df = TRUE)
    slope <- raster::extract(slope, locs, df = TRUE)
    modified <- raster::extract(HM, locs, df = TRUE)
    # pixvals18 <- raster::extract(interp_landcov18, reproj_wolf, factors = TRUE, buffer = 250, df = TRUE)
    # pixvals_df18 <- as.data.frame(pixvals18)
    # pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
    # pixvals_df19 <- as.data.frame(pixvals19)
    covs <- full_join(elevation, slope, by = "ID") %>%
      full_join(modified, by = "ID") %>%
      # full_join(pixvals18, by = "ID") %>%
      # full_join(pixvals19, by = "ID")
      transmute(
        ID = ID,
        Elev = WPPP_DEM_30m,
        Slope = round(WPPP_slope_aspect, digits = 2),
        HumanMod = WPPP_gHM
        # Landcover_2018 = interp_landcov18,
        # Landcover_2019 = interp_landcov19
      )
    return(covs)
  }
  #'  Run each species through raster extract function
  locs <- list(reproj_md, reproj_elk, reproj_wtd, reproj_coug, reproj_wolf, reproj_bob, reproj_coy)
  spp_extract <- lapply(locs, rast_extract)
  #' #'  Unlist extracted data into separate data frames for each species
  #' names(spp_extract) <- c("md_extract", "elk_extract", "wtd_extract", "coug_extract", "wolf_extract", "bob_extract", "coy_extract")
  #' invisible(lapply(names(spp_extract),function(x) assign(x,spp_extract[[x]],.GlobalEnv)))
  
  
  
  
  
  #'  Function to extract landcover value from each pixel within 250m radius of 
  #'  animal locations using interpolated landcover rasters derived from 
  #'  Cascadia landcover
  perc_landcov <- function(locs) {
    pixvals18 <- raster::extract(interp_landcov18, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df18 <- as.data.frame(pixvals18)
    pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df19 <- as.data.frame(pixvals19)
    #'  Merge together
    ID <- as.data.frame(as.numeric(seq(1:nrow(locs))))
    Animal <- cbind(ID, AnimalID)
    colnames(Animal) <- c("ID", "AnimalID")
    landcover_250m <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) %>%
      full_join(Animal, by = "ID")
    colnames(landcover_250m) <- c("ID", "landcover_2018", "landcover_2019", "AnimalID")
    
    #'  Count the number of cells in each landcover category by CameraLocation
    tbl_landcover18 <- landcover_250m %>%
      group_by(AnimalID) %>%
      count(landcover_2018) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2018, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        sumPixels = rowSums(.[2:13])
      )
    #'  Drop landcover data from 2019
    tbl_landcover18 <- cbind(Year, tbl_landcover18) %>%
      filter(Year == "Year1")
    #'  Reogranize so it's easier to keep track of each category
    tbl_landcover18 <- tbl_landcover18[, order(colnames(tbl_landcover18), decreasing = TRUE)] %>%
      relocate(sumPixels, .after = last_col()) 
    colnames(tbl_landcover18) <- c("Year", "CameraLocation", "Residential",  
                                   "Commercial", "Agriculture", "Forest",  
                                   "XericShrub", "MesicShrub", "XericGrass", 
                                   "MesicGrass", "WoodyWetland", "EmergentWetland", 
                                   "Barren", "Water", "sumPixels")
    tbl_landcover19 <- landcover_250m %>%
      group_by(CameraLocation) %>%
      count(landcover_2019) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2019, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        sumPixels = rowSums(.[2:13])
      )
    #'  Drop landcover data from 2018
    tbl_landcover19 <- cbind(Year, tbl_landcover19) %>%
      filter(Year == "Year2")
    #'  Reogranize so it's easier to keep track of each category
    tbl_landcover19 <- tbl_landcover19[, order(colnames(tbl_landcover19), decreasing = TRUE)] %>%
      relocate(sumPixels, .after = last_col())
    colnames(tbl_landcover19) <- c("Year", "CameraLocation", "Residential",  
                                   "Commercial", "Agriculture", "Forest",  
                                   "XericShrub", "MesicShrub", "XericGrass", 
                                   "MesicGrass", "WoodyWetland", "EmergentWetland", 
                                   "Barren", "Water", "sumPixels")
    #'  Merge annual landcover values together
    tbl_landcover <- rbind(tbl_landcover18, tbl_landcover19) %>%
      #'  Consolidate categories
      #'  Keeping Water to include in percent calculation but will not use for analyses
      mutate(
        Forest =  Forest + WoodyWetland + EmergentWetland,
        MesicGrass = MesicGrass + Barren,
        Developed = Residential + Commercial + Agriculture,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      dplyr::select(-c(WoodyWetland, EmergentWetland, Barren, Residential, Commercial, Agriculture)) %>%
      relocate(sumPixels, .after = last_col()) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        PercXericShrub = round(XericShrub/sumPixels, 2),
        PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2),      # Cannot be used in conjunction with other Mesic landcover types
        PercWater = round(Water/sumPixels, 2),            # Don't use for analyses (I have better data for H2o)
        PercDeveloped = round(Developed/sumPixels, 2)
      )
    return(tbl_landcover)
  }
  
  
  
  
  
  
  # cams_reproj <- st_transform(cams, crs(crs(interp_landcov18)))
  pixvals18 <- raster::extract(interp_landcov18, reproj_wolf, factors = TRUE, buffer = 250, df = TRUE)
  pixvals_df18 <- as.data.frame(pixvals18)
  pixvals19 <- raster::extract(interp_landcov19, reproj_wolf, factors = TRUE, buffer = 250, df = TRUE)
  pixvals_df19 <- as.data.frame(pixvals19)
  #'  Merge together
  ID <- as.data.frame(as.numeric(seq(1:nrow(reproj_wolf))))
  Animal <- cbind(ID, Animal)
  colnames(Animal) <- c("ID", "AnimalID")
  landcover_250m <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) %>%
    full_join(Animal, by = "ID")
  colnames(landcover_250m) <- c("ID", "landcover_2018", "landcover_2019", "AnimalID")
  
  #'  Count the number of cells in each landcover category by CameraLocation
  tbl_landcover18 <- landcover_250m %>%
    group_by(CameraLocation) %>%
    count(landcover_2018) %>%
    ungroup() %>%
    pivot_wider(names_from = landcover_2018, values_from = n) %>%
    replace(is.na(.), 0) %>% 
    mutate(
      sumPixels = rowSums(.[2:13])
    )
  #'  Drop landcover data from 2019
  tbl_landcover18 <- cbind(Year, tbl_landcover18) %>%
    filter(Year == "Year1")
  #'  Reogranize so it's easier to keep track of each category
  tbl_landcover18 <- tbl_landcover18[, order(colnames(tbl_landcover18), decreasing = TRUE)] %>%
    relocate(sumPixels, .after = last_col()) 
  colnames(tbl_landcover18) <- c("Year", "CameraLocation", "Residential",  
                                 "Commercial", "Agriculture", "Forest",  
                                 "XericShrub", "MesicShrub", "XericGrass", 
                                 "MesicGrass", "WoodyWetland", "EmergentWetland", 
                                 "Barren", "Water", "sumPixels")
  tbl_landcover19 <- landcover_250m %>%
    group_by(CameraLocation) %>%
    count(landcover_2019) %>%
    ungroup() %>%
    pivot_wider(names_from = landcover_2019, values_from = n) %>%
    replace(is.na(.), 0) %>% 
    mutate(
      sumPixels = rowSums(.[2:13])
    )
  #'  Drop landcover data from 2018
  tbl_landcover19 <- cbind(Year, tbl_landcover19) %>%
    filter(Year == "Year2")
  #'  Reogranize so it's easier to keep track of each category
  tbl_landcover19 <- tbl_landcover19[, order(colnames(tbl_landcover19), decreasing = TRUE)] %>%
    relocate(sumPixels, .after = last_col())
  colnames(tbl_landcover19) <- c("Year", "CameraLocation", "Residential",  
                                 "Commercial", "Agriculture", "Forest",  
                                 "XericShrub", "MesicShrub", "XericGrass", 
                                 "MesicGrass", "WoodyWetland", "EmergentWetland", 
                                 "Barren", "Water", "sumPixels")
  #'  Merge annual landcover values together
  tbl_landcover <- rbind(tbl_landcover18, tbl_landcover19) %>%
    #'  Consolidate categories
    #'  Keeping Water to include in percent calculation but will not use for analyses
    mutate(
      Forest =  Forest + WoodyWetland + EmergentWetland,
      MesicGrass = MesicGrass + Barren,
      Developed = Residential + Commercial + Agriculture,
      MesicMix = MesicShrub + MesicGrass,
      ForestMix = Forest + MesicMix,
      ForestMix2 = Forest + MesicShrub
    ) %>%
    dplyr::select(-c(WoodyWetland, EmergentWetland, Barren, Residential, Commercial, Agriculture)) %>%
    relocate(sumPixels, .after = last_col()) %>%
    #'  Calculate percent landcover type within 250m of each camera site
    mutate(
      PercForest = round(Forest/sumPixels, 2),
      PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
      PercForestMix2 = round(ForestMix2/sumPixels, 2),
      PercXericShrub = round(XericShrub/sumPixels, 2),
      PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
      PercXericGrass = round(XericGrass/sumPixels, 2),
      PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
      PercMesicMix = round(MesicMix/sumPixels, 2),      # Cannot be used in conjunction with other Mesic landcover types
      PercWater = round(Water/sumPixels, 2),            # Don't use for analyses (I have better data for H2o)
      PercDeveloped = round(Developed/sumPixels, 2)
    )
  

  
  
  
  
  