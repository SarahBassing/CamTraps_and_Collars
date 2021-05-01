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
  #'    -Human Modified Landscape (1km res)
  #'    -Percent Mixed Forest (within 250m of point)
  #'    -Percent Xeric Grass (within 250m of point)
  #'    -Percent Xeric Shrub (within 205m of point)
  #'    -Distance to nearest road
  #'    -Season
  #'    -Study Area... not yet actually
  #'  ============================================
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(sf)
  library(stars)
  library(rgeos)
  library(raster)
  library(tidyverse)
  
  #' #'  Read in covariate data extracted from other sources
  #' dist2road <- read.csv("./Output/dist2road18-20.csv") %>%   ### NEED TO GENERATE THIS ON LAB COMPUTER
  #'   mutate(
  #'     km2road = dist2road/1000
  #'   ) %>%
  #'   dplyr::select(-X)
  
  #'  Read in spatial data
  wppp_bound <- st_read("./Shapefiles/WPPP_CovariateBoundary", layer = "WPPP_CovariateBoundary")
  #'  Terrain rasters
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m.tif")
  Slope <- raster("./Shapefiles/WA DEM rasters/WPPP_slope_aspect.tif", band = 1)
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
  projection(Slope)
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
  dist2road$ID <- wolf_locs$ID[1]
  
  #'  Function to iterate over all points per species (this will take awhile!)
  RoadDist <- function(locs) {
    dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
    dist2road <- as.data.frame(dist2road)
    dist2road$ID <- locs$ID
    dist2road$time <- locs$time
    dist2road$obs <- c(1:nrow(locs))
    return(dist2road)
  }
  #'  Run each species through function to calculate distance to nearest road
  #'  THIS WILL TAKE FOREVER!!! Started 4/24/21 1230
  # locs <- list(md_locs, elk_locs, wtd_locs, coug_locs, wolf_locs, bob_locs, coy_locs)
  # spp_nearestRd <- lapply(locs, RoadDist)
  
  md_tst <- md_locs[md_locs$ID == "3959MD17_1" | md_locs$ID == "96MD18_279",]
  md_rd_tst <- RoadDist(md_tst)
  
  
  md_rd <- RoadDist(md_locs)
  elk_rd <- RoadDist(elk_locs)
  wtd_rd <- RoadDist(wtd_locs)
  coug_rd <- RoadDist(coug_locs)
  wolf_rd <- RoadDist(wolf_locs)
  bob_rd <- RoadDist(bob_locs)
  coy_rd <- RoadDist(coy_locs)

  save(wolf_rd, file = "./Outputs/Telemetry_covs/wolf_rd.RData")
  
  

  
  
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
    slope <- raster::extract(Slope, locs, df = TRUE)
    modified <- raster::extract(HM, locs, df = TRUE)
    join_covs <- full_join(elevation, slope, by = "ID") %>%
      full_join(modified, by = "ID") %>%
      transmute(
        obs = ID,
        Elev = WPPP_DEM_30m,
        Slope = round(WPPP_slope_aspect, digits = 2),
        HumanMod = WPPP_gHM
      )
    #'  Pull out unique animal/time information
    animal <- as.data.frame(locs) %>%
      select(c(ID, time))
    #'  Merge animal/time information with covariates
    covs <- as.data.frame(cbind(animal, join_covs))
    return(covs)
  }
  
  #'  Run each species through raster extract function
  locs <- list(reproj_md, reproj_elk, reproj_wtd, reproj_coug, reproj_wolf, reproj_bob, reproj_coy)
  spp_extract <- lapply(locs, rast_extract)
  
  #'  Unlist extracted data into separate data frames for each species
  names(spp_extract) <- c("md_extract", "elk_extract", "wtd_extract", "coug_extract", "wolf_extract", "bob_extract", "coy_extract")
  invisible(lapply(names(spp_extract),function(x) assign(x,spp_extract[[x]],.GlobalEnv)))
  # save(spp_extract, file = "./Outputs/Telemetry_covs/spp_extract.RData")  
  
  md_tst <- reproj_md[reproj_md$ID == "3959MD17_1",]
  md_cov_tst <- spp_extract(md_tst)
  
 
  
  
  #'  Function to extract landcover value from each pixel within 250m radius of 
  #'  animal locations using interpolated landcover rasters derived from 
  #'  Cascadia landcover
  landcov250 <- function(locs) {
    pixvals18 <- raster::extract(interp_landcov18, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df18 <- as.data.frame(pixvals18)
    pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
    pixvals_df19 <- as.data.frame(pixvals19)
    #'  Merge together and rename variables
    landcov <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) 
    colnames(landcov) <- c("obs", "landcover_2018", "landcover_2019")
    landcover <- landcov %>%
      mutate(
        landcover_2018 = ifelse(landcover_2018 == "101", "Water", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "121", "Barren", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "201", "EmergentWetland", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "202", "WoodyWetland", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "211", "MesicGrass", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "212", "XericGrass", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "221", "MesicShrub", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "222", "XericShrub", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "230", "Forest", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "310", "Agriculture", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "331", "Commercial", landcover_2018),
        landcover_2018 = ifelse(landcover_2018 == "332", "Residential", landcover_2018),
        landcover_2019 = ifelse(landcover_2019 == "101", "Water", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "121", "Barren", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "201", "EmergentWetland", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "202", "WoodyWetland", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "211", "MesicGrass", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "212", "XericGrass", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "221", "MesicShrub", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "222", "XericShrub", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "230", "Forest", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "310", "Agriculture", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "331", "Commercial", landcover_2019),
        landcover_2019 = ifelse(landcover_2019 == "332", "Residential", landcover_2019),
      )
    animal <- locs %>% 
      select(ID, time) %>%
      mutate(obs = 1:nrow(.))
    landcover_250m <- full_join(animal, landcover, by = "obs")
    
    return(landcover_250m)
  }
 
  # wolf_tst <- wolf_locs[wolf_locs$ID == "W61M_1" | wolf_locs$ID == "W48F_12",]
  # tst <- landcov250(wolf_tst)
  wolf_landcov <- landcov250(wolf_locs)
  save(wolf_landcov, file = "./Outputs/Telemetry_covs/wolf_landcov.RData")
  
  md_landcov_tst <- landcov250(md_tst)
  unique(md_landcov_tst$landcover_2018); unique(md_landcov_tst$landcover_2019)
  
  
  
  
  
  #'  Function to calculate the percent of each landcover type within 250m of 
  #'  each telemetry location
  perc_landcov <- function(landcover, locs) {
    #'  Create base dataframe from animal ID and location
    animal <- as.data.frame(locs) %>%
      select(c(ID, time)) %>%
      mutate(
        obs = 1:nrow(.)
      )
    #'  Count the number of cells in each landcover category by CameraLocation
    tbl_landcover18 <- as.data.frame(landcover) %>%
      select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2018) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2018, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        Forest =  Forest + WoodyWetland + EmergentWetland,
        MesicGrass = MesicGrass + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        PercXericShrub = round(XericShrub/sumPixels, 2),
        PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Join % habitat data to animal location data
      full_join(animal, by = "obs") %>%
      #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix2),
        PercXericShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericShrub),
        PercMesicShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicShrub),
        PercXericGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicMix) 
      ) %>%
      #'  Only retain relevent columns
      select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
                 PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
      #'  Filter out rows with NAs
      filter(!is.na(PercForest))

    
    tbl_landcover19 <- as.data.frame(landcover) %>%
      select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2019) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2019, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        Forest =  Forest + WoodyWetland + EmergentWetland,
        MesicGrass = MesicGrass + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        PercXericShrub = round(XericShrub/sumPixels, 2),
        PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Join % habitat data to animal location data
      full_join(animal, by = "obs") %>%
      #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix2),
        PercXericShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericShrub),
        PercMesicShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicShrub),
        PercXericGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicMix)
      ) %>%
      #'  Only retain relevant columns
      select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
             PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, time, Season) %>%
      #'  Filter out rows with NAs
      filter(!is.na(PercForest))
    
    #'  Merge annual data together so no duplicates
    tbl_landcover <- rbind(tbl_landcover18, tbl_landcover19)

    
    return(tbl_landcover)
  }
  
  #'  Run extracted landcover pixels within 250m of each point AND location data 
  #'  through function to calculate percent habitat type at each location
  # wolf_percHab <- perc_landcov(wolf_landcov, wolf_locs)
  # save(wolf_percHab, file = "./Outputs/Telemetry_covs/wolf_percHab.RData")
  
  md_percHab_tst <- perc_landcov(md_landcov_tst, md_tst)

  
  
  #'  Combine all covariate data together and save for HMM analyses
  full_cov <- md_tst %>%
    select(c(ID, time)) %>%
    full_join(md_cov_tst, by = c("ID", "time")) %>%
    full_join(md_rd_tst, by = c("ID", "time")) %>%
    full_join(md_perHab_tst, by = c("ID", "time")) %>%
    transmute(
      ID = ID,
      AnimalID = AnimalID,
      time = time,
      Season = Season,
      Elev = Elev,
      Slope = Slope,
      NearstRd = dist2road, 
      PercForMix = PercForestMix2,
      PercXGrass = PercXericGrass,
      PercXShrub = PercXericShrub,
      obs = obs
    )
  
  
  
  
  
  animal <- as.data.frame(md_tst) %>%
    select(c(ID, time)) %>%
    mutate(
      obs = 1:nrow(.)
    )
    #'  Count the number of cells in each landcover category by CameraLocation
    tbl_landcover18 <- as.data.frame(md_landcov_tst) %>%
      select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2018) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2018, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        # Forest =  Forest + WoodyWetland + EmergentWetland,
        Forest =  Forest,
        MesicGrass = MesicGrass + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        # PercXericShrub = round(XericShrub/sumPixels, 2),
        # PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        # PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Only retain relevent columns
      # select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub, 
      #        PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix)
      select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, 
             PercMesicGrass, PercMesicMix) %>%
      full_join(animal, by = "obs") %>%
       #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix2),
        # PercXericShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericShrub),
        # PercMesicShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicShrub),
        # PercXericGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicMix) 
      )  %>%
      filter(!is.na(PercForest))
     # colnames(tbl_landcover18) <- c("obs", "sumPixels", "PercForest.18", 
    #                                "PercForestMix.18", "PercForestMix2.18", 
    #                                "PercXericShrub.18", "PercMesicShrub.18", 
    #                                "PercXericGrass.18", "PercMesicGrass.18", 
    #                                "PercMesicMix.18")
    
    tbl_landcover19 <- as.data.frame(md_landcov_tst) %>%
      select(-geometry) %>%
      group_by(obs) %>%
      count(landcover_2019) %>%
      ungroup() %>%
      pivot_wider(names_from = landcover_2019, values_from = n) %>%
      replace(is.na(.), 0) %>% 
      mutate(
        #'  Count number of pixels within 250m of each location
        sumPixels = rowSums(.[2:ncol(.)]),
        #'  Combine similar habitat types
        # Forest =  Forest + WoodyWetland + EmergentWetland,
        Forest =  Forest,
        MesicGrass = MesicGrass + Barren,
        MesicMix = MesicShrub + MesicGrass,
        ForestMix = Forest + MesicMix,
        ForestMix2 = Forest + MesicShrub
      ) %>%
      #'  Calculate percent landcover type within 250m of each camera site
      mutate(
        PercForest = round(Forest/sumPixels, 2),
        PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
        PercForestMix2 = round(ForestMix2/sumPixels, 2),
        # PercXericShrub = round(XericShrub/sumPixels, 2),
        # PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        # PercXericGrass = round(XericGrass/sumPixels, 2),
        PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
        PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
      ) %>%
      #'  Only retain relevent columns
      # select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub, 
      #        PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix)
      select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, 
             PercMesicGrass, PercMesicMix) %>%
      full_join(animal, by = "obs") %>%
      #'  Drop data for year camera was NOT present
      mutate(
        Year = lubridate::year(time),
        Month = lubridate::month(time),
        Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
        Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
        Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
        Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
        Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
        PercForest = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForest),
        PercForestMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix),
        PercForestMix2 = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix2),
        # PercXericShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericShrub),
        # PercMesicShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicShrub),
        # PercXericGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericGrass),
        PercMesicGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicGrass),
        PercMesicMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicMix)
      ) %>%
      filter(!is.na(PercForest))
    # colnames(tbl_landcover19) <- c("obs", "sumPixels", "PercForest.19", 
    #                                "PercForestMix.19", "PercForestMix2.19", 
    #                                "PercXericShrub.19", "PercMesicShrub.19", 
    #                                "PercXericGrass.19", "PercMesicGrass.19", 
    #                                "PercMesicMix.19")
  
    
    #'  Merge landcover and animal ID info
    tbl_landcover <- rbind(tbl_landcover18, tbl_landcover19) 
  
  
  
  
  
  
  #' # cams_reproj <- st_transform(cams, crs(crs(interp_landcov18)))
  #' pixvals18 <- raster::extract(interp_landcov18, tst, factors = TRUE, buffer = 250, df = TRUE)
  #' pixvals_df18 <- as.data.frame(pixvals18)
  #' pixvals19 <- raster::extract(interp_landcov19, tst, factors = TRUE, buffer = 250, df = TRUE)
  #' pixvals_df19 <- as.data.frame(pixvals19)
  #' #'  Merge together and rename variables
  #' landcov <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) 
  #' colnames(landcov) <- c("obs", "landcover_2018", "landcover_2019")
  #' landcover <- landcov %>%
  #'   mutate(
  #'     landcover_2018 = ifelse(landcover_2018 == "101", "Water", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "121", "Barren", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "201", "EmergentWetland", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "202", "WoodyWetland", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "211", "MesicGrass", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "212", "XericGrass", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "221", "MesicShrub", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "222", "XericShrub", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "230", "Forest", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "310", "Agriculture", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "331", "Commercial", landcover_2018),
  #'     landcover_2018 = ifelse(landcover_2018 == "332", "Residential", landcover_2018),
  #'     landcover_2019 = ifelse(landcover_2019 == "101", "Water", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "121", "Barren", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "201", "EmergentWetland", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "202", "WoodyWetland", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "211", "MesicGrass", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "212", "XericGrass", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "221", "MesicShrub", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "222", "XericShrub", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "230", "Forest", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "310", "Agriculture", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "331", "Commercial", landcover_2019),
  #'     landcover_2019 = ifelse(landcover_2019 == "332", "Residential", landcover_2019),
  #'   )
  #' 
  #' animal <- tst %>% 
  #'   select(ID, time) %>%
  #'   mutate(obs = 1:nrow(.))
  #' landcover_250m <- full_join(animal, landcover, by = "obs")
  #' 
  #' 
  #' #'  Count the number of cells in each landcover category by CameraLocation
  #' tbl_landcover18 <- as.data.frame(landcover_250m) %>%
  #'   select(-geometry) %>%
  #'   group_by(obs) %>%
  #'   count(landcover_2018) %>%
  #'   ungroup() %>%
  #'   pivot_wider(names_from = landcover_2018, values_from = n) %>%
  #'   replace(is.na(.), 0) %>% 
  #'   mutate(
  #'     #'  Count number of pixels within 250m of each location
  #'     sumPixels = rowSums(.[2:ncol(.)]),
  #'     #'  Combine similar habitat types
  #'     Forest =  Forest + WoodyWetland + EmergentWetland,
  #'     MesicGrass = MesicGrass + Barren,
  #'     MesicMix = MesicShrub + MesicGrass,
  #'     ForestMix = Forest + MesicMix,
  #'     ForestMix2 = Forest + MesicShrub
  #'   ) %>%
  #'   #'  Calculate percent landcover type within 250m of each camera site
  #'   mutate(
  #'     PercForest = round(Forest/sumPixels, 2),
  #'     PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
  #'     PercForestMix2 = round(ForestMix2/sumPixels, 2),
  #'     PercXericShrub = round(XericShrub/sumPixels, 2),
  #'     PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'     PercXericGrass = round(XericGrass/sumPixels, 2),
  #'     PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'     PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
  #'   ) %>%
  #'   select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub, 
  #'          PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix)
  #' colnames(tbl_landcover18) <- c("obs", "sumPixels", "PercForest.18", 
  #'                                "PercForestMix.18", "PercForestMix2.18", 
  #'                                "PercXericShrub.18", "PercMesicShrub.18", 
  #'                                "PercXericGrass.18", "PercMesicGrass.18", 
  #'                                "PercMesicMix.18")
  #' 
  #' tbl_landcover19 <- as.data.frame(landcover_250m) %>%
  #'   select(-geometry) %>%
  #'   group_by(obs) %>%
  #'   count(landcover_2019) %>%
  #'   ungroup() %>%
  #'   pivot_wider(names_from = landcover_2019, values_from = n) %>%
  #'   replace(is.na(.), 0) %>% 
  #'   mutate(
  #'     #'  Count number of pixels within 250m of each location
  #'     sumPixels = rowSums(.[2:ncol(.)]),
  #'     #'  Combine similar habitat types
  #'     Forest =  Forest + WoodyWetland + EmergentWetland,
  #'     MesicGrass = MesicGrass + Barren,
  #'     MesicMix = MesicShrub + MesicGrass,
  #'     ForestMix = Forest + MesicMix,
  #'     ForestMix2 = Forest + MesicShrub
  #'   ) %>%
  #'   #'  Calculate percent landcover type within 250m of each camera site
  #'   mutate(
  #'     PercForest = round(Forest/sumPixels, 2),
  #'     PercForestMix = round(ForestMix/sumPixels,2),     # Cannot be used in conjunction with Forest or any Mesic landcover types
  #'     PercForestMix2 = round(ForestMix2/sumPixels, 2),
  #'     PercXericShrub = round(XericShrub/sumPixels, 2),
  #'     PercMesicShrub = round(MesicShrub/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'     PercXericGrass = round(XericGrass/sumPixels, 2),
  #'     PercMesicGrass = round(MesicGrass/sumPixels, 2),  # Cannot be used in conjunction with MesicMix
  #'     PercMesicMix = round(MesicMix/sumPixels, 2)      # Cannot be used in conjunction with other Mesic landcover types
  #'   ) %>%
  #'   select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub, 
  #'          PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix)
  #' colnames(tbl_landcover19) <- c("obs", "sumPixels", "PercForest.19", 
  #'                                "PercForestMix.19", "PercForestMix2.19", 
  #'                                "PercXericShrub.19", "PercMesicShrub.19", 
  #'                                "PercXericGrass.19", "PercMesicGrass.19", 
  #'                                "PercMesicMix.19")
  #' 
  #' tbl_landcover <- full_join(animal, tbl_landcover18, by = "obs") %>%
  #'   full_join(tbl_landcover19, by = c("obs", "sumPixels"))

  
  
  
  
  