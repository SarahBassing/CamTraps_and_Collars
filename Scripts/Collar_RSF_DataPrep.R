  #'  ============================================
  #'  Resource Selection Functions (cam vs collar analysis)
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  ============================================
  #'  Script to randomly sample "available" points from the home range of each 
  #'  collared individual and build resource selection function models for each
  #'  species. This will be compared to the occupancy and HMM models to evaluate
  #'  habitat associations derived from different types of data and whether those
  #'  associations vary with animal behavior.
  #'  
  #'  Cleaned telemetry and covariate data were prepared for with the
  #'  Collar_Movement_DataPrep.R and Collar_Covariate_Extraction.R scripts 
  #'  which take FOREVER to run. 
  #'  ============================================
  
  #'  Load packages for selecting available points
  library(tidyverse)
  library(sp)
  library(lme4)
  library(adehabitatHR)
  
  #'  Load telemetry data
  # load("./Outputs/Telemetry_tracks/MD_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/MD_wtr_track.RData")
  # load("./Outputs/Telemetry_tracks/ELK_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/ELK_wtr_track.RData")
  # load("./Outputs/Telemetry_tracks/WTD_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/WTD_wtr_track.RData")
  # load("./Outputs/Telemetry_tracks/COUG_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/COUG_wtr_track.RData")
  # load("./Outputs/Telemetry_tracks/WOLF_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/WOLF_wtr_track.RData")
  # load("./Outputs/Telemetry_tracks/BOB_smr_track.RData")
  # load("./Outputs/Telemetry_tracks/BOB_wtr_track.RData")
  load("./Outputs/Telemetry_tracks/COY_smr_track.RData")
  load("./Outputs/Telemetry_tracks/COY_wtr_track.RData")
  
  
  #'  Generate random "Available" locations for each individual
  #'  =========================================================
  #'  Function to pull out unique animal IDs
  unq_id <- function(locs) {
    animal <- as.data.frame(locs$AnimalID) %>%
      unique()
    colnames(animal) <- "AnimalID"
    return(animal)
  }
  coy_ID_smr <- unq_id(COY_smr_track)
  coy_ID_wtr <- unq_id(COY_wtr_track)
  
  #'  Function to split data into list of dataframes by unique animal ID and year
  #'  (Use AnimalID if ignoring year aspect of data)
  split_animal <- function(locs) {
    ind_animal <- group_split(locs, locs$FullID) #FullID so it's by year too (important for landcover data extraction)
    return(ind_animal)
  }
  coy_split_smr <- split_animal(COY_smr_track)
  coy_split_wtr <- split_animal(COY_wtr_track)
  
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
  spp_list <- list(COY_smr_track, COY_wtr_track)
  mean_locs <- lapply(spp_list, mean_obs)
  #'  Calculate mean number of used locations for all species
  mean_used <- mean(unlist(mean_locs))
  #'  RSF literature suggests 1:20 ration used:available
  navailable <- mean_used*20
  
  #'  Set projection for spatial analyses
  sa_proj <- crs("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Function to randomly select "Available" points within each animal's home
  #'  range per season. 
  avail_pts <- function(locs, plotit = F) {
    #'  1. Make each animal's locations spatial
    #'  ---------------------------------------------------------
    locs_sp <- SpatialPoints(locs[,c("x", "y")], proj4string = sa_proj)
    #'  Plot to make sure step 1 worked
    if(plotit) {
      plot(locs_sp, col = "blue", pch = 19)
    }
    
    #'  2. Create KUDs for each animal following Bivariate normal utilization distributions
    #'  ----------------------------------------------------------
    UD <- kernelUD(locs_sp, kern = "bivnorm")
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
  coy_smr_df <- lapply(coy_split_smr, avail_pts, F)
  coy_wtr_df <- lapply(coy_split_wtr, avail_pts, F)
  #'  Convert to dataframes instead of lists
  coy_smr_df <- do.call(rbind.data.frame, coy_smr_df)
  coy_wtr_df <- do.call(rbind.data.frame, coy_wtr_df)
  
  #'  Gather into one big list
  coy_available <- list(coy_smr_df, coy_wtr_df)
  
  #'  Save
  save(coy_available, file = paste0("./Outputs/RSF_available_pts/coy_available_", Sys.Date(), ".RData"))
 

  
  #'  Extract covariate data for each available location
  #'  ==================================================
  #'  This will take awhile!

  #'  Load packages for covariate extraction
  library(sf)
  library(stars)
  library(rgeos)
  library(raster)
  library(parallel)
  library(doParallel)
  library(future.apply)
  library(tidyverse)
  
  #'  Load location data
  load("./Outputs/RSF_available_pts/coy_available_2021-05-12.RData")
  
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
  formix2prop18 <- raster("./Shapefiles/Cascadia_layers/forestmix2prop_18.tif")
  formix2prop19 <- raster("./Shapefiles/Cascadia_layers/forestmix2prop_19.tif")
  xgrassprop18 <- raster("./Shapefiles/Cascadia_layers/xgrassprop_18.tif")
  xgrassprop19 <- raster("./Shapefiles/Cascadia_layers/xgrassprop_19.tif")
  xshrubprop18 <- raster("./Shapefiles/Cascadia_layers/xshrubprop_18.tif")
  xshrubprop19 <- raster("./Shapefiles/Cascadia_layers/xshrubprop_19.tif")
  roads <- st_read("./Shapefiles/Cascadia_layers/roadsForTaylor", layer = "roadsForTaylor")
  #'  Human density and human modified landscapes
  HM <- raster("./Shapefiles/Additional_WPPP_Layers/WPPP_gHM.tif")
  
  #'  Create raster stacks of 2018 and  2019 data
  perc_stack18 <- stack(formix2prop18, xgrassprop18, xshrubprop18)
  perc_stack19 <- stack(formix2prop19, xgrassprop19, xshrubprop19)
  
  #'  Identify projections & resolutions of relevant features
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Reproject road shapefile to match animal location projection
  road_reproj <- st_transform(roads, crs = st_crs(sa_proj))
  projection(road_reproj)

  
  #'  Function to make available points data a spatial sf object
  spatial_locs <- function(locs) {
    sf_locs <- st_as_sf(locs, coords = c("x", "y"), crs = sa_proj)
    return(sf_locs)
  }
  sf_locs <- lapply(coy_available, spatial_locs)
  # coy_smr_locs <- spatial_locs(coy_smr_df)
  # coy_wrt_locs <- spatial_locs(coy_wtr_df)
  

  
  #'  COVARIATE EXTRACTION & CALCULATIONS  
  #'  ===========================================
  #'  Takes forever but running in parallel helps 
  
  #'  Monitor time
  start.time <- Sys.time()
  
  #'  Setup script to run in parallel
  #'  Extract covariates for each species at once
  #'  Identify how many cores I want to use
  detectCores(logical = FALSE)
  cl <- parallel::makeCluster(20)
  #'  Run in parallel on local computer with specified number of cores
  plan(cluster, workers = cl)
  
  #'  Master function to extract and manipulate covaraite data for each species
  cov_extract <- function(locs) {
    
    #'  1. Extract data from unprojected rasters
    #'  ----------------------------------------
    #'  Reproject location data to match rasters (WGS84)
    reproj_locs <- st_transform(locs, crs = st_crs(wgs84))
    #'  Extract covariates for each location
    elevation <- raster::extract(dem, reproj_locs, df = TRUE)
    slope <- raster::extract(Slope, reproj_locs, df = TRUE)
    modified <- raster::extract(HM, reproj_locs, df = TRUE)
    #'  Merge into a single data frame of covariates
    join_covs <- full_join(elevation, slope, by = "ID") %>%
      full_join(modified, by = "ID") %>%
      transmute(
        obs = ID,
        Elev = WPPP_DEM_30m,
        Slope = round(WPPP_slope_aspect, digits = 2),
        HumanMod = WPPP_gHM
      )
    #'  Make location and animal ID information non-spatial
    animal <- as.data.frame(locs) %>%
      dplyr::select(-geometry)
    #'  Merge animal/time information with covariates
    covs <- as.data.frame(cbind(animal, join_covs))

    
    #'  2. Extract data from roads shapefile & calculate distance to nearest road
    #'     for each location
    #'  ------------------------------------------------------------------------
    dist2road <- sapply(1:nrow(locs), function(x) min(st_distance(road_reproj, locs[x, ])))
    dist2road <- as.data.frame(dist2road)
    dist2road$ID <- locs$ID
    dist2road$obs <- c(1:nrow(locs))
    #'  Append to covariate data frame
    covs <- covs %>%
      full_join(dist2road, by = c("obs", "ID"))
    
    
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
    #' #'  Extract landcover value from each pixel within 250m radius of locations
    #' #'  using interpolated landcover rasters derived from Cascadia landcover
    #' pixvals18 <- raster::extract(interp_landcov18, locs, factors = TRUE, buffer = 250, df = TRUE)
    #' pixvals_df18 <- as.data.frame(pixvals18)
    #' pixvals19 <- raster::extract(interp_landcov19, locs, factors = TRUE, buffer = 250, df = TRUE)
    #' pixvals_df19 <- as.data.frame(pixvals19)
    #' #'  Merge together and rename variables
    #' landcov <- cbind(pixvals_df18, pixvals_df19$interpolated_landcover_2019) 
    #' colnames(landcov) <- c("obs", "landcover_2018", "landcover_2019")
    #' #'  Rename categories so they're easier to work with
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
    #' #'  Add animal ID info to landcover data for further manipulation
    #' animal <- mutate(animal, obs = 1:nrow(.))
    #' landcover_250m <- full_join(animal, landcover, by = "obs")
    #' #'  Count the number of cells in each landcover category per location
    #' #'  2018 landcover and telemetry locations only
    #' tbl_landcover18 <- as.data.frame(landcover_250m) %>%
    #'   # select(-geometry) %>%
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
    #'     MesicGrass = MesicGrass,# + Barren,
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
    #'   #'  Join % habitat data to animal location data
    #'   full_join(animal, by = "obs") %>%
    #'   #'  Drop data for year camera was NOT present
    #'   mutate(
    #'     Year = lubridate::year(time),
    #'     Month = lubridate::month(time),
    #'     Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
    #'     Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
    #'     Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
    #'     Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
    #'     Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
    #'     Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
    #'     PercForest = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForest),
    #'     PercForestMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix),
    #'     PercForestMix2 = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercForestMix2),
    #'     PercXericShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericShrub),
    #'     PercMesicShrub = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicShrub),
    #'     PercXericGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercXericGrass),
    #'     PercMesicGrass = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicGrass),
    #'     PercMesicMix = ifelse(Season == "Summer19" | Season == "Winter1920", NA, PercMesicMix) 
    #'   ) %>%
    #'   #'  Only retain relevant columns
    #'   dplyr::select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
    #'                 PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, Season) %>%
    #'   #'  Filter out rows with NAs
    #'   filter(!is.na(PercForest))
    #' #'  Repeat for 2019 landcover data and telemetry locations
    #' tbl_landcover19 <- as.data.frame(landcover_250m) %>%
    #'   # select(-geometry) %>%
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
    #'     MesicGrass = MesicGrass, # + Barren,
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
    #'   #'  Join % habitat data to animal location data
    #'   full_join(animal, by = "obs") %>%
    #'   #'  Drop data for year camera was NOT present
    #'   mutate(
    #'     Year = lubridate::year(time),
    #'     Month = lubridate::month(time),
    #'     Season = ifelse(Year == 2018 & Month < 10, "Summer18", NA),
    #'     Season = ifelse(Year == 2018 & Month > 11, "Winter1819", Season),
    #'     Season = ifelse(Year == 2019 & Month < 4, "Winter1819", Season),
    #'     Season = ifelse(Year == 2019 & Month > 6 & Month < 10, "Summer19", Season),
    #'     Season = ifelse(Year == 2019 & Month > 11, "Winter1920", Season),
    #'     Season = ifelse(Year == 2020 & Month < 4, "Winter1920", Season),
    #'     PercForest = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForest),
    #'     PercForestMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix),
    #'     PercForestMix2 = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercForestMix2),
    #'     PercXericShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericShrub),
    #'     PercMesicShrub = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicShrub),
    #'     PercXericGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercXericGrass),
    #'     PercMesicGrass = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicGrass),
    #'     PercMesicMix = ifelse(Season == "Summer18" | Season == "Winter1819", NA, PercMesicMix)
    #'   ) %>%
    #'   #'  Only retain relevant columns
    #'   dplyr::select(obs, sumPixels, PercForest, PercForestMix, PercForestMix2, PercXericShrub,
    #'                 PercMesicShrub, PercXericGrass, PercMesicGrass, PercMesicMix, ID, Season) %>%
    #'   #'  Filter out rows with NAs
    #'   filter(!is.na(PercForest))
    #' #'  Merge annual landcover data together so no duplicates
    #' percHab <- rbind(tbl_landcover18, tbl_landcover19)
    
    #'  4. Join all covatiates together & clean up for inclusion in HMM
    telem_covs <- covs %>%
      full_join(percHab, by = c("obs", "ID", "Season")) %>%
      transmute(
        ID = ID,
        Season = Season,
        Year = ifelse(Season == "Summer18" | Season == "Winter1819", "Year1", "Year2"),
        Elev = Elev,
        Slope = Slope,
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
        Area = ifelse(grepl("WTD", ID), "NE", Area))
    
    return(telem_covs)
    
  }
  
  #'  Run list of species location data through function in parallel
  #'  This will take AWHILE even in parallel
  # spp_avail_covs <- lapply(sf_locs, cov_extract) # non-parallel approach
  # spp_avail_covs <- future_lapply(sf_locs, cov_extract)
  coy_avail_covs <- future_lapply(sf_locs, cov_extract)
  
  
  #'  End time keeping
  end.time <- Sys.time()
  #'  Stop running in parallel
  parallel::stopCluster(cl)
  #'  How long did this take?
  difftime(end.time, start.time, units = "hours")
  
  #' #'  Add study area to wolf data
  #' #'  No easy way of doing this because ID not associated with WPPP study areas
  #' #'  Double check lists 9 & 10 are wolf summer & winter data
  #' spp_avail_covs[[9]]$Area <- "NE"   
  #' spp_avail_covs[[9]] <- mutate(spp_avail_covs[[9]], 
  #'                               Area = ifelse(grepl("W61M", ID), "OK", Area),  #double check no "W71F" in here
  #'                               Area = ifelse(grepl("W88M", ID), "OK", Area),
  #'                               Area = ifelse(grepl("W93M", ID), "OK", Area),
  #'                               Area = ifelse(grepl("W94M", ID), "OK", Area))
  #' spp_avail_covs[[10]]$Area <- "NE"   
  #' spp_avail_covs[[10]] <- mutate(spp_avail_covs[[10]], 
  #'                                Area = ifelse(grepl("W61M", ID), "OK", Area),
  #'                                Area = ifelse(grepl("W88M", ID), "OK", Area),
  #'                                Area = ifelse(grepl("W93M", ID), "OK", Area),
  #'                                Area = ifelse(grepl("W94M", ID), "OK", Area))
  
  #'  Save and hope you never have to run this again!
  # save(spp_avail_covs, file = paste0("./Outputs/RSF_available_pts/spp_avail_covs_", Sys.Date(), ".RData"))
  save(coy_avail_covs, file = paste0("./Outputs/RSF_available_pts/coy_avail_covs_", Sys.Date(), ".RData"))
  
  
  
  
  
  
  
  