  #'  =================================================
  #'  Maps and figures for camera vs collar manuscript
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  =================================================
  #'  1. Plot study area map with camera locations and inset map showing study 
  #'     area locations in relation to Washington State
  #'  2. Plot collaring effort and MCPs based on all collars per species and 
  #'     study area
  #'  3. Plot occupancy and RSF results
  #'  =================================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(png)
  library(RCurl)
  library(RColorBrewer)
  library(rphylopic)
  library(sf)
  library(raster)
  library(tidyverse)
  
  
  #'  Get some basics pulled together to be used across most figures
  #'  -----------------------------
  ####  Spatial Data for Mapping  ####
  #'  -----------------------------
  #'  Define projections
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Read in study area data and reproject
  WA <- st_read("./Shapefiles/Washington_State_Boundary/WA_State_Geospatial_Open_Data", layer = "WA_State_Boundary") %>%
    st_transform(crs = sa_proj)
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK_SA$NAME <- "Okanogan"
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE_SA$NAME <- "Northeast"
  
  projection(WA)
  projection(OK_SA)
  extent(OK_SA)
  extent(NE_SA)
  
  #'  Reduce DEM raster resolution and prep new raster for ggplot
  # dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m_reproj.tif")
  # projection(dem)
  # dem_low <- aggregate(dem, fact = 10)
  # writeRaster(dem_low, file = "./Shapefiles/WA DEM rasters/dem_reproj_low", format = "GTiff")
  dem_low <- raster("./Shapefiles/WA DEM rasters/dem_reproj_low.tif")
  dem_p_low <- rasterToPoints(dem_low)
  dem_p_df <- as.data.frame(dem_p_low)
  colnames(dem_p_df) <- c("x", "y", "value")
  
  #'  ------------------------
  ####  Species silhouettes  ####
  #'  ------------------------
  #'  Silhouettes for each species from PhyloPic in two different formats (PNG & rastergrob)
  cougurl <- "http://phylopic.org/assets/images/submissions/3f8eff77-2868-4121-8d7d-a55ebdd49e04.64.png"
  cougimg <- readPNG(getURLContent(cougurl))
  cougg <- rasterGrob(cougimg, interpolate = TRUE)
  wolfurl <- "http://phylopic.org/assets/images/submissions/8cad2b22-30d3-4cbd-86a3-a6d2d004b201.512.png"
  wolfimg <- readPNG(getURLContent(wolfurl))
  wolfg <- rasterGrob(wolfimg, interpolate = TRUE)
  boburl <- "http://phylopic.org/assets/images/submissions/ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa.512.png"
  bobimg <- readPNG(getURLContent(boburl))
  bobg <- rasterGrob(bobimg, interpolate = TRUE)
  coyurl <- "http://phylopic.org/assets/images/submissions/5a0398e3-a455-4ca6-ba86-cf3f1b25977a.512.png"
  coyimg <- readPNG(getURLContent(coyurl)) 
  coyg <- rasterGrob(coyimg, interpolate = TRUE)
  mdurl <- "http://phylopic.org/assets/images/submissions/f889b336-9e67-4154-bc96-db4095a55be2.512.png"
  mdimg <- readPNG(getURLContent(mdurl))
  mdg <- rasterGrob(mdimg, interpolate = TRUE)
  elkmurl <- "http://phylopic.org/assets/images/submissions/72f2f997-e474-4caf-bbd5-72fc8dbcc40d.512.png"
  elkmimg <- readPNG(getURLContent(elkmurl))
  elkmg <- rasterGrob(elkmimg, interpolate = TRUE)
  elkfurl <- "http://phylopic.org/assets/images/submissions/97f83f5e-9afe-4ce8-812e-337f506ca841.512.png"
  elkfimg <- readPNG(getURLContent(elkfurl))
  elkfg <- rasterGrob(elkfimg, interpolate = TRUE)
  wtdurl <- "http://phylopic.org/assets/images/submissions/56f6fdb2-15d0-43b5-b13f-714f2cb0f5d0.512.png"
  wtdimg <- readPNG(getURLContent(wtdurl))
  wtdg <- rasterGrob(wtdimg, interpolate = TRUE)


  ####  1. Map study area and camera locations  ####
  #'  ==============================================
  #'  Read in camera locations
  station_covs <- read.csv("./Data/Camera_Station18-20_Covariates_2021-04-25.csv")
  CameraLocation <- station_covs$CameraLocation
  Year <- station_covs$Year
  
  #'  Make camera location data spatial and reproject to study area projection
  cams <- st_as_sf(station_covs[,6:8], coords = c("Longitude", "Latitude"), crs = wgs84)
  cams_reproj <- st_transform(cams, crs = sa_proj)
  cams_reproj$Year <- Year
  cams_reproj <- mutate(cams_reproj, Year = ifelse(Year == "Year1", "2018-2019", "2019-2020" ))
  
  #'  Plot state of WA with study areas
  #'  https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
  WA_SA_map <- ggplot() + 
    geom_sf(data = WA, fill = "gray95") +
    #' #'  Label map of WA with "Washington State"
    #' geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.8, vjust = 3), size = 12) +
    geom_sf(data = OK_SA, fill = "grey25", color = "grey20") +
    geom_sf(data = NE_SA, fill = "grey25", color = "grey20") +
    #'  Get rid of lines and axis labels
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  plot(WA_SA_map)

  #'  Plot study areas against DEM with camera locations
  cam_SA_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) + 
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area outlines and label with their names
    geom_sf(data = OK_SA, fill = NA, color="black", size = 0.80) +
    #'  Note the hjust & vjust need to change based on font size and coord_sf
    #'  DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
    geom_sf_text(data = OK_SA, aes(label = NAME, hjust = 1.3, vjust = -9.5), size = 4.5) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf_text(data = NE_SA, aes(label = NAME, hjust = 0.75, vjust = -6.75), size = 4.5) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = cams_reproj, aes(color = Year), shape = 16) +
    #'  Change camera data aesthetics (make sure it's colorblind friendly)
    scale_discrete_manual(aesthetics = "color", values = c("#a6611a", "#018571")) + #c("#dfc27d", "#80cdc1") #c("#601A4A", "#63ACBE")
    labs(colour = "Camera\ndeployment") +
    #'  Constrain plot to two study areas plus some room on the side & bottom
    coord_sf(xlim = c(480000.0, 810000.0), ylim = c(39000.0, 218000.0), expand = FALSE) +
    #'  Constrain map to just the two study areas only
    # coord_sf(xlim = c(504659.0, 781979.9), ylim = c(102808.3, 211000.4)) +
    #'  Get rid of lines and axis names
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank()) +
    #'  Add north arrow
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.3, "in"),
                           style = north_arrow_fancy_orienteering) +
    #'  Add scale bar (be sure to double check the scale)
    annotation_scale(location = "bl", width_hint = 0.5)
  plot(cam_SA_map)
  
  #'  Build plot with map of study areas and inset map of WA
  #'  https://upgo.lab.mcgill.ca/2019/12/13/making-beautiful-maps/
  #'  Requires "cowplot" package
  #'  Don't use png or other calls to save while plotting- formatting gets messed up
  #'  Use export option in Plot window and formatting holds
  # png(file = "./Outputs/Figures/Maps/StudyAreas_Cameras1820.png",
  #     width = 1000, height = 691)
  StudyArea_Map <- ggdraw(cam_SA_map) +
    draw_plot(
      {
        WA_SA_map +
        #'  Label map of WA with "Washington State"
        #'  hjust & vjust will depend on inset map's width/height specified below
        geom_sf_text(data = WA, aes(label = JURISDIC_3, hjust = 0.5, vjust = 2)) 
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.60,
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.19,
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.25,
      height = 0.25)
  plot(StudyArea_Map)
  # dev.off()
  
  
  
  ####  2. Maps of capture effort and MCPs  ####
  #'  ==========================================
  #'  Read in telemetry data & identify 1st location of collars included in RSFs
  #'  Track data of animals included in RSF
  load("./Outputs/Telemetry_tracks/spp_all_tracks_noDispersal.RData") #spp_all_tracks_noDispMig
  #'  Cleaned telemetry data of animals from start of collar deployment (no 
  #'  thinning, filtering, or truncating but bad locations removed)
  load("./Data/Collar_AllSpecies_AllLocations_Clean.RData")
  #'  Truncated locations 2 weeks post-capture
  load("./Data/Collar_AllSpecies_AllLocations_Truncated.RData")
  #'  Separate meso data into different species
  coy_trunk <- droplevels(clean_data[[6]][clean_data[[6]]$Species == "Coyote",])
  bob_trunk <- droplevels(clean_data[[6]][clean_data[[6]]$Species == "Bobcat",])
  
  #'  Function to identify each animal included in any of the RSFs (based on 
  #'  track data), pull out the 1st location from the truncated location data 
  #'  2 weeks post-capture, and filter those observations to just the animals
  #'  included in the RSFs.
  first_loc <- function(smrtrack, wtrtrack, spptrunk) {
    #  Identify the unique collars included in the RSFs
    smrcollars <- as.data.frame(unique(smrtrack$AnimalID)) 
    colnames(smrcollars) <- "ID"
    wtrcollars <- as.data.frame(unique(wtrtrack$AnimalID)) 
    colnames(wtrcollars) <- "ID"
    rsfcollars <- unique(rbind(smrcollars, wtrcollars))
    
    #'  Retain the first location for each collared animal in the truncated data
    collar_start <- spptrunk %>% 
      group_by(ID) %>%
      slice(1) %>%
      ungroup() %>%
      dplyr::select(c(ID, StudyArea, Longitude, Latitude, Finaldt)) %>%
    #'  Retain only location data of animals included in RSFs
    filter(ID %in% rsfcollars$ID) %>%
    #'  Make spatial
    st_as_sf(., coords = c("Longitude", "Latitude"), crs = wgs84) %>%
    st_transform(., crs = sa_proj)
    
    return(collar_start)
  }
  #'  Run seasonal track & full truncated data for each species through function
  #'  NOTE: supplying truncated data for mule deer since start of clean telemetry
  #'  locations are affected by where they were taken for capture processing
  #'  All other species use the cleaned but not truncated data since the start 
  #'  of their collar data is more representative of where they were captured.
  mdstart <- first_loc(spp_all_tracks[[1]], spp_all_tracks[[2]], trunk_data[[1]])
  elkstart <- first_loc(spp_all_tracks[[3]], spp_all_tracks[[4]], clean_data[[2]])
  wtdstart <- first_loc(spp_all_tracks[[5]], spp_all_tracks[[6]], clean_data[[3]])
  cougstart <- first_loc(spp_all_tracks[[7]], spp_all_tracks[[8]], clean_data[[4]])
  wolfstart <- first_loc(spp_all_tracks[[9]], spp_all_tracks[[10]], clean_data[[5]])
  bobstart <- first_loc(spp_all_tracks[[11]], spp_all_tracks[[12]], bob_trunk)
  coystart <- first_loc(spp_all_tracks[[13]], spp_all_tracks[[14]], coy_trunk)
  
  #'  Read in MCP polygons per species
  md_poly <- st_read("./Outputs/MCPs", layer = "md_poly")
  elk_poly <- st_read("./Outputs/MCPs", layer = "elk_poly")
  wtd_poly <- st_read("./Outputs/MCPs", layer = "wtd_poly")
  coug_NE_poly <- st_read("./Outputs/MCPs", layer = "coug_NE_poly")
  coug_OK_poly <- st_read("./Outputs/MCPs", layer = "coug_OK_poly")
  wolf_NE_poly <- st_read("./Outputs/MCPs", layer = "wolf_NE_poly")
  wolf_OK_poly <- st_read("./Outputs/MCPs", layer = "wolf_OK_poly")
  bob_NE_poly <- st_read("./Outputs/MCPs", layer = "bob_NE_poly")
  bob_OK_poly <- st_read("./Outputs/MCPs", layer = "bob_OK_poly")
  coy_NE_poly <- st_read("./Outputs/MCPs", layer = "coy_NE_poly")
  coy_OK_poly <- st_read("./Outputs/MCPs", layer = "coy_OK_poly")
  
  #'  Plot capture effort and sampling area of RSFs
  #'  ----------------------------------------
  ####  Capture Locations & MCPs per Species ####
  #'  ----------------------------------------
  md_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = md_poly, fill = NA, color = "blue", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = mdstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    #'  Change camera data aesthetics (make sure it's colorblind friendly)
    # scale_discrete_manual(aesthetics = "color", values = c("#601A4A", "#63ACBE")) +
    # labs(colour = "Camera\ndeployment") +
    #'  Constrain plot to OK study area only  
    coord_sf(xlim = c(480000.0, 600000.0), ylim = c(102000.0, 250000.0)) + #ylim = c(102000.0, 240000.0)
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(mdg, xmin = 480000.0, xmax = 510000.0, ymin = 225000, ymax = 255000) #xmin = 480000.0, xmax = 505000.0, 
  
  elk_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    # geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = elk_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = elkstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    #'  Constrain plot to NE study area only 
    coord_sf(xlim = c(680000.0, 780000.0), ylim = c(102000.0, 250000.0)) + #ylim = c(102000.0, 200000.0)
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(elkmg, xmin = 680000.0, xmax = 710000.0, ymin = 225000, ymax = 255000) #xmin = 680000.0, xmax = 695000.0, ymin = 180000, ymax = 205000
  
  wtd_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    # geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = wtd_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = wtdstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    #'  Constrain plot to NE study area only 
    coord_sf(xlim = c(680000.0, 780000.0), ylim = c(102000.0, 250000.0)) + #ylim = c(102000.0, 200000.0)
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(wtdg, xmin = 680000.0, xmax = 710000.0, ymin = 225000, ymax = 255000) #xmin = 680000.0, xmax = 695000.0, ymin = 180000, ymax = 205000
  
  coug_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = coug_OK_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf(data = coug_NE_poly, fill = NA, color = "blue", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = cougstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    #'  Constrain plot to both study areas
    coord_sf(xlim = c(490000.0, 780000.0), ylim = c(102000.0, 250000.0)) +
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(cougg, xmin = 740000.0, xmax = 780000.0, ymin = 220000, ymax = 260000)
  
  wolf_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = wolf_OK_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf(data = wolf_NE_poly, fill = NA, color = "blue", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = wolfstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    #'  Constrain plot to both study areas
    coord_sf(xlim = c(490000.0, 780000.0), ylim = c(102000.0, 250000.0)) +
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(wolfg, xmin = 750000.0, xmax = 780000.0, ymin = 220000, ymax = 260000)
  
  bob_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = bob_OK_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf(data = bob_NE_poly, fill = NA, color = "blue", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = bobstart, color = "black", shape = 16) +
    #'  Drop x and y-axis labels
    xlab("") + ylab("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    #'  Constrain plot to both study areas
    coord_sf(xlim = c(490000.0, 780000.0), ylim = c(102000.0, 250000.0)) +
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(bobg, xmin = 750000.0, xmax = 780000.0, ymin = 220000, ymax = 260000)
  
  coy_map <- ggplot() +
    geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
    #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
    scale_alpha(range = c(0.3, 0.8)) +
    #'  Change colors of the raster
    scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
    #'  Add study area and MCP polygons
    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
    geom_sf(data = coy_OK_poly, fill = NA, color = "blue", size = 0.80) +
    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
    geom_sf(data = coy_NE_poly, fill = NA, color = "blue", size = 0.80) +
    #'  Add camera locations and vary color by deployment year
    geom_sf(data = coystart, color = "black", shape = 16) +
    #'  Give x and y-axis labels
    xlab("Longitude") + ylab("Latitude") +
    #'  Constrain plot to both study areas
    coord_sf(xlim = c(490000.0, 780000.0), ylim = c(102000.0, 250000.0)) + 
    #'  Add rasterized silhouette in top left corner (min & max based on plot coordinates)
    annotation_custom(coyg, xmin = 760000.0, xmax = 780000.0, ymin = 220000, ymax = 260000)
  
  #'  Check 'em out!
  plot(md_map)
  plot(elk_map)
  plot(wtd_map)
  plot(coug_map)
  plot(wolf_map)
  plot(bob_map)
  plot(coy_map)
  
  ####  Panel of maps  ####
  capture_fig <- (md_map + elk_map + wtd_map) / (bob_map + coug_map) / (coy_map + wolf_map) 
  plot(capture_fig)
  
  #'  Save
  ggsave("./Outputs/Figures/Maps/CaptureEffort_fig.png", capture_fig)
  ggsave("./Outputs/Figures/Maps/MuleDeerCaptureEffort.png", md_map)
  ggsave("./Outputs/Figures/Maps/ElkCaptureEffort.png", elk_map)
  ggsave("./Outputs/Figures/Maps/WTDeerCaptureEffort.png", wtd_map)
  ggsave("./Outputs/Figures/Maps/CougarCaptureEffort.png", coug_map)
  ggsave("./Outputs/Figures/Maps/WolfCaptureEffort.png", wolf_map)
  ggsave("./Outputs/Figures/Maps/BobcatCaptureEffort.png", bob_map)
  ggsave("./Outputs/Figures/Maps/CoyoteCaptureEffort.png", coy_map)
  
  
  #  Also think about changing the axis- add a little space btwn plot and tick marks?
  #  Add a legend for dots and MCP?
  

 #' #'  TROUBLE SHOOTING MASSIVE MCPs  
 #'  bob_smr <- spp_all_tracks[[11]]
 #'  bob_wtr <- spp_all_tracks[[12]]
 #'  bob_all <- rbind(bob_smr, bob_wtr)
 #'  bob_sf <- st_as_sf(bob_all, coords = c("Long", "Lat"), crs = wgs84) %>%
 #'    st_transform(crs = sa_proj)
 #'  bob_all_ID <- as.data.frame(unique(bob_all$AnimalID))
 #' 
 #'  bob_trunk <- st_as_sf(bob_trunk, coords = c("Longitude", "Latitude"), crs = wgs84)
 #'  bob_skinny_sf <- st_transform(bob_trunk, crs = sa_proj)
 #'  bob_skinny_NE <- bob_skinny_sf[bob_skinny_sf$StudyArea == "NE",]
 #'  bob_skinny_OK <- bob_skinny_sf[bob_skinny_sf$StudyArea == "OK",]
 #'  bob_locs <- read.csv("./Data/BobcatData_AllLocations.csv")
 #'  bob_collars <- as.data.frame(unique(bob_locs$ID))
 #' 
 #'  bob_NE_MCP <- st_as_sf(bob_NE_mcp)
 #'  bob_NE_poly <- st_as_sf(bob_NE_poly)
 #'  bob_OK_MCP <- st_as_sf(bob_OK_mcp)
 #'  bob_OK_poly <- st_as_sf(bob_OK_poly)
 #' 
 #' ggplot() +
 #'    #' geom_raster(data = dem_p_df, aes(x = x, y = y, fill = value, alpha = value), show.legend = FALSE) +
 #'    #' #'  alpha adjusts transparency of the raster (can also just set it range = 0.7)
 #'    #' scale_alpha(range = c(0.3, 0.8)) +
 #'    #' #'  Change colors of the raster
 #'    #' scale_fill_gradient2(low = "grey95", high = "tan4") + #gray20
 #'    #'  Add study area and MCP polygons
 #'    geom_sf(data = OK_SA, fill = NA, color = "black", size = 0.80) +
 #'    # geom_sf(data = bob_OK_MCP, fill = NA, color = "red", size = 0.8) +
 #'    geom_sf(data = bob_OK_poly, fill = NA, color = "blue", size = 0.80) +
 #'    geom_sf(data = NE_SA, fill = NA, color="black", size = 0.80) +
 #'    # geom_sf(data = bob_NE_MCP, fill = NA, color = "red", size = 0.8) +
 #'    geom_sf(data = bob_NE_poly, fill = NA, color = "blue", size = 0.80) +
 #'    geom_sf(data = bob_sf, color = "black", shape = 16) +
 #'    geom_sf(data = bob_skinny_NE, color = "green", shape = 16) +
 #'    geom_sf(data = bob_skinny_OK, color = "orange", shape = 16) +
 #'    #'  Add camera locations and vary color by deployment year
 #'    geom_sf(data = bobstart, color = "black", shape = 16) +
 #'    #'  Constrain plot to both study areas
 #'    coord_sf(xlim = c(490000.0, 780000.0), ylim = c(102000.0, 250000.0))
  
 
  
  ####  3. Plot OccMod & RSF Results  ####
  #'  ====================================
  #'  Plot effect size and CI's for occ mod & RSF results for key covariates
  #'  that were generally significant in both models (elevation, % forest, % grass,
  #'  human modified landscape). Display in multiple panels by model and season.
  
  #'  Occupancy model output
  occ_out <- read.csv("./Outputs/Tables/OccMod_OccProb_Results_2021-08-10.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE
    #'  Calculate 90% confidence intervals to mirror alpha = 0.1
    mutate(
      l95 = (Estimate - (1.64 * SE)),  #### REMINDER: this is 90% CI even though column says l95/u95
      u95 = (Estimate + (1.64 * SE))   
    ) %>%
    dplyr::select(-c(X, Model))
  #'  RSF results output
  rsf_out <- read.csv("./Outputs/Tables/RSF_Results_2021-08-10.csv") %>% # MAKE SURE IT'S MOST CURRENT DATE
    #'  Calculate 95% confidence intervals to mirror alpha = 0.05
    mutate(
      l95 = (Estimate - (1.96 * SE)),  #### REMINDER: this is 95% CI
      u95 = (Estimate + (1.96 * SE))
    ) %>%
    dplyr::select(-X)
  
  #'  Isolate results for specific covariates and seasons
  #'  Elevation
  elev_occ <- filter(occ_out, Parameter == "Elev")
  elev_occ_smr <- filter(elev_occ, Season == "Summer")
  elev_occ_wtr <- filter(elev_occ, Season == "Winter")
  elev_rsf <- filter(rsf_out, Parameter == "Elev")
  elev_rsf_smr <- filter(elev_rsf, Season == "Summer")
  elev_rsf_wtr <- filter(elev_rsf, Season == "Winter")
  #'  Slope
  slope_occ <- filter(occ_out, Parameter == "Slope")
  slope_occ_smr <- filter(slope_occ, Season == "Summer")
  slope_occ_wtr <- filter(slope_occ, Season == "Winter")
  slope_rsf <- filter(rsf_out, Parameter == "Slope")
  slope_rsf_smr <- filter(slope_rsf, Season == "Summer")
  slope_rsf_wtr <- filter(slope_rsf, Season == "Winter")
  #'  Percent Forest
  for_occ <- filter(occ_out, Parameter == "PercForMix")
  for_occ_smr <- filter(for_occ, Season == "Summer")
  for_occ_wtr <- filter(for_occ, Season == "Winter")
  for_rsf <- filter(rsf_out, Parameter == "PercForMix")
  for_rsf_smr <- filter(for_rsf, Season == "Summer")
  for_rsf_wtr <- filter(for_rsf, Season == "Winter")
  #'  Percent Grass
  grass_occ <- filter(occ_out, Parameter == "PercXGrass")
  #'  Add in species that did not test for effect of percent grass in OccMods
  elk_smr <- c("Elk", "Summer", "PercXGrass", NA, NA, NA, NA, NA, NA)
  elk_wtr <- c("Elk", "Winter", "PercXGrass", NA, NA, NA, NA, NA, NA)
  wtd_smr <- c("White-tailed Deer", "Summer", "PercXGrass", NA, NA, NA, NA, NA, NA)
  wtd_wtr <- c("White-tailed Deer", "Winter", "PercXGrass", NA, NA, NA, NA, NA, NA)
  grass_occ <- rbind(grass_occ, elk_smr, elk_wtr, wtd_smr, wtd_wtr) %>%
    mutate(Estimate = as.numeric(as.character(Estimate)),
           l95 = as.numeric(as.character(l95)),
           u95 = as.numeric(as.character(u95)))
  grass_occ_smr <- filter(grass_occ, Season == "Summer")
  grass_occ_wtr <- filter(grass_occ, Season == "Winter")
  grass_rsf <- filter(rsf_out, Parameter == "PercXGrass")
  grass_rsf_smr <- filter(grass_rsf, Season == "Summer")
  grass_rsf_wtr <- filter(grass_rsf, Season == "Winter")
  #'  Percent Shrub
  shrub_occ <- filter(occ_out, Parameter == "PercXShrub")
  #'  Add in species that did not test for effect of percent shrub in OccMods
  elk_smr <- c("Elk", "Summer", "PercXShrub", NA, NA, NA, NA, NA, NA)
  elk_wtr <- c("Elk", "Winter", "PercXShrub", NA, NA, NA, NA, NA, NA)
  wtd_smr <- c("White-tailed Deer", "Summer", "PercXShrub", NA, NA, NA, NA, NA, NA)
  wtd_wtr <- c("White-tailed Deer", "Winter", "PercXShrub", NA, NA, NA, NA, NA, NA)
  wolf_smr <- c("Wolf", "Summer", "PercXShrub", NA, NA, NA, NA, NA, NA)
  wolf_wtr <- c("Wolf", "Winter", "PercXShrub", NA, NA, NA, NA, NA, NA)
  shrub_occ <- rbind(shrub_occ, elk_smr, elk_wtr, wtd_smr, wtd_wtr, wolf_smr, wolf_wtr) %>%
    mutate(Estimate = as.numeric(as.character(Estimate)),
           l95 = as.numeric(as.character(l95)),
           u95 = as.numeric(as.character(u95)))
  shrub_occ_smr <- filter(shrub_occ, Season == "Summer")
  shrub_occ_wtr <- filter(shrub_occ, Season == "Winter")
  shrub_rsf <- filter(rsf_out, Parameter == "PercXShrub")
  shrub_rsf_smr <- filter(shrub_rsf, Season == "Summer")
  shrub_rsf_wtr <- filter(shrub_rsf, Season == "Winter")
  #'  Road Density
  rdden_occ <- filter(occ_out, Parameter == "RoadDensity")
  rdden_occ_smr <- filter(rdden_occ, Season == "Summer")
  rdden_occ_wtr <- filter(rdden_occ, Season == "Winter")
  rdden_rsf <- filter(rsf_out, Parameter == "RoadDen")
  rdden_rsf_smr <- filter(rdden_rsf, Season == "Summer")
  rdden_rsf_wtr <- filter(rdden_rsf, Season == "Winter")
  #'  Human Modified
  hm_occ <- filter(occ_out, Parameter == "HumanMod")
  hm_occ_smr <- filter(hm_occ, Season == "Summer")
  hm_occ_wtr <- filter(hm_occ, Season == "Winter")
  hm_rsf <- filter(rsf_out, Parameter == "HumanMod")
  hm_rsf_smr <- filter(hm_rsf, Season == "Summer")
  hm_rsf_wtr <- filter(hm_rsf, Season == "Winter")

  
  #'  Plot effects by covariate, season, and model type for all species
  #'  --------------------------------  
  ####  PLOT OCCUPANCY MODEL RESULTS  ####
  #'  --------------------------------
  #'  Effect of ELEVATION on probability of use (on logit scale)
  #'  Summer results
  elev_occ_smr_fig <- ggplot(elev_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    # scale_color_manual(name = "Species"), 
                        # labels = c("Above Average", "Below Average"), 
                        # values = c("above"="#00ba38", "below"="#f8766d")) + 
    # geom_text(color = "black", size = 2) +
    labs(title = "Summer Elevation", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-4, 2.5) + # (-5, 5)
    coord_flip() 
  #'  Winter results
  elev_occ_wtr_fig <- ggplot(elev_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Elevation", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-5, 2) +
    coord_flip()
  
  #'  Effect of SLOPE on probability of use (on logit scale)
  #'  Summer results
  slope_occ_smr_fig <- ggplot(slope_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Slope", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-1.5, 2) + # (-5, 5)
    coord_flip() 
  #'  Winter results
  slope_occ_wtr_fig <- ggplot(slope_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Slope", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 1.5) +
    coord_flip()
  
  #'  Effect of PERCENT MIXED FOREST on Probability of Use (logit scale)
  #'  Summer results
  for_occ_smr_fig <- ggplot(for_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Percent Forest", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2.5) +
    coord_flip()
  #'  Winter results
  for_occ_wtr_fig <- ggplot(for_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Percent Forest", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.5, 3.5) +
    coord_flip()
  
  #'  Effect of PERCENT GRASS on Probability of Use (logit scale)
  #'  Warnings are OK- for wtd & elk where % grass cov was excluded from OccMods
  #'  Summer results
  grass_occ_smr_fig <- ggplot(grass_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Percent Grass", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2.5) +
    coord_flip()
  #'  Winter results
  grass_occ_wtr_fig <- ggplot(grass_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Percent Grass", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.5, 12) +
    coord_flip()
  
  #'  Effect of PERCENT SHRUB on Probability of Use (logit scale)
  #'  Warnings are OK- for wtd & elk where % shrub cov was excluded from OccMods
  #'  Summer results
  shrub_occ_smr_fig <- ggplot(shrub_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Percent Shrub", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-1.5, 1.5) +
    coord_flip()
  #'  Winter results
  shrub_occ_wtr_fig <- ggplot(shrub_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Percent Shrub", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1.5, 2) +
    coord_flip()
  
  #'  Effect of ROAD DENSITY on Probability of Use (logit scale)
  #'  Summer results
  rdden_occ_smr_fig <- ggplot(rdden_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Road Density", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2, 1.5) +
    coord_flip()
  #'  Winter results
  rdden_occ_wtr_fig <- ggplot(rdden_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Road Density", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.5, 2) +
    coord_flip()
  
  #'  Effect of PERCENT HUMAN MODIFIED LANDSCAPE on Probability of Use (logit scale)
  #'  Summer results
  hm_occ_smr_fig <- ggplot(hm_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Human Modified", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2) +
    coord_flip()
  #'  Winter results
  hm_occ_wtr_fig <- ggplot(hm_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Winter Human Modified", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.5, 2) +
    coord_flip()
  
  #'  -----------------------------------  
  ####  PLOT RESOURCE SELECTION RESULTS  ####
  #'  -----------------------------------
  #'  Effect of ELEVATION on relative probability of selection (on logit scale)
  #'  Summer results
  elev_rsf_smr_fig <- ggplot(elev_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Elevation", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.75, 0.75) +
    coord_flip()
  #'  Winter results
  elev_rsf_wtr_fig <- ggplot(elev_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", # "Elevation", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1.5, 0.5) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.3, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.3, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.3, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.3, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.3, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.3, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.3, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of SLOPE on relative probability of selection (on logit scale)
  #'  Summer results
  slope_rsf_smr_fig <- ggplot(slope_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Slope", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 0.5) +
    coord_flip()
  #'  Winter results
  slope_rsf_wtr_fig <- ggplot(slope_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", # "Slope", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.25, 0.8) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.7, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.7, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.7, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.7, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.7, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.7, ysize = 0.25, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.7, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of PERCENT MIXED FOREST on relative probability of selection (logit scale)
  #'  Summer results
  for_rsf_smr_fig <- ggplot(for_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Forest within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.5, 0.5) +
    coord_flip()
  #'  Winter results
  for_rsf_wtr_fig <- ggplot(for_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Forest within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 1.25) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 1, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 1, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 1, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 1, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 1, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 1, ysize = 0.48, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 1, ysize = 0.4, color = "black", alpha = 1)
  # add_phylopic(wolfimg, x = 7.05, y = 1.5, ysize = 2.1, color = "black", alpha = 1) +
  #   add_phylopic(wtdimg, x = 6.1, y = 1.5, ysize = 1.7, color = "black", alpha = 1) +
  #   add_phylopic(mdimg, x = 5.05, y = 1.5, ysize = 2.3, color = "black", alpha = 1) +
  #   add_phylopic(elkmimg, x = 4.05, y = 1.5, ysize = 1.85, color = "black", alpha = 1) +
  #   add_phylopic(coyimg, x = 3.05, y = 1.5, ysize = 1.6, color = "black", alpha = 1) +
  #   add_phylopic(cougimg, x = 2, y = 1.5, ysize = 2.4, color = "black", alpha = 1) +
  #   add_phylopic(bobimg, x = 1.05, y = 1.5, ysize = 1.85, color = "black", alpha = 1)
  
  #'  Effect of PERCENT GRASS on relative probability of selection (logit scale)
  #'  Summer results
  grass_rsf_smr_fig <- ggplot(grass_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Grass within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1.5, 1.5) +
    coord_flip()
  #'  Winter results
  grass_rsf_wtr_fig <- ggplot(grass_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Grass within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 1) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.8, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.8, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.8, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.8, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.8, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.8, ysize = 0.45, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.8, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of PERCENT SHRUB on relative probability of selection (logit scale)
  #'  Summer results
  shrub_rsf_smr_fig <- ggplot(shrub_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Grass within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-3, 0.5) +
    coord_flip()
  #'  Winter results
  shrub_rsf_wtr_fig <- ggplot(shrub_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Grass within 250m", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.5, 0.6) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.5, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.5, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.5, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.5, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.5, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.5, ysize = 0.25, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.5, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of ROAD DENSITY on relative probability of selection (logit scale)
  #'  Summer results
  rdden_rsf_smr_fig <- ggplot(rdden_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Road Density", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.5, 1) +
    coord_flip()
  #'  Winter results
  rdden_rsf_wtr_fig <- ggplot(rdden_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Road Density", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-0.5, 0.8) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.7, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.7, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.7, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.7, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.7, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.7, ysize = 0.3, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.7, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of PERCENT HUMAN MODIFIED LANDSCAPE on relative probability of selection (logit scale)
  #'  Summer results
  hm_rsf_smr_fig <- ggplot(hm_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Human Modified Landscape", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.0, 1) +
    coord_flip()
  #'  Winter results
  hm_rsf_wtr_fig <- ggplot(hm_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Human Modified Landscape", 
         subtitle = "Resource Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 0.5) +
    coord_flip() +
    add_phylopic(wolfimg, x = 7.05, y = 0.35, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(wtdimg, x = 6.1, y = 0.35, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.35, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.35, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.35, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.35, ysize = 0.38, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.35, ysize = 0.4, color = "black", alpha = 1)
  
  
  ####  Pair OccMod and RSF plots  ####
  #'  patchwork figures together:
  elev_fig <- elev_occ_smr_fig + plot_annotation(title = "A") + elev_rsf_smr_fig + elev_occ_wtr_fig + elev_rsf_wtr_fig + plot_layout(ncol = 4)
  slope_fig <- slope_occ_smr_fig + plot_annotation(title = "B") + slope_rsf_smr_fig + slope_occ_wtr_fig + slope_rsf_wtr_fig + plot_layout(ncol = 4)
  for_fig <- for_occ_smr_fig + plot_annotation(title = "C") + for_rsf_smr_fig + for_occ_wtr_fig + for_rsf_wtr_fig + plot_layout(ncol = 4)
  grass_fig <- grass_occ_smr_fig + plot_annotation(title = "D") + grass_rsf_smr_fig + grass_occ_wtr_fig + grass_rsf_wtr_fig + plot_layout(ncol = 4)
  shrub_fig <- shrub_occ_smr_fig + plot_annotation(title = "E") + shrub_rsf_smr_fig + shrub_occ_wtr_fig + shrub_rsf_wtr_fig + plot_layout(ncol = 4)
  rdden_fig <- rdden_occ_smr_fig + plot_annotation(title = "F") + rdden_rsf_smr_fig + rdden_occ_wtr_fig + rdden_rsf_wtr_fig + plot_layout(ncol = 4)
  hm_fig <- hm_occ_smr_fig + plot_annotation(title = "G") + hm_rsf_smr_fig + hm_occ_wtr_fig + hm_rsf_wtr_fig + plot_layout(ncol = 4)
  
    
  plot(elev_fig)
  plot(slope_fig)
  plot(for_fig)
  plot(grass_fig)
  plot(shrub_fig)
  plot(rdden_fig)
  plot(hm_fig)
  
  #'  Save different file formats
  #'  PNG
  ggsave("./Outputs/Figures/Side by Side/Elevation_Occ-RSF_plot.png", elev_fig) #, width = 15, units = "in"
  ggsave("./Outputs/Figures/Side by Side/Slope_Occ-RSF_plot.png", slope_fig)
  ggsave("./Outputs/Figures/Side by Side/Forest_Occ-RSF_plot.png", for_fig)
  ggsave("./Outputs/Figures/Side by Side/Grass_Occ-RSF_plot.png", grass_fig)
  ggsave("./Outputs/Figures/Side by Side/Shrub_Occ-RSF_plot.png", shrub_fig)
  ggsave("./Outputs/Figures/Side by Side/RoadDensity_Occ-RSF_plot.png", rdden_fig)
  ggsave("./Outputs/Figures/Side by Side/HumanMod_Occ-RSF_plot.png", hm_fig)
  #'  JPEG
  ggsave("./Outputs/Figures/Side by Side/Elevation_Occ-RSF_plot.jpeg", elev_fig)
  ggsave("./Outputs/Figures/Side by Side/Slope_Occ-RSF_plot.jpeg", slope_fig)
  ggsave("./Outputs/Figures/Side by Side/Forest_Occ-RSF_plot.jpeg", for_fig)
  ggsave("./Outputs/Figures/Side by Side/Grass_Occ-RSF_plot.jpeg", grass_fig)
  ggsave("./Outputs/Figures/Side by Side/Shrub_Occ-RSF_plot.jpeg", shrub_fig)
  ggsave("./Outputs/Figures/Side by Side/RoadDensity_Occ-RSF_plot.jpeg", rdden_fig)
  ggsave("./Outputs/Figures/Side by Side/HumanMod_Occ-RSF_plot.jpeg", hm_fig)
  
  #'  One single figure... ugly and terrible and don't do it
  # cov_fig <- elev_fig / slope_fig / for_fig / grass_fig / shrub_fig / rdden_fig / hm_fig
  # ggsave("./Outputs/Figures/OccMod_v_RSF_plot.png", cov_fig)
  # cov_fig_partial <- elev_fig / slope_fig / for_fig / hm_fig
  # ggsave("./Outputs/Figures/OccMod_v_RSF_plot_partial.png", cov_fig_partial)



  ####  4. Plot OccMod and RSF Results Against Each Other  ####
  #'  ====================================================
  #'  Plot the correlation between OccMod and RSF confidence intervals for easier 
  #'  comparison of which effects are significant and consistent across the
  #'  two different models.
  
  #'  Pull out confidence interval data by covariate and species
  occmod_90ci <- dplyr::select(occ_out, c(Species, Season, Parameter, Estimate, l95, u95)) %>%
    mutate(
      Estimate_occ = Estimate,
      l95_occ = l95,
      u95_occ = u95,
      Parameter = as.character(Parameter)
    ) %>%
    dplyr::select(-c(Estimate, l95, u95)) %>%
    filter(Parameter != "AreaOK") %>%
    filter(Parameter != "(Intercept)")
  rsf_95ci <- dplyr::select(rsf_out, c(Species, Season, Parameter, Estimate, l95, u95)) %>%
    mutate(
      Estimate_rsf = Estimate,
      l95_rsf = l95, #- 0.05, # SUBTRACTING 1 FROM L95 SO IT'S VISIBLE IN PLOT
      u95_rsf = u95, #+ 0.05, # ADDING 1 TO U95 SO IT'S VISIBLE IN PLOT
      Parameter = as.character(Parameter),
      Parameter = ifelse(Parameter == "RoadDen", "RoadDensity", Parameter)
    ) %>%
    dplyr::select(-c(Estimate, l95, u95)) %>%
    filter(Parameter != "(Intercept)")
  #'  Merge CI data from different models into single table
  combo_ci <- full_join(occmod_90ci, rsf_95ci, by = c("Species", "Season", "Parameter"))
  
  #'  Identify min and max values of confidence intervals
  min(combo_ci$l95_occ, na.rm = TRUE); max(combo_ci$u95_occ, na.rm = TRUE) # -3.89 to 10.877
  min(combo_ci$l95_rsf, na.rm = TRUE); max(combo_ci$u95_rsf, na.rm = TRUE) # -2.91 to 0.63
  
  #'  Separate CIs by Covariate
  elev_ci <- combo_ci[combo_ci$Parameter == "Elev",]
  slope_ci <- combo_ci[combo_ci$Parameter == "Slope",]
  for_ci <- combo_ci[combo_ci$Parameter == "PercForMix",]
  grass_ci <- combo_ci[combo_ci$Parameter == "PercXGrass",]
  shrub_ci <- combo_ci[combo_ci$Parameter == "PercXShrub",]
  rdden_ci <- combo_ci[combo_ci$Parameter == "RoadDensity",]
  hm_ci <- combo_ci[combo_ci$Parameter == "HumanMod",]
  #'  Separate CIs by Species
  md_ci <- combo_ci[combo_ci$Species == "Mule Deer",]
  elk_ci <- combo_ci[combo_ci$Species == "Elk",]
  wtd_ci <- combo_ci[combo_ci$Species == "White-tailed Deer",]
  coug_ci <- combo_ci[combo_ci$Species == "Cougar",]
  wolf_ci <- combo_ci[combo_ci$Species == "Wolf",]
  bob_ci <- combo_ci[combo_ci$Species == "Bobcat",]
  coy_ci <- combo_ci[combo_ci$Species == "Coyote",]
  
  #'  ---------------------------
  ####  By covariate and season  ####
  #'  ---------------------------
  #'  Plots for each covariate, season differs by point shape
  elev_ci_fig <- ggplot(elev_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Species)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Elevation") +
    xlab("RSF") + ylab("Occupancy") 
  slope_ci_fig <- ggplot(slope_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) +  
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Slope") +
    xlab("RSF") + ylab("Occupancy")
  for_ci_fig <- ggplot(for_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Percent Forest") +
    xlab("RSF") + ylab("Occupancy")
  grass_ci_fig <- ggplot(grass_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) +
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Percent Grass") +
    xlab("RSF") + ylab("Occupancy")
  shrub_ci_fig <- ggplot(shrub_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) +
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Percent Shrub") +
    xlab("RSF") + ylab("Occupancy")
  rdden_ci_fig <- ggplot(rdden_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Road Density") +
    xlab("RSF") + ylab("Occupancy")
  hm_ci_fig <- ggplot(hm_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Effect of Human Modified Landscape") +
    xlab("RSF") + ylab("Occupancy")
  
  #'  Save as PNG images
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Elevation_Occ-by-RSF_plot.png", elev_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Slope_Occ-by-RSF_plot.png", slope_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/PercentForest_Occ-by-RSF_plot.png", for_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/PercentGrass_Occ-by-RSF_plot.png", grass_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/PercentShrub_Occ-by-RSF_plot.png", shrub_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/RoadDensity_Occ-by-RSF_plot.png", rdden_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/HumanMod_Occ-by-RSF_plot.png", hm_ci_fig, width = 9.3, units = "in")
  
  #'  -------------------------
  ####  By Species and Season  ####
  #'  -------------------------
  #'  Plots for each species, season differs by point shape
  md_ci_fig <- ggplot(md_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Mule Deer") +
    xlab("RSF") + ylab("Occupancy") 
  elk_ci_fig <- ggplot(elk_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Elk") +
    xlab("RSF") + ylab("Occupancy")
  #'  NOT PLOTTING SHRUB AND GRASS RESULTS FOR RSFs
  wtd_ci_fig <- ggplot(wtd_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "White-tailed Deer") +
    xlab("RSF") + ylab("Occupancy")
  #'  NOT PLOTTING SHRUB AND GRASS RESULTS FOR RSFs
  coug_ci_fig <- ggplot(coug_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Cougar") +
    xlab("RSF") + ylab("Occupancy")
  wolf_ci_fig <- ggplot(wolf_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Wolf") +
    xlab("RSF") + ylab("Occupancy")   # NOT PLOTTING SOMETHING
  bob_ci_fig <- ggplot(bob_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Bobcat") +
    xlab("RSF") + ylab("Occupancy") 
  coy_ci_fig <- ggplot(coy_ci, aes(x = Estimate_rsf, y = Estimate_occ, col = Parameter)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, alpha = 0.5) +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.01) + 
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) + 
    geom_point(stat = 'identity', aes(shape = Season), size = 3.5) + 
    scale_shape_manual(values = c(19, 23)) +
    labs(title = "Coyote") +
    xlab("RSF") + ylab("Occupancy") 
  
  #'  Panel figure of all species
  corr_plot <- bob_ci_fig + coug_ci_fig + coy_ci_fig + elk_ci_fig + 
    md_ci_fig + wtd_ci_fig + wolf_ci_fig + guide_area() + 
    plot_layout(guides = 'collect') + plot_layout(ncol = 2) + 
    plot_annotation(title = "Correlation between Occupancy Model and RSF covariate effects") + 
    plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 10)) +
    theme(legend.box = 'horizontal')
  
  plot(corr_plot)
  
  #'  Save as PNG images
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/MuleDeer_Occ-by-RSF_plot.png", md_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Elk_Occ-by-RSF_plot.png", elk_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/WTD_Occ-by-RSF_plot.png", wtd_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Cougar_Occ-by-RSF_plot.png", coug_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Wolf_Occ-by-RSF_plot.png", wolf_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Bobcat_Occ-by-RSF_plot.png", bob_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Coyote_Occ-by-RSF_plot.png", coy_ci_fig, width = 9.3, units = "in")
  ggsave("./Outputs/Figures/Occu-RSF-Correlation/Occ-by-RSF_plot.png", corr_plot, width = 8, height = 15, units = "in")
  
  #'  ------------------------------
  ####  Stand alone seasonal plots  ####
  #'  ------------------------------
  #'  Plots for each season and covariate
  elev_smr_ci_fig <- ggplot(elev_ci[elev_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +  # "Occupancy Model vs. RSF Covariate Effects
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  elev_wtr_ci_fig <- ggplot(elev_ci[elev_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy") 
  slope_smr_ci_fig <- ggplot(slope_ci[slope_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  slope_wtr_ci_fig <- ggplot(slope_ci[slope_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  for_smr_ci_fig <- ggplot(for_ci[for_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.025) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  for_wtr_ci_fig <- ggplot(for_ci[for_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  grass_smr_ci_fig <- ggplot(grass_ci[grass_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  grass_wtr_ci_fig <- ggplot(grass_ci[grass_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  shrub_smr_ci_fig <- ggplot(shrub_ci[shrub_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none") 
  shrub_wtr_ci_fig <- ggplot(shrub_ci[shrub_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  rdden_smr_ci_fig <- ggplot(rdden_ci[rdden_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  rdden_wtr_ci_fig <- ggplot(rdden_ci[rdden_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  hm_smr_ci_fig <- ggplot(hm_ci[hm_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Summer Responses") +
    xlab("RSF") + ylab("Occupancy") +
    theme(legend.position = "none")
  hm_wtr_ci_fig <- ggplot(hm_ci[hm_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species)) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  
  hm_ci_fig <- ggplot(hm_ci, aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Species), width = 0.01) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Species)) +
    geom_point(stat = 'identity', aes(col = Species, shape = Season), size = 2.5) + #size = 3.5 
    labs(title = "Winter Responses") +
    xlab("RSF") + ylab("Occupancy")
  
  #'  Plot summer and winter figures together
  elev_ci_fig <- elev_smr_ci_fig + plot_annotation(title = "A. Elevation") + elev_wtr_ci_fig
  slope_ci_fig <- slope_smr_ci_fig + plot_annotation(title = "B. Slope") + slope_wtr_ci_fig
  for_ci_fig <- for_smr_ci_fig + plot_annotation(title = "C. Percent Forest") + for_wtr_ci_fig
  grass_ci_fig <- grass_smr_ci_fig + plot_annotation(title = "D. Percent Grass") + grass_wtr_ci_fig
  shrub_ci_fig <- shrub_smr_ci_fig + plot_annotation(title = "E. Percent Shrub") + shrub_wtr_ci_fig
  rdden_ci_fig <- rdden_smr_ci_fig + plot_annotation(title = "F. Road Density") + rdden_wtr_ci_fig
  hm_ci_fig <- hm_smr_ci_fig + plot_annotation(title = "G. Human Modification") + hm_wtr_ci_fig
  
  #'  Visualize these lovely plots!
  plot(elev_smr_ci_fig)
  plot(elev_wtr_ci_fig)
  plot(slope_smr_ci_fig)
  plot(slope_wtr_ci_fig)
  plot(for_smr_ci_fig)
  plot(for_wtr_ci_fig)
  plot(grass_smr_ci_fig)
  plot(grass_wtr_ci_fig)
  plot(shrub_smr_ci_fig)
  plot(shrub_wtr_ci_fig)
  plot(rdden_smr_ci_fig)
  plot(rdden_wtr_ci_fig)
  plot(hm_smr_ci_fig)
  plot(hm_wtr_ci_fig)
  
  plot(elev_ci_fig)
  plot(slope_ci_fig)
  plot(for_ci_fig)
  plot(grass_ci_fig)
  plot(shrub_ci_fig)
  plot(rdden_ci_fig)
  plot(hm_ci_fig)
  
  #' #'  Save as PNG images
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/Elevation_Occ-by-RSF_plot.png", elev_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/Slope_Occ-by-RSF_plot.png", slope_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/Forest_Occ-by-RSF_plot.png", for_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/Grass_Occ-by-RSF_plot.png", grass_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/Shrub_Occ-by-RSF_plot.png", shrub_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/RoadDensity_Occ-by-RSF_plot.png", rdden_ci_fig)
  #' ggsave("./Outputs/Figures/Occu-RSF-Correlation/HumanMod_Occ-by-RSF_plot.png", hm_ci_fig)
  
  
  #'  By species and season
  md_smr_ci_fig <- ggplot(md_ci[md_ci$Season == "Summer",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.025) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) +
    geom_point(stat = 'identity', aes(col = Parameter)) + #size = 3.5 
    labs(title = "Occupancy Model vs. RSF Covariate Effects",
         subtitle = "Mule Deer Summer") +
    xlab("RSF") + ylab("Occupancy") 
  # need to figure out how to make winter and summer estimates different shades so they can be plotted together
  md_wtr_ci_fig <- ggplot(md_ci[md_ci$Season == "Winter",], aes(x = Estimate_rsf, y = Estimate_occ)) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_abline(slope = 1, intercept = 0, col = "darkgray") +
    geom_errorbar(aes(ymin = l95_occ, ymax = u95_occ, col = Parameter), width = 0.025) +
    geom_errorbarh(aes(xmin = l95_rsf, xmax = u95_rsf, colour = Parameter)) +
    geom_point(stat = 'identity', aes(col = Parameter)) + #size = 3.5 
    labs(title = "Occupancy Model vs. RSF Covariate Effects",
         subtitle = "Mule Deer Winter") +
    xlab("RSF") + ylab("Occupancy")
  
  
  
  ####  5. Histograms of Covariate Values by Data Source  ####
  #'  ==================================================
  #'  Plot the frequency of covariate values for each data source to see if that
  #'  explains some of the major discrepancies between wildlife-habitat relationships
  #'  estimated by the OccMods and RSFs.
  
  #'  Occupancy input data
  source("./Scripts/CameraTrap_DetectionHistories.R")
  #'  Covariate data that went into OccMods, just not scaled here
  stations <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/CameraLocation_Covariates18-20_2021-08-10.csv") %>% #2021-05-14
    mutate(
      Study_Area = ifelse(Study_Area == "NE ", "NE", as.character(Study_Area)),
    ) %>%
    mutate(
      Monitoring = ifelse(Monitoring == "Closed road", "Dirt road", as.character(Monitoring)),
      Monitoring = ifelse(Monitoring == "Game trail", "Trail", as.character(Monitoring)),
    ) %>%
    mutate(
      Habitat_Type = ifelse(Habitat_Type == "Agriculture", "Grassland", as.character(Habitat_Type)), 
      Habitat_Type = ifelse(Habitat_Type == "Riparian", "Mixed conifer", as.character(Habitat_Type)) 
    ) %>%
    transmute(
      Year = as.factor(Year),
      Study_Area = as.factor(Study_Area),
      CameraLocation = as.factor(CameraLocation),
      Distance = as.numeric(Distance_Focal_Point),
      Height = as.numeric(Height_frm_grnd),
      Trail = as.factor(Monitoring),
      PercForestMix = as.numeric(PercForestMix2),    
      PercXericShrub = as.numeric(PercXericShrub),
      PercXericGrass = as.numeric(PercXericGrass),
      Elev = as.numeric(Elev),
      Slope = as.numeric(Slope),
      RoadDensity = as.numeric(RoadDen), 
      HumanModified = as.numeric(HumanMod)
      # Elev = as.numeric(elevation),
      # Slope = as.numeric(slope),
      # RoadDensity = as.numeric(road_density), 
      # HumanModified = as.numeric(modified)
    )
 
  #'  Format and combine detection histories and covariate data by camera station
  #'  Summer detection histories
  smr_det <- function(smr_DH) {
    DH <- as.data.frame(smr_DH)
    DH$CameraLocation <- row.names(DH)
    DH <- DH %>%
      #'  Spread detection data so each row is one sampling occasion per site
      pivot_longer(!CameraLocation, names_to = "occasion", values_to = "ndet") %>%
      #'  Remove NAs, occasions with non-detections, & DH columns
      filter(!is.na(ndet)) %>%
      filter(ndet == 1) %>%
      dplyr::select(-c(ndet, occasion)) %>%
      #'  Add season & site-specific covariate data
      mutate(Season = "Summer") %>%
      left_join(stations, by = "CameraLocation")
    return(DH)
  }
  #'  Run summer detection histories through function
  summer_DHs <- list(DH_md_smr1819, DH_elk_smr1819, DH_wtd_smr1819, DH_coug_smr1819, 
                     DH_wolf_smr1819, DH_bob_smr1819, DH_coy_smr1819)
  summer_det <- lapply(summer_DHs, smr_det)
  #'  Winter detection histories
  wtr_det <- function(wtr_DH) {
    DH <- as.data.frame(wtr_DH)
    DH$CameraLocation <- row.names(DH)
    DH <- DH %>%
      #'  Spread detection data so each row is one sampling occasion per site
      pivot_longer(!CameraLocation, names_to = "occasion", values_to = "ndet") %>%
      #'  Remove NAs, occasions with non-detections, & DH columns
      filter(!is.na(ndet)) %>%
      filter(ndet == 1) %>%
      dplyr::select(-c(ndet, occasion)) %>%
      #'  Add season & site-specific covariate data
      mutate(Season = "Winter") %>%
      left_join(stations, by = "CameraLocation")
    return(DH)
  }
  #'  Run winter detection histories through function
  winter_DHs <- list(DH_md_wtr1820, DH_elk_wtr1820, DH_wtd_wtr1820, DH_coug_wtr1820, 
                     DH_wolf_wtr1820, DH_bob_wtr1820, DH_coy_wtr1820)
  winter_det <- lapply(winter_DHs, wtr_det)
  
  #'  Combine seasonal data by species
  md_det <- rbind(summer_det[[1]], winter_det[[1]]) 
  elk_det <- rbind(summer_det[[2]], winter_det[[2]]) 
  wtd_det <- rbind(summer_det[[3]], winter_det[[3]]) 
  coug_det <- rbind(summer_det[[4]], winter_det[[4]]) 
  wolf_det <- rbind(summer_det[[5]], winter_det[[5]]) 
  bob_det <- rbind(summer_det[[6]], winter_det[[6]]) 
  coy_det <- rbind(summer_det[[7]], winter_det[[7]]) 
  
  
  #'  RSF input data
  load("./Outputs/RSF_pts/md_dat_2nd_all_2021-08-10.RData") #2021-07-22
  load("./Outputs/RSF_pts/elk_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/wtd_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/coug_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/wolf_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/bob_dat_2nd_all_2021-08-10.RData")
  load("./Outputs/RSF_pts/coy_dat_2nd_all_2021-08-10.RData")
  
  #'  Retain only the used locations and their covariate values
  md_used <- filter(md_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  elk_used <- filter(elk_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  wtd_used <- filter(wtd_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  coug_used <- filter(coug_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  wolf_used <- filter(wolf_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  bob_used <- filter(bob_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  coy_used <- filter(coy_dat_all, Used == 1) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  #'  Retain only the available locations and their covariate values
  md_avail <- filter(md_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  elk_avail <- filter(elk_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  wtd_avail <- filter(wtd_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  coug_avail <- filter(coug_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  wolf_avail <- filter(wolf_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  bob_avail <- filter(bob_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  coy_avail <- filter(coy_dat_all, Used == 0) %>%
    mutate(Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"))
  
  
  #'  Combine detection and GPS location data
  all_obs <- function(det, used) {
    cam_data <- transmute(det, Data = "Camera", 
                          Season = Season, Year = Year, Area = Study_Area, 
                          Elev = Elev, Slope = Slope, PercForMix = PercForestMix, 
                          PercXGrass = PercXericGrass, PercXShrub = PercXericShrub, 
                          RoadDen = RoadDensity, HumanMod = HumanModified)
    col_data <- mutate(used, Data = "Collar") %>%
      dplyr::select(Data, Season, Year, Area, Elev, Slope, PercForMix, PercXGrass, 
                    PercXShrub, RoadDen, HumanMod)
    combo_data <- rbind(cam_data, col_data)
    return(combo_data)
  }
  #'  Combine camera data with used data
  md_obs <- all_obs(md_det, md_used)
  elk_obs <- all_obs(elk_det, elk_used)
  wtd_obs <- all_obs(wtd_det, wtd_used)
  coug_obs <- all_obs(coug_det, coug_used)
  wolf_obs <- all_obs(wolf_det, wolf_used)
  bob_obs <- all_obs(bob_det, bob_used)
  coy_obs <- all_obs(coy_det, coy_used)
  #'  Combine camera data with available data
  md_sampled <- all_obs(md_det, md_avail)
  elk_sampled <- all_obs(elk_det, elk_avail)
  wtd_sampled <- all_obs(wtd_det, wtd_avail)
  coug_sampled <- all_obs(coug_det, coug_avail)
  wolf_sampled <- all_obs(wolf_det, wolf_avail)
  bob_sampled <- all_obs(bob_det, bob_avail)
  coy_sampled <- all_obs(coy_det, coy_avail)
  
  
  #'  Identify mean and median values for each variable by data set
  #'  List observation data
  all_data <- list(md_obs, elk_obs, wtd_obs, coug_obs, wolf_obs, bob_obs, coy_obs)
  avail_data <- list(md_sampled, elk_sampled, wtd_sampled, coug_sampled, wolf_sampled, bob_sampled, coy_sampled)
  #'  Calculate mean values per data type and species 
  mean_obs <- function(obs) {
    mu <- obs %>% 
      group_by(Data, Season) %>% 
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
      ungroup()
    return(mu)
  }
  mu_obs <- lapply(all_data, mean_obs) # ignore summarise() comment about grouping
  mu_avail <- lapply(avail_data, mean_obs)
  #'  Calculate median values per data type and species
  median_obs <- function(obs) {
    med <- obs %>% 
      group_by(Data, Season) %>% 
      summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE))) %>%
      ungroup()
    return(med)
  }
  med_obs <- lapply(all_data, median_obs)
  med_avail <- lapply(avail_data, median_obs)
  
  
  #'  -----------------------------------
  ####  Histograms of Covariate Values  ####
  #'  -----------------------------------
  #'  Camera detections by season
  ggplot(md_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Mule Deer Detections", x = "Elevation (m)", y = "Number of Detections") # technically its the number of sampling occasions with a detection across all cameras
  ggplot(elk_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity")  + labs(title = "Elk Detections", x = "Elevation (m)", y = "Number of Detections")
  ggplot(wtd_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity")  + labs(title = "White-tailed Deer Detections", x = "Elevation (m)", y = "Number of Detections")
  ggplot(coug_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Cougar Detections", x = "Elevation (m)", y = "Number of Detections")
  ggplot(wolf_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Wolf Detections", x = "Elevation (m)", y = "Number of Detections")
  ggplot(bob_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Bobcat Detections", x = "Elevation (m)", y = "Number of Detections")
  ggplot(coy_det, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Coyote Detections", x = "Elevation (m)", y = "Number of Detections")
  #'  GPS collar locations by season
  ggplot(md_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Mule Deer GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(elk_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Elk GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(wtd_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "White-tailed Deer GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(coug_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Cougar GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(wolf_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Wolf GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(bob_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Bobcat GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  ggplot(coy_used, aes(x = Elev, color = Season, fill = Season)) + geom_histogram(binwidth = 100, alpha = 0.5, position = "identity") + labs(title = "Coyote GPS Locations", x = "Elevation (m)", y = "Number of GPS Locations")
  
  #'  Compare used vs available for mule deer human modified landscape data
  md_useavail <- mutate(md_dat_all, 
                        Season = ifelse(Season == "Summer18" | Season == "Summer19", "Summer", "Winter"),
                        Used = as.factor(Used))
  ggplot(md_useavail, aes(x = HumanMod, color = Used, fill = Used)) + geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity") + labs(title = "Mule Deer GPS Locations", x = "Percent Human Modified Landscape", y = "Number of GPS Locations")
  
  
  #'  Compare camera detections and collar locations
  #'  ----------------------------------------------
  #'  NOTE: There are *SO* many more collar observations than camera detections
  #'  that the camera data aren't even visible in some of these plots. So using
  #'  the mapping = aes(y = stat(ncount)) argument to re-scale the collar data
  #'  so they are on the same scale as the camera data and proportional to the
  #'  true frequency of collar location data. I think that's what's happening here....
  #'  Could also use y = stat(density) but I'm less clear what that does.
  
  #'  ELEVATION histograms
  #'  Plotting the MEAN values for each data source as dashed lines even though 
  #'  median addresses major skews in some of these data sets. Model intercepts
  #'  are based on the mean to mean values better represent model outputs.
  #'  Summer
  md_elev_smr <- ggplot(md_obs[md_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[1]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$Elev[3], linetype = "dashed", color = "blue")
  elk_elev_smr <- ggplot(elk_obs[elk_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 1700)) +
    geom_vline(xintercept = mu_obs[[2]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$Elev[3], linetype = "dashed", color = "blue")
  wtd_elev_smr <- ggplot(wtd_obs[wtd_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "White-tailed Deer Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 1700)) +
    geom_vline(xintercept = mu_obs[[3]]$Elev[1], linetype = "dashed", color = "red") +  
    geom_vline(xintercept = mu_obs[[3]]$Elev[3], linetype = "dashed", color = "blue")
  coug_elev_smr <- ggplot(coug_obs[coug_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Cougar Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) +
    geom_vline(xintercept = mu_obs[[4]]$Elev[1], linetype = "dashed", color = "red") + #cam & collar means are identical!
    geom_vline(xintercept = mu_obs[[4]]$Elev[3], linetype = "dashed", color = "blue")
  wolf_elev_smr <- ggplot(wolf_obs[wolf_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Wolf Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) +
    geom_vline(xintercept = mu_obs[[5]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[5]]$Elev[3], linetype = "dashed", color = "blue")
  bob_elev_smr <- ggplot(bob_obs[bob_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Bobcat Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) +
    geom_vline(xintercept = mu_obs[[6]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[6]]$Elev[3], linetype = "dashed", color = "blue")
  coy_elev_smr <- ggplot(coy_obs[coy_obs$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Coyote Summer Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) +
    geom_vline(xintercept = mu_obs[[7]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[7]]$Elev[3], linetype = "dashed", color = "blue")
  #'  Winter
  md_elev_wtr <- ggplot(md_obs[md_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[1]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$Elev[4], linetype = "dashed", color = "blue")
  elk_elev_wtr <- ggplot(elk_obs[elk_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 1700)) + 
    geom_vline(xintercept = mu_obs[[2]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$Elev[4], linetype = "dashed", color = "blue")
  wtd_elev_wtr <- ggplot(wtd_obs[wtd_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "White-tailed Deer Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 1700)) + 
    geom_vline(xintercept = mu_obs[[3]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[3]]$Elev[4], linetype = "dashed", color = "blue")
  coug_elev_wtr <- ggplot(coug_obs[coug_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Cougar Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[4]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[4]]$Elev[4], linetype = "dashed", color = "blue")
  wolf_elev_wtr <- ggplot(wolf_obs[wolf_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Wolf Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[5]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[5]]$Elev[4], linetype = "dashed", color = "blue")
  bob_elev_wtr <- ggplot(bob_obs[bob_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Bobcat Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[6]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[6]]$Elev[4], linetype = "dashed", color = "blue")
  coy_elev_wtr <- ggplot(coy_obs[coy_obs$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Coyote Winter Observations", x = "Elevation (m)", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_obs[[7]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[7]]$Elev[4], linetype = "dashed", color = "blue")
  #'  Elk camera samples vs what was available for elevation
  elk_allelev_smr <- ggplot(elk_sampled[elk_sampled$Season == "Summer",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Summer Cameras vs Availble Data", x = "Elevation (m)", y = "Proportion of Observations per Data Type")  + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_avail[[2]]$Elev[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[2]]$Elev[3], linetype = "dashed", color = "blue")
  elk_allelev_wtr <- ggplot(elk_sampled[elk_sampled$Season == "Winter",], aes(x = Elev, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 100, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Winter Cameras vs Availble Data", x = "Elevation (m)", y = "Proportion of Observations per Data Type")  + 
    coord_cartesian(xlim = c(300, 2500)) + 
    geom_vline(xintercept = mu_avail[[2]]$Elev[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[2]]$Elev[4], linetype = "dashed", color = "blue")
  
  #'  PERCENT FOREST histograms
  elk_for_smr <- ggplot(elk_obs[elk_obs$Season == "Summer",], aes(x = PercForMix, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Summer Observations", x = "Percent Forest within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[2]]$PercForMix[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$PercForMix[3], linetype = "dashed", color = "blue")
  elk_for_wtr <- ggplot(elk_obs[elk_obs$Season == "Winter",], aes(x = PercForMix, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Winter Observations", x = "Percent Forest within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[2]]$PercForMix[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$PercForMix[4], linetype = "dashed", color = "blue")
  #'  Camera samples vs what was available
  elk_allfor_smr <- ggplot(elk_sampled[elk_sampled$Season == "Summer",], aes(x = PercForMix, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Summer Cameras vs Availble Data", x = "Percent Forest within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_avail[[2]]$PercForMix[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[2]]$PercForMix[3], linetype = "dashed", color = "blue")
  elk_allfor_wtr <- ggplot(elk_sampled[elk_sampled$Season == "Winter",], aes(x = PercForMix, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Winter Cameras vs Availble Data", x = "Percent Forest within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_avail[[2]]$PercForMix[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[2]]$PercForMix[4], linetype = "dashed", color = "blue")
  
  #'  Winter GRASS histograms
  coy_grass_wtr <- ggplot(coy_obs[coy_obs$Season == "Winter",], aes(x = PercXGrass, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Coyote Winter Observations", x = "Percent Grass within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[7]]$PercXGrass[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[7]]$PercXGrass[4], linetype = "dashed", color = "blue")
  
  #'  SHRUB histograms
  md_shrub_smr <- ggplot(md_obs[md_obs$Season == "Summer",], aes(x = PercXShrub, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Summer Observations", x = "Percent Shrub within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[1]]$PercXShrub[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$PercXShrub[3], linetype = "dashed", color = "blue")
  md_shrub_wtr <- ggplot(md_obs[md_obs$Season == "Winter",], aes(x = PercXShrub, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Winter Observations", x = "Percent Shrub within 250m of Observation", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[1]]$PercXShrub[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$PercXShrub[4], linetype = "dashed", color = "blue")
  
  #'  ROAD DENSITY histograms
  elk_rd_smr <- ggplot(elk_obs[elk_obs$Season == "Summer",], aes(x = RoadDen, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 400, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Summer Observations", x = "Road Density (road length/1000 m2)", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[2]]$RoadDen[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$RoadDen[3], linetype = "dashed", color = "blue")
  elk_rd_wtr <- ggplot(elk_obs[elk_obs$Season == "Winter",], aes(x = RoadDen, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 400, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Elk Winter Observations", x = "Road Density (road length/1000 m2)", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[2]]$RoadDen[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[2]]$RoadDen[4], linetype = "dashed", color = "blue")
  coug_rd_smr <- ggplot(coug_obs[coug_obs$Season == "Summer",], aes(x = RoadDen, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 400, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Cougar Summer Observations", x = "Road Density (road length/1000 m2)", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[4]]$RoadDen[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[4]]$RoadDen[3], linetype = "dashed", color = "blue")
  coug_rd_wtr <- ggplot(coug_obs[coug_obs$Season == "Winter",], aes(x = RoadDen, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 500, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Cougar Winter Observations", x = "Road Density (road length/1000 m2)", y = "Proportion of Observations per Data Type")  + 
    geom_vline(xintercept = mu_obs[[4]]$RoadDen[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[4]]$RoadDen[4], linetype = "dashed", color = "blue")
  
  #'  HUMAN MODIFIED LANDSCAPE histograms
  md_hm_smr <- ggplot(md_obs[md_obs$Season == "Summer",], aes(x = HumanMod, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Summer Observations", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1)) +
    geom_vline(xintercept = mu_obs[[1]]$HumanMod[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$HumanMod[3], linetype = "dashed", color = "blue")
  md_hm_wtr <- ggplot(md_obs[md_obs$Season == "Winter",], aes(x = HumanMod, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Winter Observations", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1))  +
    geom_vline(xintercept = mu_obs[[1]]$HumanMod[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[1]]$HumanMod[4], linetype = "dashed", color = "blue")
  coy_hm_smr <- ggplot(coy_obs[coy_obs$Season == "Summer",], aes(x = HumanMod, color = Data, fill = Data)) +
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Coyote Summer Observations", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1))  +
    geom_vline(xintercept = mu_obs[[7]]$HumanMod[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[7]]$HumanMod[3], linetype = "dashed", color = "blue")
  coy_hm_wtr <- ggplot(coy_obs[coy_obs$Season == "Winter",], aes(x = HumanMod, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Coyote Winter Observations", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1))  +
    geom_vline(xintercept = mu_obs[[7]]$HumanMod[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_obs[[7]]$HumanMod[4], linetype = "dashed", color = "blue")
  #'  Mule deer camera samples vs what was available for human modified
  md_availhm_smr <- ggplot(md_sampled[md_sampled$Season == "Summer",], aes(x = HumanMod, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Summer Cameras vs Availble Data", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1)) +
    geom_vline(xintercept = mu_avail[[1]]$HumanMod[1], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[1]]$HumanMod[3], linetype = "dashed", color = "blue")
  md_availhm_wtr <- ggplot(md_sampled[md_sampled$Season == "Winter",], aes(x = HumanMod, color = Data, fill = Data)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.5, position = "identity", mapping = aes(y = stat(ncount))) + 
    labs(title = "Mule Deer Winter Cameras vs Availble Data", x = "Percent Human Modified Landscape", y = "Proportion of Observations per Data Type") + 
    coord_cartesian(xlim = c(0, 1))  +
    geom_vline(xintercept = mu_avail[[1]]$HumanMod[2], linetype = "dashed", color = "red") +
    geom_vline(xintercept = mu_avail[[1]]$HumanMod[4], linetype = "dashed", color = "blue")
  
  
  
  #'  Save select figures as PNG images
  ggsave("./Outputs/Figures/Histograms/MuleDeer_elev_smr_plot.png", md_elev_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/MuleDeer_elev_wtr_plot.png", md_elev_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/MuleDeer_hm_smr_plot.png", md_hm_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/MuleDeer_hm_wtr_plot.png", md_hm_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Coyote_hm_smr_plot.png", coy_hm_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Coyote_hm_wtr_plot.png", coy_hm_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_elev_smr_plot.png", elk_elev_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_elev_wtr_plot.png", elk_elev_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_allelev_smr_plot.png", elk_allelev_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_allelev_wtr_plot.png", elk_allelev_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_forest_smr_plot.png", elk_for_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_forest_wtr_plot.png", elk_for_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Coyote_hm_wtr_plot.png", coy_hm_wtr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_allforest_smr_plot.png", elk_allfor_smr, width = 7.3, units = "in")
  ggsave("./Outputs/Figures/Histograms/Elk_allforest_wtr_plot.png", elk_allfor_wtr, width = 7.3, units = "in")
  
  
  
  ####  6. Mapping predicted habitat use/selection  ####
  #'  Read in extracted covariate data for entire study areas (based on 1km grids)
  OK_covs <- read.csv("./Outputs/Tables/StudyAreaWide_OK_Covariates_2021-08-10.csv") %>%
    dplyr::select(-X) %>%
    transmute(
      obs = obs,
      Elev = scale(Elev),
      Slope = scale(Slope),
      RoadDen = scale(RoadDen),
      HumanMod = scale(HumanMod),
      PercForMix = scale(PercForestMix2),
      PercXGrass = scale(PercXericGrass),
      PercXShrub = scale(PercXericShrub),
      Area = 1,
      x = x,
      y = y
    )
  NE_covs <- read.csv("./Outputs/Tables/StudyAreaWide_NE_Covariates_2021-08-09.csv") %>%
    dplyr::select(-X) %>%
    transmute(
      obs = obs,
      Elev = scale(Elev),
      Slope = scale(Slope),
      RoadDen = scale(RoadDen),
      HumanMod = scale(HumanMod),
      PercForMix = scale(PercForestMix2),
      PercXGrass = scale(PercXericGrass),
      PercXShrub = scale(PercXericShrub),
      Area = 0,
      x = x,
      y = y
    )
  all_covs <- as.data.frame(rbind(NE_covs, OK_covs))
  
  #'  Convert into SpatialPointsDataFrame
  xy <- all_covs[,c(10,11)]
  covs_spdf <- SpatialPointsDataFrame(data = all_covs, coords = xy,
                                 proj4string = CRS(sa_proj))
  #'  Extract Lat/Long coordinates and append to covariate data frame
  xy_saproj <- coordinates(covs_spdf)
  colnames(xy_saproj) <- c("x2", "y2")
  all_covs <- cbind(all_covs, xy_saproj)
  #' #'  Transform to WGS84
  #' covs_wgs84 <- spTransform(covs_spdf, wgs84)
  #' #'  Extract Lat/Long coordinates and append to covariate data frame
  #' latlong <- coordinates(covs_wgs84)
  #' colnames(latlong) <- c("Long", "Lat")
  #' all_covs <- cbind(all_covs, latlong)
  
  
  ####  Predict occupancy probability across study area  ####
  #'  Manipulate Occupancy result tables (read in above at 3. OccMod vs RSF plots)
  #'  Version 1: Include all coefficients, even if non-significant
  occ_coefs <- occ_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDensity,
      B.hm = HumanMod,
      B.area = AreaOK
    ) %>%
    #'  Change NAs to 0 (no effect) for coefficients not included in species-specific models
    mutate(
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.area = ifelse(is.na(B.area), 0, B.area)
    )
  
  #'  Manipulate Occupancy result tables (read in above at 3. OccMod vs RSF plots)
  #'  Version 2. Exclude all non-significant coefficients
  occ_coefs_signif <- occ_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Use p-values to change non-significant coefficients (alpha-level = 0.1) to 0 so there is no effect
    mutate(Estimate = ifelse(Pval > 0.1, Estimate == 0, Estimate)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Problem is many of the intercepts were reduced to 0 which we don't want 
    #'  so we need to swap this intercept column with the original intercepts 
    #'  using the occ_coef data frame above 
    cbind(occ_coefs$alpha) %>%
    dplyr::select(-Intercept) %>%
    relocate("occ_coefs$alpha", .after = "Season") %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = occ_coefs$alpha,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDensity,
      B.hm = HumanMod,
      B.area = AreaOK
    ) %>%
    #'  Change NAs to 0 (no effect) for coefficients not included in species-specific models
    mutate(
      B.grass = ifelse(is.na(B.grass), 0, B.grass),
      B.shrub = ifelse(is.na(B.shrub), 0, B.shrub),
      B.area = ifelse(is.na(B.area), 0, B.area)
    )
 
  #'  Function to predict across all grid cells based on occupancy model results
  #'  Should end up with 1 predicted value per grid cell
  predict_occ <- function(cov, coef) {
    predict_odds <- c()
    predict_prob <- c()
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$alpha + coef$B.area*cov$Area[i] + coef$B.elev*cov$Elev[i] + 
                               coef$B.slope*cov$Slope[i]+ coef$B.for*cov$PercForMix[i] + 
                               coef$B.grass*cov$PercXGrass[i] + coef$B.shrub*cov$PercXShrub[i] + 
                               coef$B.rd*cov$RoadDen[i] + coef$B.hm*cov$HumanMod[i])
      predict_prob[i] <- predict_odds[i] / (1 + predict_odds[i])
    }
    predict_prob <- as.data.frame(predict_prob) %>%
      transmute(
        Predicted_Occ = round(predict_prob, 2)
      )
    return(predict_prob)
  }
  #'  Run estimates from occupancy sub-model through function to predict probability of use
  #'  Includes ALL coefficients, even if most are non-significant
  md_smr_predict_occ <- predict_occ(OK_covs, occ_coefs[occ_coefs$Species == "Mule Deer" & occ_coefs$Season == "Summer",])
  md_wtr_predict_occ <- predict_occ(OK_covs, occ_coefs[occ_coefs$Species == "Mule Deer" & occ_coefs$Season == "Winter",])
  elk_smr_predict_occ <- predict_occ(NE_covs, occ_coefs[occ_coefs$Species == "Elk" & occ_coefs$Season == "Summer",])
  elk_wtr_predict_occ <- predict_occ(NE_covs, occ_coefs[occ_coefs$Species == "Elk" & occ_coefs$Season == "Winter",])
  wtd_smr_predict_occ <- predict_occ(NE_covs, occ_coefs[occ_coefs$Species == "White-tailed Deer" & occ_coefs$Season == "Summer",])
  wtd_wtr_predict_occ <- predict_occ(NE_covs, occ_coefs[occ_coefs$Species == "White-tailed Deer" & occ_coefs$Season == "Winter",])
  coug_smr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Cougar" & occ_coefs$Season == "Summer",])
  coug_wtr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Cougar" & occ_coefs$Season == "Winter",])
  wolf_smr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Wolf" & occ_coefs$Season == "Summer",])
  wolf_wtr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Wolf" & occ_coefs$Season == "Winter",])
  bob_smr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Bobcat" & occ_coefs$Season == "Summer",])
  bob_wtr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Bobcat" & occ_coefs$Season == "Winter",])
  coy_smr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Coyote" & occ_coefs$Season == "Summer",])
  coy_wtr_predict_occ <- predict_occ(all_covs, occ_coefs[occ_coefs$Species == "Coyote" & occ_coefs$Season == "Winter",])
  #'  Includes ONLY significant coefficients per species and season
  md_smr_predict_occ_sgnf <- predict_occ(OK_covs, occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Summer",])
  md_wtr_predict_occ_sgnf <- predict_occ(OK_covs, occ_coefs_signif[occ_coefs_signif$Species == "Mule Deer" & occ_coefs_signif$Season == "Winter",])
  elk_smr_predict_occ_sgnf <- predict_occ(NE_covs, occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Summer",])
  elk_wtr_predict_occ_sgnf <- predict_occ(NE_covs, occ_coefs_signif[occ_coefs_signif$Species == "Elk" & occ_coefs_signif$Season == "Winter",])
  wtd_smr_predict_occ_sgnf <- predict_occ(NE_covs, occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_occ_sgnf <- predict_occ(NE_covs, occ_coefs_signif[occ_coefs_signif$Species == "White-tailed Deer" & occ_coefs_signif$Season == "Winter",])
  coug_smr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Summer",])
  coug_wtr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Cougar" & occ_coefs_signif$Season == "Winter",])
  wolf_smr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Wolf" & occ_coefs_signif$Season == "Winter",])
  bob_smr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Summer",])
  bob_wtr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Bobcat" & occ_coefs_signif$Season == "Winter",])
  coy_smr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Summer",])
  coy_wtr_predict_occ_sgnf <- predict_occ(all_covs, occ_coefs_signif[occ_coefs_signif$Species == "Coyote" & occ_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predator data that spans both study areas
  Predicted_occ <- all_covs %>%
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(coug_smr_predict_occ_sgnf, coug_wtr_predict_occ_sgnf, # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o non-signif coeffs)
          wolf_smr_predict_occ_sgnf, wolf_wtr_predict_occ_sgnf, 
          bob_smr_predict_occ_sgnf, bob_wtr_predict_occ_sgnf, 
          coy_smr_predict_occ_sgnf, coy_wtr_predict_occ_sgnf)
  #'  Make sure you have the order right when you change the names!!!
  colnames(Predicted_occ) <- c("obs", "Area", "x", "y",  
                               "COUG_smr_occ", "COUG_wtr_occ", 
                               "WOLF_smr_occ", "WOLF_wtr_occ", "BOB_smr_occ", 
                               "BOB_wtr_occ", "COY_smr_occ", "COY_wtr_occ")
  #'  Okanogan-only data (mule deer)
  OK_rows <- seq(1:nrow(md_smr_predict_occ))
  Area <- rep("Okanogan", length(OK_rows))
  OK_occ <- as.data.frame(cbind(OK_rows, Area, md_smr_predict_occ_sgnf, md_wtr_predict_occ_sgnf)) # KEEP TRACK of which version of the predicted results I'm using
  colnames(OK_occ) <- c("obs", "Area", "MD_smr_occ", "MD_wtr_occ")
  #'  Northeast-only data (elk & white-tailed deer)
  NE_rows <- seq(1:nrow(elk_smr_predict_occ))
  Area <- rep("Northeast", length(NE_rows))
  NE_occ <- as.data.frame(cbind(NE_rows, Area, elk_smr_predict_occ_sgnf, elk_wtr_predict_occ_sgnf, #KEEP TRACK of which version of the predicted results I'm using
                                wtd_smr_predict_occ_sgnf, wtd_wtr_predict_occ_sgnf)) 
  colnames(NE_occ) <- c("obs", "Area", "ELK_smr_occ", "ELK_wtr_occ", "WTD_smr_occ", "WTD_wtr_occ")
  #'  Merge ungulate & predator data by study area
  Predicted_occ_OK <- Predicted_occ[Predicted_occ$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_occ = NA,
      ELK_wtr_occ = NA,
      WTD_smr_occ = NA,
      WTD_wtr_occ = NA
    ) %>%
    full_join(OK_occ, by = c("obs", "Area")) 
  Predicted_occ_NE <- Predicted_occ[Predicted_occ$Area == "Northeast",] %>%
    full_join(NE_occ, by = c("obs", "Area")) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_occ = NA,
      MD_wtr_occ = NA
    )
  #'  Merge NE and OK predictions togther
  Predicted_occ <- as.data.frame(rbind(Predicted_occ_NE, Predicted_occ_OK))
  
  
  ####  Predict resource selection across study area  ####
  #'  Manipulate RSF result tables (read in above at 3. OccMod vs RSF plots)
  #'  Version 1: Include all coefficients, even if non-significant
  rsf_coefs <- rsf_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = Intercept,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDen,
      B.hm = HumanMod
    )
  
  #'  Version 2: Include ONLY significant coefficients
  rsf_coefs_signif <- rsf_out %>%
    dplyr::select(c(Species, Season, Parameter, Estimate, Pval)) %>%
    mutate(Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter)) %>%
    #'  Use p-values to change non-significant coefficients (alpha-level = 0.05) to 0 so there is no effect
    mutate(Estimate = ifelse(Pval > 0.05, Estimate == 0, Estimate)) %>%
    #'  For some reason the mutation above changes estimates that are already 0.00 with Pval > 0.05 to equal 1!?!?!
    #'  So changing those back to 0 since the effect is non-significant
    mutate(Estimate = ifelse(Estimate == 1, Estimate == 0, Estimate)) %>%
    dplyr::select(-Pval) %>%
    #'  Spread data so each row represents model coefficients for a single season, single species model
    pivot_wider(names_from = Parameter, values_from = Estimate) %>%
    #'  No intercepts should have been changed here but just in case- swap current
    #'  intercept column with the intercept column from above
    cbind(rsf_coefs$alpha) %>%
    dplyr::select(-Intercept) %>%
    relocate("rsf_coefs$alpha", .after = "Season") %>%
    #'  Rename coefficients so they're different than covariate names
    transmute(
      Species = Species,
      Season = Season,
      alpha = rsf_coefs$alpha,
      B.elev = Elev,
      B.slope = Slope,
      B.for = PercForMix,
      B.grass = PercXGrass,
      B.shrub = PercXShrub,
      B.rd = RoadDen,
      B.hm = HumanMod
    )
  
  #'  Function to predict across all grid cells based on RSF results
  #'  Should end up with 1 predicted value per grid cell
  predict_rsf <- function(cov, coef) {
    predict_odds <- c()
    predict_prob <- c()
    for(i in 1:nrow(cov)) {
      predict_odds[i] <- exp(coef$alpha + coef$B.elev*cov$Elev[i] + coef$B.slope*cov$Slope[i] + 
                               coef$B.for*cov$PercForMix[i] + coef$B.grass*cov$PercXGrass[i] + 
                               coef$B.shrub*cov$PercXShrub[i] + coef$B.rd*cov$RoadDen[i] + 
                               coef$B.hm*cov$HumanMod[i])
      predict_prob[i] <- predict_odds[i] / (1 + predict_odds[i])
    }
    predict_prob <- as.data.frame(predict_prob) %>%
      transmute(
        Predicted_RSF = round(predict_prob, 2)
      )
    return(predict_prob)
  }
  #'  Run estimated coefficients from RSFs through function to predict relative probability of selection
  #'  Includes ALL coefficients, even if a few are non-significant
  md_smr_predict_rsf <- predict_rsf(OK_covs, rsf_coefs[rsf_coefs$Species == "Mule Deer" & rsf_coefs$Season == "Summer",])
  md_wtr_predict_rsf <- predict_rsf(OK_covs, rsf_coefs[rsf_coefs$Species == "Mule Deer" & rsf_coefs$Season == "Winter",])
  elk_smr_predict_rsf <- predict_rsf(NE_covs, rsf_coefs[rsf_coefs$Species == "Elk" & rsf_coefs$Season == "Summer",])
  elk_wtr_predict_rsf <- predict_rsf(NE_covs, rsf_coefs[rsf_coefs$Species == "Elk" & rsf_coefs$Season == "Winter",])
  wtd_smr_predict_rsf <- predict_rsf(NE_covs, rsf_coefs[rsf_coefs$Species == "White-tailed Deer" & rsf_coefs$Season == "Summer",])
  wtd_wtr_predict_rsf <- predict_rsf(NE_covs, rsf_coefs[rsf_coefs$Species == "White-tailed Deer" & rsf_coefs$Season == "Winter",])
  coug_smr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Cougar" & rsf_coefs$Season == "Summer",])
  coug_wtr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Cougar" & rsf_coefs$Season == "Winter",])
  wolf_smr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Wolf" & rsf_coefs$Season == "Summer",])
  wolf_wtr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Wolf" & rsf_coefs$Season == "Winter",])
  bob_smr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Bobcat" & rsf_coefs$Season == "Summer",])
  bob_wtr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Bobcat" & rsf_coefs$Season == "Winter",])
  coy_smr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Coyote" & rsf_coefs$Season == "Summer",])
  coy_wtr_predict_rsf <- predict_rsf(all_covs, rsf_coefs[rsf_coefs$Species == "Coyote" & rsf_coefs$Season == "Winter",])
  #'  Includes ONLY significant coefficients
  md_smr_predict_rsf_sgnf <- predict_rsf(OK_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Summer",])
  md_wtr_predict_rsf_sgnf <- predict_rsf(OK_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Mule Deer" & rsf_coefs_signif$Season == "Winter",])
  elk_smr_predict_rsf_sgnf <- predict_rsf(NE_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Summer",])
  elk_wtr_predict_rsf_sgnf <- predict_rsf(NE_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Elk" & rsf_coefs_signif$Season == "Winter",])
  wtd_smr_predict_rsf_sgnf <- predict_rsf(NE_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Summer",])
  wtd_wtr_predict_rsf_sgnf <- predict_rsf(NE_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "White-tailed Deer" & rsf_coefs_signif$Season == "Winter",])
  coug_smr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Summer",])
  coug_wtr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Cougar" & rsf_coefs_signif$Season == "Winter",])
  wolf_smr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Summer",])
  wolf_wtr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Wolf" & rsf_coefs_signif$Season == "Winter",])
  bob_smr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Summer",])
  bob_wtr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Bobcat" & rsf_coefs_signif$Season == "Winter",])
  coy_smr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Summer",])
  coy_wtr_predict_rsf_sgnf <- predict_rsf(all_covs, rsf_coefs_signif[rsf_coefs_signif$Species == "Coyote" & rsf_coefs_signif$Season == "Winter",])
  
  #'  Combine into a monster data frame
  #'  Start with predators
  Predicted_rsf <- all_covs %>%
    dplyr::select(obs, Area, x, y) %>% 
    mutate(Area = ifelse(Area == 0, "Northeast", "Okanogan")) %>%
    cbind(coug_smr_predict_rsf_sgnf, coug_wtr_predict_rsf_sgnf, # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o non-signif coeffs)
          wolf_smr_predict_rsf_sgnf, wolf_wtr_predict_rsf_sgnf, 
          bob_smr_predict_rsf_sgnf, bob_wtr_predict_rsf_sgnf, 
          coy_smr_predict_rsf_sgnf, coy_wtr_predict_rsf_sgnf)
  #'  Make sure you have the order right when you change the names!!!
  colnames(Predicted_rsf) <- c("obs", "Area", "x", "y",  
                               "COUG_smr_rsf", "COUG_wtr_rsf", 
                               "WOLF_smr_rsf", "WOLF_wtr_rsf", "BOB_smr_rsf", 
                               "BOB_wtr_rsf", "COY_smr_rsf", "COY_wtr_rsf")
  #'  Okanogan-only data (mule deer)
  OK_rows <- seq(1:nrow(md_smr_predict_rsf))
  Area <- rep("Okanogan", length(OK_rows))
  OK_rsf <- as.data.frame(cbind(OK_rows, Area, md_smr_predict_rsf_sgnf, md_wtr_predict_rsf_sgnf)) # KEEP TRACK of which version of the predicted results I'm using
  colnames(OK_rsf) <- c("obs", "Area", "MD_smr_rsf", "MD_wtr_rsf")
  #'  Northeast-only data (elk & white-tailed deer)
  NE_rows <- seq(1:nrow(elk_smr_predict_rsf))
  Area <- rep("Northeast", length(NE_rows))
  NE_rsf <- as.data.frame(cbind(NE_rows, Area, elk_smr_predict_rsf_sgnf, elk_wtr_predict_rsf_sgnf, # KEEP TRACK of which version of the predicted results I'm using
                                wtd_smr_predict_rsf_sgnf, wtd_wtr_predict_rsf_sgnf)) 
  colnames(NE_rsf) <- c("obs", "Area", "ELK_smr_rsf", "ELK_wtr_rsf", "WTD_smr_rsf", "WTD_wtr_rsf")
  #'  Merge ungulate & predator data by study area
  Predicted_rsf_OK <- Predicted_rsf[Predicted_rsf$Area == "Okanogan",] %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      ELK_smr_rsf = NA,
      ELK_wtr_rsf = NA,
      WTD_smr_rsf = NA,
      WTD_wtr_rsf = NA
    ) %>%
    full_join(OK_rsf, by = c("obs", "Area")) 
  Predicted_rsf_NE <- Predicted_rsf[Predicted_rsf$Area == "Northeast",] %>%
    full_join(NE_rsf, by = c("obs", "Area")) %>%
    #'  Need to account for columns that are present in other study area dataframe
    mutate(
      MD_smr_rsf = NA,
      MD_wtr_rsf = NA
    )
  #'  Merge NE and OK predictions togther
  Predicted_rsf <- as.data.frame(rbind(Predicted_rsf_NE, Predicted_rsf_OK))
  
  #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
  Predicted_rsf_OK_rescale <- Predicted_rsf_OK %>%
    transmute(
      obs = obs,
      Area = Area,
      x = x,
      y = y,
      COUG_smr_rsf = round(COUG_smr_rsf/(max(COUG_smr_rsf, na.rm = T)), digits = 2),
      COUG_wtr_rsf = round(COUG_wtr_rsf/(max(COUG_wtr_rsf, na.rm = T)), digits = 2),
      WOLF_smr_rsf = round(WOLF_smr_rsf/(max(WOLF_smr_rsf, na.rm = T)), digits = 2),
      WOLF_wtr_rsf = round(WOLF_wtr_rsf/(max(WOLF_wtr_rsf, na.rm = T)), digits = 2),
      BOB_smr_rsf = round(BOB_smr_rsf/(max(BOB_smr_rsf, na.rm = T)), digits = 2),
      BOB_wtr_rsf = round(BOB_wtr_rsf/(max(BOB_wtr_rsf, na.rm = T)), digits = 2),
      COY_smr_rsf = round(COY_smr_rsf/(max(COY_smr_rsf, na.rm = T)), digits = 2),
      COY_wtr_rsf = round(COY_wtr_rsf/(max(COY_wtr_rsf, na.rm = T)), digits = 2),
      ELK_smr_rsf = NA,
      ELK_wtr_rsf = NA,
      WTD_smr_rsf = NA,
      WTD_wtr_rsf = NA,
      MD_smr_rsf = round(MD_smr_rsf/(max(MD_smr_rsf, na.rm = T)), digits = 2),
      MD_wtr_rsf = round(MD_wtr_rsf/(max(MD_wtr_rsf, na.rm = T)), digits = 2)
    )
  Predicted_rsf_NE_rescale <- Predicted_rsf_NE %>%
    transmute(
      obs = obs,
      Area = Area,
      x = x,
      y = y,
      COUG_smr_rsf = round(COUG_smr_rsf/(max(COUG_smr_rsf, na.rm = T)), digits = 2),
      COUG_wtr_rsf = round(COUG_wtr_rsf/(max(COUG_wtr_rsf, na.rm = T)), digits = 2),
      WOLF_smr_rsf = round(WOLF_smr_rsf/(max(WOLF_smr_rsf, na.rm = T)), digits = 2),
      WOLF_wtr_rsf = round(WOLF_wtr_rsf/(max(WOLF_wtr_rsf, na.rm = T)), digits = 2),
      BOB_smr_rsf = round(BOB_smr_rsf/(max(BOB_smr_rsf, na.rm = T)), digits = 2),
      BOB_wtr_rsf = round(BOB_wtr_rsf/(max(BOB_wtr_rsf, na.rm = T)), digits = 2),
      COY_smr_rsf = round(COY_smr_rsf/(max(COY_smr_rsf, na.rm = T)), digits = 2),
      COY_wtr_rsf = round(COY_wtr_rsf/(max(COY_wtr_rsf, na.rm = T)), digits = 2),
      ELK_smr_rsf = round(ELK_smr_rsf/(max(ELK_smr_rsf, na.rm = T)), digits = 2),
      ELK_wtr_rsf = round(ELK_wtr_rsf/(max(ELK_wtr_rsf, na.rm = T)), digits = 2),
      WTD_smr_rsf = round(WTD_smr_rsf/(max(WTD_smr_rsf, na.rm = T)), digits = 2),
      WTD_wtr_rsf = round(WTD_wtr_rsf/(max(WTD_wtr_rsf, na.rm = T)), digits = 2),
      MD_smr_rsf = NA,
      MD_wtr_rsf = NA
    )
  Predicted_rsf_rescale <- as.data.frame(rbind(Predicted_rsf_NE_rescale, Predicted_rsf_OK_rescale))
 

  ####  Calculate Correlations between OccMod & RSF Predictions ####
  #'  Evaluate correlation between predicted space use for each paired set of models
  predict_corr <- function(predict_occu, predict_rsfs) {
    #'  Identify the maximum value in the predicted RSFs
    m <- max(predict_rsfs, na.rm = TRUE)
    #'  Re-scale RSF values so they range 0 - 1 to match occupancy predictions
    stand_rsf <- as.data.frame(predict_rsfs) %>%
      mutate(scaled_rsf = predict_rsfs/m)
    #'  Combine occupancy and re-scaled RSF values
    predicted <- as.data.frame(cbind(predict_occu, stand_rsf$scaled_rsf)) %>%
      mutate(CellID = seq(1:nrow(.)))
    colnames(predicted) <- c("Occ_predictions", "RSF_rs_predictions", "CellID")
    #'  Calculated correlation between occupancy and rsf predictions
    pred_corr <- cor(predicted$Occ_predictions, predicted$RSF_rs_predictions, use = "complete.obs")
    #' #'  Plot correlation
    #' plot(predicted$Occ_predictions, predicted$RSF_rs_predictions)
    #' abline(a = 0, b = 1, col = "red")
    return(pred_corr)
  }
  #'  Run each predictions through function
  md_smr_corr <- predict_corr(Predicted_occ_OK$MD_smr_occ, Predicted_rsf_OK$MD_smr_rsf)
  md_wtr_corr <- predict_corr(Predicted_occ_OK$MD_wtr_occ, Predicted_rsf_OK$MD_wtr_rsf)
  elk_smr_corr <- predict_corr(Predicted_occ_NE$ELK_smr_occ, Predicted_rsf_NE$ELK_smr_rsf)
  elk_wtr_corr <- predict_corr(Predicted_occ_NE$ELK_wtr_occ, Predicted_rsf_NE$ELK_wtr_rsf)
  wtd_smr_corr <- predict_corr(Predicted_occ_NE$WTD_smr_occ, Predicted_rsf_NE$WTD_smr_rsf)
  wtd_wtr_corr <- predict_corr(Predicted_occ_NE$WTD_wtr_occ, Predicted_rsf_NE$WTD_wtr_rsf)
  coug_smr_corr <- predict_corr(Predicted_occ$COUG_smr_occ, Predicted_rsf$COUG_smr_rsf)
  coug_wtr_corr <- predict_corr(Predicted_occ$COUG_wtr_occ, Predicted_rsf$COUG_wtr_rsf)
  wolf_smr_corr <- predict_corr(Predicted_occ$WOLF_smr_occ, Predicted_rsf$WOLF_smr_rsf)
  wolf_wtr_corr <- predict_corr(Predicted_occ$WOLF_wtr_occ, Predicted_rsf$WOLF_wtr_rsf)
  bob_smr_corr <- predict_corr(Predicted_occ$BOB_smr_occ, Predicted_rsf$BOB_smr_rsf)
  bob_wtr_corr <- predict_corr(Predicted_occ$BOB_wtr_occ, Predicted_rsf$BOB_wtr_rsf)
  coy_smr_corr <- predict_corr(Predicted_occ$COY_smr_occ, Predicted_rsf$COY_smr_rsf)
  coy_wtr_corr <- predict_corr(Predicted_occ$COY_wtr_occ, Predicted_rsf$COY_wtr_rsf)
  
  #'  Wrangle results into a table
  spp <- rep(c("Mule Deer", "Elk", "White-tailed Deer", "Cougar", "Wolf", "Bobcat", "Coyote"), each = 2)
  season <- rep(c("Summer", "Winter"), 7)
  corr <- c(md_smr_corr, md_wtr_corr, elk_smr_corr, elk_wtr_corr, wtd_smr_corr, 
            wtd_wtr_corr, coug_smr_corr, coug_wtr_corr, wolf_smr_corr, wolf_wtr_corr,
            bob_smr_corr, bob_wtr_corr, coy_smr_corr, coy_wtr_corr)
  corr_results <- as.data.frame(cbind(spp, season, corr)) %>%
    transmute(
      Species = spp,
      Season = season,
      Correlation = as.numeric(corr),
      Correlation = round(Correlation, digits = 2)
    ) %>%
    arrange(Species)
  
  #'  Save correlations
  write.csv(corr_results, "./Outputs/Tables/Correlation_OccMod_RSF_Predictions.csv")  # KEEP TRACK of which version of the predicted results I'm using (w/ or w/o non-signif coeffs)
  
  ####  Plot predicted estimates  ####
  #'  Keep in mind the occupancy and RSF results are on very different scales
  #'  so the coloration is going to differ just because of that.
  #'  Is there a way to weight the RSF results so they higher selected areas show
  #'  up better?
  
  ####  MULE DEER  ####
  #'  Summer Occ
  md_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + #low = "azure" #low = "floralwhite"
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  md_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_OK_rescale, aes(x = x, y = y, fill = cut(MD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_manual(breaks = c(0.25, 0.5, 0.75, 1.0), na.value = "seashell4")
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  md_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_OK, aes(x = x, y = y, fill = cut(MD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  md_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_OK_rescale, aes(x = x, y = y, fill =cut(MD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  md_smr_map <- md_smr_occ_fig + plot_annotation(title = "Predicted Summer Mule Deer Space Use") + md_smr_rsf_fig 
  md_wtr_map <- md_wtr_occ_fig + plot_annotation(title = "Predicted Winter Mule Deer Space Use") + md_wtr_rsf_fig 
  md_predicted_map <- md_smr_occ_fig + plot_annotation(title = "Predicted Mule Deer Space Use") + 
    md_wtr_occ_fig + md_smr_rsf_fig + md_wtr_rsf_fig + plot_layout(ncol = 2)
  
  #'  Visualize
  plot(md_smr_map)
  plot(md_wtr_map)
  plot(md_predicted_map)
  
  #'  How do I thin out x-axis values so they aren't so bunched in panel figure?
  
  
  
  ####  ELK  ####
  #'  Summer Occ
  elk_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  elk_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_NE_rescale, aes(x = x, y = y, fill = cut(ELK_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  elk_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(ELK_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  elk_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_NE_rescale, aes(x = x, y = y, fill = cut(ELK_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  elk_smr_map <- elk_smr_occ_fig + plot_annotation(title = "Predicted Summer Elk Space Use") + elk_smr_rsf_fig 
  elk_wtr_map <- elk_wtr_occ_fig + plot_annotation(title = "Predicted Winter Elk Space Use") + elk_wtr_rsf_fig 
  elk_predicted_map <- elk_smr_occ_fig + plot_annotation(title = "Predicted Elk Space Use") + 
    elk_wtr_occ_fig + elk_smr_rsf_fig + elk_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(elk_smr_map)
  plot(elk_wtr_map)
  plot(elk_predicted_map)
  
  
  ####  WHITE-TAILED DEER  ####
  #'  Summer Occ
  wtd_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wtd_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_NE_rescale, aes(x = x, y = y, fill = cut(WTD_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  wtd_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ_NE, aes(x = x, y = y, fill = cut(WTD_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wtd_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_NE_rescale, aes(x = x, y = y, fill = cut(WTD_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  wtd_smr_map <- wtd_smr_occ_fig + plot_annotation(title = "Predicted Summer White-tailed Deer Space Use") + wtd_smr_rsf_fig 
  wtd_wtr_map <- wtd_wtr_occ_fig + plot_annotation(title = "Predicted Winter White-tailed Deer Space Use") + wtd_wtr_rsf_fig 
  wtd_predicted_map <- wtd_smr_occ_fig + plot_annotation(title = "Predicted White-tailed Deer Space Use") + 
    wtd_wtr_occ_fig + wtd_smr_rsf_fig +  wtd_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(wtd_smr_map)
  plot(wtd_wtr_map)
  plot(wtd_predicted_map)
  
  
  ####  COUGAR  ####
  #'  Summer Occ
  coug_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coug_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  coug_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COUG_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coug_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COUG_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  coug_smr_map <- coug_smr_occ_fig + plot_annotation(title = "Predicted Summer Cougar Space Use") + coug_smr_rsf_fig + plot_layout(ncol = 1)
  coug_wtr_map <- coug_wtr_occ_fig + plot_annotation(title = "Predicted Winter Cougar Space Use") + coug_wtr_rsf_fig + plot_layout(ncol = 1)
  coug_predicted_map <- coug_smr_occ_fig + plot_annotation(title = "Predicted Cougar Space Use") + 
    coug_wtr_occ_fig + coug_smr_rsf_fig + coug_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(coug_smr_map)
  plot(coug_wtr_map)
  plot(coug_predicted_map)
  
  
  ####  WOLF  ####
  #'  Summer Occ
  wolf_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  wolf_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  wolf_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(WOLF_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  wolf_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(WOLF_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  wolf_smr_map <- wolf_smr_occ_fig + plot_annotation(title = "Predicted Summer Wolf Space Use") + wolf_smr_rsf_fig + plot_layout(ncol = 1)
  wolf_wtr_map <- wolf_wtr_occ_fig + plot_annotation(title = "Predicted Winter Wolf Space Use") + wolf_wtr_rsf_fig + plot_layout(ncol = 1)
  wolf_predicted_map <- wolf_smr_occ_fig + plot_annotation(title = "Predicted Wolf Space Use") + 
    wolf_wtr_occ_fig + wolf_smr_rsf_fig + wolf_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(wolf_smr_map)
  plot(wolf_wtr_map)
  plot(wolf_predicted_map)
  
  
  ####  BOBCAT  ####
  #'  Summer Occ
  bob_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  bob_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  bob_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(BOB_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  bob_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(BOB_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  bob_smr_map <- bob_smr_occ_fig + plot_annotation(title = "Predicted Summer Bobcat Space Use") + bob_smr_rsf_fig + plot_layout(ncol = 1)
  bob_wtr_map <- bob_wtr_occ_fig + plot_annotation(title = "Predicted Winter Bobcat Space Use") + bob_wtr_rsf_fig + plot_layout(ncol = 1)
  bob_predicted_map <- bob_smr_occ_fig + plot_annotation(title = "Predicted Bobcat Space Use") + 
    bob_wtr_occ_fig + bob_smr_rsf_fig + bob_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(bob_smr_map)
  plot(bob_wtr_map)
  plot(bob_predicted_map)
  
  
  ####  COYOTE  ####
  #'  Summer Occ
  coy_smr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_smr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4", limits = c(0, 1)) + 
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Summer Occupancy Model") 
  #'  Summer RSF
  coy_smr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_smr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "YlGn", na.translate = F) +
    # scale_fill_gradient(low = "mintcream", high = "seagreen4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Summer Resource Selection Function")
  #'  Winter Occ
  coy_wtr_occ_fig <- ggplot() +
    geom_tile(data = Predicted_occ, aes(x = x, y = y, fill = cut(COY_wtr_occ, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4", limits = c(0, 1)) +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Probability \nof Site Use')  +
    ggtitle("Winter Occupancy Model") 
  #'  Winter RSF
  coy_wtr_rsf_fig <- ggplot() +
    geom_tile(data = Predicted_rsf_rescale, aes(x = x, y = y, fill = cut(COY_wtr_rsf, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))) +
    scale_fill_brewer(type = "seq", palette = "PuBu", na.translate = F) +
    # scale_fill_gradient(low = "azure", high = "dodgerblue4", na.value = "seashell4") +
    #'  Add study area outlines for reference
    geom_sf(data = OK_SA, fill = NA, color = "grey20") +
    geom_sf(data = NE_SA, fill = NA, color = "grey20") +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    #'  Change legend, axis, & main titles
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = 'Re-Scaled Relative \nProbability of Selection')  +
    ggtitle("Winter Resource Selection Function")
  
  #'  patchwork figures together:
  coy_smr_map <- coy_smr_occ_fig + plot_annotation(title = "Predicted Summer Coyote Space Use") + coy_smr_rsf_fig + plot_layout(ncol = 1)
  coy_wtr_map <- coy_wtr_occ_fig + plot_annotation(title = "Predicted Winter Coyote Space Use") + coy_wtr_rsf_fig + plot_layout(ncol = 1)
  coy_predicted_map <- coy_smr_occ_fig + plot_annotation(title = "Predicted Coyote Space Use") + 
    coy_wtr_occ_fig + coy_smr_rsf_fig + coy_wtr_rsf_fig + plot_layout(ncol = 2) #+ plot_layout(guides = 'collect')
  
  #'  Visualize
  plot(coy_smr_map)
  plot(coy_wtr_map)
  plot(coy_predicted_map)
  
  
  #'  Save figures as PNG images
  ggsave("./Outputs/Figures/Maps/MuleDeer_predict_smr_plot.png", md_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/MuleDeer_predict_wtr_plot.png", md_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/MuleDeer_predicted_plot.png", md_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Elk_predict_smr_plot.png", elk_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/Elk_predict_wtr_plot.png", elk_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/Elk_predicted_plot.png", elk_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/WTDeer_predict_smr_plot.png", wtd_smr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/WTDeer_predict_wtr_plot.png", wtd_wtr_map, width = 14.3, units = "in")
  ggsave("./Outputs/Figures/Maps/WTDeer_predicted_plot.png", wtd_predicted_map, width = 14.3, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Cougar_predict_smr_plot.png", coug_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Cougar_predict_wtr_plot.png", coug_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Cougar_predicted_plot.png", coug_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Wolf_predict_smr_plot.png", wolf_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Wolf_predict_wtr_plot.png", wolf_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Wolf_predicted_plot.png", wolf_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Bobcat_predict_smr_plot.png", bob_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Bobcat_predict_wtr_plot.png", bob_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Bobcat_predicted_plot.png", bob_predicted_map, width = 18.5, units = "in")
  
  ggsave("./Outputs/Figures/Maps/Coyote_predict_smr_plot.png", coy_smr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Coyote_predict_wtr_plot.png", coy_wtr_map, width = 18.5, units = "in")
  ggsave("./Outputs/Figures/Maps/Coyote_predicted_plot.png", coy_predicted_map, width = 18.5, units = "in")
  
  
  
  
  
  
  