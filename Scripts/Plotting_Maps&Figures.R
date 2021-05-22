  #'  ==================================
  #'  Maps and figures
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  ==================================
  #'  Plot study area map with camera locations and inset map showing study area
  #'  locations in relation to Washington State
  

  #'  Load libraries
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(RColorBrewer)
  library(sf)
  library(raster)
  library(tidyverse)
  
  #'  Read in camera locations
  station_covs <- read.csv("./Data/Camera_Station18-20_Covariates_2021-04-25.csv")
  CameraLocation <- station_covs$CameraLocation
  Year <- station_covs$Year
  
  #'  Make camera location data spatial and reproject to study area projection
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  cams <- st_as_sf(station_covs[,6:8], coords = c("Longitude", "Latitude"), crs = wgs84)
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  cams_reproj <- st_transform(cams, crs = sa_proj)
  cams_reproj$Year <- Year
  
  #'  Read in spatial data and reproject
  WA <- st_read("./Shapefiles/Washington_State_Boundary/WA_State_Geospatial_Open_Data", layer = "WA_State_Boundary") %>%
    st_transform(crs = sa_proj)
  OK_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "METHOW_SA") %>%
    st_transform(crs = sa_proj)
  OK_SA$NAME <- "Okanogan"
  NE_SA <- st_read("./Shapefiles/fwdstudyareamaps", layer = "NE_SA") %>%
    st_transform(crs = sa_proj)
  NE_SA$NAME <- "Northeast"
  dem <- raster("./Shapefiles/WA DEM rasters/WPPP_DEM_30m_reproj.tif")
  projection(WA)
  projection(OK_SA)
  projection(dem)
  extent(OK_SA)
  extent(NE_SA)
  
  #'  Reduce raster resolution and prep new DEM raster for ggplot
  dem_low <- aggregate(dem, fact = 10)
  # writeRaster(dem_low, file = "./Shapefiles/WA DEM rasters/dem_reproj_low", format = "GTiff")
  # dem_low <- raster("./Shapefiles/WA DEM rasters/dem_reproj_low.tif")
  dem_p_low <- rasterToPoints(dem_low)
  dem_p_df <- as.data.frame(dem_p_low)
  colnames(dem_p_df) <- c("x", "y", "value")
  
  
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
    scale_discrete_manual(aesthetics = "color", values = c("#601A4A", "#63ACBE")) +
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
  # png(file = "./Outputs/Figures/StudyAreas_Cameras1820.png",
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
  
  
  
  
  
  