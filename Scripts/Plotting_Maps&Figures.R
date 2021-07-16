  #'  =================================================
  #'  Maps and figures for camera vs collar manuscript
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  May 2021
  #'  =================================================
  #'  1. Plot study area map with camera locations and inset map showing study 
  #'     area locations in relation to Washington State
  #'  2. Plot occupancy and RSF results
  #'  =================================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(ggplot2)
  library(ggspatial)
  library(cowplot)
  library(patchwork)
  library(png)
  library(RCurl)
  library(RColorBrewer)
  library(rphylopic)
  library(sf)
  library(raster)
  library(tidyverse)
  
  
  ####  1. Map study area and camera locations  ####
  #'  ==============================================
  #'  Read in camera locations
  station_covs <- read.csv("./Data/Camera_Station18-20_Covariates_2021-04-25.csv")
  CameraLocation <- station_covs$CameraLocation
  Year <- station_covs$Year
  
  #'  Define projections
  wgs84 <- projection("+proj=longlat +datum=WGS84 +no_defs")
  sa_proj <- projection("+proj=lcc +lat_1=48.73333333333333 +lat_2=47.5 +lat_0=47 +lon_0=-120.8333333333333 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs ")
  
  #'  Make camera location data spatial and reproject to study area projection
  cams <- st_as_sf(station_covs[,6:8], coords = c("Longitude", "Latitude"), crs = wgs84)
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
  # dem_low <- aggregate(dem, fact = 10)
  # writeRaster(dem_low, file = "./Shapefiles/WA DEM rasters/dem_reproj_low", format = "GTiff")
  dem_low <- raster("./Shapefiles/WA DEM rasters/dem_reproj_low.tif")
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
  
  
  
  ####  2. Maps of capture effort and MCPs  ####
  #'  ==========================================
  #'  Plot study areas with capture locations and MCPs for each species and 
  #'  display in a paneled figure
  
  
  ####  3. Plot OccMod & RSF Results  ####
  #'  ====================================
  #'  Plot effect size and CI's for occ mod & RSF results for key covariates
  #'  that were generally significant in both models (elevation, % forest, % grass,
  #'  human modified landscape). Display in multiple panels by model and season.
  
  #'  Occupancy model output
  occ_out <- read.csv("./Outputs/OccMod_OccProb_Results_2021-07-01.csv") %>%
    #'  Calculate 95% confidence intervals
    mutate(
      l95 = (Estimate - (1.95 * SE)),
      u95 = (Estimate + (1.96 * SE))
    ) %>%
    dplyr::select(-c(X, Model))
  #'  RSF results output
  rsf_out <- read.csv("./Outputs/RSF_Results_2021-07-08.csv") %>%
    #'  Calculate 95% confidence intervals
    mutate(
      l95 = (Estimate - (1.95 * SE)),
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
  #'  Human Modified
  hm_occ <- filter(occ_out, Parameter == "HumanMod")
  hm_occ_smr <- filter(hm_occ, Season == "Summer")
  hm_occ_wtr <- filter(hm_occ, Season == "Winter")
  hm_rsf <- filter(rsf_out, Parameter == "HumanMod")
  hm_rsf_smr <- filter(hm_rsf, Season == "Summer")
  hm_rsf_wtr <- filter(hm_rsf, Season == "Winter")

  #'  Get silhouettes for each species from PhyloPic
  cougurl <- "http://phylopic.org/assets/images/submissions/3f8eff77-2868-4121-8d7d-a55ebdd49e04.64.png"
  cougimg <- readPNG(getURLContent(cougurl))
  wolfurl <- "http://phylopic.org/assets/images/submissions/8cad2b22-30d3-4cbd-86a3-a6d2d004b201.512.png"
  wolfimg <- readPNG(getURLContent(wolfurl))
  boburl <- "http://phylopic.org/assets/images/submissions/ab6cfd4f-aef7-40fa-b5a5-1b79b7d112aa.512.png"
  bobimg <- readPNG(getURLContent(boburl))
  coyurl <- "http://phylopic.org/assets/images/submissions/5a0398e3-a455-4ca6-ba86-cf3f1b25977a.512.png"
  coyimg <- readPNG(getURLContent(coyurl)) 
  mdurl <- "http://phylopic.org/assets/images/submissions/f889b336-9e67-4154-bc96-db4095a55be2.512.png"
  mdimg <- readPNG(getURLContent(mdurl))
  elkmurl <- "http://phylopic.org/assets/images/submissions/72f2f997-e474-4caf-bbd5-72fc8dbcc40d.512.png"
  elkmimg <- readPNG(getURLContent(elkmurl))
  elkfurl <- "http://phylopic.org/assets/images/submissions/97f83f5e-9afe-4ce8-812e-337f506ca841.512.png"
  elkfimg <- readPNG(getURLContent(elkfurl))
  wtdurl <- "http://phylopic.org/assets/images/submissions/56f6fdb2-15d0-43b5-b13f-714f2cb0f5d0.512.png"
  wtdimg <- readPNG(getURLContent(wtdurl))

  
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
  
  #'  Effect of PERCENT MIXED FOREST on Probability of Use (logit scale)
  #'  Summer results
  for_occ_smr_fig <- ggplot(for_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Summer Percent Forest within 250m", 
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
    labs(title = "", #"Percent Forest within 250m", 
         subtitle = "Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 3.5) +
    coord_flip()
  
  #'  Effect of PERCENT GRASS on Probability of Use (logit scale)
  #'  Warnings are OK- for wtd & elk where % grass cov was excluded from OccMods
  #'  Summer results
  grass_occ_smr_fig <- ggplot(grass_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Percent Grass within 250m", 
         subtitle = "Summer Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2.5) +
    coord_flip()
  #'  Winter results
  grass_occ_wtr_fig <- ggplot(grass_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Percent Grass within 250m", 
         subtitle = "Winter Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 12) +
    coord_flip()
  
  #'  Effect of PERCENT HUMAN MODIFIED LANDSCAPE on Probability of Use (logit scale)
  #'  Summer results
  hm_occ_smr_fig <- ggplot(hm_occ_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Percent of Human Modified Landscape", 
         subtitle = "Summer Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2) +
    coord_flip()
  #'  Winter results
  hm_occ_wtr_fig <- ggplot(hm_occ_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "Percent of Human Modified Landscape", 
         subtitle = "Winter Occupancy") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none") +
    ylim(-2.5, 2) +
    coord_flip()
  
  #'  ----------------------------  
  ####  RESOURCE SELECTION RESULTS  ####
  #'  ----------------------------
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
    ylim(-0.5, 0.5) +
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
    add_phylopic(wtdimg, x = 6.1, y = 0.2, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(mdimg, x = 5.05, y = 0.3, ysize = 0.65, color = "black", alpha = 1) +
    add_phylopic(elkmimg, x = 4.05, y = 0.2, ysize = 1, color = "black", alpha = 1) +
    add_phylopic(coyimg, x = 3.05, y = 0.3, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(cougimg, x = 2, y = 0.2, ysize = 0.5, color = "black", alpha = 1) +
    add_phylopic(bobimg, x = 1.05, y = 0.3, ysize = 0.4, color = "black", alpha = 1)
    # add_phylopic(wolfimg, x = 7.05, y = 0.5, ysize = 0.5, color = "black", alpha = 1) +
    # add_phylopic(wtdimg, x = 6.1, y = 0.4, ysize = 1, color = "black", alpha = 1) +
    # add_phylopic(mdimg, x = 5.05, y = 0.5, ysize = 0.65, color = "black", alpha = 1) +
    # add_phylopic(elkmimg, x = 4.05, y = 0.4, ysize = 1, color = "black", alpha = 1) +
    # add_phylopic(coyimg, x = 3.05, y = 0.5, ysize = 0.5, color = "black", alpha = 1) +
    # add_phylopic(cougimg, x = 2, y = 0.4, ysize = 0.5, color = "black", alpha = 1) +
    # add_phylopic(bobimg, x = 1.05, y = 0.5, ysize = 0.4, color = "black", alpha = 1)
  
  #'  Effect of PERCENT MIXED FOREST on relative probability of selection (logit scale)
  #'  Summer results
  for_rsf_smr_fig <- ggplot(for_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Forest within 250m", 
         subtitle = "Summer Selection") +
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
         subtitle = "Winter Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-5, 5) +
    coord_flip()
  
  #'  Effect of PERCENT GRASS on relative probability of selection (logit scale)
  #'  Summer results
  grass_rsf_smr_fig <- ggplot(grass_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent Grass within 250m", 
         subtitle = "Summer Selection") +
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
         subtitle = "Winter Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 1) +
    coord_flip()
  
  #'  Effect of PERCENT HUMAN MODIFIED LANDSCAPE on relative probability of selection (logit scale)
  #'  Summer results
  hm_rsf_smr_fig <- ggplot(hm_rsf_smr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Human Modified Landscape", 
         subtitle = "Summer Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-2.5, 1.5) +
    coord_flip()
  #'  Winter results
  hm_rsf_wtr_fig <- ggplot(hm_rsf_wtr, aes(x = Species, y = Estimate, label = Estimate)) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    geom_errorbar(aes(ymin = l95, ymax = u95, col = Species), width = 0) +
    geom_point(stat = 'identity', aes(col = Species), size = 3.5) +
    labs(title = "", #"Percent of Human Modified Landscape", 
         subtitle = "Winter Selection") +
    xlab("") + ylab("Estimates") +
    theme(legend.position = "none", 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    ylim(-1, 0.5) +
    coord_flip()
  
  
  ####  Pair OccMod and RSF plots  ####
  #' #'  cowplot version:
  #' elev_smr_fig <- plot_grid(elev_occ_smr_fig, elev_rsf_smr_fig, align = "h")
  #' elev_wtr_fig <- plot_grid(elev_occ_wtr_fig, elev_rsf_wtr_fig, align = "h")
  #' elev_fig <- plot_grid(elev_smr_fig, elev_wtr_fig, align = "h")
  #' plot(elev_smr_fig)
  #' plot(elev_wtr_fig)
  #' plot(elev_fig)
  #' 
  #' tst <- plot_grid(elev_occ_smr_fig, elev_rsf_smr_fig, elev_occ_wtr_fig, elev_rsf_wtr_fig, align = "h")
  #' plot(tst)
  
  #'  patchwork version:
  elev_fig <- elev_occ_smr_fig + elev_rsf_smr_fig + elev_occ_wtr_fig + elev_rsf_wtr_fig + plot_layout(ncol = 4)
  plot(elev_fig)
  

  #### Combine into 8 window panel 
  #### [cov1 smr occ, cov1 smr rsf] [cov1 wtr occ, cov1 wtr rsf]
  #### [cov2 smr occ, cov2 smr rsf] [cov2 wtr occ, cov2 wtr rsf] etc.

  
  
  
  
  