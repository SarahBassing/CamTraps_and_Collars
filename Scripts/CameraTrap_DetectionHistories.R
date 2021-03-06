  #'  ============================================
  #'  Camera Trap Detection Histories
  #'  Washington Predator-Prey Project
  #'  Sarah Bassing
  #'  January 2021
  #'  ============================================
  #'  Script to combine species detection data with camera station data (mainly,
  #'  camera operation dates and problems time periods). This uses the camtrapR 
  #'  package to generate species-specific encounter histories that can be used
  #'  for occupancy models.
  #'  
  #'  Combines: 
  #'  "full_camdata_DATE.csv" from Detections_by_Camera_Station.R in 
  #'      WPPP_CameraTrapping.Rproj
  #'     -Contains ALL detections of animals, humans, & vehicles (no empties),
  #'      camera coordinates and deployment covariates (cam height, etc.)
  #'  "All_Camera_Stations_18-19_updated_DATE.csv"
  #'     -Contains camera locations (including updated names/locations when a
  #'     camera was moved), deployment & pull dates, and problem dates
  #'      
  #'  
  #'  1. Drop unnecessary fields and format dates and times correctly
  #'  2. Truncate detection data to desired date range
  #'  3. Create camera operation table (deals with cameras that were temporarily
  #'     inoperable during desired date range)
  #'  4. Create species-specific detection probabilities
  #'  
  #'  Acknowledgments: Script is based on code originally provided by Mitch 
  #'  Parsons & Michael Havrda, UW SEFS.
  #'  ============================================
  
  #'  Clean workspace & load libraries
  rm(list = ls())
  
  library(camtrapR)
  library(chron)
  library(tidyverse)
  
  #'  Read in data, format, and filter
  #'  Heads up: DON'T format the dates in the cam_stations file yet!
  cam_stationsYr1 <- read.csv("G:/My Drive/1 Predator Prey Project/Field Work/Data Entry/All_Camera_Stations_18-19_updated_1.21.21.csv")
  cam_stationsYr2 <- read.csv("G:/My Drive/1 Predator Prey Project/Field Work/Data Entry/All_Camera_Stations_19-20.csv")
  cam_stations <- rbind(cam_stationsYr1, cam_stationsYr2)
    #  Remove Probs2 because no data in these columns in Year 2
    # dplyr::select(-c(Problem2_from, Problem2_to))
   
  # megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata_2021-01-21.csv") %>%
  # megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdataYr2_2021-03-02.csv") %>%
  megadata <- read.csv("G:/My Drive/1_Repositories/WPPP_CameraTrapping/Output/full_camdata18-20_2021-03-04.csv") %>%
    dplyr::select("File", "DateTime", "Date", "Time", "CameraLocation", 
                  "Camera_Lat", "Camera_Long", "Animal", "Human", "Vehicle", 
                  "Species", "HumanActivity", "Count") %>%
    filter(!grepl("Moultrie", CameraLocation)) %>%
    #  Need to have something in the Species column for each detection
    mutate(
      Species = ifelse(Human == "TRUE" | Human == "true", "Human", Species),
      Species = ifelse(Vehicle == "TRUE" | Vehicle == "true", "Vehicle", Species),
      Species = ifelse(Species == "", "NA", Species),
      HumanActivity = ifelse(HumanActivity == "", "NA", HumanActivity)
    ) %>%
    #  Remove rows where no detection occurred but snuck into this data set somehow
    filter(!(Animal == "FALSE" & Human == "FALSE" & Vehicle == "FALSE") | (Animal == "false" & Human == "false" & Vehicle == "false")) %>%
    filter(!is.na(Species)) %>%
    mutate(
      DateTime = as.POSIXct(DateTime,
                            format="%Y-%m-%d %H:%M:%S",tz="America/Los_Angeles"),
      Date = as.Date(Date, format = "%Y-%m-%d"),
      Time = chron(times = Time)
    )
  
  #'  Number rows and add to data frame
  ID <- as.data.frame(1:nrow(megadata))
  colnames(ID) <- "ID"
  megadata <- cbind(megadata, ID)
  #'  Create unique name for each individual image file
  megadata$Image <- str_c(megadata$CameraLocation, megadata$File, megadata$ID,  sep = "-")
  megadata <- dplyr::select(megadata, -ID)
  #'  Not sure this step was actually needed anymore but oh well
  
  
   
  #'  CONSIDER RARIFYING BY INDEPENDENT DETECTION EVENTS 
  #'  SHOULD I GENERATE INDEPENDENT DETECTION EVENTS FOR MULTIPLE TIME PERIODS 
  #'  (E.G., 5 min & 30 min GAPT BTWN DETECTIONS OF SAME SPECIES)?
  #'  DEFINITELY NEED TO CHANGE THIS WHEN IT COMES TO HUMAN & VEHICLE DETECTIONS!
  
  #'  Extract independent detections
  #'  Create a column identifying whether each image is an "independent" event
  #'  If camera site is diff from previous row then give unique value. If not then...
  #'  If species detected is diff from previous row at same site then give unique value. If not then...
  #'  If DateTime is >30 min from previous DateTime at same site for same species then give unique value. If not then...
  #'  Capture value is the same as that in the previous row.
  dat <- arrange(megadata, CameraLocation, DateTime)
  caps <- c()
  caps[1] <- 1
  for (i in 2:nrow(dat)){
    if (dat$CameraLocation[i-1] != dat$CameraLocation[i]) caps[i] = i
    else (if (dat$Species[i-1] != dat$Species[i]) caps[i] = i
          else (if (difftime(dat$DateTime[i], dat$DateTime[i-1], units = c("mins")) > 30) caps[i] = i
                else caps[i] = caps[i-1]))
  }

  caps <- as.factor(caps)

  #'  Add new column to larger data set
  capdata <- cbind(as.data.frame(dat), caps)

  #'  Retain only the first image from each unique detection event
  detections <- capdata %>%
    group_by(caps) %>%
    slice(1L) %>%
    ungroup()
  #'  Keep in mind the 30 minute gap doesn't make sense for human & vehicle detections
  
  
  #'  Filter dates to specific range 
  #'  Summer 2018: 07/01/2018 - 09/29/2018 (thirteen 7-day sampling periods)
  images_summer18 <- detections %>%
    filter(Date > "2018-06-30") %>%
    filter(Date < "2018-10-01") %>%
    dplyr::select("Image", "File", "CameraLocation", "DateTime", "Date", "Time", "Species") #"DateTimeOriginal", 
  #'  Winter 2018-2019: 12/1/2018 - 03/1/2019 (thirteen 7-day sampling periods)
  images_winter1819 <- detections %>%
    filter(Date > "2018-11-30") %>%
    filter(Date < "2019-03-02") %>%
    dplyr::select("Image", "File", "CameraLocation", "DateTime", "Date", "Time", "Species") #"DateTimeOriginal", 
  #'  Summer 2019: 07/01/2019 - 09/29/2019 (thirteen 7-day sampling periods)
  images_summer19 <- detections %>%
    filter(Date > "2019-06-30") %>%
    filter(Date < "2019-10-01") %>%
    dplyr::select("Image", "File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  #'  Winter 2019-2020: 12/1/2019 - 02/29/2020 (thirteen 7-day sampling periods)... 29 days in Feb 2020
  images_winter1920 <- detections %>%
    filter(Date > "2019-11-30") %>%
    filter(Date < "2020-03-01") %>%
    dplyr::select("Image", "File", "CameraLocation", "DateTime", "Date", "Time", "Species")  
  
  
  #'  Double check these were truncated correctly
  min(images_summer18$Date); max(images_summer18$Date)
  min(images_winter1819$Date); max(images_winter1819$Date)
  
  min(images_summer19$Date); max(images_summer19$Date)
  min(images_winter1920$Date); max(images_winter1920$Date)
  
  #'  Calculate number of trap nights per camera
  #'  Doesn't work unless dates are formatted correctly above- they're not
  # trapnights <- as.numeric(cam_stations$Pull_date - cam_stations$Set_date) 
  # hist(trapnights)
  
  ####  Camera Operation Table  ####
  #'  ------------------------------ 
  #'  Creates a matrix with each camera & dates it deployed
  #'  Define date format at this step, not before!
  #'  1 = operating; 0 = not operating but deployed; NA = not deployed
  #'  
  #'  Note: consider double cameras in the NE- is each camera a unique camera
  #'  station or are there 2 cameras at the same station? Decision will affect
  #'  how I set up data and these functions.
  #'  Grid cells this applies to: Yr1 NE2048, NE5345, "b" sites?; Yr2 NE5627
  
  camop_problem <- cameraOperation(CTtable = cam_stations,
                                   stationCol = "CameraLocation",
                                   setupCol = "Setup_date",
                                   retrievalCol = "Retrieval_date",
                                   hasProblems = TRUE,
                                   dateFormat = "%m/%d/%Y", # Define date format here!
                                   writecsv = FALSE) 

  probs <- as.data.frame(camop_problem)
  head(probs)
  
  
  ####  Detection Histories  ####
  #'  ---------------------------
  #'  Function to create season-specific detection histories for each species
  #'  for each season of interest (summer 2018, winter 2018-2019)
  #'  
  #'  Key arguments:
  #'  -occasionLength: currently using 7 day sampling occasions
  #'  -day1: sampling occasion begins upon station setup date ("station"), 
  #'   first day of survey ("survey"), or a specific date ("2018-06-15")
  #'   FYI this defines start date but NOT end date so DH goes until camera pulls
  #'  -includeEffort: compute # active trap days/station/occasion- effects DH
  #'   if FALSE then when camera is not set up or malfunctioning (NA or 0 in
  #'   camOp) during part of sampling occasion DH will be 0 for that occasion
  #'  -scaleEffort: center & scale effort (I plan to do this later in model)
  #'  -dateAsOccasionNames: if day1 = "survey" then uses 1st and last day of 
  #'   occasion as name for each sampling occasion
  #'  -output: return binary detections or counts of detections; don't want to
  #'   use "count" right now b/c would count each image not independent events
  #'  -occasionStartTime: time of day (full hour) at which to begin occasions
  #'   default is midnight (0)
  #'  -unmarkedMultFrameInput: create input for multi-season occmod in unmarked
  #'   ONLY use if running multi-season models & need to add more info to camop
  #'   
  #'  FYI, cannot have any NAs in the Species column or this doesn't work
  #'  Need to remove columns that extend beyond date range of interest!
  #'  Summer 2018: July 1 - Sept 29 = 13 weeks
  #'  Winter 2018-2019: Dec 1 - Mar 1 = 13 weeks 
  #'  Keep in mind detection data stops mid-sampling occasion of the 18th week
  #'  Consider changing detection data date range to encompass full 18th week
  
  ####  BOBCATS  ####
  DH_bob_smr18 <- detectionHistory(recordTable = images_summer18,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Bobcat",
                                    occasionLength = 7,
                                    day1 = "2018-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_bob_smr18 <- DH_bob_smr18[[1]][1:125,1:13]
  # DH_bob_smr18 <- DH_bob_smr18[[1]][,1:13]
  
  DH_bob_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Bobcat",
                                      occasionLength = 7,
                                      day1 = "2018-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_bob_wtr1819 <- DH_bob_wtr1819[[1]][1:125,1:13]
  # DH_bob_wtr1819 <- DH_bob_wtr1819[[1]][,1:13]
  
  DH_bob_smr19 <- detectionHistory(recordTable = images_summer19,
                                   camOp = camop_problem,
                                   stationCol = "CameraLocation",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTime",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   species = "Bobcat",
                                   occasionLength = 7,
                                   day1 = "2019-07-01", 
                                   # datesAsOccasionNames = TRUE,
                                   # occasionStartTime = 12, # starts at noon
                                   timeZone = "America/Los_Angeles",
                                   output = "binary",
                                   includeEffort = TRUE,
                                   scaleEffort = FALSE,
                                   # writecsv = TRUE,
                                   outDir = "./Data/Detection_Histories")
  DH_bob_smr19 <- DH_bob_smr19[[1]][126:242,1:13]
  # DH_bob_smr19 <- DH_bob_smr19[[1]][,1:13]
  
  DH_bob_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                     camOp = camop_problem,
                                     stationCol = "CameraLocation",
                                     speciesCol = "Species",
                                     recordDateTimeCol = "DateTime",
                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                     species = "Bobcat",
                                     occasionLength = 7,
                                     day1 = "2019-12-01", 
                                     # datesAsOccasionNames = TRUE,
                                     # occasionStartTime = 12, # starts at noon
                                     timeZone = "America/Los_Angeles",
                                     output = "binary",
                                     includeEffort = TRUE,
                                     scaleEffort = FALSE,
                                     # writecsv = TRUE,
                                     outDir = "./Data/Detection_Histories")
  DH_bob_wtr1920 <- DH_bob_wtr1920[[1]][126:242,1:13]
  # DH_bob_wtr1920 <- DH_bob_wtr1920[[1]][,1:13]
  
  #'  Merge summer 2018 and summer 2019 data together
  DH_bob_smr1819 <- rbind(DH_bob_smr18, DH_bob_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_bob_wtr1820 <- rbind(DH_bob_wtr1819, DH_bob_wtr1920)

  
  
  ####  COUGARS  ####
  DH_coug_smr18 <- detectionHistory(recordTable = images_summer18,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Cougar",
                                    occasionLength = 7,
                                    day1 = "2018-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_coug_smr18 <- DH_coug_smr18[[1]][1:125,1:13]
  # DH_coug_smr18 <- DH_coug_smr18[[1]][,1:13]

  DH_coug_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Cougar",
                                      occasionLength = 7,
                                      day1 = "2018-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_coug_wtr1819 <- DH_coug_wtr1819[[1]][1:125,1:13]
  # DH_coug_wtr1819 <- DH_coug_wtr1819[[1]][,1:13]
  
  DH_coug_smr19 <- detectionHistory(recordTable = images_summer19,
                                   camOp = camop_problem,
                                   stationCol = "CameraLocation",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTime",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   species = "Cougar",
                                   occasionLength = 7,
                                   day1 = "2019-07-01", 
                                   # datesAsOccasionNames = TRUE,
                                   # occasionStartTime = 12, # starts at noon
                                   timeZone = "America/Los_Angeles",
                                   output = "binary",
                                   includeEffort = TRUE,
                                   scaleEffort = FALSE,
                                   # writecsv = TRUE,
                                   outDir = "./Data/Detection_Histories")
  DH_coug_smr19 <- DH_coug_smr19[[1]][126:242,1:13]
  # DH_coug_smr19 <- DH_coug_smr19[[1]][,1:13]
  
  DH_coug_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                     camOp = camop_problem,
                                     stationCol = "CameraLocation",
                                     speciesCol = "Species",
                                     recordDateTimeCol = "DateTime",
                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                     species = "Cougar",
                                     occasionLength = 7,
                                     day1 = "2019-12-01", 
                                     # datesAsOccasionNames = TRUE,
                                     # occasionStartTime = 12, # starts at noon
                                     timeZone = "America/Los_Angeles",
                                     output = "binary",
                                     includeEffort = TRUE,
                                     scaleEffort = FALSE,
                                     # writecsv = TRUE,
                                     outDir = "./Data/Detection_Histories")
  DH_coug_wtr1920 <- DH_coug_wtr1920[[1]][126:242,1:13]
  # DH_coug_wtr1920 <- DH_coug_wtr1920[[1]][,1:13]

  #'  Merge summer 2018 and summer 2019 data together
  DH_coug_smr1819 <- rbind(DH_coug_smr18, DH_coug_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_coug_wtr1820 <- rbind(DH_coug_wtr1819, DH_coug_wtr1920)

  
  ####  COYOTES  ####
  DH_coy_smr18 <- detectionHistory(recordTable = images_summer18,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Coyote",
                                    occasionLength = 7,
                                    day1 = "2018-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_coy_smr18 <- DH_coy_smr18[[1]][1:125,1:13]
  # DH_coy_smr18 <- DH_coy_smr18[[1]][,1:13]
  
  DH_coy_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Coyote",
                                      occasionLength = 7,
                                      day1 = "2018-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_coy_wtr1819 <- DH_coy_wtr1819[[1]][1:125,1:13]
  # DH_coy_wtr1819 <- DH_coy_wtr1819[[1]][,1:13]
  
  DH_coy_smr19 <- detectionHistory(recordTable = images_summer19,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Coyote",
                                    occasionLength = 7,
                                    day1 = "2019-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_coy_smr19 <- DH_coy_smr19[[1]][126:242,1:13]
  
  DH_coy_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Coyote",
                                      occasionLength = 7,
                                      day1 = "2019-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_coy_wtr1920 <- DH_coy_wtr1920[[1]][126:242,1:13]
  
  #'  Merge summer 2018 and summer 2019 data together
  DH_coy_smr1819 <- rbind(DH_coy_smr18, DH_coy_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_coy_wtr1820 <- rbind(DH_coy_wtr1819, DH_coy_wtr1920)
  
  
  ####  WOLVES  ####
  DH_wolf_smr18 <- detectionHistory(recordTable = images_summer18,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Wolf",
                                    occasionLength = 7,
                                    day1 = "2018-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_wolf_smr18 <- DH_wolf_smr18[[1]][1:125,1:13]
  # DH_wolf_smr18 <- DH_wolf_smr18[[1]][,1:13]
  
  DH_wolf_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Wolf",
                                      occasionLength = 7,
                                      day1 = "2018-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_wolf_wtr1819 <- DH_wolf_wtr1819[[1]][1:125,1:13]
  # DH_wolf_wtr1819 <- DH_wolf_wtr1819[[1]][,1:13]
  
  DH_wolf_smr19 <- detectionHistory(recordTable = images_summer19,
                                   camOp = camop_problem,
                                   stationCol = "CameraLocation",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTime",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   species = "Wolf",
                                   occasionLength = 7,
                                   day1 = "2019-07-01", 
                                   # datesAsOccasionNames = TRUE,
                                   # occasionStartTime = 12, # starts at noon
                                   timeZone = "America/Los_Angeles",
                                   output = "binary",
                                   includeEffort = TRUE,
                                   scaleEffort = FALSE,
                                   # writecsv = TRUE,
                                   outDir = "./Data/Detection_Histories")
  DH_wolf_smr19 <- DH_wolf_smr19[[1]][126:242,1:13]
  
  DH_wolf_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                     camOp = camop_problem,
                                     stationCol = "CameraLocation",
                                     speciesCol = "Species",
                                     recordDateTimeCol = "DateTime",
                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                     species = "Wolf",
                                     occasionLength = 7,
                                     day1 = "2019-12-01", 
                                     # datesAsOccasionNames = TRUE,
                                     # occasionStartTime = 12, # starts at noon
                                     timeZone = "America/Los_Angeles",
                                     output = "binary",
                                     includeEffort = TRUE,
                                     scaleEffort = FALSE,
                                     # writecsv = TRUE,
                                     outDir = "./Data/Detection_Histories")
  DH_wolf_wtr1920 <- DH_wolf_wtr1920[[1]][126:242,1:13]
  
  #'  Merge summer 2018 and summer 2019 data together
  DH_wolf_smr1819 <- rbind(DH_wolf_smr18, DH_wolf_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_wolf_wtr1820 <- rbind(DH_wolf_wtr1819, DH_wolf_wtr1920)
  
  
  ####  ELK  ####
  DH_elk_smr18 <- detectionHistory(recordTable = images_summer18,
                                   camOp = camop_problem,
                                   stationCol = "CameraLocation",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTime",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   species = "Elk",
                                   occasionLength = 7,
                                   day1 = "2018-07-01", 
                                   # datesAsOccasionNames = TRUE,
                                   # occasionStartTime = 12, # starts at noon
                                   timeZone = "America/Los_Angeles",
                                   output = "binary",
                                   includeEffort = TRUE,
                                   scaleEffort = FALSE,
                                   # writecsv = TRUE,
                                   outDir = "./Data/Detection_Histories")
  DH_elk_smr18 <- DH_elk_smr18[[1]][1:125,1:13]
  # DH_elk_smr18 <- as.data.frame(DH_elk_smr18[[1]][,1:13]
  
  DH_elk_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                     camOp = camop_problem,
                                     stationCol = "CameraLocation",
                                     speciesCol = "Species",
                                     recordDateTimeCol = "DateTime",
                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                     species = "Elk",
                                     occasionLength = 7,
                                     day1 = "2018-12-01", 
                                     # datesAsOccasionNames = TRUE,
                                     # occasionStartTime = 12, # starts at noon
                                     timeZone = "America/Los_Angeles",
                                     output = "binary",
                                     includeEffort = TRUE,
                                     scaleEffort = FALSE,
                                     # writecsv = TRUE,
                                     outDir = "./Data/Detection_Histories")
  DH_elk_wtr1819 <- DH_elk_wtr1819[[1]][1:125,1:13]
  # DH_elk_wtr1819 <- DH_elk_wtr1819[[1]][,1:13]
  
  DH_elk_smr19 <- detectionHistory(recordTable = images_summer19,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Elk",
                                    occasionLength = 7,
                                    day1 = "2019-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_elk_smr19 <- DH_elk_smr19[[1]][126:242,1:13]
  
  DH_elk_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Elk",
                                      occasionLength = 7,
                                      day1 = "2019-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_elk_wtr1920 <- DH_elk_wtr1920[[1]][126:242,1:13]
  
  #'  Merge summer 2018 and summer 2019 data together
  DH_elk_smr1819 <- rbind(DH_elk_smr18, DH_elk_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_elk_wtr1820 <- rbind(DH_elk_wtr1819, DH_elk_wtr1920)
  
  
  ####  MULE DEER  ####
  DH_md_smr18 <- detectionHistory(recordTable = images_summer18,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "Mule Deer",
                                    occasionLength = 7,
                                    day1 = "2018-07-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_md_smr18 <- DH_md_smr18[[1]][1:125,1:13]
  # DH_md_smr18 <- DH_md_smr18[[1]][,1:13]
  
  DH_md_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                      camOp = camop_problem,
                                      stationCol = "CameraLocation",
                                      speciesCol = "Species",
                                      recordDateTimeCol = "DateTime",
                                      recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                      species = "Mule Deer",
                                      occasionLength = 7,
                                      day1 = "2018-12-01", 
                                      # datesAsOccasionNames = TRUE,
                                      # occasionStartTime = 12, # starts at noon
                                      timeZone = "America/Los_Angeles",
                                      output = "binary",
                                      includeEffort = TRUE,
                                      scaleEffort = FALSE,
                                      # writecsv = TRUE,
                                      outDir = "./Data/Detection_Histories")
  DH_md_wtr1819 <- DH_md_wtr1819[[1]][1:125,1:13]
  # DH_md_wtr1819 <- DH_md_wtr1819[[1]][,1:13]
  
  DH_md_smr19 <- detectionHistory(recordTable = images_summer19,
                                   camOp = camop_problem,
                                   stationCol = "CameraLocation",
                                   speciesCol = "Species",
                                   recordDateTimeCol = "DateTime",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   species = "Mule Deer",
                                   occasionLength = 7,
                                   day1 = "2019-07-01", 
                                   # datesAsOccasionNames = TRUE,
                                   # occasionStartTime = 12, # starts at noon
                                   timeZone = "America/Los_Angeles",
                                   output = "binary",
                                   includeEffort = TRUE,
                                   scaleEffort = FALSE,
                                   # writecsv = TRUE,
                                   outDir = "./Data/Detection_Histories")
  DH_md_smr19 <- DH_md_smr19[[1]][126:242,1:13]
  
  DH_md_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                     camOp = camop_problem,
                                     stationCol = "CameraLocation",
                                     speciesCol = "Species",
                                     recordDateTimeCol = "DateTime",
                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                     species = "Mule Deer",
                                     occasionLength = 7,
                                     day1 = "2019-12-01", 
                                     # datesAsOccasionNames = TRUE,
                                     # occasionStartTime = 12, # starts at noon
                                     timeZone = "America/Los_Angeles",
                                     output = "binary",
                                     includeEffort = TRUE,
                                     scaleEffort = FALSE,
                                     # writecsv = TRUE,
                                     outDir = "./Data/Detection_Histories")
  DH_md_wtr1920 <- DH_md_wtr1920[[1]][126:242,1:13]
  
  #'  Merge summer 2018 and summer 2019 data together
  DH_md_smr1819 <- rbind(DH_md_smr18, DH_md_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_md_wtr1820 <- rbind(DH_md_wtr1819, DH_md_wtr1920)
  
  
  ####  WHITE-TAILED DEER  ####
  DH_wtd_smr18 <- detectionHistory(recordTable = images_summer18,
                                  camOp = camop_problem,
                                  stationCol = "CameraLocation",
                                  speciesCol = "Species",
                                  recordDateTimeCol = "DateTime",
                                  recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                  species = "White-tailed Deer",
                                  occasionLength = 7,
                                  day1 = "2018-07-01", 
                                  # datesAsOccasionNames = TRUE,
                                  # occasionStartTime = 12, # starts at noon
                                  timeZone = "America/Los_Angeles",
                                  output = "binary",
                                  includeEffort = TRUE,
                                  scaleEffort = FALSE,
                                  # writecsv = TRUE,
                                  outDir = "./Data/Detection_Histories")
  DH_wtd_smr18 <- DH_wtd_smr18[[1]][1:125,1:13]
  # DH_wtd_smr18 <- DH_wtd_smr18[[1]][,1:13]
  
  DH_wtd_wtr1819 <- detectionHistory(recordTable = images_winter1819,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "White-tailed Deer",
                                    occasionLength = 7,
                                    day1 = "2018-12-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_wtd_wtr1819 <- DH_wtd_wtr1819[[1]][1:125,1:13]
  # DH_wtd_wtr1819 <- DH_wtd_wtr1819[[1]][,1:13]
  
  DH_wtd_smr19 <- detectionHistory(recordTable = images_summer19,
                                  camOp = camop_problem,
                                  stationCol = "CameraLocation",
                                  speciesCol = "Species",
                                  recordDateTimeCol = "DateTime",
                                  recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                  species = "White-tailed Deer",
                                  occasionLength = 7,
                                  day1 = "2019-07-01", 
                                  # datesAsOccasionNames = TRUE,
                                  # occasionStartTime = 12, # starts at noon
                                  timeZone = "America/Los_Angeles",
                                  output = "binary",
                                  includeEffort = TRUE,
                                  scaleEffort = FALSE,
                                  # writecsv = TRUE,
                                  outDir = "./Data/Detection_Histories")
  DH_wtd_smr19 <- DH_wtd_smr19[[1]][126:242,1:13]
  
  DH_wtd_wtr1920 <- detectionHistory(recordTable = images_winter1920,
                                    camOp = camop_problem,
                                    stationCol = "CameraLocation",
                                    speciesCol = "Species",
                                    recordDateTimeCol = "DateTime",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    species = "White-tailed Deer",
                                    occasionLength = 7,
                                    day1 = "2019-12-01", 
                                    # datesAsOccasionNames = TRUE,
                                    # occasionStartTime = 12, # starts at noon
                                    timeZone = "America/Los_Angeles",
                                    output = "binary",
                                    includeEffort = TRUE,
                                    scaleEffort = FALSE,
                                    # writecsv = TRUE,
                                    outDir = "./Data/Detection_Histories")
  DH_wtd_wtr1920 <- DH_wtd_wtr1920[[1]][126:242,1:13]

  #'  Merge summer 2018 and summer 2019 data together
  DH_wtd_smr1819 <- rbind(DH_wtd_smr18, DH_wtd_smr19)
  #'  Merge winter 2018/2019 and winter 2019/2020 data together
  DH_wtd_wtr1820 <- rbind(DH_wtd_wtr1819, DH_wtd_wtr1920)


  #' #'  Is it worth keeping data from Year1 cameras collected in 2019?
  #' #'  Look at which sites have data for both years
  #' #'  How many sites, how many sampling occasions, how many detections?
  #' #'  Requires DH include all rows (don't truncate at 1:125 or 126:242)
  #' 
  #' #'  Merge summer 2018 and summer 2019 data together
  #' DH_bob_sm1819 <- rbind(DH_bob_smr18, DH_bob_smr19)
  #' #'  Reorder so repeat camera sites are grouped together 
  #' orderDH <- DH_bob_sm1819[order(row.names(DH_bob_sm1819)),]
  #' #'  Drop rows that have all NAs so it's easier to see repeat sites with data
  #' repcams <- orderDH[rowSums(is.na(orderDH)) != ncol(orderDH), ]
  #' #'  Identify which repeat camera sites have data in 2019
  #' dups <- repcams[duplicated(row.names(repcams)),]
  #' (srvys <- apply(dups, MARGIN = 1, function(x) sum(!is.na(x))))
  #' nsrvys <- table(srvys)
  #' sum(nsrvys)
  #' #'  Histogram of the number of survey occasions for repeat camera sites
  #' hist(srvys, breaks = 20, xlim = c(0,14), main = "Number of Survey Occasions, Yr1 cams in 2019")
  #' #'  Identify how many detections occurred at repeat sites in 2019
  #' (ndets <- apply(dups, MARGIN = 1, function(x) sum((x), na.rm=TRUE)))
  #' hist(ndets, breaks = 3, xlim = c(0, 5), main = "Number of Bobcat Detections, Yr1 cams in 2019")
  #' #'  Number of detections
  #' sum(ndets)
  #' #'  Number of cameras where a detection occurred
  #' sum(ndets > 0)
  
  
  