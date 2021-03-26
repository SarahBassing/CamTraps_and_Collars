  ##  Final collar cleaning steps
  ##  Washington Predator-Prey Project
  ##  Nov. 16, 2020
  ##  Sarah Bassing
  ##  =========================================================
  ##  Script takes cleaned master GPS satellite data and takes final  steps to 
  ##  truncate and filter telemetry data for analyses specific to my project.
  ##     1. Truncating 2-wks after animal was captured to ensure any movements 
  ##        affected by the capture are excluded from analyses.
  ##     2. Thinning data to only include locations on the WPPP chosen 4-hr fix
  ##        schedule. Location times should be: 2:00, 6:00, 10:00, 14:00, 18:00, 
  ##        & 22:00 for all individuals.
  ##     3. Remove any individuals that have very few locations or are missing
  ##        a lot of data due to collar malfunctions or previous filtering.
  ##  This should produce the final data set to be used with HMMs to evaluate
  ##  habitat associations under different movement states.
  ##  =========================================================
  
  #  Clean work space and load libraries
  rm(list = ls())

  library(lubridate)
  library(tidyverse)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits=15) 
  
  #  Read in data
  #  Make sure specific columns are formatted correctly for data manipulation
  md_info <- read.csv("md_info 2020-11-16.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
      ) %>%
    select(-X)
  elk_info <- read.csv("elk_info 2020-11-16.csv")%>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    select(-X)
  wtd_info <- read.csv("wtd_info 2020-11-16.csv")%>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    select(-X)
  
  md_skinny <- read.csv("md_skinny 2020-11-16.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)
  elk_skinny <- read.csv("elk_skinny 2020-11-16.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)
  wtd_skinny <- read.csv("wtd_skinny 2020-11-16.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)

  #  Function to truncate & thin telemetry data for a final data set appropriate 
  #  for HMM analyses.
  Final_telem <- function(info, telem) {
    
    #  1. Truncate telemetry data
    #  Create empty data frame to fill iteratively
    trunk <- data.frame()
    #  How many individuals are looped over?
    nrow(info)
    #  Loop over every unique individual animal and...
    for(i in 1:nrow(info)){
      #  Take the individual animal ID
      ID <- droplevels(info$IndividualIdentifier[i])
      #  Take the animal's GPS collar serial number 
      SN <- info$GPSCollarSerialNumber[i]
      #  Buffer capture date to remove locations affected by capture event
      #  Suggested to only use data from 2 weeks after the capture data
      start <- info$CaptureDate[i] + 13
      #  Exclude locations 1 day before estimated mortality date
      end <- info$EndDate[i]
      
      #  Subset telemetry data to the specific individual
      collar <- subset(telem, CollarID == SN)
      #  Add a new column to the telemetry data with the animal's individual ID
      collar$ID <- ID
      #  truncate telemetry data by new start and end dates for that individual
      collartrunk <- subset(collar, Finaldt >= start & Finaldt <= end) 
      
      #  Append each unique animal's locations to a clean dataframe
      trunk <- rbind(trunk, collartrunk)
    }
    
    #  2. Thin truncated data to only include locations on correct fix schedule
    thin_trunk <- with(trunk, trunk[hour(Floordt) == 2 | hour(Floordt) == 6 |
                                      hour(Floordt) == 10 | hour(Floordt) == 14 |
                                      hour(Floordt) == 18 | hour(Floordt) == 22,])
    thin_fix <- as.data.frame(thin_trunk) %>%
      arrange(ID, Floordt) %>%
      #  Make sure lat/long are in a numeric format
      mutate(
        Latitude = as.numeric(Latitude),
        Longitude = as.numeric(Longitude)
      ) #%>%
      # #  3. Thin data to retain only the 1st location on the hour
      # #  Important when collar goes into mortality mode but animal is still alive- 
      # #  Fix rate increases but flooring process puts all those times on the hour
      # group_by(ID) %>%
      # distinct(Floordt, .keep_all = TRUE) %>%   # .keep_all = TRUE saves all columns
      # ungroup()

    return(thin_fix)
    
    # #  Organize by individual ID and chronological order of locations
    # #  Format data fields
    # clean <- clean %>%
    #   arrange(ID, Finaldt) %>%
    #   transmute(
    #     OBJECTID = OBJECTID,
    #     PositionID = PositionID,
    #     CollarID = CollarID,
    #     Latitude = Latitude,
    #     Longitude = Longitude,
    #     ObservationDateTimePST = ObservationDateTimePST,
    #     TransmissionDateTimePST = TransmissionDateTimePST,
    #     DbLoadedDateTimePST = DbLoadedDateTimePST,
    #     ValidLocation = as.numeric(ValidLocation),
    #     ValidDate = as.numeric(ValidDate),
    #     VEC_MortalityStatus = as.factor(as.character(VEC_MortalityStatus)),
    #     VEC_FixType = as.factor(as.character(VEC_FixType)),
    #     VEC_Origin = as.factor(as.character(VEC_Origin)),
    #     VEC_DOP = as.numeric(VEC_DOP),
    #     VEC_Height = as.numeric(VEC_Height),
    #     Project = Project,
    #     CaptureID = CaptureID,
    #     DeploymentID = DeploymentID,
    #     TransmitterID = TransmitterID,
    #     Species = Species,
    #     Sex = as.factor(as.character(Sex)),
    #     IndividualName = IndividualName,
    #     Fate = Fate,
    #     SerialNumber = SerialNumber,
    #     CaptureDate = CaptureDate,
    #     FateDate = FateDate,
    #     DateMortality = DateMortality,
    #     DaysDelta = DaysDelta,
    #     daytime = daytime,
    #     UTCdt = UTCdt,
    #     Finaldt = Finaldt,
    #     Floordt = Floordt,
    #     ID = ID
    #   )
    
    #return(trunk)
  }
  
  #  Run species-specific ID and telemetry data through the function
  md_final <- Final_telem(md_info, md_skinny)
  elk_final <- Final_telem(elk_info, elk_skinny)
  wtd_final <- Final_telem(wtd_info, wtd_skinny)  

  

  