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
  md_info <- read.csv("md_info 2020-11-17.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
      ) %>%
    select(-X)
  elk_info <- read.csv("elk_info 2020-11-17.csv")%>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    select(-X)
  wtd_info <- read.csv("wtd_info 2020-11-17.csv")%>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = as_date(CaptureDate)
    ) %>%
    select(-X)
  #  Created based on data provided by L.Satterfield
  cougwolf_info <- read.csv("cougwolf_info_2021-04-01.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = mdy(CaptureDate, tz = "America/Los_Angeles"),
      EndDate = mdy(EndDate, tz = "America/Los_Angeles")
    )
  #  Created based on data provided by B.Windell
  meso_info <- read.csv("meso_info_11162020.csv") %>%
    mutate(
      IndividualIdentifier = as.factor(as.character(IndividualIdentifier)),
      CaptureDate = mdy(CaptureDate, tz = "America/Los_Angeles"),
      EndDate = mdy(EndDate, tz = "America/Los_Angeles")
    )
  
  md_skinny <- read.csv("md_skinny 2020-11-17.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)
  elk_skinny <- read.csv("elk_skinny 2020-11-17.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)
  wtd_skinny <- read.csv("wtd_skinny 2020-11-17.csv") %>%
    mutate(daytime = mdy_hms(ObservationDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           Floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-X)
  # coug_skinny <- read.csv("./Data/Cougar_Vectronic_ATS_Spring2021_All.csv") %>%
  #   mutate(daytime = as.POSIXct(LMT_DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "America/Los_Angeles"),
  #          UTCdt = with_tz(daytime, "UTC"), 
  #          Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
  #          Floordt = floor_date(Finaldt, unit = "hour"))
  # wolf_skinny <- read.csv("./Data/Wolf_Vectronic_Spring2021_All.csv")
  # meso_skinny <- read.csv()

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
      ) %>%
      
      #  3. Thin data to retain only the 1st location on the hour
      #  Important when collar goes into mortality mode but animal is still alive-
      #  Fix rate increases but flooring process puts all those times on the hour
      group_by(ID) %>%
      distinct(Floordt, .keep_all = TRUE) %>%   # .keep_all = TRUE saves all columns
      ungroup()

    return(thin_fix)
  }
  
  #  Run species-specific ID and telemetry data through the function
  md_final <- Final_telem(md_info, md_skinny)
  elk_final <- Final_telem(elk_info, elk_skinny)
  wtd_final <- Final_telem(wtd_info, wtd_skinny)  

  #  Filter data to desired date ranges that match occupancy models' primary 
  #  sampling occasion (91 days; 13 weeks)
  Seasonal_telem <- function(telem) {
    #  Summer 2018: 07/01/2018 - 09/29/2018 
    telem_summer18 <- telem %>%
      filter(Floordt > "2018-07-01 00:00:00") %>%
      filter(Floordt < "2018-09-30 00:00:00") %>%
      mutate(
        Season = "Summer18",
        Year = "Year1"
      )
    #  Summer 2019: 07/01/2019 - 09/29/2019 
    telem_summer19 <- telem %>%
      filter(Floordt > "2019-07-01 00:00:00") %>%
      filter(Floordt < "2019-09-30 00:00:00") %>%
      mutate(
        Season = "Summer19",
        Year = "Year2"
      )
    #  Winter 2018-2019: 12/1/2018 - 03/1/2019
    telem_winter1819 <- telem %>%
      filter(Floordt > "2018-12-01 00:00:00") %>%
      filter(Floordt < "2019-03-02 00:00:00") %>%
      mutate(
        Season = "Winter1819",
        Year = "Year1"
      )
    #  Winter 2019-2020: 12/1/2019 - 02/29/2020 
    telem_winter1920 <- telem %>%
      filter(Floordt > "2019-12-01 00:00:00") %>%
      filter(Floordt < "2020-03-01 00:00:00")  %>%
      mutate(
        Season = "Winter1920",
        Year = "Year2"
      )
    #  Combine into single file
    telem_smwtr <- rbind(telem_summer18, telem_winter1819, telem_summer19, telem_winter1920)
    return(telem_smwtr)
  }

  #  Run species-specific telemetry data through function to filter by data range
  md_season <- Seasonal_telem(md_final)
  elk_season <- Seasonal_telem(elk_final)
  wtd_season <- Seasonal_telem(wtd_final)
  
  #  Same thing but include data from both fix schedules
  md_season2 <- Seasonal_telem(md_skinny)
  elk_season2 <- Seasonal_telem(elk_skinny)
  wtd_season2 <- Seasonal_telem(wtd_skinny)
  
  
  ####  Summary Stats on Fix Success & Accuracy  ####
  
  #  Function to calculate number of locations per individual and season
  #  91 day sampling period with 6 fixes/day = 546 locations if no missed locations
  sum_locs <- function(telem) {
    
    Nmb_locs <- telem %>%
      group_by(ID, Season) %>%
      summarise(count = n())
    
    return(Nmb_locs)
  }
  
  #  Run species-specific data through function based on preferred fix schedule
  md_counts <- sum_locs(md_season) 
  elk_counts <- sum_locs(elk_season) 
  wtd_counts <- sum_locs(wtd_season) 
  
  #'  How much data am I losing if I stick to only 1 fix schedule?
  md_counts2 <- sum_locs(md_season2) 
  elk_counts2 <- sum_locs(elk_season2) 
  wtd_counts2 <- sum_locs(wtd_season2)
  
  md_cnt <- full_join(md_counts, md_counts2, by = c("ID", "Season"))
  elk_cnt <- full_join(elk_counts, elk_counts2, by = c("ID", "Season"))
  wtd_cnt <- full_join(wtd_counts, wtd_counts2, by = c("ID", "Season"))
  
  #  Summary stats on the seasonal locations
  summary(md_counts$count); sd(md_counts$count)
  summary(elk_counts$count); sd(elk_counts$count)
  summary(wtd_counts$count); sd(wtd_counts$count)
  
  

  