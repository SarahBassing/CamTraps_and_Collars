  ##  Telemetry Data Cleaning
  ##  April 1, 2020
  ##  WPPP -  Prugh Lab, University of Washington
  ##  Taylor Ganz, updated by Sarah Bassing Sept. 2020
  ##  =========================================================
  ##  Script to combine WPPP GPS collar location data with unique animal IDs,
  ##  and truncate by capture date and mortality date (if applicable). This
  ##  script was originally written by Taylor Ganz in the Prugh Lab to prepare
  ##  ungulate location data for the WPPP.
  ##  =========================================================

  #  Clean workspace and install libraries
  rm(list = ls())
  
  library(lubridate)
  library(zoo)
  library(tidyverse)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits=15) 
  
  
  ####  =========================================================
  ####  Read in Capture, Mortality, and Telemetry data  ####
  
  #  Capture data (latest download: 09.11.20)
  md_cap <- read.csv("./Data/Capture (MD) 091120.csv")    
  elk_cap <- read.csv("./Data/Capture (Elk) 091120.csv")  
  wtd_cap <- read.csv("./Data/Capture (WTD) 091120.csv")  
  
  str(md_cap)#; head(md_cap)
  
  #  Mortality data (latest download: 09.11.20)
  md_mort <- read.csv("./Data/Mortality (MD) 091120.csv")   
  elk_mort <- read.csv("./Data/Mortality (Elk) 091120.csv") 
  wtd_mort <- read.csv("./Data/Mortality (WTD) 091120.csv") 
  
  str(md_mort)#; head(md_mort)
  
  #  Telemetry data (latest download: 09.11.20)
  #  Add column with date/time in a useable format
  md_tel <- read.csv("./Data/telem_md_091020.csv") %>%    
    mutate(daytime = mdy_hms(ObsDateTimePST))
  elk_tel <- read.csv("./Data/telem_elk_091020.csv") %>%  
    mutate(daytime = mdy_hms(ObsDateTimePST))
  wtd_tel <- read.csv() %>%
    mutate(daytime = mdy_hms(ObsDateTimePST))
  
  str(md_tel)#; head(md_tel)
  
  #  Note: mdy_hms automatically assigns date/time to UTC unless otherwise  
  #  changed with force_tz()
  
  ####  =========================================================
  ####  Combine Capture & Mortality data for each animal  ####
  
  #  Combine capture and mortality data by unique animal ID
  md_info <- full_join(md_cap, md_mort, by = "IndividualIdentifier") %>%
    transmute(
      IndividualIdentifier = as.factor(IndividualIdentifier),
      IndividualSpecies = IndividualSpecies.x,
      IndividualSex = IndividualSex.x,
      IndividualID = IndividualID,
      CaptureID =  CaptureID,
      LifeStage = LifeStage.x,
      CaptureDate = CaptureDate,
      GPSCollarSerialNumber = GPSCollarSerialNumber,
      MortalityID = MortalityID,
      MortalityLifeStage = LifeStage.y,
      LastLiveObservation = LastLiveObservation,
      EstimatedMortalityDate = EstimatedMortalityDate
    )
  #  Make sure you got them all
  length(unique(md_cap$IndividualIdentifier))
  length(unique(md_info$IndividualIdentifier))
  
  elk_info <- full_join(elk_cap, elk_mort, by = "IndividualIdentifier") %>%
    transmute(
      IndividualIdentifier = as.factor(IndividualIdentifier),
      IndividualSpecies = IndividualSpecies.x,
      IndividualSex = IndividualSex.x,
      IndividualID = IndividualID,
      CaptureID =  CaptureID,
      LifeStage = LifeStage.x,
      CaptureDate = CaptureDate,
      GPSCollarSerialNumber = GPSCollarSerialNumber,
      MortalityID = MortalityID,
      MortalityLifeStage = LifeStage.y,
      LastLiveObservation = LastLiveObservation,
      EstimatedMortalityDate = EstimatedMortalityDate
    )
  length(unique(elk_cap$IndividualIdentifier))
  length(unique(elk_info$IndividualIdentifier))
  
  wtd_info <- full_join(wtd_cap, wtd_mort, by = "IndividualIdentifier") %>%
    transmute(
      IndividualIdentifier = as.factor(IndividualIdentifier),
      IndividualSpecies = IndividualSpecies.x,
      IndividualSex = IndividualSex.x,
      IndividualID = IndividualID,
      CaptureID =  CaptureID,
      LifeStage = LifeStage.x,
      CaptureDate = CaptureDate,
      GPSCollarSerialNumber = GPSCollarSerialNumber,
      MortalityID = MortalityID,
      MortalityLifeStage = LifeStage.y,
      LastLiveObservation = LastLiveObservation,
      EstimatedMortalityDate = EstimatedMortalityDate
    )
  length(unique(wtd_cap$IndividualIdentifier))
  length(unique(wtd_info$IndividualIdentifier))
  
  #  Add end date to truncate telemetry data
  #  Choose an end date if there is no mortality- today
  lastend <- ymd(Sys.Date())
  #  Assign mortality date as end date if animal died
  md_info$enddate <- mdy(md_info$EstimatedMortalityDate)
  elk_info$enddate <- mdy(elk_info$EstimatedMortalityDate)
  wtd_info$enddate <- mdy(wtd_info$EstimatedMortalityDate)
  #  Fill in other end dates with chosen date if animal did not die
  md_info$enddate[is.na(md_info$enddate)] <- lastend
  elk_info$enddate[is.na(elk_info$enddate)] <- lastend
  wtd_info$enddate[is.na(wtd_info$enddate)] <- lastend
  

  #  Remove individuals that don't have location data
  #  Mule deer: 3969MD17, 3958MD17, R89ND18
  md_info <- md_info[md_info$IndividualIdentifier != "3969MD17" & 
                       md_info$IndividualIdentifier != "3958MD17" &  
                       md_info$IndividualIdentifier != "R89ND18",]
  
  ####  =========================================================
  ####  Connect location data to individual animal IDs  ####
  
  #  Create empty dataframe to fill iteritively
  MDclean <- data.frame()
  #  How many individuals are looped over?
  nrow(md_info)
  #  Loop over every unique individual animal and...
  for(i in 1:nrow(md_info)){
    #  Take the individual animal ID
    mdID <- droplevels(md_info$IndividualIdentifier[i])
    #  Take the animal's GPS collar serial number 
    mdSN <- md_info$GPSCollarSerialNumber[i]
    #  Buffer capture date to remove locations affected by capture event
    #  Suggested to only use data from 2 weeks after the capture data (some papers suggest 1 month)
    start <- mdy(md_info$CaptureDate[i]) + 14
    #  Exclude locations 1 day before estimated mortality date
    end <- md_info$enddate[i] - 1
    
    #  Subset telemetry data to the specific individual
    md <- subset(md_tel, CollarID == DoeSN)
    #  Add a new column to the telemetry data with the animal's individual ID
    md$ID <- mdID
    #  truncate telemetry data by new start and end dates for taht individual
    mdlive <- subset(md, daytime >= start & daytime <= end)

    #  Append each unique animal's locations to a clean dataframe
    MDclean <- rbind(MDclean, mdlive)
  }
  
  #  Organize by individual ID and chronological order of locations
  MDclean <- MDclean %>%
    arrange(ID, daytime) %>%
    #  Filter out aberrant locations
    filter(flgLocation != 1) %>%
    filter(flgDate != 1) %>%
    filter(flgActive != 0)
  # flgLocation == 1 indicate inaccurate fixes
  # flgDate == 1 indicate dates in the future (only issue for Telonics collars)
  # flgJurisdiction == 1 indicates collar outside WA State jurisdiction (e.g., Tribal land, Canada)
  # flgActive == 0 indicates locations when the collar was not active (not deployed)
  
  
  

  ####==== Extract the doe data for the first deer===####
  # test it out for 23066 = 12MD18
  DoeID <- droplevels(md_info$IndividualIdentifier[1])
  #DoeID <- droplevels(MDids$IndividualIdentifier[1])  # Not sure where this i.. comes from but DON'T fix this in the orginal csv
  DoeSN <- md_info$GPSCollarSerialNumber[1]
  start <- mdy(md_info$CaptureDate[1]) + 7 # Here you can set what kind of buffer you want on your capture date
  # only use data from 2 weeks after the capture data (some papers suggest 1 month)
  end <- md_info$enddate[1] - 1 #for movement data, stop it 1 day before estimated mortality date
  # note enddate is already converted into date format, dont need mdy()
  
  doe <- subset(md_tel, CollarID == DoeSN) #location data for the doe you want
  doe$ID <- DoeID
  doelive <- subset(doe, daytime >= start & daytime <= end) #the timespan for the data you want
  #View(doelive)
  
  #set up a new dateframe to store the clean location data
  MDclean <- doelive
  
  #  Now iterate through the rest of the mule deer
  (n <- nrow(md_info))
  for(i in 2:n){
    DoeID <- droplevels(md_info$IndividualIdentifier[i])  
    #DoeID <- droplevels(MDids$IndividualIdentifier[i]) 
    DoeSN <- md_info$GPSCollarSerialNumber[i]
    start <- mdy(md_info$CaptureDate[i]) + 7
    end <- md_info$enddate[i] - 1
    
    doe <- subset(md_tel, CollarID == DoeSN) 
    doe$ID <- DoeID
    doelive <- subset(doe, daytime >= start & daytime <= end)
    
    #add the data onto the clean dataset - this might be bad practice but it works...
    MDclean <- rbind(MDclean, doelive)
  }
  
  #  Take stock of the data set
  length(unique(MDclean$ID))  # unique animal IDs
  length(unique(MDclean$CollarID))  # unique collar IDs; should be less than # of animals if any were refurbished
  
  MDclean <- MDclean %>%
    arrange(ID, OBJECTID) %>%  # this doesn't seem to be arrange all rows properly...
    filter(flgLocation != 1) %>%
    filter(flgDate != 1) %>%
    filter(flgActive != 0)
  # flgLocation == 1 indicate inaccurate fixes
  # flgDate == 1 indicate dates in the future (only issue for Telonics collars)
  # flgJurisdiction == 1 indicates collar outside WA State jurisdiction (e.g., Tribal land, Canada)
  # flgActive == 0 indicates locations when the collar was not active (not deployed)
  
  #  Save data!
  write.csv(MDclean, paste0('MDclean ', Sys.Date(), '.csv'))
