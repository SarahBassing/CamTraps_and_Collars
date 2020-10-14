  ##  Ungulate Telemetry Data Cleaning
  ##  Washington Predator-Prey Project
  ##  Oct. 14, 2020
  ##  Sarah Bassing
  ##  =========================================================
  ##  Script to combine WPPP GPS collar location data with unique animal IDs,
  ##  and truncate by capture date and mortality date (if applicable). This
  ##  Script was originally written by Taylor Ganz in the Prugh Lab to prepare
  ##  mule deer location data for the WPPP and I have expanded it.
  ##  =========================================================

  #  Clean workspace and install libraries
  rm(list = ls())
  
  library(lubridate)
  library(zoo)
  library(purrr)
  library(tidyverse)
  
  #  Turn off scientific notation
  options(scipen = 999) 
  #  Set digits to 15 to ensure GPS coordinates aren't truncated
  options(digits=15) 
  
  
  ####  =========================================================
  ####  Read in Capture, Mortality, and Telemetry data  ####
  
  #  Capture data (latest download: 09.11.20)
  md_cap <- read.csv("./Data/Capture (MD) 091120.csv", stringsAsFactors = FALSE) 
  elk_cap <- read.csv("./Data/Capture (Elk) 091120.csv", stringsAsFactors = FALSE)  
  wtd_cap <- read.csv("./Data/Capture (WTD) 091120.csv", stringsAsFactors = FALSE)  
  
  #str(md_cap)#; head(md_cap)
  
  #  Mortality data (latest download: 09.11.20)
  md_mort <- read.csv("./Data/Mortality (MD) 091120.csv", stringsAsFactors = FALSE) 
  elk_mort <- read.csv("./Data/Mortality (Elk) 091120.csv", stringsAsFactors = FALSE) 
  wtd_mort <- read.csv("./Data/Mortality (WTD) 091120.csv", stringsAsFactors = FALSE) 
  
  #str(md_mort)#; head(md_mort)
  
  #  Taylor combined capture & mortality data (last downloaded 10.09.20)
  md_capmort <- read.csv("./Data/WPPP_md_2020-09-28.csv", stringsAsFactors = FALSE)
  elk_capmort <- read.csv("./Data/WPPP_elk_2020-10-5.csv", stringsAsFactors = FALSE)
  wtd_capmort <- read.csv("./Data/WPPP_wtd_2020-09-30.csv", stringsAsFactors = FALSE)
  
  #  Telemetry data (latest download: 09.11.20)
  #  Add column with date/time in usable format & set timezone that data downloaded in (includes daylight savings)
  #  Add column that converts times to UTC
  #  Add column that converts times to Pacific Standard Time (add 8 hours), ignoring daylight savings
  #  Add column that floors the time to the nearest hour not ahead
  md_tel <- read.csv("./Data/telem_md_091020.csv") %>%    
    mutate(daytime = mdy_hms(ObsDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-PositionID)
  elk_tel <- read.csv("./Data/telem_elk_091020.csv") %>%  
    mutate(daytime = mdy_hms(ObsDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-PositionID)
  wtd_tel <- read.csv("./Data/telem_wtd_091420.csv") %>%
    mutate(daytime = mdy_hms(ObsDateTimePST, tz = "America/Los_Angeles"),
           UTCdt = with_tz(daytime, "UTC"),
           Finaldt = with_tz(UTCdt, tzone = "Etc/GMT+8"),
           floordt = floor_date(Finaldt, unit = "hour")) %>%
    select(-PositionID)
  
  # #  Check out available timezones
  # OlsonNames()
  
  #  Check out the time data
  #  What's the timezone?
  tz(elk_tel$daytime); tz(elk_tel$UTCdt); tz(elk_tel$Finaldt)
  #  Does it account for daylight savings time?
  head(dst(elk_tel$daytime))
  head(dst(elk_tel$Finaldt))  # better say FALSE
  
  #str(elk_tel)#; head(elk_tel)
  
  ####  =========================================================
  ####  Combine Capture & Mortality data for each animal  ####
  
  spp_info <- function(cap, mort, capmort, tel) {
  #  Choose an end date if there is no mortality (chooses today's date)
    lastend <- ymd(Sys.Date())
  #  Combine capture and mortality data by unique animal ID               
    info <- full_join(cap, mort, by = "IndividualIdentifier") %>%
      #  Add in additional info about mortality & collar status
      right_join(capmort, by = c("IndividualIdentifier")) %>%
      transmute(
        IndividualIdentifier = as.factor(IndividualIdentifier),
        IndividualSpecies = IndividualSpecies.x,
        IndividualSex = IndividualSex.x,
        IndividualID = IndividualID,
        CaptureID =  CaptureID,
        LifeStage = LifeStage.x,
        CaptureDate = mdy(CaptureDate.x),
        GPSCollarSerialNumber = GPSCollarSerialNumber,
        MortalityID = MortalityID,
        MortalityLifeStage = LifeStage.y,
        EndDate = mdy(EndDate),
        EndCause = EndCause,
        MortalityType = MortalityType.x,
        MortalitySubType = MortalitySubType.x,
        LastTransmission = mdy(LastTransmission),
        Notes = Notes
        # LastLiveObservation = LastLiveObservation,
        # EstimatedMortalityDate = EstimatedMortalityDate
      )
  #  Fill in EndDate with chosen date if animal did not die/collar did not fail
    info$EndDate[is.na(info$EndDate)] <- lastend
  #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
    which(is.na(info$GPSCollarSerialNumber))
    info <- droplevels(info[!is.na(info$GPSCollarSerialNumber),])
  #  Remove individuals that died on capture date
    deadcap <- as.character(info$IndividualIdentifier[which(info$CaptureDate == info$EndDate)])
    if(is_empty(deadcap) != TRUE) {
      info <- info[!(info$IndividualIdentifier %in% deadcap),]
    } else {
      info <- info
    }
  #  Identify any GPS collars on captured deer that never generated telemetry data
    notel <- info$GPSCollarSerialNumber[!(info$GPSCollarSerialNumber %in% tel$CollarID)]
    if(is_empty(notel) != TRUE) {
      info <- droplevels(info[info$GPSCollarSerialNumber != notel,])
    } else {
      info <- info
    }
    
  }
  
  #  Run function to merge capture and mortality data for individual animals
  md_info <- spp_info(md_cap, md_mort, md_capmort, md_tel)
  elk_info <- spp_info(elk_cap, elk_mort, elk_capmort, elk_tel)  
  wtd_info <- spp_info(wtd_cap, wtd_mort, wtd_capmort, wtd_tel)  
  
  #  Remove this one rude mule deer that for some reason ruins the function... bad coding practice, I know!
  #  89MD18 and R89ND18 have the same collar ID b/c 89MD18 died 4 days after capture
  # md_info <- md_info[md_info$IndividualIdentifier != "89MD18",]   # ... seems to be working now
  

  ####  Running list of collars that got nixed in this stage  ####
  #  Mule deer
  #  3969MD17, 3958MD17, 89MD18
  
  #  Elk
  #  4836ELK20
  
  #  White-tailed deer
  #  019WTD17, 023WTD17, 70WTD18, 85WTD19, 90WTD19
  
  ####  =========================================================
  ####  Connect location data to individual animal IDs  ####
  
  #  IDtelem function attaches animal IDs to their respective telemetry locations.
  #  This function truncates the data so that day of capture and day of mortality 
  #  (or today's locations) are excluded from the dataset. Further truncating can 
  #  happen for individual analyses.

  #  flg columns can be used to filter out some locations
  #  flgLocation == 1 indicate inaccurate fixes
  #  flgDate == 1 indicate dates in the future (only issue for Telonics collars)
  #  flgJurisdiction == 1 indicates collar outside WA State jurisdiction (e.g., Tribal land, Canada)
  #  flgActive == 0 indicates locations where the animal that generated those locations 
  #  is no longer alive (do NOT filter out 0's here)
    
  IDtelem <- function(info, telem) {
    #  Create empty dataframe to fill iteratively
    clean <- data.frame()
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
      #  Currently not buffering but can change the +/- values in the future
      start <- info$CaptureDate[i] #+ 1
      #  Exclude locations 1 day before estimated mortality date
      end <- info$EndDate[i] #- 1
      
      #  Subset telemetry data to the specific individual
      collar <- subset(telem, CollarID == SN)
      #  Add a new column to the telemetry data with the animal's individual ID
      collar$ID <- ID
      #  truncate telemetry data by new start and end dates for that individual
      collarlive <- subset(collar, Finaldt >= start & Finaldt <= end)
      
      #  Append each unique animal's locations to a clean dataframe
      clean <- rbind(clean, collarlive)
    }
    
    #  Organize by individual ID and chronological order of locations
    clean <- clean %>%
      arrange(ID, Finaldt) %>%
      #  Filter out aberrant locations
      filter(flgLocation != 1) %>%
      filter(flgDate != 1)
    
    return(clean)
  }
  
  #  Run species-specific ID and telemetry data through the function
  md_clean <- IDtelem(md_info, md_tel)
  elk_clean <- IDtelem(elk_info, elk_tel)
  wtd_clean <- IDtelem(wtd_info, wtd_tel)  

  #  Pacific Standard Time for location data!

  #  Save data!
  write.csv(md_clean, paste0('md_clean ', Sys.Date(), '.csv'))
  write.csv(elk_clean, paste0('elk_clean ', Sys.Date(), '.csv'))
  write.csv(wtd_clean, paste0('wtd_clean ', Sys.Date(), '.csv'))
  
  
  #  Fin
  
  
  # ##  =====================================================  
  # ####  Broken down species by species  ####
  # 
  # #  This version is a little easier for trouble-shooting when the two functions
  # #  aren't working thanks to weird random collar-telemetry mismatches, etc.
  # 
  # 
  # #  Choose an end date if there is no mortality- today
  # lastend <- ymd(Sys.Date())
  # 
  # 
  # ####  Mule deer info ####
  # #  Combine capture and mortality data by unique animal ID
  # md_info <- full_join(md_cap, md_mort, by = "IndividualIdentifier") %>%
  #   right_join(md_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(md_cap$IndividualIdentifier))
  # length(unique(md_info$IndividualIdentifier))
  # 
  # #  Assign mortality date as end date if animal died
  # #md_info$enddate <- mdy(md_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # #md_info$enddate[is.na(md_info$enddate)] <- lastend
  # md_info$EndDate[is.na(md_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(md_info$IndividualIdentifier)
  # dat <- cbind(dat, md_info$GPSCollarSerialNumber)
  # view(dat)  # look for NAs and remove those individuals
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # #  Mule deer: 3969MD17 , 3958MD17
  # which(is.na(md_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(md_info$GPSCollarSerialNumber))))  # number of deer w/o collars
  # md_info <- droplevels(md_info[!is.na(md_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # deadcap <- as.character(md_info$IndividualIdentifier[which(md_info$CaptureDate == md_info$EstimatedMortalityDate)])
  # print(deadcap)
  # if(is_empty(deadcap) != TRUE) {
  #   md_info <- md_info[md_info$IndividualIdentifier != deadcap,]
  # } else {
  #   md_info <- md_info
  # }
  # #  Remove this one individual that for some reason ruins the function... bad coding practice
  # #  89MD18 and R89ND18 have the same collar ID b/c 89MD18 died 4 days after capture
  # md_info <- md_info[md_info$IndividualIdentifier != "89MD18",]
  # 
  # #  Identify any GPS collars on captured deer that never generated telemetry data
  # notel <- md_info$GPSCollarSerialNumber[!(md_info$GPSCollarSerialNumber %in% md_tel$CollarID)]
  # print(notel)
  # if(is_empty(notel) != TRUE) {
  #   md_info <- droplevels(md_info[md_info$GPSCollarSerialNumber != notel,])
  # } else {
  #   md_info <- md_info
  # }
  # 
  # 
  # ####  Elk info  ####
  # # Combine capture and mortality data by unique animal ID
  # elk_info <- full_join(elk_cap, elk_mort, by = "IndividualIdentifier") %>%
  #   right_join(elk_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(elk_cap$IndividualIdentifier))
  # length(unique(elk_info$IndividualIdentifier))
  # head(elk_info)
  # 
  # #  Assign mortality date as end date if animal died
  # #elk_info$enddate <- mdy(elk_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # lastend <- ymd(Sys.Date())
  # #elk_info$enddate[is.na(elk_info$enddate)] <- lastend
  # elk_info$EndDate[is.na(elk_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(elk_info$IndividualIdentifier)
  # dat <- cbind(dat, elk_info$GPSCollarSerialNumber)
  # view(dat)
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # #  Elk: 4836ELK20
  # which(is.na(elk_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(elk_info$GPSCollarSerialNumber))))
  # elk_info <- droplevels(elk_info[!is.na(elk_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # deadcap <- as.character(elk_info$IndividualIdentifier[which(elk_info$CaptureDate == elk_info$EstimatedMortalityDate)])
  # print(deadcap)
  # if(is_empty(deadcap) != TRUE) {
  #   elk_info <- elk_info[elk_info$IndividualIdentifier != deadcap,]
  # } else {
  #   elk_info <- elk_info
  # }
  # 
  # #  Identify any GPS collars on captured deer that never generated telemetry data
  # notel <- elk_info$GPSCollarSerialNumber[!(elk_info$GPSCollarSerialNumber %in% elk_tel$CollarID)]
  # print(notel)
  # if(is_empty(notel) != TRUE) {
  #   elk_info <- droplevels(elk_info[elk_info$GPSCollarSerialNumber != notel,])
  # } else {
  #   elk_info <- elk_info
  # }
  # 
  # 
  # ####  White-tailed deer info  ####
  # #  Combine capture and mortality data by unique animal ID
  # wtd_info <- full_join(wtd_cap, wtd_mort, by = "IndividualIdentifier") %>%
  #   right_join(wtd_capmort, by = c("IndividualIdentifier")) %>%
  #   transmute(
  #     IndividualIdentifier = as.factor(IndividualIdentifier),
  #     IndividualSpecies = IndividualSpecies.x,
  #     IndividualSex = IndividualSex.x,
  #     IndividualID = IndividualID,
  #     CaptureID =  CaptureID,
  #     LifeStage = LifeStage.x,
  #     CaptureDate = mdy(CaptureDate.x),
  #     GPSCollarSerialNumber = GPSCollarSerialNumber,
  #     MortalityID = MortalityID,
  #     MortalityLifeStage = LifeStage.y,
  #     EndDate = mdy(EndDate),
  #     EndCause = EndCause,
  #     MortalityType = MortalityType.x,
  #     MortalitySubType = MortalitySubType.x,
  #     LastTransmission = mdy(LastTransmission),
  #     Notes = Notes
  #   )
  # #  Make sure you got them all
  # length(unique(wtd_cap$IndividualIdentifier))
  # length(unique(wtd_info$IndividualIdentifier))
  # 
  # #  Assign mortality date as end date if animal died
  # #wtd_info$enddate <- mdy(wtd_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # lastend <- ymd(Sys.Date())
  # #wtd_info$enddate[is.na(wtd_info$enddate)] <- lastend
  # wtd_info$EndDate[is.na(wtd_info$EndDate)] <- lastend
  # 
  # #  Check to make sure each unique individual has location data
  # dat <- as.data.frame(wtd_info$IndividualIdentifier)
  # dat <- cbind(dat, wtd_info$GPSCollarSerialNumber)
  # view(dat)
  # 
  # #  Remove individuals that don't have a corresponding GPSCollarSerialNumber
  # which(is.na(wtd_info$GPSCollarSerialNumber))
  # length(unique(which(is.na(wtd_info$GPSCollarSerialNumber))))  # number of deer w/o collars
  # wtd_info <- droplevels(wtd_info[!is.na(wtd_info$GPSCollarSerialNumber),])
  # 
  # #  Remove individuals that died on capture date
  # #  019WTD17, 70WTD18, etc.
  # deadcap <- as.character(wtd_info$IndividualIdentifier[which(wtd_info$CaptureDate == wtd_info$EndDate)])
  # if(is_empty(deadcap) != TRUE) {
  #   wtd_info <- wtd_info[!(wtd_info$IndividualIdentifier %in% deadcap),]
  #   #wtd_info <- wtd_info[wtd_info$IndividualIdentifier != deadcap,]
  # } else {
  #   wtd_info <- wtd_info
  # }
  # 
  # #  Remove individuals with collars that are not in the telemetry data
  # #  Identify mismatches between GPS collars in wtd_info vs wtd_tel
  # notel <- wtd_info$GPSCollarSerialNumber[!(wtd_info$GPSCollarSerialNumber %in% wtd_tel$CollarID)]
  # print(notel)
  # #  Collars 24833 (90WTD19), 24867 (85WTD19)
  # wtd_info <- droplevels(wtd_info[wtd_info$GPSCollarSerialNumber != notel,])
  # 
  # 
  # #### Merge info & telem for 1 individual  ####
  # # If the IDtelem function is too much, try it for a single individual!
  # 
  # #  Create empty dataframe to fill iteritively
  # MDclean <- data.frame()
  # #  How many individuals are looped over?
  # nrow(elk_info) #md_info
  # #  Loop over every unique individual animal and...
  # for(i in 1:nrow(elk_info)){ #md_info      # DON'T FOR LOOP IT IF ONLY TESTING 1 INDIVIDUAL
  #   #  Take the individual animal ID
  #   mdID <- droplevels(elk_info$IndividualIdentifier[1]) #md_info
  #   #  Take the animal's GPS collar serial number
  #   mdSN <- elk_info$GPSCollarSerialNumber[1] #md_info
  #   #  Buffer capture date to remove locations affected by capture event
  #   #  Suggested to only use data from 2 weeks after the capture data (some papers suggest 1 month)
  #   start <- elk_info$CaptureID[1] + 14
  #   #  Exclude locations 1 day before estimated mortality date
  #   end <- elk_info$EndDate[i] - 1
  # 
  #   #  Subset telemetry data to the specific individual
  #   md <- subset(elk_tel, CollarID == mdSN) #md_tel
  #   #  Add a new column to the telemetry data with the animal's individual ID
  #   md$ID <- mdID
  #   #  Truncate telemetry data by new start and end dates for that individual
  #   #  Important for collars that are redeployed- ensures locations generated
  #   #  by that specific animal are included, even if collar generates more locations
  #   #  on another animal
  #   mdlive <- subset(md, Finaldt >= start & Finaldt <= end)
  # 
  #   #  Append each unique animal's locations to a clean dataframe
  #   MDclean <- rbind(MDclean, mdlive)
  # }
  # 
  # length(unique(wtd_cap$IndividualIdentifier))
  # length(unique(wtd_info$IndividualIdentifier))
  # length(unique(MDclean$ID))
  # #  FYI 89MD18 & 24MD18 are dropped in when locations are truncated b/c
  # #  animals died within 2 weeks of capture
  # 
  # #  Organize by individual ID and chronological order of locations
  # MDclean <- MDclean %>%
  #   arrange(ID, daytime) %>%
  #   #  Filter out aberrant locations
  #   filter(flgLocation != 1) %>%
  #   filter(flgDate != 1) %>%
  # # flgLocation == 1 indicate inaccurate fixes
  # # flgDate == 1 indicate dates in the future (only issue for Telonics collars)
  # # flgJurisdiction == 1 indicates collar outside WA State jurisdiction (e.g., Tribal land, Canada)
  # # flgActive == 0 indicates locations where the animal that generated those locations is no longer alive (do NOT filter out 0's here)
  