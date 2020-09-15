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
  # md_cap <- read.csv("G:/My Drive/1_Repositories/CamTraps_and_Collars/Data/Capture (MD) 091120.csv")
  elk_cap <- read.csv("./Data/Capture (Elk) 091120.csv")  
  wtd_cap <- read.csv("./Data/Capture (WTD) 091120.csv")  
  
  str(md_cap)#; head(md_cap)
  
  #  Mortality data (latest download: 09.11.20)
  md_mort <- read.csv("./Data/Mortality (MD) 091120.csv") 
  # md_mort <- read.csv("G:/My Drive/1_Repositories/CamTraps_and_Collars/Data/Mortality (MD) 091120.csv")
  elk_mort <- read.csv("./Data/Mortality (Elk) 091120.csv") 
  wtd_mort <- read.csv("./Data/Mortality (WTD) 091120.csv") 
  
  str(md_mort)#; head(md_mort)
  
  #  Telemetry data (latest download: 09.11.20)
  #  Add column with date/time in a useable format
  md_tel <- read.csv("./Data/telem_md_091020.csv") %>%    
    mutate(daytime = mdy_hms(ObsDateTimePST))
  # md_tel <- read.csv("G:/My Drive/1_Repositories/CamTraps_and_Collars/Data/telem_md_091020.csv") %>%    
  #   mutate(daytime = mdy_hms(ObsDateTimePST))
  elk_tel <- read.csv("./Data/telem_elk_091020.csv") %>%  
    mutate(daytime = mdy_hms(ObsDateTimePST))
  wtd_tel <- read.csv("./Data/telem_wtd_091420.csv") %>%
    mutate(daytime = mdy_hms(ObsDateTimePST))
  
  str(md_tel)#; head(md_tel)
  
  #  Note: mdy_hms automatically assigns date/time to UTC unless otherwise  
  #  changed with force_tz()
  
  ####  =========================================================
  ####  Combine Capture & Mortality data for each animal  ####
  
  #  Choose an end date if there is no mortality- today
  lastend <- ymd(Sys.Date())
  
  
  ####  Mule deer  ####
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
  
  #  Assign mortality date as end date if animal died
  md_info$enddate <- mdy(md_info$EstimatedMortalityDate)
  #  Fill in other end dates with chosen date if animal did not die
  md_info$enddate[is.na(md_info$enddate)] <- lastend
  
  #  Check to make sure each unique individual has location data
  dat <- as.data.frame(md_info$IndividualIdentifier)
  dat <- cbind(dat, md_info$GPSCollarSerialNumber)
  view(dat)  # look for NAs and remove those individuals
  
  #  Remove individuals that don't have collars/location data
  #  Mule deer: 3969MD17, 3958MD17, R89ND18 (or should it be 89MD18?)
  which(is.na(md_info$GPSCollarSerialNumber))
  md_info <- md_info[!is.na(md_info$GPSCollarSerialNumber),]
  #  89MD18 and R89ND18 have the same collar ID b/c 89MD18 died 4 days after capture
  md_info <- md_info[md_info$IndividualIdentifier != "89MD18",]
  # md_info <- md_info[md_info$IndividualIdentifier != "3969MD17" &
  #                      md_info$IndividualIdentifier != "3958MD17" &
  #                      md_info$IndividualIdentifier != "89MD18",] 
  
  
  ####  Elk  ####
  #  Combine capture and mortality data by unique animal ID
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
  #  Make sure you got them all
  length(unique(elk_cap$IndividualIdentifier))
  length(unique(elk_info$IndividualIdentifier))
  head(elk_info)
  
  #  Assign mortality date as end date if animal died
  elk_info$enddate <- mdy(elk_info$EstimatedMortalityDate)
  #  Fill in other end dates with chosen date if animal did not die
  elk_info$enddate[is.na(elk_info$enddate)] <- lastend
  
  #  Check to make sure each unique individual has location data
  dat <- as.data.frame(elk_info$IndividualIdentifier)
  dat <- cbind(dat, elk_info$GPSCollarSerialNumber)
  view(dat)
  
  #  Remove individuals that don't have collars/location data
  #  Elk: 4836ELK20
  which(is.na(elk_info$GPSCollarSerialNumber))
  elk_info <- elk_info[!is.na(elk_info$GPSCollarSerialNumber),]
  #elk_info <- elk_info[elk_info$IndividualIdentifier != "4836ELK20",]
  
  
  ####  White-tailed deer  ####
  #  Combine capture and mortality data by unique animal ID
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
  #  Make sure you got them all
  length(unique(wtd_cap$IndividualIdentifier))
  length(unique(wtd_info$IndividualIdentifier))
  
  #  Assign mortality date as end date if animal died
  wtd_info$enddate <- mdy(wtd_info$EstimatedMortalityDate)
  #  Fill in other end dates with chosen date if animal did not die
  wtd_info$enddate[is.na(wtd_info$enddate)] <- lastend
  
  #  What about animals that were last observed alive on 1 day but died way later? 
  #  Does that mean the collar stopped working then?
  #  How should the locations for these animals be truncated?
  #  e.g., 013WTD17, 023WTD17, etc.
  #  What about animals where their "LastLiveObservation" makes no sense?
  #  047WTD18 observed alive on 3/22/18 but died 3/5/18???
  #  4803WTD20 observed alive on 1/3/20 but captured 1/23/20, then died 1/31/20???
  
  #  Check to make sure each unique individual has location data
  dat <- as.data.frame(wtd_info$IndividualIdentifier)
  dat <- cbind(dat, wtd_info$GPSCollarSerialNumber)
  view(dat)
  
  #  Remove individuals that don't have collar/location data
  #  White-tailed deer: 002WTD17, 003WTD17, 004WTD17, 005WTD17, 006WTD17... like 13 of them
  which(is.na(wtd_info$GPSCollarSerialNumber))
  length(unique(which(is.na(wtd_info$GPSCollarSerialNumber))))  # number of deer w/o collars
  wtd_info <- wtd_info[!is.na(wtd_info$GPSCollarSerialNumber),]
  #  Remove individuals that died on capture date
  #  019WTD17, 70WTD18
  # wtd_info <- wtd_info[wtd_info$IndividualIdentifier != "019WTD17" &
  #                        wtd_info$IndividualIdentifier != "70WTD18",]
  
  
  
  # #  Add end date to truncate telemetry data
  # #  Choose an end date if there is no mortality- today
  # lastend <- ymd(Sys.Date())
  # #  Assign mortality date as end date if animal died
  # md_info$enddate <- mdy(md_info$EstimatedMortalityDate)
  # elk_info$enddate <- mdy(elk_info$EstimatedMortalityDate)
  # wtd_info$enddate <- mdy(wtd_info$EstimatedMortalityDate)
  # #  Fill in other end dates with chosen date if animal did not die
  # md_info$enddate[is.na(md_info$enddate)] <- lastend
  # elk_info$enddate[is.na(elk_info$enddate)] <- lastend
  # wtd_info$enddate[is.na(wtd_info$enddate)] <- lastend
  # 
  # #  Trouble shoot when "Error in `$<-.data.frame`(`*tmp*`, "ID", value = 1L) : 
  # #  replacement has 1 row, data has 0" with function below.
  # #  Typically means mismatch in number of unique collar & animal IDs
  # # dat1 <- as.data.frame(unique(elk_info$IndividualIdentifier))
  # # dat2 <- as.data.frame(unique(elk_info$GPSCollarSerialNumber))
  # # dat3 <- as.data.frame(unique(elk_tel$CollarID))
  # # id <- cbind(dat1, dat2)
  # dat <- as.data.frame(elk_info$IndividualIdentifier)
  # dat <- cbind(dat, elk_info$GPSCollarSerialNumber)
  # 
  # #  Remove individuals that don't have location data
  # #  Mule deer: 3969MD17, 3958MD17, R89ND18   WHAT ABOUT 89MD18 & 24MD18
  # md_info <- md_info[md_info$IndividualIdentifier != "3969MD17" &
  #                      md_info$IndividualIdentifier != "3958MD17" &
  #                      md_info$IndividualIdentifier != "R89ND18",]
  # #  Elk: 4836ELK20
  # elk_info <- elk_info[elk_info$IndividualIdentifier != "4836ELK20",]

  
  ####  =========================================================
  ####  Connect location data to individual animal IDs  ####
  
  # #  Create empty dataframe to fill iteritively
  # MDclean <- data.frame()
  # #  How many individuals are looped over?
  # nrow(elk_info) #md_info
  # #  Loop over every unique individual animal and...
  # for(i in 1:nrow(elk_info)){ #md_info
  #   #  Take the individual animal ID
  #   mdID <- droplevels(elk_info$IndividualIdentifier[1]) #md_info
  #   #  Take the animal's GPS collar serial number
  #   mdSN <- elk_info$GPSCollarSerialNumber[1] #md_info
  #   #  Buffer capture date to remove locations affected by capture event
  #   #  Suggested to only use data from 2 weeks after the capture data (some papers suggest 1 month)
  #   start <- mdy(elk_info$CaptureDate[1]) + 14 #md_info
  #   #  Exclude locations 1 day before estimated mortality date
  #   end <- elk_info$enddate[i] - 1 #md_info
  # 
  #   #  Subset telemetry data to the specific individual
  #   md <- subset(elk_tel, CollarID == mdSN) #md_tel
  #   #  Add a new column to the telemetry data with the animal's individual ID
  #   md$ID <- mdID
  #   #  truncate telemetry data by new start and end dates for taht individual
  #   mdlive <- subset(md, daytime >= start & daytime <= end)
  # 
  #   #  Append each unique animal's locations to a clean dataframe
  #   MDclean <- rbind(MDclean, mdlive)
  # }
  # 
  # length(unique(md_cap$IndividualIdentifier))
  # length(unique(md_info$IndividualIdentifier))
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
  
  
  IDtelem <- function(info, telem) {
    #  Create empty dataframe to fill iteritively
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
      #  Suggested to only use data from 2 weeks after the capture data (some papers suggest 1 month)
      start <- mdy(info$CaptureDate[i]) + 14
      #  Exclude locations 1 day before estimated mortality date
      end <- info$enddate[i] - 1
      
      #  Subset telemetry data to the specific individual
      collar <- subset(telem, CollarID == SN)
      #  Add a new column to the telemetry data with the animal's individual ID
      collar$ID <- ID
      #  truncate telemetry data by new start and end dates for taht individual
      collarlive <- subset(collar, daytime >= start & daytime <= end)
      
      #  Append each unique animal's locations to a clean dataframe
      clean <- rbind(clean, collarlive)
    }
    
    #  Organize by individual ID and chronological order of locations
    # clean <- clean %>%
    #   arrange(ID, daytime) %>%
    #   #  Filter out aberrant locations
    #   filter(flgLocation != 1) %>%
    #   filter(flgDate != 1) %>%
    #   filter(flgActive != 0)
    
    return(clean)
  }
  
  md_clean <- IDtelem(md_info, md_tel)
  elk_clean <- IDtelem(elk_info, elk_tel)
  wtd_clean <- IDtelem(wtd_info, wtd_tel)  # not working yet
    

  #  Save data!
  write.csv(MDclean, paste0('MDclean ', Sys.Date(), '.csv'))
  
  #  Trouble shooting when I get this annoying error
  #  Error in `$<-.data.frame`(`*tmp*`, "ID", value = 1L) : 
  #  replacement has 1 row, data has 0
  a <- as.data.frame(wtd_info[,c(1,8)])
  # a <- unique(as.data.frame(wtd_info$IndividualIdentifier)) %>%
  #   mutate(a = "a")
  colnames(a) <- c("animal ID", "collar ID")
  # b <- unique(as.data.frame(wtd_info$GPSCollarSerialNumber)) %>%
  #   mutate(b = "b")
  # colnames(b) <- c("ID", "b")
  c <- unique(as.data.frame(wtd_tel$CollarID)) %>%
    mutate(c = "c")
  colnames(c) <- c("collar ID", "c")
  diff <- full_join(a, c, by = "collar ID") 
