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
  
  
  ####  =====  Import 2 data files  ====  ####
  #  Spreadsheet summarizing each individ ID, collar SN, and their start and end date
  MDids <- read.csv("./Data/MDids 20200401.csv")  # mule deer IDs
  head(MDids)
  names(MDids)
  
  #  Spreadsheet downloaded from the WebAPP containing all telemetry locations
  telem <- read.csv("./Data/telem_md_091020.csv")  #  mule deer locations
  #  Convert date and time to usable format
  #  Note mdy_hms() automatically assigns date/time to UTC unless otherwise changed with force_tz()
  telem$daytime <- mdy_hms(telem$ObsDateTimePST) 
  #  Inspect date/time stamps
  head(telem$daytime)  #  Note these are in UTC
  str(telem$daytime) 
  head(telem)
  
  #  All animals need an end date - either date or death or what you specify
  MDids$enddate <- mdy(MDids$EstimatedMortalityDate) #create a new column (enddate), and fill it with the mortality date if there is one
  lastend <- ymd(Sys.Date()) # Choose an end date if there is no mortality. For now make it today
  MDids$enddate[is.na(MDids$enddate)] <- lastend #give an end date to deer than have not died
  #View(MDids) # confirm it worked
  
  ####  ====  Extract the doe data for the first deer  ===  ####
  # test it out for 23066 = 12MD18
  DoeID <- droplevels(MDids$?..IndividualIdentifier[1])
  #DoeID <- droplevels(MDids$IndividualIdentifier[1])  # Not sure where this i.. comes from but DON'T fix this in the orginal csv
  DoeSN <- MDids$GPSCollarSerialNumber[1]
  start <- mdy(MDids$CaptureDate[1]) + 7 # Here you can set what kind of buffer you want on your capture date
    # Only use data from 2 weeks after the capture data (some papers suggest 1 month)
  end <- MDids$enddate[1] - 1 #for movement data, stop it 1 day before estimated mortality date
    # Note enddate is already converted into date format, dont need mdy()
  
  doe <- subset(telem, CollarID == DoeSN) #location data for the doe you want
  doe$ID <- DoeID
  doelive <- subset(doe, daytime >= start & daytime <= end) #the timespan for the data you want
  #View(doelive)
  
  #  Set up a new dateframe to store the clean location data
  MDclean <- doelive
  
  #  Now iterate through the rest of the mule deer
  (n <- nrow(MDids))
  for(i in 2:n){
    DoeID <- droplevels(MDids$?..IndividualIdentifier[i])  
    #DoeID <- droplevels(MDids$IndividualIdentifier[i]) 
    DoeSN <- MDids$GPSCollarSerialNumber[i]
    start <- mdy(MDids$CaptureDate[i]) + 7
    end <- MDids$enddate[i] - 1
    
    doe <- subset(telem, CollarID == DoeSN) 
    doe$ID <- DoeID
    doelive <- subset(doe, daytime >= start & daytime <= end)
    
    #add the data onto the clean dataset - this might be bad practice but it works...
    #MDclean <- rbind(MDclean, doelive)
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
