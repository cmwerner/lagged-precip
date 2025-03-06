## Forever indebted to R wizard Jane Foster for developing this code for our cold air pooling project 2022
#adapted by K Rand 2023


#this document requires the raw files from TMS probes as well as a reference document saved as a csv
# with columns in this order "tms4Id" (the serial code on the sensor)	
#               "plot" (plot number 1-10)
#               "treatment" (control or treatment)	
#               "startdatetime" (in the format mm/dd/yyyy)

## Load libraries from R
library(plyr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(cowplot)

options(stringsAsFactors=F)

## Set the working directory to the address on your local computer
setwd("/Users/klrand/Desktop/R/snowmelt_tms")

## Read in plot location .csv data file
plot.info <- read.csv("/Users/klrand/Desktop/R/snowmelt_tms/snowmelt_tms_reference.csv", header = T)
#make date into numeric values
#plot.info$startdatetime <- strptime(plot.info$startdatetime, "%m/%d/%Y")

## Look at the first few lines
print(head(plot.info))


###----First sensor download----
#it is easiest to bring in data from one field download session at a time. Then tack on the second
# download session to the end by repeating this chunk.
#Make sure the startdatetime in the "snowmelt_tms_refernce.csv" document is correct for this dataset

## List .csv files in the TMS4 directory
csvDirList <- list.files(path="Snowmelt_TMS_08_04_23", 
                         pattern='*.csv', full.names=TRUE)
## Look at the result
print(csvDirList)
## Repeat but just get the filenames with no directory path
csvList <- list.files(path="Snowmelt_TMS_08_04_23", 
                      pattern='*.csv', full.names=F)
## Look at the result
print(csvList)

## Strip out the serial numbers from filenames. Name them "tms4Id", like in plot.char table, 
#the extraction character may change based on file name
tms4Id <- lapply(strsplit(csvList, "_"), '[[', 2) %>% unlist()
## sort the ids
tms4Id <- sort(tms4Id)
## print the result
print(tms4Id)
sites <- sort(unique(plot.info$plot))

## Now loop through each file in the list of csv files from the data directory,
## Read in each file, reformat the data and extract some time variables, then
## combine all the formatted tables into one large combined file...
## NOTE: with loops you can also run each line individually to see if they work before running the 
#       whole loop. This can help isolate issues if it's not working.

for (i in 1:length(csvList)) {
  
  ## Read raw tms4 data file in. Use read.table because you want to specify the delimiter (;, not comma or tab)
  ## Also specify that there is no header information, specify that the separater is ";" with sep = ";"
  ## Specify that the decimal is a comma as in European format with dec = ","
  data_i <- read.table(csvDirList[i], header = F, sep = ";", dec = ",")
  #get rid of the last column which is created with no data in it
  data_i <- data_i[ -c(10) ]
  ## look at the result
  print(head(data_i))
  ## str(data_i) ## Temperature columns 4:6 should be numbers now
  # convert raw soil moisture count to VMC using TMS calibration curve for sand from Tomst website.
  # see TMS calibr tool on https://tomst.com/web/en/systems/tms/software/
  data_i$V7<- ((-0.000000003*(data_i$V7)^2) + 
                              0.000161192*(data_i$V7) - 0.109956505)*100
  ## Assign names to the columns. Column headers from https://tomst.com/web/en/systems/tms/software/
  names(data_i) <- c("row","dateTimeText","timeZone","temperature1", "temperature2",
                     "temperature3","volumetricsoilmoisture","shake","errFlag")
  #temperature values are read in as character, change to numberic
  data_i<- data_i %>% mutate(temperature1= as.numeric(temperature1),
                             temperature2= as.numeric(temperature2),
                             temperature3= as.numeric(temperature3))
  
  ##DateTime field is read in as character and has decimals, we need it to be recognized as a date
  data_i$dateTime <- strptime(data_i$dateTimeText, "%Y.%m.%d %H:%M", tz = "UTC")
  
  # Determine the serial number of tms4 csv (i)
       tms4Id_i <- tms4Id[i]
       data_i$tms4Id <- tms4Id_i
  
  ## Populate table with year, month, day, doy, hour, minute, decimal hour columns
  data_i <- data_i %>% mutate(Date = as.Date(dateTime),
                              year = as.numeric(format(dateTime, '%Y')),
                              month = as.numeric(format(dateTime, '%m')),
                              WOY = as.numeric(format(dateTime, '%W')),
                              DOY = as.numeric(format(dateTime, '%j')),
                              hour = as.numeric(format(dateTime, '%H')),
                              minute = as.numeric(format(dateTime, '%M')),
                              decHour = hour + (minute/60))
  
  ## Now join modified csv table with plot id table info for that plot
  ## Check if location specified in ibKey table matches site name saved in iButton file name
         plot.info.i <- plot.info %>% dplyr::filter(tms4Id == tms4Id_i)
         ploti <- plot.info.i$Plot
         treatmenti <- plot.info.i$Treatment
         Snowfreei <- plot.info.i$Snowfree
  
  ## Now, try to trim the rows that measured temperature before deployment date.
       #data_i <- data_i %>% dplyr::filter(dateTime > plot.info.i$startdatetime)
  
  ## Include the other location information for now
        data_i <- data_i %>% bind_cols(dplyr::select(plot.info.i, Treatment, Plot, Snowfree))
  
  ## Now merge this specific tms4Id file with the others that have been formatted.
        #**adjust the df title to refer just to this chunk of data.
  if (i == 1) {
    dataAllearly <- data_i
    rm(data_i)
  }
  if (exists("dataAllearly") & i > 1) {
    dataAllearly <- dataAllearly %>% bind_rows(data_i)
    rm(data_i, Treatmenti, tms4Id_i)
  }
} 

###----Second sensor download----
#it is easiest to bring in data from one field download session at a time. Then tack on the second
# download session to the end by repeating this chunk.

## List .csv files in the TMS4 directory
csvDirList <- list.files(path="Snowmelt_TMS_09_21_23", 
                         pattern='*.csv', full.names=TRUE)
## Look at the result
print(csvDirList)
## Repeat but just get the filenames with no directory path
csvList <- list.files(path="Snowmelt_TMS_09_21_23", 
                      pattern='*.csv', full.names=F)
## Look at the result
print(csvList)

## Strip out the serial numbers from filenames. Name them "tms4Id", like in plot.char table, 
#the extraction character may change based on file name
tms4Id <- lapply(strsplit(csvList, "_"), '[[', 2) %>% unlist()
## sort the ids
tms4Id <- sort(tms4Id)
## print the result
print(tms4Id)
sites <- sort(unique(plot.info$plot))

## Now loop through each file in the list of csv files from the data directory,
## Read in each file, reformat the data and extract some time variables, then
## combine all the formatted tables into one large combined file...
## NOTE: with loops you can also run each line individually to see if they work before running the 
#       whole loop. This can help isolate issues if it's not working.

for (i in 1:length(csvList)) {
  
  ## Read raw tms4 data file in. Use read.table because you want to specify the delimiter (;, not comma or tab)
  ## Also specify that there is no header information, specify that the separater is ";" with sep = ";"
  ## Specify that the decimal is a comma as in European format with dec = ","
  data_i <- read.table(csvDirList[i], header = F, sep = ";", dec = ",")
  #get rid of the last column which is created with no data in it
  data_i <- data_i[ -c(10) ]
  ## look at the result
  print(head(data_i))
  ## str(data_i) ## Temperature columns 4:6 should be numbers now
  # convert raw soil moisture count to VMC using TMS calibration curve for sand from Tomst website.
  # see TMS calibr tool on https://tomst.com/web/en/systems/tms/software/
  data_i$V7<- ((-0.000000003*(data_i$V7)^2) + 
                 0.000161192*(data_i$V7) - 0.109956505)*100
  ## Assign names to the columns. This information come from https://tomst.com/web/en/systems/tms/software/
  names(data_i) <- c("row","dateTimeText","timeZone","temperature1", "temperature2",
                     "temperature3","volumetricsoilmoisture","shake","errFlag")
  #temperature values are read in as character, change to numberic
  data_i<- data_i %>% mutate(temperature1= as.numeric(temperature1),
                             temperature2= as.numeric(temperature2),
                             temperature3= as.numeric(temperature3))
  ## DateTime field is read in as character and has decimals where we might use dashes
  data_i$dateTime <- strptime(data_i$dateTimeText, "%Y.%m.%d %H:%M", tz = "UTC")
  
  # Determine the serial number of tms4 csv (i)
  tms4Id_i <- tms4Id[i]
  data_i$tms4Id <- tms4Id_i
  
  ## Populate table with year, month, day, doy, hour, minute, decimal hour?
  data_i <- data_i %>% mutate(Date = as.Date(dateTime),
                              year = as.numeric(format(dateTime, '%Y')),
                              month = as.numeric(format(dateTime, '%m')),
                              WOY = as.numeric(format(dateTime, '%W')),
                              DOY = as.numeric(format(dateTime, '%j')),
                              hour = as.numeric(format(dateTime, '%H')),
                              minute = as.numeric(format(dateTime, '%M')),
                              decHour = hour + (minute/60))
  
  ## Now join modified csv table with plot id table info for that plot
  ## Check if location specified in ibKey table matches site name saved in iButton file name
  plot.info.i <- plot.info %>% dplyr::filter(tms4Id == tms4Id_i)
  ploti <- plot.info.i$Plot
  treatmenti <- plot.info.i$Treatment
  Snowfreei <- plot.info.i$Snowfree
  
  ## Now, try to trim the rows that measured temperature before deployment date.
  #data_i <- data_i %>% dplyr::filter(dateTime > plot.info.i$startdatetime)
  
  ## Include the other location information for now
  data_i <- data_i %>% bind_cols(dplyr::select(plot.info.i, Treatment, Plot, Snowfree))
  
  ## Now merge this specific tms4Id file with the others that have been formatted.
  #**adjust the df title to refer just to this chunk of data.
  if (i == 1) {
    dataAlllate <- data_i
    rm(data_i)
  }
  if (exists("dataAlllate") & i > 1) {
    dataAlllate <- dataAlllate %>% bind_rows(data_i)
    rm(data_i, treatmenti, tms4Id_i)
  }
} 

##----merge early and late data downloads and save csv----

#take all the rows from the first and second datasets and make a new df
dataAll<-rbind(dataAllearly,dataAlllate)
#after doing simple QAQC on dataset, remove columns that we will not use for analyses
dataAll<-dataAll[-c(1:3,8:9)]


#save files as a csv
write.csv(dataAll,file='/Users/klrand/Desktop/R/snowmelt_tms/23DOE_SM_TMS.csv', row.names=FALSE)

##----make some graphs!----

#parse the data to just the date range you want to plot
startdatetime<-"2022-02-01" #set the start date for what data you want to view
dataAllparsed<-dataAll %>% dplyr::filter(dateTime > startdatetime)

#find the row in the parsed data table which matches the snowfree date in treatment plots
snowfreetreatment<-which(dataAllparsed$dateTime== "2023-05-16 00:00:00")

#do the same for the control plot
snowfreecontrol<-which(dataAllparsed$dateTime== "2023-05-25 00:00:00")

#make individual plot for each temperature and for soil moisture

soiltemp<-ggplot(dataAllparsed, aes(x=dateTime, group=treatment))+
  stat_smooth(aes(y=temperature1, linetype=treatment))+
  labs(x = "Time",
       y = "Soil Temperature (C)")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreecontrol]),col="darkgray")+ #adds a vertical line at the snowfree date
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreetreatment]), linetype= "dashed", col="darkgray") #adds a vertical line at snowfree date

surfacetemp<-ggplot(dataAllparsed, aes(x=dateTime, group=treatment))+
  stat_smooth(aes(y=temperature2, linetype=treatment))+
  labs(x = "Time", 
       y = "Soil Surface Temperature (C)")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreecontrol]),col="darkgray")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreetreatment]), linetype= "dashed", col="darkgray")

airtemp<-ggplot(dataAllparsed, aes(x=dateTime, group=treatment))+
  stat_smooth(aes(y=temperature3, linetype=treatment))+
  labs(x = "Time", 
       y = "Air Temperature (C)")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreecontrol]),col="darkgray")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreetreatment]), linetype= "dashed", col="darkgray")


moisture<-ggplot(dataAllparsed, aes(x=dateTime, group=treatment))+
  stat_smooth(aes(y=volumetricsoilmoisture, linetype=treatment))+
  labs(x = "Time", 
       y = "Volumetric Soil Moisture (%)")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreecontrol]),col="darkgray")+
  geom_vline(xintercept = as.numeric(dataAllparsed$dateTime[snowfreetreatment]), linetype= "dashed", col="darkgray")

plot_grid(airtemp, surfacetemp, soiltemp, moisture)

#make plot with a line for each temperature and treatments have different line type
ggplot(dataAll, aes(x=dateTime, group=Treatment))+
  stat_smooth(aes(y=temperature1, color= "green", linetype=Treatment))
  stat_smooth(aes(y=temperature2, color= "blue", linetype=Treatment))+
  stat_smooth(aes(y=temperature3, color= "pink", linetype=Treatment))+
  labs(x = "Time", 
       y = "Temperature (C)")



ggplot(dataAll, aes(x=treatment, y=temperature3)) + 
  geom_boxplot()


