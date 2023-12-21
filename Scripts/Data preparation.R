# This code is used to link the sampling locations with environmental data,
# (sea-surface temperature - SST and surface chlorophyll - Chl-a)
# The first part of this code deals with the dataset of predator-prey body masses,
# while the second part deals with the larger dataset of fish body masses and abundances
# Last updated 16/12/2023

#### Loading packages ####
library(openxlsx)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(tidyverse)

## Fish Data -----------------------------------------------------------------------------------
data <- read.xlsx("Data/Inputs/Stomach_Contents_2023_12_16.xlsx",detectDates = TRUE)
# Refining the data #
data$Date <- as.character(data$Date) #making sure dates are in format comparable to environmental data
data <- subset(data,data$Cruise!="JR100") #subsetting to remove JR100 from the data

### extracting SST data ####
#loading in environmental data
data_2006 <- nc_open("Data/Inputs/Environmental data/2006_fish_SST.nc")
data_2008 <- nc_open("Data/Inputs/Environmental data/2008_fish_SST.nc")
data_2009 <- nc_open("Data/Inputs/Environmental data/2009_fish_SST.nc")

#extracting the spatial coordinates for the data (should be the same across years)
lat <- ncvar_get(data_2009,"latitude") 
lon <- ncvar_get(data_2009,"longitude")

#extracting the dates for each measurement
time_2006 <- ncvar_get(data_2006,"time")
time_2006 <- as.POSIXct(time_2006*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2006 <- as.Date(time_2006) #converting to real time

time_2008 <- ncvar_get(data_2008,"time")
time_2008 <- as.POSIXct(time_2008*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2008 <- as.Date(time_2008) #converting to real time

time_2009 <- ncvar_get(data_2009,"time")
time_2009 <- as.POSIXct(time_2009*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2009 <- as.Date(time_2009) #converting to real time

#Extracting temperature data
T_array_2006 <- ncvar_get(data_2006,"thetao")
T_array_2008 <- ncvar_get(data_2008,"thetao")
T_array_2009 <- ncvar_get(data_2009,"thetao")

# identifying the fill values used
fillvalue_T_2006 <- ncatt_get(data_2006, "thetao", "_FillValue")
fillvalue_T_2008 <- ncatt_get(data_2008, "thetao", "_FillValue")
fillvalue_T_2009 <- ncatt_get(data_2009, "thetao", "_FillValue")

#closing the connections
nc_close(data_2006)
nc_close(data_2008)
nc_close(data_2009)

#replacing fill values with NA
T_array_2006[T_array_2006 == fillvalue_T_2006$value] <- NA
T_array_2008[T_array_2008 == fillvalue_T_2008$value] <- NA
T_array_2009[T_array_2009 == fillvalue_T_2009$value] <- NA

#identifying the unique sampling stations
data_info <- data %>%
  distinct(across(c(Cruise,Event,Date,Year,Lat,Lon)))

data_info$SST <- NA #setting up SST column

data_monthly_averages <- data #making a new dataset for the 30-day average SST for each sampling station

daily_sst <- as.data.frame(matrix(nrow=30,ncol=1)) #empty vector to fill with daily SST values for each sampling location

#running a loop to iterate through each sampling location, extract the temperatures for the 30 days preceding
#and average these
for(i in 1:nrow(data_info)){
   if(data_info$Year[i]=="2006"){
    index <- which(data_info$Date[i]==time_2006)
    for(g in 1:30){
      T_slice_2006 <- T_array_2006[, ,(index-30)+g] 
      T_raster_2006 <- raster(t(as.matrix(T_slice_2006)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2006 <- flip(T_raster_2006, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2006, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (data_info$Year[i]=="2008"){
    index <- which(data_info$Date[i]==time_2008)
    for(g in 1:30){
      T_slice_2008 <- T_array_2008[, ,(index-30)+g] 
      T_raster_2008 <- raster(t(as.matrix(T_slice_2008)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2008 <- flip(T_raster_2008, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2008, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (data_info$Year[i]=="2009"){
    index <- which(data_info$Date[i]==time_2009)
    for(g in 1:30){
      T_slice_2009 <- T_array_2009[, ,(index-30)+g] 
      T_raster_2009 <- raster(t(as.matrix(T_slice_2009)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2009 <- flip(T_raster_2009, direction='y')
      
      daily_sst$V1[g] <-  (raster::extract(T_raster_2009, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$SST[i] <-  mean(daily_sst$V1)
    }
  }
}

#linking the average temperature at each location back to the original dataset:
for(i in 1:nrow(data_monthly_averages)){
  index <- which(data_monthly_averages$Event[i]==data_info$Event & data_monthly_averages$Lat[i]==data_info$Lat & data_monthly_averages$Lon[i]==data_info$Lon)
  data_monthly_averages$SST[i] <- data_info$SST[index]
}

### extracting temperature at depth data ####
#loading in environmental data - temperature at depth, with 36 depth slices and a maximum of 1062m (the depth limit for myctophid core range)
data_2006 <- nc_open("Data/Inputs/Environmental data/2006_temperature_at_depth.nc")
data_2008 <- nc_open("Data/Inputs/Environmental data/2008_temperature_at_depth.nc")
data_2009 <- nc_open("Data/Inputs/Environmental data/2009_temperature_at_depth.nc")

#extracting the spatial coordinates for the data (should be the same across years)
lat <- ncvar_get(data_2009,"latitude") 
lon <- ncvar_get(data_2009,"longitude")

#extracting the dates for each measurement
time_2006 <- ncvar_get(data_2006,"time")
time_2006 <- as.POSIXct(time_2006*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2006 <- as.Date(time_2006) #converting to real time

time_2008 <- ncvar_get(data_2008,"time")
time_2008 <- as.POSIXct(time_2008*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2008 <- as.Date(time_2008) #converting to real time

time_2009 <- ncvar_get(data_2009,"time")
time_2009 <- as.POSIXct(time_2009*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2009 <- as.Date(time_2009) #converting to real time

#Extracting temperature data
T_array_2006 <- ncvar_get(data_2006,"thetao")
T_array_2008 <- ncvar_get(data_2008,"thetao")
T_array_2009 <- ncvar_get(data_2009,"thetao")

# identifying the fill values used
fillvalue_T_2006 <- ncatt_get(data_2006, "thetao", "_FillValue")
fillvalue_T_2008 <- ncatt_get(data_2008, "thetao", "_FillValue")
fillvalue_T_2009 <- ncatt_get(data_2009, "thetao", "_FillValue")

#closing the connections
nc_close(data_2006)
nc_close(data_2008)
nc_close(data_2009)

#replacing fill values with NA
T_array_2006[T_array_2006 == fillvalue_T_2006$value] <- NA
T_array_2008[T_array_2008 == fillvalue_T_2008$value] <- NA
T_array_2009[T_array_2009 == fillvalue_T_2009$value] <- NA

#replacing the data_info dataframe with the new dataset with average SSTs
data_info <- data_monthly_averages %>%
  distinct(across(c(Cruise,Event,Date,Year,Lat,Lon)))

data_info$temp_depth <- NA #setting up temperature at depth column

data_monthly_averages <- data_monthly_averages #making a new dataset for the 30-day average temperature for each sampling station

daily_temp <- as.data.frame(matrix(nrow=30,ncol=1)) #empty vector to fill with daily SST values for each sampling location

#running a loop to iterate through each sampling location, extract the temperatures at 1062m for the 30 days preceding
#and average these. 
for(i in 1:nrow(data_info)){
  if (data_info$Year[i]=="2006"){
    index <- which(data_info$Date[i]==time_2006)
    for(g in 1:30){
      T_slice_2006 <- T_array_2006[, ,,(index-30)+g] 
      T_raster_2006_1062m <- raster(t(as.matrix(T_slice_2006[,,36])), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2006_1062m <- flip(T_raster_2006_1062m, direction='y')
      
      daily_temp$V1[g] <-  (raster::extract(T_raster_2006_1062m, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
     data_info$temp_depth[i] <-  mean(daily_temp$V1)
    }
  } else if (data_info$Year[i]=="2008"){
    index <- which(data_info$Date[i]==time_2008)
    for(g in 1:30){
      T_slice_2008 <- T_array_2008[, ,,(index-30)+g] 
      T_raster_2008_1062m <- raster(t(as.matrix(T_slice_2008[,,36])), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2008_1062m <- flip(T_raster_2008_1062m, direction='y')
      
      daily_temp$V1[g] <-  (raster::extract(T_raster_2008_1062m, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$temp_depth[i] <-  mean(daily_temp$V1)
    }
  } else if (data_info$Year[i]=="2009"){
    index <- which(data_info$Date[i]==time_2009)
    for(g in 1:30){
      T_slice_2009 <- T_array_2009[, ,,(index-30)+g] 
      T_raster_2009_1062m <- raster(t(as.matrix(T_slice_2009[,,36])), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                    crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2009_1062m <- flip(T_raster_2009_1062m, direction='y')
      
      daily_temp$V1[g] <-  (raster::extract(T_raster_2009_1062m, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$temp_depth[i] <-  mean(daily_temp$V1)
    }
  }
}

# For three points in JR161 near the south shetlands the shelf is shallower than 1062m,
# and we therefore take the 1062m temperature value as close to this location as possible
missing_data <- data_info[is.na(data_info$temp_depth)==TRUE,]$Event

for(i in missing_data){
  index <- which(data_info$Event==i)
  index2 <- which(data_info$Date[index]==time_2006)
  for(g in 1:30){
    T_slice_2006 <- T_array_2006[, ,,(index2-30)+g] 
    T_raster_2006_1062m <- raster(t(as.matrix(T_slice_2006[,,36])), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                  crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    T_raster_2006_1062m <- flip(T_raster_2006_1062m, direction='y')
    
    if(data_info$Event[index]==118){
      daily_temp$V1[g] <-  (raster::extract(T_raster_2006_1062m, SpatialPoints(cbind(data_info$Lon[index]+0.03,data_info$Lat[index]))))
    } else if (data_info$Event[index]==114){
      daily_temp$V1[g] <-  (raster::extract(T_raster_2006_1062m, SpatialPoints(cbind(data_info$Lon[index],data_info$Lat[index]+0.05))))
    } else if (data_info$Event[index]==106){
      daily_temp$V1[g] <-  (raster::extract(T_raster_2006_1062m, SpatialPoints(cbind(data_info$Lon[index]-0.1,data_info$Lat[index]+0.08))))
    }
    data_info$temp_depth[index] <-  mean(daily_temp$V1)
  }
}

#linking the average temperature at each location back to the original dataset:
for(i in 1:nrow(data_monthly_averages)){
  index <- which(data_monthly_averages$Event[i]==data_info$Event & data_monthly_averages$Lat[i]==data_info$Lat & data_monthly_averages$Lon[i]==data_info$Lon)
  data_monthly_averages$temp_depth[i] <- data_info$temp_depth[index]
}

### Extracting CHL data ####
#extracting the original fine-scale CHL data:
data_2006 <- nc_open("Data/Inputs/Environmental data/2006_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc")
data_2008 <- nc_open("Data/Inputs/Environmental data/2008_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc")
data_2009 <- nc_open("Data/Inputs/Environmental data/2009_cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D.nc")

#extracting the spatial coordinates for the data (should be the same across years)
lat <- ncvar_get(data_2009,"lat")
lon <- ncvar_get(data_2009,"lon")

#extracting the dates for each measurement
time_2006 <- ncvar_get(data_2006,"time")
time_2006 <- as.POSIXct(time_2006*24*60*60,origin='1900-01-01 00:00:00') #converting to real time
time_2006 <- as.Date(time_2006) #converting to date

time_2008 <- ncvar_get(data_2008,"time")
time_2008 <- as.POSIXct(time_2008*24*60*60,origin='1900-01-01 00:00:00') #converting to real time
time_2008 <- as.Date(time_2008) #converting to date

time_2009 <- ncvar_get(data_2009,"time")
time_2009 <- as.POSIXct(time_2009*24*60*60,origin='1900-01-01 00:00:00') #converting to real time
time_2009 <- as.Date(time_2009) #converting to date

#Extracting Chl-a data
CHL_array_2006 <- ncvar_get(data_2006,"CHL")
CHL_array_2008 <- ncvar_get(data_2008,"CHL")
CHL_array_2009 <- ncvar_get(data_2009,"CHL")

# identifying the fill values used
fillvalue_CHL_2006 <- ncatt_get(data_2006, "CHL", "_FillValue")
fillvalue_CHL_2008 <- ncatt_get(data_2008, "CHL", "_FillValue")
fillvalue_CHL_2009 <- ncatt_get(data_2009, "CHL", "_FillValue")

#closing the connections
nc_close(data_2006)
nc_close(data_2008)
nc_close(data_2009)

#replacing fill values with NA
CHL_array_2006[CHL_array_2006 == fillvalue_CHL_2006$value] <- NA
CHL_array_2008[CHL_array_2008 == fillvalue_CHL_2008$value] <- NA
CHL_array_2009[CHL_array_2009 == fillvalue_CHL_2009$value] <- NA

#replacing the data_info dataframe with the new dataset with average SSTs
data_info <- data_monthly_averages %>%
  distinct(across(c(Cruise,Event,Date,Year,Lat,Lon)))

data_info$CHL <- NA #setting up empty column for CHL

daily_chl <- as.data.frame(matrix(nrow=30,ncol=1))

#running a loop to iterate through each sampling location, extract the Chl-a for the 30 days preceding
#and average these
for(i in 1:nrow(data_info)){
   if (data_info$Year[i]=="2006"){
    index <- which(data_info$Date[i]==time_2006)
    for(g in 1:30){
      CHL_slice_2006 <- CHL_array_2006[, ,(index-30)+g] 
      CHL_raster_2006 <- raster(t(as.matrix(CHL_slice_2006)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      CHL_raster_2006 <- flip(CHL_raster_2006, direction='y')
      daily_chl$V1[g] <-  (raster::extract(CHL_raster_2006, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$CHL[i] <-  mean(daily_chl$V1)
    }
  } else if (data_info$Year[i]=="2008"){
    index <- which(data_info$Date[i]==time_2008)
    for(g in 1:30){
      CHL_slice_2008 <- CHL_array_2008[, ,(index-30)+g] 
      CHL_raster_2008 <- raster(t(as.matrix(CHL_slice_2008)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      CHL_raster_2008 <- flip(CHL_raster_2008, direction='y')
      daily_chl$V1[g] <-  (raster::extract(CHL_raster_2008, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$CHL[i] <-  mean(daily_chl$V1)
    }
  } else if (data_info$Year[i]=="2009"){
    index <- which(data_info$Date[i]==time_2009)
    for(g in 1:30){
      CHL_slice_2009 <- CHL_array_2009[, ,(index-30)+g] 
      CHL_raster_2009 <- raster(t(as.matrix(CHL_slice_2009)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                                crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      CHL_raster_2009 <- flip(CHL_raster_2009, direction='y')
      
      daily_chl$V1[g] <-  (raster::extract(CHL_raster_2009, SpatialPoints(cbind(data_info$Lon[i],data_info$Lat[i]))))
      data_info$CHL[i] <-  mean(daily_chl$V1)
    }
  }
}

#linking the average chl-a at each location back to the original dataset:
for(i in 1:nrow(data_monthly_averages)){
  index <- which(data_monthly_averages$Event[i]==data_info$Event & data_monthly_averages$Lat[i]==data_info$Lat & data_monthly_averages$Lon[i]==data_info$Lon)
  data_monthly_averages$CHL[i] <- data_info$CHL[index]
}

### Now subsetting species and refining ####
table(data_monthly_averages$Species.Code)
#subsetting to exclude two species with very few records covering a small latitudinal range:
data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$Species.Code!="GYO"& data_monthly_averages$Species.Code!="GYP") 

# Also removing stomach records which are empty, unidentifiable or fish.
# Fish removed to avoid potential net-feeding and also to remove any feedback effects of changes in fish size being reflected in prey sizes
# Also removing records for which no weight was recorded
nrow(data_monthly_averages) #3834

data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$Prey!="Jelly mass")
data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$Prey!="Fish indet")
data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$Prey!="Empty")
data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$Prey!="Unidentified")

data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$`Prey.Mass.(g)`!="NA")
data_monthly_averages <- subset(data_monthly_averages,data_monthly_averages$`Prey.Mass.(g)`!="0")

nrow(data_monthly_averages) #removed 127 prey records

# updating the unique ID column to make sure it has a unique value for each individual fish
data_monthly_averages$ID_unique <- paste(data_monthly_averages$Cruise,data_monthly_averages$Event,data_monthly_averages$Species.Code,data_monthly_averages$ID,sep="_")

#generating a new column with log10 body mass values for each fish
data_monthly_averages$fish_mass_log10 <- log10(data_monthly_averages$weight_discovery_g)

#adding in grouped locations to match to zooplankton data for preference analysis
#JR177 events 357, 378 and 379 do not have associated plankton data so will be excluded
data_monthly_averages$location_group <- NA

data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==42|data_monthly_averages$Event==43|
                                                                                      data_monthly_averages$Event==56|data_monthly_averages$Event==59))] <- 1
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==73|data_monthly_averages$Event==84|
                                                                                      data_monthly_averages$Event==91))] <- 2
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==106|data_monthly_averages$Event==114|
                                                                                      data_monthly_averages$Event==118))] <- 3
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==134|data_monthly_averages$Event==136|
                                                                                      data_monthly_averages$Event==142))] <- 4
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==157|data_monthly_averages$Event==159|
                                                                                      data_monthly_averages$Event==160))] <- 5
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==199|data_monthly_averages$Event==201|
                                                                                      data_monthly_averages$Event==214|data_monthly_averages$Event==217|
                                                                                      data_monthly_averages$Event==218))] <- 6
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==253|data_monthly_averages$Event==265|
                                                                                      data_monthly_averages$Event==267))] <- 7
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR161" & (data_monthly_averages$Event==273|data_monthly_averages$Event==275|
                                                                                      data_monthly_averages$Event==283))] <- 8

data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==75|data_monthly_averages$Event==78))] <- 9
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==108|data_monthly_averages$Event==123|
                                                                                      data_monthly_averages$Event==124))] <- 10
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==158|data_monthly_averages$Event==161|
                                                                                      data_monthly_averages$Event==165))] <- 11
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==198|data_monthly_averages$Event==199|
                                                                                      data_monthly_averages$Event==205))] <- 12
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==254|data_monthly_averages$Event==255))] <- 13
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==295|data_monthly_averages$Event==300|
                                                                                      data_monthly_averages$Event==301|data_monthly_averages$Event==305))] <- 14
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR177" & (data_monthly_averages$Event==328|data_monthly_averages$Event==329|
                                                                                      data_monthly_averages$Event==334|data_monthly_averages$Event==335))] <- 15

data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==17))] <- 16
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==42))] <- 17
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==55))] <- 18
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==81|data_monthly_averages$Event==100|
                                                                                      data_monthly_averages$Event==101))] <- 19
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==115|data_monthly_averages$Event==127|
                                                                                      data_monthly_averages$Event==128))] <- 20
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==141|data_monthly_averages$Event==142))] <- 21
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==185))] <- 22
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==225|data_monthly_averages$Event==226|
                                                                                      data_monthly_averages$Event==227))] <- 23
data_monthly_averages$location_group[which(data_monthly_averages$Cruise=="JR200" & (data_monthly_averages$Event==235|data_monthly_averages$Event==236|
                                                                                      data_monthly_averages$Event==237))] <- 24
write.csv(data_monthly_averages,"Data/Processed/FINAL_FISH_DATA.csv")


## Zooplankton data -------------------------------------------------------------

# First I need to tidy up the data. Starting with the macrozooplankton before moving on to the
# mesozooplankton

# refining the macrozooplankton data #
macro_jr161 <- read.csv("Data/Inputs/Zooplankton data/Macro_JR161.csv")
macro_jr177 <- read.csv("Data/Inputs/Zooplankton data/Macro_JR177.csv")
macro_jr200 <- read.csv("Data/Inputs/Zooplankton data/Macro_JR200.csv")

#combining these three datasets
macro_data <- rbind(macro_jr161,macro_jr177,macro_jr200)

macro_data$Date <- format(as.Date(macro_data$Date, format = "%d/%m/%Y"), "%Y-%m-%d")

#removing fish and other misc prey
macro_data <- macro_data[-(which(macro_data$Group=="Fish"|macro_data$Group=="Unknown"|macro_data$Group=="egg")),]

macro_data <- macro_data[-(which(macro_data$Taxa=="Siphonophora (Nectophores)"|macro_data$Taxa=="Siphonophora nectophores"|
                                   macro_data$Taxa=="Siphonophore bracts and nectophores"|macro_data$Taxa=="Siphonophore gonads"|
                                   macro_data$Taxa=="Siphonophore top of colony")),]

# Fixing spelling mistakes, as these otherwise prevent correct averaging across nets in some cases
macro_data$Taxa[which(macro_data$Taxa=="Acanthephyra")] <- "Acanthephyra sp."
macro_data$Taxa[which(macro_data$Taxa=="Acanthephyrasp.")] <- "Acanthephyra sp."
macro_data$Taxa[which(macro_data$Taxa=="Acathephyra sp.")] <- "Acanthephyra sp."
macro_data$Taxa[which(macro_data$Taxa=="Allurotheuthis antarcticus")] <- "Alluroteuthis antarcticus"
macro_data$Taxa[c(which(macro_data$Taxa=="Amphipoda sp."))] <- "Amphipod sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Amphipoda sp. "))] <- "Amphipod sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Amphopoda sp."))] <- "Amphipod sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Atolla sp. "))] <- "Atolla sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Atollidae"))] <- "Atolla sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Brown clear jellyfish (unknown)"))] <- "Brown jellyfish"
macro_data$Taxa[c(which(macro_data$Taxa=="Brown jellyfish (unknown)"))] <- "Brown jellyfish"
macro_data$Taxa[c(which(macro_data$Taxa=="Clear jellyfish (brown fringe) unknown "))] <- "Clear jellyfish (brown fringe) unknown"
macro_data$Taxa[c(which(macro_data$Taxa=="Clear Jellyfish (unknown)"))] <- "Clear jellyfish (unknown)"
macro_data$Taxa[c(which(macro_data$Taxa=="Clio pyrimidata"))] <- "Clio pyramidata"
macro_data$Taxa[c(which(macro_data$Taxa=="Clio pyrimidata forma sulcata"))] <- "Clio pyramidata"
macro_data$Taxa[c(which(macro_data$Taxa=="Clio recurva or pyramidata"))] <- "Clio sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Clio sp. "))] <- "Clio sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Cyllopus sp. 1"))] <- "Cyllopus sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Cyphocaris faueri"))] <- "Cyphocaris faurei"
macro_data$Taxa[c(which(macro_data$Taxa=="Decapod crustacea"))] <- "Decapoda sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Diphyes sp. 1"))] <- "Diphyes sp.1"
macro_data$Taxa[c(which(macro_data$Taxa=="Diphyes sp. 2"))] <- "Diphyes sp.2"
macro_data$Taxa[c(which(macro_data$Taxa=="Diphyes sp. 3"))] <- "Diphyes sp.3"
macro_data$Taxa[c(which(macro_data$Taxa=="Euphausidae"))] <- "Euphausiidae"
macro_data$Taxa[c(which(macro_data$Taxa=="Euphausidae spp."))] <- "Euphausiidae"
macro_data$Taxa[c(which(macro_data$Taxa=="Krill"))] <- "Euphausiidae"
macro_data$Taxa[c(which(macro_data$Taxa=="Hippopodiidae sp."))] <- "Hippopodiidae"
macro_data$Taxa[c(which(macro_data$Taxa=="Hydromedusae (indet.)"))] <- "Hydromedusae indet."
macro_data$Taxa[c(which(macro_data$Taxa=="Hyydromedusae (large)"))] <- "Hydromedusae indet."
macro_data$Taxa[c(which(macro_data$Taxa=="Paradania boecki"))] <- "Parandania boecki"
macro_data$Taxa[c(which(macro_data$Taxa=="Galiteuthis glacialis "))] <- "Galiteuthis glacialis"
macro_data$Taxa[c(which(macro_data$Taxa=="Medusea (red)"))] <- "Medusae (red)"
macro_data$Taxa[c(which(macro_data$Taxa=="Pasiphaea"))] <- "Pasiphaea sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Pasiphae sp."))] <- "Pasiphaea sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Pasiphae scotiae"))] <- "Pasiphaea scotiae"
macro_data$Taxa[c(which(macro_data$Taxa=="Pasiphaeea sp."))] <- "Pasiphaea sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Periphylla small"))] <- "Periphylla periphylla"
macro_data$Taxa[c(which(macro_data$Taxa=="Salpa thompsoni (colonial)"))] <- "Salpa thompsoni"
macro_data$Taxa[c(which(macro_data$Taxa=="Salpa thompsoni (solitary)"))] <- "Salpa thompsoni"
macro_data$Taxa[c(which(macro_data$Taxa=="Scyphomedusae (Pete Ward)"))] <- "Scyphomedusae"
macro_data$Taxa[c(which(macro_data$Taxa=="Scyphomedusae (small)"))] <- "Scyphomedusae"
macro_data$Taxa[c(which(macro_data$Taxa=="Scyphomedusae indet."))] <- "Scyphomedusae"
macro_data$Taxa[c(which(macro_data$Taxa=="Sergia sp. "))] <- "Sergia sp."
macro_data$Taxa[c(which(macro_data$Taxa=="Sibogita"))] <- "Sibogita sp."

#averaging the abundance of each taxa across nets if multiple records present
macro_data <- macro_data%>%
  group_by(Cruise,Event,Date,Lat,Lon,Group,Taxa)%>%
  dplyr::summarise(mean_wet_mass=mean(mean.Individual.wet.mass..g.),mean_density_m2=mean(total.density..ind.m.2.),mean_density_m3=mean(total.density..ind.m.3.))

macro_data$log10_mass <- log10(macro_data$mean_wet_mass) 
macro_data$biomass_m2 <- macro_data$mean_wet_mass*macro_data$mean_density_m2

# refining the mesozooplankton data #
meso_jr161 <- read.csv("Data/Inputs/Zooplankton data/Meso_JR161.csv")
meso_jr177 <- read.csv("Data/Inputs/Zooplankton data/Meso_JR177.csv")
meso_jr200 <- read.csv("Data/Inputs/Zooplankton data/Meso_JR200.csv")

# combining the datasets into a single dataframe the raw abundance data in long format
meso_161_long <- meso_jr161 %>%
  pivot_longer("E110":"E181")

meso_177_long <- meso_jr177 %>%
  pivot_longer("E53":"E203")

meso_200_long <- meso_jr200 %>%
  pivot_longer("E51":"E240")

meso_161_long$Year <- 2006
meso_177_long$Year <- 2008
meso_200_long$Year <- 2009

meso_161_long$Lat <- NA
meso_161_long$Lon <- NA
meso_177_long$Lat <- NA
meso_177_long$Lon <- NA
meso_200_long$Lat <- NA
meso_200_long$Lon <- NA

meso_161_long$Date <- NA
meso_177_long$Date <- NA
meso_200_long$Date <- NA

names(meso_161_long)[11] <- "Event"
names(meso_177_long)[11] <- "Event"
names(meso_200_long)[11] <- "Event"

# Manually matching the coordinates and date of each sampling Event, taken from the relevant cruise reports:
for(i in 1:nrow(meso_161_long)){
  if(meso_161_long$Event[i]=="E28"){
    meso_161_long$Lat[i] <- -57.74695
    meso_161_long$Lon[i] <- -50.38033
    meso_161_long$Date[i] <- "2006-10-27"
  } else if(meso_161_long$Event[i]=="E54"){
    meso_161_long$Lat[i] <- -57.74303
    meso_161_long$Lon[i] <- -50.43559
    meso_161_long$Date[i] <- "2006-10-29"
  } else if(meso_161_long$Event[i]=="E87"){
    meso_161_long$Lat[i] <- -60.64752
    meso_161_long$Lon[i] <- -48.70069
    meso_161_long$Date[i] <- "2006-11-02"
  } else if(meso_161_long$Event[i]=="E110"){
    meso_161_long$Lat[i] <- -60.43905
    meso_161_long$Lon[i] <- -44.61749
    meso_161_long$Date[i] <- "2006-11-06"
  } else if(meso_161_long$Event[i]=="E138"){
    meso_161_long$Lat[i] <- -59.67942
    meso_161_long$Lon[i] <- -44.05925
    meso_161_long$Date[i] <- "2006-11-09"
  } else if(meso_161_long$Event[i]=="E162"){
    meso_161_long$Lat[i] <- -57.43815
    meso_161_long$Lon[i] <- -42.61729
    meso_161_long$Date[i] <- "2006-11-17"
  } else if(meso_161_long$Event[i]=="E181"){
    meso_161_long$Lat[i] <- -55.20924
    meso_161_long$Lon[i] <- -41.24121
    meso_161_long$Date[i] <- "2006-11-18"
  } else if(meso_161_long$Event[i]=="E260"){
    meso_161_long$Lat[i] <- -52.86457
    meso_161_long$Lon[i] <- -40.0895
    meso_161_long$Date[i] <- "2006-11-25"
  } else if(meso_161_long$Event[i]=="E287"){
    meso_161_long$Lat[i] <- -50.00168
    meso_161_long$Lon[i] <- -38.00003
    meso_161_long$Date[i] <- "2006-11-29"
  }
}

for(i in 1:nrow(meso_177_long)){
  if(meso_177_long$Event[i]=="E53"){
    meso_177_long$Lat[i] <- -60.4982
    meso_177_long$Lon[i] <- -48.1924
    meso_177_long$Date[i] <- "2008-01-04"
  } else if(meso_177_long$Event[i]=="E91"){
    meso_177_long$Lat[i] <- -60.2088
    meso_177_long$Lon[i] <- -44.4183
    meso_177_long$Date[i] <- "2008-01-08"
  } else if(meso_177_long$Event[i]=="E149"){
    meso_177_long$Lat[i] <- -59.6888
    meso_177_long$Lon[i] <- -44.0544
    meso_177_long$Date[i] <- "2008-01-15"
  } else if(meso_177_long$Event[i]=="E188"){
    meso_177_long$Lat[i] <- -58.0228
    meso_177_long$Lon[i] <- -42.9809
    meso_177_long$Date[i] <- "2008-01-19"
  } else if(meso_177_long$Event[i]=="E203"){
    meso_177_long$Lat[i] <- -58.0271
    meso_177_long$Lon[i] <- -42.9719
    meso_177_long$Date[i] <- "2008-01-20"
  } else if(meso_177_long$Event[i]=="E219"){
    meso_177_long$Lat[i] <- -57.1408
    meso_177_long$Lon[i] <- -42.4327
    meso_177_long$Date[i] <- "2008-01-22"
  } else if(meso_177_long$Event[i]=="E230"){
    meso_177_long$Lat[i] <- -55.2136
    meso_177_long$Lon[i] <- -41.244
    meso_177_long$Date[i] <- "2008-01-25"
  } else if(meso_177_long$Event[i]=="E283"){
    meso_177_long$Lat[i] <- -52.8589
    meso_177_long$Lon[i] <- -40.0969
    meso_177_long$Date[i] <- "2008-02-01"
  } else if(meso_177_long$Event[i]=="E321"){
    meso_177_long$Lat[i] <- -52.6272
    meso_177_long$Lon[i] <- -39.1152
    meso_177_long$Date[i] <- "2008-02-04"
  } else if(meso_177_long$Event[i]=="E371"){
    meso_177_long$Lat[i] <- -53.715
    meso_177_long$Lon[i] <- -37.9645
    meso_177_long$Date[i] <- "2008-02-09"
  }
}

for(i in 1:nrow(meso_200_long)){
  if(meso_200_long$Event[i]=="E25"){
    meso_200_long$Lat[i] <- -60.4982
    meso_200_long$Lon[i] <- -48.19196
    meso_200_long$Date[i] <- "2009-03-16"
  } else if(meso_200_long$Event[i]=="E51"){
    meso_200_long$Lat[i] <- -60.20816
    meso_200_long$Lon[i] <- -44.40802
    meso_200_long$Date[i] <- "2009-03-19"
  } else if(meso_200_long$Event[i]=="E59"){
    meso_200_long$Lat[i] <- -59.68792
    meso_200_long$Lon[i] <- -44.05423
    meso_200_long$Date[i] <- "2009-03-20"
  } else if(meso_200_long$Event[i]=="E85"){
    meso_200_long$Lat[i] <- -58.03071
    meso_200_long$Lon[i] <- -42.97145
    meso_200_long$Date[i] <- "2009-03-23"
  } else if(meso_200_long$Event[i]=="E123"){
    meso_200_long$Lat[i] <- -56.76345
    meso_200_long$Lon[i] <- -42.21855
    meso_200_long$Date[i] <- "2009-03-27"
  } else if(meso_200_long$Event[i]=="E138"){
    meso_200_long$Lat[i] <- -55.25906
    meso_200_long$Lon[i] <- -41.3581
    meso_200_long$Date[i] <- "2009-03-29"
  } else if(meso_200_long$Event[i]=="E179"){
    meso_200_long$Lat[i] <- -52.71758
    meso_200_long$Lon[i] <- -40.15279
    meso_200_long$Date[i] <- "2009-04-03"
  } else if(meso_200_long$Event[i]=="E190"){
    meso_200_long$Lat[i] <- -52.85069
    meso_200_long$Lon[i] <- -37.09819
    meso_200_long$Date[i] <- "2009-04-06"
  } else if(meso_200_long$Event[i]=="E230"){
    meso_200_long$Lat[i] <- -50.05917
    meso_200_long$Lon[i] <- -33.8717
    meso_200_long$Date[i] <- "2009-04-10"
  } else if(meso_200_long$Event[i]=="E240"){
    meso_200_long$Lat[i] <- -50.59804
    meso_200_long$Lon[i] <- -33.80711
    meso_200_long$Date[i] <- "2009-04-11"
  }
}

meso_161_long$Cruise <- "JR161"
meso_177_long$Cruise <- "JR177"
meso_200_long$Cruise <- "JR200"

#combining the datasets:
meso_data_long <- rbind(meso_161_long,meso_177_long,meso_200_long)

meso_data_long$Event <- as.integer(sub('.', '', meso_data_long$Event)) #removing the E from the event names

### Matching SST data to events for both macrozooplankton and mesozooplankton ####
data_2006 <- nc_open("Data/Inputs/Environmental data/2006_zoo_SST.nc")
data_2008 <- nc_open("Data/Inputs/Environmental data/2008_zoo_SST.nc")
data_2009 <- nc_open("Data/Inputs/Environmental data/2009_zoo_SST.nc")

#extracting the spatial coordinates for the data (should be the same across years)
lat <- ncvar_get(data_2009,"latitude") 
lon <- ncvar_get(data_2009,"longitude")

#extracting the dates for each measurement
time_2006 <- ncvar_get(data_2006,"time")
time_2006 <- as.POSIXct(time_2006*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2006 <- as.Date(time_2006) #converting to real time

time_2008 <- ncvar_get(data_2008,"time")
time_2008 <- as.POSIXct(time_2008*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2008 <- as.Date(time_2008) #converting to real time

time_2009 <- ncvar_get(data_2009,"time")
time_2009 <- as.POSIXct(time_2009*3600,origin='1950-01-01 00:00:00') #converting to real time
time_2009 <- as.Date(time_2009) #converting to real time

#Extracting temperature data
T_array_2006 <- ncvar_get(data_2006,"thetao")
T_array_2008 <- ncvar_get(data_2008,"thetao")
T_array_2009 <- ncvar_get(data_2009,"thetao")

# identifying the fill values used
fillvalue_T_2006 <- ncatt_get(data_2006, "thetao", "_FillValue")
fillvalue_T_2008 <- ncatt_get(data_2008, "thetao", "_FillValue")
fillvalue_T_2009 <- ncatt_get(data_2009, "thetao", "_FillValue")

#closing the connections
nc_close(data_2006)
nc_close(data_2008)
nc_close(data_2009)

#replacing fill values with NA
T_array_2006[T_array_2006 == fillvalue_T_2006$value] <- NA
T_array_2008[T_array_2008 == fillvalue_T_2008$value] <- NA
T_array_2009[T_array_2009 == fillvalue_T_2009$value] <- NA

# FOR MACROZOOPLANKTON

#creating a new dataframe to populate with SST values for haul locations,
#to then add to the main dataset
plankton_info <- macro_data %>%
  distinct(across(c(Cruise,Event,Date,Lat,Lon)))

macro_data$SST <- NA
plankton_info$SST <- NA
daily_sst <- as.data.frame(matrix(nrow=30,ncol=1))

for(i in 1:nrow(plankton_info)){
  if (plankton_info$Cruise[i]=="JR161"){
    index <- which(plankton_info$Date[i]==time_2006)
    for(g in 1:30){
      T_slice_2006 <- T_array_2006[, ,(index-30)+g] 
      T_raster_2006 <- raster(t(as.matrix(T_slice_2006)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2006 <- flip(T_raster_2006, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2006, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (plankton_info$Cruise[i]=="JR177"){
    index <- which(plankton_info$Date[i]==time_2008)
    for(g in 1:30){
      T_slice_2008 <- T_array_2008[, ,(index-30)+g] 
      T_raster_2008 <- raster(t(as.matrix(T_slice_2008)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2008 <- flip(T_raster_2008, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2008, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (plankton_info$Cruise[i]=="JR200"){
    index <- which(plankton_info$Date[i]==time_2009)
    for(g in 1:30){
      T_slice_2009 <- T_array_2009[, ,(index-30)+g] 
      T_raster_2009 <- raster(t(as.matrix(T_slice_2009)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2009 <- flip(T_raster_2009, direction='y')
      
      daily_sst$V1[g] <-  (raster::extract(T_raster_2009, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  }
}

#linking back to original dataset:
for(i in 1:nrow(macro_data)){
  index <- which(macro_data$Event[i]==plankton_info$Event & macro_data$Lat[i]==plankton_info$Lat & macro_data$Lon[i]==plankton_info$Lon)
  macro_data$SST[i] <- plankton_info$SST[index]
}

# FOR MESOZOOPLANKTON
#creating a new dataframe to populate with SST values for haul locations,
#to then add to the main dataset

plankton_info <- meso_data_long %>%
  distinct(across(c(Cruise,Event,Date,Lat,Lon)))

meso_data_long$SST <- NA
plankton_info$SST <- NA
daily_sst <- as.data.frame(matrix(nrow=30,ncol=1))

for(i in 1:nrow(plankton_info)){
  if (plankton_info$Cruise[i]=="JR161"){
    index <- which(plankton_info$Date[i]==time_2006)
    for(g in 1:30){
      T_slice_2006 <- T_array_2006[, ,(index-30)+g] 
      T_raster_2006 <- raster(t(as.matrix(T_slice_2006)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2006 <- flip(T_raster_2006, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2006, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (plankton_info$Cruise[i]=="JR177"){
    index <- which(plankton_info$Date[i]==time_2008)
    for(g in 1:30){
      T_slice_2008 <- T_array_2008[, ,(index-30)+g] 
      T_raster_2008 <- raster(t(as.matrix(T_slice_2008)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2008 <- flip(T_raster_2008, direction='y')
      daily_sst$V1[g] <-  (raster::extract(T_raster_2008, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  } else if (plankton_info$Cruise[i]=="JR200"){
    index <- which(plankton_info$Date[i]==time_2009)
    for(g in 1:30){
      T_slice_2009 <- T_array_2009[, ,(index-30)+g] 
      T_raster_2009 <- raster(t(as.matrix(T_slice_2009)), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat),
                              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
      T_raster_2009 <- flip(T_raster_2009, direction='y')
      
      daily_sst$V1[g] <-  (raster::extract(T_raster_2009, SpatialPoints(cbind(plankton_info$Lon[i],plankton_info$Lat[i]))))
      plankton_info$SST[i] <-  mean(daily_sst$V1)
    }
  }
}

## matching the mesozooplankton sites with SSTs ##
for(i in 1:nrow(meso_data_long)){
  index <- which(meso_data_long$Event[i]==plankton_info$Event & meso_data_long$Lat[i]==plankton_info$Lat & meso_data_long$Lon[i]==plankton_info$Lon)
  meso_data_long$SST[i] <- plankton_info$SST[index]
}

### Combining the macro and meso data into a single dataset ####
macro_data$origin <- "macrozooplankton"
meso_data_long$origin <- "mesozooplankton"

macro_refined <- subset(macro_data,select=c(Taxa,Cruise,Date,Event,Group,mean_wet_mass,log10_mass,mean_density_m2,SST,Lat,Lon,origin))
colnames(macro_refined) <- c("Taxa","Cruise","Date","Event","Group","mean_wet_mass","log10_mass","density_m2", "SST","Lat","Lon","origin")

meso_data_long$log10_mass <- log10(meso_data_long$WM_g)
meso_refined <- subset(meso_data_long,select=c(Taxa,Cruise,Date,Event,WM_g,log10_mass,value,Category,SST,Lat,Lon,origin))
colnames(meso_refined) <- c("Taxa","Cruise","Date","Event","mean_wet_mass","log10_mass","density_m2", "Group", "SST","Lat","Lon","origin")

combined_data <- rbind(macro_refined,meso_refined)
combined_data <- combined_data[-which(is.na(combined_data$density_m2)==TRUE),] #removing missing densities

### adding a location group identifier to the zooplankton hauls to match to the fish data for selectivity analysis ####
combined_data$location_group <- NA
#NOTES:
#2008 event 219 only has a mesozooplankton sample, so will leave this blank
#2009 event 190 only has mesozooplankton sample
#2008 event 371 only has mesozooplankton sample

combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==28|combined_data$Event==42|combined_data$Event==43|combined_data$Event==56))] <- 1
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==73|combined_data$Event==84|combined_data$Event==87))] <- 2
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==106|combined_data$Event==110|combined_data$Event==118))] <- 3
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==134|combined_data$Event==136|combined_data$Event==138))] <- 4
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==157|combined_data$Event==159|combined_data$Event==162))] <- 5
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==181|combined_data$Event==199|combined_data$Event==214))] <- 6
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==253|combined_data$Event==260))] <- 7
combined_data$location_group[which(combined_data$Cruise=="JR161" & (combined_data$Event==273|combined_data$Event==275|combined_data$Event==287))] <- 8

combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==53|combined_data$Event==74|combined_data$Event==75|combined_data$Event==78))] <- 9
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==91|combined_data$Event==102|combined_data$Event==123|combined_data$Event==124))] <- 10
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==149|combined_data$Event==158|combined_data$Event==161))] <- 11
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==188|combined_data$Event==198|combined_data$Event==199|combined_data$Event==203))] <- 12
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==230|combined_data$Event==254|combined_data$Event==255))] <- 13
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==283|combined_data$Event==295|combined_data$Event==305))] <- 14
combined_data$location_group[which(combined_data$Cruise=="JR177" & (combined_data$Event==321|combined_data$Event==328|combined_data$Event==329))] <- 15

combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==17|combined_data$Event==18|combined_data$Event==25))] <- 16
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==42|combined_data$Event==43|combined_data$Event==51))] <- 17
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==55|combined_data$Event==56|combined_data$Event==59))] <- 18
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==82|combined_data$Event==85|combined_data$Event==100))] <- 19
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==115|combined_data$Event==123|combined_data$Event==127))] <- 20
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==138|combined_data$Event==141|combined_data$Event==142))] <- 21
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==179|combined_data$Event==184|combined_data$Event==185))] <- 22
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==225|combined_data$Event==226|combined_data$Event==230))] <- 23
combined_data$location_group[which(combined_data$Cruise=="JR200" & (combined_data$Event==235|combined_data$Event==236|combined_data$Event==240))] <- 24

#now removing events with no location assigned
combined_data <- na.omit(combined_data) 

# I originally averaged the body sizes of each taxon in the macrozooplankton across locations,
# I think to make the data more comparable to the mesozooplankton. However I'm not convinced this was the right thing to do.
# I will try the analyses without doing that, and see what the result is.

# Because the mesozooplankton body sizes are fixed averages for each taxon, I need to calculate an average density
# at each location. In fact, only one location group (12) has more than one mesozooplankton sample: JR177 events 188 and 203.
# I need to combine these to ensure I have comparable data across locations
# 
# data_macro <- subset(combined_data,combined_data$origin=="macrozooplankton")
# 
# data_macro_means <- data_macro %>%
#   dplyr::group_by(Cruise,Taxa,location_group) %>%
#   summarise(log10_mass=mean(log10_mass),mean_mass_g = mean(mean_wet_mass ),density_m2=mean(density_m2),SST=mean(SST))

#subsetting the data to extract events 188 and 203 for averaging
data_meso_subset <- subset(combined_data,combined_data$Event==188 | combined_data$Event==203)

data_meso_subset_means <- data_meso_subset %>%
  dplyr::group_by(Cruise,Taxa,location_group) %>%
  summarise(mean_wet_mass = mean(mean_wet_mass),log10_mass=mean(log10_mass),density_m2=mean(density_m2),SST=mean(SST),
            Lat=mean(Lat),Lon=mean(Lon))

data_meso_subset_means$origin <- "mesozooplankton"

# I now need to remove events 188 and 203 from the overall data before adding in my new averaged data
combined_data_reduced <- subset(combined_data,combined_data$Event!=188 & combined_data$Event!=203)
combined_data_reduced <- subset(combined_data_reduced,select = c(Cruise,Taxa,location_group,mean_wet_mass,log10_mass,density_m2,SST,Lat,Lon,origin))

combined_data_2 <- rbind(combined_data_reduced,data_meso_subset_means)
combined_data_2 <- combined_data_2[order(combined_data_2$location_group),]

write.csv(combined_data_2,"Data/Processed/PLANKTON_DATA.csv")




