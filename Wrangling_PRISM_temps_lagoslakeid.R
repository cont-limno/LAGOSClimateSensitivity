######################## Wrangling PRISM climate data ###########################################
# Date: 9-11-17
# updated: 6-21-18
# Author: Ian McCullough, immccull@gmail.com
#############################################################################################

# Creates and saves individual climate files for each lagoslakeid; (~51000)
# takes a few min
# may have later issue of replication, but this is just writing data, not modeling with it

#### R libraries ####
library(LAGOSNE)
library(raster)
library(data.table)

#### define constants ####
# input PRISM data
varname = 'tmin' #tmin, tmax, tmean, ppt

first_year = 1970
last_year = 2011

##### input data #####
# Climate data (aggregated to HUC12 watersheds) can be downloaded from: 
# Collins, S. M., et al. 2018. LAGOS-NE Annual, seasonal, and monthly climate data 
# for lakes and watersheds in a 17-state region of the U.S.. Environmental Data Initiative. 
# http://dx.doi:10.6073/pasta/4abe86a2c00dc9a628924aa149d7bf34. 
# Dataset accessed 6/19/2018.
climate = read.csv(paste0("C:/Users/FWL/Dropbox/ClimateDataforIan/Annual_monthly_calculated/hu12_",varname,"_annual.csv"))

dt = lagos_load(version = '1.087.1') #returns list of data.frame objects
names(dt)

# bring in lake shapefile to get lagoslakeid (climate data don't have that column)
lakes_4ha = shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/LAGOS_NE_All_Lakes_4ha/LAGOS_NE_All_Lakes_4ha.shp")
HUC12_lagoslakeid = cbind.data.frame(lakes_4ha@data$HU12_ZoneI, lakes_4ha@data$lagoslakei)
names(HUC12_lagoslakeid) = c('ZoneID','lagoslakeid')

######## d-fine functions #######
#shift_up <- function(x, n){
#  c(tail(x, -n), rep(NA, n))
#}

########################### main program ###############################
# To reproduce data for analysis, need to adjust directories below

#### TMEAN ####
# with help from: https://stackoverflow.com/questions/18527051/split-a-large-dataframe-into-a-list-of-data-frames-based-on-common-value-in-colu
climate[1] = NULL #get rid of useless first column (auto-index)
climate_merged = merge(climate, HUC12_lagoslakeid, by='ZoneID')
climate_x = split(climate_merged, f=climate_merged$lagoslakeid) #creates separate data file for each lagoslakeid
huc_name_vector = paste0('lagoslakeid_', names(climate_x))
#metadata for PRISM: http://www.prism.oregonstate.edu/documents/PRISM_datasets.pdf

# failed attempts...may be useful later?
#climate_test = inner_join(climate, HUC12_lagoslakeid, by='ZoneID')
#climate_test = climate_test %>%
#  distinct(ZoneID, lagoslakeid, .keep_all = TRUE)

# loop writes csv for each HUC 12 for specified climate variable
for (i in 1:length(climate_x)){
  outpathname = paste0("C:/Ian_GIS/lagoslakeid_climate/",varname,"/",huc_name_vector[i],"_",varname,".csv")
  x = as.data.frame(climate_x[i])[,1:14] #eliminating seasonal columns; will recalculate
  #row_1 = rep(NA, length(xx))
  #xx = rbind.data.frame(row_1, x)
  xx=x
  xx <- xx[order(xx[,2]),] 
  xx[,12] = data.table::shift(xx[,12],n=1, fill=NA) #loses 2012 wy oct-dec; NBD for now
  xx[,13] = data.table::shift(xx[,13],n=1, fill=NA)
  xx[,14] = data.table::shift(xx[,14],n=1, fill=NA)
  
  # calculate seasonal indices
  xx$fall_tmean = (xx[,11]+xx[,12]+xx[,13])/3 # previous year sep, oct, nov
  xx$winter_tmean = (xx[,14]+xx[,3]+xx[,4])/3 # previous year dec, current year jan and feb
  xx$spring_tmean = (xx[,5]+xx[,6]+xx[,6])/3 # current year mar, apr, may
  xx$summer_tmean = (xx[,7]+xx[,8]+xx[,9])/3 # current year jun, jul, aug
  xx$wytmean = (xx$fall_tmean + xx$winter_tmean + xx$spring_tmean + xx$summer_tmean)/4
  xx = xx[,c(1,2,12,13,14,3,4,5,6,7,8,9,10,11,15,16,17,18,19)] #rearranging months in water year order
  colnames(xx) = c('ZoneID','WY','octtmean','novtmean','dectmean','jantmean','febtmean','martmean','aprtmean',
                   'maytmean','juntmean','jultmean','augtmean','septmean','fall_tmean','winter_tmean',
                   'spring_tmean','summer_tmean','wytmean')
  # divide by 100 to correct units (ppt in mm, temp in deg C)
  xx = cbind.data.frame(xx[,1:2], (xx[,3:19]*0.01))
  write.csv(xx, file = outpathname)
}

######### TMIN #########
#climate = read.csv(paste0("C:/Users/FWL/Dropbox/ClimateDataforIan/Annual_monthly_calculated/hu12_",varname,"_annual.csv"))
climate[1] = NULL #get rid of useless first column (auto-index)
climate_merged = merge(climate, HUC12_lagoslakeid, by='ZoneID')
climate_x = split(climate_merged, f=climate_merged$lagoslakeid) #creates separate data file for each lagoslakeid
huc_name_vector = paste0('lagoslakeid_', names(climate_x))
#metadata for PRISM: http://www.prism.oregonstate.edu/documents/PRISM_datasets.pdf

# loop writes csv for each HUC 12 for specified climate variable
for (i in 1:length(climate_x)){
  outpathname = paste0("C:/Ian_GIS/lagoslakeid_climate/",varname,"/",huc_name_vector[i],"_",varname,".csv")
  x = as.data.frame(climate_x[i])[,1:14] #eliminating seasonal columns; will recalculate
  #row_1 = rep(NA, length(xx))
  #xx = rbind.data.frame(row_1, x)
  xx=x
  xx <- xx[order(xx[,2]),] 
  xx[,12] = data.table::shift(xx[,12],n=1, fill=NA) #loses 2012 wy oct-dec; NBD for now
  xx[,13] = data.table::shift(xx[,13],n=1, fill=NA)
  xx[,14] = data.table::shift(xx[,14],n=1, fill=NA)
  
  # calculate seasonal indices
  xx$fall_tmin = (xx[,11]+xx[,12]+xx[,13])/3 # previous year sep, oct, nov
  xx$winter_tmin = (xx[,14]+xx[,3]+xx[,4])/3 # previous year dec, current year jan and feb
  xx$spring_tmin = (xx[,5]+xx[,6]+xx[,6])/3 # current year mar, apr, may
  xx$summer_tmin = (xx[,7]+xx[,8]+xx[,9])/3 # current year jun, jul, aug
  xx$wytmin = (xx$fall_tmin + xx$winter_tmin + xx$spring_tmin + xx$summer_tmin)/4
  xx = xx[,c(1,2,12,13,14,3,4,5,6,7,8,9,10,11,15,16,17,18,19)] #rearranging months in water year order
  colnames(xx) = c('ZoneID','WY','octtmin','novtmin','dectmin','jantmin','febtmin','martmin','aprtmin',
                   'maytmin','juntmin','jultmin','augtmin','septmin','fall_tmin','winter_tmin',
                   'spring_tmin','summer_tmin','wytmin')
  # divide by 100 to correct units (ppt in mm, temp in deg C)
  xx = cbind.data.frame(xx[,1:2], (xx[,3:19]*0.01))
  write.csv(xx, file = outpathname)
}

######### tmax #########
#climate = read.csv(paste0("C:/Users/FWL/Dropbox/ClimateDataforIan/Annual_monthly_calculated/hu12_",varname,"_annual.csv"))
climate[1] = NULL #get rid of useless first column (auto-index)
climate_merged = merge(climate, HUC12_lagoslakeid, by='ZoneID')
climate_x = split(climate_merged, f=climate_merged$lagoslakeid) #creates separate data file for each lagoslakeid
huc_name_vector = paste0('lagoslakeid_', names(climate_x))
#metadata for PRISM: http://www.prism.oregonstate.edu/documents/PRISM_datasets.pdf

# loop writes csv for each HUC 12 for specified climate variable
for (i in 1:length(climate_x)){
  outpathname = paste0("C:/Ian_GIS/lagoslakeid_climate/",varname,"/",huc_name_vector[i],"_",varname,".csv")
  x = as.data.frame(climate_x[i])[,1:14] #eliminating seasonal columns; will recalculate
  #row_1 = rep(NA, length(xx))
  #xx = rbind.data.frame(row_1, x)
  xx=x
  xx <- xx[order(xx[,2]),] 
  xx[,12] = data.table::shift(xx[,12],n=1, fill=NA) #loses 2012 wy oct-dec; NBD for now
  xx[,13] = data.table::shift(xx[,13],n=1, fill=NA)
  xx[,14] = data.table::shift(xx[,14],n=1, fill=NA)
  
  # calculate seasonal indices
  xx$fall_tmax = (xx[,11]+xx[,12]+xx[,13])/3 # previous year sep, oct, nov
  xx$winter_tmax = (xx[,14]+xx[,3]+xx[,4])/3 # previous year dec, current year jan and feb
  xx$spring_tmax = (xx[,5]+xx[,6]+xx[,6])/3 # current year mar, apr, may
  xx$summer_tmax = (xx[,7]+xx[,8]+xx[,9])/3 # current year jun, jul, aug
  xx$wytmax = (xx$fall_tmax + xx$winter_tmax + xx$spring_tmax + xx$summer_tmax)/4
  xx = xx[,c(1,2,12,13,14,3,4,5,6,7,8,9,10,11,15,16,17,18,19)] #rearranging months in water year order
  colnames(xx) = c('ZoneID','WY','octtmax','novtmax','dectmax','jantmax','febtmax','martmax','aprtmax',
                   'maytmax','juntmax','jultmax','augtmax','septmax','fall_tmax','winter_tmax',
                   'spring_tmax','summer_tmax','wytmax')
  # divide by 100 to correct units (ppt in mm, temp in deg C)
  xx = cbind.data.frame(xx[,1:2], (xx[,3:19]*0.01))
  write.csv(xx, file = outpathname)
}