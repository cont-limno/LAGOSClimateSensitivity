######################## Wrangling PRISM ppt data ###########################################
# Date: 9-8-17
# updated:6-21-18
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
varname = 'ppt' #tmin, tmax, tmean, ppt

first_year = 1970
last_year = 2011

##### input data #####
# Climate data (aggregated to HUC12 watersheds) can be downloaded from: 
# Collins, S. M., et al. 2018. LAGOS-NE Annual, seasonal, and monthly climate data 
# for lakes and watersheds in a 17-state region of the U.S.. Environmental Data Initiative. 
# http://dx.doi:10.6073/pasta/4abe86a2c00dc9a628924aa149d7bf34. 
# Dataset accessed 6/19/2018.

climate = read.csv(paste0("C:/Users/FWL/Dropbox/ClimateDataforIan/Annual_monthly_calculated/hu12_",varname,"_seasonal_updated.csv"))

dt = lagos_load(version = '1.087.1') #returns list of data.frame objects
names(dt)

lakes_4ha = shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/LAGOS_NE_All_Lakes_4ha/LAGOS_NE_All_Lakes_4ha.shp")
HUC12_lagoslakeid = cbind.data.frame(lakes_4ha@data$HU12_ZoneI, lakes_4ha@data$lagoslakei)
names(HUC12_lagoslakeid) = c('ZoneID','lagoslakeid')

######## d-fine functions #######
#shift_up <- function(x, n){
#  c(tail(x, -n), rep(NA, n))
#}

########################### main program ###############################
# To reproduce data for analysis, need to adjust directories below
# with help from: https://stackoverflow.com/questions/18527051/split-a-large-dataframe-into-a-list-of-data-frames-based-on-common-value-in-colu
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
  xx$fall_ppt = (xx[,11]+xx[,12]+xx[,13]) # previous year sep, oct, nov
  xx$winter_ppt = (xx[,14]+xx[,3]+xx[,4]) # previous year dec, current year jan and feb
  xx$spring_ppt = (xx[,5]+xx[,6]+xx[,6]) # current year mar, apr, may
  xx$summer_ppt = (xx[,7]+xx[,8]+xx[,9]) # current year jun, jul, aug
  xx$wyppt = (xx$fall_ppt + xx$winter_ppt + xx$spring_ppt + xx$summer_ppt)
  xx = xx[,c(1,2,12,13,14,3,4,5,6,7,8,9,10,11,15,16,17,18,19)] #rearranging months in water year order
  colnames(xx) = c('ZoneID','WY','octppt','novppt','decppt','janppt','febppt','marppt','aprppt',
                   'mayppt','junppt','julppt','augppt','sepppt','fall_ppt','winter_ppt',
                   'spring_ppt','summer_ppt','wyppt')
  # divide by 100 to correct units (ppt in mm, temp in deg C)
  xx = cbind.data.frame(xx[,1:2], (xx[,3:19]*0.01))
  write.csv(xx, file = outpathname)
}
