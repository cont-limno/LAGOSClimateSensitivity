######################## Exploring lake color in LAGOSNE #################################################
# Date: 6-4-18
# updated: 6-4-18
# Author: Ian McCullough, immccull@gmail.com
# Note: NE (Northeast) refers to forested region, UM (Upper Midwest) refers to mixed agricultural region
##########################################################################################################

#### R libraries ####
library(LAGOSNE)
library(raster)
library(rgdal)

#### define constants ####
first_year = 1981 # full record is 1970-2011
last_year = 2010
first_day = '0615' #e.g., '0615' for Jun 15
last_day = '0915'

limno_var = 'colort'

##### input data from LAGOS NE #####
# About and links to public datasets: https://lagoslakes.org/
dt <- lagosne_load(version = '1.087.1') #returns list of data.frame objects in LAGOS NE
epi_nutr <- dt$epi_nutr #type ?epi_nutr for information about this data frame (Epilimnion Water Quality)

# GIS data downloaded and stored locally from: 
# Soranno P., K. Cheruvelil. (2017c). LAGOS-NE-GIS v1.0: A module for LAGOS-NE, 
# a multi-scaled geospatial and temporal database of lake ecological context and water 
# quality for thousands of U.S. Lakes: 2013-1925. Environmental Data Initiative. 
# Package ID: edi.98.1
# http://dx.doi.org/10.6073/pasta/fb4f5687339bec467ce0ed1ea0b5f0ca. Dataset accessed 9/26/2017.

# LAGOS_NE_All_Lakes_4ha_POINTS.zip (5.8 MB)
lakes_4ha_points <- shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/LAGOS_NE_All_Lakes_4ha_POINTS/LAGOS_NE_All_Lakes_4ha_POINTS.shp")
lakes_4ha_df <- lakes_4ha_points@data

# STATE.zip (1.8 MB)
LAGOS_NE_states <- shapefile("C:/Ian_GIS/LAGOS-NE-GISv1.0/STATE/STATE.shp")

########################### main program ###############################
# first part is wrangling data
# remove lakes without limno data by merge
sample_lakes <- merge(lakes_4ha_df, epi_nutr, by.x="lagoslakei", by.y="lagoslakeid", all.x=F)

# convert date factor to date format
sample_lakes$sampledate <- as.Date(sample_lakes$sampledate, format="%m/%d/%Y")
sample_lakes$monthday <- format(sample_lakes$sampledate, format="%m%d")

# subset by sample date cutoffs specified above
sample_lakes <- sample_lakes[sample_lakes$monthday >= first_day & sample_lakes$monthday <= last_day,]

# subset by sample months and years
sample_lakes <- subset(sample_lakes, sampleyear >= first_year & sampleyear <= last_year)

# calculate annual means for each lake (by lagoslakeid and sampleyear)
sample_data <- aggregate(sample_lakes[,limno_var], by=list(sample_lakes$lagoslakei, sample_lakes$sampleyear), 
                         FUN='mean')
colnames(sample_data) <- c('lagoslakeid','sampleyear',limno_var)
sample_data <- sample_data[!(is.na(sample_data[,limno_var]) | sample_data[,limno_var]==""), ] #remove rows with NA for limno var
sample_data_byLake <- aggregate(sample_data[,limno_var], by=list(sample_data$lagoslakeid), FUN='mean')
colnames(sample_data_byLake) <- c('lagoslakeid',limno_var)

# find number of sample years for each lake, merge to sample data frame
number_years <- as.data.frame(table(sample_data$lagoslakeid))
colnames(number_years) <- c("lagoslakeid","nYears")
number_years$nYears <- as.integer(number_years$nYears) #subset won't work unless no columns are factors
number_years$lagoslakeid <- as.numeric(levels(number_years$lagoslakeid))[number_years$lagoslakeid]

sample_data_nYears <- merge(sample_data_byLake, number_years, by='lagoslakeid')

# sort this df by lagoslakeid and sampleyear
sample_data_nYears <- sample_data_nYears[order(sample_data_nYears$lagoslakeid),]  

# keep only rows from states of interest by creating state data frame and joining
# create state_df for joining ease, which will have lagsolakeid, STATE and region
state_df <- lakes_4ha_df[,c('lagoslakei','STATE')]
state_df <- subset(state_df, STATE %in% c('ME','VT','NY','MN','WI','MI'))
sample_data_nYears <- merge(sample_data_nYears, state_df, by.x='lagoslakeid',by.y='lagoslakei', all.x=F) 
#sample_data_nYears$STATE <- NULL

# compare across regions
sample_data_nYears$Region <- NA
sample_data_nYears$Region <- ifelse(sample_data_nYears$STATE %in% c('ME','VT','NY'), 'NE','UM')

par(mfrow=c(1,1))
boxplot(colort ~ Region, sample_data_nYears, main=limno_var)
boxplot(colort ~ Region, sample_data_nYears, ylim=c(0,100), main='True color', ylab='PCU', las=1)

t.test(sample_data_nYears$colort ~ sample_data_nYears$Region)

# sample sizes?
nrow(subset(sample_data_nYears, Region == 'UM'))
nrow(subset(sample_data_nYears, Region == 'NE'))

# map color across lakes in study region
lakes_color_shp <- merge(lakes_4ha_points, sample_data_nYears, by.x='lagoslakei', by.y='lagoslakeid', all.x=F)

# simply plots where lakes are
plot(LAGOS_NE_states)
plot(lakes_color_shp, add=T, pch=20)

# export shapefile to map in ArcGIS
#writeOGR(lakes_color_shp, dsn="C:/Ian_GIS/GeographicPatterns/Color", layer="lakes_color_shp", driver="ESRI Shapefile", overwrite_layer = T)



