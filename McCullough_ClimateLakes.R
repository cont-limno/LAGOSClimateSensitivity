######################## Compare PRISM to limno data #####################################################
# Date: 9-11-17
# updated: 6-4-18
# Author: Ian McCullough, immccull@gmail.com
# Note: NE (Northeast) refers to forested region, UM (Upper Midwest) refers to mixed agricultural region
##########################################################################################################

#### R libraries ####
library(LAGOSNE)
library(raster)
library(dplyr)
library(pgirmess)
library(Hmisc)
library(randomForest)
library(rgdal)

#### define constants ####
first_year = 1981 # full record is 1970-2011
last_year = 2010
first_day = '0615' #e.g., '0615' for Jun 15
last_day = '0915'
min_years = 25 #minimum number of years of data allowed in analysis
pvalue_cutoff = 0.05 #significance level

limno_var = 'secchi'

# monthly and annual climate variables of interest
vars_of_interest = c('wyppt','fall_ppt','winter_ppt','spring_ppt','summer_ppt',
                     'wytmax','fall_tmax','winter_tmax','spring_tmax','summer_tmax',
                     'wytmin','fall_tmin','winter_tmin','spring_tmin','summer_tmin',
                     'wytmean','fall_tmean','winter_tmean','spring_tmean','summer_tmean')

##### input data from LAGOS NE #####
# About and links to public datasets: https://lagoslakes.org/
dt <- lagosne_load(version = '1.087.1') #returns list of data.frame objects in LAGOS NE
secchi <- dt$secchi #type ?secchi for information about this data frame (Secchi/Water Clarity)
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
sample_lakes <- merge(lakes_4ha_df, secchi, by.x="lagoslakei", by.y="lagoslakeid", all.x=F)

# convert date factor to date format
sample_lakes$sampledate <- as.Date(sample_lakes$sampledate, format="%m/%d/%Y")
sample_lakes$monthday <- format(sample_lakes$sampledate, format="%m%d")

# subset by sample date cutoffs specified above
sample_lakes <- sample_lakes[sample_lakes$monthday >= first_day & sample_lakes$monthday <= last_day,]

# subset by sample months and years
sample_lakes <- subset(sample_lakes, sampleyear >= first_year & sampleyear <= last_year)

# create data frame of just lagoslakeid and HUC 12 zoneIDs for merging later
HUC12_lagos_lakeIDs <- data.frame(HUC12_ZoneID = sample_lakes$HU12_ZoneI, lagoslakeid=sample_lakes$lagoslakei)
# remove duplicates of HUC12 IDs...this is the variable to join to climate data
HUC12_lagos_lakeIDs <- unique(HUC12_lagos_lakeIDs) #now have one row per unique lagos/HUC12 ID combo

# calculate annual means for each lake (by lagoslakeid and sampleyear)
sample_data <- aggregate(sample_lakes[,limno_var], by=list(sample_lakes$lagoslakei, sample_lakes$sampleyear), 
                 FUN='mean')
colnames(sample_data) <- c('lagoslakeid','sampleyear',limno_var)
sample_data <- sample_data[!(is.na(sample_data[,limno_var]) | sample_data[,limno_var]==""), ] #remove rows with NA for limno var

# find number of sample years for each lake
number_years <- as.data.frame(table(sample_data$lagoslakeid))
colnames(number_years) <- c("lagoslakeid","nYears")
number_years$nYears <- as.integer(number_years$nYears) #subset won't work unless no columns are factors
number_years$lagoslakeid <- as.numeric(levels(number_years$lagoslakeid))[number_years$lagoslakeid]
# identify lakes by lagoslakeid with minimum number of years or more
number_years_cut <- subset(number_years, nYears >= min_years)

# remove lakes without enough data years
sample_data_nYears <- merge(sample_data, number_years_cut, by='lagoslakeid')

# merge back with df with lagoslakeid and HUC 12 zone ID to get HUC 12 zone id column back
sample_data_nYears <- merge(sample_data_nYears, HUC12_lagos_lakeIDs, by='lagoslakeid')

# sort this df by lagoslakeid and sampleyear
sample_data_nYears <- sample_data_nYears[order(sample_data_nYears$lagoslakeid, sample_data_nYears$sampleyear),]  

# eliminate rows that didn't associate with a HUC 12
sample_data_nYears <- sample_data_nYears[!grepl("OUT_OF_HU12", sample_data_nYears$HUC12_ZoneID),]

# keep only rows from states of interest by creating state data frame and joining
# create state_df for joining ease, which will have lagsolakeid, STATE and region
state_df <- lakes_4ha_df[,c('lagoslakei','STATE')]
state_df <- subset(state_df, STATE %in% c('ME','VT','NY','MN','WI','MI'))
sample_data_nYears <- merge(sample_data_nYears, state_df, by.x='lagoslakeid',by.y='lagoslakei', all.x=F) 
sample_data_nYears$STATE <- NULL

# split up data by lagoslakeid to get list of lake-specific time series
sample_data_nYears_split <- split(sample_data_nYears, f=sample_data_nYears$lagoslakeid)

#### loop wrangled data through correlation with climate data (read from disk)
# output is data frame that can be saved or joined to the lake shapefile attribute table

#### ppt variables ####
split_list <- names(sample_data_nYears_split)
clim_var = 'ppt'

out_list <- list()
for (i in 1:length(split_list)){
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column

  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')
  
  cormat <- as.data.frame(suppressWarnings(cor(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], use='pairwise.complete.obs')))
  cor_one_lake <- cormat[limno_var]
  cor_one_lake[nrow(cor_one_lake),] <- clim_limno_merged$nYears[1] #make nYears in cor matrix the number of years in analysis, borrowed from another data frame
  colnames(cor_one_lake) <- split_list[i]
  cor_one_lake <- t(cor_one_lake)
  out_list[[i]] <- cor_one_lake
}

# store output in data frame
big_Daddy_ppt <- do.call(rbind.data.frame, out_list)
big_Daddy_ppt$lagoslakeid <- rownames(big_Daddy_ppt)
big_Daddy_ppt <- big_Daddy_ppt[,c(length(big_Daddy_ppt),1,3:length(big_Daddy_ppt)-1)] #reorder columns
big_Daddy_ppt$WY <- NULL

# which are significant? 
out_list_p <- list() #empty list to store for loop output
for (i in 1:length(split_list)){
  # read in and ready climate data
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  # get limno data, create correlation matrix based on numeric columns in merged climate and limno data df
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')#join by water year to match annual climate with annual limno data
  cordat <- as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)]) #create matrix of only numeric columns
  cormat <- rcorr(cordat, type="pearson")
  
  cormat <- rcorr(clim_limno_merged[,limno_var], as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], type='pearson'))
  cormat_p <- as.data.frame(cormat$P) #pull out pvalues from rcorr output
  cormat_p[,1] <- NULL # delete useless first column and row
  cormat_p <- cormat_p[-1,]
  
  cormat_df <- data.frame(p=cormat_p[,limno_var]) #create data frame, pulling out only data for limno_var
  colnames(cormat_df) <- c(split_list[i]) #name column by lagoslakeid
  out_list_p[[i]] <- t(cormat_df) #store in list created above, using transposition so each is row in combined data frame below
}

nameos <- names(clim_limno_merged[sapply(clim_limno_merged, is.numeric)])

pvalue_ppt_df <- do.call(rbind.data.frame, out_list_p)
colnames(pvalue_ppt_df) <- nameos
pvalue_ppt_df$lagoslakeid <- rownames(pvalue_ppt_df)
pvalue_ppt_df <- pvalue_ppt_df[,c(length(pvalue_ppt_df),1,3:length(pvalue_ppt_df)-1)] #reorder columns
pvalue_ppt_df$WY <- NULL

#### tmean variables ####
split_list <- names(sample_data_nYears_split)
clim_var = 'tmean'

out_list <- list()
for (i in 1:length(split_list)){
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')
  
  cormat <- as.data.frame(suppressWarnings(cor(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], use='pairwise.complete.obs')))
  cor_one_lake <- cormat[limno_var]
  cor_one_lake[nrow(cor_one_lake),] <- clim_limno_merged$nYears[1] #make nYears in cor matrix the number of years in analysis, borrowed from another data frame
  colnames(cor_one_lake) <- split_list[i]
  cor_one_lake <- t(cor_one_lake)
  out_list[[i]] <- cor_one_lake
}

# store output in data frame
big_Daddy_tmean <- do.call(rbind.data.frame, out_list)
big_Daddy_tmean$nYears <- NULL
big_Daddy_tmean[,limno_var] <- NULL
big_Daddy_tmean$WY <- NULL

# which are significant? 
out_list_p <- list() #empty list to store for loop output
for (i in 1:length(split_list)){
  # read in and ready climate data
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  # get limno data, create correlation matrix based on numeric columns in merged climate and limno data df
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')#join by water year to match annual climate with annual limno data
  cordat <- as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)]) #create matrix of only numeric columns
  cormat <- rcorr(cordat, type="pearson")
  
  cormat <- rcorr(clim_limno_merged[,limno_var], as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], type='pearson'))
  cormat_p <- as.data.frame(cormat$P) #pull out pvalues from rcorr output
  cormat_p[,1] <- NULL # delete useless first column and row
  cormat_p <- cormat_p[-1,]
  
  cormat_df <- data.frame(p=cormat_p[,limno_var]) #create data frame, pulling out only data for limno_var
  colnames(cormat_df) <- c(split_list[i]) #name column by lagoslakeid
  out_list_p[[i]] <- t(cormat_df) #store in list created above, using transposition so each is row in combined data frame below
}

nameos <- names(clim_limno_merged[sapply(clim_limno_merged, is.numeric)])

pvalue_tmean_df <- do.call(rbind.data.frame, out_list_p)
colnames(pvalue_tmean_df) <- nameos
pvalue_tmean_df$nYears <- NULL
pvalue_tmean_df[,limno_var] <- NULL
pvalue_tmean_df$WY <- NULL

#### tmin variables ####
split_list <- names(sample_data_nYears_split)
clim_var = 'tmin'

out_list <- list()
for (i in 1:length(split_list)){
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')
  
  cormat <- as.data.frame(suppressWarnings(cor(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], use='pairwise.complete.obs')))
  cor_one_lake <- cormat[limno_var]
  cor_one_lake[nrow(cor_one_lake),] <- clim_limno_merged$nYears[1] #make nYears in cor matrix the number of years in analysis, borrowed from another data frame
  colnames(cor_one_lake) <- split_list[i]
  cor_one_lake <- t(cor_one_lake)
  out_list[[i]] <- cor_one_lake
}

# store output in data frame
big_Daddy_tmin <- do.call(rbind.data.frame, out_list)
big_Daddy_tmin$nYears <- NULL
big_Daddy_tmin[,limno_var] <- NULL
big_Daddy_tmin$WY <- NULL

# which are significant? 
out_list_p <- list() #empty list to store for loop output
for (i in 1:length(split_list)){
  # read in and ready climate data
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  # get limno data, create correlation matrix based on numeric columns in merged climate and limno data df
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')#join by water year to match annual climate with annual limno data
  cordat <- as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)]) #create matrix of only numeric columns
  cormat <- rcorr(cordat, type="pearson")
  
  cormat <- rcorr(clim_limno_merged[,limno_var], as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], type='pearson'))
  cormat_p <- as.data.frame(cormat$P) #pull out pvalues from rcorr output
  cormat_p[,1] <- NULL # delete useless first column and row
  cormat_p <- cormat_p[-1,]
  
  cormat_df <- data.frame(p=cormat_p[,limno_var]) #create data frame, pulling out only data for limno_var
  colnames(cormat_df) <- c(split_list[i]) #name column by lagoslakeid
  out_list_p[[i]] <- t(cormat_df) #store in list created above, using transposition so each is row in combined data frame below
}

nameos <- names(clim_limno_merged[sapply(clim_limno_merged, is.numeric)])

pvalue_tmin_df <- do.call(rbind.data.frame, out_list_p)
colnames(pvalue_tmin_df) <- nameos
pvalue_tmin_df$nYears <- NULL
pvalue_tmin_df[,limno_var] <- NULL
pvalue_tmin_df$WY <- NULL

#### tmax variables ####
split_list <- names(sample_data_nYears_split)
clim_var = 'tmax'

out_list <- list()
for (i in 1:length(split_list)){
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')
  
  cormat <- as.data.frame(suppressWarnings(cor(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], use='pairwise.complete.obs')))
  cor_one_lake <- cormat[limno_var]
  cor_one_lake[nrow(cor_one_lake),] <- clim_limno_merged$nYears[1] #make nYears in cor matrix the number of years in analysis, borrowed from another data frame
  colnames(cor_one_lake) <- split_list[i]
  cor_one_lake <- t(cor_one_lake)
  out_list[[i]] <- cor_one_lake
}

# store output in data frame
big_Daddy_tmax <- do.call(rbind.data.frame, out_list)
big_Daddy_tmax$nYears <- NULL
big_Daddy_tmax[,limno_var] <- NULL
big_Daddy_tmax$WY <- NULL

# which are significant? 
out_list_p <- list() #empty list to store for loop output
for (i in 1:length(split_list)){
  # read in and ready climate data
  clim_path <- paste0("C:/Ian_GIS/lagoslakeid_climate/",clim_var,"/lagoslakeid_",split_list[i],"_",clim_var,".csv")
  clim_data <- read.csv(clim_path)
  clim_data <- clim_data[,-1] #delete useless auto-index column
  
  # get limno data, create correlation matrix based on numeric columns in merged climate and limno data df
  limno_data <- as.data.frame(sample_data_nYears_split[i])
  colnames(limno_data) <- c('lagoslakeid','WY',limno_var,'nYears','HUC12_ZoneID')
  clim_limno_merged <- inner_join(clim_data, limno_data, by='WY')#join by water year to match annual climate with annual limno data
  cordat <- as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)]) #create matrix of only numeric columns
  cormat <- rcorr(cordat, type="pearson")
  
  cormat <- rcorr(clim_limno_merged[,limno_var], as.matrix(clim_limno_merged[sapply(clim_limno_merged, is.numeric)], type='pearson'))
  cormat_p <- as.data.frame(cormat$P) #pull out pvalues from rcorr output
  cormat_p[,1] <- NULL # delete useless first column and row
  cormat_p <- cormat_p[-1,]
  
  cormat_df <- data.frame(p=cormat_p[,limno_var]) #create data frame, pulling out only data for limno_var
  colnames(cormat_df) <- c(split_list[i]) #name column by lagoslakeid
  out_list_p[[i]] <- t(cormat_df) #store in list created above, using transposition so each is row in combined data frame below
}

nameos <- names(clim_limno_merged[sapply(clim_limno_merged, is.numeric)])

pvalue_tmax_df <- do.call(rbind.data.frame, out_list_p)
colnames(pvalue_tmax_df) <- nameos
pvalue_tmax_df$nYears <- NULL
pvalue_tmax_df[,limno_var] <- NULL
pvalue_tmax_df$WY <- NULL

##### combine ppt, tmean, tmax, tmin into single DF to join to lake shapefile for mapping
big_Mama <- cbind.data.frame(big_Daddy_ppt, big_Daddy_tmean, big_Daddy_tmax, big_Daddy_tmin)
big_Mama_p <- cbind.data.frame(pvalue_ppt_df, pvalue_tmean_df, pvalue_tmax_df, pvalue_tmin_df)

# get data frame with only significant p values remaining (based on pvalue_cutoff)
big_Mama_pp <- as.data.frame(apply(big_Mama_p[,2:ncol(big_Mama_p)], 2, function(x) ifelse(x > pvalue_cutoff, NA, x)))
big_Mama_pp$lagoslakeid <- big_Mama_p$lagoslakeid

# counts of significant climate-WQ relationships by variable based on pvalue_cutoff
signif_counts <- apply(big_Mama_pp, 2, function(x) length(which(!is.na(x))))
signif_counts <- as.data.frame(sort(signif_counts, decreasing = T))
colnames(signif_counts) <- 'LakeCount'
signif_counts$pct_signif <- signif_counts$LakeCount/nrow(big_Mama)
signif_counts <- signif_counts[-1,]

order.pct_signif <- order(signif_counts$pct_signif, decreasing = T)

signif_counts$VarRank <- NA
signif_counts$VarRank[order.pct_signif] <- 1:nrow(signif_counts)
signif_counts$Var <- rownames(signif_counts)

# get stats on r values across climate variables
rmin <- sapply(big_Mama, min, na.rm=T) #for some reason, if do apply(big_Mama, 2, min, na.rm=T), get converted to factors and numbers messed up!! FML
rmax <- sapply(big_Mama, max, na.rm=T)
rmedian <- sapply(big_Mama, median, na.rm=T)
rstats_df <- cbind.data.frame(rmin[2:length(rmin)], rmax[2:length(rmax)], rmedian[2:length(rmedian)])
colnames(rstats_df) <- c('rmin','rmax','rmedian')
rstats_df$Var <- rownames(rstats_df)

# make big table to look at climate variable correlations with each other and WQ
ClimCor_stats_table <- merge(signif_counts, rstats_df, by='Var') 
ClimCor_stats_table <- subset(ClimCor_stats_table, Var!=c("nYears","secchi"))

# convert factors to numeric
w <- which( sapply(ClimCor_stats_table, class ) == 'factor' )
ClimCor_stats_table[w] <- lapply(ClimCor_stats_table[w], function(x) as.numeric(as.character(x)) )

# round off numeric columns to 3 decimals (dplyr)
ClimCor_stats_table <- ClimCor_stats_table %>% mutate_if(is.numeric, funs(round(., 3)))

#### join output to lake shapefile for mapping ####
lake_output <- merge(lakes_4ha_points, big_Mama, by.x='lagoslakei',by.y='lagoslakeid', all.x=F)

# can export this one and query p values below cutoff and map underneath colored values
lake_output_p <- merge(lakes_4ha_points, big_Mama_p, by.x='lagoslakei',by.y='lagoslakeid', all.x=F)

#### export shp for mapping in ArcGIS ####
dsnname <- paste0("C:/Ian_GIS/GeographicPatterns/",limno_var,"_climate_cor_lakes")
layername <- paste0(limno_var, "_climate_cor_lakes_POINTS25")
#writeOGR(lake_output, dsn=dsnname, layer=layername, driver="ESRI Shapefile", overwrite_layer = T)

dsnname <- paste0("C:/Ian_GIS/GeographicPatterns/",limno_var,"_climate_cor_lakes")
layername <- paste0(limno_var, "_climate_cor_lakes_POINTS_pval25")
#writeOGR(lake_output_p, dsn=dsnname, layer=layername, driver="ESRI Shapefile", overwrite_layer = T)

############# comparing mapped climate sensitivites by region ########################
lake_output_coords <- data.frame(lagoslakeid = lake_output@data$lagoslakei, xCor = lake_output@coords[,1], yCor = lake_output@coords[,2])
# subset by state (for regions) and merge with coords to get xy coordinates in non-spatial data frame
# NE: "northeast" (ME, NY, VT), UM: "upper midwest" (MI, MN, WI)
lake_output_NE <- subset(lake_output@data, STATE %in% c('ME','NY','VT'))
lake_output_NE <- merge(lake_output_NE, lake_output_coords, by.x='lagoslakei', by.y='lagoslakeid', all.X=F)
lake_output_UM <- subset(lake_output@data, STATE %in% c('MI','WI','MN','MO'))
lake_output_UM <- merge(lake_output_UM, lake_output_coords, by.x='lagoslakei', by.y='lagoslakeid', all.X=F)

##### compare climate-WQ correlations between regions
climateWQ_cor_NE <- lake_output_NE[,vars_of_interest]
climateWQ_cor_UM <- lake_output_UM[,vars_of_interest]

##### Northeast
# 5th percentile
cor_NE_5list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_NE[,vars_of_interest[i]], probs=c(0.05), na.rm=T)
  cor_NE_5list[i] <- dd
}

cor_NE_5 <- as.data.frame(t(rbind.data.frame(cor_NE_5list)))
rownames(cor_NE_5) <- vars_of_interest
colnames(cor_NE_5) <- 'NE_5pct'

# median
cor_NE_50list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_NE[,vars_of_interest[i]], probs=c(0.5), na.rm=T)
  cor_NE_50list[i] <- dd
}

cor_NE_50 <- as.data.frame(t(rbind.data.frame(cor_NE_50list)))
rownames(cor_NE_50) <- vars_of_interest
colnames(cor_NE_50) <- 'NE_50pct'

# 95th percentile
cor_NE_95list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_NE[,vars_of_interest[i]], probs=c(0.95), na.rm=T)
  cor_NE_95list[i] <- dd
}

cor_NE_95 <- as.data.frame(t(rbind.data.frame(cor_NE_95list)))
rownames(cor_NE_95) <- vars_of_interest
colnames(cor_NE_95) <- 'NE_95pct'

climateWQ_cor_NE_summary <- cbind.data.frame(cor_NE_5, cor_NE_50, cor_NE_95)

#### Upper Midwest
# 5th percentile
cor_UM_5list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_UM[,vars_of_interest[i]], probs=c(0.05), na.rm=T)
  cor_UM_5list[i] <- dd
}

cor_UM_5 <- as.data.frame(t(rbind.data.frame(cor_UM_5list)))
rownames(cor_UM_5) <- vars_of_interest
colnames(cor_UM_5) <- 'UM_5pct'

# median
cor_UM_50list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_UM[,vars_of_interest[i]], probs=c(0.5), na.rm=T)
  cor_UM_50list[i] <- dd
}

cor_UM_50 <- as.data.frame(t(rbind.data.frame(cor_UM_50list)))
rownames(cor_UM_50) <- vars_of_interest
colnames(cor_UM_50) <- 'UM_50pct'

# 95th percentile
cor_UM_95list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(climateWQ_cor_UM[,vars_of_interest[i]], probs=c(0.95), na.rm=T)
  cor_UM_95list[i] <- dd
}

cor_UM_95 <- as.data.frame(t(rbind.data.frame(cor_UM_95list)))
rownames(cor_UM_95) <- vars_of_interest
colnames(cor_UM_95) <- 'UM_95pct'

climateWQ_cor_UM_summary <- cbind.data.frame(cor_UM_5, cor_UM_50, cor_UM_95)

##### Combined UM and NE
combined_UM_NE_WQcor <- rbind.data.frame(climateWQ_cor_UM, climateWQ_cor_NE)
# 5th percentile
cor_All_5list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(combined_UM_NE_WQcor[,vars_of_interest[i]], probs=c(0.05), na.rm=T)
  cor_All_5list[i] <- dd
}

cor_All_5 <- as.data.frame(t(rbind.data.frame(cor_All_5list)))
rownames(cor_All_5) <- vars_of_interest
colnames(cor_All_5) <- 'All_5pct'

# median
cor_All_50list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(combined_UM_NE_WQcor[,vars_of_interest[i]], probs=c(0.5), na.rm=T)
  cor_All_50list[i] <- dd
}

cor_All_50 <- as.data.frame(t(rbind.data.frame(cor_All_50list)))
rownames(cor_All_50) <- vars_of_interest
colnames(cor_All_50) <- 'All_50pct'

# 95th percentile
cor_All_95list <- list()
for (i in 1:length(vars_of_interest)){
  dd <- quantile(combined_UM_NE_WQcor[,vars_of_interest[i]], probs=c(0.95), na.rm=T)
  cor_All_95list[i] <- dd
}

cor_All_95 <- as.data.frame(t(rbind.data.frame(cor_All_95list)))
rownames(cor_All_95) <- vars_of_interest
colnames(cor_All_95) <- 'All_95pct'

combined_UM_NE_WQcor_summary <- cbind.data.frame(cor_All_5, cor_All_50, cor_All_95)

####### get number of significant lakes by region
state_df$Region[state_df$STATE %in% c("VT","NY","ME")] <- "NE" #Northeast
state_df$Region[state_df$STATE %in% c("MI","WI","MN")] <- "UM" #Upper Midwest
signif_lakes <- left_join(big_Mama_pp, state_df, by=c('lagoslakeid'='lagoslakei')) #merge table with only signif p values with state/region df by lagoslakeid
signif_lakes_NE <- subset(signif_lakes, Region=='NE')
signif_lakes_UM <- subset(signif_lakes, Region=='UM')
signif_lakes_NE <- signif_lakes_NE[,vars_of_interest] #subset to only seasonal variables
signif_lakes_UM <- signif_lakes_UM[,vars_of_interest]
signif_lakes_all <- signif_lakes[,vars_of_interest]

# count number of non-NAs...this is the number of significant lakes, given that nonsignif p values were removed already
signif_lakes_NE_count <- apply(signif_lakes_NE, 2, function(x) length(which(!is.na(x))))
signif_lakes_UM_count <- apply(signif_lakes_UM, 2, function(x) length(which(!is.na(x))))
signif_lakes_all_count <- apply(signif_lakes_all, 2, function(x) length(which(!is.na(x))))

signif_lake_count_summary <- data.frame(UM_Signif_n= signif_lakes_UM_count, NE_Signif_n=signif_lakes_NE_count,
                                       all_Signif_n = signif_lakes_all_count)
signif_lake_count_summary$UM_pctSignif <- signif_lake_count_summary$UM_Signif_n/nrow(signif_lakes_UM)
signif_lake_count_summary$NE_pctSignif <- signif_lake_count_summary$NE_Signif_n/nrow(signif_lakes_NE)
signif_lake_count_summary$all_pctSignif <- signif_lake_count_summary$all_Signif_n/nrow(signif_lakes_all)
signif_lake_count_summary$UM_n <- nrow(signif_lakes_UM)
signif_lake_count_summary$NE_n <- nrow(signif_lakes_NE)
signif_lake_count_summary$all_n <- nrow(signif_lakes_all)

# get proportion of lakes not significantly correlated with any climate variable
length(which(!rowSums(!is.na(signif_lakes_NE))))/nrow(signif_lakes_NE)
length(which(!rowSums(!is.na(signif_lakes_UM))))/nrow(signif_lakes_UM)

########### final summary of climate-WQ correlations by region and for both regions combined and save to disk, if desired
final_clim_WQ_cor_summary <- cbind.data.frame(climateWQ_cor_UM_summary, climateWQ_cor_NE_summary, combined_UM_NE_WQcor_summary, signif_lake_count_summary)
#write.csv(final_clim_WQ_cor_summary, "C:/Ian_GIS/GeographicPatterns/Tables/ClimateWQCorrelationSummary.csv")

################ paneled histogram of climate-water clarity correlations, comparing regions ############
png('C:/Ian_GIS/GeographicPatterns/Figures/paneled_hist.png',width = 2.5,height = 3.87,units = 'in',res=600) #was 3.88 by 9 in
par(mfrow=c(2,1))
par(mar=c(2,3,0.5,0.5)) #bot,left,top,right
# first plot
plot_var = 'summer_ppt'
hist(lake_output_NE[,plot_var], xlim=c(-0.9,0.9), ylim=c(0,55), col="orange", breaks=seq(-0.9,0.9,0.1),
     xaxt='n', las=1, xlab='', ylab='', main='')
hist(lake_output_UM[,plot_var], add=T, xlim=c(-0.9,0.9), ylim=c(0,55), breaks=seq(-0.9,0.9,0.1), col=rgb(0.5,0.5,0.5, 0.5),
     xlab='', xaxt='n', yaxt='n', main='')
abline(v=0, lty=2, lwd=2)
axis(side=1, at=seq(-0.9,0.9,0.2), labels=seq(-0.9,0.9,0.2), cex.axis=0.8)
legend('topright',legend=c('NE','MW'), col=c('orange','gray'), pch=c(15,15), bty='n')
# second plot
plot_var = 'summer_tmax'
hist(lake_output_NE[,plot_var], xlim=c(-0.9,0.9), ylim=c(0,55), col="orange", breaks=seq(-0.9,0.9,0.1),
     xaxt='n', las=1, xlab='', ylab='', main='')
hist(lake_output_UM[,plot_var], add=T, xlim=c(-0.9,0.9), ylim=c(0,55), breaks=seq(-0.9,0.9,0.1), col=rgb(0.5,0.5,0.5, 0.5),
     xlab='', xaxt='n', yaxt='n', main='')
abline(v=0, lty=2, lwd=2)
axis(side=1, at=seq(-0.9,0.9,0.2), labels=seq(-0.9,0.9,0.2), cex.axis=0.8)
#mtext("Correlation coefficient (r)", side=3, at=c(0,50), line=-2, cex=0.7)
# # third plot
# par(mar=c(2.5,3,0.5,0.5)) #bot,left,top,right
# plot_var = 'fall_tmin'
# hist(lake_output_NE[,plot_var], xlim=c(-0.9,0.9), ylim=c(0,55), col="orange", breaks=seq(-0.9,0.9,0.1),
#      xaxt='n', las=1, xlab='Correlation coefficient (r)', ylab='', main='')
# hist(lake_output_UM[,plot_var], add=T, xlim=c(-0.9,0.9), ylim=c(0,55), breaks=seq(-0.9,0.9,0.1), col=rgb(0.5,0.5,0.5, 0.5),
#      xlab='', xaxt='n', yaxt='n', main='')
# axis(side=1, at=seq(-0.9,0.9,0.2), labels=seq(-0.9,0.9,0.2), cex.axis=0.8)
dev.off()

################################ gradient analysis ####################################
# relationships between climate sensitivities and ecological context variables
# e.g., what is the correlation between clarity-summer precip correlation and lake depth?

# create data frames of coordinates and climate sensitivities and p values of those correlations
x_df = as.data.frame(lake_output@data)
x_df$xCor = lake_output@coords[,1]
x_df$yCor = lake_output@coords[,2]

pval_df = as.data.frame(lake_output_p@data)
pval_df$xCor = lake_output_p@coords[,1]
pval_df$yCor = lake_output_p@coords[,2]

############ depth gradients from LAGOS data (mean, max depth; no predicted values) ###########
# warning: not all lakes necessarily have associated depth data
focal_var_vector <- names(x_df)[33:102]
depth2_gradient <- merge(big_Mama, dt$lakes_limno, by='lagoslakeid')#lose lakes without observed or predicted depth2; oh well

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  depth2_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], dt$lakes_limno, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(depth2_gradient[sapply(depth2_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

depth2_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(depth2_cor_table) <- rowname_vec
colnames(depth2_cor_table) <- focal_var_vector
depth2_cor_table <- as.data.frame(t(depth2_cor_table))

fall <- depth2_cor_table[grep("fall", rownames(depth2_cor_table)), ]
winter <- depth2_cor_table[grep("winter", rownames(depth2_cor_table)), ]
spring <- depth2_cor_table[grep("spring", rownames(depth2_cor_table)), ]
summer <- depth2_cor_table[grep("summer", rownames(depth2_cor_table)), ]
wy <- depth2_cor_table[grep("wy", rownames(depth2_cor_table)), ]

LAGOS_depth_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
depth2_gradient_df <- merge(big_Mama, dt$lakes_limno, by='lagoslakeid')
depth2_gradient_df <- merge(depth2_gradient_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
depth2_gradient_df$StateFac <- as.factor(depth2_gradient_df$STATE)

state_abbr <- levels(depth2_gradient_df$StateFac)

### divide depth gradient analysis by region 
depth_gradient2_df_UM <- subset(depth2_gradient_df, Region =='UM')
depth_gradient2_df_NE <- subset(depth2_gradient_df, Region =='NE')

# loop through gradient variable to plot correlations between them and select climate variable
depth_gradient2_variables <- names(depth2_gradient_df)[c(74,76)] #only plot mean and max depths
plot_var <- 'summer_ppt'
for (i in 1:length(depth_gradient2_variables)) {
  gradient_var <- depth_gradient2_variables[i]
  par(mfrow=c(1,2))
  plot(depth_gradient2_df_UM[,plot_var]~depth_gradient2_df_UM[,depth_gradient2_variables[i]], 
       main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
       xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
       xlab=depth_gradient2_variables[i])
  abline(0,0, lty=2)
  # create linear model to calculate coef (slope)
  lm <- lm(depth_gradient2_df_UM[,plot_var]~ depth_gradient2_df_UM[,depth_gradient2_variables[i]])
  slope <- round(lm$coefficients[2], digits=3)
  abline(lm)
  cortest <- cor.test(depth_gradient2_df_UM[,plot_var], depth_gradient2_df_UM[,depth_gradient2_variables[i]], alternative = 'two.sided',
                     method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
  plot_pval <- cortest$p.value
  legend('topright', bty='n', legend=paste0('r = ',round(cor(depth_gradient2_df_UM[,plot_var], depth_gradient2_df_UM[,depth_gradient2_variables[i]], 
                                                             use='pairwise.complete.obs'), digits=3)))
  legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
  legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
  legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(depth_gradient2_df_UM[,gradient_var]))))
  mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)
  
  plot(depth_gradient2_df_NE[,plot_var]~depth_gradient2_df_NE[,depth_gradient2_variables[i]], 
       main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
       xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
       xlab=depth_gradient2_variables[i])
  abline(0,0, lty=2)
  lm <- lm(depth_gradient2_df_NE[,plot_var]~ depth_gradient2_df_NE[,depth_gradient2_variables[i]])
  slope <- round(lm$coefficients[2], digits=3)
  abline(lm)
  cortest <- cor.test(depth_gradient2_df_NE[,plot_var], depth_gradient2_df_NE[,depth_gradient2_variables[i]], alternative = 'two.sided',
                     method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
  plot_pval <- cortest$p.value
  legend('topright', bty='n', legend=paste0('r = ',round(cor(depth_gradient2_df_NE[,plot_var], depth_gradient2_df_NE[,depth_gradient2_variables[i]], 
                                                             use='pairwise.complete.obs'), digits=3)))
  legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
  legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
  legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(depth_gradient2_df_NE[,gradient_var]))))
  mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)
}

####### iws variables (watershed, lake area) ########
focal_var_vector <- names(x_df)[33:102]
iws_gradient <- merge(big_Mama, dt$iws, by='lagoslakeid')
dt$iws$iws_lake_ratio <- dt$iws$iws_ha/dt$iws$iws_lakeareaha #calculate iws area ratio/lake area ratio

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  iws_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], dt$iws, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(iws_gradient[sapply(iws_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

iws_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(iws_cor_table) <- rowname_vec
colnames(iws_cor_table) <- focal_var_vector
iws_cor_table <- as.data.frame(t(iws_cor_table))

fall <- iws_cor_table[grep("fall", rownames(iws_cor_table)), ]
winter <- iws_cor_table[grep("winter", rownames(iws_cor_table)), ]
spring <- iws_cor_table[grep("spring", rownames(iws_cor_table)), ]
summer <- iws_cor_table[grep("summer", rownames(iws_cor_table)), ]
wy <- iws_cor_table[grep("wy", rownames(iws_cor_table)), ]

iws_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
iws_gradient_df <- merge(big_Mama, dt$iws, by='lagoslakeid')
iws_gradient_df <- merge(iws_gradient_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
iws_gradient_df$StateFac <- as.factor(iws_gradient_df$STATE)

state_abbr <- levels(iws_gradient_df$StateFac)

### divide iws gradient analysis by region 
iws_gradient_df_UM <- subset(iws_gradient_df, Region =='UM')
iws_gradient_df_NE <- subset(iws_gradient_df, Region =='NE')

## single variable plot of regions side by side
# deal with lake area outliers
iws_gradient_df_NE_areaclean <- subset(iws_gradient_df_NE, iws_lakeareaha < 1.2e4)
iws_gradient_df_UM_areaclean <- subset(iws_gradient_df_UM, iws_lakeareaha < 10000)

# plot these without the outliers
plot_var = 'summer_ppt'
gradient_var = 'iws_lake_ratio' #'iws_lake_ratio' ,'iws_lakeareaha'
par(mfrow=c(1,2))
plot(iws_gradient_df_UM_areaclean[,plot_var]~iws_gradient_df_UM_areaclean[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(0,20000), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_gradient_df_UM_areaclean[,plot_var]~ iws_gradient_df_UM_areaclean[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_gradient_df_UM_areaclean[,plot_var], iws_gradient_df_UM_areaclean[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_gradient_df_UM_areaclean[,plot_var], iws_gradient_df_UM_areaclean[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_gradient_df_UM_areaclean[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(iws_gradient_df_NE_areaclean[,plot_var]~iws_gradient_df_NE_areaclean[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(0,20000), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_gradient_df_NE_areaclean[,plot_var]~ iws_gradient_df_NE_areaclean[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_gradient_df_NE_areaclean[,plot_var], iws_gradient_df_NE_areaclean[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_gradient_df_NE_areaclean[,plot_var], iws_gradient_df_NE_areaclean[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_gradient_df_NE_areaclean[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

# size distributions of study lakes and watersheds?
#hist(iws_gradient_df_NE_areaclean$iws_lakeareaha, main='NE', xlab='ha', xlim=c(0,7000))

################ Land use/land cover (LULC) gradients in watersheds (iws) #################
iws_LULC_gradient <- merge(big_Mama, dt$iws.lulc, by='lagoslakeid')

# calculate total pct forest in 1992 NLCD (decid + conifer + mixed)
dt$iws.lulc$total_forest_pct_1992 <- (dt$iws.lulc$iws_nlcd1992_pct_41 + dt$iws.lulc$iws_nlcd1992_pct_42 +dt$iws.lulc$iws_nlcd1992_pct_43)

# calculate total pct ag in 1992 NLCD (pasture/hay + row crops + small grains)
dt$iws.lulc$total_ag_pct_1992 <- (dt$iws.lulc$iws_nlcd1992_pct_81 + dt$iws.lulc$iws_nlcd1992_pct_82 +dt$iws.lulc$iws_nlcd1992_pct_83)

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  iws_LULC_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], dt$iws.lulc, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(iws_LULC_gradient[sapply(iws_LULC_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

iws_LULC_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(iws_LULC_cor_table) <- rowname_vec
colnames(iws_LULC_cor_table) <- focal_var_vector
iws_LULC_cor_table <- as.data.frame(t(iws_LULC_cor_table))
# get rid of columns with "_ha_" (plain hectares, prefer percent of IWS)
iws_LULC_cor_table <- iws_LULC_cor_table[, -grep("_ha_", colnames(iws_LULC_cor_table))]
# delete columns with all NA
iws_LULC_cor_table <- iws_LULC_cor_table[,colSums(is.na(iws_LULC_cor_table))<nrow(iws_LULC_cor_table)]

fall <- iws_LULC_cor_table[grep("fall", rownames(iws_LULC_cor_table)), ]
winter <- iws_LULC_cor_table[grep("winter", rownames(iws_LULC_cor_table)), ]
spring <- iws_LULC_cor_table[grep("spring", rownames(iws_LULC_cor_table)), ]
summer <- iws_LULC_cor_table[grep("summer", rownames(iws_LULC_cor_table)), ]
wy <- iws_LULC_cor_table[grep("wy", rownames(iws_LULC_cor_table)), ]

iws_LULC_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
iws_gradient_lulc_df <- merge(big_Mama, dt$iws.lulc, by='lagoslakeid')
iws_gradient_lulc_df <- merge(iws_gradient_lulc_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
iws_gradient_lulc_df$StateFac <- as.factor(iws_gradient_lulc_df$STATE)

state_abbr <- levels(iws_gradient_df$StateFac)

### divide gradient analysis by region 
iws_gradient_lulc_df_UM <- subset(iws_gradient_lulc_df, Region =='UM')
iws_gradient_lulc_df_NE <- subset(iws_gradient_lulc_df, Region =='NE')

# single variable plot of regions side by side
plot_var = 'summer_tmax'
gradient_var = 'total_ag_pct_1992' #'total_ag_pct_1992', 'iws_slope_mean', 'total_forest_pct_1992'
par(mfrow=c(1,2))
plot(iws_gradient_lulc_df_UM[,plot_var]~iws_gradient_lulc_df_UM[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_gradient_lulc_df_UM[,plot_var]~ iws_gradient_lulc_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_gradient_lulc_df_UM[,plot_var], iws_gradient_lulc_df_UM[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_gradient_lulc_df_UM[,plot_var], iws_gradient_lulc_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_gradient_lulc_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(iws_gradient_lulc_df_NE[,plot_var]~iws_gradient_lulc_df_NE[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_gradient_lulc_df_NE[,plot_var]~ iws_gradient_lulc_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_gradient_lulc_df_NE[,plot_var], iws_gradient_lulc_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_gradient_lulc_df_NE[,plot_var], iws_gradient_lulc_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_gradient_lulc_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

### how do forested lakes in agriculture-dominated Upper Midwest respond?  
# subset Upper Midwest to isolate lakes with high % forest/low % ag in watershed
# create subset with watersheds with > 85% forest and <8% agriculture (median values in NE region)
UM_forested_iws_subset <- subset(iws_gradient_lulc_df_UM, total_forest_pct_1992 >= 85 & total_ag_pct_1992 <= 8)
UM_remaining_iws_subset <- subset(iws_gradient_lulc_df_UM, !(lagoslakeid %in% UM_forested_iws_subset$lagoslakeid))

# put into single boxplot
UM_forested_iws_subset$GroupID <- 'Subset'
UM_remaining_iws_subset$GroupID <- 'MainMW'
#iws_gradient_lulc_df_UM$GroupID <- 'MixedAg'
iws_gradient_lulc_df_NE$GroupID <- 'Forested'
UM_forested_iws_melted <- rbind.data.frame(UM_forested_iws_subset, UM_remaining_iws_subset)
#dev.off()
UM_NE_forested_iws_melted <- rbind.data.frame(UM_forested_iws_melted, iws_gradient_lulc_df_NE)

#png('C:/Ian_GIS/GeographicPatterns/Figures/NE_UM_Boxplot.png',width = 6,height = 6, units = 'in',res=600)  
boxplot(UM_NE_forested_iws_melted$summer_ppt ~ UM_NE_forested_iws_melted$GroupID, ylim=c(-1,1), las=1,
        ylab='Correlation coefficient (r)', main='Sensitivity to summer precipitation', col=c('orange','gray','gray'),
        names=c('NE','MW subset 1','MW subset 2'))
#dev.off()

# pairwise comparisons among groups (2 methods) (Forested region=northeastern region)
pairwise.t.test(UM_NE_forested_iws_melted$summer_ppt, UM_NE_forested_iws_melted$GroupID, p.adj='holm')
TukeyHSD(aov(summer_ppt ~ GroupID, data=UM_NE_forested_iws_melted))

################### what's the effect of a forest buffer right around the lake? ###################
buffer100_lulc <- dt$buffer100m.lulc
buffer100_lulc$total_forest_pct_1992 <- (buffer100_lulc$buffer100m_nlcd1992_pct_41 + 
                                                   buffer100_lulc$buffer100m_nlcd1992_pct_42 +buffer100_lulc$buffer100m_nlcd1992_pct_43)
buffer100_lulc_gradient <- merge(big_Mama, buffer100_lulc, by='lagoslakeid')

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  buffer100_lulc_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], buffer100_lulc, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(buffer100_lulc_gradient[sapply(buffer100_lulc_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

buffer100_lulc_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(buffer100_lulc_cor_table) <- rowname_vec
colnames(buffer100_lulc_cor_table) <- focal_var_vector
buffer100_lulc_cor_table <- as.data.frame(t(buffer100_lulc_cor_table))
# select out columns with pct and mperha in column names
buffer100_lulc_cor_table <- buffer100_lulc_cor_table %>% dplyr:: select(grep("pct", names(buffer100_lulc_cor_table)), grep("mperha", names(buffer100_lulc_cor_table)))
# delete columns with all NA
buffer100_lulc_cor_table <- buffer100_lulc_cor_table[,colSums(is.na(buffer100_lulc_cor_table))<nrow(buffer100_lulc_cor_table)]

fall <- buffer100_lulc_cor_table[grep("fall", rownames(buffer100_lulc_cor_table)), ]
winter <- buffer100_lulc_cor_table[grep("winter", rownames(buffer100_lulc_cor_table)), ]
spring <- buffer100_lulc_cor_table[grep("spring", rownames(buffer100_lulc_cor_table)), ]
summer <- buffer100_lulc_cor_table[grep("summer", rownames(buffer100_lulc_cor_table)), ]
wy <- buffer100_lulc_cor_table[grep("wy", rownames(buffer100_lulc_cor_table)), ]

buffer100_lulc_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
buffer100_gradient_df <- merge(big_Mama, buffer100_lulc, by='lagoslakeid')
buffer100_gradient_df <- merge(buffer100_gradient_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
buffer100_gradient_df$StateFac <- as.factor(buffer100_gradient_df$STATE)

state_abbr <- levels(buffer100_gradient_df$StateFac)

### divide buffer100 gradient analysis by region 
buffer100_gradient_df_UM <- subset(buffer100_gradient_df, Region =='UM')
buffer100_gradient_df_NE <- subset(buffer100_gradient_df, Region =='NE')

############### freshwater connectivity gradients ##################
iws_conn_gradient <- merge(big_Mama, dt$iws.conn, by='lagoslakeid')

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  iws_conn_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], dt$iws.conn, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(iws_conn_gradient[sapply(iws_conn_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

iws_conn_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(iws_conn_cor_table) <- rowname_vec
colnames(iws_conn_cor_table) <- focal_var_vector
iws_conn_cor_table <- as.data.frame(t(iws_conn_cor_table))
# select out columns with pct and mperha in column names
iws_conn_cor_table <- iws_conn_cor_table %>% dplyr:: select(grep("pct", names(iws_conn_cor_table)), grep("mperha", names(iws_conn_cor_table)))
# delete columns with all NA
iws_conn_cor_table <- iws_conn_cor_table[,colSums(is.na(iws_conn_cor_table))<nrow(iws_conn_cor_table)]

fall <- iws_conn_cor_table[grep("fall", rownames(iws_conn_cor_table)), ]
winter <- iws_conn_cor_table[grep("winter", rownames(iws_conn_cor_table)), ]
spring <- iws_conn_cor_table[grep("spring", rownames(iws_conn_cor_table)), ]
summer <- iws_conn_cor_table[grep("summer", rownames(iws_conn_cor_table)), ]
wy <- iws_conn_cor_table[grep("wy", rownames(iws_conn_cor_table)), ]

iws_conn_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
iws_conn_gradient_df <- merge(big_Mama, dt$iws.conn, by='lagoslakeid')
iws_conn_gradient_df <- merge(iws_conn_gradient_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
iws_conn_gradient_df$StateFac <- as.factor(iws_conn_gradient_df$STATE)

state_abbr <- levels(iws_conn_gradient_df$StateFac)

### divide iws conn gradient analysis by region 
iws_conn_gradient_df_UM <- subset(iws_conn_gradient_df, Region =='UM')
iws_conn_gradient_df_NE <- subset(iws_conn_gradient_df, Region =='NE')

# single variable plot of regions side by side
plot_var = 'summer_ppt'
gradient_var = 'iws_wl_allwetlandsdissolved_overlapping_area_pct' #iws_streamdensity_streams_density_mperha, iws_lakes_overlapping_area_pct, iws_wl_allwetlandsdissolved_overlapping_area_pct
par(mfrow=c(1,2))
plot(iws_conn_gradient_df_UM[,plot_var]~iws_conn_gradient_df_UM[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(0,40), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_conn_gradient_df_UM[,plot_var]~ iws_conn_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_conn_gradient_df_UM[,plot_var], iws_conn_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_conn_gradient_df_UM[,plot_var], iws_conn_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_conn_gradient_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(iws_conn_gradient_df_NE[,plot_var]~iws_conn_gradient_df_NE[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(0,40), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_conn_gradient_df_NE[,plot_var]~ iws_conn_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_conn_gradient_df_NE[,plot_var], iws_conn_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_conn_gradient_df_NE[,plot_var], iws_conn_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_conn_gradient_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

############### climate gradients (30-year norms in LAGOS) ##################
# includes runoff, groundwater by HUC 12
climate_gradient_df <- merge(dt$hu12.chag, HUC12_lagos_lakeIDs, by.x='hu12_zoneid',
                             by.y='HUC12_ZoneID', all.x=T, all.y=F) # need lagoslakeid to link to limno data
climate_gradient_df <- merge(climate_gradient_df, state_df, by.x='lagoslakeid', by.y='lagoslakei',all.x=T) #add state column #climate_gradient_df = merge(big_Mama, climate_gradient_df, by='lagoslakeid')
# 
cor_list <- list()
for (i in 1:length(focal_var_vector)){
  hu12_climate_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], climate_gradient_df, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(hu12_climate_gradient[sapply(hu12_climate_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

hu12_climate_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(hu12_climate_cor_table) <- rowname_vec
colnames(hu12_climate_cor_table) <- focal_var_vector
hu12_climate_cor_table <- as.data.frame(t(hu12_climate_cor_table))
# get rid of columns with "_ha_" (plain hectares, prefer percent of hu12)
hu12_climate_cor_table <- hu12_climate_cor_table[, -grep("_ha", colnames(hu12_climate_cor_table))]
# delete columns with all NA
hu12_climate_cor_table <- hu12_climate_cor_table[,colSums(is.na(hu12_climate_cor_table))<nrow(hu12_climate_cor_table)]

fall <- hu12_climate_cor_table[grep("fall", rownames(hu12_climate_cor_table)), ]
winter <- hu12_climate_cor_table[grep("winter", rownames(hu12_climate_cor_table)), ]
spring <- hu12_climate_cor_table[grep("spring", rownames(hu12_climate_cor_table)), ]
summer <- hu12_climate_cor_table[grep("summer", rownames(hu12_climate_cor_table)), ]
wy <- hu12_climate_cor_table[grep("wy", rownames(hu12_climate_cor_table)), ]

hu12_climate_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
hu12_gradient_df <- merge(big_Mama, climate_gradient_df, by='lagoslakeid')
hu12_gradient_df$StateFac <- as.factor(hu12_gradient_df$STATE)

state_abbr <- levels(hu12_gradient_df$StateFac)

#### divide hu12 climate gradient analysis by region 
hu12_gradient_df_UM <- subset(hu12_gradient_df, Region =='UM')
hu12_gradient_df_NE <- subset(hu12_gradient_df, Region=='NE')

############## climate gradient analysis with Ian-extracted PRISM ####################
# warning: slow to load large csv (contains PRISM climate normals for all lakes)
# big_climate_df from: 
# Collins, S. M., et al. 2018. LAGOS-NE Annual, seasonal, and monthly climate data for lakes 
# and watersheds in a 17-state region of the U.S.. 
# Environmental Data Initiative. http://dx.doi:10.6073/pasta/4abe86a2c00dc9a628924aa149d7bf34. 
# Dataset accessed 6/19/2018.
big_climate_df <- read.csv("C:/Ian_GIS/lagoslakeid_PRISM_Normals_1981_2010.csv")
big_climate_df[,1] <- NULL #delete useless first column
clim_var_colnames <- big_climate_df$Var #create vector of clim variable names

big_climate_df <- as.data.frame(t(big_climate_df)) #transpose to get lagoslakeid as a column
colnames(big_climate_df) <- clim_var_colnames #rename columns by clim var name vector above
big_climate_df <- big_climate_df[-1,] #delete useless first row

#create column of lagoslakeid
lagoslakeids <- rownames(big_climate_df)
lagoslakeids <- gsub("X","", lagoslakeids)
big_climate_df$lagoslakeid <- lagoslakeids
big_climate_df <- big_climate_df[,c(length(big_climate_df),1,3:length(big_climate_df)-1)]

# convert factors to numeric
w <- which( sapply(big_climate_df, class ) == 'factor' )
big_climate_df[w] <- lapply(big_climate_df[w], function(x) as.numeric(as.character(x)) )

focal_var_vector <- names(x_df)[33:102]

cor_list <- list()
for (i in 1:length(focal_var_vector)){
  PRISM_clim_gradient <- merge(big_Mama[,c(focal_var_vector[i],'lagoslakeid')], big_climate_df, by='lagoslakeid')
  whee <- as.data.frame(suppressWarnings(cor(PRISM_clim_gradient[sapply(PRISM_clim_gradient, is.numeric)], use='pairwise.complete.obs')))
  whee <- as.data.frame(whee[1])
  rowname_vec <- rownames(whee)[2:nrow(whee)]
  whee <- as.data.frame(whee[-1,]) #delete first row; correlation with itself
  names(whee) <- focal_var_vector[i]
  rownames(whee) <- rowname_vec
  cor_list[i] <- whee
}

PRISM_clim_cor_table <- do.call(cbind.data.frame, cor_list)
rownames(PRISM_clim_cor_table) <- rowname_vec
colnames(PRISM_clim_cor_table) <- focal_var_vector
PRISM_clim_cor_table <- as.data.frame(t(PRISM_clim_cor_table))
# delete columns with all NA
PRISM_clim_cor_table <- PRISM_clim_cor_table[,colSums(is.na(PRISM_clim_cor_table))<nrow(PRISM_clim_cor_table)]

fall <- PRISM_clim_cor_table[grep("fall", rownames(PRISM_clim_cor_table)), ]
winter <- PRISM_clim_cor_table[grep("winter", rownames(PRISM_clim_cor_table)), ]
spring <- PRISM_clim_cor_table[grep("spring", rownames(PRISM_clim_cor_table)), ]
summer <- PRISM_clim_cor_table[grep("summer", rownames(PRISM_clim_cor_table)), ]
wy <- PRISM_clim_cor_table[grep("wy", rownames(PRISM_clim_cor_table)), ]

PRISM_climate_cor_seasonal <- rbind.data.frame(fall,winter,spring,summer,wy)

# data frames for scatter plots
PRISM_gradient_df <- merge(big_Mama, big_climate_df, by='lagoslakeid')
PRISM_gradient_df <- merge(PRISM_gradient_df, state_df, by.x='lagoslakeid',by.y='lagoslakei')
PRISM_gradient_df$StateFac <- as.factor(PRISM_gradient_df$STATE)

state_abbr <- levels(PRISM_gradient_df$StateFac)

# adjust names manually to avoid confusion between climate variables and climate sensitivities
first_half_names <- names(PRISM_gradient_df)[1:71]
first_half_names <- gsub("ppt.x", "ppt", first_half_names)
first_half_names <- gsub("tmean.x", "tmean", first_half_names)
first_half_names <- gsub("tmin.x", "tmin", first_half_names)
first_half_names <- gsub("tmax.x", "tmax", first_half_names)

second_half_names <- names(PRISM_gradient_df[72:ncol(PRISM_gradient_df)])
second_half_names <- gsub("tmn", "tmin", second_half_names)
second_half_names <- gsub("tmx", "tmax", second_half_names)
second_half_names <- gsub("ppt.y", "ppt_PRISM", second_half_names)
second_half_names <- gsub("tmax.y", "tmax_PRISM", second_half_names)
second_half_names <- gsub("tmin.y", "tmin_PRISM", second_half_names)
second_half_names <- gsub("tmean.y", "tmean_PRISM", second_half_names)
second_half_names[69:71] <- c('wytmin_PRISM','wytmax_PRISM','wytmean_PRISM')

full_names <- c(first_half_names, second_half_names)
colnames(PRISM_gradient_df) <- full_names

### divide PRISM climate gradient analysis by region 
PRISM_gradient_df_UM <- subset(PRISM_gradient_df, Region =='UM')
PRISM_gradient_df_NE <- subset(PRISM_gradient_df, Region =='NE')

# single variable plot of regions side by side
plot_var = 'summer_tmax'
gradient_var = 'summer_tmax_PRISM'
par(mfrow=c(1,2))
plot(PRISM_gradient_df_UM[,plot_var]~PRISM_gradient_df_UM[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(20,28), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_UM[,plot_var]~ PRISM_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(PRISM_gradient_df_NE[,plot_var]~PRISM_gradient_df_NE[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(20,28), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_NE[,plot_var]~ PRISM_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

# single variable plot of regions side by side
plot_var = 'summer_ppt'
gradient_var = 'summer_ppt_PRISM'
par(mfrow=c(1,2))
plot(PRISM_gradient_df_UM[,plot_var]~PRISM_gradient_df_UM[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(200,450), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_UM[,plot_var]~ PRISM_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_UM[,plot_var], PRISM_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(PRISM_gradient_df_NE[,plot_var]~PRISM_gradient_df_NE[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(200,450), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(PRISM_gradient_df_NE[,plot_var]~ PRISM_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(PRISM_gradient_df_NE[,plot_var], PRISM_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(PRISM_gradient_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

############## calculate gradients of TP, chla and true color ##########
WQ_gradients_df <- data.frame(lagoslakeid=epi_nutr$lagoslakeid, colort = epi_nutr$colort, TP = epi_nutr$tp,
                             chla = epi_nutr$chla, DOC=epi_nutr$doc, sampledate= epi_nutr$sampledate, sampleyear=epi_nutr$sampleyear)

WQ_gradients_df$sampledate <- as.Date(WQ_gradients_df$sampledate, format="%m/%d/%Y")
WQ_gradients_df$monthday <- format(WQ_gradients_df$sampledate, format="%m%d")

# subset by sample date cutoffs specified above
WQ_gradients_df <- WQ_gradients_df[WQ_gradients_df$monthday >= first_day & WQ_gradients_df$monthday <= last_day,]
WQ_gradients_df <- subset(WQ_gradients_df, sampleyear >= first_year & sampleyear <= last_year)

## DOC 
# calculate annual means for each lake (by lagoslakeid and sampleyear)
DOC_mean <- aggregate(WQ_gradients_df$DOC, by=list(WQ_gradients_df$lagoslakeid, WQ_gradients_df$sampleyear), 
                        FUN='mean')
colnames(DOC_mean) <- c('lagoslakeid','sampleyear','DOC')

DOC_mean <- DOC_mean[!(is.na(DOC_mean$DOC) | DOC_mean$DOC==""), ] #remove rows with NA for limno var
DOC_mean_n <- DOC_mean %>% count(lagoslakeid)

# get single mean value by lagoslakeid
DOC_mean <- aggregate(DOC_mean$DOC, by=list(DOC_mean$lagoslakeid), FUN='mean')
colnames(DOC_mean) <- c('lagoslakeid','DOC')

# similar to color, few lakes in midwest with DOC data
DOC_gradient <- merge(big_Mama, DOC_mean, by='lagoslakeid')
DOC_gradient <- merge(DOC_gradient, state_df, by.x='lagoslakeid', by.y='lagoslakei')
DOC_gradient_df_UM <- subset(DOC_gradient, Region =='UM')
DOC_gradient_df_NE <- subset(DOC_gradient, Region =='NE')

boxplot(DOC_gradient$DOC ~ DOC_gradient$Region, col=c('orange','gray'), main='DOC')

map_DOC_lakes <- DOC_gradient$lagoslakeid
map_DOC_lakes <- subset(lakes_4ha_points, lagoslakei %in% map_DOC_lakes)
plot(LAGOS_NE_states)
plot(map_DOC_lakes, add=T, pch=20, col='dodgerblue')

# does DOC concentration affect summer ppt sensitivity?
plot(summer_ppt ~ DOC, DOC_gradient_df_NE, col='dodgerblue', pch=20)
cor(DOC_gradient_df_NE$DOC, DOC_gradient_df_NE$summer_ppt)

## color
# calculate annual means for each lake (by lagoslakeid and sampleyear)
color_mean <- aggregate(WQ_gradients_df$colort, by=list(WQ_gradients_df$lagoslakeid, WQ_gradients_df$sampleyear), 
                        FUN='mean')
colnames(color_mean) <- c('lagoslakeid','sampleyear','colort')

color_mean <- color_mean[!(is.na(color_mean$colort) | color_mean$colort==""), ] #remove rows with NA for limno var

# get single mean value by lagoslakeid
color_mean <- aggregate(color_mean$color, by=list(color_mean$lagoslakeid), FUN='mean')
colnames(color_mean) <- c('lagoslakeid','colort')

## chla
# calculate annual means for each lake (by lagoslakeid and sampleyear)
chla_mean <- aggregate(WQ_gradients_df$chla, by=list(WQ_gradients_df$lagoslakeid, WQ_gradients_df$sampleyear), 
                       FUN='mean')
colnames(chla_mean) <- c('lagoslakeid','sampleyear','chla')

chla_mean <- chla_mean[!(is.na(chla_mean$chla) | chla_mean$chla==""), ] #remove rows with NA for limno var

# get single mean value by lagoslakeid
chla_mean <- aggregate(chla_mean$chla, by=list(chla_mean$lagoslakeid), FUN='mean')
colnames(chla_mean) <- c('lagoslakeid','chla')

## TP
# calculate annual means for each lake (by lagoslakeid and sampleyear)
TP_mean <- aggregate(WQ_gradients_df$TP, by=list(WQ_gradients_df$lagoslakeid, WQ_gradients_df$sampleyear), 
                     FUN='mean')
colnames(TP_mean) <- c('lagoslakeid','sampleyear','TP')

TP_mean <- TP_mean[!(is.na(TP_mean$TP) | TP_mean$TP==""), ] #remove rows with NA for limno var

# get single mean value by lagoslakeid
TP_mean <- aggregate(TP_mean$TP, by=list(TP_mean$lagoslakeid), FUN='mean')
colnames(TP_mean) <- c('lagoslakeid','TP')

### next need to join these to the regional climate sensitity data frames, compare as gradients as for other vars
# COLORT (Log transforming doesn't help)
color_gradient <- merge(big_Mama, color_mean, by='lagoslakeid')
color_gradient <- merge(color_gradient, state_df, by.x='lagoslakeid', by.y='lagoslakei')
color_gradient_df_UM <- subset(color_gradient, Region =='UM')
color_gradient_df_NE <- subset(color_gradient, Region =='NE')

# side by side plots
plot_var = 'summer_ppt'
gradient_var = 'colort' #iws_streamdensity_streams_density_mperha
par(mfrow=c(1,2))
plot(color_gradient_df_UM[,plot_var]~color_gradient_df_UM[,gradient_var], 
     main=paste0(limno_var, ' Upper Midwest'), pch=20, col='red',
     xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(color_gradient_df_UM[,plot_var]~ color_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(color_gradient_df_UM[,plot_var], color_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(color_gradient_df_UM[,plot_var], color_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(color_gradient_df_UM[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

plot(color_gradient_df_NE[,plot_var]~color_gradient_df_NE[,gradient_var], 
     main=paste0(limno_var, ' Northeast'), pch=20, col='dodgerblue',
     xlim=c(), ylim=c(-1,1), ylab=paste0('sensitivity of ', limno_var, ' to ', plot_var),
     xlab=gradient_var)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(color_gradient_df_NE[,plot_var]~ color_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(color_gradient_df_NE[,plot_var], color_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval = cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(color_gradient_df_NE[,plot_var], color_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('topleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(color_gradient_df_NE[,gradient_var]))))
mtext(side=3, paste0('gradient variable = ',gradient_var, ', climate variable = ', plot_var), cex=0.75)

### TP (log transforming makes very minor difference)
TP_gradient <- merge(big_Mama, TP_mean, by='lagoslakeid')
TP_gradient <- merge(TP_gradient, state_df, by.x='lagoslakeid', by.y='lagoslakei')
TP_gradient_df_UM <- subset(TP_gradient, Region =='UM')
TP_gradient_df_NE <- subset(TP_gradient, Region =='NE')

### chla
chla_gradient <- merge(big_Mama, chla_mean, by='lagoslakeid')
chla_gradient <- merge(chla_gradient, state_df, by.x='lagoslakeid', by.y='lagoslakei')
chla_gradient_df_UM <- subset(chla_gradient, Region =='UM')
chla_gradient_df_NE <- subset(chla_gradient, Region =='NE')

#### determine number of lakes as eutrophic, mesotrophic, oligotrophic based on Dodds et al. (2006)
# based on TP: E: >30 ug/L, M: 10-30 ug/L, O: <10 ug/L
nrow(subset(TP_gradient_df_NE, TP < 10)) #oligo
nrow(subset(TP_gradient_df_NE, TP < 30 & TP >= 10)) #meso
nrow(subset(TP_gradient_df_NE, TP >= 30)) #eutro

nrow(subset(TP_gradient_df_UM, TP < 10)) #oligo
nrow(subset(TP_gradient_df_UM, TP < 30 & TP >= 10)) #meso
nrow(subset(TP_gradient_df_UM, TP >= 30)) #eutro

########## paneled gradient plots for paper ############
# summer ppt sensitivity
#png('C:/Ian_GIS/GeographicPatterns/Figures/CorPlot_summer_ppt_panel.png',width = 6,height = 6, units = 'in',res=600)  
# first plot
plot_var = 'summer_ppt'
gradient_var = 'chla'
par(mfrow=c(2,2))
par(mar=c(5,4,1,1))
plot(chla_gradient_df_NE[,plot_var]~chla_gradient_df_NE[,gradient_var], 
     pch=20, col='orange', main='',
     xlim=c(0,40), ylim=c(-1,1), ylab='correlation coefficient (r)',
     xlab='Chlorophyll-a (µg/L)', las=1)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(chla_gradient_df_NE[,plot_var]~ chla_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(chla_gradient_df_NE[,plot_var], chla_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(chla_gradient_df_NE[,plot_var], chla_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(chla_gradient_df_NE[,gradient_var]))))

# second plot
par(mar=c(5,1,1,4))
gradient_var = 'TP'
plot_var = 'summer_ppt'
plot(TP_gradient_df_NE[,plot_var]~TP_gradient_df_NE[,gradient_var], 
     main='', pch=20, col='orange',
     xlim=c(0,50), ylim=c(-1,1), ylab='', yaxt='n',
     xlab='Total phosphorus (µg/L)', las=1)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(TP_gradient_df_NE[,plot_var]~ TP_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(TP_gradient_df_NE[,plot_var], TP_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(TP_gradient_df_NE[,plot_var], TP_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(TP_gradient_df_NE[,gradient_var]))))

# third plot
plot_var = 'summer_ppt'
gradient_var = 'hu12_runoff_mean'
par(mar=c(5,4,1,1))
plot(hu12_gradient_df_NE[,plot_var]~hu12_gradient_df_NE[,gradient_var], 
     pch=20, col='orange', main='',
     xlim=c(18,32), ylim=c(-1,1), ylab='correlation coefficient (r)',
     xlab='Mean annual runoff (in/yr)')
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(hu12_gradient_df_NE[,plot_var]~ hu12_gradient_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(hu12_gradient_df_NE[,plot_var], hu12_gradient_df_NE[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(hu12_gradient_df_NE[,plot_var], hu12_gradient_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(hu12_gradient_df_NE[,gradient_var]))))

# fourth plot
plot_var = 'summer_ppt'
gradient_var = 'maxdepth'
par(mar=c(5,1,1,4))
plot(depth_gradient2_df_NE[,plot_var]~depth_gradient2_df_NE[,gradient_var], 
     main='',pch=20, col='orange',
     xlim=c(0,50), ylim=c(-1,1), ylab='', yaxt='n',
     xlab='Maximum depth (m)')
abline(0,0, lty=2)
lm <- lm(depth_gradient2_df_NE[,plot_var]~ depth_gradient2_df_NE[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm)
cortest <- cor.test(depth_gradient2_df_NE[,plot_var], depth_gradient2_df_NE[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
legend('topright', bty='n', legend=paste0('r = ',round(cor(depth_gradient2_df_NE[,plot_var], depth_gradient2_df_NE[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(depth_gradient2_df_NE[,gradient_var]))))
#dev.off()

## summer tmax sensitivity paneled plot
#png('C:/Ian_GIS/GeographicPatterns/Figures/CorPlot_summer_tmax_panel.png',width = 6,height = 6, units = 'in',res=600)  
par(mfrow=c(2,2))
# first plot
par(mar=c(5,4,1,1))
gradient_var = 'TP'
plot_var = 'summer_tmax'
plot(TP_gradient_df_UM[,plot_var]~TP_gradient_df_UM[,gradient_var], 
     main='', pch=20, col='gray66',
     xlim=c(0,50), ylim=c(-1,1), ylab='correlation coefficient (r)',
     xlab='Total phosphorus (µg/L)', las=1)
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(TP_gradient_df_UM[,plot_var]~ TP_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(TP_gradient_df_UM[,plot_var], TP_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                   method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(TP_gradient_df_UM[,plot_var], TP_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(TP_gradient_df_UM[,gradient_var]))))
#mtext(side=3, 'Average summer tmax (°C)', cex=0.75)

# second plot
plot_var = 'summer_tmax'
gradient_var = 'total_forest_pct_1992'
par(mar=c(5,1,1,4))
plot(buffer100_gradient_df_UM[,plot_var]~buffer100_gradient_df_UM[,gradient_var], 
     main='', pch=20, col='gray66', las=1,
     xlim=c(0,100), ylim=c(-1,1), ylab='', yaxt='n',
     xlab='Percent forest within 100 m buffer')
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(buffer100_gradient_df_UM[,plot_var]~ buffer100_gradient_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(buffer100_gradient_df_UM[,plot_var], buffer100_gradient_df_UM[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(buffer100_gradient_df_UM[,plot_var], buffer100_gradient_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(buffer100_gradient_df_UM[,gradient_var]))))

# third plot
plot_var = 'summer_tmax'
gradient_var = 'iws_lakeareaha'
par(mar=c(5,4,1,1))
plot(iws_gradient_df_UM_areaclean[,plot_var]~iws_gradient_df_UM_areaclean[,gradient_var], 
     pch=20, col='gray66', main='',
     xlim=c(0,2000), ylim=c(-1,1), ylab='correlation coefficient (r)',
     xlab='Lake area (ha)')
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(iws_gradient_df_UM_areaclean[,plot_var]~ iws_gradient_df_UM_areaclean[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(iws_gradient_df_UM_areaclean[,plot_var], iws_gradient_df_UM_areaclean[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(iws_gradient_df_UM_areaclean[,plot_var], iws_gradient_df_UM_areaclean[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(iws_gradient_df_UM_areaclean[,gradient_var]))))

# fourth plot
plot_var = 'summer_tmax'
gradient_var = 'maxdepth'
par(mar=c(5,1,1,4))
plot(depth_gradient2_df_UM[,plot_var]~depth_gradient2_df_UM[,gradient_var], 
     pch=20, col='gray66', main='',
     xlim=c(0,50), ylim=c(-1,1), ylab='', yaxt='n',
     xlab='Maximum depth (m)')
abline(0,0, lty=2)
# create linear model to calculate coef (slope)
lm <- lm(depth_gradient2_df_UM[,plot_var]~ depth_gradient2_df_UM[,gradient_var])
slope <- round(lm$coefficients[2], digits=3)
abline(lm, lty=1)
cortest <- cor.test(depth_gradient2_df_UM[,plot_var], depth_gradient2_df_UM[,gradient_var], alternative = 'two.sided',
                    method='pearson',conf.level=(1-pvalue_cutoff)) #get pvalue from correlation so it can be plotted
plot_pval <- cortest$p.value
legend('topright', bty='n', legend=paste0('r = ',round(cor(depth_gradient2_df_UM[,plot_var], depth_gradient2_df_UM[,gradient_var], 
                                                           use='pairwise.complete.obs'), digits=3)))
legend('bottomleft', bty='n', legend=paste0('p = ', round(plot_pval, digits=3)))
#legend('bottomleft', bty='n', legend=paste0('coef = ', slope))
legend('bottomright', bty='n', legend=paste0('n = ', length(na.omit(depth_gradient2_df_UM[,gradient_var]))))
#dev.off()

##### merge together gradient analyses from seasonal clim variables #####
seasonal_gradients_df <- cbind.data.frame(LAGOS_depth_cor_seasonal,
                                         iws_cor_seasonal, iws_LULC_cor_seasonal,
                                         iws_conn_cor_seasonal, hu12_climate_cor_seasonal,
                                         PRISM_climate_cor_seasonal)

#write.csv(seasonal_gradients_df, "C:/Ian_GIS/GeographicPatterns/SeasonalClimSensitivities_gradients.csv")

######### Summary table of basic lake characteristics by region ############
NE_TP <- quantile(TP_gradient_df_NE$TP, c(0.05,0.50,0.95), na.rm=T)
UM_TP <- quantile(TP_gradient_df_UM$TP, c(0.05,0.50,0.95), na.rm=T)

NE_chla <- quantile(chla_gradient_df_NE$chla, c(0.05,0.50,0.95), na.rm=T)
UM_chla <- quantile(chla_gradient_df_UM$chla, c(0.05,0.50,0.95), na.rm=T)

NE_color <- quantile(color_gradient_df_NE$color, c(0.05,0.50,0.95), na.rm=T)
UM_color <- quantile(color_gradient_df_UM$color, c(0.05,0.50,0.95), na.rm=T)

NE_depth <- quantile(depth_gradient2_df_NE$maxdepth, c(0.05,0.50,0.95), na.rm=T)
UM_depth <- quantile(depth_gradient2_df_UM$maxdepth, c(0.05,0.50,0.95), na.rm=T)

NE_area <- quantile(iws_gradient_df_NE_areaclean$iws_lakeareaha, c(0.05,0.50,0.95), na.rm=T)
UM_area <- quantile(iws_gradient_df_UM_areaclean$iws_lakeareaha, c(0.05,0.50,0.95), na.rm=T)

NE_area_lake_ratio <- quantile(iws_gradient_df_NE_areaclean$iws_lake_ratio, c(0.05,0.50,0.95), na.rm=T)
UM_area_lake_ratio <- quantile(iws_gradient_df_UM_areaclean$iws_lake_ratio, c(0.05,0.50,0.95), na.rm=T)

NE_slope <- quantile(iws_gradient_lulc_df_NE$iws_slope_mean, c(0.05,0.50,0.95), na.rm=T)
UM_slope <- quantile(iws_gradient_lulc_df_UM$iws_slope_mean, c(0.05,0.50,0.95), na.rm=T)

NE_forest <- quantile(iws_gradient_lulc_df_NE$total_forest_pct_1992, c(0.05,0.50,0.95), na.rm=T)
UM_forest <- quantile(iws_gradient_lulc_df_UM$total_forest_pct_1992, c(0.05,0.50,0.95), na.rm=T)

NE_ag <- quantile(iws_gradient_lulc_df_NE$total_ag_pct_1992, c(0.05,0.50,0.95), na.rm=T)
UM_ag <- quantile(iws_gradient_lulc_df_UM$total_ag_pct_1992, c(0.05,0.50,0.95), na.rm=T)

NE_stream <- quantile(iws_conn_gradient_df_NE$iws_streamdensity_streams_density_mperha, c(0.05,0.50,0.95), na.rm=T)
UM_stream <- quantile(iws_conn_gradient_df_UM$iws_streamdensity_streams_density_mperha, c(0.05,0.50,0.95), na.rm=T)

NE_wetland <- quantile(iws_conn_gradient_df_NE$iws_wl_allwetlandsdissolved_overlapping_area_pct, c(0.05,0.50,0.95), na.rm=T)
UM_wetland <- quantile(iws_conn_gradient_df_UM$iws_wl_allwetlandsdissolved_overlapping_area_pct, c(0.05,0.50,0.95), na.rm=T)

NE_runoff <- quantile(hu12_gradient_df_NE$hu12_runoff_mean, c(0.05,0.50,0.95), na.rm=T)
UM_runoff <- quantile(hu12_gradient_df_UM$hu12_runoff_mean, c(0.05,0.50,0.95), na.rm=T)

NE_GW <- quantile(hu12_gradient_df_NE$hu12_groundwaterrecharge_mean, c(0.05,0.50,0.95), na.rm=T)
UM_GW <- quantile(hu12_gradient_df_UM$hu12_groundwaterrecharge_mean, c(0.05,0.50,0.95), na.rm=T)

lake_summary_table_NE <- as.data.frame(t(data.frame(lakearea = NE_area, maxdepth=NE_depth, lake_iws_ratio = NE_area_lake_ratio,
                                   slope=NE_slope, forest=NE_forest, ag=NE_ag, wetland=NE_wetland,
                                   stream_density=NE_stream, runoff=NE_runoff, GW=NE_GW, TP=NE_TP,
                                   chla=NE_chla, color=NE_color)))
#write.csv(lake_summary_table_NE, file = "C:/Ian_GIS/GeographicPatterns/Tables/lake_summary_table_NE.csv")

lake_summary_table_UM <- as.data.frame(t(data.frame(lakearea = UM_area, maxdepth=UM_depth, lake_iws_ratio = UM_area_lake_ratio,
                                                   slope=UM_slope, forest=UM_forest, ag=UM_ag, wetland=UM_wetland,
                                                   stream_density=UM_stream, runoff=UM_runoff, GW=UM_GW, TP=UM_TP,
                                                   chla=UM_chla, color=UM_color)))
#write.csv(lake_summary_table_UM, file = "C:/Ian_GIS/GeographicPatterns/Tables/lake_summary_table_UM.csv")

######## Random forest #########
par(mfrow=c(1,1))
## UM - summer tmax sensitvity predictors
# Create data frame for UM with variables of interest
# summer_tmax sensitivity, lake area, lake max depth, forest 100 m buffer, mean summer tmax

t1_df <- data.frame(lagoslakeid=as.character(PRISM_gradient_df_UM$lagoslakeid), summer_tmax = PRISM_gradient_df_UM$summer_tmax,
                   summer_tmax_PRISM = PRISM_gradient_df_UM$summer_tmax_PRISM)
t2_df <- data.frame(buffer100_gradient_df_UM$total_forest_pct_1992, lagoslakeid=as.character(buffer100_gradient_df_UM$lagoslakeid))
t3_df <- data.frame(depth_gradient2_df_UM$maxdepth, lagoslakeid=as.character(depth_gradient2_df_UM$lagoslakeid))
t4_df <- data.frame(iws_gradient_df_UM_areaclean$iws_lakeareaha, lagoslakeid=as.character(iws_gradient_df_UM_areaclean$lagoslakeid))
#t5_df <- data.frame(color_mean$colort, lagoslakeid=as.character(color_mean$lagoslakeid)) #color omitted due to low sample size
t6_df <- data.frame(TP_mean$TP, lagoslakeid=as.character(TP_mean$lagoslakeid))
t7_df <- data.frame(chla_mean$chla, lagoslakeid=as.character(chla_mean$lagoslakeid))
RF_UM_df <- left_join(t1_df, t2_df, by='lagoslakeid')     
RF_UM_df <- left_join(RF_UM_df, t3_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t4_df, by='lagoslakeid')
#RF_UM_df <- left_join(RF_UM_df, t5_df, by='lagoslakeid') #not using color in UM due to low number of lakes with data
RF_UM_df <- left_join(RF_UM_df, t6_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t7_df, by='lagoslakeid')

colnames(RF_UM_df) <- c('lagoslakeid','summer_tmax','summer_tmax_PRISM','forest100m','maxdepth','lakeareaha','TP','chla')
RF_UM_df <- RF_UM_df[complete.cases(RF_UM_df), ]
set.seed(999)
RF_UM_tmax <- randomForest(summer_tmax ~ forest100m + summer_tmax_PRISM + TP + chla + maxdepth + lakeareaha, 
                        data=RF_UM_df, ntree=500, importance=T)
print(RF_UM_tmax)
importance(RF_UM_tmax, type=1)
importance(RF_UM_tmax, type=2) #MSE for regression

# plot1 is %Increase in MSE, with larger increases for a given variable indicating that variable is relatively more important
# high %IncMSE means without that variable, that much predictive capacity is lost 
# a neg value therefore means you are better off not including that variable
varImpPlot(RF_UM_tmax, main='Variable importance, summer tmax sensitivity - UM')

# diagnostic for how many trees to run? Processing is done in seconds, so I don't see major downside to large # of trees
plot(RF_UM_tmax, xlim=c(0,500))

## NE - summer ppt sensitivity predictors
# Create data frame for NE with variables of interest
# summer_ppt sensitivity, lake area, lake max depth, watershed mean slope, summer_ppt_PRISM, % IWS forest,
# % IWS ag, % IWS wetlands, % IWS lakes, stream density_mperha, runoff, GW, watershed/lake area
t1_df <- data.frame(lagoslakeid=PRISM_gradient_df_NE$lagoslakeid, summer_ppt = PRISM_gradient_df_NE$summer_ppt,
                   summer_ppt_PRISM = PRISM_gradient_df_NE$summer_ppt_PRISM)
t2_df <- data.frame(wetland_pct = iws_conn_gradient_df_NE$iws_wl_allwetlandsdissolved_overlapping_area_pct, 
                   lake_pct = iws_conn_gradient_df_NE$iws_lakes_overlapping_area_pct,
                   stream_density = iws_conn_gradient_df_NE$iws_streamdensity_streams_density_mperha,
                   lagoslakeid=iws_conn_gradient_df_NE$lagoslakeid)
t3_df <- data.frame(maxdepth = depth_gradient2_df_NE$maxdepth, lagoslakeid=depth_gradient2_df_NE$lagoslakeid)
t4_df <- data.frame(lakeareaha = iws_gradient_df_NE_areaclean$iws_lakeareaha, lagoslakeid=iws_gradient_df_NE_areaclean$lagoslakeid,
                   iws_lake_ratio = iws_gradient_df_NE_areaclean$iws_lake_ratio)
t5_df <- data.frame(lagoslakeid= iws_gradient_lulc_df_NE$lagoslakeid, total_forest_pct = iws_gradient_lulc_df_NE$total_forest_pct_1992,
                   total_ag_pct = iws_gradient_lulc_df_NE$total_ag_pct_1992,
                   slope_mean = iws_gradient_lulc_df_NE$iws_slope_mean)
t6_df <- data.frame(lagoslakeid = hu12_gradient_df_NE$lagoslakeid, runoff=hu12_gradient_df_NE$hu12_runoff_mean,
                   baseflow = hu12_gradient_df_NE$hu12_baseflowindex_mean,
                   groundwater = hu12_gradient_df_NE$hu12_groundwaterrecharge_mean)
#t7_df <- data.frame(colort=color_mean$colort, lagoslakeid=as.character(color_mean$lagoslakeid)) #color omitted due to low sample size
t8_df <- data.frame(TP=TP_mean$TP, lagoslakeid=as.character(TP_mean$lagoslakeid))
t9_df <- data.frame(chla=chla_mean$chla, lagoslakeid=as.character(chla_mean$lagoslakeid))
RF_NE_df <- left_join(t1_df, t2_df, by='lagoslakeid')     
RF_NE_df <- left_join(RF_NE_df, t3_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t4_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t5_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t6_df, by='lagoslakeid')
#RF_NE_df <- left_join(RF_NE_df, t7_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t8_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t9_df, by='lagoslakeid')
RF_NE_df <- RF_NE_df[complete.cases(RF_NE_df), ]
# not using lake pct in RF analysis because has so many zeros
set.seed(888)
RF_NE_ppt <- randomForest(summer_ppt ~ lakeareaha + summer_ppt_PRISM + maxdepth + wetland_pct +
                       stream_density + iws_lake_ratio + total_forest_pct + total_ag_pct + slope_mean +
                       runoff + groundwater + chla + TP , data=RF_NE_df, ntree=500, importance=T)
print(RF_NE_ppt)
importance(RF_NE_ppt, type=1)
importance(RF_NE_ppt, type=2) #MSE for regression

# plot1 is %Increase in MSE, with larger increases for a given variable indicating that variable is relatively more important
# high %IncMSE means without that variable, that much predictive capacity is lost 
# a neg value therefore means you are better off not including that variable
varImpPlot(RF_NE_ppt, main='Variable importance, summer ppt sensitivity - NE')

# diagnostic for how many trees to run? Processing is done in seconds, so I don't see major downside to large # of trees unless error increases
plot(RF_NE_ppt, xlim=c(0,500))

#### reciprocal random forests for comparison across regions
# NE - summer tmax
t1_df <- data.frame(lagoslakeid=as.character(PRISM_gradient_df_NE$lagoslakeid), summer_tmax = PRISM_gradient_df_NE$summer_tmax,
                   summer_tmax_PRISM = PRISM_gradient_df_NE$summer_tmax_PRISM)
t2_df <- data.frame(buffer100_gradient_df_NE$total_forest_pct_1992, lagoslakeid=as.character(buffer100_gradient_df_NE$lagoslakeid))
t3_df <- data.frame(depth_gradient2_df_NE$maxdepth, lagoslakeid=as.character(depth_gradient2_df_NE$lagoslakeid))
t4_df <- data.frame(iws_gradient_df_NE_areaclean$iws_lakeareaha, lagoslakeid=as.character(iws_gradient_df_NE_areaclean$lagoslakeid))
#t5_df <- data.frame(color_mean$colort, lagoslakeid=as.character(color_mean$lagoslakeid))
t6_df <- data.frame(TP_mean$TP, lagoslakeid=as.character(TP_mean$lagoslakeid))
t7_df <- data.frame(chla_mean$chla, lagoslakeid=as.character(chla_mean$lagoslakeid))
RF_NE_df <- left_join(t1_df, t2_df, by='lagoslakeid')     
RF_NE_df <- left_join(RF_NE_df, t3_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t4_df, by='lagoslakeid')
#RF_NE_df <- left_join(RF_NE_df, t5_df, by='lagoslakeid') #not using color in NE due to low number of lakes with data
RF_NE_df <- left_join(RF_NE_df, t6_df, by='lagoslakeid')
RF_NE_df <- left_join(RF_NE_df, t7_df, by='lagoslakeid')

colnames(RF_NE_df) <- c('lagoslakeid','summer_tmax','summer_tmax_PRISM','forest100m','maxdepth','lakeareaha','TP','chla')
RF_NE_df <- RF_NE_df[complete.cases(RF_NE_df), ]
set.seed(777)
RF_NE_tmax <- randomForest(summer_tmax ~ lakeareaha + summer_tmax_PRISM + maxdepth + TP + chla,
                        data=RF_NE_df, ntree=500, importance=T) #forest omitted due to neg variance explained
print(RF_NE_tmax)
importance(RF_NE_tmax, type=1)
importance(RF_NE_tmax, type=2) #MSE for regression

# plot1 is %Increase in MSE, with larger increases for a given variable indicating that variable is relatively more important
# high %IncMSE means without that variable, that much predictive capacity is lost 
# a neg value therefore means you are better off not including that variable
varImpPlot(RF_NE_tmax, main='Variable importance, summer tmax sensitivity - NE')

# diagnostic for how many trees to run? Processing is done in seconds, so I don't see major downside to large # of trees
plot(RF_NE_tmax, xlim=c(0,5000))

## UM - summer ppt
# Create data frame for UM with variables of interest
# summer_ppt sensitivity, lake area, lake max depth, watershed mean slope, summer_ppt_PRISM, % IWS forest,
# % IWS ag, % IWS wetlands, % IWS lakes, stream density_mperha, runoff, GW, watershed/lake area
t1_df <- data.frame(lagoslakeid=PRISM_gradient_df_UM$lagoslakeid, summer_ppt = PRISM_gradient_df_UM$summer_ppt,
                   summer_ppt_PRISM = PRISM_gradient_df_UM$summer_ppt_PRISM)
t2_df <- data.frame(wetland_pct = iws_conn_gradient_df_UM$iws_wl_allwetlandsdissolved_overlapping_area_pct, 
                   lake_pct = iws_conn_gradient_df_UM$iws_lakes_overlapping_area_pct,
                   stream_density = iws_conn_gradient_df_UM$iws_streamdensity_streams_density_mperha,
                   lagoslakeid=iws_conn_gradient_df_UM$lagoslakeid)
t3_df <- data.frame(maxdepth = depth_gradient2_df_UM$maxdepth, lagoslakeid=depth_gradient2_df_UM$lagoslakeid)
t4_df <- data.frame(lakeareaha = iws_gradient_df_UM_areaclean$iws_lakeareaha, lagoslakeid=iws_gradient_df_UM_areaclean$lagoslakeid,
                   iws_lake_ratio = iws_gradient_df_UM_areaclean$iws_lake_ratio)
t5_df <- data.frame(lagoslakeid= iws_gradient_lulc_df_UM$lagoslakeid, total_forest_pct = iws_gradient_lulc_df_UM$total_forest_pct_1992,
                   total_ag_pct = iws_gradient_lulc_df_UM$total_ag_pct_1992,
                   slope_mean = iws_gradient_lulc_df_UM$iws_slope_mean)
t6_df <- data.frame(lagoslakeid = hu12_gradient_df_UM$lagoslakeid, runoff=hu12_gradient_df_UM$hu12_runoff_mean,
                   baseflow = hu12_gradient_df_UM$hu12_baseflowindex_mean,
                   groundwater = hu12_gradient_df_UM$hu12_groundwaterrecharge_mean)
#t7_df <- data.frame(colort=color_mean$colort, lagoslakeid=as.character(color_mean$lagoslakeid))
t8_df <- data.frame(TP=TP_mean$TP, lagoslakeid=as.character(TP_mean$lagoslakeid))
t9_df <- data.frame(chla=chla_mean$chla, lagoslakeid=as.character(chla_mean$lagoslakeid))
RF_UM_df <- left_join(t1_df, t2_df, by='lagoslakeid')     
RF_UM_df <- left_join(RF_UM_df, t3_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t4_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t5_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t6_df, by='lagoslakeid')
#RF_UM_df <- left_join(RF_UM_df, t7_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t8_df, by='lagoslakeid')
RF_UM_df <- left_join(RF_UM_df, t9_df, by='lagoslakeid')
RF_UM_df <- RF_UM_df[complete.cases(RF_UM_df), ]
# omitted slope, stream density, ag, iws_lake_ratio due to neg influence
set.seed(666)
RF_UM_ppt <- randomForest(summer_ppt ~ lakeareaha + summer_ppt_PRISM + maxdepth + wetland_pct +
                       total_forest_pct + iws_lake_ratio + stream_density+ total_ag_pct + 
                       runoff + groundwater + chla + TP , data=RF_UM_df, ntree=500, importance=T)
print(RF_UM_ppt)
importance(RF_UM_ppt, type=1)
importance(RF_UM_ppt, type=2) #MSE for regression

# plot1 is %Increase in MSE, with larger increases for a given variable indicating that variable is relatively more important
# high %IncMSE means without that variable, that much predictive capacity is lost 
# a neg value therefore means you are better off not including that variable
varImpPlot(RF_UM_ppt, main='Variable importance, summer ppt sensitivity - UM')

# diagnostic for how many trees to run? Processing is done in seconds, so I don't see major downside to large # of trees unless error increases
plot(RF_UM_ppt, xlim=c(0,500))
######################## end ###########################