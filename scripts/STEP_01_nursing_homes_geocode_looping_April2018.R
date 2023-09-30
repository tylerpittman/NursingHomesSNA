## Takes Medicare nursing home SAS data and finds geocoordinates from zipcodes
## Tyler Pittman, April 2018
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_01_nursing_homes_geocode_looping_April2018.R

##########################################################
#### DON'T RUN THIS SCRIPT ANYMORE, GOT DATA PREPARED ####
##########################################################

setwd("/Users/tylerpittman/GitHub/NursingHomesSNA/data");

library(dplyr); 	# Gives inner_join function for data merges like SAS, don't load if using plyr package
library(sqldf);		# Does SQL joins similar to SAS
library(TeachingDemos);
library(sgeostat);
#library(car);
library(sp);
library(shapefiles);
library(PBSmapping);
library(mapdata);
library(MASS);
library(maptools);
library(spBayes); 	## important gives Bayesian spatial modeling functions
library(RColorBrewer);
library(classInt);     	# finds class intervals for continuous variables
library(MBA);		# make sure to have `boost headers for C++' installed on linux system	
#library(gpclib);
library(rgeos);
library(raster);
library(fossil);
library(akima);
library(mapproj);
library(rgdal);
library(fields);
gpclibPermit() 		##need to use this to permit open license for unionSpatialPolygons feature
#library(RSurvey);	# requires x11 to be installed in ubuntu
#library(rgdal);
library(spatstat); 	## important gives point process modeling functions for kernel density
library(gstat); 	### Must load this after spatstat or will receive error message!!
#library(RMySQL);
#library(spTimer);
library(imputation);
library(xtable);
library(spam);
library(haven); 	# Read in SAS or Stata files directly
library(foreign);	# Write SAS or Stata files directly
library(RJSONIO);	# Use Google API for lat and lon
#library(plyr);		# Gives a nice join function for merges, BUT USE dpylr instead!
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
library(ggmap); 	# Google Maps geocode function alterative
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function


ll = "+proj=longlat";
usCounties <- readShapeSpatial(paste("cb_2016_us_county_20m", ".shp", sep=""), proj4string=CRS(ll));
usCounties.dbf <- read.dbf(paste("cb_2016_us_county_20m", ".dbf", sep=""), header=TRUE);
#proj4string(usCounties) <- CRS("+proj=longlat");

#projected <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
projected <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs";
usCountiesProj <- spTransform(usCounties, CRS(projected));

geocodeAdddress <- function(address) {
  require(RJSONIO)
  url <- "http://maps.google.com/maps/api/geocode/json?address="
  url <- URLencode(paste(url, address, "&sensor=false", sep = ""))
  x <- fromJSON(url, simplify = FALSE)
  if (x$status == "OK") {
    out <- c(x$results[[1]]$geometry$location$lng,
             x$results[[1]]$geometry$location$lat)
  } else {
    out <- NA
  }
  Sys.sleep(0.2)  # API only allows 5 requests per second
  #Sys.sleep(5)  # API only allows 5 requests per second
  out
};


df <- read_sas("nhs_combined.sas7bdat"); #uses tibbles automatically; 77,673 NHs; #df should be 77,669 rows x 61 columns;
#write_sas(nhs_combined.sas7bdat);
#df$lat;
#df$lon;
nhNA <- df[which(is.na(df$lat)),]; #has 61 variables, including PERIOD;
#nhNA10 <- nhNA[c(1:10),];
#nhNA <- nhNA10;
#ex <- geocodeAdddress("Time Square, New York City");
#[1] -73.98722  40.7575

nhLLmiss <- read.csv("nhLLmiss_TMP.csv", header = TRUE);
nhLLmiss$Federal_Provider_Number <- as.character(nhLLmiss$Federal_Provider_Number);
nhLLmiss$Provider_Address <- as.character(nhLLmiss$Provider_Address);
nhLLmiss$Provider_City <- as.character(nhLLmiss$Provider_City);
nhLLmiss$Provider_State <- as.character(nhLLmiss$Provider_State);
nhLLmiss$Provider_Zip_Code <- as.character(nhLLmiss$Provider_Zip_Code);
nhLLmiss$PERIOD <- as.character(nhLLmiss$PERIOD);
#nhLLmiss[duplicated(nhLLmiss[, which(colnames(nhLLmiss) %in% c("Federal_Provider_Number", "PERIOD"))]),]; #shows duplicates;
nhLLmiss <- nhLLmiss[!duplicated(nhLLmiss[, which(colnames(nhLLmiss) %in% c("Federal_Provider_Number", "PERIOD"))]),]; #remove duplicates for Federal_Provider_Number and PERIOD;

##df[which(df$Federal_Provider_Number %in% nhLLmiss$Federal_Provider_Number), ]$lat <- nhLLmiss$lat;
##df[which(df$Federal_Provider_Number %in% nhLLmiss$Federal_Provider_Number), ]$lon <- nhLLmiss$lon;
##test <- join(df, nhLLmiss, by=c("Federal_Provider_Number"), type="full", match="first");

df$Federal_Provider_Number <- as.character(df$Federal_Provider_Number);
nhLLmiss <- nhLLmiss[, which(colnames(nhLLmiss) %in% c("Federal_Provider_Number", "PERIOD", "lat", "lon"))];
#str(df$Federal_Provider_Number);
#str(nhLLmiss$Federal_Provider_Number);
#str(df$lat);
#str(df$lon);
#str(nhLLmiss);



###----###----###----###----###----###----###----###----###----###----###----###----###----###----###
###----###----###----###----###----###----###----###----###----###----###----###----###----###----###
#################### SUBMIT THIS SECTION RECURSIVELY TO UPDATE GEOCOORDINATES #######################
#---------- Use library(sqldf) in this section for merging and rename test to df to resubmit nhNA below in cycles ------------#
#test <- full_join(df, nhLLmiss, copy=TRUE); #doesn't work, puts duplicate entires each entry and can't fix this;
#test <- sqldf("select df.*, nhLLmiss.* from (df left join nhLLmiss on df.Federal_Provider_Number = nhLLmiss.Federal_Provider_Number)");
#test <- merge(df, nhLLmiss, by="Federal_Provider_Number", all.x=TRUE);
#--------------------------------#

##df <- read.csv("nhLLcomplete_TMP.csv", header = TRUE); #must load this file each time when restarting program for more complete lat and lon;
## Might want to load below instead for first merge;
#nhLLmiss <- read.csv("nhLLcomplete_TMP.csv", header = TRUE); #must load this file each time when restarting program for more complete lat and lon;
##nhLLmiss$Federal_Provider_Number <- as.character(nhLLmiss$Federal_Provider_Number);
##nhLLmiss$PERIOD <- as.character(nhLLmiss$PERIOD);
##nhLLmiss$lat <- as.numeric(as.character(nhLLmiss$lat));
##nhLLmiss$lon <- as.numeric(as.character(nhLLmiss$lon));

#geocode("801 E GRANT         MEADE             KS");
#geocode("801 E GRANT, MEADE, KS", force = TRUE);

nhLLmiss <- nhLLmiss[, which(colnames(nhLLmiss) %in% c("Federal_Provider_Number", "PERIOD", "lat", "lon"))];
nhLLmiss <- nhLLmiss[!duplicated(nhLLmiss[, which(colnames(nhLLmiss) %in% c("Federal_Provider_Number", "PERIOD"))]),]; #remove duplicates for Federal_Provider_Number and PERIOD;
test <- left_join(df, nhLLmiss, by=c("Federal_Provider_Number", "PERIOD"), copy=TRUE);
length(df$lat); #Should be 77673;
length(test$lat.x); #Should be 77673;
#colnames(test);
test$lat = test$lat.x;
test[which(is.na(test$lat)), ]$lat <- test[which(is.na(test$lat)), ]$lat.y;
test$lon = test$lon.x;
test[which(is.na(test$lon)), ]$lon <- test[which(is.na(test$lon)), ]$lon.y;
test <- test[, which(colnames(test) %in% colnames(df))];
length(which(is.na(df$lat)));
length(which(is.na(test$lat)));
length(df);
length(test);
length(df$lat); #Should be 77673;
length(test$lat); #Should be 77673;
df[which(df$Federal_Provider_Number == "04A329"),]$lat;
nhLLmiss[which(nhLLmiss$Federal_Provider_Number == "04A329"),]$lat;
test[which(test$Federal_Provider_Number == "04A329"),]$lat;
test[which(test$Federal_Provider_Number == "04A329"),]$Federal_Provider_Number;
###Delete whole numbers because they are errors;
#test[which(test$lat%%1==0), ]$lat <- NA;
#test[which(test$lon%%1==0), ]$lon <- NA;
df <- test;
write.table(df, file = "nhLLcomplete_TMP.csv",row.names=FALSE, na="",col.names=TRUE, sep=",");
nhNA <- df[which(is.na(df$lat)),];
nhNA$Provider_Address <- as.character(nhNA$Provider_Address);
nhNA$Provider_City <- as.character(nhNA$Provider_City);
nhNA$Provider_State <- as.character(nhNA$Provider_State);
nhNA$Provider_Zip_Code <- as.character(nhNA$Provider_Zip_Code);
#-----------------------------------------------------------------------------------------------------------------------------#

outputSM <- matrix(NA, 1, length(c("Federal_Provider_Number", "PERIOD", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Location", "lat", "lon")));
outputSM <- as.data.frame(outputSM);
colnames(outputSM) <- c("Federal_Provider_Number", "PERIOD", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Location", "lat", "lon");
outputSM <- outputSM[-1,];
for(i in 1:length(unique(nhNA$Federal_Provider_Number))){
  	if (nhNA$Location[i] == "") {
    	#submitted <- paste(nhNA$Provider_Address[i], ", ", nhNA$Provider_City[i], ", ", nhNA$Provider_State[i], ", ", nhNA$Provider_Zip_Code[i])
    	#submitted <- paste(nhNA$Provider_City[i], ", ", nhNA$Provider_State[i], ", ", nhNA$Provider_Zip_Code[i])
    	submitted <- paste(nhNA$Provider_Address[i], ", ", nhNA$Provider_City[i], ", ", nhNA$Provider_State[i])
    	#submitted <- paste(nhNA$Provider_Name[i], ", ", nhNA$Provider_City[i], ", ", nhNA$Provider_State[i])
    	#submitted <- paste(nhNA$Provider_Address[i],  nhNA$Provider_City[i], nhNA$Provider_State[i])
    	#submitted <- paste(nhNA$Provider_City[i], ", ", nhNA$Provider_State[i])
  	} else {
    	submitted <- nhNA$Location[i]
  	}
	tmpJ3 <- geocodeAdddress(submitted);
	#tmpJ3 <- geocode(submitted);
	tmpSM <- cbind(nhNA$Federal_Provider_Number[i], nhNA$PERIOD[i], nhNA$Provider_Address[i], nhNA$Provider_City[i], nhNA$Provider_State[i], nhNA$Provider_Zip_Code[i], nhNA$Location[i], tmpJ3[2], tmpJ3[1]);
	outputSM <- rbind(outputSM, tmpSM);
}
colnames(outputSM) <- c("Federal_Provider_Number", "PERIOD", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Location", "lat", "lon");
nhLLmiss <- outputSM;
nhLLmiss$lat <- as.numeric(as.character(nhLLmiss$lat));
nhLLmiss$lon <- as.numeric(as.character(nhLLmiss$lon));
#### save temp ####
write.table(nhLLmiss, file = "nhLLmiss_TMP.csv",row.names=FALSE, na="",col.names=TRUE, sep=",");
#nhLLmiss <- read.csv("nhLLmiss_TMP.csv", header = TRUE);
###----###----###----###----###----###----###----###----###----###----###----###----###----###----###
###----###----###----###----###----###----###----###----###----###----###----###----###----###----###

df[which(df$Federal_Provider_Number %in% "555904"), ]$lat <- 34.679361;
df[which(df$Federal_Provider_Number %in% "555904"), ]$lon <- -118.146952;
df[which(df$Federal_Provider_Number %in% "676427"), ]$lat <- 31.686175;
df[which(df$Federal_Provider_Number %in% "676427"), ]$lon <- -96.481509;

##The below four NHs do not have name or any address (i.e. city or state) information;
df <- df[-which(df$Federal_Provider_Number %in% "37E628"), ];
df <- df[-which(df$Federal_Provider_Number %in% "45F841"), ];
df <- df[-which(df$Federal_Provider_Number %in% "45F849"), ];
df <- df[-which(df$Federal_Provider_Number %in% "45F851"), ];

#df should be 77,669 rows x 61 columns;

#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#
##This from STEP_03_nursing_homes_map.R previously to re-geocode, shouldn't have to redo more than once;
df[which(df$Federal_Provider_Number %in% "475058"), ]$lat <- 42.817580;
df[which(df$Federal_Provider_Number %in% "475058"), ]$lon <- -96.551201;
df[which(df$Federal_Provider_Number %in% "375199"), ]$lat <- 45.320215;
df[which(df$Federal_Provider_Number %in% "375199"), ]$lon <- -96.445103;
df[which(df$Federal_Provider_Number %in% "245451"), ]$lat <- 35.100607;
df[which(df$Federal_Provider_Number %in% "245451"), ]$lon <- -96.405755;
df[which(df$Federal_Provider_Number %in% "165595"), ]$lat <- 43.943257;
df[which(df$Federal_Provider_Number %in% "165595"), ]$lon <- -72.607451;

df[which(df$Federal_Provider_Number %in% "225222"), ]$lat <- 33.575846;
df[which(df$Federal_Provider_Number %in% "225222"), ]$lon <- -86.414153;
df[which(df$Federal_Provider_Number %in% "195498"), ]$lat <- 40.162671;
df[which(df$Federal_Provider_Number %in% "195498"), ]$lon <- -103.213240;
df[which(df$Federal_Provider_Number %in% "115644"), ]$lat <- 33.251500;
df[which(df$Federal_Provider_Number %in% "115644"), ]$lon <- -82.234953;
df[which(df$Federal_Provider_Number %in% "065309"), ]$lat <- 29.959535;
df[which(df$Federal_Provider_Number %in% "065309"), ]$lon <- -91.039173;
df[which(df$Federal_Provider_Number %in% "015400"), ]$lat <- 42.305501;
df[which(df$Federal_Provider_Number %in% "015400"), ]$lon <- -71.296822;
df[which(df$Federal_Provider_Number %in% "655000"), ]$lat <- 13.504824;
df[which(df$Federal_Provider_Number %in% "655000"), ]$lon <- 144.774517;
#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#

# write out text datafile and an SAS program to read it
write.foreign(df, "nhsll_combined.txt", "nhsll_combined.sas", package="SAS");

o2 <- as.data.frame(df)
#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS("+proj=longlat")
o2new <- spTransform(o2, CRS(projected))

#par(mfrow=c(1,1), mar=c(1,1,1,11), oma=c(1,1,4,1), xpd=TRUE);
#plot(o2new);
#identify(x=coordinates(o2new)[,1], y=coordinates(o2new)[,2], lab=o2new$Federal_Provider_Number, tolerance = 0.05);
which(df$Federal_Provider_Number %in% c("225222", "195498", "115644", "065309", "015400", "655000"));
regeocode <- df[which(df$Federal_Provider_Number %in% c("225222", "195498", "115644", "065309", "015400", "655000")), ];
#as.data.frame(regeocode); #note the NAs, have to revise above to remerge properly...;
#MEADOWVIEW NURSING CENTER 7300 OLD HIGHWAY PELL CITY, AL 35128 "015400" 33.575846, -86.414153
#WASHINGTON COUNTY NURSING HOME 676 W GREENHOUSE DRIVE AKRON, CO 80720 "065309" 40.162671, -103.213240
#KEYSVILLE NURSING HOME & REHAB CENTER 1005 HWY KEYSVILLE, GA 30816  "115644"  33.251500, -82.234953
#ASSUMPTION HEALTHCARE AND REHABILITATION 252 HWY 402 NAPOLEONVILLE, LA 70390  "195498" 29.959535, -91.039173
#NEWTON WELLESLEY CENTER FOR ALZHEIMERS CARE 694 WORCESTER RD WELLESLEY FMS, MA 02181 "225222" 42.305501, -71.296822
#GUAM MEMORIAL HOSPITAL AUTHORITY 499 NORTH SABANA DRIVE BARRIGADA, GU 96913 "655000" 13.504824, 144.774517


