## Takes Medicare nursing home SAS data with geocoordinates and merges to 2010 Census at county and census tract levels, 2015 ACS at county level
## Tyler Pittman, April 2018
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_02_nursing_homes_US_Census_ACS_looping_April2018.R

##########################################################
#### DON'T RUN THIS SCRIPT ANYMORE, GOT DATA PREPARED ####
##########################################################

setwd("/Users/tylerpittman/GitHub/NursingHomesSNA/data");

library(tibble); 	## Must load before ggplot or haven;
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
library(RMySQL);
#library(spTimer);
library(imputation);
library(xtable);
library(spam);
library(haven); 	# Read in SAS or Stata files directly
library(foreign);	# Write SAS or Stata files directly
library(RJSONIO);	# Use Google API for lat and lon
library(plyr);		# Gives a nice join function for merges
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
library(spatialEco); 	#point.in.poly function for attribute data
library(SASxport);
library(data.table);	#gives fast fwrite and fread functions for .csv files
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function


ll = "+proj=longlat";
usCounties <- readShapeSpatial(paste("County_2010Census_DP1/County_2010Census_DP1", ".shp", sep=""), proj4string=CRS(ll));
usCounties.dbf <- read.dbf(paste("County_2010Census_DP1/County_2010Census_DP1", ".dbf", sep=""), header=TRUE);
usCensusTracts <- readShapeSpatial(paste("Tract_2010Census_DP1/Tract_2010Census_DP1", ".shp", sep=""), proj4string=CRS(ll));
usCensusTracts.dbf <- read.dbf(paste("Tract_2010Census_DP1/Tract_2010Census_DP1", ".dbf", sep=""), header=TRUE);
#proj4string(usCounties) <- CRS("+proj=longlat");

#projected <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
projected <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs";
usCountiesProj <- spTransform(usCounties, CRS(projected));
usCensusTractsProj <- spTransform(usCensusTracts, CRS(projected));
#plot(usCountiesProj);
#plot(usCensusTractsProj);

df <- read_sas("nhsLL.sas7bdat"); #uses tibbles automatically
df$lat;
df$lon;
nhNA <- df[which(is.na(df$lat)),];
o2 <- as.data.frame(df)
#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS("+proj=longlat");
o2new <- spTransform(o2, CRS(projected));


require(rgdal);
# The input file geodatabase
fgdb = "ACS_2015_5YR_COUNTY.gdb";
# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name));
fc_list = ogrListLayers(fgdb);
print(fc_list);
ogrInfo(dsn=fgdb, layer="ACS_2015_5YR_COUNTY");
# Read the feature class
fc = readOGR(dsn=fgdb,layer="ACS_2015_5YR_COUNTY", dropNULLGeometries=FALSE);
# Determine the FC extent, projection, and attribute information
summary(fc);
# View the feature class
#plot(fc);
colnames(fc@data);

##Have to use this method to read attribute data for geodatabase;
##which ogr2ogr true false  ##This gives location path of ogr2ogr on OSX;
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X17_POVERTY.csv ACS_2015_5YR_COUNTY.gdb X17_POVERTY");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X00_COUNTS.csv ACS_2015_5YR_COUNTY.gdb X00_COUNTS");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X19_INCOME.csv ACS_2015_5YR_COUNTY.gdb X19_INCOME");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X01_AGE_AND_SEX.csv ACS_2015_5YR_COUNTY.gdb X01_AGE_AND_SEX");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X02_RACE.csv ACS_2015_5YR_COUNTY.gdb X02_RACE");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X27_HEALTH_INSURANCE.csv ACS_2015_5YR_COUNTY.gdb X27_HEALTH_INSURANCE");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X10_GRANDPARENTS_GRANDCHILDREN.csv ACS_2015_5YR_COUNTY.gdb X10_GRANDPARENTS_GRANDCHILDREN");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X99_IMPUTATION.csv ACS_2015_5YR_COUNTY.gdb X99_IMPUTATION");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X18_DISABILITY.csv ACS_2015_5YR_COUNTY.gdb X18_DISABILITY");
#system("/Library/Frameworks/GDAL.framework/Programs/ogr2ogr -f CSV X07_MIGRATION.csv ACS_2015_5YR_COUNTY.gdb X07_MIGRATION");
mutext1 <- fread("X17_POVERTY.csv");
mutext2 <- fread("X00_COUNTS.csv");
mutext3 <- fread("X19_INCOME.csv");
mutext4 <- fread("X01_AGE_AND_SEX.csv");
mutext5 <- fread("X02_RACE.csv");
mutext6 <- fread("X27_HEALTH_INSURANCE.csv");
mutext7 <- fread("X10_GRANDPARENTS_GRANDCHILDREN.csv");
mutext8 <- fread("X99_IMPUTATION.csv");
mutext9 <- fread("X18_DISABILITY.csv");
mutext10 <- fread("X07_MIGRATION.csv");
#mutext1 <- read.csv("X17_POVERTY.csv");
#mutext2 <- read.csv("X00_COUNTS.csv");
#mutext3 <- read.csv("X19_INCOME.csv");
#mutext4 <- read.csv("X01_AGE_AND_SEX.csv");
#mutext5 <- read.csv("X02_RACE.csv");
#mutext6 <- read.csv("X27_HEALTH_INSURANCE.csv");
#mutext7 <- read.csv("X10_GRANDPARENTS_GRANDCHILDREN.csv");
#mutext8 <- read.csv("X99_IMPUTATION.csv");
#mutext9 <- read.csv("X18_DISABILITY.csv");
#mutext10 <- read.csv("X07_MIGRATION.csv");
row.names(mutext1) <- mutext1$GEOID;
row.names(mutext2) <- mutext2$GEOID;
row.names(mutext3) <- mutext3$GEOID;
row.names(mutext4) <- mutext4$GEOID;
row.names(mutext5) <- mutext5$GEOID;
row.names(mutext6) <- mutext6$GEOID;
row.names(mutext7) <- mutext7$GEOID;
row.names(mutext8) <- mutext8$GEOID;
row.names(mutext9) <- mutext9$GEOID;
row.names(mutext10) <- mutext10$GEOID;
acs <- join_all(list(mutext1, mutext2, mutext3, mutext4, mutext5, mutext6, mutext7, mutext8, mutext9, mutext10), by = 'GEOID', type = 'full');
#help.search("tibble");
acs <- as_tibble(acs);
acs$GEOID; #3220 counties;
colnames(usCountiesProj@data);
usCountiesProj@data$GEOID10; #3221 counties;
usCensusTractsProj@data$GEOID10; #74002 census tracts;
#sub("^[^S]*", "", acs$GEOID);
acs$GEOID10 <- substring(acs$GEOID, 8);
colnames(usCensusTractsProj@data) <- paste("CT", colnames(usCensusTractsProj@data), sep = "_");
usCountiesProj@data <- join_all(list(usCountiesProj@data, acs), by = 'GEOID10', type = "left");
#writePolyShape(usCountiesProj, "UScounties_ACS2015_Census2010_combined.shp");

nh.counties <- over(o2new, usCountiesProj);  
nh.comb.counties <- spCbind(o2new, nh.counties);
nh.countiesCTs <- over(nh.comb.counties, usCensusTractsProj);  
nh.comb.countiesCTs <- spCbind(nh.comb.counties, nh.countiesCTs);
str(nh.comb.countiesCTs@data); #Can't write "POSIXct", "POSIXt" to shapefile dbf, convert to Date class below;
nh.comb.countiesCTs@data[] <- lapply(nh.comb.countiesCTs@data, function(x) if(inherits(x, "POSIXct")) as.Date(x) else x);
##Below commands don't save formats, or are slow to read back into R;
##writePointsShape(nh.comb.countiesCTs, "Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.shp", factor2char=TRUE, max_nchar=254);
#write.csv(cbind(coordinates(nh.comb.countiesCTs), nh.comb.countiesCTs@data), file = "Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.csv", row.names = FALSE);

##Fwrite is much faster than write.csv as of 2016;
set.seed(1);
system.time(fwrite(cbind(coordinates(nh.comb.countiesCTs), nh.comb.countiesCTs@data), file = "Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.csv"));
#system.time(trialDF <- fread("Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.csv"));

##Use the SQLite format which stores long field names and can be read into QGIS or ArcGIS, but still doesn't work..;
#writeOGR(obj=nh.comb.countiesCTs, dsn="nhACScensusCT.sqlite", layer="nh", overwrite_layer=TRUE, driver="SQLite"); 
save(nh.comb.countiesCTs, file = "nhACScensusCT.RData");
#unlink("nhACScensusCT.RData");

## write out text datafile and an SAS program to read it, but too many formats to work export for 15,000 columns;
#write.foreign(nh.comb.countiesCTs@data, "Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.txt", "Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.sas", package="SAS");