## Takes Medicare nursing home SAS data owner groups;
## Tyler Pittman, May 2019
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_05_owner_groups_6May2019.R

#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#
## Order to run MedicareUS_scripts ##
## 1_nursing_homes_SNA_May2017.sas
## STEP_01_nursing_homes_geocode.R
## 1_nursing_homes_SNA_May2017.sas #Yes, run this again;
## 2_nursing_homes_SNA_May2017.sas
## STEP_02_nursing_homes_US_Census_ACS.R
## STEP_03_nursing_homes_map.R
#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#


#toBibtex(citation("ergm"));
#toBibtex(citation("RSiena"));

setwd("/Users/tylerpittman/GitHub/NursingHomesSNA/data");

library(tibble);	# load before ggplot or haven;
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
library(spatstat); 	## important gives point process modeling functions for kernel density
library(gstat); 	### Must load this after spatstat or will receive error message!!
library(RMySQL);
#library(spTimer);
###library(imputation);
library(xtable);
library(spam);
library(haven); 	# Read in SAS or Stata files directly
library(foreign);	# Write SAS or Stata files directly
library(RJSONIO);	# Use Google API for lat and lon
library(plyr);		# Gives a nice join function for merges
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
###library(spatialEco); 	#point.in.poly function for attribute data
library(SASxport);
library(data.table);	#gives fast fwrite and fread functions for .csv files
library(devEMF); 	# Allows plots to be saved as .emf nicely in Windows
library(ggplot2);
library(ggrepel);
#library(plotly);	# Makes nice pie charts
library(igraphdata);
library(igraph);
library(networkD3);
library(png); 
library(plyr); 		#gives multicolumn frequencies for tables with ddply();
library(extrafont); #Need to load this for png files to have correct font size on Cedar;
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function

##Only have to do this once;
#font_import();
#fonts();

##Do this to load fonts for plotting with MacOS;
#quartzFonts(avenir=c("Avenir Bookfamily="FreeSans"", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"));
#quartzFonts(ariel=c("Ariel Regular", "Ariel Italic", "Ariel Bold", "Ariel Bold Italic"));
#quartzFonts(helvetica=c("Helvetica Regular", "Helvetica Oblique", "Helvetica Bold", "Helvetica Bold Obilque"));
#quartzFonts(DejaVuSansMono=c("DejaVuSansMono", "DejaVuSansMono-Oblique", "DejaVuSansMono-Bold", "DejaVuSansMono-BoldOblique"));
#quartzFonts(garamond=c("GaramondNo8-Regular", "GaramondNo8-Italic", "GaramondNo8-Bold", "GaramondNo8-Bold-Italic"));
quartzFonts(FreeSans=c("FreeSans", "FreeSansBold", "FreeSansOblique", "FreeSansBoldOblique")); #must do this order on MacOS;
#par(family="FreeSans"); #do this after or during plot function call, the call to ‘quartzFonts’ above has to be executed first;

ll = "+proj=longlat";
usStates <- readShapeSpatial(paste("cb_2016_us_state_20m/cb_2016_us_state_20m", ".shp", sep=""), proj4string=CRS(ll));
usStates.dbf <- read.dbf(paste("cb_2016_us_state_20m/cb_2016_us_state_20m", ".dbf", sep=""), header=TRUE);
#usCounties <- readShapeSpatial(paste("County_2010Census_DP1/County_2010Census_DP1", ".shp", sep=""), proj4string=CRS(ll));
#usCounties.dbf <- read.dbf(paste("County_2010Census_DP1/County_2010Census_DP1", ".dbf", sep=""), header=TRUE);
#usCensusTracts <- readShapeSpatial(paste("Tract_2010Census_DP1/Tract_2010Census_DP1", ".shp", sep=""), proj4string=CRS(ll));
#usCensusTracts.dbf <- read.dbf(paste("Tract_2010Census_DP1/Tract_2010Census_DP1", ".dbf", sep=""), header=TRUE);
##proj4string(usCounties) <- CRS("+proj=longlat");
#projected <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs";
projected <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs";
#usCountiesProj <- spTransform(usCounties, CRS(projected));
#usCensusTractsProj <- spTransform(usCensusTracts, CRS(projected));
load("nhACScensusCT.RData");
length(nh.comb.countiesCTs); #should be 77760;
length(colnames(nh.comb.countiesCTs@data)); #should be 15305;
#nh.comb.countiesCTs@data$Processing_Date;
#nh.comb.countiesCTs$Processing_Date;

hrrBoundary <- readShapeSpatial(paste("hrr_bdry/hrr", ".shp", sep=""), proj4string=CRS("+proj=longlat +ellps=WGS84"));
hrrBoundary.dbf <- read.dbf(paste("hrr_bdry/hrr", ".dbf", sep=""), header=TRUE);
hrrBoundary <- spTransform(hrrBoundary, CRS(ll));

##Below is for cleaning Boundary file with regards to Erie, PA having correct hrrnum of 351 to lookup table;
#hrrBoundary@data[which(hrrBoundary@data$hrr_num == 305), ]$hrr_num <- 351;
#writePolyShape(hrrBoundary, paste("hrr_bdry/hrr", sep=""));

### This gives CMS HSA, PCSA and HHR look up tables based on zip code as of 2016 ###; 
lookup <- fread(paste("ZipHsaHrr16", ".csv", sep=""));
#colnames(lookup);
#str(lookup$zipcode16); #have to add leading 0 for zip code to make 5 digit length to merge;
lookup$zipcode <- sprintf("%05s", lookup$zipcode16); # fix to 5 characters;  

##Below is for cleaning lookup table with regards to Victoria, TX and Corpus Christi, TX in same boundary file;
#lookup[which(lookup$hrrnum == 417), ]$hrrcity <- "Corpus Christi";
#lookup[which(lookup$hrrnum == 417), ]$hrrnum <- 390;
#lookup[which(lookup$hrrnum == 390), ];
#fwrite(lookup, paste("ZipHsaHrr16", ".csv", sep=""));


##
## Check, all were correct;
#nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$PERIOD=="March2016"),]$Processing_Date <- "2016-03-01";
#nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$PERIOD=="September2016"),]$Processing_Date <- "2016-09-01";
#nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$PERIOD=="March2017"),]$Processing_Date <- "2017-03-01";
#nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$PERIOD=="September2017"),]$Processing_Date <- "2017-09-01";
#nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$PERIOD=="March2018"),]$Processing_Date <- "2018-03-01";
#save(nh.comb.countiesCTs, file = "nhACScensusCT.RData");
##

###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###;
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
#nh.comb.countiesCTs <- nh.comb.countiesCTs[c(which(nh.comb.countiesCTs@data$Provider_State == "MD")), ];
#nh.comb.countiesCTs@data$Provider_State;
#set.seed(1);
#system.time(fwrite(cbind(coordinates(nh.comb.countiesCTs), nh.comb.countiesCTs@data), file = "AnyLogicNhsAllCycles_MD.csv"));
###
###
###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###;


keyPeriod <- matrix(c( 
"March2016","2016-03-01",
"September2016","2016-09-01",
"March2017","2017-03-01",
"September2017","2017-09-01",
"March2018","2018-03-01"
), ncol=2, byrow=T);


####################### All Looping starts below here ###############################;

####################### Looping starts below here ###############################;
counterPeriod <- 5;
	periodLoop <- function(z) {
	k <- z;
	
	#k <- 1;
	keyDate = keyPeriod[k,1];

#df <- read_sas("nhACScensusCT.sas7bdat"); #uses tibbles automatically, SAS Studio can't export files this big properly so use .RData instead;
#colnames(df@data);
#df@data$Processing_Date;
#df <- nh.comb.countiesCTs;
df <- nh.comb.countiesCTs[which(nh.comb.countiesCTs$PERIOD == keyDate), ];
#length(df); 
#rm(nh.comb.countiesCTs); #unload memory intensive;
#which(is.na(df$Ownership_Type)); #none, good check;
#which(df$Ownership_Type == ""); #1749 are blank, not good, check;
#which(is.na(df$Processing_Date)); #none;

#colnames(df@data);
#table(df@data$Continuing_Care_Retirement_Commu);
#df@data$Legal_Business_Name;
#table(df@data$With_a_Resident_and_Family_Counc);
#table(df@data$Automatic_Sprinkler_Systems_in_A);

#df$lat;
#df$lon;
nhNA <- df[which(is.na(df$lat)),];
o2 <- as.data.frame(df)
#colnames(o2);
#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS(projected);

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
### This gives frequency of QOL Deficiency Severity Codes ### 
def <- fread(paste("Deficiencies_", keyDate, "_fixed.csv", sep=""));
#colnames(def);
#str(def$provnum); #have to add leading 0 for Federal_Provider_Number to make 6 digit length to merge;
def$Federal_Provider_Number <- sprintf("%06s", def$provnum); # fix to 6 characters; Some might have E+ and can't be linked; 
#def$filedate;
#def$scope;
facCount <- length(unique(def$Federal_Provider_Number));

### These 3 methods all get same freq count ###
countsM1 <- ddply(def, c("scope", "defstat"), nrow);
sumDef <- sum(countsM1$V1);
def$scopestat <- paste(def$scope, def$defstat, sep="_");
countsM2 <- ddply(def, .(def$scopestat), nrow);
countsM3 <- aggregate(rep(1, length(def$Federal_Provider_Number)), by=list(Category=def$scopestat), FUN=sum);
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 

#sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);
#usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#####################################################################

#nh.dat = read.csv("Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.csv", header = TRUE);
#nh.spdf <- SpatialPointsDataFrame(coords=c(which(colnames(nh.dat) %in% "lon"), which(colnames(nh.dat) %in% "lat")), data=data.frame(nh.dat));
##coordinates(nh.spdf) <- c(which(colnames(nh.spdf) %in% "lon"), which(colnames(nh.spdf) %in% "lat"));
#proj4string(nh.spdf) <- CRS("+proj=longlat");

## ESRI shapefile fields in dbf can only be 10 characters, below will truncate data;
##nh <- readOGR(".","Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined"); #this driver is slow
#nh <- readShapeSpatial(paste("Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined", ".shp", sep=""), proj4string=CRS(projected));
#nh.dbf <- read.dbf(paste("Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined", ".dbf", sep=""), header=TRUE);

###########################;

##http://www.shizukalab.com/toolkits/sna/node-level-calculations
reach2=function(x){
r=vector(length=vcount(x))
for (i in 1:vcount(x)){
n=neighborhood(x,2,nodes=i)
ni=unlist(n)
l=length(ni)
r[i]=(l)/vcount(x)}
r}

reach3=function(x){
r=vector(length=vcount(x))
for (i in 1:vcount(x)){
n=neighborhood(x,3,nodes=i)
ni=unlist(n)
l=length(ni)
r[i]=(l)/vcount(x)}
r}

dwreach=function(x){
distances=shortest.paths(x) #create matrix of geodesic distances
diag(distances)=1 # replace the diagonal with 1s
weights=1/distances # take the reciprocal of distances
apply(weights,1,sum) # sum for each node (row)
}

## Create adjacency matrix of Owner_Name by Federal_Provider_Number for nursing homes;
own <- fread(paste("Ownership_", keyDate, "_fixed.csv", sep=""));
#colnames(own);
#str(own$Federal_Provider_Number); #have to add leading 0 for Federal_Provider_Number to make 6 digit length to merge;
own$Federal_Provider_Number <- sprintf("%06s", own$Federal_Provider_Number); # fix to 6 characters; Some might have E+ and can't be linked; 
#own$Association_Date;
#own$Processing_Date;
#which(is.na(own$Processing_Date)); #none;
unique(own$Role_played_by_Owner_or_Manager);
# [1] "5% OR GREATER DIRECT OWNERSHIP INTEREST"  
# [2] "5% OR GREATER INDIRECT OWNERSHIP INTEREST"
# [3] "OPERATIONAL/MANAGERIAL CONTROL"           
# [4] "OFFICER"                                  
# [5] "MANAGING EMPLOYEE"                        
# [6] "5% OR GREATER SECURITY INTEREST"          
# [7] "DIRECTOR"                                 
# [8] "Ownership Data Not Available"             
# [9] "5% OR GREATER MORTGAGE INTEREST"          
#[10] "PARTNERSHIP INTEREST"  
#own$Owner_Name;
#own$Federal_Provider_Number;
fpntable <- table(own$Federal_Provider_Number);
otable <- table(own$Owner_Name);
ownSmall <- own;
#ownSmall <- own[c(1:200),];
#ownSmall  <- ownSmall[-c(which(is.na(ownSmall$Owner_Name))), ];
ownSmall  <- ownSmall[-c(which(ownSmall$Owner_Name == "")), ];
#ownSmall <- unique(ownSmall[ ,c("Owner_Name","Federal_Provider_Number")]);
##two_way <- as.data.frame(table(subset(ownSmall,select=c("Owner_Name","Federal_Provider_Number"))));
two_way_count <- count(ownSmall,c("Federal_Provider_Number", "Owner_Name"));
#which(is.na(two_way_count$Federal_Provider_Number)); #none;
#which(is.na(two_way_count$Owner_Name)); #none;
#which(is.na(two_way_count$freq)); #none;
set.seed(1);
system.time(fwrite(two_way_count, file = paste("Federal_Provider_Number_Owner_Name_2WAY_", keyDate, ".csv", sep="")));

##Processing_Date and other variables are merged twice below;
################
#edgelist <- cbind(two_way_count$Federal_Provider_Number, two_way_count$Owner_Name, two_way_count$freq);
#colnames(edgelist) <- c("Federal_Provider_Number", "Owner_Name", "freq");
#edgelist <- as.data.frame(edgelist);
#str(df@data$Federal_Provider_Number);
#str(two_way_count$Federal_Provider_Number);
#linkedDataNH <- df@data;
#linkedDataNH <- as.data.frame(linkedDataNH);
#edgelist0 <- join_all(list(edgelist, linkedDataNH), by = 'Federal_Provider_Number', type = "left", match = "first");
#edgelist00 <- join_all(list(edgelist0, own), by = 'Owner_Name', type = "left", match = "first");
#colnames(edgelist00); #check for duplicate column names;
#edgelist00 <- edgelist00[, -c(15308:15313, 15319)]; #This may change for each PERIOD but didn't for files used so far;
#colnames(edgelist00);
#edgelist00_tibble <- as_tibble(edgelist00);
#save(edgelist00, file = paste("edgelist00_", keyDate, ".RData", sep=""));
################
load(paste("edgelist00_", keyDate, ".RData", sep=""));

#Total_Fine_Amount
#Total_Amount_of_Fines_in_Dollars
edgelist <- edgelist00[,colnames(edgelist00) %in% c("Federal_Provider_Number", "Provider_Zip_Code", "Owner_Name", "Owner_Type", "freq", "Provider_State", "Total_Fine_Amount", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Association_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Processing_Date")]; 
edgelist[which(is.na(edgelist$Total_Fine_Amount)), ]$Total_Fine_Amount <- 0;

#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
edgelist <- edgelist[c(which(edgelist$Owner_Type == "Organization")), ];
edgelist <- edgelist[-c(which(edgelist$Provider_State %in% c(NA, "PR", "GU"))), ];
edgelist_count <- count(edgelist$Federal_Provider_Number);
number_NHs <- length(unique(edgelist$Federal_Provider_Number)); #11384;
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 255334),];
#edgelist[32173,];
#edgelist[which(edgelist$Provider_State %in% "AZ"),];
#edgelist[which(edgelist$Provider_State %in% "CO"),];
#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#


#edgelist$Provider_Zip_Code; #looks like all are 5 digit characters. Be certain below;
#colnames(lookup);
#str(edgelist$Provider_Zip_Code); 
#which(is.na(edgelist$Provider_Zip_Code)); #none;
edgelist$zipcode <- sprintf("%05s", edgelist$Provider_Zip_Code); # fix to 5 characters;  
#str(edgelist$zipcode); 


#tryCatch({
############################;
### 11 unlinked NHs, Google search and update with correct zip code to link;
###035216 075330 105178 105733 225222 225509 225674 225695 245622 365515 676196
#edgelist[which(edgelist$Federal_Provider_Number == "035216"), ]$zipcode <- "85147"; #85247 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "075330"), ]$zipcode <- "06828"; #06432 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "105178"), ]$zipcode <- "34102"; #33940 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "105733"), ]$zipcode <- "33770"; #34641 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "225222"), ]$zipcode <- "02481"; #02181 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "225509"), ]$zipcode <- "02445"; #02146 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "225674"), ]$zipcode <- "02492"; #02194 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "225695"), ]$zipcode <- "02492"; #02194 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "245622"), ]$zipcode <- "55084"; #55092 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "365515"), ]$zipcode <- "45402"; #45418 zipcode;
#edgelist[which(edgelist$Federal_Provider_Number == "676196"), ]$zipcode <- "76502"; #76505 zipcode;
####nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$Federal_Provider_Number == "676196"), ]$Provider_City;
#############################;
#}, error=function(e){})

#str(lookup$zipcode);
length(edgelist$zipcode); #44323;
#edgelist <- join_all(list(edgelist, lookup), by = 'zipcode', type = "left", match = "first"); 
edgelist <- join_all(list(edgelist, lookup), by = 'zipcode', type = "left"); 
length(edgelist$zipcode); #44323;

tryCatch({
#length(unique(edgelist[which(is.na(edgelist$hrrnum)), ]$Federal_Provider_Number)); #11 NHs can't be linked to HRRs;
edgelist <- edgelist[-c(which(is.na(edgelist$hrrnum))), ];
#which(is.na(edgelist$zipcode)); #none;
#edgelist$hrrnum;
#which(is.na(edgelist$hrrnum));
}, error=function(e){})

length(edgelist$zipcode); #44273;
#which(is.na(edgelist$zipcode)); #none;
#length(which(is.na(edgelist$hrrnum))); #0;
#unique(edgelist[which(is.na(edgelist$hrrnum)), ]$Federal_Provider_Number);
#unique(edgelist[which(is.na(edgelist$hrrnum)), ]$Provider_State);
#length(unique(edgelist[which(is.na(edgelist$hrrnum)), ]$Federal_Provider_Number)); #should be 0 now;
#edgelist[which(edgelist$hrrnum %in% "344"),]; 
#colnames(edgelist);

#save(edgelist, file = "edgelist.RData");
#load("edgelist.RData");
#edgelist <- as_tibble(edgelist);
#colnames(edgelist);
#length(edgelist$freq);


nh_states <- unique(edgelist$Provider_State);
length(nh_states); #51;

if( length(which(edgelist$Association_Date == "NO DATE PROVIDED")) > 0 ) {
    edgelist <- edgelist[-c(which(edgelist$Association_Date == "NO DATE PROVIDED")), ];
}

length(unique(edgelist$Federal_Provider_Number)); #11362 NHs;
length(edgelist$zipcode); #44262;
edgelist2 <- edgelist; #keep this for second loop for HRRs;


                
####################### Looping starts below here ###############################;
#counterState <- 1;
counterState <- length(nh_states);
stateLoop <- function(y) {
	j <- y;
	#j <- 10;
	#j <- 51;	

###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###
tryCatch({

st = nh_states[j];
stateNH = st; 
edgelist <- edgelist[which(edgelist$Provider_State == st), ];
#edgelist <- edgelist[which(edgelist$Provider_State == "MD"), ];

#edgelist$Owner_Name[grep("ONTARIO", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("BAY", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("REVERA", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("CHARTWELL", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("CHG", edgelist$Owner_Name)];

#edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING, LLC"),];
#edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING MANAGEMENT LLC"),];
#edgelist[which(edgelist$Owner_Name=="STRAND ADVISORS INC"),];
#edgelist[which(edgelist$Owner_Name=="COWAY"),];

#sel_fac <- edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING, LLC"),]$Federal_Provider_Number;
#sel_own <- edgelist[which(edgelist$Federal_Provider_Number %in% sel_fac),]$Owner_Name;
#edgelist[which(edgelist$Owner_Name %in% sel_own),];

#edgelist[which(edgelist$Owner_Name=="BAY BRIDGE CAPITAL PARTNERS, LLC"),];
#edgelist[which(edgelist$Federal_Provider_Number=="035068"),];
#edgelist[which(edgelist$Federal_Provider_Number=="055064"),];

#save(edgelist, file = "edgelist.RData");
#load("edgelist.RData");
#edgelist <- as_tibble(edgelist);
#nh_states <- nh_states[-c(which(is.na(nh_states)), which(nh_states == c("PR", "GU")))];
#
#if( length(which(nh_states == c("PR")) & which(nh_states == c("GU"))) > 0 ) {
#    nh_states <- nh_states[-c(which(nh_states == c("PR")), which(nh_states == c("GU")))];
#}
#
#if( length(which(edgelist$Association_Date == "NO DATE PROVIDED")) > 0 ) {
#    edgelist <- edgelist[-c(which(edgelist$Association_Date == "NO DATE PROVIDED")), ];
#}
#
#edgelist <- edgelist[-c(which(is.na(edgelist$Number_of_Residents_in_Certified))), ];
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
#which(is.na(edgelist$Ownership_Type)); #none, good check;


###
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
#edgelist$Provider_State;
igraph <- graph.data.frame(edgelist, directed=FALSE);
###
###
###


#igraph;
#igraph <- delete.vertices(igraph, which(is.na(V(igraph)$name)));

###let's remove multi-edges and loops; Don't use this need multi-edges;
##igraph <- simplify(igraph);
##E(igraph)$Total_Fine_Amount;
#classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="fisher");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="sd");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="quantile");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="hclust", method="complete");
brks <- c(0.0, 34457.0, 126503.5, 323213.0, 933453.0, 1508727.0);

##This works for whole US, but not with some states;
#E(igraph)$sdTFA <- cut(E(igraph)$Total_Fine_Amount, brks, right=FALSE);
##Below stopped working as of April 2018;
#E(igraph)$sdTFA <- log(E(igraph)$Total_Fine_Amount);
##E(igraph)$sdTFA <- cut(E(igraph)$Total_Fine_Amount, classTFA$brks, right=FALSE);
##quantile(edgelist$Total_Fine_Amount);

'%!in%' <- function(x,y)!('%in%'(x,y)) ;
V(igraph)$type <- V(igraph)$name %!in% edgelist[,2]; #either use %in% or %!in% operatote to get correct order of projections in bipartite.projection() below;
V(igraph)$color <- V(igraph)$type + 1;
V(igraph)$color = gsub("1", "red", V(igraph)$color) #Owners will be red;
V(igraph)$color = gsub("2", "blue", V(igraph)$color) #Nursing Homes will be blue;
#V(igraph)$size = 1;
V(igraph)$size = degree(igraph)*1 #because 1 is a small size for a node, multiply it by 5;
E(igraph)$weight <- as.numeric(edgelist[,4]);
#igraph;

tryCatch({
igraph.bi <- bipartite.projection(igraph, type=V(igraph)$type); #make sure NHs are proj1, and owners are proj2;
}, error=function(e){})

#-------------------------#;
key_fac <- matrix(c( 
"For profit - Corporation","#FEE5D9","circle",
"For profit - Individual","#FCAE91","circle",
"For profit - Limited Lia","#FB6A4A","circle",
"For profit - Partnership","#CB181D","circle",
"Government - City","#E0E8FF","square",
"Government - City/county","#C6DBEF","square",
"Government - County","#9ECAE1","square",
"Government - Federal","#6BAED6","square",
"Government - Hospital di","#3182BD","square",
"Government - State","#08519C","square",
"Non profit - Church rela","#E5F5E0","circle",
"Non profit - Corporation","#A1D99B","circle",
"Non profit - Other","#31A354","circle",
NA,"#FFFFFF","square"
), ncol=3, byrow=T);
key_fac <- as.data.frame(key_fac);
colnames(key_fac) <- c("Ownership_Type", "color", "shape");
key_own <- matrix(c( 
"Individual","#FFF44F","circle",
"Organization","#FF69B4","square"
), ncol=3, byrow=T);
key_own <- as.data.frame(key_own);
colnames(key_own) <- c("Owner_Type", "color", "shape");
edgelist$matchRetVal_fac = match(edgelist$Ownership_Type, key_fac$Ownership_Type);
edgelist$matchRetVal_own = match(edgelist$Owner_Type, key_own$Owner_Type);

edgelist.bi.proj2 <- edgelist[!duplicated(edgelist$Owner_Name), ];
length(edgelist.bi.proj2[,1]); #49632 unique owner names;
#edgelist.bi.proj2$Owner_Name[[49630]]; #should be LEADERSTAT LTD;
#V(igraph.bi$proj2)$name[[49630]]; #should be LEADERSTAT LTD;
V(igraph.bi$proj2)$name;
#set_vertex_attr(igraph.bi$proj2, name="type", index=V(igraph.bi$proj2), value=as.character(edgelist.bi.proj2$Owner_Type));
V(igraph.bi$proj2)$type <- as.character(edgelist.bi.proj2$Owner_Type);
V(igraph.bi$proj2)$type;
nodes_proj2 <- vcount(igraph.bi$proj2);
edges_proj2 <- ecount(igraph.bi$proj2);
stateNH <- edgelist.bi.proj2$Provider_State[[1]]; 
processingDate <- edgelist.bi.proj2$Processing_Date[[1]]; 
##
#edgelist.bi.proj2$Association_Date;
#which(is.na(edgelist.bi.proj2$Association_Date));
unlist(strsplit("since 01/01/2000", " "))[2];
unlist(strsplit("since 01/01/2000", " "))[c(FALSE, TRUE)];
ownAD <- sapply(strsplit(edgelist.bi.proj2$Association_Date, " "), `[`, 2);
ndatesAD <- as.Date(ownAD, "%m/%d/%Y"); 
procDate <- as.Date(processingDate, "%m/%d/%Y");
daysAD <- procDate - ndatesAD;
yearsAD <- daysAD/365.25;
sumyearsAD <- as.integer(sum(as.numeric(yearsAD), na.rm=TRUE));
V(igraph.bi$proj2)$size <- as.integer(yearsAD);
V(igraph.bi$proj2)$size;
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#igraph.bi$proj2 <- delete.vertices(igraph.bi$proj2, which(is.na(V(igraph.bi$proj2)$size)));
smallest_proj2 <- min(V(igraph.bi$proj2)$size);
largest_proj2 <- max(V(igraph.bi$proj2)$size);
tmpType2 = factor(edgelist.bi.proj2$Owner_Type, 
	levels = c("Individual", "Organization", NA));
V(igraph.bi$proj2)$color <- as.character(key_own$color[tmpType2]);
V(igraph.bi$proj2)$color;
V(igraph.bi$proj2)$shape <- as.character(key_own$shape[tmpType2]);
V(igraph.bi$proj2)$shape;
##
#V(igraph.bi$proj2)$color <- as.integer(factor(V(igraph.bi$proj2)$type));
#V(igraph.bi$proj2)$color <- key_own$color %*% edgelist.bi.proj2$matchRetVal_own;
#I(brewer.pal(nlevels(factor(V(igraph.bi$proj2)$type)), name = 'Dark2'));
#display.brewer.all(n=12)
#brewer.pal(1, "Greys");
#brewer.pal(1, "Oranges");

edgelist.bi.proj1 <- edgelist[!duplicated(edgelist$Federal_Provider_Number), ];
length(edgelist.bi.proj1[,1]); #14711 facilities;
#edgelist.bi.proj1$Federal_Provider_Number[[14705]]; #should be 675687;
#V(igraph.bi$proj1)$name[[14705]]; #should be 675687;
V(igraph.bi$proj1)$name;
#set_vertex_attr(igraph.bi$proj1, name="type", index=V(igraph.bi$proj1), value=as.character(edgelist.bi.proj1$Ownership_Type));
#V(igraph.bi$proj1)$type <- as.character(edgelist.bi.proj1$Ownership_Type);
tmpType1 = factor(edgelist.bi.proj1$Ownership_Type, 
	levels = c("For profit - Corporation", "For profit - Individual", 
	"For profit - Limited Lia", "For profit - Partnership",
	"Government - City", "Government - City/county",
	"Government - County", "Government - Federal",    
	"Government - Hospital di", "Government - State",      
	"Non profit - Church rela", "Non profit - Corporation",
	"Non profit - Other", NA));
V(igraph.bi$proj1)$color <- as.character(key_fac$color[tmpType1]);
V(igraph.bi$proj1)$color;
V(igraph.bi$proj1)$shape <- as.character(key_fac$shape[tmpType1]);
V(igraph.bi$proj1)$shape;
sumBeds_proj1 <- sum(edgelist.bi.proj1$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_proj1 <- sum(edgelist.bi.proj1$Number_of_Residents_in_Certified, na.rm=TRUE);
#V(igraph.bi$proj1)$color <- as.integer(factor(V(igraph.bi$proj1)$type));
#V(igraph.bi$proj1)$color <- I(brewer.pal(nlevels(factor(V(igraph.bi$proj1)$type)), name = 'Dark2'));

#igraph.bi$proj2 <- delete.vertices(igraph.bi$proj2, which(is.na(V(igraph.bi$proj2)$name)));


#igraph.bi$proj1[ order(-igraph.bi$proj1[,4], igraph.bi$proj1[,1]), ];
#display.brewer.all(n=6,type="seq",exact.n=TRUE);
#brewer.pal(4, "Reds");
#brewer.pal(6, "Blues");
#brewer.pal(3, "Greens");
##FFFFFF

#unique(V(igraph.bi$proj2)$type);
#[1] "Individual"   "Organization"
#unique(V(igraph.bi$proj1)$type);
# [1] "For profit - Corporation" "For profit - Individual" 
# [3] "For profit - Limited Lia" "For profit - Partnership"
# [5] "Government - City"        "Government - City/county"
# [7] "Government - County"      "Government - Federal"    
# [9] "Government - Hospital di" "Government - State"      
#[11] "Non profit - Church rela" "Non profit - Corporation"
#[13] "Non profit - Other"       NA    
#NHshape = c("circle", "circle", "circle", "circle", "square", "square", "square", "square", "square", "square", "circle", "circle", "circle", "square");
#NHcolor = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354", "#FFFFFF");


#V(igraph.bi$proj1)$size = degree(igraph.bi$proj1)*1;
#V(igraph.bi$proj1)$size = degree(igraph.bi$proj1)*1;
V(igraph.bi$proj1)$size <- edgelist.bi.proj1$Number_of_Residents_in_Certified;
V(igraph.bi$proj1)$size;
smallest_proj1 <- min(V(igraph.bi$proj1)$size);
largest_proj1 <- max(V(igraph.bi$proj1)$size);
#V(igraph.bi$proj1)$type %*% key_own$color;
nodes_proj1 <- vcount(igraph.bi$proj1);
edges_proj1 <- ecount(igraph.bi$proj1);
#edgelist.bi.proj1$Total_Fine_Amount;
#classTFA.proj1 <- classIntervals(edgelist.bi.proj1$Total_Fine_Amount, n=5, style="fisher");
#set_edge_attr(igraph.bi$proj1, "esize", index = E(igraph.bi$proj1), edgelist.bi.proj1$Number_of_Residents_in_Certified);
#E(igraph.bi$proj1)$sdTFA <- cut(edgelist.bi.proj1$Total_Fine_Amount, classTFA.proj1$brks, right=FALSE);
#save(igraph.bi, file = "igraph.bi.RData");
#load("igraph.bi.RData");
legend1_text <- c("For-profit -- Corporation", "For-profit -- Individual", 
	"For-profit -- Limited Liability company", "For-profit -- Partnership",
	"Government -- City", "Government -- City/county",
	"Government -- County", "Government -- Federal",    
	"Government -- Hospital district", "Government -- State",      
	"Non-profit -- Church related", "Non-profit -- Corporation",
	"Non-profit -- Other");
legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354");
legend1_symbols <- c(21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 21, 21, 21);
legend2_text <- c("Individual", "Organization");
legend2_colors <- c("#FFF44F", "#FF69B4");
legend2_symbols <- c(21, 22);
#-------------------------#;

igraph.edgelist <- get.edgelist(igraph);

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;

##https://rpubs.com/pjmurphy/323120
##Girvan-Newman edge betweenness;
#eb1 <- cluster_edge_betweenness(igraph.bi$proj1);
##export_pajek(eb1, igraph.bi$proj1, file=paste("owners_SNA_proj1_cluster_edge_betweenness_", keyDate, ".paj", sep=''));
#V(igraph.bi$proj1)$group <- cluster_edge_betweenness(igraph.bi$proj1) %>% membership();
#eb_comms1 <- cluster_edge_betweenness(igraph.bi$proj1);
###mark.col sets color of grouping boxes in plot options below;
###col sets color of vertex nodes like vertex.color did before;
##eb_comms1$membership = paste(stateNH, '_', eb_comms1$membership, sep='');
###sizes(eb_comms1);
##clusters_eb <- data.frame(eb_comms1$names, eb_comms1$membership);
##colnames(clusters_eb) <- c("Federal_Provider_Number", "Subgroup");
##clusters_eb;
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_cluster_edge_betweenness_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
##par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
#plot(eb_comms1, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##title("Predicted Response\nStandard Deviation", outer=F);
#par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n\nGirvan-Newman edge betweenness\nsubgroups"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
#options(scipen=999);
#options(scipen=999);
#mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
##mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
#dev.off();

##Spinglass clusters;
#sg1 <- cluster_spinglass(igraph.bi$proj1); #doesn't work with unaffiliated graph;

#Walktrap clusters;
#wc <- cluster_walktrap(igraph.bi$proj1);
#wc <- cluster_walktrap(igraph.bi$proj1, weights = E(igraph.bi$proj1)$weight, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE);
#modularity(wc);
#membership(wc);
#plot(wc, igraph.bi$proj1);
#plot(wc, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);


#----------------------------------#;
#16 October 2018 work;
##This places nodes into "affiliated" and "unaffiliated" organization ownership groups;
V(igraph.bi$proj1)$group <- cluster_edge_betweenness(igraph.bi$proj1) %>% membership();
V(igraph.bi$proj1)$group;

#must use [duplicated(vec) | duplicated(vec, fromLast=TRUE)] to get right groups;
#https://stackoverflow.com/questions/7854433/finding-all-duplicate-rows-including-elements-with-smaller-subscripts
#which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE));
#which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)));

V(igraph.bi$proj1)[which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE))]$Subgroup <- "Multiple";
V(igraph.bi$proj1)[which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)))]$Subgroup <- "Single";
V(igraph.bi$proj1)[which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE))]$membership <- 1;
V(igraph.bi$proj1)[which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)))]$membership <- 2;
#V(igraph.bi$proj1)$Subgroup;
#V(igraph.bi$proj1)$membership;
#names(V(igraph.bi$proj1));

clusters_simple <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup);
colnames(clusters_simple) <- c("Federal_Provider_Number", "Subgroup_Simple");
fwrite(clusters_simple, file = paste("NH_simple_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));

#--- November 2018 added;
V(igraph.bi$proj1)$degree <- degree(igraph.bi$proj1)*1;
degreeNotZero <- V(igraph.bi$proj1)$degree[which(V(igraph.bi$proj1)$degree != 0)];
meanDegreeMultiple <- sum(degreeNotZero) / length(degreeNotZero);
proportionMultiple <- length(degreeNotZero) / length(V(igraph.bi$proj1)$degree);

V(igraph.bi$proj1)$Subgroup = paste(stateNH, '_', V(igraph.bi$proj1)$Subgroup, sep='');
clusters_simple_state <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup, meanDegreeMultiple, proportionMultiple);
colnames(clusters_simple_state) <- c("Federal_Provider_Number", "Subgroup_Simple_State", "meanDegreeMultiple_State", "proportionMultiple_State");
#clusters_simple_state;
fwrite(clusters_simple_state, file = paste("NH_simple_state_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));
#---;

#----- May 2019 added for HHI calculation;
HHI_tmp1 <- (edgelist.bi.proj1$Number_of_Certified_Beds/sumBeds_proj1)^2;
HHI <- sum(HHI_tmp1, na.rm=TRUE);

V(igraph.bi$proj1)$Number_of_Certified_Beds <- edgelist.bi.proj1$Number_of_Certified_Beds;
clusters_simple <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup, V(igraph.bi$proj1)$group, V(igraph.bi$proj1)$Number_of_Certified_Beds);
colnames(clusters_simple) <- c("Federal_Provider_Number", "Subgroup_Simple", "Group", "Number_of_Certified_Beds");

HHI_tmp2 <- aggregate(clusters_simple$Number_of_Certified_Beds, by=list(Category=clusters_simple$Group), FUN=sum);
HHI_tmp3 <- (HHI_tmp2$x/sumBeds_proj1)^2;
affiliation_HHI <- sum(HHI_tmp3, na.rm=TRUE);

HHI <- round(HHI, digits=4);
affiliation_HHI <- round(affiliation_HHI, digits=4);
#-----;

#https://rdrr.io/cran/igraph/man/make_clusters.html
#simple1 <- make_clusters(igraph.bi$proj1, membership=V(igraph.bi$proj1)$membership, algorithm=V(igraph.bi$proj1)$membership, merges=NULL, modularity=FALSE);
simple1 <- make_clusters(igraph.bi$proj1, membership=V(igraph.bi$proj1)$membership, modularity=TRUE);
#length(membership(simple1)); 

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_cluster_simple_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,5.0,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(simple1, mark.col=c("#FFA700", "#DCD0FF"), mark.border=c("#FFA700", "#DCD0FF"), igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds\nMultiple and single class subgroups"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
options(scipen=999);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners\nwho are organizations"), adj=0);
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Herfindahl-Hirschman Index (HHI): ", HHI, "\n", "Affiliation Accounted HHI: ", affiliation_HHI, sep=""), cex=1.20, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();
#----------------------------------#;


##Louvain clusters; #uses modularity to detect clusters; https://arxiv.org/abs/0803.0476
lvc1 <- cluster_louvain(igraph.bi$proj1); 
#https://www.rdocumentation.org/packages/igraph/versions/1.2.2/topics/membership
#membership(lvc1);
#is_hierarchical(lvc1);
##mark.col sets color of grouping boxes in plot options below;
##col sets color of vertex nodes like vertex.color did before;
#modularity(lvc1);
#sizes(lvc1);

lvc1$membership = paste(stateNH, '_', lvc1$membership, sep='');
clusters_lvc1 <- data.frame(lvc1$names, lvc1$membership);
colnames(clusters_lvc1) <- c("Federal_Provider_Number", "Subgroup_Louvain");
#clusters_lvc1;
fwrite(clusters_lvc1, file = paste("NH_louvain_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_cluster_louvain_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,5.0,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(lvc1, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds\nLouvain cluster subgroups"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
options(scipen=999);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners\nwho are organizations"), adj=0);
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

##http://cneurocvs.rmki.kfki.hu/igraph/doc/R/cohesive.blocks.html
## Find cohesive blocks:
#gBlocks1 <- cohesive.blocks(igraph.bi$proj1);
##export_pajek(gBlocks1, igraph.bi$proj1, file=paste("owners_SNA_proj1_gBlocks_", keyDate, ".paj", sep=''));  
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_gBlocks_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
##par(xpd = T, mar = par()$mar + c(0,0,4,7))
##plot(igraph.bi$proj1);
#plot(gBlocks1, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##title("Predicted Response\nStandard Deviation", outer=F);
#par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n\nCohesive block grouping"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
#options(scipen=999);
#options(scipen=999);
#mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
##mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
#dev.off();

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;

##https://rpubs.com/pjmurphy/323120
##Girvan-Newman edge betweenness;
#eb <- cluster_edge_betweenness(igraph.bi$proj2);
#export_pajek(eb, igraph.bi$proj2, file=paste("owners_SNA_cluster_edge_betweenness_", keyDate, ".paj", sep=''));
#V(igraph.bi$proj2)$group <- cluster_edge_betweenness(igraph.bi$proj2) %>% membership();
#eb_comms <- cluster_edge_betweenness(igraph.bi$proj2);
##plot(eb_comms, igraph.bi$proj2, vertex.label = V(igraph.bi$proj2)$id, layout = layout_with_kk, main="Max Modularity Solution");
##mark.col sets color of grouping boxes in plot options below;
##col sets color of vertex nodes like vertex.color did before;
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_cluster_edge_betweenness_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
##par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
##plot(igraph.bi$proj2);
#plot(eb_comms, igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, col=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##title("Predicted Response\nStandard Deviation", outer=F);
#par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership\n\nGirvan-Newman edge betweenness\nsubgroups"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
#options(scipen=999);
#mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", stateNH, "facilities accepting Medicare or Medicaid"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
##mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
#dev.off();

##Spinglass clusters;
#sg <- cluster_spinglass(igraph.bi$proj2); #doesn't work with unaffiliated graph;

##Louvain clusters;
lvc <- cluster_louvain(igraph.bi$proj2); 
#https://www.rdocumentation.org/packages/igraph/versions/1.2.2/topics/membership
#membership(lvc);
#is_hierarchical(lvc);
#mark.col sets color of grouping boxes in plot options below;
#col sets color of vertex nodes like vertex.color did before;
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_cluster_louvain_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,5.0,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
plot(lvc, igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, col=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership\n\nLouvain cluster subgroups"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", stateNH, "facilities accepting Medicare or Medicaid\nwho are organizations"), adj=0);
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

##http://cneurocvs.rmki.kfki.hu/igraph/doc/R/cohesive.blocks.html
## Find cohesive blocks:
#gBlocks <- cohesive.blocks(igraph.bi$proj2);
##export_pajek(gBlocks, igraph.bi$proj2, file=paste("owners_SNA_gBlocks_", keyDate, ".paj", sep=''));
##blocks(gBlocks);
##cohesion(gBlocks);
##plot(gBlocks, igraph.bi$proj2, vertex.size=7, layout=layout.kamada.kawai);   
##mark.col sets color of grouping boxes in plot options below;
##col sets color of vertex nodes like vertex.color did before;
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_gBlocks_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
##par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
##plot(igraph.bi$proj2);
#plot(gBlocks, igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, col=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##title("Predicted Response\nStandard Deviation", outer=F);
#par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership\n\nCohesive block grouping"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
#options(scipen=999);
#mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", stateNH, "facilities accepting Medicare or Medicaid"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
##mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
#dev.off();

}, error=function(e){})
###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;
#-------------------------#;
#######################################################################
print(j);
j <- j + 1;

}
mclapply(1:counterState, stateLoop);
#mclapply(10, stateLoop);
#mclapply(5, stateLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;
#stopCluster(cl);




edgelist <- edgelist2; 
#edgelist[which(edgelist$hrrcity == "Minot"), ]$hrrnum; #none;
#edgelist[which(edgelist$hrrcity == "Albuquerque"), ]$hrrnum; #293;
#which(unique(edgelist$hrrnum) == 421); #287 in counter;
nh_HRRs <- unique(edgelist$hrrnum);
counterHRR <- length(nh_HRRs); #306 in latest 2016 file;

####################### Looping starts below here ###############################;

HRRLoop <- function(y) {
	c <- y;
	#c <- 10;
	#c <- 287;	

###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###
tryCatch({

hrr = nh_HRRs[c];
hrrNH = hrr; 
edgelist <- edgelist[which(edgelist$hrrnum == hrr), ];
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
#which(is.na(edgelist$Ownership_Type)); #none, good check;

#colnames(edgelist);
hrrCity <- unique(edgelist[which(edgelist$hrrnum == hrr), ]$hrrcity);
hrrState <- unique(edgelist[which(edgelist$hrrnum == hrr), ]$hrrstate);
#paste(hrr, ": ", hrrCity, ", ", hrrState, sep="");

###
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
#edgelist$Provider_State;
igraph <- graph.data.frame(edgelist, directed=FALSE);
###
###
###


#igraph;
#igraph <- delete.vertices(igraph, which(is.na(V(igraph)$name)));

###let's remove multi-edges and loops; Don't use this need multi-edges;
##igraph <- simplify(igraph);
##E(igraph)$Total_Fine_Amount;
#classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="fisher");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="sd");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="quantile");
##classTFA <- classIntervals(E(igraph)$Total_Fine_Amount, n=5, style="hclust", method="complete");
brks <- c(0.0, 34457.0, 126503.5, 323213.0, 933453.0, 1508727.0);

##This works for whole US, but not with some states;
#E(igraph)$sdTFA <- cut(E(igraph)$Total_Fine_Amount, brks, right=FALSE);
##Below stopped working as of April 2018;
#E(igraph)$sdTFA <- log(E(igraph)$Total_Fine_Amount);
##E(igraph)$sdTFA <- cut(E(igraph)$Total_Fine_Amount, classTFA$brks, right=FALSE);
##quantile(edgelist$Total_Fine_Amount);

'%!in%' <- function(x,y)!('%in%'(x,y)) ;
V(igraph)$type <- V(igraph)$name %!in% edgelist[,2]; #either use %in% or %!in% operatote to get correct order of projections in bipartite.projection() below;
V(igraph)$color <- V(igraph)$type + 1;
V(igraph)$color = gsub("1", "red", V(igraph)$color) #Owners will be red;
V(igraph)$color = gsub("2", "blue", V(igraph)$color) #Nursing Homes will be blue;
#V(igraph)$size = 1;
V(igraph)$size = degree(igraph)*1 #because 1 is a small size for a node, multiply it by 5;
E(igraph)$weight <- as.numeric(edgelist[,4]);
#igraph;

tryCatch({
igraph.bi <- bipartite.projection(igraph); #make sure NHs are proj1, and owners are proj2;
}, error=function(e){})

#-------------------------#;
key_fac <- matrix(c( 
"For profit - Corporation","#FEE5D9","circle",
"For profit - Individual","#FCAE91","circle",
"For profit - Limited Lia","#FB6A4A","circle",
"For profit - Partnership","#CB181D","circle",
"Government - City","#E0E8FF","square",
"Government - City/county","#C6DBEF","square",
"Government - County","#9ECAE1","square",
"Government - Federal","#6BAED6","square",
"Government - Hospital di","#3182BD","square",
"Government - State","#08519C","square",
"Non profit - Church rela","#E5F5E0","circle",
"Non profit - Corporation","#A1D99B","circle",
"Non profit - Other","#31A354","circle",
NA,"#FFFFFF","square"
), ncol=3, byrow=T);
key_fac <- as.data.frame(key_fac);
colnames(key_fac) <- c("Ownership_Type", "color", "shape");
key_own <- matrix(c( 
"Individual","#FFF44F","circle",
"Organization","#FF69B4","square"
), ncol=3, byrow=T);
key_own <- as.data.frame(key_own);
colnames(key_own) <- c("Owner_Type", "color", "shape");
edgelist$matchRetVal_fac = match(edgelist$Ownership_Type, key_fac$Ownership_Type);
edgelist$matchRetVal_own = match(edgelist$Owner_Type, key_own$Owner_Type);

edgelist.bi.proj2 <- edgelist[!duplicated(edgelist$Owner_Name), ];
length(edgelist.bi.proj2[,1]); #49632 unique owner names;
#edgelist.bi.proj2$Owner_Name[[49630]]; #should be LEADERSTAT LTD;
#V(igraph.bi$proj2)$name[[49630]]; #should be LEADERSTAT LTD;
V(igraph.bi$proj2)$name;
#set_vertex_attr(igraph.bi$proj2, name="type", index=V(igraph.bi$proj2), value=as.character(edgelist.bi.proj2$Owner_Type));
V(igraph.bi$proj2)$type <- as.character(edgelist.bi.proj2$Owner_Type);
V(igraph.bi$proj2)$type;
nodes_proj2 <- vcount(igraph.bi$proj2);
edges_proj2 <- ecount(igraph.bi$proj2);
hrrNH <- edgelist.bi.proj2$hrrnum[[1]]; 
processingDate <- edgelist.bi.proj2$Processing_Date[[1]]; 
##
#edgelist.bi.proj2$Association_Date;
#which(is.na(edgelist.bi.proj2$Association_Date));
unlist(strsplit("since 01/01/2000", " "))[2];
unlist(strsplit("since 01/01/2000", " "))[c(FALSE, TRUE)];
ownAD <- sapply(strsplit(edgelist.bi.proj2$Association_Date, " "), `[`, 2);
ndatesAD <- as.Date(ownAD, "%m/%d/%Y"); 
procDate <- as.Date(processingDate, "%m/%d/%Y");
daysAD <- procDate - ndatesAD;
yearsAD <- daysAD/365.25;
sumyearsAD <- as.integer(sum(as.numeric(yearsAD), na.rm=TRUE));
V(igraph.bi$proj2)$size <- as.integer(yearsAD);
V(igraph.bi$proj2)$size;
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#igraph.bi$proj2 <- delete.vertices(igraph.bi$proj2, which(is.na(V(igraph.bi$proj2)$size)));
smallest_proj2 <- min(V(igraph.bi$proj2)$size);
largest_proj2 <- max(V(igraph.bi$proj2)$size);
tmpType2 = factor(edgelist.bi.proj2$Owner_Type, 
	levels = c("Individual", "Organization", NA));
V(igraph.bi$proj2)$color <- as.character(key_own$color[tmpType2]);
V(igraph.bi$proj2)$color;
V(igraph.bi$proj2)$shape <- as.character(key_own$shape[tmpType2]);
V(igraph.bi$proj2)$shape;
##
#V(igraph.bi$proj2)$color <- as.integer(factor(V(igraph.bi$proj2)$type));
#V(igraph.bi$proj2)$color <- key_own$color %*% edgelist.bi.proj2$matchRetVal_own;
#I(brewer.pal(nlevels(factor(V(igraph.bi$proj2)$type)), name = 'Dark2'));
#display.brewer.all(n=12)
#brewer.pal(1, "Greys");
#brewer.pal(1, "Oranges");

edgelist.bi.proj1 <- edgelist[!duplicated(edgelist$Federal_Provider_Number), ];
length(edgelist.bi.proj1[,1]); #14711 facilities;
#edgelist.bi.proj1$Federal_Provider_Number[[14705]]; #should be 675687;
#V(igraph.bi$proj1)$name[[14705]]; #should be 675687;
V(igraph.bi$proj1)$name;
#set_vertex_attr(igraph.bi$proj1, name="type", index=V(igraph.bi$proj1), value=as.character(edgelist.bi.proj1$Ownership_Type));
#V(igraph.bi$proj1)$type <- as.character(edgelist.bi.proj1$Ownership_Type);
tmpType1 = factor(edgelist.bi.proj1$Ownership_Type, 
	levels = c("For profit - Corporation", "For profit - Individual", 
	"For profit - Limited Lia", "For profit - Partnership",
	"Government - City", "Government - City/county",
	"Government - County", "Government - Federal",    
	"Government - Hospital di", "Government - State",      
	"Non profit - Church rela", "Non profit - Corporation",
	"Non profit - Other", NA));
V(igraph.bi$proj1)$color <- as.character(key_fac$color[tmpType1]);
V(igraph.bi$proj1)$color;
V(igraph.bi$proj1)$shape <- as.character(key_fac$shape[tmpType1]);
V(igraph.bi$proj1)$shape;
sumBeds_proj1 <- sum(edgelist.bi.proj1$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_proj1 <- sum(edgelist.bi.proj1$Number_of_Residents_in_Certified, na.rm=TRUE);
#V(igraph.bi$proj1)$color <- as.integer(factor(V(igraph.bi$proj1)$type));
#V(igraph.bi$proj1)$color <- I(brewer.pal(nlevels(factor(V(igraph.bi$proj1)$type)), name = 'Dark2'));

#igraph.bi$proj2 <- delete.vertices(igraph.bi$proj2, which(is.na(V(igraph.bi$proj2)$name)));


#igraph.bi$proj1[ order(-igraph.bi$proj1[,4], igraph.bi$proj1[,1]), ];
#display.brewer.all(n=6,type="seq",exact.n=TRUE);
#brewer.pal(4, "Reds");
#brewer.pal(6, "Blues");
#brewer.pal(3, "Greens");
##FFFFFF

#unique(V(igraph.bi$proj2)$type);
#[1] "Individual"   "Organization"
#unique(V(igraph.bi$proj1)$type);
# [1] "For profit - Corporation" "For profit - Individual" 
# [3] "For profit - Limited Lia" "For profit - Partnership"
# [5] "Government - City"        "Government - City/county"
# [7] "Government - County"      "Government - Federal"    
# [9] "Government - Hospital di" "Government - State"      
#[11] "Non profit - Church rela" "Non profit - Corporation"
#[13] "Non profit - Other"       NA    
#NHshape = c("circle", "circle", "circle", "circle", "square", "square", "square", "square", "square", "square", "circle", "circle", "circle", "square");
#NHcolor = c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354", "#FFFFFF");


#V(igraph.bi$proj1)$size = degree(igraph.bi$proj1)*1;
#V(igraph.bi$proj1)$size = degree(igraph.bi$proj1)*1;
V(igraph.bi$proj1)$size <- edgelist.bi.proj1$Number_of_Residents_in_Certified;
V(igraph.bi$proj1)$size;
smallest_proj1 <- min(V(igraph.bi$proj1)$size);
largest_proj1 <- max(V(igraph.bi$proj1)$size);
#V(igraph.bi$proj1)$type %*% key_own$color;
nodes_proj1 <- vcount(igraph.bi$proj1);
edges_proj1 <- ecount(igraph.bi$proj1);
#edgelist.bi.proj1$Total_Fine_Amount;
#classTFA.proj1 <- classIntervals(edgelist.bi.proj1$Total_Fine_Amount, n=5, style="fisher");
#set_edge_attr(igraph.bi$proj1, "esize", index = E(igraph.bi$proj1), edgelist.bi.proj1$Number_of_Residents_in_Certified);
#E(igraph.bi$proj1)$sdTFA <- cut(edgelist.bi.proj1$Total_Fine_Amount, classTFA.proj1$brks, right=FALSE);
#save(igraph.bi, file = "igraph.bi.RData");
#load("igraph.bi.RData");
legend1_text <- c("For-profit -- Corporation", "For-profit -- Individual", 
	"For-profit -- Limited Liability company", "For-profit -- Partnership",
	"Government -- City", "Government -- City/county",
	"Government -- County", "Government -- Federal",    
	"Government -- Hospital district", "Government -- State",      
	"Non-profit -- Church related", "Non-profit -- Corporation",
	"Non-profit -- Other");
legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354");
legend1_symbols <- c(21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 21, 21, 21);
legend2_text <- c("Individual", "Organization");
legend2_colors <- c("#FFF44F", "#FF69B4");
legend2_symbols <- c(21, 22);
#-------------------------#;

igraph.edgelist <- get.edgelist(igraph);

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;

##https://rpubs.com/pjmurphy/323120
##Girvan-Newman edge betweenness;
#eb1 <- cluster_edge_betweenness(igraph.bi$proj1);
##export_pajek(eb1, igraph.bi$proj1, file=paste("owners_SNA_proj1_cluster_edge_betweenness_", keyDate, ".paj", sep=''));
#V(igraph.bi$proj1)$group <- cluster_edge_betweenness(igraph.bi$proj1) %>% membership();
#eb_comms1 <- cluster_edge_betweenness(igraph.bi$proj1);
###mark.col sets color of grouping boxes in plot options below;
###col sets color of vertex nodes like vertex.color did before;
##eb_comms1$membership = paste(stateNH, '_', eb_comms1$membership, sep='');
###sizes(eb_comms1);
##clusters_eb <- data.frame(eb_comms1$names, eb_comms1$membership);
##colnames(clusters_eb) <- c("Federal_Provider_Number", "Subgroup");
##clusters_eb;
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_cluster_edge_betweenness_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
##par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
#plot(eb_comms1, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##title("Predicted Response\nStandard Deviation", outer=F);
#par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n\nGirvan-Newman edge betweenness\nsubgroups"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
#options(scipen=999);
#options(scipen=999);
#mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
##mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
#dev.off();

##Spinglass clusters;
#sg1 <- cluster_spinglass(igraph.bi$proj1); #doesn't work with unaffiliated graph;

#Walktrap clusters;
#wc <- cluster_walktrap(igraph.bi$proj1);
#wc <- cluster_walktrap(igraph.bi$proj1, weights = E(igraph.bi$proj1)$weight, steps = 4, merges = TRUE, modularity = TRUE, membership = TRUE);
#modularity(wc);
#membership(wc);
#plot(wc, igraph.bi$proj1);
#plot(wc, igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);


#----------------------------------#;
#16 October 2018 work;
##This places nodes into "affiliated" and "unaffiliated" organization ownership groups;
V(igraph.bi$proj1)$group <- cluster_edge_betweenness(igraph.bi$proj1) %>% membership();
V(igraph.bi$proj1)$group;

#must use [duplicated(vec) | duplicated(vec, fromLast=TRUE)] to get right groups;
#https://stackoverflow.com/questions/7854433/finding-all-duplicate-rows-including-elements-with-smaller-subscripts
#which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE));
#which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)));

V(igraph.bi$proj1)[which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE))]$Subgroup <- "Multiple";
V(igraph.bi$proj1)[which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)))]$Subgroup <- "Single";
V(igraph.bi$proj1)[which(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE))]$membership <- 1;
V(igraph.bi$proj1)[which(!(duplicated(V(igraph.bi$proj1)$group) | duplicated(V(igraph.bi$proj1)$group, fromLast=TRUE)))]$membership <- 2;
#V(igraph.bi$proj1)$Subgroup;
#V(igraph.bi$proj1)$membership;
#names(V(igraph.bi$proj1));


clusters_simple <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup);
colnames(clusters_simple) <- c("Federal_Provider_Number", "Subgroup_Simple");
#fwrite(clusters_simple, file = paste("NH_simple_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));

#--- November 2018 added;
V(igraph.bi$proj1)$degree <- degree(igraph.bi$proj1)*1;
degreeNotZero <- V(igraph.bi$proj1)$degree[which(V(igraph.bi$proj1)$degree != 0)];
meanDegreeMultiple <- sum(degreeNotZero) / length(degreeNotZero);
proportionMultiple <- length(degreeNotZero) / length(V(igraph.bi$proj1)$degree);

V(igraph.bi$proj1)$Subgroup = paste(hrrNH, '_', V(igraph.bi$proj1)$Subgroup, sep='');
clusters_simple_hrr <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup, meanDegreeMultiple, proportionMultiple);
colnames(clusters_simple_hrr) <- c("Federal_Provider_Number", "Subgroup_Simple_HRR", "meanDegreeMultiple_HRR", "proportionMultiple_HRR");
#clusters_simple_hrr;
fwrite(clusters_simple_hrr, file = paste("NH_simple_HRR_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));
#---;

#----- 5 May 2019 added for HHI calculation;
HHI_tmp1 <- (edgelist.bi.proj1$Number_of_Certified_Beds/sumBeds_proj1)^2;
HHI <- sum(HHI_tmp1, na.rm=TRUE);

#edgelist$Owner_Name; #gives nice list of organization owner names;
V(igraph.bi$proj1)$Number_of_Certified_Beds <- edgelist.bi.proj1$Number_of_Certified_Beds;
clusters_simple <- data.frame(names(V(igraph.bi$proj1)), V(igraph.bi$proj1)$Subgroup, V(igraph.bi$proj1)$group, V(igraph.bi$proj1)$Number_of_Certified_Beds);
colnames(clusters_simple) <- c("Federal_Provider_Number", "Subgroup_Simple", "Group", "Number_of_Certified_Beds");

HHI_tmp2 <- aggregate(clusters_simple$Number_of_Certified_Beds, by=list(Category=clusters_simple$Group), FUN=sum);
HHI_tmp3 <- (HHI_tmp2$x/sumBeds_proj1)^2;
affiliation_HHI <- sum(HHI_tmp3, na.rm=TRUE);

HHI <- round(HHI, digits=4);
affiliation_HHI <- round(affiliation_HHI, digits=4);

clusters_simple_hrr <- cbind(clusters_simple_hrr, HHI, affiliation_HHI);
#clusters_simple_hrr;
fwrite(clusters_simple_hrr, file = paste("NH_simple_HRR_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));
#-----;

#---- 6 May 2019 added for Louvain clusters;
##Louvain clusters; #uses modularity to detect clusters; https://arxiv.org/abs/0803.0476
lvc1 <- cluster_louvain(igraph.bi$proj1); 

#-*-# 7 May 2019 added for number of affiliated organizations per HRR ;
sizeTmp <- as.numeric(as.character(V(igraph.bi$proj2)$size));
nameTmp <- as.character(V(igraph.bi$proj2)$name);
affiliatedOwners <- 1;
affiliatedOwners <- cbind(nameTmp, sizeTmp);
affiliatedOwners <- data.frame(affiliatedOwners);
colnames(affiliatedOwners) <- c("name", "size");
affiliatedOwners$size <- as.numeric(as.character(affiliatedOwners$size));
affiliatedOwners$name <- as.character(affiliatedOwners$name);
str(affiliatedOwners);
afOwners <- affiliatedOwners[which(affiliatedOwners$size > 1), ];
numberAffiliatedOwners <- length(unique(afOwners$name));
numberTiesAffiliatedOwners <- sum(afOwners$size);
#-*-#;

#https://www.rdocumentation.org/packages/igraph/versions/1.2.2/topics/membership
#membership(lvc1);
#is_hierarchical(lvc1);
##mark.col sets color of grouping boxes in plot options below;
##col sets color of vertex nodes like vertex.color did before;
#modularity(lvc1);
#sizes(lvc1);
lvc1$membership = paste(hrrNH, '_', lvc1$membership, sep='');
clusters_lvc1 <- data.frame(lvc1$names, lvc1$membership, numberAffiliatedOwners);
colnames(clusters_lvc1) <- c("Federal_Provider_Number", "Subgroup_Louvain_HRR", "Num_Affiliated_Owners");
#clusters_lvc1;

### Added 25 June 2019 ###;
smallData <- cbind(names(V(igraph.bi$proj1)), degree(igraph.bi$proj1)*1);
smallData <- as.data.frame(smallData);
colnames(smallData) <- c("Federal_Provider_Number", "Degree_Fac");
### figure out Degree_Own next.... can get this from NHC Ownership spreadsheet though and not of interest;
clusters_lvc1 <- join_all(list(clusters_lvc1, smallData), by="Federal_Provider_Number", type="left", match="all"); 
##########################;
fwrite(clusters_lvc1, file = paste("NH_louvain_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));
#-----;


#https://rdrr.io/cran/igraph/man/make_clusters.html
#simple1 <- make_clusters(igraph.bi$proj1, membership=V(igraph.bi$proj1)$membership, algorithm=V(igraph.bi$proj1)$membership, merges=NULL, modularity=FALSE);
simple1 <- make_clusters(igraph.bi$proj1, membership=V(igraph.bi$proj1)$membership, modularity=TRUE);
#length(membership(simple1)); 

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_HRR_Owners_proj1_SNA_cluster_simple_", hrrNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,5.0,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(simple1, mark.col=c("#FFA700", "#DCD0FF"), mark.border=c("#FFA700", "#DCD0FF"), igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, col=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds\nMultiple and single class subgroups"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n", y.intersp=0.0);
options(scipen=999);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Nursing homes affiliated to other facilities accepting\nMedicare or Medicaid through owners who are organizations", "\nin Hospital Referral Region of ", hrr, ": ", hrrCity, ", ", hrrState, ".", sep=""), adj=0);
mtext(text=paste("Herfindahl-Hirschman Index (HHI): ", HHI, "\n", "Affiliation Accounted HHI: ", affiliation_HHI, sep=""), cex=1.20, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();
#----------------------------------#;

}, error=function(e){})
###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;
#-------------------------#;
#######################################################################
print(c);
c <- c + 1;

}
mclapply(1:counterHRR, HRRLoop);
#mclapply(10, HRRLoop);
#mclapply(5, HRRLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;
#stopCluster(cl);


datalist = list();
for (i in 1:counterState) {
	stateNH = nh_states[i];
	dat <- fread(paste("NH_louvain_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));
    datalist[[i]] <- dat; # add it to your list;
    print(stateNH);
}
combined_lvc1 = do.call(rbind, datalist);
fwrite(combined_lvc1, file = paste("NH_louvain_ownership_clusters_", keyDate, ".csv", sep=""));
#11362 NHs;

datalist = list();
for (i in 1:counterState) {
	stateNH = nh_states[i];
	dat <- fread(paste("NH_simple_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));
    datalist[[i]] <- dat; # add it to your list;
    print(stateNH);
}
combined_simple1 = do.call(rbind, datalist);
fwrite(combined_simple1, file = paste("NH_simple_ownership_clusters_", keyDate, ".csv", sep=""));
#11362 NHs;

datalist = list();
for (i in 1:counterState) {
	stateNH = nh_states[i];
	dat <- fread(paste("NH_simple_state_ownership_clusters_", stateNH, "_", keyDate, ".csv", sep=""));
    datalist[[i]] <- dat; # add it to your list;
    print(stateNH);
}
combined_simple_state1 = do.call(rbind, datalist);
fwrite(combined_simple_state1, file = paste("NH_simple_state_ownership_clusters_", keyDate, ".csv", sep=""));
#11362 NHs; 

datalist = list();
for (i in 1:counterHRR) {
	hrrNH = nh_HRRs[i];
	dat <- fread(paste("NH_simple_HRR_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));
    datalist[[i]] <- dat; # add it to your list;
    print(hrrNH);
}
combined_simple_hrr1 = do.call(rbind, datalist);
fwrite(combined_simple_hrr1, file = paste("NH_simple_HRR_ownership_clusters_", keyDate, ".csv", sep=""));
#11362 NHs;

datalist = list();
for (i in 1:counterHRR) {
	hrrNH = nh_HRRs[i];
	dat <- fread(paste("NH_louvain_ownership_clusters_", hrrNH, "_", keyDate, ".csv", sep=""));
    datalist[[i]] <- dat; # add it to your list;
    print(hrrNH);
}
combined_lvc_hrr1 = do.call(rbind, datalist);
fwrite(combined_lvc_hrr1, file = paste("NH_louvain_ownership_clusters_HRR_", keyDate, ".csv", sep=""));
#11362 NHs;


#####################################################################################;
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#---------------------- Merge Owner Groups to Data ------------------------------#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#

##Convert to lat and lon from projected coordinates for spatial models;
#colnames(o2@data);
#o2@data$DP0150001; #Number of county households with individuals 65 years and over;
#coordinates(o2);
#proj4string(o2);
us_NHs = spTransform(o2, CRS("+proj=longlat"));
#colnames(us_NHs@data);
#coordinates(us_NHs);
us_NHs_small = us_NHs[, which(colnames(us_NHs@data) %in% c("Federal_Provider_Number", "Provider_Name", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Processing_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Total_Fine_Amount", "Association_Date", "Provider_SSA_County_Code", "Number_of_Certified_Beds", "Health_Inspection_Rating", "QM_Rating", "Staffing_Rating", "RN_Staffing_Rating", "Reported_CNA_Staffing_Hours_per", "Reported_LPN_Staffing_Hours_per", "Reported_RN_Staffing_Hours_per_R", "Reported_Licensed_Staffing_Hours", "Reported_Total_Nurse_Staffing_Ho", "Reported_Physical_Therapist_Staf", "Expected_CNA_Staffing_Hours_per", "Expected_LPN_Staffing_Hours_per", "Expected_RN_Staffing_Hours_per_R", "Expected_Total_Nurse_Staffing_Ho", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "Total_Weighted_Health_Survey_Sco", "Number_of_Facility_Reported_Inci", "Number_of_Substantiated_Complain", "Number_of_Fines", "Total_Amount_of_Fines_in_Dollars", "Number_of_Payment_Denials", "Total_Number_of_Penalties", "Special_Focus_Facility", "Provider_Changed_Ownership_in_La", "Provider_Resides_in_Hospital", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc"))];
colnames(us_NHs_small@data);
#us_NHs_small@data$Federal_Provider_Number;
#"Federal_Provider_Number", "Provider_Zip_Code", "Provider_State", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Processing_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Total_Fine_Amount", "Association_Date"
#"Provider_SSA_County_Code", "Number_of_Certified_Beds", "Health_Inspection_Rating", "QM_Rating", "Staffing_Rating", "RN_Staffing_Rating", "Reported_CNA_Staffing_Hours_per", "Reported_LPN_Staffing_Hours_per", "Reported_RN_Staffing_Hours_per_R", "Reported_Licensed_Staffing_Hours", "Reported_Total_Nurse_Staffing_Ho", "Reported_Physical_Therapist_Staf", "Expected_CNA_Staffing_Hours_per", "Expected_LPN_Staffing_Hours_per", "Expected_RN_Staffing_Hours_per_R", "Expected_Total_Nurse_Staffing_Ho", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "Total_Weighted_Health_Survey_Sco", "Number_of_Facility_Reported_Inci", "Number_of_Substantiated_Complain", "Number_of_Fines", "Total_Amount_of_Fines_in_Dollars", "Number_of_Payment_Denials", "Total_Number_of_Penalties", "Special_Focus_Facility", "Provider_Changed_Ownership_in_La", "Provider_Resides_in_Hospital", "Continuing_Care_Retirement_Commu"
#us_NHs_small[which(us_NHs_small$Provider_State %in% "AZ"),]; #148 NHs;
#us_NHs_small[which(us_NHs_small$Provider_State %in% "CO"),]; #222 NHs;

##############;
##must add leading 0's to five digit edgelist$Federal_Provider_Number to make 6 digit length to merge with us_NHs_small@data$Federal_Provider_Number;
#str(def$provnum); #have to add leading 0 for Federal_Provider_Number to make 6 digit length to merge;
combined_lvc1$Federal_Provider_Number <- sprintf("%06s", combined_lvc1$Federal_Provider_Number); # fix to 6 characters; 
combined_simple1$Federal_Provider_Number <- sprintf("%06s", combined_simple1$Federal_Provider_Number); # fix to 6 characters; 
combined_simple_state1$Federal_Provider_Number <- sprintf("%06s", combined_simple_state1$Federal_Provider_Number); # fix to 6 characters;   
combined_simple_hrr1$Federal_Provider_Number <- sprintf("%06s", combined_simple_hrr1$Federal_Provider_Number); # fix to 6 characters; 
combined_lvc_hrr1$Federal_Provider_Number <- sprintf("%06s", combined_lvc_hrr1$Federal_Provider_Number); # fix to 6 characters; 
##############;
#str(combined_lvc1);
#str(us_NHs_small@data);

selected_NHs <- join_all(list(us_NHs_small@data, combined_lvc1), by = 'Federal_Provider_Number', type = "right", match = "first"); 
#selected_NHs[which(selected_NHs$Provider_State %in% "AZ"),]; 
selected_NHs <- join_all(list(selected_NHs, combined_simple1), by = 'Federal_Provider_Number', type = "right", match = "first"); 
selected_NHs <- join_all(list(selected_NHs, combined_simple_state1), by = 'Federal_Provider_Number', type = "right", match = "first");
selected_NHs <- join_all(list(selected_NHs, combined_simple_hrr1), by = 'Federal_Provider_Number', type = "right", match = "first");
selected_NHs <- join_all(list(selected_NHs, combined_lvc_hrr1), by = 'Federal_Provider_Number', type = "right", match = "first");

#
#selected_NHs$zipcode <- sprintf("%05s", selected_NHs$Provider_Zip_Code); # fix to 5 characters;  
###selected_NHs$zipcode; 
#############################;
### 11 unlinked NHs, Google search and update with correct zip code to link;
###035216 075330 105178 105733 225222 225509 225674 225695 245622 365515 676196
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "035216"), ]$zipcode <- "85147"; #85247 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "075330"), ]$zipcode <- "06828"; #06432 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "105178"), ]$zipcode <- "34102"; #33940 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "105733"), ]$zipcode <- "33770"; #34641 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "225222"), ]$zipcode <- "02481"; #02181 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "225509"), ]$zipcode <- "02445"; #02146 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "225674"), ]$zipcode <- "02492"; #02194 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "225695"), ]$zipcode <- "02492"; #02194 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "245622"), ]$zipcode <- "55084"; #55092 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "365515"), ]$zipcode <- "45402"; #45418 zipcode;
#selected_NHs[which(selected_NHs$Federal_Provider_Number == "676196"), ]$zipcode <- "76502"; #76505 zipcode;
###nh.comb.countiesCTs@data[which(nh.comb.countiesCTs@data$Federal_Provider_Number == "676196"), ]$Provider_City;
#############################;
####which(is.na(selected_NHs$zipcode)); #none;
#selected_NHs <- join_all(list(selected_NHs, lookup), by = 'zipcode', type = "left", match = "first"); 
####length(selected_NHs$Federal_Provider_Number);
####which(is.na(selected_NHs$zipcode)); #none;
#

###selected_NHs[which(selected_NHs$Provider_State %in% "AZ"),]; 
length(selected_NHs$Federal_Provider_Number); #15664 NHs, should be 11373;
###selected_NHs$Federal_Provider_Number;
###selected_NHs[which(selected_NHs$Federal_Provider_Number == "01A193"), ]; #no Louvain cluster Subgroup, no organization owners?;
###selected_NHs[which(selected_NHs$Federal_Provider_Number == "075332"), ]; 
selected_NHs[which(is.na(selected_NHs$Total_Fine_Amount)), ]$Total_Fine_Amount <- 0;


#----------- ONLY KEEP NON-MISSING DATA, need this for statistical modeling part ahead ------------#;
#FOR EXPLORATIVE BOXPLOTS AND INFERENTIAL MODELING SECTIONS;
selected_NHs <- na.omit(selected_NHs);  #Reduces from 11072 organization owner NHs to 10728 NHs for March 2016, 15 November 2018;
#--------------------------------------------------------------------------------------------------#;

fwrite(selected_NHs, file = paste("NH_base_hierarchial_analysis_data_", keyDate, ".csv", sep=""));
  
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#


#################### HAVE TO EXECUTE ABOVE EACH TIME RUNNING THIS FILE ###################

print(k);
k <- k + 1;

}
mclapply(1:counterPeriod, periodLoop);
#mclapply(1, periodLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl1);