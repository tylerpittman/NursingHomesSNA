## Takes Medicare nursing home SAS data spatial hierarchical analysis;
## Tyler Pittman, June 2019
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_08_nursing_homes_longitudinal_15June2019.R

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

#Sys.setenv(R_MAX_NUM_DLLS = 512);
#Sys.getenv();

library(colorspace); #gives scale_fill_continuous_diverging() function, look at manual;
library(tibble);	# load before ggplot or haven;
#library(TeachingDemos);
library(sgeostat);
#library(car);
library(sp);
library(shapefiles);
#library(PBSmapping);
library(mapdata);
library(MASS);
library(maptools);
#library(spBayes); 	## important gives Bayesian spatial modeling functions
library(RColorBrewer);
#library(classInt);     	# finds class intervals for continuous variables
library(MBA);		# make sure to have `boost headers for C++' installed on linux system	
#library(gpclib);
library(rgeos);
#library(raster);
#library(fossil);
library(akima);
library(mapproj);
library(rgdal);
library(fields);
gpclibPermit() 		##need to use this to permit open license for unionSpatialPolygons feature
#library(RSurvey);	# requires x11 to be installed in ubuntu
#library(spatstat); 	## important gives point process modeling functions for kernel density
#library(gstat); 	### Must load this after spatstat or will receive error message!!
#library(RMySQL);
#library(spTimer);
###library(imputation);
#library(xtable);
#library(spam);
#library(haven); 	# Read in SAS or Stata files directly
#library(foreign);	# Write SAS or Stata files directly
library(RJSONIO);	# Use Google API for lat and lon
library(plyr);		# Gives a nice join function for merges
library(parallel);	# Does multicore processing when using mclapply WORKS on Linux!
###library(spatialEco); 	#point.in.poly function for attribute data
library(SASxport);
library(data.table);	#gives fast fwrite and fread functions for .csv files
#library(devEMF); 	# Allows plots to be saved as .emf nicely in Windows
#library(ggplot2);
library(ggrepel);
#library(plotly);	# Makes nice pie charts
#library(igraphdata);
#library(igraph);
#library(networkD3);
#library(png); 
library(plyr); 		#gives multicolumn frequencies for tables with ddply();
library(extrafont); #Need to load this for png files to have correct font size on Cedar;
#library(RSiena);
library(lme4);
library(lattice);
#library(VGAM);
library(MCMCglmm);
#library(btergm);
#library(networkDynamic);
library(viridis); 	#better color scales for ggplot
library(dplyr);  	#need dplyr for the bind_rows function
library(stargazer); 	#use to create nice summary table
library(tidyr);
###library(reshape2);	#gives melt() function
library(postMCMCglmm); #need for ranef() function, may not install on cedar;
library(reporttools);	#gives tableContinous() summary function;
library(gridExtra);		#gives grid.arrange() to stack ggplots;
library(doBy);			#gives summaryBy() function for descriptive statistics;
library(cowplot);		#gives fancy ways to arrange ggplots in one image;
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
##GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]
##epsg:4269 for most US Federal Agencies;
##use this online calculator http://prj2epsg.org/search
##https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf
usStates <- readShapeSpatial(paste("cb_2016_us_state_20m/cb_2016_us_state_20m", ".shp", sep=""), proj4string=CRS("+init=epsg:4269"));
#usStates <- readShapeSpatial(paste("cb_2016_us_state_20m/cb_2016_us_state_20m", ".shp", sep=""), proj4string=CRS(ll));
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

##https://data.cms.gov/dataset/Hospital-Referral-Region/ia25-mrsk
##click export data from above link and click shapefile;
##GEOGCS["WGS84(DD)", DATUM["WGS84", SPHEROID["WGS84", 6378137.0, 298.257223563]], PRIMEM["Greenwich", 0.0], UNIT["degree", 0.017453292519943295], AXIS["Geodetic longitude", EAST], AXIS["Geodetic latitude", NORTH]]
##epsg:4326 for most Global Agencies;
##use this online calculator http://prj2epsg.org/search
##https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/OverviewCoordinateReferenceSystems.pdf
#ogr2ogr -f "ESRI Shapefile" -t_srs EPSG:4269 -s_srs EPSG:4326 hrr_4269.shp hrr.shp
#https://github.com/nancyorgan/Hospital-Referral-Region-graphic/tree/master/hrr_bdry
#hrrBoundary <- readShapeSpatial(paste("hrr_bdry/hrr_4269", ".shp", sep=""), proj4string=CRS("+init=epsg:4269"));
hrrBoundary <- readShapeSpatial(paste("hrr_bdry/hrr", ".shp", sep=""), proj4string=CRS("+init=epsg:4326"));
#hrrBoundary <- readShapeSpatial(paste("hrr_bdry/hrr", ".shp", sep=""), proj4string=CRS("+proj=longlat +ellps=WGS84"));
hrrBoundary.dbf <- read.dbf(paste("hrr_bdry/hrr", ".dbf", sep=""), header=TRUE);

##Below is for cleaning Boundary file with regards to Erie, PA having correct hrrnum of 351 to lookup table;
#hrrBoundary@data[which(hrrBoundary@data$hrr_num == 305), ]$hrr_num <- 351;
#writePolyShape(hrrBoundary, paste("hrr_bdry/hrr", sep=""));


#### IMPORTANT FOR SHIFTING LATITUDE OR LONGITUDE !!!!
###An R option using package maptools and its elide function:
##shift y coordinates down 0.192 of a degree latitude to align;
##https://gis.stackexchange.com/questions/21256/how-to-easily-shift-all-features-in-a-vector-dataset/73162
#shift.xy <- c(0, -0.192);
### update the geometry with elide arguments
#shifted <- elide(hrrBoundary, shift = shift.xy);
### write out a new shapfile
#writeSpatialShape(shifted, paste("hrr_bdry/hrr", sep = ""));
###should align;
#plot(usStates, xlim=c(-110, -80), ylim=c(30, 50));
###plot(shifted, xlim=c(-110, -80), ylim=c(30, 50), border="red", add=TRUE);
#plot(hrrBoundary, xlim=c(-110, -80), ylim=c(30, 50), border="red", add=TRUE);



### This gives CMS HSA, PCSA and HHR look up tables based on zip code as of 2016 ###; 
lookup <- fread(paste("ZipHsaHrr16", ".csv", sep=""));
#colnames(lookup);
#str(lookup$zipcode16); #have to add leading 0 for zip code to make 5 digit length to merge;
lookup$zipcode <- sprintf("%05s", lookup$zipcode16); # fix to 5 characters;  

#http://www.gastonsanchez.com/visually-enforced/how-to/2014/01/15/Center-data-in-R/
# center with 'apply()'
# centering with 'scale()'
center_scale <- function(x) {
    scale(x, scale = FALSE)
}


###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###---###;

#-------------------------#;
## This does quick merge to see how many NHs linked longitudinally over all periods;
mdata1 <- fread(paste("ProviderInfo_", "March2016", "_fixed.csv", sep=""));
mdata2 <- fread(paste("ProviderInfo_", "September2016", "_fixed.csv", sep=""));
mdata3 <- fread(paste("ProviderInfo_", "March2017", "_fixed.csv", sep=""));
mdata4 <- fread(paste("ProviderInfo_", "September2017", "_fixed.csv", sep=""));
mdata5 <- fread(paste("ProviderInfo_", "March2018", "_fixed.csv", sep=""));
mdata1 <- mdata1[,c(1,81)];
mdata2 <- mdata2[,c(1,81)];
mdata3 <- mdata3[,c(1,81)];
mdata4 <- mdata4[,c(1,81)];
mdata5 <- mdata5[,c(1,81)];
colnames(mdata1) <- c("Federal_Provider_Number", "Processing_Date1");
colnames(mdata2) <- c("Federal_Provider_Number", "Processing_Date2");
colnames(mdata3) <- c("Federal_Provider_Number", "Processing_Date3");
colnames(mdata4) <- c("Federal_Provider_Number", "Processing_Date4");
colnames(mdata5) <- c("Federal_Provider_Number", "Processing_Date5");
mdata1$Federal_Provider_Number <- sprintf("%06s", mdata1$Federal_Provider_Number);  
mdata2$Federal_Provider_Number <- sprintf("%06s", mdata2$Federal_Provider_Number);  
mdata3$Federal_Provider_Number <- sprintf("%06s", mdata3$Federal_Provider_Number);  
mdata4$Federal_Provider_Number <- sprintf("%06s", mdata4$Federal_Provider_Number);  
mdata5$Federal_Provider_Number <- sprintf("%06s", mdata5$Federal_Provider_Number); 
#mdataAll <- join_all(list(mdata1, mdata2, mdata3, mdata4, mdata5), by="Federal_Provider_Number", type="inner", match="first"); 
mdataAll <- join_all(list(mdata1, mdata2, mdata3, mdata4, mdata5), by="Federal_Provider_Number", type="inner", match="all"); 
##15,264 NHs with data in each period;
#which(is.na(mdataAll$Processing_Date5)); #they all matched!
#-------------------------#;


#-------------------------#;
## This does quick merge to see how many NHs linked longitudinally over all periods with organization owners;
mydata1 <- fread(paste("NH_base_hierarchial_analysis_data_mds_", "March2016", ".csv", sep=""));
mydata2 <- fread(paste("NH_base_hierarchial_analysis_data_mds_", "September2016", ".csv", sep=""));
mydata3 <- fread(paste("NH_base_hierarchial_analysis_data_mds_", "March2017", ".csv", sep=""));
mydata4 <- fread(paste("NH_base_hierarchial_analysis_data_mds_", "September2017", ".csv", sep=""));
mydata5 <- fread(paste("NH_base_hierarchial_analysis_data_mds_", "March2018", ".csv", sep=""));
mydata1$Federal_Provider_Number <- sprintf("%06s", mydata1$Federal_Provider_Number);  
mydata2$Federal_Provider_Number <- sprintf("%06s", mydata2$Federal_Provider_Number);  
mydata3$Federal_Provider_Number <- sprintf("%06s", mydata3$Federal_Provider_Number);  
mydata4$Federal_Provider_Number <- sprintf("%06s", mydata4$Federal_Provider_Number);  
mydata5$Federal_Provider_Number <- sprintf("%06s", mydata5$Federal_Provider_Number); 
mydata1a <- mydata1[,c(1,14)];
mydata2a <- mydata2[,c(1,14)];
mydata3a <- mydata3[,c(1,14)];
mydata4a <- mydata4[,c(1,14)];
mydata5a <- mydata5[,c(1,14)];
colnames(mydata1a) <- c("Federal_Provider_Number", "Processing_Date1");
colnames(mydata2a) <- c("Federal_Provider_Number", "Processing_Date2");
colnames(mydata3a) <- c("Federal_Provider_Number", "Processing_Date3");
colnames(mydata4a) <- c("Federal_Provider_Number", "Processing_Date4");
colnames(mydata5a) <- c("Federal_Provider_Number", "Processing_Date5");
##mydataAlla <- join_all(list(mydata1a, mydata2a, mydata3a, mydata4a, mydata5a), by="Federal_Provider_Number", type="inner", match="first"); 
#mydataAlla <- join_all(list(mydata1a, mydata2a), by="Federal_Provider_Number", type="left"); 
mydataAlla <- join_all(list(mydata1a, mydata2a, mydata3a, mydata4a, mydata5a), by="Federal_Provider_Number", type="inner", match="all"); 
###9,001 NHs with complete data in each period with organization owner, most NA are due to missing MDS quality measures;
##which(is.na(mydataAlla$Processing_Date5)); #they all matched!
#-------------------------#;


#length(unique(mydata$HRR)); #295 HRRs with non missing data for NHs;
#mydata$QM_Rating;


#colnames(mydata1);
mydata <- rbind(mydata1, mydata2, mydata3, mydata4, mydata5); #54,291 nursing homes;
keyDate <- "March2016toMarch2018";
processingDate <- "2016-03-01 to 2018-03-01";

mydata$HRR <- sub("\\_.*", "", mydata$Subgroup_Simple_HRR);
mydata$Subgroup_Simple_Overall_by_State <- mydata$Subgroup_Simple;
mydata$Subgroup_Simple_Overall_by_HRR <- sub('.*_', '', mydata$Subgroup_Simple_HRR);
#str(mydata);

mydata$Occupancy_Rate <- (mydata$Number_of_Residents_in_Certified / mydata$Number_of_Certified_Beds); 
mydata$Occupancy_Rate <- round(mydata$Occupancy_Rate, digits = 3);
mydata[which(mydata$Occupancy_Rate > 1.0),]$Occupancy_Rate <- 1.0 ;
mydata$Occupancy_Rate <- mydata$Occupancy_Rate * 100;
summary(mydata$Occupancy_Rate); #mean is 80.73%;

mydata$Years_Business <- (as.Date(mydata$Processing_Date, "%Y-%m-%d") - as.Date(mydata$Date_First_Approved_to_Provide_M, "%Y-%m-%d"))/365.25; 

mydata[which(mydata$Subgroup_Simple_Overall_by_HRR == "Single"), ]$Degree_Fac <- 0; #just in case some with single affiliation are not zero;
#---------------------------#;


#----------- CONVERT TO INTEGER TO USE POISSON ------------#;
#mydata$Total_Weighted_Health_Survey_Sco <- as.integer(mydata$Total_Weighted_Health_Survey_Sco);
mydata$Total_Weighted_Health_Survey_Sco <- round(mydata$Total_Weighted_Health_Survey_Sco, digits=0);
#str(mydata$Total_Weighted_Health_Survey_Sco);
#---------------------------------------------------#;


summary(mydata$Total_Weighted_Health_Survey_Sco); #mean is total weighted deficiency score of 60.89;
sd(mydata$Total_Weighted_Health_Survey_Sco); #sd is 74.34;

summary(mydata$QM_Rating);
mydata$QM_Rating <- as.factor(mydata$QM_Rating);

summary(mydata$Subgroup_Simple);
mydata$Subgroup_Simple <- as.factor(mydata$Subgroup_Simple);
legend1_text <- levels(mydata$Subgroup_Simple);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,160), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR, Overall (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();


legend2_text <- c("Multiple  ", "Single  ", "Multiple  ", "Single  ", "Multiple  ", "Single  ", "Multiple  ", "Single  ", "Multiple  ", "Single  ");
#legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Processing_Date, sep=" - ", lex.order=FALSE)));
#fourCol <- c("grey", "white", "grey", "white", "grey", "white", "grey", "white", "grey", "white");
fourCol <- c('#b3e2cd', '#b3e2cd', '#fdcdac', '#fdcdac', '#cbd5e8', '#cbd5e8', '#f4cae4', '#f4cae4', '#e6f5c9', '#e6f5c9');
angles <- c(30, 0, 30, 0, 30, 0, 30, 0, 30, 0);
densities <- c(30, 0, 30, 0, 30, 0, 30, 0, 30, 0);

#unique(mydata$Processing_Date)[1];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_timeperiod_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple*Processing_Date, data=mydata, ylab="Total Weighted Health Survey Score", cex=1.5, xaxt="n", xlab="", ylim=c(-10,160), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=1.0);
text(x=1.5, y=-10, labels=unique(mydata$Processing_Date)[1]);
text(x=3.5, y=-10, labels=unique(mydata$Processing_Date)[2]);
text(x=5.5, y=-10, labels=unique(mydata$Processing_Date)[3]);
text(x=7.5, y=-10, labels=unique(mydata$Processing_Date)[4]);
text(x=9.5, y=-10, labels=unique(mydata$Processing_Date)[5]);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score in Nursing Homes by \nMultiple or Single Affiliation Class per HRR and Processing Date (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", " Processing date:", processingDate), cex=1.3, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;


#--------#;

# MUST MEAN CENTER BY TIME PERIOD!!!;

#--------#;



#--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--#;
#DESCRIPTIVE TABLE PART;
#--@--@--@--@--@--@--@--@--@;

#str(mydata);
#mydata$Years_Business;

#colnames(mydata);
#mydata$Subgroup_Simple_Overall_by_HRR;
#unique(mydata$Processing_Date);


#### Below works nicely for categorical variables, used for Excel input ####
#-------;
neat.table <- function(x, name){
  xx <- data.frame(x)
  names(xx) <- c("Value", "Count")
  xx$Fraction <- with(xx, Count/sum(Count))
  data.frame(Variable = name, xx)
}

x <- lapply(mydata[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata20160301 <- mydata[which(mydata$Processing_Date == "2016-03-01"), ];
x <- lapply(mydata20160301[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata20160901 <- mydata[which(mydata$Processing_Date == "2016-09-01"), ];
x <- lapply(mydata20160901[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata20170301 <- mydata[which(mydata$Processing_Date == "2017-03-01"), ];
x <- lapply(mydata20170301[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata20170901 <- mydata[which(mydata$Processing_Date == "2017-09-01"), ];
x <- lapply(mydata20170901[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydata20180301 <- mydata[which(mydata$Processing_Date == "2018-03-01"), ];
x <- lapply(mydata20180301[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "Processing_Date")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))


#x <- lapply(mydata, table)
#do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))
#-------;

#stargazer(ftable(mydata$Number_of_Certified_Beds_Cat, mydata$Provider_Resides_in_Hospital), type="text");
#stargazer(ftable(mydata$Number_of_Certified_Beds_Cat, mydata$Provider_Resides_in_Hospital));
## display default statistics, only use a subset of observations, grouped analysis
#tableContinuous(vars=mydata, prec=1, cap="Table of 'len','dose' by 'supp' ", lab="tab: descr stat");
#tableNominal(vars=mydata$Ownership_Type, cap="Table of nominal variables.", lab="tab: nominal", longtable=TRUE);

#stargazer(subset(mydata, Subgroup_Simple_Overall_by_HRR=="Single"), type="text", digits=2, title="Nursing Home Characteristics by Single Organization Ownership Class per HRR");
#stargazer(subset(mydata, Subgroup_Simple_Overall_by_HRR=="Multiple"), type="text", digits=2, title="Nursing Home Characteristics by Multiple Organization Ownership Class per HRR");
#stargazer(subset(mydata[, c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR")], mydata$Subgroup_Simple_Overall_by_HRR=="Multiple"), type="text", digits=2, title="Nursing Home Characteristics by Multiple Organization Ownership Class per HRR", order=c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR"));

#colnames(mydata);
mydata2 <- mydata[, c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Years_Business", "Subgroup_Simple_Overall_by_HRR", "Degree_Fac", "Processing_Date")];
#mydata2 <- mydata2[c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Subgroup_Simple_Overall_by_HRR")];
#colnames(mydata2);
#mydata2;
  
mydata2 %>%
    group_by(Subgroup_Simple_Overall_by_HRR, Processing_Date) %>%
    mutate(id = 1:n()) %>%
    ungroup() %>%
    gather(temp, val, Total_Weighted_Health_Survey_Sco, Number_of_Residents_in_Certified, Number_of_Certified_Beds, Occupancy_Rate, Adjusted_CNA_Staffing_Hours_per, Adjusted_LPN_Staffing_Hours_per, Adjusted_RN_Staffing_Hours_per_R, Adjusted_Total_Nurse_Staffing_Ho, Years_Business, proportionMultiple_HRR, meanDegreeMultiple_HRR, Degree_Fac) %>%
    unite(temp1, temp, Subgroup_Simple_Overall_by_HRR, Processing_Date, sep = '_') %>%
    spread(temp1, val) %>%
    select(-id) %>%
    as.data.frame() %>%
    stargazer(digits=2, type='text')
    

#https://stackoverflow.com/questions/43840086/total-of-group-by-summarize-values

# install.packages("reporttools")  #Use this to install it, do this only once
#require(reporttools);
vars <- mydata[, c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Years_Business", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Years_Business", "Degree_Fac")];
group <- interaction(as.factor(mydata$Subgroup_Simple_Overall_by_HRR), as.factor(mydata$Processing_Date)); #use this code for multiplying factor levels to use for grouping in tables;
#levels(group);
## display default statistics, only use a subset of observations, grouped analysis
result <- tableContinuous(vars=vars, group=group, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "n"), longtable=TRUE);
#cat(gsub("\\\\hline\n[^\n]+& all &[^\n]+\n", "", result)); #get rid of all category;
#cat(gsub("\\hline", "", result)); #get rid of \hline;


vars2 <- mydata[, c("mds_401", "mds_402", "mds_403", "mds_404", "mds_405", "mds_406", "mds_407", "mds_408", "mds_409", "mds_410", "mds_411", "mds_415", "mds_419", "mds_424", "mds_425", "mds_426", "mds_430", "mds_434", "mds_451", "mds_452", "mds_471")];
group2 <- interaction(as.factor(mydata$Subgroup_Simple_Overall_by_HRR), as.factor(mydata$Processing_Date)); 
## display default statistics, only use a subset of observations, grouped analysis
result2 <- tableContinuous(vars=vars2, group=group2, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);


#############################################################################################################
#### Might add this in later for descriptive tables of Nursing HRD and American State, but inconclusive #####
varsNurse <- mydata[, c("Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R")]; 
groupPS <- interaction(as.factor(mydata$Provider_State), as.factor(mydata$Processing_Date)); 
result3  <- tableContinuous(vars=varsNurse, group=groupPS, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);

sumHRDPS <- summaryBy(Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R ~ Provider_State + Subgroup_Simple_Overall_by_HRR, data=mydata, FUN = function(x) { c(m = mean(x), s = sd(x)) } );
#############################################################################################################


#-------------------------------#;
##### This does for HRR-level explanatory variables
#str(mydata);
#length(unique(mydata$HRR));

mydataHRR <- unique(mydata[, c("HRR", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "HHI", "affiliation_HHI", "Processing_Date")]); #should be 294 HRRs;
mydataHRR$delta_HHI <- mydataHRR$affiliation_HHI - mydataHRR$HHI; 
mydataHRR$proportionMultiple_HRR <- mydataHRR$proportionMultiple_HRR * 100;
vars3 <- mydataHRR[, c("proportionMultiple_HRR", "meanDegreeMultiple_HRR", "HHI", "affiliation_HHI", "delta_HHI")];
group3 <- as.factor(mydataHRR$Processing_Date);
result3 <- tableContinuous(vars=vars3, group=group3, prec=2, cap="Table of Nursing Home Random Effect Characteristics", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);
#####;
#-------------------------------#;




#-*-!-*-!-*-!-*-!-*-!-*-#;
### Takes this from part below for HRR-level descriptives of 9,001 NHs, done recursively 15June2019;
load(paste("NHs_Model_1_", keyDate, ".RData", sep=""));
mydataAlla1sample <- join_all(list(mydataAlla, mydata), by="Federal_Provider_Number", type="left"); 
mydataHRRsample <- unique(mydata[, c("HRR", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "HHI", "affiliation_HHI", "Processing_Date")]); #should be 294 HRRs;
mydataHRRsample$delta_HHI <- mydataHRRsample$affiliation_HHI - mydataHRRsample$HHI; 
mydataHRRsample$proportionMultiple_HRR <- mydataHRRsample$proportionMultiple_HRR * 100;
vars3sample <- mydataHRRsample[, c("proportionMultiple_HRR", "meanDegreeMultiple_HRR", "HHI", "affiliation_HHI", "delta_HHI")];
group3sample <- as.factor(mydataHRRsample$Processing_Date);
result3sample <- tableContinuous(vars=vars3sample, group=group3sample, prec=2, cap="Table of Nursing Home Random Effect Characteristics", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);

##Gives Total_Weighted_Health_Survey_Sco by Processing_Date;
mydataAlla1sample <- join_all(list(mydataAlla, mydata), by="Federal_Provider_Number", type="left"); 
#colnames(mydataAlla1sample);
varsTWHSsample <- mydataAlla1sample[, c("Degree_Fac", "Total_Weighted_Health_Survey_Sco")];
groupTWHSsample <- as.factor(mydataAlla1sample$Processing_Date);
resultTWHSsample <- tableContinuous(vars=varsTWHSsample, group=groupTWHSsample, prec=2, cap="Table of Nursing Home Outcome Characteristics", lab="tab2_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);
#-*-!-*-!-*-!-*-!-*-!-*-#;




### Create HRR-Level Explanatory Variables and link to NH-Level mydata ###
#mydataHRR$proportionMultiple_HRR;
#mydataHRR$meanDegreeMultiple_HRR;

## MUST CENTER SCALE (Mean-center) BY PROCESSING DATE!!!!!;


##http://r.789695.n4.nabble.com/How-to-convert-the-output-of-tapply-so-that-it-has-the-same-order-as-the-input-td4643238.html
##tapply(variabletoscale, list(groupvar1, groupvar2), scale)
##tmp1 <- ave(mydataHRR$proportionMultiple_HRR, list(as.factor(mydataHRR$Processing_Date)), center_scale); #this applies a Function Over A Ragged Array and works;
## BELOW DOES MEAN CENTERING CORRECTLY BY GROUP TO CHECK with group.center function further below ##
#tmp1 <- tapply(mydataHRR$proportionMultiple_HRR, list(as.factor(mydataHRR$Processing_Date)), center_scale); #this Apply A Function Over A Ragged Array and works;
##list(tmp1);
##unlist(tmp1[1]);
##https://stackoverflow.com/questions/35775696/trying-to-use-dplyr-to-group-by-and-apply-scale
#scaled_data <-  mydataHRR %>% group_by(Processing_Date, HRR) %>% mutate(C_proportionMultiple_HRR = center_scale(proportionMultiple_HRR));

#https://stats.stackexchange.com/questions/15305/how-to-group-center-standardize-variables-in-r
group.center <- function(var,grp) {
    return(var-tapply(var,grp,mean,na.rm=T)[grp])
}

#mydataHRR$proportionMultiple_HRRSc <- center_scale(mydataHRR$proportionMultiple_HRR); #old way of mean centering for one processing date;
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, proportionMultiple_HRRSc = center_scale(proportionMultiple_HRR)); ##This works for group using ddply function!!!;
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, meanDegreeMultiple_HRRSc = center_scale(meanDegreeMultiple_HRR)); ##This works for group using ddply function!!!;
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, HHI_HRRSc = center_scale(HHI)); ##This works for group using ddply function!!!;
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, affiliation_HHI_HRRSc = center_scale(affiliation_HHI)); ##This works for group using ddply function!!!;
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, delta_HHI_HRRSc = center_scale(delta_HHI)); ##This works for group using ddply function!!!;
#mydataHRR$proportionMultiple_HRRSc <- as.vector(group.center(mydataHRR$proportionMultiple_HRR, mydataHRR$Processing_Date)); #works for mean centering by group!;
#mydataHRR$meanDegreeMultiple_HRRSc <- as.vector(group.center(mydataHRR$meanDegreeMultiple_HRR, mydataHRR$Processing_Date)); #works for mean centering by group!;


mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, proportionMultiple_HRR_Quin = cut(proportionMultiple_HRR, breaks=c(quantile(proportionMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE));
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, meanDegreeMultiple_HRR_Quin = cut(meanDegreeMultiple_HRR, breaks=c(quantile(meanDegreeMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE));

mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, proportionMultiple_HRR_Bi = cut(proportionMultiple_HRR, breaks=c(quantile(proportionMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE));
mydataHRR <- ddply(mydataHRR, .(Processing_Date), transform, meanDegreeMultiple_HRR_Bi = cut(meanDegreeMultiple_HRR, breaks=c(quantile(meanDegreeMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE));

mydataHRR <- mydataHRR[, c("HRR", "Processing_Date", "proportionMultiple_HRRSc", "meanDegreeMultiple_HRRSc", "HHI_HRRSc", "affiliation_HHI_HRRSc", "delta_HHI_HRRSc", "proportionMultiple_HRR_Quin", "meanDegreeMultiple_HRR_Quin", "proportionMultiple_HRR_Bi", "meanDegreeMultiple_HRR_Bi")];
mydata <- join_all(list(mydata, mydataHRR), by = c('HRR', 'Processing_Date'), type = "left"); 
mydata$proportionMultiple_HRR_Quin <- as.factor(mydata$proportionMultiple_HRR_Quin);
mydata$meanDegreeMultiple_HRR_Quin <- as.factor(mydata$meanDegreeMultiple_HRR_Quin);
mydata$proportionMultiple_HRR_Bi <- as.factor(mydata$proportionMultiple_HRR_Bi);
mydata$meanDegreeMultiple_HRR_Bi <- as.factor(mydata$meanDegreeMultiple_HRR_Bi);
##################;


#col.tit=c("Mean", "Median", "S.D.", "Missing", "No.")

### These 3 methods all get same freq count ###
#countsM1 <- ddply(mydata, c("Ownership_Type", "Subgroup_Simple_Overall_by_HRR"), nrow, .drop=FALSE);

##### Added 6 May 2019;
#mydata$Subgroup_Louvain_HRR;
#mydata$HRR;


##mydataAffiliated <- mydata[which(duplicated(mydata$Subgroup_Louvain_HRR) == TRUE), ];
#mydataAffiliated <- ddply(mydata, .(Processing_Date), transform, mydata[which(duplicated(mydata$Subgroup_Louvain_HRR) == TRUE), ]);
#ddply(mydataAffiliated, .(Processing_Date), transform, length(unique(mydataAffiliated$Subgroup_Louvain_HRR)));

#unique(mydataAffiliated[,c('Subgroup_Louvain_HRR','Processing_Date')]);
#unique(mydataAffiliated$Processing_Date);

#### BE MINDFUL THAT NEW FACILITIES IN AFFILIATION GROUPS HAVE MISSING DATA FROM OSCAR CYCLES AND ARE NOT CAPTURED BELOW ####;
###"2016-03-01";
mydataAffiliatedPeriod <- mydata[which(mydata$Processing_Date == "2016-03-01"), ];
length(mydataAffiliatedPeriod$Federal_Provider_Number); #10,728 NHs;
length(mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$Subgroup_Simple_Overall_by_HRR == "Multiple"), ]$Federal_Provider_Number); #6,940 NHs;
mydataAffiliatedPeriod <- mydataAffiliatedPeriod[which(duplicated(mydataAffiliatedPeriod$Subgroup_Louvain_HRR) == TRUE), ];
mydataPeriod <- mydata[which(mydata$Processing_Date == "2016-03-01"), ];
length(unique(mydataAffiliatedPeriod$Subgroup_Louvain_HRR)); #1,870 removing commas and periods and multiple spaces;
#mydataPeriod$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydataPeriod[which(duplicated(mydataPeriod$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #17,202 removing commas and periods and multiple spaces;
checkAff <- mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #7 this is correct, removing commas and periods and multiple spaces;
###;
###"2016-09-01";
mydataAffiliatedPeriod <- mydata[which(mydata$Processing_Date == "2016-09-01"), ];
length(mydataAffiliatedPeriod$Federal_Provider_Number); #10,875 NHs;
length(mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$Subgroup_Simple_Overall_by_HRR == "Multiple"), ]$Federal_Provider_Number); #7,083 NHs;
mydataAffiliatedPeriod <- mydataAffiliatedPeriod[which(duplicated(mydataAffiliatedPeriod$Subgroup_Louvain_HRR) == TRUE), ];
mydataPeriod <- mydata[which(mydata$Processing_Date == "2016-09-01"), ];
length(unique(mydataAffiliatedPeriod$Subgroup_Louvain_HRR)); #1,896 removing commas and periods and multiple spaces;
#mydataPeriod$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydataPeriod[which(duplicated(mydataPeriod$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #17,289 removing commas and periods and multiple spaces;
checkAff <- mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #9 this is correct, removing commas and periods and multiple spaces;
###;
###"2017-03-01";
mydataAffiliatedPeriod <- mydata[which(mydata$Processing_Date == "2017-03-01"), ];
length(mydataAffiliatedPeriod$Federal_Provider_Number); #11,009 NHs;
length(mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$Subgroup_Simple_Overall_by_HRR == "Multiple"), ]$Federal_Provider_Number); #7,191 NHs;
mydataAffiliatedPeriod <- mydataAffiliatedPeriod[which(duplicated(mydataAffiliatedPeriod$Subgroup_Louvain_HRR) == TRUE), ];
mydataPeriod <- mydata[which(mydata$Processing_Date == "2017-03-01"), ];
length(unique(mydataAffiliatedPeriod$Subgroup_Louvain_HRR)); #1,905 removing commas and periods and multiple spaces;
#mydataPeriod$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydataPeriod[which(duplicated(mydataPeriod$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #18,847 removing commas and periods and multiple spaces;
checkAff <- mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #9 this is correct, removing commas and periods and multiple spaces;
###;
###"2017-09-01";
mydataAffiliatedPeriod <- mydata[which(mydata$Processing_Date == "2017-09-01"), ];
length(mydataAffiliatedPeriod$Federal_Provider_Number); #10,736 NHs;
length(mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$Subgroup_Simple_Overall_by_HRR == "Multiple"), ]$Federal_Provider_Number); #7,058 NHs;
mydataAffiliatedPeriod <- mydataAffiliatedPeriod[which(duplicated(mydataAffiliatedPeriod$Subgroup_Louvain_HRR) == TRUE), ];
mydataPeriod <- mydata[which(mydata$Processing_Date == "2017-09-01"), ];
length(unique(mydataAffiliatedPeriod$Subgroup_Louvain_HRR)); #1,865 removing commas and periods and multiple spaces;
#mydataPeriod$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydataPeriod[which(duplicated(mydataPeriod$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #18,271 removing commas and periods and multiple spaces;
checkAff <- mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #10 this is correct, removing commas and periods and multiple spaces;
###;
###"2018-03-01";
mydataAffiliatedPeriod <- mydata[which(mydata$Processing_Date == "2018-03-01"), ];
length(mydataAffiliatedPeriod$Federal_Provider_Number); #10,943 NHs;
length(mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$Subgroup_Simple_Overall_by_HRR == "Multiple"), ]$Federal_Provider_Number); #7,146 NHs;
mydataAffiliatedPeriod <- mydataAffiliatedPeriod[which(duplicated(mydataAffiliatedPeriod$Subgroup_Louvain_HRR) == TRUE), ];
mydataPeriod <- mydata[which(mydata$Processing_Date == "2018-03-01"), ];
length(unique(mydataAffiliatedPeriod$Subgroup_Louvain_HRR)); #1,878 removing commas and periods and multiple spaces;
#mydataPeriod$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydataPeriod[which(duplicated(mydataPeriod$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #18,839 removing commas and periods and multiple spaces;
checkAff <- mydataAffiliatedPeriod[which(mydataAffiliatedPeriod$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #8 this is correct, removing commas and periods and multiple spaces;
###;


#--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--#;
#ANALYTICAL TABLE PART;
#--@--@--@--@--@--@--@--@--@;

#unique(mydata$Ownership_Type);
#summary(mydata$Ownership_Type);
legend1_text <- c("For-profit -- \nCorporation", "For-profit -- \nIndividual", 
	"For-profit -- \nLimited Liability company", "For-profit -- \nPartnership",
	"Government -- \nCity", "Government -- \nCity/county",
	"Government -- \nCounty", "Government -- \nFederal",    
	"Government -- \nHospital district", "Government -- \nState",      
	"Non-profit -- \nChurch related", "Non-profit -- \nCorporation",
	"Non-profit -- \nOther");
legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354");

mydata$Ownership_Type <- as.factor(mydata$Ownership_Type);
levels(mydata$Ownership_Type) <- legend1_text;

#### Collapse categories;
legend1_text <- c("For-profit", "Government", "Non-profit");
mydata$Ownership_Type <- as.character(mydata$Ownership_Type);
mydata[which(mydata$Ownership_Type == "For-profit -- \nCorporation"), ]$Ownership_Type <- "For-profit";
mydata[which(mydata$Ownership_Type == "For-profit -- \nIndividual"), ]$Ownership_Type <- "For-profit";
mydata[which(mydata$Ownership_Type == "For-profit -- \nLimited Liability company"), ]$Ownership_Type <- "For-profit";
mydata[which(mydata$Ownership_Type == "For-profit -- \nPartnership"), ]$Ownership_Type <- "For-profit";
mydata[which(mydata$Ownership_Type == "Government -- \nCity"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Government -- \nCity/county"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Government -- \nCounty"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Government -- \nFederal"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Government -- \nHospital district"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Government -- \nState"), ]$Ownership_Type <- "Government";
mydata[which(mydata$Ownership_Type == "Non-profit -- \nChurch related"), ]$Ownership_Type <- "Non-profit";
mydata[which(mydata$Ownership_Type == "Non-profit -- \nCorporation"), ]$Ownership_Type <- "Non-profit";
mydata[which(mydata$Ownership_Type == "Non-profit -- \nOther"), ]$Ownership_Type <- "Non-profit";
mydata$Ownership_Type <- as.factor(mydata$Ownership_Type);
levels(mydata$Ownership_Type) <- legend1_text;
summary(mydata$Ownership_Type);

legend1_text <- c("No", "Yes");
mydata$With_a_Resident_and_Family_Counc <- as.character(mydata$With_a_Resident_and_Family_Counc);
mydata[which(mydata$With_a_Resident_and_Family_Counc == "Both"), ]$With_a_Resident_and_Family_Counc <- "Yes";
mydata[which(mydata$With_a_Resident_and_Family_Counc == "None"), ]$With_a_Resident_and_Family_Counc <- "No";
mydata[which(mydata$With_a_Resident_and_Family_Counc == "Resident"), ]$With_a_Resident_and_Family_Counc <- "Yes";
mydata[which(mydata$With_a_Resident_and_Family_Counc == "Family"), ]$With_a_Resident_and_Family_Counc <- "Yes";
mydata$With_a_Resident_and_Family_Counc <- as.factor(mydata$With_a_Resident_and_Family_Counc);
levels(mydata$With_a_Resident_and_Family_Counc) <- legend1_text;
summary(mydata$With_a_Resident_and_Family_Counc);

####; 


#http://www.gastonsanchez.com/visually-enforced/how-to/2014/01/15/Center-data-in-R/
# center with 'apply()'
# centering with 'scale()'
center_scale <- function(x) {
    scale(x, scale = FALSE)
}

#Adjusted_RN_Staffing_Hours_per_R_mean <- mean(mydata$Adjusted_RN_Staffing_Hours_per_R);
mydata <- ddply(mydata, .(Processing_Date), transform, Adjusted_RN_Staffing_Hours_per_R = center_scale(Adjusted_RN_Staffing_Hours_per_R)); 
mydata <- ddply(mydata, .(Processing_Date), transform, Adjusted_LPN_Staffing_Hours_per = center_scale(Adjusted_LPN_Staffing_Hours_per)); 
mydata <- ddply(mydata, .(Processing_Date), transform, Adjusted_CNA_Staffing_Hours_per = center_scale(Adjusted_CNA_Staffing_Hours_per)); 

mydata <- ddply(mydata, .(Processing_Date), transform, Occupancy_Rate = center_scale(Occupancy_Rate)); 

mydata <- ddply(mydata, .(Processing_Date), transform, mds_401 = center_scale(mds_401)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_402 = center_scale(mds_402)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_403 = center_scale(mds_403)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_404 = center_scale(mds_404)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_405 = center_scale(mds_405)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_406 = center_scale(mds_406)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_407 = center_scale(mds_407)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_408 = center_scale(mds_408)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_409 = center_scale(mds_409)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_410 = center_scale(mds_410)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_411 = center_scale(mds_411)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_415 = center_scale(mds_415)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_419 = center_scale(mds_419)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_424 = center_scale(mds_424)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_425 = center_scale(mds_425)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_426 = center_scale(mds_426)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_430 = center_scale(mds_430)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_434 = center_scale(mds_434)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_451 = center_scale(mds_451)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_452 = center_scale(mds_452)); 
mydata <- ddply(mydata, .(Processing_Date), transform, mds_471 = center_scale(mds_471)); 


###mydata$proportionMultiple_HRR <- mydata$proportionMultiple_HRR * 100;
### Don't mean transform here, must do per HRR in above section;
##mydata$proportionMultiple_HRR <- center_scale(mydata$proportionMultiple_HRR); #64.45131 was mean;
##mydata$meanDegreeMultiple_HRR <- center_scale(mydata$meanDegreeMultiple_HRR); #4.562907 was mean;

# below MCMC function is not yet implemented for `glmer'
#jpeg("manitobaAsthmaLP15posteriorExample.jpeg", bg="white", res=150, units="in", width=11, height=11);
#pvals.fnc(manitobaAsthmaLP15.lmer, nsim=10000);
#dev.off()

# below MCMC function is not yet implemented for `glmer'
#manitobaAsthmaLP15.lmerExact <- lmer(asth_child ~ 1 + fh_asthma + LF_PAR_15P + (1|NEIGH/CTUID),data=manitobaAsthmaHLM1996);
#postscript("manitobaAsthmaLP15posteriorExample.eps", paper="special", width=10, height=10, horizontal=FALSE);
#jpeg("manitobaAsthmaLP15posteriorExample.jpeg", bg="white", res=150, units="in", width=11, height=11);
#pvals.fnc(manitobaAsthmaLP15.lmerExact, nsim=40000);
#dev.off()

#fit <- glmer(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Adjusted_RN_Staffing_Hours_per_R + (1|Provider_State/Subgroup_Simple_State), family=negative.binomial(2), data=mydata); #theta = 2;
#summary(fit);
#pvals.fnc(fit, nsim=40000);

#prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)));
summary(mydata$Total_Weighted_Health_Survey_Sco); #overall mean in total weighted deficiency score is 60.89;
sd(mydata$Total_Weighted_Health_Survey_Sco); #overall standard deviation is 74.33, variance is 5526.05;
#must use negative binomial if data is overdispersed;
#str(mydata);

#fit <- glmer(Total_Weighted_Health_Survey_Sco ~ Adjusted_RN_Staffing_Hours_per_R + (1|Provider_State/Subgroup), family="poisson", data=mydata);
#summary(fit);
#transform(mydata, obs=factor(seq(nrow(mydata))); #check for overdispersion https://stackoverflow.com/questions/19414336/using-glmer-for-nested-data
#update(fit, .~. + (1|obs));

#look at page 64 of Hadfield's Course Notes to specify random effects;
#us variance structure, of idh variance structure?;
citation();
citation("MCMCglmm");
citation("igraph");
packageVersion("MCMCglmm");
packageVersion("igraph");

#https://gkhajduk.github.io/2017-10-25-cleanMCMCglmm/
clean.MCMC <- function(x) {
    sols <- summary(x)$solutions  ## pull out relevant info from model summary
    Gcovs <- summary(x)$Gcovariances
    Rcovs <- summary(x)$Rcovariances
    fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
    random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
    residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
    names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
    names(random)[names(random) == "row.names.Gcovs."] <- "variable"
    names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
    fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
    random$effect <- "random"
    residual$effect <- "residual"
    modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # add the model name

#vignette("CourseNotes",package="MCMCglmm"); #look at pg 100 for zero inflated Poisson model;

#fit.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR + Subgroup_Simple_Overall_by_HRR), family="poisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#fit.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per, random=~us(1):(HRR + Subgroup_Simple_Overall_by_HRR), family="poisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR + Subgroup_Simple_Overall_by_HRR), family="poisson", data=mydata, saveX=TRUE, verbose=TRUE);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR:Subgroup_Simple_HRR), family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR), family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR + meanDegreeMultiple_HRR + proportionMultiple_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);

#-----------;
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q3/004433.html
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per, random=~HRR:Subgroup_Simple_Overall_by_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + Subgroup_Simple_HRR + HRR + Subgroup_Simple_HRR:HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE);


#--------------------------------#;
###### THIS MODEL #####
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR + meanDegreeMultiple_HRR + proportionMultiple_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=60000, nitt=1060000, thin=1000);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Years_Business + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + Subgroup_Simple_Overall_by_HRR + meanDegreeMultiple_HRR + proportionMultiple_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + proportionMultiple_HRR + meanDegreeMultiple_HRR + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);

#odel1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model1);
#mydataMDS <- na.omit(mydata);
#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model2);

#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR + Subgroup_Simple_HRR + HRR:Subgroup_Simple_HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model1);
#mydataMDS <- na.omit(mydata);
#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~us(proportionMultiple_HRR):(HRR), family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model2);

#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~us(proportionMultiple_HRR):(Provider_State) + us(proportionMultiple_HRR):(HRR), family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model1);
#mydataMDS <- na.omit(mydata);
#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~us(proportionMultiple_HRR):(Provider_State) + us(proportionMultiple_HRR):(HRR), family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model2);


#http://www.maths.bath.ac.uk/~jjf23/mixchange/nested.html

#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019896.html FIT RANDOM AS random=~Provider_State + HRR AS THESE EFFECTS ARE CROSS-CLASSIFIED IMPORTANT
#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model1);
#mydataMDS <- na.omit(mydata);
#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10);
#summary(model2);

#str(mydata);
#(1000000 - 150000)/85; #sample size of 10000;

#---------#;
#rcov=~proportionMultiple_HRR:units + meanDegreeMultiple_HRR:units
#sd(mydata$proportionMultiple_HRR); #15 is sd;
#sd(mydata$meanDegreeMultiple_HRR); #3 is sd;


                 
prior1 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*0^2),
         G2=list(V        = diag(3),
                 n        = 3,
                 alpha.mu = rep(0, 3),
                 alpha.V  = diag(3)*0^2)));
    
prior2 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(1)^2),
         G2=list(V        = diag(3),
                 n        = 3,
                 alpha.mu = rep(0, 3),
                 alpha.V  = diag(3)*0^2)));

data.V <- c(0,0,0,0,0,0,0,0,0);
matrix.V <- matrix(data.V, nrow = 3, ncol = 3, byrow = TRUE);
matrix.V;
prior3 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(3),
                 n        = 3,
                 alpha.mu = rep(0, 3),
                 alpha.V  = matrix.V)));
                 
prior4 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(3),
                 n        = 3,
                 alpha.mu = rep(0, 3),
                 alpha.V  = matrix.V),
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(1)^2)        
                 ));
                 
prior5 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2)      
                 ));
                 
prior6 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2), 
         G3=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = diag(2)*0^2)      
                 ));
                 
prior7 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)        
                 ));

#data.V8 <- c(25,0,0,0,100,0,0,0,100);
data.V8 <- c(25,0,0,100);
matrix.V8 <- matrix(data.V8, nrow = 2, ncol = 2, byrow = TRUE);
matrix.V8;                 
prior8 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(2),
                 n        = 2,
                 alpha.mu = rep(0, 2),
                 alpha.V  = matrix.V8)));
                 
prior9 <- list(
  R=list(V=1, nu = 0.02),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)      
                 ));
                 
priorLT1 <- list(
  R=list(V=1, nu = 0),
  G=list(G1=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G2=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2),
         G3=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(10)^2),      
         G4=list(V        = diag(1),
                 n        = 1,
                 alpha.mu = rep(0, 1),
                 alpha.V  = diag(1)*(0)^2)));

##how to specify fixed effect as random slope;            
##https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019896.html 

#cb_tmp <- mydata[which(mydata$Processing_Date == "2016-03-01"), ];
#cb_tmp <- cb_tmp[which(cb_tmp$HRR == 218), ];
#cb_tmp$meanDegreeMultiple_HRRSc;
#cb_tmp$Degree_Fac;

#length(mydata[which(mydata$Special_Focus_Facility == TRUE), ]$Federal_Provider_Number); #only 296 Special Focus Facilities over all cycles, not enough for interation;

###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=25000, nitt=100000, thin=10, pr=TRUE);
###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10, pr=TRUE);
###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10, pr=TRUE);
###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior3);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR + Provider_State:HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior4);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior5);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRR_Bi + meanDegreeMultiple_HRR_Bi, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior5);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior5);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR + idh(proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior6);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + HRR + proportionMultiple_HRR + meanDegreeMultiple_HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior7);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=100000, nitt=475000, thin=75, pr=TRUE, prior=prior2);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior2);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);

######################;
#### re-arrange reference category for ownership type to government, 24 October 2019;
#### this part doesn't work, as MCMCglmm ignores contrasts for a multivariate model 
#### https://stat.ethz.ch/pipermail/r-sig-mixed-models/2015q2/023699.html
#mydata1 <- mydata;
#levels(mydata1$Ownership_Type);
#mydata1$intvar <- interaction(mydata1$Ownership_Type, mydata1$Subgroup_Simple_Overall_by_HRR); 
#levels(mydata1$intvar);
##[1] "For-profit.Multiple" "Government.Multiple" "Non-profit.Multiple" "For-profit.Single"   "Government.Single"   "Non-profit.Single"  
#mydata1$intvar <- factor(mydata1$intvar, levels(mydata1$intvar)[c(2,5,1,3,4,6)]);
#
#
##https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/
##look at user defined coding in above link at bottom;
#mat <- matrix(c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1, -1, 0, 0, 0, 0), ncol=2);
#mat;
#my.contrasts <- mat[, 2];
#contrasts(mydata1$intvar) <- my.contrasts;
#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + intvar, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", rcov=~intvar:units, data=mydata1, saveX=TRUE, verbose=TRUE, burnin=500, nitt=1000, thin=10, pr=TRUE, prior=priorLT1, singular.ok=TRUE);
#summary(model1);
#str(model1);
######################;


#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=15000, nitt=115000, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type/Subgroup_Simple_Overall_by_HRR + Special_Focus_Facility/Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type/meanDegreeMultiple_HRRSc + Special_Focus_Facility/meanDegreeMultiple_HRRSc, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
#####model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Degree_Fac + Ownership_Type/Degree_Fac + Special_Focus_Facility/Degree_Fac, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type/meanDegreeMultiple_HRRSc + Special_Focus_Facility/meanDegreeMultiple_HRRSc, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type/Subgroup_Simple_Overall_by_HRR + Special_Focus_Facility/Subgroup_Simple_Overall_by_HRR + Special_Focus_Facility/meanDegreeMultiple_HRRSc, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type*Subgroup_Simple_Overall_by_HRR + Special_Focus_Facility*Provider_Changed_Ownership_in_La + Special_Focus_Facility*Subgroup_Simple_Overall_by_HRR + Special_Focus_Facility*meanDegreeMultiple_HRRSc, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
##model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + HHI_HRRSc + affiliation_HHI_HRRSc + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydata, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=priorLT1);
###model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);
#save(model1, file=paste("PREM_LONG_summary_Model_1_", keyDate, ".RData", sep=""), compress="xz");
##summary(model1);
#colnames(mydata);

mydataModel1Check <- mydata[, c("Total_Weighted_Health_Survey_Sco", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Number_of_Residents_in_Certified", "Occupancy_Rate", "Years_Business", "Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR", "proportionMultiple_HRRSc", "meanDegreeMultiple_HRRSc", "delta_HHI_HRRSc", "Federal_Provider_Number", "Provider_State", "HRR", "Processing_Date")]; 
ch1 <- mydataModel1Check[which(mydataModel1Check$Processing_Date == "2016-03-01"), ];
ch2 <- mydataModel1Check[which(mydataModel1Check$Processing_Date == "2016-09-01"), ];
ch3 <- mydataModel1Check[which(mydataModel1Check$Processing_Date == "2017-03-01"), ];
ch4 <- mydataModel1Check[which(mydataModel1Check$Processing_Date == "2017-09-01"), ];
ch5 <- mydataModel1Check[which(mydataModel1Check$Processing_Date == "2018-03-01"), ];
length(ch1$Processing_Date); #10,728;
length(ch2$Processing_Date); #10,875;
length(ch3$Processing_Date); #11,009;
length(ch4$Processing_Date); #10,736;
length(ch5$Processing_Date); #10,943;
ch1 <- na.omit(ch1);
ch2 <- na.omit(ch2);
ch3 <- na.omit(ch3);
ch4 <- na.omit(ch4);
ch5 <- na.omit(ch5);
length(ch1$Processing_Date); #10,728;
length(ch2$Processing_Date); #10,875;
length(ch3$Processing_Date); #11,009;
length(ch4$Processing_Date); #10,736;
length(ch5$Processing_Date); #10,943;
colnames(mydataModel1Check);
mydata1a <- ch1[,c(19,21,22)];
mydata2a <- ch2[,c(19,21,22)];
mydata3a <- ch3[,c(19,21,22)];
mydata4a <- ch4[,c(19,21,22)];
mydata5a <- ch5[,c(19,21,22)];
colnames(mydata1a) <- c("Federal_Provider_Number", "HRR", "Processing_Date1");
colnames(mydata2a) <- c("Federal_Provider_Number", "HRR", "Processing_Date2");
colnames(mydata3a) <- c("Federal_Provider_Number", "HRR", "Processing_Date3");
colnames(mydata4a) <- c("Federal_Provider_Number", "HRR", "Processing_Date4");
colnames(mydata5a) <- c("Federal_Provider_Number", "HRR", "Processing_Date5");
#mydataAlla <- join_all(list(mydata1a, mydata2a), by="Federal_Provider_Number", type="left"); 
#mydataAlla <- join_all(list(mydataAlla, mydata3a), by="Federal_Provider_Number", type="left"); 
#mydataAlla <- join_all(list(mydataAlla, mydata4a), by="Federal_Provider_Number", type="left"); 
#mydataAlla <- join_all(list(mydataAlla, mydata5a), by="Federal_Provider_Number", type="left"); 
#mydataAlla <- na.omit(mydataAlla);
#length(mydataAlla$Federal_Provider_Number); #9,001 NHs, same as below;
mydataAlla <- join_all(list(mydata1a, mydata2a, mydata3a, mydata4a, mydata5a), by="Federal_Provider_Number", type="inner", match="all"); 
###9,001 NHs with complete data in each period with organization owner matched to first baseline cohort of 10,728 from March 2016;
###length(unique(mydataAlla$HRR)); #294 HRRs represented;
##which(is.na(mydataAlla$Processing_Date5)); #they all matched!

save(mydataAlla, file=paste("NHs_Model_1_", keyDate, ".RData", sep=""), compress="xz");
#load(paste("NHs_Model_1_", keyDate, ".RData", sep=""));

#-------------------------#;



########## -------------- Start with baseline cohort and link ahead for location -------------------;
df <- nh.comb.countiesCTs[which(nh.comb.countiesCTs$PERIOD == "March2016"), ];
#length(df); 
#rm(nh.comb.countiesCTs); #unload memory intensive;
#which(is.na(df$Ownership_Type)); #none, good check;
#which(df$Ownership_Type == ""); #1749 are blank, not good, check;
#which(is.na(df$Processing_Date)); #none;

#df$lat;
#df$lon;
nhNA <- df[which(is.na(df$lat)),];
o2 <- as.data.frame(df)

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;
##### This does for NHs owned by organizations and those without missing values in statistical model #####;
o2 = o2[which(o2$Federal_Provider_Number %in% mydataAlla$Federal_Provider_Number),];
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;

#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS(projected);
########## -------------- Start with baseline cohort and link ahead -------------------;


mydataMDS <- na.omit(mydata);

### count nursing homes with complete MDS data;
mydataModel2Check <- mydataMDS;
cha1 <- mydataModel2Check[which(mydataModel2Check$Processing_Date == "2016-03-01"), ];
cha2 <- mydataModel2Check[which(mydataModel2Check$Processing_Date == "2016-09-01"), ];
cha3 <- mydataModel2Check[which(mydataModel2Check$Processing_Date == "2017-03-01"), ];
cha4 <- mydataModel2Check[which(mydataModel2Check$Processing_Date == "2017-09-01"), ];
cha5 <- mydataModel2Check[which(mydataModel2Check$Processing_Date == "2018-03-01"), ];
length(cha1$Processing_Date); #7,993;
length(cha2$Processing_Date); #8,996;
length(cha3$Processing_Date); #9,148;
length(cha4$Processing_Date); #9,005;
length(cha5$Processing_Date); #9,224;
cha1 <- na.omit(cha1);
cha2 <- na.omit(cha2);
cha3 <- na.omit(cha3);
cha4 <- na.omit(cha4);
cha5 <- na.omit(cha5);
length(cha1$Processing_Date); #10,728;
length(cha2$Processing_Date); #10,875;
length(cha3$Processing_Date); #11,009;
length(cha4$Processing_Date); #10,736;
length(cha5$Processing_Date); #10,943;
colnames(mydataModel2Check);
mydata1aa <- cha1[,c(1,78,14)];
mydata2aa <- cha2[,c(1,78,14)];
mydata3aa <- cha3[,c(1,78,14)];
mydata4aa <- cha4[,c(1,78,14)];
mydata5aa <- cha5[,c(1,78,14)];
colnames(mydata1aa) <- c("Federal_Provider_Number", "HRR", "Processing_Date1");
colnames(mydata2aa) <- c("Federal_Provider_Number", "HRR", "Processing_Date2");
colnames(mydata3aa) <- c("Federal_Provider_Number", "HRR", "Processing_Date3");
colnames(mydata4aa) <- c("Federal_Provider_Number", "HRR", "Processing_Date4");
colnames(mydata5aa) <- c("Federal_Provider_Number", "HRR", "Processing_Date5");
mydataAllMDSa <- join_all(list(mydata1aa, mydata2aa, mydata3aa, mydata4aa, mydata5aa), by="Federal_Provider_Number", type="inner", match="all"); 
#length(mydataAllMDSa$Federal_Provider_Number); #6,693 NHs, same as below;
###length(unique(mydataAllMDSa$HRR)); #293 HRRs represented;
###


###length(mydataMDS$Federal_Provider_Number);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=25000, nitt=100000, thin=10, pr=TRUE);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital+ Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_406 + mds_407 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10, pr=TRUE);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior1);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior5);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior1);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=100000, nitt=450000, thin=75, pr=TRUE, prior=prior2);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior2);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);

#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471 + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=15000, nitt=115000, thin=10, pr=TRUE, prior=priorLT1);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + delta_HHI_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471 + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1500, nitt=11500, thin=10, pr=TRUE, prior=priorLT1);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + HHI_HRRSc + affiliation_HHI_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471 + Ownership_Type*Subgroup_Simple_Overall_by_HRR, random=~Federal_Provider_Number + Provider_State + HRR + Processing_Date, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=priorLT1);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);
#save(model2, file=paste("PREM_LONG_summary_Model_2_", keyDate, ".RData", sep=""), compress="xz");
##summary(model2);

load(paste("HSRmodel_published/PREM_LONG_summary_Model_1_", keyDate, ".RData", sep=""));
load(paste("HSRmodel_published/PREM_LONG_summary_Model_2_", keyDate, ".RData", sep=""));
##load(paste("Chapter4ModelBackup_published/PREM_summary_Model_1_", keyDate, ".RData", sep=""));
##load(paste("Chapter4ModelBackup_published/PREM_summary_Model_2_", keyDate, ".RData", sep=""));
#summary(model1);
#summary(model2);

#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV);
#autocorr.diag(model2$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;
#---------#;



#https://stats.stackexchange.com/questions/307021/why-does-mcmcglmm-result-in-small-effective-sample-sizes-for-logistic-regression/310291
#autocorr.diag(model1$VCV);
#autocorr.diag(model2$VCV); #use this to choose thinning interval, autocorrelation less than 0.01 ideally;
#---------#;

fe1 <- summary(model1);
#str(fe1$solutions);
fe1$solutions[,1] <- exp(fe1$solutions[,1]);
fe1$solutions[,2] <- exp(fe1$solutions[,2]);
fe1$solutions[,3] <- exp(fe1$solutions[,3]);
fe2 <- summary(model2);
#str(fe1$solutions);
fe2$solutions[,1] <- exp(fe2$solutions[,1]);
fe2$solutions[,2] <- exp(fe2$solutions[,2]);
fe2$solutions[,3] <- exp(fe2$solutions[,3]);

writeLines(capture.output(fe1), paste("PREM_LONG_RRs_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(fe2), paste("PREM_LONG_RRs_Model_2_", keyDate, ".txt", sep=""));
writeLines(capture.output(summary(model1)), paste("PREM_LONG_summary_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(summary(model2)), paste("PREM_LONG_summary_Model_2_", keyDate, ".txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html
#colnames(model1$Sol); #1:24 are explanatory variables;
#colnames(model2$Sol); #1:41 are explanatory variables;

EVs1 <- model1$Sol[, 1:24];
EVs2 <- model2$Sol[, 1:41];
intPSs1 <- model1$Sol[, 12008:12058];
intPSs2 <- model2$Sol[, 10371:10421];
intHRRs1 <- model1$Sol[, 12059:12358];
intHRRs2 <- model2$Sol[, 10422:10721];
#slopeProportionMultHRRs1 <- model1$Sol[, 368:662];
#slopeProportionMultHRRs2 <- model2$Sol[, 384:677];
#slopeDegreeMultHRRs1 <- model1$Sol[, 663:957];
#slopeDegreeMultHRRs2 <- model2$Sol[, 678:971];
##plot(slopePropotionMultHRRs1);

length(posterior.mode(intHRRs1)); #300;  #mean intercept;
#length(posterior.mode(slopeProportionMultHRRs1)); #295; #these slopes are not much different from 0;
#length(posterior.mode(slopeDegreeMultHRRs1)); #295;  #these slopes are not much different from 0;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_LONG_VCV_", keyDate, ".pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_LONG_sol_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_LONG_provider_state_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPSs1)),2), mar=c(2,2,1,0));
plot(intPSs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_LONG_hrr_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=200, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intHRRs1)),2), mar=c(2,2,1,0));
plot(intHRRs1, auto.layout=F);
dev.off()

#---;

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_LONG_VCV_", keyDate, ".pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model2$VCV)),2), mar=c(2,2,1,0));
plot(model2$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_LONG_sol_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs2)),2), mar=c(2,2,1,0));
plot(EVs2, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_LONG_provider_state_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPSs2)),2), mar=c(2,2,1,0));
plot(intPSs2, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_LONG_hrr_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=200, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intHRRs2)),2), mar=c(2,2,1,0));
plot(intHRRs2, auto.layout=F);
dev.off()

#plot(model2);
#############################################################;

##mpred1 <- predict(model1, interval="confidence"); #This gives predicted value for each NH;
##mpredHRRs1 <-predict(model1, marginal=model1$Random$formula, type="response", interval="confidence", level=0.95, posterior="all", verbose=FALSE, approx="numerical");
#mpred1 <- predict(model1); #This gives predicted value for each NH;
#
#datapred1 <- unique(data.frame(cbind(fit=mpred1, hrr=mydata$HRR, proportionMultiple_HRR=mydata$proportionMultiple_HRR)));
#colnames(datapred1);
#length(datapred1$V1);
#
#
##datapred1  <- unique(data.frame(cbind(predlogit = predlogit, comm = mydata$comm, mage =  mydata$mage)));
#x#yplot(V1 ~ proportionMultiple_HRR, data=datapred1, groups=hrr, type= c("p", "l", "g"), col="blue", xlim=c(9, 51), ylim=c(-4, 4));
#xyplot(V1 ~ proportionMultiple_HRR, data=datapred1, groups=hrr, type= c("p", "l", "g"), col="blue");
#
#pdf.options(encoding='CP1250');
#cairo_pdf(file=paste("model1_random_slope_", keyDate, ".pdf", sep=''), bg="transparent", width=200, height=11, pointsize=6, family="FreeSans");
#xyplot(V1 ~ proportionMultiple_HRR, data=datapred1, groups=hrr, type= c("p", "l", "g"), col="blue");
#dev.off()


## BELOW LINK IS IMPORTANT!!!
##https://tomhouslay.files.wordpress.com/2017/02/indivvar_plasticity_tutorial_mcmcglmm1.pdf



#https://github.com/tmalsburg/MCMCglmm-intro
## need to calculate ICC for mcmcglmm.....;

#model1G <- lapply(model1, function(m) m$Sol)
#model1G <- do.call(mcmc.list, model1$Sol)
#
#gelman.diag(model1G);
#gelman.diag(model2G);

#---------------------------------
#str(model2);
#colnames(model1$Sol); #select the random effects manually from here;
#colnames(model2$Sol); #select the random effects manually from here;
#colnames(model1$Sol[, 73:367]); #(Intercept).HRR...;
#colnames(model2$Sol[, 90:383]); #(Intercept).HRR...;

#HPDinterval(intHRRs1 + slopeProportionMultHRRs1 + slopeDegreeMultHRRs1, prob=0.95);

ranefHRR1 <- cbind(B = posterior.mode(intHRRs1), CI = HPDinterval(intHRRs1, prob=0.95));
ranefHRR2 <- cbind(B = posterior.mode(intHRRs2), CI = HPDinterval(intHRRs2, prob=0.95));
ranefHRR1 <- as.data.frame(ranefHRR1);
ranefHRR2 <- as.data.frame(ranefHRR2);
#colnames(ranefHRR1);
ranefHRR1$B <- exp(ranefHRR1$B);
ranefHRR1$lower <- exp(ranefHRR1$lower);
ranefHRR1$upper <- exp(ranefHRR1$upper);
ranefHRR2$B <- exp(ranefHRR2$B);
ranefHRR2$lower <- exp(ranefHRR2$lower);
ranefHRR2$upper <- exp(ranefHRR2$upper);


qfun <- function(x, lev) unname(quantile(x, lev));

rsumPS1 <- as.data.frame(t(apply(t(intPSs1), 1, function(x) c(est=mean(x), min=qfun(x, 0.025),max=qfun(x, 0.975)))));
rownames(rsumPS1) <- sub('Provider_State.', '', rownames(t(intPSs1)));
#colnames(rsumPS1);
rsumPS1$est <- exp(rsumPS1$est);
rsumPS1$min <- exp(rsumPS1$min);
rsumPS1$max <- exp(rsumPS1$max);
rsumPS1$term <- reorder(factor(rownames(rsumPS1)), rsumPS1$est);

rsumPS2 <- as.data.frame(t(apply(t(intPSs2), 1, function(x) c(est=mean(x), min=qfun(x, 0.025),max=qfun(x, 0.975)))));
rownames(rsumPS2) <- sub('Provider_State.', '', rownames(t(intPSs2)));
rsumPS2$est <- exp(rsumPS2$est);
rsumPS2$min <- exp(rsumPS2$min);
rsumPS2$max <- exp(rsumPS2$max);
rsumPS2$term <- reorder(factor(rownames(rsumPS2)), rsumPS2$est);

rownames(ranefHRR1) <- gsub("[^0-9]", "", rownames(ranefHRR1)) 
ranefHRR1$hrrnum <- rownames(ranefHRR1);
ranefHRR1$term <- reorder(factor(rownames(ranefHRR1)), ranefHRR1$B);

rownames(ranefHRR2) <- gsub("[^0-9]", "", rownames(ranefHRR2)) 
ranefHRR2$hrrnum <- rownames(ranefHRR2);
ranefHRR2$term <- reorder(factor(rownames(ranefHRR2)), ranefHRR2$B);


pREStates1 <- 
	ggplot(rsumPS1,aes(term,est)) +
	geom_hline(yintercept=1, color="grey") +
	theme_bw() +
  	theme(panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
    geom_pointrange(aes(ymin=min,ymax=max)) +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_flip() + 
    ggtitle("Caterpillar Plot of Prevalence Ratio for\nAmerican State with 95% Highest Posterior\nDensity Interval from Model 1") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) + 
 	labs(y="Prevalence Ratio", x="American State") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model1_LONG_Provider_State_", keyDate, ".pdf", sep=""), pREStates1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

pREStates2 <- 
	ggplot(rsumPS2,aes(term,est)) +
	geom_hline(yintercept=1, color="grey") +
	theme_bw() +
  	theme(panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
    geom_pointrange(aes(ymin=min,ymax=max)) +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_flip() + 
    ggtitle("Caterpillar Plot of Prevalence Ratio for\nAmerican State with 95% Highest Posterior\nDensity Interval from Model 2") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) + 
 	labs(y="Prevalence Ratio", x="American State") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model2_LONG_Provider_State_", keyDate, ".pdf", sep=""), pREStates2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

pREHRRs1 <- 
	ggplot(ranefHRR1, aes(term,B)) +
	geom_hline(yintercept=1, color="grey") +
	theme_bw() +
  	theme(panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
    geom_pointrange(aes(ymin=lower,ymax=upper)) +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_flip() + 
    ggtitle("Caterpillar Plot of Prevalence Ratio\nfor HRR with 95% Highest Posterior\nDensity Interval from Model 1") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.text.y=element_text(family="FreeSans", size=3)) + 
 	labs(y="Prevalence Ratio", x="HRR") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model1_LONG_HRR_", keyDate, ".pdf", sep=""), pREHRRs1, width = 10, height = 14, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

pREHRRs2 <- 
	ggplot(ranefHRR2, aes(term,B)) +
	geom_hline(yintercept=1, color="grey") +
	theme_bw() +
  	theme(panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
    geom_pointrange(aes(ymin=lower,ymax=upper)) +
    scale_y_continuous(breaks=seq(0,10,0.25)) +
    coord_flip() + 
    ggtitle("Caterpillar Plot of Prevalence Ratio\nfor HRR with 95% Highest Posterior\nDensity Interval from Model 2") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.text.y=element_text(family="FreeSans", size=3)) + 
 	labs(y="Prevalence Ratio", x="HRR") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model2_LONG_HRR_", keyDate, ".pdf", sep=""), pREHRRs2, width = 10, height = 14, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

#---------------------------------

colnames(model1$VCV);
ICC_Provider_State.1 <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
ICC_HRR.1 <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
#ICC_proportionMultiple_HRR.1 <- model1$VCV[, 3]/(rowSums(model1$VCV)); 
#ICC_meanDegreeMultiple_HRR.1 <- model1$VCV[, 4]/(rowSums(model1$VCV)); 
posterior.mode(ICC_Provider_State.1);
posterior.mode(ICC_HRR.1);
#posterior.mode(ICC_proportionMultiple_HRR.1);
#posterior.mode(ICC_meanDegreeMultiple_HRR.1);

dft1.1 <- cbind(ICC = posterior.mode(ICC_Provider_State.1), CI = HPDinterval(ICC_Provider_State.1));
dft2.1 <- cbind(ICC = posterior.mode(ICC_HRR.1), CI = HPDinterval(ICC_HRR.1));
#dft3.1 <- cbind(ICC = posterior.mode(ICC_proportionMultiple_HRR.1), CI = HPDinterval(ICC_proportionMultiple_HRR.1));
#dft4.1 <- cbind(ICC = posterior.mode(ICC_meanDegreeMultiple_HRR.1), CI = HPDinterval(ICC_meanDegreeMultiple_HRR.1));
#dft1 <- rbind(dft1.1, dft2.1, dft3.1, dft4.1);
#row.names(dft1) <- c("Provider_State", "(Intercept).HRR", "proportionMultiple_HRR.HRR", "meanDegreeMultiple_HRR.HRR");
dft1 <- rbind(dft1.1, dft2.1);
row.names(dft1) <- c("Provider_State", "(Intercept).HRR");

colnames(model2$VCV);
ICC_Provider_State.2 <- model2$VCV[, 2]/(rowSums(model2$VCV)); 
ICC_HRR.2 <- model2$VCV[, 3]/(rowSums(model2$VCV)); 
#ICC_proportionMultiple_HRR.2 <- model2$VCV[, 3]/(rowSums(model2$VCV)); 
#ICC_meanDegreeMultiple_HRR.2 <- model2$VCV[, 4]/(rowSums(model2$VCV)); 
posterior.mode(ICC_Provider_State.2);
posterior.mode(ICC_HRR.2);
#posterior.mode(ICC_proportionMultiple_HRR.2);
#posterior.mode(ICC_meanDegreeMultiple_HRR.2);

dft1.2 <- cbind(ICC = posterior.mode(ICC_Provider_State.2), CI = HPDinterval(ICC_Provider_State.2));
dft2.2 <- cbind(ICC = posterior.mode(ICC_HRR.2), CI = HPDinterval(ICC_HRR.2));
#dft3.2 <- cbind(ICC = posterior.mode(ICC_proportionMultiple_HRR.2), CI = HPDinterval(ICC_proportionMultiple_HRR.2));
#dft4.2 <- cbind(ICC = posterior.mode(ICC_meanDegreeMultiple_HRR.2), CI = HPDinterval(ICC_meanDegreeMultiple_HRR.2));
#dft2 <- rbind(dft1.2, dft2.2, dft3.2, dft4.2);
#row.names(dft2) <- c("Provider_State", "(Intercept).HRR", "proportionMultiple_HRR.HRR", "meanDegreeMultiple_HRR.HRR");
dft2 <- rbind(dft1.2, dft2.2);
row.names(dft2) <- c("Provider_State", "(Intercept).HRR");

writeLines(capture.output(dft1), paste("PREM_LONG_ICCs_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(dft2), paste("PREM_LONG_ICCs_Model_2_", keyDate, ".txt", sep=""));


##https://jrosen48.github.io/blog/explorations-in-markov-chain-monte-carlo-mcmc/
#above tells how to interpret post.mean, it is variance for fixed effects;

#int.slope.cor <- model2$VCV[, 2]/sqrt(model2$VCV[,1] * model2$VCV[, 4]);
#posterior.mode(int.slope.cor);
#
#
#ref <- data.frame(model2$VCV);
#rdf <- melt(ref, value.name="Total_Weighted_Health_Survey_Sco");
#ggplot(rdf, aes(x=sqrt(Total_Weighted_Health_Survey_Sco),color=variable)) + geom_density()
#
#ggplot(ref, aes(x=sqrt(Provider_State),y=sqrt(HRR)))+geom_density2d()+geom_abline(int=0,slope=1)
#

oneModel <- clean.MCMC(model1);  # get all the info from summary(modelName)
oneModel$modelName <- getName.MCMC(model1);  # add the model's name in a new column
oneModel;  # check out the created dataframe

dataList <- list(model1, model2);
dataListNames <- list("model1", "model2");
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F);
mcmcOutputs <- as.data.frame(do.call(rbind, readyList), stringsAsFactors = FALSE)
stargazer(mcmcOutputs, type="text", summary=FALSE);
writeLines(capture.output(stargazer(mcmcOutputs, type="text", summary=FALSE)), paste("PREM_LONG_stargazer_Models_1_and2_", keyDate, ".txt", sep=""));

#######################
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2011q1/012403.html

#-~-~-~-~-~-~-~-~-#;
#https://stackoverflow.com/questions/47598123/how-do-i-extract-random-effects-from-mcmcglmm
#if (!require("postMCMCglmm")) {
#    devtools::install_github("JWiley/postMCMCglmm", force=TRUE)
#    library("postMCMCglmm")
#}
#-~-~-~-~-~-~-~-~-#;

#--------------------------------#;


#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR), family="poisson", data=mydata, saveX=TRUE, verbose=TRUE);
#summary(fit1.MCMC);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per, random=~us(1):(HRR + Subgroup_Simple_Overall_by_HRR), family="poisson", data=mydata, verbose=TRUE);
#summary(fit1.MCMC);
#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_CNA_Staffing_Hours_per + meanDegreeMultiple_HRR + proportionMultiple_HRR + Subgroup_Simple_Overall_by_HRR, random=~Subgroup_Simple_HRR, family="poisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#summary(fit1.MCMC);

#plot(mcmc.list(fit1.MCMC$VCV[, 1], fit1.MCMC$VCV[, 1]));
#
#nz <- 1:1000;
#oz <- sum(mydata$Total_Weighted_Health_Survey_Sco == 0);
#for (i in 1:1000) {
#	pred.l <- rnorm(915, (fit1.MCMC$X %*% fit1.MCMC$Sol[i])@x, sqrt(fit1.MCMC$VCV[i]))
#	nz[i] <- sum(rpois(915, exp(pred.l)) == 0)
#}
#-----------;

#fit1.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per + Total_Adjusted_Staffing_HRD, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR + Subgroup_Simple_Overall_by_HRR), family="zipoisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#summary(fit1.MCMC);
#plot(fit1.MCMC);

#offvec <- c(1,1,2,rep(1,5));
#prior_overdisp  <- list(R=list(V=diag(c(1,1)),nu=0.002,fix=2), G=list(list(V=diag(c(1,1e-6)),nu=0.002,fix=2)));
#prior_overdisp_broodoff <- c(prior_overdisp, list(B=list(mu=c(0,1)[offvec], V=diag(c(1e8,1e-6)[offvec]))));
#fit2.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Number_of_Certified_Beds_Cat + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per + Total_Adjusted_Staffing_HRD, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR):(HRR + Subgroup_Simple_Overall_by_HRR), family="zipoisson", data=mydata, prior=prior_overdisp_broodoff, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#summary(fit2.MCMC);

#fit2.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per, random=~us(1 + meanDegreeMultiple_HRR + proportionMultiple_HRR + Subgroup_Simple_Overall_by_HRR):(HRR), family="poisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#summary(fit2.MCMC);

#fit3.MCMC <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Occupancy_Rate + Adjusted_RN_Staffing_Hours_per_R + Adjusted_LPN_Staffing_Hours_per + Adjusted_CNA_Staffing_Hours_per, random=~HRR:Federal_Provider_Number, family="poisson", data=mydata, verbose=TRUE, nitt=6000, burnin=1000, thin=10);
#summary(fit3.MCMC);

#fit <- glmer(Total_Weighted_Health_Survey_Sco ~ Adjusted_RN_Staffing_Hours_per_R + (1|HRR/Subgroup_Simple_HRR), family=negative.binomial(2), data=mydata2);
#summary(fit);
#oddsratios.lmer <- exp(,,1)^fixef(fit);

#-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-$-#;




#--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--#;
#--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--#;


#--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--#;
#--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--#;




#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#


################################ SHAPEFILE MAPPING ################################;
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

###########################;
remove.territories = function(.df) {
    subset(.df, 
        .df$id != "AS" &
        .df$id != "MP" &
        .df$id != "GU" & 
        .df$id != "PR" &
        .df$id != "VI" 
    )
}
remove.fake = function(.df) {
    subset(.df, 
        .df$Federal_Provider_Number != "fake"
    )
}
#save("remove.territories", file = "helpers/remove.territories")

plain_theme = theme(axis.text=element_blank()) + 
    theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.ticks = element_blank())

no_ylab = ylab("") 
no_xlab = xlab("")


us_aea = spTransform(usStates, CRS(projected))
#colnames(us_aea@data);
us_aea = us_aea[, which(colnames(us_aea@data) %in% "STUSPS")]
us_aea@data$id = us_aea@data$STUSPS

alaska = us_aea[us_aea$STUSPS=="AK",]
	#bbox_a1 <- bbox(alaska)
	#center_a1 <- coordinates(alaska)
alaska = elide(alaska, rotate=-50)
	#bbox_a2 <- bbox(alaska)
	#center_a2 <- coordinates(alaska)
alaska = elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska = elide(alaska, shift=c(-2100000, -2500000))
proj4string(alaska) = CRS(projected)

hawaii = us_aea[us_aea$STUSPS=="HI",]
hawaii = elide(hawaii, rotate=-35)
hawaii = elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) = CRS(projected)

us_aea = us_aea[!us_aea$STUSPS %in% c("AK", "HI"),]
us_aea = rbind(us_aea, alaska, hawaii)

us50 <- fortify(us_aea, region="STUSPS")
us50 = remove.territories(us50)
#save('helpers/us50')

#------------------------------------------------------------------;

us_hrr = spTransform(hrrBoundary, CRS(projected))
us_hrr@data$STUSPS = us_hrr@data$hrr_num
us_hrr = us_hrr[, which(colnames(us_hrr@data) %in% "STUSPS")]
us_hrr@data$id = us_hrr@data$STUSPS

alaska_hrr = alaska
hawaii_hrr = hawaii

us_hrr = us_hrr[!us_hrr$STUSPS %in% c("10", "150"),]
us_hrr = rbind(us_hrr, alaska_hrr, hawaii_hrr) #add alaska hawaii from other map;
#plot(us_hrr);

us_hrr@data[which(us_hrr@data$STUSPS %in% "AK"), ] <- 10; #hrrnum is 10; 
us_hrr@data[which(us_hrr@data$STUSPS %in% "HI"), ] <- 150; #hrrnum is 150;
us_hrr@data[which(us_hrr@data$id %in% "AK"), ] <- 10; #hrrnum is 10; 
us_hrr@data[which(us_hrr@data$id %in% "HI"), ] <- 150; #hrrnum is 150;
row.names(us_hrr@data)[which(us_hrr@data$STUSPS %in% "10")] <- "AK - Anchorage";
row.names(us_hrr@data)[which(us_hrr@data$STUSPS %in% "150")] <- "HI - Honolulu";
colnames(us_hrr@data) <- c("hrrnum", "id");

#plot(us_hrr);
#writePolyShape(us_hrr, paste("HRRs_projected.shp", sep=""));

######
mydata$zipcode <- sprintf("%05s", mydata$Provider_Zip_Code); # fix to 5 characters;  
#str(lookup$zipcode);
length(mydata$zipcode); #10793;
mydata <- join_all(list(mydata, lookup), by = 'zipcode', type = "left"); 
length(mydata$zipcode); #10793;

tryCatch({
mydata <- edgelist[-c(which(is.na(mydata$hrrnum))), ];
}, error=function(e){})

length(mydata$zipcode); #10793;

dt <- data.table(mydata);
summaryHRR <- dt[ , list(mean=mean(Total_Weighted_Health_Survey_Sco), sd=sd(Total_Weighted_Health_Survey_Sco)), by=hrrnum];

### Rounds mean and SD to integer for nice plotting in maps;
summaryHRR$mean <- round(summaryHRR$mean, digits=0);
summaryHRR$sd <- round(summaryHRR$sd, digits=0);


#-@-@-@-@-@-
#---
#eval((ranefHRR1$lower < 1) & (ranefHRR1$upper < 1));
#eval((ranefHRR1$lower > 1) & (ranefHRR1$upper > 1));
ranefHRR1$sig <- 0;
ranefHRR1[which((ranefHRR1$lower < 1) & (ranefHRR1$upper < 1)), ]$sig <- 1;
ranefHRR1[which((ranefHRR1$lower > 1) & (ranefHRR1$upper > 1)), ]$sig <- 1;
ranefHRR1 <- ranefHRR1[ , -c(5)];
colnames(ranefHRR1) <- c("model1_RR", "model1_lower", "model1_upper", "hrrnum", "model1_sig");

ranefHRR2$sig <- 0;
ranefHRR2[which((ranefHRR2$lower < 1) & (ranefHRR2$upper < 1)), ]$sig <- 1;
ranefHRR2[which((ranefHRR2$lower > 1) & (ranefHRR2$upper > 1)), ]$sig <- 1;
ranefHRR2 <- ranefHRR2[ , -c(5)];
colnames(ranefHRR2) <- c("model2_RR", "model2_lower", "model2_upper", "hrrnum", "model2_sig");
#---

us_hrr@data <- join_all(list(us_hrr@data, summaryHRR), by = 'hrrnum', type = "left");
summaryHRR <- join_all(list(summaryHRR, ranefHRR1), by = 'hrrnum', type = "left");
summaryHRR <- join_all(list(summaryHRR, ranefHRR2), by = 'hrrnum', type = "left");
#writePolyShape(us_hrr, paste("HRRs_projected.shp", sep=""));

#------------------#;
###check, must address this in STEP_05 file prior and check here after ###
#us_hrr@data[which(is.na(us_hrr@data$mean)), ]; #305 hrrnum;
#lookup[which(lookup$hrrnum == 305), ];
#
#sort(unique(lookup$hrrnum)); #no 305 in lookup table, but in boundary file;
#length(unique(lookup$hrrnum)); #there are 306 unique hrrnums in lookup table;
#length(unique(us_hrr@data$hrrnum)); #there are 305 unique hrrnums in boundary file;
#
#'%!in%' <- function(x,y)!('%in%'(x,y))
#which(unique(us_hrr@data$hrrnum) %!in% unique(lookup$hrrnum)); #257;
#unique(us_hrr@data$hrrnum)[257]; # HRR 305;
#
#which(unique(lookup$hrrnum) %!in% unique(us_hrr@data$hrrnum)); #35, 248;
#unique(lookup$hrrnum)[35]; #HRR 351;
#unique(lookup$hrrnum)[248]; #HRR 417;
#lookup[which(lookup$hrrnum == 351), ]; #Erie, PA this one is miscoded as HRR 305 in boundary file;
#lookup[which(lookup$hrrnum == 417), ]; #Victoria, TX. Not in boundary file;
#lookup[which(lookup$hrrcity %in% "Corpus Christi"), ]; #hrrnum is 390;
#
#mydata[which(mydata$hrrnum == 351), ]; #Yes, in there;
#mydata[which(mydata$hrrnum == 417), ]; #Must change to hrrnum = 390 for Corpus Christi;
#------------------#;

us50_hrr <- fortify(us_hrr, region="hrrnum");
us50_hrr = remove.territories(us50_hrr);
summaryHRR$id <- summaryHRR$hrrnum;
us50_hrr <- join_all(list(us50_hrr, summaryHRR), by = 'id', type = "left");
colnames(us50_hrr);
#us50_hrr$mean;
#-@-@-@-@-@-


#-@-@-@-@-@-
#---
rsumPS1$id <- rsumPS1$term;
rsumPS1$sig <- 0;
rsumPS1[which((rsumPS1$min < 1) & (rsumPS1$max < 1)), ]$sig <- 1;
rsumPS1[which((rsumPS1$min > 1) & (rsumPS1$max > 1)), ]$sig <- 1;
rsumPS1 <- rsumPS1[ , -c(4)];
colnames(rsumPS1) <- c("model1_RR", "model1_lower", "model1_upper", "id", "model1_sig");

rsumPS2$id <- rsumPS2$term;
rsumPS2$sig <- 0;
rsumPS2[which((rsumPS2$min < 1) & (rsumPS2$max < 1)), ]$sig <- 1;
rsumPS2[which((rsumPS2$min > 1) & (rsumPS2$max > 1)), ]$sig <- 1;
rsumPS2 <- rsumPS2[ , -c(4)];
colnames(rsumPS2) <- c("model2_RR", "model2_lower", "model2_upper", "id", "model2_sig");
#---

us50 <- fortify(us50, region="id");
us50 = remove.territories(us50);
summaryPS <- join_all(list(rsumPS1, rsumPS2), by = 'id', type = "left");
us50 <- join_all(list(us50, summaryPS), by = 'id', type = "left");
#writePolyShape(us50, paste("PSs_projected.shp", sep=""));
colnames(us50);
#us50$mean;
#us50$group;
#-@-@-@-@-@-

#-------------------------------------------------------------;


p = ggplot(data=us50) + 
    geom_map(map=us50, aes(map_id=id, group=group), ,fill="white", color="dark grey", size=0.15) + 
    no_ylab + 
    no_xlab + 
    plain_theme
    
p_hrr = ggplot(data=us50_hrr) + 
    geom_map(map=us50_hrr, aes(map_id=id, group=group), ,fill="white", color="dark grey", size=0.15) + 
    no_ylab + 
    no_xlab + 
    plain_theme
#p
#us50$id 
###################################
#Below does Alaska Hawaii inset for Nursing home points like for US states above;
#alaskaCorner1 <- c(52.222162, -174.207006);
#alaskaCorner2 <- c(69.642491, -141.189882);
#alaskaCorner3 <- c(54.874814, -130.834586);
options(digits=9);
lon <- c(-174.207006, -141.189882, -130.834586, 172.471398, -167.684893, -166.058917, -156.676022,
			-160.184347, -154.854820, -155.680197, -157.998298, -160.542870);
lat <- c(52.222162, 69.642491, 54.874814, 52.922414, 65.604794, 68.812759, 71.279249,
			21.892327, 19.500730, 18.956699, 21.697147, 21.653282);
lon <- as.numeric(lon);
lat <- as.numeric(lat);
id <- c(rep("AK", 7), rep("HI", 5));
Federal_Provider_Number <- c(rep("fake", 12));
o2af <- cbind(lon, lat, id, Federal_Provider_Number);
o2af <- as.data.frame(o2af);
o2af$lon <- as.numeric(as.character(o2af$lon));
o2af$lat <- as.numeric(as.character(o2af$lat));
coordinates(o2af) <- c(which(colnames(o2af) %in% "lon"), which(colnames(o2af) %in% "lat"))
proj4string(o2af) <- CRS("+proj=longlat");
o2aft = spTransform(o2af, CRS(projected));
o2aft@data$lon <- coordinates(o2aft)[,1];
o2aft@data$lat <- coordinates(o2aft)[,2];

o2@data$id <- o2@data$Provider_State;
us_aeaNH  = spTransform(o2, CRS(projected));

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;
##### This does for NHs owned by organizations and those without missing values in statistical model #####;
##us_aeaNH = us_aeaNH[which(us_aeaNH@data$Federal_Provider_Number %in% mydata$Federal_Provider_Number),];

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;

alaskaNH = us_aeaNH[us_aeaNH$id=="AK",]
alaskaNH@data$lon <- coordinates(alaskaNH)[,1];
alaskaNH@data$lat <- coordinates(alaskaNH)[,2];
alaskaNHm1 <- join_all(list(alaskaNH@data, o2aft@data[which(o2aft@data$id == "AK"), ]), by = 'Federal_Provider_Number', type = "full", match = "first");
alaskaNHm1$Federal_Provider_Number;
alaskaNHm1$lon;
alaskaNHm1$lat;
coordinates(alaskaNHm1) <- c(which(colnames(alaskaNHm1) %in% "lon"), which(colnames(alaskaNHm1) %in% "lat"));
proj4string(alaskaNHm1) <- CRS(projected);
alaskaNH <- alaskaNHm1;

hawaiiNH = us_aeaNH[us_aeaNH$id=="HI",]
hawaiiNH@data$lon <- coordinates(hawaiiNH)[,1];
hawaiiNH@data$lat <- coordinates(hawaiiNH)[,2];
hawaiiNHm1 <- join_all(list(hawaiiNH@data, o2aft@data[which(o2aft@data$id == "HI"), ]), by = 'Federal_Provider_Number', type = "full", match = "first");
hawaiiNHm1$Federal_Provider_Number;
hawaiiNHm1$lon;
hawaiiNHm1$lat;
coordinates(hawaiiNHm1) <- c(which(colnames(hawaiiNHm1) %in% "lon"), which(colnames(hawaiiNHm1) %in% "lat"));
proj4string(hawaiiNHm1) <- CRS(projected);
hawaiiNH <- hawaiiNHm1;

alaskaNH = elide(alaskaNH, rotate=-50);
alaskaNH = elide(alaskaNH, scale=max(apply(bbox(alaskaNH), 1, diff))/2.3); 
alaskaNH = elide(alaskaNH, shift=c(-2100000, -2500000));
proj4string(alaskaNH) = CRS(projected)

hawaiiNH = elide(hawaiiNH, rotate=-35)
hawaiiNH = elide(hawaiiNH, shift=c(5400000, -1400000))
proj4string(hawaiiNH) = CRS(projected)

us_aeaNH = us_aeaNH[!us_aeaNH$id %in% c("AK", "HI"),]
us_aeaNH = rbind(us_aeaNH, alaskaNH, hawaiiNH)
us_aeaNH = remove.territories(us_aeaNH)
us_aeaNH = remove.fake(us_aeaNH)
#which(is.na(us_aeaNH@data$Ownership_Type)); #none, good check;
#which(us_aeaNH@data$Ownership_Type == ""); #1749 are blank, not good, check;
###################################
#plot(us_aea);
#plot(us_aeaNH, add=TRUE);


#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,11), oma=c(1,1,4,1), xpd=TRUE);
#plot(us_aeaNH);
#identify(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], lab=us_aeaNH$Federal_Provider_Number, tolerance = 0.05);

#identify(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], lab=us_aeaNH$Provider_State, tolerance = 0.05);
## Nursing homes with Federal_Provider_Number in (475058, 375199, 245451, 165595) need to be re-geocoded;
#str(df@data$Federal_Provider_Number)
which(df@data$Federal_Provider_Number %in% c("475058", "375199", "245451", "165595"));
regeocode <- df[which(df@data$Federal_Provider_Number %in% c("475058", "375199", "245451", "165595")), ];
#regeocode@data; #note the NAs, have to revise STEP_02_nursing_homes_US_Census_ACS.R script to remerge properly...;
#4680 	AKRON  CARE CENTER, INC         991 HIGHWAY AKRON, IA 51001   	"475058"   42.817580, -96.551201        
#7047	FAIRWAY VIEW NEIGHBORHOODS     	201 MARK DRIVE ORTONVILLE, MN 56278    "375199" 45.320215, -96.445103         
#10926  HERITAGE VILLAGE NURSING HOME  	801 HWY HOLDENVILLE, OK 74848    "245451"   35.100607, -96.405755      
#13139  MENIG NURSING HOME215 			TOM WICKER LANE RANDOLPH CENTER, VT 05061 	"165595" 43.943257, -72.607451
#df[which(df$Federal_Provider_Number %in% "475058"), ]$lat <- 42.817580;
#df[which(df$Federal_Provider_Number %in% "475058"), ]$lon <- -96.551201;
#df[which(df$Federal_Provider_Number %in% "375199"), ]$lat <- 45.320215;
#df[which(df$Federal_Provider_Number %in% "375199"), ]$lon <- -96.445103;
#df[which(df$Federal_Provider_Number %in% "245451"), ]$lat <- 35.100607;
#df[which(df$Federal_Provider_Number %in% "245451"), ]$lon <- -96.405755;
#df[which(df$Federal_Provider_Number %in% "165595"), ]$lat <- 43.943257;
#df[which(df$Federal_Provider_Number %in% "165595"), ]$lon <- -72.607451;
#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#-*-#

##----------##;
us_aeaNH <- us_aeaNH[order(us_aeaNH@data$Ownership_Type), ]; #must do this to have in right order for plots;
unique(us_aeaNH@data$Ownership_Type);
##----------##;

#--------#;
mydataAlla <- join_all(list(mydata1a, mydata2a, mydata3a, mydata4a, mydata5a), by="Federal_Provider_Number", type="inner", match="all"); 
mydataSmall <- mydata[, c("Federal_Provider_Number", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Processing_Date")];
mydataSmall_ch1 <- mydataSmall[which(mydataSmall$Processing_Date == "2016-03-01"), ];
mydataSmall_ch2 <- mydataSmall[which(mydataSmall$Processing_Date == "2016-09-01"), ];
mydataSmall_ch3 <- mydataSmall[which(mydataSmall$Processing_Date == "2017-03-01"), ];
mydataSmall_ch4 <- mydataSmall[which(mydataSmall$Processing_Date == "2017-09-01"), ];
mydataSmall_ch5 <- mydataSmall[which(mydataSmall$Processing_Date == "2018-03-01"), ];

mydataSmallAll_ch1 <- join_all(list(mydataAlla, mydataSmall_ch1), by="Federal_Provider_Number", type="left", match="all"); 
mydataSmallAll_ch2 <- join_all(list(mydataAlla, mydataSmall_ch2), by="Federal_Provider_Number", type="left", match="all"); 
mydataSmallAll_ch3 <- join_all(list(mydataAlla, mydataSmall_ch3), by="Federal_Provider_Number", type="left", match="all"); 
mydataSmallAll_ch4 <- join_all(list(mydataAlla, mydataSmall_ch4), by="Federal_Provider_Number", type="left", match="all"); 
mydataSmallAll_ch5 <- join_all(list(mydataAlla, mydataSmall_ch5), by="Federal_Provider_Number", type="left", match="all"); 

mydataSmallAlla <- rbind(mydataSmallAll_ch1, mydataSmallAll_ch2, mydataSmallAll_ch3, mydataSmallAll_ch4, mydataSmallAll_ch5); 
length(mydataSmallAlla$Federal_Provider_Number)/5; #9,001 NHs;

### More efficient alternative method below yields same result;
#mydataSmallAlla <- join_all(list(mydataSmall, mydataAlla), by="Federal_Provider_Number", type="inner", match="all"); 
#length(mydataSmallAlla$Federal_Provider_Number)/5; #9,001 NHs;

shapefileData <- cbind(us_aeaNH@data$Federal_Provider_Number, us_aeaNH@data$Ownership_Type);
shapefileData <- as.data.frame(shapefileData);
colnames(shapefileData) <- c("Federal_Provider_Number", "Ownership_Type");
shapefileDataAll <- join_all(list(shapefileData, mydataSmallAlla), by="Federal_Provider_Number", type="inner", match="all"); 
#--------#;

nodes_US <- length(shapefileDataAll$Number_of_Certified_Beds)/5;
sumBeds_US <- sum(shapefileDataAll$Number_of_Certified_Beds, na.rm=TRUE)/5;
sumRes_US <- sum(shapefileDataAll$Number_of_Residents_in_Certified, na.rm=TRUE)/5;
#nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
#sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
#sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(shapefileDataAll$Number_of_Residents_in_Certified/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usBeds <- aggregate(shapefileDataAll$Number_of_Certified_Beds/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

sumBeds_US <- floor(sumBeds_US);
sumRes_US <- floor(sumRes_US);
usRes$x <- floor(usRes$x);
usBeds$x <- floor(usBeds$x);
#--------#;


counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
 
 	#theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
legend1_symbols2 <- c(16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 16, 16, 16);
#img = readPNG(system.file("img", "Rlogo.png", package="png"));
p1 <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
#embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_LONG_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);


p2 <-
	ggplot(data=us50_hrr) +
	geom_map(map=us50_hrr, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
#embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_LONG_HRR_", keyDate, ".pdf", sep=""), p2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);


p3 <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group ), colour="red", fill=NA) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

##p1;
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
###embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_LONG_state_HRR_", keyDate, ".pdf", sep=""), p3, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
###ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);



p4 <-
	ggplot(data=us50) +
	#geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group ), fill="black", color="white",) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("CMS Hospital Referral Regions, 1998 to Present") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=28, hjust=0), plot.subtitle=element_text(family="FreeSans", size=16)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Minimum population of 120,000 per market (305 regions; Victoria, TX in with Corpus Christi, TX)\nAggregated from hospital service areas derived by zip code\n[Dartmouth Atlas of Health Care, 1998]", sep=""), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

##p1;
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
###embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("HRR_boundaries_US_", keyDate, ".pdf", sep=""), p4, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
###ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);





#colnames(us50_hrr);
#us50_hrr$RR1 <- factor(
#    cut(us50_hrr$model1_RR, c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, Inf)),
#    #labels = c("0 to 0.25", "0.25 to 0.5", "0.5 to 0.75", "0.75 to 1", "1 to 1.25", "1.25 to 1.5", "1.5 to 1.75", "1.75+")
#);

#-----!-----!-----!
us50$RRs1 <- us50$model1_RR * us50$model1_sig;
us50[which(us50$RRs1 == 0), ]$RRs1 <- 1;
us50$RRs2 <- us50$model2_RR * us50$model2_sig;
us50[which(us50$RRs2 == 0), ]$RRs2 <- 1;
us50_hrr$RRs1 <- us50_hrr$model1_RR * us50_hrr$model1_sig;
us50_hrr[which(us50_hrr$RRs1 == 0), ]$RRs1 <- 1;
us50_hrr$RRs2 <- us50_hrr$model2_RR * us50_hrr$model2_sig;
us50_hrr[which(us50_hrr$RRs2 == 0), ]$RRs2 <- 1;

## set non-significant to NA, so can set color to white in Viridis;
#us50[which(us50$RRs1 == 1), ]$RRs1 <- NA;
#us50[which(us50$RRs2 == 1), ]$RRs2 <- NA;
#us50_hrr[which(us50_hrr$RRs1 == 1), ]$RRs1 <- NA;
#us50_hrr[which(us50_hrr$RRs2 == 1), ]$RRs2 <- NA;

#below sets limit for colorscale;
maxValuePRR <- max(us50$RRs1, us50$RRs2, us50_hrr$RRs1, us50_hrr$RRs2, na.rm=TRUE);
minValuePRR <- min(us50$RRs1, us50$RRs2, us50_hrr$RRs1, us50_hrr$RRs2, na.rm=TRUE);
maxValuePRR <- round(maxValuePRR, digits = 2) + 0.01;
minValuePRR <- round(minValuePRR, digits = 2) - 0.01;
#us50$RRs1;

#hclplot(diverging_hcl(7, h = c(260, 0), c = 80, l = c(35, 95), power = 1))
#colorspace::diverging_hcl(n = 40, h = c(240, 15), c = c(60, 80), l = c(75, 5), power = c(1.2, 1.5), register = "Custom-Palette");


p7.1 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs1), colour="white", size=0.5) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	#geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(fill="Prevalence\nRatio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_continuous_diverging("Berlin", mid=1, h1=240, h2=15, c1=60, cmax=80, l1=75, l2=5, p1=1.2, p2=1.5, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_fill_continuous_diverging("Berlin", mid=1, c1=100, l1=100, l2=0, p1=0.01, p2=1.0, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_fill_continuous_diverging("Vik", mid=1, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	scale_fill_continuous_diverging("Berlin", mid=1, h1=200, h2=45, c1=100, l1=100, l2=0, p1=0.01, p2=1.0, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_fill_continuous_diverging("Custom-Palette", mid=1, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_fill_viridis(option="plasma", na.value="black", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 1: Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by\nAmerican State of Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50_hrr, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="white", size=0.1, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model1_LONG_RR_PS_US_", keyDate, ".pdf", sep=""), p7.1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;




#us50_hrr$RRs1;
p7.2 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs1), colour="white", size=0.1) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	#geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(fill="Prevalence\nRatio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	scale_fill_continuous_diverging("Berlin", mid=1, h1=200, h2=45, c1=100, l1=100, l2=0, p1=0.01, p2=1.0, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 1: Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by HRR\nof Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="white", size=0.5, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model1_LONG_RR_HRR_US_", keyDate, ".pdf", sep=""), p7.2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;


p7.3 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs1), colour="white", size=0.5) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	labs(fill="Prevalence\nRatio") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
 	scale_fill_continuous_diverging("Berlin", mid=1, h1=305, h2=45, c1=125, l1=110, l2=0, p1=0.01, p2=0.5, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
 	#ggtitle("Model 1: Additive Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState and HRR of Nursing Homes Accepting Medicare or\nMedicaid Funding in the United States Owned by\nOrganizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=14)) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#labs(subtitle=paste("Model 1"), caption="") +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="white", size=0.1, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	theme(legend.position="bottom") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

p7.4 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs1), colour="white", size=0.1) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	labs(fill="Prevalence\nRatio") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	scale_fill_continuous_diverging("Berlin", mid=1, h1=305, h2=45, c1=125, l1=110, l2=0, p1=0.01, p2=0.5, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#ggtitle("Model 1: Significant Prevalence Ratios Different from Unity of\nTotal Weighted Health Survey Score by HRR of Nursing\nHomes Accepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#labs(subtitle=paste("Model 1"), caption="") +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="white", size=0.5, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	theme(legend.position="bottom") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
title1 <- textGrob("Model 1: Additive Significant Prevalence Ratios\nDifferent from Unity of Total Weighted Health Survey\nScore by American State and HRR of Nursing Homes\nAccepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle1 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=12));
ggsave(paste("Model1_LONG_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p7.3, p7.4, ncol=1, top=title1, bottom=subtitle1), width=10, height=16, device=cairo_pdf); 
#https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
#get.gpar();
#-----!-----!-----!

#-----!-----!-----!
#us50$RRs2;
p8.1 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs2), colour="white", size=0.5) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	#geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(fill="Prevalence\nRatio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	scale_fill_continuous_diverging("Berlin", mid=1, h1=200, h2=45, c1=100, l1=100, l2=0, p1=0.01, p2=1.0, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 2: Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState of Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50_hrr, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="white", size=0.1, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model2_LONG_RR_PS_US_", keyDate, ".pdf", sep=""), p8.1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;




#us50_hrr$RRs2;
p8.2 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs2), colour="white", size=0.1) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	#geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(fill="Prevalence\nRatio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	scale_fill_continuous_diverging("Berlin", mid=1, h1=200, h2=45, c1=100, l1=100, l2=0, p1=0.01, p2=1.0, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 2: Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by HRR\nof Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="white", size=0.5, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model2_LONG_RR_HRR_US_", keyDate, ".pdf", sep=""), p8.2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

p8.3 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs2), colour="white", size=0.5) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	labs(fill="Prevalence\nRatio") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
 	scale_fill_continuous_diverging("Berlin", mid=1, h1=305, h2=45, c1=125, l1=110, l2=0, p1=0.01, p2=0.5, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
 	#ggtitle("Model 1: Additive Significant Prevalence Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState and HRR of Nursing Homes Accepting Medicare or\nMedicaid Funding in the United States Owned by\nOrganizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#labs(subtitle=paste("Model 2"), caption="") +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="white", size=0.1, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	theme(legend.position="bottom") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

p8.4 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs2), colour="white", size=0.1) +
	#coord_equal() + 
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	labs(fill="Prevalence\nRatio") +
	#scale_fill_viridis(option="inferno", na.value="grey50", begin=0, end=0.85, limits=c(minValuePRR, maxValuePRR)) +
	#scale_fill_gradient2(low="#130B82", mid="white", high="#80490B", midpoint=1, limit=c(minValuePRR, maxValuePRR), space="Lab", na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") +
	scale_fill_continuous_diverging("Berlin", mid=1, h1=305, h2=45, c1=125, l1=110, l2=0, p1=0.01, p2=0.5, limit=c(minValuePRR, maxValuePRR), na.value="grey50", guide="legend", breaks=rev(c(minValuePRR, 0.5, 0.75, 1, 1.25, 1.50, 1.75, 2, 2.25, maxValuePRR)), aesthetics="fill") + 
	#ggtitle("Model 1: Significant Prevalence Ratios Different from Unity of\nTotal Weighted Health Survey Score by HRR of Nursing\nHomes Accepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#labs(subtitle=paste("Model 2"), caption="") +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="white", size=0.5, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	theme(legend.position="bottom") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
title2 <- textGrob("Model 2: Additive Significant Prevalence Ratios\nDifferent from Unity of Total Weighted Health Survey\nScore by American State and HRR of Nursing Homes\nAccepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle2 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=12));
ggsave(paste("Model2_LONG_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p8.3, p8.4, ncol=1, top=title2, bottom=subtitle2), width=10, height=16, device=cairo_pdf); 
#https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
#get.gpar();
#-----!-----!-----!

title4 <- textGrob("Models 1 (Left) and 2 (Right): Additive Significant Prevalence Ratios Different from\n1.00 of Total Weighted Health Survey Score by American State (Above) and HRR (Below)\nof Nursing Homes Accepting Medicare or Medicaid Funding Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle4 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=14));
ggsave(paste("Model1andModel2_LONG_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p7.3, p8.3, p7.4, p8.4, ncol=2, top=title4, bottom=subtitle4), width=14, height=14, device=cairo_pdf); 



####### 27 June 2019 image;
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

myLegendOne <- g_legend(p8.4);
#myRectOne <- grid.rect(.5,.5,width=unit(.99,"npc"), height=unit(0.99,"npc"), gp=gpar(lwd=3, fill=NA, col="black"));
myLineOne <- linesGrob(x=unit(c(0.5, 0.5), "npc"), y=unit(c(0.1, 0.9), "npc"), gp=gpar(), vp=NULL);

pCombinedOneLegend <- grid.arrange(arrangeGrob(p7.3 + theme(legend.position="none"),
                         p8.3 + theme(legend.position="none"),
                         p7.4 + theme(legend.position="none"),
                         p8.4 + theme(legend.position="none"), 
                         nrow=2, top=title4),
             		myLegendOne, nrow=2, heights=c(10, 1));

tmpQ <- ggdraw(pCombinedOneLegend) + draw_plot(myLineOne);
ggsave(paste("Model1andModel2_OneLegend_LONG_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), tmpQ, width=14, height=14, device=cairo_pdf);            		
#ggsave(paste("Model1andModel2_OneLegend_LONG_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), pCombinedOneLegend, width=14, height=14, device=cairo_pdf); 
########;



#---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---#;
us_aeaNH_backup <- us_aeaNH;

#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("For profit - Corporation", "For profit - Partnership", "For profit - Individual", "For profit - Limited Lia")), ];

#--------#;
shapefileData <- cbind(us_aeaNH@data$Federal_Provider_Number, us_aeaNH@data$Ownership_Type);
shapefileData <- as.data.frame(shapefileData);
colnames(shapefileData) <- c("Federal_Provider_Number", "Ownership_Type");
shapefileDataAll <- join_all(list(shapefileData, mydataSmallAlla), by="Federal_Provider_Number", type="inner", match="all"); 

nodes_US <- length(shapefileDataAll$Number_of_Certified_Beds)/5;
sumBeds_US <- sum(shapefileDataAll$Number_of_Certified_Beds, na.rm=TRUE)/5;
sumRes_US <- sum(shapefileDataAll$Number_of_Residents_in_Certified, na.rm=TRUE)/5;
#nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
#sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
#sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(shapefileDataAll$Number_of_Residents_in_Certified/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usBeds <- aggregate(shapefileDataAll$Number_of_Certified_Beds/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

sumBeds_US <- floor(sumBeds_US);
sumRes_US <- floor(sumRes_US);
usRes$x <- floor(usRes$x);
usBeds$x <- floor(usBeds$x);
#--------#;


legend1_text <- c("For-profit -- Corporation", "For-profit -- Individual", 
	"For-profit -- Limited Liability company", "For-profit -- Partnership");
	
counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
 
 	#theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
legend1_symbols2 <- c(16, 16, 16, 16);
legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D");
#us_aeaNH@data$Ownership_Type

pfp <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("For-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_organizations_forprofit_statistical_", keyDate, ".pdf", sep=""), pfp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;



#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("Non profit - Other", "Non profit - Corporation", "Non profit - Church rela")), ];

#--------#;
shapefileData <- cbind(us_aeaNH@data$Federal_Provider_Number, us_aeaNH@data$Ownership_Type);
shapefileData <- as.data.frame(shapefileData);
colnames(shapefileData) <- c("Federal_Provider_Number", "Ownership_Type");
shapefileDataAll <- join_all(list(shapefileData, mydataSmallAlla), by="Federal_Provider_Number", type="inner", match="all"); 

nodes_US <- length(shapefileDataAll$Number_of_Certified_Beds)/5;
sumBeds_US <- sum(shapefileDataAll$Number_of_Certified_Beds, na.rm=TRUE)/5;
sumRes_US <- sum(shapefileDataAll$Number_of_Residents_in_Certified, na.rm=TRUE)/5;
#nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
#sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
#sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(shapefileDataAll$Number_of_Residents_in_Certified/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usBeds <- aggregate(shapefileDataAll$Number_of_Certified_Beds/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

sumBeds_US <- floor(sumBeds_US);
sumRes_US <- floor(sumRes_US);
usRes$x <- floor(usRes$x);
usBeds$x <- floor(usBeds$x);
#--------#;

legend1_text <- c("Non-profit -- Church related", "Non-profit -- Corporation",
	"Non-profit -- Other");
	
counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
 
 	#theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
legend1_symbols2 <- c(16, 16, 16);
legend1_colors <- c("#E5F5E0", "#A1D99B", "#31A354");

pnp <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Non-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_organizations_nonprofit_statistical_", keyDate, ".pdf", sep=""), pnp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;

#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("Government - County", "Government - City/county", "Government - City", "Government - Federal", "Government - State", "Government - Hospital di")), ];

#--------#;
shapefileData <- cbind(us_aeaNH@data$Federal_Provider_Number, us_aeaNH@data$Ownership_Type);
shapefileData <- as.data.frame(shapefileData);
colnames(shapefileData) <- c("Federal_Provider_Number", "Ownership_Type");
shapefileDataAll <- join_all(list(shapefileData, mydataSmallAlla), by="Federal_Provider_Number", type="inner", match="all"); 

nodes_US <- length(shapefileDataAll$Number_of_Certified_Beds)/5;
sumBeds_US <- sum(shapefileDataAll$Number_of_Certified_Beds, na.rm=TRUE)/5;
sumRes_US <- sum(shapefileDataAll$Number_of_Residents_in_Certified, na.rm=TRUE)/5;
#nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
#sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
#sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(shapefileDataAll$Number_of_Residents_in_Certified/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usBeds <- aggregate(shapefileDataAll$Number_of_Certified_Beds/5, by=list(Category=shapefileDataAll$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
#usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

sumBeds_US <- floor(sumBeds_US);
sumRes_US <- floor(sumRes_US);
usRes$x <- floor(usRes$x);
usBeds$x <- floor(usBeds$x);
#--------#;

legend1_text <- c("Government -- City", "Government -- City/county",
	"Government -- County", "Government -- Federal",    
	"Government -- Hospital district", "Government -- State");
	
counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
 
 	#theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
legend1_symbols2 <- c(15, 15, 15, 15, 15, 15);
legend1_colors <- c("#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C");

pgo <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="black", color="white", size=0.15) +
	no_ylab +
	no_xlab +
	theme_bw() +
  	theme(axis.line = element_blank(),
  		axis.ticks = element_blank(),
  		axis.text = element_blank(),
    	panel.grid.major = element_blank(),
    	panel.grid.minor = element_blank(),
    	panel.border = element_blank(),
    	panel.background = element_blank()) +
	geom_point(data=us_aeaNH@data, mapping=aes(x=coordinates(us_aeaNH)[,1], y=coordinates(us_aeaNH)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Average Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type (Average Number)", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type (Average Number)", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Government Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nAverage Number of Residents in Certified Beds:", sumRes_US, "\nAverage Number of Certified Beds:", sumBeds_US, "\nProcessing Date (Semi-Annually):", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_organizations_government_statistical_", keyDate, ".pdf", sep=""), pgo, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;

legend1_text <- c("For-profit -- Corporation", "For-profit -- Individual", 
	"For-profit -- Limited Liability company", "For-profit -- Partnership",
	"Government -- City", "Government -- City/county",
	"Government -- County", "Government -- Federal",    
	"Government -- Hospital district", "Government -- State",      
	"Non-profit -- Church related", "Non-profit -- Corporation",
	"Non-profit -- Other");
legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354");
#---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---#;


#summary(us50_hrr$mean);
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
# 13.0000  40.0000  56.0000  64.0945  78.0000 226.0000     2659 

#us50_hrr[which(us50_hrr$mean %in% 13), ]; #HRR 282;
#hrrBoundary@data[which(hrrBoundary$hrr_num %in% 282), ]; #Manchester, NH;

#us50_hrr[which(us50_hrr$mean %in% 226), ]; #HRR 388;
#hrrBoundary@data[which(hrrBoundary$hrr_num %in% 388), ]; #Bryan, TX;

#hrrBoundary@data[which(hrrBoundary$hrr_name %in% "CA - Chico"), ]; #HRR 31;
#hrrBoundary@data[which(hrrBoundary$hrr_num %in% 31), ]; #Chico, CA;
#us50_hrr[which(us50_hrr$hrrnum %in% 31), ]; #mean 225;