## Takes Medicare nursing home SAS data spatial hierarchical analysis;
## Tyler Pittman, May 2019
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_07_nursing_homes_spatial_hierarchical_6May2019.R

#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#
## Order to run MedicareUS_scripts ##
## 1_nursing_homes_SNA_May2017.sas
## STEP_01_nursing_homes_geocode.R
## 1_nursing_homes_SNA_May2017.sas #Yes, run this again;
## 2_nursing_homes_SNA_May2017.sas
## STEP_02_nursing_homes_US_Census_ACS.R
## STEP_03_nursing_homes_map.R
#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#!-@-!#

------------------------------#

#toBibtex(citation("ergm"));
#toBibtex(citation("RSiena"));

setwd("/Users/tylerpittman/GitHub/NursingHomesSNA/data");

#Sys.setenv(R_MAX_NUM_DLLS = 512);
#Sys.getenv();

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
library(igraph);
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


keyPeriod <- matrix(c( 
"March2016","2016-03-01",
"September2016","2016-09-01",
"March2017","2017-03-01",
"September2017","2017-09-01",
"March2018","2018-03-01"
), ncol=2, byrow=T);


####################### All Looping starts below here ###############################;

####################### Looping starts below here ###############################;
#---
#counterPeriod <- 5;
#	periodLoop <- function(z) {
#	k <- z;
	k <- 1;
#---	
	keyDate = keyPeriod[k,1];
	processingDate <- keyPeriod[k,2];
	

#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#--------------------------- Statistical Modeling -------------------------------#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#

mydata <- fread(paste("NH_base_hierarchial_analysis_data_mds_", keyDate, ".csv", sep=""));
##############;
##must add leading 0's to five digit bd$Federal_Provider_Number to make 6 digit length;
mydata$Federal_Provider_Number <- sprintf("%06s", mydata$Federal_Provider_Number);  
##############;
mydata$HRR <- sub("\\_.*", "", mydata$Subgroup_Simple_HRR);
mydata$Subgroup_Simple_Overall_by_State <- mydata$Subgroup_Simple;
mydata$Subgroup_Simple_Overall_by_HRR <- sub('.*_', '', mydata$Subgroup_Simple_HRR);
#str(mydata);

#length(unique(mydata$HRR)); #295 HRRs with non missing data for NHs;
#mydata$QM_Rating;

#----------- DON'T DO BELOW, MISSING for key part base removed in STEP_05 ------------#;
#mydata <- na.omit(mydata); #7993 NHs for March 2016, 15 November 2018;
#-------------------------------------------------------------------------------------#;

#----------- CONVERT TO INTEGER TO USE POISSON ------------#;
#mydata$Total_Weighted_Health_Survey_Sco <- as.integer(mydata$Total_Weighted_Health_Survey_Sco);
mydata$Total_Weighted_Health_Survey_Sco <- round(mydata$Total_Weighted_Health_Survey_Sco, digits=0);
#str(mydata$Total_Weighted_Health_Survey_Sco);
#---------------------------------------------------#;


	
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

#df$lat;
#df$lon;
nhNA <- df[which(is.na(df$lat)),];
o2 <- as.data.frame(df)

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;
##### This does for NHs owned by organizations and those without missing values in statistical model #####;
o2 = o2[which(o2$Federal_Provider_Number %in% mydata$Federal_Provider_Number),];
#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;

#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS(projected);

#mydata$Subgroup_Louvain;

#summary(mydata$Total_Number_of_Penalties); #Max is 11 penalties;
#mydata$Health_Inspection_Rating; 
#summary(mydata$Health_Inspection_Rating); #Min is 1, max is 5 and mean is 2.764 but 1337 NAs;

summary(mydata$Total_Weighted_Health_Survey_Sco); #this is the outcome variable, 109 have missing values;
##Have to remove NHs with missing values of outcome variable to do proper histogram;
mydata <- mydata[which(!is.na(mydata$Total_Weighted_Health_Survey_Sco)), ];
summary(mydata$Total_Weighted_Health_Survey_Sco); #mean is total weighted deficiency score of 56.93;
sd(mydata$Total_Weighted_Health_Survey_Sco); #sd is 66.49;
overall_mean <- mean(mydata$Total_Weighted_Health_Survey_Sco); #mean is 56.93;
overall_mean <- round(overall_mean, digits=2);

histogram <-ggplot(data=mydata, aes(mydata$Total_Weighted_Health_Survey_Sco)) + 
  		theme_bw() +
  		geom_histogram(aes(y =..count..), breaks=seq(0, 1350, by = 10), col="black", fill="white") + 
  		theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
  		theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13), axis.text=element_text(family="FreeSans", size=12), axis.title=element_text(family="FreeSans", size=12)) +
  		ggtitle("Histogram of Total Weighted Health Survey Score") +
  		labs(x="Total Weighted Health Survey Score", y="Count") +
  		labs(subtitle=paste("For nursing homes with organizations as owners.", "\nProcessing Date:", processingDate), caption="");
  		
ggsave(paste("TWHSS_proj1_histogram_", keyDate, ".pdf", sep=""), histogram, width = 11, height = 6, dpi = 300, units = c("in"), scale = 1, device=cairo_pdf);
#ggsave(paste("TWHSS_proj1_histogram_", keyDate, ".pdf", sep=""), histogram, width = 11, height = 6, device=cairo_pdf);
#theme_get();


#---------------------------#;

summary(mydata$Ownership_Type);
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
#levels(mydata$Ownership_Type);
#order(mydata$Ownership_Type);


pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Ownership_Type_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(7.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Ownership_Type, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nNursing Home Ownership Type"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();


pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Ownership_Type_boxplot_cutoff_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(7.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Ownership_Type, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190));
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nNursing Home Ownership Type (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Subgroup_Simple_Overall_by_HRR);
mydata$Subgroup_Simple_Overall_by_HRR<- as.factor(mydata$Subgroup_Simple_Overall_by_HRR);
legend1_text <- levels(mydata$Subgroup_Simple_Overall_by_HRR);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_HRR_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple_Overall_by_HRR, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR, Overall by HRR (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;


summary(mydata$Subgroup_Simple);
mydata$Subgroup_Simple <- as.factor(mydata$Subgroup_Simple);
legend1_text <- levels(mydata$Subgroup_Simple);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR, Overall (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Subgroup_Simple_State);
mydata$Subgroup_Simple_State <- as.factor(mydata$Subgroup_Simple_State);
legend1_text <- levels(mydata$Subgroup_Simple_State);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_State_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple_State, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR, by Provider State (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Subgroup_Simple_HRR);
mydata$Subgroup_Simple_HRR <- as.factor(mydata$Subgroup_Simple_HRR);
legend1_text <- levels(mydata$Subgroup_Simple_HRR);
twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_Subgroup_Simple_HRR_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Subgroup_Simple_HRR, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#"1967-07-24";
mydata$Years_Business <- (as.Date(mydata$Processing_Date, "%Y-%m-%d") - as.Date(mydata$Date_First_Approved_to_Provide_M, "%Y-%m-%d"))/365.25; 
#plot(mydata$Years_Business, mydata$Total_Weighted_Health_Survey_Sco, ylim = c(0,180));
#plot(mydata$Years_Business, mydata$Total_Weighted_Health_Survey_Sco);
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
legend1_text <- levels(mydata$Subgroup_Simple_State);
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_YB_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Years_Business, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Years of CMS Funding", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Years of CMS Funding"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

twoCol <- c("grey", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("YB_proj1_Subgroup_Simple_State_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Years_Business ~ Subgroup_Simple_State, data=mydata, ylab="Years of CMS Funding", xaxt="n", xlab="", ylim=c(0,55), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Years of CMS Funding per Nursing Home by \nMultiple or Single Ownership Class per HRR, by Provider State (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();


#---------------------------#;

#mydata$Number_of_Certified_Beds;
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_number_beds_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Number_of_Certified_Beds, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Number of Certified Beds", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Number of Certified Beds"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#str(mydata);
#mydata$Number_of_Certified_Beds;
mydata$Occupancy_Rate <- (mydata$Number_of_Residents_in_Certified / mydata$Number_of_Certified_Beds); 
mydata$Occupancy_Rate <- round(mydata$Occupancy_Rate, digits = 3);
mydata[which(mydata$Occupancy_Rate > 1.0),]$Occupancy_Rate <- 1.0 ;
mydata$Occupancy_Rate <- mydata$Occupancy_Rate * 100;
summary(mydata$Occupancy_Rate); #mean is 81.6%;

#mydata[which(mydata$Occupancy_Rate > 1.0),]; #53 NHs, possible to have overcrowded NHs over capacity;
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_occupancy_rate_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Occupancy_Rate, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Occupancy Rate", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Occupancy Rate"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$QM_Rating);
mydata$QM_Rating <- as.factor(mydata$QM_Rating);
legend1_text <- levels(mydata$QM_Rating);
twoCol <- c("white", "white", "white", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_fivestar_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ QM_Rating, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,1200), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nFive-Star Quality Measure"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, QM_Rating, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "grey", "grey", "grey", "white", "white", "white", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_fivestar_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ QM_Rating*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,1200), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Five-Star Quality Measure"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Special_Focus_Facility);
mydata$Special_Focus_Facility <- as.character(mydata$Special_Focus_Facility);
mydata[which(mydata$Special_Focus_Facility == "TRUE"), ]$Special_Focus_Facility <- "Yes";
mydata[which(mydata$Special_Focus_Facility == "FALSE"), ]$Special_Focus_Facility <- "No";
mydata$Special_Focus_Facility <- as.factor(mydata$Special_Focus_Facility);
legend1_text <- levels(mydata$Special_Focus_Facility);
twoCol <- c("white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_sff_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Special_Focus_Facility, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,1200), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nSpecial Focus Facility"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Special_Focus_Facility, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_sff_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Special_Focus_Facility*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,1200), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Special Focus Facility"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Provider_Changed_Ownership_in_La);
mydata$Provider_Changed_Ownership_in_La <- as.character(mydata$Provider_Changed_Ownership_in_La);
mydata[which(mydata$Provider_Changed_Ownership_in_La == "TRUE"), ]$Provider_Changed_Ownership_in_La <- "Yes";
mydata[which(mydata$Provider_Changed_Ownership_in_La == "FALSE"), ]$Provider_Changed_Ownership_in_La <- "No";
mydata$Provider_Changed_Ownership_in_La <- as.factor(mydata$Provider_Changed_Ownership_in_La);
legend1_text <- levels(mydata$Provider_Changed_Ownership_in_La);
twoCol <- c("white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_ownership_change_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Provider_Changed_Ownership_in_La, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nProvider Changed Ownership in Past Year (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Provider_Changed_Ownership_in_La, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_ownership_change_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Provider_Changed_Ownership_in_La*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR per Ownership Change in Past Year (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Provider_Resides_in_Hospital);
mydata$Provider_Resides_in_Hospital <- as.character(mydata$Provider_Resides_in_Hospital);
mydata[which(mydata$Provider_Resides_in_Hospital == "TRUE"), ]$Provider_Resides_in_Hospital <- "Yes";
mydata[which(mydata$Provider_Resides_in_Hospital == "FALSE"), ]$Provider_Resides_in_Hospital <- "No";
mydata$Provider_Resides_in_Hospital <- as.factor(mydata$Provider_Resides_in_Hospital);
legend1_text <- levels(mydata$Provider_Resides_in_Hospital);
twoCol <- c("white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_hospital_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Provider_Resides_in_Hospital, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nProvider Reside in Hospital (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Provider_Resides_in_Hospital, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_hospital_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Provider_Resides_in_Hospital*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Provider Reside in Hospital (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$Continuing_Care_Retirement_Commu);
mydata$Continuing_Care_Retirement_Commu <- as.character(mydata$Continuing_Care_Retirement_Commu);
mydata[which(mydata$Continuing_Care_Retirement_Commu == "TRUE"), ]$Continuing_Care_Retirement_Commu <- "Yes";
mydata[which(mydata$Continuing_Care_Retirement_Commu == "FALSE"), ]$Continuing_Care_Retirement_Commu <- "No";
mydata$Continuing_Care_Retirement_Commu <- as.factor(mydata$Continuing_Care_Retirement_Commu);
legend1_text <- levels(mydata$Continuing_Care_Retirement_Commu);
twoCol <- c("white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_cont_care_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Continuing_Care_Retirement_Commu, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nContinuing Care Community (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Continuing_Care_Retirement_Commu, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_cont_care_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Continuing_Care_Retirement_Commu*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Continuing Care Community (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

summary(mydata$With_a_Resident_and_Family_Counc);
mydata$With_a_Resident_and_Family_Counc <- as.factor(mydata$With_a_Resident_and_Family_Counc);
legend1_text <- levels(mydata$With_a_Resident_and_Family_Counc);
twoCol <- c("white", "white","white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_res_council_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ With_a_Resident_and_Family_Counc, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=twoCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nResident or Family Council (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, With_a_Resident_and_Family_Counc, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "grey", "grey", "white", "white", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_res_council_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ With_a_Resident_and_Family_Counc*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.9);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Resident or Family Council (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#str(mydata);
#mydata$Adjusted_RN_Staffing_Hours_per_R;
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_adj_rn_hrd_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Adjusted_RN_Staffing_Hours_per_R, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Adjusted RN Staffing HRD", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Adjusted RN Staffing Hours per Resident Day (HRD)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#str(mydata);
#mydata$Adjusted_LPN_Staffing_Hours_per;
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_adj_lpn_hrd_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Adjusted_LPN_Staffing_Hours_per, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Adjusted LPN Staffing HRD", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Adjusted LPN Staffing Hours per Resident Day (HRD)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#str(mydata);
#mydata$Adjusted_CNA_Staffing_Hours_per;
twoCol <- c("grey", "black");
twoSymbol <- c(1, 3);
legend_text <- levels(mydata$Subgroup_Simple);
legend_symbols <- twoSymbol[mydata$Subgroup_Simple];
legend_color <- twoCol[mydata$Subgroup_Simple];
#twoCol[mydata$Subgroup_Simple];

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_adj_cna_hrd_scatterplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 0), oma=c(1,0,4,8), xpd=TRUE);
plot(mydata$Adjusted_CNA_Staffing_Hours_per, mydata$Total_Weighted_Health_Survey_Sco, ylab="Total Weighted Health Survey Score", xlab="Adjusted CNA Staffing HRD", ylim=c(0,1200), col=legend_color, pch=legend_symbols, , bty='L');
#legend("topright", legend=legend_text, title="Ownership Class per HRR", pch=twoSymbol, col=twoCol, cex=1.0, xjust=0, box.lwd=0, box.col="white", bty="o", bg="white");
par(lwd=1); legend(xpd=NA, inset=c(-0.17,0.3), "topright", box.lwd=1, title="Ownership Class per HRR", legend=legend_text, pch=twoSymbol, col=twoCol, cex=1.0, bg="transparent", xjust=0, bty="n");
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
#text(x=seq_along(legend1_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend1_text, xpd=TRUE, cex=0.4);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Scatter Plot of Years of Total Weighted Health Survey Score per Nursing Home \nby Adjusted CNA Staffing Hours per Resident Day (HRD)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#---------------------------#;

#str(mydata);
#mydata$HRR;
#mydata$Subgroup_Simple_State;
#mydata$Subgroup_Simple_HRR;
#mydata$Subgroup_Simple;

#mydata$Provider_Resides_in_Hospital;
mydata$Provider_Resides_in_Hospital <- relevel(mydata$Provider_Resides_in_Hospital, ref="No");
#cutoffs as used by Harrington, 2000;
mydata$Number_of_Certified_Beds_Cat <- cut(mydata$Number_of_Certified_Beds, c(0, 60, 120, 160, Inf));

#---------------------------#;
summary(mydata$Number_of_Certified_Beds_Cat);
mydata$Number_of_Certified_Beds_Cat <- as.character(mydata$Number_of_Certified_Beds_Cat);
mydata$Number_of_Certified_Beds_Cat <- as.factor(mydata$Number_of_Certified_Beds_Cat);
mydata$Number_of_Certified_Beds_Cat <- factor(mydata$Number_of_Certified_Beds_Cat, levels=c("(0,60]", "(60,120]", "(120,160]", "(160,Inf]"));
legend1_text <- levels(mydata$Number_of_Certified_Beds_Cat);
twoCol <- c("white", "white");

legend2_text <- levels(with(mydata, interaction(Subgroup_Simple, Number_of_Certified_Beds_Cat, sep=" - ", lex.order=TRUE)));
fourCol <- c("grey", "grey", "grey", "grey", "white", "white", "white", "white");

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("TWHSS_proj1_number_beds_category_cluster_boxplot_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=6, pointsize=12, family="FreeSans");
#par(family="FreeSans", mfrow=c(1,1), mar=c(6, 4.1, 4.1, 2.1), oma=c(1,1,3.2,16), xpd=TRUE);
par(family="FreeSans", mfrow=c(1,1), mar=c(4.5, 4.1, 0.0, 2.1), oma=c(1,0,4,0), xpd=TRUE);
boxplot(Total_Weighted_Health_Survey_Sco ~ Number_of_Certified_Beds_Cat*Subgroup_Simple, data=mydata, ylab="Total Weighted Health Survey Score", xaxt="n", xlab="", ylim=c(0,190), col=fourCol);
# x axis with ticks but without labels
#axis(1, labels = FALSE);
# Plot x labs at default x position;
text(x=seq_along(legend2_text), y=par("usr")[3]-1, srt=45, adj=1, labels=legend2_text, xpd=TRUE, cex=0.7);
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.5, line=0.5, text=paste("Boxplot of Total Weighted Health Survey Score per Nursing Home by \nMultiple or Single Ownership Class per HRR by Number of Beds by Category (Delimited)"));
#mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("For nursing homes with organizations as owners.", "\nProcessing date:", processingDate), cex=1.0, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();
#---------------------------#;


#--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--#;
#DESCRIPTIVE TABLE PART;
#--@--@--@--@--@--@--@--@--@;

#str(mydata);
#mydata$Years_Business;

#### Below works nicely for categorical variables, used for Excel input ####
#-------;
neat.table <- function(x, name){
  xx <- data.frame(x)
  names(xx) <- c("Value", "Count")
  xx$Fraction <- with(xx, Count/sum(Count))
  data.frame(Variable = name, xx)
}

x <- lapply(mydata[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydataMultiple <- mydata[which(mydata$Subgroup_Simple_Overall_by_HRR == "Multiple"), ];
x <- lapply(mydataMultiple[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR")], table)
do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

mydataSingle <- mydata[which(mydata$Subgroup_Simple_Overall_by_HRR == "Single"), ];
x <- lapply(mydataSingle[, c("Ownership_Type", "Provider_Resides_in_Hospital", "Special_Focus_Facility", "Continuing_Care_Retirement_Commu", "With_a_Resident_and_Family_Counc", "Provider_Changed_Ownership_in_La", "QM_Rating", "Subgroup_Simple_Overall_by_HRR")], table)
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

mydata2 <- mydata[, c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Years_Business", "Subgroup_Simple_Overall_by_HRR")];
#mydata2 <- mydata2[c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Subgroup_Simple_Overall_by_HRR")];
#colnames(mydata2);
#mydata2;
  
mydata2 %>%
    group_by(Subgroup_Simple_Overall_by_HRR) %>%
    mutate(id = 1:n()) %>%
    ungroup() %>%
    gather(temp, val, Total_Weighted_Health_Survey_Sco, Number_of_Residents_in_Certified, Number_of_Certified_Beds, Occupancy_Rate, Adjusted_CNA_Staffing_Hours_per, Adjusted_LPN_Staffing_Hours_per, Adjusted_RN_Staffing_Hours_per_R, Adjusted_Total_Nurse_Staffing_Ho, Years_Business, proportionMultiple_HRR, meanDegreeMultiple_HRR) %>%
    unite(temp1, temp, Subgroup_Simple_Overall_by_HRR, sep = '_') %>%
    spread(temp1, val) %>%
    select(-id) %>%
    as.data.frame() %>%
    stargazer(digits=2, type='text')
    

#https://stackoverflow.com/questions/43840086/total-of-group-by-summarize-values
#----------; 
dataGrouped <- mydata2 %>%    
   group_by(Subgroup_Simple_Overall_by_HRR) %>%   
   summarise(Count = length(Subgroup_Simple_Overall_by_HRR),
             Total = n(),
             Percent = Count/length(mydata2$Subgroup_Simple_Overall_by_HRR))
             
mydata2 %>%     
   summarise(Count = length(Subgroup_Simple_Overall_by_HRR),
             Total = n(),
             Percent = Count/length(mydata2$Subgroup_Simple_Overall_by_HRR)) %>%
   mutate( Subgroup_Simple_Overall_by_HRR = 'Total' ) %>%
   ungroup() %>%
   rbind( ungroup(dataGrouped) ) %>%
   mutate( Subgroup_Simple_Overall_by_HRR = as.factor(Subgroup_Simple_Overall_by_HRR) )
#----------;  
             

# install.packages("reporttools")  #Use this to install it, do this only once
#require(reporttools);
vars <- mydata[, c("Total_Weighted_Health_Survey_Sco", "Number_of_Residents_in_Certified", "Number_of_Certified_Beds", "Occupancy_Rate", "Years_Business", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "proportionMultiple_HRR", "meanDegreeMultiple_HRR", "Years_Business")];
group <- as.factor(mydata$Subgroup_Simple_Overall_by_HRR);
## display default statistics, only use a subset of observations, grouped analysis
result <- tableContinuous(vars=vars, group=group, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "n"), longtable=TRUE);
#cat(gsub("\\\\hline\n[^\n]+& all &[^\n]+\n", "", result)); #get rid of all category;
#cat(gsub("\\hline", "", result)); #get rid of \hline;


vars2 <- mydata[, c("mds_401", "mds_402", "mds_403", "mds_404", "mds_405", "mds_406", "mds_407", "mds_408", "mds_409", "mds_410", "mds_411", "mds_415", "mds_419", "mds_424", "mds_425", "mds_426", "mds_430", "mds_434", "mds_451", "mds_452", "mds_471")];
group2 <- as.factor(mydata$Subgroup_Simple_Overall_by_HRR);
## display default statistics, only use a subset of observations, grouped analysis
result2 <- tableContinuous(vars=vars2, group=group2, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);


#############################################################################################################
#### Might add this in later for descriptive tables of Nursing HRD and Provider State, but inconclusive #####
varsNurse <- mydata[, c("Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R")]; 
groupPS <- as.factor(mydata$Provider_State);
result3  <- tableContinuous(vars=varsNurse, group=groupPS, prec=2, cap="Table of Nursing Home Characteristics by Organization Ownership Class per HRR in HRR", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);

sumHRDPS <- summaryBy(Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R ~ Provider_State + Subgroup_Simple_Overall_by_HRR, data=mydata, FUN = function(x) { c(m = mean(x), s = sd(x)) } );
#############################################################################################################



##### This does for HRR-level explanatory variables
#str(mydata);
#length(unique(mydata$HRR));
mydataHRR <- unique(mydata[, c("HRR", "proportionMultiple_HRR", "meanDegreeMultiple_HRR")]); #should be 295 HRRs;
mydataHRR$proportionMultiple_HRR <- mydataHRR$proportionMultiple_HRR * 100;
vars3 <- mydataHRR[, c("proportionMultiple_HRR", "meanDegreeMultiple_HRR")];
result3 <- tableContinuous(vars=vars3, prec=2, cap="Table of Nursing Home Random Effect Characteristics", lab="tab1_descr_stat", stats=c("mean", "s", "na"), longtable=TRUE);
#####;

### Create HRR-Level Explanatory Variables and link to NH-Level mydata ###
#mydataHRR$proportionMultiple_HRR;
#mydataHRR$meanDegreeMultiple_HRR;
mydataHRR$proportionMultiple_HRRSc <- center_scale(mydataHRR$proportionMultiple_HRR); 
mydataHRR$meanDegreeMultiple_HRRSc <- center_scale(mydataHRR$meanDegreeMultiple_HRR);

mydataHRR$proportionMultiple_HRR_Quin <- cut(mydataHRR$proportionMultiple_HRR, breaks=c(quantile(mydataHRR$proportionMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE);
mydataHRR$meanDegreeMultiple_HRR_Quin <- cut(mydataHRR$meanDegreeMultiple_HRR, breaks=c(quantile(mydataHRR$meanDegreeMultiple_HRR, probs = seq(0, 1, by = 0.20))), labels=c("1","2","3","4","5"), include.lowest=TRUE);

mydataHRR$proportionMultiple_HRR_Bi <- cut(mydataHRR$proportionMultiple_HRR, breaks=c(quantile(mydataHRR$proportionMultiple_HRR, probs = seq(0, 1, by = 0.50))), labels=c("1","2"), include.lowest=TRUE);
mydataHRR$meanDegreeMultiple_HRR_Bi <- cut(mydataHRR$meanDegreeMultiple_HRR, breaks=c(quantile(mydataHRR$meanDegreeMultiple_HRR, probs = seq(0, 1, by = 0.50))), labels=c("1","2"), include.lowest=TRUE);


mydataHRR <- mydataHRR[, c("HRR", "proportionMultiple_HRRSc", "meanDegreeMultiple_HRRSc", "proportionMultiple_HRR_Quin", "meanDegreeMultiple_HRR_Quin", "proportionMultiple_HRR_Bi", "meanDegreeMultiple_HRR_Bi")];
mydata <- join_all(list(mydata, mydataHRR), by = 'HRR', type = "left"); 
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
#countsAffiliation_Groups <- ddply(mydata, c("HRR", "Subgroup_Louvain_HRR"), nrow, .drop=FALSE);
mydataAffiliated <- mydata[which(duplicated(mydata$Subgroup_Louvain_HRR) == TRUE), ];
length(unique(mydataAffiliated$Subgroup_Louvain_HRR)); #there are 1,868 affiliation groups; 1,870 removing commas and periods and multiple spaces;

checkAff <- mydataAffiliated[which(mydataAffiliated$HRR == 19), ]; #Little Rock, AR;
length(unique(checkAff$Subgroup_Louvain_HRR)); #7 this is correct, removing commas and periods and multiple spaces;

checkAff <- mydataAffiliated[which(mydataAffiliated$HRR == 428), ]; #Lynchburg, VA;
length(unique(checkAff$Subgroup_Louvain_HRR)); #1 this is correct, removing commas and periods and multiple spaces;

checkAff <- mydataAffiliated[which(mydataAffiliated$HRR == 119), ]; #Fort Myers, FL;
length(unique(checkAff$Subgroup_Louvain_HRR)); #5 this is correct, removing commas and periods and multiple spaces;


colnames(mydata);
mydata$Num_Affiliated_Owners;
mydataAffiliatedOwnerCount <- mydata[which(duplicated(mydata$HRR) == FALSE), ];
sum(mydataAffiliatedOwnerCount$Num_Affiliated_Owners); #there are 17,238 organization owners; 17,202 removing commas and periods and multiple spaces;
#####

#--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--@--#;
#ANALYTICAL TABLE PART;
#--@--@--@--@--@--@--@--@--@;

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
mydata$Adjusted_RN_Staffing_Hours_per_R <- center_scale(mydata$Adjusted_RN_Staffing_Hours_per_R); #0.5689037 was mean;
mydata$Adjusted_LPN_Staffing_Hours_per <- center_scale(mydata$Adjusted_LPN_Staffing_Hours_per); #1.05767 was mean;
mydata$Adjusted_CNA_Staffing_Hours_per <- center_scale(mydata$Adjusted_CNA_Staffing_Hours_per); #2.431182 was mean;

mydata$Occupancy_Rate <- center_scale(mydata$Occupancy_Rate); #81.6 was mean;

mydata$mds_401 <- center_scale(mydata$mds_401); #15.60717 was mean;
mydata$mds_402 <- center_scale(mydata$mds_402); #8.528339 was mean;
mydata$mds_403 <- center_scale(mydata$mds_403); #5.782952 was mean;
mydata$mds_404 <- center_scale(mydata$mds_404); #7.188552 was mean;
mydata$mds_405 <- center_scale(mydata$mds_405); #47.02208 was mean;
mydata$mds_406 <- center_scale(mydata$mds_406); #3.140269 was mean;
mydata$mds_407 <- center_scale(mydata$mds_407); #4.754765 was mean;
mydata$mds_408 <- center_scale(mydata$mds_408); #5.022968 was mean;
mydata$mds_409 <- center_scale(mydata$mds_409); #0.7385239 was mean;
mydata$mds_410 <- center_scale(mydata$mds_410); #3.318802 was mean;
mydata$mds_411 <- center_scale(mydata$mds_411); #93.90044 was mean;
mydata$mds_415 <- center_scale(mydata$mds_415); #92.90735 was mean;
mydata$mds_419 <- center_scale(mydata$mds_419); #17.21389 was mean;
mydata$mds_424 <- center_scale(mydata$mds_424); #17.16885 was mean;
mydata$mds_425 <- center_scale(mydata$mds_425); #1.206925 was mean;
mydata$mds_426 <- center_scale(mydata$mds_426); #79.8 was mean;
mydata$mds_430 <- center_scale(mydata$mds_430); #80.69524 was mean;
mydata$mds_434 <- center_scale(mydata$mds_434); #2.196039 was mean;
mydata$mds_451 <- center_scale(mydata$mds_451); #18.60322 was mean;
mydata$mds_452 <- center_scale(mydata$mds_452); #23.82713 was mean;
mydata$mds_471 <- center_scale(mydata$mds_471); #63.42185 was mean;

###mydata$proportionMultiple_HRR <- mydata$proportionMultiple_HRR * 100;
### Don't mean transform here, must do per HRR in above section;
##mydata$proportionMultiple_HRR <- center_scale(mydata$proportionMultiple_HRR); #64.45131 was mean;
##mydata$meanDegreeMultiple_HRR <- center_scale(mydata$meanDegreeMultiple_HRR); #4.562907 was mean;


#fit <- glmer(Total_Weighted_Health_Survey_Sco ~ Ownership_Type + Special_Focus_Facility + Provider_Changed_Ownership_in_La + Provider_Resides_in_Hospital + Adjusted_RN_Staffing_Hours_per_R + (1|Provider_State/Subgroup_Simple_State), family=negative.binomial(2), data=mydata); #theta = 2;
#summary(fit);
#pvals.fnc(fit, nsim=40000);

#prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)));
summary(mydata$Total_Weighted_Health_Survey_Sco); #mean is total weighted deficiency score of 56.97;
sd(mydata$Total_Weighted_Health_Survey_Sco); #standard deviation is 66.49, variance is 4420.88;
#must use negative binomial, data is overdispersed;
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

##how to specify fixed effect as random slope;            
##https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019896.html 

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

#model1 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + QM_Rating + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc, random=~Provider_State + HRR, family="poisson", data=mydata, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);
#save(model1, file=paste("PREM_summary_Model_1_", keyDate, ".RData", sep=""), compress="xz");
##summary(model1);
#colnames(mydata);
mydataMDS <- na.omit(mydata);

###length(mydataMDS$Federal_Provider_Number);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_Total_Nurse_Staffing_Ho + Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Certified_Beds + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=25000, nitt=100000, thin=10, pr=TRUE);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital+ Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_406 + mds_407 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=TRUE, burnin=1000, nitt=6000, thin=10, pr=TRUE);
###model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior1);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=1000, nitt=6000, thin=10, pr=TRUE, prior=prior5);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior1);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=100000, nitt=450000, thin=75, pr=TRUE, prior=prior2);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + idh(1 + proportionMultiple_HRR + meanDegreeMultiple_HRR):HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior2);
##model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);

#model2 <- MCMCglmm(Total_Weighted_Health_Survey_Sco ~ Adjusted_CNA_Staffing_Hours_per + Adjusted_LPN_Staffing_Hours_per + Adjusted_RN_Staffing_Hours_per_R + Number_of_Residents_in_Certified + Occupancy_Rate + Years_Business + Ownership_Type + Provider_Resides_in_Hospital + Special_Focus_Facility + Continuing_Care_Retirement_Commu + With_a_Resident_and_Family_Counc + Provider_Changed_Ownership_in_La + Subgroup_Simple_Overall_by_HRR + proportionMultiple_HRRSc + meanDegreeMultiple_HRRSc + mds_401 + mds_402 + mds_403 + mds_404 + mds_405 + mds_406 + mds_407 + mds_408 + mds_409 + mds_410 + mds_411 + mds_415 + mds_419 + mds_424 + mds_425 + mds_426 + mds_430 + mds_434 + mds_451 + mds_452 + mds_471, random=~Provider_State + HRR, family="poisson", data=mydataMDS, saveX=TRUE, verbose=FALSE, burnin=150000, nitt=1000000, thin=85, pr=TRUE, prior=prior5);
#save(model2, file=paste("PREM_summary_Model_2_", keyDate, ".RData", sep=""), compress="xz");
##summary(model2);

load(paste("Chapter4ModelBackup_published/PREM_summary_Model_1_", keyDate, ".RData", sep=""));
load(paste("Chapter4ModelBackup_published/PREM_summary_Model_2_", keyDate, ".RData", sep=""));
#load(paste("PREM_summary_Model_1_", keyDate, ".RData", sep=""));
#load(paste("PREM_summary_Model_2_", keyDate, ".RData", sep=""));
#summary(model1);
#summary(model2);

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

writeLines(capture.output(fe1), paste("PREM_RRs_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(fe2), paste("PREM_RRs_Model_2_", keyDate, ".txt", sep=""));
writeLines(capture.output(summary(model1)), paste("PREM_summary_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(summary(model2)), paste("PREM_summary_Model_2_", keyDate, ".txt", sep=""));

### DO THIS FOR DIAGNOSTICS OF HOW CHAIN FIT with mcmcglmm;
#http://www.maths.bath.ac.uk/~jjf23/mixchange/split.html
#colnames(model1$Sol); #1:19 are explanatory variables;
#colnames(model2$Sol); #1:38 are explanatory variables;

EVs1 <- model1$Sol[, 1:21];
EVs2 <- model2$Sol[, 1:38];
intPSs1 <- model1$Sol[, 22:72];
intPSs2 <- model2$Sol[, 39:89];
intHRRs1 <- model1$Sol[, 73:367];
intHRRs2 <- model2$Sol[, 90:383];
#slopeProportionMultHRRs1 <- model1$Sol[, 368:662];
#slopeProportionMultHRRs2 <- model2$Sol[, 384:677];
#slopeDegreeMultHRRs1 <- model1$Sol[, 663:957];
#slopeDegreeMultHRRs2 <- model2$Sol[, 678:971];
##plot(slopePropotionMultHRRs1);

length(posterior.mode(intHRRs1)); #295;  #mean intercept;
#length(posterior.mode(slopeProportionMultHRRs1)); #295; #these slopes are not much different from 0;
#length(posterior.mode(slopeDegreeMultHRRs1)); #295;  #these slopes are not much different from 0;

### https://github.com/tmalsburg/MCMCglmm-intro
###

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_VCV_", keyDate, ".pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model1$VCV)),2), mar=c(2,2,1,0));
plot(model1$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_sol_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs1)),2), mar=c(2,2,1,0));
plot(EVs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_provider_state_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPSs1)),2), mar=c(2,2,1,0));
plot(intPSs1, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model1_hrr_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=200, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intHRRs1)),2), mar=c(2,2,1,0));
plot(intHRRs1, auto.layout=F);
dev.off()

#---;

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_VCV_", keyDate, ".pdf", sep=''), bg="transparent", width=10, height=8, pointsize=12, family="FreeSans");
par(mfrow=c(length(colnames(model2$VCV)),2), mar=c(2,2,1,0));
plot(model2$VCV, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_sol_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=12, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(EVs2)),2), mar=c(2,2,1,0));
plot(EVs2, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_provider_state_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=30, pointsize=6, family="FreeSans");
par(mfrow=c(length(colnames(intPSs2)),2), mar=c(2,2,1,0));
plot(intPSs2, auto.layout=F);
dev.off()

pdf.options(encoding='CP1250');
cairo_pdf(file=paste("model2_hrr_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=200, pointsize=6, family="FreeSans");
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
    ggtitle("Caterpillar Plot of Prevalence Rate Ratio for\nProvider State with 95% Highest Posterior\nDensity Interval from Model 1") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) + 
 	labs(y="Prevalence Rate Ratio", x="Provider State") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model1_Provider_State_", keyDate, ".pdf", sep=""), pREStates1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

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
    ggtitle("Caterpillar Plot of Prevalence Rate Ratio for\nProvider State with 95% Highest Posterior\nDensity Interval from Model 2") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) + 
 	labs(y="Prevalence Rate Ratio", x="Provider State") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model2_Provider_State_", keyDate, ".pdf", sep=""), pREStates2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

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
    ggtitle("Caterpillar Plot of Prevalence Rate Ratio\nfor HRR with 95% Highest Posterior\nDensity Interval from Model 1") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.text.y=element_text(family="FreeSans", size=3)) + 
 	labs(y="Prevalence Rate Ratio", x="HRR") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model1_HRR_", keyDate, ".pdf", sep=""), pREHRRs1, width = 10, height = 14, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

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
    ggtitle("Caterpillar Plot of Prevalence Rate Ratio\nfor HRR with 95% Highest Posterior\nDensity Interval from Model 2") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.text.y=element_text(family="FreeSans", size=3)) + 
 	labs(y="Prevalence Rate Ratio", x="HRR") +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="");
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Random_intercept_Model2_HRR_", keyDate, ".pdf", sep=""), pREHRRs2, width = 10, height = 14, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

#---------------------------------

colnames(model1$VCV);
ICC_Provider_State.1 <- model1$VCV[, 1]/(rowSums(model1$VCV)); 
ICC_HRR.1 <- model1$VCV[, 2]/(rowSums(model1$VCV)); 
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
ICC_Provider_State.2 <- model2$VCV[, 1]/(rowSums(model2$VCV)); 
ICC_HRR.2 <- model2$VCV[, 2]/(rowSums(model2$VCV)); 
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

writeLines(capture.output(dft1), paste("PREM_ICCs_Model_1_", keyDate, ".txt", sep=""));
writeLines(capture.output(dft2), paste("PREM_ICCs_Model_2_", keyDate, ".txt", sep=""));


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
writeLines(capture.output(stargazer(mcmcOutputs, type="text", summary=FALSE)), paste("PREM_stargazer_Models_1_and2_", keyDate, ".txt", sep=""));

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
#--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--@@--#


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

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

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
	geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
#embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);


p2 <-
	ggplot(data=us50_hrr) +
	geom_map(map=us50_hrr, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
#embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_HRR_", keyDate, ".pdf", sep=""), p2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);


p3 <-
	ggplot(data=us50) +
	geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

##p1;
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
###embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_statistical_state_HRR_", keyDate, ".pdf", sep=""), p3, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
###ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
###ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);



p4 <-
	ggplot(data=us50) +
	#geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group ), colour="black", fill="lightblue2") +
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



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Makes piechart of nursing home residents by nursing home type
#legend1_colors;
#usLegend;

plabel <- sapply(strsplit(usLegend, ":"), "[", 1);
ptmp <- sapply(strsplit(usLegend, "\n"), "[", 2);
ptmp2 <- sapply(strsplit(ptmp, "r"), "[", 1);
pcount <- as.numeric(ptmp2);
ptmp3 <- paste(round(pcount/sum(pcount)*100, digits=1),"%", sep='');
#ptmp3 <- paste(round(pcount/sum(pcount)*100, digits=1), sep='');
#ptmp3 <- paste(pcount, ", ", round(pcount/sum(pcount)*100, digits=1),"%", sep='');
pcolor <- legend1_colors;
pieChart <- cbind(plabel, pcount, pcolor);
pieChart <- as.data.frame(pieChart);
pieChart$pcount <- as.numeric(ptmp2);
pieChart$pos = (cumsum(c(0, pieChart$pcount)) + c(pieChart$pcount / 2, .01))[1:nrow(pieChart)];

pie <- 	ggplot(pieChart, aes(1, pcount, fill = plabel)) +
		coord_polar("y", start=0) +
		scale_fill_manual(values=pcolor) +
		#geom_text(aes(label = ptmp3), position = position_stack(vjust = 0.5)) +
		geom_col(position = position_stack(reverse = TRUE)) +
		no_ylab +
		no_xlab +
		theme_bw() +
  		theme(text=element_text(family="FreeSans"),
  			axis.line = element_blank(),
  			axis.ticks = element_blank(),
  			axis.text = element_blank(),
    		panel.grid.major = element_blank(),
    		panel.grid.minor = element_blank(),
    		panel.border = element_blank(),
    		panel.background = element_blank()) +
    	ggtitle("Resident Composition of Nursing Homes \nin the United States Owned by Organizations") +
		theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 		theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 		guides(fill=guide_legend(title="Ownership Type")) + 
 		geom_text_repel(aes(x = 1.4, y = pieChart$pos, label = ptmp3), size = 6, nudge_x = .3, segment.size = .7, family="FreeSans") +
 		labs(subtitle=paste("Number of Residents in Certified Beds:", sumRes_US, "\nProcessing Date:", processingDate), caption="");
#pie;

#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, width = 10, height = 8);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, device="png", width = 10, height = 8, dpi = 100, units = c("in"), scale = 1);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".svg", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1);
ggsave(paste("Nursing_home_pie_US_organizations_statistical_", keyDate, ".pdf", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1, device=cairo_pdf);

#----------------------------------------------------------------------------------------------------------------------#;

p5 <-
	ggplot(data=us50_hrr) +
	#geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=mean), colour="black") +
	#geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=mean), colour="black") +
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
	labs(fill="Mean") +
	scale_fill_viridis(option="inferno", na.value="darkgrey") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Average Total Weighted Health Survey Score by HRR of\nNursing Homes Accepting Medicare or Medicaid Funding in\nthe United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("TotalDef_Mean_HRR_US_", keyDate, ".pdf", sep=""), p5, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;



overall_max <- max(summaryHRR$mean);
overall_min <- min(summaryHRR$mean);

p6 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=mean), colour="black", size=0.5) +
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
	labs(fill="Mean of Nursing Homes\nin Hospital Referral\nRegion (HRR)") +
	#scale_fill_viridis(option="inferno", na.value="grey50") +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=overall_mean, space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_fill_gradientn(colours=heat.colors(9, alpha = 1), na.value="grey50", breaks=c(0,56.93,226), labels=c("Minimum",56.93,"Maximum"), limits=c(0,226)) + 
	#scale_fill_gradientn(colours=rev(brewer.pal(n = 9, name = "RdBu")), na.value="grey50", breaks=c(0,overall_mean,226), labels=c("Minimum",overall_mean,"Maximum"), limits=c(0,226)) + 
	#scale_fill_gradientn(colours=rev(brewer.pal(n = 9, name = "RdBu")), na.value="grey50", breaks=c(overall_min,overall_mean,overall_max), labels=c(overall_min,overall_mean,overall_max), limits=c(overall_min,overall_max)) + 
	#guides(fill= guide_colorbar(barheight=15)) + 
	#scale_fill_gradientn(colours = c("cyan", "white", "red"), values=c(Inf, 60, 50, 40, 0), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Mean Total Weighted Health Survey Score by HRR of\nNursing Homes Accepting Medicare or Medicaid Funding in\nthe United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Overall Mean of Total Weighted Health Survey Score for Nursing Homes:", overall_mean, "\nProcessing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("TotalDef_Mean_HRR_US_", keyDate, ".pdf", sep=""), p6, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;


#colnames(us50_hrr);
#us50_hrr$RR1 <- factor(
#    cut(us50_hrr$model1_RR, c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, Inf)),
#    #labels = c("0 to 0.25", "0.25 to 0.5", "0.5 to 0.75", "0.75 to 1", "1 to 1.25", "1.25 to 1.5", "1.5 to 1.75", "1.75+")
#);

#-----!-----!-----!
us50$RRs1 <- us50$model1_RR * us50$model1_sig;
us50[which(us50$RRs1 == 0), ]$RRs1 <- 1;
#us50$RRs1;
p7.1 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs1), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="viridis", na.value="grey50") +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 1: Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by\nProvider State of Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50_hrr, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model1_RR_PS_US_", keyDate, ".pdf", sep=""), p7.1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;




us50_hrr$RRs1 <- us50_hrr$model1_RR * us50_hrr$model1_sig;
us50_hrr[which(us50_hrr$RRs1 == 0), ]$RRs1 <- 1;
#us50_hrr$RRs1;
p7.2 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs1), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="viridis", na.value="darkgrey") +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 1: Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by HRR\nof Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model1_RR_HRR_US_", keyDate, ".pdf", sep=""), p7.2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;


p7.3 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs1), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio\n(Model 1)") +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
 	#ggtitle("Model 1: Additive Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState and HRR of Nursing Homes Accepting Medicare or\nMedicaid Funding in the United States Owned by\nOrganizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

p7.4 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs1), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio\n(Model 1)") +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#ggtitle("Model 1: Significant Prevalence Rate Ratios Different from Unity of\nTotal Weighted Health Survey Score by HRR of Nursing\nHomes Accepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
title1 <- textGrob("Model 1: Additive Significant Prevalence Rate Ratios\nDifferent from Unity of Total Weighted Health Survey\nScore by Provider State and HRR of Nursing Homes\nAccepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle1 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=12));
ggsave(paste("Model1_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p7.3, p7.4, ncol=1, top=title1, bottom=subtitle1), width=10, height=16, device=cairo_pdf); 
#https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
#get.gpar();
#-----!-----!-----!

#-----!-----!-----!
us50$RRs2 <- us50$model2_RR * us50$model2_sig;
us50[which(us50$RRs2 == 0), ]$RRs2 <- 1;
#us50$RRs2;
p8.1 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs2), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="viridis", na.value="darkgrey") +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 2: Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState of Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50_hrr, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model2_RR_PS_US_", keyDate, ".pdf", sep=""), p8.1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;




us50_hrr$RRs2 <- us50_hrr$model2_RR * us50_hrr$model2_sig;
us50_hrr[which(us50_hrr$RRs2 == 0), ]$RRs2 <- 1;
#us50_hrr$RRs2;
p8.2 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs2), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio") +
	#viridis::scale_fill_viridis(option = "viridis", discrete = TRUE) +
	#scale_color_brewer(type = 'div', palette = 2, direction = 1) +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#scale_color_brewer(type="div", palette="RdBu") +
	#scale_fill_viridis(option="viridis", na.value="darkgrey") +
	#scale_fill_distiller(palette = "RdBu", limits = c(-0.5,2.5)) +
	#scale_fill_distiller(type = "seq", palette = 1, direction = -1, values = NULL, space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "fill") +
	#scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	#scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Model 2: Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by HRR\nof Nursing Homes Accepting Medicare or Medicaid\nFunding in the United States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	#guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	#guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="red") +
 	#geom_map(map=us50, aes(map_id=id, group=group), fill=NA, color="grey50", size=0.25) +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
ggsave(paste("Model2_RR_HRR_US_", keyDate, ".pdf", sep=""), p8.2, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;

p8.3 <-
	ggplot(data=us50) +
	geom_polygon(data=us50, aes(x=long, y=lat, group=group, fill=RRs2), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio\n(Model 2)") +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
 	#ggtitle("Model 1: Additive Significant Prevalence Rate Ratios Different from\nUnity of Total Weighted Health Survey Score by Provider\nState and HRR of Nursing Homes Accepting Medicare or\nMedicaid Funding in the United States Owned by\nOrganizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

p8.4 <-
	ggplot(data=us50) +
	geom_polygon(data=us50_hrr, aes(x=long, y=lat, group=group, fill=RRs2), colour="black", size=0.5) +
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
	labs(fill="Prevalence\nRate Ratio\n(Model 2)") +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1, limit=c(0, 4), space="Lab", na.value="grey50", guide="colourbar", aesthetics="fill") +
	#ggtitle("Model 1: Significant Prevalence Rate Ratios Different from Unity of\nTotal Weighted Health Survey Score by HRR of Nursing\nHomes Accepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	#labs(subtitle=paste("Processing Date:", processingDate), caption="") +
 	geom_polygon(data=us50, aes(x=long, y=lat, group=group), colour="grey50", size=0.25, alpha=0) +
 	theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;
title2 <- textGrob("Model 2: Additive Significant Prevalence Rate Ratios\nDifferent from Unity of Total Weighted Health Survey\nScore by Provider State and HRR of Nursing Homes\nAccepting Medicare or Medicaid Funding in the\nUnited States Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle2 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=12));
ggsave(paste("Model2_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p8.3, p8.4, ncol=1, top=title2, bottom=subtitle2), width=10, height=16, device=cairo_pdf); 
#https://stat.ethz.ch/R-manual/R-devel/library/grid/html/grid.text.html
#get.gpar();
#-----!-----!-----!

title4 <- textGrob("Models 1 and 2: Additive Significant Prevalence Rate Ratios Different from Unity of\nTotal Weighted Health Survey Score by Provider State and HRR of Nursing Homes\nAccepting Medicare or Medicaid Funding in the United States Owned by Organizations", gp=gpar(fontfamily="FreeSans", fontface="bold", fontsize=24, lineheight=1.0), x=0, hjust=0);
subtitle4 <- textGrob(paste("Processing Date:", processingDate), x=0, hjust=0, gp=gpar(fontfamily="FreeSans", fontsize=14));
ggsave(paste("Model1andModel2_RR_HRRandPS_US_", keyDate, ".pdf", sep=""), grid.arrange(p7.3, p8.3, p7.4, p8.4, ncol=2, top=title4, bottom=subtitle4), width=14, height=14, device=cairo_pdf); 




#---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---#;
us_aeaNH_backup <- us_aeaNH;

#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("For profit - Corporation", "For profit - Partnership", "For profit - Individual", "For profit - Limited Lia")), ];

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

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
	geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("For-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_organizations_forprofit_statistical_", keyDate, ".pdf", sep=""), pfp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;

#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("Non profit - Other", "Non profit - Corporation", "Non profit - Church rela")), ];

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

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
	geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Non-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_organizations_nonprofit_statistical_", keyDate, ".pdf", sep=""), pnp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;

#-----------------;
#unique(us_aeaNH@data$Ownership_Type);
us_aeaNH <- us_aeaNH_backup[which(us_aeaNH_backup@data$Ownership_Type %in% c("Government - County", "Government - City/county", "Government - City", "Government - Federal", "Government - State", "Government - Hospital di")), ];

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);

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
	geom_map(map=us50, aes(map_id=id, group=group), fill="white", color="black", size=0.15) +
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
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Government Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
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

us50_hrr[which(us50_hrr$mean %in% 13), ]; #HRR 282;
hrrBoundary@data[which(hrrBoundary$hrr_num %in% 282), ]; #Manchester, NH;

us50_hrr[which(us50_hrr$mean %in% 226), ]; #HRR 388;
hrrBoundary@data[which(hrrBoundary$hrr_num %in% 388), ]; #Bryan, TX;

hrrBoundary@data[which(hrrBoundary$hrr_name %in% "CA - Chico"), ]; #HRR 31;
hrrBoundary@data[which(hrrBoundary$hrr_num %in% 31), ]; #Chico, CA;
us50_hrr[which(us50_hrr$hrrnum %in% 31), ]; #mean 225;



#######################################;
print(k);
k <- k + 1;

}
#mclapply(1:counterPeriod, periodLoop);
mclapply(1, periodLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl1);