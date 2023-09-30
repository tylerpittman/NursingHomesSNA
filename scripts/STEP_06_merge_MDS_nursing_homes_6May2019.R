## Merges CMS MDS data to base statistical nursing home data;
## Tyler Pittman, May 2019
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_06_merge_MDS_nursing_homes_6May2019.R

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
library(dplyr);		#gives left_join function by multiple keys;
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
	processingDate <- keyPeriod[k,2];
	

#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#---------------------------- Merge to Base Data --------------------------------#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#


mdst <- read_sas(paste("mdst_fixed_", k, ".sas7bdat", sep="")); #uses tibbles automatically;
mdst$Federal_Provider_Number <- sprintf("%06s", mdst$Federal_Provider_Number); 
#mdst$Federal_Provider_Number1 <- mdst$Federal_Provider_Number;
mdst$Provider_Name1 <- mdst$Provider_Name;
mdst$Provider_Address1 <- mdst$Provider_Address;
mdst$Provider_City1 <- mdst$Provider_City;
mdst$Provider_State1 <- mdst$Provider_State;
mdst$Provider_Zip_Code1 <- mdst$Provider_Zip_Code;
mdst <- mdst[,-c(which(colnames(mdst) %in% c("Provider_Name", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code")))];
#colnames(mdst);
#str(mdst);
#mdst$Federal_Provider_Number;
#mdst$Provider_Name;
length(mdst$Federal_Provider_Number); #15655;


mydata <- fread(paste("NH_base_hierarchial_analysis_data_", keyDate, ".csv", sep=""));
##############;
##must add leading 0's to five digit bd$Federal_Provider_Number to make 6 digit length;
mydata$Federal_Provider_Number <- sprintf("%06s", mydata$Federal_Provider_Number);  
##must add leading 0's to five digit bd$Provider_Zip_Code to make 5 digit length;
mydata$Provider_Zip_Code <- sprintf("%05s", mydata$Provider_Zip_Code);  
##############;
#colnames(mydata);
str(mydata);
#mydata$Provider_Name;
length(mydata$Federal_Provider_Number); #10728 organization owner NHs with complete base statistical data;


#"Federal_Provider_Number", "Provider_Name", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code"
mds_NHs <- join_all(list(mydata, mdst), by="Federal_Provider_Number", type="left", match="first"); 
#which(is.na(mds_NHs$Provider_Name1)); #they all matched!
#mds_NHs[1,];

mds_NHs$mds_401 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score1));
mds_NHs$mds_402 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score2));
mds_NHs$mds_403 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score3));
mds_NHs$mds_404 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score4));
mds_NHs$mds_405 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score5));
mds_NHs$mds_406 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score6));
mds_NHs$mds_407 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score7));
mds_NHs$mds_408 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score8));
mds_NHs$mds_409 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score9));
mds_NHs$mds_410 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score10));
mds_NHs$mds_411 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score11));
mds_NHs$mds_415 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score12));
mds_NHs$mds_419 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score13));
mds_NHs$mds_424 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score14));
mds_NHs$mds_425 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score15));
mds_NHs$mds_426 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score16));
mds_NHs$mds_430 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score17));
mds_NHs$mds_434 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score18));
mds_NHs$mds_451 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score19));
mds_NHs$mds_452 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score20));
mds_NHs$mds_471 <- as.numeric(gsub("%", "", mds_NHs$Four_Quarter_Average_Score21));

#colnames(mds_NHs);
#"Federal_Provider_Number", "Provider_Name", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Ownership_Type", "Provider_Resides_in_Hospital", "Continuing_Care_Retirement_Commu", "Special_Focus_Facility", "Provider_Changed_Ownership_in_La", "With_a_Resident_and_Family_Counc", "Date_First_Approved_to_Provide_M", "Processing_Date", "Provider_SSA_County_Code", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Health_Inspection_Rating", "QM_Rating", "Staffing_Rating", "RN_Staffing_Rating", "Reported_CNA_Staffing_Hours_per", "Reported_LPN_Staffing_Hours_per", "Reported_RN_Staffing_Hours_per_R", "Reported_Licensed_Staffing_Hours", "Reported_Total_Nurse_Staffing_Ho", "Reported_Physical_Therapist_Staf", "Expected_CNA_Staffing_Hours_per", "Expected_LPN_Staffing_Hours_per", "Expected_RN_Staffing_Hours_per_R", "Expected_Total_Nurse_Staffing_Ho", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "Total_Weighted_Health_Survey_Sco", "Number_of_Facility_Reported_Inci", "Number_of_Substantiated_Complain", "Number_of_Fines", "Total_Amount_of_Fines_in_Dollars", "Number_of_Payment_Denials", "Total_Number_of_Penalties", "Total_Fine_Amount", "Subgroup_Louvain", "Subgroup_Simple", "Subgroup_Simple_State", "meanDegreeMultiple_State", "proportionMultiple_State", "Subgroup_Simple_HRR", "meanDegreeMultiple_HRR", "mds_401", "mds_402", "mds_403", "mds_404", "mds_405", "mds_406", "mds_407", "mds_408", "mds_409", "mds_410", "mds_411", "mds_415", "mds_419", "mds_424", "mds_425", "mds_426", "mds_430", "mds_434", "mds_451", "mds_452", "mds_471"
#mds_NHs2 <- left_join(mydataMis, mdst2, by=c("Provider_Name"="Provider_Name1", "Provider_Zip_Code"="Provider_Zip_Code1"));

#mds_NHs <- mds_NHs[, c(which(colnames(mds_NHs) %in% c("Federal_Provider_Number", "Provider_Name", "Provider_Address", "Provider_City", "Provider_State", "Provider_Zip_Code", "Ownership_Type", "Provider_Resides_in_Hospital", "Continuing_Care_Retirement_Commu", "Special_Focus_Facility", "Provider_Changed_Ownership_in_La", "With_a_Resident_and_Family_Counc", "Date_First_Approved_to_Provide_M", "Processing_Date", "Provider_SSA_County_Code", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Health_Inspection_Rating", "QM_Rating", "Staffing_Rating", "RN_Staffing_Rating", "Reported_CNA_Staffing_Hours_per", "Reported_LPN_Staffing_Hours_per", "Reported_RN_Staffing_Hours_per_R", "Reported_Licensed_Staffing_Hours", "Reported_Total_Nurse_Staffing_Ho", "Reported_Physical_Therapist_Staf", "Expected_CNA_Staffing_Hours_per", "Expected_LPN_Staffing_Hours_per", "Expected_RN_Staffing_Hours_per_R", "Expected_Total_Nurse_Staffing_Ho", "Adjusted_CNA_Staffing_Hours_per", "Adjusted_LPN_Staffing_Hours_per", "Adjusted_RN_Staffing_Hours_per_R", "Adjusted_Total_Nurse_Staffing_Ho", "Total_Weighted_Health_Survey_Sco", "Number_of_Facility_Reported_Inci", "Number_of_Substantiated_Complain", "Number_of_Fines", "Total_Amount_of_Fines_in_Dollars", "Number_of_Payment_Denials", "Total_Number_of_Penalties", "Total_Fine_Amount", "Subgroup_Louvain", "Subgroup_Simple", "Subgroup_Simple_State", "meanDegreeMultiple_State", "proportionMultiple_State", "Subgroup_Simple_HRR", "meanDegreeMultiple_HRR", "mds_401", "mds_402", "mds_403", "mds_404", "mds_405", "mds_406", "mds_407", "mds_408", "mds_409", "mds_410", "mds_411", "mds_415", "mds_419", "mds_424", "mds_425", "mds_426", "mds_430", "mds_434", "mds_451", "mds_452", "mds_471")))];
mds_NHs <- mds_NHs[, c(1:57, 147:167)];
#colnames(mds_NHs);

fwrite(mds_NHs, file = paste("NH_base_hierarchial_analysis_data_mds_", keyDate, ".csv", sep=""));



#######################################;
print(k);
k <- k + 1;

}
mclapply(1:counterPeriod, periodLoop);
#mclapply(1, periodLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl1);