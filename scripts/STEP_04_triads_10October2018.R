## Takes Medicare nursing home SAS data triads;
## Tyler Pittman, October 2018
# Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_04_triads_10October2018.R

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
	
	#k <- 3;
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

#df$lat;
#df$lon;
nhNA <- df[which(is.na(df$lat)),];
o2 <- as.data.frame(df)
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
#own <- fread(paste("Ownership_", keyDate, "_fixed.csv", sep=""));
own <- fread(paste("Ownership_", keyDate, ".csv", sep="")); #DO NOT USE _fixed.csv, IT TRUNCATES Owner_Name TO 50 CHARACTERS;
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
edgelist <- edgelist00[,colnames(edgelist00) %in% c("Federal_Provider_Number", "Owner_Name", "Owner_Type", "freq", "Provider_State", "Total_Fine_Amount", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Association_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Processing_Date")]; 

#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
edgelist <- edgelist[c(which(edgelist$Owner_Type == "Organization")), ];
edgelist_count <- count(edgelist$Federal_Provider_Number);
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 255334),];
#edgelist[32173,];
#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#

#edgelist$Owner_Name[grep("ONTARIO", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("BAY", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("REVERA", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("CHARTWELL", edgelist$Owner_Name)];
#edgelist$Owner_Name[grep("CHG", edgelist$Owner_Name)];

#edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING, LLC"),];
#edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING MANAGEMENT LLC"),];
#edgelist[which(edgelist$Owner_Name=="STRAND ADVISORS INC"),];
#edgelist[which(edgelist$Owner_Name=="COWAY"),];

sel_fac <- edgelist[which(edgelist$Owner_Name=="CHG SENIOR LIVING, LLC"),]$Federal_Provider_Number;
sel_own <- edgelist[which(edgelist$Federal_Provider_Number %in% sel_fac),]$Owner_Name;
edgelist[which(edgelist$Owner_Name %in% sel_own),];

#edgelist[which(edgelist$Owner_Name=="BAY BRIDGE CAPITAL PARTNERS, LLC"),];
#edgelist[which(edgelist$Federal_Provider_Number=="035068"),];
#edgelist[which(edgelist$Federal_Provider_Number=="055064"),];

edgelist[which(is.na(edgelist$Total_Fine_Amount)), ]$Total_Fine_Amount <- 0;
#save(edgelist, file = "edgelist.RData");
#load("edgelist.RData");
#edgelist <- as_tibble(edgelist);
nh_state_count <- count(edgelist$Provider_State)
nh_states <- unique(edgelist$Provider_State);
#nh_states <- nh_states[-c(which(is.na(nh_states)), which(nh_states == c("PR", "GU")))];
nh_states <- nh_states[-c(which(nh_states == c("PR")), which(nh_states == c("GU")))];
edgelist <- edgelist[-c(which(edgelist$Association_Date == "NO DATE PROVIDED")), ];
#edgelist <- edgelist[-c(which(is.na(edgelist$Number_of_Residents_in_Certified))), ];
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
#which(is.na(edgelist$Ownership_Type)); #none, good check;


###
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
#edgelist <- edgelist[c(which(edgelist$Provider_State == "MD")), ];
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
E(igraph)$sdTFA <- log(E(igraph)$Total_Fine_Amount);
##E(igraph)$sdTFA <- cut(E(igraph)$Total_Fine_Amount, classTFA$brks, right=FALSE);
##quantile(edgelist$Total_Fine_Amount);

V(igraph)$type <- V(igraph)$name %in% edgelist[,1];
V(igraph)$color <- V(igraph)$type + 1;
V(igraph)$color = gsub("1", "red", V(igraph)$color) #Owners will be red;
V(igraph)$color = gsub("2", "blue", V(igraph)$color) #Nursing Homes will be blue;
#V(igraph)$size = 1;
V(igraph)$size = degree(igraph)*1 #because 1 is a small size for a node, multiply it by 5;
E(igraph)$weight <- as.numeric(edgelist[,3]);
#igraph;
#igraph.bi <- bipartite.projection(igraph);
igraph.edgelist <- get.edgelist(igraph);
#plot.igraph(igraph, vertex.label=NA);
#plot(igraph, vertex.label=NA, layout=layout.fruchterman.reingold, vertex.color=V(igraph)$color);
#plot(igraph, vertex.label=NA, vertex.color=V(igraph)$color);
#plot(igraph, vertex.label=NA, edge.width=E(igraph)$weight);
#plot.igraph(igraph, layout=layout.circle);
#write.graph(igraph, file='NH_owner_name_graph.dl', format="pajek");
#write.graph(igraph, file='NH_owner_name_graph.txt', format="edgelist");
processingDate <- edgelist$Processing_Date[[1]]; 
#processingDate <- "2017-03-01"

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

#Backup highlight yellow is #FFF44F;

#V(igraph)$name[grep("REG", V(igraph)$name)];
V(igraph)$name[grep("CHG", V(igraph)$name)];
#"1397225 ONTARIO LIMITED"
#"1749877 ONTARIO LTD"

#Date_First_Approved_to_Provide_M may be variable of interest for Medicare facility;

##From lab_3.R;
#friend_comm_wt <- walktrap.community(igraph, steps=200,modularity=TRUE);
#friend_comm_wt;
#friend_comm_eb <- edge.betweenness.community(igraph);
#friend_comm_eb;
#plot(as.dendrogram(friend_comm_eb));

## let's see if we have communities here using the 
## Grivan-Newman algorithm
## 1st we calculate the edge betweenness, merges, etc...
#igraph.ebc <- edge.betweenness.community(igraph, directed=F);
#save(igraph.ebc, file = "igraphebc.RData");
##load("igraphebc.RData"); 

#save(igraph, file = "igraph.RData");
#load("igraph.RData");


#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#
#::::::::#::::::::# Overall Descriptive Section for Chapter 3 #::::::::#::::::::#

##::::::::#::::::::# This ran and executed; 16 hours, 256 GB ram, 5 cpus (comment out) #::::::::#::::::::#
igraph.bi <- bipartite.projection(igraph, remove.type = FALSE);
##igraph.bi$proj1; #facilities connected to other facilities by owner;
##igraph.bi$proj2; #owners connected to other owners by facility;

save(igraph.bi, file=paste("igraph_bi_", keyDate, ".RData", sep=''));
#load(paste("igraph_bi_", keyDate, ".RData", sep=''));

#---!--!--!--- COMMENT BLOCK ---!--!--!---#
.f = function() {
V(igraph.bi$proj1)$name[degree(igraph.bi$proj1)==max(degree(igraph.bi$proj1))];
V(igraph.bi$proj2)$name[degree(igraph.bi$proj2)==max(degree(igraph.bi$proj2))]; #"GENERAL ELECTRIC CAPITAL CORPORATION";

vcount(igraph.bi$proj1); #14620;
ecount(igraph.bi$proj1); #517560;
vcount(igraph.bi$proj2); #50114;
ecount(igraph.bi$proj2); #420441;

#count.cl.tri <- count_max_cliques(igraph.bi$proj2, min=3); #13,111 cliques with 3 or more owners;
#count.cl.tri.num <- clique_num(igraph.bi$proj2); #90 owners is clique of largest size;
#cl.tri <- cliques(igraph.bi$proj2,min=3);
cl.tri <- cliques(igraph.bi$proj2,min=3,max=3);
df <- lapply(cl.tri,function(x){V(igraph.bi$proj2)$name[x]});
df2 <- data.frame(matrix(unlist(df),ncol=3,byrow=T));
save(df2, file=paste("df_cliques_3_", keyDate, ".Rda", sep=''));
#load(paste("df_cliques_3_", keyDate, ".Rda", sep=''));
#write.graph(cl.tri, paste("cliques_3_", keyDate, sep=''), "pajek"); #Not supported for bipartite graphs;
save(cl.tri, file=paste("cliques_3_", keyDate, ".RData", sep=''));
#load(paste("cliques_3_", keyDate, ".RData", sep=''));
set.seed(1);
system.time(fwrite(df2, file = paste("df_cliques_3_", keyDate, ".csv", sep="")));
}
#---!--!--!--- COMMENT BLOCK ---!--!--!---#
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#


#transitivity(igraph.bi$proj2, type = "average"); #gives clustering coefficient;
#transitivity(igraph.bi$proj2, type = "local"); #gives clustering coefficient;
#mean(transitivity(igraph.bi$proj2, type = "local"), na.rm = TRUE);

#clos.proj2 <- closeness(igraph.bi$proj2); 
#betw.proj2 <- betweeness(igraph.bi$proj2); 

#count.cl.tri <- count_max_cliques(igraph.bi$proj2, min=3); #13,111 cliques with 3 or more owners;
#count.cl.tri.num <- clique_num(igraph.bi$proj2); #90 owners is clique of largest size;
#cl.tri <- cliques(igraph.bi$proj2,min=3);
##cl.tri <- cliques(igraph.bi$proj2,min=3,max=3);
#df <- lapply(cl.tri,function(x){V(igraph.bi$proj2)$name[x]});
#df2 <- data.frame(matrix(unlist(df),ncol=3,byrow=T));

#triangle <- graph.full(3);
##below takes a long time to run;
##has to be divided by 6 because there are six possible ways to map a triangle to itself (by permuting the vertices)
##hence each triangle will be counted six times
#igraph.triangles.proj1 <- graph.count.subisomorphisms.vf2(igraph.bi$proj1, triangle) / 6; #56,999,435 triangles;
#igraph.triangles.proj2 <- graph.count.subisomorphisms.vf2(igraph.bi$proj2, triangle) / 6; #3,203,789 triangles;

##Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
##This is sometimes also called the clustering coefficient.
#igraph.trans.proj1 <- transitivity(igraph.bi$proj1, type = c("globalundirected"), isolates = c("NaN")); #0.8922108;
#igraph.trans.proj2 <- transitivity(igraph.bi$proj2, type = c("globalundirected"), isolates = c("NaN")); #0.4355276;

##can only use below for directed graphs, junk results otherwise; 
#tri.proj1 <- triad.census(igraph.bi$proj1); 
#tri.proj2 <- triad.census(igraph.bi$proj2); 


#######

##if we wanted to use the fastgreedy.community algorithm we would do
#igraph.fc <- fastgreedy.community(igraph);
#igraph.wc <- cluster_walktrap(igraph);
#igraph.mod <- modularity(igraph.wc);
#igraph.mem <- membership(igraph.wc); #2428 groups of nursing home owners;

## Convert to object suitable for networkD3
#igraph.d3 <- igraph_to_networkD3(igraph, group = igraph.mem);

## Create force directed network plot
#forceNetwork(Links = igraph.d3$links, Nodes = igraph.d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group');

## if we wanted to use the fastgreedy.community agorithm we would do
#igraph.fc <- fastgreedy.community(igraph);
#igraph.fc.mem <- membership(igraph.fc); #2428 groups of nursing home owners;
#str(igraph.fc.mem);

#igraph.fc.d3 <- igraph_to_networkD3(igraph, group = igraph.fc.mem);
## Create force directed network plot
#forceNetwork(Links = igraph.fc.d3$links, Nodes = igraph.fc.d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group');
             
#igraph.res <- simplify(contract(igraph, membership(igraph.fc))); 
#pdf.options(encoding='CP1250');
#pdf(file=paste("NH_state_Owners_communities_SNA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#plot(igraph.res, vertex.label=NA, vertex.label.color="black", vertex.size=V(igraph)$size/3,  vertex.color=V(igraph)$color, cex=1, edge.color="black", edge.width=E(igraph)$sdTFA/3, edge.arrow.size=1.5);
#dev.off()

##com <- community.to.membership(igraph, igraph.fc$merges, steps= which.max(igraph.fc$modularity)-1);
#com <- membership(igraph.fc); #2428 groups of nursing home owners;

#igraph.adj <- get.adjacency(igraph, attr='weight');
##igraph.shortpath <- shortest.paths(igraph, mode="out");
##igraph.shortpath <- shortest.paths(igraph, mode="in");

#get.vertex.attribute(igraph, 'Provider_State');

#igraph.deg <- degree(igraph);
##igraph.clos <- closeness(igraph);
#igraph.betw <- betweenness(igraph);
#igraph.eig <- evcent(igraph);
#igraph.dwreach <- dwreach(igraph);

####################################################;
print(k);
k <- k + 1;

}
mclapply(1:counterPeriod, periodLoop);
#mclapply(1, periodLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl1);