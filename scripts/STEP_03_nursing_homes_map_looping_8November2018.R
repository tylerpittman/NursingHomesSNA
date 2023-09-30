## Takes Medicare nursing home SAS data and creates adjacency matrix of owners, makes ggplot map;
## Tyler Pittman, November 2018
# R_MAX_VSIZE=30Gb Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_03_nursing_homes_map_looping_8November2018.R

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
library(stringr);	#Removes characters from strings easily;
source("shape2poly.R"); # Reads in shape2poly function and others
source("polygonizer.R"); # Reads in polygonizer function

Sys.setenv('R_MAX_VSIZE'=60000000000); #60 Gb virtual size;
Sys.getenv('R_MAX_VSIZE'); #need this to set memory to run jobs on Macbook with newer R versions;

##Only have to do this once;
#font_import();
#fonts();

##Do this to load fonts for plotting with MacOS;
#quartzFonts(avenir=c("Avenir Bookfamily="FreeSans"", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"));
#quartzFonts(ariel=c("Ariel Regular", "Ariel Italic", "Ariel Bold", "Ariel Bold Italic"));
#quartzFonts(helvetica=c("Helvetica Regular", "Helvetica Oblique", "Helvetica Bold", "Helvetica Bold Obilque"));
#quartzFonts(DejaVuSansMono=c("DejaVuSansMono", "DejaVuSansMono-Oblique", "DejaVuSansMono-Bold", "DejaVuSansMono-BoldOblique"));
#quartzFonts(garamond=c("GaramondNo8-Regular", "GaramondNo8-Italic", "GaramondNo8-Bold", "GaramondNo8-Bold-Italic"));

quartzFonts(FreeSans=c("FreeSans", "FreeSansOblique", "FreeSansBold", "FreeSansBoldOblique")); #must do this order on other computers;
#quartzFonts(FreeSans=c("FreeSans", "FreeSansBold", "FreeSansOblique", "FreeSansBoldOblique")); #must do this order on Macbook;
#par(family="FreeSans"); #do this after or during plot function call, the call to ‘quartzFonts’ above has to be executed first;
#quartzFonts();

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
own <- fread(paste("Ownership_", keyDate, ".csv", sep="")); ## DO NOT USE _fixed.csv file, TRUNCATES Owner_Name to 50 characters;
#colnames(own);
#own$Owner_Name;

## ---- 7 May 2019;
tmp0 <- str_remove(own$Owner_Name, "\\."); #removes periods;
tmp1 <- str_remove(tmp0, ","); #removes commas;
tmp2 <- str_squish(tmp1); #removes multiple spaces between characters;
own$Owner_Name <- tmp2;
## ----;

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
unique(edgelist$Owner_Type);
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 555237),];
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

sel_fac <- edgelist[which(edgelist$Owner_Name=="1749877 ONTARIO LTD"),]$Federal_Provider_Number;
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
igraph <- graph.data.frame(edgelist, directed=FALSE);
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
V(igraph)$name[grep("ONTARIO", V(igraph)$name)];
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
#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#::::::::#


#################### HAVE TO EXECUTE ABOVE EACH TIME RUNNING THIS FILE ###################


# if we wanted to use the fastgreedy.community algorithm we would do
igraph.fc <- fastgreedy.community(igraph);

#karate <- make_graph("Zachary")
#wc <- cluster_walktrap(karate)
#modularity(wc)
#membership(wc)
#plot(wc, karate)

igraph.wc <- cluster_walktrap(igraph);
igraph.mod <- modularity(igraph.wc);
igraph.mem <- membership(igraph.wc); #2428 groups of nursing home owners;
#pdf("NH_state_Owners_cluster_SNA.pdf");
#png(file=paste("NH_state_Owners_cluster_SNA_", keyDate, ".png", sep=""), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_cluster_SNA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
plot(igraph.wc, igraph, vertex.label=NA, vertex.label.color="black", vertex.size=V(igraph)$size/3, vertex.color=V(igraph)$color, cex=1, edge.color="black", edge.width=E(igraph)$sdTFA/3, edge.arrow.size=1.5);
dev.off()

## Convert to object suitable for networkD3
#igraph.d3 <- igraph_to_networkD3(igraph, group = igraph.mem);

## Create force directed network plot
#forceNetwork(Links = igraph.d3$links, Nodes = igraph.d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group');


# if we wanted to use the fastgreedy.community agorithm we would do
igraph.fc <- fastgreedy.community(igraph);
igraph.fc.mem <- membership(igraph.fc); #2428 groups of nursing home owners;
str(igraph.fc.mem);

#igraph.fc.d3 <- igraph_to_networkD3(igraph, group = igraph.fc.mem);
## Create force directed network plot
#forceNetwork(Links = igraph.fc.d3$links, Nodes = igraph.fc.d3$nodes, Source = 'source', Target = 'target', NodeID = 'name', Group = 'group');
             
igraph.res <- simplify(contract(igraph, membership(igraph.fc))); 
#pdf("NH_state_Owners_communities_SNA.pdf");
#png(file=paste("NH_state_Owners_communities_SNA_", keyDate, ".png", sep=""), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_communities_SNA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
plot(igraph.res, vertex.label=NA, vertex.label.color="black", vertex.size=V(igraph)$size/3,  vertex.color=V(igraph)$color, cex=1, edge.color="black", edge.width=E(igraph)$sdTFA/3, edge.arrow.size=1.5);
dev.off()

#com <- community.to.membership(igraph, igraph.fc$merges, steps= which.max(igraph.fc$modularity)-1);
com <- membership(igraph.fc); #2428 groups of nursing home owners;

igraph.adj <- get.adjacency(igraph, attr='weight');
#igraph.shortpath <- shortest.paths(igraph, mode="out");
#igraph.shortpath <- shortest.paths(igraph, mode="in");

get.vertex.attribute(igraph, 'Provider_State');

igraph.deg <- degree(igraph);
#igraph.clos <- closeness(igraph);
igraph.betw <- betweenness(igraph);
igraph.eig <- evcent(igraph);
#igraph.dwreach <- dwreach(igraph);

#pdf("NH_state_Owners_SNA.pdf");
#png(file=paste("NH_state_Owners_SNA_", keyDate, ".png", sep=""), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_SNA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
#plot(igraph, vertex.label=NA, layout=layout.fruchterman.reingold, vertex.size=V(igraph)$size, vertex.color=V(igraph)$color, edge.width=E(igraph)$sdTFA, edge.arrow.size=.5);
#plot(igraph, vertex.label=NA, layout=layout.fruchterman.reingold, vertex.color=V(igraph)$color, edge.width=E(igraph)$sdTFA);
plot.igraph(igraph, vertex.label=NA, vertex.label.color="black", vertex.size=V(igraph)$size/3,  vertex.color=V(igraph)$color, cex=1, edge.color="black", edge.width=E(igraph)$sdTFA/3, edge.arrow.size=1.5)
dev.off()


################################ SHAPEFILE MAPPING ################################;

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
us_aea@data$id = rownames(us_aea@data)

alaska = us_aea[us_aea$STATEFP=="02",]
	#bbox_a1 <- bbox(alaska)
	#center_a1 <- coordinates(alaska)
alaska = elide(alaska, rotate=-50)
	#bbox_a2 <- bbox(alaska)
	#center_a2 <- coordinates(alaska)
alaska = elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska = elide(alaska, shift=c(-2100000, -2500000))
proj4string(alaska) = CRS(projected)

hawaii = us_aea[us_aea$STATEFP=="15",]
hawaii = elide(hawaii, rotate=-35)
hawaii = elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) = CRS(projected)

us_aea = us_aea[!us_aea$STATEFP %in% c("02", "15"),]
us_aea = rbind(us_aea, alaska, hawaii)

us50 <- fortify(us_aea, region="STUSPS")
us50 = remove.territories(us50)
#save('helpers/us50')

p = ggplot(data=us50) + 
    geom_map(map=us50, aes(map_id=id, group=group), ,fill="white", color="dark grey", size=0.15) + 
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
img = readPNG(system.file("img", "Rlogo.png", package="png"));
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
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States") +
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
ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);


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
    	ggtitle("Resident Composition of Nursing Homes \nin the United States") +
		theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 		theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 		guides(fill=guide_legend(title="Ownership Type")) + 
 		geom_text_repel(aes(x = 1.4, y = pieChart$pos, label = ptmp3), size = 6, nudge_x = .3, segment.size = .7, family="FreeSans") +
 		labs(subtitle=paste("Number of Residents in Certified Beds:", sumRes_US, "\nProcessing Date:", processingDate), caption="");
#pie;

#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, width = 10, height = 8);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, device="png", width = 10, height = 8, dpi = 100, units = c("in"), scale = 1);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".svg", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1);
ggsave(paste("Nursing_home_pie_US_", keyDate, ".pdf", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1, device=cairo_pdf);
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#



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
 	ggtitle("For-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_forprofit_", keyDate, ".pdf", sep=""), pfp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
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
 	ggtitle("Non-profit Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_nonprofit_", keyDate, ".pdf", sep=""), pnp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
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
 	ggtitle("Government Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
ggsave(paste("Nursing_home_locations_US_government_", keyDate, ".pdf", sep=""), pgo, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
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






#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;
##### This does for NHs owned by organizations #####;
us_aeaNH <- us_aeaNH_backup; ##note that this is right here!!!!!!!!!!
us_aeaNH_org = us_aeaNH[which(us_aeaNH@data$Federal_Provider_Number %in% edgelist$Federal_Provider_Number),];

#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#;

nodes_US_org <- length(us_aeaNH_org@data$Number_of_Certified_Beds);
sumBeds_US_org <- sum(us_aeaNH_org@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US_org <- sum(us_aeaNH_org@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes_org <- aggregate(us_aeaNH_org@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH_org@data$Ownership_Type), FUN=sum);
usBeds_org <- aggregate(us_aeaNH_org@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH_org@data$Ownership_Type), FUN=sum);
usNHs_org <- aggregate(rep(1, length(us_aeaNH_org@data$Ownership_Type)), by=list(Category=us_aeaNH_org@data$Ownership_Type), FUN=sum);

counter <- 1;
while(counter < length(usNHs_org$Category)+1){
	usNHs_org$legend[counter] <- paste(legend1_text[counter], ": ", usNHs_org$x[counter], "\n  ", usRes_org$x[counter], " residents, ", usBeds_org$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend_org <- usNHs_org$legend;
 
 	#theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
legend1_symbols2 <- c(16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 16, 16, 16);
img = readPNG(system.file("img", "Rlogo.png", package="png"));
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
	geom_point(data=us_aeaNH_org@data, mapping=aes(x=coordinates(us_aeaNH_org)[,1], y=coordinates(us_aeaNH_org)[,2], color=Ownership_Type, shape=Ownership_Type, size=Number_of_Residents_in_Certified)) +  
	labs(size="Number of Residents in \nCertified Beds") +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend_org) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend_org) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle("Nursing Homes Accepting Medicare or Medicaid \nFunding in the United States Co-Owned by Organizations") +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US_org, "\nNumber of Residents in Certified Beds:", sumRes_US_org, "\nNumber of Certified Beds:", sumBeds_US_org, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".png", sep=""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8); #doesn't embed fonts right on Macbook or Cedar;
#embed_fonts(paste("Nursing_home_locations_US_", keyDate, ".pdf", sep=""));
ggsave(paste("Nursing_home_locations_US_organizations_", keyDate, ".pdf", sep=""), p1, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#ggsave(paste("Nursing_home_locations_US", ".eps", sep = ""), p1, width = 10, height = 8);
#ggsave(paste("Nursing_home_locations_US", ".jpg", sep = ""), p1, width = 10, height = 8);



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
## Makes piechart of nursing home residents by nursing home type
plabel_org <- sapply(strsplit(usLegend_org, ":"), "[", 1);
ptmp_org <- sapply(strsplit(usLegend_org, "\n"), "[", 2);
ptmp2_org <- sapply(strsplit(ptmp_org, "r"), "[", 1);
pcount_org <- as.numeric(ptmp2_org);
ptmp3_org <- paste(round(pcount_org/sum(pcount_org)*100, digits=1),"%", sep='');
pcolor <- legend1_colors;
pieChart_org <- cbind(plabel_org, pcount_org, pcolor);
pieChart_org <- as.data.frame(pieChart_org);
pieChart_org$pcount_org <- as.numeric(ptmp2_org);
pieChart_org$pos_org = (cumsum(c(0, pieChart_org$pcount_org)) + c(pieChart_org$pcount_org / 2, .01))[1:nrow(pieChart_org)];

pie <- 	ggplot(pieChart_org, aes(1, pcount_org, fill = plabel_org)) +
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
    	ggtitle("Resident Composition of Nursing Homes \nin the United States Co-Owned by Organizations") +
		theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12)) +
 		theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 		guides(fill=guide_legend(title="Ownership Type")) + 
 		geom_text_repel(aes(x = 1.4, y = pieChart_org$pos_org, label = ptmp3_org), size = 6, nudge_x = .3, segment.size = .7, family="FreeSans") +
 		labs(subtitle=paste("Number of Residents in Certified Beds:", sumRes_US_org, "\nProcessing Date:", processingDate), caption="");
#pie;

#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, width = 10, height = 8);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, device="png", width = 10, height = 8, dpi = 100, units = c("in"), scale = 1);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".svg", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1);
ggsave(paste("Nursing_home_pie_US_organizations_", keyDate, ".pdf", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1, device=cairo_pdf);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, width = 10, height = 8, dpi = 300, units = c("in"), scale = 1);
#ggsave(paste("Nursing_home_pie_US_", keyDate, ".png", sep=""), pie, device="png", width = 254, height = 203.2, dpi = 300, units = c("mm"), scale = 1);
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#devtools::session_info();
#theme_get();
#?grid::unit
#?grid::convertX
#?grid::convertY
#ggplot2:::.pt
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#





#---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---#;
us_aeaNH_backup <- us_aeaNH_org; ##note that this is right here!!!!!!!!!!

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
ggsave(paste("Nursing_home_locations_US_organizations_forprofit_", keyDate, ".pdf", sep=""), pfp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
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
ggsave(paste("Nursing_home_locations_US_organizations_nonprofit_", keyDate, ".pdf", sep=""), pnp, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
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
ggsave(paste("Nursing_home_locations_US_organizations_government_", keyDate, ".pdf", sep=""), pgo, width = 10, height = 8, device=cairo_pdf); ## this will embed fonts right on Macbook and Cedar;
#-----------------;


legend1_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D", "#E0E8FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#E5F5E0", "#A1D99B", "#31A354");
#---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---!!!---#;


##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##
## THIS PART BELOW ADDED FOR HYTPOTHETICAL REGION FOR CONCEPTUAL FRAMEWORK PLOTS USED IN DISSERTATION PROPOSAL, 28 MARCH 2018


edgelist <- edgelist[which(edgelist$Provider_State == "MS"), ];
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
igraph <- graph.data.frame(edgelist, directed=FALSE);
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

#bipartite.mapping(igraph)$type;
igraph.bi <- bipartite.projection(igraph, remove.type = FALSE);

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

##pdf("NH_state_Owners_proj2_SNA.pdf", width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#pdf(paste("NH_state_Owners_proj2_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,4,1));
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#legend("right", legend=legend2_text, pch=legend2_symbols, col="black", pt.bg=legend2_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - Owners affiliated to other owners through ", stateNH, " facilities accepting Medicare or Medicaid \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, " Cumulative years of ownership: ",sumyearsAD, " Processing date: ", processingDate), outer=TRUE);
#dev.off()
#
##pdf("NH_state_Owners_proj1_SNA.pdf", width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#pdf(paste("NH_state_Owners_proj1_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,5,1));
#plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.label=NA, edge.color="black", edge.arrow.size=1.5);
##plot(igraph.bi$proj1);
#legend("right", legend=legend1_text, pch=legend1_symbols, col="black", pt.bg=legend1_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - ", stateNH, " Facilities affiliated to other ", stateNH, " facilities \n accepting Medicare or Medicaid through owners \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, "Number of residents in certified beds: ", sumRes_proj1, ", Number of certified beds: ", sumBeds_proj1, "\n Processing date: ", processingDate), outer=TRUE);
#dev.off()

#which(is.na(V(igraph.bi$proj2)$size)); #if any of these have NAs plot below won't work (xlim error);
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#which(is.na(V(igraph.bi$proj2)$size));  
#which(is.na(V(igraph.bi$proj2)$shape));
#which(is.na(V(igraph.bi$proj2)$color));

which(is.na(V(igraph.bi$proj2)$shape));
which(is.na(V(igraph.bi$proj2)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj2)$color));
igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(is.na(V(igraph.bi$proj2)$size))]; 
igraph.bi$proj2.small2 <- igraph.bi$proj2.small - V(igraph.bi$proj2.small)[which(is.na(V(igraph.bi$proj2.small)$color))];
igraph.bi$proj2.small3 <- igraph.bi$proj2.small2 - V(igraph.bi$proj2.small2)[which(is.na(V(igraph.bi$proj2.small2)$shape))];
igraph.bi$proj2 <-igraph.bi$proj2.small3;

#postscript(paste("NH_state_Owners_proj2_SNA_", stateNH, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj2_SNA_", "hypothetical_region_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_", "hypothetical_region_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Hypothetical region owners affiliated to each other \nthrough", "facilities accepting Medicare or Medicaid"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

#postscript(paste("NH_state_Owners_proj1_SNA_", stateNH, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj1_SNA_", "hypothetical_region_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_", "hypothetical_region_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Hypothetical region facilities accepting Medicare or \nMedicaid affiliated to each other through owners"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();
#-------------------------#;
##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##


####################### All Looping starts below here ###############################;
edgelist <- edgelist00[,colnames(edgelist00) %in% c("Federal_Provider_Number", "Owner_Name", "Owner_Type", "freq", "Provider_State", "Total_Fine_Amount", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Association_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Processing_Date")]; 

#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
edgelist <- edgelist[c(which(edgelist$Owner_Type == "Organization")), ];
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 555237),];
#edgelist[32173,];
#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#

####################### Looping starts below here ###############################;
#counterState <- 1;
counterState <- length(unique(edgelist$Provider_State));
stateLoop <- function(y) {
	j <- y;
	#j <- 10;
	#j <- 1;	
	
###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###
tryCatch({

st = nh_states[j];
edgelist <- edgelist[which(edgelist$Provider_State == st), ];
#edgelist <- edgelist[which(edgelist$Provider_State == "WY"), ];
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
igraph <- graph.data.frame(edgelist, directed=FALSE);
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

#E(igraph)$sdTFA <- log(E(igraph)$Total_Fine_Amount);

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

#bipartite.mapping(igraph)$type;
#igraph.bi <- bipartite.projection(igraph, remove.type = FALSE);
igraph.bi <- bipartite.projection(igraph, multiplicity = TRUE, remove.type = TRUE);  #Multiplicity gives weight of ties, might use for triad join part;

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

##pdf("NH_state_Owners_proj2_SNA.pdf", width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#pdf(paste("NH_state_Owners_proj2_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,4,1));
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#legend("right", legend=legend2_text, pch=legend2_symbols, col="black", pt.bg=legend2_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - Owners affiliated to other owners through ", stateNH, " facilities accepting Medicare or Medicaid \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, " Cumulative years of ownership: ",sumyearsAD, " Processing date: ", processingDate), outer=TRUE);
#dev.off()
#
##pdf("NH_state_Owners_proj1_SNA.pdf", width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#pdf(paste("NH_state_Owners_proj1_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,5,1));
#plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.label=NA, edge.color="black", edge.arrow.size=1.5);
##plot(igraph.bi$proj1);
#legend("right", legend=legend1_text, pch=legend1_symbols, col="black", pt.bg=legend1_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - ", stateNH, " Facilities affiliated to other ", stateNH, " facilities \n accepting Medicare or Medicaid through owners \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, "Number of residents in certified beds: ", sumRes_proj1, ", Number of certified beds: ", sumBeds_proj1, "\n Processing date: ", processingDate), outer=TRUE);
#dev.off()

#which(is.na(V(igraph.bi$proj2)$size)); #if any of these have NAs plot below won't work (xlim error);
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#which(is.na(V(igraph.bi$proj2)$size));  
#which(is.na(V(igraph.bi$proj2)$shape));
#which(is.na(V(igraph.bi$proj2)$color));


#length(V(igraph.bi$proj2));
#length(E(igraph.bi$proj2));
#which(E(igraph.bi$proj2)$weight > 2);

#ig.layout <- layout.fruchterman.reingold(igraph.bi$proj2);
#which(is.na(ig.layout));
#layout <- layout.norm(ig.layout, -1, 1, -1, 1);
#l <- layout_with_fr(igraph.bi$proj2); 
#plot(igraph.bi$proj2, layout=l);
#plot(igraph.bi$proj2, layout=layout, xlim = c(-1, 1));
#is_igraph(igraph.bi$proj2);


#igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(V(igraph.bi$proj2)$size < 15)]; # another way to delete vertices for less than 15 NHs ties;
#V(igraph.bi$proj2.m)$color;
#is.na(E(igraph.bi$proj2.m));
#is.na(igraph.bi$proj2.m);
#plot(igraph.bi$proj2.m);

#str(V(igraph.bi$proj2.m));
#V(igraph.bi$proj2.m)$size;
#V(igraph.bi$proj2.m)$color;
#V(igraph.bi$proj2.m)$shape;


#ERROR : need finite 'xlim' values 
#Error in igraph.check.shapes(params("vertex", "shape")) : 
#  Bad vertex shape(s): NA.
#Calls: plot -> plot -> plot.igraph -> igraph.check.shapes
#Execution halted
#igraph.check.shapes() ##this is error function for proj2;

which(is.na(V(igraph.bi$proj2)$shape));
which(is.na(V(igraph.bi$proj2)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj2)$color));
igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(is.na(V(igraph.bi$proj2)$size))]; 
igraph.bi$proj2.small2 <- igraph.bi$proj2.small - V(igraph.bi$proj2.small)[which(is.na(V(igraph.bi$proj2.small)$color))];
igraph.bi$proj2.small3 <- igraph.bi$proj2.small2 - V(igraph.bi$proj2.small2)[which(is.na(V(igraph.bi$proj2.small2)$shape))];
igraph.bi$proj2 <-igraph.bi$proj2.small3;

}, error=function(e){})
###@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!@@@!!!!###

tryCatch({	

#postscript(paste("NH_state_Owners_proj2_SNA_", stateNH, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj2_SNA_", stateNH, "_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", stateNH, "facilities accepting Medicare or Medicaid"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #allows job to contiune loop with error;


#postscript(paste("NH_state_Owners_proj1_SNA_", stateNH, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj1_SNA_", stateNH, "_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_", stateNH, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", stateNH, "facilities affiliated to other", stateNH, "\nfacilities accepting Medicare or Medicaid through owners"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();
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



#-----------------------# FOR ENTIRE USA #-----------------------#
edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
igraph <- graph.data.frame(edgelist, directed=FALSE);
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

#bipartite.mapping(igraph)$type;
#igraph.bi <- bipartite.projection(igraph, remove.type = FALSE);
igraph.bi <- bipartite.projection(igraph, remove.type=TRUE);

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

##pdf("NH_state_Owners_proj2_SNA.pdf", width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#pdf(paste("NH_state_Owners_proj2_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,4,1));
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#legend("right", legend=legend2_text, pch=legend2_symbols, col="black", pt.bg=legend2_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - Owners affiliated to other owners through ", stateNH, " facilities accepting Medicare or Medicaid \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, " Cumulative years of ownership: ",sumyearsAD, " Processing date: ", processingDate), outer=TRUE);
#dev.off()
#
##pdf("NH_state_Owners_proj1_SNA.pdf", width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#pdf(paste("NH_state_Owners_proj1_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,5,1));
#plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.label=NA, edge.color="black", edge.arrow.size=1.5);
##plot(igraph.bi$proj1);
#legend("right", legend=legend1_text, pch=legend1_symbols, col="black", pt.bg=legend1_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - ", stateNH, " Facilities affiliated to other ", stateNH, " facilities \n accepting Medicare or Medicaid through owners \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, "Number of residents in certified beds: ", sumRes_proj1, ", Number of certified beds: ", sumBeds_proj1, "\n Processing date: ", processingDate), outer=TRUE);
#dev.off()

#which(is.na(V(igraph.bi$proj2)$size)); #if any of these have NAs plot below won't work (xlim error);
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#which(is.na(V(igraph.bi$proj2)$size));  
#which(is.na(V(igraph.bi$proj2)$shape));
#which(is.na(V(igraph.bi$proj2)$color));
#which(is.na(V(igraph.bi$proj1)$size));  
#which(is.na(V(igraph.bi$proj1)$shape));
#which(is.na(V(igraph.bi$proj1)$color));
#edgelist.bi.proj1$Number_of_Residents_in_Certified[which(is.na(V(igraph.bi$proj1)$size))];

which(is.na(V(igraph.bi$proj2)$shape));
which(is.na(V(igraph.bi$proj2)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj2)$color));
igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(is.na(V(igraph.bi$proj2)$size))]; 
igraph.bi$proj2.small2 <- igraph.bi$proj2.small - V(igraph.bi$proj2.small)[which(is.na(V(igraph.bi$proj2.small)$color))];
igraph.bi$proj2.small3 <- igraph.bi$proj2.small2 - V(igraph.bi$proj2.small2)[which(is.na(V(igraph.bi$proj2.small2)$shape))];
igraph.bi$proj2 <-igraph.bi$proj2.small3;

which(is.na(V(igraph.bi$proj1)$shape));
which(is.na(V(igraph.bi$proj1)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj1)$color));
igraph.bi$proj1.small <- igraph.bi$proj1 - V(igraph.bi$proj1)[which(is.na(V(igraph.bi$proj1)$size))]; 
igraph.bi$proj1.small2 <- igraph.bi$proj1.small - V(igraph.bi$proj1.small)[which(is.na(V(igraph.bi$proj1.small)$color))];
igraph.bi$proj1.small3 <- igraph.bi$proj1.small2 - V(igraph.bi$proj1.small2)[which(is.na(V(igraph.bi$proj1.small2)$shape))];
igraph.bi$proj1 <-igraph.bi$proj1.small3;

tryCatch({

#postscript(paste("NH_state_Owners_proj2_SNA_USA", ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj2_SNA_", "USA_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_", "USA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership"), pch=c(legend2_symbols, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", "USA", "facilities accepting Medicare or Medicaid"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #allows job to contiune loop with error;

tryCatch({

#postscript(paste("NH_state_Owners_proj1_SNA_", "USA", ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj1_SNA_", "USA_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_", "USA_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds"), pch=c(legend1_symbols, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", "USA", "facilities affiliated to other", "USA", "\nfacilities accepting Medicare or Medicaid through owners"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
dev.off();

}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) #allows job to contiune loop with error;
#-----------------------# FOR ENTIRE USA #-----------------------#



#-----------------------# FOR pri_own #-----------------------#
edgelist <- edgelist00[,colnames(edgelist00) %in% c("Federal_Provider_Number", "Owner_Name", "Owner_Type", "freq", "Provider_State", "Total_Fine_Amount", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Association_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Processing_Date")]; 

#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
edgelist <- edgelist[c(which(edgelist$Owner_Type == "Organization")), ];
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 555237),];
#edgelist[32173,];
#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#

edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
pri_own <- "BAY BRIDGE CAPITAL PARTNERS LLC";
#pri_own <- "1397225 ONTARIO LIMITED";
sel_fac <- edgelist[which(edgelist$Owner_Name %in% pri_own),]$Federal_Provider_Number;
sel_own <- edgelist[which(edgelist$Federal_Provider_Number %in% sel_fac),]$Owner_Name;
edgelist <- edgelist[which(edgelist$Owner_Name %in% sel_own),];

igraph <- graph.data.frame(edgelist, directed=FALSE);
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
#E(igraph)$matchRetVal_own <- as.numeric(edgelist$matchRetVal_own);
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

#bipartite.mapping(igraph)$type;
### NOTE THAT BIPARTITE DOESN'T KEEP EDGE ATTRIBUTES EXCEPT WEIGHT, ONLY VERTEX ATTRIBUTES;
#igraph.bi <- bipartite.projection(igraph, remove.type = FALSE);
igraph.bi <- bipartite.projection(igraph, remove.type=TRUE);



#-------------------------#;
key_fac <- matrix(c( 
"For profit - Corporation","#FEE5D9","circle",16,
"For profit - Individual","#FCAE91","circle",16,
"For profit - Limited Lia","#FB6A4A","circle",16,
"For profit - Partnership","#CB181D","circle",16,
"Government - City","#E0E8FF","square",15,
"Government - City/county","#C6DBEF","square",15,
"Government - County","#9ECAE1","square",15,
"Government - Federal","#6BAED6","square",15,
"Government - Hospital di","#3182BD","square",15,
"Government - State","#08519C","square",15,
"Non profit - Church rela","#E5F5E0","circle",16,
"Non profit - Corporation","#A1D99B","circle",16,
"Non profit - Other","#31A354","circle",16,
NA,"#FFFFFF","square",15
), ncol=4, byrow=T);
key_fac <- as.data.frame(key_fac);
colnames(key_fac) <- c("Ownership_Type", "color", "shape", "shape1");
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

##pdf("NH_state_Owners_proj2_SNA.pdf", width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#pdf(paste("NH_state_Owners_proj2_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,4,1));
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#legend("right", legend=legend2_text, pch=legend2_symbols, col="black", pt.bg=legend2_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - Owners affiliated to other owners through ", stateNH, " facilities accepting Medicare or Medicaid \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, " Cumulative years of ownership: ",sumyearsAD, " Processing date: ", processingDate), outer=TRUE);
#dev.off()
#
##pdf("NH_state_Owners_proj1_SNA.pdf", width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#pdf(paste("NH_state_Owners_proj1_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,5,1));
#plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.label=NA, edge.color="black", edge.arrow.size=1.5);
##plot(igraph.bi$proj1);
#legend("right", legend=legend1_text, pch=legend1_symbols, col="black", pt.bg=legend1_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - ", stateNH, " Facilities affiliated to other ", stateNH, " facilities \n accepting Medicare or Medicaid through owners \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, "Number of residents in certified beds: ", sumRes_proj1, ", Number of certified beds: ", sumBeds_proj1, "\n Processing date: ", processingDate), outer=TRUE);
#dev.off()

#which(is.na(V(igraph.bi$proj2)$size)); #if any of these have NAs plot below won't work (xlim error);
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#which(is.na(V(igraph.bi$proj2)$size));  
#which(is.na(V(igraph.bi$proj2)$shape));
#which(is.na(V(igraph.bi$proj2)$color));
#which(is.na(V(igraph.bi$proj1)$size));  
#which(is.na(V(igraph.bi$proj1)$shape));
#which(is.na(V(igraph.bi$proj1)$color));
#edgelist.bi.proj1$Number_of_Residents_in_Certified[which(is.na(V(igraph.bi$proj1)$size))];

colHighlight <- "#450045";
V(igraph.bi$proj2)$label <- NA;
V(igraph.bi$proj1)$label <- NA;
V(igraph.bi$proj2)$label[which(V(igraph.bi$proj2)$name %in% pri_own)] <- "\u2020";
V(igraph.bi$proj1)$label[which(V(igraph.bi$proj1)$name %in% sel_fac)] <- "\u2020";
#* is "\u002A";
#"\u2020" is a dagger;


#get.vertex.attribute(igraph.bi$proj2);
#get.edge.attribute(igraph.bi$proj2);

which(is.na(V(igraph.bi$proj2)$shape));
which(is.na(V(igraph.bi$proj2)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj2)$color));
igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(is.na(V(igraph.bi$proj2)$size))]; 
igraph.bi$proj2.small2 <- igraph.bi$proj2.small - V(igraph.bi$proj2.small)[which(is.na(V(igraph.bi$proj2.small)$color))];
igraph.bi$proj2.small3 <- igraph.bi$proj2.small2 - V(igraph.bi$proj2.small2)[which(is.na(V(igraph.bi$proj2.small2)$shape))];
igraph.bi$proj2 <-igraph.bi$proj2.small3;


#postscript(paste("NH_state_Owners_proj2_SNA", pri_own, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj2_SNA_", pri_own, "_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_", pri_own, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=V(igraph.bi$proj2)$label, vertex.label.family="FreeSans", vertex.label.cex=3, vertex.label.color=colHighlight, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership", paste("\u2020 Denotes ", length(pri_own), " owner(s) of interest", sep='')), pch=c(legend2_symbols, NA, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA, colHighlight), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", "facilities accepting Medicare or Medicaid"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
mtext(text=paste("Owner of interest:", pri_own), cex=1.2, line=-2, side=3, adj=1, outer=TRUE);
dev.off();

#postscript(paste("NH_state_Owners_proj1_SNA_", pri_own, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj1_SNA_", pri_own, "_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_", pri_own, "_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=V(igraph.bi$proj1)$label, vertex.label.family="FreeSans", vertex.label.cex=3, vertex.label.color=colHighlight, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviationcex=1.2, line=-2,", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds", NA, paste("\u2020 Denotes ", length(sel_fac), " owner of interest\n   NH(s)", sep='')), pch=c(legend1_symbols, NA, NA, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", "facilities affiliated to other", "\nfacilities accepting Medicare or Medicaid through owners"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
mtext(text=paste("Owner of interest:", pri_own), cex=1.2, line=-2, side=3, adj=1, outer=TRUE);
dev.off();
#-----------------------# FOR pri_own #-----------------------#

################################ SHAPEFILE MAPPING FOR pri_own ################################;
sel_pri_own <- unique(edgelist$Federal_Provider_Number);
df1 <- df[which(df$Federal_Provider_Number %in% sel_pri_own),];
#df1$Ownership_Type; #sometimes is missing some for some reason;
#df1$Federal_Provider_Number;
#edgelist.bi.proj1$Federal_Provider_Number;
key <- edgelist.bi.proj1;
key <- key[sort.list(key$Federal_Provider_Number),];
#df1$Federal_Provider_Number;
#key$Federal_Provider_Number;

lat <- df1$lat;
lon <- df1$lon;
tmpDF <- cbind(key, lat, lon); #have to do this to get correct Node attributes (i.e. Ownership_Type)
#edgelist.bi.proj1$Number_of_Residents_in_Certified;

o2 <- as.data.frame(tmpDF);
#o2 <- as.data.frame(df1);

#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS(projected);


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
us_aea@data$id = rownames(us_aea@data)

alaska = us_aea[us_aea$STATEFP=="02",]
	#bbox_a1 <- bbox(alaska)
	#center_a1 <- coordinates(alaska)
alaska = elide(alaska, rotate=-50)
	#bbox_a2 <- bbox(alaska)
	#center_a2 <- coordinates(alaska)
alaska = elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska = elide(alaska, shift=c(-2100000, -2500000))
proj4string(alaska) = CRS(projected)

hawaii = us_aea[us_aea$STATEFP=="15",]
hawaii = elide(hawaii, rotate=-35)
hawaii = elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) = CRS(projected)

us_aea = us_aea[!us_aea$STATEFP %in% c("02", "15"),]
us_aea = rbind(us_aea, alaska, hawaii)

us50 <- fortify(us_aea, region="STUSPS")
us50 = remove.territories(us50)
#save('helpers/us50')

p = ggplot(data=us50) + 
    geom_map(map=us50, aes(map_id=id, group=group), ,fill="white", color="dark grey", size=0.15) + 
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
#unique(us_aeaNH@data$Ownership_Type);
##----------##;

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);


tmpType1 = factor(usNHs$Category, 
	levels = c("For profit - Corporation", "For profit - Individual", 
	"For profit - Limited Lia", "For profit - Partnership",
	"Government - City", "Government - City/county",
	"Government - County", "Government - Federal",    
	"Government - Hospital di", "Government - State",      
	"Non profit - Church rela", "Non profit - Corporation",
	"Non profit - Other", NA));
legend1_colors <- as.character(key_fac$color[tmpType1]);
legend1_symbols2 <- as.numeric(as.character(key_fac$shape1[tmpType1]));
legend1_sel_text <- legend1_text[tmpType1];

counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_sel_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
facLegend <- "\u2020 Denotes owner of interest";

 	#theme(plot.title=element_text(face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

singledata <- us_aeaNH[which(us_aeaNH$Federal_Provider_Number %in% sel_fac),];

	#scale_colour_manual("Federal Provider Number", values=colHighlight, labels=facLegend) +
 	#scale_shape_manual("Federal Provider Number", values=8, labels=facLegend) +
#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
#legend1_symbols2 <- c(16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 16, 16, 16);
img = readPNG(system.file("img", "Rlogo.png", package="png"));
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
	labs(size="Number of Residents in \nCertified Beds", x=paste("\u2020 Denotes ", length(sel_fac), " owner of interest NH(s)", sep='')) +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle(paste("Nursing Homes Accepting Medicare or Medicaid \nAffiliated by Ownership Association to\n", pri_own, sep='')) +
 	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.title=element_text(family="FreeSans", size=12)) +
 	#theme(plot.title=element_text(face="bold", size=24, hjust=0), plot.subtitle=element_text(size=12), axis.title=element_text(size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	geom_point(data=singledata@data, mapping=aes(x=coordinates(singledata)[,1], y=coordinates(singledata)[,2]), color=colHighlight, shape="†", size=6) + 
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US", pri_own, "_", keyDate, ".png", sep = ""), p1, width = 10, height = 8);
ggsave(paste("Nursing_home_locations_US", pri_own, "_", keyDate, ".pdf", sep = ""), p1, width = 10, height = 8, device=cairo_pdf);



##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##
## THIS PART BELOW ADDED FOR HYPOTHETICAL OWNER FOR CONCEPTUAL FRAMEWORK PLOTS USED IN DISSERTATION PROPOSAL, 28 MARCH 2018

#-----------------------# FOR "Hypothetical Owner" #-----------------------#
edgelist <- edgelist00[,colnames(edgelist00) %in% c("Federal_Provider_Number", "Owner_Name", "Owner_Type", "freq", "Provider_State", "Total_Fine_Amount", "Ownership_Type", "Date_First_Approved_to_Provide_M", "Association_Date", "Number_of_Certified_Beds", "Number_of_Residents_in_Certified", "Overall_Rating", "Processing_Date")]; 

#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#
### DO THIS FOR A SIMPLER DATASET TO CHECK CODE;
###
edgelist <- edgelist[c(which(edgelist$Owner_Type == "Organization")), ];
#edgelist$freq;
#edgelist[which(edgelist$Federal_Provider_Number %in% 555237),];
#edgelist[32173,];
#-------------------------------#-----------------------------------#
#-------------------------------#-----------------------------------#

edgelist <- edgelist[order(edgelist$Ownership_Type, edgelist$Owner_Type),];
#pri_own <- "BAY BRIDGE CAPITAL PARTNERS, LLC";
#pri_own <- "1397225 ONTARIO LIMITED";
#pri_own <- "OAKLEAF HOLDING LLC";
#pri_own <- "EXLEY, TOMMY";
pri_own <- "LONG TERM CARE MANAGEMENT LLC";
sel_fac <- edgelist[which(edgelist$Owner_Name %in% pri_own),]$Federal_Provider_Number;
sel_own <- edgelist[which(edgelist$Federal_Provider_Number %in% sel_fac),]$Owner_Name;
edgelist <- edgelist[which(edgelist$Owner_Name %in% sel_own),];
length(edgelist$Federal_Provider_Number);

igraph <- graph.data.frame(edgelist, directed=FALSE);
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

#bipartite.mapping(igraph)$type;
E(igraph)$Ownership_Type;

igraph.bi <- bipartite.projection(igraph, remove.type=TRUE);

#-------------------------#;
key_fac <- matrix(c( 
"For profit - Corporation","#FEE5D9","circle",16,
"For profit - Individual","#FCAE91","circle",16,
"For profit - Limited Lia","#FB6A4A","circle",16,
"For profit - Partnership","#CB181D","circle",16,
"Government - City","#E0E8FF","square",15,
"Government - City/county","#C6DBEF","square",15,
"Government - County","#9ECAE1","square",15,
"Government - Federal","#6BAED6","square",15,
"Government - Hospital di","#3182BD","square",15,
"Government - State","#08519C","square",15,
"Non profit - Church rela","#E5F5E0","circle",16,
"Non profit - Corporation","#A1D99B","circle",16,
"Non profit - Other","#31A354","circle",16,
NA,"#FFFFFF","square",15
), ncol=4, byrow=T);
key_fac <- as.data.frame(key_fac);
colnames(key_fac) <- c("Ownership_Type", "color", "shape", "shape1");
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

##pdf("NH_state_Owners_proj2_SNA.pdf", width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#pdf(paste("NH_state_Owners_proj2_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some owners are affiliated to other owners through facilities;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,4,1));
#plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
#legend("right", legend=legend2_text, pch=legend2_symbols, col="black", pt.bg=legend2_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - Owners affiliated to other owners through ", stateNH, " facilities accepting Medicare or Medicaid \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, " Cumulative years of ownership: ",sumyearsAD, " Processing date: ", processingDate), outer=TRUE);
#dev.off()
#
##pdf("NH_state_Owners_proj1_SNA.pdf", width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#pdf(paste("NH_state_Owners_proj1_SNA_", stateNH, ".pdf", sep=''), width=10, height=7.5); ##some facilities are affiliated to other facilities through owners;
#par(xpd=F, mar=c(0,0,0,1), oma=c(0,0,5,1));
#plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=NA, edge.color="grey", edge.arrow.size=1.5);
##plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.label=NA, edge.color="black", edge.arrow.size=1.5);
##plot(igraph.bi$proj1);
#legend("right", legend=legend1_text, pch=legend1_symbols, col="black", pt.bg=legend1_colors, cex=1.0, bty = "n", xjust=0,);
#title(paste("Bipartite Projection - ", stateNH, " Facilities affiliated to other ", stateNH, " facilities \n accepting Medicare or Medicaid through owners \n", "Nodes: ", nodes_proj1, ", Edges: ", edges_proj1, "Number of residents in certified beds: ", sumRes_proj1, ", Number of certified beds: ", sumBeds_proj1, "\n Processing date: ", processingDate), outer=TRUE);
#dev.off()

#which(is.na(V(igraph.bi$proj2)$size)); #if any of these have NAs plot below won't work (xlim error);
#edgelist.bi.proj2$Association_Date[which(is.na(V(igraph.bi$proj2)$size))];
#which(is.na(V(igraph.bi$proj2)$size));  
#which(is.na(V(igraph.bi$proj2)$shape));
#which(is.na(V(igraph.bi$proj2)$color));
#which(is.na(V(igraph.bi$proj1)$size));  
#which(is.na(V(igraph.bi$proj1)$shape));
#which(is.na(V(igraph.bi$proj1)$color));
#edgelist.bi.proj1$Number_of_Residents_in_Certified[which(is.na(V(igraph.bi$proj1)$size))];

colHighlight <- "#450045";
#colHighlight_back <- "#FFF44F";
V(igraph.bi$proj2)$label <- NA;
V(igraph.bi$proj1)$label <- NA;
V(igraph.bi$proj2)$label[which(V(igraph.bi$proj2)$name %in% pri_own)] <- "\u2020";
V(igraph.bi$proj1)$label[which(V(igraph.bi$proj1)$name %in% sel_fac)] <- "\u2020";
#* is "\u002A";
#\u019a is ƚ; 
#\u0197 is Ɨ;
#\u02d0 is ː;
#\u205d is ⁝;
#\u00A6 is ¦;
#\u25CA is ◊;
#\u007C is |;
#\u2020 is a dagger †;

which(is.na(V(igraph.bi$proj2)$shape));
which(is.na(V(igraph.bi$proj2)$size)); #sometimes this is not assigned for all vertices and gives problem in plotting;
which(is.na(V(igraph.bi$proj2)$color));
igraph.bi$proj2.small <- igraph.bi$proj2 - V(igraph.bi$proj2)[which(is.na(V(igraph.bi$proj2)$size))]; 
igraph.bi$proj2.small2 <- igraph.bi$proj2.small - V(igraph.bi$proj2.small)[which(is.na(V(igraph.bi$proj2.small)$color))];
igraph.bi$proj2.small3 <- igraph.bi$proj2.small2 - V(igraph.bi$proj2.small2)[which(is.na(V(igraph.bi$proj2.small2)$shape))];
igraph.bi$proj2 <-igraph.bi$proj2.small3;

#postscript(paste("NH_state_Owners_proj2_SNA", pri_own, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj2_SNA_", "Hypothetical_Owner_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj2_SNA_", "Hypothetical_Owner_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(family="FreeSans", mfrow=c(1,1), mar=c(1,1,1,1), oma=c(1,1,4,1));
#plot(igraph.bi$proj2);
plot(igraph.bi$proj2, vertex.size=V(igraph.bi$proj2)$size/3, vertex.shape=V(igraph.bi$proj2)$shape, vertex.color=V(igraph.bi$proj2)$color, vertex.label=V(igraph.bi$proj2)$label, vertex.label.family="FreeSans", vertex.label.cex=3, vertex.label.color=colHighlight, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend2_text, NA, "Size based on years of ownership", paste("\u2020 Denotes ", length(pri_own), " owner(s) of interest", sep='')), pch=c(legend2_symbols, NA, NA, NA), col="black", pt.bg=c(legend2_colors, NA, NA, colHighlight), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(outer=TRUE, cex=1.8, text=paste("Bipartite projection -- owners affiliated to other owners \nthrough", "facilities accepting Medicare or Medicaid"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj2, "owners", "\nEdges:", edges_proj2, "associations", "\nShortest ownership duration:", smallest_proj2, "years", "\nLongest ownership duration:", largest_proj2, "years", "\nCumulative years of ownership:",sumyearsAD, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
mtext(text=paste("Owner of interest:", "Hypothetical Owner"), cex=1.2, line=-2, side=3, adj=1, outer=TRUE);
dev.off();


#postscript(paste("NH_state_Owners_proj1_SNA_", pri_own, ".eps", sep=''), paper="letter", horizontal=TRUE, fonts=c("serif", "Palatino"));
#png(file=paste("NH_state_Owners_proj1_SNA_", "Hypothetical_Owner_", keyDate, ".png", sep=''), bg="transparent", units="in", width=11, height=8.5, res=300, pointsize=12, type="cairo");
pdf.options(encoding='CP1250');
cairo_pdf(file=paste("NH_state_Owners_proj1_SNA_", "Hypothetical_Owner_", keyDate, ".pdf", sep=''), bg="transparent", width=11, height=8.5, pointsize=12, family="FreeSans");
par(family="FreeSans", mfrow=c(1,1), mar=c(0,4,0,0), oma=c(1,0,3.2,16), xpd=TRUE);
#par(xpd = T, mar = par()$mar + c(0,0,4,7))
#plot(igraph.bi$proj1);
plot(igraph.bi$proj1, vertex.size=V(igraph.bi$proj1)$size/10, vertex.shape=V(igraph.bi$proj1)$shape, vertex.color=V(igraph.bi$proj1)$color, vertex.label=V(igraph.bi$proj1)$label, vertex.label.family="FreeSans", vertex.label.cex=3, vertex.label.color=colHighlight, edge.color="grey", edge.arrow.size=1.5);
#title("Predicted Response\nStandard Deviation", outer=F);
par(lwd=1); legend(xpd=NA, inset=c(-0.47,0), "right", box.lwd=1, legend=c(legend1_text, NA, "Size based on number of residents\n     in certified beds", NA, paste("\u2020 Denotes ", length(sel_fac), " owner of interest\n   NH(s)", sep='')), pch=c(legend1_symbols, NA, NA, NA, NA), col="black", pt.bg=c(legend1_colors, NA, NA, NA, NA), pt.cex=2, cex=1.2, bg="transparent", xjust=0, bty="n");
options(scipen=999);
mtext(xpd=NA, outer=TRUE, cex=1.8, text=paste("Bipartite projection --", "facilities affiliated to other", "\nfacilities accepting Medicare or Medicaid through owners"));
mtext(text="", cex=0.90, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
#mtext(text=paste("Processing date:", processingDate), cex=1.10, line=0, side=SOUTH<-1, adj=1.0, outer=TRUE);
mtext(text=paste("Nodes:", nodes_proj1, "nursing homes", "\nEdges:", edges_proj1, "associations", "\nSmallest facility:", smallest_proj1, "residents", "\nLargest facility:", largest_proj1, "residents", "\nNumber of residents in certified beds:", sumRes_proj1, "\nNumber of certified beds:", sumBeds_proj1, "\nProcessing date:", processingDate), cex=1.2, line=0, side=SOUTH<-1, adj=0.0, outer=TRUE);
mtext(text=paste("Owner of interest:", "Hypothetical Owner"), cex=1.2, line=-2, side=3, adj=1, outer=TRUE);
dev.off();
#-----------------------# FOR "Hypothetical Owner" #-----------------------#

################################ SHAPEFILE MAPPING FOR "Hypothetical Owner" ################################;
sel_pri_own <- unique(edgelist$Federal_Provider_Number);
df1 <- df[which(df$Federal_Provider_Number %in% sel_pri_own),];
#df1$Ownership_Type; #sometimes is missing some for some reason;
#df1$Federal_Provider_Number;
#edgelist.bi.proj1$Federal_Provider_Number;
key <- edgelist.bi.proj1;
key <- key[sort.list(key$Federal_Provider_Number),];
#df1$Federal_Provider_Number;
#key$Federal_Provider_Number;

lat <- df1$lat;
lon <- df1$lon;
tmpDF <- cbind(key, lat, lon); #have to do this to get correct Node attributes (i.e. Ownership_Type)
#edgelist.bi.proj1$Number_of_Residents_in_Certified;

o2 <- as.data.frame(tmpDF);
#o2 <- as.data.frame(df1);
#o2 <- o2[complete.cases(o2),]
coordinates(o2) <- c(which(colnames(o2) %in% "lon"), which(colnames(o2) %in% "lat"))
proj4string(o2) <- CRS(projected);


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
#us_aeaNH@data$Number_of_Certified_Beds
us_aea@data$id = rownames(us_aea@data)

alaska = us_aea[us_aea$STATEFP=="02",]
	#bbox_a1 <- bbox(alaska)
	#center_a1 <- coordinates(alaska)
alaska = elide(alaska, rotate=-50)
	#bbox_a2 <- bbox(alaska)
	#center_a2 <- coordinates(alaska)
alaska = elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska = elide(alaska, shift=c(-2100000, -2500000))
proj4string(alaska) = CRS(projected)

hawaii = us_aea[us_aea$STATEFP=="15",]
hawaii = elide(hawaii, rotate=-35)
hawaii = elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) = CRS(projected)

us_aea = us_aea[!us_aea$STATEFP %in% c("02", "15"),]
us_aea = rbind(us_aea, alaska, hawaii)

us50 <- fortify(us_aea, region="STUSPS")
us50 = remove.territories(us50)
#save('helpers/us50')

p = ggplot(data=us50) + 
    geom_map(map=us50, aes(map_id=id, group=group), ,fill="white", color="dark grey", size=0.15) + 
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
#us_aeaNH@data$Number_of_Certified_Beds;


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

###################################
#plot(us_aea);
#plot(us_aeaNH, add=TRUE);

#plot(us50); #doesn't plot right;
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
#unique(us_aeaNH@data$Ownership_Type);
##----------##;

nodes_US <- length(us_aeaNH@data$Number_of_Certified_Beds);
sumBeds_US <- sum(us_aeaNH@data$Number_of_Certified_Beds, na.rm=TRUE);
sumRes_US <- sum(us_aeaNH@data$Number_of_Residents_in_Certified, na.rm=TRUE);

usRes <- aggregate(us_aeaNH@data$Number_of_Residents_in_Certified, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usBeds <- aggregate(us_aeaNH@data$Number_of_Certified_Beds, by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);
usNHs <- aggregate(rep(1, length(us_aeaNH@data$Ownership_Type)), by=list(Category=us_aeaNH@data$Ownership_Type), FUN=sum);


tmpType1 = factor(usNHs$Category, 
	levels = c("For profit - Corporation", "For profit - Individual", 
	"For profit - Limited Lia", "For profit - Partnership",
	"Government - City", "Government - City/county",
	"Government - County", "Government - Federal",    
	"Government - Hospital di", "Government - State",      
	"Non profit - Church rela", "Non profit - Corporation",
	"Non profit - Other", NA));
legend1_colors <- as.character(key_fac$color[tmpType1]);
legend1_symbols2 <- as.numeric(as.character(key_fac$shape1[tmpType1]));
legend1_sel_text <- legend1_text[tmpType1];

counter <- 1;
while(counter < length(usNHs$Category)+1){
	usNHs$legend[counter] <- paste(legend1_sel_text[counter], ": ", usNHs$x[counter], "\n  ", usRes$x[counter], " residents, ", usBeds$x[counter], " beds", sep="");
	counter <- counter+1;
};
usLegend <- usNHs$legend;
facLegend <- "\u2020 Denotes owner of interest";

 	#theme(plot.title=element_text(face="bold", size=24, hjust=0), legend.key.size=unit(1.6, 'lines')) +

	#guides(colour=guide_legend(override.aes=list(size=12))) +

singledata <- us_aeaNH[which(us_aeaNH$Federal_Provider_Number %in% sel_fac),];

	#scale_colour_manual("Federal Provider Number", values=colHighlight, labels=facLegend) +
 	#scale_shape_manual("Federal Provider Number", values=8, labels=facLegend) +
#p + geom_point(aes(coordinates(o2)[,1], coordinates(o2)[,2], size=Number_of_Residents_in_Certified), data=o2@data, colour=alpha("red",0.2)) +  facet_wrap(~Ownership_Type, ncol=2)
#p1 <- p + geom_point(aes(coordinates(us_aeaNH)[,1], coordinates(us_aeaNH)[,2]), size=1, colour=alpha("black", 0.2), data=us_aeaNH@data);
#legend1_symbols2 <- c(16, 16, 16, 16, 15, 15, 15, 15, 15, 15, 16, 16, 16);
img = readPNG(system.file("img", "Rlogo.png", package="png"));
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
	labs(size="Number of Residents in \nCertified Beds", x=paste("\u2020 Denotes ", length(sel_fac), " owner of interest NH(s)", sep='')) +
	scale_colour_manual("Ownership Type", values=legend1_colors, labels=usLegend) +
 	scale_shape_manual("Ownership Type", values=legend1_symbols2, labels=usLegend) +
 	##guides(colour=guide_legend(override.aes=list(size=12))) +
 	ggtitle(paste("Nursing Homes Accepting Medicare or Medicaid \nAffiliated by Ownership Association to\n", "Hypothetical Owner", sep='')) +
	theme(plot.title=element_text(family="FreeSans", face="bold", size=24, hjust=0), plot.subtitle=element_text(family="FreeSans", size=12), axis.title=element_text(family="FreeSans", size=12)) +
 	guides(size=guide_legend(keyheight=unit(1.0, 'lines'), override.aes = list(color="darkgrey"))) +
 	guides(colour=guide_legend(keyheight=unit(1.7, 'lines'), override.aes = list(size=4))) +
 	theme(legend.text=element_text(family="FreeSans", size=12), legend.title=element_text(family="FreeSans", face="italic", size=13)) +
 	geom_point(data=singledata@data, mapping=aes(x=coordinates(singledata)[,1], y=coordinates(singledata)[,2]), color=colHighlight, shape="†", size=6) + 
 	labs(subtitle=paste("Number of Nursing Homes:", nodes_US, "\nNumber of Residents in Certified Beds:", sumRes_US, "\nNumber of Certified Beds:", sumBeds_US, "\nProcessing Date:", processingDate), caption="") +
 	expand_limits(x=us50$long, y=us50$lat);
###########################;
#theme_get(); #check this to see theme options for ggplot with element_text() object!;

#p1;
#ggsave(paste("Nursing_home_locations_US", "Hypothetical_Owner_", keyDate, ".png", sep = ""), p1, width = 10, height = 8);
ggsave(paste("Nursing_home_locations_US", "Hypothetical_Owner_", keyDate, ".pdf", sep = ""), p1, width = 10, height = 8, device=cairo_pdf);
##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##-~-##

print(k);
k <- k + 1;

}
mclapply(1:counterPeriod, periodLoop);
#mclapply(1, periodLoop);
####################### Looping ends here ###############################;

####################### All Looping ends here ###############################;

#stopCluster(cl1);