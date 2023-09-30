/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;
*Medicare Nursing Home Open-Data combine annual cycles;
*April 2018, Tyler Pittman;

*************************************************************;
******DON'T EXECUTE THIS SCRIPT, GOT DATA LINKED RIGHT*******;
*************************************************************;

proc datasets library=work kill memtype=data;

libname mylib "/folders/myshortcuts/tylerpittman/NursingHomesSNA/data";
%let dir = /folders/myshortcuts/tylerpittman/NursingHomesSNA/data;

data WORK.NHGEO    ;
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "&dir./nursing_homes_geocoordinates_COMPLETE.csv" delimiter = ',' MISSOVER DSD  firstobs=2 ;
informat Federal_Provider_Number $35. ;
format Federal_Provider_Number $35. ;
input
Federal_Provider_Number $
lat
lon
;
if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;

data period;
length a $19.;
input a$ ;
datalines;   
March2016 
September2016 
March2017 
September2017 
March2018 
;
run;

/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/
proc sql noprint;
   select count(*) into : nobs
   from period;
quit;
*%put &nobs;

%MACRO Get_data(myDataset=,myLine=,myColumn=,myMVar=);
%GLOBAL &myMVar.;
data _null_;
set &myDataset.;
if _N_ = &myLine. then do;
call symput(symget('myMVar'),&myColumn.);
end;
run;
%MEND Get_data;

/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/

/*---------- TEST ----------*/;
%macro Assemble_data;
	%let count = 1;
	%do count=1 %to &nobs;

%put &count;


%Get_data(myDataset=period, myLine=&count, myColumn=a, myMVar=t);
%let t1 = %qsysfunc(compress(&t,%str(%")));
%let time =%superq(t1);
%put &time;

/*
%let time = March2016;
%let count = 1;
*/;

PROC SQL;
CREATE TABLE WORK.MDST_&count AS
SELECT PERIOD , Federal_Provider_Number , Provider_Name , Provider_Address , Provider_City , Provider_State , 
Provider_Zip_Code , Provider_Phone_Number , Provider_SSA_County_Code , Provider_County_Name , Ownership_Type , 
Number_of_Certified_Beds , Number_of_Residents_in_Certified , Provider_Type , Provider_Resides_in_Hospital , Legal_Business_Name ,
Date_First_Approved_to_Provide_M , Continuing_Care_Retirement_Commu , Special_Focus_Facility , Most_Recent_Health_Inspection_Mo , 
Provider_Changed_Ownership_in_La , With_a_Resident_and_Family_Counc , Automatic_Sprinkler_Systems_in_A , Overall_Rating , 
Overall_Rating_Footnote , Health_Inspection_Rating , Health_Inspection_Rating_Footnot , QM_Rating , QM_Rating_Footnote , Staffing_Rating , 
Staffing_Rating_Footnote , RN_Staffing_Rating , RN_Staffing_Rating_Footnote , Reported_Staffing_Footnote , Physical_Therapist_Staffing_Foot , 
Reported_CNA_Staffing_Hours_per , Reported_LPN_Staffing_Hours_per , Reported_RN_Staffing_Hours_per_R , Reported_Licensed_Staffing_Hours , 
Reported_Total_Nurse_Staffing_Ho , Reported_Physical_Therapist_Staf , Expected_CNA_Staffing_Hours_per , Expected_LPN_Staffing_Hours_per , 
Expected_RN_Staffing_Hours_per_R , Expected_Total_Nurse_Staffing_Ho , Adjusted_CNA_Staffing_Hours_per , Adjusted_LPN_Staffing_Hours_per , 
Adjusted_RN_Staffing_Hours_per_R , Adjusted_Total_Nurse_Staffing_Ho , Total_Weighted_Health_Survey_Sco , Number_of_Facility_Reported_Inci , 
Number_of_Substantiated_Complain , Number_of_Fines , Total_Amount_of_Fines_in_Dollars , Number_of_Payment_Denials , Total_Number_of_Penalties , 
'Location'n , Processing_Date , Total_Fine_Amount , lat , lon FROM MYLIB.MDST_&count;
RUN;
QUIT;

/*
data MYLIB.MDST; set MYLIB.MDST;
PERIOD = "March2017";
run;
*/;

PROC SQL;
CREATE TABLE WORK.MDST_3 AS
SELECT PERIOD , Federal_Provider_Number , Provider_Name , Provider_Address , Provider_City , Provider_State , 
Provider_Zip_Code , Provider_Phone_Number , Provider_SSA_County_Code , Provider_County_Name , Ownership_Type , 
Number_of_Certified_Beds , Number_of_Residents_in_Certified , Provider_Type , Provider_Resides_in_Hospital , Legal_Business_Name ,
Date_First_Approved_to_Provide_M , Continuing_Care_Retirement_Commu , Special_Focus_Facility , Most_Recent_Health_Inspection_Mo , 
Provider_Changed_Ownership_in_La , With_a_Resident_and_Family_Counc , Automatic_Sprinkler_Systems_in_A , Overall_Rating , 
Overall_Rating_Footnote , Health_Inspection_Rating , Health_Inspection_Rating_Footnot , QM_Rating , QM_Rating_Footnote , Staffing_Rating , 
Staffing_Rating_Footnote , RN_Staffing_Rating , RN_Staffing_Rating_Footnote , Reported_Staffing_Footnote , Physical_Therapist_Staffing_Foot , 
Reported_CNA_Staffing_Hours_per , Reported_LPN_Staffing_Hours_per , Reported_RN_Staffing_Hours_per_R , Reported_Licensed_Staffing_Hours , 
Reported_Total_Nurse_Staffing_Ho , Reported_Physical_Therapist_Staf , Expected_CNA_Staffing_Hours_per , Expected_LPN_Staffing_Hours_per , 
Expected_RN_Staffing_Hours_per_R , Expected_Total_Nurse_Staffing_Ho , Adjusted_CNA_Staffing_Hours_per , Adjusted_LPN_Staffing_Hours_per , 
Adjusted_RN_Staffing_Hours_per_R , Adjusted_Total_Nurse_Staffing_Ho , Total_Weighted_Health_Survey_Sco , Number_of_Facility_Reported_Inci , 
Number_of_Substantiated_Complain , Number_of_Fines , Total_Amount_of_Fines_in_Dollars , Number_of_Payment_Denials , Total_Number_of_Penalties , 
'Location'n , Processing_Date , Total_Fine_Amount , lat , lon FROM MYLIB.MDST;
RUN;
QUIT;

/*
proc freq data=MYLIB.MDST_3;
table Ownership_type; *missing 1866;
run;

proc freq data=MDST_3;
table Ownership_type; *missing none;
run;

proc freq data=MYLIB.MDST;
table Ownership_type; *missing none;
run;

proc freq data=MYLIB.MDST_1;
table Ownership_type; *missing 15;
run;

proc freq data=MYLIB.MDST_2;
table Ownership_type; *missing none;
run;

proc freq data=MYLIB.MDST_4;
table Ownership_type; *missing none;
run;

proc freq data=MYLIB.MDST_5;
table Ownership_type; *missing 15;
run;
*/;

/*
PROC DATASETS NOLIST NODETAILS;
CONTENTS DATA=WORK.MDST_&count OUT=WORK.details;
RUN;

PROC PRINT DATA=WORK.details;
RUN;
*/;

data WORK.MDST_&count; set WORK.MDST_&count;
if missing(Ownership_Type) then delete;
run;

data checkEnd; set WORK.MDST_&count;
where Federal_Provider_Number = "1.40E+148" ;
run;

data MDST_&count; set MDST_&count;
if find(Federal_Provider_Number,"+") then delete;
run;

data checkEnd; set MDST_&count;
where Federal_Provider_Number = "1.40E+148" ;
run;

data checkEnd; set MDST_&count;
where Federal_Provider_Number = "15009" ; *must make this 015009 to merge with nhGeo;
run;

data MDST_&count; set MDST_&count;
Federal_Provider_Number2 = REVERSE(SUBSTR(TRIM(LEFT(REVERSE(Federal_Provider_Number))) || "000000", 1, 6));
run;

data MDST_&count; set MDST_&count;
Federal_Provider_Number = Federal_Provider_Number2;
drop Federal_Provider_Number2;
run;

data checkEnd; set MDST_&count;
where Federal_Provider_Number = "015009" ; *must make this 015009 to merge with nhGeo;
run;

data MDST_&count; set MDST_&count;
if missing(Federal_Provider_Number) then delete;
drop lat lon;
run;

data nhGeo; set nhGeo;
if missing(Federal_Provider_Number) then delete;
run;

proc sql;
create table nhs_&count as
select A.*, B.lat, B.lon
from MDST_&count as A left join nhGeo as B on (A.Federal_Provider_Number=B.Federal_Provider_Number);
quit;


data mylib.nhs_&count; set nhs_&count;
run;

proc export data=nhs_&count
    outfile="&dir./cms_nursing_homes_&time..csv"
    dbms=csv
    replace;
run;


proc datasets lib=work;
delete 
mdsll_&count mdst_&count;
run;
quit;



	%end;
%exit: %mend Assemble_data;
%Assemble_data
/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/


data nhs_combined;
set nhs_1 nhs_2 nhs_3 nhs_4 nhs_5;
run;


proc freq data=nhs_combined; 
tables PERIOD;
run;


data mylib.nhs_combined; set nhs_combined;
run;

proc export data=nhs_combined
    outfile="&dir./nhs_combined.csv"
    dbms=csv
    replace;
run;
*/;
/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/

/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;