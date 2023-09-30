/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;
*Medicare Nursing Home Open-Data combine annual cycles;
*April 2018, Tyler Pittman;

proc datasets library=work kill memtype=data;

libname mylib "/folders/myshortcuts/tylerpittman/NursingHomesSNA/data";
%let dir = /folders/myshortcuts/tylerpittman/NursingHomesSNA/data;


/****************************************************************************************/;
/*
Next run Rscript --no-save /Users/tylerpittman/GitHub/NursingHomesSNA/scripts/STEP_01_nursing_homes_geocode_looping_April2018.R;
*/;

/*--------------------------*/;
* Written by R;
*  write.foreign(df, "nhsll_combined.txt", "nhsll_combined.sas",  ;

DATA nhsLL;
LENGTH
 PERIOD $ 13
 Federal_Provider_Number $ 6
 Provider_Name $ 48
 Provider_Address $ 42
 Provider_City $ 13
 Provider_State $ 2
 Provider_Zip_Code $ 5
 Provider_County_Name $ 22
 Ownership_Type $ 24
 Provider_Type $ 21
 Provider_Resides_in_Hospital $ 5
 Legal_Business_Name $ 44
 Continuing_Care_Retirement_Commu $ 5
 Special_Focus_Facility $ 5
 Most_Recent_Health_Inspection_Mo $ 5
 Provider_Changed_Ownership_in_La $ 5
 With_a_Resident_and_Family_Counc $ 8
 Automatic_Sprinkler_Systems_in_A $ 3
 Overall_Rating_Footnote $ 1
 Health_Inspection_Rating_Footnot $ 1
 QM_Rating_Footnote $ 1
 Staffing_Rating_Footnote $ 1
 RN_Staffing_Rating_Footnote $ 1
 Reported_Staffing_Footnote $ 1
 Physical_Therapist_Staffing_Foot $ 1
 Location $ 58
;

INFORMAT
 Date_First_Approved_to_Provide_M
 Processing_Date
 YYMMDD10.
;

INFILE  "&dir./nhsll_combined.txt" 
     DSD 
     LRECL= 600 ;
INPUT
 PERIOD
 Federal_Provider_Number
 Provider_Name
 Provider_Address
 Provider_City
 Provider_State
 Provider_Zip_Code
 Provider_Phone_Number
 Provider_SSA_County_Code
 Provider_County_Name
 Ownership_Type
 Number_of_Certified_Beds
 Number_of_Residents_in_Certified
 Provider_Type
 Provider_Resides_in_Hospital
 Legal_Business_Name
 Date_First_Approved_to_Provide_M
 Continuing_Care_Retirement_Commu
 Special_Focus_Facility
 Most_Recent_Health_Inspection_Mo
 Provider_Changed_Ownership_in_La
 With_a_Resident_and_Family_Counc
 Automatic_Sprinkler_Systems_in_A
 Overall_Rating
 Overall_Rating_Footnote
 Health_Inspection_Rating
 Health_Inspection_Rating_Footnot
 QM_Rating
 QM_Rating_Footnote
 Staffing_Rating
 Staffing_Rating_Footnote
 RN_Staffing_Rating
 RN_Staffing_Rating_Footnote
 Reported_Staffing_Footnote
 Physical_Therapist_Staffing_Foot
 Reported_CNA_Staffing_Hours_per
 Reported_LPN_Staffing_Hours_per
 Reported_RN_Staffing_Hours_per_R
 Reported_Licensed_Staffing_Hours
 Reported_Total_Nurse_Staffing_Ho
 Reported_Physical_Therapist_Staf
 Expected_CNA_Staffing_Hours_per
 Expected_LPN_Staffing_Hours_per
 Expected_RN_Staffing_Hours_per_R
 Expected_Total_Nurse_Staffing_Ho
 Adjusted_CNA_Staffing_Hours_per
 Adjusted_LPN_Staffing_Hours_per
 Adjusted_RN_Staffing_Hours_per_R
 Adjusted_Total_Nurse_Staffing_Ho
 Total_Weighted_Health_Survey_Sco
 Number_of_Facility_Reported_Inci
 Number_of_Substantiated_Complain
 Number_of_Fines
 Total_Amount_of_Fines_in_Dollars
 Number_of_Payment_Denials
 Total_Number_of_Penalties
 Location
 Processing_Date
 Total_Fine_Amount
 lat
 lon
;
FORMAT Date_First_Approved_to_Provide_M yymmdd10.;
FORMAT Processing_Date yymmdd10.;
RUN;
/****************************************************************************************/;

/*
proc freq data=nhsLL;
tables lat;
*tables Measure_Description lat lon;
run; *most lats have frequency 5 for 5 timepoints;
*/;

data mylib.nhsLL; set work.nhsLL;
run;




/*--------------------*/;
/*
%web_drop_table(MYLIB.nhACScensusCT);
FILENAME REFFILE "&dir./Medicare_nursinghomes_ACS2015_Census2010countyCTs_combined.csv";
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=MYLIB.nhACScensusCT;
	GETNAMES=YES;
RUN;
PROC CONTENTS DATA=MYLIB.nhACScensusCT; RUN;
%web_open_table(MYLIB.nhACScensusCT);

data nh; set MYLIB.nhACScensusCT; *only 9765 varibles as opposed to 15625, but has 15664 obs;
run; *15408 nursing homes, missing 256 as treats Federal_Provider_Number as int as opposed to character in read-in;
*/;

data nhsLL; set mylib.nhsLL; *61 variables;
run;

/*
data nhGeo; set nhsLL;
keep Federal_Provider_Number lat lon;
run;

proc export data=nhGeo
    outfile="&dir./nursing_homes_geocoordinates_COMPLETE.csv"
    dbms=csv
    replace;
run;
*/;

proc freq data=nhsLL;
tables Federal_Provider_Number Ownership_Type Number_of_Certified_Beds Number_of_Residents_in_Certified Provider_Type Provider_Resides_in_Hospital 
Total_Fine_Amount Provider_Name Legal_Business_Name; 
run;
*77760 rows, 0 missing for all except the below;
*1 missing Number_of_Certified_Beds;
*58092 missing Total_Fine_Amount;



proc univariate data=nhsLL normal;
var Total_Fine_Amount;
class Ownership_Type;
run;

/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;