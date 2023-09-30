/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;
*Medicare Nursing Home Open-Data combine annual cycles;
*Using SAS University Edition;
*April 2018, Tyler Pittman;

**************************************************************************************************;
******DON'T EXECUTE THIS SCRIPT, UNIVERSITY SAS RUNS OUT OF SPACE AND GOT DATA LINKED RIGHT*******;
**************************************************************************************************;

proc datasets library=work kill memtype=data;

libname mylib "/folders/myshortcuts/tylerpittman/NursingHomesSNA/data";
%let dir = /folders/myshortcuts/tylerpittman/NursingHomesSNA/data;

/*
libname xptfile xport "&dir./MMSA2015.xpt" access=readonly;
*/;

/*
proc copy inlib=xptfile outlib=mylib;
run;
libname xptfile2 xport "&dir./CNTY10.xpt" access=readonly;
proc copy inlib=xptfile2 outlib=mylib;
run;
*/;

/*
data brfss2015; set mylib.TRNSPORT;
run;

proc freq data=brfss2015; *BRFSS2015 data;
tables _CNTYNAM;
run;
*/;

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


/* Use this to clean csv files carriage returns from Windows to Unix */;
/*
filename old "&dir./Deficiencies_September2018.csv";
filename new "&dir./Deficiencies_September2018_fixed.csv";

filename old "&dir./Penalties_September2018.csv";
filename new "&dir./Penalties_September2018_fixed.csv";

filename old "&dir./ProviderInfo_September2018.csv";
filename new "&dir./ProviderInfo_September2018_fixed.csv";

filename old "&dir./Ownership_September2018.csv";
filename new "&dir./Ownership_September2018_fixed.csv";

filename old "&dir./QualityMsrMDS_September2018.csv";
filename new "&dir./QualityMsrMDS_September2018_fixed.csv";

*/;
/*
filename old "&dir./QualityMsrMDS_September2018.csv";
filename new "&dir./QualityMsrMDS_September2018_fixed.csv";

data _null_ ;
  if eof then put 'NOTE: Records read=' newn 'Records with missing quotes=' missq ;
  infile old lrecl=10000 end=eof ;
  file new lrecl=10000;
  nq=0;
  do until (mod(nq,2)=0 or eof );
     input;
     newn+1;
     nq = nq + countc(_infile_,'"');
     put _infile_ @;
     if mod(nq,2) then do;
       missq+1;
       put ' ' @;
     end;
  end;
  put;
run;
*/;

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
%let time = March2018;
%let count = 5;
%let time = March2016;
%let count = 1;
*/;

/*
%let time = March2018;
%put &time;
data tmp; 
length PERIOD $20;
set NHGEO;
PERIOD = "&time";
run;

	%end;
%exit: %mend Assemble_data;
%Assemble_data
*/;
/*---------- TEST ----------*/;


data WORK.OWN_&count    ;
%let _EFIERR_ = 0; 
infile "&dir./Ownership_&time._fixed.csv" delimiter = ',' MISSOVER DSD  firstobs=2 ;
informat Federal_Provider_Number $35. ;
informat Provider_Name $47. ;
informat Provider_Address $50. ;
informat Provider_City $50. ;
informat Provider_State $2. ;
informat Provider_Zip_Code $10. ;
informat Role_played_by_Owner_or_Manager $41. ;
informat Owner_Type $12. ;
informat Owner_Name $38. ;
informat Ownership_Percentage $22. ;
informat Association_Date $16. ;
informat Location $22. ;
informat Processing_Date mmddyy10. ;
format Federal_Provider_Number $35. ;
format Provider_Name $47. ;
format Provider_Address $50. ;
format Provider_City $50. ;
format Provider_State $2. ;
format Provider_Zip_Code $10. ;
format Role_played_by_Owner_or_Manager $41. ;
format Owner_Type $12. ;
format Owner_Name $38. ;
format Ownership_Percentage $22. ;
format Association_Date $16. ;
format Location $22. ;
format Processing_Date mmddyy10. ;
input
Federal_Provider_Number $
Provider_Name $
Provider_Address $
Provider_City $
Provider_State $
Provider_Zip_Code $
Role_played_by_Owner_or_Manager $
Owner_Type $
Owner_Name $
Ownership_Percentage $
Association_Date $
Location $
Processing_Date
;
if _ERROR_ then call symputx('_EFIERR_',1);  
run;


/*
proc import datafile="&dir./Ownership_&time._fixed.xls"
   out=WORK.OWN_&count
   dbms=xls replace;
   getnames=yes;
   sheet="Ownership_&time._fixed";
run;

data OWN_&count; set OWN_&count;
Federal_Provider_Number2 = input(strip(Federal_Provider_Number),$35.);
format Federal_Provider_Number2 $35.;
drop Federal_Provider_Number;
run;

data OWN_&count; set OWN_&count;
Federal_Provider_Number = Federal_Provider_Number2;
drop Federal_Provider_Number2;
run;
*/;

/*
data OWNb_&count;
retain Federal_Provider_Number Provider_Name Provider_Address Provider_City Provider_State Provider_Zip_Code 
Role_played_by_Owner_or_Manager Owner_Type Owner_Name Ownership_Percentage Association_Date 'Location'n 
Processing_Date  ;
set OWN_&count;
run;

proc export data=OWNb_&count
    outfile="&dir./Ownership_&time._fixed.csv"
    dbms=csv
    replace;
run;
*/;

/*
proc import datafile="&dir./Ownership_&time._fixed.csv"
     out=WORK.OWN_&count 
     dbms=csv
     replace;
     getnames=yes;
run;
*/;

/*
proc freq data=OWN_&count; 
tables Owner_Name Federal_Provider_Number;
run;
*/;

data WORK.PEN_&count    ;
%let _EFIERR_ = 0; 
infile "&dir./Penalties_&time._fixed.csv" delimiter = ',' MISSOVER DSD  firstobs=2 ;
informat Federal_Provider_Number $35. ;
informat Provider_Name $47. ;
informat Provider_Address $28. ;
informat Provider_City $11. ;
informat Provider_State $2. ;
informat Provider_Zip_Code best32. ;
informat Penalty_Date mmddyy10. ;
informat Penalty_Type $14. ;
informat Fine_Amount dollar10. ;
informat Payment_Denial_Start_Date mmddyy10. ;
informat Payment_Denial_Length_in_Days best32. ;
informat Location $20. ;
format Federal_Provider_Number $35. ;
format Provider_Name $47. ;
format Provider_Address $28. ;
format Provider_City $11. ;
format Provider_State $2. ;
format Provider_Zip_Code best12. ;
format Penalty_Date mmddyy10. ;
format Penalty_Type $14. ;
format Fine_Amount dollar10. ;
format Payment_Denial_Start_Date mmddyy10. ;
format Payment_Denial_Length_in_Days best12. ;
format Location $20. ;
input
Federal_Provider_Number $
Provider_Name $
Provider_Address $
Provider_City $
Provider_State $
Provider_Zip_Code
Penalty_Date
Penalty_Type $
Fine_Amount 
Payment_Denial_Start_Date
Payment_Denial_Length_in_Days
Location $
;
if _ERROR_ then call symputx('_EFIERR_',1);  
run;

data WORK.MDS_&count ;
%let _EFIERR_ = 0; 
infile "&dir./QualityMsrMDS_&time._fixed.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
informat Federal_Provider_Number $35. ;
informat Provider_Name $49. ;
informat Provider_Address $21. ;
informat Provider_City $16. ;
informat Provider_State $2. ;
informat Provider_Zip_Code $10. ;
informat Measure_Code best32. ;
informat Measure_Description $140. ;
informat Resident_type $10. ;
informat Q1_Measure_Score $10. ;
informat Footnote_for_Q1_Measure_Score $1. ;
informat Q2_Measure_Score $9. ;
informat Footnote_for_Q2_Measure_Score $1. ;
informat Q3_Measure_Score $10. ;
informat Footnote_for_Q3_Measure_Score $1. ;
informat Q4_Measure_Score $10. ;
informat Footnote_for_Q4_Measure_Score $1. ;
informat Four_Quarter_Average_Score $10. ;
informat Footnote_for_Four_Quarter_Averag $1. ;
informat Used_in_Quality_Measure_Five_Sta $5. ;
informat Q1_quarter anydtdtm40. ;
informat Q2_quarter anydtdtm40. ;
informat Q3_quarter anydtdtm40. ;
informat Q4_quarter anydtdtm40. ;
informat Location $140. ;
informat Processing_Date anydtdtm40. ;
format Federal_Provider_Number $35. ;
format Provider_Name $49. ;
format Provider_Address $21. ;
format Provider_City $16. ;
format Provider_State $2. ;
format Provider_Zip_Code $10. ;
format Measure_Code best12. ;
format Measure_Description $140. ;
format Resident_type $10. ;
format Q1_Measure_Score $10. ;
format Footnote_for_Q1_Measure_Score $1. ;
format Q2_Measure_Score $9. ;
format Footnote_for_Q2_Measure_Score $1. ;
format Q3_Measure_Score $10. ;
format Footnote_for_Q3_Measure_Score $1. ;
format Q4_Measure_Score $10. ;
format Footnote_for_Q4_Measure_Score $1. ;
format Four_Quarter_Average_Score $10. ;
format Footnote_for_Four_Quarter_Averag $1. ;
format Used_in_Quality_Measure_Five_Sta $5. ;
format Q1_quarter datetime. ;
format Q2_quarter datetime. ;
format Q3_quarter datetime. ;
format Q4_quarter datetime. ;
format Location $140. ;
format Processing_Date datetime. ;
input
Federal_Provider_Number $
Provider_Name $
Provider_Address $
Provider_City $
Provider_State $
Provider_Zip_Code $
Measure_Code
Measure_Description $
Resident_type $
Q1_Measure_Score $
Footnote_for_Q1_Measure_Score $
Q2_Measure_Score $
Footnote_for_Q2_Measure_Score $
Q3_Measure_Score $
Footnote_for_Q3_Measure_Score $
Q4_Measure_Score $
Footnote_for_Q4_Measure_Score $
Four_Quarter_Average_Score $
Footnote_for_Four_Quarter_Averag $
Used_in_Quality_Measure_Five_Sta $
Q1_quarter
Q2_quarter
Q3_quarter
Q4_quarter
Location $
Processing_Date
;
if _ERROR_ then call symputx('_EFIERR_',1);  
run;

/******* Put Federal_Provider_Number into regular number and non scientific notation for correct merge ***********/;
/*
data checkEnd; set mds_1;
where Federal_Provider_Number = "1.40E+148" ;
run;

data mds_tmp; set mds_1;
if find(Federal_Provider_Number,"+") then delete;
run;

data checkEnd; set mds_tmp;
where Federal_Provider_Number = "1.40E+148" ;
run;

proc sql;
create table mds_tmp1 as
select A.*, C.Federal_Provider_Number
from own_1 as C Left Join mds_tmp as A on (C.Federal_Provider_Number = A.Federal_Provider_Number)
;
quit;

data mds_tmp1; set mds_tmp1;
if missing(Federal_Provider_Number) then delete;
run;

proc export data=mds_tmp1
    outfile="&dir./QualityMsrMDS_March2016_fixed.csv"
    dbms=csv
    replace;
run;
*/;

/*
data mds_tmp1; set mds_&count;
format name $47.;
zip = Provider_Zip_Code*1;
name = Provider_Name;
drop Federal_Provider_Number;
run;

data pen_tmp1; set pen_&count;
format Address $21.;
Address = Provider_Address;
run;

proc sql;
create table mds_tmp2 as
select A.*, C.Federal_Provider_Number, C.Provider_Name, C.Provider_Address, C.Provider_City, C.Provider_State, 
C.Provider_Zip_Code
from mds_tmp1 as A Left Join pen_tmp1 as C on (A.Provider_Address=C.Address and A.zip=C.Provider_Zip_Code)
;
quit;

data mds_tmp2; 
retain Federal_Provider_Number Provider_Name Provider_Address Provider_City Provider_State Provider_Zip_Code Measure_Code Measure_Description Resident_type 
Q1_Measure_Score Footnote_for_Q1_Measure_Score Q2_Measure_Score Footnote_for_Q2_Measure_Score Q3_Measure_Score Footnote_for_Q3_Measure_Score 
Q4_Measure_Score Footnote_for_Q4_Measure_Score Four_Quarter_Average_Score Footnote_for_Four_Quarter_Averag Used_in_Quality_Measure_Five_Sta 
Q1_quarter Q2_quarter Q3_quarter Q4_quarter 'Location'n Processing_Date;
set mds_tmp2;
if missing (Federal_Provider_Number) then delete;
drop Address name zip;
run;

data checkEnd; set mds_tmp2;
where Federal_Provider_Number = "1.40E+148" ;
run;

proc export data=mds_tmp2
    outfile="&dir./QualityMsrMDS_March2016_fixed.csv"
    dbms=csv
    replace;
run;
*/;
/****************************************************************************************************************/;


data mds_&count; set mds_&count;
latitude = SCAN(SUBSTR(Location,INDEX(Location,'(')+1),1,',');
longitude = scan(substr(Location, 1, index(Location, ')') - 1),-1,',');
lat = latitude * 1;
lon = longitude * 1;
drop latitude longitude;
run;

data mds_&count; set mds_&count; 
  len=length(lat); 
   do i=1 to len;                                                                                                                       
     str=substr(lat,i,1);     
     if anyalpha(str) = 1 then do;
        lat = .;
        lon = .;
     end; 
   end;                                                                                                                               
run; 

proc sort data=MDS_&count;
by Federal_Provider_Number; 
run;

proc sort data=nhGeo;
by Federal_Provider_Number; 
run;

/*
DATA mdsll_&count; *does one-to-many merge;
MERGE mds_&count nhGeo; 
BY Federal_Provider_Number; 
RUN;
*/;

data mds_&count; set mds_&count;
if missing(Federal_Provider_Number) then delete;
drop lat lon;
run;

data nhGeo; set nhGeo;
if missing(Federal_Provider_Number) then delete;
run;

proc sql;
create table mdsll_&count as
select A.*, B.lat, B.lon
from mds_&count as A left join nhGeo as B on (A.Federal_Provider_Number=B.Federal_Provider_Number);
quit;

data MDS_&count; set mdsll_&count;
run;


data WORK.PI_&count;
%let _EFIERR_ = 0; 
infile "&dir./ProviderInfo_&time._fixed.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
informat Federal_Provider_Number $35. ;
informat Provider_Name $48. ;
informat Provider_Address $42. ;
informat Provider_City $13. ;
informat Provider_State $2. ;
informat Provider_Zip_Code $10. ;
informat Provider_Phone_Number best32. ;
informat Provider_SSA_County_Code best32. ;
informat Provider_County_Name $22. ;
informat Ownership_Type $24. ;
informat Number_of_Certified_Beds best32. ;
informat Number_of_Residents_in_Certified best32. ;
informat Provider_Type $21. ;
informat Provider_Resides_in_Hospital $5. ;
informat Legal_Business_Name $44. ;
informat Date_First_Approved_to_Provide_M mmddyy10. ;
informat Continuing_Care_Retirement_Commu $5. ;
informat Special_Focus_Facility $5. ;
informat Most_Recent_Health_Inspection_Mo $5. ;
informat Provider_Changed_Ownership_in_La $5. ;
informat With_a_Resident_and_Family_Counc $8. ;
informat Automatic_Sprinkler_Systems_in_A $3. ;
informat Overall_Rating best32. ;
informat Overall_Rating_Footnote $1. ;
informat Health_Inspection_Rating best32. ;
informat Health_Inspection_Rating_Footnot $1. ;
informat QM_Rating best32. ;
informat QM_Rating_Footnote $1. ;
informat Staffing_Rating best32. ;
informat Staffing_Rating_Footnote $1. ;
informat RN_Staffing_Rating best32. ;
informat RN_Staffing_Rating_Footnote $1. ;
informat Reported_Staffing_Footnote $1. ;
informat Physical_Therapist_Staffing_Foot $1. ;
informat Reported_CNA_Staffing_Hours_per best32. ;
informat Reported_LPN_Staffing_Hours_per best32. ;
informat Reported_RN_Staffing_Hours_per_R best32. ;
informat Reported_Licensed_Staffing_Hours best32. ;
informat Reported_Total_Nurse_Staffing_Ho best32. ;
informat Reported_Physical_Therapist_Staf best32. ;
informat Expected_CNA_Staffing_Hours_per best32. ;
informat Expected_LPN_Staffing_Hours_per best32. ;
informat Expected_RN_Staffing_Hours_per_R best32. ;
informat Expected_Total_Nurse_Staffing_Ho best32. ;
informat Adjusted_CNA_Staffing_Hours_per best32. ;
informat Adjusted_LPN_Staffing_Hours_per best32. ;
informat Adjusted_RN_Staffing_Hours_per_R best32. ;
informat Adjusted_Total_Nurse_Staffing_Ho best32. ;
informat Cycle_1_Total_Number_of_Health_D best32. ;
informat Cycle_1_Number_of_Standard_Healt best32. ;
informat Cycle_1_Number_of_Complaint_Heal best32. ;
informat Cycle_1_Health_Deficiency_Score best32. ;
informat Cycle_1_Standard_Survey_Health_D mmddyy10. ;
informat Cycle_1_Number_of_Health_Revisit best32. ;
informat Cycle_1_Health_Revisit_Score best32. ;
informat Cycle_1_Total_Health_Score best32. ;
informat Cycle_2_Total_Number_of_Health_D best32. ;
informat Cycle_2_Number_of_Standard_Healt best32. ;
informat Cycle_2_Number_of_Complaint_Heal best32. ;
informat Cycle_2_Health_Deficiency_Score best32. ;
informat Cycle_2_Standard_Health_Survey_D mmddyy10. ;
informat Cycle_2_Number_of_Health_Revisit best32. ;
informat Cycle_2_Health_Revisit_Score best32. ;
informat Cycle_2_Total_Health_Score best32. ;
informat Cycle_3_Total_Number_of_Health_D best32. ;
informat Cycle_3_Number_of_Standard_Healt best32. ;
informat Cycle_3_Number_of_Complaint_Heal best32. ;
informat Cycle_3_Health_Deficiency_Score best32. ;
informat Cycle_3_Standard_Health_Survey_D mmddyy10. ;
informat Cycle_3_Number_of_Health_Revisit best32. ;
informat Cycle_3_Health_Revisit_Score best32. ;
informat Cycle_3_Total_Health_Score best32. ;
informat Total_Weighted_Health_Survey_Sco best32. ;
informat Number_of_Facility_Reported_Inci best32. ;
informat Number_of_Substantiated_Complain best32. ;
informat Number_of_Fines best32. ;
informat Total_Amount_of_Fines_in_Dollars nlnum32. ;
informat Number_of_Payment_Denials best32. ;
informat Total_Number_of_Penalties best32. ;
informat Location $58. ;
informat Processing_Date mmddyy10. ;
format Federal_Provider_Number $35. ;
format Provider_Name $48. ;
format Provider_Address $42. ;
format Provider_City $13. ;
format Provider_State $2. ;
format Provider_Zip_Code $10. ;
format Provider_Phone_Number best12. ;
format Provider_SSA_County_Code best12. ;
format Provider_County_Name $22. ;
format Ownership_Type $24. ;
format Number_of_Certified_Beds best12. ;
format Number_of_Residents_in_Certified best12. ;
format Provider_Type $21. ;
format Provider_Resides_in_Hospital $5. ;
format Legal_Business_Name $44. ;
format Date_First_Approved_to_Provide_M mmddyy10. ;
format Continuing_Care_Retirement_Commu $5. ;
format Special_Focus_Facility $5. ;
format Most_Recent_Health_Inspection_Mo $5. ;
format Provider_Changed_Ownership_in_La $5. ;
format With_a_Resident_and_Family_Counc $8. ;
format Automatic_Sprinkler_Systems_in_A $3. ;
format Overall_Rating best12. ;
format Overall_Rating_Footnote $1. ;
format Health_Inspection_Rating best12. ;
format Health_Inspection_Rating_Footnot $1. ;
format QM_Rating best12. ;
format QM_Rating_Footnote $1. ;
format Staffing_Rating best12. ;
format Staffing_Rating_Footnote $1. ;
format RN_Staffing_Rating best12. ;
format RN_Staffing_Rating_Footnote $1. ;
format Reported_Staffing_Footnote $1. ;
format Physical_Therapist_Staffing_Foot $1. ;
format Reported_CNA_Staffing_Hours_per best12. ;
format Reported_LPN_Staffing_Hours_per best12. ;
format Reported_RN_Staffing_Hours_per_R best12. ;
format Reported_Licensed_Staffing_Hours best12. ;
format Reported_Total_Nurse_Staffing_Ho best12. ;
format Reported_Physical_Therapist_Staf best12. ;
format Expected_CNA_Staffing_Hours_per best12. ;
format Expected_LPN_Staffing_Hours_per best12. ;
format Expected_RN_Staffing_Hours_per_R best12. ;
format Expected_Total_Nurse_Staffing_Ho best12. ;
format Adjusted_CNA_Staffing_Hours_per best12. ;
format Adjusted_LPN_Staffing_Hours_per best12. ;
format Adjusted_RN_Staffing_Hours_per_R best12. ;
format Adjusted_Total_Nurse_Staffing_Ho best12. ;
format Cycle_1_Total_Number_of_Health_D best12. ;
format Cycle_1_Number_of_Standard_Healt best12. ;
format Cycle_1_Number_of_Complaint_Heal best12. ;
format Cycle_1_Health_Deficiency_Score best12. ;
format Cycle_1_Standard_Survey_Health_D mmddyy10. ;
format Cycle_1_Number_of_Health_Revisit best12. ;
format Cycle_1_Health_Revisit_Score best12. ;
format Cycle_1_Total_Health_Score best12. ;
format Cycle_2_Total_Number_of_Health_D best12. ;
format Cycle_2_Number_of_Standard_Healt best12. ;
format Cycle_2_Number_of_Complaint_Heal best12. ;
format Cycle_2_Health_Deficiency_Score best12. ;
format Cycle_2_Standard_Health_Survey_D mmddyy10. ;
format Cycle_2_Number_of_Health_Revisit best12. ;
format Cycle_2_Health_Revisit_Score best12. ;
format Cycle_2_Total_Health_Score best12. ;
format Cycle_3_Total_Number_of_Health_D best12. ;
format Cycle_3_Number_of_Standard_Healt best12. ;
format Cycle_3_Number_of_Complaint_Heal best12. ;
format Cycle_3_Health_Deficiency_Score best12. ;
format Cycle_3_Standard_Health_Survey_D mmddyy10. ;
format Cycle_3_Number_of_Health_Revisit best12. ;
format Cycle_3_Health_Revisit_Score best12. ;
format Cycle_3_Total_Health_Score best12. ;
format Total_Weighted_Health_Survey_Sco best12. ;
format Number_of_Facility_Reported_Inci best12. ;
format Number_of_Substantiated_Complain best12. ;
format Number_of_Fines best12. ;
format Total_Amount_of_Fines_in_Dollars nlnum12. ;
format Number_of_Payment_Denials best12. ;
format Total_Number_of_Penalties best12. ;
format Location $58. ;
format Processing_Date mmddyy10. ;
input
Federal_Provider_Number $
Provider_Name $
Provider_Address $
Provider_City $
Provider_State $
Provider_Zip_Code $
Provider_Phone_Number
Provider_SSA_County_Code
Provider_County_Name $
Ownership_Type $
Number_of_Certified_Beds
Number_of_Residents_in_Certified
Provider_Type $
Provider_Resides_in_Hospital $
Legal_Business_Name $
Date_First_Approved_to_Provide_M
Continuing_Care_Retirement_Commu $
Special_Focus_Facility $
Most_Recent_Health_Inspection_Mo $
Provider_Changed_Ownership_in_La $
With_a_Resident_and_Family_Counc $
Automatic_Sprinkler_Systems_in_A $
Overall_Rating
Overall_Rating_Footnote $
Health_Inspection_Rating
Health_Inspection_Rating_Footnot $
QM_Rating
QM_Rating_Footnote $
Staffing_Rating
Staffing_Rating_Footnote $
RN_Staffing_Rating
RN_Staffing_Rating_Footnote $
Reported_Staffing_Footnote $
Physical_Therapist_Staffing_Foot $
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
Cycle_1_Total_Number_of_Health_D
Cycle_1_Number_of_Standard_Healt
Cycle_1_Number_of_Complaint_Heal
Cycle_1_Health_Deficiency_Score
Cycle_1_Standard_Survey_Health_D
Cycle_1_Number_of_Health_Revisit
Cycle_1_Health_Revisit_Score
Cycle_1_Total_Health_Score
Cycle_2_Total_Number_of_Health_D
Cycle_2_Number_of_Standard_Healt
Cycle_2_Number_of_Complaint_Heal
Cycle_2_Health_Deficiency_Score
Cycle_2_Standard_Health_Survey_D
Cycle_2_Number_of_Health_Revisit
Cycle_2_Health_Revisit_Score
Cycle_2_Total_Health_Score
Cycle_3_Total_Number_of_Health_D
Cycle_3_Number_of_Standard_Healt
Cycle_3_Number_of_Complaint_Heal
Cycle_3_Health_Deficiency_Score
Cycle_3_Standard_Health_Survey_D
Cycle_3_Number_of_Health_Revisit
Cycle_3_Health_Revisit_Score
Cycle_3_Total_Health_Score
Total_Weighted_Health_Survey_Sco
Number_of_Facility_Reported_Inci
Number_of_Substantiated_Complain
Number_of_Fines
Total_Amount_of_Fines_in_Dollars
Number_of_Payment_Denials
Total_Number_of_Penalties
Location $
Processing_Date
;
if _ERROR_ then call symputx('_EFIERR_',1); 
run;

proc sort data=pen_&count out=pen_&count;
by Federal_Provider_Number Penalty_Date;
run;

proc transpose data=pen_&count prefix=Penalty_Date out=pen1_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Penalty_Date;
run;
proc transpose data=pen_&count prefix=Penalty_Type out=pen2_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Penalty_Type;
run;
proc transpose data=pen_&count prefix=Fine_Amount out=pen3_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Fine_Amount;
run;
proc transpose data=pen_&count prefix=Payment_Denial_Start_Date out=pen4_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Payment_Denial_Start_Date;
run;
proc transpose data=pen_&count prefix=Payment_Denial_Length_in_Days out=pen5_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Payment_Denial_Length_in_Days;
run;

data pen0_&count;
merge pen1_&count pen2_&count pen3_&count pen4_&count pen5_&count;
by Federal_Provider_Number;
run;
proc sort data=pen0_&count out=pen0_&count nodupkey;
 	by Federal_Provider_Number;
run;
data pen0_&count; set pen0_&count;
Total_Fine_Amount = sum(of Fine_Amount1-Fine_Amount8);
run;

proc sort data=mds_&count out=mds_&count;
by Federal_Provider_Number Measure_Code;
run;

data mds_&count; set mds_&count;
if missing (Federal_Provider_Number) then delete;
run;

proc transpose data=mds_&count prefix=Measure_Code out=mds1_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Measure_Code;
run;
proc transpose data=mds_&count prefix=Measure_Description out=mds2_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Measure_Description;
run;
proc transpose data=mds_&count prefix=Resident_type out=mds3_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Resident_type;
run;
proc transpose data=mds_&count prefix=Q1_Measure_Score out=mds4_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Q1_Measure_Score;
run;
proc transpose data=mds_&count prefix=Footnote_for_Q1_Measure_Score out=mds5_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Footnote_for_Q1_Measure_Score;
run;
proc transpose data=mds_&count prefix=Q2_Measure_Score out=mds6_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Q2_Measure_Score;
run;
proc transpose data=mds_&count prefix=Footnote_for_Q2_Measure_Scoree out=mds7_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Footnote_for_Q2_Measure_Score;
run;
proc transpose data=mds_&count prefix=Q3_Measure_Score out=mds8_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Q3_Measure_Score;
run;
proc transpose data=mds_&count prefix=Footnote_for_Q3_Measure_Score out=mds9_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Footnote_for_Q3_Measure_Score;
run;
proc transpose data=mds_&count prefix=Q4_Measure_Score out=mds10_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Q4_Measure_Score;
run;
proc transpose data=mds_&count prefix=Footnote_for_Q4_Measure_Score out=mds11_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Footnote_for_Q4_Measure_Score;
run;
proc transpose data=mds_&count prefix=Four_Quarter_Average_Score out=mds12_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Four_Quarter_Average_Score;
run;
proc transpose data=mds_&count prefix=Footnote_for_Four_Quarter_Averag out=mds13_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Footnote_for_Four_Quarter_Averag;
run;
proc transpose data=mds_&count prefix=Used_in_Quality_Measure_Five_Sta out=mds14_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Used_in_Quality_Measure_Five_Sta;
run;

/******************/;
proc datasets lib=work;
delete own_&count 
pen_&count pen0__&count pen1__&count pen2__&count pen3__&count pen4__&count pen5__&count
mdsll_&count 
run;
quit;
/******************/;

data mdst0_&count;
merge mds1_&count mds2_&count mds3_&count mds4_&count mds5_&count mds6_&count mds7_&count mds8_&count mds9_&count mds10_&count mds11_&count mds12_&count mds13_&count mds14_&count;
by Federal_Provider_Number;
run;

/******************/;
proc datasets lib=work;
delete own_&count 
pen_&count pen0__&count pen1__&count pen2__&count pen3__&count pen4__&count pen5__&count
mdsll_&count 
mds1_&count mds2_&count mds3_&count mds4_&count mds5_&count mds6_&count mds7_&count mds8_&count mds9_&count 
mds10_&count mds11_&count mds12_&count mds13_&count mds14_&count;
run;
quit;
/******************/;

proc sort data=mdst0_&count out=mdst0_&count nodupkey;
 	by Federal_Provider_Number;
run;

proc sort data=pi_&count out=piB_&count nodupkey;
 	by Federal_Provider_Number;
run;

data pi_&count; set piB_&count;
run;

proc sort data=mds_&count out=mdsB_&count nodupkey;
 	by Federal_Provider_Number;
run;

proc sql;
create table mdst_&count as
select A.*,  B.*, D.*, C.Federal_Provider_Number, C.Provider_Name, C.Provider_Address, C.Provider_City, C.Provider_State, 
C.Provider_Zip_Code, C.Q1_quarter, C.Q2_quarter, C.Q3_quarter, C.Q4_quarter, C.lat, C.lon, C.Processing_Date
from mdst0_&count as A Left Join pen0_&count as D on (A.Federal_Provider_Number=D.Federal_Provider_Number)
				Left Join mdsB_&count as C on (A.Federal_Provider_Number=C.Federal_Provider_Number)
				Left Join pi_&count as B on (A.Federal_Provider_Number=B.Federal_Provider_Number)
;
quit;

/******************/;
proc datasets lib=work;
delete own_&count pi_&count pib_&count
pen_&count pen0__&count pen1__&count pen2__&count pen3__&count pen4__&count pen5__&count
mdsll_&count mdst0_&count mdsb_&count
mds_&count 
mds1_&count mds2_&count mds3_&count mds4_&count mds5_&count mds6_&count mds7_&count mds8_&count mds9_&count 
mds10_&count mds11_&count mds12_&count mds13_&count mds14_&count;
run;
quit;
/******************/;

data mdst_&count; 
length PERIOD $20;
set mdst_&count;
PERIOD = "&time";
run;

data mylib.mdst_&count; set mdst_&count;
run;

proc export data=mdst_&count
    outfile="&dir./medicare_nursing_homes_&time..csv"
    dbms=csv
    replace;
run;


proc datasets lib=work;
delete own_&count pi_&count pib_&count
pen_&count pen0__&count pen1__&count pen2__&count pen3__&count pen4__&count pen5__&count
mdsll_&count mdst0_&count mdsb_&count
mds_&count mdst_&count
mds1_&count mds2_&count mds3_&count mds4_&count mds5_&count mds6_&count mds7_&count mds8_&count mds9_&count 
mds10_&count mds11_&count mds12_&count mds13_&count mds14_&count;
run;
quit;


	%end;
%exit: %mend Assemble_data;
%Assemble_data
/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/

/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;