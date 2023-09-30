/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;
*Transpose CMS MDS data fixed ;
*November 2018, Tyler Pittman;

proc datasets library=work kill memtype=data;

libname mylib "/folders/myshortcuts/tylerpittman/NursingHomesSNA/data";
%let dir = /folders/myshortcuts/tylerpittman/NursingHomesSNA/data;


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


data MDSadd_&count; set WORK.MDS_&count ;
keep Federal_Provider_Number Provider_Name Provider_Address Provider_City Provider_State Provider_Zip_Code ;
run;

proc sort data=MDSadd_&count out=MDSadd_&count;
by Federal_Provider_Number;
run;


proc sort data=mds_&count out=mds_&count nodupkey;
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
proc transpose data=mds_&count prefix=Four_Quarter_Average_Score out=mds12_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Four_Quarter_Average_Score;
run;
proc transpose data=mds_&count prefix=Used_in_Quality_Measure_Five_Sta out=mds14_&count (drop= _name_ _label_);
by Federal_Provider_Number;
var Used_in_Quality_Measure_Five_Sta;
run;



data mdst0_&count;
merge MDSadd_&count mds1_&count mds2_&count mds12_&count mds14_&count;
by Federal_Provider_Number;
run;

/******************/;
proc datasets lib=work;
delete  
MDSadd_&count mds1_&count mds2_&count mds12_&count mds14_&count;
run;
quit;
/******************/;

proc sort data=mdst0_&count out=mdst0_&count nodupkey;
 	by Federal_Provider_Number;
run;




data mylib.mdst_fixed_&count; set mdst0_&count;
run;

proc export data=mdst0_&count
    outfile="&dir./medicare_nursing_homes_fixed_&time..csv"
    dbms=csv
    replace;
run;


proc datasets lib=work;
delete mdst0_&count 
mds_&count MDSadd_&count
mds1_&count mds2_&count mds12_&count mds14_&count;
run;
quit;


	%end;
%exit: %mend Assemble_data;
%Assemble_data
/*-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^-------^^^^^^^^----------*/

/***---***---***---***---***---***---***---***---***---***---***---***---***---***---***/;