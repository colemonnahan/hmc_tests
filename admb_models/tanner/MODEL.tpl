// Updates of note:
// 2015-05-30:
// Added a fifth column to the tagging data (AIGKC)
// 2015-07-27
//  Calculate F35 based on mid-year SSB


GLOBALS_SECTION
  #define AGEMODEL 1
  #define LENMODEL 2
  #define AGELENMODEL 3
  #include <admodel.h>
  #include <time.h>
  #include "statsLib.h"
  ofstream CheckFile,OutFile1;
  ifstream GradFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
  
  adstring_array fleet_names;
  adstring hello;


// ===========================================================================
// ===========================================================================

 double CheckBounds(const prevariable& xx, double lower, double upper)
  {
   RETURN_ARRAYS_INCREMENT();  
   int Status;
   double Range;
   
   Status = 0;
   Range = upper - lower;
   if (xx < lower+Range*0.01) Status = 1;
   if (xx > upper-Range*0.01) Status = 2;
   OutFile1 << lower << " " << upper << " ";
   if (Status == 1) OutFile1 << "*l" << " ";
   if (Status == 2) OutFile1 << "*u" << " ";

   RETURN_ARRAYS_DECREMENT();  
   return (Status);
  }

 dvar_vector logistic(const dvector& x,  const prevariable& a,  const prevariable& b)
  {
   RETURN_ARRAYS_INCREMENT();  
   dvar_vector y;
   int n1,n2,Isize; n1 = x.indexmin(); n2 = x.indexmax(); 
   y=elem_div(mfexp(a*(x-b)), 1.0+mfexp(a*(x-b)));
   //for (Isize=n1;Isize<=n2;Isize++) y(Isize) /= y(n2);
 
   RETURN_ARRAYS_DECREMENT();  
   return (y);
  }

 dvar_vector logistic2(const dvector& x,  const prevariable& a,  const prevariable& b, const prevariable& c)
  {
   RETURN_ARRAYS_INCREMENT();  
   dvar_vector y;
   
   y=c*elem_div(mfexp(a*(x-b)), 1.0+mfexp(a*(x-b)));
 
   RETURN_ARRAYS_DECREMENT();  
   return (y);
  }

 dvar_vector logistic3(const dvector& x,  const prevariable& a,  const prevariable& b,  const prevariable& c)
  {
   RETURN_ARRAYS_INCREMENT();  
   dvar_vector y;
   int n1,n2,Isize; n1 = x.indexmin(); n2 = x.indexmax(); 
   y=elem_div(mfexp(a*(x-b)), 1.0+mfexp(a*(x-b)));
   for (Isize=n1;Isize<=n2;Isize++) y(Isize) = c*y(Isize)/y(n2);
 
   RETURN_ARRAYS_DECREMENT();  
   return (y);
  }

 dvar_vector doublelogistic(const dvector& MidLen,  const dvar_vector& sp,const double& binwidth2)
  {
   RETURN_ARRAYS_INCREMENT();  
   int Isize;
   int n1,n2; n1 = MidLen.indexmin(); n2 = MidLen.indexmax(); 
   dvar_vector y(n1,n2);
   dvariable peak,upselex,downselex,final,point1,point2,peak2;
   dvariable join1,join2,t1,t2;
   dvariable t1min,t2min;
   dvariable asc,dsc;

   peak = sp(1);
   upselex = exp(sp(3));
   downselex = exp(sp(4));
   final = sp(6);
   point1 = 1.0/(1.0+mfexp(-sp(5)));
   t1min = mfexp(-(square(MidLen(n1)-peak)/upselex)); 
   peak2 = peak + binwidth2+ (0.99*MidLen(n2)-peak-binwidth2)/(1.+mfexp(-sp(2)));
   point2 = 1.0/(1.0+mfexp(-final));
   t2min = mfexp(-(square(MidLen(n2)-peak2)/downselex));
   for (Isize=n1;Isize<=n2;Isize++)
    {
     t1 = MidLen(Isize)-peak;  t2 = MidLen(Isize)-peak2;
     join1 = 1.0/(1.0+mfexp(-(20.0*t1/(1.0+fabs(t1)))));                                           //  note the logit transform on t1 causes range of mfexp to be over -20 to 20
     join2 = 1.0/(1.0+mfexp(-(20.0*t2/(1.0+fabs(t2)))));
     if (sp(5) > -999)
      {asc = point1+(1.0-point1)*(mfexp(-square(t1)/upselex)-t1min)/(1.0-t1min);}
     else
      {asc = mfexp(-square(t1)/upselex);}
     if(sp(6) > -999)
      {dsc = 1.0+(point2-1.0)*(mfexp(-square(t2)/downselex)-1.0)/(t2min-1.0);}
     else
      {dsc = mfexp(-square(t2)/downselex);}
      y(Isize) = asc*(1.0-join1)+join1*(1.0-join2+dsc*join2);
     }
   RETURN_ARRAYS_DECREMENT();  
   return (y);
  } 


 dvar_vector splineAEP(const dvector& x, const dvar_vector& V,const dvector& MidLen)
  {
   RETURN_ARRAYS_INCREMENT();  
   int Isize;
   int n1,n2; n1 = MidLen.indexmin(); n2 = MidLen.indexmax(); 
   dvar_vector y(n1,n2);
   
   dvar_vector ders = spline(x,V,0,0);
   for (Isize=n1;Isize<=n2;Isize++) y(Isize) = sqrt(square(splint(x,V,ders,MidLen(Isize))));  
   
   RETURN_ARRAYS_DECREMENT();  
   return (y);
  }

 dvariable GenPrior(const prevariable& xx, double mean, double SD, const int Type)
  {
   RETURN_ARRAYS_INCREMENT();  
   dvariable Outcome;
   
   Outcome = 0;
   if (Type == 1) Outcome = log(SD) + square(xx-mean)/(2.0*SD*SD);

   RETURN_ARRAYS_DECREMENT();  
   return (Outcome);
  }



// ===========================================================================
// ===========================================================================

TOP_OF_MAIN_SECTION
  arrmblsize = 50000000;
  // Cole set this from 450000000 to 4500000
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(4500000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(30000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
//  arrmblsize = 50000000;
//  gradient_structure::set_GRADSTACK_BUFFER_SIZE(900000000);
//  gradient_structure::set_CMPDIF_BUFFER_SIZE(30000000);
//  gradient_structure::set_MAX_NVAR_OFFSET(2000);
//  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  time(&start);

  CheckFile.open("Check.Out");

// ======================================================================================================================
// ======================================================================================================================
// ======================================================================================================================

DATA_SECTION
 int IcntA;
 int IfleetA;
 int IyearA;
 int IsexA;
 int ItypeA;
 int IsizeA;
 int IageA;
 int Ilen1A;
 int Ilen2A;
 int III;
 int JJJ;
 int IparCnt;
 int BlockPnt;
 int KKK;
 int IsDataError;
 !! IsDataError = 0;
 int IsCTLError;
 !! IsCTLError = 0;
 number TotalA;

 init_int Model_Type;                                                              // Select a type of model
 init_int First_Year;
 init_int Last_Year;
 init_int Nsex;
 init_int Nage;
 int MaxAge
 !! MaxAge = Nage*5;
 !! if (MaxAge < 50) MaxAge = 50;
 int NageEst
 !! NageEst = Nage; if (NageEst < 20) NageEst = 20;                               // Trick to output growth curves
 init_int NsizeDym;
 int NsizeMax;
 int NsizeAge;
 !! if (Model_Type == AGEMODEL) { NsizeAge = NsizeDym; NsizeDym = 1; }
 !! if (Model_Type == LENMODEL | Model_Type == AGELENMODEL) { NsizeAge = 0; }
 !! if (Model_Type == LENMODEL) { Nage = 1; }
 !! if (NsizeAge > NsizeDym) NsizeMax = NsizeAge; else NsizeMax = NsizeDym;
 !! CheckFile << "NsizeDym NsizeAge NsizeMax" << endl << NsizeDym << " " << NsizeAge << " " << NsizeMax << endl;
 init_int Nplatoon;
 init_int NfishFleet;
 init_int NsurveyFleet;
 int Nfleet;
 !! Nfleet = NfishFleet + NsurveyFleet;
 
 // Read the timing of the survey fleets
 vector SurveyTime(1,Nfleet);
 !!for (IfleetA=NfishFleet+1;IfleetA<=Nfleet;IfleetA++)
 !! *(ad_comm::global_datafile) >> SurveyTime(IfleetA);

 init_adstring name_read_fleet;
 imatrix iname_fleet(1,Nfleet,1,2);
 !! CheckFile << name_read_fleet << endl;
 
 LOCAL_CALCS
   int k;
   adstring_array CRLF;   // Blank to terminate lines (not sure why this is needed...)
   CRLF+="";
   // Set whole array equal to 1 in case not enough names are read:
   for(k=1;k<=Nfleet;k++) 
   { iname_fleet(k,1)=1;  iname_fleet(k,2)=1; }    
   k=1;
   for(III=1;III<=strlen(name_read_fleet);III++)
    {
     if(adstring(name_read_fleet(III))==adstring(":")) 
      { iname_fleet(k,2)=III-1; k++;  iname_fleet(k,1)=III+1; }
    }
   iname_fleet(Nfleet,2)=strlen(name_read_fleet);
   for(k=1;k<=Nfleet;k++)
   {
     fleet_names += name_read_fleet(iname_fleet(k,1),iname_fleet(k,2))+CRLF(1);
   }
   CheckFile << "fleet_names: " << fleet_names << endl;
   CheckFile << "iname_fleet: " << iname_fleet << endl;
 END_CALCS 

 init_int LengthClassType;
 number Length1;                                                                     // First Length 
 number LengthInc;                                                                   // Length increment
 vector InputLengthClass(1,NsizeMax+1);
 !!if (LengthClassType == 1) *(ad_comm::global_datafile) >>  Length1 >> LengthInc;
 !!if (LengthClassType == 2)
 !! for (IsizeA=1;IsizeA<=NsizeMax+1;IsizeA++) *(ad_comm::global_datafile) >>  InputLengthClass(IsizeA);

 // Treatment of missing data
 init_int TreatMissF;

 // Catch data specifications
 init_matrix Catch_Specs(1,NfishFleet,1,3);
 !! CheckFile << "Catch Specifications" << endl << Catch_Specs << endl;
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) 
 !!  if (Catch_Specs(IfleetA,3) != 1 & Catch_Specs(IfleetA,3) != 0) 
 !!   { cout << "Retained / discarded for fleet " << IfleetA << " must be 0 or 1" << endl; IsDataError = 1; }

 // Initial catch data
 init_matrix InitialCinp(1,NfishFleet,0,3);
 vector InitialC(1,NfishFleet);
 vector SigmaIntC(1,NfishFleet);
 ivector PhaseInitC(1,NfishFleet);
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) 
 !!  {
 !!   InitialC(IfleetA) = InitialCinp(IfleetA,1);
 !!   SigmaIntC(IfleetA) = InitialCinp(IfleetA,2);
 !!   PhaseInitC(IfleetA) = int(InitialCinp(IfleetA,3));
 !!  }
 
 // Catch data input
 int NestCatch;
 !! NestCatch = 0;
 matrix CatchFix(1,NfishFleet,First_Year,Last_Year);
 imatrix CatchType(1,NfishFleet,First_Year,Last_Year);
 !! CatchFix.initialize();
 !! CatchType.initialize();
 init_int Ncatches;                                                                // Number of catches
 init_matrix CatchInput(1,Ncatches,1,4);
 int PhaseMissingF; 
 ivector FleetMissF(1,NfishFleet);
 !! CheckFile << "Catch Data" << endl;
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) FleetMissF(IfleetA) = -1;
 !! for (IcntA=1;IcntA<=Ncatches;IcntA++)
 !!  {
 !!   IyearA = int(CatchInput(IcntA,1)); IfleetA = int(CatchInput(IcntA,2)); 
 !!   CatchFix(IfleetA,IyearA) = CatchInput(IcntA,4)*Catch_Specs(IfleetA,2);
 !!   CatchType(IfleetA,IyearA) = int(CatchInput(IcntA,3));
 !!   CheckFile << IcntA << " " << CatchInput(IcntA) << " " << Catch_Specs(IfleetA,2) << " " << CatchFix(IfleetA,IyearA) << " " << CatchType(IfleetA,IyearA) << endl;
 !!   if (CatchInput(IcntA,4) < 0) { NestCatch += 1; if (TreatMissF==2) FleetMissF(IfleetA) = 1; }
 !!  }
 !! CheckFile << "Number of missing Fs" << endl << NestCatch << endl;
 !! CheckFile << "Number of fleets with missing Fs" << endl << FleetMissF << " " << sum(FleetMissF) << endl ;
 !! PhaseMissingF = -1;
 !! if (TreatMissF==1 || TreatMissF==3) PhaseMissingF = 1;
 !! if (TreatMissF==4) PhaseMissingF = -1;
 // Index Timing
 init_imatrix IndexTiming(1,Nfleet,1,2);                                           // Index timing
 
 // Index data input
 init_int NindexData;                                                              // Number of data points
 init_matrix IndexInput(1,NindexData,1,4);
 ivector IndexYear(1,NindexData);
 ivector IndexFleet(1,NindexData);
 matrix IndexData(1,NindexData,1,2);
 !! CheckFile << "Index Data" << endl;
 !! for (IcntA=1;IcntA<=NindexData;IcntA++)
 !!  {
 !!   IyearA = int(IndexInput(IcntA,1)); IfleetA = int(IndexInput(IcntA,2));
 !!   IndexYear(IcntA) = IyearA; IndexFleet(IcntA) = IfleetA;
 !!   IndexData(IcntA,1) = IndexInput(IcntA,3); IndexData(IcntA,2) = IndexInput(IcntA,4);
 !!   CheckFile << IcntA << " " << IndexInput(IcntA,1) << " " << IndexData(IcntA,1) << " " << IndexData(IcntA,2) << endl;
 !!  }
 !!  

 // Discard data specifications
 init_matrix Discard_Specs(1,NfishFleet,1,4);
 !! CheckFile << "Discard Specifications" << endl << Discard_Specs << endl;
 ivector Discard_units(1,NfishFleet);
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) Discard_units(IfleetA) = int(Discard_Specs(IfleetA,2));
 ivector Discard_type(1,NfishFleet);
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) Discard_type(IfleetA) = int(Discard_Specs(IfleetA,4));

 // Discard data input
 init_int NdiscardData;                                                             // Number of data points
 init_matrix DiscardInput(1,NdiscardData,1,5);                                      // Basic data
 ivector DiscardYear(1,NdiscardData);
 ivector DiscardFleet(1,NdiscardData);
 ivector DiscardSex(1,NdiscardData);
 matrix DiscardData(1,NdiscardData,1,2);
 !! CheckFile << "Discard Data" << endl;
 !! for (IcntA=1;IcntA<=NdiscardData;IcntA++)
 !!  {
 !!   IyearA = int(DiscardInput(IcntA,1)); IfleetA = int(DiscardInput(IcntA,2)); IsexA = int(DiscardInput(IcntA,3));
 !!   DiscardYear(IcntA) = IyearA;  DiscardFleet(IcntA) = IfleetA; DiscardSex(IcntA) = IsexA;
 !!   DiscardData(IcntA,1) = DiscardInput(IcntA,4)*Discard_Specs(IfleetA,3);
 !!   DiscardData(IcntA,2) = DiscardInput(IcntA,5);
 !!   CheckFile << IcntA << " " << DiscardInput(IcntA) << " " << Discard_Specs(IfleetA,3) << " " << Discard_Specs(IfleetA,4) << " " << DiscardData(IcntA,1) << " " << DiscardData(IcntA,2) << endl;
 !!  }
  
 // Effort data input
 init_int NeffortData;                                                              // Number of data points
 init_matrix EffortInput(1,NeffortData,1,4);                                        // Basic data
 ivector EffortYear(1,NeffortData);
 ivector EffortFleet(1,NeffortData);
 matrix EffortData(1,NeffortData,1,2);
 imatrix EffortByYear(1,NfishFleet,First_Year,Last_Year);
 vector MaxEffortByFleet(1,NfishFleet);                                             // Maximum effort for this fleet 
 !! CheckFile << "Effort Data" << endl;
 !! for (IcntA=1;IcntA<=NeffortData;IcntA++)
 !!  {
 !!   IyearA = int(EffortInput(IcntA,1)); IfleetA = int(EffortInput(IcntA,2));
 !!   EffortYear(IcntA) = IyearA; EffortFleet(IcntA) = IfleetA;
 !!   EffortData(IcntA,1) = EffortInput(IcntA,3); EffortData(IcntA,2) = EffortInput(IcntA,4);
 !!   EffortByYear(IfleetA,IyearA) = EffortInput(IcntA,3);
 !!   CheckFile << IcntA << " " << EffortData(IcntA,1) << " " << EffortData(IcntA,2) << endl;
 !!   if (EffortInput(IcntA,3) > MaxEffortByFleet(IfleetA)) MaxEffortByFleet(IfleetA) = EffortInput(IcntA,3);
 !!   if (EffortData(IcntA,1) <=0 || EffortData(IcntA,2) <=0) { cout << "Effort or its CV cannot be zero" << endl; exit(1); }
 !!  }

 // Length comp data
 init_int NlengthData                                                               // Number of data points
 init_int LengthCompLike;
 init_matrix LengthData(1,NlengthData,-4,Nsex*NsizeMax)
 !! CheckFile << "Length data" << endl << LengthData << endl;
 ivector LengthDataYear(1,NlengthData);
 ivector LengthDataFleet(1,NlengthData);
 ivector LengthDataSex(1,NlengthData);
 ivector LengthDataType(1,NlengthData);
 vector LengthDataSS(1,NlengthData);
 !! for (IcntA=1;IcntA<=NlengthData;IcntA++)
 !!  {
 !!   IyearA = int(LengthData(IcntA,-4)); IfleetA = int(LengthData(IcntA,-3));
 !!   IsexA = int(LengthData(IcntA,-2)); ItypeA = int(LengthData(IcntA,-1));
 !!   if (IfleetA > NfishFleet & ItypeA != 0) { CheckFile << "Error Survey data for fleet " << IfleetA << " In year " << IyearA << " should be 0" << endl; IsDataError = 1; }
 !!   LengthDataYear(IcntA) = IyearA; LengthDataFleet(IcntA) = IfleetA;
 !!   LengthDataSex(IcntA) = IsexA; LengthDataType(IcntA) = ItypeA;
 !!   LengthDataSS(IcntA) = LengthData(IcntA,0);
 !!   TotalA = 0;
 !!   for (IsizeA=1;IsizeA<=Nsex*NsizeMax;IsizeA++) TotalA += LengthData(IcntA,IsizeA);
 !!   for (IsizeA=1;IsizeA<=Nsex*NsizeMax;IsizeA++) LengthData(IcntA,IsizeA) /= TotalA;
 !!  }
 !! CheckFile << "Length data" << endl << LengthData << endl;
 imatrix LenMinMax(1,NlengthData,1,4);

 init_number PlusMinus;


 // Age-length comp data
 init_int NagelengthData;
 !! CheckFile << "NagelengthData: " << NagelengthData << endl;
 init_int AgeCompLike;
 init_matrix AgeLengthData(1,NagelengthData,-6,Nsex*Nage+Nsex)
 !! CheckFile << "AgeLength data" << endl << AgeLengthData << endl;
 ivector AgeLengthDataYear(1,NagelengthData);
 ivector AgeLengthDataFleet(1,NagelengthData);
 ivector AgeLengthDataSex(1,NagelengthData);
 ivector AgeLengthDataType(1,NagelengthData);
 ivector AgeLengthDataLen1(1,NagelengthData);
 ivector AgeLengthDataLen2(1,NagelengthData);
 vector AgeLengthDataSS(1,NagelengthData);
 !! for (IcntA=1;IcntA<=NagelengthData;IcntA++)
 !!  {
 !!   IyearA = int(AgeLengthData(IcntA,-6)); IfleetA = int(AgeLengthData(IcntA,-5));
 !!   IsexA = int(AgeLengthData(IcntA,-4)); ItypeA = int(AgeLengthData(IcntA,-3));
 !!   Ilen1A = int(AgeLengthData(IcntA,-2)); Ilen2A = int(AgeLengthData(IcntA,-1));
 !!   if (IfleetA > NfishFleet & ItypeA != 0) { CheckFile << "Error Survey data for fleet " << IfleetA << " In year " << IyearA << " should be 0" << endl; IsDataError = 1; }
 !!   AgeLengthDataYear(IcntA) = IyearA; AgeLengthDataFleet(IcntA) = IfleetA;
 !!   AgeLengthDataSex(IcntA) = IsexA; AgeLengthDataType(IcntA) = ItypeA;
 !!   AgeLengthDataLen1(IcntA) = Ilen1A; AgeLengthDataLen2(IcntA) = Ilen2A;
 !!   AgeLengthDataSS(IcntA) = AgeLengthData(IcntA,0);
 !!   TotalA = 0;
 !!   for (IageA=1;IageA<=Nsex*(Nage+1);IageA++) TotalA += AgeLengthData(IcntA,IageA);
 !!   for (IageA=1;IageA<=Nsex*(Nage+1);IageA++) AgeLengthData(IcntA,IageA) /= TotalA;
 !!  }
 !! CheckFile << "AgeLength data" << endl << AgeLengthData << endl;
 !! CheckFile << "AgeLength data complete" << endl;

 // Mean size-at-age data
 init_int NmeansizeData;
 ivector MeanSizeDataYear(1,NmeansizeData);
 ivector MeanSizeDataFleet(1,NmeansizeData);
 ivector MeanSizeDataSex(1,NmeansizeData);
 matrix MeanSizeData(1,NmeansizeData,1,Nsex*Nage+Nsex);
 matrix MeanSizeDataSS(1,NmeansizeData,1,Nsex*Nage+Nsex);
 !! for (IcntA=1;IcntA<=NmeansizeData;IcntA++)
 !! {
 !!  *(ad_comm::global_datafile) >> MeanSizeDataYear(IcntA) >> MeanSizeDataFleet(IcntA) >> MeanSizeDataSex(IcntA);
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   for (IageA=0;IageA<=Nage;IageA++) *(ad_comm::global_datafile) >> MeanSizeData(IcntA,(IsexA-1)*(Nage+1)+IageA+1);
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   for (IageA=0;IageA<=Nage;IageA++) *(ad_comm::global_datafile) >> MeanSizeDataSS(IcntA,(IsexA-1)*(Nage+1)+IageA+1);
 !! }
 !! CheckFile << "MeanSize data" << endl << MeanSizeData << endl;
 !! CheckFile << "MeanSize data complete" << endl;
 
 // Tagging data
 init_int NtagData;
 init_int TagYear;
 int NTagDiag2;
 init_matrix TagData(1,NtagData,1,5);
 vector TagOffset(1,NtagData);
 !! CheckFile << "Tagging data" << endl << TagData << endl;
 

 init_int CheckSum;
 !! if (CheckSum != 12345678) { cout << "CheckSum Error Dat File " << CheckSum << endl; exit(1); }
 !! if (IsDataError != 0) { cout << "Errors detected in DAT file; stopping; look at Check.Out" << endl; exit(1); }
 
 !! cout << "Completed Data Read-in" << endl;
  
 // ========================================================================================================================
 // ========================== CTL FILE ====================================================================================
 // ========================================================================================================================
 
 !! ad_comm::change_datafile_name("model.ctl");
 
 init_int Nblock;                                                                  // Number of blocks 
 imatrix Blocks(0,999,First_Year,Last_Year+1);
 !! for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++) Blocks(0,IyearA) = 1;
 !! for (IcntA=1;IcntA<=Nblock;IcntA++) for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)*(ad_comm::global_datafile) >> Blocks(IcntA,IyearA);
 ivector BlocksCnt(0,999);
 !! BlocksCnt.initialize();
 !! for (IcntA=1;IcntA<=Nblock;IcntA++)
 !!  {
 !!   for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!    if (Blocks(IcntA,IyearA) > BlocksCnt(IcntA)) BlocksCnt(IcntA) = Blocks(IcntA,IyearA);
 !!   CheckFile << "Blocks" << Blocks(IcntA) << endl;
 !!   for (III=1;III<=BlocksCnt(IcntA);III++)
 !!    {
 !!     JJJ=0;
 !!     for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!     if (Blocks(IcntA,IyearA) == III) JJJ = 1;
 !!     if (JJJ ==0) { cout << "Expected " << BlocksCnt(IcntA) << " years for block " << IcntA << endl; exit(1); }
 !!    }
 !!   BlocksCnt(IcntA) -= 1;
 !!  }
 !! CheckFile << "BlocksCnt" << endl << BlocksCnt << endl;
 
 init_int MatOpt;                                                                  // Maturity Option
 init_int MatAgeLen;                                                               // Maturity Option
 init_int GrowthOpt;                                                               // Growth Option
 init_int BiolOffset;                                                              // Biological parameters are offsets (1)
 init_int SigmaOrCV;                                                               // Parameterized in terms of a CV (1) or sigma (2)?
 init_int Jump;                                                                    // Max size classes 
 //!! Jump = 15;
 init_matrix RecPlat(1,Nsex,1,Nplatoon);                                           // Allocation of recruits to platoons 
 init_number PropWith;                                                             // Proportion of growth within and between
 !! if (Nplatoon==1) PropWith = 0;
 
 // Spline specifications for maturity
 int xknotsMat;                                                                    // Knots for Spline
 !! if (MatOpt == 2) *(ad_comm::global_datafile) >>  xknotsMat; else xknotsMat = 1;
 ivector xvalsMatL(1,xknotsMat); 
 !! if (MatOpt == 2) *(ad_comm::global_datafile) >>  xvalsMatL;

 // Set the base parameters
 int NBiolPar                                                                        // Number of parameters (total)
 !! NBiolPar = 2;                                                                    // Sex Ratio and Phi
 !! NBiolPar += 2*Nsex;                                                              // Length-weight 
 !! NBiolPar += 2*Nsex;                                                              // Natural mortality
 !! if (GrowthOpt == 1) NBiolPar += 3*Nsex;                                          // Vonbert Growth
 !! if (GrowthOpt == 2) NBiolPar += 7*Nsex;                                          // Linear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (GrowthOpt == 3) NBiolPar += 7*Nsex;                                          // Curvalinear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (GrowthOpt == 4) NBiolPar += 7*Nsex;                                          // Linear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (GrowthOpt == 5) NBiolPar += 7*Nsex;                                          // Curvalinear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (GrowthOpt == 6) NBiolPar += 7*Nsex;                                          // Linear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (GrowthOpt == 7) NBiolPar += 7*Nsex;                                          // Curvalinear (with molt prob) {len1; len2; molt1; molt2; growth1; growth2}
 !! if (MatOpt == 1) NBiolPar += 3*Nsex;
 !! if (MatOpt == 2) NBiolPar += Nsex*xknotsMat;
 !!  CheckFile << "Initial number of Biological Parameters = " << NBiolPar << endl;
 
 
 imatrix BiolInt(1,1000,1,6);
 matrix BiolReal(1,1000,1,6);
 
 !! JJJ = 0;
 !! for (IcntA=1;IcntA<=NBiolPar;IcntA++) 
 !!  {
 !!   JJJ += 1;
 !!   KKK = JJJ;
 !!   *(ad_comm::global_datafile) >> BiolReal(KKK,1) >> BiolReal(KKK,2) >> BiolReal(KKK,3) >> BiolInt(KKK,1);
 !!   *(ad_comm::global_datafile) >> BiolInt(KKK,6) >> BiolReal(KKK,5) >> BiolReal(KKK,6);
 !!   *(ad_comm::global_datafile) >> BiolInt(KKK,2) >> BiolInt(KKK,3) >> BiolInt(KKK,4) >> BiolInt(KKK,5) >> BiolReal(KKK,4);
 !!   if (BiolInt(KKK,3) != 0)
 !!    {
 !!     if (BiolInt(KKK,4) < First_Year+1) { cout << "WARNING: A biological dev vector starts before " << First_Year+1 << endl; exit(1); }
 !!     if (BiolInt(KKK,4) > Last_Year) { cout << "WARNING: A biological dev vector ends after " << Last_Year << endl; exit(1); }
 !!    }
 !!   BlockPnt = BiolInt(KKK,2);
 !!   if (BlockPnt > 0)
 !!    for (III=1;III<=BlocksCnt(BlockPnt);III++)
 !!     {
 !!      JJJ += 1;
 !!      *(ad_comm::global_datafile) >> BiolReal(JJJ,1) >> BiolReal(JJJ,2) >> BiolReal(JJJ,3) >> BiolInt(JJJ,1);
 !!      *(ad_comm::global_datafile) >> BiolInt(JJJ,6) >> BiolReal(JJJ,5) >> BiolReal(JJJ,6);
 !!     }
 !!   BlockPnt = BiolInt(KKK,3);
 !!   if (BlockPnt > 0)
 !!    for (III=BiolInt(KKK,4);III<=BiolInt(KKK,5);III++)
 !!      {
 !!       JJJ += 1;
 !!       BiolReal(JJJ,1) = -10; BiolReal(JJJ,2) = 10; BiolReal(JJJ,3) = 0; BiolInt(JJJ,1) = BiolInt(KKK,3);
 !!       BiolInt(JJJ,6) = 1; BiolReal(JJJ,5) = 0; BiolReal(JJJ,6) = BiolReal(KKK,4);
 !!      }
 !!  }
 !! NBiolPar = JJJ;
 !!  CheckFile << "Adjusted number of Biological Parameters = " << NBiolPar << endl;
 !!   CheckFile << "Refer to Section B.2 (of the .ctl file) for Biological Parameter Names" << endl;
 !!   CheckFile << "Param Num" << " " << "Min" << " " << "Max" << " " << "Init Num" << " " << "Phase" << endl;
 
 
 ivector BiolPhases(1,NBiolPar);                                                   // Phases
 vector BiolInit(1,NBiolPar);                                                      // Initial values
 vector BiolParMin(1,NBiolPar);                                                    // Lower bounds
 vector BiolParMax(1,NBiolPar);                                                    // Upper bounds
 !! for (IcntA=1;IcntA<=NBiolPar;IcntA++) 
 !!  {
 !!   BiolInit(IcntA) = BiolReal(IcntA,3); BiolPhases(IcntA) = BiolInt(IcntA,1);
 !!   BiolParMin(IcntA) = BiolReal(IcntA,1); BiolParMax(IcntA) = BiolReal(IcntA,2);
 !!   CheckFile << IcntA << " " << BiolParMin(IcntA) << " " << BiolParMax(IcntA) << " " << BiolInit(IcntA) << " " << BiolPhases(IcntA) << endl;
 !!  }
 
 int NSRPar; 
 !! NSRPar = 6;                                                                    //  SR pars 
 ivector SRPhases(1,NSRPar);                                                       // Phases
 vector SRInit(1,NSRPar);                                                          // initial values
 vector SRParMin(1,NSRPar);                                                        // Lower bounds
 vector SRParMax(1,NSRPar);                                                        // Upper bounds
 init_matrix SRParSet(1,NSRPar,1,4)                                                // Stock-recruit parameter set up
 init_int Last_Rec_size_Class;                                                     // Last recruitment size-class
 !! if (Last_Rec_size_Class > NsizeDym) Last_Rec_size_Class = NsizeDym;
 !! for (IcntA=1;IcntA<=NSRPar;IcntA++) 
 !!  {
 !!   SRInit(IcntA) = SRParSet(IcntA,3); SRPhases(IcntA) = int(SRParSet(IcntA,4));
 !!   SRParMin(IcntA) = SRParSet(IcntA,1); SRParMax(IcntA) = SRParSet(IcntA,2);
 !!  }
 !! CheckFile << "Stock and Recruitment Parameter Set" << endl << "R0 ; R1; Steepness ; SigmaR ; Alpha ; Beta" << endl << "Min Max Init Phase" << endl << SRParSet << endl;
 !! CheckFile << "last class to which animals recruit: " << Last_Rec_size_Class << endl;

 init_int RecEstYr1;                                                               // First year rec  
 int FirstProj_Yr;                                                                 // First year of model projection 
 !! if (RecEstYr1 < First_Year) FirstProj_Yr = RecEstYr1-1; else FirstProj_Yr = First_Year-1;
 !! FirstProj_Yr -= 2*NageEst;
 init_int RecEstYr2;                                                               // Last year rec
 init_int RecDevPhase;                                                             // Phase for rec_devs
 !! CheckFile << "Recruitment estimation years: " << RecEstYr1 << " : " << RecEstYr2 << endl;
 !! CheckFile << "First model projection year: " << FirstProj_Yr << endl;
 init_int PhaseYr1;                                                                // Phase in          
 init_int PhaseYr2;
 !! if (PhaseYr2 <= PhaseYr1) { cout << "PhaseY2 must be > PhaseY1: " << PhaseYr1 << " " << PhaseYr2 << endl; exit(1); }
 init_int PhaseYr3;
 init_int PhaseYr4;                                                                // Phase  
 !! if (PhaseYr4 <= PhaseYr3) { cout << "PhaseY4 must be > PhaseY3: " << PhaseYr3 << " " << PhaseYr4 << endl; exit(1); }
 vector SigPhi(RecEstYr1,RecEstYr2);
 init_int EarlyDevEnd;                                                             // When do the early devs end 
 init_number SigmaRecDev;                                                          // What is the sigma
 
 init_int F_tune;                                                                       // Number of times the hybrid method is applied 
 init_number max_harvest_rate;                                                          // Maximum harvest rate in the hybrid method
 
 // Read in Initial Fs
 init_matrix InitHrateParSet(1,NfishFleet,1,4)                                          // Initial F Par Setup
 !! CheckFile << "Initial Fs " << endl << InitHrateParSet << endl;
 vector InitHrateInit(1,NfishFleet);
 vector InitHrateMin(1,NfishFleet);
 vector InitHrateMax(1,NfishFleet);
 ivector InitHratePhases(1,NfishFleet);
 !! for (IcntA=1;IcntA<=NfishFleet;IcntA++) 
 !!  {
 !!   InitHrateInit(IcntA) = InitHrateParSet(IcntA,3); InitHratePhases(IcntA) = int(InitHrateParSet(IcntA,4));
 !!   InitHrateMin(IcntA) = InitHrateParSet(IcntA,1); InitHrateMax(IcntA) = InitHrateParSet(IcntA,2);
 !!  }

 // Selectivity parameterizarion
 init_int SelexParStyle;
 !! CheckFile << "SelexParStyle " << SelexParStyle << endl;

 // Read in age-selectivity specifications
 3darray SelexASpex(1,Nfleet,1,Nsex,1,7);                                             
 matrix SelexLSpexR(1,Nfleet,1,Nsex)                                      // Lambda 
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!  {
 !!   for (IcntA=1;IcntA<=7; IcntA++) *(ad_comm::global_datafile) >> SelexASpex(IfleetA,IsexA,IcntA);
 !!   *(ad_comm::global_datafile) >> SelexLSpexR(IfleetA,IsexA);
 !!  }
 !! CheckFile << "SelexASpex = " << endl << SelexASpex << endl;
 
 // Read in size-selectivity specifications
 3darray SelexLSpex(1,Nfleet,1,Nsex,1,7);                                             
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   for (IcntA=1;IcntA<=7; IcntA++) *(ad_comm::global_datafile) >> SelexLSpex(IfleetA,IsexA,IcntA);
 !! CheckFile << "SelexLSpex" << endl << SelexLSpex << endl;
  
 int BaseRetainAPars;                                                                 // Basic retention parameters
 int NRetainAPatterns;                                                                // Retention parameters
 int BlockRetainAPars;                                                                // Retention parameters (blocks)
 imatrix RetainAPnt(1,Nfleet*Nsex,First_Year,Last_Year+1);                           // Pointers between fleets and actual retain patterns  
 !! BaseRetainAPars = 0; BlockRetainAPars = 0; NRetainAPatterns = 0; IcntA = 0;                                         
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   {
 !!    IparCnt = 0;
 !!    if (SelexASpex(IfleetA,IsexA,6) == 0) IparCnt = 0;                            // Flat Selex
 !!    if (SelexASpex(IfleetA,IsexA,6) == 1) IparCnt = 3;                            // Logistic
 !!    if (SelexASpex(IfleetA,IsexA,6) > 0) IcntA += 1;
 !!    if (Catch_Specs(IfleetA,3) == 0 & SelexASpex(IfleetA,IsexA,6) >=0) 
 !!     { cout << "Age retained spec for fleet " << IfleetA << "must be -1 as all catch is discarded" << endl; IsCTLError = 1; }
 !!    BaseRetainAPars += IparCnt;
 !!    NRetainAPatterns += 1;
 !!    if (SelexASpex(IfleetA,IsexA,7) > 0) BlockRetainAPars += BlocksCnt(SelexASpex(IfleetA,IsexA,7))*IparCnt;
 !!    if (SelexASpex(IfleetA,IsexA,7) > 0) NRetainAPatterns += BlocksCnt(SelexASpex(IfleetA,IsexA,7));
 !!    for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!     {
 !!      if (SelexASpex(IfleetA,IsexA,7) > 0)
 !!        RetainAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = IcntA + Blocks(SelexASpex(IfleetA,IsexA,7),IyearA)-1;
 !!       else
 !!        RetainAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = IcntA;
 !!      }
 !!     if (SelexASpex(IfleetA,IsexA,7) > 0) IcntA += BlocksCnt(SelexASpex(IfleetA,IsexA,7));
 !!   }
 !! CheckFile << "Base RetainPars (age) = " << BaseRetainAPars << endl;
 !! CheckFile << "Number of retention (age) patterns =  " << NRetainAPatterns << endl;
 !! CheckFile << "Block RetainPars (age) = " << BlockRetainAPars << endl;
 !! CheckFile << RetainAPnt << endl;
 
 int BaseRetainLPars;                                                                 // Basic retention parameters
 int NRetainLPatterns;                                                                // Retention parameters
 int BlockRetainLPars;                                                                // Retention parameters (blocks)
 imatrix RetainLPnt(1,Nfleet*Nsex,First_Year,Last_Year+1);                           // Pointers between fleets and actual retain patterns  
 !! BaseRetainLPars = 0; BlockRetainLPars = 0; NRetainLPatterns = 0; IcntA = 0;                                         
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   {
 !!    IparCnt = 0;
 !!    if (SelexLSpex(IfleetA,IsexA,6) == 0) IparCnt = 0;                            // Flat Selex
 !!    if (SelexLSpex(IfleetA,IsexA,6) == 1) IparCnt = 3;                            // Logistic
 !!    if (SelexLSpex(IfleetA,IsexA,6) > 0) IcntA += 1;
 !!    if (Catch_Specs(IfleetA,3) == 0 & SelexLSpex(IfleetA,IsexA,6) >=0) 
 !!     { cout << "Length retained spec for fleet " << IfleetA << "must be -1 as all catch is discarded" << endl; IsCTLError = 1; }
 !!    BaseRetainLPars += IparCnt;
 !!    NRetainLPatterns += 1;
 !!    if (SelexLSpex(IfleetA,IsexA,7) > 0) BlockRetainLPars += BlocksCnt(SelexLSpex(IfleetA,IsexA,7))*IparCnt;
 !!    if (SelexLSpex(IfleetA,IsexA,7) > 0) NRetainLPatterns += BlocksCnt(SelexLSpex(IfleetA,IsexA,7));
 !!    for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!     {
 !!      if (SelexLSpex(IfleetA,IsexA,7) > 0)
 !!        RetainLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = IcntA + Blocks(SelexLSpex(IfleetA,IsexA,7),IyearA)-1;
 !!       else
 !!        RetainLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = IcntA;
 !!      }
 !!     if (SelexLSpex(IfleetA,IsexA,7) > 0) IcntA += BlocksCnt(SelexLSpex(IfleetA,IsexA,7));
 !!   }
 !! CheckFile << "Base RetainPars (length)= " << BaseRetainLPars << endl;
 !! CheckFile << "Number of retention (length) patterns  " << NRetainLPatterns << endl;
 !! CheckFile << "Block RetainPars (length) " << BlockRetainLPars << endl;
 !! CheckFile << RetainLPnt << endl;

 // Specifications for selectivity knots (age)
 matrix xvalsSelexA(1,Nfleet,0,Nage)                                                 // Knots
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   if (SelexASpex(IfleetA,IsexA,1) == 3)
 !!    for (IcntA=1;IcntA<=nknotsSelexA(IfleetA);IcntA++)
 !!    *(ad_comm::global_datafile) >> xvalsSelexA(IfleetA,IcntA);
 !! CheckFile << "Knots for age-specific selctivity" << endl;
 !! CheckFile << "nknotsSelexA" << endl << nknotsSelexA << endl;
 !! CheckFile << "xvalsSelexA" << endl << xvalsSelexA << endl;

 // Specifications for selectivity knots (length)
 matrix xvalsSelexL(1,Nfleet,1,NsizeMax)                                             // Knots
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   if (SelexLSpex(IfleetA,IsexA,1) == 3)
 !!    for (IcntA=1;IcntA<=nknotsSelexL(IfleetA);IcntA++)
 !!    *(ad_comm::global_datafile) >> xvalsSelexL(IfleetA,IcntA);
 !! CheckFile << "Knots for length-specific selctivity" << endl;
 !! CheckFile << "nknotsSelexL" << endl << nknotsSelexL << endl;
 !! CheckFile << "xvalsSelexL" << endl << xvalsSelexL << endl;

 int LLL;
 int MMM;
 int IFound;
 ivector Itest(1,10);
 ivector Jtest(1,10);

 // Read in selectivity parameters (ageh) (base)
 imatrix SelexAInt(1,1000,1,6);
 matrix SelexAReal(1,1000,1,6);
  
 int BaseSelexAPars;                                                                // Basic selectivity parameters 
 int NSelexAPatterns;                                                               // Selectivity parameters
 imatrix SelexAPnt(1,Nfleet*Nsex,First_Year,Last_Year+1);                           // Pointers between fleets and actual selex patterns  
 ivector nknotsSelexA(1,Nfleet)                                                     // Number of knots
 !! nknotsSelexA.initialize();
 !! SelexAReal.initialize();
 !! NSelexAPatterns = 0; JJJ = 0;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!    for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!      SelexAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = -1;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   {
 !!    if (SelexASpex(IfleetA,IsexA,1) == 0) IparCnt = 0;                                 // Flat Selex
 !!    if (SelexASpex(IfleetA,IsexA,1) == 1) IparCnt = 2;                                 // Logistic
 !!    if (SelexASpex(IfleetA,IsexA,1) == 2) IparCnt = 6;                                 // Double normal 
 !!    if (SelexASpex(IfleetA,IsexA,1) == 3) { IparCnt = SelexASpex(IfleetA,IsexA,2); nknotsSelexA(IfleetA) = IparCnt; }  // Spline
 !!    if (SelexASpex(IfleetA,IsexA,1) == 4) IparCnt = 3;                                 // Logistic with extra parameter
 !!    if (SelexASpex(IfleetA,IsexA,1) == 5) IparCnt = 0;                                 // Mirrored selectivity
 !!    if (SelexASpex(IfleetA,IsexA,1) == 5)
 !!     for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!      {
 !!       if (SelexASpex(IfleetA,IsexA,5) == 0 | SelexASpex(IfleetA,IsexA,5) == IsexA)
 !!        SelexAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = SelexAPnt(SelexASpex(IfleetA,IsexA,3),IyearA);
 !!       else
 !!        SelexAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = 0;
 !!      }
 !!    SelexASpex(IfleetA,IsexA,4) = IparCnt;
 !!    LLL = JJJ;
 !!    for (IcntA=1;IcntA<=IparCnt;IcntA++) 
 !!     {
 !!      JJJ += 1;
 !!      KKK = JJJ;
 !!      Jtest(IcntA) = KKK;
 !!      CheckFile << "Age-selex reading for fleet " << IfleetA << " and sex " << IsexA;
 !!      CheckFile << "; parameter is " << IcntA << " of " << IparCnt << endl;
 !!      *(ad_comm::global_datafile) >> SelexAReal(KKK,1) >> SelexAReal(KKK,2) >> SelexAReal(KKK,3) >> SelexAInt(KKK,1);
 !!      *(ad_comm::global_datafile) >> SelexAInt(KKK,6) >> SelexAReal(KKK,5) >> SelexAReal(KKK,6);
 !!      *(ad_comm::global_datafile) >> SelexAInt(KKK,2) >> SelexAInt(KKK,3) >> SelexAInt(KKK,4) >> SelexAInt(KKK,5) >> SelexAReal(KKK,4);
 !!      if (SelexAInt(KKK,3) != 0)
 !!        {
 !!         if (SelexAInt(KKK,4) < First_Year+1) { cout << "WARNING: An age-selectivity dev vector starts before " << First_Year+1 << "for " << IfleetA << " " << IsexA << endl; exit(1); }
 !!         if (SelexAInt(KKK,4) > Last_Year) { cout << "WARNING: An age-selectivity dev vector ends after " << Last_Year << endl; exit(1); }
 !!        }
 !!      CheckFile << SelexAReal(KKK) << " " << SelexAInt(KKK) << endl;
 !!      BlockPnt = SelexAInt(KKK,2);
 !!      if (BlockPnt > 0) CheckFile << "Need a block input here: " << BlockPnt << endl;
 !!      if (BlockPnt > 0)
 !!      for (III=1;III<=BlocksCnt(BlockPnt);III++)
 !!       {
 !!        JJJ += 1;
 !!        *(ad_comm::global_datafile) >> SelexAReal(JJJ,1) >> SelexAReal(JJJ,2) >> SelexAReal(JJJ,3) >> SelexAInt(JJJ,1);
 !!        CheckFile << "reading block offset" << III << endl;
 !!       }
 !!      BlockPnt = SelexAInt(KKK,3);
 !!      if (BlockPnt > 0)
 !!       {
 !!        Nblock = Nblock + 1;
 !!        SelexAInt(KKK,3) = Nblock;
 !!        SelexAInt(KKK,2) = Nblock;
 !!        MMM = 1;
 !!        for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!         {
 !!          if (IyearA < SelexAInt(KKK,4)) Blocks(Nblock,IyearA) = 1;
 !!          if (IyearA >=SelexAInt(KKK,4) & IyearA <=SelexAInt(KKK,5))
 !!           {
 !!            MMM += 1; JJJ += 1;
 !!            Blocks(Nblock,IyearA) = MMM;
 !!            SelexAReal(JJJ,1) = -10; SelexAReal(JJJ,2) = 10; SelexAReal(JJJ,3) = 0; SelexAInt(JJJ,1) = SelexAInt(JJJ,3);
 !!            SelexAReal(JJJ,4) = SelexAReal(KKK,4); SelexAInt(JJJ,1) = SelexAInt(KKK,1)+1;
 !!           }
 !!          if (IyearA > SelexAInt(KKK,5)) Blocks(Nblock,IyearA) = 1;
 !!         }
 !!        BlocksCnt(Nblock) = MMM-1;
 !!       }
 !!      }
 //    AEP FIX
 !!    if (SelexASpex(IfleetA,IsexA,1) != 5) NSelexAPatterns += 1;
 !!    for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),First_Year);
 !!    for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!     {
 !!      IFound = 0; 
 !!      for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexAInt(Jtest(IcntA),2),IyearA)) IFound = 1;
 !!      if (IFound == 1) NSelexAPatterns += 1;
 !!      if (SelexASpex(IfleetA,IsexA,1) != 5)
 !!       SelexAPnt((IfleetA-1)*Nsex+IsexA,IyearA) = NSelexAPatterns;
 !!      for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),IyearA);
 !!     }
 !!   }
 !! BaseSelexAPars = JJJ;
 !! CheckFile << "Number of age selectivity patterns =  " << NSelexAPatterns << endl;
 !! CheckFile << "Base SelexPars (age) = " << BaseSelexAPars << endl;
 !! CheckFile << SelexAPnt << endl;

 vector SelexAInit(1,BaseSelexAPars);
 vector SelexAParMin(1,BaseSelexAPars);
 vector SelexAParMax(1,BaseSelexAPars);
 ivector SelexAPhases(1,BaseSelexAPars);
 !! for (IcntA=1;IcntA<=BaseSelexAPars;IcntA++) 
 !!  {
 !!   SelexAInit(IcntA) = SelexAReal(IcntA,3); SelexAPhases(IcntA) = SelexAInt(IcntA,1);
 !!   SelexAParMin(IcntA) = SelexAReal(IcntA,1); SelexAParMax(IcntA) = SelexAReal(IcntA,2);
 !!   CheckFile << IcntA << " " << SelexAParMin(IcntA) << " " << SelexAParMax(IcntA) << " " << SelexAInit(IcntA) << " " << SelexAPhases(IcntA) << endl;
 !!  }

 // Read in selectivity parameters (length) (base)
 imatrix SelexLInt(1,1000,1,6);
 matrix SelexLReal(1,1000,1,6);
 
 // Read in selectivity parameters (length) (base)
 int BaseSelexLPars;                                                                // Basic selectivity parameters 
 int NSelexLPatterns;                                                               // Selectivity parameters
 imatrix SelexLPnt(1,Nfleet*Nsex,First_Year,Last_Year+1);                           // Pointers between fleets and actual selex patterns  
 ivector nknotsSelexL(1,Nfleet)                                                     // Number of knots
 !! nknotsSelexL.initialize();
 !! SelexLReal.initialize();
 !! NSelexLPatterns = 0; JJJ = 0;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!    for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!      SelexLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = -1;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++)
 !!  for (IsexA=1;IsexA<=Nsex;IsexA++)
 !!   {
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 0) IparCnt = 0;                                 // Flat Selex
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 1) IparCnt = 2;                                 // Logistic
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 2) IparCnt = 6;                                 // Double normal 
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 3) { IparCnt = SelexLSpex(IfleetA,IsexA,2); nknotsSelexL(IfleetA) = IparCnt; }  // Spline
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 4) IparCnt = 3;                                 // Logistic with extra parameter
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 5) IparCnt = 0;                                 // Mirrored selectivity
 !!    if (SelexLSpex(IfleetA,IsexA,1) == 5)
 !!     for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!      {
 !!       if (SelexLSpex(IfleetA,IsexA,5) == 0 | SelexLSpex(IfleetA,IsexA,5) == IsexA)
 !!        SelexLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = SelexLPnt(SelexLSpex(IfleetA,IsexA,3),IyearA);
 !!       else
 !!        SelexLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = 0;
 !!      }
 !!    SelexLSpex(IfleetA,IsexA,4) = IparCnt;
 !!    LLL = JJJ;
 !!    for (IcntA=1;IcntA<=IparCnt;IcntA++) 
 !!     {
 !!      JJJ += 1;
 !!      KKK = JJJ;
 !!      Jtest(IcntA) = KKK;
 !!      CheckFile << "length-selex Reading for fleet " << IfleetA << " and sex " << IsexA;
 !!      CheckFile << "; parameter is " << IcntA << " of " << IparCnt << endl;
 !!      *(ad_comm::global_datafile) >> SelexLReal(KKK,1) >> SelexLReal(KKK,2) >> SelexLReal(KKK,3) >> SelexLInt(KKK,1);
 !!      *(ad_comm::global_datafile) >> SelexLInt(KKK,6) >> SelexLReal(KKK,5) >> SelexLReal(KKK,6);
 !!      *(ad_comm::global_datafile) >> SelexLInt(KKK,2) >> SelexLInt(KKK,3) >> SelexLInt(KKK,4) >> SelexLInt(KKK,5) >> SelexLReal(KKK,4);
 !!      if (SelexLInt(KKK,3) != 0)
 !!        {
 !!         if (SelexLInt(KKK,4) < First_Year+1) { cout << "WARNING: An length-selectivity dev vector starts before " << First_Year+1 << endl; exit(1); }
 !!         if (SelexLInt(KKK,4) > Last_Year) { cout << "WARNING: An length-selectivity dev vector ends after " << Last_Year << endl; exit(1); }
 !!        }
 !!      CheckFile << SelexLReal(KKK) << " " << SelexLInt(KKK) << endl;
 !!      BlockPnt = SelexLInt(KKK,2);
 !!      if (BlockPnt > 0) CheckFile << "Need a block input here: " << BlockPnt << endl;
 !!      if (BlockPnt > 0)
 !!      for (III=1;III<=BlocksCnt(BlockPnt);III++)
 !!       {
 !!        JJJ += 1;
 !!        *(ad_comm::global_datafile) >> SelexLReal(JJJ,1) >> SelexLReal(JJJ,2) >> SelexLReal(JJJ,3) >> SelexLInt(JJJ,1);
 !!        CheckFile << "reading block offset" << III << endl;
 !!       }
 !!      BlockPnt = SelexLInt(KKK,3);
 !!      if (BlockPnt > 0)
 !!       {
 !!        Nblock = Nblock + 1;
 !!        SelexLInt(KKK,3) = Nblock;
 !!        SelexLInt(KKK,2) = Nblock;
 !!        MMM = 1;
 !!        for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!         {
 !!          if (IyearA < SelexLInt(KKK,4)) Blocks(Nblock,IyearA) = 1;
 !!          if (IyearA >=SelexLInt(KKK,4) & IyearA <=SelexLInt(KKK,5))
 !!           {
 !!            MMM += 1; JJJ += 1;
 !!            Blocks(Nblock,IyearA) = MMM;
 !!            SelexLReal(JJJ,1) = -10; SelexLReal(JJJ,2) = 10; SelexLReal(JJJ,3) = 0; SelexLInt(JJJ,1) = SelexLInt(JJJ,3);
 !!            SelexLReal(JJJ,4) = SelexLReal(KKK,4); SelexLInt(JJJ,1) = SelexLInt(KKK,1)+1;
 !!           }
 !!          if (IyearA > SelexLInt(KKK,5)) Blocks(Nblock,IyearA) = 1;
 !!         }
 !!        BlocksCnt(Nblock) = MMM-1;
 !!       }
 !!     }
 //    AEP FIX
 !!    if (SelexLSpex(IfleetA,IsexA,1) != 5)NSelexLPatterns += 1;
 !!    for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),First_Year);
 !!    CheckFile << IfleetA << " " << IsexA << " " << IparCnt << " "; for (IcntA=1;IcntA<=IparCnt;IcntA++) CheckFile << Itest(IcntA) << " "; CheckFile << endl;
 !!     for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!      {
 !!       IFound = 0; 
 !!       for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexLInt(Jtest(IcntA),2),IyearA)) IFound = 1;
 !!       if (IFound == 1) NSelexLPatterns += 1;
 !!       if (SelexLSpex(IfleetA,IsexA,1) != 5)
 !!        SelexLPnt((IfleetA-1)*Nsex+IsexA,IyearA) = NSelexLPatterns;
 !!       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),IyearA);
 !!      }
 !!   }
 !! BaseSelexLPars = JJJ;
 !! CheckFile << "Number of length selectivity patterns =  " << NSelexLPatterns << endl;
 !! CheckFile << "Base SelexPars (length) = " << BaseSelexLPars << endl;
 !! CheckFile << SelexLPnt << endl;
 !! CheckFile << "Updated blocks" << endl;
 !! for (IcntA=1;IcntA<=Nblock;IcntA++) CheckFile << Blocks(IcntA) << endl;

 vector SelexLInit(1,BaseSelexLPars);
 vector SelexLParMin(1,BaseSelexLPars);
 vector SelexLParMax(1,BaseSelexLPars);
 ivector SelexLPhases(1,BaseSelexLPars);
 !!   CheckFile << "Selectivity Parameters, Section B.5" << endl;
 !! for (IcntA=1;IcntA<=BaseSelexLPars;IcntA++) 
 !!  {
 !!   SelexLInit(IcntA) = SelexLReal(IcntA,3); SelexLPhases(IcntA) = SelexLInt(IcntA,1);
 !!   SelexLParMin(IcntA) = SelexLReal(IcntA,1); SelexLParMax(IcntA) = SelexLReal(IcntA,2);
 !!   CheckFile << IcntA << " " << SelexLParMin(IcntA) << " " << SelexLParMax(IcntA) << " " << SelexLInit(IcntA) << " " << SelexLPhases(IcntA) << endl;
 !!  }

 // Read in Retention parameters (age) (base)
  init_matrix RetainAParSet(1,BaseRetainAPars,1,4)                                        // Selex Par Setup
  !! CheckFile << "RetainAParSet = " << endl << RetainAParSet << endl;
  vector RetainAInit(1,BaseRetainAPars);
  vector RetainAParMin(1,BaseRetainAPars);
  vector RetainAParMax(1,BaseRetainAPars);
  ivector RetainAPhases(1,BaseRetainAPars);
  !! for (IcntA=1;IcntA<=BaseRetainAPars;IcntA++) 
  !!  {
  !!   RetainAInit(IcntA) = RetainAParSet(IcntA,3); RetainAPhases(IcntA) = int(RetainAParSet(IcntA,4));
  !!   RetainAParMin(IcntA) = RetainAParSet(IcntA,1); RetainAParMax(IcntA) = RetainAParSet(IcntA,2);
  !!  }
  init_matrix RetainABlkParSet(1,BlockRetainAPars,1,4)                                        // Selex Par Setup
  !! CheckFile << "RetainABlkParSet = " << endl << RetainABlkParSet << endl;
  vector RetainABlkInit(1,BlockRetainAPars);
  vector RetainABlkParMin(1,BlockRetainAPars);
  vector RetainABlkParMax(1,BlockRetainAPars);
  ivector RetainABlkPhases(1,BlockRetainAPars);
  !! for (IcntA=1;IcntA<=BlockRetainAPars;IcntA++) 
  !!  {
  !!   RetainABlkInit(IcntA) = RetainABlkParSet(IcntA,3); RetainABlkPhases(IcntA) = int(RetainABlkParSet(IcntA,4));
  !!   RetainABlkParMin(IcntA) = RetainABlkParSet(IcntA,1); RetainABlkParMax(IcntA) = RetainABlkParSet(IcntA,2);
  !!  }

 // Read in Retention parameters (length) (base)
 init_matrix RetainLParSet(1,BaseRetainLPars,1,4)                                        // Selex Par Setup
 !! CheckFile << "RetainLParSet = " << endl << RetainLParSet << endl;
 vector RetainLInit(1,BaseRetainLPars);
 vector RetainLParMin(1,BaseRetainLPars);
 vector RetainLParMax(1,BaseRetainLPars);
 ivector RetainLPhases(1,BaseRetainLPars);
 !! for (IcntA=1;IcntA<=BaseRetainLPars;IcntA++) 
 !!  {
 !!   RetainLInit(IcntA) = RetainLParSet(IcntA,3); RetainLPhases(IcntA) = int(RetainLParSet(IcntA,4));
 !!   RetainLParMin(IcntA) = RetainLParSet(IcntA,1); RetainLParMax(IcntA) = RetainLParSet(IcntA,2);
 !!  }
 init_matrix RetainLBlkParSet(1,BlockRetainLPars,1,4)                                        // Selex Par Setup
 !! CheckFile << "RetainLBlkParSet = " << endl << RetainLBlkParSet << endl;
 vector RetainLBlkInit(1,BlockRetainLPars);
 vector RetainLBlkParMin(1,BlockRetainLPars);
 vector RetainLBlkParMax(1,BlockRetainLPars);
 ivector RetainLBlkPhases(1,BlockRetainLPars);
 !! for (IcntA=1;IcntA<=BlockRetainLPars;IcntA++) 
 !!  {
 !!   RetainLBlkInit(IcntA) = RetainLBlkParSet(IcntA,3); RetainLBlkPhases(IcntA) = int(RetainLBlkParSet(IcntA,4));
 !!   RetainLBlkParMin(IcntA) = RetainLBlkParSet(IcntA,1); RetainLBlkParMax(IcntA) = RetainLBlkParSet(IcntA,2);
 !!  }

 // Q prior specifications
 ivector QpriorType(1,Nfleet);
 imatrix QpriorBlock(1,Nfleet,1,2);
 init_matrix QpriorSpex(1,Nfleet,0,3);
 !! CheckFile << "Catchability, Section B.6" << endl;
 !! CheckFile << "QpriorSpex" << endl << QpriorSpex << endl;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++) QpriorType(IfleetA) = int(QpriorSpex(IfleetA,0));
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++) 
 !!  {
 !!   QpriorBlock(IfleetA,2) = int(QpriorSpex(IfleetA,3));
 !!   QpriorBlock(IfleetA,1) = 1;
 !!   for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
 !!    if (Blocks(QpriorBlock(IfleetA,2),IyearA) > QpriorBlock(IfleetA,1)) QpriorBlock(IfleetA,1) =Blocks(QpriorBlock(IfleetA,2),IyearA);
 !!   CheckFile << "Qs for fleet " << IfleetA << " is " << QpriorBlock(IfleetA,1) << " " << QpriorBlock(IfleetA,2) << endl;
 !!  } 

 // Additional variance specifications
 ivector AddVarType(1,Nfleet);
 init_matrix AddVarSpex(1,Nfleet,0,4);
 vector AddVarParInit(1,Nfleet);
 vector AddVarParMin(1,Nfleet);
 vector AddVarParMax(1,Nfleet);
 ivector AddVarParPhases(1,Nfleet);
 !! CheckFile << "Additional Variance" << endl;
 !! CheckFile << "AddVarSpex" << endl << AddVarSpex << endl;
 !! for (IfleetA=1;IfleetA<=Nfleet;IfleetA++) 
 !! {
 !!  AddVarType(IfleetA) = int(AddVarSpex(IfleetA,0));
 !!  AddVarParMin(IfleetA) = AddVarSpex(IfleetA,1);
 !!  AddVarParMax(IfleetA) = AddVarSpex(IfleetA,2);
 !!  AddVarParInit(IfleetA) = AddVarSpex(IfleetA,3);
 !!  AddVarParPhases(IfleetA) = int(AddVarSpex(IfleetA,4));
 !!  if (AddVarType(IfleetA) != 1) AddVarParPhases(IfleetA) = -1;
 !! }
 !!
  
 // ========================================================================================================================

 // SigmaCatch
 init_number SigmaCatch;


 // Lambdas
 vector Catch_Lambda(1,Nfleet);
 !! Catch_Lambda.initialize();
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) *(ad_comm::global_datafile) >> Catch_Lambda(IfleetA);
 vector Discard_Lambda(1,Nfleet);
 !! Discard_Lambda.initialize();
 !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) *(ad_comm::global_datafile) >> Discard_Lambda(IfleetA);
 init_vector Index_Lambda(1,Nfleet);
 vector Effort_Lambda(1,Nfleet);
  !! Effort_Lambda.initialize();
  !! for (IfleetA=1;IfleetA<=NfishFleet;IfleetA++) *(ad_comm::global_datafile) >> Effort_Lambda(IfleetA);
 init_vector Length_Lambda(1,Nfleet);
 init_vector AgeLength_Lambda(1,Nfleet);
 init_vector MeanSize_Lambda(1,Nfleet);
 init_number Tag_Lambda;
 !! CheckFile << "Weights, Section B.7" << endl;
 !! CheckFile << "Catch_Lambda" << endl << Catch_Lambda << endl;
 !! CheckFile << "Discard_Lambda" << endl << Discard_Lambda << endl;
 !! CheckFile << "Index_Lambda" << endl << Index_Lambda << endl;
 !! CheckFile << "Effort_Lambda" << endl << Effort_Lambda << endl;
 !! CheckFile << "Length_Lambda" << endl << Length_Lambda << endl;
 !! CheckFile << "AgeLength_Lambda" << endl << AgeLength_Lambda << endl;
 !! CheckFile << "MeanSize_Lambda" << endl << MeanSize_Lambda << endl;

 init_int DiagLevel;
 init_int Stop_criterion;
 
 init_int CheckSum2;
 !! if (CheckSum2 != 12345678) { cout << "CheckSum Error CTL File " << CheckSum2 << endl; exit(1); }

 int IsFinalPhase;
 !! IsFinalPhase = 0;
 
 // ========================================================================================================================
 // ========================== MISC ====================================================================================
 // ========================================================================================================================

 int Year;                                                                         // Year counter
 int YearPass;                                                                     // Year to pass to biomass computation
 int FleetPass;                                                                    // Fleet to pass to biomass computation
 int TypePass;                                                                     // General Pass variable
 int DoTagDiag;                                                                    // Diagnostics for tagging datas

 int Proj_SR;                                                                      // Projection SR 
 
 int IFuncCallCount;
 !! IFuncCallCount = 0;
 vector StartLen(1,NsizeMax);
 vector MidLen(1,NsizeMax);
 vector EndLen(1,NsizeMax);
 vector Ages(0,Nage);
 
 ivector GrowthLow(1,NsizeDym);                                                   // Which size classes can growth into this size-class

 int NVarPar;
 !! NVarPar = 0;
 !! for (III=1;III<=NBiolPar;III++) if (BiolPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "BIOL " << NVarPar << endl;
 !! CheckFile << "Phases " << BiolPhases << endl;
 !! for (III=1;III<=NSRPar;III++) if (SRPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "SR " << NVarPar << endl;
 !! CheckFile << "Phases " << SRPhases << endl;
 !! for (III=1;III<=BaseSelexAPars;III++) if (SelexAPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "SA " << NVarPar << endl;
 !! CheckFile << "Phases " << SelexAPhases << endl;
 !! for (III=1;III<=BaseSelexLPars;III++) if (SelexLPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "SL " << NVarPar << endl;
 !! CheckFile << "Phases " << SelexLPhases << endl;
 !! for (III=1;III<=BaseRetainAPars;III++) if (RetainAPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "RA " << NVarPar << endl;
 !! CheckFile << "Phases " << RetainAPhases << endl;
 !! for (III=1;III<=BlockRetainAPars;III++) if (RetainABlkPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "Phases " << RetainABlkPhases << endl;
 !! CheckFile << "RA " << NVarPar << endl;
 !! for (III=1;III<=BaseRetainLPars;III++) if (RetainLPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "RL " << NVarPar << endl;
 !! CheckFile << "Phases " << RetainLPhases << endl;
 !! for (III=1;III<=BlockRetainLPars;III++) if (RetainLBlkPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "RL " << NVarPar << endl;
 !! CheckFile << "Phases " << RetainLBlkPhases << endl;
 !! if (PhaseMissingF > 0) {NVarPar += NestCatch; }
 !! CheckFile << "FM " << NVarPar << endl;
 !! CheckFile << "Phases " <<PhaseMissingF << endl;
 !! for (III=RecEstYr1;III<=RecEstYr2;III++) if (RecDevPhase > 0) {NVarPar += 1; }
 !! CheckFile << "REC " << NVarPar << endl;
 !! CheckFile << "Phases " <<RecDevPhase << endl;
 !! for (III=1;III<=NfishFleet;III++) if (InitHratePhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "INIT " << NVarPar << endl;
 !! CheckFile << "Phases " <<InitHratePhases << endl;
 !! for (III=1;III<=Nfleet;III++) if (AddVarParPhases(III) > 0) {NVarPar += 1; }
 !! CheckFile << "ADDVAR " << NVarPar << endl;
 !! CheckFile << "Phases " <<AddVarParPhases << endl;
 
 !! CheckFile << "No. Estimated Parameters " << NVarPar << endl;


 !! cout << "Completed Parameter Specification Read-in" << endl;

// ======================================================================================================================

PARAMETER_SECTION
 init_bounded_number_vector BiolPars(1,NBiolPar,BiolParMin,BiolParMax,BiolPhases);                  // Biological Parameters  
 init_bounded_number_vector SRPars(1,NSRPar,SRParMin,SRParMax,SRPhases);                            // Stock-recruit Parameters  
 init_bounded_number_vector SelexAPars(1,BaseSelexAPars,SelexAParMin,SelexAParMax,SelexAPhases);      // Selectivity parameters
 init_bounded_number_vector SelexLPars(1,BaseSelexLPars,SelexLParMin,SelexLParMax,SelexLPhases);      // Selectivity parameters
 init_bounded_number_vector RetainAPars(1,BaseRetainAPars,RetainAParMin,RetainAParMax,RetainAPhases); // Retension parameters
 init_bounded_number_vector RetainABlkPars(1,BlockRetainAPars,RetainABlkParMin,RetainABlkParMax,RetainABlkPhases);             // Selectivity parameters
 init_bounded_number_vector RetainLPars(1,BaseRetainLPars,RetainLParMin,RetainLParMax,RetainLPhases); // Retension parameters
 init_bounded_number_vector RetainLBlkPars(1,BlockRetainLPars,RetainLBlkParMin,RetainLBlkParMax,RetainLBlkPhases);             // Selectivity parameters
 init_bounded_vector MissingLogF(1,NestCatch,-10,0,PhaseMissingF);
 init_bounded_number_vector MissingLogq(1,NfishFleet,-10,0,FleetMissF);                             // Q for missing fleets
 init_bounded_vector RecDev(RecEstYr1,RecEstYr2,-10,10,RecDevPhase);                                // Recruit devs
 init_bounded_number_vector InitHrate(1,NfishFleet,InitHrateMin,InitHrateMax,InitHratePhases);      // Initial harvest rate
 init_bounded_number_vector AddVarPar(1,Nfleet,AddVarParMin,AddVarParMax,AddVarParPhases);          // Additional variance parameters
 init_number dummy;
 
 vector PriorVAPars(1,BaseSelexAPars);                                             // Prior component for age-specific sel pars 
 vector PriorVLPars(1,BaseSelexLPars);                                             // Prior component for len-specific sel pars 
 vector BiolPrior(1,NBiolPar);                                                     // Prior component for biological parameters
 
 3darray NatMBase(1,2,1,Nsex,First_Year,Last_Year);                                // Natural mortality (base level)
 5darray NatM(1,Nsex,First_Year,Last_Year,1,Nplatoon,0,Nage,1,NsizeDym);           // Natural mortality
 vector Linf(1,Nsex);                                                              // Linf
 vector Kappa(1,Nsex);                                                             // Kappa
 vector SigmaL(1,Nsex);                                                            // Sigma
 matrix RefLenGrow(1,Nsex,1,2);                                                    // Reference length
 matrix GrowInc(1,Nsex,1,2);                                                       // growth increment 
 matrix MoltProb(1,Nsex,1,2);                                                      // Molt probability
 
 vector Mat50(1,2);                                                                // Length-at-50% maturity
 vector MatSlope(1,2);                                                             // Maturity slope
 vector MatA(1,2);                                                                 // Maturity asymptotic
 vector NatSplinePars(1,2*xknotsMat);                                              // Maturity spline
 number phi;                                                                       // Timing of MMB output
 
 number Steepness;                                                                 // Steepness
 number R0;                                                                        // R0
 number R1Scalar;                                                                  // Multiplier for the first year
 number SigmaR;                                                                    // SigmaR

 5darray N(First_Year-1,Last_Year+1,1,Nsex,0,Nage,1,Nplatoon,1,NsizeDym);          // Numbers-at-age
 4darray Ntemp(1,Nsex,0,Nage,1,Nplatoon,1,NsizeDym);                               // Numbers-at-age after mortality
 4darray NtempA(1,Nsex,0,Nage,1,Nplatoon,1,NsizeDym);                              // Numbers-at-age after mortality
 vector RecDevMult(FirstProj_Yr,Last_Year+1);                                      // Recruitment multiplication

 3darray fec(1,Nsex,1,NsizeMax,0,Nage);                                            // Fecundity-at-length
 3darray fecAge(1,Nsex,1,Nplatoon,0,Nage);                                         // Fecundity-at-age 
 matrix mat(1,Nsex,1,NsizeMax);                                                    // Maturity-at-length
 matrix matA(1,Nsex,0,Nage);                                                       // Maturity-at-length
 3darray matAge(1,Nsex,1,Nplatoon,0,Nage);                                         // Maturity-at-age 
 4darray TransX(1,Nsex,1,Nplatoon,1,NsizeDym,1,NsizeDym);                          // Transition matrix 
 4darray PhiGrow(1,Nsex,1,Nplatoon,0,Nage,1,NsizeAge);                             // Phigrowth

 4darray Z_rate(1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                              // Temp Z
 4darray Z_rate2(1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                             // Temp Z
 matrix Hrate(1,NfishFleet,First_Year,Last_Year);                                  // Harvest rate
 matrix wghtL(1,Nsex,1,NsizeMax);                                                  // Weight-at-size   
 3darray selexL(1,Nfleet,1,Nsex,1,NsizeMax);                                       // Selectivity-at -sex and -size
 3darray retainL(1,Nfleet,1,Nsex,1,NsizeMax);                                      // Retention-at -sex and -size
 matrix wghtA(1,Nsex,0,Nage);                                                      // Weight-at-size   
 3darray selexA(1,Nfleet,1,Nsex,0,Nage);                                           // Selectivity-at -sex and -size
 3darray retainA(1,Nfleet,1,Nsex,0,Nage);                                          // Retention-at -sex and -size
 5darray selretwght1(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                // Retention-at -sex, -age and -size (numbers)
 5darray selretwght2(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                // Retention-at -sex, -age and -size (weights)
 5darray selretwghtUse(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);              // Retention-at -sex, -age and -size (used for computing F)
 5darray selwght1(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                   // Captured-at -sex, -age and -size (numbers)
 5darray selwght2(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                   // Captured-at -sex, -age and -size (weights)
 5darray selwghtUse(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                 // Captured-at -sex, -age and -size (used for computing F)
 5darray selexS(1,Nfleet,1,Nsex,1,Nplatoon,0,Nage,1,NsizeDym);                     // Selectivity accounting for discard survival
 
 matrix CatchPred(1,NfishFleet,First_Year,Last_Year);                              // Predicted catch
 matrix CatchPredN(1,NfishFleet,First_Year,Last_Year);                             // Predicted catch (numbers)
 matrix CatchPredB(1,NfishFleet,First_Year,Last_Year);                             // Predicted catch (biomass)
 6darray CatchPredAgeSize(1,3,1,Nfleet,First_Year,Last_Year+1,1,Nsex,0,Nage,1,NsizeMax);  // Predicted catch by component 
 vector InitCatch(1,NfishFleet);                                                   // Initial catch
 matrix MidExpBio(1,Nfleet,First_Year,Last_Year);                                  // Mid-year biomass
 matrix MidExpNum(1,Nfleet,First_Year,Last_Year);                                  // Mid-year numbers
 matrix SurveyBio(1,Nfleet,First_Year,Last_Year+1);                                // Specify predicted biomass  
 matrix BaseSelexA(0,NSelexAPatterns,0,Nage);                                      // Selectivity as a function of age
 matrix BaseRetainA(-1,NSelexAPatterns,0,Nage);                                     // Retention as a function of age
 matrix BaseSelexL(0,NSelexLPatterns,1,NsizeMax);                                  // Selectivity as a function of length
 matrix BaseRetainL(-1,NSelexLPatterns,1,NsizeMax);                                 // Retention as a function of length
 
 number MeanSexRatio;                                                              // Sex ratio in logit space
 vector SexRatioDev(FirstProj_Yr,Last_Year+1);                                     // Sex ratio devs
 
 3darray RecLen(1,Nsex,1,Nplatoon,1,NsizeDym);                                     // Distribution to recruits by size
 vector SSB(FirstProj_Yr,Last_Year+1);                                             // SSB  
 vector SSBMid(FirstProj_Yr,Last_Year);                                            // SSB  less phi

 matrix Catch(1,NfishFleet,First_Year,Last_Year);                                  // Catches (differentiable)
 vector TotalCatch(First_Year,Last_Year);                                          // Total catch   

 vector CatchLikeCompFleet(1,Nfleet);                                              // Catch by fleet
 matrix CatchLikeCompData(1,Nfleet,First_Year,Last_Year);                          // Catch by year and fleet
 vector LikeInitC(1,Nfleet);                                                       // Vector of likelihoods 

 vector DiscardPred(1,NdiscardData);                                               // Predicted discard
 vector DiscardLikeCompFleet(1,Nfleet);                                            // Discard by fleet
 vector DiscardLikeCompData(1,NdiscardData);                                       // Discard by data point

 matrix PredTag(1,NtagData,0,6);                                                   // Predicted stuff
 matrix PredProb(1,NtagData,1,NsizeMax)
 matrix TagDiag2(1,1000,1,9)                                                       // Place to store tag output          
 4darray PredTagTable(1,2,1,10,1,2,1,NsizeMax);                                    // Tabular output
 number TagLike;                                                                   // Tag likelihood
 3darray TagProps(1,Nsex,1,NsizeDym,1,Nplatoon)                                    // Proportion by year in each platoon

 matrix qest(1,Nfleet,First_Year,Last_Year+1);                                     // Catchability
 vector IndexLikeCompFleet(1,Nfleet);                                              // Index by fleet
 vector IndexLikeCompData(1,NindexData);                                           // Index by datum
 matrix IndexPred(1,2,1,NindexData);                                               // Predicted indexes
 vector AddVarParUse(1,Nfleet);                                                    // Additional variance 
 
 4darray LenMat(1,2,1,Nplatoon,1,NsizeMax,0,NageEst);                                 // Output of length matrix

 vector qestE(1,Nfleet);                                                           // Catchability
 vector EffortLikeCompFleet(1,Nfleet);                                             // Effort by fleet
 vector EffortLikeCompData(1,NeffortData);                                         // Effort by datum
 vector EffortPred(1,NeffortData);                                                 // Predicted effort
 
 matrix PredLengthComp(1,NlengthData,1,2*NsizeMax);                                // Predicted length comp  
 vector LengthLikeCompFleet(1,Nfleet);                                             // Length data by fleet
 vector LengthLikeCompData(1,NlengthData);                                         // Length data by datum

 
 matrix PredAgeLengthComp(1,NagelengthData,1,Nsex*Nage+Nsex);                            // Predicted age-length comp  
 vector AgeLengthLikeCompFleet(1,Nfleet);                                          // Age-length data by fleet
 vector AgeLengthLikeCompData(1,NagelengthData);                                   // Age-length data by datum
 5darray TransX2(1,2,1,6,1,Nplatoon,1,NsizeMax,1,NsizeMax);
 3darray TempX(1,Nplatoon,1,NsizeMax,1,NsizeMax);
 3darray TempY(1,Nplatoon,1,NsizeMax,1,NsizeMax);
 3darray MeanLP(1,Nsex,1,Nplatoon,0,MaxAge)
 3darray SDLP(1,Nsex,1,Nplatoon,0,MaxAge);

 matrix PredMeanSize(1,NmeansizeData,1,Nsex*Nage+Nsex);                            // Predicted Meansize 
 matrix PredMeanSizeSQ(1,NmeansizeData,1,Nsex*Nage+Nsex);                          // Predicted Meansize SD
 vector MeanSizeLikeCompFleet(1,Nfleet);                                           // Mean Size data by fleet
 vector MeanSizeLikeCompData(1,NmeansizeData);                                     // Mean Size data by datum

 number SSBPass;                                                                   // Pass for SSB
 vector HrateProj(1,NfishFleet);                                                   // Set the harvest rate
 vector OFL(1,NfishFleet);                                                         // OFL
   
 number DevPenal;                                                                  // Dev pelanlty 
 vector Penal(1,5);                                                                // Penalties

 sdreport_vector ParsOut(1,NVarPar);
 sdreport_vector SSBOut(FirstProj_Yr,Last_Year+1);
 sdreport_matrix RecOut(FirstProj_Yr,Last_Year+1,1,Nsex)                           // Recruitment
 objective_function_value Obj;
  
 !! cout << "Completed Parameter Specification" << endl;

// ======================================================================================================================
  
PRELIMINARY_CALCS_SECTION  
 int Ilen,II,IY,Ifleet,Isize,Jsize,Iyear,Iage,Icnt,IOK,Itag,Isex,Jtag,TimeAtLib;
 int MinSize,MaxSize;
 float Slope, Int, RelTot, Total, Returns;
 
 for (Ilen=1;Ilen<=NsizeMax;Ilen++)
  if (LengthClassType==1)
   {
    StartLen(Ilen) = (Ilen-1)*LengthInc+Length1;  
    MidLen(Ilen) = (Ilen-0.5)*LengthInc+Length1;  
    EndLen(Ilen) = Ilen*LengthInc+Length1;  
   }
  else 
   {
    StartLen(Ilen) = InputLengthClass(Ilen);  
    MidLen(Ilen) =  (InputLengthClass(Ilen)+InputLengthClass(Ilen+1))/2.0;  
    EndLen(Ilen) = InputLengthClass(Ilen+1); 
   }
 for (Iage=0;Iage<=Nage;Iage++) Ages(Iage) = float(Iage);
 CheckFile << "Lengths specified" << endl;
 CheckFile << MidLen << endl;
  
 // Set the initial  values of parameters
 for (II=1;II<=NBiolPar;II++) BiolPars(II) = BiolInit(II);
 for (II=1;II<=NSRPar;II++) SRPars(II) = SRInit(II);
 for (II=1;II<=BaseSelexAPars;II++) SelexAPars(II) = SelexAInit(II);
 for (II=1;II<=BaseRetainAPars;II++) RetainAPars(II) = RetainAInit(II);
 for (II=1;II<=BlockRetainAPars;II++)  RetainABlkPars(II) = RetainABlkInit(II);
 for (II=1;II<=BaseSelexLPars;II++) SelexLPars(II) = SelexLInit(II);
 for (II=1;II<=BaseRetainLPars;II++) RetainLPars(II) = RetainLInit(II);
 for (II=1;II<=BlockRetainLPars;II++)  RetainLBlkPars(II) = RetainLBlkInit(II);
 for (II=1;II<=NfishFleet;II++)  InitHrate(II) = InitHrateInit(II);
 for (II=1;II<=Nfleet;II++) AddVarPar(II) = AddVarParInit(II);
 if (DiagLevel > 4) CheckFile << "Initial values set" << endl;
  
 // Set the catches to the values that are read in
 for (IY=First_Year;IY<=Last_Year;IY++)
  for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
   Catch(Ifleet,IY) = CatchFix(Ifleet,IY);
 if (DiagLevel > 4) CheckFile << "Catches set" << endl;
 
 // Find limits of growth increments
 for (Isize=1;Isize<=NsizeDym;Isize++)
  {
   Jsize = Isize-Jump+1; 
   if (Jsize < 1) Jsize = 1;
   GrowthLow(Isize) = Jsize;
  }
 CheckFile << "GrowthLow" << endl << GrowthLow << endl;

 // Phase in 
 for (Iyear=RecEstYr1;Iyear<=PhaseYr1;Iyear++) SigPhi(Iyear) = 0;
 Slope = 1.0/float(PhaseYr2-PhaseYr1);
 if (PhaseYr1 <= RecEstYr2 & PhaseYr2 >= RecEstYr1 & PhaseYr2 <= RecEstYr2)
  for (Iyear=max(PhaseYr1,RecEstYr1);Iyear<=PhaseYr2;Iyear++) SigPhi(Iyear) = Slope*float(Iyear-PhaseYr1);
 if (PhaseYr2 >= RecEstYr1 & PhaseYr2 <= RecEstYr2 & PhaseYr3 >= RecEstYr1 & PhaseYr3 <= RecEstYr2)
  for (Iyear=PhaseYr2;Iyear<=PhaseYr3;Iyear++) SigPhi(Iyear) = 1;
 Slope = 1.0/float(PhaseYr4-PhaseYr3);
 if (PhaseYr3 >= RecEstYr1 & PhaseYr3 <= RecEstYr2 & PhaseYr4 >= RecEstYr1)
  for (Iyear=PhaseYr3;Iyear<=min(PhaseYr4,RecEstYr2);Iyear++) SigPhi(Iyear) = 1-Slope*float(Iyear-PhaseYr3);
 for (Iyear=PhaseYr4;Iyear<=RecEstYr2;Iyear++) SigPhi(Iyear) = 0;
 CheckFile << "Rec Phase In" << endl << SigPhi << endl;
 if (DiagLevel > 4) cout << "Phase In Set" << endl;
 
 for (Icnt=1;Icnt<=NlengthData;Icnt++)
  {
   if (LengthDataSex(Icnt) == 0) { MinSize = 1;MaxSize=Nsex*NsizeMax; }
   if (LengthDataSex(Icnt) == 1) { MinSize = 1;MaxSize=NsizeMax; }
   if (LengthDataSex(Icnt) == 2) { MinSize = NsizeMax+1;MaxSize=Nsex*NsizeMax; }

   RelTot = 0; IOK = 1;
   for (Isize=1;Isize<=NsizeMax & IOK == 1;Isize++) { RelTot += LengthData(Icnt,Isize); if (RelTot > PlusMinus) IOK = 0; }
   for (Jsize=1;Jsize<=Isize-1;Jsize++) { LengthData(Icnt,Isize) += LengthData(Icnt,Jsize);LengthData(Icnt,Jsize) = 0; }
   LenMinMax(Icnt,1) = Isize;
   
   RelTot = 0; IOK = 1;
   for (Isize=NsizeMax;Isize>=1 & IOK == 1;Isize--) { RelTot += LengthData(Icnt,Isize); if (RelTot > PlusMinus) IOK = 0; }
   for (Jsize=NsizeMax;Jsize>=Isize+1;Jsize--) { LengthData(Icnt,Isize) += LengthData(Icnt,Jsize); LengthData(Icnt,Jsize) = 0; }
   LenMinMax(Icnt,2) = Isize;

   RelTot = 0; IOK = 1;
   for (Isize=NsizeMax+1;Isize<=Nsex*NsizeMax & IOK == 1;Isize++) { RelTot += LengthData(Icnt,Isize); if (RelTot > PlusMinus) IOK = 0; }
   if (Isize>2*NsizeMax) Isize = Nsex*NsizeMax;
   for (Jsize=NsizeMax+1;Jsize<=Isize-1;Jsize++) { LengthData(Icnt,Isize) += LengthData(Icnt,Jsize);LengthData(Icnt,Jsize) = 0; }
   LenMinMax(Icnt,3) = Isize;
   
   RelTot = 0; IOK = 1;
   for (Isize=Nsex*NsizeMax;Isize>=NsizeMax+1 & IOK == 1;Isize--) { RelTot += LengthData(Icnt,Isize); if (RelTot > PlusMinus) IOK = 0; }
   for (Jsize=Nsex*NsizeMax;Jsize>=Isize+1;Jsize--) { LengthData(Icnt,Isize) += LengthData(Icnt,Jsize); LengthData(Icnt,Jsize) = 0; }
   LenMinMax(Icnt,4) = Isize;
  }
 CheckFile << "Length min-max" << endl << LenMinMax << endl;
 CheckFile << "Modified lengths" << endl << LengthData << endl; 
 CheckFile << "CheckSum" << CheckSum << CheckSum2 << endl;
   
 // Calculate offsets
 for (Itag=1;Itag<=NtagData;Itag++)
  {
   Isex = TagData(Itag,1); 
   Isize = TagData(Itag,2); 
   TimeAtLib = TagData(Itag,4); 
   Returns = float(TagData(Itag,5)); 
   Total = 0;
   for (Jtag=1;Jtag<=NtagData;Jtag++)
    if (TagData(Jtag,1)==Isex & TagData(Jtag,2)==Isize & TagData(Jtag,4)==TimeAtLib) Total += float(TagData(Jtag,5)); 
   TagOffset(Itag) = Returns/Total;
  }

 cout << "Preliminary Calcs Section Completed" << endl;
  
// ======================================================================================================================
PROCEDURE_SECTION
 if (DiagLevel > 4) cout << "Begin Procedure Section" << endl;
 int Isex,Iplatoon,Isize,Iage,Ifleet,IY,Ipnt,II,f,ImissF,Iyear;
 dvariable Fbar,NF;
 dvar_vector Flast(1,NfishFleet);

 if (DiagLevel > 4) cout << "VarParameters extract start" << endl;
 Ipnt = 0;
 for (II=1;II<=NBiolPar;II++) if (BiolPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = BiolPars(II); }
 for (II=1;II<=NSRPar;II++) if (SRPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = SRPars(II); }
 for (II=1;II<=BaseSelexAPars;II++) if (SelexAPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = SelexAPars(II); }
 for (II=1;II<=BaseSelexLPars;II++) if (SelexLPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = SelexLPars(II); }
 for (II=1;II<=BaseRetainAPars;II++) if (RetainAPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = RetainAPars(II); }
 for (II=1;II<=BlockRetainAPars;II++) if (RetainABlkPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = RetainABlkPars(II); }
 for (II=1;II<=BaseRetainLPars;II++) if (RetainLPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = RetainLPars(II); }
 for (II=1;II<=BlockRetainLPars;II++) if (RetainLBlkPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = RetainLBlkPars(II); }
 for (II=1;II<=NfishFleet;II++) if (FleetMissF(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = MissingLogq(II); }
 for (II=RecEstYr1;II<=RecEstYr2;II++) if (RecDevPhase > 0) {Ipnt +=1; ParsOut(Ipnt) = RecDev(II); }
 for (II=1;II<=NfishFleet;II++) if (InitHratePhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) =InitHrate(II); }
 for (II=1;II<=Nfleet;II++) if (AddVarParPhases(II) > 0) {Ipnt +=1; ParsOut(Ipnt) = AddVarPar(II); }
 if (DiagLevel > 4) cout << "VarParameters extracted" << endl;

 // Update the function call
 IFuncCallCount += 1;
 
 // First extract the data
 Extract_Pars();
 if (DiagLevel > 4) cout << "Parameters extracted" << endl;

 // Compute the total catch
 TotalCatch.initialize();
 Ipnt = 0;
 for (IY=First_Year;IY<=Last_Year;IY++)
  for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
   {
    if (CatchFix(Ifleet,IY) < 0)
     {
      if (TreatMissF == 1 || TreatMissF == 3)
       { Ipnt += 1; Hrate(Ifleet,IY) = exp(MissingLogF(Ipnt)); }
      if (TreatMissF == 2)
        Hrate(Ifleet,IY) = exp(MissingLogq(Ifleet))*EffortByYear(Ifleet,IY)/MaxEffortByFleet(Ifleet);
      if (TreatMissF == 4) Hrate(Ifleet,IY) = 0.01;
     } 
    else 
     TotalCatch(IY) += Catch(Ifleet,IY);
   }
 //if (DiagLevel > 4) cout << "Hrate extracted" << endl;
 //if (DiagLevel > 4) cout <<  Hrate << endl;

 // Create the growth matrix
 if (Model_Type == LENMODEL | Model_Type == AGELENMODEL) Growth_matrix();
 if (Model_Type == AGEMODEL) LenAge_matrix();
 if (DiagLevel > 4) cout << "Growth defined" << endl;

 // Create Maturity at length
 Get_Maturity();
 if (DiagLevel > 4) cout << "Maturity defined" << endl;

 // Define Selectivity
 SetLengthSelexL(); 
 if (DiagLevel > 4) cout << "Length-selectivity defined" << endl;
 SetLengthSelexA(); 
 if (DiagLevel > 4) cout << "Age-selectivity defined" << endl;

 // Define Retention
 SetLengthRetainL(); 
 if (DiagLevel > 4) cout << "Length Retention defined" << endl;
 SetLengthRetainA(); 
 if (DiagLevel > 4) cout << "Age Retention defined" << endl;

 N.initialize();
 CatchPred.initialize();
 CatchPredN.initialize();
 CatchPredB.initialize();
 DiscardPred.initialize();
 MidExpBio.initialize();
 MidExpNum.initialize();
 SurveyBio.initialize();
 CatchPredAgeSize.initialize();
 if (DiagLevel > 4) cout << "variables Initialized" << endl;

 // Set up the initial conditions
 Year = First_Year;
 Get_SelRet();
 if (DiagLevel > 4) cout << "SetRet for Initial" << endl;
 Initial_State();
 if (DiagLevel > 4) cout << "Initial State defined" << endl;

 // Do the annual update
 for (Year=First_Year;Year<=Last_Year;Year++)
  {
   Get_SelRet();
   Update_Dynamics();
  }
 if (DiagLevel > 4) cout << "Done Projection" << endl; 
 Flast = 0.01;
 if (TreatMissF == 4)
  for (ImissF=1;ImissF<=10;ImissF++)
   {
    for (f=1;f<=NfishFleet;f++)
     {
      Fbar = 0; NF = 0;
      for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
       if (CatchFix(f,Iyear) > 0) { Fbar += Hrate(f,Iyear); NF += 1; } 
      Fbar /= NF; 
      for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
       if (CatchFix(f,Iyear) < 0) Hrate(f,Iyear)= sqrt(Fbar*Flast(f)); 
      Flast(f) = Fbar; 
     }
   // Set up the initial conditions
   N.initialize();
   CatchPred.initialize();
   CatchPredN.initialize();
   CatchPredB.initialize();
   DiscardPred.initialize();
   MidExpBio.initialize();
   MidExpNum.initialize();
   SurveyBio.initialize();
   CatchPredAgeSize.initialize();
   Year = First_Year;
   Get_SelRet(); Initial_State();
    for (Year=First_Year;Year<=Last_Year;Year++)
     { Get_SelRet(); Update_Dynamics(); }
    if (DiagLevel > 4) cout << "Done Projection" << endl; 
   }
 
 //  compute the likelihood contributions for the catch data
 InitCLikelihood();
 CatchLikelihood();
 if (DiagLevel > 4)  cout << "Done Catch Likelihood" << endl; 
 //  compute the likelihood contributions for the discard data
 DiscardLikelihood();
 if (DiagLevel > 4)  cout << "Done Discard Likelihood" << endl; 
 //  compute the likelihood contributions for the index data
 IndexLikelihood();
 if (DiagLevel > 4)  cout << "Done Index Likelihood" << endl; 
 //  compute the likelihood contributions for the effort data
 EffortLikelihood();
 if (DiagLevel > 4)  cout << "Done Effort Likelihood" << endl; 
 //  compute the likelihood contributions for the length data 
 LengthLikelihood();
 if (DiagLevel > 4)  cout << "Done Length Likelihood" << endl; 
 ConditionalLikelihood();
 if (DiagLevel > 4)  cout << "Done Conditional age-at-length data" << endl; 
 MeanSizeLikelihood();
 if (DiagLevel > 4)  cout << "Done Conditional meansize data" << endl; 
 TaggingLikelihood();
 if (DiagLevel > 4)  cout << "Done tagging data" << endl; 
 //  compute the penalties
 Penalties(); 
 if (DiagLevel > 4)  cout << "Done Penalties" << endl; 


 Obj = dummy*dummy;
 Obj += sum(elem_prod(CatchLikeCompFleet,Catch_Lambda));
 Obj += sum(LikeInitC);
 Obj += sum(elem_prod(DiscardLikeCompFleet,Discard_Lambda));
 Obj += sum(elem_prod(IndexLikeCompFleet,Index_Lambda));
 Obj += sum(elem_prod(EffortLikeCompFleet,Effort_Lambda));
 Obj += sum(elem_prod(LengthLikeCompFleet,Length_Lambda));
 Obj += sum(elem_prod(AgeLengthLikeCompFleet,AgeLength_Lambda));
 Obj += sum(elem_prod(MeanSizeLikeCompFleet,MeanSize_Lambda));
 Obj += TagLike*Tag_Lambda;
 Obj += sum(Penal);
 if (DiagLevel > 1) cout << "CatchLikeCompFleet: " << CatchLikeCompFleet << " ";
 if (DiagLevel > 1) cout << " ; DiscardLikeCompFleet: " << DiscardLikeCompFleet << " ";
 if (DiagLevel > 1) cout << " ; IndexLikeCompFleet: " << IndexLikeCompFleet << " ";
 if (DiagLevel > 1)  cout << "; EffortLikeCompFleet : " << EffortLikeCompFleet << " ";
 if (DiagLevel > 1)  cout << " ; LengthLikeCompFleet: " << LengthLikeCompFleet << " ";
 if (DiagLevel > 1)  cout << " ; AgeLengthLikeCompFleet: " << AgeLengthLikeCompFleet << endl;
 if (DiagLevel > 1)  cout << " ; MeanSizeLikeCompFleet: " << MeanSizeLikeCompFleet << endl;
 if (DiagLevel > 1)  cout << " ; TagLike: " << TagLike << endl;
 if (DiagLevel > 1)  cout << " ; Penal: " << Penal << endl;
 if (DiagLevel > -1)  cout << "Phase: " << current_phase() << endl << " Func Call Count: " << IFuncCallCount << endl << " Obj Fun: " << Obj << endl;

 if (mceval_phase())
  {
    cout << Obj << " ";
    cout << CatchLikeCompFleet << " ";
    cout << DiscardLikeCompFleet << " ";
    cout << IndexLikeCompFleet << " ";
    cout << EffortLikeCompFleet << " ";
    cout << LengthLikeCompFleet << " ";
    cout << AgeLengthLikeCompFleet << " ";
    cout << HrateProj << " ";
    cout << OFL << " ";
    cout << Penal << " ";
    cout << Fbar << " ";
    cout << SSB(Last_Year+1) << " ";
    cout << RecOut(Last_Year+1,1) << " ";
    cout << RecOut(Last_Year+1,2) << endl;
    cout << "SSB" << endl;
    cout << SSB << endl;
    cout << "Pred Length Comp" << endl;
    cout << PredLengthComp << endl;
    //cout << FrancisOut << endl;                 //output 4 cols: obsv len, pred len, var pred len, resid
    cout << "Pred Age" << endl;
    cout << PredAgeLengthComp << endl;
    cout << "Effort" << endl;
    cout << EffortByYear << endl;
    cout << "Survey" << endl;
    cout << SurveyBio << endl;
  }

 SSBOut = SSB;
 if (IFuncCallCount == Stop_criterion) { CreateOutput();exit(1); }
 
 
// ================================================================================================================
// ================================================================================================================

FUNCTION Extract_Pars
 // Extract the biological parameters
 int Ilen,II,Ipar,Iyear,Isex,Iplatoon,Isize,NoBlock,Ifleet,IparStore;
 dvariable NatMM,NatMM2,Mval,AlphaRec,BetaRec,MeanRecLen,Total,TempV;
 dvar_vector Alpha(1,2);
 dvar_vector Beta(1,2);

 // Extract Biological parameters
 Ipar = 1; 
 DevPenal = 0;
 BiolPrior.initialize();
 NatMBase.initialize();
 for (Isex=1;Isex<=Nsex;Isex++)
  {
   if (DiagLevel > 4) cout << "extracting parameters for sex=" << Isex << endl;
   NoBlock = 0;
   IparStore = Ipar;
   if (Isex==1) NatMM = BiolPars(Ipar); else  NatMM = NatMM*exp(BiolPars(Ipar));
   if (BiolInt(Ipar,2)>0)
    {
     for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
      { if (Blocks(BiolInt(Ipar,2),Iyear)==1) Mval = NatMM; else Mval = NatMM*exp(BiolPars(Ipar+Blocks(BiolInt(Ipar,2),Iyear)-1)); NatMBase(1,Isex,Iyear) = Mval; }
     NoBlock = BlocksCnt(BiolInt(Ipar,2));
    }  
   else  
    for (Iyear=First_Year;Iyear<=Last_Year;Iyear++) NatMBase(1,Isex,Iyear) = NatMM;
   Ipar = Ipar + NoBlock + 1;
   if (BiolInt(IparStore,3)>0)
    {
     for (Iyear=BiolInt(IparStore,4);Iyear<=BiolInt(IparStore,5);Iyear++)
      { 
       NatMBase(1,Isex,Iyear) = NatMBase(1,Isex,Iyear) * exp(BiolPars(Ipar)); 
       //BiolPrior(Ipar) += square(BiolPars(Ipar))/(2.0*square(BiolReal(IparStore,6)));
       Ipar += 1;
      } 
    }
   NoBlock = 0;
   IparStore = Ipar;
   NatMM2 =  NatMBase(1,Isex,First_Year)*exp(BiolPars(Ipar));
   if (BiolInt(Ipar,2)>0)
    {
     for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
      { if (Blocks(BiolInt(Ipar,2),Iyear)==1) Mval = NatMM2; else Mval = NatMBase(1,Isex,Iyear)*exp(BiolPars(Ipar+Blocks(BiolInt(Ipar,2),Iyear)-1)); NatMBase(2,Isex,Iyear) = Mval; }
     NoBlock = BlocksCnt(BiolInt(Ipar,2));
    }  
   else  
    for (Iyear=First_Year;Iyear<=Last_Year;Iyear++) NatMBase(2,Isex,Iyear) = NatMM2;
   Ipar = Ipar + NoBlock + 1;
   if (BiolInt(IparStore,3)>0)
    {
     for (Iyear=BiolInt(IparStore,4);Iyear<=BiolInt(IparStore,5);Iyear++)
      { 
       NatMBase(2,Isex,Iyear) = NatMBase(1,Isex,Iyear) * exp(BiolPars(Ipar)); 
       DevPenal += square(BiolPars(Ipar))/(2.0*square(BiolReal(IparStore,6)));
       Ipar += 1;
      } 
    }
   
   // Female growth
   if (Isex==1)
    {
     if (GrowthOpt==1)
      { Linf(1) = BiolPars(Ipar); Kappa(1) = BiolPars(Ipar+1); SigmaL(1) = BiolPars(Ipar+2); }
     if (GrowthOpt==2 || GrowthOpt==3 || GrowthOpt==4 || GrowthOpt==5 || GrowthOpt==6 || GrowthOpt==7)
      { RefLenGrow(1,1) = BiolPars(Ipar); RefLenGrow(1,2) = BiolPars(Ipar+1);
        GrowInc(1,1) = BiolPars(Ipar+2); GrowInc(1,2) = BiolPars(Ipar+3); 
        MoltProb(1,1) = BiolPars(Ipar+4); MoltProb(1,2) = BiolPars(Ipar+5);
        SigmaL(1) = BiolPars(Ipar+6); }
    }
   else 
    {
     if (BiolOffset == 1)
      {
       if (GrowthOpt==1)
        { Linf(2) = Linf(1)*exp(BiolPars(Ipar)); Kappa(2) = Kappa(1)*exp(BiolPars(Ipar+1)); SigmaL(2) = SigmaL(1)*exp(BiolPars(Ipar+2));}
       if (GrowthOpt==2 || GrowthOpt==3 || GrowthOpt==4 || GrowthOpt==5 || GrowthOpt==6 || GrowthOpt==7)
        { RefLenGrow(2,1) = BiolPars(Ipar); RefLenGrow(2,2) = BiolPars(Ipar+1);
         GrowInc(2,1) = GrowInc(1,1)*exp(BiolPars(Ipar+2)); GrowInc(2,2) = GrowInc(1,2)*exp(BiolPars(Ipar+3)); 
         TempV = log(1.0/MoltProb(1,1)-1.0); MoltProb(2,1) = 1.0/(1+exp(TempV+BiolPars(Ipar+4)));
         TempV = log(1.0/MoltProb(1,2)-1.0); MoltProb(2,2) = 1.0/(1+exp(TempV+BiolPars(Ipar+5)));
         //MoltProb(2,1) = MoltProb(1,1)*exp(BiolPars(Ipar+4)); MoltProb(2,2) = MoltProb(1,2)*exp(BiolPars(Ipar+5)); 
         SigmaL(2) = SigmaL(1)*exp(BiolPars(Ipar+6)); }
      }
     else
      {
       if (GrowthOpt==1)
        { Linf(2) = BiolPars(Ipar); Kappa(2) = BiolPars(Ipar+1); SigmaL(2) = BiolPars(Ipar+2); }
       if (GrowthOpt==2 || GrowthOpt==3 || GrowthOpt==4 || GrowthOpt==5 || GrowthOpt==6 || GrowthOpt==7)
        { RefLenGrow(2,1) = BiolPars(Ipar); RefLenGrow(2,2) = BiolPars(Ipar+1);
          GrowInc(2,1) = BiolPars(Ipar+2); GrowInc(2,2) = BiolPars(Ipar+3); 
          MoltProb(2,1) = BiolPars(Ipar+4); MoltProb(2,2) = BiolPars(Ipar+5);
          SigmaL(2) = BiolPars(Ipar+6); }
      }
    }
   if (GrowthOpt==1) Ipar = Ipar+3;
   if (GrowthOpt==2) Ipar = Ipar+7;
   if (GrowthOpt==3) Ipar = Ipar+7;
   if (GrowthOpt==4) Ipar = Ipar+7;
   if (GrowthOpt==5) Ipar = Ipar+7;
   if (GrowthOpt==6) Ipar = Ipar+7;
   if (GrowthOpt==7) Ipar = Ipar+7;

   // Weight length
   Alpha(Isex) = BiolPars(Ipar);
   Beta(Isex) = BiolPars(Ipar+1);
   Ipar += 2;
   
   // Maturity
   if (MatOpt == 1) 
    { 
     Mat50(Isex) = BiolPars(Ipar);
     MatSlope(Isex) = BiolPars(Ipar+1);
     MatA(Isex) = BiolPars(Ipar+2);
     Ipar += 3;
    }
   if (MatOpt == 2) 
    {
     for (II = 1; II<=xknotsMat;II++) NatSplinePars(II+(Isex-1)*xknotsMat) = BiolPars(Ipar+II-1);
     Ipar = Ipar +xknotsMat;
    } 
   if (DiagLevel > 4) cout << "Parameters extract dor sex " << Isex << endl; 
   
  } 
 
 // Phi
 phi = BiolPars(Ipar);
 Ipar = Ipar + 1;
 
 // Sex ratio
 MeanSexRatio = BiolPars(Ipar);
 // Set up the sex-ratio multiplier
 SexRatioDev.initialize();
 NoBlock = 0;
 if (BiolInt(Ipar,2)>0)
  {
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    { 
     if (Blocks(BiolInt(Ipar,2),Iyear)>1) { SexRatioDev(Iyear) = BiolPars(Ipar+Blocks(BiolInt(Ipar,2),Iyear)-1); }
     NoBlock = BlocksCnt(BiolInt(Ipar,2));
    }  
  }  
 Ipar += (NoBlock + 1);
 if (DiagLevel > 4) cout << "Biological parameters extracted" << endl;

 // Extract stock-recruit parameters
 R0 = exp(SRPars(1));
 R1Scalar = exp(SRPars(2));
 Steepness = SRPars(3);
 SigmaR = SRPars(4);
 MeanRecLen = SRPars(5);
 BetaRec = SRPars(6);

 // Set up the recruitment multiplier
 for (Iyear=FirstProj_Yr;Iyear<=Last_Year+1;Iyear++) RecDevMult(Iyear) = 1;
 for (Iyear=RecEstYr1;Iyear<=RecEstYr2;Iyear++) RecDevMult(Iyear) = exp(RecDev(Iyear));
 
 // Weight-at-length
 for (Ilen=1;Ilen<=NsizeMax;Ilen++)
  for (Isex=1;Isex<=Nsex;Isex++)
   wghtL(Isex,Ilen) = Alpha(Isex)*pow(MidLen(Ilen),Beta(Isex));
 
 // Additional variance
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   if (AddVarType(Ifleet) == 1)
    AddVarParUse(Ifleet) = exp(AddVarPar(Ifleet));
   else
    AddVarParUse(Ifleet) = 0;
  } 
 if (DiagLevel > 4) cout << "Additional variance extracted" << endl; 
 
 // Recruitment
 RecLen.initialize();
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    AlphaRec = MeanRecLen/BetaRec;
    for (Isize=1;Isize<=Last_Rec_size_Class;Isize++)
     RecLen(Isex,Iplatoon,Isize) = pow(MidLen(Isize),(AlphaRec-1.0))*exp(-1*MidLen(Isize)/BetaRec);
    Total = sum(RecLen(Isex,Iplatoon));
    for (Isize=1;Isize<=Last_Rec_size_Class;Isize++) RecLen(Isex,Iplatoon,Isize) = RecPlat(Isex,Iplatoon)*RecLen(Isex,Iplatoon,Isize)/Total;
   }
 
// ================================================================================================================
// ================================ Biological and Fishery Functions ==============================================
// ================================================================================================================

FUNCTION Get_Maturity
 int II,Isex,Ioffset,Iage,Isize,Iplatoon,Iyear;
 dvector x(1,xknotsMat);
 dvar_vector VyI(1,xknotsMat);
 
 // Reset matyrity
 for (Isex=1;Isex<=Nsex;Isex++)
  {
   matA(Isex,0) = 0;
   for (Isize=1;Isize<=NsizeMax;Isize++) mat(Isex,Isize) = 1;
   for (Iage=1;Iage<=Nage;Iage++) matA(Isex,Iage) = 1;
  } 
 
 // Logistic curve (age)
 if (MatOpt==1 & MatAgeLen == 1) 
  for (Isex=1;Isex<=Nsex;Isex++)
   matA(Isex) = logistic2(Ages,  MatSlope(Isex),  Mat50(Isex), MatA(Isex));
   
 // Spline curve (age)
 if (MatOpt==2 & MatAgeLen == 1) 
   for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ioffset = (Isex-1)*xknotsMat;
    for (II = 1; II<=xknotsMat;II++) 
     {
      x(II) = Ages(xvalsMatL(II));
      VyI(II) = exp(NatSplinePars(II+Ioffset))/(1.0+exp(NatSplinePars(II+Ioffset)));
     } 
    matA(Isex) = splineAEP(x, VyI, MidLen);
   } 

 // Logistic curve (length)
 if (MatOpt==1 & MatAgeLen == 2) 
  for (Isex=1;Isex<=Nsex;Isex++)
   mat(Isex) = logistic2(MidLen,  MatSlope(Isex),  Mat50(Isex), MatA(Isex));
   
 // Spline curve (length)
 if (MatOpt==2 & MatAgeLen == 2) 
   for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ioffset = (Isex-1)*xknotsMat;
    for (II = 1; II<=xknotsMat;II++) 
     {
      x(II) = MidLen(xvalsMatL(II));
      VyI(II) = exp(NatSplinePars(II+Ioffset))/(1.0+exp(NatSplinePars(II+Ioffset)));
     } 
    mat(Isex) = splineAEP(x, VyI, MidLen);
   } 
   
 // Multiply proportion mature   
 if (Model_Type == LENMODEL | Model_Type == AGELENMODEL)
  for (Isex=1;Isex<=Nsex;Isex++) 
   for (Isize=1;Isize<=NsizeMax;Isize++) 
    for (Iage=1;Iage<=Nage;Iage++) 
     fec(Isex,Isize,Iage) = matA(Isex,Iage)*mat(Isex,Isize)*wghtL(Isex,Isize); 
 if (Model_Type == AGEMODEL)
  {
   fecAge.initialize();
   matAge.initialize();
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Iage=0;Iage<=Nage;Iage++)
      for (Isize=1;Isize<=NsizeAge;Isize++)
       {
        fec(Isex,Isize,Iage) = matA(Isex,Iage)*mat(Isex,Isize)*wghtL(Isex,Isize); 
        matAge(Isex,Iplatoon,Iage) += matA(Isex,Iage)*PhiGrow(Isex,Iplatoon,Iage,Isize)*mat(Isex,Isize);
        fecAge(Isex,Iplatoon,Iage) += PhiGrow(Isex,Iplatoon,Iage,Isize)*fec(Isex,Isize,Iage);
       } 
  }   
   
 // Set natural mortality
 for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    for (Iage=0;Iage<=Nage;Iage++)
    {
     if (Model_Type == LENMODEL | Model_Type == AGELENMODEL)
      for (Isize=1;Isize<=NsizeDym;Isize++)
       NatM(Isex,Iyear,Iplatoon,Iage,Isize) = NatMBase(1,Isex,Iyear)*matA(Isex,Iage)*mat(Isex,Isize)+(1-matA(Isex,Iage)*mat(Isex,Isize))*NatMBase(2,Isex,Iyear);
     if (Model_Type == AGEMODEL)
      for (Isize=1;Isize<=NsizeDym;Isize++)
       NatM(Isex,Iyear,Iplatoon,Iage,Isize) = NatMBase(1,Isex,Iyear)*matAge(Isex,Iplatoon,Iage)+(1-matAge(Isex,Iplatoon,Iage))*NatMBase(2,Isex,Iyear);
    }
   
   
// ================================================================================================================

FUNCTION Growth_matrix                                                       // Equations C.1-C.3
 // Specify Growth
 int Isex,Ilen,Jlen,Iplatoon,MaxLen2,Iplatoff;
 dvariable LenMid2,Len,Total,Zstat,Value,LenUpp,scale,shape;
 dvariable AlphaA,BetaA,Zed1,Zed2,AlphaB,BetaB,Probb;
 dvariable SigBetween,SigWithin,SigmaV;
 
 // set the platoon offset
 Iplatoff = (Nplatoon+1)/2;
 
 // Type 1 (Vonbert)
 if (GrowthOpt==1)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2+(Linf(Isex) - LenMid2)*(1-mfexp(-Kappa(Isex)));
        if (SigmaOrCV == 1) SigmaV = Len*SigmaL(Isex); else SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len + (Iplatoon-Iplatoff)*SigBetween;

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2-1;Jlen++)
         {
          LenUpp = EndLen(Jlen);
          Zstat = (LenUpp-Len)/SigWithin;
          Value = cumd_norm(Zstat);
          TransX(Isex,Iplatoon,Jlen,Ilen) = Value - Total;
          Total = Value;
         }
        TransX(Isex,Iplatoon,MaxLen2,Ilen) = 1.0 - Total; 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

 // Type 2 (Linear; normal)
 if (GrowthOpt==2)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (GrowInc(Isex,2) - GrowInc(Isex,1))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  GrowInc(Isex,1) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2 + (AlphaA + BetaA*LenMid2);
        if (SigmaOrCV == 1) SigmaV = Len*SigmaL(Isex); else SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len + (Iplatoon-Iplatoff)*SigBetween;
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2-1;Jlen++)
         {
          LenUpp = EndLen(Jlen);
          Zstat = (LenUpp-Len)/SigWithin;
          Value = cumd_norm(Zstat);
          if (Jlen==Ilen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb) + Probb*(Value - Total);
          else 
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*(Value - Total);
          Total = Value;
         }
        TransX(Isex,Iplatoon,MaxLen2,Ilen) = Probb*(1.0 - Total); 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   
 // Type 4 (Linear; log-normal)
 if (GrowthOpt==4)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (GrowInc(Isex,2) - GrowInc(Isex,1))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  GrowInc(Isex,1) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2 + (AlphaA + BetaA*LenMid2);
        SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len*exp((Iplatoon-Iplatoff)*SigBetween);
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2-1;Jlen++)
         {
          LenUpp = EndLen(Jlen);
          Zstat = (log(LenUpp)-log(Len))/SigWithin;
          Value = cumd_norm(Zstat);
          if (Jlen==Ilen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb) + Probb*(Value - Total);
          else 
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*(Value - Total);
          Total = Value;
         }
        TransX(Isex,Iplatoon,MaxLen2,Ilen) = Probb*(1.0 - Total); 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

 // Type 6 (Linear; gamma)
 if (GrowthOpt==6)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (GrowInc(Isex,2) - GrowInc(Isex,1))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  GrowInc(Isex,1) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2 + (AlphaA + BetaA*LenMid2);
        SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len*exp((Iplatoon-Iplatoff)*SigBetween);
        scale = 1.0/square(SigWithin);
        shape = Len/scale;
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2;Jlen++)
         {
          LenUpp = MidLen(Jlen);
          Value = pow(LenUpp,shape-1)*exp(-LenUpp/scale);
          TransX(Isex,Iplatoon,Jlen,Ilen) = Value;
          Total += Value;
         }
        for (Jlen=Ilen;Jlen<=MaxLen2;Jlen++)
         {
          if (Ilen == Jlen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb)+Probb*TransX(Isex,Iplatoon,Jlen,Ilen)/Total; 
          else
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*TransX(Isex,Iplatoon,Jlen,Ilen)/Total; 
         } 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

 // Type 3 (Log-Linear; normal)
 if (GrowthOpt==3)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (log(GrowInc(Isex,2)) - log(GrowInc(Isex,1)))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  log(GrowInc(Isex,1)) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2+ exp((AlphaA + BetaA*LenMid2));
        if (SigmaOrCV == 1) SigmaV = Len*SigmaL(Isex); else SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len + (Iplatoon-Iplatoff)*SigBetween;
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2-1;Jlen++)
         {
          LenUpp = EndLen(Jlen);
          Zstat = (LenUpp-Len)/SigWithin;
          Value = cumd_norm(Zstat);
          if (Jlen==Ilen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb) + Probb*(Value - Total);
          else 
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*(Value - Total);
          Total = Value;
         }
        TransX(Isex,Iplatoon,MaxLen2,Ilen) = Probb*(1.0 - Total); 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

 // Type 5 (Log-Linear; log-normal)
 if (GrowthOpt==5)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (log(GrowInc(Isex,2)) - log(GrowInc(Isex,1)))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  log(GrowInc(Isex,1)) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2 + exp((AlphaA + BetaA*LenMid2));
        SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len*exp((Iplatoon-Iplatoff)*SigBetween);
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2-1;Jlen++)
         {
          LenUpp = EndLen(Jlen);
          Zstat = (LenUpp-Len)/SigWithin;
          Value = cumd_norm(Zstat);
          if (Jlen==Ilen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb) + Probb*(Value - Total);
          else 
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*(Value - Total);
          Total = Value;
         }
        TransX(Isex,Iplatoon,MaxLen2,Ilen) = Probb*(1.0 - Total); 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

 // Type 7 (Log-Linear; gamma)
 if (GrowthOpt==7)
  {
   // Loop over sex, plattoon, etc
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      BetaA = (log(GrowInc(Isex,2)) - log(GrowInc(Isex,1)))/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      AlphaA =  log(GrowInc(Isex,1)) - BetaA*RefLenGrow(Isex,1);
      Zed1 = log(1.0/MoltProb(Isex,1)-1.0);
      Zed2 = log(1.0/MoltProb(Isex,2)-1.0);
      AlphaB = (Zed2-Zed1)/(RefLenGrow(Isex,2)-RefLenGrow(Isex,1));
      BetaB = -1.0*(Zed1/AlphaB-RefLenGrow(Isex,1));
      scale = 1.0/square(SigWithin);
      for (Ilen=1;Ilen<=NsizeDym-1;Ilen++)
       {
        LenMid2 = MidLen(Ilen);
        Len = LenMid2 + exp((AlphaA + BetaA*LenMid2));
        SigmaV = SigmaL(Isex);
        SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
        SigBetween = SigWithin * PropWith;
        Len = Len*exp((Iplatoon-Iplatoff)*SigBetween);
        Probb = 1.0/(1.0+exp(AlphaB*(LenMid2-BetaB)));
        scale = 1.0/square(SigWithin);
        shape = Len/scale;

        Total = 0;
        MaxLen2 = Ilen+Jump-1;
        if (MaxLen2 > NsizeDym) MaxLen2 = NsizeDym;
        for (Jlen=Ilen;Jlen<=MaxLen2;Jlen++)
         {
          LenUpp = MidLen(Jlen);
          Value = pow(LenUpp,shape-1)*exp(-LenUpp/scale);
          TransX(Isex,Iplatoon,Jlen,Ilen) = Value;
          Total += Value;
         }
        for (Jlen=Ilen;Jlen<=MaxLen2;Jlen++)
         {
          if (Ilen == Jlen)
           TransX(Isex,Iplatoon,Jlen,Ilen) = (1-Probb)+Probb*TransX(Isex,Iplatoon,Jlen,Ilen)/Total+1.0e-10; 
          else
           TransX(Isex,Iplatoon,Jlen,Ilen) = Probb*TransX(Isex,Iplatoon,Jlen,Ilen)/Total+1.0e-10; 
         } 
       }
      TransX(Isex,Iplatoon,NsizeDym,NsizeDym) = 1.0; 
     }
  }   

   
// ================================================================================================================

FUNCTION LenAge_matrix                                                       // Equations C.4-C.6
 // Specify Growth
 int Isex,Iage,Iplatoon,Ilen,Iplatoff;
 dvariable Len,Total,Zstat,Value;
 dvariable SigBetween,SigWithin,SigmaV;
 
 PhiGrow.initialize();
 
 // Loop over sex, plattoon, etc
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    PhiGrow(Isex,Iplatoon,0,1) = 1;
    for (Iage=1;Iage<=Nage;Iage++)
     {
      Len = Linf(Isex)*(1.0 - exp(-Kappa(Isex)*float(Iage)));
      if (SigmaOrCV == 1) SigmaV = Len*SigmaL(Isex); else SigmaV = SigmaL(Isex);
      SigWithin = SigmaV / sqrt(PropWith*PropWith+1.0);
      SigBetween = SigWithin * PropWith;
      Len = Len + (Iplatoon-Iplatoff)*SigBetween;
      Total = 0;
      for (Ilen=1;Ilen<=NsizeAge-1;Ilen++)
       {
        Zstat = (EndLen(Ilen)-Len)/SigWithin;
        Value = cumd_norm( Zstat );
        PhiGrow(Isex,Iplatoon,Iage,Ilen) = Value - Total;
        Total = Value;
       } 
      PhiGrow(Isex,Iplatoon,Iage,NsizeAge) = 1.0 - Total;
     } 
   }  

// ================================================================================================================
// ================================= Population dynamics parameters ===============================================
// ================================================================================================================

FUNCTION Get_SelRet
  int Ifleet,Isex,Ipnt,Isize,Iage,Iplatoon;
  dvariable sum1,sum2,sum3,sum4,sum5;
  dvar_vector temp1(1,NsizeMax);
  dvar_vector temp2(0,Nage);

  // Extract Age selectivity
  Ipnt = 0;
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++) 
   for (Isex=1;Isex<=Nsex;Isex++)
    { Ipnt += 1; selexA(Ifleet,Isex) = BaseSelexA(SelexAPnt(Ipnt,Year)); }
  // Extract size selectivity
  Ipnt = 0;
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++) 
   for (Isex=1;Isex<=Nsex;Isex++)
    { Ipnt += 1; selexL(Ifleet,Isex) = BaseSelexL(SelexLPnt(Ipnt,Year)); }
  if (DiagLevel > 4) cout << "Selectivity extracted" << endl;
  
  // Extract age retention
  Ipnt = 0;
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++) 
   for (Isex=1;Isex<=Nsex;Isex++)
    { Ipnt += 1;
      if (SelexASpex(Ifleet,Isex,6) > 0)
       retainA(Ifleet,Isex) = BaseRetainA(RetainAPnt(Ipnt,Year)); 
      else 
       if (SelexASpex(Ifleet,Isex,6) == 0)
        retainA(Ifleet,Isex) = BaseRetainA(0);
       else
        retainA(Ifleet,Isex) = BaseRetainA(-1);
    } 
  // Extract size retention
  Ipnt = 0;
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++) 
   for (Isex=1;Isex<=Nsex;Isex++)
    { Ipnt += 1; 
      if (SelexLSpex(Ifleet,Isex,6) > 0)
       retainL(Ifleet,Isex) = BaseRetainL(RetainLPnt(Ipnt,Year)); 
      else
       if (SelexLSpex(Ifleet,Isex,6) == 0)
        retainL(Ifleet,Isex) = BaseRetainL(0);
       else 
        retainL(Ifleet,Isex) = BaseRetainL(-1);
    } 
  if (DiagLevel > 4) cout << "Retention extracted" << endl;
    
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++) 
   for (Isex=1;Isex<=Nsex;Isex++)
    { 
     temp1 = retainL(Ifleet,Isex)+(1-SelexLSpexR(Ifleet,Isex))*(1-retainL(Ifleet,Isex));
     temp2 = retainA(Ifleet,Isex)+(1-SelexLSpexR(Ifleet,Isex))*(1-retainA(Ifleet,Isex));
     if (Model_Type==AGEMODEL)                                 // Equation A.3b
      {
       for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
        for (Iage=0;Iage<=Nage;Iage++)
         {
          sum1 = 0; sum2 = 0; sum3 = 0; sum4 = 0; sum5 = 0;
          for (Isize=1;Isize<=NsizeAge;Isize++)
           {
            sum1 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*retainL(Ifleet,Isex,Isize)*wghtL(Isex,Isize);
            sum2 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*temp1(Isize);
            sum3 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*retainL(Ifleet,Isex,Isize);
            sum4 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*wghtL(Isex,Isize);
            sum5 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize);
           }
           if (Ifleet <= NfishFleet)
            {
             selretwght1(Ifleet,Isex,Iplatoon,Iage,1) = sum3*selexA(Ifleet,Isex,Iage)*retainA(Ifleet,Isex,Iage); 
             selretwght2(Ifleet,Isex,Iplatoon,Iage,1) = sum1*selexA(Ifleet,Isex,Iage)*retainA(Ifleet,Isex,Iage); 
             if (CatchType(Ifleet,Year) == 1) selretwghtUse(Ifleet,Isex,Iplatoon,Iage,1) = selretwght1(Ifleet,Isex,Iplatoon,Iage,1);
             if (CatchType(Ifleet,Year) == 2) selretwghtUse(Ifleet,Isex,Iplatoon,Iage,1) = selretwght2(Ifleet,Isex,Iplatoon,Iage,1);
             selwght1(Ifleet,Isex,Iplatoon,Iage,1) = sum4*selexA(Ifleet,Isex,Iage); 
             selwght2(Ifleet,Isex,Iplatoon,Iage,1) = sum5*selexA(Ifleet,Isex,Iage); 
             if (CatchType(Ifleet,Year) == 1) selwghtUse(Ifleet,Isex,Iplatoon,Iage,1) = selwght1(Ifleet,Isex,Iplatoon,Iage,1);
             if (CatchType(Ifleet,Year) == 2) selwghtUse(Ifleet,Isex,Iplatoon,Iage,1) = selwght2(Ifleet,Isex,Iplatoon,Iage,1);
            } 
          selexS(Ifleet,Isex,Iplatoon,Iage,1) = sum2*selexA(Ifleet,Isex,Iage)*temp2(Iage);
         } 
      }
     if (Model_Type==LENMODEL || Model_Type == AGELENMODEL)      // Equation A.3a
      {
       for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
        for (Iage=0;Iage<=Nage;Iage++)
         for (Isize=1;Isize<=NsizeDym;Isize++)
          {
           if (Ifleet <= NfishFleet)
            {
             selretwght1(Ifleet,Isex,Iplatoon,Iage,Isize) = selexA(Ifleet,Isex,Iage)*selexL(Ifleet,Isex,Isize)*retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage); 
             selretwght2(Ifleet,Isex,Iplatoon,Iage,Isize) = selexA(Ifleet,Isex,Iage)*selexL(Ifleet,Isex,Isize)*retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage)*wghtL(Isex,Isize); 
             if (CatchType(Ifleet,Year) == 1) selretwghtUse(Ifleet,Isex,Iplatoon,Iage,Isize) = selretwght1(Ifleet,Isex,Iplatoon,Iage,Isize);
             if (CatchType(Ifleet,Year) == 2) selretwghtUse(Ifleet,Isex,Iplatoon,Iage,Isize) = selretwght2(Ifleet,Isex,Iplatoon,Iage,Isize);
             selwght1(Ifleet,Isex,Iplatoon,Iage,Isize) = selexA(Ifleet,Isex,Iage)*selexL(Ifleet,Isex,Isize); 
             selwght2(Ifleet,Isex,Iplatoon,Iage,Isize) = selexA(Ifleet,Isex,Iage)*selexL(Ifleet,Isex,Isize)*wghtL(Isex,Isize); 
             if (CatchType(Ifleet,Year) == 1) selwghtUse(Ifleet,Isex,Iplatoon,Iage,Isize) = selwght1(Ifleet,Isex,Iplatoon,Iage,Isize);
             if (CatchType(Ifleet,Year) == 2) selwghtUse(Ifleet,Isex,Iplatoon,Iage,Isize) = selwght2(Ifleet,Isex,Iplatoon,Iage,Isize);
            }  
           selexS(Ifleet,Isex,Iplatoon,Iage,Isize) = selexA(Ifleet,Isex,Iage)*temp2(Iage)*selexL(Ifleet,Isex,Isize)*temp1(Isize);
          } 
      }
    }
  if (DiagLevel > 4) cout << "Sel Ret defined" << endl;
 

// ================================================================================================================

FUNCTION Initial_State
 // Set up the initial conditions
 
 int Isex,Iplatoon,Isize,Jsize,Iage,II,Iyear,Ifleet;
 dvar_vector NextL(1,NsizeDym),OldL(1,NsizeDym);
 dvar_matrix NNext(0,Nage,1,NsizeDym);
 dvariable Test,Recr,Surv,SSBTemp,Numer,Denom,Recruit,BioAtLength;
 dvar_vector SexRatio(1,Nsex);
 
 // Determine total mortality
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Isize=1;Isize<=NsizeDym;Isize++)
     Z_rate(Isex,Iplatoon,Iage,Isize) = NatM(Isex,First_Year,Iplatoon,Iage,Isize);
 
 SexRatio(1) = 2.0/(1+exp(MeanSexRatio));
 if (Nsex == 2) SexRatio(2) = 2.0 - SexRatio(1);

 // Specify year FirstYear-1 age and size structure
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
   
    // Specify Age 0
    for (Isize=1;Isize<=NsizeDym;Isize++)
     {
      Recr =RecLen(Isex,Iplatoon,Isize)*R0*SexRatio(Isex);
      N(First_Year-1,Isex,0,Iplatoon,Isize) = Recr;
     } 
     
    // Specify Ages 1-Nage
    for (Iage=1;Iage<=Nage;Iage++)
     {
      if (Model_Type==AGEMODEL)
       {
         Surv = exp(-NatM(Isex,First_Year,Iplatoon,Iage-1,1));
         N(First_Year-1,Isex,Iage,Iplatoon,1) = N(First_Year-1,Isex,Iage-1,Iplatoon,1)*Surv;
       }
      else
       {
        NextL.initialize();     
        for (Jsize=1;Jsize<=NsizeDym;Jsize++)
         for (Isize=GrowthLow(Jsize);Isize<=Jsize;Isize++)
          NextL(Jsize) += N(First_Year-1,Isex,Iage-1,Iplatoon,Isize)*TransX(Isex,Iplatoon,Jsize,Isize);
        for (Isize=1;Isize<=NsizeDym;Isize++)
         {
          Surv = exp(-NatM(Isex,First_Year,Iplatoon,Iage-1,Isize));
          N(First_Year-1,Isex,Iage,Iplatoon,Isize) = NextL(Isize)*Surv;
         } 
       }  
     }
    
    // Plus Group
    for (Isize=1;Isize<=NsizeDym;Isize++) OldL(Isize) = N(First_Year-1,Isex,Nage,Iplatoon,Isize);
    for (Iage=1;Iage<=2*MaxAge;Iage++)
     if (Model_Type==AGEMODEL)
      {
       Surv = exp(-NatM(Isex,First_Year,Iplatoon,Nage,1));
       OldL(1) = OldL(1)*Surv;
       N(First_Year-1,Isex,Nage,Iplatoon,1) += OldL(1);
      }
     else 
      {
       NextL.initialize();     
       for (Jsize=1;Jsize<=NsizeDym;Jsize++)
        for (Isize=GrowthLow(Jsize);Isize<=Jsize;Isize++)
         NextL(Jsize) += OldL(Isize)*TransX(Isex,Iplatoon,Jsize,Isize);
       for (Isize=1;Isize<=NsizeDym;Isize++)
        {
         Surv = exp(-NatM(Isex,First_Year,Iplatoon,Nage,Isize));
         OldL(Isize) = NextL(Isize)*Surv;
         N(First_Year-1,Isex,Nage,Iplatoon,Isize) += OldL(Isize);
        } 
      }   
     
   }

 // Determine total mortality
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Isize=1;Isize<=NsizeDym;Isize++)
     {
      Z_rate(Isex,Iplatoon,Iage,Isize) = NatM(Isex,First_Year,Iplatoon,Iage,Isize);
      for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
       Z_rate(Isex,Iplatoon,Iage,Isize) += selexS(Ifleet,Isex,Iplatoon,Iage,Isize) * InitHrate(Ifleet);
      Z_rate2(Isex,Iplatoon,Iage,Isize) = (1-mfexp(-Z_rate(Isex,Iplatoon,Iage,Isize)))/Z_rate(Isex,Iplatoon,Iage,Isize);  
     }  
   
 // SSB (females only in this version)
 SSB(FirstProj_Yr) = 0;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=1;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
     SSB(FirstProj_Yr) +=  N(First_Year-1,Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage);   
    else
     for (Isize=1;Isize<=NsizeDym;Isize++)
      SSB(FirstProj_Yr) +=  N(First_Year-1,Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage);   
     
 // Initialize Ntemp   
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Isize=1;Isize<=NsizeDym;Isize++)
     Ntemp(Isex,Iage,Iplatoon,Isize) = N(First_Year-1,Isex,Iage,Iplatoon,Isize);

 RecOut(FirstProj_Yr,1) = R0*SexRatio(1);
 if (Nsex == 2) RecOut(FirstProj_Yr,2) = R0*SexRatio(2);

 // Adjust by R1scalar
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=0;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
     {
      Ntemp(Isex,Iage,Iplatoon,1) *= R1Scalar;   
       }
      else
       { 
        for (Isize=1;Isize<=NsizeDym;Isize++)
         Ntemp(Isex,Iage,Iplatoon,Isize) *= R1Scalar;   
       }  

 // Project forward (at least one year)
 for (Iyear=FirstProj_Yr;Iyear<=First_Year-1;Iyear++)
  {
   
   // Find the mid-year SSB 
   SSBMid(Iyear) = 0;
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iage=1;Iage<=Nage;Iage++)
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      if (Model_Type==AGEMODEL)
       {
        SSBMid(Iyear) +=  Ntemp(Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,1));   
       }
      else
       { 
        for (Isize=1;Isize<=NsizeDym;Isize++)
         SSBMid(Iyear) +=  Ntemp(Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,Isize));   
       }  
   
   // Multiply by transition matrix and project
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     {
      NNext.initialize();
      // Ages 1-N-1
      for (Iage=1;Iage<=Nage-1;Iage++)
       if (Model_Type==AGEMODEL)
        NNext(Iage,1) += Ntemp(Isex,Iage-1,Iplatoon,1)*exp(-Z_rate(Isex,Iplatoon,Iage-1,1));
       else 
        {
         for (Isize=1;Isize<=NsizeDym;Isize++)
          for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
           NNext(Iage,Isize) += Ntemp(Isex,Iage-1,Iplatoon,Jsize)*TransX(Isex,Iplatoon,Isize,Jsize)*exp(-Z_rate(Isex,Iplatoon,Iage-1,Jsize));
        }
      // Plus-group 
      if (Model_Type==AGEMODEL)
       NNext(Nage,1) += Ntemp(Isex,Nage,Iplatoon,1)*exp(-Z_rate(Isex,Iplatoon,Nage,1))+Ntemp(Isex,Nage-1,Iplatoon,1)*exp(-Z_rate(Isex,Iplatoon,Nage-1,1));
      else
       for (Isize=1;Isize<=NsizeDym;Isize++)
        for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
         NNext(Nage,Isize) += (Ntemp(Isex,Nage,Iplatoon,Jsize)*exp(-Z_rate(Isex,Iplatoon,Nage,Jsize))+Ntemp(Isex,Nage-1,Iplatoon,Jsize)*exp(-Z_rate(Isex,Iplatoon,Nage-1,Jsize)))*TransX(Isex,Iplatoon,Isize,Jsize);

      // Copy back  
      for (Iage=0;Iage<=Nage;Iage++)
       for (Isize=1;Isize<=NsizeDym;Isize++)
        Ntemp(Isex,Iage,Iplatoon,Isize) = NNext(Iage,Isize);
     }
   
   // SSB
   SSBTemp = 0;
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iage=1;Iage<=Nage;Iage++)
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      if (Model_Type==AGEMODEL)
       SSBTemp +=  Ntemp(Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage);   
      else
       for (Isize=1;Isize<=NsizeDym;Isize++)
        SSBTemp +=  Ntemp(Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage);   
   SSB(Iyear+1) = SSBTemp;       

   // Recruitment
   Numer = 4.0*Steepness*R0*SSBTemp/SSB(FirstProj_Yr);
   Denom = (1.0-Steepness)+(5*Steepness-1)*SSBTemp/SSB(FirstProj_Yr);
   Recruit = Numer/Denom;
   Recruit = R1Scalar*Recruit*RecDevMult(Iyear+1); 
   
   // sex ratio
   SexRatio(1) = 2.0/(1+exp(MeanSexRatio+SexRatioDev(Iyear+1)));
   if (Nsex == 2) SexRatio(2) = 2.0 - SexRatio(1);
   RecOut(Iyear+1,1) = Recruit*SexRatio(1);
   if (Nsex == 2) RecOut(Iyear+1,2) = Recruit*SexRatio(2);

   // Add in recruitment   
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Isize=1;Isize<=NsizeDym;Isize++)  
      Ntemp(Isex,0,Iplatoon,Isize) = RecLen(Isex,Iplatoon,Isize)*Recruit*SexRatio(Isex);
  }
  
 // Copy back
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Isize=1;Isize<=NsizeDym;Isize++)
     N(First_Year,Isex,Iage,Iplatoon,Isize) = Ntemp(Isex,Iage,Iplatoon,Isize); 
 
 // Compute the initial catch    
 InitCatch.initialize();
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Isize=1;Isize<=NsizeDym;Isize++)
      {
       BioAtLength = N(First_Year,Isex,Iage,Iplatoon,Isize)*selretwght2(Ifleet,Isex,Iplatoon,Iage,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
       InitCatch(Ifleet) += InitHrate(Ifleet)*BioAtLength;
      }   
    
 // Compute biomasses
 YearPass = First_Year;
 for (FleetPass=NfishFleet+1;FleetPass<=Nfleet;FleetPass++)
  Create_Survey_Pred();
 
// ================================================================================================================
FUNCTION Solve_HR
 // Apply the Hybrid method to solve for F by fleet

 int f,Isex,Iage,Iplatoon,Isize,tune_F;
 dvariable vbio,temp,temp1,join1,Z_adjuster2,Z_adjuster,BioAtLength,sum1,Surv;

 // Get initial Hrate estimate
 for (f=1;f<=NfishFleet;f++)
  if (Catch(f,Year) > 0)
   {
    vbio = 0;
    for (Isex=1;Isex<=Nsex;Isex++)
     for (Iage=0;Iage<=Nage;Iage++)
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
       if (Catch_Specs(f,3) == 1)
        vbio += sum(elem_prod(N(Year,Isex,Iage,Iplatoon),selretwghtUse(f,Isex,Iplatoon,Iage)));     
       else
        vbio += sum(elem_prod(N(Year,Isex,Iage,Iplatoon),selwghtUse(f,Isex,Iplatoon,Iage)));     
    temp = Catch(f,Year)/(vbio + Catch(f,Year));    
    join1=1.0/(1.0+mfexp(30.*(temp-0.95)));
    temp1=join1*temp + (1.0-join1)*0.95;
    Hrate(f,Year) = -log(1.-temp1);
   }
  else
    if (Catch(f,Year) == 0) Hrate(f,Year) = 0;
  //if (DiagLevel > 4) cout << "Hrate calculated" << endl;
  //if (DiagLevel > 4) cout <<  Hrate << endl;
 
 // Tune 
 for (tune_F=1;tune_F<=F_tune;tune_F++)
  {
   // Compute Z given F and M
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Iage=0;Iage<=Nage;Iage++)
      for (Isize=1;Isize<=NsizeDym;Isize++)
       {
        Z_rate(Isex,Iplatoon,Iage,Isize) = NatM(Isex,Year,Iplatoon,Iage,Isize);
        for (f=1;f<=NfishFleet;f++)
         Z_rate(Isex,Iplatoon,Iage,Isize) += selexS(f,Isex,Iplatoon,Iage,Isize) * Hrate(f,Year);
        Z_rate2(Isex,Iplatoon,Iage,Isize) = (1-mfexp(-Z_rate(Isex,Iplatoon,Iage,Isize)))/Z_rate(Isex,Iplatoon,Iage,Isize);  
       }
      
   // Now tune
   if (tune_F < F_tune)
    {
     Z_adjuster2 = 0;
     for (f=1;f<=NfishFleet;f++)
      if (Catch(f,Year) > 0)
       {
        for (Isex=1;Isex<=Nsex;Isex++)
         for (Iage=0;Iage<=Nage;Iage++)
          for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
           for (Isize=1;Isize<=NsizeDym;Isize++)
            if (Catch_Specs(f,3) == 1)
             Z_adjuster2 += Hrate(f,Year)*N(Year,Isex,Iage,Iplatoon,Isize)*selretwghtUse(f,Isex,Iplatoon,Iage,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
            else 
             Z_adjuster2 += Hrate(f,Year)*N(Year,Isex,Iage,Iplatoon,Isize)*selwghtUse(f,Isex,Iplatoon,Iage,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
       }
     Z_adjuster = TotalCatch(Year)/(Z_adjuster2+0.0001);
     
     // Adjust total Z
     for (Isex=1;Isex<=Nsex;Isex++)
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
       for (Iage=0;Iage<=Nage;Iage++)
        for (Isize=1;Isize<=NsizeDym;Isize++)
         {
          Z_rate(Isex,Iplatoon,Iage,Isize)  = NatM(Isex,Year,Iplatoon,Iage,Isize) + Z_adjuster*(Z_rate(Isex,Iplatoon,Iage,Isize)-NatM(Isex,Year,Iplatoon,Iage,Isize));
          Z_rate2(Isex,Iplatoon,Iage,Isize) = (1-mfexp(-Z_rate(Isex,Iplatoon,Iage,Isize)))/Z_rate(Isex,Iplatoon,Iage,Isize);  
         }
     
     // Adjust total exploitable biomass
     for (f=1;f<=NfishFleet;f++)
      if (Catch(f,Year) > 0)
       {
        Z_adjuster2 = 0;
        for (Isex=1;Isex<=Nsex;Isex++)
         for (Iage=0;Iage<=Nage;Iage++)
          for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
           for (Isize=1;Isize<=NsizeDym;Isize++)
            if (Catch_Specs(f,3) == 1)
             Z_adjuster2 += N(Year,Isex,Iage,Iplatoon,Isize)*selretwghtUse(f,Isex,Iplatoon,Iage,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
            else 
             Z_adjuster2 += N(Year,Isex,Iage,Iplatoon,Isize)*selwghtUse(f,Isex,Iplatoon,Iage,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
        temp = Catch(f,Year)/(Z_adjuster2 + 0.00001);    
        join1=1.0/(1.0+mfexp(30.*(temp-0.95*max_harvest_rate)));
        Hrate(f,Year)=join1*temp + (1.0-join1)*max_harvest_rate;
        //if (DiagLevel > 4) cout << "Hrate Z-Adjusted" << endl;
        //if (DiagLevel > 4) cout <<  Hrate << endl; 
       }
    }
   else
    {
     for (f=1;f<=NfishFleet;f++)
      {
       for (Isex=1;Isex<=Nsex;Isex++)
        for (Iage=0;Iage<=Nage;Iage++)
         for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
          for (Isize=1;Isize<=NsizeDym;Isize++)
           {
            // Note selretwght will not be biomass when Catch_Specs=1
            BioAtLength = N(Year,Isex,Iage,Iplatoon,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
            if (Catch_Specs(f,3) == 1)
             {
              CatchPred(f,Year) += Hrate(f,Year)*BioAtLength*selretwghtUse(f,Isex,Iplatoon,Iage,Isize);
              CatchPredN(f,Year) += Hrate(f,Year)*BioAtLength*selretwght1(f,Isex,Iplatoon,Iage,Isize);
              CatchPredB(f,Year) += Hrate(f,Year)*BioAtLength*selretwght2(f,Isex,Iplatoon,Iage,Isize);
              MidExpNum(f,Year) += BioAtLength*selretwght1(f,Isex,Iplatoon,Iage,Isize);
              MidExpBio(f,Year) += BioAtLength*selretwght2(f,Isex,Iplatoon,Iage,Isize);
             } 
            else
             {
              CatchPred(f,Year) += Hrate(f,Year)*BioAtLength*selwghtUse(f,Isex,Iplatoon,Iage,Isize);
              CatchPredN(f,Year) += Hrate(f,Year)*BioAtLength*selwght1(f,Isex,Iplatoon,Iage,Isize);
              CatchPredB(f,Year) += Hrate(f,Year)*BioAtLength*selwght2(f,Isex,Iplatoon,Iage,Isize);
              MidExpNum(f,Year) += BioAtLength*selwght1(f,Isex,Iplatoon,Iage,Isize);
              MidExpBio(f,Year) += BioAtLength*selwght2(f,Isex,Iplatoon,Iage,Isize);
             } 
           } 
       }   
    }
  } // Tune
  
  //Other fleets
  for (f=NfishFleet+1;f<=Nfleet;f++)
   {
    for (Isex=1;Isex<=Nsex;Isex++)
     for (Iage=0;Iage<=Nage;Iage++)
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
       if (Model_Type == AGEMODEL)
        {
         Surv = mfexp(-1*SurveyTime(f)*Z_rate(Isex,Iplatoon,Iage,1));
         sum1 = 0;
         for (Isize=1;Isize<=NsizeAge;Isize++)
          sum1 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(f,Isex,Isize);
         MidExpNum(f,Year) += N(Year,Isex,Iage,Iplatoon,1)*selexA(f,Isex,Iage)*sum1*Surv;
         sum1 = 0;
         for (Isize=1;Isize<=NsizeAge;Isize++)
          sum1 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(f,Isex,Isize)*wghtL(Isex,Isize);
         MidExpBio(f,Year) += N(Year,Isex,Iage,Iplatoon,1)*selexA(f,Isex,Iage)*sum1*Surv;
        }
       else
        {
         for (Isize=1;Isize<=NsizeDym;Isize++)
          {
           Surv = mfexp(-1*SurveyTime(f)*Z_rate(Isex,Iplatoon,Iage,Isize));
           MidExpNum(f,Year) += N(Year,Isex,Iage,Iplatoon,Isize)*selexL(f,Isex,Isize)*selexA(f,Isex,Iage)*Surv;
          } 
         for (Isize=1;Isize<=NsizeDym;Isize++)
          {
           Surv = mfexp(-1*SurveyTime(f)*Z_rate(Isex,Iplatoon,Iage,Isize));
           MidExpBio(f,Year) += N(Year,Isex,Iage,Iplatoon,Isize)*selexL(f,Isex,Isize)*selexA(f,Isex,Iage)*Surv*wghtL(Isex,Isize);
          } 
        }
    }   
  

// ================================================================================================================

FUNCTION Update_Dynamics
 // Annual update routine 
 int Isex,Iage,Iplatoon,Isize,Jsize,Ifleet;
 dvariable Test,Numer,Denom,Recruit,CatchAtLength,Temp;
 dvar_vector SexRatio(1,Nsex);

 // First compute the harvest rates
 Solve_HR();

 // Find the mid-year SSB 
 Create_Mid_SSB();
   
 // Remove catch and natural mortality (record the catch - both retained and total)
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=0;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equations D1.b1 & D.1b2
     {
      Ntemp(Isex,Iage,Iplatoon,1) = N(Year,Isex,Iage,Iplatoon,1)*exp(-Z_rate(Isex,Iplatoon,Iage,1));
      for (Isize=1;Isize<=NsizeAge;Isize++)
       for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
        {
         Temp = selexA(Ifleet,Isex,Iage)*PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize);
         CatchAtLength = Hrate(Ifleet,Year)*Temp*N(Year,Isex,Iage,Iplatoon,1)*Z_rate2(Isex,Iplatoon,Iage,1);
         CatchPredAgeSize(1,Ifleet,Year,Isex,Iage,Isize) += retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage)*CatchAtLength;
         CatchPredAgeSize(2,Ifleet,Year,Isex,Iage,Isize) += CatchAtLength;
         CatchPredAgeSize(3,Ifleet,Year,Isex,Iage,Isize) += (1.0-retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage))*CatchAtLength;;
        }
     }
    else                                                                     // Equations D1.a1 & D.1a2  
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       {
        Ntemp(Isex,Iage,Iplatoon,Isize) = N(Year,Isex,Iage,Iplatoon,Isize)*exp(-Z_rate(Isex,Iplatoon,Iage,Isize));
        for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
         {
          Temp = selexA(Ifleet,Isex,Iage)*selexL(Ifleet,Isex,Isize);
          CatchAtLength = Hrate(Ifleet,Year)*Temp*N(Year,Isex,Iage,Iplatoon,Isize)*Z_rate2(Isex,Iplatoon,Iage,Isize);
          CatchPredAgeSize(1,Ifleet,Year,Isex,Iage,Isize) += retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage)*CatchAtLength;
          CatchPredAgeSize(2,Ifleet,Year,Isex,Iage,Isize) += CatchAtLength;
          CatchPredAgeSize(3,Ifleet,Year,Isex,Iage,Isize) += (1.0-retainL(Ifleet,Isex,Isize)*retainA(Ifleet,Isex,Iage))*CatchAtLength;
         } 
       }  
     } 
  
 // Multiply by transition matrix
 for (Iage=1;Iage<=Nage-1;Iage++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equation A1.b                                   
     {
      N(Year+1,Isex,Iage,Iplatoon,1) = Ntemp(Isex,Iage-1,Iplatoon,1);
     }
    else                                                                     // Equation A1.a
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
       N(Year+1,Isex,Iage,Iplatoon,Isize) += Ntemp(Isex,Iage-1,Iplatoon,Jsize)*TransX(Isex,Iplatoon,Isize,Jsize);
     }
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equation A1.b
     {
      N(Year+1,Isex,Nage,Iplatoon,1) = (Ntemp(Isex,Nage,Iplatoon,1)+Ntemp(Isex,Nage-1,Iplatoon,1));
     }
    else                                                                     // Equation A1.a  
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
        N(Year+1,Isex,Nage,Iplatoon,Isize) += (Ntemp(Isex,Nage,Iplatoon,Jsize)+Ntemp(Isex,Nage-1,Iplatoon,Jsize))*TransX(Isex,Iplatoon,Isize,Jsize);
     }
   
 // SSB                                                                        Equation A.5
 SSB(Year+1) = 0;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=1;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
      SSB(Year+1) +=  N(Year+1,Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage);   
    else
     for (Isize=1;Isize<=NsizeDym;Isize++)
      SSB(Year+1) +=  N(Year+1,Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage);   
 //CheckFile << "SSB " << Year+1 << " " << SSB(Year+1) << endl;   
    
 // Recruitment                                                                Equation A.4
 Numer = 4.0*Steepness*R0*SSB(Year+1)/SSB(FirstProj_Yr);
 Denom = (1.0-Steepness)+(5*Steepness-1)*SSB(Year+1)/SSB(FirstProj_Yr);
 Recruit = Numer/Denom;
 Recruit = Recruit*RecDevMult(Year+1); 

 // sex ratio
 SexRatio(1) = 2.0/(1+exp(MeanSexRatio+SexRatioDev(Year+1)));
 if (Nsex == 2) SexRatio(2) = 2.0 - SexRatio(1);
 
 RecOut(Year+1,1) = Recruit*SexRatio(1);
 if (Nsex == 2) RecOut(Year+1,2) = Recruit*SexRatio(2);

 // Add in recruitment   
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Isize=1;Isize<=NsizeDym;Isize++)  
    N(Year+1,Isex,0,Iplatoon,Isize) = RecLen(Isex,Iplatoon,Isize)*Recruit*SexRatio(Isex);

 // Create biomass
 YearPass = Year;
 for (FleetPass=1;FleetPass<=NfishFleet;FleetPass++) Create_Survey_Pred();
 for (FleetPass=NfishFleet+1;FleetPass<=Nfleet;FleetPass++)
  {
   if (SurveyTime(FleetPass) > 0)YearPass = Year; else YearPass = Year + 1;
   Create_Survey_Pred();
  } 

// ================================================================================================================

FUNCTION Create_Mid_SSB
 int Iage,Isex,Isize,Iplatoon;

 SSBMid(Year) = 0;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=1;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
     {
      SSBMid(Year) +=  N(Year,Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,1));   
     }
    else
     { 
      for (Isize=1;Isize<=NsizeDym;Isize++)
       SSBMid(Year) +=  N(Year,Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,Isize));   
     }  
  
// ================================================================================================================
// ============================ Selectivity and retention =========================================================
// ================================================================================================================

FUNCTION SetLengthSelexA   
 int Ifleet,Isex,Iage,Ioffset1,Ioffset2,Ipnt,Iblock,II;
 int IcntA,JcntA,IparCnt,III;
 dvariable L50, Slope, Asympt, L50U, SlopeU, AsymptU; ;
 dvar_matrix Pars(1,100,1,6);
 dvar_vector sp(1,6),sp2(1,6);
 dvariable peak,upselex,downselex,final,point1,point2,peak2;
 dvariable join1,join2,t1,t2;
 dvariable t1min,t2min;
 dvariable asc,dsc;
 dvariable binwidth2;
 ivector Itest(1,10), Jtest(1,10);
  
 PriorVAPars.initialize();
 BaseSelexA.initialize();
 
 Ioffset1 = 0; Ipnt = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     dvector x(1,nknotsSelexA(Ifleet));
     dvar_vector VyI(1,nknotsSelexA(Ifleet));
     dvar_vector VyIU(1,nknotsSelexA(Ifleet));
     dvar_vector VyIU2(1,nknotsSelexA(Ifleet));
     IparCnt = SelexASpex(Ifleet,Isex,4);

     // Extract parameters 
     for (IcntA=1;IcntA<=IparCnt;IcntA++)  
      {
       Ioffset1 += 1;
       Jtest(IcntA) = Ioffset1;
       Pars(1,IcntA) = SelexAPars(Ioffset1);
       PriorVAPars(Ioffset1) = GenPrior(SelexAPars(Ioffset1), SelexAReal(Ioffset1,5), SelexAReal(Ioffset1,6), SelexAInt(Ioffset1,6));
       Ioffset2 = Ioffset1;
       JcntA = 1;
       if (SelexAInt(Ioffset2,2) > 0)
        for (III=1;III<=BlocksCnt(SelexAInt(Ioffset2,2));III++)
         {
          Ioffset1 += 1; JcntA += 1;
          Pars(JcntA,IcntA) = SelexAPars(Ioffset1);
          if (SelexAReal(Ioffset2,4) > 0) PriorVAPars(Ioffset1) += square(Pars(JcntA,IcntA))/(2.0*square(SelexAReal(Ioffset2,4)));
         }
      } 
     
     // Flat selectivity
     if (SelexASpex(Ifleet,Isex,1) == 0)   
      {
       Ipnt += 1;
       if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
        for (Iage=0;Iage<=Nage;Iage++) BaseSelexA(Ipnt,Iage) = 1;
       else 
        BaseSelexA(Ipnt).initialize();
       }
     
     // Logistic selectivity 
     if (SelexASpex(Ifleet,Isex,1) == 1)   
      {
       Ipnt += 1;
       if (Isex == 1) { L50 = Pars(1,1); Slope = Pars(1,2);}
       if (Isex == 2 & SelexParStyle==1) { L50 *= exp(Pars(1,1)); Slope *= exp(Pars(1,2)); }
       if (Isex == 2 & SelexParStyle==2) { L50 = Pars(1,1); Slope = Pars(1,2); }
       L50U = L50; SlopeU = Slope;
      
       if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
        BaseSelexA(Ipnt) = logistic(Ages,  SlopeU,  L50U);
       else 
        BaseSelexA(Ipnt).initialize();
     
       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexAInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           if (Itest(1)==1) L50U = L50; else L50U = L50*exp(Pars(Itest(1),1)); 
           if (Itest(2)==1) SlopeU = Slope; else SlopeU = Slope*exp(Pars(Itest(2),2)); 
           if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
            BaseSelexA(Ipnt) = logistic(Ages,  SlopeU,  L50U);
           else 
            BaseSelexA(Ipnt).initialize();
          }
        }

      } 
    
     // Double logistic selectivity 
     if (SelexASpex(Ifleet,Isex,1) == 2)   
      {
       Ipnt += 1;
       for (II=1;II<=6;II++)
        {
         if (Isex == 1) sp(II) = Pars(1,II);
         if (Isex == 2 & SelexParStyle==1) sp(II) = sp(II) + Pars(1,II);
         if (Isex == 2 & SelexParStyle==2) sp(II) = Pars(1,II);
        } 
       sp2 = sp;
       if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
        BaseSelexA(Ipnt) = doublelogistic(Ages,  sp2, 1.0);
       else 
        BaseSelexA(Ipnt).initialize();

       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexAInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           for (II=1;II<=6;II++) if (Itest(II)==1) sp2(II) = sp(II); else sp2(II) = sp(II)*exp(Pars(Itest(II),II));
           if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
            BaseSelexA(Ipnt) = doublelogistic(Ages,  sp2, 1.0);
           else 
            BaseSelexA(Ipnt).initialize();
          }
        }
      }
      
     // Spline selectivity 
     if (SelexASpex(Ifleet,Isex,1) == 3)   
      {
       Ipnt += 1;
       for (II = 1; II<=nknotsSelexA(Ifleet);II++) 
        {
         x(II) = Ages(xvalsSelexA(Ifleet,II));
         if (Isex==1) VyI(II) = SelexAPars(Ioffset1+II); 
         if (Isex==2) VyI(II) = VyI(II)+ SelexAPars(Ioffset1+II); 
        } 
       VyIU = elem_div(exp(VyI),1.0+exp(VyI));
       if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
        BaseSelexA(Ipnt) = splineAEP(x, VyIU,Ages);
       else 
        BaseSelexA(Ipnt).initialize();
       Ioffset1 += nknotsSelexA(Ifleet);

       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexAInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           for (II = 1; II<=nknotsSelexA(Ifleet);II++) if (Itest(II)==1) VyIU(II) = VyI(II); else VyIU(II) = VyI(II) + Pars(Itest(II),2); 
           if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
            BaseSelexA(Ipnt) = splineAEP(x, VyIU2,Ages);
           else 
            BaseSelexA(Ipnt).initialize();
          }
        }
      }  

     // Logistic selectivity with asymptote
     if (SelexASpex(Ifleet,Isex,1) == 4)   
      {
       Ipnt += 1;
       if (Isex == 1) { L50 = Pars(1,1); Slope = Pars(1,2); Asympt = Pars(1,3);}
       if (Isex == 2 & SelexParStyle==1) { L50 *= exp(Pars(1,1)); Slope *= exp(Pars(1,2)); Asympt *= exp(Pars(1,3));}
       if (Isex == 2 & SelexParStyle==2) { L50 = Pars(1,1); Slope = Pars(1,2); Asympt = Pars(1,3); }
       L50U = L50; SlopeU = Slope; AsymptU = Asympt;
      
       if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
        BaseSelexA(Ipnt) = logistic3(Ages,  SlopeU,  L50U, AsymptU);
       else 
        BaseSelexA(Ipnt).initialize();
     
       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexAInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexAInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           if (Itest(1)==1) L50U = L50; else L50U = L50*exp(Pars(Itest(1),1)); 
           if (Itest(2)==1) SlopeU = Slope; else SlopeU = Slope*exp(Pars(Itest(2),2)); 
           if (Itest(3)==1) AsymptU = Asympt; else AsymptU = Asympt*exp(Pars(Itest(3),3)); 
           if (SelexASpex(Ifleet,Isex,5) == Isex | SelexASpex(Ifleet,Isex,5) == 0) 
            BaseSelexA(Ipnt) = logistic3(Ages,  SlopeU,  L50U, AsymptU);
           else 
            BaseSelexA(Ipnt).initialize();
          }
        }

      } 
      
    } // Sex
     
   } // Fleet
 
// ================================================================================================================
FUNCTION SetLengthSelexL   
 int Ifleet,Isex,Isize,Ioffset1,Ioffset2,Ipnt,Iblock,II;
 int IcntA,JcntA,IparCnt,III;
 dvariable L50, Slope, Asympt, L50U, SlopeU, AsympU;
 dvar_matrix Pars(1,100,1,6);
 dvar_vector sp(1,6),sp2(1,6);
 dvariable peak,upselex,downselex,final,point1,point2,peak2;
 dvariable join1,join2,t1,t2;
 dvariable t1min,t2min;
 dvariable asc,dsc;
 dvariable binwidth2;
 ivector Itest(1,10), Jtest(1,10);
 
 PriorVLPars.initialize();
 BaseSelexL.initialize();
 
 Ioffset1 = 0; Ipnt = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     dvector x(1,nknotsSelexL(Ifleet));
     dvar_vector VyI(1,nknotsSelexL(Ifleet));
     dvar_vector VyIU(1,nknotsSelexL(Ifleet));
     dvar_vector VyIU2(1,nknotsSelexL(Ifleet));
     IparCnt = SelexLSpex(Ifleet,Isex,4);

     // Extract parameters 
     for (IcntA=1;IcntA<=IparCnt;IcntA++)  
      {
       Ioffset1 += 1;
       Jtest(IcntA) = Ioffset1;
       Pars(1,IcntA) = SelexLPars(Ioffset1);
       //cout << Ifleet << " " << IcntA << " " << Ioffset1 << " " << Pars(1,IcntA) << endl;
       PriorVLPars(Ioffset1) = GenPrior(SelexLPars(Ioffset1), SelexLReal(Ioffset1,5), SelexLReal(Ioffset1,6), SelexLInt(Ioffset1,6));
       Ioffset2 = Ioffset1;
       JcntA = 1;
       if (SelexLInt(Ioffset2,2) > 0)
        for (III=1;III<=BlocksCnt(SelexLInt(Ioffset2,2));III++)
         {
          Ioffset1 += 1; JcntA += 1;
          Pars(JcntA,IcntA) = SelexLPars(Ioffset1);
          //cout << "O" << Ifleet << " " << JcntA << " " << Ioffset1 << " " << Pars(JcntA,IcntA) << endl;
          if (SelexLReal(Ioffset2,4) > 0) PriorVLPars(Ioffset1) += square(Pars(JcntA,IcntA))/(2.0*square(SelexLReal(Ioffset2,4)));
         }
      } 
     
     // Flat selectivity
     if (SelexLSpex(Ifleet,Isex,1) == 0)   
      {
       Ipnt += 1;
       if (SelexLSpex(Ifleet,Isex,5) == Isex | SelexLSpex(Ifleet,Isex,5) == 0) 
        for (Isize=1;Isize<=NsizeMax;Isize++) BaseSelexL(Ipnt,Isize) = 1;
       else 
        BaseSelexL(Ipnt).initialize();
       }
     
     // Logistic selectivity 
     if (SelexLSpex(Ifleet,Isex,1) == 1)   
      {
       Ipnt += 1;
       if (Isex == 1) { L50 = Pars(1,1); Slope = Pars(1,2);}
       if (Isex == 2 & SelexParStyle==1) { L50 *= exp(Pars(1,1)); Slope *= exp(Pars(1,2)); }
       if (Isex == 2 & SelexParStyle==2) { L50 = Pars(1,1); Slope = Pars(1,2); }
       L50U = L50; SlopeU = Slope;
      
       if (SelexLSpex(Ifleet,Isex,5) == Isex  | SelexLSpex(Ifleet,Isex,5) == 0) 
        BaseSelexL(Ipnt) = logistic(MidLen,  SlopeU,  L50U);
       else 
        BaseSelexL(Ipnt).initialize();
     
       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexLInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           if (Itest(1)==1) L50U = L50; else L50U = L50*exp(Pars(Itest(1),1)); 
           if (Itest(2)==1) SlopeU = Slope; else SlopeU = Slope*exp(Pars(Itest(2),2)); 
           if (SelexLSpex(Ifleet,Isex,5) == Isex | SelexLSpex(Ifleet,Isex,5) == 0) 
            BaseSelexL(Ipnt) = logistic(MidLen,  SlopeU,  L50U);
           else 
            BaseSelexL(Ipnt).initialize();
          }
        }
      } 
    
     // Double logistic selectivity 
     if (SelexLSpex(Ifleet,Isex,1) == 2)   
      {
       Ipnt += 1;
       for (II=1;II<=6;II++)
        {
         if (Isex == 1) sp(II) = Pars(1,II);
         if (Isex == 2 & SelexParStyle==1) sp(II) = sp(II) + Pars(1,II);
         if (Isex == 2 & SelexParStyle==2) sp(II) = Pars(1,II);
        } 
       sp2 = sp;
       if (SelexLSpex(Ifleet,Isex,5) == Isex  | SelexLSpex(Ifleet,Isex,5) == 0) 
        BaseSelexL(Ipnt) = doublelogistic(MidLen,  sp2, LengthInc);
       else 
        BaseSelexL(Ipnt).initialize();

       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexLInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           for (II=1;II<=6;II++) if (Itest(II)==1) sp2(II) = sp(II); else sp2(II) = sp(II)*exp(Pars(Itest(II),II));
           if (SelexLSpex(Ifleet,Isex,5) == Isex  | SelexLSpex(Ifleet,Isex,5) == 0) 
            BaseSelexL(Ipnt) = doublelogistic(MidLen,  sp2, LengthInc);
           else 
            BaseSelexL(Ipnt).initialize();
          }
        }
      }
      
     // Spline selectivity 
     if (SelexLSpex(Ifleet,Isex,1) == 3)   
      {
       Ipnt += 1;
       for (II = 1; II<=nknotsSelexL(Ifleet);II++) 
        {
         x(II) = MidLen(xvalsSelexL(Ifleet,II));
         if (Isex==1) VyI(II) = SelexLPars(Ioffset1+II); 
         if (Isex==2) VyI(II) = VyI(II)+ SelexLPars(Ioffset1+II); 
        } 
       VyIU = elem_div(exp(VyI),1.0+exp(VyI));
       if (SelexLSpex(Ifleet,Isex,5) == Isex | SelexLSpex(Ifleet,Isex,5) == 0) 
        BaseSelexL(Ipnt) = splineAEP(x, VyIU,MidLen);
       else 
        BaseSelexL(Ipnt).initialize();
       Ioffset1 += nknotsSelexL(Ifleet);

       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexLInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           for (II = 1; II<=nknotsSelexL(Ifleet);II++) if (Itest(II)==1) VyIU(II) = VyI(II); else VyIU(II) = VyI(II) + Pars(Itest(II),2); 
           if (SelexLSpex(Ifleet,Isex,5) == Isex | SelexLSpex(Ifleet,Isex,5) == 0) 
            BaseSelexL(Ipnt) = splineAEP(x, VyIU2,MidLen);
           else 
            BaseSelexL(Ipnt).initialize();
          }
        }
      }  
      
     // Logistic selectivity with offset
     if (SelexLSpex(Ifleet,Isex,1) == 4)   
      {
       Ipnt += 1;
       if (Isex == 1) { L50 = Pars(1,1); Slope = Pars(1,2); Asympt = Pars(1,3); }
       if (Isex == 2 & SelexParStyle==1) { L50 *= exp(Pars(1,1)); Slope *= exp(Pars(1,2)); Asympt *= exp(Pars(1,3)); }
       if (Isex == 2 & SelexParStyle==2) { L50 = Pars(1,1); Slope = Pars(1,2); Asympt = Pars(1,3); }
       L50U = L50; SlopeU = Slope; AsympU = Asympt;
      
       if (SelexLSpex(Ifleet,Isex,5) == Isex  | SelexLSpex(Ifleet,Isex,5) == 0) 
        BaseSelexL(Ipnt) = logistic3(MidLen,  SlopeU,  L50U, AsympU);
       else 
        BaseSelexL(Ipnt).initialize();
     
       for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),First_Year);
       for (IyearA=First_Year;IyearA<=Last_Year+1;IyearA++)
        {
         IFound = 0; 
         for (IcntA=1;IcntA<=IparCnt;IcntA++) if (Itest(IcntA)!= Blocks(SelexLInt(Jtest(IcntA),2),IyearA)) IFound = 1;
         if (IFound == 1) 
          {
           for (IcntA=1;IcntA<=IparCnt;IcntA++) Itest(IcntA) = Blocks(SelexLInt(Jtest(IcntA),2),IyearA);
           Ipnt += 1;
           if (Itest(1)==1) L50U = L50; else L50U = L50*exp(Pars(Itest(1),1)); 
           if (Itest(2)==1) SlopeU = Slope; else SlopeU = Slope*exp(Pars(Itest(2),2)); 
           if (Itest(3)==1) AsympU = Asympt; else AsympU = Asympt*exp(Pars(Itest(3),3)); 
           if (SelexLSpex(Ifleet,Isex,5) == Isex | SelexLSpex(Ifleet,Isex,5) == 0) 
            BaseSelexL(Ipnt) = logistic3(MidLen,  SlopeU,  L50U, AsympU);
           else 
            BaseSelexL(Ipnt).initialize();
          }
        }
      } 
      
    } // Sex
     
   } // Fleet

// ================================================================================================================

FUNCTION SetLengthRetainA
 int Ifleet,Isex,Iage,Ioffset1,Ioffset2,Iblock,Ipnt;
 dvariable L50, Slope, Asymp;
 dvariable L50U, SlopeU, AsympU;
 
 // All discarded / retained
 for (Iage=0;Iage<=Nage;Iage++) BaseRetainA(-1,Iage) = 0;
 for (Iage=0;Iage<=Nage;Iage++) BaseRetainA(0,Iage) = 1;
  
 Ioffset1 = 0; Ioffset2 = 0; Ipnt = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     
     // Flat selectivity
     if (SelexASpex(Ifleet,Isex,6) == 0)  {}
     
     // Logistic selectivity 
     if (SelexASpex(Ifleet,Isex,6) == 1)   
      {
       Ipnt += 1;
       // Logistic selectivity (parameters)
       L50 = RetainAPars(Ioffset1+1); Slope = RetainAPars(Ioffset1+2); Asymp = RetainAPars(Ioffset1+3);
       L50U = L50; SlopeU = Slope; AsympU = Asymp;
       Ioffset1 += 3; 
       BaseRetainA(Ipnt) = logistic2(Ages,  Slope,  L50, Asymp);

       if (SelexASpex(Ifleet,Isex,7) > 0)
        for (Iblock=1;Iblock<=BlocksCnt(SelexASpex(Ifleet,Isex,7));Iblock++)
         {
          Ipnt += 1;
          L50U = L50*exp(RetainABlkPars(Ioffset2+1)); 
          SlopeU = Slope*exp(RetainABlkPars(Ioffset2+2)); 
          AsympU = Asymp*exp(RetainABlkPars(Ioffset2+3)); 
          Ioffset2 += 3;
          BaseRetainA(Ipnt) = logistic2(Ages,  SlopeU,  L50U, AsympU);
         } 
      } 
    } // Sex
   } // Fleet

// ================================================================================================================

FUNCTION SetLengthRetainL
 int Ifleet,Isex,Isize,Ioffset1,Ioffset2,Iblock,Ipnt;
 dvariable L50, Slope, Asymp;
 dvariable L50U, SlopeU, AsympU;

 // All discarded / retained
 for (Isize=1;Isize<=NsizeMax;Isize++) BaseRetainL(-1,Isize) = 0;
 for (Isize=1;Isize<=NsizeMax;Isize++) BaseRetainL(0,Isize) = 1;
 
 Ioffset1 = 0; Ioffset2 = 0; Ipnt = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     
     // Flat selectivity
     if (SelexLSpex(Ifleet,Isex,6) == 0) {}  
     
     // Logistic selectivity 
     if (SelexLSpex(Ifleet,Isex,6) == 1)   
      {
       Ipnt += 1;
       // Logistic selectivity (parameters)
       L50 = RetainLPars(Ioffset1+1); Slope = RetainLPars(Ioffset1+2); Asymp = RetainLPars(Ioffset1+3);
       L50U = L50; SlopeU = Slope; AsympU = Asymp;
       Ioffset1 += 3; 
       BaseRetainL(Ipnt) = logistic2(MidLen,  Slope,  L50, Asymp);

       if (SelexLSpex(Ifleet,Isex,7) > 0)
        for (Iblock=1;Iblock<=BlocksCnt(SelexLSpex(Ifleet,Isex,7));Iblock++)
         {
          Ipnt += 1;
          L50U = L50*exp(RetainLBlkPars(Ioffset2+1)); 
          SlopeU = Slope*exp(RetainLBlkPars(Ioffset2+2)); 
          AsympU = Asymp*exp(RetainLBlkPars(Ioffset2+3)); 
          Ioffset2 += 3;
          BaseRetainL(Ipnt) = logistic2(MidLen,  SlopeU,  L50U, AsympU);
         } 
      } 
    } // Sex
   } // Fleet

// ================================================================================================================

FUNCTION Create_Survey_Pred
 int Ifleet,Iage,Isex,Isize,Ipnt,Iplatoon;
 dvariable NumbersAtAge,NumbersAtAgeLen,Surv;
 
 // Survey biomass
 Ifleet = FleetPass;                                            // Equation F.2
  {
   if (Ifleet <= NfishFleet)
    if (YearPass <= Last_Year) 
     if (IndexTiming(Ifleet,2) == 1)
      SurveyBio(Ifleet,YearPass) = MidExpNum(Ifleet,YearPass);
     else 
      SurveyBio(Ifleet,YearPass) = MidExpBio(Ifleet,YearPass);
    
   if (Ifleet > NfishFleet)
    {
     for (Isex=1;Isex<=Nsex;Isex++)
      {
       Ipnt = (Ifleet-1)*Nsex+Isex;
       for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
        for (Iage=0;Iage<=Nage;Iage++)
         if (Model_Type == AGEMODEL)
          {
           Surv = mfexp(-1*SurveyTime(Ifleet)*Z_rate(Isex,Iplatoon,Iage,1));
           NumbersAtAge = N(YearPass,Isex,Iage,Iplatoon,1)*BaseSelexA(SelexAPnt(Ipnt,YearPass),Iage)*Surv;
           for (Isize=1;Isize<=NsizeMax;Isize++)
            {
             NumbersAtAgeLen = NumbersAtAge*PhiGrow(Isex,Iplatoon,Iage,Isize)*BaseSelexL(SelexLPnt(Ipnt,YearPass),Isize);
             CatchPredAgeSize(2,Ifleet,YearPass,Isex,Iage,Isize) += NumbersAtAgeLen;
             if (IndexTiming(Ifleet,2) == 1) SurveyBio(Ifleet,YearPass) += NumbersAtAgeLen;
             if (IndexTiming(Ifleet,2) == 2) SurveyBio(Ifleet,YearPass) += NumbersAtAgeLen*wghtL(Isex,Isize);
            }
          }
         else
          {
           for (Isize=1;Isize<=NsizeDym;Isize++)
            {
             Surv = mfexp(-1*SurveyTime(Ifleet)*Z_rate(Isex,Iplatoon,Iage,Isize));
             NumbersAtAge = N(YearPass,Isex,Iage,Iplatoon,Isize)*BaseSelexA(SelexAPnt(Ipnt,YearPass),Iage)*Surv;
             NumbersAtAgeLen = NumbersAtAge*BaseSelexL(SelexLPnt(Ipnt,YearPass),Isize);
             CatchPredAgeSize(2,Ifleet,YearPass,Isex,Iage,Isize) += NumbersAtAgeLen;
             if (IndexTiming(Ifleet,2) == 1) SurveyBio(Ifleet,YearPass) += NumbersAtAgeLen;
             if (IndexTiming(Ifleet,2) == 2) SurveyBio(Ifleet,YearPass) += NumbersAtAgeLen*wghtL(Isex,Isize);
            }
          } // Model type 
      }  // Sex     
    } // If statement
  } 
  
// ====================================================================================================================
// ================================= Data Fitting Routines ============================================================
// ====================================================================================================================

FUNCTION CatchLikelihood
 int Ifleet,Iyear;
 dvariable error;

 CatchLikeCompFleet.initialize();
 CatchLikeCompData.initialize();

 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
   if (Catch(Ifleet,Iyear) > 0)
    {
     error = log(Catch(Ifleet,Iyear)/CatchPred(Ifleet,Iyear))/SigmaCatch; 
     CatchLikeCompData(Ifleet,Iyear) = 0.5*square(error);
     CatchLikeCompFleet(Ifleet) += square(error);
    } 
 
FUNCTION InitCLikelihood
 int Ifleet;
 dvariable error;
 
 LikeInitC.initialize();
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  if (InitialC(Ifleet) > 0)
   {
    error = log(InitialC(Ifleet)/InitCatch(Ifleet))/SigmaIntC(Ifleet);
    if (current_phase() >= PhaseInitC(Ifleet)) LikeInitC(Ifleet) = square(error);    
   }

FUNCTION DiscardLikelihood
 // Equations D.1 and D.2 as well as F.3 and F.4
 int Ifleet,Ipnt,Iyear,Isex,Iplatoon,Iage,Isize;
 dvariable error,NumbersAtAge,LikeComp,Residual,Discard;

 DiscardLikeCompFleet.initialize();
 DiscardLikeCompData.initialize();

 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=0;Isex<=Nsex;Isex++)
   {
    for (Ipnt=1;Ipnt<=NdiscardData;Ipnt++)
     if (DiscardFleet(Ipnt) == Ifleet  & DiscardSex(Ipnt)==Isex)
      {
       Iyear = DiscardYear(Ipnt);
       DiscardPred(Ipnt) = 0;
       for (Iage=0;Iage<=Nage;Iage++)
        for (Isize=1;Isize<=NsizeMax;Isize++)
         {
          if (Discard_type(Ifleet) == 1)
           {
            if (Discard_Specs(Ifleet,2) == 1)
             {
              if (Isex==0) { NumbersAtAge = CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize)+ CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize); }
              if (Isex==1) NumbersAtAge = CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize);
              if (Isex==2) NumbersAtAge = CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize);
             } 
            if (Discard_Specs(Ifleet,2) == 2)
             {
              if (Isex==0) { NumbersAtAge = wghtL(1,Isize)*CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize)+ wghtL(2,Isize)*CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize); }
              if (Isex==1) NumbersAtAge = wghtL(1,Isize)*CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize);
              if (Isex==2) NumbersAtAge = wghtL(2,Isize)*CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize);
             }
           }  
          if (Discard_type(Ifleet) == 2)
           {
            if (Discard_Specs(Ifleet,2) == 1)
             {
              NumbersAtAge = 0;
              if (Isex==1 || Isex==0) 
               {
                NumbersAtAge += CatchPredAgeSize(2,Ifleet,Iyear,1,Iage,Isize);
               } 
              if (Isex==2 || Isex==0) 
               {
                NumbersAtAge += CatchPredAgeSize(2,Ifleet,Iyear,2,Iage,Isize);
               } 
             } 
            if (Discard_Specs(Ifleet,2) == 2)
             {
              NumbersAtAge = 0;
              if (Isex==1 || Isex==0) 
               {
                NumbersAtAge += wghtL(1,Isex)*CatchPredAgeSize(2,Ifleet,Iyear,Isex,Iage,Isize);
               } 
              if (Isex==2 || Isex==0) 
               {
                NumbersAtAge += wghtL(2,Isex)*CatchPredAgeSize(2,Ifleet,Iyear,Isex,Iage,Isize);
               } 
             }
            
           }  
          if (NumbersAtAge <0) 
           {
            cout << Ifleet << " " << Isex << " " << Iage << " " << Isize << endl;
            cout << CatchPredAgeSize(2,Ifleet,Iyear,1,Iage,Isize)  << " " << CatchPredAgeSize(1,Ifleet,Iyear,1,Iage,Isize) << endl;
            if (Nsex > 1)
             cout << CatchPredAgeSize(2,Ifleet,Iyear,2,Iage,Isize)  << " " << CatchPredAgeSize(1,Ifleet,Iyear,2,Iage,Isize) << endl;
            cout << "CRAP" << endl;
            //exit(1);
           }
          // AEP FIX
          //DiscardPred(Ipnt) += Discard_Specs(Ifleet,5)*NumbersAtAge;
          DiscardPred(Ipnt) += NumbersAtAge;
         }
       Residual = log(DiscardData(Ipnt,1)/DiscardPred(Ipnt));
       LikeComp = log(DiscardData(Ipnt,2)) + 0.5*square(Residual)/square(DiscardData(Ipnt,2));
       DiscardLikeCompData(Ipnt) = LikeComp;
       DiscardLikeCompFleet(Ifleet) += LikeComp;
      }
   }
 

FUNCTION IndexLikelihood                                                     // Equations F.1 and F.2  
 int Ifleet,Ipnt,Iyear,Iblk,IblkPnt;
 dvariable Residual,LikeComp,SigmaUse2,Top,Bot,qestV;
 
 // Solve for Q
 qest.initialize();
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   if (QpriorType(Ifleet) == 0) 
    {
     for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++) qest(Ifleet,Iyear) = exp(QpriorSpex(Ifleet,1));
    }
   else
    {
     IblkPnt = QpriorBlock(Ifleet,2);
     for (Iblk=1;Iblk<=BlocksCnt(QpriorBlock(Ifleet,2))+1;Iblk++)
      {
       Top = 0; Bot = 0;
       if (QpriorType(Ifleet) == 1)
        { Top = QpriorSpex(Ifleet,1)/square(QpriorSpex(Ifleet,2)); Bot = 1.0/square(QpriorSpex(Ifleet,2)); }
       for (Ipnt=1;Ipnt<=NindexData;Ipnt++)
        if (IndexFleet(Ipnt) == Ifleet & Blocks(IblkPnt,IndexYear(Ipnt))==Iblk)
         {
          Iyear = IndexYear(Ipnt);
          SigmaUse2 = square(IndexData(Ipnt,2))+AddVarParUse(Ifleet);
          Residual = log(IndexData(Ipnt,1)/SurveyBio(Ifleet,Iyear));
          Top += Residual/SigmaUse2;
          Bot += 1.0/SigmaUse2;
         }
       if (Bot > 0)
        qestV = exp(Top/Bot); 
       for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++) 
        if (Blocks(QpriorBlock(Ifleet,2),Iyear)==Iblk)  qest(Ifleet,Iyear) = qestV;
      }  
    }  
  }  

 // Initialize likelihoods
 IndexLikeCompFleet.initialize();
 IndexLikeCompData.initialize();

 // Compute the contributions to the likelihood function  
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {  
   for (Ipnt=1;Ipnt<=NindexData;Ipnt++)
    if (IndexFleet(Ipnt) == Ifleet)
     {
      Iyear = IndexYear(Ipnt);
      IndexPred(1,Ipnt) = SurveyBio(Ifleet,Iyear);
      IndexPred(2,Ipnt) = qest(Ifleet,Iyear)*SurveyBio(Ifleet,Iyear);
      SigmaUse2 = square(IndexData(Ipnt,2))+AddVarParUse(Ifleet);
      //cout << "M " << Ifleet << " " << Iyear << " " << SurveyBio(Ifleet,Iyear) << " " << IndexPred(2,Ipnt) << endl;
      Residual = log(IndexData(Ipnt,1)/(qest(Ifleet,Iyear)*SurveyBio(Ifleet,Iyear)));
      LikeComp = 0.5*log(SigmaUse2) + 0.5*square(Residual)/SigmaUse2;
      IndexLikeCompData(Ipnt) = LikeComp;
      IndexLikeCompFleet(Ifleet) += LikeComp;
     }
  } 
   
FUNCTION EffortLikelihood                                                    // Equations F.5 and F.6
 int Ifleet,Ipnt,Iyear;
 dvariable Residual,LikeComp;
 dvar_vector npnts(1,Nfleet);
 
 // Solve for Q
 qestE.initialize();
 npnts.initialize();
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NeffortData;Ipnt++)
    if (EffortFleet(Ipnt) == Ifleet)
     {
      Iyear = EffortYear(Ipnt);
      npnts(Ifleet) += 1.0/square(EffortData(Ipnt,2));
      Residual = log(EffortData(Ipnt,1)/Hrate(Ifleet,Iyear));
      qestE(Ifleet) += Residual/square(EffortData(Ipnt,2));
     }
   if (npnts(Ifleet) > 0)  
    qestE(Ifleet) = exp(qestE(Ifleet)/npnts(Ifleet)); 
  }  

 // Initialize likelihoods
 EffortLikeCompFleet.initialize();
 EffortLikeCompData.initialize();

 // Compute the contributions to the likelihood function  
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NeffortData;Ipnt++)
    if (EffortFleet(Ipnt) == Ifleet)
     {
      Iyear = EffortYear(Ipnt);
      EffortPred(Ipnt) = qestE(Ifleet)*Hrate(Ifleet,Iyear);
      Residual = log(EffortData(Ipnt,1)/(qestE(Ifleet)*Hrate(Ifleet,Iyear)));
      LikeComp = log(EffortData(Ipnt,2)) + 0.5*square(Residual)/square(EffortData(Ipnt,2));
      EffortLikeCompData(Ipnt) = LikeComp;
      EffortLikeCompFleet(Ifleet) += LikeComp;
     }
  } 

FUNCTION LengthLikelihood                                                    // Equations F.7 and F.8
 int Ifleet,Ipnt,Isex,Iage,Iyear,Isize,Itype;
 dvariable error,totalA,LikeComp,Temp1,Sigma2,ExtraPen;
 dvar_vector temp(1,2*NsizeMax);

 // Initialize likelihoods
 LengthLikeCompFleet.initialize();
 LengthLikeCompData.initialize();

 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NlengthData;Ipnt++)
    if (LengthDataFleet(Ipnt) == Ifleet)
     {
      ExtraPen = 0;
      Isex = LengthDataSex(Ipnt);
      Itype = LengthDataType(Ipnt);
      Iyear = LengthDataYear(Ipnt);
    
      // Find and store the predictions (0=total;1=discard;2=landed)
      temp.initialize();
      if (Itype == 2)
       for (Isize=1;Isize<=NsizeMax;Isize++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Iage=0;Iage<=Nage;Iage++) temp(Isize) += CatchPredAgeSize(1,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Iage=0;Iage<=Nage;Iage++) temp(NsizeMax+Isize) += CatchPredAgeSize(1,Ifleet,Iyear,2,Iage,Isize);
        }  
      if (Itype == 1)
       for (Isize=1;Isize<=NsizeMax;Isize++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Iage=0;Iage<=Nage;Iage++) temp(Isize) += CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Iage=0;Iage<=Nage;Iage++) temp(NsizeMax+Isize) += CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize);
        } 
      if (Itype == 0)
       for (Isize=1;Isize<=NsizeMax;Isize++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Iage=0;Iage<=Nage;Iage++) temp(Isize) += CatchPredAgeSize(2,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Iage=0;Iage<=Nage;Iage++) temp(NsizeMax+Isize) += CatchPredAgeSize(2,Ifleet,Iyear,2,Iage,Isize);
        } 
      if (LenMinMax(Ipnt,1) < LenMinMax(Ipnt,2))
       {
        for (Isize=1;Isize<=LenMinMax(Ipnt,1)-1;Isize++)
         { temp(LenMinMax(Ipnt,1)) += temp(Isize); temp(Isize) = 0; }
        for (Isize=NsizeMax;Isize>LenMinMax(Ipnt,2);Isize--)
         { temp(LenMinMax(Ipnt,2)) += temp(Isize); temp(Isize) = 0; }
       }
      if (LenMinMax(Ipnt,3) < LenMinMax(Ipnt,4)) 
       {
        for (Isize=NsizeMax+1;Isize<=LenMinMax(Ipnt,3)-1;Isize++)
         { temp(LenMinMax(Ipnt,3)) += temp(Isize); temp(Isize) = 0; }
        for (Isize=2*NsizeMax;Isize>LenMinMax(Ipnt,4);Isize--)
         { temp(LenMinMax(Ipnt,4)) += temp(Isize); temp(Isize) = 0; }
       }  
        
      temp = temp/sum(temp);         
      PredLengthComp(Ipnt) = temp; 
      
      // Likelihood function (multinomial)
      LikeComp = 0;
      if (LengthCompLike==1)
       {
        for (Isize=1;Isize<=Nsex*NsizeMax;Isize++)
         { 
          if (LengthData(Ipnt,Isize) > 0)
           {
            temp(Isize) = posfun(temp(Isize),1.0e-10,ExtraPen);
            error = LengthData(Ipnt,Isize) * log((0.0e-20+temp(Isize))/LengthData(Ipnt,Isize));
            // Brute force
            //if (temp(Isize) > 0)
            // error = LengthData(Ipnt,Isize) * log((0.0e-20+temp(Isize))/LengthData(Ipnt,Isize));
            //else
            // error = -1000;
           }  
          else
           error = 0;
          LikeComp += LengthDataSS(Ipnt)*error; 
         } 
        LikeComp = -1*LikeComp; 
       } 
      // Likelihood function (robust normal)
      if (LengthCompLike==2)
       {
        for (Isize=1;Isize<=Nsex*NsizeMax;Isize++)
         { 
          Sigma2 = LengthData(Ipnt,Isize)*(1.0-LengthData(Ipnt,Isize));
          Sigma2 = (Sigma2 + 0.1/float(Nsex*NsizeMax))/LengthDataSS(Ipnt);
          Temp1 = log(mfexp(-1.0*square(LengthData(Ipnt,Isize)-temp(Isize))/(2*Sigma2))+0.01);
          LikeComp += 0.5*log(Sigma2)-Temp1;
         } 
       } 
      LengthLikeCompData(Ipnt) = LikeComp+ExtraPen;
      LengthLikeCompFleet(Ifleet) += LikeComp+ExtraPen;
     }

  }  

FUNCTION ConditionalLikelihood                                               // Equations F.9 and F.10
 int Ifleet,Ipnt,Jpnt,Isex,Iage,Iyear,Isize,Itype,Len1,Len2;
 int MaxAgeLike;
 dvariable error,totalA,LikeComp,Temp1,Sigma2;
 dvar_vector temp(1,2*(Nage+1));

 // Initialize likelihoods
 AgeLengthLikeCompFleet.initialize();
 AgeLengthLikeCompData.initialize();

 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NagelengthData;Ipnt++)
    if (AgeLengthDataFleet(Ipnt) == Ifleet)
     {
      Isex = AgeLengthDataSex(Ipnt);
      Itype = AgeLengthDataType(Ipnt);
      Iyear = AgeLengthDataYear(Ipnt);
      Len1 = AgeLengthDataLen1(Ipnt); 
      Len2 = AgeLengthDataLen2(Ipnt); 
    
      // Find and store the predictions (0=total;1=discard;2=landed)
      temp.initialize();
      if (Itype == 2)
       for (Iage=0;Iage<=Nage;Iage++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Iage+1) += CatchPredAgeSize(1,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Nage+2+Iage) += CatchPredAgeSize(1,Ifleet,Iyear,2,Iage,Isize);
        }  
      if (Itype == 1)
       for (Iage=0;Iage<=Nage;Iage++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Iage+1) += CatchPredAgeSize(3,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Nage+2+Iage) += CatchPredAgeSize(3,Ifleet,Iyear,2,Iage,Isize);
        } 
      if (Itype == 0)
       for (Iage=0;Iage<=Nage;Iage++)
        {
         if (Isex == 1 | Isex == 0 | Isex == 3)
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Iage+1) += CatchPredAgeSize(2,Ifleet,Iyear,1,Iage,Isize);
         if (Nsex == 2 & (Isex == 2 | Isex == 0))
          for (Isize=Len1;Isize<=Len2;Isize++) temp(Nage+2+Iage) += CatchPredAgeSize(2,Ifleet,Iyear,2,Iage,Isize);
        } 
      temp = temp/sum(temp);         
      for (Jpnt=1;Jpnt<=(Nsex*Nage+Nsex);Jpnt++) PredAgeLengthComp(Ipnt,Jpnt) = temp(Jpnt); 
      
      // Likelihood function (multinomial)
      LikeComp = 0;
      if (Nsex==1) MaxAgeLike = Nage+1; else MaxAgeLike = 2*Nage+2;
      if (AgeCompLike==1)
       {
        for (Iage=1;Iage<=MaxAgeLike;Iage++)
         { 
          if (AgeLengthData(Ipnt,Iage) > 0)
           error = AgeLengthData(Ipnt,Iage) * log(temp(Iage)/AgeLengthData(Ipnt,Iage));
          else
           error = 0;
          LikeComp += AgeLengthDataSS(Ipnt)*error; 
         } 
        LikeComp = -1*LikeComp; 
       } 
      // Likelihood function (robust normal)
      if (AgeCompLike==2)
       {
        for (Iage=1;Iage<=MaxAgeLike;Iage++)
         { 
          Sigma2 = AgeLengthData(Ipnt,Iage)*(1.0-AgeLengthData(Ipnt,Iage));
          Sigma2 = (Sigma2 + 0.1/float(Nsex*(Nage+1)))/AgeLengthDataSS(Ipnt);
          Temp1 = log(mfexp(-1.0*square(AgeLengthData(Ipnt,Iage)-temp(Iage))/(2*Sigma2))+0.01);
          LikeComp += 0.5*log(Sigma2)-Temp1;
         } 
       } 
      AgeLengthLikeCompData(Ipnt) = LikeComp;
      AgeLengthLikeCompFleet(Ifleet) += LikeComp;
     }
  }  

FUNCTION MeanSizeLikelihood
 int Ifleet,Ipnt,Iqnt,Isex,Iage,Iyear,Isize,Iplatoon,Offset;
 dvariable MeanL,SDL,Denom,SEM,LikeComp,Resid;
 dvar_vector Weights(1,NsizeMax);
 
 MeanSizeLikeCompFleet.initialize();
 MeanSizeLikeCompData.initialize();

 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NmeansizeData;Ipnt++)
    if (MeanSizeDataFleet(Ipnt) == Ifleet)
     {
      Iyear = MeanSizeDataYear(Ipnt);
      for (Isex=1;Isex<=Nsex;Isex++)
       {
        //cout << Ifleet << " " << Isex << " " << Iyear << " " << Ipnt << endl;
        Iqnt = (Ifleet-1)*Nsex+Isex;
        if (Model_Type==AGEMODEL)
         {
          Offset = (Isex-1)*(Nage+1)+1;
          for (Iage=0;Iage<=Nage;Iage++)
           {
            // Find the weights
            Weights.initialize();  MeanL = 0; Denom = 0;
            for (Isize=1;Isize<=NsizeMax;Isize++)
             {
              for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
               Weights(Isize) += N(Iyear,Isex,Iage,Iplatoon,1)*PhiGrow(Isex,Iplatoon,Iage,Isize)*BaseSelexL(SelexLPnt(Iqnt,Iyear),Isize);
              MeanL += MidLen(Isize)*Weights(Isize); Denom += Weights(Isize);
             }
            MeanL = MeanL / Denom;
            //cout << Isex << " " << Iage << " " << MeanL << endl;
            SDL = 1.0e-20;
            for (Isize=1;Isize<=NsizeMax;Isize++)
             SDL += square(MidLen(Isize)-MeanL)*Weights(Isize)/Denom;
            SDL = sqrt(SDL);
            PredMeanSize(Ipnt,Offset+Iage) = MeanL;
            PredMeanSizeSQ(Ipnt,Offset+Iage) = SDL;
           } 
         }
        else
         {
          Offset = (Isex-1)*(Nage+1)+1;
          for (Iage=0;Iage<=Nage;Iage++)
           {
            // Find the weights
            Weights.initialize();  MeanL = 0; Denom = 0;
            for (Isize=1;Isize<=NsizeMax;Isize++)
             {
              for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
               Weights(Isize) += N(Iyear,Isex,Iage,Iplatoon,Isize)*BaseSelexL(SelexLPnt(Iqnt,Iyear),Isize);
              MeanL += MidLen(Isize)*Weights(Isize); Denom += Weights(Isize);
             }
            MeanL = MeanL / Denom;
            //cout << Isex << " " << Iage << " " << MeanL << endl;
            SDL = 1.0e-20;
            for (Isize=1;Isize<=NsizeMax;Isize++)
             SDL += square(MidLen(Isize)-MeanL)*Weights(Isize)/Denom;
            SDL = sqrt(SDL);
            PredMeanSize(Ipnt,Offset+Iage) = MeanL;
            PredMeanSizeSQ(Ipnt,Offset+Iage) = SDL;
           }
         } // model type
         
        for (Iage=0;Iage<=Nage;Iage++) 
         if (MeanSizeData(Ipnt,Offset+Iage) > 0 && MeanSizeDataSS(Ipnt,Offset+Iage) > 0) 
          {
           SEM = PredMeanSizeSQ(Ipnt,Offset+Iage)/sqrt(MeanSizeDataSS(Ipnt,Offset+Iage));
           Resid = (PredMeanSize(Ipnt,Offset+Iage) -MeanSizeData(Ipnt,Offset+Iage))/SEM;
           LikeComp = log(SEM) + 0.5*square(Resid);
           MeanSizeLikeCompData(Ipnt) += LikeComp;
           MeanSizeLikeCompFleet(Ifleet) += LikeComp;
          }
       } // Sex
     } // if  
    } // Fleet  

FUNCTION TaggingLikelihood
 int Isex,Iplatoon,Isize,Jsize,Ksize,Iyear,Ipnt,Iage;
 int RelSize,RecSize,TimeAtLib,Itag;
 dvariable MeanLenRec,Total,Total1,Slope,Intc,Pred,VarLenRec,Cnt;
 dvariable Bot,Top,MeanPred,MeanObs,VarPred,EffN;
 dvar_vector Preds(1,NsizeMax);
 dvar_vector Obs(1,NsizeMax);
 
 TransX2.initialize();
 if (NtagData > 0)
  {
   if (Model_Type==AGEMODEL)
    {
    }
   else
    {
     for (Isex=1;Isex<=Nsex;Isex++)
      {
      // Set up (one-year-ahead)
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
       {
        for (Isize=1;Isize<=NsizeMax;Isize++)
         for (Jsize=1;Jsize<=NsizeMax;Jsize++)
         {
//                       To    From        
          TempX(Iplatoon,Isize,Jsize) = TransX(Isex,Iplatoon,Isize,Jsize);
          TransX2(Isex,1,Iplatoon,Isize,Jsize) = TempX(Iplatoon,Isize,Jsize);
         } 
       }   
      
      // Project ahead
      for (Iyear=2;Iyear<=6;Iyear++)
       {
        TempY.initialize();
        for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
         {
          for (Isize=1;Isize<=NsizeMax;Isize++)
           for (Jsize=1;Jsize<=NsizeMax;Jsize++)
            for (Ksize=1;Ksize<=NsizeMax;Ksize++)
//                          To    From                            To    From                  To    From
            TempY(Iplatoon,Isize,Jsize) += TransX(Isex,Iplatoon,Isize,Ksize)*TempX(Iplatoon,Ksize,Jsize);
         }   
        TempX = TempY;  
        for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
         for (Isize=1;Isize<=NsizeMax;Isize++)
          for (Jsize=1;Jsize<=NsizeMax;Jsize++)
           TransX2(Isex,Iyear,Iplatoon,Isize,Jsize) = TempX(Iplatoon,Isize,Jsize);
       }
     }
   }
  }
  
 // Find the proportions
 if (NtagData > 0)
  {
   TagProps.initialize();
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Isize=1;Isize<=NsizeMax;Isize++)
     {
      Total = 0;
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
       {      
        for (Iage=1;Iage<=Nage;Iage++) TagProps(Isex,Isize,Iplatoon) += N(TagYear,Isex,Iage,Iplatoon,Isize);
        Total += TagProps(Isex,Isize,Iplatoon);
       }  
      for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++) TagProps(Isex,Isize,Iplatoon) /= Total;
     }
   }  

 // Likelihood
 TagLike = 0;
 for (Itag=1;Itag<=NtagData;Itag++)
  {
   Isex = TagData(Itag,1);
   RelSize = TagData(Itag,2);
   RecSize = TagData(Itag,3);
   TimeAtLib = TagData(Itag,4);
   Total = 0;
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    for (Isize=1;Isize<=NsizeMax;Isize++)
     Total += TagProps(Isex,RelSize,Iplatoon)*selexL(1,Isex,Isize)*TransX2(Isex,TimeAtLib,Iplatoon,Isize,RelSize);
   Pred = 0; 
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    Pred += TagProps(Isex,RelSize,Iplatoon)*selexL(1,Isex,RecSize)*TransX2(Isex,TimeAtLib,Iplatoon,RecSize,RelSize)/Total+0.00001;
   //cout << Itag << " " << Pred << " " << TagOffset(Itag) << endl;
   TagLike = TagLike + -1.0*float(TagData(Itag,5))*log(Pred/TagOffset(Itag));
  }
  
 if (DoTagDiag == 1)
  {
   PredTagTable.initialize();
   PredProb.initialize();
   for (Itag=1;Itag<=NtagData;Itag++)
    {
     Isex = TagData(Itag,1);
     RelSize = TagData(Itag,2);
     RecSize = TagData(Itag,3);
     TimeAtLib = TagData(Itag,4);
     Preds.initialize();
     Total = 0;
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      for (Isize=RelSize;Isize<=NsizeMax;Isize++)
       Total += TagProps(Isex,RelSize,Iplatoon)*selexL(1,Isex,Isize)*TransX2(Isex,TimeAtLib,Iplatoon,Isize,RelSize);
     for (Isize=RelSize;Isize<=NsizeMax;Isize++)
      {
       Preds(Isize) = 0;
       for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
        Preds(Isize) += TagProps(Isex,RelSize,Iplatoon)*selexL(1,Isex,Isize)*TransX2(Isex,TimeAtLib,Iplatoon,Isize,RelSize)/Total;
       PredProb(Itag,Isize) = Preds(Isize);
      } 
      
     PredTag(Itag,0) = MidLen(RelSize);
     PredTag(Itag,1) = Preds(RecSize);
     PredTag(Itag,2) = MidLen(RecSize);
     MeanLenRec = 0;
     for (Isize=1;Isize<=NsizeMax;Isize++)
      MeanLenRec += MidLen(Isize)*Preds(Isize);
     PredTag(Itag,3) = MeanLenRec;
     PredTagTable(Isex,TimeAtLib,1,RecSize) += TagData(Itag,5);
     for (Isize=1;Isize<=NsizeMax;Isize++)
      PredTagTable(Isex,TimeAtLib,2,Isize) += TagData(Itag,5)*Preds(Isize);
    }   
   
   Ipnt = 0;
   for (TimeAtLib=1;TimeAtLib<=10;TimeAtLib++)
    {
     for (RelSize=1;RelSize<NsizeMax;RelSize++)
      {
       Preds.initialize();
       Obs.initialize();
       Cnt = 0;
       for (Itag=1;Itag<=NtagData;Itag++)
        if (TagData(Itag,4) == TimeAtLib & TagData(Itag,2) == RelSize)
         {
          Cnt += float(TagData(Itag,5));
          Preds = PredProb(Itag);
          Obs(TagData(Itag,3)) += float(TagData(Itag,5));
         }
       if (Cnt > 0)
        {
         Ipnt += 1;
         Bot = 0; Top = 0; MeanPred = 0; MeanObs = 0; VarPred = 0;
         for (RecSize=RelSize;RecSize<=NsizeMax;RecSize++)
          {
           Obs(RecSize) /= Cnt;
           Top += Preds(RecSize)*(1.0-Preds(RecSize));
           Bot += square(Obs(RecSize) - Preds(RecSize));
           MeanPred += Preds(RecSize)*MidLen(RecSize);
           MeanObs += Obs(RecSize)*MidLen(RecSize);
           VarPred += Preds(RecSize)*MidLen(RecSize)*MidLen(RecSize);
          }
         VarPred = (VarPred-MeanPred*MeanPred)/Cnt; 
         EffN = Top/Bot; 
         TagDiag2(Ipnt,1) = TimeAtLib;
         TagDiag2(Ipnt,2) = RelSize;
         TagDiag2(Ipnt,3) = MidLen(RelSize);
         TagDiag2(Ipnt,4) = Cnt;
         TagDiag2(Ipnt,5) = MeanPred;
         TagDiag2(Ipnt,6) = MeanObs;
         TagDiag2(Ipnt,7) = sqrt(VarPred);
         TagDiag2(Ipnt,8) = EffN/Cnt;
         TagDiag2(Ipnt,9) = (MeanObs-MeanPred)/sqrt(VarPred);
        } 
      }    
     }
    NTagDiag2 = Ipnt; 
  }

FUNCTION Penalties
 int Iyear,Ipar,f;
 dvariable Fbar,NF,SS;

 // Initialize penalties
 Penal.initialize();

 // Recruitment deviation penalty (major devs)
 for (Iyear=max(EarlyDevEnd,RecEstYr1);Iyear<=RecEstYr2;Iyear++) 
  Penal(1) += square(RecDev(Iyear)+SigPhi(Iyear)*square(SigmaR)/2.0)/(2*square(SigmaR));

 if (EarlyDevEnd > RecEstYr1)
  for (Iyear=min(EarlyDevEnd,RecEstYr1);Iyear<=max(EarlyDevEnd,RecEstYr1)-1;Iyear++) 
   Penal(1) += square(RecDev(Iyear)-RecDev(Iyear+1))/(2*square(SigmaRecDev));
  
 for (Ipar=1;Ipar<=NBiolPar;Ipar++)
  {
   //cout << BiolInt(Ipar,6) << endl;
   if (BiolInt(Ipar,6)>0)
    BiolPrior(Ipar) += GenPrior(BiolPars(Ipar), BiolReal(Ipar,5), BiolReal(Ipar,6), BiolInt(Ipar,6));
  }
 Penal(2) = sum(BiolPrior);
 Penal(4) = DevPenal; 
 Penal(5) =sum(PriorVLPars) + sum(PriorVAPars);

 if (TreatMissF==3) 
  {
   SS = 0;
   for (f=1;f<=NfishFleet;f++)
    {
     Fbar = 0; NF = 0;
     for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
      if (CatchFix(f,Iyear) > 0) { Fbar += Hrate(f,Iyear); NF += 1; } 
     Fbar /= NF; 
     for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
      if (CatchFix(f,Iyear) < 0) { SS += 1000000.0* square(Hrate(f,Iyear)-Fbar); } 
    }
   Penal(3) = SS;
  }

// ====================================================================================================================
// ====================================================================================================================
// ================================================================================================================

FUNCTION SolveF35
 dvariable fbar,SSBZero,fmin,fmax,Temp,CatchAtLength;
 dvariable ZZ;
 int Icnt,Iyear,f,Isex,Iage,Iplatoon,Isize;
 
 HrateProj.initialize();
 Project();
 SSBZero = SSBPass;
 cout << "SSBZ " << SSBZero << endl;
 cout << " HrateProj(1)"  << "SSBPass" << "SSBPass/SSBZero" << endl;
 
 for (f=2;f<=NfishFleet;f++)
  {
   fbar = 0;
   for (Iyear=Last_Year-4;Iyear<=Last_Year;Iyear++) fbar += Hrate(f,Iyear);
   HrateProj(f) = fbar/5.0;
  }

 fmin = 0;
 fmax = 10;
 for (Icnt=1;Icnt<=20;Icnt++)
  {
   HrateProj(1) = (fmin+fmax)/2.0;
   Project();
   cout << " " << HrateProj(1) << " " << SSBPass << " " << SSBPass/SSBZero << endl; 
   if (SSBPass > 0.35*SSBZero)
    fmin = HrateProj(1);
   else
    fmax = HrateProj(1);
  }
 cout << " " << HrateProj(1) << " " << SSBPass << " " << SSBPass/SSBZero << endl; 
  

 // Remove catch and natural mortality (record the catch - bot retained and total)
 OFL.initialize();
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=0;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equations D1.b1 & D.1b2
     {
      for (Isize=1;Isize<=NsizeAge;Isize++)
       for (f=1;f<=NfishFleet;f++)
        {
         Temp = selexA(f,Isex,Iage)*PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(f,Isex,Isize);
         ZZ = Z_rate(Isex,Iplatoon,Iage,1);
         CatchAtLength = HrateProj(f)*Temp*N(Last_Year+1,Isex,Iage,Iplatoon,1)*(1-exp(-ZZ))/ZZ;
         OFL(f) += wghtL(Isex,Isize)*retainL(f,Isex,Isize)*retainA(f,Isex,Iage)*CatchAtLength;
        }
     }
    else                                                                     // Equations D1.a1 & D.1a2  
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       {
        for (f=1;f<=NfishFleet;f++)
         {
          Temp = selexA(f,Isex,Iage)*selexL(f,Isex,Isize);
          ZZ = Z_rate(Isex,Iplatoon,Iage,Isize);
          CatchAtLength = HrateProj(f)*Temp*N(Last_Year+1,Isex,Iage,Iplatoon,Isize)*(1-exp(-ZZ))/ZZ;
          OFL(f) += wghtL(Isex,Isize)*retainL(f,Isex,Isize)*retainA(f,Isex,Iage)*CatchAtLength;
         } 
       }  
     } 
   cout << "OFL per Fishery Fleet" << OFL << endl;  
  
FUNCTION Project
 int Iyear,Isex,Iplatoon,Iage,Isize,f;

 // Selex is equal to that for the last year
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Iage=0;Iage<=Nage;Iage++)
    for (Isize=1;Isize<=NsizeDym;Isize++)
     {
      NtempA(Isex,Iage,Iplatoon,Isize) = 100.0;
      Z_rate(Isex,Iplatoon,Iage,Isize) = NatM(Isex,Last_Year,Iplatoon,Iage,Isize);
      for (f=1;f<=NfishFleet;f++)
       Z_rate(Isex,Iplatoon,Iage,Isize) += selexS(f,Isex,Iplatoon,Iage,Isize) * HrateProj(f);
     }

 // Constant recruitment
 Proj_SR = 0;

 for (Iyear=1;Iyear<=100;Iyear++) Project_update();

FUNCTION Project_update
 // Annual update routine 
 int Isex,Iage,Iplatoon,Isize,Jsize,Ifleet;
 dvariable Test,Numer,Denom,Recruit,Temp,SSBO,SSBP;
 dvar_vector SexRatio(1,Nsex);

 // Set the sex-ratio
 SexRatio(1) = 2.0/(1+exp(MeanSexRatio));
 if (Nsex == 2) SexRatio(2) = 2.0 - SexRatio(1);

 // Remove catch and natural mortality (record the catch - bot retained and total)
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=0;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equations D1.b1 & D.1b2
     {
      Ntemp(Isex,Iage,Iplatoon,1) = NtempA(Isex,Iage,Iplatoon,1)*exp(-Z_rate(Isex,Iplatoon,Iage,1));
     }
    else                                                                     // Equations D1.a1 & D.1a2  
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       Ntemp(Isex,Iage,Iplatoon,Isize) = NtempA(Isex,Iage,Iplatoon,Isize)*exp(-Z_rate(Isex,Iplatoon,Iage,Isize));
     } 

 // SSB (Mid-year)                                                                        Equation X.X
 SSBP = 0;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=1;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
      SSBP +=  NtempA(Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,1));   
    else
     for (Isize=1;Isize<=NsizeDym;Isize++)
      SSBP +=  NtempA(Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage)*exp(-phi*Z_rate(Isex,Iplatoon,Iage,Isize));;   
  
 // Multiply by transition matrix
 NtempA.initialize();
 for (Iage=1;Iage<=Nage-1;Iage++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equation A1.b                                   
     {
      NtempA(Isex,Iage,Iplatoon,1) = Ntemp(Isex,Iage-1,Iplatoon,1);
     }
    else                                                                     // Equation A1.a
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
       NtempA(Isex,Iage,Iplatoon,Isize) += Ntemp(Isex,Iage-1,Iplatoon,Jsize)*TransX(Isex,Iplatoon,Isize,Jsize);
     }
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type == AGEMODEL)                                              // Equation A1.b
     {
      NtempA(Isex,Nage,Iplatoon,1) = (Ntemp(Isex,Nage,Iplatoon,1)+Ntemp(Isex,Nage-1,Iplatoon,1));
     }
    else                                                                     // Equation A1.a  
     {
      for (Isize=1;Isize<=NsizeDym;Isize++)
       for (Jsize=GrowthLow(Isize);Jsize<=Isize;Jsize++)
        NtempA(Isex,Nage,Iplatoon,Isize) += (Ntemp(Isex,Nage,Iplatoon,Jsize)+Ntemp(Isex,Nage-1,Iplatoon,Jsize))*TransX(Isex,Iplatoon,Isize,Jsize);
     }
   
 // SSB                                                                        Equation A.5
 SSBO = 0;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iage=1;Iage<=Nage;Iage++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    if (Model_Type==AGEMODEL)
      SSBO +=  NtempA(Isex,Iage,Iplatoon,1)*fecAge(Isex,Iplatoon,Iage);   
    else
     for (Isize=1;Isize<=NsizeDym;Isize++)
      SSBO +=  NtempA(Isex,Iage,Iplatoon,Isize)*fec(Isex,Isize,Iage);   
 //CheckFile << "SSB " << SSBO << endl;   
    
 // Recruitment Equation A.4
 if (Proj_SR==1)
  {
   Numer = 4.0*Steepness*R0*SSB(Year+1)/SSB(FirstProj_Yr);
   Denom = (1.0-Steepness)+(5*Steepness-1)*SSB(Year+1)/SSB(FirstProj_Yr);
   Recruit = Numer/Denom;
   Recruit = Recruit*RecDevMult(Year+1); 
  }
 else
  Recruit = 100.0;

 // Add in recruitment   
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   for (Isize=1;Isize<=NsizeDym;Isize++)  
    NtempA(Isex,0,Iplatoon,Isize) = RecLen(Isex,Iplatoon,Isize)*Recruit*SexRatio(Isex);

 SSBPass = SSBP;

// ================================================================================================================
// ================================================================================================================

FUNCTION CreateOutput

 int Isex,Iplatoon,Ilen,Jlen,Iyear,Ifleet,Ipnt,Iage,Isize,IcntB;
 int Npar, NparEst, Ipar, Iselex, Iblk, Itag, Nout, Itype,Ipnt2;
 int TimeAtLib,Jsize;
 dvariable temp,sum1,sum2,Top,Bot,McAllister;
 dvariable ObsLen, PredLen, VarPredLen, Residual, VarN, VarO;
 dvariable MeanLenRec,Term1,Term2,Total;
 dvar_matrix MeanL(1,Nsex,0,MaxAge),SDL(1,Nsex,0,MaxAge);
 dvar_vector Prob1(1,NsizeMax),Prob2(1,NsizeMax);
 dvar_vector temp1(1,NsizeMax);
 dvar_matrix FrancisOut(1,NlengthData,1,4);
 
 cout << "Calling Output" << endl;
 SolveF35();
 cout << "Done F35%" << endl;
 
 OutFile1.close();
 OutFile1.open("OutputFile1.Out");


 // Likelihood components
 // =====================

 OutFile1 << "Penalties" << endl;
 OutFile1 << "Recruitment penalty : " << Penal(1) << endl;
 OutFile1 << "Biological parameters priors penalty : " << Penal(2) << endl;
 OutFile1 << "Fishing mortality penalties : " << Penal(3) << endl;
 OutFile1 << "Devs penalty : " << Penal(4) << endl;
 OutFile1 << "Selectivity parameters prior : " << Penal(5) << endl;
 OutFile1 << "Total penalties : " << sum(Penal) << endl;
 OutFile1 << "Likelihood components" << endl; 
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << "Initial catch " << Ifleet << " : " << fleet_names[Ifleet] << " : " << LikeInitC(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << "Catch components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << CatchLikeCompFleet(Ifleet) << " " << Catch_Lambda(Ifleet) << " " << CatchLikeCompFleet(Ifleet)*Catch_Lambda(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << "Discard components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << DiscardLikeCompFleet(Ifleet)  << " " << Discard_Lambda(Ifleet) << " " << DiscardLikeCompFleet(Ifleet)*Discard_Lambda(Ifleet) << endl;  
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << "Index components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << IndexLikeCompFleet(Ifleet)  << " " << Index_Lambda(Ifleet) << " " << IndexLikeCompFleet(Ifleet)*Index_Lambda(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << "Effort components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << EffortLikeCompFleet(Ifleet)  << " " << Effort_Lambda(Ifleet) << " " << EffortLikeCompFleet(Ifleet)*Effort_Lambda(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << "Length components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << LengthLikeCompFleet(Ifleet)  << " " << Length_Lambda(Ifleet) << " " << LengthLikeCompFleet(Ifleet)*Length_Lambda(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << "AgeLength components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << AgeLengthLikeCompFleet(Ifleet)  << " " << AgeLength_Lambda(Ifleet) << " " << AgeLengthLikeCompFleet(Ifleet)*AgeLength_Lambda(Ifleet) << endl; 
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << "MeanSize components " << Ifleet << " : " << fleet_names[Ifleet] << " : " << MeanSizeLikeCompFleet(Ifleet)  << " " << MeanSize_Lambda(Ifleet) << " " << MeanSizeLikeCompFleet(Ifleet)*MeanSize_Lambda(Ifleet) << endl; 
 OutFile1 << "Tagging : " << TagLike << " " << Tag_Lambda << " " << TagLike*Tag_Lambda << endl;
 OutFile1 << "Overall objective function " << Obj << " Phase/function call " << current_phase() << " " << IFuncCallCount << endl;
 OutFile1 << endl;

 OutFile1 << endl << "Structure (#sex, #platoons, #age #length)" << endl;
 OutFile1 << Nsex << " " << Nplatoon << " " << Nage << " " << NsizeMax << " " << Nfleet << endl;
 OutFile1 << FirstProj_Yr << " " << First_Year << " " << Last_Year << endl;
 OutFile1 << MidLen << endl;
 
 OutFile1 << endl << "OFL outputs" << endl;
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << HrateProj(Ifleet) << " " << OFL(Ifleet) << endl;
 OutFile1 << endl; 


 // Estimated parameters
 // ====================
 Npar = 0; NparEst = 0;
 for (Ipar=1;Ipar<=NBiolPar;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : Biological Par " << Ipar << " : " << BiolPars(Ipar) << " ";
   OutFile1 << BiolPrior(Ipar) << " ";
   if (BiolPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(BiolPars(Ipar),BiolParMin(Ipar),BiolParMax(Ipar)); } 
   OutFile1 << endl;
  }

 for (Ipar=1;Ipar<=NSRPar;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : Stock_recruit Par " << Ipar << " : " << SRPars(Ipar) << " ";
   if (SRPhases(Ipar) > 0) { NparEst +=1; OutFile1  << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(SRPars(Ipar),SRParMin(Ipar),SRParMax(Ipar)); } 
   OutFile1 << endl;
  }

 Ipar = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    if (SelexASpex(Ifleet,Isex,1) == 0) IparCnt = 0;                                 // Flat Selex
    if (SelexASpex(Ifleet,Isex,1) == 1) IparCnt = 2;                                 // Logistic
    if (SelexASpex(Ifleet,Isex,1) == 2) IparCnt = 6;                                 // Double normal 
    if (SelexASpex(Ifleet,Isex,1) == 3) IparCnt = SelexASpex(Ifleet,Isex,2);         // Spline
    if (SelexASpex(Ifleet,Isex,1) == 4) IparCnt = 3;                                 // Logistic with fixed asymptote
    if (SelexASpex(Ifleet,Isex,1) == 5) IparCnt = 0;                                 // Mirrored selectivity
    if (IparCnt > 0)
     for (IcntA=1;IcntA<=IparCnt;IcntA++)  
      {
       Npar += 1;  Ipar += 1;
       OutFile1 << "AgeSel_" << fleet_names[Ifleet] << "_" << Isex << "_A" << IcntA << "_" << SelexASpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << SelexAPars(Ipar) << " ";
       OutFile1 << PriorVAPars(Ipar) << " ";
       if (SelexAPhases(Ipar) > 0) { NparEst +=1;  OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(SelexAPars(Ipar),SelexAParMin(Ipar),SelexAParMax(Ipar)); } 
       OutFile1 << endl;
       IcntB = Ipar;
       if (SelexAInt(IcntB,2) > 0)
        for (III=1;III<=BlocksCnt(SelexAInt(IcntB,2));III++)
         {
          Npar += 1; Ipar += 1;
          OutFile1 << "AgeSel_" << fleet_names[Ifleet] << "_" << Isex << "_B" << IcntA << "_" << III << "_"<< SelexASpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << SelexAPars(Ipar) << " ";
          OutFile1 << PriorVAPars(Ipar) << " ";
          if (SelexAPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(SelexAPars(Ipar),SelexAParMin(Ipar),SelexAParMax(Ipar));  } 
          OutFile1 << endl;
         }
      } 
   }    

 Ipar = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    if (SelexLSpex(Ifleet,Isex,1) == 0) IparCnt = 0;                                 // Flat Selex
    if (SelexLSpex(Ifleet,Isex,1) == 1) IparCnt = 2;                                 // Logistic
    if (SelexLSpex(Ifleet,Isex,1) == 2) IparCnt = 6;                                 // Double normal 
    if (SelexLSpex(Ifleet,Isex,1) == 3) IparCnt = SelexLSpex(Ifleet,Isex,2);         // Spline
    if (SelexLSpex(Ifleet,Isex,1) == 4) IparCnt = 3;                                 // Logistic with fixed asymptote
    if (SelexLSpex(Ifleet,Isex,1) == 5) IparCnt = 0;                                 // Mirrored selectivity
    if (IparCnt > 0)
     for (IcntA=1;IcntA<=IparCnt;IcntA++)  
      {
       Npar +=1;
       Ipar += 1;
       OutFile1 << "SizeSel_" << fleet_names[Ifleet] << "_" << Isex << "_A" << IcntA << "_" <<SelexLSpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << SelexLPars(Ipar) << " ";
       OutFile1 << PriorVLPars(Ipar) << " ";
       if (SelexLPhases(Ipar) > 0) { NparEst +=1;   OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(SelexLPars(Ipar),SelexLParMin(Ipar),SelexLParMax(Ipar));} 
       OutFile1 << endl;
       IcntB = Ipar;
       if (SelexLInt(IcntB,2) > 0)
        for (III=1;III<=BlocksCnt(SelexLInt(IcntB,2));III++)
         {
          Npar +=1;
          Ipar += 1;
          OutFile1 << "SizeSel_" << fleet_names[Ifleet] << "_" << Isex << "_B" << IcntA << "_" << III << "_"<< SelexLSpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << SelexLPars(Ipar) << " ";
          OutFile1 << PriorVLPars(Ipar) << " ";
          if (SelexLPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(SelexLPars(Ipar),SelexLParMin(Ipar),SelexLParMax(Ipar)); } 
          OutFile1 << endl;
         }
      } 
   }    

 IcntA = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    IparCnt = 0;
    if (SelexASpex(Ifleet,Isex,6) == 0) IparCnt = 0;                                 // Uniform
    if (SelexASpex(Ifleet,Isex,6) == 1) IparCnt = 3;                                 // Logistic
    if (IparCnt > 0)
     for (Ipar=IcntA+1;Ipar<=IcntA+IparCnt;Ipar++)
      {
       Npar +=1; 
       OutFile1 << "AgeRet_" << fleet_names[Ifleet] << "_" << Isex << "_" << SelexASpex(Ifleet,Isex,6) << " " << Ipar-IcntA << " " << Ipar << " " << RetainAPars(Ipar) << " ";
       if (RetainAPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(RetainAPars(Ipar),RetainAParMin(Ipar),RetainAParMax(Ipar)); } 
       OutFile1 << endl;
      } 
    IcntA += IparCnt;  
   }    

 IcntA = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    if (SelexASpex(Ifleet,Isex,6) != 5)
     {
      if (SelexASpex(Ifleet,Isex,7) > 0)
       for (Iblk=1;Iblk<=BlocksCnt(SelexASpex(Ifleet,Isex,7));Iblk++)
        {
         if (SelexASpex(Ifleet,Isex,6) == 1) IparCnt = 3;                                 // Logistic
         if (IparCnt > 0)
          for (Ipar=IcntA+1;Ipar<=IcntA+IparCnt;Ipar++)
           {
            Npar +=1; 
            OutFile1 << "AgeRetBlk_" << fleet_names[Ifleet] << "_" << Isex << "_" << SelexASpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << RetainABlkPars(Ipar) << " ";
            if (RetainABlkPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(RetainABlkPars(Ipar),RetainABlkParMin(Ipar),RetainABlkParMax(Ipar)); } 
            OutFile1 << endl;
           } 
         IcntA += IparCnt;  
        }
     }
   }  

 IcntA = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    IparCnt = 0;
    if (SelexLSpex(Ifleet,Isex,6) == 0) IparCnt = 0;                                 // Uniform
    if (SelexLSpex(Ifleet,Isex,6) == 1) IparCnt = 3;                                 // Logistic
    if (IparCnt > 0)
     for (Ipar=IcntA+1;Ipar<=IcntA+IparCnt;Ipar++)
      {
       Npar +=1; 
       OutFile1 << "LenRet_" << fleet_names[Ifleet] << "_" << Isex << "_" << SelexLSpex(Ifleet,Isex,6) << " " << Ipar-IcntA << " " << Ipar << " " << RetainLPars(Ipar) << " ";
       if (RetainLPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(RetainLPars(Ipar),RetainLParMin(Ipar),RetainLParMax(Ipar));} 
       OutFile1 << endl;
      } 
    IcntA += IparCnt;  
   }    

 IcntA = 0; Npar = 0;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    if (SelexLSpex(Ifleet,Isex,6) != 5)
     {
      if (SelexLSpex(Ifleet,Isex,7) > 0)
       for (Iblk=1;Iblk<=BlocksCnt(SelexLSpex(Ifleet,Isex,7));Iblk++)
        {
         if (SelexLSpex(Ifleet,Isex,6) == 1) IparCnt = 3;                                 // Logistic
         if (IparCnt > 0)
          for (Ipar=IcntA+1;Ipar<=IcntA+IparCnt;Ipar++)
           {
            Npar +=1; 
            OutFile1 << "LenRetBlk_" << fleet_names[Ifleet] << "_" << Isex << "_" << SelexLSpex(Ifleet,Isex,1) << " " << Ipar-IcntA << " " << Ipar << " " << RetainLBlkPars(Ipar) << " ";
            if (RetainLBlkPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst << " "; CheckBounds(RetainLBlkPars(Ipar),RetainLBlkParMin(Ipar),RetainLBlkParMax(Ipar)); } 
            OutFile1 << endl;
           } 
         IcntA += IparCnt;  
        }
     }
   }  


 for (Ipar=1;Ipar<=NestCatch;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : MissinglogF " << Ipar << " : " << MissingLogF(Ipar) << " ";
   if (PhaseMissingF > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst; } 
   OutFile1 << endl;
  }

 for (Ipar=1;Ipar<=NfishFleet;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : MissingLogq Par " << Ipar << " : " << MissingLogq(Ipar) << " ";
   if (FleetMissF(Ipar)> 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst; } 
   OutFile1 << endl;
  }

 for (Ipar=RecEstYr1;Ipar<=RecEstYr2;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : RecDev " << Ipar << " : " << RecDev(Ipar) << " ";
   { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst; } 
   OutFile1 << endl;
  }

 for (Ipar=1;Ipar<=NfishFleet;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : Initial_Fs " << Ipar << " : " << InitHrate(Ipar) << " ";
   if (InitHratePhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst; } 
   OutFile1 << endl;
  }

 // Additional variance
 for (Ipar=1;Ipar<=Nfleet;Ipar++)
  {
   Npar +=1; 
   OutFile1 << Npar << " : Additional variance " << Ipar << " : " << AddVarPar(Ipar) << " ";
   if (AddVarParPhases(Ipar) > 0) { NparEst +=1; OutFile1 << ParsOut.sd(NparEst) << " " << NparEst; } 
   OutFile1 << endl;
  }

 OutFile1 << endl;
 OutFile1 << "SSB" << endl;
 for (Iyear=FirstProj_Yr;Iyear<=Last_Year+1;Iyear++)
  if (Iyear <= Last_Year)
   OutFile1 << Iyear << " " << SSB(Iyear) << " " << SSBOut.sd(Iyear) << " " << SSBMid(Iyear) << endl;
  else 
   OutFile1 << Iyear << " " << SSB(Iyear) << " " << SSBOut.sd(Iyear) << endl;
   

 OutFile1 << endl;
 OutFile1 << "Recruitment" << endl;
 for (Iyear=FirstProj_Yr;Iyear<=Last_Year+1;Iyear++)
  {
   if (Nsex==1)
    OutFile1 << Iyear << " " << RecOut(Iyear,1) << " " << RecOut.sd(Iyear,1);
   else 
    OutFile1 << Iyear << " " << RecOut(Iyear,1) << " " << RecOut.sd(Iyear,1)<< " " << RecOut(Iyear,2) << " " << RecOut.sd(Iyear,2);
   OutFile1 << endl;
  } 

 // Fits to Data 
 // ============

 OutFile1 << endl << "Initial_catch" << endl;
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  OutFile1 << Ifleet << " " << InitialC(Ifleet) << " " << InitCatch(Ifleet) << " " << SigmaIntC(Ifleet) << " " << LikeInitC(Ifleet) << endl;
 
 OutFile1 << endl << "Mortality_and_catch" << endl;
 OutFile1 << NfishFleet*(Last_Year-First_Year+1) << endl;
 OutFile1 << "# Fleet Year Type Obs_Catch Harvest_rate Pred_Catch Pred_CatchN  Pred_CatchW Log_likelihood" << endl;
 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
   OutFile1 << Ifleet << " " << Iyear << " " << CatchType(Ifleet,Iyear) << " " << Catch(Ifleet,Iyear) << " " << Hrate(Ifleet,Iyear) << " " << CatchPred(Ifleet,Iyear) << " " << CatchPredN(Ifleet,Iyear) << " " << CatchPredB(Ifleet,Iyear) << " " << CatchLikeCompData(Ifleet,Iyear) << endl; 

 OutFile1 << endl << "Index_fit" << endl;
 OutFile1 <<  NindexData << endl;
 OutFile1 << "# Fleet Year Type Obs_Index Input_CV Total_CV Abundance Pred_Index log-likelihod" << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NindexData;Ipnt++)
    if (IndexFleet(Ipnt) == Ifleet)
     {
      Iyear = IndexYear(Ipnt);
      OutFile1 << Ifleet << " " << Iyear << " " << IndexData(Ipnt,1) << " " << IndexData(Ipnt,2) << " " << sqrt(square(IndexData(Ipnt,2))+AddVarParUse(Ifleet)) << " " << IndexPred(1,Ipnt) << " " << IndexPred(2,Ipnt) << " " << IndexLikeCompData(Ipnt) << endl;
     } 
  } 
 OutFile1 << endl << "Catchability_values" << endl;
 OutFile1 << qest << endl;
  

 OutFile1 << endl << "Discard_fit" << endl;
 OutFile1 <<  NdiscardData << endl;
 OutFile1 << "# Fleet Year Sex Obs_Discard Input_CV Pred_Discard log-likelihod" << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NdiscardData;Ipnt++)
    if (DiscardFleet(Ipnt) == Ifleet)
     {
      Iyear = DiscardYear(Ipnt);
      Isex = DiscardSex(Ipnt);
      OutFile1 << Ifleet << " " << Iyear << " " << Isex << " " << DiscardData(Ipnt,1) << " " << DiscardData(Ipnt,2) << " " << DiscardPred(Ipnt) << " " << DiscardLikeCompData(Ipnt) << endl;
     } 
  } 
  
 OutFile1 << endl << "Effort_fit" << endl;
 OutFile1 <<  NeffortData << endl;
 OutFile1 << "# Fleet Year S Obs_Effort Input_CV Pred_Effort log-likelihod" << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NeffortData;Ipnt++)
    if (EffortFleet(Ipnt) == Ifleet)
     {
      Iyear = EffortYear(Ipnt);
      OutFile1 << Ifleet << " " << Iyear << " " << EffortData(Ipnt,1) << " " << EffortData(Ipnt,2) << " " << EffortPred(Ipnt) << " " << EffortLikeCompData(Ipnt) << endl;
     } 
  } 

 OutFile1 << endl << "Length_fit" << endl;
 OutFile1 <<  NlengthData << endl;
 OutFile1 << "# Fleet Year Fleet Type Sex Eff_N" << " " << MidLen << endl;
 Ipnt2 = 0;
 for (Itype=0;Itype<=2;Itype++)
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
   {
    VarN = 0; Nout = 0;
    for (Ipnt=1;Ipnt<=NlengthData;Ipnt++)
     if (LengthData(Ipnt,-3) == Ifleet & LengthData(Ipnt,-1) == Itype)
      {
       Top = 0; Bot = 0; 
       ObsLen = 0; PredLen = 0; VarPredLen = 0;
       for (Ilen = 1;Ilen<=Nsex*NsizeMax;Ilen++)
        {
         Top += PredLengthComp(Ipnt,Ilen)*(1.0-PredLengthComp(Ipnt,Ilen));
         Bot += square(PredLengthComp(Ipnt,Ilen)- LengthData(Ipnt,Ilen));
         if (Ilen <= NsizeMax)
          {
           ObsLen += MidLen(Ilen) * LengthData(Ipnt,Ilen);
           PredLen += MidLen(Ilen) * PredLengthComp(Ipnt,Ilen);
           VarPredLen += square(MidLen(Ilen)) * PredLengthComp(Ipnt,Ilen);
          }
         else 
          {
           ObsLen += MidLen(Ilen-NsizeMax) * LengthData(Ipnt,Ilen);
           PredLen += MidLen(Ilen-NsizeMax) * PredLengthComp(Ipnt,Ilen);
           VarPredLen += square(MidLen(Ilen-NsizeMax)) * PredLengthComp(Ipnt,Ilen);
          }
        }
       if (VarPredLen - square(PredLen) > 0) 
        VarPredLen = sqrt(VarPredLen - square(PredLen))/sqrt(LengthData(Ipnt,0));
       Residual = (ObsLen - PredLen)/VarPredLen;
       Ipnt2 += 1;
       FrancisOut(Ipnt2,1) = ObsLen;
       FrancisOut(Ipnt2,2) = PredLen;
       FrancisOut(Ipnt2,3) = VarPredLen;
       FrancisOut(Ipnt2,4) = Residual;
       McAllister = Top/Bot;
       OutFile1 << "Obs" << " " << LengthData(Ipnt) << " " << McAllister << endl;
       OutFile1 << "Pred" << " " << LengthData(Ipnt,-4) << " " <<  LengthData(Ipnt,-3) << " " << LengthData(Ipnt,-2) << " "  << LengthData(Ipnt,-1) << " " << LengthData(Ipnt,0) << " ";
       OutFile1 << PredLengthComp(Ipnt) << " " << McAllister << endl;
      }
  } 
 
 OutFile1 << endl << "Francis_diags" << endl;
 OutFile1 <<  NlengthData << endl;
 OutFile1 <<  "Year Fleet Type Sex EffN Obs_mean_len Pred_mean_len SD_mean_len Residual" << endl;
 Ipnt2 = 0;
 for (Itype=0;Itype<=2;Itype++)
  for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
   {
    VarN = 0; Nout = 0; VarO = 0;
    for (Ipnt=1;Ipnt<=NlengthData;Ipnt++)
     if (LengthData(Ipnt,-3) == Ifleet & LengthData(Ipnt,-1) == Itype)
      {
       Ipnt2 += 1;
       Nout += 1;
       OutFile1 << LengthData(Ipnt,-4) << " " <<  LengthData(Ipnt,-3) << " " << LengthData(Ipnt,-2) << " "  << LengthData(Ipnt,-1) << " " << LengthData(Ipnt,0) << " ";
       OutFile1 << FrancisOut(Ipnt2) << endl; 
       VarN += square(FrancisOut(Ipnt2,4));
       VarO += FrancisOut(Ipnt2,4);
      }
    if (Nout > 1) VarN = VarN - VarO*VarO/Nout;
    if (Nout > 1) OutFile1 << Ifleet << " " << Itype << " " << 1.0/(VarN/float(Nout-1)) << endl;  
   }   
      
 OutFile1 << endl << "AgeLength_fit" << endl;
 OutFile1 <<  NagelengthData << endl;
 OutFile1 << "# Fleet Year Fleet Len1 Len2 Type ZZ Eff_N" << " " << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NagelengthData;Ipnt++)
    if (AgeLengthData(Ipnt,-5) == Ifleet)
     {
      OutFile1 << "Obs" << " " << AgeLengthData(Ipnt) << endl;
      OutFile1 << "Pred" << " " << AgeLengthData(Ipnt,-6) << " " <<  AgeLengthData(Ipnt,-5) << " " << AgeLengthData(Ipnt,-4) << " " <<  AgeLengthData(Ipnt,-3) << " " << AgeLengthData(Ipnt,-2) << " "  << AgeLengthData(Ipnt,-1) << " " << AgeLengthData(Ipnt,0) << " ";
      OutFile1 << PredAgeLengthComp(Ipnt) << endl;
     }
  } 
 OutFile1 << endl;
 OutFile1 << "AgeLength_fit_diags" << endl;
 for (Ipnt=1;Ipnt<=NagelengthData;Ipnt++)
  OutFile1 << Ipnt << " " <<  AgeLengthDataSex(Ipnt) << " " <<  AgeLengthDataType(Ipnt) << " " <<  AgeLengthDataYear(Ipnt) << " " << AgeLengthDataLen1(Ipnt) << " " << AgeLengthDataLen2(Ipnt) << " " << AgeLengthLikeCompData(Ipnt) << endl;

 OutFile1 << endl << "Meansize_fit" << endl;
 OutFile1 <<  NmeansizeData << endl;
 OutFile1 << "# Fleet Year Sex " << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  {
   for (Ipnt=1;Ipnt<=NmeansizeData;Ipnt++)
    if (MeanSizeDataFleet(Ipnt) == Ifleet)
     {
      OutFile1 << "Obs" << " " << MeanSizeDataFleet(Ipnt) << " " << MeanSizeDataYear(Ipnt) << " " <<  MeanSizeDataSex(Ipnt) << " " << MeanSizeData(Ipnt) << endl;
      OutFile1 << "Pred" << " " << MeanSizeDataFleet(Ipnt) << " " << MeanSizeDataYear(Ipnt) << " " <<  MeanSizeDataSex(Ipnt) << " " << PredMeanSize(Ipnt) << endl;
      OutFile1 << "sq" << " " << MeanSizeDataFleet(Ipnt) << " " << MeanSizeDataYear(Ipnt) << " " <<  MeanSizeDataSex(Ipnt) << " " << PredMeanSizeSQ(Ipnt) << endl;
      OutFile1 << "N" << " " << MeanSizeDataFleet(Ipnt) << " " << MeanSizeDataYear(Ipnt) << " " <<  MeanSizeDataSex(Ipnt) << " " << MeanSizeDataSS(Ipnt) << endl;
     }
  } 
 OutFile1 << endl;
 OutFile1 << "MeanSize_fit_diags" << endl;
 for (Ipnt=1;Ipnt<=NmeansizeData;Ipnt++)
  OutFile1 << Ipnt << " " <<  MeanSizeDataFleet(Ipnt) << " " << MeanSizeDataYear(Ipnt) << " " <<  MeanSizeDataSex(Ipnt) << " "  << MeanSizeLikeCompData(Ipnt) << endl;

 OutFile1 << endl << "Tagging_fit" << endl;
 DoTagDiag = 1;
 TaggingLikelihood();
 OutFile1 << NtagData << endl;
 OutFile1 << "# To Come " << endl;
 for (Isex=1;Isex<=2;Isex++)
  for (Isize=1;Isize<=NsizeMax;Isize++)
   for (TimeAtLib=1;TimeAtLib<=6;TimeAtLib++)
    for (Jsize=1;Jsize<=NsizeMax;Jsize++)
     for (Itag=1;Itag<=NtagData;Itag++)
      if (TagData(Itag,1)==Isex&TagData(Itag,2)==Isize&TagData(Itag,3)==Jsize&TagData(Itag,4)==TimeAtLib)
       OutFile1 << TagData(Itag) << " " << PredTag(Itag,0) << " " << TagOffset(Itag) << " " << PredTag(Itag,1) << " " << PredTag(Itag,2)<< " " << PredTag(Itag,3)<< " " << PredTag(Itag,4)<< " " << PredTag(Itag,5)<< " " << PredTag(Itag,6) << endl;
 OutFile1 << endl << "Tagging_fit_table" << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Ipnt=1;Ipnt<=10;Ipnt++)
   {
    OutFile1 << Isex << " " << Ipnt << " " << PredTagTable(Isex,Ipnt,1) << " " << PredTagTable(Isex,Ipnt,2) << endl;
   }
 OutFile1 << endl << "Tagging_EffN_table" << endl;
 OutFile1 << NTagDiag2 << endl;
 for (Itag=1;Itag<=NTagDiag2;Itag++)
  OutFile1 << TagDiag2(Itag) << endl;
 DoTagDiag = 0;

 // Natural mortality
 OutFile1 << endl << "Natural_mortality" << endl;
 OutFile1 << "# Sex Year Platon Age Size" << " " << MidLen << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
   for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    for (Iage=0;Iage<=Nage;Iage++)
     //for (Iyear=First_Year;Iyear<=Last_Year;Iyear++) 
      {
       OutFile1 << Isex << " " << Iyear << " " << Iplatoon << " " << Iage << " ";
       OutFile1 << NatM(Isex,Iyear,Iplatoon,Iage) << " " << endl;
      } 

 // Partial Recruitment
 OutFile1 << endl << "Partial_Recruitment" << endl;
 OutFile1 << "Sex Platoon " << MidLen << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   OutFile1 << Isex << " " << Iplatoon << " " << RecLen(Isex,Iplatoon) << endl;
 
 // Age-specific selectivity
 OutFile1 << endl << "Age_specific_selex_base" << endl;
 OutFile1 <<  NSelexAPatterns << endl;
 OutFile1 << "SelexNo" << Ages << endl;
 for (Iselex=1;Iselex<=NSelexAPatterns;Iselex++)
  OutFile1 << Iselex << " " << BaseSelexA(Iselex) << endl;

 OutFile1 << endl << "Age_specific_selex_base_points" << endl;
 OutFile1 <<  Nfleet*Nsex << endl; Ipnt = 0;
 OutFile1 << "Fleet Sex ";
 for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  OutFile1 << Iyear << " ";
 OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   { Ipnt += 1; OutFile1 << Ifleet << " " << Isex << " " << SelexAPnt(Ipnt) << endl; }

 OutFile1 << endl << "Age_specific_selex_by_year" << endl;
 Ipnt = 0;
 OutFile1 << "Fleet Sex Year" << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ipnt += 1;
    for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  
     OutFile1 << Ifleet << " " << Isex << " " << Iyear << BaseSelexA(SelexAPnt(Ipnt,Iyear)) << endl;
   }  

 // Length-specific selectivity
 OutFile1 << endl << "Length_specific_selex_base" << endl;
 OutFile1 <<  NSelexLPatterns << endl;
 OutFile1 << "SelexNo" << MidLen << endl;
 for (Iselex=1;Iselex<=NSelexLPatterns;Iselex++)
  OutFile1 << Iselex << " " << BaseSelexL(Iselex) << endl;

 OutFile1 << endl << "Length_specific_selex_base_point" << endl;
 OutFile1 <<  Nfleet*Nsex << endl; Ipnt = 0;
 OutFile1 << "Fleet Sex ";
 for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  OutFile1 << Iyear << " ";
 OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   { Ipnt += 1; OutFile1 << Ifleet << " " << Isex << " " << SelexLPnt(Ipnt) << endl; }

 OutFile1 << endl << "Length_specific_selex_by_year" << endl;
 Ipnt = 0;
 OutFile1 << "Fleet Sex Year " << MidLen << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ipnt += 1;
    for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  
     OutFile1 << Ifleet << " " << Isex << " " << Iyear << BaseSelexL(SelexLPnt(Ipnt,Iyear)) << endl;
   }  


 // Retention
 OutFile1 << endl << "Age_specific_retain_base" << endl;
 OutFile1 <<  NRetainAPatterns << endl;
 OutFile1 << "RetainNo" << Ages << endl;
 for (Iselex=1;Iselex<=NRetainAPatterns;Iselex++)
  OutFile1 << Iselex << " " << BaseRetainA(Iselex) << endl;

 OutFile1 << endl << "Age_specific_retain_base_point" << endl;
 OutFile1 <<  Nfleet*Nsex << endl; Ipnt = 0;
 OutFile1 << "Fleet Sex ";
 for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  OutFile1 << Iyear << " ";
 OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   { Ipnt += 1; OutFile1 << Ifleet << " " << Isex << " " << RetainAPnt(Ipnt) << endl; }

 OutFile1 << endl << "Age_specific_selex_by_year" << endl;
 Ipnt = 0;
 OutFile1 << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ipnt += 1;
    if (SelexASpex(Ifleet,Isex,6) > 0)
     for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  
      OutFile1 << Ifleet << " " << Isex << " " << Iyear << BaseRetainA(RetainAPnt(Ipnt,Iyear)) << endl;
   }  

// Retention
 OutFile1 << endl << "Length_specific_retain_base" << endl;
 OutFile1 <<  NRetainLPatterns << endl;
 OutFile1 << "RetainNo" << MidLen << endl;
 for (Iselex=1;Iselex<=NRetainLPatterns;Iselex++)
  OutFile1 << Iselex << " " << BaseRetainL(Iselex) << endl;

 OutFile1 << endl << "Length_specific_retain_base_point" << endl;
 OutFile1 <<  Nfleet*Nsex << endl; Ipnt = 0;
 OutFile1 << "Fleet Sex ";
 for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  OutFile1 << Iyear << " ";
 OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   { Ipnt += 1; OutFile1 << Ifleet << " " << Isex << " " << RetainLPnt(Ipnt) << endl; }

 OutFile1 << endl << "Length_specific_retain_by_year" << endl;
 Ipnt = 0;
 OutFile1 << MidLen << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   {
    Ipnt += 1;
    if (SelexLSpex(Ifleet,Isex,6) > 0)
     for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  
      OutFile1 << Ifleet << " " << Isex << " " << Iyear << BaseRetainL(RetainLPnt(Ipnt,Iyear)) << endl;
   }  

 if (Model_Type == AGEMODEL)
  {
   OutFile1 << endl;
   OutFile1 << "Length_age_selectivity" << endl;
   Ipnt = 0;
   OutFile1 << Ages << endl;
   for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
    for (Isex=1;Isex<=Nsex;Isex++)
     {
      Ipnt += 1; 
      for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)  
       {
        selexL(Ifleet,Isex) = BaseSelexL(SelexLPnt(Ipnt,Iyear)); 
        if (SelexLSpex(Ifleet,Isex,6) > 0)
         retainL(Ifleet,Isex) = BaseRetainL(RetainLPnt(Ipnt,Iyear)); 
        else
         retainL(Ifleet,Isex) = BaseRetainL(0);
      
        temp1 = retainL(Ifleet,Isex)+(1-SelexLSpexR(Ifleet,Isex))*(1-retainL(Ifleet,Isex));
        for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
         {
          OutFile1 << Ifleet << " " << Isex << " " << Iplatoon << " " << Iyear << " ";
          for (Iage=0;Iage<=Nage;Iage++)
           {
            sum1 = 0; 
            for (Isize=1;Isize<=NsizeAge;Isize++)
             sum1 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*retainL(Ifleet,Isex,Isize)*wghtL(Isex,Isize);
            OutFile1 << sum1 << " "; 
           }
          OutFile1 << endl;
          OutFile1 << Ifleet << " " << Isex << " " << Iplatoon << " " << Iyear << " ";
          for (Iage=0;Iage<=Nage;Iage++)
           {
            sum2 = 0;
            for (Isize=1;Isize<=NsizeAge;Isize++)
             sum2 += PhiGrow(Isex,Iplatoon,Iage,Isize)*selexL(Ifleet,Isex,Isize)*temp1(Isize);
            OutFile1 << sum2 << " "; 
           } 
          OutFile1 << endl;
         } 
       } // Year  
     } // Fleet/Sex  
   }  

 if (Model_Type == AGEMODEL)
  {
   OutFile1 << endl;
   OutFile1 << "Age-Length Transition_matrix" << endl;
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Iage= 0;Iage<=Nage;Iage++)
      OutFile1 << Isex << " " << Iplatoon << " " << Iage << " " << PhiGrow(Isex,Iplatoon,Iage) << endl;

   MeanL.initialize();
   SDL.initialize();
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      {
       // Find the mean length and the variance of length-at-age
       for (Iage=0;Iage<=Nage;Iage++)
        {
         Term1 = 0; Term2 = 0;
         for (Ilen=1;Ilen<=NsizeAge;Ilen++)
          {
           Term1 += PhiGrow(Isex,Iplatoon,Iage,Ilen)*MidLen(Ilen);
           Term2 += PhiGrow(Isex,Iplatoon,Iage,Ilen)*square(MidLen(Ilen));
          }
         MeanLP(Isex,Iplatoon,Iage) = Term1;
         if (Term2-Term1*Term1 > 0)
          SDLP(Isex,Iplatoon,Iage) = sqrt(Term2-Term1*Term1);
         else 
          SDLP(Isex,Iplatoon,Iage) = 0;
         MeanL(Isex,Iage) += RecPlat(Isex,Iplatoon)*MeanLP(Isex,Iplatoon,Iage);
        }
      }  
    }  
       
   OutFile1 << endl;
   OutFile1 << "LAA_Transition_matrix" << endl;
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      {
       OutFile1 << Isex << " " << Iplatoon << " ";   
       for (Iage=0;Iage<=Nage;Iage++) OutFile1 << MeanLP(Isex,Iplatoon,Iage) << " "; OutFile1 << endl;
       OutFile1 << Isex << " " << Iplatoon << " ";   
       for (Iage=0;Iage<=Nage;Iage++) OutFile1 << SDLP(Isex,Iplatoon,Iage) << " "; OutFile1 << endl;
       for (Iage=0;Iage<=Nage;Iage++)
        SDL(Isex,Iage) += RecPlat(Isex,Iplatoon)*(square(SDLP(Isex,Iplatoon,Iage))+square(MeanLP(Isex,Iplatoon,Iage)-MeanL(Isex,Iage)));
      } 
     OutFile1 << "# Mean_length_at-age" << endl; 
     OutFile1 << Isex << " ";   
     for (Iage=0;Iage<=Nage;Iage++) OutFile1 << MeanL(Isex,Iage) << " "; OutFile1 << endl;
     OutFile1 << Isex << " ";   
     for (Iage=0;Iage<=Nage;Iage++) OutFile1 << sqrt(SDL(Isex,Iage)+0.00000001) << " "; OutFile1 << endl;
    }  
  }
  
 if (Model_Type == LENMODEL | Model_Type == AGELENMODEL)
  {
   OutFile1 << endl;
   OutFile1 << "Size_transition_matrix" << endl;
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
     for (Ilen=1;Ilen<=NsizeDym;Ilen++)
      {
       OutFile1 << Isex << " " << Iplatoon << " " << Ilen << " ";
       for (Jlen=1;Jlen<=NsizeDym;Jlen++) OutFile1 << " " << TransX(Isex,Iplatoon,Jlen,Ilen) << " ";
       OutFile1 << endl;
      } 

   MeanL.initialize();
   SDL.initialize();
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      {
       // Size-structure for age-0
       Total = 0; for (Ilen=1;Ilen<=NsizeDym;Ilen++)Total+= RecLen(Isex,Iplatoon,Ilen);
       for (Ilen=1;Ilen<=NsizeDym;Ilen++)
        LenMat(Isex,Iplatoon,Ilen,0) = RecLen(Isex,Iplatoon,Ilen)/Total;
       
       // Update by age 
       for (Iage=1;Iage<=NageEst;Iage++)
        {
         for (Ilen=1;Ilen<=NsizeDym;Ilen++)
          { Prob1(Ilen) = LenMat(Isex,Iplatoon,Ilen,Iage-1); Prob2(Ilen) = 0; }
         for (Ilen=1;Ilen<=NsizeDym;Ilen++)
          for (Jlen=1;Jlen<=NsizeDym;Jlen++)
           Prob2(Jlen) += TransX(Isex,Iplatoon,Jlen,Ilen)*Prob1(Ilen);
         for (Ilen=1;Ilen<=NsizeDym;Ilen++)
          LenMat(Isex,Iplatoon,Ilen,Iage) = Prob2(Ilen);
        }  

       // Find the mean length and the variance of length-at-age
       for (Iage=0;Iage<=NageEst;Iage++)
        {
         Term1 = 0; Term2 = 0;
         for (Ilen=1;Ilen<=NsizeDym;Ilen++)
          {
           Term1 += LenMat(Isex,Iplatoon,Ilen,Iage)*MidLen(Ilen);
           Term2 += LenMat(Isex,Iplatoon,Ilen,Iage)*square(MidLen(Ilen));
          }
         MeanLP(Isex,Iplatoon,Iage) = Term1;
         if (Term2-Term1*Term1 > 0)
          SDLP(Isex,Iplatoon,Iage) = sqrt(Term2-Term1*Term1);
         else 
          SDLP(Isex,Iplatoon,Iage) = 0;
         MeanL(Isex,Iage) += RecPlat(Isex,Iplatoon)*MeanLP(Isex,Iplatoon,Iage);
        }
      }  
    }  
       
   OutFile1 << endl;
   OutFile1 << "LAA_Transition_matrix" << endl;
   for (Isex=1;Isex<=Nsex;Isex++)
    {
     for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
      {
       for (Ilen=1;Ilen<=NsizeDym;Ilen++)
        {
         OutFile1 << Isex << " " << Iplatoon << " " << MidLen(Ilen) << " ";
         for (Iage=0;Iage<=NageEst;Iage++) OutFile1 << LenMat(Isex,Iplatoon,Ilen,Iage) << " "; OutFile1 << endl;
        }
       OutFile1 << Isex << " " << Iplatoon << " ";   
       for (Iage=0;Iage<=NageEst;Iage++) OutFile1 << MeanLP(Isex,Iplatoon,Iage) << " "; OutFile1 << endl;
       OutFile1 << Isex << " " << Iplatoon << " ";   
       for (Iage=0;Iage<=NageEst;Iage++) OutFile1 << SDLP(Isex,Iplatoon,Iage) << " "; OutFile1 << endl;
       for (Iage=0;Iage<=NageEst;Iage++)
        SDL(Isex,Iage) += RecPlat(Isex,Iplatoon)*(square(SDLP(Isex,Iplatoon,Iage))+square(MeanLP(Isex,Iplatoon,Iage)-MeanL(Isex,Iage)));
      } 
     OutFile1 << "# Mean_length_at-age" << endl; 
     OutFile1 << Isex << " ";   
     for (Iage=0;Iage<=NageEst;Iage++) OutFile1 << MeanL(Isex,Iage) << " "; OutFile1 << endl;
     OutFile1 << Isex << " ";   
     for (Iage=0;Iage<=NageEst;Iage++) OutFile1 << sqrt(SDL(Isex,Iage)+0.00000001) << " "; OutFile1 << endl;
    }  
  }
  
 OutFile1 << endl;
 OutFile1 << "Biological parameters" << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  {
   if (Model_Type == LENMODEL | Model_Type == AGELENMODEL)
    for (Jlen=1;Jlen<=NsizeDym;Jlen++) 
     {
      OutFile1 << " " << Jlen << " " << MidLen(Jlen) << " " << mat(Isex,Jlen) << " " << wghtL(Isex,Jlen) << " ";
      for (Iage=0;Iage<=Nage;Iage++) OutFile1 << matA(Isex,Iage) << " " << fec(Isex,Jlen,Iage) << " ";
      OutFile1 << endl;
     } 
   else 
    for (Jlen=1;Jlen<=NsizeAge;Jlen++) 
     {
      OutFile1 << " " << Jlen << " " << MidLen(Jlen) << " " << mat(Isex,Jlen) << " " << wghtL(Isex,Jlen) << " ";
      for (Iage=0;Iage<=Nage;Iage++) OutFile1 << matA(Isex,Iage) << " " << fec(Isex,Jlen,Iage) << " " << wghtL(Isex,Jlen) << " ";
      OutFile1 << endl;
     } 
  } 
 
 if (Model_Type == AGEMODEL)
  {
   OutFile1 << endl;
   OutFile1 << "Fecundity-at-age" << endl;
   for (Isex=1;Isex<=Nsex;Isex++)
    for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
    {
     OutFile1 << Isex << " " << Iplatoon;
     for (Iage=0;Iage<=Nage;Iage++) OutFile1 << " " << matA(Isex,Iage) << " " << fecAge(Isex,Iplatoon,Iage) << " ";
     OutFile1 << endl;
    } 
  }  
 
 // N matrix (initial and final)
 OutFile1 << endl;
 OutFile1 << "Numbers_at_age_matrix(init year)" << endl;
 OutFile1 << "Sex Platoon Age ";
 for (Iage= 0;Iage<=Nage;Iage++) OutFile1 << Iage << " ";
 OutFile1 << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    for (Iage= 0;Iage<=Nage;Iage++)
     {
      OutFile1 << Isex << " " << Iplatoon << " "<< Iage << " ";
      for (Isize=1;Isize<=NsizeDym;Isize++) OutFile1 << N(First_Year-1,Isex,Iage,Iplatoon,Isize) << " ";
      OutFile1 << endl; 
     }
   } 
 OutFile1 << endl;
 OutFile1 << "Numbers_at_age_matrix(final year)" << endl;
 OutFile1 << "Sex Platoon Age ";
 for (Iage= 0;Iage<=Nage;Iage++) OutFile1 << Iage << " ";
 OutFile1 << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    for (Iage= 0;Iage<=Nage;Iage++)
     {
      OutFile1 << Isex << " " << Iplatoon << " "<< Iage << " ";
      for (Isize=1;Isize<=NsizeDym;Isize++) OutFile1 << N(Last_Year+1,Isex,Iage,Iplatoon,Isize) << " ";
      OutFile1 << endl; 
     }
   } 
    
 // N matrix (integrated over age)
 OutFile1 << endl;
 OutFile1 << "Numbers_at_size_matrix" << endl;
 OutFile1 << "Year Sex Platoon " << MidLen << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    OutFile1 << "Eqn" << " " << Isex << " " << Iplatoon << " ";
    for (Isize=1;Isize<=NsizeDym;Isize++)
     { temp = 0; for (Iage= 0;Iage<=Nage;Iage++) temp += N(First_Year-1,Isex,Iage,Iplatoon,Isize); OutFile1 << temp << " ";}
    OutFile1 << endl; 
    for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)
     {
      OutFile1 << Iyear << " " << Isex << " " << Iplatoon << " ";
      for (Isize=1;Isize<=NsizeDym;Isize++)
       { temp = 0; for (Iage= 0;Iage<=Nage;Iage++) temp += N(Iyear,Isex,Iage,Iplatoon,Isize); OutFile1 << temp << " ";}
      OutFile1 << endl; 
     }
   } 


 // N matrix (integrated over size)
 OutFile1 << endl;
 OutFile1 << "Numbers_at_age_matrix" << endl;
 OutFile1 << "Year Sex Platoon " << Ages << endl;;
 for (Iage= 0;Iage<=Nage;Iage++)OutFile1 << Iage << " ";
 OutFile1 << endl;
 for (Isex=1;Isex<=Nsex;Isex++)
  for (Iplatoon=1;Iplatoon<=Nplatoon;Iplatoon++)
   {
    OutFile1 << "Eqn" << " " << Isex << " " << Iplatoon << " ";
    for (Iage= 0;Iage<=Nage;Iage++)
     { temp = 0; for (Isize=1;Isize<=NsizeDym;Isize++) temp += N(First_Year-1,Isex,Iage,Iplatoon,Isize); OutFile1 << temp << " ";}
     OutFile1 << endl; 
    for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++)
     {
      OutFile1 << Iyear << " " << Isex << " " << Iplatoon << " ";
      for (Iage= 0;Iage<=Nage;Iage++)
       { temp = 0; for (Isize=1;Isize<=NsizeDym;Isize++) temp += N(Iyear,Isex,Iage,Iplatoon,Isize); OutFile1 << temp << " ";}
      OutFile1 << endl; 
     } 
   } 

 OutFile1 << endl;
 OutFile1 << "Fishing_mortality" << endl;
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    {
     OutFile1 << Iyear << " ";
     for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)OutFile1 << Hrate(Ifleet,Iyear) << " ";
     OutFile1 << endl;
    }

// Composition matrix (integrated over age)
 OutFile1 << endl;
 OutFile1 << "Catch_at_size_matrix (retained)" << endl;
 OutFile1 << "Fleet Year Sex " << MidLen << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    {
     OutFile1 << Ifleet << " " << Iyear << " " << Isex << " ";
     for (Isize=1;Isize<=NsizeMax;Isize++)
      { temp = 0; for (Iage= 0;Iage<=Nage;Iage++) temp += CatchPredAgeSize(1,Ifleet,Iyear,Isex,Iage,Isize); OutFile1 << temp << " ";}
     OutFile1 << endl; 
    }

// Composition matrix (integrated over age)
 OutFile1 << endl;
 OutFile1 << "Catch_at_size_matrix (total)" << endl;
 OutFile1 << "Fleet Year Sex " << MidLen << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    {
     OutFile1 << Ifleet << " " << Iyear << " " << Isex << " ";
     for (Isize=1;Isize<=NsizeMax;Isize++)
      { temp = 0; for (Iage= 0;Iage<=Nage;Iage++) temp += CatchPredAgeSize(2,Ifleet,Iyear,Isex,Iage,Isize); OutFile1 << temp << " ";}
     OutFile1 << endl; 
    }

// Composition matrix (integrated over size)
 OutFile1 << endl;
 OutFile1 << "Catch_at_age_matrix (retained)" << endl;
 OutFile1 << "Fleet Year Sex " << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    {
     OutFile1 << Ifleet << " " << Iyear << " " << Isex << " ";
     for (Iage= 0;Iage<=Nage;Iage++)
      { temp = 0; for (Isize=1;Isize<=NsizeMax;Isize++) temp += CatchPredAgeSize(1,Ifleet,Iyear,Isex,Iage,Isize); OutFile1 << temp << " ";}
     OutFile1 << endl; 
    } 

// Composition matrix (integrated over size)
 OutFile1 << endl;
 OutFile1 << "Catch_at_age_matrix (total)" << endl;
 OutFile1 << "Fleet Year Sex " << Ages << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  for (Isex=1;Isex<=Nsex;Isex++)
   for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    {
     OutFile1 << Ifleet << " " << Iyear << " " << Isex << " ";
     for (Iage=0;Iage<=Nage;Iage++){ temp = 0; for (Isize=1;Isize<=NsizeMax;Isize++) temp += CatchPredAgeSize(2,Ifleet,Iyear,Isex,Iage,Isize); OutFile1 << temp << " ";}
     OutFile1 << endl; 
    } 

 // Indices relative to biomass
 OutFile1 << endl;
 OutFile1 << "Survey biomass indices" << endl;
 OutFile1 << "Fleet" << " ";  for (Iyear=First_Year;Iyear<=Last_Year+1;Iyear++) OutFile1 << Iyear << " "; OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << Ifleet << " " << SurveyBio(Ifleet) << endl;
 
 // Mid-year numbers and biomass
 OutFile1 << endl << "Mid-year biomass" << endl;
 OutFile1 << "Fleet" << " ";  for (Iyear=First_Year;Iyear<=Last_Year;Iyear++) OutFile1 << Iyear << " "; OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << Ifleet << " " << MidExpBio(Ifleet) << endl;
 OutFile1 << endl << "Mid-year numbers" << endl;
 OutFile1 << "Fleet" << " "; for (Iyear=First_Year;Iyear<=Last_Year;Iyear++) OutFile1 << Iyear << " "; OutFile1 << endl;
 for (Ifleet=1;Ifleet<=Nfleet;Ifleet++)
  OutFile1 << Ifleet << " " << MidExpNum(Ifleet) << endl;
 
 OutFile1 << endl;
 OutFile1 << "SSB_final" << endl;
 for (Iyear=FirstProj_Yr;Iyear<=Last_Year+1;Iyear++)
  {
   if (Iyear < First_Year)
    OutFile1 << Iyear << " PRE_FISHERY " << SSB(Iyear) << " " << SSBOut.sd(Iyear) << " " << SSBMid(Iyear) << endl;
   else 
    if (Iyear <= Last_Year)
     OutFile1 << Iyear << " MAIN " << SSB(Iyear) << " " << SSBOut.sd(Iyear) << " " << SSBMid(Iyear) << endl;
    else 
     OutFile1 << Iyear << " MAIN " << SSB(Iyear) << " " << SSBOut.sd(Iyear) << endl;
   }  
 OutFile1 << endl;
 OutFile1 << "Recruitment_final" << endl;
 for (Iyear=FirstProj_Yr;Iyear<=Last_Year+1;Iyear++)
  if (Iyear < First_Year)
   {
    if (Nsex==1)
     OutFile1 << Iyear << " PRE_FISHERY " << RecOut(Iyear,1);
    else 
     OutFile1 << Iyear << " PRE_FISHERY " << RecOut(Iyear,1) << " " << RecOut(Iyear,2);
    OutFile1 << endl;
   } 
  else 
   {
    if (Nsex==1)
     OutFile1 << Iyear << " MAIN " << RecOut(Iyear,1);
    else 
     OutFile1 << Iyear << " MAIN " << RecOut(Iyear,1) << " " << RecOut(Iyear,2);
    OutFile1 << endl;
   } 
 OutFile1 << endl; 
  
 for (Ipnt=1;Ipnt<=NindexData;Ipnt++)
  if (IndexPred(2,Ipnt) < 0) { cout << "Index pred < 0" << endl; exit(1); }

 for (Ifleet=1;Ifleet<=NfishFleet;Ifleet++)
  for (Iyear=First_Year;Iyear<=Last_Year;Iyear++)
    if (MidExpBio(Ifleet,Iyear) < 0) { cout << "MidExpBio < 0" << endl;  exit(1); }
   
 //if (sum(CatchLikeCompFleet) > 5) { cout << "CatchLike" << endl; exit(1); }

//***********************************************************************

RUNTIME_SECTION
 maximum_function_evaluations 3000,3000,3000,50000   
 convergence_criteria 00001,.00001,1e-9 


FINAL_SECTION
 IsFinalPhase = 1;
 CreateOutput();


 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;

//***********************************************************************

REPORT_SECTION
 
 
 save_gradients(gradients);
 GradFile.close();
 GradFile.open("gradient.1");
 GradFile >> hello >> hello >> hello;
 //cout << hello << endl;
 //exit(1);
 CreateOutput();
