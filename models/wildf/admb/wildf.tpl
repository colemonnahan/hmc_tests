
GLOBALS_SECTION
 #include <admodel.h>
 #include "statsLib.h"
 #include <df1b2fun.h>
 #include <adrndeff.h>
 // Negative log-likelihood for the Cauchy distribution; scalar arguments.
 dvariable dcauchy(const prevariable& mu, const double& x, const double& scale)
 {
  #ifdef M_PI
  double pi = M_PI;
  #else
  double pi = 3.14159265358979323844;
  #endif
  dvariable logres=0.0;
  logres+=log(pi);
  logres+= log(scale);
  logres+= log(1 + pow( (x-mu)/scale, 2));
  return logres;
 }

DATA_SECTION
  init_number Ndata;
  init_number Nstage;
  init_number Nyear;
  init_number Nplant;
  init_ivector year(1,Ndata);
  init_ivector plant(1,Ndata);
  init_ivector stage(1,Ndata);
  init_vector Pods(1,Ndata);
  init_ivector toF(1,Ndata);

PARAMETER_SECTION
  init_number yearInterceptSD;
  init_number plantInterceptSD;
  init_number plantSlopeSD;
  init_vector intercept(1,Nstage);
  init_number slope;
  // random effect vectors
  init_vector yearInterceptEffect_raw(1,Nyear);
  init_vector plantInterceptEffect_raw(1,Nplant);
  init_vector plantSlopeEffect_raw(1,Nplant);
  // predictions
  number ypred;
  // exponentiated versions of variance parameters since ADMB can't do lower bounds
  number yearInterceptSD2;
  number plantInterceptSD2;
  number plantSlopeSD2;
  number nlp;
  number nll;
  objective_function_value nld;

 INITIALIZATION_SECTION
 yearInterceptSD .1;
 plantInterceptSD .1;
 plantSlopeSD .1;
 slope .001;


PROCEDURE_SECTION
 nlp=0.0;
 nll=0.0;
 double zero=0.0; double five=5.0;
 // Exponentiate sds since ADMB can't do (0,Inf) bounds. Need to adjust for
 // Jacobian below
 yearInterceptSD2=exp(yearInterceptSD);
 plantInterceptSD2=exp(plantInterceptSD);
 plantSlopeSD2=exp(plantSlopeSD);

 // priors
 nlp+= dcauchy(yearInterceptSD2, zero, five);
 nlp+= dcauchy(plantInterceptSD2, zero, five);
 nlp+= dcauchy(plantSlopeSD2, zero, five);
 nlp+= dnorm(slope, 0.0, 10.0);
 nlp+= dnorm(intercept, 0.0, 10.0);

// Model predictions
 for(int i=1; i<=Ndata; i++){
 // prediction logit scale
  ypred= intercept(stage(i)) +
   yearInterceptEffect_raw(year(i))*yearInterceptSD2 +
   plantInterceptEffect_raw(plant(i))*plantInterceptSD2+
   Pods(i) * plantSlopeEffect_raw(plant(i))*plantSlopeSD2+
   Pods(i) * slope;
 // likelihood contribution
 if(toF(i)==1){
  nll+= log(1+exp(-ypred));
   } else {
  nll+= ypred+log(1+exp(-ypred));
 }
 }
  // Jacobian adjustment for variances
  nll -= yearInterceptSD + plantInterceptSD + plantSlopeSD;

  // random effects; non-centered
  nll+=dnorm(yearInterceptEffect_raw, 0,1);
  nll+=dnorm(plantInterceptEffect_raw, 0,1);
  nll+=dnorm(plantSlopeEffect_raw, 0,1);

 // temporrary strong priors to get an invertible hessian, also set nld=nlp;
 // and it will work
 //nlp+=dnorm(plantInterceptSD2, .1, .1);
 // nlp+=dnorm(yearInterceptSD2, .1, .1); nlp+=dnorm(plantSlopeSD2, .1, .1);
 // nlp+=dnorm(slope, .1, .1); nlp+=dnorm(intercept, 0.0, 10.0);
 // nlp+=dnorm(yearInterceptEffect_raw, .1,1);
 // nlp+=dnorm(plantInterceptEffect_raw, .1,1);
 // nlp+=dnorm(plantSlopeEffect_raw, .1,1);

 nld=nll+nlp; // negative log density

// REPORT_SECTION
//  cout << ypred << endl;
//  cout << nll2 << endl;
//   cout << sum(nll2) << endl;
//  cout <<  yearInterceptSD + plantInterceptSD + plantSlopeSD << endl;
//  cout << dnorm(plantInterceptEffect_raw, 0,1) << endl;
//  cout << dnorm(plantSlopeEffect_raw, 0,1) << endl;
//  cout << dnorm(yearInterceptEffect_raw,0,1) << endl;
//  cout << "nll=" << nll << "; nlp=" << nlp << endl;
