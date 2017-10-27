GLOBALS_SECTION
 #include <admodel.h>
 #include "statsLib.h"
 #include <df1b2fun.h>
 #include <adrndeff.h> // Negative log-likelihood for the Cauchy distribution; scalar arguments.
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
  init_number Nfish;
  init_number Nobs;
  init_vector loglengths(1,Nobs);
  init_ivector fish(1,Nobs);
  init_ivector ages(1,Nobs);

PARAMETER_SECTION
  init_bounded_number delta(.5,5);
  init_number sigma_obs;
  init_number logLinf_mean;
  init_number logk_mean;
  init_number logLinf_sigma;
  init_number logk_sigma;
  init_vector logLinf_raw(1,Nfish);
  init_vector logk_raw(1,Nfish);
  vector logLinf(1,Nfish);
  vector logk(1,Nfish);
  vector ypred(1,Nobs);
  number Linf;
  number k;
  number nlp;
  number nll;
  objective_function_value nld;

INITIALIZATION_SECTION
 delta 1;
 sigma_obs .1;
 logLinf_sigma .1;
 logk_sigma .1;

PROCEDURE_SECTION
 nlp=0.0;
 nll=0.0;
 logLinf = logLinf_mean+logLinf_raw*logLinf_sigma;
 logk = logk_mean+logk_raw*logk_sigma;
 double zero=0.0; double five=5.0;
  // Calculate predicted lengths
  for(int i=1; i<=Nobs; i++){
    Linf  = exp(logLinf(fish(i)));
    k = exp(logk(fish(i)));
    ypred(i) = log( Linf*pow(1-exp(-k*(ages(i)-5)),delta) );
  }
  // delta is uniform above
  nlp += dcauchy(sigma_obs, zero, five);

  // hyperpriors
  nlp += dcauchy(logk_sigma, zero, five);
  nlp += dcauchy(logLinf_sigma, zero, five);
  // hyper means are uniform above

  // random effects; non-centered
  nll+=dnorm(logLinf_raw, 0,1);
  nll+=dnorm(logk_raw, 0,1);

  // likelihood of data
   dvariable SS=norm2(loglengths-ypred);
   dvariable tmp=Nobs*(0.5*log(2.*M_PI)+log(sigma_obs))+0.5*SS/(sigma_obs*sigma_obs);
  nll+=tmp;
  // for(int i=1; i<=Nobs; i++){	  
  //   nll+=dnorm(loglengths(i), ypred(i), sigma_obs);
  // }
  nld=nll+nlp; // negative log density

REPORT_SECTION
 cout << ypred << endl;

