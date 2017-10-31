#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <growth.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  Nfish.allocate("Nfish");
  Nobs.allocate("Nobs");
  loglengths.allocate(1,Nobs,"loglengths");
  fish.allocate(1,Nobs,"fish");
  ages.allocate(1,Nobs,"ages");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  delta.allocate(.5,5,"delta");
  sigma_obs.allocate("sigma_obs");
  logLinf_mean.allocate("logLinf_mean");
  logk_mean.allocate("logk_mean");
  logLinf_sigma.allocate("logLinf_sigma");
  logk_sigma.allocate("logk_sigma");
  logLinf_raw.allocate(1,Nfish,"logLinf_raw");
  logk_raw.allocate(1,Nfish,"logk_raw");
  logLinf.allocate(1,Nfish,"logLinf");
  #ifndef NO_AD_INITIALIZE
    logLinf.initialize();
  #endif
  logk.allocate(1,Nfish,"logk");
  #ifndef NO_AD_INITIALIZE
    logk.initialize();
  #endif
  ypred.allocate(1,Nobs,"ypred");
  #ifndef NO_AD_INITIALIZE
    ypred.initialize();
  #endif
  sigma_obs2.allocate("sigma_obs2");
  #ifndef NO_AD_INITIALIZE
  sigma_obs2.initialize();
  #endif
  logk_sigma2.allocate("logk_sigma2");
  #ifndef NO_AD_INITIALIZE
  logk_sigma2.initialize();
  #endif
  logLinf_sigma2.allocate("logLinf_sigma2");
  #ifndef NO_AD_INITIALIZE
  logLinf_sigma2.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  nlp.allocate("nlp");
  #ifndef NO_AD_INITIALIZE
  nlp.initialize();
  #endif
  nll.allocate("nll");
  #ifndef NO_AD_INITIALIZE
  nll.initialize();
  #endif
  nld.allocate("nld");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::initializationfunction(void)
{
  delta.set_initial_value(1);
  sigma_obs.set_initial_value(-2.3);
  logLinf_sigma.set_initial_value(-2.3);
  logk_sigma.set_initial_value(-2.3);
  logk_mean.set_initial_value(-2.3);
  logLinf_mean.set_initial_value(3.9);
}

void model_parameters::userfunction(void)
{
  nld =0.0;
 nlp=0.0;
 nll=0.0;
 // Exponentiate sds since ADMB can't do (0,Inf) bounds. Need to adjust for
 // Jacobian below
 sigma_obs2=exp(sigma_obs);
 logLinf_sigma2=exp(logLinf_sigma);
 logk_sigma2=exp(logk_sigma);
 logLinf = logLinf_mean+logLinf_raw*logLinf_sigma2;
 logk = logk_mean+logk_raw*logk_sigma2;
 double zero=0.0; double five=5.0;
  // Calculate predicted lengths
  for(int i=1; i<=Nobs; i++){
    Linf  = exp(logLinf(fish(i)));
    k = exp(logk(fish(i)));
    ypred(i) = log( Linf*pow(1-exp(-k*(ages(i)-5)),delta) );
  }
  // delta is uniform above
  nlp += dnorm(delta, 1, 0.25);
  nlp += dcauchy(sigma_obs2, zero, five);
  // hyperpriors
  nlp += dcauchy(logk_sigma2, zero, five);
  nlp += dcauchy(logLinf_sigma2, zero, five);
  // hyper means are uniform above
  // Jacobian adjustment for variances
  nll -= sigma_obs + logk_sigma + logLinf_sigma;
  // random effects; non-centered
  nll+=dnorm(logLinf_raw, 0,1);
  nll+=dnorm(logk_raw, 0,1);
  // likelihood of data
   dvariable SS=norm2(loglengths-ypred);
   dvariable tmp=Nobs*(0.5*log(2.*M_PI)+log(sigma_obs2))+0.5*SS/(sigma_obs2*sigma_obs2);
  nll+=tmp;
  // for(int i=1; i<=Nobs; i++){
  //   nll+=dnorm(loglengths(i), ypred(i), sigma_obs);
  // }
  nld=nll+nlp; // negative log density
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
 cout << ypred << endl;
}

void model_parameters::preliminary_calculations(void){
#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
