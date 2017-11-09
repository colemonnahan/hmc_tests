#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
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
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <wildf.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  Ndata.allocate("Ndata");
  Nstage.allocate("Nstage");
  Nyear.allocate("Nyear");
  Nplant.allocate("Nplant");
  year.allocate(1,Ndata,"year");
  plant.allocate(1,Ndata,"plant");
  stage.allocate(1,Ndata,"stage");
  Pods.allocate(1,Ndata,"Pods");
  toF.allocate(1,Ndata,"toF");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  yearInterceptSD.allocate("yearInterceptSD");
  plantInterceptSD.allocate("plantInterceptSD");
  plantSlopeSD.allocate("plantSlopeSD");
  intercept.allocate(1,Nstage,"intercept");
  slope.allocate("slope");
  yearInterceptEffect_raw.allocate(1,Nyear,"yearInterceptEffect_raw");
  plantInterceptEffect_raw.allocate(1,Nplant,"plantInterceptEffect_raw");
  plantSlopeEffect_raw.allocate(1,Nplant,"plantSlopeEffect_raw");
  ypred.allocate("ypred");
  #ifndef NO_AD_INITIALIZE
  ypred.initialize();
  #endif
  yearInterceptSD2.allocate("yearInterceptSD2");
  #ifndef NO_AD_INITIALIZE
  yearInterceptSD2.initialize();
  #endif
  plantInterceptSD2.allocate("plantInterceptSD2");
  #ifndef NO_AD_INITIALIZE
  plantInterceptSD2.initialize();
  #endif
  plantSlopeSD2.allocate("plantSlopeSD2");
  #ifndef NO_AD_INITIALIZE
  plantSlopeSD2.initialize();
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
  yearInterceptSD.set_initial_value(.1);
  plantInterceptSD.set_initial_value(.1);
  plantSlopeSD.set_initial_value(.1);
  slope.set_initial_value(.001);
}

void model_parameters::userfunction(void)
{
  nld =0.0;
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

void model_parameters::report(const dvector& gradients){}

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
