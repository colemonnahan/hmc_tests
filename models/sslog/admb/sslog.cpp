#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <sslog.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  N.allocate("N");
  catches.allocate(1,N,"catches");
  logcpue.allocate(1,N,"logcpue");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  logK.allocate("logK");
  logr.allocate("logr");
  iq.allocate(1,10,"iq");
  isigma2.allocate("isigma2");
  itau2.allocate("itau2");
  u_raw.allocate(1,N,"u_raw");
  u.allocate(1,N,"u");
  #ifndef NO_AD_INITIALIZE
    u.initialize();
  #endif
  ypred.allocate(1,N,"ypred");
  #ifndef NO_AD_INITIALIZE
    ypred.initialize();
  #endif
  B.allocate(1,N,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  temp.allocate("temp");
  #ifndef NO_AD_INITIALIZE
  temp.initialize();
  #endif
  sigma2.allocate("sigma2");
  #ifndef NO_AD_INITIALIZE
  sigma2.initialize();
  #endif
  sigma.allocate("sigma");
  #ifndef NO_AD_INITIALIZE
  sigma.initialize();
  #endif
  tau2.allocate("tau2");
  #ifndef NO_AD_INITIALIZE
  tau2.initialize();
  #endif
  q.allocate("q");
  #ifndef NO_AD_INITIALIZE
  q.initialize();
  #endif
  K.allocate("K");
  #ifndef NO_AD_INITIALIZE
  K.initialize();
  #endif
  r.allocate("r");
  #ifndef NO_AD_INITIALIZE
  r.initialize();
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
  logK.set_initial_value(6.5);
  logr.set_initial_value(-1);
  iq.set_initial_value(5);
  isigma2.set_initial_value(400);
  itau2.set_initial_value(92);
}

void model_parameters::userfunction(void)
{
  nld =0.0;
 nlp=0.0;
 nll=0.0;
  sigma2 = 1/isigma2;
  sigma=sqrt(sigma2);
  tau2 = 1/itau2;
  q = 1/iq;
  K = exp(logK);
  r = exp(logr);
  // non-centered random effects used below; implies u~N(0,sigma)
  u=u_raw*sigma;
   // priors; note the parameterization is different between stan, where
  // scale is 1/scale here.
  nlp+=dnorm(logr, -1.38, 0.51);
  nlp+=dgamma(iq, 0.001, 0.001);
  nlp+=dgamma(isigma2, 3.785518, 0.010223);
  nlp+=dgamma(itau2, 1.708603, 0.008613854);
  // project dynamics
 B(1) = K;
 ypred(1) = log(B(1)) +log(q);
 for(int i=2; i<=N; i++){
   temp = (B(i-1)+r*B(i-1)*(1.0-B(i-1)/K)-catches(i-1))*exp(u(i));
   if(temp<.001){
     //increment_log_prob(-1*(temp-1)^2);
     B(i) = 1/(2-temp/.001);
   } else {
     B(i) = temp;
   }
    B(i) = temp;
   ypred(i) = log(B(i)) +log(q);
 }
 // likelihoods
 nll+=dnorm(u_raw, 0.0, 1.0);
 // The likelihood
 // nll-=dnorm(ypred, logcpue, sqrt(tau2));
 dvariable SS=norm2(logcpue-ypred);
 dvariable tmp= N*(0.5*log(2.*M_PI)+log(sqrt(tau2)))+0.5*SS/tau2;
 nll+= tmp;
 nld=nll+nlp;
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
