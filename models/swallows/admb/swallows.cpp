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
 dvariable inv_logit(const dvariable& x){
  dvariable y= 1/(1+mfexp(-x));
  return y;
 }
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <swallows.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  I.allocate("I");
  K.allocate("K");
  nfam.allocate("nfam");
  CH0.allocate(1,K,1,I,"CH0");
  CH.allocate(1,I,1,K);
 CH = trans(CH0);
  carez.allocate(1,I,"carez");
  year.allocate(1,I,"year");
  agec.allocate(1,K,"agec");
  family.allocate(1,I,"family");
  last.allocate(1,I,"last");
  ones.allocate(1,I,1,K,"ones");
  ones2.allocate(1,I,"ones2");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  sigmayearphi.allocate("sigmayearphi");
  sigmaphi.allocate("sigmaphi");
  sigmap.allocate("sigmap");
  a.allocate(1,17,"a");
  a1.allocate("a1");
  b0.allocate(1,4,"b0");
  b1.allocate(1,4,"b1");
  fameffphi_raw.allocate(1,72,"fameffphi_raw");
  fameffp_raw.allocate(1,72,"fameffp_raw");
  yeareffphi_raw.allocate(1,4,"yeareffphi_raw");
  sigmayearphi2.allocate("sigmayearphi2");
  #ifndef NO_AD_INITIALIZE
  sigmayearphi2.initialize();
  #endif
  sigmaphi2.allocate("sigmaphi2");
  #ifndef NO_AD_INITIALIZE
  sigmaphi2.initialize();
  #endif
  sigmap2.allocate("sigmap2");
  #ifndef NO_AD_INITIALIZE
  sigmap2.initialize();
  #endif
  x.allocate("x");
  #ifndef NO_AD_INITIALIZE
  x.initialize();
  #endif
  nlp.allocate("nlp");
  #ifndef NO_AD_INITIALIZE
  nlp.initialize();
  #endif
  nll.allocate("nll");
  #ifndef NO_AD_INITIALIZE
  nll.initialize();
  #endif
  p.allocate(1,I,1,K,"p");
  #ifndef NO_AD_INITIALIZE
    p.initialize();
  #endif
  phi.allocate(1,I,1,K-1,"phi");
  #ifndef NO_AD_INITIALIZE
    phi.initialize();
  #endif
  chi.allocate(1,I,1,K+1,"chi");
  #ifndef NO_AD_INITIALIZE
    chi.initialize();
  #endif
  nld.allocate("nld");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::initializationfunction(void)
{
  sigmaphi.set_initial_value(1);
  sigmayearphi.set_initial_value(1);
  sigmap.set_initial_value(1);
  a.set_initial_value(.1);
  a1.set_initial_value(.1);
  b0.set_initial_value(.1);
  b1.set_initial_value(.1);
  fameffphi_raw.set_initial_value(.1);
  fameffp_raw.set_initial_value(.1);
  yeareffphi_raw.set_initial_value(.1);
}

void model_parameters::userfunction(void)
{
  nld =0.0;
 nlp=0.0; // negative log prior
 nll=0.0; // negative log likelihood
 // Jacobian below
 sigmayearphi2=exp(sigmayearphi);
 sigmaphi2=exp(sigmaphi);
 sigmap2=exp(sigmap);
 int k;
 for(int i=1; i<=I; i++){ // loop over each individual
  // calculate phi as a function of fixed and random effects
  for(int t=1; t<=(K-1); t++) {
    x=a(t)+ a1*carez(i)+
      sigmayearphi2*yeareffphi_raw(year(i))+
      sigmaphi2*fameffphi_raw(family(i));
    phi(i,t) = inv_logit(x);
  }
  // calculate p as a function of fixed and random effects
  p(i,1) = 1;  // first occasion is marking occasion
  for(int t=2; t<=K; t++){
    x=b0(year(i))+ b1(year(i))*agec(t)+
      sigmap2*fameffp_raw(family(i));
    p(i,t) = inv_logit(x);
  }
  // probabilitiy of never being seen after last observation. ind here is
  // a reverse index so this loop goes from K:2, recursively calculating
  // backward.
  chi(i,K+1) = 1.0;
  k=K;
  while (k > 1) {
    chi(i,k) = (1-phi(i,k-1)) + phi(i,k-1) * (1-p(i,k)) * chi(i,k+1);
    k=k-1;
  }
  chi(i,1) = (1-p(i,1)) * chi(i,2);
 }
 // Jacobian adjustment for variances
 nll -= sigmaphi + sigmayearphi + sigmap;
 // priors
 nlp+= dcauchy(sigmaphi2, 0.0, 1.0);
 nlp+= dnorm(sigmayearphi2, 0.0, 3);
 nlp+= dcauchy(sigmap2, 0.0, 1.0);
 nlp+= dnorm(a, 0.0, 1.5);
 nlp+= dnorm(a1, 0.0, 5.0);
 nlp+= dnorm(b0, 0.0, 5.0);
 nlp+= dnorm(b1, 0.0, 5.0);
 // random effects; non-centered
 nll+= dnorm(fameffphi_raw, 0,1);
 nll+= dnorm(fameffp_raw, 0,1);
 nll+= dnorm(yeareffphi_raw, 0,1);
 // // likelihood
 for(int i=1; i<=I; i++){ // loop over each individual
    // probability of survival, known alive since k<last
    for (int t=2; t<=last(i); t++) {
      // ones[i,t]~bernoulli(phi[i, t-1]);
    	nll-= log(phi(i,t-1));
    }
    // // probability of observation given known alive
    for(int t=1; t<=last(i); t++){
      // CH[i,t]~bernoulli(p[i,t]);
      if(CH(i,t)==1){
    	nll-= log(p(i,t));
      } else {
    	nll-= log(1-p(i,t));
      }
    }
    // probability of no observations after time period last
    // ones2[i]~bernoulli(chi[i,last[i]+1]);
     nll-= log(chi(i,last(i)+1));
  }
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
