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
  sigma_obs.allocate(0,5,"sigma_obs");
  logLinf_mean.allocate("logLinf_mean");
  logk_mean.allocate("logk_mean");
  logLinf_sigma.allocate(0,5,"logLinf_sigma");
  logk_sigma.allocate(0,5,"logk_sigma");
  logLinf.allocate(1,Nfish,"logLinf");
  logk.allocate(1,Nfish,"logk");
  jnll.allocate("jnll");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::initializationfunction(void)
{
  sigma_obs.set_initial_value(.11);
  logLinf_mean.set_initial_value(3.1);
  logk_mean.set_initial_value(-2.1);
  logLinf_sigma.set_initial_value(.11);
  logk_sigma.set_initial_value(.11);
  logLinf.set_initial_value(3);
  logk.set_initial_value(-2);
}

void model_parameters::userfunction(void)
{
  jnll =0.0;
  jnll=0;
  dvar_vector ypred(1,Nobs);
  dvariable Linf;
  dvariable k;
  double delta=1.0;
  // Priors
  jnll+= pow((sigma_obs-.1)/.1, 2);
  jnll+= pow((logLinf_sigma-.1)/.1, 2);
  jnll+= pow((logk_sigma-.2)/.1, 2);
  // Random effects
  for(int i=1; i<=Nfish; i++){
     jnll+= pow((logLinf_mean-logLinf(fish(i)))/logLinf_sigma, 2);
     jnll+= pow((logk_mean-logk(fish(i)))/logk_sigma, 2);
  }
  // Predict fish length and add likelihood
  for(int i=1; i<=Nobs; i++){
     Linf=exp(logLinf(fish(i)));
     k=exp(logk(fish(i)));
     ypred(i)=pow(log( Linf*(1-exp(-k*(ages(i)-5.0)))), delta);
     jnll+=pow( (loglengths(i)-ypred(i))/sigma_obs,2);
  }
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
