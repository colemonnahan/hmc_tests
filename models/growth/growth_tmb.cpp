#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nfish);
  DATA_INTEGER(Nobs);
  DATA_VECTOR(loglengths);
  DATA_IVECTOR(fish);
  DATA_IVECTOR(ages);
  // PARAMETER(delta);
  PARAMETER(sigma_obs);
  PARAMETER(logLinf_mean);
  PARAMETER(logk_mean);
  PARAMETER(logLinf_sigma);
  PARAMETER(logk_sigma);
  PARAMETER_VECTOR(logLinf);
  PARAMETER_VECTOR(logk);

  // Start of model
  vector<Type> ypred(Nobs); // predictions
  Type Linf;
  Type k;
  Type nll=0; 			// negative log likelihood

  // Priors
  nll-= dnorm(sigma_obs, Type(0), Type(5), true);
  nll-= dnorm(logLinf_sigma, Type(0), Type(5), true);
  nll-= dnorm(logk_sigma, Type(0), Type(5), true);

 // Random effects
nll-= dnorm(logLinf_mean, logLinf, logLinf_sigma, true).sum();
nll-= dnorm(logk_mean, logk, logk_sigma, true).sum();

  // Calculate likelihood
  Type delta=1;
   for(int i=0; i<Nobs; i++){
    Linf = exp(logLinf(fish(i)-1));
    k = exp(logk(fish(i)-1));
    ypred(i) = log( Linf*pow(1-exp(-k*(ages(i)-Type(5))),delta));
	nll-=dnorm(loglengths(i), ypred(i), sigma_obs, true);
   }
  return nll;
}

