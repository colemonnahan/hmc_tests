#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Nfish);
  DATA_INTEGER(Nobs);
  DATA_VECTOR(loglengths);
  DATA_IVECTOR(fish);
  DATA_IVECTOR(ages);

  // fixed effects
  PARAMETER(delta);
  PARAMETER(sigma_obs);

  // hyperparameters with bounds added in R
  PARAMETER(logLinf_mean);
  PARAMETER(logk_mean);
  PARAMETER(logLinf_sigma);
  PARAMETER(logk_sigma);

  // non-centered random effects
  PARAMETER_VECTOR(logLinf_raw);
  PARAMETER_VECTOR(logk_raw);
  vector<Type> logk_raw(Nfish);

  // transform random effects
  vector<Type> logk_raw(Nfish);
  vector<Type> logLinf_raw(Nfish);
  logLinf = logLinf_mean+logLinf_raw*logLinf_sigma;
  logk = logk_mean+logk_raw*logk_sigma;

  // Calculate predicted lengths
  vector<Type> ypred;
  Type Linf;
  Type k;
  for(i in 1:Nobs){
    Linf  = exp(logLinf(fish(i)));
    k = exp(logk(fish(i)));
    ypred(i) = log( Linf*(1-exp(-k*(ages(i)-5)))^delta );
  }

  Type logprior=0.0;
  Type loglikelihood=0.0;

  // delta is uniform above
  sigma_obs~cauchy(0,5);

  // hyperpriors
  logLinf_sigma~cauchy(0,5);
  logk_sigma~cauchy(0,5);
  // hyper means are uniform above

  // random effects; non-centered
  logLinf_raw~normal(0, 1);
  logk_raw~normal(0, 1);

  // calculate likelihood of data

  loglengths~normal(ypred, sigma_obs);

    Type nll= -1*dnorm(mu, x, Type(1.0), true).sum();
  return(nll);
}
