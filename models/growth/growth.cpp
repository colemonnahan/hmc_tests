#include <TMB.hpp>

// Hand-coded Cauchy distribution
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

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

  // transform random effects
  vector<Type> logk(Nfish);
  vector<Type> logLinf(Nfish);
  logLinf = logLinf_mean+logLinf_raw*logLinf_sigma;
  logk = logk_mean+logk_raw*logk_sigma;

  // Calculate predicted lengths
  vector<Type> ypred(Nobs);
  Type Linf; Type k;
  for(int i=0; i<Nobs; i++){
    Linf  = exp(logLinf(fish(i)-1));
    k = exp(logk(fish(i)-1));
    ypred(i) = log( Linf*pow(1-exp(-k*(ages(i)-5)),delta) );
  }

  Type nlp=0.0; // negative log prior
  Type nll=0.0; // negative log likelihood

  // delta is uniform above
  nlp -= dcauchy(sigma_obs, Type(0), Type(5), true);

  // hyperpriors
  nlp -= dcauchy(logk_sigma, Type(0), Type(5), true);
  nlp -= dcauchy(logLinf_sigma, Type(0), Type(5), true);
  // hyper means are uniform above

  // random effects; non-centered
  nll-=dnorm(logLinf_raw, Type(0), Type(1), true).sum();
  nll-=dnorm(logk_raw, Type(0), Type(1), true).sum();

  // likelihood of data
  nll-=dnorm(ypred, loglengths, sigma_obs, true).sum();
  //Type nll= -1*dnorm(mu, x, Type(1.0), true).sum();

  Type nld=nll+nlp; // negative log density

  REPORT(ypred);
  return(nld);
}
