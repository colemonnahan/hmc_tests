#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_SCALAR(covar);
  DATA_VECTOR(bounds); // lower and upper bounds
  PARAMETER(mu); // parameters
  // vector<Type> mu_bounded(mu.size());
  Type mu_bounded;
  Type nll=0;
  mu_bounded=bounds(0)+(bounds(1)-bounds(0))/(1+exp(-mu));
  // jacobian adjustment by hand
  Type log_scales =
    log((bounds(1)-bounds(0))*(1/(1+exp(-mu)))*(1-1/(1+exp(-mu))));
    // log((bounds(1)-bounds(0))) -mu -2*log(1+exp(-mu)) ;
   nll -= -mu_bounded*mu_bounded;
    nll -= log_scales;
   REPORT(mu_bounded);
   REPORT(log_scales);
  return(nll);
}
