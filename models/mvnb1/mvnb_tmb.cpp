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
  Type log_scales = log( (bounds(1)-bounds(0))*exp(-mu)/pow(1+exp(-mu),2) );
   // for(int i=0; i<=mu.size()-1; i++){
  //   mu_bounded(i)=bounds(0)+(bounds(1)-bounds(0))/(1+exp(-mu(i)));
  //   // jacobian adjustment by hand
  //   nll -= log( (bounds(1)-bounds(0))*exp(-mu(i))/pow(1+exp(-mu(i)),2) );
  // }
   using namespace density;
  // MVNORM_t<Type> neg_log_density(covar);
  // nll += neg_log_density(mu_bounded);
   nll -= dnorm(mu_bounded, Type(0), covar, true);
   nll -= log_scales;
   REPORT(mu_bounded);
   REPORT(log_scales);
  return(nll);
}
