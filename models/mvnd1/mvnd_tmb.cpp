#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_SCALAR(covar);
  PARAMETER(mu); // parameters
  using namespace density;
  // MVNORM_t<Type> neg_log_density(covar);
  // Type nll= neg_log_density(mu);
  Type nll= mu*mu; //dnorm(mu_bounded, Type(0), covar, true);
  return(nll);
}
