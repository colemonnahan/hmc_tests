#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_MATRIX(covar);
  PARAMETER_VECTOR(mu); // parameters
  vector<Type> mu_bounded(2);
  Type nll=0;
  for(int i=0; i<=1; i++){
    mu_bounded(i)=-2+4/(1+exp(-mu(i)));
    // jacobian adjustment by hand
    nll -= log( (4.0*exp(-mu(i)))/pow(1+exp(-mu(i)),2) );
  }
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  nll -= neg_log_density(mu_bounded);
  REPORT(mu_bounded);
  return(nll);
}
