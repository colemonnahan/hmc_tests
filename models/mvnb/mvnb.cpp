#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(covar);
  PARAMETER(mu1); // parameters
  PARAMETER(mu2); // parameters
  PARAMETER(mu3); // parameters

  Type nll=0;
  vector<Type> mu(3);
  mu(0)=mu1;
  mu(1)=mu2;
  mu(2)=mu3;
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  nll += neg_log_density(mu);
  return(nll);
}
