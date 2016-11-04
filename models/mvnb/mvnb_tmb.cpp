#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_MATRIX(covar);
  PARAMETER_VECTOR(mu); // parameters
  vector<Type> mu_bounded(2);
  for(int i=0; i<=1; i++)
    mu_bounded(i)=-2+(2- (-2))/(1+exp(-mu(i)));
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  Type nll= neg_log_density(mu_bounded);
  REPORT(nll);
  return(nll);
}
