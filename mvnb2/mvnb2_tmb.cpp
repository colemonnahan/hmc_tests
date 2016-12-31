#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_MATRIX(covar);
  PARAMETER_VECTOR(mu); // parameters
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  Type nll= neg_log_density(mu);
  REPORT(nll);
  return(nll);
}
