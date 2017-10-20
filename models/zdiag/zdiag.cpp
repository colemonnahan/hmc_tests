#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n); // dimension
  DATA_VECTOR(sds);
  PARAMETER_VECTOR(mu);
  Type nll= -1*dnorm(mu, Type(0.0), sds, true).sum();
  return(nll);
}
