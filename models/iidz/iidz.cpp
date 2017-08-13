#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n); // dimension
  DATA_VECTOR(x);
  PARAMETER_VECTOR(mu);
  Type nll= -1*dnorm(mu, x, Type(1.0), true).sum();
  return(nll);
}
