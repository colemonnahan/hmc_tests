#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(v);
  PARAMETER(theta);
  Type nll= 0;
  // Add arbitrary prior on x1 to keep it a little tigher
  nll-=dnorm(v, Type(0.0), Type(1.0), true);
  nll-=dnorm(theta, Type(0.0), exp(v), true);
  return nll;
}

