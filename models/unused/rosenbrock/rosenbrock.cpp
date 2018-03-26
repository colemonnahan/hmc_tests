#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(x1);
  PARAMETER(x2);
  Type nll= (pow(100*(x2-x1*x1),2)+pow(1-x1,2))/20;
  // Add arbitrary prior on x1 to keep it a little tigher
  nll+=dnorm(x1, Type(0), Type(2), true);
  return nll;
}

