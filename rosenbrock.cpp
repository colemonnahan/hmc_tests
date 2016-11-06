#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(x1);
  PARAMETER(x2);
  Type nll= (pow(100*(x2-x1*x1),2)+pow(1-x1,2))/20;
  return nll;
}

