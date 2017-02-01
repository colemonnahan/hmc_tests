#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  PARAMETER(x1);
  PARAMETER(x2);
  Type nll=  x1*x1/2+x2*x2/2;
  return(nll);
}
