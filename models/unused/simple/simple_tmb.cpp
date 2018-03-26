


// init_int nobs
//   init_vector Y(1,nobs)
//   init_vector x(1,nobs)
// PARAMETER_SECTION
//   init_number a
//   init_number b
//   vector pred_Y(1,nobs)
//   sdreport_number aa
//   objective_function_value f
// PROCEDURE_SECTION
//  aa=a;
//   pred_Y=a*x+b;
//   f=(norm2(pred_Y-Y));
//   f=nobs/2.*log(f);    // make it a likelihood function so that
//                        // covariance matrix is correct

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(x);
  DATA_INTEGER(nobs);
  PARAMETER(a);
  PARAMETER(b);
  vector<Type> pred_Y(nobs);
  Type nll = -sum(dnorm(Y, a*x+b, Type(1.0), true));
  return nll;
}
