// the centered version of the 8 schools model

#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_VECTOR(sigma);
  PARAMETER(mu);
  PARAMETER(logtau);
  PARAMETER_VECTOR(theta);
  Type tau=exp(logtau);
  //ADREPORT(2*b);
  Type nll=0;
  //vector<Type> theta(8);
  for(int i=0;i<=7;i++){
   // theta(i)=mu+tau*eta(i);
  // likelihood
      nll-=dnorm(Y(i), theta(i),sigma(i),true);
  // hyperdistribution
     nll-=dnorm( theta(i),mu, tau,true);
    //nll-=dnorm(theta(i), Type(0), Type(1), true);
  }
  REPORT(nll);
  REPORT(theta);
  return nll;
}

