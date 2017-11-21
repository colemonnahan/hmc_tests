
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N); // number of years
  DATA_VECTOR(catches);
  DATA_VECTOR(logcpue);

 // bounded parameters provide uniform priors
  PARAMETER(logK);
  PARAMETER(logr);
  PARAMETER(iq);
  PARAMETER(isigma2);
  PARAMETER(itau2);
  PARAMETER_VECTOR(u_raw);

  Type sigma2;
  Type sigma;
  Type tau2;
  Type q;
  Type K;
  Type r;
  vector<Type> u(N);
  sigma2 = 1/isigma2;
  sigma=sqrt(sigma2);
  tau2 = 1/itau2;
  q = 1/iq;
  K = exp(logK);
  r = exp(logr);
  // non-centered random effects used below; implies u~N(0,sigma)
  u=u_raw*sigma;

  Type nll=0.0;
  Type nlp=0.0;

  vector<Type> B(N);
  vector<Type> ypred(N);
  Type temp;
  // priors; note the parameterization is different between stan, where
  // scale is 1/scale here.
  nlp-=dnorm(logr, Type(-1.38), Type(0.51), true);
  nlp-=dgamma(iq, Type(0.001), Type(1/0.001), true);
  nlp-=dgamma(isigma2, Type(3.785518), Type(1/0.010223), true);
  nlp-=dgamma(itau2, Type(1.708603), Type(1/0.008613854), true);
  // project dynamics
 B(0) = K;
 ypred(0) = log(B(0)) +log(q);


 for(int i=1; i<=(N-1); i++){
   temp = (B(i-1)+r*B(i-1)*(Type(1.)-B(i-1)/K)-catches(i-1))*exp(u(i));
   if(temp<.001){
     //increment_log_prob(-1*(temp-1)^2);
     B(i) = Type(1)/(Type(2)-temp/Type(.001));
   } else {
     B(i) = temp;
   }
    B(i) = temp;
   ypred(i) = log(B(i)) +log(q);
 }

 // likelihoods
 nll-=dnorm(u_raw, Type(0.0), Type(1.0), true).sum();
 // The likelihood
 nll-=dnorm(ypred, logcpue, sqrt(tau2), true).sum();

 Type nld=nll+nlp;
 return(nld);
}
