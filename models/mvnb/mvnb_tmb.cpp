#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(Npar); // dimension of matrix
  DATA_MATRIX(covar);
  DATA_VECTOR(bounds);
  PARAMETER_VECTOR(mu); // parameters
  vector<Type> mu_bounded(2);
  Type nll=0;
  vector<Type> log_scales(2);
  for(int i=0; i<=1; i++){
    mu_bounded(i)=bounds(0)+(bounds(1)-bounds(0))/(1+exp(-mu(i)));
    // jacobian adjustment by hand
	log_scales(i) = log( (bounds(1)-bounds(0)) ) -mu(i) -2*log(1+exp(-mu(i))) ;
	nll-=log_scales(i);
  }
  using namespace density;
  MVNORM_t<Type> neg_log_density(covar);
  nll += neg_log_density(mu_bounded);
  REPORT(mu_bounded);
  REPORT(log_scales);
  return(nll);
}
