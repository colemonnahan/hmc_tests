// Spatial poisson GLMM on a grid, with exponentially decaying correlation function
#include <TMB.hpp>
// Hand-coded Cauchy distribution
template<class Type>
Type dcauchy(Type x, Type mean, Type shape, int give_log=0){
  Type logres = 0.0;
  logres-= log(M_PI);
  logres-= log(shape);
  // Note, this is unstable and should switch to log1p formulation
  logres-= log(1 + pow( (x-mean)/shape ,2));
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_VECTOR(y);
  DATA_MATRIX(X)
  DATA_MATRIX(dd)
  PARAMETER_VECTOR(b);
  PARAMETER(a);
  PARAMETER(sigma);
  PARAMETER_VECTOR(u);

  using namespace density;
  int i,j;
  Type res=0;

  vector<Type> eta(n);
  eta = X*b + sigma*u;

  //
  matrix<Type> cov(n,n);
  for (i=0;i<n;i++)
  {
    cov(i,i)=Type(1);
    for ( j=0;j<i;j++)
    {
      cov(i,j)=exp(-a*dd(i,j));			// Exponentially decaying correlation
      cov(j,i)=cov(i,j);
    }
  }

  MVNORM_t<Type> neg_log_density(cov);
  res+=neg_log_density(u);

  // logdpois = N log lam - lam
  for(i=0;i<n;i++) res -= y[i]*eta[i]-exp(eta[i]);

  // priors; a is U(0,2) bounded externally but putting a normal here too
  res-= dnorm(a, Type(1.0), Type(.35), true);
  res-= dnorm(b, Type(0.0), Type(10.0), true).sum();
  res-= dcauchy(sigma, Type(0.0), Type(1.0), true);

  return res;

}
