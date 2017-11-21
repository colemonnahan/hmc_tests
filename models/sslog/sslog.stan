data {
  int<lower=0> N; // number of years
  real catches[N];
  real logcpue[N];
}
parameters {
 // bounded parameters provide uniform priors
  real logK;
  real logr;
  real iq;
  real isigma2;
  real itau2;
  vector[N] u_raw;
}

transformed parameters {
  real sigma2;
  real sigma;
  real tau2;
  real q;
  real K;
  real r;
  vector[N] u;
  sigma2 = 1/isigma2;
  sigma=sqrt(sigma2);
  tau2 = 1/itau2;
  q = 1/iq;
  K = exp(logK);
  r = exp(logr);
  // non-centered random effects used below; implies u~N(0,sigma)
  u=u_raw*sigma;
}

model {
 real B[N];
 real ypred[N];
 real temp;
 // priors
 logr~normal(-1.38, 0.51);
 iq~gamma(0.001, 0.001);
 isigma2~gamma(3.785518, 0.010223);
 itau2~gamma(1.708603, 0.008613854);
 // project dynamics
 B[1] = K;
 ypred[1] = log(B[1]) +log(q);


for(i in 2:N){
   temp = (B[i-1]+r*B[i-1]*(1-B[i-1]/K)-catches[i-1])*exp(u[i]);
   if(temp<.001){
  //increment_log_prob(-1*(temp-1)^2);
   B[i] = 1/(2-temp/.001);
   } else {
      B[i] = temp;
   }
   ypred[i] = log(B[i]) +log(q);
 }
 // hyper prior
 u_raw~normal(0, 1);
 // The likelihood
 logcpue~normal(ypred, sqrt(tau2));
}
