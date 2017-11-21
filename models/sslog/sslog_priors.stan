
// this is a hack version I started to change priors then gave up
data {
  int<lower=0> N; // number of years
  real catches[N];
  real logcpue[N];
}
parameters {
  real<lower=0> logK;
  real<lower=0> logr;
  real<lower=0> q;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[N] u_raw;
}

transformed parameters {
  vector[N] u;
  // non-centered random effects; implies u~N(0,sigma)
  u = u_raw*sigma;
}

model {
 real B[N];
 real ypred[N];
 real temp;

// priors
 logr~normal(-1.38, 0.51);
 iq~gamma(0.001, 0.001);
 sigma~lognormal(-2.8, 0.25);
 tau~lognormal(-2.5,0.75);
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
