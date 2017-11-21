
// this is a hack version I started to change priors then gave up

// par(mfrow=c(2,3))
// isigma2 <- rgamma(1e5, 3.78, .01)
// hist(1/sqrt(isigma2), breaks=750, xlim=c(0,.2))
// sigma <- rlnorm(1e5, -2.81, .25)
// hist(sigma, breaks=750, xlim=c(0,.2))
// qqplot(1/sqrt(isigma2), sigma);  abline(0,1)
// itau2 <- rgamma(1e5, 1.708603, 0.008613854)
// hist(1/sqrt(itau2), breaks=750, xlim=c(0,1))
// tau <- rlnorm(1e5, -2.5, .75)
// hist(tau, breaks=750, xlim=c(0,1))
// qqplot(1/sqrt(itau2), tau);  abline(0,1)
// iq2 <- rgamma(1e5, .001, 0.001)
// hist(1/sqrt(iq2), breaks=750, xlim=c(0,1000))
// q <- rlnorm(1e5, -2.5, .75)
// hist(q, breaks=750, xlim=c(0,1))
// qqplot(1/sqrt(iq2), q);  abline(0,1)

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
