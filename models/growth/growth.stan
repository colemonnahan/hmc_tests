data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  vector[Nobs] loglengths;  // observed log cpue
  int fish[Nobs];      // vector to index fish
  int ages[Nobs];      // observed ages
}
parameters {
  // fixed effects
  // real<lower=0, upper=5> delta;
  real<lower=0, upper=5> sigma_obs; // data on log scale

  // hyperparameters with bounds
  real logLinf_mean;
  real logk_mean;
  real<lower=0, upper=5> logLinf_sigma;
  real<lower=0, upper=5> logk_sigma;

  // random effects
  vector[Nfish] logLinf;
  vector[Nfish] logk;
}

model {
  vector[Nobs] ypred;
  real Linf;
  real k;
  real delta=1.0;

  // priors
  // delta is uniform above
  // sigma_obs~cauchy(0,5);
  sigma_obs~normal(.1,.1);
  // hyperpriors
  // logLinf_sigma~cauchy(0,5);
  // logk_sigma~cauchy(0,5);
  logLinf_sigma~normal(.1, .1);
  logk_sigma~normal(.2, .1);
  // hyper means are uniform above

  // random effects
  logLinf~normal(logLinf_mean, logLinf_sigma);
  logk~normal(logk_mean, logk_sigma);

  // calculate likelihood of data
   for(i in 1:Nobs){
    Linf = exp(logLinf[fish[i]]);
    k = exp(logk[fish[i]]);
    ypred[i] = log( Linf*(1-exp(-k*(ages[i]-5)))^delta );
   }
  loglengths~normal(ypred, sigma_obs);
}
