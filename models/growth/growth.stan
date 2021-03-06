	data {
  int<lower=0> Nfish; // number of groups
  int<lower=0> Nobs;   // number of observations
  vector[Nobs] loglengths;  // observed log cpue
  int fish[Nobs];      // vector to index fish
  int ages[Nobs];      // observed ages
}
parameters {
  // fixed effects
  real<lower=0> delta;
  real<lower=0> sigma_obs; // data on log scale

  // hyperparameters with bounds
  real logLinf_mean;
  real logk_mean;
  real<lower=0> logLinf_sigma;
  real<lower=0> logk_sigma;

  // non-centered random effects
  vector[Nfish] logLinf_raw;
  vector[Nfish] logk_raw;
}

transformed parameters {
  // non-centered random effects, implies logk~N(logk_mean, logk_sigma) etc.
  vector[Nfish] logLinf;
  vector[Nfish] logk;
  logLinf = logLinf_mean+logLinf_raw*logLinf_sigma;
  logk = logk_mean+logk_raw*logk_sigma;
}

model {
  vector[Nobs] ypred;
  real k;

  // priors
  delta~normal(1,.25);
  sigma_obs~cauchy(0,5);

  // hyperpriors
  logLinf_sigma~cauchy(0,5);
  logk_sigma~cauchy(0,5);
  // hyper means are uniform above

  // random effects; non-centered
  logLinf_raw~normal(0, 1);
  logk_raw~normal(0, 1);

  // calculate likelihood of data
  for(i in 1:Nobs){
    k = exp(logk[fish[i]]);
    ypred[i] = logLinf[fish[i]] + delta*log(1-exp(-k*(ages[i]-5)));
  }
  loglengths~normal(ypred, sigma_obs);
}
