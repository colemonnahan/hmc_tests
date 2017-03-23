data {
  matrix[3,3] covar;
  vector[3] x;
}
parameters {
  real<lower=-2, upper=2> mu1;
  real<lower=0> mu2;
  real mu3;
}

transformed parameters{
  vector[3] mu;
  mu[1]=mu1;
  mu[2]=mu2;
  mu[3]=mu3;
}

model {
  x~multi_normal(mu, covar);
}
