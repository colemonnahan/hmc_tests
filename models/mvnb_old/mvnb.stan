data {
  int<lower=0> Npar;
  matrix[Npar,Npar] covar;
  vector[Npar] x;
}
parameters {
  vector<lower=-2, upper=2>[Npar] mu;
}

model {
  x~multi_normal(mu, covar);
}
