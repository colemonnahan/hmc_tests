data {
  int<lower=0> n;
  vector[n] x;
  vector[n] sds;
}
parameters {
  vector[n] mu;
}

model {
  x~normal(mu,sds);
}

