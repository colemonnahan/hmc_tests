data {
  int<lower=0> n;
  vector[n] x;
}
parameters {
  vector[n] mu;
}

model {
  x~normal(mu,1.0);
}
