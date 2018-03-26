// Neal's funnel, copied from Betancourt and Girolami (2015) HMC for hierarchical models
data {
  // no data
}
parameters {
  real v; 			//  log variance
  real theta;			// random effect
}

model {
  theta~normal(0, exp(v));
  v~normal(0, 1);
}
