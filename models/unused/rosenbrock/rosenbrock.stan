
data {
  // no data
}
parameters {
  real x1;
  real x2;
}

model {
  // Calculate surface directly
  target+= -(pow(100*(x2-x1*x1),2)+pow(1-x1,2))/20;
  x1~normal(0, 2);
}
