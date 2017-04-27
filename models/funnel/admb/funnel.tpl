
GLOBALS_SECTION
 #include "statsLib.h"
DATA_SECTION

PARAMETER_SECTION
  init_number v;
  init_number theta;
  objective_function_value jnll;

PROCEDURE_SECTION
 double v2=value(exp(v));
  jnll=0;
  jnll+= dnorm(v,0,1);
  jnll+= dnorm(theta,0, v2);

