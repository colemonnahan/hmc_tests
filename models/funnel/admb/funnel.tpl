
DATA_SECTION

PARAMETER_SECTION
  init_number v;
  init_number theta;
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0;
  jnll+= (v/3)*(v/3)/2;
  jnll+= (theta/exp(v))*(theta/exp(v))/2;

