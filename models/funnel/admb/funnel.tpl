
DATA_SECTION

PARAMETER_SECTION
  init_number v;
  init_number theta;
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0;
  jnll+= pow(v/3, 2);
  jnll+= pow(theta/(exp(v)), 2);

