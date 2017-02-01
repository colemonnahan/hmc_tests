DATA_SECTION
  init_number dim;

PARAMETER_SECTION
  init_number x1;
  init_number x2;
  objective_function_value jnll;

PROCEDURE_SECTION
   jnll= x1*x1/2+x2*x2/2;

