DATA_SECTION

PARAMETER_SECTION
  init_bounded_number x1(-2,1);
  init_bounded_number x2(-1,2);
  objective_function_value jnll;

PROCEDURE_SECTION
   jnll= x1*x1/2+x2*x2/2;

