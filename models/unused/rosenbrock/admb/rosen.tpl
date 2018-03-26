
DATA_SECTION

PARAMETER_SECTION
  init_number x1;
  init_number x2;
  objective_function_value jnll;

PROCEDURE_SECTION
  jnll=0;
  jnll+= (pow(100*(x2-x1*x1),2)+pow(1-x1,2))/20;
  // Prior
  jnll+= pow(x1/2, 2)/2;

