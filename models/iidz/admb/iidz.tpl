DATA_SECTION
  init_number n;

PARAMETER_SECTION
  init_vector mu(1,n);
  objective_function_value jnll;

PROCEDURE_SECTION
 jnll=dnorm(mu, 0, 1);

