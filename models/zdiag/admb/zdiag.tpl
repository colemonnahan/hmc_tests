DATA_SECTION
  init_number n;
  init_vector x(1,n);
  init_vector sds(1,n);

PARAMETER_SECTION
  init_vector mu(1,n);
  objective_function_value jnll;

PROCEDURE_SECTION
 jnll=0.0;
 for(int i=1; i<=n; i++){
   jnll+=dnorm(mu(i), 0.0, sds(i));
 }


