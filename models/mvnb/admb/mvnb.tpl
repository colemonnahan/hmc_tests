GLOBALS_SECTION
  #include <mvnorm.h>
  MVNORM_t neg_log_density;

DATA_SECTION
  init_matrix covar(1,3,1,3);

PARAMETER_SECTION
  init_bounded_number mu1(-2,2);
  init_number mu2;
  init_number mu3;
  objective_function_value jnll;

PRELIMINARY_CALCS_SECTION
  neg_log_density.setSigma(covar);
  // dvector x(1,3);
  // for(int i=1; i<=3; i++) x(i)=0;

PROCEDURE_SECTION
  dvar_vector mu(1,3);
  mu(1)=mu1;
  // Manually transform paramter
  mu(2)=exp(mu2);
  mu(3)=mu3;
  jnll=neg_log_density(mu);
  jnll-=mu2;

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
