GLOBALS_SECTION
  #include <mvnorm.h>
  MVNORM_t neg_log_density;

DATA_SECTION
  init_number dim;
  init_matrix covar(1,dim,1,dim);
  init_vector x(1,dim);
  init_number lwr;
  init_number upr;

PARAMETER_SECTION
 init_vector mu(1,dim);
  //init_bounded_vector mu(1,dim,lwr,upr,1);
  objective_function_value jnll;

PRELIMINARY_CALCS_SECTION
  neg_log_density.setSigma(covar);

PROCEDURE_SECTION
    jnll=neg_log_density(x-mu);

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
