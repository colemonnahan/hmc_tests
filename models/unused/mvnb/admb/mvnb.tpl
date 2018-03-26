GLOBALS_SECTION
  #include <mvnorm.h>
  MVNORM_t neg_log_density;

DATA_SECTION
  init_matrix covar(1,3,1,3);
  !! cout << covar << endl;

PARAMETER_SECTION
  init_bounded_number mu1(-2,2);
  init_bounded_number mu2(0,5);
  init_number mu3;
  objective_function_value jnll;

PRELIMINARY_CALCS_SECTION
  neg_log_density.setSigma(covar);
  // dvector x(1,3);
  // for(int i=1; i<=3; i++) x(i)=0;

PROCEDURE_SECTION
  dvar_vector mu(1,3);
  mu(1)=mu1-.15;
  // Manually transform paramter
  mu(2)=mu2+.1;
  mu(3)=mu3+.98;
  jnll=neg_log_density(mu);
 // jnll-=mu2;

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;

// // nuts hbf 0
// old scale= 0.937504 1 1
// current scale= 0.937504 1 1
// S before= 0.995545 0.266671 0.399997
//  0.266671 0.500024 0.250006
//  0.399997 0.250006 0.875
// Starting from theta= 0.510808 -2.49733e-005 0.499988
// Starting from z= 0.509669 0.136506 0.613011
// Starting from nll=3.34832
// Starting from chd= 0.99777 0 0
//  0.267267 0.65467 0
//  0.400891 0.218218 0.816497

// // nuts hbf 1
// old scale= 0.937504 1 1
// current scale= 0.937504 1 1
// S before= 0.995545 0.266671 0.399997
//  0.266671 0.500024 0.250006
//  0.399997 0.250006 0.875
// Starting from theta= 0.510808 -2.49733e-005 0.499988
// Starting from z= 0.509669 0.136506 0.613011
// Starting from nll=3.34832
// Starting from chd= 0.99777 0 0
//  0.267267 0.65467 0
//  0.400891 0.218218 0.816497
