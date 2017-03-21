
DATA_SECTION
  init_number Nfish;
  init_number Nobs;
  init_vector loglengths(1,Nobs);
  init_ivector fish(1,Nobs);
  init_ivector ages(1,Nobs);

PARAMETER_SECTION
  //init_bounded_number delta(0,5);
  init_bounded_number sigma_obs(0,5);
  init_number logLinf_mean;
  init_number logk_mean;
  init_bounded_number logLinf_sigma(0,5);
  init_bounded_number logk_sigma(0,5);
  init_vector logLinf(1,Nfish);
  init_vector logk(1,Nfish);
  objective_function_value jnll;

INITIALIZATION_SECTION
 //delta 1.1;
 sigma_obs .11;
 logLinf_mean 3.1;
 logk_mean -2.1;
 logLinf_sigma .11;
 logk_sigma .11;
 logLinf 3;
 logk -2;

PROCEDURE_SECTION
  jnll=0;
  dvar_vector ypred(1,Nobs);
  dvariable Linf;
  dvariable k;
  double delta=1.0;
  // Priors
  jnll+= pow((sigma_obs-.1)/.1, 2);
  jnll+= pow((logLinf_sigma-.1)/.1, 2);
  jnll+= pow((logk_sigma-.2)/.1, 2);

  // Random effects
  for(int i=1; i<=Nfish; i++){
     jnll+= pow((logLinf_mean-logLinf(fish(i)))/logLinf_sigma, 2);
     jnll+= pow((logk_mean-logk(fish(i)))/logk_sigma, 2);
  }

  // Predict fish length and add likelihood
  for(int i=1; i<=Nobs; i++){
     Linf=exp(logLinf(fish(i)));
     k=exp(logk(fish(i)));
     ypred(i)=pow(log( Linf*(1-exp(-k*(ages(i)-5.0)))), delta);
     jnll+=pow( (loglengths(i)-ypred(i))/sigma_obs,2);
  }

REPORT_SECTION
