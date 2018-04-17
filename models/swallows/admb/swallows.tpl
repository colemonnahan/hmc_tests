GLOBALS_SECTION
 #include <admodel.h>
 #include "statsLib.h"
 #include <df1b2fun.h>
 #include <adrndeff.h>
 // Negative log-likelihood for the Cauchy distribution; scalar arguments.
 dvariable dcauchy(const prevariable& mu, const double& x, const double& scale)
 {
  #ifdef M_PI
  double pi = M_PI;
  #else
  double pi = 3.14159265358979323844;
  #endif
  dvariable logres=0.0;
  logres+=log(pi);
  logres+= log(scale);
  logres+= log(1 + pow( (x-mu)/scale, 2));
  return logres;
 }
 dvariable inv_logit(const dvariable& x){
  dvariable y= 1/(1+mfexp(-x));
  return y;
 }


DATA_SECTION
  init_number I;
  init_number K;
  init_int nfam;
  // For some reason ADMB reads this in the wrong order. So just go with it
  // and transpose
  init_matrix CH0(1,K,1,I);
  matrix CH(1,I,1,K);
  !! CH = trans(CH0);
  init_vector carez(1,I);
  init_ivector year(1,I);
  init_vector agec(1,K);
  init_ivector family(1,I);
  init_ivector last(1,I);
  init_imatrix ones(1,I,1,K);
  init_ivector ones2(1,I);

PARAMETER_SECTION
  // fixed effects
  init_number sigmayearphi;
  init_number sigmaphi;
  init_number sigmap;
  init_vector a(1,17);
  init_number a1;
  init_vector b0(1,4);
  init_vector b1(1,4);
  // non-centered random effects
  init_vector fameffphi_raw(1,72);
  init_vector fameffp_raw(1,72);
  init_vector yeareffphi_raw(1,4);

  // exponentiated versions of variance parameters since ADMB can't do lower bounds
  number sigmayearphi2;
  number sigmaphi2;
  number sigmap2;
  number x;
  number nlp;
  number nll;
  matrix p(1,I,1,K);
  matrix phi(1,I,1,K-1);
  matrix chi(1,I,1,K+1);
  objective_function_value nld; // negative log density

INITIALIZATION_SECTION
 sigmaphi 1
 sigmayearphi 1
 sigmap 1
 a .1
 a1 .1
 b0 .1
 b1 .1
 fameffphi_raw .1
 fameffp_raw .1
 yeareffphi_raw .1

PROCEDURE_SECTION
 nlp=0.0; // negative log prior
 nll=0.0; // negative log likelihood

// Exponentiate sds since ADMB can't do (0,Inf) bounds. Need to adjust for
 // Jacobian below
 sigmayearphi2=exp(sigmayearphi);
 sigmaphi2=exp(sigmaphi);
 sigmap2=exp(sigmap);

 int k;
 for(int i=1; i<=I; i++){ // loop over each individual
  // calculate phi as a function of fixed and random effects
  for(int t=1; t<=(K-1); t++) {
    x=a(t)+ a1*carez(i)+
      sigmayearphi2*yeareffphi_raw(year(i))+
      sigmaphi2*fameffphi_raw(family(i));
    phi(i,t) = inv_logit(x);
  }
  // calculate p as a function of fixed and random effects
  p(i,1) = 1;  // first occasion is marking occasion
  for(int t=2; t<=K; t++){
    x=b0(year(i))+ b1(year(i))*agec(t)+
      sigmap2*fameffp_raw(family(i));
    p(i,t) = inv_logit(x);
  }
  // probabilitiy of never being seen after last observation. ind here is
  // a reverse index so this loop goes from K:2, recursively calculating
  // backward.
  chi(i,K+1) = 1.0;
  k=K;
  while (k > 1) {
    chi(i,k) = (1-phi(i,k-1)) + phi(i,k-1) * (1-p(i,k)) * chi(i,k+1);
    k=k-1;
  }
  chi(i,1) = (1-p(i,1)) * chi(i,2);
 }

 // Jacobian adjustment for variances
 nll -= sigmaphi + sigmayearphi + sigmap;

 // priors
 nlp+= dcauchy(sigmaphi2, 0.0, 1.0);
 nlp+= dnorm(sigmayearphi2, 0.0, 3);
 nlp+= dcauchy(sigmap2, 0.0, 1.0);
 nlp+= dnorm(a, 0.0, 1.5);
 nlp+= dnorm(a1, 0.0, 5.0);
 nlp+= dnorm(b0, 0.0, 5.0);
 nlp+= dnorm(b1, 0.0, 5.0);

 // random effects; non-centered
 nll+= dnorm(fameffphi_raw, 0,1);
 nll+= dnorm(fameffp_raw, 0,1);
 nll+= dnorm(yeareffphi_raw, 0,1);

 // // likelihood
 for(int i=1; i<=I; i++){ // loop over each individual
    // probability of survival, known alive since k<last
    for (int t=2; t<=last(i); t++) {
      // ones[i,t]~bernoulli(phi[i, t-1]);
    	nll-= log(phi(i,t-1));
    }
    // // probability of observation given known alive
    for(int t=1; t<=last(i); t++){
      // CH[i,t]~bernoulli(p[i,t]);
      if(CH(i,t)==1){
    	nll-= log(p(i,t));
      } else {
    	nll-= log(1-p(i,t));
      }
    }
    // probability of no observations after time period last
    // ones2[i]~bernoulli(chi[i,last[i]+1]);
     nll-= log(chi(i,last(i)+1));
  }
 nld=nll+nlp; // negative log density

