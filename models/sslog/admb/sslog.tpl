
DATA_SECTION
  init_number N;
  init_vector catches(1,N);
  init_vector logcpue(1,N);
  
PARAMETER_SECTION
  init_number logK;
  init_number logr;
  init_number iq;
  init_number isigma2;
  init_number itau2;
  init_vector u_raw(1,N);
  vector u(1,N);
  vector ypred(1,N);
  vector B(1,N);
  number temp;
 // exponentiated versions of variance parameters since ADMB can't do lower bounds
  number sigma2;
  number sigma;
  number tau2;
  number q;
  number K;
  number r;
  number nlp;
  number nll;
  
  objective_function_value nld;

INITIALIZATION_SECTION
  logK 6.5;
  logr -1;
  iq 5;
  isigma2 400;
  itau2 92;
  
PROCEDURE_SECTION
 nlp=0.0;
 nll=0.0;
  sigma2 = 1/isigma2;
  sigma=sqrt(sigma2);
  tau2 = 1/itau2;
  q = 1/iq;
  K = exp(logK);
  r = exp(logr);
  // non-centered random effects used below; implies u~N(0,sigma)
  u=u_raw*sigma;

   // priors; note the parameterization is different between stan, where
  // scale is 1/scale here.
  nlp+=dnorm(logr, -1.38, 0.51);
  nlp+=dgamma(iq, 0.001, 0.001);
  nlp+=dgamma(isigma2, 3.785518, 0.010223);
  nlp+=dgamma(itau2, 1.708603, 0.008613854);
  // project dynamics
 B(1) = K;
 ypred(1) = log(B(1)) +log(q);
 for(int i=2; i<=N; i++){
   temp = (B(i-1)+r*B(i-1)*(1.0-B(i-1)/K)-catches(i-1))*exp(u(i));
   if(temp<.001){
     //increment_log_prob(-1*(temp-1)^2);
     B(i) = 1/(2-temp/.001);
   } else {
     B(i) = temp;
   }
    B(i) = temp;
   ypred(i) = log(B(i)) +log(q);
 }

 // likelihoods
 nll+=dnorm(u_raw, 0.0, 1.0);
 // The likelihood
 // nll-=dnorm(ypred, logcpue, sqrt(tau2));
 dvariable SS=norm2(logcpue-ypred);
 dvariable tmp= N*(0.5*log(2.*M_PI)+log(sqrt(tau2)))+0.5*SS/tau2;
 nll+= tmp;

 nld=nll+nlp;
 
