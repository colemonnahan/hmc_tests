DATA_SECTION
  int n
 !! n=500;
  matrix M(1,n,1,n)
  matrix M1(1,n,1,n)
 LOC_CALCS
  random_number_generator rng(501);
  M.initialize();
  M1.initialize();
  M(1,1)=1.0;
  for (int i=2;i<=n;i++)
  {
    M(i,i)=5.0;
    M(i)(1,i-1).fill_randn(rng);
    M(i)/=norm(M(i));
  }
  M=M*trans(M);
  dvector ev=eigenvalues(M);
  cout << "M condition number " << max(ev)/min(ev) << endl;
  M1(1,1)=1.0;
  for (int i=2;i<=n;i++)
  {
    M1(i,i)=5.0;
    M1(i)(1,i-1).fill_randn(rng);
    M1(i)/=norm(M1(i));
  }
  M1=M1*trans(M1);
  dvector ev1=eigenvalues(M1);
  cout << "M1 condition number " << max(ev1)/min(ev1) << endl;
  //cout << eigenvalues(M) << endl;
  //ad_exit(1);
PARAMETER_SECTION
  init_vector x(1,n)
  sdreport_vector sx(1,n)
  objective_function_value f
PROCEDURE_SECTION
  dvariable n2x=x*(M*x);
  dvariable n4x=square(x*(M1*x));
  f+=0.5*n2x + 0.1*n4x;
  sx=x;
  if (mceval_phase())
   cout << sqrt(value(n2x)) << endl;
  
