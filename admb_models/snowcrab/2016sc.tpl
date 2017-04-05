DATA_SECTION
//Bering Sea snow crab model
  !!CLASS ofstream post("eval.csv")
  int call_no;
  !! call_no = 0;
  number spmo																						// spawning month
  !! spmo=0.;																							// spmo=deviation in fraction of year from time of fishery to mating 

  init_int styr  																							//start year of the model
  init_int endyr 																							//end year of the model

  // Data stuff only from here
  init_int nirec     																				    	 //number of intial recruitments to estimate 
  init_int nlenm      																				     //number of length bins for males in the model
  init_int nobs_fish 																					 //number of years of fishery retained length data
  init_ivector yrs_fish(1,nobs_fish) 															//years when have fishery retained length data
  init_matrix nsamples_fish(1,2,1,nobs_fish)											//nsamples weight for fish length comps needmatrix each year,new/old shell
  init_int nobs_fish_discf  																		//number of years of fishery female discard length data
  init_ivector yrs_fish_discf(1,nobs_fish_discf) 										//years when have fishery discard length data
  init_vector nsamples_fish_discf(1,nobs_fish_discf)							//nsamples weight for fish length comps needmatrix each year and sex
  init_int nobs_fish_discm 																		//number of years of fishery male discard length data
  init_ivector yrs_fish_discm(1,nobs_fish_discm)									//years when have fishery discard length data
  init_matrix nsamples_fish_discm(1,2,1,nobs_fish_discm)				 	//nsamples weight for fish length comps needmatrix each year and sex
  init_number nobs_trawl
  init_vector yrs_trawl(1,nobs_trawl)
  init_matrix nsamples_trawl(1,2,1,nobs_trawl)
  init_number nobs_srv1
  init_ivector yrs_srv1(1,nobs_srv1)            								             //years when have biomass estimates
  init_int nobs_srv1_length              												         //number of years of survey length data
  init_ivector yrs_srv1_length(1,nobs_srv1_length) 						       	//years when have length data
  init_4darray nsamples_srv1_length(1,2,1,2,1,2,1,nobs_srv1_length)	 //number of samples for each length comp by immat,mat,new/old,sex,year
   
//extra survey in 2009
  init_int yrs_srv2                     																 //years when have biomass estimates for extra survey
  init_3darray nsamples_srv2_length(1,2,1,2,1,2)									 //number of samples for each length comp by sruvey,sex,immat,mat
  init_4darray obs_p_srv2_lend(1,2,1,2,1,2,1,nlenm) 							 //numbers by length two surveys, sex, immat, mat
  init_vector  obs_srv2(1,2)
  init_matrix  obs_srv2_cv(1,2,1,2)
  
//extra survey in 2010
  init_int yrs_srv10                          															  //years when have biomass estimates for extra survey
  init_3darray nsamples_srv10_length(1,2,1,2,1,2)								 //number of samples for each length comp by sruvey,sex,immat,mat
  init_4darray obs_p_srv10_lend(1,2,1,2,1,2,1,nlenm)  							//numbers by length two surveys, sex, immat, mat
  init_vector  obs_srv10(1,2)
  init_matrix  obs_srv10_cv(1,2,1,2)
  
 //standard survey for length data
 //first index,1 immat, 2 mature,1 new shell, 2 old shell, then female 1 male 2
  init_5darray obs_p_srv1_lend(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm)  //immat,mat,new, old survey length data,female,male,year then bin
  init_3darray obs_p_fish_retd(1,2,1,nobs_fish,1,nlenm)  							 //new, old,male retained fishery length data
  init_matrix obs_p_fish_discfd(1,nobs_fish_discf,1,nlenm)     					 //female,discard length data
  init_3darray obs_p_fish_discmd(1,2,1,nobs_fish_discm,1,nlenm)           //male,discard length data new-old shel
  init_3darray obs_p_trawld(1,2,1,nobs_trawl,1,nlenm)        
  init_vector catch_numbers(styr,endyr)
  !! catch_numbers /= 1000.;                      												 //retained catch number of crab 
  init_vector catch_ret(styr,endyr)   															//retained catch lbs of crab 
  !! catch_ret /= 1000.;
  init_matrix catch_odisc(1,2,styr,endyr) 													 //estimated discard catch numbers female,male
  !! catch_odisc /= 1000.;
  init_vector catch_trawl(styr,endyr)                     										 //trawl bycatch numbers sex combined need to apply mort 80%
  !! catch_trawl /= 1000.;
  init_vector obs_srv1(1,nobs_srv1)  															  //survey numbers
  !! obs_srv1 /= 1000.;
  init_matrix cv_srv1o(1,2,1,nobs_srv1)   													 //survey cv
 
  init_vector maturity_logistic(1,nlenm)   													 //logistic maturity probability curve for new shell immature males
  init_matrix maturity_average(1,2,1,nlenm)  												 //probability mature for new immature females, avearge proportion mature by length for new males
  init_matrix maturity_old_average(1,2,1,nlenm) 										 //average proportion mature by length for old shell females, males
  init_3darray maturity(1,2,styr,endyr,1,nlenm)  											 //proportion mature by length for new shell females,males by year
  init_3darray maturity_old(1,2,styr,endyr,1,nlenm) 									 //proportion mature by length for old shell females, males by year
  init_matrix cv_mean_length_obs(1,2,1,2)  												 //cv of mean length for female and male min and max age
  init_vector length_bins(1,nlenm)
  init_vector catch_midptIn(styr,endyr)
  init_vector meanlength(1,13)  																	//mean length at age for F40 calc
  init_vector catch_fracmature(1,22)  														 //fraction morphometrically mature males in the pot fishery catch
  init_vector cpue(styr,endyr) 																		 //cpue fishery retained males only numbers/potlift
  init_vector catch_ghl(styr,endyr)
  init_int nobs_growf
  init_vector femalegrowdatx(1,nobs_growf)
  init_vector femalegrowdaty(1,nobs_growf)
  init_int nobs_growm
  init_vector malegrowdatx(1,nobs_growm)
  init_vector malegrowdaty(1,nobs_growm)

  // End of reading normal data file 
 // Open control file....
 !! ad_comm::change_datafile_name("2016sc.ctl");

  init_int styr_fut   																							  //start year of future projections
  init_int endyr_fut   																							  //end year of future projections
  init_int nsellen 																								 //selectivity is set to the selectivity at nselages-1 after age nselages
  init_int nsellen_srv1   																					    //same as above for survey selectivities
  init_number p_const
  init_number q1       																						   // Q  mult by pop biomass to get survey biomass
 
  init_vector M_in(1,2)   																					  	//natural mortality females then males
  init_vector M_matn_in(1,2)																				//natural mortality mature new shell female/male
  init_vector M_mato_in(1,2) 																				//natural mortality mature old shell female/male
  init_number m_disc        																					//fraction of discards that die (e.g 0.25)
  init_number m_trawl         																				//fraction of trawl discards that die(.8)

  init_vector median_rec(1,2) 																				//median recruitment value to use for last years in model
  init_int median_rec_yrs																					 //median recruitment fixed for endyr to endyr-median_rec_yrs+1
  init_number var_rec_obs
  init_number sd_var_rec

  init_vector sel_som(1,5)  																					 // parameters for somerton 2009 study selectivity curve
  init_number sel_avg_Nyrs
  
  init_number alpha_wt_imm_f
  init_number beta_wt_imm_f
  init_number alpha_wt_mat_f
  init_number beta_wt_mat_f
  init_number alpha_wt_m
  init_number beta_wt_m

  init_number linff_obs
  init_number sd_linff
  init_number linfm_obs
  init_number sd_linfm
  init_number growthkf_obs
  init_number sd_growthkf
  init_number growthkm_obs
  init_number sd_growthkm
  init_number af_obs
  init_number sd_af
  init_number am_obs
  init_number sd_am
  init_number bf_obs
  init_number sd_bf
  init_number bm_obs
  init_number sd_bm
  init_number a1_obs
  init_number sd_a1
  init_number b1_obs
  init_number sd_b1
  init_number sd_meetpt
  init_number var_last_obs
  init_number sd_var_last
  init_number mate_ratio

  init_number fraction_new_error
  init_int nages     																									  //number of ages to track for mature old shell 
  init_int matest_n
  init_int matestm_n

  init_vector wt_like(1,8)																							 //weights for selectivity likelihoods 1 fishery female, 2 survey female, 3 fishery male, 4 survey male
  init_vector like_wght(1,7) 																						 //likelihood weights for fishery length data, survey length, age data, catch likelihood, survey biomass likelihood,growth like
  init_number like_wght_mbio 																				 //likelihood weight for male biomass fit
  init_number like_wght_rec
  init_number like_wght_recf 
  init_number like_wght_sexr
  init_number like_wght_sel50
  init_number like_wght_fph1
  init_number like_wght_fph2
  init_number like_wght_fdev
  init_number wght_total_catch																				  //weight for total catch biomass
  init_number wght_female_potcatch  																	 //weight for female pot bycatch
  init_number cpue_cv          																					 //cv for fit to fishery pot cpue
  init_number wt_lmlike

  init_number old_shell_constraint
  init_vector disclen_mult(1,2)
  init_number selsmo_wght
  init_number femSel_wght
  
  init_number init_yr_len_smooth_f
  init_number init_yr_len_smooth_m 
  init_number natm_mult_wght 
  init_number natm_mult_var 
  init_number natm_Immult_wght 
  init_number natm_Immult_var 
  init_number smooth_mat_wght 
  init_number mat_est_wght
  init_number mat_est_vs_obs_sd 
  init_number smooth_mat_wght_f 
  init_number mat_est_vs_obs_sd_f 
  init_number growth_data_wght_m 
  init_number growth_data_wght_f 
  init_number extra_wght_ind_m 
  init_number smooth_disc_catch 
  init_number disc_catch_wght_fem 
  
  init_number fmort_phase
  init_number rec_phase
  init_number growth_phase
  init_number growth_phase2
  init_number maturity_phase
  init_number natM_phase
  
  init_int phase_moltingp																								 //phase to estimate molting prob for mature males
  init_int phase_fishsel																								 //phase to estimate dome shape parameters for fishery selectivities
  init_int survsel_phase  																								 //switch for which survey selectivty to use for 1989 to present - positive estimated negative fixed at somerton and otto
  init_int survsel1_phase																								 //switch for fixing all survey sel to somerton and otto - <0 fix, >0 estimate
  init_int phase_fut  																										 //phase to do F40% and future projection calculations
  init_int phase_logistic_sel 																						//phase to estimate selectivities using logistic function
  init_int phase_selcoffs																								//phase to estimate smooth selectivities

  init_int growth_switch 																								 //switch for which growth function to use
  init_int somertonsel
  init_int monot_sel 																										 //switch for monotonically increasing selectivities (1 on 0 off)
  init_int monot_sel_srv1																								 //sames as above for survey
  init_number maturity_switch
  init_int f_penalties
  init_int retro_years

   int styr_rec;   
   //year,age,sex
  int i;
  int j; 
  int k;   

// counters
   int iy;
   int ii;
   int m;
   int iage;
   int ilen;
   int nyl;
   int ll;
   int l;
   int il2;
   int is;
   int what;
   int ipass;
   int Nproj;
   !! Nproj = 101;
   int IsB0;
 //==============================================================================  
 LOCAL_CALCS
   cout<<"to local calcs"<<endl;
   styr_rec=styr-nirec; 																									  //year to start estimating recruits to get initial age comp
   if(nsellen>nlenm) nsellen=nlenm; 																				 //make sure nselages not greater than nages
   if(nsellen_srv1>nlenm) nsellen_srv1=nlenm; 															 //same as above for survey
   obs_srv1=obs_srv1*1000000;              																		 //survey numbers read in are millions of crab
   obs_srv2=obs_srv2*1000000;              																		  //survey numbers read in are millions of crab
   catch_ret=catch_ret/2204.6; 
   //   cout<<"end of local calcs"<<endl;
 END_CALCS

//==============================================================================
INITIALIZATION_SECTION
  mean_log_rec1 13.1
  rec_dev_mean  0.0
  log_avg_fmort -1.0
  log_avg_fmortt  -7.0
  log_avg_sel50_mn  4.4427  //this is 85.0 mm
  fmort_dev 0.00001
  moltp_af 2.0
  moltp_bf 150.
  moltp_am 0.02
  moltp_bm 300.
  moltp_ammat 0.05
  moltp_bmmat 105.0
  srv1_q 1.0
  srv2_q 1.0
  srv3_q 1.0
  srv1_sel95 60
  srv1_sel50 40 
  fish_fit_sel50_mn 95.1
  fish_fit_sel50_mo 100.0
  log_sel50_dev_mn  0.0000
  discard_mult 1.0
  mateste -0.7
  matestfe -0.7
  cpueq 0.001
  fish_disc_sel50_f 4.2
 //==============================================================================
PARAMETER_SECTION
 //growth pars
  init_bounded_number af(-100,0.0,growth_phase)
  init_bounded_number am(-50.0,0.0,growth_phase)
  init_bounded_number bf(1.0,10.0,growth_phase)
  init_bounded_number bm(1.0,5.0,growth_phase)
  init_bounded_number b1(1.0,1.5,growth_phase2)
  init_bounded_number bf1(1.0,2.0,growth_phase2)
  init_bounded_number deltam(10.0,50.0,growth_phase2)
  init_bounded_number deltaf(5.0,50.0,growth_phase2)
  init_bounded_number st_gr(0.5,0.5,-growth_phase2)
  init_bounded_vector growth_beta(1,2,0.749,0.751,-2)
  
 //Maturity parameters
   init_bounded_number matest50f(30,60,-4)      
   init_bounded_number matestslpf(0.05,0.5,-4)         
   init_bounded_number matest50m(60,110,-4)      
   init_bounded_number matestslpm(0.02,0.5,-4)    
   init_bounded_vector mateste(1,matestm_n,-6.0,-0.000000001,maturity_phase)
   init_bounded_vector matestfe(1,matest_n,-6.0,-0.000000001,maturity_phase)
   
   //Molting parameters
   //Females molt to maturity, therefore the prob of molting for mature females is 0
   init_bounded_number moltp_af(0.01,3.0,-6) 											     //paramters for logistic function molting
   init_bounded_number moltp_bf(20,200,-6)    											     //female
  init_bounded_number moltp_am(0.04,3.0,-5)  										    //paramters for logistic function molting
  init_bounded_number moltp_bm(130.0,300.0,-5)   									    //immature males
  init_bounded_number moltp_ammat(.0025,3.0,phase_moltingp) 			    //logistic molting prob for mature males
  init_bounded_number moltp_bmmat(1,120,phase_moltingp)  				   //logistic molting prob for mature males

    //recruitment parameters  (CSS: why two mean log recs?)
    init_number mean_log_rec1(1)
    init_bounded_number rec_dev_mean(0,0,-1)
    vector mean_log_rec(1,2)
    init_bounded_dev_vector rec_devf(styr,endyr,-15,15,rec_phase)
    init_bounded_dev_vector rec_devm(styr,endyr,0,0,-rec_phase)
    init_bounded_number alpha1_rec(11.49,11.51,-8)
    init_bounded_number beta_rec(3.99,4.01,-8)
	
	//initial numbers at length 
    init_bounded_matrix mnatlen_styr(1,2,1,nlenm,-3,15,1)
    init_bounded_matrix fnatlen_styr(1,2,1,12,-10,15,1)
	
	//Fishing mortality parameters
  init_number log_avg_fmort(1)
  init_bounded_dev_vector fmort_dev(styr,endyr-1,-5,5,fmort_phase)
  init_bounded_number log_avg_fmortdf(-8.0,-0.0001,fmort_phase)
  init_bounded_dev_vector fmortdf_dev(styr,endyr-1,-15,15,fmort_phase)
  init_bounded_number log_avg_fmortt(-8.0,-0.0001,fmort_phase)
  init_bounded_dev_vector fmortt_dev_era1(styr,1991,-15,15,fmort_phase)
  init_bounded_dev_vector fmortt_dev_era2(1992,endyr-1,-15,15,fmort_phase)
  init_bounded_number discard_mult(0.999,1.01,-fmort_phase)
  init_bounded_number Fem_F_prop_constant(0.0000001,1,fmort_phase)
  
  //Fishing selectivity parameters 
  //fish_slope, fish_fit, fish_disc, mn vs. mo, mo vs. mo2????
  init_bounded_number log_avg_sel50_mn(4,5.0,phase_logistic_sel)
  init_bounded_dev_vector log_sel50_dev_mn(1978,endyr,-5,5,-phase_logistic_sel+2)
  init_bounded_number log_avg_sel50_mo(4,5.0,-phase_logistic_sel)
  init_bounded_dev_vector log_sel50_dev_mo(1978,endyr,-5,5,-phase_logistic_sel+2)
  init_bounded_number fish_slope_mn(0.1,0.5,phase_logistic_sel)
  init_bounded_number fish_slope_mo(0.1,0.8,-phase_logistic_sel)

  init_bounded_number fish_fit_slope_mn(.05,.5,phase_logistic_sel)
  init_bounded_number fish_fit_sel50_mn(85.0,120.0,phase_logistic_sel)
  init_bounded_number fish_fit_slope_mo(.05,.5,-phase_logistic_sel)
  init_bounded_number fish_fit_sel50_mo(90.0,300.0,-phase_logistic_sel)

  init_bounded_number fish_slope_mo2(1.9,2.0,phase_fishsel)
  init_bounded_number fish_sel50_mo2(159.0,160.0,phase_fishsel)
  init_bounded_number fish_slope_mn2(.01,2.0,phase_fishsel)
  init_bounded_number fish_sel50_mn2(100.0,160.0,phase_fishsel)

  init_bounded_number fish_disc_slope_f(0.10,0.70,phase_logistic_sel)
  // init_bounded_number fish_disc_sel50_f(4.20,4.201,-phase_logistic_sel)
  init_bounded_number fish_disc_sel50_f(1,5,phase_logistic_sel)
  init_bounded_dev_vector log_dev_50f(1995,endyr-1,-5,5,-phase_logistic_sel)
  init_bounded_number fish_disc_slope_tf(.01,.3,phase_logistic_sel)
  init_bounded_number fish_disc_sel50_tf(30.,120.0,phase_logistic_sel)

  //Survey selectivity parameters
 // three survey eras:
  //srv1: 1978-1982
   init_bounded_number srv1_q(0.2,1.000,survsel1_phase+1)
   init_bounded_number srv1_q_f(0.2,1.000,survsel1_phase+1)
   init_bounded_number srv1_sel95(30.0,150.0,survsel1_phase)
   init_bounded_number srv1_sel50(0.0,150.0,survsel1_phase)
   
   //srv2: 1983-1988
   init_bounded_number srv2_q(0.2,1.000,survsel1_phase+1)
   init_bounded_number srv2_q_f(0.2,1.000,survsel1_phase+1)
   init_bounded_number srv2_sel95(50.0,160.0,survsel1_phase)
   init_bounded_number srv2_sel50(0.0,80.0,survsel1_phase)
   
   //srv3: 1989-present
   //males and females split
   init_bounded_number srv3_q(0.20,1.0000,survsel1_phase)
   init_bounded_number srv3_sel95(40.0,200.0,survsel1_phase)
   init_bounded_number srv3_sel50(25.0,90.0,survsel1_phase)
   init_bounded_number srv3_q_f(0.20,1.0000,survsel1_phase)
   init_bounded_number srv3_sel95_f(40.0,150.0,survsel1_phase+1)
   init_bounded_number srv3_sel50_f(0.0,90.0,survsel1_phase+1)
   
   //2009 side by side trawls
  //industry
   init_bounded_number srvind_q(0.10,1.000,survsel1_phase+2)
   init_bounded_number srvind_q_f(0.01,1.000,survsel1_phase+2)
   init_bounded_number srvind_sel95_f(55.0,120.0,survsel1_phase)
   init_bounded_number srvind_sel50_f(-50.0,55.0,survsel1_phase)

 //2010 bsfrf study area double logistic
  //industry
   init_bounded_number srv10ind_q(0.1,1.000,survsel1_phase+2)
   init_bounded_number srv10ind_q_f(0.01,1.00,survsel1_phase+2)
   init_bounded_number srv10ind_sel95_f(0.05,50.0,-survsel1_phase)
   init_bounded_number srv10ind_sel50_f(0.0,50.0,-survsel1_phase)
   
   //smooth selectivity for the industry survey 
   init_bounded_vector selsmo10ind(1,nlenm,-4.0,-0.0010,survsel1_phase)
   init_bounded_vector selsmo09ind(1,nlenm,-4.0,-0.0010,survsel1_phase)

//model 5 estimate M
    init_bounded_number Mmult_imat(0.2000,2.000001,natM_phase)
    init_bounded_number Mmult(0.20000,2.0000001,natM_phase)
    init_bounded_number Mmultf(1.0000,1.000001,-5)
    
//pot fishery cpue q
    // init_bounded_number cpueq(0.000877,0.000877,-5)
    init_bounded_number cpueq(0.0000877,0.00877,5)
	
// proportion new shell in recruits
  init_bounded_number proprecn(1.0,1.0,-2)

////end of estimated parameters///////////////
  //Weight at length vectors
  matrix wtf(1,2,1,nlenm)   																	 //weight at length juvenile and mature females (from kodiak program)
  vector wtm(1,nlenm)  
  
 //Selectivity related variables
  vector fish_sel50_mn(styr,endyr+Nproj)
  vector fish_sel50_mo(styr,endyr+Nproj)
  vector log_sel_fish(1,nlenm)
  matrix log_sel_fish_disc(1,2,1,nlenm)
  matrix log_sel_srv1(1,2,1,nlenm)
  3darray sel(1,2,styr,endyr+Nproj,1,nlenm)
  matrix sel_discf(1995,endyr+Nproj,1,nlenm)
  vector sel_discf_e(1,nlenm)
  3darray sel_fit(1,2,styr,endyr+Nproj,1,nlenm)
  3darray sel_ret(1,2,styr,endyr+Nproj,1,nlenm)
  matrix sel_trawl(1,2,1,nlenm)
  matrix sel_srv1(1,2,1,nlenm)
  matrix sel_srv2(1,2,1,nlenm)
  matrix sel_srv3(1,2,1,nlenm)
  matrix sel_srvind(1,2,1,nlenm)
  matrix sel_srvnmfs(1,2,1,nlenm)  
  matrix sel_srv10ind(1,2,1,nlenm)
  matrix sel_srv10nmfs(1,2,1,nlenm)  
  number avgsel_fish
  vector avgsel_fish_disc(1,2)
  vector avgsel_srv1(1,2)
  
  //summary population quantities (year and sex)
  vector popn_lmale(styr,endyr+Nproj)
  matrix popn_disc(1,2,styr,endyr+Nproj)
  matrix totn_srv1(1,2,styr,endyr+Nproj)
  matrix totn_srv2(1,2,1,2)							//2009 survey numbers
  matrix totn_srv10(1,2,1,2)
  matrix totn_trawl(1,2,styr,endyr+Nproj)
  vector explbiom(styr,endyr+Nproj)
  vector pred_bio(styr,endyr+Nproj)
  vector fspbio(styr,endyr+Nproj)
  vector mspbio(styr,endyr+Nproj)
  
  //population quantities by length, year, and sex
  3darray pred_srv1(1,2,styr,endyr+Nproj,1,nlenm)
  3darray pred_p_fish(1,2,styr,endyr+Nproj,1,nlenm)
  3darray pred_p_fish_fit(1,2,styr,endyr+Nproj,1,nlenm)
  3darray pred_p_fish_discm(1,2,styr,endyr+Nproj,1,nlenm)
  matrix  pred_p_fish_discf(styr,endyr+Nproj,1,nlenm)
  3darray pred_p_trawl(1,2,styr,endyr+Nproj,1,nlenm)
  4darray pred_p_srv1_len(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  4darray pred_p_srv1_len_new(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  4darray pred_p_srv1_len_old(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  3darray pred_p_srv2_len(1,2,1,2,1,nlenm)
  matrix  pred_p_srv2_len_ind(1,2,1,nlenm)
  matrix  pred_p_srv2_len_nmfs(1,2,1,nlenm)
  matrix  pred_p_srv10_len_ind(1,2,1,nlenm)
  matrix  pred_p_srv10_len_nmfs(1,2,1,nlenm)
  
  //summary catch statistics (year and sex)
  vector pred_catch(styr,endyr+Nproj)
  vector pred_catch_ret(styr,endyr+Nproj)
  matrix pred_catch_disc(1,2,styr,endyr+Nproj)
  vector pred_catch_trawl(styr,endyr+Nproj)
  matrix catch_male_ret(styr,endyr+Nproj,1,nlenm)
  matrix catch_lmale(styr,endyr+Nproj,1,nlenm)
  
  //catch quantities by length, year, and sex
  3darray catch_discp(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_iold(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_inew(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_mold(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_mnew(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_old(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_new(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_i(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natlength_mat(1,2,styr,endyr+Nproj,1,nlenm)
  
  //growth and molting related variables
  3darray len_len(1,2,1,nlenm,1,nlenm)          	  	 //length to length growth array
  matrix moltp(1,2,1,nlenm)                      				 //molting probabilities for female, male by length bin 
  matrix moltp_mat(1,2,1,nlenm)                 		    //molting probs for mature female, male by length bin
 
 //fishing mortality related variables
 4darray Ftot(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  3darray F(1,2,styr,endyr+Nproj,1,nlenm)
  3darray F_ret(1,2,styr,endyr+Nproj,1,nlenm)
  3darray Fmat(1,2,styr,endyr+Nproj,1,nlenm)
  3darray Fmat_ret(1,2,styr,endyr+Nproj,1,nlenm)
  3darray Fimat(1,2,styr,endyr+Nproj,1,nlenm)
  3darray Fimat_ret(1,2,styr,endyr+Nproj,1,nlenm)
  3darray Fdiscm(1,2,styr,endyr+Nproj,1,nlenm)
  matrix  Fdiscf(styr,endyr+Nproj,1,nlenm)
  3darray Fdisct(1,2,styr,endyr+Nproj,1,nlenm)  
  4darray S(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  4darray Simat(1,2,1,2,styr,endyr+Nproj,1,nlenm)
  4darray Smat(1,2,1,2,styr,endyr+Nproj,1,nlenm)

  //how is fmort, fmortret, different from F, F_ret??
  vector fmort(styr,endyr+Nproj)
  vector fmortret(styr,endyr+Nproj)
  vector fmortdf(styr,endyr+Nproj)
  vector fmortt(styr,endyr+Nproj)
  vector fmort_disc(styr,endyr+Nproj)
    
  number rbar
  vector surv(1,2)
  vector offset(1,5)
  vector matest(1,nlenm)  
  vector matestf(1,nlenm)
  
  //likelihood components
  number rec_like
  number catch_like1
  number catch_like2
  number catch_liket
  number catch_likef
  number growth_like
  4darray len_likeyr(1,2,1,2,1,2,1,nobs_srv1_length)
  vector len_like(1,7)
  vector sel_like(1,4)
  number sel_like_50m
  number like_natm
  number fpen    
  number surv_like
  number surv2_like
  number surv3_like
  number surv10_like
  number sexr_like
  number like_initsmo
  number like_q
  
  number like_linff
  number like_linfm
  number like_growthkf
  number like_growthkm
  number like_af
  number like_am
  number like_bf
  number like_bf1
  number like_bm
  number like_a1
  number like_b1
  number like_meetpt
  
  number like_srvsel
  number selsmo_like
  number discf_like
  number largemale_like
  number like_mat
  number like_initnum
  number len_like10_ind
  number len_like10_nmfs
  number surv10nmfs_like
    
  
  //quantities to report in SD file
  sdreport_vector fspbios(styr,endyr+Nproj)
  sdreport_vector mspbios(styr,endyr+Nproj)
  sdreport_vector legal_malesd(styr,endyr+Nproj)
  sdreport_vector recf_sd(styr,endyr+Nproj-1)
  sdreport_vector recm_sd(styr,endyr+Nproj-1)
  sdreport_number depletion
  
  objective_function_value f
  
  number tmp
  vector tmpn(1,nlenm)
  vector pred_sexr(styr,endyr+Nproj)
  vector preds_sexr(styr,endyr+Nproj)
  vector predpop_sexr(styr,endyr+Nproj)
  
  number maxsel_fish
  number maxsel_srv1
  number maxseld
  matrix mean_length(1,2,1,nlenm)
  matrix sd_mean_length(1,2,1,nlenm)
  matrix sum_len(1,2,1,nlenm)
  number tmp1
  number tmp2
  number tmp3
  vector sumfish(1,nobs_fish)
  vector sumfishret(1,nobs_fish)
  vector sumfishtot(1,nobs_fish_discm)
  vector sumfishdiscm(1,nobs_fish_discm)
  vector sumfishdiscf(1,nobs_fish_discf)
  vector sumtrawl(1,nobs_trawl)
  vector sumsrv(1,nobs_srv1_length)
  
  //observed proportion at length for surveys
  4darray obs_p_srv2_len(1,2,1,2,1,2,1,nlenm) 
  4darray obs_p_srv10_len(1,2,1,2,1,2,1,nlenm) 
  5darray obs_p_srv1_len(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm) 
  5darray obs_p_srv1_len1(1,2,1,2,1,2,1,nobs_srv1_length,1,nlenm) 
  3darray obs_p_srv1_lenc(1,2,1,nobs_srv1_length,1,nlenm) 											 //observed proportion survey length female, males,combined new/old
  4darray obs_p_fish(1,2,1,2,1,nobs_fish,1,nlenm) 															 //observed proportion fishery length-combined ret/disc
  vector obs_sexr(1,nobs_fish)  																					        //observed fraction female fishery length data
  vector obs_sexr_srv1_l(1,nobs_srv1_length)  																 //observed fraction female survey length data
  
  //observed proportion at length fishery
  3darray obs_p_fish_ret(1,2,1,nobs_fish,1,nlenm)
  3darray obs_p_fish_tot(1,2,1,nobs_fish,1,nlenm)
  3darray obs_p_fish_discm(1,2,1,nobs_fish_discm,1,nlenm)
  matrix obs_p_fish_discf(1,nobs_fish_discf,1,nlenm)
  3darray obs_p_trawl(1,2,1,nobs_trawl,1,nlenm)
  vector pf(1,nlenm)
  vector pb(1,nlenm)

   //management related quanities
  vector catch_tot(styr,endyr+Nproj) 
  vector legal_males(styr,endyr+Nproj)
  vector legal_srv_males(styr,endyr+Nproj)
  vector popn(styr,endyr+Nproj)
  vector popn_fit(styr,endyr+Nproj)
  3darray obs_srv1_num(1,2,styr,endyr+Nproj,1,nlenm)
  matrix tmpo(1,2,styr,endyr+Nproj)
  matrix tmpp(1,2,styr,endyr+Nproj)
  matrix obs_srv1_bioms(1,2,styr,endyr+Nproj)
  vector obs_srv1_biom(styr,endyr+Nproj)
  matrix obs_srv1_spbiom(1,2,styr,endyr+Nproj)
  matrix obs_srv2_spbiom(1,2,1,2)
  matrix obs_srv10_spbiom(1,2,1,2)
  3darray obs_srv1_spnum(1,2,1,2,styr,endyr+Nproj)
  vector pred_srv1_biom(styr,endyr+Nproj)
  matrix pred_srv1_bioms(1,2,styr,endyr+Nproj)
  vector obs_srv1t(styr,endyr+Nproj)
  matrix cv_srv1(1,2,styr,endyr+Nproj)
  matrix cv_srv1_nowt(1,2,styr,endyr+Nproj)
  vector obs_catchdm_biom(styr,endyr+Nproj)
  vector obs_catchdf_biom(styr,endyr+Nproj)
  vector obs_catcht_biom(styr,endyr+Nproj)
  vector obs_catchtot_biom(styr,endyr+Nproj)
  vector avgp(1,nlenm)
  vector avgpf(1,nlenm)
  vector avgpm(1,nlenm)
  matrix biom_tmp(1,2,styr,endyr+Nproj)
  matrix discardc(1,2,1,nlenm)
  vector rec_len(1,nlenm)
  vector avg_rec(1,2)
  vector obs_lmales(1,nobs_srv1_length)
  matrix rec_dev(1,2,styr,endyr)
  matrix avgsel(1,2,1,nlenm)
  matrix avgsel_fit(1,2,1,nlenm)
  number sel_avg_like
  4darray effn_srv1(1,2,1,2,1,2,styr,endyr+Nproj)
  matrix effn_fish_ret(1,2,styr,endyr+Nproj)
  matrix effn_fish_tot(1,2,styr,endyr+Nproj)

  //population broken up by new and old shell reported at different time
  vector legal_males_bio(styr,endyr+Nproj)
  vector legal_srv_males_n(styr,endyr+Nproj)
  vector legal_srv_males_o(styr,endyr+Nproj)
  vector legal_srv_males_bio(styr,endyr+Nproj)
  vector obs_lmales_bio(1,nobs_srv1_length)
  3darray natl_new_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natl_old_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natl_inew_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natl_iold_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natl_mnew_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  3darray natl_mold_fishtime(1,2,styr,endyr+Nproj,1,nlenm)
  vector popn_lmale_new(styr,endyr+Nproj)
  vector popn_lmale_old(styr,endyr+Nproj)
  
  //catch in different fleets of different shell conditions
  matrix catch_lmale_new(styr,endyr+Nproj,1,nlenm)
  matrix catch_lmale_old(styr,endyr+Nproj,1,nlenm)
  matrix catch_male_ret_new(styr,endyr+Nproj,1,nlenm)
  matrix catch_male_ret_old(styr,endyr+Nproj,1,nlenm)
  vector popn_fit_new(styr,endyr+Nproj)
  vector popn_fit_old(styr,endyr+Nproj)
  vector popn_lmale_bio(styr,endyr+Nproj)
  number like_mmat
  number sumrecf
  number sumrecm
  vector cv_rec(1,2)
  vector fspbio_srv1(styr,endyr+Nproj)
  vector mspbio_srv1(styr,endyr+Nproj)
  number fspbio_srv2_ind
  number mspbio_srv2_ind
  number fspbio_srv2_nmfs
  number mspbio_srv2_nmfs
  number fspbio_srv10_ind
  number mspbio_srv10_ind
  number fspbio_srv10_nmfs
  number mspbio_srv10_nmfs
  matrix fspbio_srv1_num(1,2,styr,endyr+Nproj)
  matrix mspbio_srv1_num(1,2,styr,endyr+Nproj)
  3darray offset_srv(1,2,1,2,1,2)
  3darray len_like_srv(1,2,1,2,1,2)
  number growinc_90
  number growinc_67
  number hrate
  number ghl
  number ghl_number
  vector emspbio_matetime(styr,endyr+Nproj)
  vector efspbio_matetime(styr,endyr+Nproj)
  vector mspbio_matetime(styr,endyr+Nproj)
  vector mspbio_old_matetime(styr,endyr+Nproj)
  vector fspbio_new_matetime(styr,endyr+Nproj)
  vector efspbio_new_matetime(styr,endyr+Nproj)
  vector fspnum_new_matetime(styr,endyr+Nproj)
  vector fspbio_matetime(styr,endyr+Nproj)
  vector mspbio_fishtime(styr,endyr+Nproj)
  vector fspbio_fishtime(styr,endyr+Nproj)
  vector emspnum_old_matetime(styr,endyr+Nproj)
  vector mspnum_matetime(styr,endyr+Nproj)
  vector efspnum_matetime(styr,endyr+Nproj)
  vector pred_catch_gt101(styr,endyr+Nproj)
  vector pred_catch_no_gt101(styr,endyr+Nproj)
  vector num_males_gt101(styr,endyr+Nproj)
  vector bio_males_gt101(styr,endyr+Nproj)
  vector obs_tmp(styr,endyr+Nproj)
  vector catch_midpt(styr,endyr+Nproj)
  vector alpha(1,2)
  number alpha_rec
  number avg_beta
  number devia
  number recsum
  vector tmps(1,nlenm)
  number tmpi
  matrix catch_disc(1,2,styr,endyr+Nproj)
  number old_mult
  matrix maturity_est(1,2,1,nlenm)
  vector cpue_pred(styr,endyr+Nproj)
  number cpue_like
  vector ftarget(styr,endyr+Nproj)
  vector pred_catch_target(styr,endyr+Nproj)
  number ghl_like
  number wt_lmlike2
  number multi
  vector tmpmnew(1,nlenm)
  vector tmpinew(1,nlenm)
  vector tmpmold(1,nlenm)
  vector offset_srv2(1,2)
  vector offset_srv10(1,2)
   vector selsmo2_n(1,nlenm)
   vector selsmof2_n(1,nlenm)
   vector Fout(1,30)
   number call_no
    vector M_matn(1,2)
    vector M_mato(1,2)
    vector M(1,2)
    vector tmpp1(1,nlenm)
    vector tmpp2(1,nlenm)
    vector tmpp3(1,nlenm)
    vector tmpp4(1,nlenm)
    number len_like_immat
    number len_like_mat
    number len_like_ex
    number sumobs
    number sumobs2
    number sumobs10
    number sumobs102

   likeprof_number lp_q;
  
// For projection functions
   number FutRec
   number FutMort
   number Bzero
   number Target
   number F35
   number SBPRF35
   number Ratio
   number Btest
   number OFL
   number FOFL
   number Bmsy
 
 //Francis weighting stuff 
  4darray Lbar(1,2,1,2,1,2,1,nobs_srv1_length) 
  3darray Lbar_hat_old(1,2,1,2,1,nobs_srv1_length) 
  3darray Lbar_hat_new(1,2,1,2,1,nobs_srv1_length)
  3darray SE_Lbar_hat_old(1,2,1,2,1,nobs_srv1_length) 
  3darray SE_Lbar_hat_new(1,2,1,2,1,nobs_srv1_length)
  3darray Francis_var_temp_new(1,2,1,2,1,nobs_srv1_length)
  3darray Francis_var_temp_old(1,2,1,2,1,nobs_srv1_length)
  number Francis_weight_m
  number Francis_weight_f
  number countFem
  number countMal
  number totalFem
  number totalMal
  number FemMeanVarTerm
  number MaleMeanVarTerm
  
  //==============================================================================
PRELIMINARY_CALCS_SECTION
   int mat;
   int m;

  // Calculate weight at length
  for(i=1;i<=nlenm;i++)
  {
   wtm(i) = alpha_wt_m*pow(length_bins(i),beta_wt_m);
   wtf(1,i) = alpha_wt_imm_f*pow(length_bins(i),beta_wt_imm_f);
   wtf(2,i) = alpha_wt_mat_f*pow(length_bins(i),beta_wt_mat_f);
  }
    
  //change weight to tons
   wtm=wtm*.000001;
   wtf=wtf*.000001; 
   
   for(m=styr;m<=endyr;m++)
    catch_midpt(m) = catch_midptIn(m);

   for(m=endyr+1;m<=endyr+Nproj;m++)
    catch_midpt(m) = catch_midptIn(endyr);

   mnatlen_styr(1) = log(10*(obs_p_srv1_lend(1,1,2,1)+obs_p_srv1_lend(2,1,2,1)+1e-02));
   mnatlen_styr(2) = log(10*(obs_p_srv1_lend(1,2,2,1)+obs_p_srv1_lend(2,2,2,1)+1e-02));
   
   for(j=1;j<=12;j++){
   fnatlen_styr(1,j) = log(10*(obs_p_srv1_lend(1,1,1,1,j)+obs_p_srv1_lend(2,1,1,1,j)+1e-02));
   fnatlen_styr(2,j) = log(10*(obs_p_srv1_lend(1,2,1,1,j)+obs_p_srv1_lend(2,2,1,1,j)+1e-02));
   }

//use logistic maturity curve for new shell males instead of fractions by year if switch>0
//this would be for initial population not probability of moving to mature
   if(maturity_switch>0){
     maturity_average(2) = maturity_logistic;
      for(i=styr;i<=endyr;i++){
         maturity(2,i)= maturity_logistic;
        }
     }

  sumfishdiscm.initialize();
  sumfishret.initialize();
  sumfishtot.initialize();
  sumfishdiscf.initialize();
  sumtrawl.initialize();

  //get total catch - sum up catches mult by assumed mortality
  for(i=styr; i<=endyr; i++)
  {
    catch_trawl(i)=catch_trawl(i)*m_trawl;
    catch_disc(1,i)=catch_odisc(1,i)*m_disc;
    catch_disc(2,i)=catch_odisc(2,i)*m_disc;      
    catch_tot(i)=catch_numbers(i)+catch_disc(2,i);

  }
  //sum each  fishery data
  for(i=1; i<=nobs_fish_discm; i++)
  { 
     sumfishdiscm(i)+=sum(obs_p_fish_discmd(1,i));
     sumfishdiscm(i)+=sum(obs_p_fish_discmd(2,i));
  } 
  
  for(i=1; i<=nobs_fish; i++){
    for(j=1; j<=2; j++){
      sumfishret(i)+=sum(obs_p_fish_retd(j,i));
    }
  } 
  
  for(i=1; i<=nobs_fish_discf; i++){
    sumfishdiscf(i)+=sum(obs_p_fish_discfd(i));
  }
  
  for(i=1; i<=nobs_trawl; i++){
    for(k=1;k<=2;k++){
      sumtrawl(i)+=sum(obs_p_trawld(k,i));
    }
  }

//length obs sex ratio in survey
  for(i=1; i<=nobs_srv1_length;i++)
  {
    obs_p_srv1_lenc(1,i)=obs_p_srv1_lend(1,1,1,i)+obs_p_srv1_lend(1,2,1,i)+obs_p_srv1_lend(2,1,1,i)+obs_p_srv1_lend(2,2,1,i);
    obs_p_srv1_lenc(2,i)=obs_p_srv1_lend(1,1,2,i)+obs_p_srv1_lend(1,2,2,i)+obs_p_srv1_lend(2,1,2,i)+obs_p_srv1_lend(2,2,2,i);
    obs_sexr_srv1_l(i)=sum(obs_p_srv1_lenc(1,i))/(sum(obs_p_srv1_lenc(1,i))+sum(obs_p_srv1_lenc(2,i)));
  }
  
  offset.initialize();
  offset_srv2.initialize(); 

   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
      if(k<2){   //new shell
        obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)*fraction_new_error)/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
         }
        else{ 
          obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)+obs_p_fish_discmd(1,i,j)*(1.-fraction_new_error))/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
           }
    }

  for(i=1;i<=nobs_fish;i++)
   {
    obs_p_fish_ret(1,i) = obs_p_fish_retd(1,i)*fraction_new_error;
    obs_p_fish_ret(2,i) = obs_p_fish_retd(2,i)+(1.-fraction_new_error)*obs_p_fish_retd(1,i);
   }

    //make observations proportions by year      
         //fishery offset
   for (i=1; i <= nobs_fish; i++)
    {
        
      for (j=1; j<=nlenm; j++)
       {
        for(k=1; k<=2; k++)
           obs_p_fish_ret(k,i,j)=(obs_p_fish_ret(k,i,j)/sumfishret(i));

               offset(1)-=nsamples_fish(1,i)*(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j))*log(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j)+p_const);
        }
	}

   for(k=1; k<=2; k++)
    {
       for (i=1; i <= nobs_fish_discm; i++)
       {
         //make observations proportions by year      
         //fishery offset
           for (j=1; j<=nlenm; j++)
           {
               obs_p_fish_tot(k,i,j)=((obs_p_fish_ret(k,i+yrs_fish_discm(1)-1978,j)*catch_numbers(yrs_fish_discm(i)))/catch_tot(yrs_fish_discm(i)))+obs_p_fish_discm(k,i,j);
             //old and new shell together
             if (k<2)
                offset(2)-=nsamples_fish(k,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j) )*log(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j)+p_const);
           }

       }
    }
  
  //make observations proportions by year      
        //fishery offset
  for (i=1; i <= nobs_fish_discf; i++)
   for (j=1; j<=nlenm; j++)
    {
        obs_p_fish_discf(i,j)=((obs_p_fish_discfd(i,j))/sumfishdiscf(i));
        offset(3)-=nsamples_fish_discf(i)*obs_p_fish_discf(i,j)*log(obs_p_fish_discf(i,j)+p_const);
	}        

//trawl length freq
  for(k=1; k<=2; k++)
   for (i=1; i <= nobs_trawl; i++)
    {
      //make observations proportions by year      
      //fishery offset
      for (j=1; j<=nlenm; j++)
      {
        obs_p_trawl(k,i,j)=((obs_p_trawld(k,i,j))/sumtrawl(i));
        offset(5)-=nsamples_trawl(k,i)*obs_p_trawl(k,i,j)*log(obs_p_trawl(k,i,j)+p_const);
      }
      obs_catcht_biom(yrs_trawl(i))=(obs_p_trawl(1,i)*catch_trawl(yrs_trawl(i)))*wtf(2)+(obs_p_trawl(2,i)*catch_trawl(yrs_trawl(i)))*wtm;
    }

  sumsrv.initialize();
//survey length offset
 for(ll=1; ll<=nobs_srv1_length; ll++)
  for(mat=1;mat<=2;mat++)
   for(k=1; k<=2; k++)
    for(j=1; j<=2; j++)
    sumsrv(ll)+=sum(obs_p_srv1_lend(mat,k,j,ll));


  for(mat=1; mat<=2; mat++) //maturity
   for(l=1; l<=2; l++) //shell condition
    for(k=1; k<=2;k++) //sex
    for (i=1; i <= nobs_srv1_length; i++)
    for (j=1; j<=nlenm; j++)
    {
 //only do new/old shell correction for mature crab
       if(mat<2){
                obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j))/sumsrv(i);
                }
        else{
        if(l<2){
                obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j)*fraction_new_error)/sumsrv(i);
         }
        else{
               obs_p_srv1_len1(mat,l,k,i,j)=(obs_p_srv1_lend(mat,l,k,i,j)+obs_p_srv1_lend(mat,1,k,i,j)*(1.-fraction_new_error))/sumsrv(i);
            }
         }
    }




//
//use logistic maturity curve for new and old shell male survey data if switch>0 instead of yearly samples
// old shell already uses ok maturity curve
   if(maturity_switch>0)
   {
    for(i=1; i <= nobs_srv1_length; i++)
    {
     //oldshell
        tmps = (obs_p_srv1_len1(1,2,2,i)+obs_p_srv1_len1(2,2,2,i));
       obs_p_srv1_len1(2,2,2,i) = elem_prod(maturity_old_average(2),tmps);
       obs_p_srv1_len1(1,2,2,i) = elem_prod(1.0-maturity_old_average(2),tmps);
    }
   }

   avgpf=0;
   avgpm=0;
   for(i=1; i <= nobs_trawl; i++)
   {
    avgpf+=obs_p_trawl(1,i);
    avgpm+=obs_p_trawl(2,i);
   }
   avgpf=avgpf/nobs_trawl;
   avgpm=avgpm/nobs_trawl;

   for(i=styr;i<=(yrs_trawl(1)-1);i++)
    obs_catcht_biom(i)=avgpf*catch_trawl(i)*wtf(2)+avgpm*catch_trawl(i)*wtm;

  obs_catchdm_biom.initialize();
  obs_catchdf_biom.initialize();
  avgp.initialize();

  for(i=1;i<=nobs_fish_discm;i++)
  {
   //obs_p_fish_discm are proportional to total catch (retained+discard)
    obs_catchdm_biom(yrs_fish_discm(i)) = catch_tot(yrs_fish_discm(i)) * (obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i))*wtm;
    avgp += obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i);
  }
  avgp=avgp/nobs_fish_discm;
  avgpf = mean(obs_p_fish_discf);

  // Historical approximation to observed catch biomass
  for(i=styr;i<=(yrs_fish_discm(1)-1);i++)
  {
   // cout<<obs_catchdm_biom<<endl;
    obs_catchdm_biom(i)=avgp*catch_tot(i)*wtm;
  } 
  obs_catchdm_biom(endyr)=avgp*catch_tot(endyr)*wtm;

  for(i=styr;i<=endyr;i++)
  {
    obs_catchtot_biom(i)=obs_catchdm_biom(i)+catch_ret(i);
    obs_catchdf_biom(i)=(avgpf*catch_disc(1,i))*wtf(2);
  }
  
  // Compute the moulting probabilities
  get_moltingp();
  // estimate growth function
  get_growth();
  // Set maturity
  get_maturity();

  what=1;

  

 
  
   cout<<"end prelim calcs"<<endl;

//==============================================================================
PROCEDURE_SECTION

//calculating natural mortality by maturity state, sex, and shell condition
    M(1)=M_in(1)*Mmult_imat;    			 //natural mortality immature females then males
    M(2)=M_in(2)*Mmult_imat;    			 //natural mortality immature females then males
    M_matn(2)=M_matn_in(2)*Mmult;  //natural mortality mature new shell female/male
    M_mato(2)=M_mato_in(2)*Mmult;  //natural mortality mature old shell female/male
    M_matn(1)=M_matn_in(1)*Mmultf;  //natural mortality mature new shell female/male
    M_mato(1)=M_mato_in(1)*Mmultf;  //natural mortality mature old shell female/male
	
// do this when doing likelihood profile
  lp_q=srv3_q;
  get_maturity();
  // cout<<"maturity"<<endl;
  if(what<2)
   {
    get_obs_survey();
    what=what+1;
   }

   get_moltingp();
  // cout<<"molting"<<endl;
 
   if(active(am) || active(af))
    get_growth();
   // cout<<" end of growth "<<endl;

   get_selectivity();
   // cout<<" end of sel "<<endl;

   get_mortality();
   // cout<<" end of mort "<<endl;

   get_numbers_at_len();
   // cout<<" end of numbers at len "<<endl;

   get_catch_at_len();
   // cout<<" end of catch at len "<<endl;

   evaluate_the_objective_function();

   
   if(mceval_phase())
   {
	   Find_F35();
	   Find_OFL();
	   
	   post<<f<<" "<<Bmsy<<" "<< F35 << " " << FOFL << " " << OFL << endl;
	   post<<fspbio_srv1<<endl;
	   post<<mspbio_srv1<<endl;
	   post<<mspbio_matetime<<endl;

	   post<<fspbio_srv2_ind<<endl;
	   post<<mspbio_srv2_ind<<endl;
	   post<<fspbio_srv2_nmfs<<endl;
	   post<<mspbio_srv2_nmfs<<endl;

	   post<<fspbio_srv10_ind<<endl;
	   post<<mspbio_srv10_ind<<endl;
	   post<<fspbio_srv10_nmfs<<endl;
	   post<<mspbio_srv10_nmfs<<endl;
	   
	   post<<pred_catch<<endl;
	   post<<pred_catch_disc<<endl;
	   post<<pred_catch_ret<<endl;
	   post<<pred_catch_trawl<<endl;

      //Don't need these now because they can be taken from the .psv file, but I'm scared to get rid of them yet
	   // post<<af<< " "<<am<<" "<<bf<<" "<<bm<<" "<<b1<< " " <<bf1 << " " << deltam << " " << deltaf<< " " <<st_gr<< " "<<growth_beta<<" ";
	   // post<<matest50f<< " "<<matestslpf<<" "<<matest50m<<" "<<matestslpm<<" "<<mateste<< " " <<matestfe << " " ;
	   // post<<moltp_af<< " "<<moltp_bf<<" "<<moltp_am<<" "<<moltp_bm<<" "<<moltp_ammat<< " " <<moltp_bmmat << " " ;  
 	   // post<<mean_log_rec1<< " "<<rec_dev_mean<<" "<<mean_log_rec<<" "<<rec_devf<<" "<<alpha1_rec<< " " <<beta_rec << " " ;  
	   // post<<mnatlen_styr<<" "<< fnatlen_styr << " ";
	   // post<<log_avg_fmort<< " "<<fmort_dev<<" "<<log_avg_fmortdf<<" "<<fmortdf_dev<<" "<<log_avg_fmortt<<" " <<fmortt_dev<< " " <<discard_mult << " ";
	   // post<<log_avg_sel50_mn<< " "<<log_sel50_dev_mn<<" "<<log_avg_sel50_mo<<" "<<log_sel50_dev_mo<<" "<<fish_slope_mn<<" " <<fish_slope_mo << " ";
 	   // post<<fish_fit_slope_mn<< " "<<fish_fit_sel50_mn<<" "<<fish_fit_slope_mo<<" "<<fish_fit_sel50_mo << " " ;  	
 	   // post<<fish_slope_mo2<< " "<<fish_sel50_mo2<<" "<<fish_slope_mn2<<" "<<fish_sel50_mn2 << " " ;  	
 	   // post<<fish_disc_slope_f<< " "<<fish_disc_sel50_f<<" "<<log_dev_50f<<" "<<fish_disc_slope_tf << " " <<fish_disc_sel50_tf<<" " ;
	   // post<<srv1_q<< " "<<srv1_sel95<<" "<<srv1_sel50<<" "<<srv2_q<<" "<<srv2_sel95<<" " <<srv2_sel50 << " "<<srv3_q<<" "<<srv3_sel95<<" " <<srv3_sel50 << " ";   
	   // post<<srv3_sel95_f<<" "<< srv3_sel50_f << " "<<" ";   
 	   // post<<srvind_q<< " "<<srvind_q_f<<" "<<srvind_sel95<<" "<<srvind_sel50<<" "<<srvind_sel95_f<< " " <<srvind_sel50_f << " " ;  
 	   // post<<srvnmfs_sel95<< " "<<srvnmfs_sel50<<" "<<srvnmfs_sel95_f<<" "<<srvnmfs_sel50_f << " " ;  	
 	   // post<<srv10ind_q<< " "<<srv10ind_q_f<<" "<<srv10ind_sel95<<" "<<srv10ind_sel50<<" "<<srv10ind_sel95_f<< " " <<srv10ind_sel50_f << " " ;  
 	   // post<<srv10nmfs_sel95<< " "<<srv10nmfs_sel50<<" "<<srv10nmfs_sel95_f<<" "<<srv10nmfs_sel50_f << " " ;  	
 	   // post<<selsmo10ind<< " "<<selsmo09ind<<" "<<Mmult_imat<<" "<<Mmult<<" "<<Mmultf << " " ;  
	   // post<<cpueq<<" "<< proprecn << " "<<endl;
   }

//==============================================================================
FUNCTION get_maturity


  if(active(matestfe))
   {
     for(j=1;j<=nlenm;j++)
      {
	      if(j<(matest_n+1)){
	      matestf(j)=matestfe(j);
          }
          else{
	          matestf(j)=0.0;
            }
           maturity_est(1,j) = mfexp(matestf(j));
       }
      
   }
   else
    {    
	  maturity_est(1) = maturity_average(1);
    }

  if(active(mateste))
   {
     for(j=1;j<=nlenm;j++)
      {
	      if(j<(matestm_n+1)){
	      matest(j)=mateste(j);
          }
          else{
	          matest(j)=0.0;
            }
          maturity_est(2,j) = mfexp(matest(j));
      }
    }
    else{    
       maturity_est(2) = maturity_average(2);
    }

//==============================================================================
FUNCTION get_catch_tot

  for(i=styr;i<=endyr;i++)
  {
    catch_disc(1,i)=catch_odisc(1,i)*discard_mult*m_disc;
    catch_disc(2,i)=catch_odisc(2,i)*discard_mult*m_disc;      
    catch_tot(i)=catch_numbers(i)+catch_disc(2,i);
  }
 
//make observations proportions by year     
   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
      if(k<2)
	   {
        obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j))/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
        }
        else
		{
          obs_p_fish_discm(k,i,j)=((obs_p_fish_discmd(k,i,j)*old_mult)/sumfishdiscm(i))*(catch_disc(2,yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)));
        }
    }
   
   offset(2)=0.0;

  //make observations proportions by year; fishery offset
   for(k=1; k<=2; k++)
    for (i=1; i <= nobs_fish_discm; i++)
    for (j=1; j<=nlenm; j++)
    {
        obs_p_fish_tot(k,i,j)=(obs_p_fish_retd(k,i+yrs_fish_discm(1)-1978,j)/sumfishret(i+yrs_fish_discm(1)-1978))*(catch_numbers(yrs_fish_discm(i))/catch_tot(yrs_fish_discm(i)))+obs_p_fish_discm(k,i,j);
   //old and new shell together
       if (k<2){
         offset(2)-=nsamples_fish(k,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j) )*log(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j)+p_const);
        }
    }

  obs_catchdm_biom.initialize();
  obs_catchdf_biom.initialize();
  avgp.initialize();

  for(i=1;i<=nobs_fish_discm;i++)
  {
//obs_p_fish_discm are proportional to total catch (retained+discard)
    obs_catchdm_biom(yrs_fish_discm(i)) = catch_tot(yrs_fish_discm(i)) * (obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i))*wtm;
    avgp += obs_p_fish_discm(1,i)+obs_p_fish_discm(2,i);
  }
  avgp=avgp/nobs_fish_discm;
  avgpf = (obs_p_fish_discf(1)+obs_p_fish_discf(2)+obs_p_fish_discf(3)+obs_p_fish_discf(4)+obs_p_fish_discf(5))/5.;
 
  for(i=styr;i<=(yrs_fish_discm(1)-1);i++)
  {
    obs_catchdm_biom(i)=avgp*catch_tot(i)*wtm;
  } 
  obs_catchdm_biom(endyr)=avgp*catch_tot(endyr)*wtm;

  // cout <<obs_catchdm_biom<<endl;exit(1);

  for(i=styr;i<=endyr;i++)
  {
    obs_catchtot_biom(i)=obs_catchdm_biom(i)+catch_ret(i);
    obs_catchdf_biom(i)=(avgpf*catch_disc(1,i))*wtf(2);
  }
//==============================================================================
FUNCTION get_obs_survey
   //don't do the move with fraction_new_error -adjustment made in prelimn calcs for obs_p_srv1_len1
  //move new shells to old shells for mature only due to error in shell condition age
         for(k=1;k<=2;k++){
          obs_p_srv1_len(2,1,k)=obs_p_srv1_len1(2,1,k);
          obs_p_srv1_len(2,2,k)=obs_p_srv1_len1(2,2,k);
         }
          obs_p_srv1_len(1)=obs_p_srv1_len1(1);

 int mat;
  offset(4)=0.0;
  offset_srv2.initialize();

//for maturity and shell condition together in survey length comp fits
  for(k=1; k<=2;k++) //sex
   for (i=1; i <= nobs_srv1_length; i++)
    for (j=1; j<=nlenm; j++)
    {
          offset(4)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(1,1,k,i,j)+
                       obs_p_srv1_len(1,2,k,i,j)+
                       obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j))*
                            log(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j)+
                                   obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j)+p_const);
          
    }
//extra 2009 survey
//make proportions - first index is survey BSFRF 1, NMFS 2
    sumobs=sum(obs_p_srv2_lend(1,1,2)(5,nlenm)+obs_p_srv2_lend(1,2,2)(5,nlenm)+obs_p_srv2_lend(1,1,1)(5,nlenm)+obs_p_srv2_lend(1,2,1)(5,nlenm));
    sumobs2=sum(obs_p_srv2_lend(2,1,2)(5,nlenm)+obs_p_srv2_lend(2,2,2)(5,nlenm)+obs_p_srv2_lend(2,1,1)(5,nlenm)+obs_p_srv2_lend(2,2,1)(5,nlenm));

//
  for(il2=1;il2<=2;il2++) 
    for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
      if(j>4){
        obs_p_srv2_len(1,il2,mat,j)=obs_p_srv2_lend(1,il2,mat,j)/sumobs;
        obs_p_srv2_len(2,il2,mat,j)=obs_p_srv2_lend(2,il2,mat,j)/sumobs2;
        }
        else{                    
        obs_p_srv2_len(1,il2,mat,j)=0.0;
        obs_p_srv2_len(2,il2,mat,j)=0.0;
        }
   }

  //offset
 for(k=1; k<=2; k++)
  for(il2=1; il2<=2 ;il2++)
   for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
      offset_srv2(k) -= nsamples_srv2_length(k,il2,mat)* obs_p_srv2_len(k,il2,mat,j)*log(obs_p_srv2_len(k,il2,mat,j)+p_const);
    }

//extra 2010 survey
//make proportions - first index is survey BSFRF 1, NMFS 2
    sumobs10=sum(obs_p_srv10_lend(1,1,2)(1,nlenm)+obs_p_srv10_lend(1,2,2)(1,nlenm)+obs_p_srv10_lend(1,1,1)(1,nlenm)+obs_p_srv10_lend(1,2,1)(1,nlenm));
    sumobs102=sum(obs_p_srv10_lend(2,1,2)(1,nlenm)+obs_p_srv10_lend(2,2,2)(1,nlenm)+obs_p_srv10_lend(2,1,1)(1,nlenm)+obs_p_srv10_lend(2,2,1)(1,nlenm));

 for(il2=1;il2<=2;il2++) 
   for(mat=1;mat<=2;mat++)
    for(j=1;j<=nlenm;j++)
    {
       obs_p_srv10_len(1,il2,mat,j)=obs_p_srv10_lend(1,il2,mat,j)/sumobs10;
       obs_p_srv10_len(2,il2,mat,j)=obs_p_srv10_lend(2,il2,mat,j)/sumobs102;
    }
	
  //offset
   offset_srv10.initialize();
//
 for(k=1; k<=2; k++)
  for(il2=1; il2<=2 ;il2++)
   for(mat=1;mat<=2;mat++)
	for(j=1;j<=nlenm;j++)
	 offset_srv10(k) -= nsamples_srv10_length(k,il2,mat)* obs_p_srv10_len(k,il2,mat,j)*log(obs_p_srv10_len(k,il2,mat,j)+p_const);


  obs_srv1_num.initialize();
  obs_srv1_biom.initialize();
  obs_srv1_bioms.initialize();
  obs_srv1_spbiom.initialize();
  obs_srv1_spnum.initialize();
  obs_srv2_spbiom.initialize();
  obs_srv10_spbiom.initialize();

  for(i=1;i<=nobs_srv1;i++)
    obs_srv1t(yrs_srv1(i))=obs_srv1(i);

  for(mat=1;mat<=2;mat++)  //maturity status
   for(l=1;l<=2;l++)  //shell condition
    for(k=1;k<=2;k++)  //sex
    for (i=1; i <= nobs_srv1_length; i++)
    {
         obs_srv1_num(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i));
        if(k<2){
          obs_srv1_bioms(k,yrs_srv1_length(i))+=obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
         obs_srv1_biom(yrs_srv1_length(i))+= obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
         }
        else{
         obs_srv1_bioms(k,yrs_srv1_length(i))+=obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
         obs_srv1_biom(yrs_srv1_length(i))+= obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
        }        

//  sum to get mature biomass by sex
       if(mat>1)
       {
        if(k<2){
        obs_srv1_spbiom(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtf(mat);
        }
        else{
        obs_srv1_spbiom(k,yrs_srv1_length(i)) += obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i))*wtm;
        }
        obs_srv1_spnum(l,k,yrs_srv1_length(i)) += sum(obs_p_srv1_len(mat,l,k,i)*obs_srv1t(yrs_srv1_length(i)));
       }
    }

//industry area survey 2009
//first index bsfrf or nmfs, 2nd index sex female, male, 3rd immat, mature	this is 1000s of tons        
   obs_srv2_spbiom(1,1) = (obs_p_srv2_len(1,1,2)*sumobs*wtf(2))*1000.;
   obs_srv2_spbiom(2,1) = (obs_p_srv2_len(2,1,2)*sumobs2*wtf(2))*1000.;
   obs_srv2_spbiom(1,2) = (obs_p_srv2_len(1,2,2)*sumobs*wtm)*1000.;
   obs_srv2_spbiom(2,2) = (obs_p_srv2_len(2,2,2)*sumobs2*wtm)*1000.;

//industry area survey 2010
//first index bsfrf or nmfs, 2nd index sex female, male, 3rd immat, mature	this is 1000s of tons        
   obs_srv10_spbiom(1,1) = (obs_p_srv10_len(1,1,2)*sumobs10*wtf(2))*1000.;
   obs_srv10_spbiom(2,1) = (obs_p_srv10_len(2,1,2)*sumobs102*wtf(2))*1000.;
   obs_srv10_spbiom(1,2) = (obs_p_srv10_len(1,2,2)*sumobs10*wtm)*1000.;
   obs_srv10_spbiom(2,2) = (obs_p_srv10_len(2,2,2)*sumobs102*wtm)*1000.;


//==============================================================================
FUNCTION dvector bootstrap(dvector& data,int& sample_size, long int seed)

 //data contains the data from which the bootstrap sample is to be taken
 //sample_size is the size of the desired bootstrap sample
 //seed is a seed for the random number generator
 //use like this
 //  dvector boot_sample=bootstrap(data,100,1231);  takes a sample with replacement of size
 // 100 from data using seed 1231 and puts it into the vector boot_sample.

 dvector tmp(1,sample_size);  //declare vector tmp of length sample_size
 tmp.fill_randu(seed);       // gives uniform random numbers from 0 to 1 size sample_size
 tmp=(tmp*data.size())+1.0;      // makes tmp random numbers from 0 to length of data+1
 ivector iselect(tmp);       // takes the integer part of tmp and puts it into a vector iselect
 tmp=data(iselect);          // takes the bootstrap sample from data puts into tmp

 return(tmp);                // return the bootstrap sample in tmp
//==============================================================================
FUNCTION get_growth
  int ilen,il2,sex;
  double devia ;
  dvariable lensum;
  len_len.initialize();
  rec_len.initialize();

  for(ilen=1;ilen<=nlenm;ilen++)
  { 
//linear growth curve
  if(growth_switch==1)
   {
//growth slope different males and females first two bins different
     mean_length(1,ilen) = ((af+bf*length_bins(ilen))*(1-cumd_norm((length_bins(ilen)-deltaf)/st_gr))+((af+(bf-bf1)*deltaf)+bf1*length_bins(ilen))*(cumd_norm((length_bins(ilen)-deltaf)/st_gr)));
     mean_length(2,ilen) = ((am+bm*length_bins(ilen))*(1-cumd_norm((length_bins(ilen)-deltam)/st_gr))+((am+(bm-b1)*deltam)+b1*length_bins(ilen))*(cumd_norm((length_bins(ilen)-deltam)/st_gr)));
   }
  }

   for (sex=1;sex<=2;sex++)
  {
    for(ilen=1;ilen<=nlenm;ilen++)
    {
      // subract the 2.5 from the midpoint of the length bin to get the lower bound
      alpha(sex) = (mean_length(sex,ilen)-(length_bins(ilen)-2.5))/growth_beta(sex);
      lensum = 0;
      for(il2=ilen;il2<=ilen+min(8,nlenm-ilen);il2++)
      {
          devia = length_bins(il2)+2.5-length_bins(ilen);
           len_len(sex,ilen,il2) = pow(devia,(alpha(sex)-1.))*mfexp(-devia/growth_beta(sex));
          lensum += len_len(sex,ilen,il2);
          // cout << devia<<" "<<alpha<<" "<<growth_beta(sex)<<endl;

   //standardize so each row sums to 1.0
      }  
      len_len(sex,ilen) /= lensum;
    }
  }
 // Fraction recruiting
  recsum=0.0;
  alpha_rec=alpha1_rec/beta_rec;
  
 //do only first 6 bins
  for(ilen=1;ilen<=6;ilen++)
  {
    devia = length_bins(ilen)+2.5-length_bins(1);
    rec_len(ilen) = pow(devia,alpha_rec-1.)*mfexp(-devia/beta_rec);
    recsum += rec_len(ilen);
  }

  //standardize so each row sums to 1.0
  for(ilen=1;ilen<=nlenm;ilen++)
   rec_len(ilen) = rec_len(ilen)/recsum;
 
//==============================================================================
FUNCTION get_moltingp
//assuming a declining logistic function
  for(j=1;j<=nlenm;j++)
  {
 //  these next two lines are logistic molting females then males
     moltp(1,j)=1-(1./(1.+mfexp(-1.*moltp_af*(length_bins(j)-moltp_bf))));
      moltp(2,j)=1-(1./(1.+mfexp(-1.*moltp_am*(length_bins(j)-moltp_bm))));
     //set molting prob for mature females at 0.0
      moltp_mat(1,j)=0.0;
     if(phase_moltingp>0)
       {
       moltp_mat(2,j)=1-(1./(1.+mfexp(-1.*moltp_ammat*(length_bins(j)-moltp_bmmat))));
        }
        else
        {
 //mature male molting probability equal zero - terminal molt--I ithought we all agreed that there is a terminal molt?
          moltp_mat(2,j)=0.0;
        }
    }
//  cout<<moltp<<endl;
//==============================================================================
FUNCTION get_selectivity
  dvariable tempSeln, tempSelo;
  int x;
//logistic selectivity curves
 for(iy=styr;iy<=endyr+Nproj;iy++)
   {
    if(iy<1978)
    {
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn);
     if(active(log_avg_sel50_mo))
      fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo);
    }
    if(iy>=1978 & iy<endyr)
    {      
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn+log_sel50_dev_mn(iy));
     if(active(log_avg_sel50_mo))
       fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo+log_sel50_dev_mo(iy));
    }
    if(iy>=endyr)  
    {  
     tempSeln = 0;
	 tempSelo = 0;
     for(x=(endyr-sel_avg_Nyrs);x<=(endyr-1);x++)
	 {
      tempSeln += log_sel50_dev_mn(x);
	  tempSelo += log_sel50_dev_mo(x);
	 }
	 tempSeln = tempSeln/sel_avg_Nyrs;
	 tempSelo = tempSelo/sel_avg_Nyrs;
	 
     fish_sel50_mn(iy)=mfexp(log_avg_sel50_mn+tempSeln);
     if(active(log_avg_sel50_mo))
       fish_sel50_mo(iy)=mfexp(log_avg_sel50_mo+tempSelo);
    }	
	
  //logistic selectivity curve
 for (j=1;j<=nlenm;j++)
  { 
    sel(1,iy,j)=1./(1.+mfexp(-1.*fish_slope_mn*(length_bins(j)-fish_sel50_mn(iy))));
   //set new and old sel same
    sel(2,iy,j)= sel(1,iy,j);
	//for dome shaped add this part
    if(phase_fishsel>0)
     {
      tmp2=1./(1.+mfexp(fish_slope_mn2*(length_bins(j)-fish_sel50_mn2)));
      tmp3=1./(1.+mfexp(fish_slope_mo2*(length_bins(j)-fish_sel50_mo2)));
      sel(1,iy,j)=sel(1,iy,j)*tmp2;
      sel(2,iy,j)=sel(1,iy,j);
      }
      sel_ret(1,iy,j)=1./(1.+mfexp(-1.*fish_fit_slope_mn*(length_bins(j)-fish_fit_sel50_mn)));
      sel_fit(1,iy,j)=sel_ret(1,iy,j)*sel(1,iy,j);
      sel_ret(2,iy,j)=sel_ret(1,iy,j);
      sel_fit(2,iy,j)=sel_ret(2,iy,j)*sel(2,iy,j);
   }
  }  //end of year loop

 //female discards ascending logistic curve 
  for (j=1;j<=nlenm;j++)
  {
   for(iy=yrs_fish_discf(1);iy<(yrs_fish_discf(nobs_fish_discf)+.5);iy++)
    sel_discf(iy,j)=1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins(j)-mfexp(fish_disc_sel50_f+log_dev_50f(iy)))));

   for(iy=yrs_fish_discf(nobs_fish_discf);iy<endyr+Nproj;iy++)
     sel_discf(iy,j)=1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins(j)-mfexp(fish_disc_sel50_f+log_dev_50f(endyr-1)))));
	
	
// trawl fishery selectivity
    sel_trawl(1,j)=1./(1.+mfexp(-1.*fish_disc_slope_tf*(length_bins(j)-fish_disc_sel50_tf)));
//set male = female for trawl bycatch
    sel_trawl(2,j)=sel_trawl(1,j); 

 if(j<=nsellen_srv1)
   {
//somerton and otto curve for survey selectivities
    if(somertonsel<0)
    { 
       sel_srv3(1,j)=1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50_f)/(srv3_sel95_f-srv3_sel50_f)));
       sel_srv3(2,j)=sel_som(1)/(1.+mfexp(-1.*(sel_som(2)+(sel_som(3)*length_bins(j)))));           
    }
    else
    { 
    sel_srv3(2,j)=srv3_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50)/(srv3_sel95-srv3_sel50)));
//female survey 3
    sel_srv3(1,j)=srv3_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv3_sel50_f)/(srv3_sel95_f-srv3_sel50_f)));
    }

// surv sel 2009 study area                
   sel_srvind(2,j)=srvind_q*mfexp(selsmo09ind(j));
   sel_srvind(1,j)= srvind_q_f* 1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srvind_sel50_f)/(srvind_sel95_f-srvind_sel50_f)));
   sel_srvnmfs(1,j)= sel_srvind(1,j)*sel_srv3(1,j);
   sel_srvnmfs(2,j)= sel_srvind(2,j)*sel_srv3(2,j);

   // surv sel 2010 study area
   sel_srv10ind(2,j)=srv10ind_q*mfexp(selsmo10ind(j));
   sel_srv10ind(1,j)= srv10ind_q_f* 1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv10ind_sel50_f)/(srv10ind_sel95_f-srv10ind_sel50_f)));
   sel_srv10nmfs(1,j)=sel_srv10ind(1,j)*sel_srv3(1,j);
   sel_srv10nmfs(2,j)=sel_srv10ind(2,j)*sel_srv3(2,j); 

   // this sets time periods 1 and 2 survey selectivities to somerton otto as well
  if(survsel1_phase<0){
    sel_srv1(2,j)=sel_srv3(2,j);
    sel_srv2(2,j)=sel_srv3(2,j);
    sel_srv1(1,j)=sel_srv3(1,j);
    sel_srv2(1,j)=sel_srv3(1,j);
    }
    else
    { 
//logistic curve if estimating selectivity parameters
    sel_srv1(2,j)=srv1_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv1_sel50)/(srv1_sel95-srv1_sel50)));
    sel_srv2(2,j)=srv2_q*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv2_sel50)/(srv2_sel95-srv2_sel50)));
    sel_srv1(1,j)=srv1_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv1_sel50)/(srv1_sel95-srv1_sel50)));
    sel_srv2(1,j)=srv2_q_f*1./(1.+mfexp(-1.*log(19.)*(length_bins(j)-srv2_sel50)/(srv2_sel95-srv2_sel50)));
    }
    }
    else
   {
    sel_srv1(1,j)=sel_srv1(1,j-1);
    sel_srv1(2,j)=sel_srv1(2,j-1);
    }
  }
  
  for(iy=styr;iy<=endyr+Nproj;iy++)
    {
     maxsel_fish=max(sel(1,iy));
     if(maxsel_fish<max(sel(2,iy)))
	   maxsel_fish=max(sel(2,iy));
     
	 sel(1,iy)=sel(1,iy)/maxsel_fish;
     sel(2,iy)=sel(2,iy)/maxsel_fish;

     sel_fit(1,iy)=elem_prod(sel_ret(1,iy),sel(1,iy));
     sel_fit(2,iy)=elem_prod(sel_ret(2,iy),sel(2,iy));
    }
  
//==============================================================================
FUNCTION get_mortality
 int x;
 
 fmort(styr,endyr-1) = mfexp(log_avg_fmort+fmort_dev);
 fmort(endyr)=mfexp(log_avg_fmort);  // last year not fit to

 fmortdf(styr,endyr-1)=mfexp(log_avg_fmortdf+fmortdf_dev);
 fmortdf(endyr)=mfexp(log_avg_fmortdf);

 for(x=styr;x<=1991;x++)
  fmortt(x)=mfexp(log_avg_fmortt+fmortt_dev_era1(x));
 for(x=1992;x<=endyr-1;x++)
  fmortt(x)=mfexp(log_avg_fmortt+fmortt_dev_era2(x));

  fmortt(endyr)=mfexp(log_avg_fmortt);
 
  for (i=styr;i<=endyr;i++)
 {
//have discard mort for females and males fishing F for males only
   if(i>=yrs_fish_discf(1) && i<endyr)
   {
    Fdiscf(i)=sel_discf(i)*fmortdf(i);
   }
   else
   {
    sel_discf_e=(1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins-mfexp(fish_disc_sel50_f)))));
    Fdiscf(i)=sel_discf_e*fmortdf(i);
    }
   Fdisct(1,i)=sel_trawl(1)*fmortt(i);   
   Fdisct(2,i)=sel_trawl(2)*fmortt(i);   
   
   for(k=1;k<=2;k++) //over new (k=1) and old (k=2) shell...
   { 
     F(k,i) = sel(k,i)*fmort(i);       
     F_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fmat(k,i) = sel(k,i)*fmort(i);       
     Fmat_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fimat(k,i) = sel(k,i)*fmort(i);      // Fishing mort on immature males new or old shell  
     Fimat_ret(k,i)= sel_fit(k,i)*fmort(i);
     Ftot(1,k,i)=Fdiscf(i) + Fdisct(1,i);
     S(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Simat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(2,k,i)=mfexp(-1.0*(Fmat(k,i)+Fdisct(2,i)));
     Simat(2,k,i)=mfexp(-1.0*(Fimat(k,i)+Fdisct(2,i)));

   } 
 }
//==============================================================================
FUNCTION get_numbers_at_len
  int itmp;
  natlength_new.initialize();
  natlength_old.initialize();
  natlength_inew.initialize();
  natlength_iold.initialize();
  natlength_mnew.initialize();
  natlength_mold.initialize();
  natlength_i.initialize();
  natlength_mat.initialize();
  natlength.initialize();
  efspbio_matetime.initialize();
  emspbio_matetime.initialize();
  mspbio_old_matetime.initialize();
  fspbio_new_matetime.initialize();
  efspbio_new_matetime.initialize();
  fspnum_new_matetime.initialize();
  efspnum_matetime.initialize();
  emspnum_old_matetime.initialize();
  mspnum_matetime.initialize();
  bio_males_gt101.initialize();
  num_males_gt101.initialize();
   popn.initialize();
   explbiom.initialize();
   pred_bio.initialize();
   pred_srv1.initialize();
   pred_srv1_biom.initialize();
   pred_srv1_bioms.initialize();
   fspbio.initialize();
   mspbio.initialize(); 

  mean_log_rec(1)				=mean_log_rec1*mfexp(rec_dev_mean);
  mean_log_rec(2)				=mean_log_rec1;
  rec_dev(1)					=rec_devf;
  rec_dev(2)					=rec_devf;
 
//use average maturity curve for males and females in beginning years
  dvar_matrix mat_not_average(1,2,1,nlenm);
   int jk;
  for(k=1;k<=2;k++)
   for(jk=1;jk<=nlenm;jk++)
    mat_not_average(k,jk)	=1.0-maturity_average(k,jk);

//initial numbers at length
// maturity, shell cond, sex, year, length
//split initial numbers into maturity, shell condition by obs proportions for first year
  for(j=1;j<=12;j++)
   {
      natlength_mnew(1,styr,j)	 	= maturity_est(1,j)*mfexp(fnatlen_styr(1,j));
      natlength_inew(1,styr,j) 		= (1.0-maturity_est(1,j))*mfexp(fnatlen_styr(1,j));
      natlength_mold(1,styr,j) 		= mfexp(fnatlen_styr(2,j));
    }

   for(j=13;j<=nlenm;j++)
	{
      natlength_mnew(1,styr,j)	= 0.0;
      natlength_inew(1,styr,j) 	= 0.0;
      natlength_mold(1,styr,j) 	= 0.0;
    }
   natlength_iold(1,styr) 		= 0.0;
   

  //males
   natlength_inew(2,styr) 		= mfexp(mnatlen_styr(1));
   natlength_mnew(2,styr)		= mfexp(mnatlen_styr(2));
   natlength_mold(2,styr) 		= 0.0;
   natlength_iold(2,styr) 		= 0.0;


  for(k=1;k<=2;k++)  //k is sex 
   {
	 natlength_new(k,styr)		=natlength_inew(k,styr)+natlength_mnew(k,styr);
      natlength_old(k,styr)		=natlength_iold(k,styr)+natlength_mold(k,styr);
      natlength_mat(k,styr)		=natlength_mnew(k,styr)+natlength_mold(k,styr);
      natlength_i(k,styr)		=natlength_inew(k,styr)+natlength_iold(k,styr);
      natlength(k,styr)			=natlength_inew(k,styr)+natlength_iold(k,styr) + natlength_mnew(k,styr)+natlength_mold(k,styr);
    }

 for (ipass=styr;ipass<=endyr;ipass++)
    get_num_at_len_yr();

//===========================================================
FUNCTION get_num_at_len_yr	
 i = ipass;
//numbers at length from styr to endyr
 // if(i < endyr)
 // {
 for(k=1;k<=2;k++)
  {
      // Numbers advancing to new shell...
      dvar_vector tmp 			= elem_prod(moltp(k)*mfexp(-(1-catch_midpt(i))*M(k)),elem_prod(Simat(k,1,i),mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i)));
      //is this is the same as:
	  // dvar_vector tmp 				= elem_prod(moltp(k),elem_prod(Simat(k,1,i),mfexp(*M(k))*natlength_inew(k,i)));
	  
      natlength_new(k,i+1)		=  tmp * len_len(k);
	  
      dvar_vector tmpo 			= elem_prod(moltp(k)*mfexp(-(1-catch_midpt(i))*M(k)), elem_prod(Simat(k,2,i),mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i)));

      natlength_new(k,i+1)		+=  tmpo * len_len(k);

      natlength_iold(k,i+1) 	= mfexp(-(1-catch_midpt(i))*M(k))*(elem_prod(Simat(k,1,i),mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i)) +  elem_prod(Simat(k,2,i),mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i))) - tmp-tmpo;
 
      dvar_vector tmpm 			= elem_prod(moltp_mat(k)*mfexp(-(1-catch_midpt(i))*M_matn(k)), elem_prod(Smat(k,1,i),mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i)));

      natlength_mnew(k,i+1)		=   tmpm * len_len(k);

      dvar_vector tmpmo			= elem_prod(moltp_mat(k)*mfexp(-(1-catch_midpt(i))*M_mato(k)), elem_prod(Smat(k,2,i),mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i)));

      natlength_mnew(k,i+1) 	+=   tmpmo * len_len(k);
 
      natlength_mold(k,i+1)		= (mfexp(-(1-catch_midpt(i))*M_matn(k)) * elem_prod(Smat(k,1,i),mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i))) +  (mfexp(-(1-catch_midpt(i))*M_mato(k)) * elem_prod(Smat(k,2,i),mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i))) - tmpm-tmpmo;

 //this is for estimating the fractino of new shell that move to old shell to fit
//the survey data that is split by immature and mature
      natlength_mnew(k,i+1) 	+= elem_prod(maturity_est(k),natlength_new(k,i+1));
      natlength_inew(k,i+1) 	= elem_prod(1.0-maturity_est(k),natlength_new(k,i+1));

//add recruits for next year to new shell immature
 if(i<=endyr)
      natlength_inew(k,i+1) 	+= mfexp(mean_log_rec(k)+rec_devf(i))*rec_len*proprecn;
 else
      natlength_inew(k,i+1) 	+= FutRec*rec_len*proprecn;
  
//using estimated maturity curve for new shells
      natlength_mat(k,i+1)    	 = natlength_mnew(k,i+1) + natlength_mold(k,i+1);
      natlength_i(k,i+1)      	 = natlength_inew(k,i+1) + natlength_iold(k,i+1);
      natlength(k,i+1)        	 = natlength_mat(k,i+1)  + natlength_i(k,i+1);
      natlength_old(k,i+1)    	 = natlength_mold(k,i+1) + natlength_iold(k,i+1);
      natlength_new(k,i+1)    	 = natlength_inew(k,i+1) + natlength_mnew(k,i+1);
      
      natl_inew_fishtime(k,i)	 = mfexp(-catch_midpt(i)*M(k))*natlength_inew(k,i);
      natl_iold_fishtime(k,i) 	= mfexp(-catch_midpt(i)*M(k))*natlength_iold(k,i);
      natl_mnew_fishtime(k,i)	=  mfexp(-catch_midpt(i)*M_matn(k))*natlength_mnew(k,i);
      natl_mold_fishtime(k,i) 	=  mfexp(-catch_midpt(i)*M_mato(k))*natlength_mold(k,i);
      natl_new_fishtime(k,i)	= natl_inew_fishtime(k,i)+natl_mnew_fishtime(k,i); 
      natl_old_fishtime(k,i) 	= natl_iold_fishtime(k,i)+natl_mold_fishtime(k,i);

      mspbio_matetime(i) 		= (elem_prod(Smat(2,1,i)*mfexp(-(spmo)*M_matn(2)),mfexp(-catch_midpt(i)*M_matn(2))*natlength_mnew(2,i))+elem_prod(Smat(2,2,i)*mfexp(-(spmo)*M_mato(2)),mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i)))*wtm;
      fspbio_matetime(i) 		= (elem_prod(Smat(1,1,i)*mfexp(-(spmo)*M_matn(1)),mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i))+elem_prod(Smat(1,2,i)*mfexp(-(spmo)*M_mato(1)),mfexp(-catch_midpt(i)*M_mato(1))*natlength_mold(1,i)))*wtf(2);
      mspbio_fishtime(i) 		= (natl_mnew_fishtime(2,i)+natl_mold_fishtime(2,i))*wtm;
      fspbio_fishtime(i) 		= (natl_mnew_fishtime(1,i)+natl_mold_fishtime(1,i))*wtf(2);

    if(k>1)
	{
//only want to do this once not for each k
     for(j=1;j<=nlenm;j++)
       {
        emspnum_old_matetime(i) 		+= Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j);
        mspnum_matetime(i) 				+= Smat(2,1,i,j)*mfexp(-(spmo)*M_matn(2))*mfexp(-catch_midpt(i)*M_matn(2))*natlength_mnew(2,i,j) + Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j);
        mspbio_old_matetime(i) 			+= (Smat(2,2,i,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(i)*M_mato(2))*natlength_mold(2,i,j))*wtm(j);
        fspnum_new_matetime(i) 			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j));
        fspbio_new_matetime(i)			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j))*wtf(2,j);
        efspnum_matetime(i) 			+= (Smat(1,1,i,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(i)*M_matn(1))*natlength_mnew(1,i,j)+Smat(1,2,i,j)*mfexp(-(spmo)*M_mato(1))*mfexp(-catch_midpt(i)*M_mato(1))*natlength_mold(1,i,j));
        if(j>15){
          num_males_gt101(i)			+= natl_inew_fishtime(k,i,j) + natl_iold_fishtime(k,i,j) + natl_mnew_fishtime(k,i,j) + natl_mold_fishtime(k,i,j);
          bio_males_gt101(i)			+= (natl_inew_fishtime(k,i,j) + natl_iold_fishtime(k,i,j) + natl_mnew_fishtime(k,i,j) + natl_mold_fishtime(k,i,j))*wtm(j);
        if(j<17){
              num_males_gt101(i)		=num_males_gt101(i)*0.5;
              bio_males_gt101(i)		=bio_males_gt101(i)*0.5;
            }
          }
       }
//effective sp numbers
      efspbio_matetime(i)				=fspbio_matetime(i);
      emspbio_matetime(i)				=mspbio_old_matetime(i);

//for male old shell mating only
      if(emspnum_old_matetime(i)<(efspnum_matetime(i)/mate_ratio))
        efspbio_matetime(i)				=fspbio_matetime(i)*((emspnum_old_matetime(i)*mate_ratio)/efspnum_matetime(i));

//effective sp numbers for new shell females
      efspbio_new_matetime(i)			=fspbio_new_matetime(i);

//for male old shell mating only
      if(emspnum_old_matetime(i)<fspnum_new_matetime(i)/mate_ratio)
        efspbio_new_matetime(i)			=fspbio_new_matetime(i)*((emspnum_old_matetime(i)*mate_ratio)/fspnum_new_matetime(i));

       popn_lmale_new(i) 				= natl_inew_fishtime(2,i)*sel(1,i)+natl_mnew_fishtime(2,i)*sel(1,i);
       popn_lmale_old(i) 				=  natl_iold_fishtime(2,i)*sel(2,i)+natl_mold_fishtime(2,i)*sel(2,i);
       popn_lmale_bio(i)   			 	= elem_prod(natl_new_fishtime(2,i),sel(1,i))*wtm+elem_prod(natl_old_fishtime(2,i),sel(2,i))*wtm;
       popn_lmale(i)       				= popn_lmale_new(i)+popn_lmale_old(i);
       popn_fit_new(i)  				= natl_inew_fishtime(2,i)*sel_fit(1,i)+natl_mnew_fishtime(2,i)*sel_fit(1,i);
       popn_fit_old(i)  				= natl_iold_fishtime(2,i)*sel_fit(2,i)+natl_mold_fishtime(2,i)*sel_fit(2,i);
       popn_fit(i)          			= popn_fit_new(i) + popn_fit_old(i);
      } //end of if(k>1) loop

     if(k<2)
         {
         if(i>=yrs_fish_discf(1) && i<endyr)
		  {
		    popn_disc(1,i)       = (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*sel_discf(i);
          }
         else
		   {
             popn_disc(1,i)       = (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*sel_discf_e;
            }
         
         }
    }



  //isn't this all repeated code?
 if(i ==endyr)
  {
  natl_inew_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M(1))*natlength_inew(1,endyr));
  natl_iold_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M(1))*natlength_iold(1,endyr));
  natl_inew_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M(2))*natlength_inew(2,endyr));
  natl_iold_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M(2))*natlength_iold(2,endyr));
  natl_mnew_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr));
  natl_mold_fishtime(1,endyr) = (mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr));
  natl_mnew_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr));
  natl_mold_fishtime(2,endyr) = (mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr));
  natl_new_fishtime(1,endyr)  = natl_inew_fishtime(1,endyr)+natl_mnew_fishtime(1,endyr);
  natl_old_fishtime(1,endyr)  = natl_iold_fishtime(1,endyr)+natl_mold_fishtime(1,endyr);
  natl_new_fishtime(2,endyr)  = natl_inew_fishtime(2,endyr)+natl_mnew_fishtime(2,endyr);
  natl_old_fishtime(2,endyr)  = natl_iold_fishtime(2,endyr)+natl_mold_fishtime(2,endyr);

  mspbio_matetime(endyr) = (elem_prod(Smat(2,1,endyr)*mfexp(-(spmo)*M_matn(2)),mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr))+elem_prod(Smat(2,2,endyr)*mfexp(-(spmo)*M_mato(2)),mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr)))*wtm;
  fspbio_matetime(endyr) = (elem_prod(Smat(1,1,endyr)*mfexp(-(spmo)*M_matn(1)),mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr))+elem_prod(Smat(1,2,endyr)*mfexp(-(spmo)*M_mato(1)),mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr)))*wtf(2);
  mspbio_fishtime(endyr) = (natl_mnew_fishtime(2,endyr)+natl_mold_fishtime(2,endyr))*wtm;
  fspbio_fishtime(endyr) = (natl_mnew_fishtime(1,endyr)+natl_mold_fishtime(1,endyr))*wtf(2);

  for(j=1;j<=nlenm;j++)
   {
      emspnum_old_matetime(endyr) 	+= Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j);
      mspnum_matetime(endyr) 		+= Smat(2,1,endyr,j)*mfexp(-(spmo)*M_matn(2))*mfexp(-catch_midpt(endyr)*M_matn(2))*natlength_mnew(2,endyr,j) + Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j);
        mspbio_old_matetime(endyr) 	+= (Smat(2,2,endyr,j)*mfexp(-(spmo)*M_mato(2))*mfexp(-catch_midpt(endyr)*M_mato(2))*natlength_mold(2,endyr,j))*wtm(j);
        fspnum_new_matetime(endyr) 	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j));
        fspbio_new_matetime(endyr)	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j))*wtf(2,j);
        efspnum_matetime(endyr) 	+= (Smat(1,1,endyr,j)*mfexp(-(spmo)*M_matn(1))*mfexp(-catch_midpt(endyr)*M_matn(1))*natlength_mnew(1,endyr,j)+Smat(1,2,endyr,j)*mfexp(-(spmo)*M_mato(1))*mfexp(-catch_midpt(endyr)*M_mato(1))*natlength_mold(1,endyr,j));
        if(j>15){
          num_males_gt101(endyr)	+= natl_inew_fishtime(2,endyr,j) + natl_iold_fishtime(2,endyr,j) + natl_mnew_fishtime(2,endyr,j) + natl_mold_fishtime(2,endyr,j);
          bio_males_gt101(endyr)	+= (natl_inew_fishtime(2,endyr,j) + natl_iold_fishtime(2,endyr,j) + natl_mnew_fishtime(2,endyr,j) + natl_mold_fishtime(2,endyr,j))*wtm(j);
        if(j<17){
              num_males_gt101(endyr)=num_males_gt101(endyr)*0.5;
              bio_males_gt101(endyr)=bio_males_gt101(endyr)*0.5;
            }
         }

    }
//effective sp numbers
    efspbio_matetime(endyr)			=fspbio_matetime(endyr);
    emspbio_matetime(endyr)			=mspbio_old_matetime(endyr);

//for male old shell mating only
   if(emspnum_old_matetime(endyr)	>=efspnum_matetime(endyr)/mate_ratio)
    efspbio_matetime(endyr)			=fspbio_matetime(endyr)*((emspnum_old_matetime(endyr)*mate_ratio)/efspnum_matetime(endyr));

//effective sp numbers for new shell females
   efspbio_new_matetime(endyr)		=fspbio_new_matetime(endyr);

//for male old shell mating only
   if(emspnum_old_matetime(endyr)<fspnum_new_matetime(endyr)/mate_ratio)
    efspbio_new_matetime(endyr)		=fspbio_new_matetime(endyr)*((emspnum_old_matetime(endyr)*mate_ratio)/fspnum_new_matetime(endyr));
 
  popn_lmale_new(endyr)  			= natl_inew_fishtime(2,endyr)*sel(1,endyr)+natl_mnew_fishtime(2,endyr)*sel(1,endyr);
  popn_lmale_old(endyr) 			= natl_iold_fishtime(2,endyr)*sel(2,endyr)+natl_mold_fishtime(2,endyr)*sel(2,endyr);
  popn_lmale_bio(endyr)   			= elem_prod(natl_new_fishtime(2,endyr),sel(1,endyr))*wtm+elem_prod(natl_old_fishtime(2,endyr),sel(2,endyr))*wtm;
  popn_lmale(endyr)    				= popn_lmale_new(endyr) + popn_lmale_old(endyr);
  popn_fit_new(endyr)  				= natl_inew_fishtime(2,endyr)*sel_fit(1,endyr)+natl_mnew_fishtime(2,endyr)*sel_fit(1,endyr);
  popn_fit_old(endyr)  				= natl_iold_fishtime(2,endyr)*sel_fit(2,endyr)+natl_mold_fishtime(2,endyr)*sel_fit(2,endyr);
  popn_fit(endyr)    				= popn_fit_new(endyr)+ popn_fit_old(endyr);
  popn_disc(1,endyr) 				= (natl_new_fishtime(1,endyr)+natl_old_fishtime(1,endyr))*sel_discf(yrs_fish_discf(nobs_fish_discf));
  }
  //predicted survey values 
// cout<<"natlength 2"<<endl;
  for(k=1;k<=2;k++)
      {
       if(i<1982){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv1(k));
        }
       if(i>1981 && i<1989){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv2(k));
        }
       if(i>1988){
        totn_srv1(k,i)=(natlength(k,i)*sel_srv3(k));
        if(i==2009){
	        		totn_srv2(1,k)=(natlength(k,i)*sel_srvind(k));
       				totn_srv2(2,k)=(natlength(k,i)*sel_srvnmfs(k));
           		   }
           if(i==2010){
	        		totn_srv10(1,k)=(natlength(k,i)*sel_srv10ind(k));
       				totn_srv10(2,k)=(natlength(k,i)*sel_srv10nmfs(k));
           		   }
            }

        totn_trawl(k,i)= (natl_new_fishtime(k,i)+natl_old_fishtime(k,i))*sel_trawl(k);

   }

    predpop_sexr(i)=sum(natlength(1,i))/(sum(natlength(1,i))+sum(natlength(2,i)));
    fspbio(i)=natlength_mat(1,i)*wtf(2);
    mspbio(i)=natlength_mat(2,i)*wtm;
    if(i<1982)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv1(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv1(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv1(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv1(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv1(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv1(2);
    }
    if(i>1981 && i<1989)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv2(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv2(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv2(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv2(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv2(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv2(2);
    }
    if(i>1988)
	{
     fspbio_srv1(i) = q1*natlength_mat(1,i)*elem_prod(wtf(2),sel_srv3(1));
     mspbio_srv1(i) = q1*natlength_mat(2,i)*elem_prod(wtm,sel_srv3(2));
     fspbio_srv1_num(1,i) = q1*natlength_mnew(1,i)*sel_srv3(1);
     mspbio_srv1_num(1,i) = q1*natlength_mnew(2,i)*sel_srv3(2);
     fspbio_srv1_num(2,i) = q1*natlength_mold(1,i)*sel_srv3(1);
     mspbio_srv1_num(2,i) = q1*natlength_mold(2,i)*sel_srv3(2);
      if(i==2009)
        {
	     fspbio_srv2_ind = natlength_mat(1,i)*elem_prod(wtf(2),sel_srvind(1));
   		 mspbio_srv2_ind = natlength_mat(2,i)*elem_prod(wtm,sel_srvind(2));
		 fspbio_srv2_nmfs = natlength_mat(1,i)*elem_prod(wtf(2),sel_srvnmfs(1));
   		 mspbio_srv2_nmfs = natlength_mat(2,i)*elem_prod(wtm,sel_srvnmfs(2));
         pred_p_srv2_len_ind(1)=elem_prod(sel_srvind(1),natlength(1,i))/(totn_srv2(1,1)+totn_srv2(1,2));
         pred_p_srv2_len_nmfs(1)=elem_prod(sel_srvnmfs(1),natlength(1,i))/(totn_srv2(2,1)+totn_srv2(2,2));
         pred_p_srv2_len_ind(2)=elem_prod(sel_srvind(2),natlength(2,i))/(totn_srv2(1,1)+totn_srv2(1,2));
         pred_p_srv2_len_nmfs(2)=elem_prod(sel_srvnmfs(2),natlength(2,i))/(totn_srv2(2,1)+totn_srv2(2,2));
        }

      if(i==2010)
       {
	     fspbio_srv10_ind = natlength_mat(1,i)*elem_prod(wtf(2),sel_srv10ind(1));
   		 mspbio_srv10_ind = natlength_mat(2,i)*elem_prod(wtm,sel_srv10ind(2));
		 fspbio_srv10_nmfs = natlength_mat(1,i)*elem_prod(wtf(2),sel_srv10nmfs(1));
   		 mspbio_srv10_nmfs = natlength_mat(2,i)*elem_prod(wtm,sel_srv10nmfs(2));
         pred_p_srv10_len_ind(1)=elem_prod(sel_srv10ind(1),natlength(1,i))/(totn_srv10(1,1)+totn_srv10(1,2));
         pred_p_srv10_len_nmfs(1)=elem_prod(sel_srv10nmfs(1),natlength(1,i))/(totn_srv10(2,1)+totn_srv10(2,2));
         pred_p_srv10_len_ind(2)=elem_prod(sel_srv10ind(2),natlength(2,i))/(totn_srv10(1,1)+totn_srv10(1,2));
         pred_p_srv10_len_nmfs(2)=elem_prod(sel_srv10nmfs(2),natlength(2,i))/(totn_srv10(2,1)+totn_srv10(2,2));
	   }
    }

 for(k=1;k<=2;k++)
  {
// this is predicted survey in numbers not biomass-don't adjust by max selectivity 
    if(i<1982)
	{
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv1(k));
     if(k<2)
	 {
	   pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv1(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv1(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv1(k),wtm));
      }

//immature
      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv1(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv1(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
//mature
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv1(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv1(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));

    }

  if(i>1981 && i<1989)
  {
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv2(k));

     if(k<2)
	 {
      pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv2(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv2(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv2(k),wtm));
      }

      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv2(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv2(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv2(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv2(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
    }
	
   if(i>1988)
   {
      pred_srv1(k,i) = q1*elem_prod(natlength(k,i),sel_srv3(k));
     if(k<2)
	 {
      pred_srv1_bioms(k,i) = q1*((natlength_inew(k,i)*elem_prod(sel_srv3(k),wtf(1)))+((natlength_mnew(k,i)+natlength_mold(k,i))*elem_prod(sel_srv3(k),wtf(2))));
      }
     else
	 {
      pred_srv1_bioms(k,i) = q1*(natlength(k,i)*elem_prod(sel_srv3(k),wtm));
      }

      pred_p_srv1_len_new(1,k,i)=elem_prod(sel_srv3(k),natlength_inew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(1,k,i)=elem_prod(sel_srv3(k),natlength_iold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_new(2,k,i)=elem_prod(sel_srv3(k),natlength_mnew(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
      pred_p_srv1_len_old(2,k,i)=elem_prod(sel_srv3(k),natlength_mold(k,i))/(totn_srv1(1,i)+totn_srv1(2,i));
    }
    //next line used to fix q1 to 1.0 - problem is if you start from a bin file, even if the bounds
    // are set different in the tpl file the program will take to value from the bin file and use that 
    //   pred_srv1(i)=1.0*(natage(i)*elem_prod(sel_srv1,wt));
    explbiom(i)=((natlength_old(2,i)*elem_prod(sel(2,i),wtm))+(natlength_new(2,i)*elem_prod(sel(1,i),wtm)));
    popn(i)+=sum(natlength(k,i));
   }  //end of k loop

    pred_bio(i)+=natlength_inew(1,i)*wtf(1)+(natlength_mnew(1,i)+natlength_mold(1,i))*wtf(2)+(natlength_inew(2,i)+natlength_mnew(2,i)+natlength_mold(2,i))*wtm;

 
  depletion = pred_bio(endyr) / pred_bio(styr);
  fspbios=fspbio;
  mspbios=mspbio_matetime;

    legal_males(i)=0.;
    legal_srv_males(i)=0.;
    legal_srv_males_n(i)=0.0;
    legal_srv_males_o(i)=0.0;
// legal is >102mm take half the numbers in the 100-105 bin
        legal_males(i)=0.5*natlength(2,i,16);
        legal_srv_males_n(i)=0.5*natlength_new(2,i,16);
        legal_srv_males_o(i)=0.5*natlength_old(2,i,16);
        legal_males_bio(i)=legal_males(i)*wtm(16);

        if(i<1982)
		{
        legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv1(2,16);
        legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv1(2,16);
        legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv1(2,16);
        }
        if(i>1981 && i<1989)
		{
         legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv2(2,16);
         legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv2(2,16);
         legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv2(2,16);
        }
        if(i>1988)
		{
         legal_srv_males(i)=0.5*natlength(2,i,16)*sel_srv3(2,16);
         legal_srv_males_n(i)=0.5*natlength_new(2,i,16)*sel_srv3(2,16);
         legal_srv_males_o(i)=0.5*natlength_old(2,i,16)*sel_srv3(2,16);
        }
 
        legal_srv_males_bio(i)=legal_srv_males(i)*wtm(16);
    
	for(j=17;j<=nlenm;j++)
    {
        legal_males(i)+=natlength(2,i,j);
        legal_males_bio(i)+=natlength(2,i,j)*wtm(j);
        if(i<1982)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv1(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv1(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv1(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv1(2,j)*wtm(j);
         }
        if(i>1981 && i<1989)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv2(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv2(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv2(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv2(2,j)*wtm(j);
          }
        if(i>1988)
		{
          legal_srv_males(i)+=natlength(2,i,j)*sel_srv3(2,j);
          legal_srv_males_n(i)+=natlength_new(2,i,j)*sel_srv3(2,j);
          legal_srv_males_o(i)+=natlength_old(2,i,j)*sel_srv3(2,j);
          legal_srv_males_bio(i)+=natlength(2,i,j)*sel_srv3(2,j)*wtm(j);
          }

    }
    // legal_malesd=legal_males;
    // recf_sd=mfexp(mean_log_rec(1)+rec_dev(1));
    // recm_sd=mfexp(mean_log_rec(2)+rec_dev(2));

//==============================================================================
FUNCTION get_catch_at_len
 //cout<<" begin catch at len"<<endl;

    pred_catch.initialize();
    pred_catch_ret.initialize();
    pred_catch_gt101.initialize();
    pred_catch_no_gt101.initialize();

	// take out catch all at once - survey is start, catch_midpt has midpoint of the
    // catch(fishing season weighted by the catch)
    // then rest of nat mort and growth.
  for (i=styr;i<=endyr;i++)
  { 
   for(k=1;k<=2;k++)
    {      
    for (j = 1 ; j<= nlenm; j++)
    {
        if(k>1)
        {
        if(j>15)
		{
             pred_catch_no_gt101(i) += (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + 
                                         (Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j))+ 
                                          (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+
                                           (Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
             pred_catch_gt101(i)+= (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + 
                                         (Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j))+ 
                                          (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+
                                           (Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j)) * wtm(j);
           if(j<17)
		   {
               pred_catch_gt101(i)=pred_catch_gt101(i)*0.5;
               pred_catch_no_gt101(i)=pred_catch_no_gt101(i)*0.5;
            }
        }
        catch_lmale_new(i,j) = (Fimat(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j))+(Fmat(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j)); 
        catch_lmale_old(i,j) = (Fimat(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))+(Fmat(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
        catch_lmale(i,j)= catch_lmale_new(i,j)+catch_lmale_old(i,j);
        pred_catch(i) += catch_lmale(i,j)*wtm(j);
        catch_male_ret_new(i,j) = (Fimat_ret(1,i,j)/(Fimat(1,i,j)+Fdisct(2,i,j)))*natl_inew_fishtime(2,i,j)*(1-Simat(2,1,i,j)) + (Fmat_ret(1,i,j)/(Fmat(1,i,j)+Fdisct(2,i,j)))*natl_mnew_fishtime(2,i,j)*(1-Smat(2,1,i,j)); 
        catch_male_ret_old(i,j) = (Fimat_ret(2,i,j)/(Fimat(2,i,j)+Fdisct(2,i,j)))*natl_iold_fishtime(2,i,j)*(1-Simat(2,2,i,j))  + (Fmat_ret(2,i,j)/(Fmat(2,i,j)+Fdisct(2,i,j)))*natl_mold_fishtime(2,i,j)*(1-Smat(2,2,i,j));
        catch_male_ret(i,j)=catch_male_ret_new(i,j)+catch_male_ret_old(i,j);
        pred_catch_ret(i) += catch_male_ret(i,j)*wtm(j);
        
       } //end k
    catch_discp(1,i,j) = (Fdiscf(i,j)/(Fdiscf(i,j)+Fdisct(1,i,j))) *(natl_new_fishtime(1,i,j)+natl_old_fishtime(1,i,j))*(1-mfexp(-1.0*(Fdiscf(i,j)+Fdisct(1,i,j))));
    // pred_catch_trawl
    } 
    
	if(k<2)
	{
      pred_catch_disc(k,i)=catch_discp(k,i)*wtf(2);
    }
    else
	{
      pred_catch_disc(k,i)=catch_discp(k,i)*wtm;
    }
    //  cout<<" 3 catlen"<<endl;
   }
   pred_catch_trawl(i)= (natl_new_fishtime(1,i)+natl_old_fishtime(1,i))*elem_prod(1-mfexp(-1.0*Fdisct(1,i)),wtf(2))+ (natl_new_fishtime(2,i)+natl_old_fishtime(2,i))*elem_prod(1-mfexp(-1.0*Fdisct(2,i)),wtm);
   pred_p_fish_fit(1,i)=(elem_prod(sel_fit(1,i),natl_new_fishtime(2,i)))/(popn_fit(i));
   pred_p_fish_fit(2,i)=(elem_prod(sel_fit(2,i),natl_old_fishtime(2,i)))/(popn_fit(i));
   pred_p_fish(1,i)=(elem_prod(sel(1,i),natl_new_fishtime(2,i))/(popn_lmale(i)));
   pred_p_fish(2,i)=(elem_prod(sel(2,i),natl_old_fishtime(2,i))/(popn_lmale(i)));

   if(i>=yrs_fish_discf(1) && i<endyr)
   {
   pred_p_fish_discf(i)=elem_prod(sel_discf(i),(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/popn_disc(1,i);
   }
   else
   {
   pred_p_fish_discf(i)=elem_prod(sel_discf_e,(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/popn_disc(1,i);
   }
   pred_p_trawl(1,i)=elem_prod(sel_trawl(1),(natl_new_fishtime(1,i)+natl_old_fishtime(1,i)))/(totn_trawl(1,i)+totn_trawl(2,i));
   pred_p_trawl(2,i)=elem_prod(sel_trawl(2),(natl_new_fishtime(2,i)+natl_old_fishtime(2,i)))/(totn_trawl(1,i)+totn_trawl(2,i));
   preds_sexr(i)=totn_srv1(1,i)/(totn_srv1(1,i)+totn_srv1(2,i));
  }
  

//==============================================================================
FUNCTION evaluate_the_objective_function
  len_likeyr.initialize();
  len_like.initialize();
  len_like_srv.initialize();
  len_like10_ind.initialize();
  len_like10_nmfs.initialize();
  sel_like=0.;
  fpen=.0;
  like_initsmo=0.0;
  ghl_like=0.0;
  rec_like=.0;
  surv_like=.0;
  surv2_like=0.0;
  surv3_like=0.0;
  catch_like1=.0;
  catch_like2=.0;
  catch_likef=.0;
  catch_liket=.0;
  like_q=0.0;
  like_natm.initialize();
  like_af.initialize();
  like_bf.initialize();
  like_bf1.initialize();
  like_mmat.initialize();
  like_am.initialize();
  like_bm.initialize();
  sexr_like.initialize();
  like_natm.initialize();
  sumrecf.initialize();
  sumrecm.initialize();
  sel_like_50m.initialize();
   sel_avg_like.initialize();
  like_srvsel.initialize();
  like_initnum.initialize();
  Fout.initialize();
  len_like_ex=0.0;
  selsmo_like=0.0;
  //  catch_likef=0.;
  growth_like=0.0;
  f.initialize();
  largemale_like=0.0;

 if (active(rec_devf))
   {
    rec_like = like_wght_recf*norm2(rec_devf)+like_wght_rec*norm2(rec_devm);
    Fout(1)=rec_like;
    f += rec_like;
   if(active(rec_devm))
	{
      for(i=styr;i<endyr-retro_years;i++)
        sexr_like += like_wght_sexr*square((mean_log_rec(1)+rec_devf(i))-(mean_log_rec(2)+rec_devm(i)));
    }

    }  //end of active rec_devf

 //constraint on intial numbers in small length bins for old shell males
   for(i=1;i<=5;i++)
    {
     like_initnum += old_shell_constraint  * square(mfexp(mnatlen_styr(2,i)));
     }
     f+= like_initnum;
    Fout(2)=like_initnum;

 //jim said take the log--CSS if only adding the likecomp when active...why is it commented out?
  if(active(log_sel50_dev_mn))
   {
    sel_like_50m=like_wght_sel50*norm2(first_difference((log_avg_sel50_mn+log_sel50_dev_mn)(styr,endyr-1)));
    sel_like_50m+=like_wght_sel50*norm2(first_difference((log_avg_sel50_mo+log_sel50_dev_mo)(styr,endyr-1)));
 //   f +=sel_like_50m;
  }

  int ii;
  int ij;
  int ik;
  int mat;
  
    //==========Retained fishery lengths likelihood component===============
  for (i=1; i <= (nobs_fish - retro_years); i++)
  {
    ii=yrs_fish(i);
   for (j=1; j<=nlenm; j++)
    len_like(1)-=nsamples_fish(1,i)*(obs_p_fish_ret(1,i,j)+obs_p_fish_ret(2,i,j))*log(pred_p_fish_fit(2,ii,j)+pred_p_fish_fit(1,ii,j)+p_const); //old and new together

   for(k=1;k<=2;k++)
    effn_fish_ret(k,ii)=1/(norm2(pred_p_fish_fit(k,ii)-obs_p_fish_ret(k,i))/(pred_p_fish_fit(1,ii)*(1-pred_p_fish_fit(1,ii))+pred_p_fish_fit(2,ii)*(1-pred_p_fish_fit(2,ii))));
  }

 for (i=1; i <= (nobs_fish_discm - retro_years); i++)
  {
    ij=yrs_fish_discm(i);
   for (j=1; j<=nlenm; j++)
    len_like(2)-=disclen_mult(1)*nsamples_fish_discm(1,i)*(obs_p_fish_tot(1,i,j)+obs_p_fish_tot(2,i,j))*log(pred_p_fish(1,ij,j)+pred_p_fish(2,ij,j)+p_const);
   for(k=1;k<=2;k++)
    effn_fish_tot(k,ij)=1/(norm2(pred_p_fish(k,ij)-obs_p_fish_tot(k,i))/(pred_p_fish(1,ij)*(1-pred_p_fish(1,ij))+pred_p_fish(2,ij)*(1-pred_p_fish(2,ij))));
  }
 
  //==========Discards female lengths likelihood component=============== 
  for (i=1; i <= (nobs_fish_discf - retro_years); i++)
  {
    ik=yrs_fish_discf(i);
   for (j=1; j<=nlenm; j++)
	len_like(3)-=nsamples_fish_discf(i)*(obs_p_fish_discf(i,j))*log(pred_p_fish_discf(ik,j)+p_const);
  }
  
  //==========Trawl fishery lengths likelihood component===============
  for (i=1; i <= (nobs_trawl - retro_years); i++)
  {
    ij=yrs_trawl(i);
  //trawlfishery length likelihood sur
    for(k=1;k<=2;k++)
      for (j=1; j<=nlenm; j++)
       len_like(5)-=nsamples_trawl(k,i)*(obs_p_trawl(k,i,j))*log(pred_p_trawl(k,ij,j)+p_const);
  }

  //add the offset to the likelihood   
  len_like(1)-=offset(1);
  len_like(2)-=offset(2);
  len_like(3)-=offset(3);
  len_like(5)-=offset(5);

 for (j=1; j<=nlenm; j++)
   {
// mature and immature together
// don't fit the first four bins of BSFRF survey length freq data - need to rescale obs proportions when leaving out first four bins
   	len_like(6)-=1.*nsamples_srv2_length(1,1,1)*((obs_p_srv2_len(1,1,1,j)+obs_p_srv2_len(1,1,2,j))) *log(pred_p_srv2_len_ind(1,j)+p_const);   			  	
    len_like(6)-=1.*nsamples_srv2_length(1,2,1)*((obs_p_srv2_len(1,2,1,j)+obs_p_srv2_len(1,2,2,j))) *log(pred_p_srv2_len_ind(2,j)+p_const);
    len_like(7)-=1.*nsamples_srv2_length(2,1,1)*(obs_p_srv2_len(2,1,1,j)+obs_p_srv2_len(2,1,2,j)) *log(pred_p_srv2_len_nmfs(1,j)+p_const);
    len_like(7)-=1.*nsamples_srv2_length(2,2,1)*(obs_p_srv2_len(2,2,1,j)+obs_p_srv2_len(2,2,2,j)) *log(pred_p_srv2_len_nmfs(2,j)+p_const);
//2010 bsfrf survey length data
   	len_like10_ind-=1.*nsamples_srv10_length(1,1,1)*((obs_p_srv10_len(1,1,1,j)+obs_p_srv10_len(1,1,2,j))) *log(pred_p_srv10_len_ind(1,j)+p_const);   			  	
    len_like10_ind-=1.*nsamples_srv10_length(1,2,1)*((obs_p_srv10_len(1,2,1,j)+obs_p_srv10_len(1,2,2,j))) *log(pred_p_srv10_len_ind(2,j)+p_const);
    len_like10_nmfs-=1.*nsamples_srv10_length(2,1,1)*(obs_p_srv10_len(2,1,1,j)+obs_p_srv10_len(2,1,2,j)) *log(pred_p_srv10_len_nmfs(1,j)+p_const);
    len_like10_nmfs-=1.*nsamples_srv10_length(2,2,1)*(obs_p_srv10_len(2,2,1,j)+obs_p_srv10_len(2,2,2,j)) *log(pred_p_srv10_len_nmfs(2,j)+p_const);
    }               
  
   len_like(6)-=offset_srv2(1);
   len_like(7)-=offset_srv2(2);
   len_like10_ind-=offset_srv10(1);
   len_like10_nmfs-=offset_srv10(2);

 //===========survey lengths likelihood components==================
   for(k=1;k<=2;k++)  //sex
   {
    for (i=1; i <=(nobs_srv1_length - retro_years); i++)
    {
      ii=yrs_srv1_length(i);
      //survey likelihood 
      for (j=1; j<=nlenm; j++)
      {
 //  obs(maturity, SC, sex, year), pred(maturity,sex, year)
 // this is for mature new and old shell together
       len_like(4)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j)+obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+pred_p_srv1_len_old(2,k,ii,j)+p_const);
 //immature new and old together
        len_like(4)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
 //immature new and old together in likelihood indices are (mat,shell,sex,year,length)
        len_likeyr(1,1,k,i)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
        len_like_srv(1,2,k)=0.0;
        len_like_srv(1,1,k)-=nsamples_srv1_length(1,1,k,i)*(obs_p_srv1_len(1,1,k,i,j)+obs_p_srv1_len(1,2,k,i,j))*log(pred_p_srv1_len_new(1,k,ii,j)+pred_p_srv1_len_old(1,k,ii,j)+p_const);
     //mature
        len_likeyr(2,1,k,i)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+p_const);
        len_likeyr(2,2,k,i)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_old(2,k,ii,j)+p_const);
        len_like_srv(2,1,k)-=nsamples_srv1_length(2,1,k,i)*(obs_p_srv1_len(2,1,k,i,j))*log(pred_p_srv1_len_new(2,k,ii,j)+p_const);
        len_like_srv(2,2,k)-=nsamples_srv1_length(2,2,k,i)*(obs_p_srv1_len(2,2,k,i,j))*log(pred_p_srv1_len_old(2,k,ii,j)+p_const);
    }  //j loop     
    effn_srv1(1,1,k,ii)=1./(norm2(pred_p_srv1_len_new(1,k,ii)-obs_p_srv1_len(1,1,k,i))/(pred_p_srv1_len_new(1,1,ii)*(1-pred_p_srv1_len_new(1,1,ii))+pred_p_srv1_len_new(1,2,ii)*(1-pred_p_srv1_len_new(1,2,ii))+pred_p_srv1_len_old(1,1,ii)*(1-pred_p_srv1_len_old(1,1,ii))+pred_p_srv1_len_old(1,2,ii)*(1-pred_p_srv1_len_old(1,2,ii))));
 
  if(k>1) //what is this if statement about?? if and else have the same statement inside...
	{
    effn_srv1(1,2,k,ii)=0.0;
    }
    else
	{
    //effn is zero for old shell immature females
    effn_srv1(1,2,k,ii)=0.0;
    }

    effn_srv1(2,1,k,ii)=1./(norm2(pred_p_srv1_len_new(2,k,ii)-(obs_p_srv1_len(2,1,k,i)))/(pred_p_srv1_len_new(2,1,ii)*(1-pred_p_srv1_len_new(2,1,ii))+pred_p_srv1_len_new(2,2,ii)*(1-pred_p_srv1_len_new(2,2,ii))+pred_p_srv1_len_old(2,1,ii)*(1-pred_p_srv1_len_old(2,1,ii))+pred_p_srv1_len_old(2,2,ii)*(1-pred_p_srv1_len_old(2,2,ii))));
    effn_srv1(2,2,k,ii)=1./(norm2(pred_p_srv1_len_old(2,k,ii)-(obs_p_srv1_len(2,2,k,i)))/(pred_p_srv1_len_new(2,1,ii)*(1-pred_p_srv1_len_new(2,1,ii))+pred_p_srv1_len_new(2,2,ii)*(1-pred_p_srv1_len_new(2,2,ii))+pred_p_srv1_len_old(2,1,ii)*(1-pred_p_srv1_len_old(2,1,ii))+pred_p_srv1_len_old(2,2,ii)*(1-pred_p_srv1_len_old(2,2,ii))));
 
    } // year loop

  }  //k (sex) loop

 
//=============added likelihood for smoothness of initial year length comp=============
   like_initsmo = init_yr_len_smooth_f *norm2(first_difference(mnatlen_styr(1)))+init_yr_len_smooth_m *norm2(first_difference(first_difference(mnatlen_styr(2))))+norm2(first_difference(fnatlen_styr(1)))+norm2(first_difference(fnatlen_styr(2)));
   Fout(30)=like_initsmo;
   f += like_initsmo;
   
  //===================extra weight for start year length comp.==========================   
  //==immature males of both shel conditions
   len_like(4)-=nsamples_srv1_length(1,1,2,1)*(obs_p_srv1_len(1,1,2,1))*log(pred_p_srv1_len_new(1,2,styr)+p_const);
   if(moltp(2,nlenm)<0.99)
    len_like(4)-=nsamples_srv1_length(1,2,2,1)*(obs_p_srv1_len(1,2,2,1))*log(pred_p_srv1_len_old(1,2,styr)+p_const);

 //  obs(maturity, SC, sex, year), pred(maturity,sex, year)

    multi=1.0;
	//==mature both sexes and shell conditions
    len_like_ex-=multi*nsamples_srv1_length(2,1,2,1)*(obs_p_srv1_len(2,1,2,1))*log(pred_p_srv1_len_new(2,2,styr)+p_const);
    len_like_ex-= multi*nsamples_srv1_length(2,2,2,1)*(obs_p_srv1_len(2,2,2,1))*log(pred_p_srv1_len_old(2,2,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(2,1,1,1)*(obs_p_srv1_len(2,1,1,1))*log(pred_p_srv1_len_new(2,1,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(2,2,1,1)*(obs_p_srv1_len(2,2,1,1))*log(pred_p_srv1_len_old(2,1,styr)+p_const);
    len_like_ex-=multi*nsamples_srv1_length(1,1,1,1)*(obs_p_srv1_len(1,1,1,1)+obs_p_srv1_len(1,2,1,1))*log(pred_p_srv1_len_new(1,1,styr)+pred_p_srv1_len_old(1,1,styr)+p_const);
   
   Fout(25)=len_like_ex;
   f+=len_like_ex;
   
//put this back in when doing fit by length proportions
  len_like(4)-=offset(4);

 for(mat=1;mat<=2;mat++)
  for(k=1;k<=1;k++)
    len_like_srv(mat,k,1)-=offset_srv(mat,k,1);

  if(active(selsmo10ind))
  {
    selsmo_like= selsmo_wght *norm2(first_difference(first_difference(selsmo10ind)));
    f+=selsmo_like;
    Fout(28)+=selsmo_like;
   }
  if(active(selsmo09ind))
  {
    selsmo_like= selsmo_wght *norm2(first_difference(first_difference(selsmo09ind)));
    f+=selsmo_like;
    Fout(28)+=selsmo_like;
   }
   if(active(log_dev_50f))
   {
    discf_like= femSel_wght *norm2(first_difference(first_difference(log_dev_50f)));
    f+=discf_like;
    Fout(29)=discf_like;
    }
 
  f+=like_wght(1)*len_like(1);     // Retained fishery
  Fout(3)=like_wght(1)*len_like(1);
  f+=like_wght(2)*len_like(2);     // Total (ret+disc)
  Fout(4)=like_wght(2)*len_like(2);

  f+=like_wght(3)*len_like(3);     // Female
  Fout(5)=like_wght(3)*len_like(3);
  f+=like_wght(4)*len_like(4);     // Survey
  Fout(6)=like_wght(4)*len_like(4);
  f+=like_wght(7)*len_like(5);     // Trawl
  Fout(7)=like_wght(7)*len_like(5);
  f+=len_like(6);                 // Industry Survey BSFRF length
  Fout(8)=len_like(6);
  f+=len_like(7);                // Industry survey NMFS length
  Fout(9)=len_like(7);
  f+=len_like10_ind;                 // Industry 2010 Survey BSFRF length
  Fout(26)=len_like10_ind;
  f+=len_like10_nmfs;                // Industry 2010 survey NMFS length
  Fout(27)=len_like10_nmfs;
  
 like_natm=0.0;
 if(active(Mmult))
 {  
  like_natm   = natm_mult_wght  * square((Mmult    - 1.0)    / natm_mult_var);
  f += like_natm;
  Fout(10)=like_natm;
 }
 
 if(active(Mmult_imat))
 {  
  like_natm   += natm_Immult_wght  * square((Mmult_imat    - 1.0)    / natm_Immult_var);
  f += like_natm;
  Fout(10)+=like_natm;
 }

   like_mat=0.0;
  if(active(mateste))
   {
     like_mat = smooth_mat_wght *norm2(first_difference(first_difference(mateste)));
     // like_mat += mat_est_wght *norm2(maturity_est(2)-maturity_logistic)/(mat_est_vs_obs_sd *mat_est_vs_obs_sd );
     f += like_mat;
     Fout(11)+=like_mat;
    }
 
 if(active(matestfe))
    {
     like_mat += smooth_mat_wght_f *norm2(first_difference(first_difference(matestfe)));
     // like_mat += 1.0*norm2(maturity_est(1)-maturity_average(1))/(mat_est_vs_obs_sd_f *mat_est_vs_obs_sd_f );
     f += like_mat;
     Fout(11)+=like_mat;    
    }

// use growth data to estimate parameters 
   if(active(am))
   {
	 for(i=1;i<=nobs_growm ;i++)
	  like_am+= growth_data_wght_m *pow(malegrowdaty(i)-((am+bm*malegrowdatx(i))*(1-cumd_norm((malegrowdatx(i)-deltam)/st_gr))+((am+(bm-b1)*deltam)+b1*malegrowdatx(i))*(cumd_norm((malegrowdatx(i)-deltam)/st_gr))),2);       
      f += like_am;
      Fout(12)=like_am;
    }
	
   if(active(af))
   {
	for(i=1;i<=nobs_growf ;i++)
     like_af+= growth_data_wght_f *pow(femalegrowdaty(i)-((af+bf*femalegrowdatx(i))*(1-cumd_norm((femalegrowdatx(i)-deltaf)/st_gr))+((af+(bf-bf1)*deltaf)+bf1*femalegrowdatx(i))*(cumd_norm((femalegrowdatx(i)-deltaf)/st_gr))),2);
	 f += like_af;
     Fout(13)=like_af;
   }
  
  // ==================Fit to indices===============================
  // (lognormal) weight each years estimate by 1/(2*variance) - use cv of biomass in sqrt(log(cv^2+1)) as sd of log(biomass) 
 for(i=1;i<=nobs_srv1;i++)
  {
    biom_tmp(1,yrs_srv1(i))=fspbio_srv1(yrs_srv1(i));
    biom_tmp(2,yrs_srv1(i))=mspbio_srv1(yrs_srv1(i));
  }
//likelihood for survey biomass by sex
//  This loop allows for increasing the weight onthe survey biomass through a constant specified in the control file (should be set to 1)
  for(i=1;i<=nobs_srv1;i++)
    {
     cv_srv1(1,yrs_srv1(i))=like_wght(5)*cv_srv1o(1,i);
     cv_srv1(2,yrs_srv1(i))=like_wght_mbio*cv_srv1o(2,i);
     cv_srv1_nowt(1,yrs_srv1(i))=cv_srv1o(1,i);
     cv_srv1_nowt(2,yrs_srv1(i))=cv_srv1o(2,i);
    }

//this fits mature biomass separate male and female

    // surv_like += norm2(elem_div( log(obs_srv1_spbiom(1)(yrs_srv1)+p_const )-log(biom_tmp(1)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(1)(yrs_srv1),cv_srv1(1)(yrs_srv1))+1.0))));
    // surv_like += norm2(elem_div( log(obs_srv1_spbiom(2)(yrs_srv1)+p_const )-log(biom_tmp(2)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(2)(yrs_srv1),cv_srv1(2)(yrs_srv1))+1.0))));

  for(i=1;i<=(nobs_srv1 - retro_years);i++)
  {
   surv_like += pow((log(obs_srv1_spbiom(1)(yrs_srv1(i))+p_const )-log(biom_tmp(1)(yrs_srv1(i))+p_const ))/(sqrt(2)*sqrt(log((cv_srv1(1)(yrs_srv1(i))*cv_srv1(1)(yrs_srv1(i)))+1.0))),2);
   surv_like += pow((log(obs_srv1_spbiom(2)(yrs_srv1(i))+p_const )-log(biom_tmp(2)(yrs_srv1(i))+p_const ))/(sqrt(2)*sqrt(log((cv_srv1(2)(yrs_srv1(i))*cv_srv1(2)(yrs_srv1(i)))+1.0))),2);
  }
	
//industry survey
//2009
    surv2_like = pow(( log(obs_srv2_spbiom(1,1)+p_const)-log(fspbio_srv2_ind+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(1,1)*obs_srv2_cv(1,1))+1.0))),2.0);
    surv2_like += extra_wght_ind_m *pow((log(obs_srv2_spbiom(1,2)+p_const)-log(mspbio_srv2_ind+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(1,2)*obs_srv2_cv(1,2))+1.0))),2.0);
    surv3_like = pow(( log(obs_srv2_spbiom(2,1)+p_const)-log(fspbio_srv2_nmfs+p_const))/ (sqrt(2)*sqrt(log((obs_srv2_cv(2,1)*obs_srv2_cv(2,1))+1.0))),2.0);
    surv3_like += pow(( log(obs_srv2_spbiom(2,2)+p_const)-log(mspbio_srv2_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv2_cv(2,2)*obs_srv2_cv(2,2))+1.0))),2.0);
   f+=1.0*surv2_like;
   f+=1.0*surv3_like;			  	 
   Fout(14)=1.0*surv2_like;
   Fout(15)=1.0*surv3_like;

//industry survey
//2010
   surv10_like = pow(( log(obs_srv10_spbiom(1,1)+p_const)-log(fspbio_srv10_ind+ p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(1,1)*obs_srv10_cv(1,1))+1.0))),2.0);
   surv10_like += extra_wght_ind_m *pow((log(obs_srv10_spbiom(1,2)+p_const)-log(mspbio_srv10_ind+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(1,2)*obs_srv10_cv(1,2))+1.0))),2.0);
   surv10nmfs_like = pow(( log(obs_srv10_spbiom(2,1)+p_const)-log(fspbio_srv10_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(2,1)*obs_srv10_cv(2,1))+1.0))),2.0);
   surv10nmfs_like += pow(( log(obs_srv10_spbiom(2,2)+p_const)-log(mspbio_srv10_nmfs+p_const))/(sqrt(2)*sqrt(log((obs_srv10_cv(2,2)*obs_srv10_cv(2,2))+1.0))),2.0);
   f+=1.0*surv10_like;
   f+=1.0*surv10nmfs_like;			  	 
   Fout(23)=1.0*surv10_like;
   Fout(24)=1.0*surv10nmfs_like;

   int yrc;
   cpue_pred = cpueq * legal_males;
 
   // calculate the observed large males
   obs_lmales.initialize();
   obs_lmales_bio.initialize();
   for(i=1;i<=nobs_srv1_length;i++)
    {
    //take 1/2 of the 100-104 bin, 
          obs_lmales(i)=0.5*obs_srv1_num(2,yrs_srv1_length(i),16);
          obs_lmales_bio(i)=obs_lmales(i)*wtm(16);
      for(j=17;j<=nlenm;j++)
        {
          obs_lmales(i)+=obs_srv1_num(2,yrs_srv1_length(i),j);
          obs_lmales_bio(i)+=obs_srv1_num(2,yrs_srv1_length(i),j)*wtm(j);
        }
     }
 
  //fishery cpue likelihood
    cpue_like=0.0;
   if(active(cpueq))
   { 
   for(yrc=styr;yrc<=(endyr-1- retro_years);yrc++)
    cpue_like += pow(((log(cpue(yrc+1)+1e-9)-log(cpue_pred(yrc)+1e-9))/(sqrt(2)* cpue_cv)),2.0);
    //cpue_like += norm2(elem_div( log(obs_srv1_spbiom(2)(yrs_srv1)+p_const )-log(biom_tmp(2)(yrs_srv1)+p_const ),sqrt(2)*sqrt(log(elem_prod(cv_srv1(2)(yrs_srv1),cv_srv1(2)(yrs_srv1))+1.0))));

    f+=cpue_like;
    Fout(16)=cpue_like;
   }
   
//catch likelihoods
 //don't include last year as that would be endyr+1 fishery season
   catch_like1 = norm2((obs_catchtot_biom(1992,endyr-1-retro_years)-catch_ret(1992,endyr-1-retro_years)+p_const)-(pred_catch(1992,endyr-1-retro_years)-pred_catch_ret(1992,endyr-1- retro_years)+p_const));
//total catch likelihood
//female catch
    catch_likef = norm2((obs_catchdf_biom)-(pred_catch_disc(1)));
    catch_likef += smooth_disc_catch *norm2(first_difference(pred_catch_disc(1)));
//retained catch
    catch_like2 = norm2((catch_ret(styr,endyr-1-retro_years))-(pred_catch_ret(styr,endyr-1-retro_years)));
//groundfish catch
    catch_liket = norm2((obs_catcht_biom)-(pred_catch_trawl));

  f += like_wght(6)*1.*catch_like2;
  Fout(17)= like_wght(6)*1.*catch_like2;
  f += wght_total_catch*1.*catch_like1;
  Fout(18)= wght_total_catch*1.*catch_like1;
   
  f += like_wght(6)*catch_liket;
  Fout(19)= like_wght(6)*catch_liket;
  f += like_wght(6)*disc_catch_wght_fem *catch_likef;
  Fout(20)= like_wght(6)*disc_catch_wght_fem *catch_likef;
  f += surv_like;
  Fout(21)= surv_like;

//  f += ghl_like;
  fpen.initialize();
 if(f_penalties>0)
 {
// f of 1.15 is about an exploitation rate of 0.68(.58 increased by nat. mort. 7 months).
  if (current_phase()<2)
    fpen = like_wght_fph1*norm2(mfexp(fmort_dev+log_avg_fmort)-1.15);
  else
    fpen = like_wght_fph2*norm2(mfexp(fmort_dev+log_avg_fmort)-1.15);
 }
 
  if (active(fmortt_dev_era1))
  {
   fpen += norm2(fmortt_dev_era1(styr,1991));
   fpen += norm2(fmortt_dev_era2(1992,endyr-1));
  }
  
  if (active(fmort_dev))
  {
 // this keeps fmorts from going too high if like_wght_fdev is high...THIS ACTUALLY CONSTRAINS F TO ZERO
    fpen += like_wght_fdev*norm2(mfexp(fmort_dev(styr,endyr-1)+log_avg_fmort));
    fpen += 1.0*norm2(fmortdf_dev(styr,endyr-1));
   }
  Fout(22)= fpen;
  f+=fpen;

  call_no += 1;
  //cout <<"Likes = "<< Fout << endl;
  //cout <<"phase = "<< current_phase() << " call = " << call_no << " Total Like = " << f << endl;

// ========================y==================================================   
FUNCTION get_fut_mortality  
// specify all of the fishing mortality sources into the future
// allow for another variable (FutMort) to be passed when this function is called in Find F35 and FindOFL
// set non-directed fishing to the average, set directed fishing to an input value
 for(i=ipass;i<=endyr+Nproj-1;i++)
 {
	if(IsB0==1)
	{
      fmort(i)=0.0; 
	  fmortdf(i)=0.0;
      fmortt(i)=0.0;
	}
	else
	{
	fmort(i)=FutMort;  
	fmortdf(i)=mfexp(log_avg_fmortdf);
    fmortt(i)=mfexp(log_avg_fmortt);
	}	
 }
 
  for (i=ipass;i<=endyr+Nproj;i++)
 {
   if(i>=yrs_fish_discf(1) && i<endyr)
   {
    Fdiscf(i)=sel_discf(i)*fmortdf(i);
   }
   else
   {
    sel_discf_e=(1./(1.+mfexp(-1.*fish_disc_slope_f*(length_bins-mfexp(fish_disc_sel50_f)))));
    Fdiscf(i)=sel_discf_e*fmortdf(i);
    }
   Fdisct(1,i)=sel_trawl(1)*fmortt(i);   
   Fdisct(2,i)=sel_trawl(2)*fmortt(i);   
   
   for(k=1;k<=2;k++) //over new (k=1) and old (k=2) shell...
   { 
     F(k,i) = sel(k,i)*fmort(i);       
     F_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fmat(k,i) = sel(k,i)*fmort(i);       
     Fmat_ret(k,i)=sel_fit(k,i)*fmort(i);
     Fimat(k,i) = sel(k,i)*fmort(i);      
     Fimat_ret(k,i)= sel_fit(k,i)*fmort(i);
     Ftot(1,k,i)=Fdiscf(i) + Fdisct(1,i);
     S(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Simat(1,k,i)=mfexp(-1.0*Ftot(1,k,i));
     Smat(2,k,i)=mfexp(-1.0*(Fmat(k,i)+Fdisct(2,i)));
     Simat(2,k,i)=mfexp(-1.0*(Fimat(k,i)+Fdisct(2,i)));

   } 
 }

  //===============================================================================
FUNCTION Find_F35  
 int icnt;

//Find Bzero
  IsB0 = 1;
  FutRec = 100000;
  FutMort = 0;
  ipass = endyr+1;
  get_fut_mortality();
  
  for (ipass=endyr+1;ipass<=endyr+Nproj-1;ipass++) 
	get_num_at_len_yr(); 

  Bzero =   mspbio_matetime(endyr+Nproj-2);  

  // Find F35%  
  Target = 0.35;
  IsB0 = 0;
  FutMort = 0.3;
  for (icnt=1;icnt<=20;icnt++)
   {
    ipass = endyr+1;
    get_fut_mortality();
   for (ipass=endyr+1;ipass<=endyr+Nproj-1;ipass++)
    get_num_at_len_yr(); 
 
    Btest = mspbio_matetime(endyr+Nproj-2) ;   

    Ratio = Btest/Bzero;
    cout << FutMort << " " << Ratio << endl;
    FutMort = FutMort * Ratio / Target;
   }
  F35 = FutMort; 
  SBPRF35 = Btest/FutRec;

// ========================y==================================================   
FUNCTION Find_OFL  
  dvariable Fmsy,Rbar,nn,alpha,beta;
  int BMSY_Yr1, BMSY_Yr2,ii,Iyr,kk,jj;
  
 //Define time period for BMSY  
  BMSY_Yr1 = 1979;BMSY_Yr2 = endyr-1;
  alpha = 0.1;
  beta = 0.25;
  Fmsy = F35;
  
 // Find Rbar 
  Rbar = 0; nn= 0;
  for (Iyr=BMSY_Yr1;Iyr<=BMSY_Yr2;Iyr++)
   {
    Rbar += mfexp(mean_log_rec(2)+rec_devf(Iyr));
    nn += 1;
   }
  Rbar = Rbar / nn;
  
 // Specify the BMSY proxy
  Bmsy = SBPRF35 * Rbar;

  //Begin projection
  // This code adjusts F so that the MMB at mate time next year is equal to the BMSY proxy
  // The resulting F is the FOFL
  ipass = endyr;
  if (ipass > endyr) 
	FutRec = Rbar;

  FutMort = Fmsy;
  get_fut_mortality();
  get_num_at_len_yr();
  
   if (mspbio_matetime(ipass) < Bmsy)
   {
     FutMort = 0;
     get_fut_mortality();
     get_num_at_len_yr();
     if (mspbio_matetime(ipass) > beta*Bmsy)
     {
       FutMort = Fmsy/2;
       get_fut_mortality();
       get_num_at_len_yr();
      for (ii=1;ii<=15;ii++)
       {
        FutMort = Fmsy*(mspbio_matetime(ipass)/Bmsy-alpha)/(1-alpha);
        get_fut_mortality();
        get_num_at_len_yr();
       }
      }
    }

   FOFL = FutMort;
   cout<<"Bzero"<<Bzero<<endl;     
   cout<<"SBPRF35"<<SBPRF35<<endl;   
   cout<<"FOFL"<<FOFL<<endl;
   
  //project under FOFL
   get_fut_mortality();
   get_num_at_len_yr();
   
   // find OFL from applying the FOFL
   get_catch_at_len();
   OFL = 0;
   OFL += pred_catch_ret(endyr);
   OFL += pred_catch_trawl(endyr);

   for(j=1;j<=2;j++)
    OFL += pred_catch_disc(j,endyr);

   cout<<"OFL"<<OFL<<endl; 
 
// ===========================================================================
FUNCTION Francis_weights  
  int w,x,y,z;

  Lbar.initialize();
  Lbar_hat_new.initialize();
  Lbar_hat_old.initialize();
  SE_Lbar_hat_new.initialize();
  SE_Lbar_hat_old.initialize();
  Francis_var_temp_new.initialize();
  Francis_var_temp_old.initialize();
  
  // obs(maturity, SC, sex, year,length), pred(maturity,sex, year,length)
  // calculate observed and predicted average lengths by year (z) 
  for(w=1;w<=2;w++)
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {
	  Lbar(w,x,y,z) 		=sum(elem_prod(obs_p_srv1_len(w,x,y,z)/sum(obs_p_srv1_len(w,x,y,z)),length_bins));
	  Lbar_hat_new(x,y,z) 	=sum(elem_prod(pred_p_srv1_len_new(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_new(x,y,yrs_srv1_length(z))),length_bins));
	  Lbar_hat_old(x,y,z) 	=sum(elem_prod(pred_p_srv1_len_old(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_old(x,y,yrs_srv1_length(z))),length_bins));	
     }

  // calculate the predicted standard error of the mean length of the catch by year (z)	 
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {	
      SE_Lbar_hat_new(x,y,z) = sqrt(sum(elem_prod((pred_p_srv1_len_new(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_new(x,y,yrs_srv1_length(z)))),pow((length_bins-Lbar_hat_new(x,y,z)),2))))/sqrt(nsamples_srv1_length(x,1,y,z)); 
      SE_Lbar_hat_old(x,y,z) = sqrt(sum(elem_prod((pred_p_srv1_len_old(x,y,yrs_srv1_length(z))/sum(pred_p_srv1_len_old(x,y,yrs_srv1_length(z)))),pow((length_bins-Lbar_hat_old(x,y,z)),2))))/sqrt(nsamples_srv1_length(x,2,y,z)); 
	  }
  // calculate the term for which the inverse variance is used to modify the input effective N for weighting len comps	
   for(x=1;x<=2;x++)
	for(y=1;y<=2;y++) 
     for(z=1;z<=nobs_srv1_length;z++)
	 {	
	  Francis_var_temp_new(x,y,z) = (Lbar(x,1,y,z)-Lbar_hat_new(x,y,z))/SE_Lbar_hat_new(x,y,z);
	  Francis_var_temp_old(x,y,z) = (Lbar(x,2,y,z)-Lbar_hat_old(x,y,z))/SE_Lbar_hat_old(x,y,z);
	 }

   // find the mena for females and males separately to calculate the variance
   // why isn't there a 'var()' function?! ugh.
   countFem = 0;
   countMal = 0;
   totalFem = 0;
   totalMal = 0;
  for(x=1;x<=2;x++)
   for(z=1;z<=nobs_srv1_length;z++)
   {
    // the if statements are because some of these are 'NaN' e.g. immature old shell crab
	// should find a better way to do this
	//females
	 if(Francis_var_temp_new(x,1,z)>0)
	 {
	  totalFem += Francis_var_temp_new(x,1,z);
	  countFem += 1;
	 }
	 if(Francis_var_temp_old(x,1,z)>0)
	 {
	  totalFem += Francis_var_temp_old(x,1,z);
	  countFem += 1;
	 }

	//males
	 if( Francis_var_temp_new(x,2,z)>0)
	 {
	  totalMal +=  Francis_var_temp_new(x,2,z);
	  countMal += 1;
	 }
 	 if( Francis_var_temp_old(x,2,z)>0)
	 {
	  totalMal +=  Francis_var_temp_old(x,2,z);
	  countMal += 1;
	 }
    }

	FemMeanVarTerm = 0.0;
	MaleMeanVarTerm = 0.0;
  // find the variance 
  for(x=1;x<=2;x++)
   for(z=1;z<=nobs_srv1_length;z++)
   {
	if(Francis_var_temp_new(x,1,z)>0)
     FemMeanVarTerm += square(totalFem/countFem - Francis_var_temp_new(x,1,z));
 	if(Francis_var_temp_old(x,1,z)>0)
	 FemMeanVarTerm += square(totalFem/countFem - Francis_var_temp_old(x,1,z));
 	if(Francis_var_temp_new(x,2,z)>0)
	 MaleMeanVarTerm += square(totalMal/countMal - Francis_var_temp_new(x,2,z));
 	if(Francis_var_temp_old(x,2,z)>0)
	 MaleMeanVarTerm += square(totalMal/countMal - Francis_var_temp_old(x,2,z));
   }
   
   Francis_weight_m = 1/(MaleMeanVarTerm/(countMal-1));
   Francis_weight_f = 1/(FemMeanVarTerm/(countFem-1));
// ========================y==================================================   
   
REPORT_SECTION
  Find_F35();
  Find_OFL();
  Francis_weights();
  
  report<<f<<" "<<Bmsy<<" "<< mspbio_srv1(endyr) << " " << mspbio_fishtime(endyr)<< " "<< mspbio_matetime(endyr-1)<< " " << F35 << " " << FOFL << " " << OFL << endl;
  report<<Fout<< " "<<endl;
  report<<"# mature male biomass at survey"<<endl;
  report<<mspbio_srv1<<endl;
  report<<"# Francis weight F"<<endl;
  report<<Francis_weight_f<<endl;
  report<<"# Francis weight M"<<endl;
  report<<Francis_weight_m<<endl;
  
  // this writes gradients for each parameter at each phase and at end "gradients.dat"
   save_gradients(gradients);

   //The report file holds all of the quantities used in Jack's projection program
   //The "Rout" below this writes to a file that is read for plotting in the safe document
  report<<"#number of length bins"<<endl;
  report<<nlenm<<endl;
  report<<"#Nat mort immature female/male"<<endl;
  report<<M<<endl;
  report<<"#nat mort mature new shell female/male"<<endl;
  report<<M_matn<<endl;
  report<<"#nat mort mature old shell female/male"<<endl;
  report<<M_mato<<endl;
  report<<"#constant recruitment"<<endl;
  report<<"1000000"<<endl;
  report<<"#average of last 4 years sel total male new old shell"<<endl;
  report<<(sel(1,endyr-4)+sel(1,endyr-3)+sel(1,endyr-2)+sel(1,endyr-1))/4.0<<endl;
  report<<(sel(1,endyr-4)+sel(2,endyr-3)+sel(2,endyr-2)+sel(2,endyr-1))/4.0<<endl;
  report<<"#average of last 4 years sel retained curve male new old shell"<<endl;
  report<<(sel_fit(1,endyr-3)+sel_fit(1,endyr-2)+sel_fit(1,endyr-1))/3.0<<endl;
  report<<(sel_fit(2,endyr-3)+sel_fit(2,endyr-2)+sel_fit(2,endyr-1))/3.0<<endl;
  report<<"#trawl selectivity female male"<<endl;
  report<<sel_trawl<<endl;
  report<<"#female pot discard selectivity"<<endl;
  report<<sel_discf(2009)<<endl;
  report<<"#maturity curve new shell female male"<<endl;
  report<<maturity_est(1)<<endl;
  report<<maturity_est(2)<<endl;
//  report<<maturity_average(1)<<endl;
//  report<<maturity_logistic<<endl;
  report<<"#maturity curve old shell female male"<<endl;
  report<<maturity_old_average<<endl;
  report<<"#molting probability immature female male"<<endl;
  report<<moltp<<endl;
  report<<"#molting probability mature female male"<<endl;
  report<<moltp_mat<<endl;
  report<<"#prop recruits to new shell"<<endl;
  report<<proprecn<<endl;
  report<<"#distribution of recruits to length bins"<<endl;
  report<<rec_len<<endl;
  report<<"#time of catch in fraction of year from survey - 7 months"<<endl;
  report<<catch_midpt(endyr)<<endl;
  report<<"#number at length new shell females males at time of fishery endyr from model"<<endl;
  report<<natl_new_fishtime(1,endyr)<<endl;
  report<<natl_new_fishtime(2,endyr)<<endl;
  report<<"#number at length old shell females males at time of fishery endyr from model"<<endl;
  report<<natl_old_fishtime(1,endyr)<<endl;
  report<<natl_old_fishtime(2,endyr)<<endl;
  report<<"#last year male spawning biomass"<<endl;
  report<<mspbio(endyr)<<endl;
  report<<"#last year female spawning biomass"<<endl;
  report<<fspbio(endyr)<<endl;
  report<<"#last year male spawning biomass at matingtime"<<endl;
  report<<mspbio_matetime(endyr-1)<<endl;
  report<<"#last year female spawning biomass at matingtime"<<endl;
  report<<fspbio_matetime(endyr-1)<<endl;
  report<<"#numbers at length immature new shell female male last year"<<endl;
  report<<natlength_inew(1,endyr)<<endl;
  report<<natlength_inew(2,endyr)<<endl;
  report<<"#numbers at length immature old shell female male last year"<<endl;
  report<<natlength_iold(1,endyr)<<endl;
  report<<natlength_iold(2,endyr)<<endl;
  report<<"#numbers at length mature new shell female male last year"<<endl;
  report<<natlength_mnew(1,endyr)<<endl;
  report<<natlength_mnew(2,endyr)<<endl;
  report<<"#numbers at length mature old shell female male last year"<<endl;
  report<<natlength_mold(1,endyr)<<endl;
  report<<natlength_mold(2,endyr)<<endl;
  report<<"#weight at length female juvenile"<<endl;
  report<<wtf(1)<<endl;
  report<<"#weight at length female mature"<<endl;
  report<<wtf(2)<<endl;
  report<<"#weight at length male"<<endl;
  report<<wtm<<endl;
  report<<"#length-length transition matrix"<<endl;
  report<<len_len<<endl;
  report<<"#female discard pot fishing F"<<endl;
  report<<fmortdf(endyr-1)<<endl;
  report<<"#trawl fishing F female male"<<endl;
  report<<fmortt(endyr-1)<<endl;
  report<<"#number of recruits from the model styr to endyr-1"<<endl;
  report<<endyr-styr<<endl;
  report <<"#recruitments female, male start year to endyr-1 from model" << endl;
  for(i=styr; i<endyr; i++)
  {
    report << mfexp(mean_log_rec(1)+rec_dev(1,i))<<" ";
  }
  report <<endl<< "#recruitments male, male start year+1 to endyr-1 from model" << endl;
  for(i=styr; i<endyr; i++)
  {
    report << mfexp(mean_log_rec(2)+rec_dev(2,i))<<" ";
  }
  report<<endl;
  report<<"#male spawning biomass at matetime for endyr-5 to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(endyr-5,endyr-1)<<endl;
  report<<"#male spawning biomass at matetime for str year to endyr-1 for spawner recruit curve to estimate recruitments"<<endl;
  report<<mspbio_matetime(styr,endyr-1)<<endl;
  report <<"#selectivity survey males 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  report << sel_srv3(2) << endl;
   
//Rout section 
    if (last_phase())
  {
//  R_out << "Estimated numbers of immature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
//
    R_out << "$Estimated numbers of immature new shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_inew(1,i) << endl;
       }
//  R_out << "Estimated numbers of immature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of immature old shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_iold(1,i) << endl;
       }
//  R_out << "Estimated numbers of mature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of mature new shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_mnew(1,i) << endl;
       }
//  R_out << "Estimated numbers of mature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of mature old shell female crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out <<  i<<" "<<natlength_mold(1,i) << endl;
       }

//  R_out << "Estimated numbers of immature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of immature new shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_inew(2,i) << endl;
      }
// R_out << "Estimated numbers of immature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of immature old shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_iold(2,i) << endl;
      }
// R_out << "Estimated numbers of mature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of mature new shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_mnew(2,i) << endl;
      }
 //  R_out << "Estimated numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Estimated numbers of mature old shell male crab" << endl;
     for(i=styr;i<=endyr;i++)
      {
        R_out << i<<" "<<natlength_mold(2,i) << endl;
      }
 //  R_out << "Observed numbers of immature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature new shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,1,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature new shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature new shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,1,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature old shell female crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature old shell female crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,2,1,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of immature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature new shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,1,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of immature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of immature old shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(1,2,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature new shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature new shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,1,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }
 //  R_out << "Observed numbers of mature old shell male crab by length: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed numbers of mature old shell male crab" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
           R_out<<yrs_srv1_length(i)<<" "<<obs_p_srv1_len(2,2,2,i)*obs_srv1t(yrs_srv1_length(i))<<endl;
      }

  //  R_out << "Observed Survey Numbers by length females:  'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed survey numbers female" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
        R_out<<yrs_srv1_length(i)<<" " << obs_srv1_num(1,yrs_srv1_length(i)) << endl;
      }

  //  R_out << "Observed Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed survey numbers male" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
        R_out<<yrs_srv1_length(i)<<" " << obs_srv1_num(2,yrs_srv1_length(i))<< endl;
      }
  //  R_out << "Predicted Survey Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted survey numbers female" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
       R_out<<yrs_srv1_length(i)<<" "  << pred_srv1(1,yrs_srv1_length(i)) << endl;
      }
  //  R_out << "Predicted Survey Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted survey numbers male" << endl;
      for (i=1; i <= nobs_srv1_length; i++)
      {
         R_out<<yrs_srv1_length(i)<<" "  << pred_srv1(2,yrs_srv1_length(i)) << endl;
       }
  //  R_out << "Predicted pop Numbers by length females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted population numbers female" << endl;
     for(i=styr;i<=endyr;i++)
      {
       R_out<<i<<" "<< natlength(1,i)<< endl;
      }
  //  R_out << "Predicted pop Numbers by length males: 'year', '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted population numbers male" << endl;
     for(i=styr;i<=endyr;i++)
      {
         R_out<<i<<" "<< natlength(2,i)<< endl;
       }

  //  R_out<<"observed number of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed number males greater than 101mm" << endl;
  R_out<<obs_lmales<<endl;
  //  R_out<<"observed biomass of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed biomass males greater than 101mm" << endl;
  R_out<<obs_lmales_bio<<endl;
  //  R_out<<"Predicted number of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted number males greater than 101mm" << endl;
  R_out<<num_males_gt101<<endl;
  //  R_out<<"Predicted biomass of males greater than 101 mm: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted biomass males greater than 101mm" << endl;
  R_out<<bio_males_gt101<<endl;
  
  //  R_out<<"pop estimate numbers of males >101: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Population numbers male" << endl;
        R_out<<legal_males<<endl;
      //  R_out<<"estimated population biomass of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Population biomass male" << endl;
      R_out<<legal_males_bio<<endl;
      //  R_out<<"estimated survey numbers of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Estimated survey numbers male" << endl;
      R_out<<legal_srv_males<<endl;
      //  R_out<<"estimated survey biomass of males > 101: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Estimated survey biomass male" << endl;
      R_out<<legal_srv_males_bio<<endl;
  //  R_out << "Observed survey biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey biomass" << endl;
  R_out << obs_srv1_biom(styr,endyr)<<endl;
  //  R_out << "predicted survey biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey biomass" << endl;
  R_out << pred_srv1_bioms(1)+pred_srv1_bioms(2)<<endl;
  //survey numbers
    for(k=1;k<=2;k++)
    {
     for(i=styr;i<=endyr;i++)
      {
       tmpo(k,i)=sum(obs_srv1_num(k,i));
       tmpp(k,i)=sum(pred_srv1(k,i));
      }
     }
  //  R_out << "Observed survey numbers female: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey numbers female" << endl;
  R_out << tmpo(1)<<endl;
  //  R_out << "Observed survey numbers male: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey numbers male" << endl;
  R_out << tmpo(2)<<endl;
  //  R_out << "predicted survey numbers female: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey numbers female" << endl;
  R_out << tmpp(1)<<endl;
  //  R_out << "predicted survey numbers male: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Predicted survey numbers male" << endl;
  R_out << tmpp(2)<<endl;
  //  R_out << "Observed survey female spawning biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female spawning biomass" << endl;
  R_out << obs_srv1_spbiom(1)<<endl;
  //  R_out << "Observed survey male spawning biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male spawning biomass" << endl;
  R_out << obs_srv1_spbiom(2)<<endl;
  //  R_out << "Observed survey female new spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$bserved survey female new spawning numbers" << endl;
  R_out << obs_srv1_spnum(1,1)<<endl;
  //  R_out << "Observed survey female old spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female old spawning numbers" << endl;
  R_out << obs_srv1_spnum(2,1)<<endl;
  //  R_out << "Observed survey male new spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male new spawning numbers" << endl;
  R_out << obs_srv1_spnum(1,2)<<endl;
  //  R_out << "Observed survey male old spawning numbers: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male old spawning numbers" << endl;
  R_out << obs_srv1_spnum(2,2)<<endl;
  //  R_out << "Observed survey female biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey female biomass" << endl;
  R_out << obs_srv1_bioms(1)<<endl;
  //  R_out << "Observed survey male biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$Observed survey male biomass" << endl;
  R_out << obs_srv1_bioms(2)<<endl;
  //  R_out << "natural mortality immature females, males: 'FemM','MaleM'" << endl;
  R_out << "$natural mortality immature" << endl;
  R_out << M << endl;
  //  R_out << "natural mortality mature females, males: 'FemMm','MaleMm'" << endl;
  R_out << "$atural mortality mature" << endl;
  R_out << M_matn << endl;
  //  R_out << "natural mortality mature old shell females, males: 'FemMmo','MaleMmo'" << endl;
  R_out << "$natural mortality mature old shell" << endl;
  R_out << M_mato << endl;
  //  R_out << "Predicted Biomass: seq(1978,"<<endyr<<")" << endl;
  R_out << "$Predicted Biomass" << endl;
  R_out << pred_bio << endl;
  //  R_out << "Predicted total population numbers: seq(1978,"<<endyr<<") "<<endl;
  R_out << "$Predicted total population numbers" << endl;
  R_out <<popn<<endl;
//  //  R_out << "Predicted Exploitable Biomass: seq(1978,"<<endyr<<") " << endl;
//  R_out << explbiom << endl;
//  R_out << "Female Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Female Spawning Biomass" << endl;
  R_out << fspbio << endl;
  //  R_out << "Male Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature Male Biomass" << endl;
  R_out << mspbio << endl;
  //  R_out << "Total Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total Spawning Biomass" << endl;
  R_out << fspbio+mspbio << endl;
  //  R_out << "Female Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Female Spawning Biomass at fish time" << endl;
  R_out << fspbio_fishtime << endl;
  //  R_out << "Male Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature Male Biomass at fish time" << endl;
  R_out << mspbio_fishtime << endl;
  //  R_out << "Total Spawning Biomass at fish time: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total Spawning Biomass at fish time" << endl;
  R_out << fspbio_fishtime+mspbio_fishtime << endl;
  //  R_out << "Mating time Female Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Female Spawning Biomass" << endl;
  R_out << fspbio_matetime << endl;
  //  R_out << "Mating time Male Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mature male biomass at mating" << endl;
  R_out << mspbio_matetime << endl;
  //  R_out << "Mating time Male old shell Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Old shell mature male biomass at mating" << endl;
  R_out << mspbio_old_matetime << endl;
  //  R_out << "Mating time female new shell Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$New shell female spawning biomass at mating time" << endl;
  R_out << fspbio_new_matetime << endl;
  //  R_out << "Mating time Total Spawning Biomass : seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total mature biomass at mating time" << endl;
  R_out << fspbio_matetime+mspbio_matetime << endl;
  //  R_out << "Mating time effective Female Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Effective female spawning biomass at mating time" << endl;
  R_out << efspbio_matetime << endl;
  //  R_out << "Mating time effective Male Spawning Biomass(old shell only): seq(1978,"<<endyr<<") " << endl;
  R_out << "$Effective mature male biomass at mating time" << endl;
  R_out << emspbio_matetime << endl;
  //  R_out << "Mating time Total effective Spawning Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Total effective spawning biomass at mating time" << endl;
  R_out << efspbio_matetime+emspbio_matetime << endl;
  //  R_out << "Mating time male Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Spawning numbers male at mating time" << endl;
  R_out << mspnum_matetime << endl;
  //  R_out << "Mating time Female Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Spawning number female at mating time" << endl;
  R_out << efspnum_matetime << endl;
  //  R_out << "Mating time Male Spawning numbers(old shell only): seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Male Spawning numbers old " << endl;
  R_out << emspnum_old_matetime << endl;
  //  R_out << "ratio Mating time Female Spawning numbers to male old shell mature numbers : seq(1978,"<<endyr<<") " << endl;
  R_out << "$ratio Mating time Female Spawning numbers to male old shell mature numbers" << endl;
  R_out << elem_div(efspnum_matetime,emspnum_old_matetime) << endl;
  //  R_out << "Mating time effective Female new shell Spawning biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time effective Female new shell Spawning biomass" << endl;
  R_out <<efspbio_new_matetime << endl;
  //  R_out << "Mating time Female new shell Spawning numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Mating time Female new shell Spawning numbers" << endl;
  R_out << fspnum_new_matetime << endl;
  //  R_out << "ratio Mating time Female new shell Spawning numbers to male old shell mature numbers : seq(1978,"<<endyr<<") " << endl;
  R_out << "$ratio Mating time Female new shell Spawning numbers to male old shell mature numbers" << endl;
  R_out << elem_div(fspnum_new_matetime,emspnum_old_matetime) << endl;
  //  R_out << "Predicted Female survey Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey Biomass" << endl;
  R_out << pred_srv1_bioms(1) << endl;
  //  R_out << "Predicted Male survey Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey Biomass" << endl;
  R_out << pred_srv1_bioms(2)<< endl;
  //  R_out << "Predicted Female survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey mature Biomass" << endl;
  R_out << fspbio_srv1 << endl;
  //  R_out << "Predicted Male survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey mature Biomass" << endl;
  R_out << mspbio_srv1<< endl;
  //  R_out << "Predicted total survey mature Biomass: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted total survey mature Biomass" << endl;
  R_out << fspbio_srv1+mspbio_srv1<< endl;
  //  R_out << "Predicted Female survey new mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Female survey new mature numbers" << endl;
  R_out << fspbio_srv1_num(1) << endl;
  //  R_out << "Predicted Female survey old mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$redicted Female survey old mature numbers" << endl;
  R_out << fspbio_srv1_num(2) << endl;
  //  R_out << "Predicted Male survey new mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey new mature numbers" << endl;
  R_out << mspbio_srv1_num(1)<< endl;
  //  R_out << "Predicted Male survey old mature numbers: seq(1978,"<<endyr<<") " << endl;
  R_out << "$Predicted Male survey old mature numbers" << endl;
  R_out << mspbio_srv1_num(2)<< endl;
//2009 bsfrf study
  //  R_out << "Observed industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Observed industry survey mature biomass" << endl;
  R_out << obs_srv2_spbiom(1,1)<<" "<<obs_srv2_spbiom(1,2)<<" "<<obs_srv2_spbiom(2,1)<<" "<<obs_srv2_spbiom(2,2)<<endl;
  //  R_out << "Predicted industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Predicted industry survey mature biomass" << endl;
  R_out << fspbio_srv2_ind<<" "<<mspbio_srv2_ind<<" "<<fspbio_srv2_nmfs<<" "<<mspbio_srv2_nmfs<<endl;
  //  R_out << "Observed Prop industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$bserved Prop industry survey female" << endl;
  R_out <<(obs_p_srv2_len(1,1,1)+obs_p_srv2_len(1,1,2))<<endl;
  //  R_out << "Observed Prop industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry survey male" << endl;
  R_out <<(obs_p_srv2_len(1,2,1)+obs_p_srv2_len(1,2,2))<<endl;
  //  R_out << "Observed Prop industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry nmfs survey female" << endl;
  R_out <<obs_p_srv2_len(2,1,1)+obs_p_srv2_len(2,1,2)<<endl;
  //  R_out << "Observed Prop industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop industry nmfs survey male" << endl;
  R_out <<obs_p_srv2_len(2,2,1)+obs_p_srv2_len(2,2,2)<<endl;
  //  R_out << "Predicted Prop industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry survey female" << endl;
  R_out <<pred_p_srv2_len_ind(1)<<endl;
  //  R_out << "Predicted Prop industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry survey male" << endl;
  R_out <<pred_p_srv2_len_ind(2)<<endl;
  //  R_out << "Predicted Prop industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop industry nmfs survey female" << endl;
  R_out <<pred_p_srv2_len_nmfs(1)<<endl;
  //  R_out << "Predicted Prop industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$redicted Prop industry nmfs survey male" << endl;
  R_out <<pred_p_srv2_len_nmfs(2)<<endl;
//2010 bsfrf study   
  //  R_out << "Observed 2010 industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Observed 2010 industry survey mature biomass" << endl;
  R_out << obs_srv10_spbiom(1,1)<<" "<<obs_srv10_spbiom(1,2)<<" "<<obs_srv10_spbiom(2,1)<<" "<<obs_srv10_spbiom(2,2)<<endl;
  //  R_out << "Predicted 2010 industry survey mature biomass: seq(1,4) " << endl;
  R_out << "$Predicted 2010 industry survey mature biomass" << endl;
  R_out << fspbio_srv10_ind<<" "<<mspbio_srv10_ind<<" "<<fspbio_srv10_nmfs<<" "<<mspbio_srv10_nmfs<<endl;
  //  R_out << "Observed Prop 2010 industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry survey female" << endl;
  R_out <<(obs_p_srv10_len(1,1,1)+obs_p_srv10_len(1,1,2))<<endl;
  //  R_out << "Observed Prop 2010 industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry survey male" << endl;
  R_out <<(obs_p_srv10_len(1,2,1)+obs_p_srv10_len(1,2,2))<<endl;
  //  R_out << "Observed Prop 2010 industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry nmfs survey female" << endl;
  R_out <<obs_p_srv10_len(2,1,1)+obs_p_srv10_len(2,1,2)<<endl;
  //  R_out << "Observed Prop 2010 industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed Prop 2010 industry nmfs survey male" << endl;
  R_out <<obs_p_srv10_len(2,2,1)+obs_p_srv10_len(2,2,2)<<endl;
  //  R_out << "Predicted Prop 2010 industry survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry survey female" << endl;
  R_out <<pred_p_srv10_len_ind(1)<<endl;
  //  R_out << "Predicted Prop 2010 industry survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry survey male" << endl;
  R_out <<pred_p_srv10_len_ind(2)<<endl;
  //  R_out << "Predicted Prop 2010 industry nmfs survey female:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry nmfs survey female" << endl;
  R_out <<pred_p_srv10_len_nmfs(1)<<endl;
  //  R_out << "Predicted Prop 2010 industry nmfs survey male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted Prop 2010 industry nmfs survey male" << endl;
  R_out <<pred_p_srv10_len_nmfs(2)<<endl;



   //  R_out << "Observed Prop fishery ret new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery retained new male" << endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i) << " " << obs_p_fish_ret(1,i)<< endl;
         }
    //  R_out << "Predicted length prop fishery ret new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery retained new male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish_fit(1,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery ret old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery retained old male" << endl;
      for (i=1; i<=nobs_fish; i++)
        {
          R_out << yrs_fish(i) << " " << obs_p_fish_ret(2,i)<< endl;
         }
    //  R_out << "Predicted length prop fishery ret old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery retained old male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish; i++)
        {
          ii=yrs_fish(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish_fit(2,ii)  << endl;
         }

  //  R_out << "Observed Prop fishery total new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery total new male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_tot(1,i) << endl;
         }
    //  R_out << "Predicted length prop fishery total new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery total new male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish_discm; i++)
        {
          ii=yrs_fish_discm(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish(1,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery total old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery total old male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_tot(2,i) << endl;
         }
    //  R_out << "Predicted length prop fishery total old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery total old male" << endl;
   //R_out<<pred_p_fish<<endl;
       for (i=1; i<=nobs_fish_discm; i++)
        {
          ii=yrs_fish_discm(i);  
          R_out <<  ii  <<  " "  <<  pred_p_fish(2,ii)  << endl;
         }
  //  R_out << "Observed Prop fishery discard new males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery discard new male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_discm(1,i) << endl;
         }

 // R_out << "Observed Prop fishery discard old males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Observed proportion fishery discard old male" << endl;
      for (i=1; i<=nobs_fish_discm; i++)
        {
          R_out << yrs_fish_discm(i) << " " << obs_p_fish_discm(2,i)<< endl;
         }

    //  R_out << "Observed length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion fishery discard all female" << endl;
    for (i=1; i<=nobs_fish_discf; i++)
    {
      R_out <<  yrs_fish_discf(i)  <<  " "  <<  obs_p_fish_discf(i)  << endl;
    }
    //  R_out << "Predicted length prop fishery discard all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion fishery discard all female" << endl;
    for (i=1; i<=nobs_fish_discf; i++)
    {
      ii=yrs_fish_discf(i);  
      R_out <<  ii  <<  " "  <<  pred_p_fish_discf(ii)  << endl;
    }

    //  R_out << "Predicted length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion trawl female" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {
          ii=yrs_trawl(i);  
          R_out <<  ii  <<  " "  <<  pred_p_trawl(1,ii)  << endl;
         }
    //  R_out << "Observed length prop trawl females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion trawl female" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {  
           R_out <<  yrs_trawl(i)  <<  " "  <<  obs_p_trawl(1,i)  << endl;
         }
    //  R_out << "Predicted length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion trawl male" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {
          ii=yrs_trawl(i);  
          R_out <<  ii  <<  " "  <<  pred_p_trawl(2,ii)  << endl;
         }
    //  R_out << "Observed length prop trawl males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion trawl male" << endl;
       for (i=1; i<=nobs_trawl; i++)
        {  
           R_out <<  yrs_trawl(i)  <<  " "  <<  obs_p_trawl(2,i)  << endl;
         }

//  R_out << "Observed length Prop fishery males:'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  //    for (i=1; i<=nobs_fish; i++)
  //      {
  //        R_out << yrs_fish(i) << " " << obs_p_fish(2,i) << endl;
  //       }
//    R_out << "Predicted length prop fishery males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
   //R_out<<pred_p_fish<<endl;
   //    for (i=1; i<=nobs_fish; i++)
  //      {
    //      ii=yrs_fish(i);  
    //      R_out <<  ii  <<  " "  <<  pred_p_fish(2,ii)  << endl;
    //     }
 //  R_out << "Observed Length Prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature new female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey immature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature new female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(1,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature old female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,2,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey immature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature old female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_old(1,1,ii) << endl;
          }
 
  //  R_out << "Observed Length Prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature new male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(1,1,2,i) << endl;
          }
  //  R_out << "Predicted length prop survey immature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature new male" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(1,2,ii) << endl;
         }
 //  R_out << "Observed Length Prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey immature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(1,2,2,i) << endl;
 }
 //  R_out << "Predicted length prop survey immature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey immature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_old(1,2,ii) << endl;
 }
 //  R_out << "Observed Length Prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature new female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(2,1,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey mature new females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature new female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(2,1,ii) << endl;
          }
  //  R_out << "Observed Length Prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature old female" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(2,2,1,i) << endl;
           }
  //  R_out << "Predicted length prop survey mature old females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature old female" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_old(2,1,ii) << endl;
          }
 
  //  R_out << "Observed Length Prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature new male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(2,1,2,i) << endl;
          }
  //  R_out << "Predicted length prop survey mature new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature new male" << endl;
       for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);  
             R_out << ii << " " << pred_p_srv1_len_new(2,2,ii) << endl;
         }
 //  R_out << "Observed Length Prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey mature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   R_out << yrs_srv1_length(i) <<" " <<obs_p_srv1_len(2,2,2,i) << endl;
 }
 //  R_out << "Predicted length prop survey mature old males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey mature old male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_old(2,2,ii) << endl;
 }
 //  R_out << "Observed Length Prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey all females" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,1,i)+obs_p_srv1_len(2,1,1,i)+obs_p_srv1_len(2,2,1,i)<< endl;
              tmpp4+=obs_p_srv1_len(1,1,1,i)+obs_p_srv1_len(2,1,1,i)+obs_p_srv1_len(2,2,1,i);
                         }
 //  R_out << "Predicted length prop survey all females: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey all female" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_new(1,1,ii)+pred_p_srv1_len_new(2,1,ii)+pred_p_srv1_len_old(2,1,ii) << endl;
    tmpp1+=pred_p_srv1_len_new(1,1,ii)+pred_p_srv1_len_new(2,1,ii)+pred_p_srv1_len_old(2,1,ii);
    }
 //  R_out << "Observed Length Prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Observed proportion survey all male" << endl;
         for (i=1; i<=nobs_srv1_length; i++)
          {
             ii=yrs_srv1_length(i);
              R_out << ii <<" " <<obs_p_srv1_len(1,1,2,i)+obs_p_srv1_len(1,2,2,i)+obs_p_srv1_len(2,1,2,i)+obs_p_srv1_len(2,2,2,i)<< endl;
         tmpp2+=obs_p_srv1_len(1,1,2,i)+obs_p_srv1_len(1,2,2,i)+obs_p_srv1_len(2,1,2,i)+obs_p_srv1_len(2,2,2,i);
                        }
 //  R_out << "Predicted length prop survey all males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Predicted proportion survey all male" << endl;
 for (i=1; i<=nobs_srv1_length; i++)
 {
   ii=yrs_srv1_length(i);  
   R_out << ii << " " << pred_p_srv1_len_new(1,2,ii)+pred_p_srv1_len_new(2,2,ii)+pred_p_srv1_len_old(2,2,ii) << endl;
  tmpp3+=pred_p_srv1_len_new(1,2,ii)+pred_p_srv1_len_new(2,2,ii)+pred_p_srv1_len_old(2,2,ii);
    }
  //  R_out << "Sum of predicted prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of predicted proportion survey all female" << endl;
            R_out <<tmpp1<<endl;
  //  R_out << "Sum of predicted prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of predicted proportion survey all male" << endl;
            R_out <<tmpp3<<endl;
  //  R_out << "Sum of Observed prop survey all females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of observed proportion survey all female" << endl;
            R_out <<tmpp4<<endl;
  //  R_out << "Sum of Observed prop survey all males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
  R_out << "$Sum of observed proportion survey all male" << endl;
            R_out <<tmpp2<<endl;
                

// R_out << "Observed Length Prop fishery new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
//         for (i=1; i<=nobs_fish_length; i++)
//          {
//             R_out << yrs_fish_length(i) <<" " <<obs_p_fish_len(1,i) << endl;
//          }
//  R_out << "Predicted length prop fishery new males: 'year','27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'" << endl;
//       for (i=1; i<=nobs_fish_length; i++)
//          {
//             ii=yrs_fish_length(i);  
//             R_out << ii << " " << pred_p_fish_discm(1,ii) << endl;
//         }
//          
//  R_out << "Observed mean length at age females:  'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
//        for(i=1;i<=nyrs_mlen; i++)
//  {
//                   ii=yrs_mlen(i);
//                 R_out << ii << " " << mean_length_obs(1,i) << endl;
//  }
//  R_out << "Observed mean length at age males:  'year','3','4','5','6','7','8','9','10','11','12','13','14','15'" << endl;
//        for(i=1;i<=nyrs_mlen; i++)
//  {
//                  ii=yrs_mlen(i); 
//                R_out << ii << " " << mean_length_obs(2,i) << endl;
//  }
  //  R_out << "Predicted mean postmolt length females:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted mean postmolt length female" << endl;
  R_out << mean_length(1) << endl;
  //  R_out << "Predicted mean postmolt length males:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted mean postmolt length male" << endl;
  R_out << mean_length(2)<<endl; 
  //  R_out<< "sd mean length females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$sd mean length female" << endl;
  R_out<<sd_mean_length(1)<<endl;
  //  R_out<< "sd mean length males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$sd mean length male" << endl;
  R_out<<sd_mean_length(2)<<endl;
  //  R_out << "af: 'females'" << endl;
  R_out << "$af " << endl;
  R_out << af << endl;
  //  R_out << "am: 'males'" << endl;
  R_out << "$am" << endl;
  R_out << am << endl;
  R_out << "$af2" << endl;
  R_out << af +(bf-bf1)*deltaf << endl;
  R_out << "$am2" << endl;
  R_out << am+(bm-b1)*deltam << endl;
  //  R_out << "bf: 'females'" << endl;
  R_out << "$bf " << endl;
  R_out << bf << endl;
  R_out << "$bf1" << endl;
  R_out << bf1 << endl;
  R_out << "$EN4" << endl;
    R_out << deltaf << endl;
  R_out << "$bm " << endl;
  R_out << bm << endl;
  R_out << "$b1" << endl;
  R_out << b1 << endl;
  R_out << "$delta male" << endl;
    R_out << deltam << endl;
  R_out << "$st_gr" << endl;
  R_out << st_gr << endl;
  R_out << "$st_gr" << endl;
  R_out << st_gr << endl;
    //  R_out << "Predicted probability of maturing females:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted probability of maturing female" << endl;
  R_out << maturity_est(1)<<endl; 
  //  R_out << "Predicted probability of maturing males:  '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Predicted probability of maturing male" << endl;
  R_out << maturity_est(2)<<endl; 
  //  R_out<<"molting probs female: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$Molting probability female" << endl;
  R_out<<moltp(1)<<endl;
  //  R_out<<"molting probs male:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Molting probability male" << endl;
  R_out<<moltp(2)<<endl;
  //  R_out <<"Molting probability mature males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$Molting probability mature males" << endl;
  R_out <<moltp_mat(2)<<endl;
  //  R_out << "observed pot fishery cpue 1979 fishery to endyr fishery: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed pot fishery cpue" << endl;
  R_out <<cpue(1979,endyr)<<endl;
  //  R_out << "predicted pot fishery cpue 1978 to endyr-1 survey: seq(1978,"<<endyr-1<<")" << endl;
  R_out << "$predicted pot fishery cpue" << endl;
  R_out <<cpue_pred(1978,endyr-1)<<endl;
  //  R_out << "observed retained catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$Observed retained catch biomass" << endl;
  R_out << catch_ret(styr,endyr-1) << endl;
  //  R_out << "predicted retained catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$Predicted retained catch biomas" << endl;
  R_out << pred_catch_ret(styr,endyr-1)<<endl;
  //  R_out << "predicted retained new catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained new catch biomass" << endl;
  R_out << (catch_male_ret_new*wtm)(styr,endyr-1)<<endl;
  //  R_out << "predicted retained old catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained old catch biomass" << endl;
  R_out << (catch_male_ret_old*wtm)(styr,endyr-1)<<endl;
  //  R_out << "observed retained+discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed retained plus discard male catch biomass" << endl;
  R_out << obs_catchtot_biom(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard male catch biomass" << endl;
  R_out << pred_catch(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard new male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard new male catch biomass" << endl;
  R_out << (catch_lmale_new*wtm)(styr,endyr-1) << endl;
  //  R_out << "predicted retained+discard old male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted retained plus discard old male catch biomass" << endl;
  R_out << (catch_lmale_old*wtm)(styr,endyr-1) << endl;
  //  R_out << "observed discard male mortality biomass: seq(1979,"<<endyr<<")"<<endl;
  // R_out << "$FD" << endl;
  // R_out << (obs_catchtot_biom-catch_ret)(styr,endyr-1) <<endl;
  //  R_out << "predicted discard male catch biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted discard male catch biomass" << endl;
  R_out << pred_catch(styr,endyr-1) -pred_catch_ret(styr,endyr-1)<< endl;
  //  R_out << "observed female discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed female discard mortality biomass" << endl;
  R_out << obs_catchdf_biom(styr,endyr-1) << endl;
  //  R_out << "predicted female discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$predicted female discard mortality biomass" << endl;
  R_out << pred_catch_disc(1)(styr,endyr-1) << endl;
  //  R_out << "observed male discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
  R_out << "$observed male discard mortality biomass" << endl;
  R_out << obs_catchdm_biom(styr,endyr-1) << endl;
//  R_out << "predicted male discard mortality biomass: seq(1979,"<<endyr<<")" << endl;
//  R_out << pred_catch_disc(2)(styr,endyr-1) << endl;
  //  R_out << "observed trawl catch biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$observed trawl catch biomass" << endl;
  R_out << obs_catcht_biom<<endl;
  //  R_out << "predicted trawl catch biomass: seq(1978,"<<endyr<<")"<<endl;
  R_out << "$predicted trawl catch biomass" << endl;
  R_out <<pred_catch_trawl<<endl;
  //  R_out << "estimated retained catch div. by male spawning biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated retained catch div by male spawning biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch_ret,mspbio_fishtime)(styr,endyr-1) << endl;
  //  R_out << "estimated total catch div. by male spawning biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch div by male spawning biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch,mspbio_fishtime)(styr,endyr-1) << endl;
  //  R_out << "estimated total catch of males >101 div. by males >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch of males 101 div by males 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_gt101(styr,endyr-1),bio_males_gt101(styr,endyr-1)) << endl;
  //  R_out << "estimated total catch numbers of males >101 div. by males numbers >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch numbers of males 101 div by males numbers 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_no_gt101(styr,endyr-1),num_males_gt101(styr,endyr-1)) << endl;
  //  R_out << "estimated total catch numbers of males >101 div. by survey estimate males numbers >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch numbers of males 101 div by survey estimate males numbers 101 at fishtime" << endl;
     for(i=styr;i<endyr;i++)
        {
         obs_tmp(i) = obs_lmales(i-(styr-1));
        }
  R_out <<elem_div(pred_catch_no_gt101(styr,endyr-1),obs_tmp(styr,endyr-1)*mfexp(-M_matn(2)*(7/12))) << endl;
     for(i=styr;i<endyr;i++)
        {
         obs_tmp(i) = obs_lmales_bio(i-(styr-1));
        }

  //  R_out << "estimated total catch biomass of males >101 div. by survey estimate male biomass >101 at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch biomass of males 101 div by survey estimate male biomass 101 at fishtime" << endl;
  R_out <<elem_div(pred_catch_gt101(styr,endyr-1),obs_tmp(styr,endyr-1)*mfexp(-M_matn(2)*(7/12)) ) << endl;

  //  R_out << "estimated total catch biomass div. by survey estimate male mature biomass at fishtime: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated total catch biomass div. by survey estimate male mature biomass at fishtime" << endl;
  R_out <<elem_div(pred_catch(styr,endyr-1),((obs_srv1_spbiom(2))(styr,endyr-1))*mfexp(-M_matn(2)*(7/12))) << endl;
  //  R_out << "estimated annual total fishing mortality: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual total fishing mortality" << endl;
  R_out << mfexp(log_avg_fmort+fmort_dev)(styr,endyr-1) << endl;
//  R_out << "estimated target annual fishing mortality rate: seq(1978,"<<endyr<<")" << endl;
//  R_out <<ftarget<<endl;
  //  R_out <<"retained F: seq(1978,"<<endyr<<")" << endl;
  R_out << "$retained F" << endl;
         for(i=styr;i<=endyr;i++){
          R_out <<F_ret(1,i)(22)<<" ";
         }
  R_out<<endl;
//  R_out <<"predicted ghl: seq(1978,"<<endyr<<")" << endl;
//  R_out <<pred_catch_target<<endl;
  //  R_out <<"ghl: seq(1978,"<<endyr<<")" << endl;
  R_out << "$ghl" << endl;
  R_out <<catch_ghl/2.2<<endl;
    //  R_out << "estimated annual fishing mortality females pot: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual fishing mortality females pot" << endl;
  R_out << fmortdf(styr,endyr-1) <<endl;
  //  R_out << "estimated annual fishing mortality trawl bycatch: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated annual fishing mortality trawl bycatch" << endl;
  R_out << fmortt(styr,endyr-1) <<endl;
//  R_out << "initial number of recruitments female: seq(1963,1977)" << endl;
//  for(i=styr-nirec; i<=styr-1; i++)
//  {
//    R_out << mfexp(mean_log_rec(1)+rec_dev(1,i))<<" ";
//  }
//  R_out <<endl<< "initial number of recruitments male: seq(1963,1977)" << endl;
//  for(i=styr-nirec; i<=styr-1; i++)
//  {
//    R_out << mfexp(mean_log_rec(2)+rec_dev(2,i))<<" ";
//  } 
//recruits in the model are 1978 to 2004, the 1978 recruits are those that enter the population
//in spring of 1979, before the 1979 survey - since using the survey as the start of the year
// in the model spring 1979 is stil 1978.  the last recruits are 2003 that come in spring 2004
//  R_out << "estimated number of recruitments female: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated number of recruits female" << endl;
  for(i=styr; i<=endyr-1; i++)
  {
    R_out << mfexp(mean_log_rec(1)+rec_dev(1,i))<<" ";
  }
  R_out<<endl;
  //  R_out <<endl<< "estimated number of recruitments male: seq(1979,"<<endyr<<")" << endl;
  R_out << "$estimated number of recruits male" << endl;
  for(i=styr; i<=endyr-1; i++)
  {
    R_out << mfexp(mean_log_rec(2)+rec_dev(2,i))<<" ";
  }
  for(i=1;i<=median_rec_yrs;i++)R_out<<2*median_rec<<" ";
    R_out<<endl;
 // R_out<<"distribution of recruits to length bins: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<<endl;
  R_out << "$distribution of recruits to length bins" << endl;
  R_out<<rec_len<<endl;
 // R_out<<"fishery total selectivity new shell 50% parameter: seq(1979,"<<endyr<<")"<<endl;
  R_out << "$fishery total selectivity new shell 50 parameter" << endl;
  R_out <<mfexp(log_avg_sel50_mn+log_sel50_dev_mn)(styr,endyr-1)<<endl;
 // R_out <<"fishery total selectivity old shell 50% parameter: seq(1979,"<<endyr<<")"<<endl;
  R_out << "$fishery total selectivity old shell 50 parameter" << endl;
  R_out <<mfexp(log_avg_sel50_mo+log_sel50_dev_mo)(styr,endyr-1)<<endl;
//  R_out << "selectivity fishery total new males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery total new male" << endl;
  R_out << sel(1) << endl;
 // R_out << "selectivity fishery total old males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery total old male" << endl;
  R_out << sel(2) << endl;
 // R_out << "selectivity fishery ret new males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery retained new male" << endl;
  R_out << sel_fit(1) << endl;
 // R_out << "selectivity fishery ret old males: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity fishery retained old male" << endl;
  R_out << sel_fit(2) << endl;
 // R_out <<"retention curve males new: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$retention curve new male" << endl;
  R_out <<sel_ret(1,endyr-1)<<endl;
 // R_out <<"retention curve males old: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$retention curve old male" << endl;
  R_out <<sel_ret(2,endyr-1)<<endl;
 // R_out << "selectivity discard females: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity discard female" << endl;
  R_out <<sel_discf<<endl;
 // R_out << "selectivity trawl females:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity trawl female" << endl;
  R_out <<sel_trawl(1)<<endl;
 // R_out << "selectivity trawl males:'27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity trawl male" << endl;
  R_out <<sel_trawl(2)<<endl;
 // R_out << "selectivity survey females 1978 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 1" << endl;
  R_out << sel_srv1(1) << endl;
 // R_out << "selectivity survey males 1978 1981: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 1" << endl;
  R_out << sel_srv1(2) << endl;
 // R_out << "selectivity survey females 1982 to 1988: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 2" << endl;
  R_out << sel_srv2(1) << endl;
 // R_out << "selectivity survey males 1982 to 1988: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 2" << endl;
  R_out << sel_srv2(2) << endl;
 // R_out << "selectivity survey females 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey female Era 3" << endl;
  R_out << sel_srv3(1) << endl;
 // R_out << "selectivity survey males 1989 to endyr: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity survey male Era 3" << endl;
  R_out << sel_srv3(2) << endl;
 // R_out << "selectivity industry survey females 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey females 2009" << endl;
  R_out << sel_srvind(1) << endl;
 // R_out << "selectivity industry survey males 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey males 2009" << endl;
  R_out << sel_srvind(2) << endl;
 // R_out << "selectivity nmfs industry survey females 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey females 2009" << endl;
  R_out << sel_srvnmfs(1) << endl;
 // R_out << "selectivity nmfs industry survey males 2009: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey males 2009" << endl;
  R_out << sel_srvnmfs(2) << endl;
 // R_out << "selectivity industry survey females 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey females 2010" << endl;
  R_out << sel_srv10ind(1) << endl;
 // R_out << "selectivity industry survey males 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity industry survey males 2010" << endl;
  R_out << sel_srv10ind(2) << endl;
 // R_out << "selectivity nmfs industry survey females 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey females 2010" << endl;
  R_out << sel_srv10nmfs(1) << endl;
 // R_out << "selectivity nmfs industry survey males 2010: '27.5','32.5','37.5','42.5','47.5','52.5','57.5','62.5','67.5','72.5','77.5','82.5','87.5','92.5','97.5','102.5','107.5','112.5','117.5','122.5','127.5','132.5'"<< endl;
  R_out << "$selectivity nmfs industry survey males 2010" << endl;
  R_out << sel_srv10nmfs(2) << endl;

  R_out << "$survey CV" << endl;
  R_out << cv_srv1o << endl;
  
  R_out << "$Observed Lbar" <<endl;
  R_out << Lbar <<endl;

  R_out << "$Predicted Lbar new shell" <<endl;
  R_out << Lbar_hat_new <<endl;
  
  R_out << "$Predicted Lbar old shell" <<endl;
  R_out << Lbar_hat_old <<endl;
  
  
 }


GLOBALS_SECTION
 #include <math.h>
 #include <admodel.h>
  #include <time.h>
 ofstream CheckFile;
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;
 ofstream R_out;

// EVALUATE THE GAMMA FUNCTION
 dvariable gammln(dvariable& xx)
 {
   RETURN_ARRAYS_INCREMENT();	
        dvariable x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
               0.1208650973866179e-2,-0.5395239384953e-5}; 
        int j;
        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) 
          {
            ser += cof[j]/(y+1.);
          }
         dvariable value_=-tmp+log(2.5066282746310005*ser/x);
    RETURN_ARRAYS_DECREMENT();         
         return(value_); 
 }

RUNTIME_SECTION
//one number for each phase, if more phases then uses the last number
  maximum_function_evaluations 500,1000,1000,1000,1000,1000,3000
  convergence_criteria 1,1,1,1,.01,.001,1e-3,1e-3
TOP_OF_MAIN_SECTION
  arrmblsize = 3500000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000); // this may be incorrect in
  // the AUTODIF manual.
  gradient_structure::set_CMPDIF_BUFFER_SIZE(10000000);
  time(&start);
  CheckFile.open("Check.Out");
  R_out.open("R_input.txt");

FINAL_SECTION
 time(&finish); 
 elapsed_time = difftime(finish,start);
 hour = long(elapsed_time)/3600;
 minute = long(elapsed_time)%3600/60;
 second = (long(elapsed_time)%3600)%60;
 cout << endl << endl << "Starting time: " << ctime(&start);
 cout << "Finishing time: " << ctime(&finish);
 cout << "This run took: " << hour << " hours, " << minute << " minutes, " << second << " seconds." << endl << endl;



