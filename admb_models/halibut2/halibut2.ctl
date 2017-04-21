##################################
# Halibut control file for 2015
# Coastwide short time-series model
##################################

1 # Number of growth morphs
1 # Number of sub-morphs 
0 # Number of time-varying block patterns

### Growth and mortality setup section ###
0.5 # Fraction female at birth
0   # M type: 0=1Parm; 1=N_breakpoints; 2=Lorenzen, 3=age specific, 4=agespec withseasinterpolate

###### Bypass all growth and fecundity calculations with empirical weight-at-age and spawning output at age ######
1   # Growth model: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=not implemented; 4=not implemented
5   # Growth_Age_for_L1
18  # Growth_Age_for_L2 (999 to use as Linf)
0   # SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0   # CV_Growth_Pattern: 0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)

######
# Invoke empirical weight-at-age
5   # Maturity option: 1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=read fec and wt from wtatage.ss
######

### More bypassed parameters
8   # First Mature Age
1   # Fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0   # Hermaphroditism option: 0=none; 1=age-specific fxn
1   # Parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
1   # Env/block/dev_adjust_method (1=standard; 2=logistic transform keeps in base parm bounds; 3=standard w/ no bound check)

### Natural mortality ###
# Lo	 Hi	 	Init	 Prior	 Prior	 Prior	Param	Env	Use	Dev		Dev		Dev	Block	block
# bnd	 bnd 	value	 mean	 type	 SD		phase	var	dev	minyr	maxyr	SD	design	switch
# Female natural mortality
  0.05   0.25 	0.15 	 0.15    0       0.02   -8    	0   0   0   	0 		0 	0 		0		# M females

#######################
### Growth bypassed ###
  25 	80 	 49   65 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # Len Age-6 Females
  70    200  120  120 	-1 	99 	-50		0 	0 	0 	0 	0 	0 	0 # Len Age-18 Females
  0.01 	0.4  0.1  0.2 	-1 	99 	-50 	0 	0 	0 	0 	0  	0 	0 # K Females
  0.05 	0.3  0.2  0.1 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # CV young Females
  0.05 	0.35 0.25 0.1 	-1 	99 	-50 	0 	0 	0	0 	0 	0 	0 # CV old Females
#######################

# Lo	 Hi	 	Init	 Prior	 Prior	 Prior	Param	Env	Use	Dev		Dev		Dev	Block	block
# bnd	 bnd 	value	 mean	 type	 SD		phase	var	dev	minyr	maxyr	SD	design	switch
# Male natural mortality
  0.05   0.25 	0.15 	 0.15    -1      99 	5    	0   0   0   	0 		0 	0 		0          # M males

#######################
### Growth and biological parameters bypassed ###
  25 	80 	52 	65 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # Len Age-6 males
  50    160 	98 	90 	-1 	99 	-50	0 	0 	0 	0 	0 	0 	0 # Len Age-18 males
  0.01 	0.4 	0.07 	0.2 	-1 	99 	-50 	0 	0 	0 	0 	0  	0 	0 # K males
  0.05 	0.3 	0.2 	0.1 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # CV young males
  0.05 	0.35 	0.15 	0.1 	-1 	99 	-50 	0 	0 	0	0 	0 	0 	0 # CV old males
# Weight-length, maturity, and fecundity ignored (empirical vector read in is used)
  -1 	1   0.00000691 	0.01 	-1 	99	-50 	0	0 	0 	0 	0 	0 	0 # W-L a
  2 	4   3.24019356  3.24 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # W-L b
# Maturity ignored (empirical vector read in is used)
  50 60 55 55 -1 0.8 -3 0 0 0 0 0 0 0 # Mat50%_Fem
  -3 3 -0.25 -0.25 -1 0.8 -3 0 0 0 0 0 0 0 # Mat_slope_Fem
  -3 	3 	1 	1 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # Eggs/kg_inter_Fem
  -3 	3 	0 	0 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # Eggs/kg_slope_wt_Fem
  -1 	1   0.00000691 	0.01 	-1 	99	-50 	0	0 	0 	0 	0 	0 	0 # W-L a
  2 	4   3.24019356  3.24 	-1 	99 	-50 	0 	0 	0 	0 	0 	0 	0 # W-L b
# Misc unused vectors
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_GP_1
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Area_1
 0 0 0 0 -1 0 -4 0 0 0 0 0 0 0 # RecrDist_Seas_1
 0 2 1 0 -1 0 -4 0 0 0 0 0 0 0 # CohortGrowDev
 0 0 0 0 0 0 0 0 0 0 #_seasonal_effects_on_biology_parms
########################

########################
### Spawner-Recruit section ###
3 # S-R functions: 2=Ricker, 3=B-H, 4=SCAA, 5=Hockey, 6=B-H_flattop, 7=3-par survival
# LO	HI	INIT	PRIOR	PR_type	SD		PHASE
  8     15  12.96   10.3    -1      99      1		# Log(R0)
  0.2   1.0 0.75    0.9     1       0.05    -50     # Steepness
  0.05  1.1 0.9     0.6     -1      0.8     -50     # SigmaR
  -5    5   0.1     0       -1      1       -50     # Environmental link coefficient
  -1    1   0       0       -1      1       6       # R1 offset
  0     0   0       0       -1      0       -99     # Autocorrelation

0 # Environmental index for S-R parameters
0 # Environmental link to: 0=none, 1=rec devs, 2=log(R0), 3=Steepness

# Recruitment deviation section
2		# Deviation type do_recdev:  0=none; 1=devvector; 2=simple deviations
1977 	# Start main rec devs
2007 	# End main rec devs
2		# Phase for main rec devs 
1		# Read advanced options
1976	# Start early rec devs
-9		# Phase for early rec devs_recdev_early_phase
9		# Phase for forecast rec devs
1		# Lambda for Fcast_recr_like occurring before endyr+1
1975 	# Last_early_yr_nobias_adj_in_MPD
1976 	# First_yr_fullbias_adj_in_MPD
2003 	# Last_yr_fullbias_adj_in_MPD
2007 	# First_recent_yr_nobias_adj_in_MPD
-1		# Max_bias_adj_in_MPD (-1 to override ramp and set to 1.0 for all deviations)
0		# Period of cycles in recruitment
-4		# Recruitment deviation lower bound 
6		# Recruitment deviation upper bound 
0		# Read initial values for recruitment deviations: 0 = no, 1 = yes

### Fishing mortality section ###
0.3		# F target for tuning early phases
-2013   # F target year (negative value disables this feature)
1 	    # F Method: 1=Popes; 2=instan. F; 3=hybrid
0.8		# Max F/harvest rate, depends on F_Method

# Initial equilibrium F parameters
# LO  HI INIT PRIOR PR_type SD  PHASE
  -1  1  0    0.0   -1      99  -1	  # Commercial
  -1  1  0    0.0   -1      99  -1    # Discards
  -1  1  0    0.0   -1      99  -1    # Bycatch
  -1  1  0    0.0   -1      99  -1    # Sport
  -1  1  0    0.0   -1      99  -1    # Personal use

### Catchability section ###
# A) Density dependence  
# B) Environmental variable  
# C) Extra standard error 
# D) Type: <0=mirror, 0=median_float, 1=mean_float, 2=parameter, 3=parm_w_random_dev, 4=parm_w_randwalk, 5=mean_unbiased_float_assign_to_parm
# A  B  C  D
  0  0  1  4    # 1 Commercial
  0  0  0  1    # 2 Discards
  0  0  0  1    # 3 Bycatch
  0  0  0  1    # 4 Sport
  0  0  0  1    # 5 Personal use
  0  0  1  1    # 6 Survey
  
1 # Q parameter detail for time-varying parameters (1=one line for each)

# Q parameters
# LO   HI    INIT  PRIOR  PR_type   SD     PHASE
  0    0.5   0.0001  0.05 -1        50     -8      # 1 Commercial index extra SD
  0    0.5   0.0001  0.05 -1        50     -8      # 6 Survey index extra SD
  -12  0     -7.5  0.05   -1        50     3       # 1 Fishery index log(1996 base parameter)
  -2   2     0     0.0    0         0.003  4       # 1997 deviation
  -2   2     0     0.0    0         0.003  4       # 1998 deviation
  -2   2     0     0.0    0         0.003  4       # 1999 deviation
  -2   2     0     0.0    0         0.003  4       # 2000 deviation
  -2   2     0     0.0    0         0.003  4       # 2001 deviation
  -2   2     0     0.0    0         0.003  4       # 2002 deviation
  -2   2     0     0.0    0         0.003  4       # 2003 deviation
  -2   2     0     0.0    0         0.003  4       # 2004 deviation
  -2   2     0     0.0    0         0.003  4       # 2005 deviation
  -2   2     0     0.0    0         0.003  4       # 2006 deviation
  -2   2     0     0.0    0         0.003  4       # 2007 deviation
  -2   2     0     0.0    0         0.003  4       # 2008 deviation
  -2   2     0     0.0    0         0.003  4       # 2009 deviation
  -2   2     0     0.0    0         0.003  4       # 2010 deviation
  -2   2     0     0.0    0         0.003  4       # 2011 deviation
  -2   2     0     0.0    0         0.003  4       # 2012 deviation
  -2   2     0     0.0    0         0.003  4       # 2013 deviation
  -2   2     0     0.0    0         0.003  4       # 2014 deviation

### Length-based selectivity ###
# A) Pattern, B) Discard, C) Male, D) Special
# A  B  C  D 
  0  0  0  0   # 1 Commercial
  0  2  0  0   # 2 Discards
  0  0  0  0   # 3 Bycatch
  0  0  0  0   # 4 Sport
  0  0  0  0   # 5 Personal use
  0  0  0  0   # 6 Survey

### Age-based selectivity ###
# A) Pattern, B) ---, C) Male, D) Special
# A  B  C  D 
  20 0  3  0   # 1 Commercial
  20 0  4  0   # 2 Discards
  20 0  0  0   # 3 Bycatch
  20 0  0  0   # 4 Sport
  15 0  0  4   # 5 Personal use
  20 0  3  0   # 6 Survey

#############################################################################################################################
### Length-based selectivity parameters ###
#############################################################################################################################
### Discards: Retention and mortality ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
# Discard retention logistic (this corrects the input and sets the scale)
  0	    2	 1      100    -1	 99		-50	    0   0   0	 	0	  0		0		0		# Retention inflection
  0  	1	 0.1    100    -1    99 	-50   	0   0   0 	    0 	  0 	0		0    	# Retention slope
  0	    2	 0.0016 100    -1    99 	-50	    0   0   0	 	0	  0		0		0    	# Asymptotic retention 
  -1  	1    0      9      -1    99 	-50 	0   0   0 	    0 	  0     0		0    	# Male offset
# Discard mortality logistic (interpreted as DMR, but 0.15864 = 16%)
  0	    2	 1      100    -1	 99		-50	    0   0   0	 	0	  0		0		0		# Mortality inflection
  0  	1	 0.1    100    -1    99 	-50   	0   0   0 	    0 	  0 	0		0    	# Mortality slope
  0.01  0.99 0.15864 0.15864 0   0.03 	-10	    0   0   0	 	0	  0		0		0    	# Asymptotic mortality 
  -1  	1    0      9      -1    99 	-50 	0   0   0 	    0 	  0     0		0    	# Male offset
#############################################################################################################################
### Age-based selectivity parameters ###
#############################################################################################################################
### Fishery: double normal ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  0	    25	 15     100    -1	 99		4 	    0   3   1997 	2014  0.06	0		0		# Age at peak selectivity
  1  	30   25     9      -1    99 	-50   	0   0   0 	    0 	  0 	0		0    	# Top width
  0	    25	 2.9    100    -1    99 	4 	    0   3   1997 	2014  0.09	0		0    	# Ascending width
  0.1  	100  99     9      -1    99 	-50 	0   0   0 	    0 	  0     0		0    	# Descending width
  -1000	0	 -999   100    -1    99 	-50 	0   0   0 	    0 	  0		0 	    0       # Initial (-999 to ignore)
  -1000 25   -999   9      -1    99     -50 	0   0   0 	    0 	  0     0		0    	# Final
# Fishery: male parameter offset to females
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  -8	8    0		10     -1    99     5       0   0   0       0     0 	0 	    0		# Additive to female peak
  -4  	5    0		9      -1    99     5       0   0   0       0     0 	0		0       # Additive to ascending width
  -10	0    0      9      -1    99     -50     0   0   0 	    0 	  0		0 	    0       # Additive to descending width
  -1	1    0      9      -1    99 	-50     0   0   0 	    0 	  0	    0 	    0       # Additive to final
  0.1	1.00 0.7    9      -1    99 	5       0   3   1997 	2014  0.03	0 	    0       # Asymptote for males 
#############################################################################################################################
### Discards: double normal ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  7	    14	 12		100    -1    99 	5 	    0   0   0 	    0     0	    0 	    0		# Age at peak selectivity
  -8    8	 -2     9      -1    99 	-5   	0   0   0 	    0     0 	0  	    0       # Top width
  -4	5	 2		100    -1    99 	5 	    0   0   0 	    0 	  0		0	 	0    	# Ascending width
  -4  	5    3		9      -1    99     5 	    0   0   0 	    0 	  0		0       0    	# Descending width
  -1003	0	 -1002  100    -1    99 	-50 	0   0   0 	    0 	  0	    0	 	0       # Init (-999 to ignore)
  -8	8	 4		9      -1    99     5	 	0   0   0 	    0 	  0     0       0    	# Final
# Discards: female parameter offset to males
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  -8	8    -1		10     -1    99     6       0   0   0       0     0 	0 	    0		# Additive to male peak
  -5  	5    -1		9      -1    99     6       0   0   0 	    0 	  0		0		0       # Additive to ascending width
  -10	0    0      9      -1    99     -50     0   0   0 	    0 	  0		0 	    0       # Additive to descending width
  -1	1    0      9      -1    99 	-50     0   0   0 	    0 	  0	    0 	    0       # Additive to final
  0.1	1.0  0.8    9      -1    99 	6       0   0   0 	    0 	  0		0 	    0       # Asymptote for females 
#############################################################################################################################

### XXXX !!!!!!!!!!!!! XXXXX
## Cole fixed these parameters at MLEs by hand and turned of phase. 
### Bycatch: Double normal ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  2	    8	 3.6354	10	   -1    99 	-5 	    0   0   0 	    0     0	    0 	    0		# Age at peak selectivity
  -8    8	 -4     9      -1    99 	-5   	0   0   0 	    0     0 	0  	    0       # Top width
  -8	8	 0		100    -1    99 	-5 	    0   0   0 	    0 	  0		0	 	0    	# Ascending width
  -6  	6    2.5805	9      -1    99     -5 	    0   0   0 	    0 	  0		0       0    	# Descending width
  -1003	0	 -1002  100    -1    99 	-50 	0   0   0 	    0 	  0	    0	 	0       # Init (-999 to ignore)
  -8	8	-0.2808 9      -1    99     -5	 	0   0   0 	    0 	  0     0       0    	# Final
#############################################################################################################################
### Sport: Double normal ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  2	    14	 8.485	10	   -1    99 	-5 	    0   0   0 	    0     0	    0 	    0		# Age at peak selectivity
  -8    8	 -4     9      -1    99 	-5   	0   0   0 	    0     0 	0  	    0       # Top width
  -8	8	 1.6645	100    -1    99 	-5 	    0   0   0 	    0 	  0		0	 	0    	# Ascending width
  -6  	6   1.40205 9      -1    99     -5 	    0   0   0 	    0 	  0		0       0    	# Descending width
  -1003	0	 -1002  100    -1    99 	-50 	0   0   0 	    0 	  0	    0	 	0       # Init (-999 to ignore)
  -8	8	0.2571	9      -1    99     -5	 	0   0   0 	    0 	  0     0       0    	# Final
#############################################################################################################################
## END of changes by Cole


### Personal use: mirror sport ###
#############################################################################################################################
### Survey: Simple double normal  (option 20) ###
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  0		25	 12     100    -1    99 	4 	    0   3   1998 	2014  0.55	0 	    0    	# Age at peak selectivity
  1		30   25     9      -1    99 	-50   	0   0   0 	    0 	  0 	0  	  	0    	# Top width
  0		25	 2.9    100    -1    99 	4 	    0   3   1998 	2014  0.1	0 	    0    	# Ascending width
  0.1	100  99     9      -1    99 	-50 	0   0   0 	    0 	  0     0  	  	0    	# Descending width
  -1000	0	 -999   100    -1    99 	-50 	0   0   0 	    0 	  0		0 	    0       # Init (-999 to ignore)
  -1000	0    -999   9      -1    99     -50 	0   0   0 	    0 	  0     0       0    	# Final
# Survey male parameter offset to females
# Lo	Hi	 Init	Prior  Prior Prior	Param	Env	Use	Dev	  	Dev	  Dev	Block	Block
# bnd	bnd  value  mean   type  SD		phase	var	dev	minyr	maxyr SD	design	switch
  0	    8	 2      10     -1    99     5       0   0   0       0     0 	0 	    0       # Additive to female peak
  -1  	5    0      9      -1    99     5       0   0   0       0     0 	0  	   	0       # Additive to ascending width
  -10	0	 0      9      -1    99     -50     0   0   0 	    0 	  0	   	0 	    0       # Additive to descending width
  -1	1	 0      9      -1    99 	-50     0   0   0 	    0 	  0	    0 	    0       # Additive to final
  0.1	1.00 0.8    9      -1    99 	5       0   3   1998 	2014  0.06	0 	    0       # Asymptote for males
#############################################################################################################################

### Extra parameters for selectivity deviations
5 # Phase for selectivity parameter deviations
1 # Selectivity parameter adjustment method

0  # No tagging parameters

### Variance adjustments ###
1 # Variance_adjustments by fleet
# 1		2     3	     4	   5	6    
# Comm. Disc. Bycat. Sport Pers Survey 
  0		0     0      0     0    0       # Additive survey SE constant
  0		0     0      0     0    0       # Additive discard SE constant
  0		0     0      0	   0    0       # Additive mean body weight SE constant
  0		0     0      0     0    0       # Length comp. N multiplier
  0.1	0.007 0.1	 0.1   1    0.3		# Age comp N multiplier
  0		0     0      0     0    0       # Size-at-age N multiplier
#############################################################################################################################

### Lambdas ###
1 # Maximum phase for lambda implementation
1 # Include SD offset in likelihoods: 0=no, 1=yes

0 # Number of changes to make to default (value=1.0) Lambdas
# Type: 1=index, 2=discard, 3=mnwt, 4=length, 5=age, 6=SizeFreq, 7=sizeage, 8=catch, 9=init_equ_catch,
#       10=recrdev, 11=parm_prior, 12=parm_dev, 13=CrashPen, 14=Morphcomp, 15=Tag-comp, 16=Tag-negbin.
# Type  Fleet  Phase  Lambda  Sizemthd

# Extra SD reporting section
0 # 0=none, 1=read specs below

999 # End of file marker