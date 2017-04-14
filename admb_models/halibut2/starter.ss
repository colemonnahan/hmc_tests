############################
# Halibut model settings for 2015
# Coastwide short time-series model
############################

# Input file names
halibut.dat
halibut.ctl

0 		# Initial values: 0 = control file, 1 = ss3.par
0 		# Run display detail (0,1,2)
1 		# Include age-structured reports: 0 = yes, 1 = no 
0 		# Write detailed checkup.sso file (0,1) 
0 		# Write parm values to ParmTrace.sso (0=no,1=good,active; 2=good,all; 3=every_iter,all_parms; 4=every,active)
0 		# Write to cumreport.sso (0=no,1=like&timeseries; 2=add survey fits)
1 		# Include prior likelihood for non-estimated parameters (0,1) 
0 		# Use "soft bounds" to aid convergence (0,1)
0 		# Number datafiles to produce: 1st is input, 2nd is estimates, 3rd and higher are bootstrap
50     	# Turn off estimation for parameters entering after this phase
1		# MCeval burn interval
1 		# MCeval thin interval
0 		# Jitter initial parm value by this fraction
-1 		# min yr for sdreport outputs (-1 for styr)
-2 		# max yr for sdreport outputs (-1 for endyr; -2 for endyr+Nforecastyrs)
0 		# Number of individual STD years 
0.0001  # final convergence criteria (e.g. 1.0e-04) 
0 		# Retrospective year relative to end year (e.g. -4)
8 		# Minimum summary biomass age
1 		# Depletion basis: denom is: 0=skip; 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr
0.3		# Fraction (X) for Depletion denominator (e.g. 0.4)
4 		# SPR_report_basis:  0=skip; 1=(1-SPR)/(1-SPR_tgt); 2=(1-SPR)/(1-SPR_MSY); 3=(1-SPR)/(1-SPR_Btarget); 4=rawSPR
1 		# F_report_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num); 3=sum(Frates); 4=true F for range of ages
0 		# F_report_basis: 0=raw; 1=F/Fspr; 2=F/Fmsy ; 3=F/Fbtgt
999   	# End of starter file
