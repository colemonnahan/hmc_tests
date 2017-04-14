##############################
# Halibut model settings for 2015
# Coastwide short time-series model
##############################

1		# Benchmarks: 0=skip, 1=calc F_spr,F_btgt,F_msy 
2		# MSY: 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr) 
0.3 	# SPR target (e.g. 0.40)
0.3 	# Biomass target (e.g. 0.40)
#_Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF, end_relF (enter actual year, or values of 0 or -integer to be rel. endyr)
-10 0 -3 0 -3 0
1 		# Bmark_relF_Basis: 1 = use year range; 2 = set relF same as forecast below
4 		# Forecast: 0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
4 		# Number of forecast years 
1.0 	# F scalar (only used for Do_Forecast==5)
# Fcast_years:  beg_selex, end_selex, beg_relF, end_relF  (enter actual year, or values of 0 or -integer to be rel. endyr)
 -3 0 -3 0
1 		# Control rule method (1=catch=f(SSB) west coast; 2=F=f(SSB) ) 
0.3 	# Control rule Biomass level for constant F (as frac of Bzero, e.g. 0.40); (Must be > the no F level below) 
0.2 	# Control rule Biomass level for no F (as frac of Bzero, e.g. 0.10) 
1.0 	# Control rule target as fraction of Flimit (e.g. 0.75) 
3 		# Number of forecast loops (1=OFL only; 2=ABC; 3=get F from forecast ABC catch with allocations applied)
3 		# First forecast loop with stochastic recruitment
0 		# Forecast loop control #3 (reserved for future bells&whistles) 
0 		# Forecast loop control #4 (reserved for future bells&whistles) 
0 		# Forecast loop control #5 (reserved for future bells&whistles) 
2019  	# FirstYear for caps and allocations (should be after years with fixed inputs) 
0		# SD of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error)
0 		# Do West Coast gfish rebuilder output (0/1) 
1999 	# Rebuilder: first year catch could have been set to zero (Ydecl)(-1 to set to 1999)
2002 	# Rebuilder: year for current age structure (Yinit) (-1 to set to endyear+1)
1 		# Fleet relative F:  1=use first-last alloc year; 2=read seas(row) x fleet(col) below
# Note that fleet allocation is used directly as average F if Do_Forecast=4 
2 		# Basis for fcast catch tuning and for fcast catch caps and allocation  (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum)
# Maximum total catch by fleet (-1 to have no max), must enter value for each fleet
 -1 -1 -1 -1 -1
# Maximum total catch by area (-1 to have no max), must enter value for each area 
 -1
# Fleet assignment to allocation group (enter group ID# for each fleet, 0 for not included in an allocation group)
 0 0 0 0 0
15 		# Number of forecast catch levels to input (else calc catch from forecast F) 
2 		# Basis for input Fcast catch:  2=dead catch; 3=retained catch; 99=input Hrate(F) (units are from fleetunits; note new codes in SSV3.20)
# yr season flt catch
2016 1 1 10
2016 1 2 10
2016 1 3 10
2016 1 4 10
2016 1 5 10
2017 1 1 10
2017 1 2 10
2017 1 3 10
2017 1 4 10
2017 1 5 10
2018 1 1 10
2018 1 2 10
2018 1 3 10
2018 1 4 10
2018 1 5 10

999 # End of forecast file