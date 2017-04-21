#===============input parameters====================
#==model structure==
# styr_fut : start year future projections
2010
# endyr_fut : end year future
2020
# nsellen : number of length bins female
22
# nsellen_srv1 : number of length bins male
22
# p_const : small number added to length comps
0.001

#==catchability==
# q1 : catchability
1.0

#==mortality sources==
# M_in : Nat mort immature female/male
0.23  0.23
# M_matn_in : nat mort mature new shell female/male
0.23  0.23
# M_mato_in : nat mort mature old shell female/male
0.23  0.23
# m_disc : assumed mortality for pot fishery bycatch
 0.3
# m_trawl : mortality for trawl bycatch
 0.8

#==recruitment===
# median_rec : median rec for last years and future years
692528.  542191.
# median_rec_yrs : number of years to fix recruits
0
# var_rec_obs : sd at recruitment
3.6
# sd_var_rec : 
0.1

#==selectivity==
# sel_som : parameters for somerton selectivity curve (not standard logistic)
# model:  rat.m ~ a/(1 + exp(-(b + c * width.m)))
#      a          b            c
# 0.8418 -2.6466  0.0354
0.8418
-2.6466
0.0354
7.905
0.081
# number of final years to average over for fishery sel50 par when time varying
4

#==weight at length
# juvenile female
# alpha_wt_imm_f
0.001047
# beta_wt_imm_f
2.708367
# alpha_wt_mat_f
0.001158
# beta_wt_mat_f
2.708793
# alpha_wt_m
0.000267
# beta_wt_m
3.097253

#==growth==
#growth parameters from ST. Marie with sd(for priors not estimated)
# linff_obs
5.9417
# sd_linff
0.59
# linfm_obs
6.0652
# sd_linfm
.5
# growthkf_obs
1.1067
# sd_growthkf
0.25
# growthkm_obs
1.1238
# sd_growthkM
0.25
# af_obs : female a
6.773   # this fixes intercept to be same for males and females - average of the two intercepts
# sd_af :sd a
0.1
# am_obs : male a
6.593
# sd_am : sd from SOmerton
0.73
# bf_obs : female b
1.05      # this gives approx same growth with intercept at 6.773 same for males and females
# sd_bf : sd b
0.10
# bm_obs : male b
1.17
# sd_bm : sd b
0.012
# a1_obs : a1
 -4.0
# sd_a1 : a1 sd from somerton
1.36
# b1_obs : b1
1.46
# sd_b1 : b1 sd from somerton
0.052
# sd_meetpt : sd for meet point at 36.1mm premolt width
0.2
# var_last_obs : sd at max length
10.5
# sd_var_last : 
0.2
# mate_ratio : mating ratio for estimation of effective spawning biomass
 1.72

#==maturity==
# fraction_new_error : fraction of new shell animals that are 0<1 yr from molting (error in shell condition)
1.0
# nages : number of ages to track for old shell mature animals
10
# matest_n : number of parameters (one per length bin) to estimate for maturity estimation females
 8
# matestm_n : number of parameters (one per length bin) to estimate for maturity estimation males
 17

#===============weighting===========================
# wt_like : weights for smooth sel - not used
 10. 40. 40. 40.
# second number is weight for retained catch - multiplied by 10* in later phases
 2000. 1000. 100000. 10000.
# like_wght :  weight for retained length comps
1.0
#weight tot length comps
1.0
#weight for female pot fishery discard length comps
0.2
#weight for survey length comps
1.0
#survey biomass cv factor multiplied times the obs cv's
#female
1.0
# weight for retained catch biomass
10.0
#weight for trawl bycatch length comps
0.25
# like_wght_mbio : survey biomass cv factor multiplied times the obs cv's
1.0
# like_wght_rec : weight for male recruitment devs
1.0
# like_wght_recf :weight for female rec devs
1.0
# like_wght_sexr :weight for sex ratio
0.001
# like_wght_sel50 :weight for fishery 50% selectivity deviations
100
# like_wght_fph1 :weight for fmorts pot fishery phase 1 like_wght_fph1
10.0
# like_wght_fph2 :weight for fmorts pot fishery after phase 1 like_wght_fph2
2.0
# like_wght_fdev :weight for fmort deviations pot fishery like_wght_fdev
0.1
# wght_total_catch : weight for total (retained plus discard) catch likelihood
20.0
# wght_female_potcatch : weight for female pot bycatch likelihood
10.0
# cput_cv : CV for the fit to pot fishery likelihood - larger means less weight on fit (CSS: this is not a CV, it is a SD for a normal likelihood)
5
# wt_lmlike : wt for fitting large male (>101) survey numbers
0.001

#==constants taken from objective function
# old_shell_constraint : weighting for a penalty place on initial numbers of old shell in small length bins
 0.000001
# disclen_mult : multipler on the discard at length likelihood component
 1.0 1.0
# selsmo_wght : penalty on smoothness of selsmo10 and 09ind; only active if parameter estimated
 5.
# femSel_wght : penalty on smoothness of time-varying selectivity for females
 10.
# init_yr_len_smooth_f : smoothness on female lengths in initial year
 2
# init_yr_len_smooth_m : smoothness on female lengths in initial year
 20
# natm_mult_wght : weighting on the multiplier for M for males
0.5
# natm_mult_var : variance for prior estimated M multiplier
0.054
# natm_Immult_wght : weighting on the multiplier for M for immature
0.5
# natm_Immult_var : variance for prior estimated M multiplier for immature (was 1)
0.15
# smooth_mat_wght : smoothing for estimated maturity (2)
50
# mat_est_wght : weighting for the maturity estimates
0.5
# mat_est_vs_obs_sd : standard deviation for deviations from observed maturity
0.05
# smooth_mat_wght_f : smoothness for estimated maturity for females (2)
50
# mat_est_vs_obs_sd_f : standard deviation for deviation from observed maturity curve
0.05
# growth_data_wght_m : weight on the growth data (was 1)
 10
# growth_data_wght_f : weight on the growth data (was 5)
 10
# extra_wght_ind_m : add extra weight to the male lengths from the industry survey
 2
# smooth_disc_catch : weighting onthe smoothness of the discarded female catch
0.0001
# disc_catch_wght_fem : multiply like_wght(6) by this to increase weight on female discards
30

#===============phase estimation====================
# fmort_phase : phase to estimate fishing mortality
1
# rec_phase : phase to estimate recruitment
2
# growth_phase : phase to estimate growth parameters
6
# growth_phase2 : phase to estimate growth parameters
6
# maturity_phase : phase to estimate maturity
5
# natM_phase : estimate multiplier for natural mortality
5
# phase_moltingp : molting prob mature male phase (set to -1 for terminal molt)
-1
# phase_fishsel : phase for estimating dome shaped fishery selectivity parameters (-1 for no dome shpaed)
-1
# survsel_phase : phase to estimate 1989 to present period survey selectivity(set to 4)
4
# survsel1_phase : phase to estimate 1978 to 1982 and 1983-1988 survey selectivity(set to 4) if negative then fixed at Somerton and Otto
# to fix all survey sel at somerton and otto set both to negative
4
# phase_fut : phase to estimate future projections
9
# phase_logistic_sel : phase for estimating fishery and bycatch selectivity curves
4
# phase_selcoffs : phase for estimating smooth selectivities
-1

#===============code switches=======================
# growth_switch : growth switch males 1 for st. marie 2 for st.marie then same after 90mm
1
# somertonsel : if want to use somerton sel (2009 industry survey) set < 0
1
# monot_sel : switch for smooth selectivity fishery
1
# monot_sel_srv1 : switch for smooth selectivity survey
1
# maturity_switch : if greater than 0 then use maturity logistic curve instead of fractions by year for probability
#of new shell maturity
2
# f_penalties : switch for penalties on F
0
# retro_years : number of years from endyr to chop off when fitting
0