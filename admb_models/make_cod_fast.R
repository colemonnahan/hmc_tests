## This file generates a very pared down ss3sim cod model. Few years and
## wide bins so it hsould run really fast per iteration.

## devtools::install_github("ss3sim/ss3sim")
## devtools::install_github('r4ss/r4ss')
case_folder <- 'cases'
em.paths <- list('cod_fast'='cod_fast_dummy/em')
om.paths <- list('cod_fast'='cod_fast_dummy/om')
library(ss3sim)
### ------------------------------------------------------------

## Craft a contrived but fast model (few years, big bins)
spp <- 'cod_fast'
case_files <- list(F="F", B="em_binning", I="data",
                   D=c("index","lcomp","agecomp"))
Nsim <- 1
scenarios <-
  expand_scenarios(cases=list(D=0, F=1, I=0, B=0), species=spp)
run_ss3sim(iterations=1:Nsim, scenarios=scenarios,
           parallel=FALSE, parallel_iterations=FALSE,
           case_folder=case_folder, om_dir=om.paths[spp],
           em_dir=em.paths[spp], case_files=case_files,
           admb_options= " -maxfn 1 -phase 50",
           call_change_data=TRUE, seed=21415)

## Move the EM folder into more convenient place, copy exe in and run it
## with Hessian. Then make plots
newdir <- 'cod_fast'
olddir <- file.path(scenarios[1], '1/em')
if(dir.exists(newdir)) unlink(newdir, TRUE)
dir.create(newdir)
trash <- file.copy(from=list.files(olddir, full.names=TRUE), to=newdir)
file.copy('cod_fast_dummy/ss3.exe', to=newdir)
unlink(scenarios, TRUE)
setwd(newdir)
system('ss3')
## !!!!! Manually set starter.ss to start from .par file and add .1 to
## jitter. This gets rid of an ADMB optimizer error, which I'm not sure
## whether it is important.

system('ss3')
## Make plots to check it
replist <- r4ss::SS_output(getwd(), covar=TRUE)
r4ss::SS_plots(replist)
setwd('..')
## DONE! Can now use it for testing


### May want to use this code to turn off some parameters if they are being
### difficult.
## data.old <- r4ss::SS_readdat("cod_fast_dummy/om/ss3.dat")
## change_e(ctl_file_in = 'cod_fast_dummy/om/ss3.ctl',
##          ctl_file_out = "change_e.ctl",
##          dat_list = data.old, for_file_in = "cod_fast_dummy/om/forecast.ss",
##          natM_type = "1parm", natM_val=NA,
##          par_name = c("_steep", "SizeSel_1P_1_Fishery"),
##          par_int = c(NA,NA), par_phase = c(-5, -5),
##          forecast_num = 0, run_change_e_full = TRUE )


## dat <- r4ss::SS_readdat(file='models/cod/om/ss3.dat')
## ## use ss3sim function to extend burn in for yellow model
## setwd('C:/Users/Cole/binning/models/cod_fast/om')
## ## Run this and do a copjle manaual updates
## ss3sim::change_year(year_begin=1, year_end=50, burnin=0,
##                     ctl_file_in='ss3.ctl',
##                     ctl_file_out='ss3_new.ctl',
##                     dat_file_in='ss3.dat',
##                     dat_file_out='ss3_new.dat',
##                     par_file_in='ss3.par',
##                     par_file_out='ss3_new.par',
##                     str_file_in='starter.ss',
##                     str_file_out='starter_new.ss')
## ## Run this to make wider population bins and create dummy data (to be
## ## overwritten)
## dat <- r4ss::SS_readdat(file='ss3.dat')
## dat2 <- change_data(dat, outfile='ss3_new.dat', fleets=list(1,2),
##                     years=list(c(1:5), 1:5), types=c('len', 'age'),
##             pop_binwidth = 3, pop_minimum_size =10,
##             pop_maximum_size = 190, write_file = TRUE)
## ## Rename the new files
