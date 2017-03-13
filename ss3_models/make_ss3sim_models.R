## This file runs a vvery pared down ss3sim cod model. Few years and wide
## bins so it hsould run really fast per iteration.

## devtools::install_github("ss3sim/ss3sim")
## devtools::install_github('r4ss/r4ss')
case_folder <- 'cases'
em.paths <- list('cod_fast'='cod_fast/em')
om.paths <- list('cod_fast'='cod_fast/om')
library(ss3sim)
### ------------------------------------------------------------


## Try to craft a fast model
spp <- 'cod_fast'
case_files <- list(F="F", B="em_binning", I="data",
                   D=c("index","lcomp","agecomp"))
Nsim <- 1
scenarios <-
  expand_scenarios(cases=list(D=c(2), F=1, I=0, B=0), species=spp)
run_ss3sim(iterations=1:Nsim, scenarios=scenarios,
           parallel=FALSE, parallel_iterations=FALSE,
           case_folder=case_folder, om_dir=om.paths[spp],
           em_dir=em.paths[spp], case_files=case_files,
           ## bias_adjust=TRUE, bias_nsim=10,
           ## admb_options= " -maxfn 1 -phase 50",
           call_change_data=TRUE)
get_results_all(user=scenarios, over=TRUE)
xx <- read.csv("ss3sim_scalar.csv")





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
