### I used this file to generate two exmaple ss3sim models for use in
### testing NUTS. Clone the ss3sim/binning repo, set WD to that directory
### and run this file.

## As of 3/10/17 the get_results was broke. I suspect this is an r4ss issue
## so ignoring for now. These models were run with ss3.24o, so need to make
## sure to recompile that version for use here.

cores <- 5   # parallel cores
devtools::install_github("ss3sim/ss3sim")
devtools::install_github('r4ss/r4ss')
## major cases for binning section
D.binning <- c(2,5)
B.binning <- 2#c(0:4, 11:14)
species <- c('cod','flatfish','yellow-long')[1]
source("startup.R")
### ------------------------------------------------------------


## We split up the data poor and rich b/c to have different sample sizes
## since data poor converged a lot slower for our MARE metric
case_files <- list(F="F", B="em_binning", I="data",
                   D=c("index","lcomp","agecomp","calcomp"))
for(spp in 'cod'){
    Nsim <- 1
    scenarios <- expand_scenarios(cases=list(D=c(2), F=1, I=0,
                                  B=B.binning), species=spp)
    run_ss3sim(iterations=1:Nsim, scenarios=scenarios,
               parallel=FALSE, parallel_iterations=FALSE,
               case_folder=case_folder, om_dir=om.paths[spp],
               em_dir=em.paths[spp], case_files=case_files,
               ## bias_adjust=TRUE, bias_nsim=10,
               ## admb_options= " -maxfn 1 -phase 50",
               call_change_data=TRUE)
}
for(spp in species){
    scenarios <- expand_scenarios(cases=list(D=c(5), F=1, I=0,
                                  B=B.binning), species=spp)
    run_ss3sim(iterations=1:Nsim, scenarios=scenarios,
               parallel=FALSE, parallel_iterations=FALSE,
               case_folder=case_folder, om_dir=om.paths[spp],
               em_dir=em.paths[spp], case_files=case_files,
               ## bias_adjust=TRUE, bias_nsim=20,
               ## admb_options= " -maxfn 1 -phase 50",
               call_change_data=TRUE)
}

## Read in results
scenarios.binning.all <-
    expand_scenarios(cases=list(D=D.binning, F=1, I=0, B=B.binning),
                     species=species)
get_results_all(user=scenarios.binning.all, parallel=FALSE, over=TRUE)
xx <- read.csv("ss3sim_scalar.csv")
saveRDS(xx, file="results/results_binning.sc.RData")
## xx <- read.csv("ss3sim_ts.csv")
## saveRDS(xx, file="results/results_binning.ts.RData")
## file.remove(c('ss3sim_ts.csv', 'ss3sim_scalar.csv'))
## unlink(scen.all, TRUE)
rm(xx, scenarios, Nsim, case_files)
