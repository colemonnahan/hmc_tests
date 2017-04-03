setwd('ss3_models/')
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)
cores <- 3
sfStop()
sfInit(parallel=TRUE, cpus=cores)
chains <- cores
inits <- rep(list('mle'), chains)
m <- 'cod_fast'
mle <- r4ss::read.admbFit('cod_fast/ss3')
inits <- rep(list(mle$est[1:mle$nopar]), chains)
par.names <- mle$names[1:mle$nopar]
sfExportAll()

cod <- sample_admb(model='ss3', chains=chains, dir=m, iter=10,
                       init=inits, par.names=par.names, warmup=2,
                       parallel=TRUE, cores=cores,
                       control=list(stepsize=.2, max_treedepth=3))

launch_shinystan_admb(codrich)
