library(adnuts)
library(snowfall)
library(shinystan)

setwd('cod_fast')
mle <- r4ss::read.admbFit('ss3')
par.names <- mle$names[1:mle$nopar]
setwd('..')

sfStop()
sfInit(parallel=TRUE, cpus=3)
cc <- 3
inits <- rep(list(mle$est[1:mle$nopar]), cc)
sfExportAll()

fit.nuts <-
  sample_admb('ss3', iter=50, init=inits, par.names=par.names,
              parallel=TRUE,
  chains=cc, warmup=2, dir='cod_fast', cores=3,
  control=list(max_treedepth=8, stepsize=.1))
launch_shinystan_admb(fit.nuts)

fit.nuts <-
  sample_admb('ss3', iter=50, init=inits, par.names=par.names,
              parallel=TRUE,
  chains=cc, warmup=2, dir='cod_fast', cores=3,
  control=list(algorithm='RWM'))
