setwd("ss3_models/")
library(adnuts)
library(snowfall)
library(shinystan)

## Run model if needed
m <- 'cod_fast'
setwd(m)
system("ss3")
mle <- r4ss::read.admbFit('ss3')
replist <- r4ss::SS_output(getwd(), covar=TRUE)
unlink('plots', TRUE)
r4ss::SS_plots(replist)
par.names <- mle$names[1:mle$nopar]
setwd('..')


sfStop()
reps <- 5
sfInit(parallel=TRUE, cpus=reps)
inits <- rep(list(mle$est[1:mle$nopar]), reps)
td <- 12
warmup <- 50
iter <- 1000
hh <- .1                           # hours to run
eps <- .1
mm <- diag(length(par.names))
sfExportAll()

fit.nuts <-
  sample_admb('ss3', iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=mm))
saveRDS(fit.nuts, file='fit.nuts.RDS')
stats.nuts <- data.frame(rstan::monitor(fit.nuts$samples, warmup=warmup, probs=.5, print=FALSE))
perf.nuts <- min(stats.nuts$n_eff)/sum(fit.nuts$time.total)
tt <- 1000 # thin rate
fit.rwm <-
  sample_admb('ss3', iter=iter*tt, duration=hh*60, init=inits, par.names=par.names,
              parallel=TRUE, chains=reps, warmup=warmup*tt, dir=m, cores=reps,
              thin=tt, control=list(algorithm='RWM'))
saveRDS(fit.rwm, file='fit.rwm.RDS')
stats.rwm <- data.frame(rstan::monitor(fit.rwm$samples, warmup=warmup, probs=.5, print=FALSE))
perf.rwm <- min(stats.rwm$n_eff)/sum(fit.rwm$time.total)

## How much faster is NUTS?
perf.nuts/perf.rwm

launch_shinystan_admb(fit.nuts)
launch_shinystan_admb(fit.rwm)
