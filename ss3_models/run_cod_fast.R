setwd("ss3_models/")
library(adnuts)
library(snowfall)
library(shinystan)

m <- 'cod_fast'
setwd(m)
mle <- r4ss::read.admbFit('ss3')
par.names <- mle$names[1:mle$nopar]
setwd('..')
sfStop()
reps <- 3
sfInit(parallel=TRUE, cpus=reps)
inits <- rep(list(mle$est[1:mle$nopar]), reps)
td <- 8
warmup <- 25
iter <- 225
sfExportAll()

fit.nuts <-
  sample_admb('ss3', iter=iter, init=inits, par.names=par.names,
              parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=.3))
stats.nuts <- data.frame(rstan::monitor(fit.nuts$samples, warmup=warmup, probs=.5, print=FALSE))
perf.nuts <- min(stats.nuts$n_eff)/sum(fit.nuts$time.total)
launch_shinystan_admb(fit.nuts)


fit.rwm <-
  sample_admb('ss3', iter=iter, init=inits, par.names=par.names,
              parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(algorithm='RWM'))
launch_shinystan_admb(fit.rwm)
stats.rwm <- data.frame(rstan::monitor(fit.rwm$samples, warmup=warmup, probs=.5, print=FALSE))
perf.rwm <- min(stats.rwm$n_eff)/sum(fit.rwm$time.total)
