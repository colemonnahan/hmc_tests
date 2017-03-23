setwd("ss3_models/")
library(adnuts)
library(snowfall)
library(shinystan)
m <- 'cod_fast'

## Run model if needed
setwd(m)
system('admb ss3')
system("ss3")
replist <- r4ss::SS_output(getwd(), covar=TRUE)
unlink('plots', TRUE)
r4ss::SS_plots(replist)
setwd('..')


sfStop()
mle <- r4ss::read.admbFit('cod_fast/ss3')
covar <- get.admb.cov(m)$cov.bounded
N <- mle$nopar
par.names <- mle$names[1:N]
reps <- 6                        # chains/reps to run
## Draw inits from MVN using MLE and covar
inits <- rep(list(as.vector(mvtnorm::rmvnorm(n=1, mean=mle$est[1:N], sigma=covar))),reps)
inits <- rep(list(mle$est[1:N]), reps)
td <- 8
iter <- 300000
warmup <- iter/2
hh <- 3                           # hours to run
eps <- NULL
mm <- NULL #diag(length(par.names))
sfInit(parallel=TRUE, cpus=reps)
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
fit.rwm$warmup <- fit.rwm$warmup/tt
saveRDS(fit.rwm, file='fit.rwm.RDS')
stats.rwm <- data.frame(rstan::monitor(fit.rwm$samples, warmup=warmup/tt, probs=.5, print=FALSE))
perf.rwm <- min(stats.rwm$n_eff)/sum(fit.rwm$time.total)

dimnames(fit.rwm$samples)[[3]] <- paste0(1:53, dimnames(fit.rwm$samples)[[3]])

## How much faster is NUTS?
perf.nuts/perf.rwm

launch_shinystan_admb(fit.nuts)
launch_shinystan_admb(fit.rwm)
