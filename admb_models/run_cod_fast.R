library(adnuts)
library(snowfall)
library(shinystan)
m <- 'cod_fast'

## ## Run model if needed
## setwd(m)
## system('admb ss3')
## system("ss3")
## replist <- r4ss::SS_output(getwd(), covar=TRUE)
## unlink('plots', TRUE)
## r4ss::SS_plots(replist)
## setwd('..')
## temp <- get.admb.cov(m)

sfStop()
mle <- r4ss::read.admbFit('cod_fast/ss3')
covar <- get.admb.cov(m)$cov.unbounded
N <- mle$nopar
par.names <- mle$names[1:N]
par.names <- paste0(1:N, "_", par.names)
reps <- 4                        # chains/reps to run
## Draw inits from MVN using MLE and covar
#inits <- rep(list(as.vector(mvtnorm::rmvnorm(n=1, mean=mle$est[1:N], sigma=covar))),reps)
inits <- NULL#rep(list(mle$est[1:N]), reps)
td <- 10
iter <- 500
warmup <- 250
tt <- 1 # thin rate
hh <- 10                           # hours to run
eps <- NULL
mm <- NULL #diag(length(par.names))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()

fit.nuts <-
  sample_admb('ss3', iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=mm, adapt_delta=.9))
saveRDS(fit.nuts, file='fit.nuts.RDS')
stats.nuts <- data.frame(rstan::monitor(fit.nuts$samples, warmup=warmup, probs=.5, print=FALSE))
perf.nuts <- min(stats.nuts$n_eff)/sum(fit.nuts$time.total)
fit.rwm <-
  sample_admb('ss3', iter=iter*tt, duration=hh*60, init=inits, par.names=par.names,
              parallel=FALSE, chains=reps, warmup=warmup*tt, dir=m, cores=reps,
              thin=tt, control=list(algorithm='RWM', metric=mm))
fit.rwm$warmup <- fit.rwm$warmup/tt
saveRDS(fit.rwm, file='fit.rwm.RDS')
stats.rwm <- data.frame(rstan::monitor(fit.rwm$samples, warmup=warmup/tt, probs=.5, print=FALSE))
perf.rwm <- min(stats.rwm$n_eff)/sum(fit.rwm$time.total)

## How much faster is NUTS?
perf.nuts/perf.rwm

launch_shinystan_admb(fit.nuts)
launch_shinystan_admb(fit.rwm)

fit.rwm$samples[1, 1, N+1]
fit.nuts$samples[1, 1, N+1]
as.vector(fit.rwm$samples[1, 1, -(N+1)])-inits[[1]]
fit.nuts$samples[1, 1, N]

inits=rep(list(fit.rwm$samples[100,1,1:N]+.1),reps)
