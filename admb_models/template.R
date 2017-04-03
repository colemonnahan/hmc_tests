library(ggplot2)
library(adnuts)
library(snowfall)
library(shinystan)

## Rerun model
setwd(d)
system(paste('admb', m))
system(m)
setwd('..')

sfStop()
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- mle$names[1:N]
par.names <- paste0(1:N, "_", par.names)
reps <- 6                        # chains/reps to run
## Draw inits from MVN using MLE and covar
inits <- NULL#rep(list(mle$est[1:N]), reps)
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) as.vector(mvtnorm::rmvnorm(n=1, mean=mle$est[1:N], sigma=covar)))
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, sigma=covar)))
td <- 12
iter <- 2000
warmup <- iter/2
tt <- 100 # thin rate
hh <- 10                           # hours to run
eps <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()

## Run NUTS for different mass matrices
fit.nuts.unit <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric='unit', adapt_delta=.9))
fit.nuts.mle <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=NULL, adapt_delta=.9))
covar.diag <- diag(x=diag(fit.nuts.mle$covar.est))
fit.nuts.diag <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=covar.diag, adapt_delta=.9))
covar.dense <- fit.nuts.mle$covar.est
fit.nuts.dense <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=covar.dense, adapt_delta=.9))

## Run RWM for different mass matrices
fit.rwm.unit <-
  sample_admb(m, iter=tt*iter, init=inits, par.names=par.names, thin=tt,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric='unit'),
              algorithm='RWM')
fit.rwm.mle <-
  sample_admb(m, iter=tt*iter, init=inits, par.names=par.names, thin=tt,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=NULL),
              algorithm='RWM')
covar.diag <- diag(x=diag(fit.rwm.mle$covar.est))
fit.rwm.diag <-
  sample_admb(m, iter=tt*iter, init=inits, par.names=par.names, thin=tt,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=covar.diag),
              algorithm='RWM')
covar.dense <- fit.rwm.mle$covar.est
fit.rwm.dense <-
  sample_admb(m, iter=tt*iter, init=inits, par.names=par.names, thin=tt,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=covar.dense),
              algorithm='RWM')

## Gather adaptation and performance metrics
ff <- function(labels, ...){
  fits <- list(...)
  if(length(fits)!=length(labels)) stop("aa")
  do.call(rbind,  lapply(1:length(fits), function(i){
    x <- fits[[i]]$sampler_params
    data.frame(m=labels[i], chain=1:length(x),
               eps=as.numeric(do.call(rbind, lapply(x, function(l) tail(l[,2],1)))))
  }))
}
adaptation <- ff(c("unit", "mle", "diag", "dense"), fit.nuts.unit, fit.nuts.mle,
  fit.nuts.diag, fit.nuts.dense)
ggplot(adaptation, aes(y=log10(eps), x=m)) + geom_point(alpha=.5)

stats.nuts.unit <- with(fit.nuts.unit, data.frame(alg='nuts', m='unit', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.unit <- data.frame(alg='nuts', m='unit', efficiency=min(stats.nuts.unit$n_eff)/sum(fit.nuts.unit$time.total))
stats.nuts.mle <- with(fit.nuts.mle, data.frame(alg='nuts', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.mle <- data.frame(alg='nuts', m='mle', efficiency=min(stats.nuts.mle$n_eff)/sum(fit.nuts.mle$time.total))
stats.nuts.diag <- with(fit.nuts.diag, data.frame(alg='nuts', m='diag', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.diag <- data.frame(alg='nuts', m='diag', efficiency=min(stats.nuts.diag$n_eff)/sum(fit.nuts.diag$time.total))
stats.nuts.dense <- with(fit.nuts.dense, data.frame(alg='nuts', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.dense <- data.frame(alg='nuts', m='dense', efficiency=min(stats.nuts.dense$n_eff)/sum(fit.nuts.dense$time.total))
stats.rwm.unit <- with(fit.rwm.unit, data.frame(alg='rwm', m='unit', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.unit <- data.frame(alg='rwm', m='unit', efficiency=min(stats.rwm.unit$n_eff)/sum(fit.rwm.unit$time.total))
stats.rwm.mle <- with(fit.rwm.mle, data.frame(alg='rwm', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.mle <- data.frame(alg='rwm', m='mle', efficiency=min(stats.rwm.mle$n_eff)/sum(fit.rwm.mle$time.total))
stats.rwm.diag <- with(fit.rwm.diag, data.frame(alg='rwm', m='diag', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.diag <- data.frame(alg='rwm', m='diag', efficiency=min(stats.rwm.diag$n_eff)/sum(fit.rwm.diag$time.total))
stats.rwm.dense <- with(fit.rwm.dense, data.frame(alg='rwm', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.dense <- data.frame(alg='rwm', m='dense', efficiency=min(stats.rwm.dense$n_eff)/sum(fit.rwm.dense$time.total))

stats.all <- rbind(stats.nuts.unit, stats.nuts.mle, stats.nuts.diag,
                   stats.nuts.dense, stats.rwm.unit, stats.rwm.mle,
                   stats.rwm.diag, stats.rwm.dense)
perf.all <- rbind(perf.nuts.unit, perf.nuts.mle, perf.nuts.diag,
                   perf.nuts.dense, perf.rwm.unit, perf.rwm.mle,
                   perf.rwm.diag, perf.rwm.dense)

ggplot(perf.all, aes(m, log10(efficiency), color=alg)) + geom_point()
ggplot(stats.all, aes(m, Rhat, color=alg)) + geom_point()
ggplot(stats.all, aes(m, n_eff, color=alg)) + geom_point()
ggplot(stats.all, aes(m, time.total, color=alg)) + geom_point()
