devtools::install('C:/Users/Cole/adnuts')
library(adnuts)
library(snowfall)
library(shinystan)
m <- 'catage'

## Run model if needed
setwd(m)
system('admb catage')
system("catage")
setwd('..')
temp <- get.admb.cov(m)

sfStop()
mle <- r4ss::read.admbFit('catage/catage')
N <- mle$nopar
par.names <- mle$names[1:N]
par.names <- paste0(1:N, "_", par.names)
reps <- 6                        # chains/reps to run
## Draw inits from MVN using MLE and covar
inits <- NULL#rep(list(mle$est[1:N]), reps)
covar <- get.admb.cov(m)$cov.unbounded
inits <- lapply(1:reps, function(i) as.vector(mvtnorm::rmvnorm(n=1, mean=mle$est[1:N], sigma=covar)))
td <- 10
iter <- 2000
warmup <- 1000
tt <- 1 # thin rate
hh <- NULL                           # hours to run
eps <- NULL
mm <- NULL #diag(length(par.names))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()

## Run it with MLE covar
fit.nuts.unit <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric='unit', adapt_delta=.8))
fit.nuts.mle <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=NULL, adapt_delta=.8))
covar.diag <- diag(x=diag(fit.nuts.mle$covar.est))
fit.nuts.diag <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=covar.diag, adapt_delta=.8))
covar.dense <- fit.nuts.mle$covar.est
fit.nuts.dense <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup, dir=m, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=covar.dense, adapt_delta=.95))

ff <- function(labels, ...){
  fits <- list(...)
  if(length(fits)!=length(labels)) stop("aa")
  do.call(rbind,  lapply(1:length(fits), function(i){
    x <- fits[[i]]$sampler_params
    data.frame(m=labels[i], chain=1:length(x), eps=as.numeric(do.call(rbind, lapply(x, function(l) tail(l[,2],1)))))
  }))
}
results <- ff(c("unit", "mle", "diag", "dense"), fit.nuts.unit, fit.nuts.mle,
  fit.nuts.diag, fit.nuts.dense)
ggplot(results, aes(y=eps, x=m)) + geom_violin()

launch_shinystan_admb(fit.nuts.dense)


stats.nuts <- data.frame(rstan::monitor(fit.nuts$samples, warmup=warmup, probs=.5, print=FALSE))
perf.nuts <- min(stats.nuts$n_eff)/sum(fit.nuts$time.total)
fit.rwm <-
  sample_admb(m, iter=iter*tt, duration=hh*60, init=inits, par.names=par.names,
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
