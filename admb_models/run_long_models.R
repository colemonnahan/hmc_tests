## Run long chains for individual models using RWM.
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

reps <- 4 # chains to run in parallel
hh <- 24 # hours to run

sfStop()
d <- m <- 'cod_fast'
thin <- 100
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM', extra.args=' -ainp cod_fast.par')

saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=3, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS',
              control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'hake'
thin <- 1000
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'tanner'
thin <- 100
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'halibut'
thin <- 100
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))
