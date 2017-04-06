## Run long chains for individual models using RWM.
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

reps <- 6 # chains to run in parallel
hh <- 24 # hours to run

sfStop()
d <- m <- 'cod_fast'
thin <- 1
iter <- 2000*thin
warmup <- 1000*thin
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.long <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.long, file=paste0("results/long_run_", m, ".RDS")
## launch_shinyadmb(fit.long)


sfStop()
d <- m <- 'hake'
thin <- 1
iter <- 2000*thin
warmup <- 1000*thin
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.hake <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.hake, file=paste0("results/long_run_", m, ".RDS")
 launch_shinyadmb(fit.hake)
