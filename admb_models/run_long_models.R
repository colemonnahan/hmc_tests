## Run long chains for individual models using RWM.
setwd('admb_models/')
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

reps <- 6 # chains to run in parallel
hh <- 24 # hours to run

sfStop()
d <- m <- 'cod_fast'
thin <- 100
iter <- 2000*thin
warmup <- 1000*thin
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=5, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()

fit.long <-
  sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='RWM')
launch_shinyadmb(fit.long)
