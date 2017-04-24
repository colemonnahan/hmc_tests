## Run pilot chains for individual models using RWM, and then the "fixed"
## version of the model
## devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)
library(r4ss)

reps <- 4 # chains to run in parallel

sfStop()
d <- m <- 'cod'
d <- m <- 'cod2'
thin <- 1000
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <-
  sample_admb(m, iter=iter*thin, init=inits, thin=thin, mceval=TRUE,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_45")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


d <- m <- 'hake'
thin <- 1000
iter <- 1000
warmup <- iter/4
inits <- NULL
sfStop()
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits,  thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_2013", "Bratio_2013")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=F, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))


sfStop()
d <- m <- 'tanner'
thin <- 1
iter <- 1000
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
             parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
               parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/pilot_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'halibut'
d <- m <- 'halibut2'
thin <- 1000
iter <- 1000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "SPB_45", "Bratio_45")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq


sfStop()
d <- 'snowcrab'; m <- 'snowcrab'
d <- 'snowcrab2'; m <- 'snowcrab2'
thin <- 1000
iter <- 1000
warmup <- iter/4
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits,  thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))

sfStop()
d <- m <- 'canary2'
d <- m <- 'canary'
thin <- 1000
iter <- 1000
warmup <- iter/4
mle <- read_mle_fit(paste0(d,'/',m))
N <- mle$nopar
inits <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin, mceval=TRUE,
              dir=d, cores=reps, algorithm='RWM')
## Get posterior draws of dqs to cbind onto parameter draws later
dq.names <- c("SSB_MSY", "OFLCatch_2015", "Bratio_2015")
fit.rwm$dq.post <- r4ss::SSgetMCMC(dir=m)[[1]][,dq.names]
xx <- SS_output(m, model=m, verbose=TRUE, covar=T, ncols=500)
## Get estimates for derived quantitiesd
dq <- subset(xx$derived_quants, LABEL %in% dq.names)[,1:3]
names(dq) <- c('dq','mle', 'se'); rownames(dq) <- NULL
fit.rwm$dq <- dq
saveRDS(fit.rwm, file=paste0("results/pilot_rwm_", m, ".RDS"))
launch_shinyadmb(fit.rwm)

