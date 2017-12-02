## Script to demonstrate adnuts and tmbstan functions

## Working directory needs to be set to the folder of this file before
## proceeding.
source("startup.R")
library(snowfall)
library(coda)
set.seed(121413) # ensures random inits are the same

## swallows model
data <- swallows_setup()$data
inits <- swallows_setup()$inits
setwd('models/swallows')
seeds <- 1:3
admb1 <- sample_admb(model='swallows', path='admb', seeds=seeds, init=inits,
                         iter=1000, parallel=TRUE, cores=3)
saveRDS(admb1, file='admb1.RDS')

admb1 <- readRDS("models/swallows/admb1.RDS")

## Can extract parameters as a data.frame or list object
post <- extract_samples(admb1, as.list=TRUE)
## The list can be converted to a CODA mcmc.list object for use with that
## package
postlist <- mcmc.list(lapply(post, mcmc))
par(mfrow=c(3,3))
coda::traceplot(postlist)
## Or shinystan can be used
launch_shinyadmb(admb1)
## Extract the NUTS metadata
sp <- extract_sampler_params(admb1, inc_warmup=TRUE)
sum(sp$divergent__)

seeds <- 1:3
admb2 <- sample_admb(model='swallows', path='admb', seeds=seeds, init=inits,
                         parallel=TRUE, cores=3, control=list(adapt_delta=.9))
## Now the divergences are gone.
sum(extract_sampler_params(admb2)$divergent__)
setwd('../..')

cov.est <- admb1$covar.est
admb3 <- sample_admb(model='swallows', path='admb', seeds=seeds, init=inits,
                     parallel=TRUE, cores=3,
                     control=list(adapt_delta=.9, metric=cov.est))

## admb4 <- sample_admb(model='swallows', path='admb', seeds=seeds,
##                      init=inits, iter=200000, thin=100,
##                      parallel=TRUE, chains=3, cores=3, algorithm='RWM')

## admb5 <- sample_admb(model='swallows', path='admb', seeds=seeds,
##                      init=inits, iter=200000, thin=100,
##                      parallel=TRUE, chains=3, cores=3, algorithm='RWM',
##                      control=list(metric=admb4$covar.est))




### Demonstrate tmbstan
data <- wildf_setup()$data
inits.fn <- wildf_setup()$inits
setwd('models/wildf')
compile('wildf.cpp')
dyn.load('wildf')
obj <- MakeADFun(data=data, parameters=inits.fn(), DLL='wildf')

set.seed(325234)
inits <- lapply(1:3, function(i) inits.fn())
mcmc.tmb <- tmbstan(obj=obj, seed=1, iter=2000, chains=3)

## Methods provided by 'rstan'
class(mcmc.tmb)
methods(class="stanfit")

## Pairs plot
pairs(mcmc.tmb, pars=names(obj$par))
post <- as.data.frame(mcmc.tmb) # get data.frame

## Trace plot
rstan::traceplot(mcmc.tmb, pars=names(obj$par), inc_warmup=TRUE)

## Run in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
