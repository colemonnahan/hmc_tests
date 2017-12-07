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
admb1 <- sample_admb(model='swallows', path='admb', init=inits, seeds=1:3,
                         parallel=TRUE, cores=3)
saveRDS(admb1, file='admb1.RDS')

admb1 <- readRDS("models/swallows/admb1.RDS")
mon <- rstan::monitor(admb1$samples, print=FALSE)
mon[1:4,'n_eff']
mon[1:4,'Rhat']
min(mon[,'n_eff'])
max(mon[,'Rhat'])


mon <- as.data.frame(rstan::monitor(fit$samples, print=FALSE))
mon$pars <- row.names(mon)
mon <- mon[order(mon$n_eff, decreasing=FALSE),]
post <- extract_samples(fit)
sp <- extract_sampler_params(fit)
str(post[,1:5])
str(sp)


slow <-  c("sigmayearphi", "yeareffphi_raw[3]", "yeareffphi_raw[2]",
           "yeareffphi_raw[4]", "yeareffphi_raw[1]")
png('slow_mixing.png', width=5, height=3.5, units='in', res=300)
pairs_admb(fit, pars=slow)
dev.off()





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

set.seed(23523)
admb2 <- sample_admb(model='swallows', path='admb', seeds=1:3, init=inits,
                         parallel=TRUE, cores=3, control=list(adapt_delta=.9))
## Now the divergences are gone.
sum(extract_sampler_params(admb2)$divergent__)
## And ESS and Rhats are good too
mon <- rstan::monitor(admb2$samples, print=FALSE)
min(mon[,'n_eff'])
max(mon[,'Rhat'])



setwd('../..')

cov.est <- admb2$covar.est
admb3 <- sample_admb(model='swallows', path='admb', seeds=1:3, init=inits,
                     parallel=TRUE, cores=3,
                     control=list(adapt_delta=.8, metric=cov.est))

## admb4 <- sample_admb(model='swallows', path='admb', seeds=seeds,
##                      init=inits, iter=200000, thin=100,
##                      parallel=TRUE, chains=3, cores=3, algorithm='RWM')

## admb5 <- sample_admb(model='swallows', path='admb', seeds=seeds,
##                      init=inits, iter=200000, thin=100,
##                      parallel=TRUE, chains=3, cores=3, algorithm='RWM',
##                      control=list(metric=admb4$covar.est))




### Demonstrate tmbstan
data <- wildf_setup()$data
inits <- wildf_setup()$inits
setwd('models/wildf')
compile('wildf.cpp')
dyn.load('wildf')
random <- c('yearInterceptEffect_raw', 'plantInterceptEffect_raw',
            'plantSlopeEffect_raw')
obj <- MakeADFun(data=data, parameters=inits(), random=random,
            DLL='wildf')

set.seed(325234)
inits <- lapply(1:3, function(i) inits.fn())
tmb1 <- tmbstan(obj=obj, init=inits)

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
