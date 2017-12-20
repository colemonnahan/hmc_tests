## Script to demonstrate adnuts and tmbstan functions

## Working directory needs to be set to the folder of this file before
## proceeding.
source("startup.R")
library(snowfall)
cores <- parallel::detectCores()-1 # parallel cores
set.seed(145)
seeds <- sample(1:1e4, size=cores)

## Demonstrate ADMB with the swallows model
data <- swallows_setup()$data
inits.fn <- swallows_setup()$inits
## init can be a list of lists, or a function which returns a list
inits <- lapply(1:3, function(x) inits.fn())

setwd('models/swallows')
fit <- sample_admb(model='swallows', path='admb', init=inits,
                     seeds=seeds, parallel=TRUE, cores=cores)
## Extract the NUTS info, per iteration
sp <- extract_sampler_params(fit, inc_warmup=FALSE)
str(sp)
## Check for post-warmup divergences (also shown on console)
tapply(sp$divergent__, sp$chain, sum)

## Rerun model with higher target acceptance rate to see if divergences
## disappear
fit <- sample_admb(model='swallows', path='admb', init=inits,
                   seeds=seeds, parallel=TRUE, cores=cores,
                   control=list(adapt_delta=.9))
## Now the divergences are gone.
sum(extract_sampler_params(fit)$divergent__)
## And ESS and Rhats are good too
mon <- rstan::monitor(admb2$samples, print=FALSE)
min(mon[,'n_eff'])
max(mon[,'Rhat'])

## Can extract parameters as a data.frame or list object
post <- extract_samples(fit)
quantile(post[,1], c(0.1, 0.5, 0.9))
## The list can be converted to a CODA mcmc.list object for use with that
## package
library(coda)
post <- extract_samples(fit, as.list=TRUE)
postlist <- mcmc.list(lapply(post, mcmc))
par(mfrow=c(4,4))
coda::traceplot(postlist)
## Or shinystan can be used
launch_shinyadmb(fit)
setwd('../..')


### Demonstrate tmbstan with wildf (wildflower) model.
data <- wildf_setup()$data
inits.fn <- wildf_setup()$inits
setwd('models/wildf')
compile('wildf.cpp')
dyn.load(dynload('wildf'))
random <- c('yearInterceptEffect_raw', 'plantInterceptEffect_raw',
            'plantSlopeEffect_raw')
obj <- MakeADFun(data=data, parameters=inits.fn(), random=random,
            DLL='wildf')
set.seed(325234)
inits <- lapply(1:3, function(i) inits.fn())

## Run with Stan defaults, including in parallel
rstan_options(auto_write = TRUE)
options(mc.cores=cores)
tmb1 <- tmbstan(obj=obj, chains=3, init=inits)

## Methods provided by 'rstan'
class(mcmc.tmb)
methods(class="stanfit")

## Pairs plot
pairs(mcmc.tmb, pars=names(obj$par))
post <- as.data.frame(mcmc.tmb) # get data.frame

## Trace plot
rstan::traceplot(mcmc.tmb, pars=names(obj$par), inc_warmup=TRUE)


