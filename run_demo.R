## Script to demonstrate adnuts and tmbstan functions

## Working directory needs to be set to the folder of this file before
## proceeding.
source("startup.R")
library(snowfall)
cores <- 3
set.seed(1459)
seeds <- sample(1:1e4, size=cores)

## Demonstrate ADMB with the swallows model
data <- swallows_setup()$data
inits.fn <- swallows_setup()$inits
## init can be a list of lists, or a function
inits <- lapply(1:cores, function(x) inits.fn())

#setwd('models/swallows')
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
mon <- rstan::monitor(fit$samples, print=FALSE)
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
par(mfrow=c(5,5))
coda::traceplot(postlist)
## Or shinystan can be used
launch_shinyadmb(fit)
setwd('../..')


### Demonstrate tmbstan with wildf (wildflower) model.
library(TMB)
library(tmbstan)
detach("package:snowfall", unload=TRUE)
detach("package:snow", unload=TRUE)
data <- wildf_setup()$data
inits.fn <- wildf_setup()$inits
set.seed(325234)
inits <- lapply(1:cores, function(i) inits.fn())
setwd('models/wildf')
compile('wildf.cpp')
dyn.load(dynlib('wildf'))
random <- c('yearInterceptEffect_raw', 'plantInterceptEffect_raw',
            'plantSlopeEffect_raw')
obj <- MakeADFun(data=data, parameters=inits[[1]], random=random,
            DLL='wildf')

## Run with Stan defaults, including in parallel
rstan_options(auto_write = TRUE)
options(mc.cores=cores)
fit <- tmbstan(obj=obj, chains=cores, init=inits, seed=1934)
## Equivalent Stan code. Note that even with same inits and seeds small
## rounding errors accumulate and the chains may not be identical,
## particularly for complex models.
## fit <- stan(file='wildf.stan', data=data, chains=cores,init=inits,
##             seed=1934)

## Methods provided by 'rstan'
class(fit)
methods(class="stanfit")

## Pairs plot
pairs(fit, pars=names(obj$par))
post <- as.data.frame(fit) # get data.frame

## Trace plot
rstan::traceplot(fit, pars=names(obj$par), inc_warmup=TRUE)

## Now refit with the Laplace approximation turned on
fit.la <- tmbstan(obj=obj, chains=cores, init=inits, seed=1934, laplace=TRUE)
## Notice that Stan only sampled from the fixed effects:
str(as.data.frame(fit.la))

## It's mixing better
(m1 <- min(monitor(fit, print=FALSE)[,'n_eff']))
(m2 <- min(monitor(fit.la, print=FALSE)[,'n_eff']))
## But takes much longer to run
(t1 <- get_elapsed_time(fit))
(t2 <- get_elapsed_time(fit.la))
