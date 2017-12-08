## Script to do MLE fits for the example models

## Working directory needs to be set to the folder of this file before
## proceeding.

## Note: these models include priors, and use exponentiation in the
## template with a Jacobian adjustment instead of external bounds of
## (0,Inf) for hypervariance parameters. Thus, the SD parameters are really
## in log space. See template for more details. Should we turn the priors
## off for MLE?

source("startup.R")

## swallows model
data <- swallows_setup()$data
inits <- swallows_setup()$inits
compile('models/swallows/swallows.cpp')
dyn.load('models/swallows/swallows')
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'),
            DLL='swallows')
fit.swallows <- with(obj, nlminb(par, fn, gr))
mcmc.swallows <- tmbstan(obj, iter=2000, chains=3)

## wildf model
data <- wildf_se tup()$data
inits <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'),
            DLL='wildf')
fit.wildf <- with(obj, nlminb(par, fn, gr))
mcmc.wildf <- tmbstan(obj, iter=2000, chains=1)

## Quickly rerun to plot MCMC vs MLE for fixed effects
mcmc.wildf2 <- sample_tmb(obj, iter=5000, warmup=1000, init=inits, chains=3, cores=3,
                          parallel=TRUE, path='models/wildf', control=list(adapt_delta=.59)
fit.wildf <- with(obj, nlminb(par, fn, gr))
## Have to manually modify the list so that pairs_admb can use it. It's not
## designed for tmb output.
sd <- sdreport(obj)
se <- sqrt(diag(sd$cov.fixed))
cor <- sd$cov.fixed / se %o% se
mle <- list(par.names=dimnames(mcmc.wildf2$samples)[[3]][1:9],
            est=fit.wildf$par, se=se, cor=cor)
mcmc.wildf2$mle <- mle
png('pairs_wildf.png', width=6, height=4, units='in', res=500)
pairs_admb(mcmc.wildf2, pars=mcmc.wildf2$mle$par.names)
dev.off()
