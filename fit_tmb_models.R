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
dyn.load(dynlib('models/swallows/swallows'))
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'),
            DLL='swallows')
fit.swallows <- with(obj, nlminb(par, fn, gr))
mcmc.swallows <- tmbstan(obj, iter=2000, chains=3, seed=1)

## Compare
layout(matrix(1:30,5))
s <- as.array(mcmc.swallows)
for(i in 1:length(obj$par)) {
    plot(density(s[,,i]))
    abline(v=fit.swallows$par[i])
}

## wildf model
data <- wildf_setup()$data
data$flag_prior <- 0
data$flag_jac <- 0
inits <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load(dynlib('models/wildf/wildf'))
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'),
            DLL='wildf')
fit.wildf <- with(obj, nlminb(par, fn, gr))
mcmc.wildf <- tmbstan(obj, iter=2000, chains=3, seed=1)
## Compare
layout(matrix(1:9,3))
s <- as.array(mcmc.wildf)
for(i in 1:length(obj$par)) {
    plot(density(s[,,i]), main=names(obj$par)[i])
    abline(v=fit.wildf$par[i])
}

## Try Laplace checker:
set.seed(1)
chk <- checkConsistency(obj, n = 1000)
s <- summary(chk)
s$marginal$bias
sdr <- sdreport(obj)
s$marginal$bias / summary(sdr,"fixed")[,2]

## Understand MCMC problem (Posterior moments all infinite due to nll plateau):
lpb <- obj$env$last.par.best ## Best encountered LA parameter (fixed + random)
f <- Vectorize(function(x) {
    ## Slice: "plantSlopeSD"
    lpb[3] <- lpb[3] + x
    obj$env$f(lpb)
})
layout(1)
plot(f,-20,2)
