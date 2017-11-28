## Script to do MLE fits for the example models

## Working directory needs to be set to the folder of this file before
## proceeding.

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
data <- wildf_setup()$data
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
mcmc.wildf <- tmbstan(obj, iter=2000, chains=3)
