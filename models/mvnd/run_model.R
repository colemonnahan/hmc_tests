### Template file to run a single model. Copy this over and modify below
### for the specific model (update data, inits, etc.).

## Then sourcing this file will run everything for this model, given the
## MCMC arguments are in the global workspace.
setwd(paste0('models/',m))

## Setup data, inits and pars
Npar <- 5
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
covar <- diag(Npar)
data <- list(Npar=Npar, covar=covar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
pars <- 'mu'


## Compile Stan, TMB and ADMB models
obj.stan <- stan(file= paste0(m, '.stan'), data=data, iter=5000,
                   chains=1, init=list(inits[[1]]),
                   control=list(adapt_engaged=TRUE))
## Use these samples to get an estimated covariance to use in TMB and ADMB
## since no adapation of M yet.
samples <- extract(obj.stan, permuted=FALSE)[,1,1:Npar]
covar.est <- cov(samples)               # estimated mass matrix

#covar[1,1] <- covar[1,1]*10
data <- list(Npar=Npar, covar=covar, x=rep(0, len=Npar))
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj.tmb <- MakeADFun(data=data, parameters=inits[[1]], DLL=m)
devtools::install("c:/users/cole/adnuts")
x1 <- sample_tmb(obj.tmb, iter=2000, init=inits[[1]], control=list(metric=diag(5)))
x2 <- sample_tmb(obj.tmb, iter=2000, init=inits[[1]], control=list(metric=covar))

launch_shinystan_tmb(x1)
launch_shinystan_tmb(x2)

setwd('admb')
write.table(x=unlist(data), file=paste0(m,'.dat'), row.names=FALSE,
            col.names=FALSE )
system(paste('admb',m))
system(m)
setwd('..')

## Get independent samples from each model to ensure identical
## posteriors. For now I am using covar.est since do not care about fair
## comparisons, and it should make TMB and ADMB run faster.
if(verify)
  verify.models(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m, dir='admb',
                pars=pars, inits=inits, data=data, Nout=Nout.ind,
                Nthin=Nthin.ind, covar=covar.est)
## Load initial values from those sampled above.
sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
## Setup inits for efficiency tests (fit.empirical). You need to adjust the
## named arguments for each model.
inits <- lapply(1:length(seeds), function(i) list(mu=as.numeric(sims.ind[i,])))

## Fit empirical data with no thinning for efficiency tests. Since mass
## matrix adapataion is not working for TMB/ADMB, approximate it with the
## diagonal of that estimated from Stan.
covar2 <- diag(x=diag(covar.est))
fit.empirical(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m, pars=pars, inits=inits, data=data,
              delta=delta, metric='unit', seeds=seeds, covar=covar2,
              Nout=Nout, max_treedepth=12)

## If there is a simulation component put it in this file
if(TRUE)
  source("simulation.R")

rm(obj.stan, obj.tmb, data, inits, pars, lower, upper, covar.est)
dyn.unload(m)
message(paste('Finished with model:', m))
setwd('../..')
