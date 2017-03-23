### Template file to run a single model. Copy this over and modify below
### for the specific model (update data, inits, etc.).

## Then sourcing this file will run everything for this model, given the
## MCMC arguments are in the global workspace.
setwd(paste0('models/',m))

## Setup data, inits and pars
Npar <- 3
covar <- matrix(.5, 3,3)
diag(covar) <- 1
data <- list(covar=covar, x=rep(0, len=Npar))
inits <- list(list(mu1=1, mu2=1, mu3=1))
pars <- c('mu1', 'mu2', 'mu3')

## Compile Stan, TMB and ADMB models
obj.stan <- stan(file= paste0(m, '.stan'), data=data, iter=5000,
                   chains=1, init=list(inits[[1]]),
                   control=list(adapt_engaged=TRUE))
## Use these samples to get an estimated covariance to use in TMB and ADMB
## since no adapation of M yet.
samples <- extract(obj.stan, permuted=FALSE)[,1,1:Npar]
covar.est <- cov(samples)               # estimated mass matrix
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj.tmb <- MakeADFun(data=data, parameters=inits[[1]])
lower <- c(-2,0,-Inf)
upper <- c(2,Inf, Inf)
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
                Nthin=Nthin.ind, covar=covar.est, lower=lower, upper=upper,
                admb.columns=2)
## Load initial values from those sampled above.
sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
## Setup inits for efficiency tests (fit.empirical). You need to adjust the
## named arguments for each model.
inits <- lapply(1:length(seeds), function(i)
  list(mu1=sims.ind$mu1[i], mu2=sims.ind$mu2[i], mu3=sims.ind$mu3[i]))

## Fit empirical data with no thinning for efficiency tests. Since mass
## matrix adapataion is not working for TMB/ADMB, approximate it with the
## diagonal of that estimated from Stan.
covar2 <- diag(x=diag(covar.est))
fit.empirical(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m, pars=pars, inits=inits, data=data,
              delta=delta, metric='unit', seeds=seeds, covar=covar2,
              Nout=Nout, max_treedepth=12)

## If there is a simulation component put it in this file
if(FALSE)
  source("simulation.R")

rm(obj.stan, obj.tmb, data, inits, pars, lower, upper, covar.est)
clean.TMB.files()
message(paste('Finished with model:', m))
setwd('../..')
