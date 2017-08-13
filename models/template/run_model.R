### Template file to run a single model. Copy this over and modify below
### for the specific model (update data, inits, etc.).

## Then sourcing this file will run everything for this model, given the
## MCMC arguments are in the global workspace.
setwd(paste0('models/',m))

## Setup data, inits and pars
data <- NULL
inits <- NULL
pars <- NULL

## Compile Stan, TMB and ADMB models
obj.stan <- stan(file= paste0(m, '.stan'), data=data, iter=500,
                   chains=1, init=list(inits[[1]]),
                   control=list(adapt_engaged=TRUE))
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj.tmb <- MakeADFun(data=data, parameters=inits[[1]])
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

## Fit empirical data with no thinning for efficiency tests. Using Stan
## defaults for TMB and ADMB too: estimated diagonal mass matrix. I also
## dropped RWM since it wont work for mixed effects models
fit.empirical(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m, pars=pars, inits=inits, data=data,
              delta=delta, metric='diag', seeds=seeds, covar=covar2,
              Nout=Nout, max_treedepth=12)

## If there is a simulation component put it in this file
if(FALSE)
  source("simulation.R")

rm(obj.stan, obj.tmb, data, inits, pars, lower, upper, covar.est)
dyn.unload(m)
message(paste('Finished with model:', m))
setwd('../..')
