### Template file to run a single model. Copy this over and modify below
### for the specific model (update data, inits, etc.).

## Then sourcing this file will run everything for this model, given the
## MCMC arguments are in the global workspace.
setwd(paste0('models/',m))

## Setup data, inits and pars
logLinf.mean <- log(50)
logk.mean <- log(.1)
t0 <- 6
logLinf.sigma <- .1
logk.sigma <- .2
Ntime <- 40
sigma.obs <- .1
Nfish <- 15
set.seed(115)
dat <- sample.lengths(Nfish=Nfish, n.ages=5, logLinf.mean=logLinf.mean,
                           logLinf.sigma=logLinf.sigma, logk.mean=logk.mean,
                           logk.sigma=logk.sigma, sigma.obs=sigma.obs, t0=t0)
g <- ggplot(dat, aes(ages, loglengths, group=fish, color=fish)) +
    geom_point(alpha=.5) + geom_line()
ggsave('plots/simulated_growth.png', g, width=9, height=5)
data <- list(Nfish=Nfish, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
inits <- list(list(sigma_obs=sigma.obs, logLinf_mean=logLinf.mean+.1,
                   logk_mean=logk.mean+.1, logLinf_sigma=logLinf.sigma,
                   logk_sigma=logk.sigma, logLinf=rep(logLinf.mean, len=Nfish),
                   logk=rep(logk.mean, len=Nfish) ))
pars <-
    c("logLinf_mean", "logLinf_sigma", "logk_mean", "logk_sigma", "logk", "logLinf",
      "sigma_obs", "delta")
par.names <- names(unlist(inits[[1]]))
Npar <- length(unlist(inits[[1]]))

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
## Setup the bounds
lower <- rep(-Inf, length(par.names))
names(lower) <- par.names
upper <- -1*lower
lower['sigma_obs'] <- 0
lower['logLinf_sigma'] <- 0
lower['logk_sigma'] <- 0
upper['sigma_obs'] <- 5
upper['logLinf_sigma'] <- 5
upper['logk_sigma'] <- 5
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
                Nthin=Nthin.ind, covar=covar.est, upper=upper, lower=lower)
## Load initial values from those sampled above.
sims.ind <- readRDS(file='sims.ind.RDS')
sims.ind <- sims.ind[sample(x=1:NROW(sims.ind), size=length(seeds)),]
## Setup inits for efficiency tests (fit.empirical). You need to adjust the
## named arguments for each model.
inits <- lapply(1:length(seeds), function(i)
    list(
        logLinf_mean=sims.ind$logLinf_mean[i],
        logLinf_sigma=sims.ind$logLinf_sigma[i],
        logk_mean=sims.ind$logk_mean[i],
        logk_sigma=sims.ind$logk_sigma[i],
        sigma_obs=sims.ind$sigma_obs[i],
        delta=sims.ind$delta[i],
        logk=as.numeric(sims.ind[i, grep('logk\\.', x=names(sims.ind))]),
        logLinf=as.numeric(sims.ind[i, grep('logLinf\\.', x=names(sims.ind))])))


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
dyn.unload(m)
message(paste('Finished with model:', m))
setwd('../..')
