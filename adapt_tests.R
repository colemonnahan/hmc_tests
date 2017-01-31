library(TMB)
library(plyr)
library(coda)
## devtools::install_github('bbolker/R2admb/R2admb')
library(R2admb)
library(rstan)
library(shinystan)
devtools::load_all("c:/Users/Cole/adnuts")

## Test that step size and mass matrix adaptation is consistent across TMB
## and ADMB Start with simple MVN model.
set.seed(235)
Npar <- 2
corr <- matrix(.98, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
covar2 <- diag(x=se^2)
covar.inv <- solve(covar)
samples <- mvtnorm::rmvnorm(n=1e5, sigma=covar)
apply(samples, 2, var)
setwd('C:/Users/Cole/hmc_tests/')
dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')
model.path <- 'models/mvnd'
model.name <- 'mvnd'
write.table(x=c(2, covar), file='models/mvnd/mvnd.dat', row.names=FALSE, col.names=FALSE)
setwd('models/mvnd')
system('admb mvnd')
system('mvnd')
setwd('../..')

chains <- 5
init <- lapply(1:chains, function(x) rnorm(2, sd=se*1))
iter <- 1000
tmb1 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=diag(2))
tmb2 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar2)
tmb3 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar)
adapt.tmb1 <- cbind(form=1, extract_adapt(tmb1))
adapt.tmb2 <- cbind(form=2, extract_adapt(tmb2))
adapt.tmb3 <- cbind(form=3, extract_adapt(tmb3))
adapt.tmb <- rbind(adapt.tmb1, adapt.tmb2, adapt.tmb3)
admb1 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=diag(2))
admb2 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=covar2)
admb3 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=covar)
adapt.admb1 <- cbind(form=1, extract_adapt(admb1))
adapt.admb2 <- cbind(form=2, extract_adapt(admb2))
adapt.admb3 <- cbind(form=3, extract_adapt(admb3))
adapt.admb <- rbind(adapt.admb1, adapt.admb2, adapt.admb3)

adapt <- rbind(cbind(name='tmb', adapt.tmb),
               cbind(name='admb', adapt.admb))
ggplot(adapt, aes(form, stepsize, color=name)) + geom_jitter()
ggplot(adapt, aes(form, nsteps.median, color=name)) + geom_jitter()
ggplot(adapt, aes(form, accept.prob, color=name)) + geom_jitter()
ggplot(adapt, aes(form, log(stepsize0), color=name)) + geom_jitter()
## launch_shinystan_tmb(tmb1)
## launch_shinystan_tmb(tmb2)
## launch_shinystan_tmb(tmb3)

extract_adapt <- function(fit){
  ldply(1:length(fit$sampler_params), function(i){
    ind <- -(1:fit$warmup)
    xx <- fit$sampler_params[[i]]
    cbind(chain=i, stepsize=tail(xx[ind,2],1), nsteps.median= median(xx[ind,4]),
          nsteps.mean=mean(xx[ind,4]),
          accept.prob=mean(xx[ind,1]),
          stepsize0=head(xx[ind,2],1))})
}



model.path="C:/Users/Cole/hmc_tests/models/catage"
model.name='catage'
x <- run_admb_mcmc(model.path=model.path, model.name=model.name, iter=20000,
                   chains=1, eps=.2, max_treedepth=10, thin=100)
launch_shinystan_admb(x)
setwd(model.path)
system('admb catage')
system('catage -nohess -mcmc 10 -nuts -mcseed 5')
adapt <- read.csv("adaptation.csv")
pars <- read_psv(model.name)
fit <- read_admb(model.name)
