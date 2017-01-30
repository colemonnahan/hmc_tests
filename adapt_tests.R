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
corr <- matrix(-.954, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, 5)
covar <- corr * (se %o% se)
covar2 <- diag(2); diag(covar2) <- se^2
covar2[1,2] <- covar2[2,1] <- 4
covar.inv <- solve(covar)
samples <- mvtnorm::rmvnorm(n=1e5, sigma=covar)
apply(samples, 2, var)
setwd('C:/Users/Cole/hmc_tests/')
dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')

chains <- 5
init <- lapply(1:chains, function(x) rnorm(2, sd=se*10))
iter <- 2000
tmb1 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=NULL)
tmb2 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar2)
tmb3 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar)
eps.tmb1 <- laply(tmb1$sampler_params, function(x) tail(x[,2],1))
eps.tmb2 <- laply(tmb2$sampler_params, function(x) tail(x[,2],1))
eps.tmb3 <- laply(tmb3$sampler_params, function(x) tail(x[,2],1))
boxplot(eps.tmb1, eps.tmb2, eps.tmb3)
## launch_shinystan_tmb(tmb1)
## launch_shinystan_tmb(tmb2)
## launch_shinystan_tmb(tmb3)






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
