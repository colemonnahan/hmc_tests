## Quick code to test TMB bounding
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)

## Run a cross of bounded/unbounded and covar/no covar to make sure they
## all work
## 2d example
covar <- matrix(.954, nrow=2, ncol=2)
diag(covar) <- 1
covar.inv <- solve(covar)
Npar <- 2
Niter <- 500

dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')


source("../adcomp/TMB/R/mcmc.R")

lower <- c(-20,-20); upper <- c(20,20)
init <- lapply(1:1, function(x) rnorm(2, sd=1))

fit <- run_mcmc(obj, iter=2000, init=NULL, chains=3, covar=covar)
## samples from all chains, excluding warmup samples
x <- extract_samples(fit)
## Rhat, ESS, and quantiles
diagnostics <- rstan::monitor(sims=fit$samples)
## Quick way to get reported values for each row
derived.report <- unlist(apply(x, 1, function(i) as.numeric(mvnd.obj$report(par=i))))
## A shiny tool to look at NUTS output. Requires installing dev version of
## shinystan with:
## devtools::install_github('stan-dev/shinystan', ref='feature-issue-136')
launch_shinystan_tmb(fit)

