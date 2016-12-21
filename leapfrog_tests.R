library(plyr)
library(coda)
library(R2admb)
library(rstan)
library(shinystan)
library(TMB)

setwd("models/simple/")

#dyn.unload(dynlib('simple_tmb'))
compile(file='simple_tmb.cpp')
dyn.load(dynlib('simple_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')
