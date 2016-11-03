## Quick code to test TMB bounding
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)

#' The logistic transformation function for bounding parameters in MCMC
boundp <- function(x, minb, maxb) minb+(maxb-minb)/(1+exp(-x))
#' The inverse of the transformation
boundpin <- function(y, minb, maxb) -log( (maxb-y)/(y-minb) )
#' The derivative of boundp
ndfboundp <- function(x, minb, maxb) (maxb-minb)*exp(-x)/(1+exp(-x))^2
#' The derivative of boundpin
#'
ndfboundpin <- function(y, minb, maxb) (maxb-minb)/( (y-minb)*(maxb-y))

boundp(-5, 0, 5)
boundpin(.033, 0,5)
ndfboundp(-5, 0, 5)

Npar <- 2
covar <- diag(Npar)

covar <- matrix(.5121, nrow=2, ncol=2)
diag(covar) <- c(1,1)
covar.inv <- solve(covar)

dyn.unload(dynlib('models/mvnb/mvnb_tmb'))
compile(file='models/mvnb/mvnb_tmb.cpp')
dyn.load(dynlib('models/mvnb/mvnb_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
inits <- list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2))
pars <- 'mu'
mvnb.obj <- MakeADFun(data=data, parameters=inits[[1]])


## Analytical MVN log-densities and gradients
a <- -2
b <- 2
fn <- function(x) as.vector(dmvnorm(as.vector(x), sigma=covar, log=TRUE))
gr <- function(x) -as.vector(covar.inv%*%x)
fnb <- function(x) fn(x=boundp(x, a,b))
grb <- function(x) gr(boundp(x, a,b))* ndfboundp(x,a,b)

x1.seq <- x2.seq <- seq(-5,5, len=10)
res <- ldply(x1.seq, function(x1)
  ldply(x2.seq, function(x2){
    print(c(x1,x2))
    data.frame(x1=x1, x2=x2, NLL1=mvnb.obj$fn(c(x1,x2)), NLL2=fnb(c(x1,x2)),
          grad1=mvnb.obj$gr(c(x1,x2)), grad2=t(grb(c(x1,x2))))
    }))
res
plot(0,0, xlim=range(x1.seq), ylim=range(x2.seq), type='n')
with(res, arrows(x0=x1, y0=x2, x1=x1+grad2.1, y1=x2+grad2.2, add=F))

library(shinystan)
source('mcmc.R')

set.seed(1)
out1 <- TMB:::run_mcmc(nsim=200, obj=mvnb.obj, params.init=c(0,0), L=50, eps=.5, alg='HMC')
set.seed(1)
out2 <- TMB:::run_mcmc.hmc(nsim=200, fn=fnb, gr=grb, params.init=c(0,0), L=50, eps=.5)
par(mfrow=c(2,2))
with(out1, plot(samples[,1,1], samples[,1,2]))
with(out2, plot(par[,1],par[,2]))
with(out1, plot(boundp(samples[,1,1], a,b), boundp(samples[,1,2], a,b)))
with(out2, plot(boundp(par[,1], a,b), boundp(par[,2], a,b)))
