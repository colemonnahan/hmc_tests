## Quick code to test TMB bounding
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
logit.inv <- function(y) 1/(1+exp(-y))
#' The logistic transformation function for bounding parameters in MCMC
boundp <- function(x, minb, maxb) minb+(maxb-minb)/(1+exp(-x))
#' The inverse of the transformation
boundpin <- function(y, minb, maxb) -log( (maxb-y)/(y-minb) )
#' The derivative of boundp
ndfboundp <- function(x, minb, maxb) (maxb-minb)*exp(-x)/(1+exp(-x))^2

## First test with a 1d normal. Three cases, mvnd1 is normal with bounds
## created in R, mvnb1 is bounded in the template, and analytically in R.
a <- -2
b <- 2
covar <- covar.inv <- 1
Npar <- 1

dyn.unload(dynlib('models/mvnb1/mvnb_tmb'))
compile(file='models/mvnb1/mvnb_tmb.cpp')
dyn.load(dynlib('models/mvnb1/mvnb_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar), bounds=c(a,b))
mvnb.obj <- MakeADFun(data=data, parameters=list(mu=rep(0, NROW(covar))), DLL='mvnb_tmb')
nlminb(start=0.01*c(1), objective=mvnb.obj$fn, gradient=mvnb.obj$gr)
mvnb.obj$env$beSilent()
mvnb.obj$report()

##
dyn.unload(dynlib('models/mvnd1/mvnd_tmb'))
compile(file='models/mvnd1/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd1/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=0), DLL='mvnd_tmb')
nlminb(start=0.01, objective=mvnd.obj$fn, gradient=mvnd.obj$gr)
mvnd.obj$env$beSilent()
mvnd.obj$report()

## Analytical MVN log-densities and gradients
## fn <- function(x) as.vector(dmvnorm(as.vector(x), sigma=covar, log=TRUE))
## gr <- function(x) -as.vector(covar.inv%*%x)
fn <- function(x) -x^2
gr <- function(x) -2*x
plot(x <- seq(a,b, len=1000), y=fn(x))
plot(x <- seq(a,b, len=1000), y=gr(x))
fnb <- function(y){
  scales <- ndfboundp(y,a,b)
  fn(boundp(y, a,b)) + log(scales)
  }
grb <- function(y) {
  gr(boundp(y, a,b))* ndfboundp(y,a,b)-1+2/(1+exp(y))
}
fn2 <- function(y){
 -mvnd.obj$fn(boundp(y, a,b)) + sum(log(abs(ndfboundp(y,a,b))))
}
gr2 <- function(y) -mvnd.obj$gr(boundp(y, a,b))* ndfboundp(y,a,b) -1 +2/(1+exp(y))
x.seq <- seq(-10,10, len=1000)
par(mfrow=c(2,3))
plot(x.seq, y=fnb(x.seq))
plot(x.seq, y=sapply(x.seq, function(i) fn2(i)), col='red')
plot(x.seq, y=-sapply(x.seq, function(i) mvnb.obj$fn(i)), col='red')
plot(x.seq, y=grb(x.seq)); abline(h=0)
plot(x.seq, y=sapply(x.seq, function(i) gr2(i)), col='red')
plot(x.seq, y=-sapply(x.seq, function(i) mvnb.obj$gr(i)), col='red')
x1 <- rnorm(1)
mvnb.obj$gr(x1)
gr2(x1)
grb(x1)
-mvnb.obj$fn(x1)+mvnb.obj$fn(x2)
fnb(x1)-fnb(x2)
fn2(x1)-fn2(x2)

## Quick MCMC tests
out.rwm <- TMB:::run_mcmc.rwm(nsim=20000, fn=fnb, params.init=c(0))
plot(boundp(out.rwm[,1], a,b))
## RWM works fine, so fnb must be right
out.nuts <- TMB:::run_mcmc.nuts(nsim=2000, fn=fnb, gr=grb, params.init=c(0), eps=NULL, max_doublings=6)
##with(out.nuts, plot(boundp(par[,1], a,b), boundp(par[,2], a,b)))
with(out.nuts, plot(boundp(par[,1], a,b)))
## OK it works for 1d example.

## 2d example
covar <- matrix(.354, nrow=2, ncol=2)
diag(covar) <- 1
covar.inv <- solve(covar)
Npar <- 2

## dyn.unload(dynlib('models/mvnb/mvnb_tmb'))
## compile(file='models/mvnb/mvnb_tmb.cpp')
## dyn.load(dynlib('models/mvnb/mvnb_tmb'))
## data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar), bounds=c(a,b))
## mvnb.obj <- MakeADFun(data=data, parameters=list(mu=rep(0, 2)), DLL='mvnb_tmb')
## mvnb.obj$env$beSilent()
## nlminb(start=0.01*c(1,1), objective=mvnb.obj$fn, gradient=mvnb.obj$gr)
## mvnb.obj$report(par=c(1,1))

##
dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')
nlminb(start=c(1,1), objective=mvnd.obj$fn, gradient=mvnd.obj$gr)
mvnd.obj$env$beSilent()

## Analytical MVN log-densities and gradients
fn <- function(x) as.vector(dmvnorm(as.vector(x), sigma=covar, log=TRUE))
gr <- function(x) -as.vector(covar.inv%*%x)
fnb <- function(y){
  x <- sapply(1:Npar, function(i) .transform(y[i], a[i], b[i], cases[i]))
  scales <- sapply(1:Npar, function(i) .transform.grad(y[i],a[i],b[i], cases[i]))
  fn(x) + sum(log(scales))
  }
grb <- function(y) {
  x <- sapply(1:Npar, function(i) .transform(y[i], a[i], b[i], cases[i]))
  scales <- sapply(1:Npar, function(i) .transform.grad(y[i],a[i],b[i], cases[i]))
  scales2 <- sapply(1:Npar, function(i) .transform.grad2(y[i],a[i],b[i], cases[i]))
  gr(x)*scales+scales2
}
fn2 <- function(y){
 -mvnd.obj$fn(boundp(y, a,b)) + sum(log(abs(ndfboundp(y,a,b))))
}
gr2 <- function(y) -mvnd.obj$gr(boundp(y, a,b))* ndfboundp(y,a,b) -1 +2/(1+exp(y))
x.seq <- seq(-10,10, len=1000)
par(mfrow=c(2,3))
x1 <- rnorm(2)
x2 <- rnorm(2)
mvnb.obj$gr(x1)
gr2(x1)
grb(x1)
-mvnb.obj$fn(x1)+mvnb.obj$fn(x2)
fnb(x1)-fnb(x2)
fn2(x1)-fn2(x2)

## Quick MCMC tests
a <- c(-5,-5); b <- c(Inf, Inf)
a <- b <- NULL
a <- c(-Inf, -Inf); b <- c(Inf, Inf)
.transform.cases(a,b)
out.rwm <- run_mcmc(obj=mvnd.obj, nsim=20000,params.init=c(0,0),
                    lower=a, upper=b, alg='RWM', chain=2, alpha=.1)
plot(out.rwm$samples[,1,1], out.rwm$samples[,1,2])
out.nuts <- run_mcmc(obj=mvnd.obj, nsim=2000, alg='NUTS', max_doublings=7,
                     chains=2, lower=a, upper=b)
plot(out.nuts$samples[,1,1], out.nuts$samples[,1,2])
sso.tmb <- with(out.nuts, as.shinystan(samples, burnin=warmup, max_treedepth=8,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso.tmb)

x1.seq <- x2.seq <- seq(-5,5, len=10)
res <- ldply(x1.seq, function(x1)
  ldply(x2.seq, function(x2){
    data.frame(x1=x1, x2=x2, NLL1=mvnb.obj$fn(c(x1,x2)), NLL2=fnb(c(x1,x2)),
          grad1=mvnb.obj$gr(c(x1,x2)), grad2=t(grb(c(x1,x2))))
    }))
res
plot(0,0, xlim=range(x1.seq), ylim=range(x2.seq), type='n')
eps <- .5
with(res, arrows(x0=x1, y0=x2, x1=x1+eps*grad2.1, y1=x2+eps*grad2.2))

library(shinystan)
source('mcmc.R')

set.seed(1)
out.rwm <- run_mcmc.rwm(nsim=2000, fn=fnb, params.init=c(0,0))
out1 <- run_mcmc(nsim=2000, obj=mvnb.obj, params.init=c(0,0), L=50, eps=.5, alg='HMC')
set.seed(1)
out2 <- run_mcmc.hmc(nsim=2000, fn=fnb, gr=grb, params.init=c(0,0),
                           L=50, eps=.1)
out2 <- run_mcmc.nuts(nsim=2000, fn=fnb, gr=grb, params.init=c(0,0), eps=.1)
par(mfrow=c(2,2))
with(out1, plot(samples[,1,1], samples[,1,2]))
with(out2, plot(par[,1],par[,2]))
with(out1, plot(boundp(samples[,1,1], a,b), boundp(samples[,1,2], a,b)))
with(out2, plot(boundp(par[,1], a,b), boundp(par[,2], a,b)))
plot(boundp(out.rwm[,1], a,b))
plot(boundp(out.rwm[,2], a,b))
plot(sims.stan[,-3])



## boundp(-5, 0, 5)
## boundpin(.033, 0,5)
## ndfboundp(-5, 0, 5)


## Compare TMB vs Stan
chains <- 3
iter <- 500
td <- 9
inits <- rep(list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2)),chains)
stan.fit <- stan(file= 'models/mvnb/mvnb.stan', data=data, iter=iter, par='mu',
                  chains=chains, thin=1, algorithm='NUTS',
                   init=inits, seed=1, verbose=FALSE,
                   control=list(adapt_engaged=FALSE, metric='unit_e', stepsize=.3))
sims.stan <- data.frame(extract(stan.fit, permuted=FALSE)[,1,])
sso.stan <- as.shinystan(stan.fit)
launch_shinystan(sso.stan)

tmb.fit <- run_mcmc(nsim=iter, obj=mvnb.obj, alg='NUTS', lower=c(-100,-100),
                    upper=c(100,100), eps=NULL, chains=1, max_doubling=td,
                    covar=NULL)
sims.tmb <- data.frame(tmb.fit$samples[-(1:iter/2),1,])
sso.tmb <- with(tmb.fit, as.shinystan(samples, burnin=warmup, max_treedepth=td,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso.tmb)

tmb.fit2 <- run_mcmc(nsim=iter, obj=mvnb.obj, alg='NUTS', lower=c(-2,-2),
                    upper=c(2,2), eps=NULL, chains=chains, max_doubling=td,
                    covar=covar)
sims.tmb2 <- data.frame(tmb.fit$samples[-(1:iter/2),1,])
sso.tmb2 <- with(tmb.fit2, as.shinystan(samples, burnin=warmup, max_treedepth=td,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso.tmb2)

thin <- 50
qqplot(sims.stan[seq(1,nrow(sims.stan), by=thin),1],
       sims.tmb[seq(1,nrow(sims.tmb), by=thin),1])
