## Quick code to test TMB bounding
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)

Npar <- 2
covar <- matrix(.5121, nrow=2, ncol=2)
diag(covar) <- c(1,1)
covar.inv <- solve(covar)

dyn.unload(dynlib('models/mvnb/mvnb_tmb'))
compile(file='models/mvnb/mvnb_tmb.cpp')
dyn.load(dynlib('models/mvnb/mvnb_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
pars <- 'mu'
mvnb.obj <- MakeADFun(data=data, parameters=inits[[1]])
nlminb(start=c(1,-1)*10, objective=mvnb.obj$fn, gradient=mvnb.obj$gr)

chains <- 3
iter <- 500
td <- 5
inits <- rep(list(list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar)))/2)),chains)
stan.fit <- stan(file= 'models/mvnb/mvnb.stan', data=data, iter=iter, par='mu',
                  chains=chains, thin=1, algorithm='NUTS',
                   init=inits, seed=1, verbose=FALSE,
                   control=list(adapt_engaged=FALSE, metric='unit_e', stepsize=.3))
sims.stan <- data.frame(extract(stan.fit, permuted=FALSE)[,1,])
sso.stan <- as.shinystan(stan.fit)
launch_shinystan(sso.stan)

tmb.fit <- run_mcmc(nsim=iter, obj=mvnb.obj, alg='NUTS', lower=c(-2,-2),
                    upper=c(2,2), eps=5, chains=1, max_doubling=td,
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

## Analytical MVN log-densities and gradients
a <- -2
b <- 2
fn <- function(x) as.vector(dmvnorm(as.vector(x), sigma=covar, log=TRUE))
gr <- function(x) -as.vector(covar.inv%*%x)
fnb <- function(y){
  yinv <- logit.inv(y)
  scales <- (b-a)*yinv*(1-yinv)
  fn(a+(b-a)*yinv) + sum(log(scales))
  }
grb <- function(x) gr(boundp(x, a,b))* ndfboundp(x,a,b)
plot(x <- seq(-5,5, len=1000), y=boundp(x, a,b))

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
plot(boundp(out.rwm[,1], a,b), boundp(out.rwm[,2], a,b))
plot(boundp(out.rwm[,1], a,b))
plot(boundp(out.rwm[,2], a,b))
plot(sims.stan[,-3])


## logit.inv <- function(y) 1/(1+exp(-y))
## #' The logistic transformation function for bounding parameters in MCMC
## boundp <- function(x, minb, maxb) minb+(maxb-minb)/(1+exp(-x))
## #' The inverse of the transformation
## boundpin <- function(y, minb, maxb) -log( (maxb-y)/(y-minb) )
## #' The derivative of boundp
## ndfboundp <- function(x, minb, maxb) (maxb-minb)*exp(-x)/(1+exp(-x))^2

## boundp(-5, 0, 5)
## boundpin(.033, 0,5)
## ndfboundp(-5, 0, 5)
