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
Niter <- 1000

dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')

lower <- c(-1.5,-2); upper <- c(0,1.13)
init <- lapply(1:3, function(x) rnorm(2, sd=10))
## No covar
out1 <- run_mcmc(mvnd.obj, nsim=Niter, init=init, chains=3)
out2 <- run_mcmc(mvnd.obj, nsim=Niter, init=init, lower=lower, upper=upper,
                 chains=3)
## With covar
out3 <- run_mcmc(mvnd.obj, nsim=Niter, init=init, covar=covar, chains=3)
out4 <- run_mcmc(mvnd.obj, nsim=Niter, init=init, lower=lower, upper=upper,
                 covar=covar, chains=3)

get_efficiency_label<- function(x){
  samples <- x$samples[-(1:x$warmup),1,, drop=FALSE]
  minESS <- min(monitor(samples, probs=.5, print=FALSE)[,5])
  paste0("minESS/time= ", round(minESS/x$time.total[1],2))
}
par(mfrow=c(2,2))
ylim <- xlim <- c(-4,4)
plot(out1$samples[-(1:out1$warmup),1,], xlim=xlim, ylim=ylim)
text(3,-3, labels=get_efficiency_label(out1))
plot(out2$samples[-(1:out2$warmup),1,], xlim=xlim, ylim=ylim)
abline(v=c(lower[1], upper[1]), h=c(lower[2], upper[2]))
text(3,-3, labels=get_efficiency_label(out2))
plot(out3$samples[-(1:out1$warmup),1,], xlim=xlim, ylim=ylim)
text(3,-3, labels=get_efficiency_label(out3))
plot(out4$samples[-(1:out2$warmup),1,], xlim=xlim, ylim=ylim)
text(3,-3, labels=get_efficiency_label(out4))
abline(v=c(lower[1], upper[1]), h=c(lower[2], upper[2]))



## First test with a 1d normal. Two cases: mvnb1 is bounded in the
## template, and analytically in R.
a <- -2.124
b <- 3.215
cases <- .transform.cases(a,b)
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

## Analytical MVN log-densities and gradients
fn <- function(x) -x^2
gr <- function(x) -2*x
fnb <- function(y){
  x <- .transform(y, a, b, cases)
  scales <- .transform.grad(y, a, b, cases)
  fn(x) + sum(log(abs(scales)))
}
grb <- function(y){
  x <- .transform(y, a, b, cases)
  scales <- .transform.grad(y, a, b, cases)
  scales2 <- .transform.grad2(y, a, b, cases)
  gr(x)*scales + scales2
}
x1 <- rnorm(1)
x2 <- rnorm(1)
-mvnb.obj$gr(x1)
grb(x1)
-mvnb.obj$fn(x1)+mvnb.obj$fn(x2)
fnb(x1)-fnb(x2)
x.seq <- seq(-10,10, len=1000)
par(mfrow=c(2,3))
plot(x.seq, y=sapply(x.seq, function(i) fnb(i))); abline(v=0)
plot(x.seq, y=-sapply(x.seq, function(i) mvnb.obj$fn(i)), col='red')
plot(x.seq, y=sapply(x.seq, function(i) grb(i))); abline(h=0)
plot(x.seq, y=-sapply(x.seq, function(i) mvnb.obj$gr(i)), col='red'); abline(h=0)

## Quick MCMC tests
out.rwm <- run_mcmc.rwm(nsim=20000, fn=fnb, init=c(0))
plot(sapply(1:nrow(out.rwm$par), function(i) .transform(out.rwm$par[i,1], a,b, cases)))
## RWM works fine, so fnb must be right
out.nuts <- run_mcmc.nuts(nsim=2000, fn=fnb, gr=grb, init=c(0), eps=NULL, max_doublings=6)
##with(out.nuts, plot(boundp(par[,1], a,b), boundp(par[,2], a,b)))
plot(sapply(1:nrow(out.nuts$par), function(i)
  .transform(out.nuts$par[i,1], a, b, cases)))
## OK it works for 1d example.

## 2d example
a <- c(-2.42, -2.42); b <- c(3.5, 3.5)
covar <- matrix(.954, nrow=2, ncol=2)
cases <- .transform.cases(a,b)
diag(covar) <- 1
covar.inv <- solve(covar)
Npar <- 2

dyn.unload(dynlib('models/mvnb/mvnb_tmb'))
compile(file='models/mvnb/mvnb_tmb.cpp')
dyn.load(dynlib('models/mvnb/mvnb_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar), bounds=c(a[1],b[1]))
mvnb.obj <- MakeADFun(data=data, parameters=list(mu=rep(0, 2)), DLL='mvnb_tmb')
mvnb.obj$env$beSilent()
nlminb(start=0.01*c(1,1), objective=mvnb.obj$fn, gradient=mvnb.obj$gr)
mvnb.obj$report(c(1.2,4.15))
log(.transform.grad(c(1.2, 4.15), a,b, cases))

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
x.seq <- seq(-10,10, len=1000)
par(mfrow=c(2,3))
x1 <- rnorm(2)
x2 <- rnorm(2)
-mvnb.obj$gr(x1)
grb(x1)
-mvnb.obj$fn(x1)+mvnb.obj$fn(x2)
fnb(x1)-fnb(x2)

## Quick MCMC tests
a <- c(-5,-5); b <- c(Inf, Inf)
a <- b <- NULL
a <- c(-Inf, -Inf); b <- c(Inf, Inf)

out.rwm <- run_mcmc(obj=mvnd.obj, nsim=20000,init=c(0,0),
                    lower=a, upper=b, alg='RWM', chain=1, alpha=.1,
                    thin=10)
plot(out.rwm$samples[,1,1], out.rwm$samples[,1,2])
out.nuts <- run_mcmc(obj=mvnd.obj, nsim=2000, alg='NUTS', max_doublings=10,
                     chains=3, lower=a, upper=b)
plot(out.nuts$samples[,1,1], out.nuts$samples[,1,2])
sso.tmb <- with(out.nuts, as.shinystan(samples, burnin=warmup, max_treedepth=8,
             sampler_params=sampler_params, algorithm='NUTS', model_name='test'))
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
out.rwm <- run_mcmc.rwm(nsim=2000, fn=fnb, init=c(0,0))
out1 <- run_mcmc(nsim=2000, obj=mvnb.obj, init=c(0,0), L=50, eps=.5, alg='HMC')
set.seed(1)
out2 <- run_mcmc.hmc(nsim=2000, fn=fnb, gr=grb, init=c(0,0),
                           L=50, eps=.1)
out2 <- run_mcmc.nuts(nsim=2000, fn=fnb, gr=grb, init=c(0,0), eps=.1)
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
