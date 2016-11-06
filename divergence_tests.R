## Quick code to test TMB divergences from NUTS
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)

compile(file='rosenbrock.cpp')
dyn.load(dynlib('rosenbrock'))
data <- list()
pars <- c('x1', 'x2')
obj <- MakeADFun(data=data, parameters=list(x1=0, x2=0))
nlminb(start=c(2,-1), objective=obj$fn, gradient=obj$gr)
out <- run_mcmc(obj, nsim=500, alg='NUTS', eps=.1)
sum(out$sampler_params[[1]][,5])
plot(out$samples[,1,])
sso.tmb <- with(out, as.shinystan(samples, burnin=warmup, max_treedepth=12,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso.tmb)

data <- list(Y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
## dyn.unload(dynlib('rosenbrock'))
compile(file='8schools.cpp')
dyn.load(dynlib('8schools'))
obj <- MakeADFun(data=data, parameters=list(mu=1, logtau=-5, theta=rep(0,8)))
out <- run_mcmc(obj, nsim=2000, alg='NUTS', eps=NULL, chains=3)
sso.tmb <- with(out, as.shinystan(samples, burnin=warmup, max_treedepth=12,
             sampler_params=sampler_params, algorithm='NUTS'))
launch_shinystan(sso.tmb)



plot(0,0, type='n', xlim=c(-1.5,1.5), ylim=c(-.5,3))

rm(theta.trajectory)
theta0 <- theta <- c(0,0)
r0 <- r <- rnorm(2)
v <- 1
u <- TMB:::.sample.u(theta=theta, r=r, fn=obj$fn)
j <- 12
divergent <<- 0
info <- as.environment(list(n.calls=0))
xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j, eps=.1,
                 theta0=theta0, r0=r0, fn=fn2, gr=gr2, info=info)
c(info$n.calls, divergent)
points(theta[1], theta[2], pch=16)
points(theta.trajectory[,2], theta.trajectory[,3], type='b', xlim=c(-3,3), ylim=c(-3,3))
with(xx, points(theta.prime[1], theta.prime[2], pch=16, col='red', cex=.5))


with(xx, (theta.plus-theta.minus) %*% r.minus)
drop(with(xx, crossprod(theta.plus-theta.minus, r.minus)))
cbind(step=0, t(c(1,1)))

xx <- t(sapply(1:5000, function(i) .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=.3,
                 theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)$theta.prime))
barplot(table(xx[,1]))
