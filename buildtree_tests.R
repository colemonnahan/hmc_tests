## Quick code to test that TMB and ADMB build tree functions do the same
## thing and work with bounds and arbitrary mass matrix
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
devtools::load_all("c:/Users/Cole/rnuts")

## Setup posterior surface
corr <- matrix(-.954, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, 5)
covar <- corr * (se %o% se)
covar.inv <- solve(covar)
zz <- rmvnorm(n=1e4, sigma=covar)
## Setup temporary functions
source("buildtree2.R")
setwd('mvnb2/')
dyn.unload(dynlib('mvnb2_tmb'))
compile(file='mvnb2_tmb.cpp')
dyn.load(dynlib('mvnb2_tmb'))
data <- list(covar=covar, Npar=2, x=rep(0, len=2))
obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnb2_tmb')
obj$env$beSilent()
## Manual bounding functions
fn <- function(y){
  x <- .transform(y, lower, upper, cases)
  scales <- .transform.grad(y, lower, upper, cases)
  -obj$fn(x) + sum(log(scales))
}
gr <- function(y){
  x <- .transform(y, lower, upper, cases)
  scales <- .transform.grad(y, lower, upper, cases)
  scales2 <- .transform.grad2(y, lower, upper, cases)
  -as.vector(obj$gr(x))*scales + scales2
}
## Rotated, bounded functions
fn2 <- function(theta) fn(chd %*% theta)
gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
## Quick fun to run a single buildtree trajectory and return key things,
## used to plot below
f <- function(x0, r0, u, v=1, j, eps){
  rm(theta.trajectory, pos='.GlobalEnv')
  temp <- buildtree2(x0, r=r0, u=u, v=v, j=j, eps=eps,
                  theta0=x0, r0=r0, fn=fn2, gr=gr2)
  traj.x <- rbind(x0,theta.trajectory)
  dimnames(traj.x) <- NULL
  traj.y <- t(apply(traj.x, 1, function(i) chd %*% i))
  traj.z <- t(sapply(1:nrow(traj.y), function(i)
    .transform(traj.y[i,], lower,upper, cases)))
  return(list(s=temp$s, x=traj.x, y=traj.y, z=traj.z, theta.prime=temp$theta.prime))
}
## Quick tests
z0 <- c(.5,2.1)
lower <- -3*se; upper <- 3*se
lower <- c(-Inf, -3*se[2]); upper <- c(Inf, 3*se[2])
cases <- .transform.cases(lower, upper)
y0 <- .transform.inv(z0, a=lower, b=upper, cases=cases)
## Check gradients in y space
h <- 1e-6
c(fn(y0+c(h,0))-fn(y0), fn(y0+c(0,h))-fn(y0))/h
gr(y0)
M <- diag(2); M[1,1] <- .234
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
x0 <- as.vector(chd.inv %*% y0)
## Recover z0  by calculating backwards, to test that my functions are
## working properly.
y02 <- as.vector(chd %*% x0)
z0- as.vector(.transform(y02, lower, upper, cases))
## Test gradients of x0 using finite differences and compare to gr2
c(fn2(x0+c(h,0))-fn2(x0), fn2(x0+c(0,h))-fn2(x0))/h
gr2(x0)

## Add single trajectory to surfaces, with and without a mass matrix for
## both TMB and ADMB. THIS NEEDS TO WORK!!
hh <- function(a,c, main, add=TRUE){
  if(add) plot(a, pch='.', main=main)
  points(c[1,1], c[1,2],col='red', pch=16)
  lines(c, type='l', col='red', pch=16, cex=.5, lwd=1)
}
z0 <- -c(2,-12.1)
lower <- c(-Inf, -Inf); upper <- c(Inf, Inf)
lower <- -3*se; upper <- 3*se
cases <- .transform.cases(lower, upper)
y0 <- .transform.inv(z0, lower,upper,cases)
samples.z <- zz[zz[,1] > lower[1] & zz[,1] < upper[1] &
              zz[,2] > lower[2] & zz[,2] < upper[2],]
samples.y <- t(sapply(1:nrow(samples.z), function(i)
  .transform.inv(samples.z[i,], lower,upper, cases)))
## First with unit diag M
r0 <- rnorm(2)#c(0,0)
M <- diag(2)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
x0 <- as.vector(chd.inv %*% y0)
gr2(x0)
res <- f(x0=x0, r0=r0, u=1e-5, v=1, j=15, eps=.005)
par(mfrow=c(2,3))
hh(samples.x, res$x, main='Unbounded')
hh(samples.y, res$y, main='Uounded + Rotated')
hh(samples.z, res$z, main='Model Space')
## Now with M=covar
M <- cov(samples.y)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
x0 <- as.vector(chd.inv %*% y0)
gr2(x0)
res <- f(x0=x0, r0=r0, u=1e-5, v=1, j=15, eps=.005)
hh(samples.x, res$x, main='Unbounded')
hh(samples.y, res$y, main='Unbounded + Rotated')
hh(samples.z, res$z, main='Model Space')



## Add multiple random trajectories
nsim <- 15
temp <- lapply(1:nsim, function(i){
  x0 <- rnorm(2, mean=c(-.1, .2), sd=c(.1, .1))#c(-.1,1.2)
  y0 <- as.vector(chd %*% x0)
  z0 <- .transform(y0, a=lower, b=upper, cases=cases)
  res <- f(x0, rnorm(2), u=1e-5, v=1, j=15, eps=.005)
  return(res)
})
par(mfrow=c(1,3))
tt <- sapply(1:nsim, function(i)
  hh(samples.x, temp[[i]]$x, main='Unbounded + Rotated', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  hh(samples.y, temp[[i]]$y, main='Unbounded + Rotated', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  hh(samples.z, temp[[i]]$z, main='Unbounded + Rotated', add=(i==1)))




## Check that it is drawing uniformly from the trajectory
test <- ldply(1:50, function(i){
  u <- .sample.u(theta=y0, r=r0, fn=fn)
 .buildtree(x0, r=r0, u=u, v=1, j=0, eps=5000,
            theta0=x0, r0=r0, fn=fn, gr=gr)$theta.prime})
test.bounded <- t(sapply(1:nrow(test), function(i)
  .transform(as.matrix(test)[i,], lower,upper, cases)))

unique(test[,1])
barplot(table(test[,1]))


