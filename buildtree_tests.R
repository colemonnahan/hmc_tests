## Quick code to test that TMB and ADMB build tree functions do the same
## thing and work with bounds and arbitrary mass matrix
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
devtools::load_all("c:/Users/Cole/rnuts")


covar <- matrix(-.54, nrow=2, ncol=2)
diag(covar) <- 1
covar.inv <- solve(covar)

## Make single buildtree trajectory using TMB
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
if(TRUE){
  fn2 <- function(theta) fn(chd %*% theta)
  gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
  chd <- t(chol(covar))               # lower triangular Cholesky decomp.
  chd.inv <- solve(chd)               # inverse
} else {
  fn2 <- fn; gr2 <- gr
}


lower <- c(-3.1,-3.8); upper <- c(3.4,2.9)
lower <- 10*c(-3.1,-3.8); upper <- 10*c(3.4,2.9)
## lower <- c(-Inf, -Inf); upper <- c(Inf, Inf)
cases <- .transform.cases(lower, upper)
z0 <- c(.5,2.1)
y0 <- .transform.inv(z0, a=lower, b=upper, cases=cases)
x0 <- as.vector(chd.inv %*% y0)
## ## Recover z0
## y02 <- as.vector(chd %*% x0)
## z0- as.vector(.transform(y02, lower, upper, cases))
## ## Recover gradients. Why doesnt this work??
## .transform(as.vector(chd %*% gr2(x0)), upper,lower,cases)
## -as.vector(covar.inv%*%z0)
zz <- rmvnorm(n=1e4, sigma=covar)
## filter out those outside bounds
samples.z <- zz[zz[,1] > lower[1] & zz[,1] < upper[1] &
              zz[,2] > lower[2] & zz[,2] < upper[2],]
samples.y <- t(sapply(1:nrow(samples.z), function(i)
  .transform.inv(samples.z[i,], lower,upper, cases)))
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
u <- 1e-5
r0 <- -c(1.5,-1.1)
f <- function(x0, r0, u, v=1, j, eps){
  rm(theta.trajectory, pos='.GlobalEnv')
  temp <- buildtree2(x0, r=r0, u=u, v=v, j=j, eps=eps,
                  theta0=x0, r0=r0, fn=fn2, gr=gr2)
  traj.x <- rbind(x0,theta.trajectory)
  dimnames(traj.x) <- NULL
  traj.y <- t(apply(traj.x, 1, function(i) chd %*% i))
  traj.z <- t(sapply(1:nrow(traj.y), function(i)
    .transform(traj.y[i,], lower,upper, cases)))
  return(list(x=traj.x, y=traj.y, z=traj.z, theta.prime=temp$theta.prime))
}
res <- f(x0, r0, u=1e-5, v=1, j=10, eps=.05)


## Add single trajectory to surfaces
x0 <- c(-.1,.2)
y0 <- as.vector(chd %*% x0)
z0 <- .transform(y0, a=lower, b=upper, cases=cases)
res <- f(x0, rnorm(2), u=1e-5, v=1, j=15, eps=.005)
par(mfrow=c(1,3))
hh <- function(a,c, main, add=TRUE){
  if(add) plot(a, pch='.', main=main)
  points(c[1,1], c[1,2],col='red', pch=16)
  lines(c, type='l', col='red', pch=16, cex=.5, lwd=1)
}
hh(samples.x, res$x, main='Unbounded + Rotated')
hh(samples.y, res$y, main='Unbounded + Rotated')
hh(samples.z, res$z, main='Unbounded + Rotated')

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


