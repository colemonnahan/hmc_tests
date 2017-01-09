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
se <- c(1, 25)
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
setwd('..')
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

## Quick tests
z0 <- c(.5,2.1)
lower <- c(0, -3.5*se[2]); upper <- c(Inf, 3.5*se[2])
lower <- c(-Inf, -Inf); upper <- c(Inf, Inf)
lower <- -2.5*se; upper <- 2.5*se
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

