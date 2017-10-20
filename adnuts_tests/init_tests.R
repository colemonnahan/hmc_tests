## Make sure that the init values are processed the same way for ADMB and
## TMB. This seems to be broken when using bounds.
devtools::install('C:/Users/Cole/adnuts')
library(adnuts)
library(shinystan)
library(TMB)


## Setup mvnb
m <- 'mvnb'
setwd(file.path('models',m))
Npar <- 3
covar <- matrix(.5, 3,3)
diag(covar) <- 1
data <- list(covar=covar, x=rep(0, len=Npar))
inits <- list(list(mu1=1, mu2=1, mu3=1))
pars <- c('mu1', 'mu2', 'mu3')
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj <- MakeADFun(data=data, parameters=inits[[1]])
obj$env$beSilent()
lower <- c(-2,0,-Inf)
upper <- c(2,5, Inf)


fit <- sample_tmb(obj, iter=2000, init=inits, lower=lower, upper=upper, control=list(metric=diag(3)))
pairs(extract_samples(fit))


## Setup the TMB parts manually. First the transformations.
cases <- .transform.cases(lower, upper)
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
## Rotation done using choleski decomposition of mass matrix M
M <- diag(3)
system('mvnb')
temp <- get.admb.cov()
M <- temp$cov.unbounded
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
## Redefine these functions
fn2 <- function(theta) fn(chd %*% theta)
gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )


## Need to adjust the current parameters so the chain is continuous
## initial bounded parameters =  model space

mle <- r4ss::read.admbFit(m)
z <- mle$est
z <- c(1.12,.91,10)
## unbounded
y <- .transform.inv(z, lower, upper, cases)
## unbounded + rotated = algorithm space
x <- as.numeric(chd.inv %*% y)

## TMB values
z;y;x;gr(y);gr2(x)


setwd('admb')
write.table(x=z, file='init.pin', row.names=FALSE, col.names=FALSE)
system('admb mvnb')
system('mvnb -mcmc 10 -nuts -hyeps .0001 -nohess')
adnuts:::write.admb.cov(M)
system('mvnb -mcmc 10 -nuts -hyeps .0001 -mcpin init.pin -noest')
setwd('..')
