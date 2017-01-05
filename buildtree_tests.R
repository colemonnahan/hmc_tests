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
lower <- c(0, -3.5*se[2]); upper <- c(Inf, 3.5*se[2])
lower <- c(-Inf, -Inf); upper <- c(Inf, Inf)
lower <- -3.5*se; upper <- 3.5*se
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


## Make trajectories with and without a mass matrix for
## both TMB and ADMB. THIS NEEDS TO WORK!!
plot.tree <- function(a,c,space, main=NA, add=TRUE){
  ## Function to add trajectory to surfaces.
  if(add) {
    xlim <- quantile(x=a[,1], probs=c(.025, .975))
    ylim <- quantile(x=a[,2], probs=c(.025, .975))
    plot(a, pch='.', main=NA, xlab=paste0(space, '[1]'),
         ylab=paste0(space, '[2]'), col=rgb(0,0,0,.5),
         xlim=xlim, ylim=ylim)
    mtext(main, line=.5)
    }
  points(c[1,1], c[1,2],col='red', pch=16, cex=.5)
  lines(c, type='l', col='red', pch=16, cex=.5, lwd=1)
}
## Add multiple random trajectories
nsim <- 50
eps <- .01
x.ind <- sample(1:nrow(samples.x), size=nsim)
r0 <- matrix(rnorm(nsim*2), ncol=2)
v <- sample(c(-1,1), nsim, replace=TRUE)
png('plots/tree_trajectories.png', width=7, height=5, units='in', res=500)
par(mfrow=c(3,3), mar=c(3,3,.5,.1), mgp=c(1.5,.5,.05), oma=c(0,0,1.5,0))
## First with unit diag M
M <- diag(2)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
tmp <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  return(f(x0, r0[i,], u=1e-5, v=v[i], j=15, eps=eps))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp[[i]]$x, space='x', main='Unbounded', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.y, tmp[[i]]$y, space='y', main='Unbounded + Rotated', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.z, tmp[[i]]$z, space='z',  main='Model Space', add=(i==1)))
## Now with M= estimated diagonals
M <- matrix(0, nrow=2, ncol=2)
diag(M) <- apply(samples.y, 2, var)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
tmp2 <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  return(f(x0, r0[i,], u=1e-5, v=v[i], j=15, eps=eps))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp2[[i]]$x, space='x', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.y, tmp2[[i]]$y, space='y', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.z, tmp2[[i]]$z, space='z', add=(i==1)))
## Now with M=covar
M <- cov(samples.y)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
tmp3 <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  return(f(x0, r0[i,], u=1e-5, v=v[i], j=15, eps=eps))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp3[[i]]$x, space='x', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.y, tmp3[[i]]$y, space='y', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.z, tmp3[[i]]$z, space='z', add=(i==1)))
dev.off()

any(unlist(lapply(tmp, function(x) x$s)) )
any(unlist(lapply(tmp2, function(x) x$s)))
any(unlist(lapply(tmp3, function(x) x$s)))
par(mfrow=c(1,3))
barplot(table(unlist(lapply(tmp, function(x) nrow(x$x)))))
barplot(table(unlist(lapply(tmp2, function(x) nrow(x$x)))))
barplot(table(unlist(lapply(tmp3, function(x) nrow(x$x)))))



## Check that it is drawing uniformly from the trajectory
test <- ldply(1:50, function(i){
  u <- .sample.u(theta=y0, r=r0, fn=fn)
 .buildtree(x0, r=r0, u=u, v=1, j=0, eps=5000,
            theta0=x0, r0=r0, fn=fn, gr=gr)$theta.prime})
test.bounded <- t(sapply(1:nrow(test), function(i)
  .transform(as.matrix(test)[i,], lower,upper, cases)))

unique(test[,1])
barplot(table(test[,1]))


