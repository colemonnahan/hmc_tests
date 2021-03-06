## Quick code to test that TMB and ADMB build tree functions do the same
## thing and work with bounds and arbitrary mass matrix
rm(list=ls())
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
devtools::load_all("c:/Users/Cole/adnuts")
source("buildtree2.R")

## !!!! THIS REQUIRES COMPILING A SPECIAL TESTING BRANCH OF ADMB !!!!
## Use: buildtree_tests. OTherwise this wont work at all.

## Generate single trajectory and make sure identical path and acceptance
## probabilities. Start with basic parabola as easiest case.
set.seed(6)
x0 <- c(1,.5)
r0 <- c(.27, -.63) #rnorm(2)
eps <- .5
setwd('parabola/')
compile(file='parabola_tmb.cpp')
dyn.load(dynlib('parabola_tmb'))
data <- list(Npar=2)
obj <- MakeADFun(data=data, parameters=list(x1=0, x2=0), DLL='parabola_tmb')
obj$env$beSilent()
fn2 <- function(x) -obj$fn(x)
gr2 <- function(x) -as.vector(obj$gr(x))
rm(theta.trajectory, pos='.GlobalEnv')
temp <- buildtree2(x0, r=r0, u=1e-15, v=1, j=10, eps=eps,
                   theta0=x0, r0=r0, fn=fn2, gr=gr2)
theta.tmb <- theta.trajectory[,1:2]
p.tmb <- theta.trajectory[, 5:6]
alpha.tmb <- theta.trajectory[,3:4]
dimnames(p.tmb) <- dimnames(theta.tmb) <- dimnames(alpha.tmb) <- NULL
## Repeat with admb
write.table(c(x0, r0, 10), file='input.txt', sep=' ', row.names=F,
            col.names=F)
write.table(x=2, file='parabola.dat', row.names=FALSE, col.names=FALSE)
system('admb parabola'); system('parabola')
file <- 'trajectory.txt'
pp <- if(file.exists(file)) file.remove(file)
system(paste0('parabola -noest -mcmc 1 -nuts -hyeps ', eps), ignore.stdout=F)
yy <- as.matrix(read.table(file, sep=' ')[,-1])
setwd('..')
dimnames(yy) <- NULL
alpha.admb <- yy[,7:8]
p.admb <- yy[,9:10]
theta.admb <- yy[,1:2]
alpha.tmb <- exp(alpha.tmb[,1]-alpha.tmb[,2])
alpha.admb <- exp(alpha.admb[,1]-alpha.admb[,2])
png('plots/single_trajectory.png', width=2.5, height=5, units='in', res=500)
par(mfrow=c(3,1), mar=c(3,3,.5,.5), mgp=c(1.5,.5,0), oma=c(0,0,1.5,0))
plot(alpha.tmb, (alpha.tmb-alpha.admb)/alpha.tmb,
     ylim=1e-6*c(-1,1), ylab=NA)
mtext('RE in alpha (TMB vs ADMB)', line=-1.5, cex=.8)
col1 <- 1; col2 <- 2
plot(theta.admb[,1], theta.admb[,2], col=col1, type='l', lwd=3,
      xlab='x[1]', ylab='x[2]')
lines(theta.tmb[,1], theta.tmb[,2], col=col2, type='l', lwd=1)
legend('top', legend=c('ADMB', 'TMB'), lty=1, col=c(col1, col2), ncol=2)
plot(p.admb[,1], p.admb[,2], col=col1, type='l', lwd=3,
      xlab='p[1]', ylab='p[2]')
lines(p.tmb[,1], p.tmb[,2], col=col2, type='l', lwd=1)
dev.off()


## Setup more complicated posterior surface
set.seed(235)
corr <- matrix(-.954, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, 25)
covar <- corr * (se %o% se)
covar.inv <- solve(covar)
zz <- rmvnorm(n=1e4, sigma=covar)
setwd('mvnb2/')
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
f <- function(x0, r0, u, v=1, j, eps, chd){
  rm(theta.trajectory, pos='.GlobalEnv')
  temp <- buildtree2(x0, r=r0, u=u, v=v, j=j, eps=eps,
                     theta0=x0, r0=r0, fn=fn2, gr=gr2)
  traj.x <- theta.trajectory[,1:2]
  alpha.tmb <- exp(theta.trajectory[,3]-theta.trajectory[,4])
  dimnames(traj.x) <- NULL
  traj.y <- t(apply(traj.x, 1, function(i) chd %*% i))
  traj.z <- t(sapply(1:nrow(traj.y), function(i)
    .transform(traj.y[i,], lower,upper, cases)))
  return(list(s=temp$s, x=traj.x, y=traj.y, z=traj.z,
              theta.prime=temp$theta.prime, alpha.tmb=alpha.tmb))
}
plot.tree <- function(a,c,space, main=NA, add=TRUE){
  ## Function to add trajectory to surfaces.
  if(add) {
    xlim <- quantile(x=a[,1], probs=c(2/10000, .9998))
    ylim <- quantile(x=a[,2], probs=c(2/10000, .9998))
    ## xlim <- ylim <- NULL
    plot(a, pch='.', main=NA, xlab=NA, ylab=NA,
       #  xlab=paste0(space, '[1]'), ylab=paste0(space, '[2]'),
         col=rgb(0,0,0,.5),
         xlim=xlim, ylim=ylim)
    mtext(main, line=.5)
    }
  points(c[1,1], c[1,2],col='red', pch=16, cex=.5)
  lines(c, type='l', col='red', pch=16, cex=.5, lwd=1)
}## Quick tests
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

## Make trajectories with and without a mass matrix for TMB
## Add multiple random trajectories
set.seed(3512)
nsim <- 15
eps <- .01
lower <- c(-1, -4)*se; upper <- c(1,4)*se
cases <- .transform.cases(lower,upper)
## Filter out points out of bounds
samples.z <-
  zz[zz[,1] > lower[1] & zz[,1] < upper[1] &
     zz[,2] > lower[2] & zz[,2] < upper[2],]
samples.y <- t(apply(samples.z, 1, function(i)
  .transform.inv(i, lower, upper, cases)))
r0 <- matrix(rnorm(nsim*2), ncol=2)
v <- rep(1, nsim)                       # all 1 so matches admb
png('plots/tree_trajectories.png', width=7, height=5, units='in', res=500)
par(mfrow=c(3,3), mar=c(2,2,.5,.1), mgp=c(1.5,.5,0), oma=c(0,1,1.5,.5))
## First with unit diag M
M <- diag(2)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
x.ind <- sample(1:nrow(samples.x), size=nsim)
tmp <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  return(f(x0, r0[i,], u=1e-15, v=v[i], j=15, eps=eps, chd=chd))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp[[i]]$x, space='x', main='Unbounded', add=(i==1)))
mtext("Identity M", side=2, line=1.5)
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
  return(f(x0, r0[i,], u=1e-15, v=v[i], j=15, eps=eps, chd=chd))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp2[[i]]$x, space='x', add=(i==1)))
mtext("Diagonal M", side=2, line=1.5)
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
  return(f(x0, r0[i,], u=1e-15, v=v[i], j=15, eps=eps, chd=chd))
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp3[[i]]$x, space='x', add=(i==1)))
mtext("Dense M", side=2, line=1.5)
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.y, tmp3[[i]]$y, space='y', add=(i==1)))
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.z, tmp3[[i]]$z, space='z', add=(i==1)))
dev.off()
## Check some properties of the trajectories. Why so many duplicated number
## of steps for tmp3?? Shouldn't they be random??
any(unlist(lapply(tmp, function(x) x$s)) )
any(unlist(lapply(tmp2, function(x) x$s)))
any(unlist(lapply(tmp3, function(x) x$s)))
png('plots/trajectory_lengths.png', width=7, height=3, units='in', res=500)
par(mfrow=c(1,3), mar=c(3,3,1,.1), mgp=c(1.5,.5,.05), oma=c(0,0,1,.5))
barplot(table(unlist(lapply(tmp, function(x) nrow(x$x)))))
mtext("Identity M");box()
barplot(table(unlist(lapply(tmp2, function(x) nrow(x$x)))))
mtext("Diagonal M");box()
barplot(table(unlist(lapply(tmp3, function(x) nrow(x$x)))))
mtext("Dense M");box()
mtext("Frequency", side=2, outer=TRUE, line=-1.5)
mtext("Trajectory length (steps)", side=1, outer=TRUE, line=-1.4)
dev.off()
## End of trajectories


## Make sure TMB and ADMB produce identical buildtree trajectories given
## the same model and initial values.
compare.traj <- function(x0,r0, j=15, eps=.05, seed=1, lower,upper, covar){
  setwd('mvnb2')
  on.exit(setwd('..'))
  write.table(c(x0, r0, j), file='input.txt', sep=' ', row.names=F, col.names=F)
  write.table(x=c(2, covar, lower, upper), file='mvnb2.dat', row.names=FALSE, col.names=FALSE)
  ## system('admb mvnb2')
  file <- 'trajectory.txt'
  pp <- if(file.exists(file)) file.remove(file)
  file2 <- 'theta_prime.txt'
  pp <- if(file.exists(file2)) file.remove(file2)
  system(paste0('mvnb2 -noest -mcmc 1 -nuts -hyeps ', eps, ' -mcseed ',seed), ignore.stdout=F)
  traj.admb <- as.matrix(read.table(file, sep=' ')[-1,-1])
  theta.admb <- as.matrix(read.table('theta_prime.txt'))
  alpha.admb <- exp(traj.admb[,7]-traj.admb[,8])
  set.seed(seed)
  res <- f(x0, r0, exp(-15), 1, j=j, eps=eps, chd=chd)
  traj.tmb <- cbind(res$x, res$y, res$z)
  theta.tmb <- matrix(res$theta.prime, ncol=2)
  return(list(traj.tmb=traj.tmb, traj.admb=traj.admb,
              theta.admb=theta.admb, theta.tmb=theta.tmb,
              alpha.tmb=res$alpha.tmb[-1], alpha.admb=alpha.admb))
}
## Make a single trajectory in three spaces for both platforms
set.seed(352)
M <- cov(samples.y)
chd <- t(chol(M))
write.admb.cov(cov.unbounded=M, model.path='mvnb2')
yy <- compare.traj(rnorm(2), rnorm(2), j=15, eps=.01, lower=lower,
                   upper=upper, covar=covar)
range(yy$alpha.tmb-yy$alpha.admb)
png('plots/single_trajectory_bounded.png', width=7, height=3, units='in', res=500)
par(mfrow=c(1,3), mar=c(3,3,1.5,.1), mgp=c(1.5,.5,.05), oma=c(0,0,1.5,0))
plot(yy$traj.admb[,1], yy$traj.admb[,2], col=col1, type='l', lwd=3,
     main="Unbounded, Rotated", xlab='x[1]', ylab='x[2]')
lines(yy$traj.tmb[,1], yy$traj.tmb[,2], col=col2, type='l', lwd=1)
plot(yy$traj.admb[,3], yy$traj.admb[,4], col=col1, type='l', lwd=3,
          main="Unbounded", xlab='y[1]', ylab='y[2]')
lines(yy$traj.tmb[,3], yy$traj.tmb[,4], col=col2, type='l', lwd=1)
plot(yy$traj.admb[,5], yy$traj.admb[,6], col=col1, type='l', lwd=3,
          main="Model Space", xlab='z[1]', ylab='z[2]')
lines(yy$traj.tmb[,5], yy$traj.tmb[,6], col=col2, type='l', lwd=1)
legend('bottomleft', legend=c('ADMB', 'TMB'), lwd=c(3,1), col=c(col1,col2))
dev.off()


## Check that it is drawing uniformly from the trajectory for both ADMB and
## TMB. Run the same trajectory but get different theta_primes by using a
## different seed.
## !!!! This takes a long time to run !!!!!
x0 <- c(1,2)
r0 <- c(.112, -1.147)
test <- ldply(1:5000, function(i){
  yy <- compare.traj(x0, r0, j=4, eps=.1, seed=i, lower=lower, upper=upper,
                    covar=covar)
  data.frame(tmb=yy$theta.tmb[1], admb=yy$theta.admb[1])})
x1 <- sort(as.numeric(as.factor(test[,1])))
x2 <- sort(as.numeric(as.factor(test[,2])))
## length(unique(x1))
## length(unique(x2))
## sort(unique(test[,1]))
## sort(unique(test[,2]))

png('plots/trajectory_sampling.png', width=6, height=3, units='in', res=500)
par(mfrow=c(1,2), mar=c(3,2,1.5,.1), oma=c(0,0,2,0), mgp=c(1.5,.5,0))
barplot(table(x1), main="TMB", xlab="Step # selected");box()
barplot(table(x2), main="ADMB",xlab="Step # selected");box()
mtext('5000 sampled steps from same trajectory', line=0, outer=TRUE, cex=1.5)
dev.off()

