## Quick code to test that TMB and ADMB build tree functions do the same
## thing and work with bounds and arbitrary mass matrix
rm(list=ls())
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
devtools::load_all("c:/Users/Cole/admbtools")
devtools::load_all("c:/Users/Cole/rnuts")


set.seed(235)
## Setup posterior surface
corr <- matrix(-.954, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, 25)
covar <- corr * (se %o% se)
covar.inv <- solve(covar)
zz <- rmvnorm(n=1e4, sigma=covar)
admb.traj <- function(x0,r0, j=15, eps=.05, seed=1){
  setwd('mvnb2')
  on.exit(setwd('..'))
  write.table(c(x0, r0, j), file='input.txt', sep=' ', row.names=F, col.names=F)
  write.table(x=c(2, covar, lower, upper), file='mvnb2.dat', row.names=FALSE, col.names=FALSE)
  ## system('admb mvnb2')
  file <- 'trajectory.txt'
  pp <- if(file.exists(file)) file.remove(file)
  file2 <- 'theta_prime.txt'
  pp <- if(file.exists(file2)) file.remove(file2)
  system(paste0('mvnb2 -noest -mcmc 1 -nuts -hyeps ', eps, ' -mcseed ',seed), ignore.stdout=T)
  theta.trajectory <- as.matrix(read.table(file, sep=' ')[,-1])
  traj.x <- theta.trajectory[,1:2]
  traj.y <- theta.trajectory[,3:4]
  traj.z <- theta.trajectory[,5:6]
  dimnames(traj.x) <- dimnames(traj.y) <- dimnames(traj.z) <- NULL
  return(list(x=traj.x, y=traj.y, z=traj.z))
}
## Make trajectories with and without a mass matrix for ADMB
plot.tree <- function(a,c,space, main=NA, add=TRUE){
  ## Function to add trajectory to surfaces.
  if(add) {
    xlim <- quantile(x=a[,1], probs=c(2/10000, .9998))
    ylim <- quantile(x=a[,2], probs=c(2/10000, .9998))
    ## xlim <- ylim <- NULL
    plot(a, pch='.', main=NA, xlab=paste0(space, '[1]'),
         ylab=paste0(space, '[2]'), col=rgb(0,0,0,.5),
         xlim=xlim, ylim=ylim)
    mtext(main, line=.5)
    }
  points(c[1,1], c[1,2],col='red', pch=16, cex=.5)
  lines(c, type='l', col='red', pch=16, cex=.5, lwd=1)
}

## Add multiple random trajectories in all three spaces for three mass
## matrices
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
rr <- matrix(rnorm(nsim*2), ncol=2)
png('plots/tree_trajectories_admb.png', width=7, height=5, units='in', res=500)
par(mfrow=c(3,3), mar=c(3,3,.5,.1), mgp=c(1.5,.5,.05), oma=c(0,0,1.5,0))
## First with unit diag M
M <- diag(2)
chd <- t(chol(M))               # lower triangular Cholesky decomp.
chd.inv <- solve(chd)               # inverse
write.admb.cov.unbounded(cov.unbounded=M, 'mvnb2')
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
x.ind <- sample(1:nrow(samples.x), size=nsim)
tmp <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  admb.traj(x0, rr[i,], j=15, eps=eps)
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
write.admb.cov.unbounded(cov.unbounded=M, 'mvnb2')
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
tmp2 <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  admb.traj(x0, rr[i,], j=15, eps=eps)
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
write.admb.cov.unbounded(cov.unbounded=M, 'mvnb2')
samples.x <- t(apply(samples.y, 1, function(i) chd.inv %*% i))
tmp3 <- lapply(1:nsim, function(i){
  x0 <- samples.x[x.ind[i],]
  admb.traj(x0, rr[i,], j=15, eps=eps)
})
tt <- sapply(1:nsim, function(i)
  plot.tree(samples.x, tmp3[[i]]$x, space='x', add=(i==1)))
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
png('plots/trajectory_lengths_admb.png', width=7, height=3, units='in', res=500)
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


