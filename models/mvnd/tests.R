### This file quickly demonstrates how to use the static HMC and NUTS
### algorithms provided in the mcmc.R file in this folder. Here we use
### simple multivariate normal model with analytical gradients, avoiding
### the need for automatic differentiation. See function delcarations for
### argument definitions.

library(plyr)
library(mvtnorm)
source('../../mcmc.R')
col.label <- gray(.3)
col.border <- gray(.5)
col.tick <- gray(.4)
col.contour <- gray(.5)
lty.contour <- 2
lwd.contour <- .75

a <- -2
b <- 2

## log density in parameter space
fn <- function(y, covar) -t(y)%*%covar.inv%*%y
gr <- function(y, covar) as.vector(-covar.inv%*%y)
plot.gr <- function(y, covar){
  g <- gr(y,covar)
  points(t(y), pch=16)
  arrows(x0=y[1], y0=y[2], x1=g[1], y1=g[2])
}
gr2 <- function(y, A) -solve(t(A))%*% gr(y)
plot.gr2 <- function(y, A){
  g <- gr2(y, A)
  points(t(y), pch=16)
  arrows(x0=y[1], y0=y[2], x1=g[1], y1=g[2])
}
gr3 <- function(x, A){
  xbounded <- boundp(x=x, a,b,1)
  g <-  -t(solve(A)) %*% solve(A) %*% x
  x2 <- solve(A) %*% x
  ## apply chain run from bounding function
  scale1 <- ndfboundp(x=x2[1], minb=a, maxb=b, 1)
  scale2 <- ndfboundp(x=x2[2], minb=a, maxb=b, 1)
  return(g*c(scale1,scale2))
}
plot.gr3 <- function(x, A){
  xbounded <- boundp(x=x, a,b,1)
  g <- gr3(x, A)
  points(t(x), pch=16)
  arrows(x0=x[1], y0=x[2], x1=g[1], y1=g[2])
}


sims <- read.csv('stan.sims.csv')
par(mfrow=c(1,3))
plot(sims, pch='.')
sims2 <- cbind(boundpin(sims[,1], a,b,1), boundpin(sims[,2], a,b,1))
plot(sims2, pch='.')
sims3 <- t(A %*% t(sims2))
plot(sims3, pch='.')

## 2d parabola: y=x1^2+x2^2 with static HMC leapfrog trajectories, with and
## without mass matrix
covar <- matrix(c(1,.99,.99,1), 2)
##covar <- diag(2)
covar.inv <- solve(covar)
chd <- t(chol(covar))               # lower triangular Cholesky decomp.
A <- solve(chd)               # inverse
xx <- rmvnorm(1000, sigma=covar)
## make grid for contours
ngrid <- 50
x1.seq <- x2.seq <- seq(-3, 3, len=ngrid)
mvn.grid <- ldply(x1.seq, function(x1)
    ldply(x2.seq, function(x2){
        cbind(x1=x1, x2=x2, NLL=-dmvnorm(x=c(x1,x2), sigma=covar, log=TRUE))}))
mvn.grid$NLL[is.nan(mvn.grid$NLL)] <- NA
z.mvn <- matrix(exp(-mvn.grid$NLL), nrow=ngrid, ncol=ngrid, byrow=TRUE)
y1.seq <- temp[,1]; y2.seq <- temp[,2]
mvn.grid2 <- ldply(x1.seq, function(y1)
    ldply(x2.seq, function(y2){
      cbind(y1=y1, y2=y2, NLL2= -dmvnorm(x=c(y1,y2), sigma=diag(2), log=TRUE))}))
mvn.grid2$NLL2[is.nan(mvn.grid2$NLL2)] <- NA
z2.mvn <- matrix(exp(-mvn.grid2$NLL), nrow=ngrid, ncol=ngrid, byrow=TRUE)
mvn.grid3 <- ldply(x1.seq, function(y1)
    ldply(x2.seq, function(y2){
      cbind(y1=boundp(y1,a,b,1), y2=boundp(y2,a,b,1), NLL3= -dmvnorm(x=boundp(c(y1,y2), a,b,1), sigma=diag(2), log=TRUE))}))
mvn.grid3$NLL3[is.nan(mvn.grid3$NLL3)] <- NA
z3.mvn <- matrix(exp(-mvn.grid3$NLL3), nrow=ngrid, ncol=ngrid, byrow=TRUE)
nlevels <- 10


pts2 <- matrix(rnorm(4, sd=1), ncol=2)
pts <- t(solve(A) %*% t(pts2))
par(mfrow=c(1,3))
## Original space, unbounded
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE)
trash <- apply(pts, 1, function(x) plot.gr(x, covar))
## Rotated space, calculated via linear transformation
contour(x=x1.seq, y=x2.seq, z=z2.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE)
trash <- apply(pts2, 1, function(x) plot.gr2(x, A))
## Original space, bounded
contour(x=boundp(x1.seq,a,b,1), y=boundp(x2.seq,a,b,1), z=z3.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE)
trash <- apply(pts, 1, function(x) plot.gr3(x, A))







## log density in transformed space
2 <- function(x, A) A %*% x
fn2.inv <- function(x2, A) solve(A) %*% x2

fn <- function(x) -t(x)%*%covar.inv%*%x
gr <- function(x) as.vector(-covar.inv%*%x)
fn2 <- function(x) -t(x)%*%covar.inv%*%x
gr2 <- function(x) as.vector(-covar.inv%*%x)
(xx <- theta2.fn(c(.12345,-.54321), A))
(yy <- theta2.inv.fn(xx, A))


## Rotate the likelihood surface
par(mfrow=c(1,2))
contour(x=y1.seq, y=y2.seq, z=z2.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE)
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE)

## Do trajectory in transformed space and back calculate
init <- c(0,1)
r.cur <- c(.1,.1)
zz1 <- leap(L=5, fn=fn, gr=gr, eps=.1, init=init, r.cur=r.cur)

L.big <- 50
eps.big <- .1
## without mass matrix
mvn.leap1 <- leapfrog(L=L.big, fn=fn, gr=gr, eps=eps.big,
                      init=c(2.3, -2), r.cur=c(-1,.1))
## with it
mvn.leap2 <- leapfrog(L=L.big, fn=fn, gr=gr, eps=eps.big, covar=covar,
                      init=c(2.3, -2), r.cur=c(-1,.1))
mvn.leap.list <- list(mvn.leap1, mvn.leap2)
nlevels <- 8
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
a <- mvn.leap1$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05)
points(a[1,1], a[1,2], col='red', pch=16)
a <- mvn.leap2$par; N <- nrow(a)
arrows(a[-N,1], a[-N,2], a[-1,1], a[-1,2], length=.05, col='blue')
points(a[1,1], a[1,2], col='red', pch=16)
mtext(expression(x[1]), 1, line=2.5)
mtext(expression(x[2]), 2, line=2.5)


## Run NUTS with and without mass matrix
nuts.nocovar <- nuts(nsim=500, fn=fn, gr=gr, params.init=c(2.3, -2), delta=.8,
           diagnostic=TRUE)
nuts.covar <- nuts(nsim=500, fn=fn, gr=gr, params.init=c(2.3, -2), delta=.8,
           diagnostic=TRUE, covar=covar)
## The adaptation of the step size
plot(nuts.nocovar$epsbar, ylab='Step size', type='l')
lines(nuts.covar$epsbar, col='blue')
par(mfrow=c(1,2))
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
points(nuts.nocovar$par)
contour(x=x1.seq, y=x2.seq, z=z.mvn, nlevels=nlevels, col=col.contour, lty=lty.contour,
        drawlabels=FALSE, xlim=c(-2.5,2.5), ylim=c(-2.5,2.5))
points(nuts.covar$par)
