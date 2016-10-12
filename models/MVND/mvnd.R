library(R2admb)
library(shinystan)
library(rstan)
## ------------------------------------------------------------
## The bounded transformation(s)
boundp <- function(x, minb, maxb, hbf){
    ## The internal transformations used in ADMB depending on the value of the
    ## Hybrid_bounded_flag (hbf) value.
    if(hbf==1)
        result <- minb+(maxb-minb)/(1+exp(-x))
    else if(hbf==0)
        result <- minb+(maxb-minb)*(.5*sin(x*pi/2)+.5)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}
boundpin <- function(x, minb, maxb, hbf) {
    ## The inverse of the transformation
    if(hbf==1)
        result <- -log( (maxb-x)/(x-minb) )
    else if(hbf==0)
        result <- asin(2*(x-minb)/(maxb-minb)-1)/(pi/2)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}
ndfboundp <- function(x, minb, maxb, hbf) {
    ## The derivative used to find the "scales"
    if(hbf==1)
        result <- (maxb-minb)*exp(-x)/(1+exp(-x))^2
    else if(hbf==0)
        result <- (maxb-minb)*.5*pi/2*cos(x*pi/2)
    else stop("Invalid hbf value, should be 0 or 1")
    return(result)
}
## ------------------------------------------------------------


N <- 2
covar <- diag(N)
covar <- matrix(c(1,-.95,-.95,1), nrow=2)
chd <- t(chol(covar))               # lower triangular Cholesky decomp.
A <- solve(chd)               # inverse
## covar <- matrix(rWishart(n=1, df=N, Sigma=diag(N)), nrow=N)
write.table(x=c(N, covar), file='mvnd.dat', , row.names=FALSE,
            col.names=FALSE)

## Run model in Stan
data <- list(covar=covar, Npar=N, x=rep(0, len=N))
inits <- list(list(mu=as.vector(mvtnorm::rmvnorm(1,sigma=covar))))
params.jags <- 'mu'
fit <- stan(file='mvnd.stan', data=data, iter=100000, init=inits, chains=1)
sims <- extract(fit, permuted=TRUE)$mu
eps <- tail(get_sampler_params(fit)[[1]][,'stepsize__'],1)
pairs(sims)

system('mvnd -mcmc 10000000 -mcsave 100')
rwm <- read_psv('mvnd')
system('mvnd -mcmc 100000 -hmc -hynstep 100 -hyeps .5')
hmc <- read_psv('mvnd')

## Thin them
rwm <- rwm[seq(1,nrow(rwm), len=10000),]
hmc <- hmc[seq(1,nrow(hmc), len=10000),]

pairs(rwm, pch='.')
pairs(hmc, pch='.')
chains <- array(data=NA, dim=c(nrow(rwm),2, ncol(rwm)))
chains[,1,] <- as.matrix(rwm)
chains[,2,] <- as.matrix(hmc)


## launch_shinystan(as.shinystan(chains))

zz <- hmc
par(mfrow=c(1,3))
## model space
plot(zz[,1], zz[,2], pch='.')
## unbounded model space
xx <- apply(zz, 2, function(x) boundpin(x, -2,2,1))
plot(xx[,1], xx[,2], pch='.')
## unbounded, transformed
yy <- t(A %*% t(xx))
plot(yy[,1], yy[,2], pch='.')
