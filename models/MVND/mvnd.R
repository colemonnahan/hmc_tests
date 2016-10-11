library(R2admb)
library(shinystan)
library(rstan)

N <- 2
covar <- diag(N)
covar <- matrix(c(1,-.95,-.95,1), nrow=2)
## covar <- matrix(rWishart(n=1, df=N, Sigma=diag(N)), nrow=N)
write.table(x=c(N, covar), file='mvnd.dat', , row.names=FALSE,
            col.names=FALSE)

## Setup for Stan
data <- list(covar=covar, Npar=N, x=rep(0, len=N))
inits <- list(list(mu=as.vector(mvtnorm::rmvnorm(1,sigma=covar))))
params.jags <- 'mu'
fit <- stan(file='mvnd.stan', data=data, iter=100000, init=inits, chains=1)
sims <- extract(fit, permuted=TRUE)$mu
eps <- tail(get_sampler_params(fit)[[1]][,'stepsize__'],1)
pairs(sims)

system('mvnd -mcmc 10000000 -mcsave 100')
rwm <- read_psv('mvnd')
system('mvnd -mcmc 100000 -hybrid -hynstep 100 -hyeps .5')
hmc <- read_psv('mvnd')

## Thin them
rwm <- rwm[seq(1,nrow(rwm), len=10000),]
hmc <- hmc[seq(1,nrow(hmc), len=10000),]

pairs(rwm, pch='.')
pairs(hmc, pch='.')
chains <- array(data=NA, dim=c(nrow(rwm),2, ncol(rwm)))
chains[,1,] <- as.matrix(rwm)
chains[,2,] <- as.matrix(hmc)


launch_shinystan(as.shinystan(chains))
launch_shinystan(as.shinystan(hmc2))

ans <- R2admb:::read_admbbin(model)

read.psv <- function(model){
psv <- file(paste0(model,".psv"), "rb")
nparams <- readBin(psv, "integer", n = 1)
mcmc <- matrix(readBin(psv, "numeric", n = nparams * Nout), ncol = nparams,
               byrow = TRUE)
close(psv)
