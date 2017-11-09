### ------------------------------------------------------------
## Code for testing performance of NUTS between Stan, TMB and ADMB, across
## real and simulated models.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'C:/Users/Cole/hmc_tests/'
setwd(main.dir)
source("startup.R")
Nreps <- 4                 # number of replicates
Nout.ind <- 1000            # number of independent samples if verify==TRUE
set.seed(241)
seeds <- 1:Nreps#sample(1:1e5, size=Nreps)         #
metric <- 'diag'#c('unit_e', 'diag_e', 'dense_e')[1]
## Suppress output to file via sink? Useful after debugging to clean up
## console and judge progress.
sink <- FALSE
### End of Step 1.
### ------------------------------------------------------------

### ------------------------------------------------------------


### Step 2: Run the models.
## Run multivariate normal, empirical and simulated
##Npar.vec <- c(2,4,8,16,32,64, 128)
Npar <- 16
covar <- diag(Npar)
data <- list(Npar=Npar, covar=covar, x=rep(0, len=Npar))
inits <- function() list(mu=rnorm(n=Npar, mean=0, sd=sqrt(diag(covar))))
obj.stan <- stan_model(file= 'models/mvnd/mvnd.stan')
run_model(m='mvnd', data=data, inits=inits, pars=pars, verify=FALSE)

## Run iid normal increasing in size
## Setup data, inits and pars
data <- list(n=50, x=rep(0, 50))
inits <- function() list(mu=rnorm(50))
obj.stan <- stan_model(file= 'models/iidz/iidz.stan')
Npar.vec <- 2^(4+1:4)
run_model(m='iidz', obj.stan=obj.stan, data=data, inits=inits,
          simulation=TRUE, empirical=TRUE, verify=FALSE)

## Run independent normal with variable SDs
data <- list(n=50, x=rep(0, 50), sds=1:50)
inits <- function() list(mu=rnorm(50))
obj.stan <- stan_model(file= 'models/zdiag/zdiag.stan')
Npar.vec <- 2^(4+1:4)
run_model(m='zdiag', obj.stan=obj.stan, data=data, inits=inits,
          verify=FALSE, simulation=TRUE, empirical=TRUE)

## VB growth, simulated
m <- 'growth'
temp <- growth_setup(N=30, seed=2345)
data <- temp$data; inits <- temp$inits
Npar.vec <- 2^(3+1:4)
obj.stan <- stan_model(file= 'models/growth/growth.stan')
run_model(m='growth', obj.stan=obj.stan, data=data, inits=inits, delta=0.9,
          verify=TRUE, simulation=TRUE, empirical=TRUE, Nthin.ind=3,
          exp.columns=c(1,2,5,6))

## Wildflower
m <- 'wildf'
temp <- wildf_setup()
data <- temp$data
inits <- temp$inits
obj.stan <- stan_model(file= 'models/wildf/wildf.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.95,
          verify=TRUE, simulation=FALSE, empirical=FALSE, Nthin.ind=5,
          Nout.ind=500, exp.columns=c(1,2,3))

## Swallows bird count data
m <- 'swallows'
temp <- swallows_setup()
data <- temp$data
inits <- temp$inits
lower <- abs(unlist(inits()))*-Inf
upper <- abs(unlist(inits()))*Inf
lower[c('sigmayearphi', 'sigmaphi', 'sigmap')] <- 0
obj.stan <- stan_model(file= 'models/swallows/swallows.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.9,
          verify=FALSE, simulation=FALSE, empirical=TRUE, Nthin.ind=1,
          Nout.ind=500, lower=lower, upper=upper, exp.columns=c(1,2,3))

## Simulated spatial model, TMB example
m <- 'spatial'
temp <- spatial_setup()
data <- temp$data; inits <- temp$inits
lower <- abs(unlist(inits()))*-Inf
upper <- abs(unlist(inits()))*Inf
lower['sigma'] <- 0; lower['a'] <- 0
upper['a'] <- 2
obj <- MakeADFun(data=data, parameters=inits(), DLL = "spatial", random = "u")
opt <- nlminb(obj$par, obj$fn, obj$gr,
              lower=c(-100.0, -100.0, 0.01, -3.0),
              upper=c( 100.0,  100.0, 3.00,  3.0) )
test <- tmbstan(obj, lower=lower, upper=upper, chains=3, iter=700,
                init=lapply(1:3, function(i) inits()))

### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables
setwd(main.dir)
source('load_data.R')
source('make_plots.R')
### End of Step 3.
### ------------------------------------------------------------

## options(mc.cores=parallel::detectCores()-1)
## test <- stan('models/swallows/swallows.stan', data=data, init=inits, iter=700, chains=3,
##                  pars=names(inits()))

