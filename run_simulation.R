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
data <- list(n=5, x=rep(0, 5))
inits <- function() list(mu=rnorm(5))
pars <- 'mu'
obj.stan <- stan_model(file= 'models/iidz/iidz.stan')
Npar.vec <- 2^(4+1:4)
run_model(m='iidz', obj.stan=obj.stan, data=data, inits=inits, pars=pars, verify=FALSE,
          simulation=TRUE, empirical=FALSE)

## Run independent normal with variable SDs
data <- list(n=5, x=rep(0, 5), sds=1:5)
inits <- function() list(mu=rnorm(5))
pars <- 'mu'
obj.stan <- stan_model(file= 'models/zdiag/zdiag.stan')
Npar.vec <- 2^(4+1:4)
run_model(m='zdiag', obj.stan=obj.stan, data=data, inits=inits, pars=pars,
          verify=FALSE, simulation=TRUE, empirical=FALSE)


## VB growth, simulated
## Run independent normal with variable SDs
m <- 'growth_nc'
dat <-
  sample.lengths(Nfish=30, n.ages=5)
data <- list(Nfish=30, Nobs=nrow(dat), loglengths=dat$loglengths,
                  fish=dat$fish, ages=dat$ages)
inits <- function()
  list(sigma_obs=runif(1, .01, .2), logLinf_mean=runif(1, 2, 5),
       logk_mean=runif(1,-4,0), logLinf_sigma=runif(1, .01, .4),
       logk_sigma=runif(1, .01, .4), logLinf=rnorm(30, 4, 1),
       logk=rnorm(30, -2, 1))
pars <- NULL
obj.stan <- stan_model(file= 'models/growth_nc/growth_nc.stan')
Npar.vec <- 2^(4+1:4)
run_model(m='growth_nc', obj.stan=obj.stan, data=data, inits=inits, pars=pars,
          verify=FALSE, simulation=TRUE, empirical=FALSE)

verify <- TRUE
delta <- 0.8
Nout <- 500
Nthin <- 1
Nthin.ind <- 1
Npar.vec <- c(5, 15, 25, 50, 100, 200)
source(paste0('models/',m,'/run_model.R'))


### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables
setwd(main.dir)
source('load_data.R')
source('make_plots.R')
### End of Step 3.
### ------------------------------------------------------------

