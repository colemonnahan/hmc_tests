### ------------------------------------------------------------
## Code for testing performance of NUTS between Stan, TMB and ADMB, across
## real and simulated models. This script was adapted from some earlier,
## more intensive testing, so it's unfortuantely a bit awkward to use.

## The run_model function runs combinations of verification (to make
## qqplots aftering thinning), simulation (if appropriate, using the
## Nvar.vec vector) and empirical (the real model). Its kind of a mess but
## see function for more infor. It should just run without any
## modifications as long as you have the right packages installed.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'C:/Users/Cole/hmc_tests/'
setwd(main.dir)
source("startup.R")
Nreps <- 30                 # number of replicates
Nout.ind <- 1000            # number of independent samples if verify==TRUE
set.seed(241)
seeds <- sample(1:1e5, size=Nreps)         #
metric <- 'diag'

packageVersion('rstan') # 2.17.3
packageVersion('adnuts') # 1.0.0
packageVersion('tmbstan') # 1.0.0
## ADMB version 12.0. Models compiled with optimized mode.

### End of Step 1.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 2: Run the models.
## Run independent normal with variable SDs
data <- list(n=500, x=rep(0, 500), sds=1:500)
inits <- function() list(mu=rnorm(500))
obj.stan <- stan_model(file= 'models/zdiag/zdiag.stan')
Npar.vec <- 2^(3:10)
run_model(m='zdiag', obj.stan=obj.stan, data=data, inits=inits,
          verify=FALSE, simulation=TRUE, empirical=FALSE)

## VB growth, simulated
m <- 'growth'
temp <- growth_setup(N=20, seed=2345)
data <- temp$data; inits <- temp$inits
Npar.vec <- 2^(2:9) # will double again inside model
obj.stan <- stan_model(file= 'models/growth/growth.stan')
run_model(m='growth', obj.stan=obj.stan, data=data, inits=inits, delta=0.9,
          verify=FALSE, simulation=TRUE, empirical=FALSE, Nthin.ind=3,
          exp.columns=c(1,2,5,6))

## Wildflower
m <- 'wildf'
temp <- wildf_setup()
data <- temp$data
inits <- temp$inits
obj.stan <- stan_model(file= 'models/wildf/wildf.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.8,
          verify=FALSE, empirical=TRUE, Nthin.ind=5,
          Nout.ind=500, exp.columns=c(1,2,3))

## Swallows bird count data. See files for more info.
m <- 'swallows'
temp <- swallows_setup()
data <- temp$data
inits <- temp$inits
obj.stan <- stan_model(file= 'models/swallows/swallows.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.9,
          verify=FALSE, empirical=TRUE, Nthin.ind=10,
          Nout.ind=400, exp.columns=c(1,2,3))

## State space logistic fisheries assessment model. The adapt_delta needs
## to be high since it's not very well defined (classic banana
## shape). Plus, this model was converted from BUGS with really specific
## priors and so it's funky in other ways.
m <- 'sslog'
setwd(main.dir)
temp <- sslog_setup()
data <- temp$data
inits <- temp$inits
lower <- -Inf*abs(.1+unlist(inits()))
upper <- Inf*abs(.1+unlist(inits()))
lower['iq'] <- 1; upper['iq'] <- 10
obj.stan <- stan_model(file= 'models/sslog/sslog.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.98,
          verify=FALSE, empirical=TRUE, Nthin.ind=2,
          Nout.ind=500, lower=lower, upper=upper)

### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables

### old and not used.. moved to adnuts_tests for now

## setwd(main.dir)
## source('load_data.R')
## source('make_plots.R')
### End of Step 3.
### ------------------------------------------------------------
