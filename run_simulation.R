### ------------------------------------------------------------
## Code for testing performance of NUTS between Stan, TMB and ADMB, across
## real and simulated models.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
## Set the reference directory from which everything else is relative.
main.dir <- 'C:/Users/Cole/hmc_tests/'
setwd(main.dir)
source("startup.R")
Nreps <- 5                 # number of replicates
Nout.ind <- 1000            # number of independent samples if verify==TRUE
set.seed(241)
seeds <- sample(1:1e5, size=Nreps)         #
## Which metric to use for NUTS. Paper only used estimated diagonal (Stan
## default) but others could be used.
metric <- c('unit_e', 'diag_e', 'dense_e')[1]
## Suppress output to file via sink? Useful after debugging to clean up
## console and judge progress.
sink <- TRUE
## The paper was run with these key software versions
version$version.string                  # R version 3.2.3
packageVersion('rstan')                 # 2.11.1
packageVersion('R2jags')                # 0.5.7
packageVersion('rjags')                 # 4.4
### End of Step 1.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 2: Run the models.
## Run multivariate normal, empirical and simulated
m <- 'mvnd'                             # model name
verify <- FALSE                         # whether to verify
delta <- 0.8                            # adapt_delta for Stan
Nout <- 20000                           # no. of samples out
Nthin <- 1                              # thin rate for emp/sim modes (leave at 1!)
Nthin.ind <- 100                        # thin rate for verify mode
## Settings for simulation mode. cor is a factor for independent (0) or
## from wishart (1) (see paper). Npar is how many parameters.
cor.vec <- c(0,1)
Npar.vec <- c(2, 5, 15, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run MVN with varying correlations and a fixed Npar
m <- 'mvnc'
delta <- 0.8
Npar <- 5
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
cor.vec <- c(0, .25, .5, .75, .85, .9, .95)
Npar.vec <- c(2, 25, 50)
source(paste0('models/',m,'/run_model.R'))

## Run simulated growth tests, cross between centered/noncentered
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
Npar.vec <- c(2:10, floor(10^(seq(1.25, 2.5, by=.25))))
delta <- 0.8
m <- 'growth'
source(paste0('models/',m,'/run_model.R'))
delta <- 0.8
m <- 'growth_nc'
source(paste0('models/',m,'/run_model.R'))

## State space logistic. Note had to raise delta here.
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
delta <- 0.95
m <- 'ss_logistic'
source(paste0('models/',m,'/run_model.R'))
m <- 'ss_logistic_nc'
delta <- 0.95
source(paste0('models/',m,'/run_model.R'))

## Red kite example from Kery and Schaub; 8.4 w/ informative prior
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
delta <- 0.8
m <- 'redkite'
source(paste0('models/',m,'/run_model.R'))

## swallows -- Modified from example 14.5 from Korner-Nievergelt et al
Nout <- 20000; Nthin <- 1; Nthin.ind <- 100
delta <- .8
m <- 'swallows_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'swallows'
source(paste0('models/',m,'/run_model.R'))

## Wildflower; from Bolker et al 2013
Nout <- 20000; Nthin <- 1; Nthin.ind <- 500
delta <- .8
m <- 'wildflower_nc'
source(paste0('models/',m,'/run_model.R'))
delta <- .8
m <- 'wildflower'
source(paste0('models/',m,'/run_model.R'))
### End of Step 2.
### ------------------------------------------------------------

### ------------------------------------------------------------
### Step 3: Load and prepare result data frames for plotting and tables
source('load_data.R')
source('make_plots.R')
### End of Step 3.
### ------------------------------------------------------------

