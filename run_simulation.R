### ------------------------------------------------------------
## Code for testing performance of NUTS between Stan, TMB and ADMB, across
## real and simulated models.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'C:/Users/Cole/hmc_tests/'
setwd(main.dir)
source("startup.R")
Nreps <- 3                 # number of replicates
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
m <- 'iidz'                             # model name
verify <- FALSE                         # whether to verify
delta <- 0.8                            # adapt_delta for Stan
Nout <- 500                            # no. of samples out
Nthin <- 1                              # thin rate for emp/sim modes (leave at 1!)
Nthin.ind <- 1                          # thin rate for verify mode
Npar.vec <-2^(1:8)
source(paste0('models/',m,'/run_model.R'))
## Run multivariate normal, empirical and simulated



## Run multivariate normal, empirical and simulated
m <- 'mvnd'                             # model name
verify <- FALSE                         # whether to verify
delta <- 0.8                            # adapt_delta for Stan
Nout <- 1000                            # no. of samples out
Nthin <- 1                              # thin rate for emp/sim modes (leave at 1!)
Nthin.ind <- 1                          # thin rate for verify mode
## Settings for simulation mode. cor is a factor for independent (0) or
## from wishart (1) (see paper). Npar is how many parameters.
cor.vec <- c(0,1)
## Npar.vec <- c(5, 15, 25, 50, 100, 200, 300, 400)
Npar.vec <- c(2,4,8,16,32,64, 128)
source(paste0('models/',m,'/run_model.R'))
## Run multivariate normal, empirical and simulated

## Bounded MVN to test parameter transformations. Has 3 parameters, one
## bounded above and below, one bounded only below, and one unbounded. ADMB
## has to do the second case manually. Hence some funky code.
m <- 'mvnb'
verify <- TRUE
delta <- 0.8
Nout <- 1000
Nthin <- 1
Nthin.ind <- 10
source(paste0('models/',m,'/run_model.R'))

## Run multivariate normal, empirical and simulated
m <- 'growth'
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

