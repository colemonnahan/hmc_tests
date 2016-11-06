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
Nreps <- 3                 # number of replicates
Nout.ind <- 1000            # number of independent samples if verify==TRUE
set.seed(241)
seeds <- 1:Nreps#sample(1:1e5, size=Nreps)         #
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
Nout <- 5000                            # no. of samples out
Nthin <- 1                              # thin rate for emp/sim modes (leave at 1!)
Nthin.ind <- 1                          # thin rate for verify mode
## Settings for simulation mode. cor is a factor for independent (0) or
## from wishart (1) (see paper). Npar is how many parameters.
cor.vec <- c(0,1)
Npar.vec <- c(5, 15, 25, 50, 100, 200, 300, 400)
source(paste0('models/',m,'/run_model.R'))

## Quick exploration of mvnd plots
setwd(main.dir)
adapt <- read.csv('results/mvnd_adapt_simulated.csv')[,-1]
ggplot(adapt, aes(Npar, eps.final, group=platform, color=platform)) +
  geom_point() + facet_wrap('cor')
adapt.long <- melt(adapt, id.vars=c('platform', 'seed', 'Npar', 'cor', 'Nsims'),
                   measure.vars=c('delta.mean', 'eps.final',
                                  'max_treedepths', 'nsteps.mean'))
ggplot(adapt.long, aes(Npar, value, group=platform, color=platform)) + geom_point() +
  facet_grid(variable~cor, scales='free')

perf <- read.csv('results/mvnd_perf_simulated.csv')[,-1]
ggplot(perf, aes(Npar, samples.per.time, group=platform, color=platform)) +
  geom_point() + facet_wrap('cor') + scale_y_log10()
ggplot(perf, aes(Npar, minESS/Nsims, group=platform, color=platform)) +
  geom_point() + facet_wrap('cor') + ylim(0,1)
ggplot(perf, aes(Npar, time.total, group=platform, color=platform)) +
  geom_point() + facet_wrap('cor') + scale_y_log10()

## Run multivariate normal, empirical and simulated
m <- 'growth'
verify <- FALSE
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

