### ------------------------------------------------------------
## Code for testing performance of NUTS between Stan, TMB and ADMB, across
## real and simulated models.

### ------------------------------------------------------------
### Step 1: prepare working space; load libraries, functions, and global
### variables
main.dir <- 'C:/Users/Cole/hmc_tests/'
setwd(main.dir)
source("startup.R")
Nreps <- 10                 # number of replicates
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
Npar.vec <- 2^(4+1:6)
run_model(m='iidz', obj.stan=obj.stan, data=data, inits=inits,
          simulation=TRUE, empirical=TRUE, verify=FALSE)

## Run independent normal with variable SDs
data <- list(n=50, x=rep(0, 50), sds=1:50)
inits <- function() list(mu=rnorm(50))
obj.stan <- stan_model(file= 'models/zdiag/zdiag.stan')
Npar.vec <- 2^(4+1:6)
run_model(m='zdiag', obj.stan=obj.stan, data=data, inits=inits,
          verify=FALSE, simulation=TRUE, empirical=TRUE)

## VB growth, simulated
m <- 'growth'
temp <- growth_setup(N=128, seed=2345)
data <- temp$data; inits <- temp$inits
Npar.vec <- 2^(3+1:5)
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
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.95,
          verify=FALSE, empirical=TRUE, Nthin.ind=5,
          Nout.ind=500, exp.columns=c(1,2,3))

## Swallows bird count data
m <- 'swallows'
temp <- swallows_setup()
data <- temp$data
inits <- temp$inits
obj.stan <- stan_model(file= 'models/swallows/swallows.stan')
run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.9,
          verify=TRUE, simulation=FALSE, empirical=TRUE, Nthin.ind=10,
          Nout.ind=500, exp.columns=c(1,2,3))

## Still trying to get this model working
m <- 'sslog'
setwd(main.dir)
temp <- sslog_setup()
data <- temp$data
inits <- temp$inits
obj.stan <- stan_model(file= 'models/sslog/sslog.stan')
test <- sampling(obj.stan, iter=2000, pars=names(inits()[[1]]), data=data,
                 init=inits,
                 control=list(adapt_delta=.98))
setwd(paste0('models/',m))
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj.tmb <- MakeADFun(data=data, parameters=inits(), DLL=m)
test <- tmbstan(obj.tmb, iter=2000, init=list(inits()), chains=1 )
setwd('admb')
write.table(x=unlist(data), file=paste0(m,'.dat'), row.names=FALSE,
            col.names=FALSE )
## Don't need a good hessian since using adaptation now, so just run it a
## single iteration to get the files with names
system(paste('admb',m))
system(paste(m, ' -maxfn 1 -nox -nohess'))


run_model(m=m, obj.stan=obj.stan, data=data, inits=inits, delta=.9,
          verify=FALSE, simulation=FALSE, empirical=TRUE, Nthin.ind=1,
          Nout.ind=500, exp.columns=c(1,2,3))



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

## ## Old tests for a growth model which was slow to converge and thus
## adapating poorly
## devtools::install("C:/Users/Cole/adnuts")
## NN <- 15
## set.seed(1)
## seeds0 <- sample(1:1e5, size=NN)
## ii <- ii2 <- lapply(1:NN, function(i) inits())
## for(i in 1:NN){
##   for(j in c(1,2,5,6)){
##     ii[[i]][[j]] <- log(ii[[i]][[j]])
##   }
## }
## ii <- lapply(1:NN, function(i) ii[[1]])
## ii2 <- lapply(1:NN, function(i) ii2[[1]])
## fit.tmb <-
##   sample_tmb(obj=obj.tmb, iter=1100, warmup=1000, chains=NN, seeds=1:15,
##              init=ii, control=list(metric=NULL, adapt_delta=.9,
##                           max_treedepth=8, stepsize=.001, adapt_mass=TRUE))
## sp1 <- extract_sampler_params(fit.tmb, TRUE)
## ggplot(subset(sp, chain!=11), aes(iteration, energy__, group=chain,
##                                   color=factor(chain), size=divergent__)) +
##   geom_point() + geom_line() + xlim(0,300)
## ggplot(subset(sp, chain!=11), aes(iteration, divergent__, group=chain, color=factor(chain))) +
##   geom_line() + xlim(0,25)
## plot(sp$iteration, log(sp$stepsize__), type='l'); abline(v=c(125,175,475))
## plot(sp$iteration, sp$energy__, type='l'); abline(v=c(125,175,475))
## ss <- extract_samples(fit.tmb, TRUE, TRUE)
## which(ss$lp__-sp$energy__>0)

## fit.admb <-
##   sample_admb('growth', 'admb',  iter=1100, warmup=1000, chains=15, seeds=1:15,
##              init=ii, control=list(metric=NULL, adapt_delta=.9,
##                           max_treedepth=8, adapt_mass=TRUE))
## sp2 <- extract_sampler_params(fit.admb, TRUE)
## plot(sp$iteration, log(sp$stepsize__), type='l'); abfline(v=c(125,175,475))
## plot(sp$iteration, sp$energy__, type='l'); abline(v=c(125,175,475))
## ss <- extract_samples(fit.admb, TRUE, TRUE)
## which(ss$lp__-sp$energy__>0)

## fit.stan <- sampling(obj.stan, data=data,  iter=1100, warmup=1000, chains=15, init=ii2,
##                  control=list(adapt_delta=.9, max_treedepth=8))
## sp3 <- data.frame(chain=rep(1:15, each=1100), iteration=1:1100, do.call(rbind,
##                  get_sampler_params(fit.stan, inc_warmup=TRUE)))
## sp3$energy__ <- sp3$energy__*-1

## sp.all <- rbind(cbind(alg='tmb', sp1), cbind(alg='admb', sp2),
##                 cbind(alg='stan', sp3))
## ggplot(subset(sp.all, chain!=11), aes(iteration, energy__, group=chain, color=factor(chain))) +
##   geom_line() + xlim(150,750) + facet_wrap("alg")

## ## Simulated spatial model, TMB example
## m <- 'spatial'
## temp <- spatial_setup()
## data <- temp$data; inits <- temp$inits
## lower <- abs(unlist(inits()))*-Inf
## upper <- abs(unlist(inits()))*Inf
## lower['sigma'] <- 0; lower['a'] <- 0
## upper['a'] <- 2
## setwd('models/spatial')
## compile('spatial.cpp')
## dyn.load('spatial')

## obj <- MakeADFun(data=data, parameters=inits(), DLL = "spatial", random = "u")
## opt <- nlminb(obj$par, obj$fn, obj$gr,
##               lower=c(-100.0, -100.0, 0.01, -3.0),
##               upper=c( 100.0,  100.0, 3.00,  3.0) )
## test <- tmbstan(obj,  chains=3, iter=700,
##                 init=lapply(1:3, function(i) inits()))
