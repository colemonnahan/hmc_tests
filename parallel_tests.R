## Get parallel working with TMB
## devtools::install("../adnuts")
source("startup.R")
set.seed(235)
Npar <- 2
corr <- matrix(.8, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
#covar <- diag(2)
covar2 <- diag(x=diag(covar))
covar.inv <- solve(covar)
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
compile(file='mvnd.cpp')
dyn.load(dynlib('mvnd'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd')
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')

library(snowfall)
library(adnuts)
sfInit(parallel=TRUE, cpus=6)
sfExportAll()
init <- lapply(1:6, function(i)  list(mu=c(0,0)))
admb <- sample_admb(model, iter=2000, dir=dir, parallel=TRUE, chains=6,
                    cores=3, init=inits, algorithm='RWM')

devtools::install("c:/Users/cole/adnuts")
devtools::document("c:/Users/cole/adnuts")
sfInit(parallel=TRUE, cpus=6)
sfLibrary(TMB)
sfLibrary(adnuts)
sfExportAll()
tmb <- sample_tmb(obj=mvnd.obj, iter=2000, init=inits, chains=6)


##mcmc.out <- lapply(1:2, function(i)
mcmc.out <- sfLapply(1:2, function(i)
  adnuts:::sample_tmb_parallel(parallel_number=i))

sfLibrary(TMB)
mcmc.out <- sfLapply(1:6, function(i){
  dir='C:/Users/Cole/hmc_tests/models/mvnd'
  setwd(dir)
  dyn.load(dynlib(model))
  obj <- MakeADFun(data=inputs$data,
                   parameters=list(mu=rnorm(2)))
  for(i in 1:100000) obj$gr()
  return(obj$fn())
})

## rm(list=ls())
## sample_tmb_parallel2 <-  function(parallel_number, inputs, init, dir, model){
##   setwd(dir)
##   dyn.load(dynlib(model))
##   obj <- MakeADFun(data=inputs$data,
##                    parameters=init)
##   return(data.frame(chain=parallel_number, x=obj$fn()))
## }
inits <- lapply(1:6, function(i)  list(mu=rnorm(2)))
dir <- 'C:/Users/Cole/hmc_tests/models/mvnd'
inputs <- list(data=list(covar=diag(2), Npar=2, x=c(0,0)))
model <- 'mvnd'
sfInit(parallel=TRUE, slaveOutfile='out.txt', cpus=6)
sfLibrary(TMB)
sfExportAll()
test <- sfLapply(x=1:6, fun=function(i)
  sample_tmb_parallel2(i, inputs, init= inits[[i]], dir, model, iter=200000))



