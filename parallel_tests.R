## Get parallel working with TMB
## devtools::install("../adnuts")
source("startup.R")
library(snowfall)
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
## setwd('admb')
## write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
## system('admb mvnd')
## system('mvnd')
## setwd('..')

sfInit(parallel=TRUE, cpus=2)
inits <- lapply(1:2, function(i)  list(mu=c(0,0)))
admb <- sample_admb(model, iter=2000, dir=dir, parallel=TRUE, chains=2,
                    cores=3, init=inits, algorithm='RWM')


init <- NULL
init <- function() list(mu=rnorm(2))
devtools::install("c:/Users/cole/adnuts")
tmb <- sample_tmb(obj=mvnd.obj, iter=205, init=init, chains=2,
                  parallel=TRUE, cores=2)


## devtools::install("c:/Users/cole/adnuts")
## inits <- lapply(1:3, function(i)  list(mu=rnorm(2)))
## sfInit(parallel=TRUE, slaveOutfile='out.txt', cpus=2)
## sfLibrary(TMB)
## rm(test)
## sfExportAll()
## test <- sfLapply(1:2, function(i)
## ##test <- lapply(1:2, function(i)
##   adnuts:::sample_tmb_parallel(i, obj=mvnd.obj, init= inits[[i]],
##                                dir=getwd() , algorithm='NUTS', lower=NULL,
##                                upper=NULL, iter=2000 ))



