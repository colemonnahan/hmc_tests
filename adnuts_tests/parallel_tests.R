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
setwd('admb')
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
## system('admb mvnd')
system('mvnd -hbf 1')
setwd('..')

devtools::install("c:/Users/cole/adnuts")

alg <- "RWM"; system('admb/mvnd')
alg <- "NUTS"; system('admb/mvnd -hbf 1')
control <- list(metric=covar)
## Make sure the seeds work across algorithms and parallel and serial.
inits <- lapply(1:2, function(i)  list(mu=c(0,0)))
admb.para <- sample_admb('mvnd', iter=2000, dir='admb', parallel=TRUE, chains=2,
                         cores=3, init=inits, algorithm=alg, seeds=1:2,
                         control=control)
admb.serial <- sample_admb('mvnd', iter=2000, dir='admb', parallel=FALSE, chains=2,
                           cores=3, init=inits, algorithm=alg, seeds=1:2,
                           control=control)
init <- function() list(mu=c(0,0))
tmb.para <- sample_tmb(obj=mvnd.obj, iter=200, init=init, chains=2,
                       parallel=TRUE, cores=2, seeds=2:3, algorithm=alg,
                       control=control, lower=c(-1,-1))
tmb.serial <- sample_tmb(obj=mvnd.obj, iter=200, init=init, chains=2,
                         parallel=FALSE, cores=2, seeds=2:3, algorithm=alg,
                         control=control, lower=c(-1,-1))
## Make sure samples work with same seeds; these should match
tail(extract_samples(admb.para)[,1,],1)
tail(extract_samples(admb.serial)[,1,],1)
tail(extract_samples(tmb.para)[,1,],1)
tail(extract_samples(tmb.serial)[,1,],1)

### Use this to test the function outside of the package. Useful for
### debugging
###
devtools::install("c:/Users/cole/adnuts")
inits <- lapply(1:3, function(i)  list(mu=rnorm(2)))
sfInit(parallel=TRUE, slaveOutfile='out.txt', cpus=2)
sfLibrary(TMB)
rm(test)
sfExportAll()
test <- sfLapply(1:2, function(i)
test <- lapply(1:2, function(i)
  adnuts:::sample_tmb_parallel(i, obj=mvnd.obj, init= inits[[i]], seed=1,
                               path=getwd() , algorithm='NUTS', lower=lower.bounds,
                               upper=NULL, iter=100 ))



