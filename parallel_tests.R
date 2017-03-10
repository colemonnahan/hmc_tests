library(snowfall)
set.seed(235)
Npar <- 2
corr <- matrix(.98, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
covar <- diag(2)
covar2 <- diag(x=se^2)
covar.inv <- solve(covar)
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
compile(file='mvnd.cpp')
dyn.load(dynlib('mvnd'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')

devtools::document('C:/Users/Cole/adnuts')

devtools::load_all('C:/Users/Cole/adnuts')

devtools::install('C:/Users/Cole/adnuts')
library(adnuts)
cores <- 4
sfStop()
sfInit(parallel=TRUE, cpus=cores)
sfExportAll()
chains <- 3
out <- sample_admb(dir='admb', model='mvnd', iter=2000,
                   init=rep(list(c(0,0)), chains), chains=chains,
                   parallel=TRUE, cores=3, control=list(algorithm="RWM"))


xx <- sfLapply(1:8, sample_admb_parallel, dir='admb', iter=4000,
               model='mvnd', init=rep(list(c(0,0)),1),
               control=list(algorithm='NUTS', thin=2))
yy <- combine_fits(xx)
out <- sample_admb_parallel(1, dir=dir, model='mvnd', iter=2000, init=list(rnorm(2)))
out <- sample_admb(chains=2, dir=dir, model='mvnd', mceval=TRUE, iter=2000,
                   init=rep(list(rnorm(2)),2), control=list(algorithm='RWM'))
