library(snowfall)
set.seed(235)
Npar <- 10
corr <- matrix(.98, nrow=Npar, ncol=Npar)
diag(corr) <- 1
se <- runif(Npar, .1, 5)
covar <- corr * (se %o% se)
#covar <- diag(Npar)
covar2 <- diag(x=se^2)
covar.inv <- solve(covar)
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
compile(file='mvnd.cpp')
dyn.load(dynlib('mvnd'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=rep(0,Npar)), DLL='mvnd')
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')


#rm(list=ls())
devtools::install('C:/Users/Cole/adnuts')
library(adnuts)
Npar <- 10
cores <- 3
chains <- 3
ww <- 1000 # warmup
dd <- .5 # duration
pp <- TRUE
ii <- 10000
sfStop()
sfInit(parallel=TRUE, cpus=cores)
sfExportAll()


out.nuts1 <- sample_admb(dir='admb', model='mvnd', iter=ii, duration=dd, thin=5,
                   init=rep(list(rep(0,Npar)), chains), chains=chains, warmup=ww,
                   parallel=pp,  control=list(algorithm="NUTS"))
out.rwm1 <- sample_admb(dir='admb', model='mvnd', iter=10*ii, duration=dd,
                        init=rep(list(rep(0,Npar)), chains), chains=chains, thin=50,
                        warmup=10*ww, parallel=pp,
                        control=list(algorithm="RWM"))

sfStop()
launch_shinystan_admb(out.nuts1)
launch_shinystan_admb(out.rwm1)

## lps of both
x <- as.vector(out.nuts1$samples[-(1:500),,11])
y <- as.vector(out.rwm1$samples[-(1:500),,11])
true <- mvtnorm::rmvnorm(n=length(x), sigma=covar)
z <- mvtnorm::dmvnorm(x=true, sigma=covar, log=TRUE)
par(mfrow=c(1,3))
qqplot(x,y)
qqplot(z,x)
qqplot(z,y)
