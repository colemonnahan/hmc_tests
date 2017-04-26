## Quick code to test that TMB and ADMB generate identical trajectories
## with the same model and initial theta and r values. Setup to work with
## the funnel model since it has difficult properties
rm(list=ls())
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
library(adnuts)
##!!!!! This needs to be the test_trajectory branch!!!!

##devtools::install("c:/Users/Cole/adnuts")

## !!!! THIS REQUIRES COMPILING A SPECIAL TESTING BRANCH OF ADMB !!!!
## Use: trajectory_tests. OTherwise this wont work at all. Also need to use
## the adnuts branch test_trajectory.

## Generate single trajectory and make sure identical path and acceptance
## probabilities. Start with basic parabola as easiest case.
#set.seed(200) # 200 gives a bad trajectory
x0 <- rnorm(2)
r0 <- rnorm(2)
eps <- .8
setwd('models/parabola/')
compile(file='parabola_tmb.cpp')
dyn.load(dynlib('parabola_tmb'))
data <- list(Npar=2)
obj <- MakeADFun(data=data, parameters=list(x1=0, x2=0), DLL='parabola_tmb')
obj$env$beSilent()
fn2 <- function(x) -obj$fn(x)
gr2 <- function(x) -as.vector(obj$gr(x))
write.table(c(x0, r0, 10), file='input.txt', sep=' ', row.names=F, col.names=F)
write.table(x=2, file='parabola.dat', row.names=FALSE, col.names=FALSE)
## devtools::load_all("c:/Users/Cole/adnuts")
## system('admb parabola'); system('parabola')
file <- c('trajectory.txt', 'out2.txt')
## trash <- lapply(file, function(f) if(file.exists(f)) file.remove(f))
temp <- run_mcmc.nuts.test(iter=1, warmup=1, fn=fn2, gr=gr2, init=x0, r0=r0,
                               chain=1, control=list(stepsize=eps))
system(paste0('parabola -noest -mcmc 1 -nuts -mcdiag -hyeps ', eps, ' -mcseed ',10), ignore.stdout=FALSE)
out1.admb <- read.table(file[1], sep=' ', row.names=NULL, header=TRUE)[,-1]
out2.admb <- cbind(read.table(file[2], sep=' ', row.names=NULL, header=TRUE))
out1.tmb <- data.frame(temp$out, row.names=NULL)
out2.tmb <- data.frame(temp$out2)
names(out2.tmb) <- names(out2.admb)
names(out1.tmb) <- names(out1.admb)
out1.tmb$treedepth <- c(-1, factor(rep(0:10, times=2^(0:10))))[1:nrow(out1.tmb)]
out1.admb$treedepth <- c(-1, factor(rep(0:10, times=2^(0:10))))[1:nrow(out1.admb)]
rm(out1, out2)
out1 <- rbind(cbind(m='admb', out1.admb), cbind(m='tmb', out1.tmb))
out2 <- rbind(cbind(m='admb', out2.admb), cbind(m='tmb', out2.tmb))
g <- ggplot(data=out1) + facet_wrap('m') + xlim(-3,3) + ylim(-3,3)
g+ geom_point(aes(theta1, theta2, color=factor(treedepth)))
## g+ geom_point(aes(thetaminus1, thetaminus2, color=treedepth))
## g+ geom_point(aes(thetaplus1, thetaplus2, color=treedepth))
plot(out1.tmb$alpha, out1.tmb$alpha-out1.admb$alpha)
setwd('../../')



## Generate single trajectory and make sure identical path and acceptance
## probabilities. Start with basic parabola as easiest case.
#set.seed(6)
setwd('models/funnel/')
#stan.fit <- stan(file='funnel.stan', data=list())
compile(file='funnel.cpp')
dyn.load(dynlib('funnel'))
x0 <- c(0,0)#c(2,.1)
r0 <- rnorm(2)
eps <- .1
obj <- MakeADFun(data=list(), parameters=list(v=0, theta=0), DLL='funnel')
obj$env$beSilent()
fn2 <- function(x) -obj$fn(x)
gr2 <- function(x) -as.vector(obj$gr(x))
setwd('admb')
write.table(c(x0, r0, 10), file='input.txt', sep=' ', row.names=F, col.names=F)
 system('admb funnel'); system('funnel')
file <- c('trajectory.txt', 'out2.txt')
## trash <- lapply(file, function(f) if(file.exists(f)) file.remove(f))
temp <- run_mcmc.nuts.test(iter=1, warmup=1, fn=fn2, gr=gr2, init=x0, r0=r0, chain=1, control=list(stepsize=eps))
## Be careful here with admodel.cov beine right!!
# write.admb.cov(diag(2))
system(paste0('funnel -noest -mcmc 1 -mcdiag -nuts -hyeps ', eps, ' -mcseed ',10), ignore.stdout=FALSE)
out1.admb <- cbind(read.table(file[1], sep=' ', row.names=NULL, header=TRUE)[,-1])
out2.admb <- cbind(read.table(file[2], sep=' ', row.names=NULL, header=TRUE))
out1.tmb <- data.frame(temp$out, row.names=NULL)
out2.tmb <- data.frame(temp$out2)
names(out2.tmb) <- names(out2.admb)
names(out1.tmb) <- names(out1.admb)
out1.tmb$treedepth <- c(-1, factor(rep(0:10, times=2^(0:10))))[1:nrow(out1.tmb)]
out1.admb$treedepth <- c(-1, factor(rep(0:10, times=2^(0:10))))[1:nrow(out1.admb)]
rm(out1, out2)
out1 <- rbind(cbind(m='admb', out1.admb), cbind(m='tmb', out1.tmb))
out2 <- rbind(cbind(m='admb', out2.admb), cbind(m='tmb', out2.tmb))
g <- ggplot(data=out1) + facet_wrap('m') + xlim(-3,3) + ylim(-3,3)
g+ geom_point(aes(theta1, theta2, color=factor(treedepth)))
## g+ geom_point(aes(thetaminus1, thetaminus2, color=treedepth))
## g+ geom_point(aes(thetaplus1, thetaplus2, color=treedepth))
#plot(out1.tmb$alpha, out1.tmb$alpha-out1.admb$alpha)
setwd('..')




