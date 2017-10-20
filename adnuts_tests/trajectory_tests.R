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

devtools::load_all("c:/Users/Cole/adnuts")

## !!!! This requires to use the test nuts function in ADMB -nuts_test and
## !!!! NOT -nuts. This special fucntion writes the first directory to file
## !!!! so it can be read into R. Also you need to use the test_trajectory
## !!!! branch of adnuts which has the run_mcmc.nuts.test function.

## Generate single trajectory and make sure identical path and acceptance
## probabilities. Start with basic parabola as easiest case.
## set.seed(200) # 200 gives a bad trajectory
set.seed(546)
out.list <- list()
k <- 1
for(i in 1:15){
r0 <- rnorm(2)
x0 <- rnorm(2)
for(eps in c(.01, .05)){
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
trash <- lapply(file, function(f) if(file.exists(f)) file.remove(f))
temp <- run_mcmc.nuts.test(iter=1, warmup=1, fn=fn2, gr=gr2, init=x0, r0=r0,
                               chain=1, control=list(stepsize=eps))
system(paste0('parabola -noest -mcmc 1 -nuts_test -mcdiag -hyeps ', eps, ' -mcseed ',10), ignore.stdout=FALSE)
out1.admb <- read.table(file[1], sep=' ', row.names=NULL, header=TRUE)[,-1]
out2.admb <- cbind(read.table(file[2], sep=' ', row.names=NULL, header=TRUE))
out1.tmb <- data.frame(temp$out, row.names=NULL)
out2.tmb <- data.frame(temp$out2)
names(out2.tmb) <- names(out2.admb)
names(out1.tmb) <- names(out1.admb)
td <- c(-1,rep(0:10, times=2^(0:10)))
out1.tmb$treedepth <- factor(td[1:nrow(out1.tmb)])
out1.admb$treedepth <- factor(td[1:nrow(out1.admb)])
rm(out1, out2)
out1 <- rbind(cbind(m='admb', out1.admb), cbind(m='tmb', out1.tmb))
out2 <- rbind(cbind(m='admb', out2.admb), cbind(m='tmb', out2.tmb))
## plot(out1.tmb$alpha, out1.tmb$alpha-out1.admb$alpha)
out.list[[k]] <- cbind(eps=eps, rep=i, out1)
setwd('../../')
k <- k+1
}}
out.all <- do.call(rbind, out.list)
g <- ggplot(data=out.all, aes(theta1, theta2, group=rep, color=treedepth)) + facet_grid(m~eps)
g <- g+ geom_point(size=.1, alpha=1)
ggsave('plots/paraboloa_trajectories.png', g, width=7, height=5)
g <- ggplot(data=out.all, aes(theta1, log(alpha), color=rep)) + geom_point(alpha=.5) + facet_wrap('m')
ggsave('plots/paraboloa_alpha.png', g, width=7, height=5)


## Same parabola but this time bounded tightly
system('admb parabola_bounded'); system('parabola_bounded')
compile(file='parabola_bounded_tmb.cpp')
set.seed(546)
out.list <- list()
k <- 1
for(i in 1:15){
r0 <- rnorm(2)
x0 <- rnorm(2) # unbounded!
for(eps in c(.01, .05)){
setwd('models/parabola_bounded/')
dyn.load(dynlib('parabola_bounded_tmb'))
data <- list(Npar=2)
obj <- MakeADFun(data=data, parameters=list(x1=0, x2=0), DLL='parabola_bounded_tmb')
obj$env$beSilent()
lower <- c(-2,-1); upper <- c(1,2)
cases <- .transform.cases(lower,upper)
fn2 <- function(y){
  x <- .transform(y, lower, upper, cases)
  scales <- .transform.grad(y, lower, upper, cases)
  -obj$fn(x) + sum(log(scales))
}
gr2 <- function(y){
  x <- .transform(y, lower, upper, cases)
  scales <- .transform.grad(y, lower, upper, cases)
  scales2 <- .transform.grad2(y, lower, upper, cases)
  -as.vector(obj$gr(x))*scales + scales2
}
write.table(c(x0, r0, 10), file='input.txt', sep=' ', row.names=F, col.names=F)
## devtools::load_all("c:/Users/Cole/adnuts")
file <- c('trajectory.txt', 'out2.txt')
trash <- lapply(file, function(f) if(file.exists(f)) file.remove(f))
temp <- run_mcmc.nuts.test(iter=1, warmup=1, fn=fn2, gr=gr2, init=x0, r0=r0,
                               chain=1, control=list(stepsize=eps))
system(paste0('parabola_bounded -mcscale 0 -noest -mcmc 1 -nuts_test -mcdiag -hyeps ', eps, ' -mcseed ',10), ignore.stdout=FALSE)
out1.admb <- read.table(file[1], sep=' ', row.names=NULL, header=TRUE)[,-1]
out2.admb <- cbind(read.table(file[2], sep=' ', row.names=NULL, header=TRUE))
out1.tmb <- data.frame(temp$out, row.names=NULL)
out2.tmb <- data.frame(temp$out2)
names(out2.tmb) <- names(out2.admb)
names(out1.tmb) <- names(out1.admb)
td <- c(-1,rep(0:10, times=2^(0:10)))
out1.tmb$treedepth <- factor(td[1:nrow(out1.tmb)])
out1.admb$treedepth <- factor(td[1:nrow(out1.admb)])
rm(out1, out2)
out1 <- rbind(cbind(m='admb', out1.admb), cbind(m='tmb', out1.tmb))
out2 <- rbind(cbind(m='admb', out2.admb), cbind(m='tmb', out2.tmb))
## plot(out1.tmb$alpha, out1.tmb$alpha-out1.admb$alpha)
out.list[[k]] <- cbind(eps=eps, rep=i, out1)
setwd('../../')
k <- k+1
}}
## Filter out rejected points
out.all <- do.call(rbind, out.list)
## bound them manually
for(i in 1:nrow(out.all))
  out.all[i, 4:5] <- .transform(as.numeric(out.all[i, 4:5]), lower, upper, cases)
g <- ggplot(data=out.all, aes(theta1, theta2, group=rep, color=treedepth)) + facet_grid(m~eps)
g <- g+ geom_point(size=.1, alpha=1)
ggsave('plots/paraboloa_bounded_trajectories.png', g, width=7, height=5)
g <- ggplot(data=out.all, aes(theta1, log(alpha), color=rep)) + geom_point(alpha=.5) + facet_wrap('m')
ggsave('plots/paraboloa_bounded_alpha.png', g, width=7, height=5)

## Same thing with funnel
set.seed(546)
out.list <- list()
k <- 1
for(i in 1:15){
r0 <- rnorm(2, sd=c(1,.1))#c(0,0)#c(2,.1)
x0 <- rnorm(2)
  for(eps in c(.01, .05)){
setwd('models/funnel/')
#stan.fit <- stan(file='funnel.stan', data=list())
## compile(file='funnel.cpp')
dyn.load(dynlib('funnel'))
obj <- MakeADFun(data=list(), parameters=list(v=0, theta=0), DLL='funnel')
obj$env$beSilent()
fn2 <- function(x) -obj$fn(x)
gr2 <- function(x) -as.vector(obj$gr(x))
setwd('admb')
write.table(c(x0, r0, 15), file='input.txt', sep=' ', row.names=F, col.names=F)
## system('admb funnel'); system('funnel')
file <- c('trajectory.txt', 'out2.txt')
## trash <- lapply(file, function(f) if(file.exists(f)) file.remove(f))
temp <- run_mcmc.nuts.test(iter=1, warmup=1, fn=fn2, gr=gr2, init=x0,
                           r0=r0, chain=1, control=list(stepsize=eps, max_treedepth=12))
system(paste0('funnel -max_treedepth 12 -noest -mcmc 1 -mcdiag -nuts_test -hyeps ', eps, ' -mcseed ',10), ignore.stdout=FALSE)
out1.admb <- cbind(read.table(file[1], sep=' ', row.names=NULL, header=TRUE)[,-1])
out2.admb <- cbind(read.table(file[2], sep=' ', row.names=NULL, header=TRUE))
out1.tmb <- data.frame(temp$out, row.names=NULL)
out2.tmb <- data.frame(temp$out2)
names(out2.tmb) <- names(out2.admb)
names(out1.tmb) <- names(out1.admb)
td <- c(-1,rep(0:15, times=2^(0:15)))
out1.tmb$treedepth <- factor(td[1:nrow(out1.tmb)])
out1.admb$treedepth <- factor(td[1:nrow(out1.admb)])
rm(out1, out2)
out1 <- rbind(cbind(m='admb', out1.admb), cbind(m='tmb', out1.tmb))
out2 <- rbind(cbind(m='admb', out2.admb), cbind(m='tmb', out2.tmb))
out.list[[k]] <- cbind(eps=eps, rep=i, out1)
setwd('../../..')
k <- k+1
}}
out.all <- do.call(rbind, out.list)
g <- ggplot(data=out.all, aes(theta1, theta2, group=rep, color=treedepth)) + facet_grid(m~eps)
g <- g+ geom_point(size=.1, alpha=1)
ggsave('plots/funnel_trajectories.png', g, width=7, height=5)
g <- ggplot(data=out.all, aes(theta1, log(alpha), color=rep)) + geom_point(alpha=.5) + facet_wrap('m')
ggsave('plots/funnel_alpha.png', g, width=7, height=5)




