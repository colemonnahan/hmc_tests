## Quick code to test TMB divergences from NUTS
library(TMB)
library(rstan)
library(plyr)
library(mvtnorm)
library(shinystan)
library(adnuts)


## ## Rosenbrock is a classic difficult one to fit and useful for checking
## ## divergences. Turns out this one is too difficult!! Not really that
## ## informative for these tests.

## ## ADMB wsa being weird with the psv file so I changed it's name
## ## temporarily to rosen.
## setwd('models/rosenbrock/')
## compile(file='rosenbrock.cpp')
## dyn.load(dynlib('rosenbrock'))
## data <- list()
## pars <- c('x1', 'x2')
## obj <- MakeADFun(data=data, parameters=list(x1=0, x2=0))
## ## admb model compiled manually
## temp <- list()
## k <- 1
## for(eps in c(.05, .1)){
##   iter <- 10000; warmup <- iter/2

## fit.tmb <- sample_tmb(obj, iter=iter, init=list(c(0,0)), warmup=warmup,
##                       control=list(metric=diag(2), stepsize=eps))
## fit.admb <- sample_admb(dir='admb', model='rosen', iter=iter,
##                         init=list(list(c(0,0))), warmup=warmup,
##                       control=list(metric=diag(2), stepsize=eps))
## fit.stan <- stan(file='rosenbrock.stan', iter=iter, warmup=warmup,
##                  init=list(list(x1=0, x2=0)), chains=1,
##                  control=list(metric='unit_e', stepsize=eps, adapt_engaged=FALSE))
## ## Make pairs plots with divergences
## x1 <- data.frame(m='tmb',eps=eps, extract_samples(fit.tmb), div=extract_sampler_params(fit.tmb)[,'divergent__'])
## x2 <- data.frame(m='admb',eps=eps, extract_samples(fit.admb), div=extract_sampler_params(fit.admb)[,'divergent__'])
## x3 <- data.frame(m='stan',eps=eps, extract(fit.stan), div=get_sampler_params(fit.stan)[[1]][-(1:warmup),'divergent__'])
## temp[[k]] <- rbind.fill(x1, x2,x3)
## k <- k +1
## }
## res <- do.call(rbind, temp)
## ggplot(res, aes(x1, x2, col=div==1)) + geom_point(alpha=.1) + facet_grid(m~eps)
## setwd('../..')





## Neal's funnel, copied from Betancourt and Girolami (2015) HMC for
## hierarchical models
setwd('models/funnel/')
compile(file='funnel.cpp')
dyn.load(dynlib('funnel'))
data <- list()
pars <- c('v', 'theta')
obj <- MakeADFun(data=data, parameters=list(v=0, theta=0), DLL='funnel')
## admb model compiled manually
temp <- list(); k <- 1
for(eps in c(.1, .3, .7)){
iter <- 5000; warmup <- 50
fit.tmb <- sample_tmb(obj, iter=iter, init=list(c(0,0)), warmup=warmup,
                      control=list(metric=diag(2), stepsize=eps))
fit.admb <- sample_admb(dir='admb', model='funnel', iter=iter,
                        init=list(c(0,0)),
                        warmup=warmup,
                      control=list(metric=diag(2), stepsize=eps))
fit.stan <- stan(file='funnel.stan', iter=iter, warmup=warmup,
                 init=list(list(v=0, theta=0)), chains=1,
                 control=list(metric='unit_e', max_treedepth=10, stepsize=eps, adapt_engaged=FALSE))
## Make pairs plots with divergences
x1 <- data.frame(m='tmb',eps=eps, extract_samples(fit.tmb), div=extract_sampler_params(fit.tmb)[,'divergent__'])
x2 <- data.frame(m='admb',eps=eps, extract_samples(fit.admb), div=extract_sampler_params(fit.admb)[,'divergent__'])
x3 <- data.frame(m='stan',eps=eps, extract(fit.stan), div=get_sampler_params(fit.stan, inc_warmup=FALSE)[[1]][,'divergent__'])
temp[[k]] <- rbind.fill(x1, x2,x3)
k <- k +1
}
res <- do.call(rbind, temp)
res$eps <- paste0("stepsize=",res$eps)
res <- res[order(res$div),]
setwd('../..')
g <- ggplot(res, aes(v, theta, col=div==1)) + geom_point(alpha=.5, size=.2) + facet_grid(m~eps)
ggsave('plots/divergence_tests.png',g, width=7, height=6)


### OLD code used during development of TMB. Probably not useful anymore
## #but leaving here for now
## #
## data <- list(Y = c(28,  8, -3,  7, -1,  1, 18, 12),
##                     sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
## ## dyn.unload(dynlib('rosenbrock'))
## compile(file='8schools.cpp')
## dyn.load(dynlib('8schools'))
## obj <- MakeADFun(data=data, parameters=list(mu=1, logtau=-5, theta=rep(0,8)))
## out <- run_mcmc(obj, nsim=2000, alg='NUTS', eps=NULL, chains=3)
## sso.tmb <- with(out, as.shinystan(samples, burnin=warmup, max_treedepth=12,
##              sampler_params=sampler_params, algorithm='NUTS'))
## launch_shinystan(sso.tmb)


## launch_shinystan_tmb(fit.tmb)


## plot(0,0, type='n', xlim=c(-1.5,1.5), ylim=c(-.5,3))

## rm(theta.trajectory)
## theta0 <- theta <- c(0,0)
## r0 <- r <- rnorm(2)
## v <- 1
## u <- TMB:::.sample.u(theta=theta, r=r, fn=obj$fn)
## j <- 12
## divergent <<- 0
## info <- as.environment(list(n.calls=0))
## xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j, eps=.1,
##                  theta0=theta0, r0=r0, fn=fn2, gr=gr2, info=info)
## c(info$n.calls, divergent)
## points(theta[1], theta[2], pch=16)
## points(theta.trajectory[,2], theta.trajectory[,3], type='b', xlim=c(-3,3), ylim=c(-3,3))
## with(xx, points(theta.prime[1], theta.prime[2], pch=16, col='red', cex=.5))


## with(xx, (theta.plus-theta.minus) %*% r.minus)
## drop(with(xx, crossprod(theta.plus-theta.minus, r.minus)))
## cbind(step=0, t(c(1,1)))

## xx <- t(sapply(1:5000, function(i) .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=.3,
##                  theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)$theta.prime))
## barplot(table(xx[,1]))
