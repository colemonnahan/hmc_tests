library(TMB)
library(plyr)
library(coda)
devtools::install_github('bbolker/R2admb/R2admb')
library(R2admb)
library(rstan)
library(shinystan)
devtools::load_all("c:/Users/Cole/rnuts")
## devtools::install("c:/Users/Cole/rnuts")
## library(rnuts)
## devtools::document("c:/Users/Cole/admbtools")
## devtools::install("c:/Users/Cole/admbtools")
## devtools::load_all("c:/Users/Cole/admbtools")
## library(admbtools)

## build_tree tests
d <- read.csv('build_tree.csv', header=F)[-1,]
names(d) <- c('x1', 'x2', 'y1', 'y2', 'z1', 'z2')
par(mfrow=c(1,3))
f <- function() points(0,0, pch=16)
plot(d$x1, d$x2, type='b'); f()
plot(d$y1, d$y2, type='b'); f()
plot(d$z1, d$z2, type='b'); f()



## Super quick ADMB tests.
model.path="C:/Users/Cole/hmc_tests/models/catage"
model.name='catage'
x <- run_admb_mcmc(model.path=model.path, model.name=model.name, iter=1000,
                   chains=3, eps=.2, max_treedepth=14)
launch_shinystan_admb(x)
setwd(model.path)
system('admb catage')
system('catage -nohess -mcmc 10 -nuts -mcseed 5')
adapt <- read.csv("adaptation.csv")
pars <- read_psv(model.name)
fit <- read_admb(model.name)

covar <- matrix(.954, nrow=2, ncol=2)
diag(covar) <- 1
covar.inv <- solve(covar)
Npar <- 2
Niter <- 1000

setwd('C:/Users/Cole/hmc_tests/')
dyn.unload(dynlib('models/mvnd/mvnd_tmb'))
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')

init <- lapply(1:1, function(x) rnorm(2, sd=10))
Niter <- 1000
nuts <- run_mcmc(mvnd.obj, iter=2*Niter, init=init, chains=1, alg="NUTS",
                 covar=covar, thin=2)
hmc <- run_mcmc(mvnd.obj, iter=2*Niter, init=init, L=10, chains=3, alg="HMC",
                covar=covar, thin=2)
rwm <- run_mcmc(mvnd.obj, iter=10*Niter, init=init, chains=3, alg="RWM",
                thin=10, alpha=.2, covar=covar)
sso.nuts <- as.shinystan.tmb(nuts)
launch_shinystan(sso.nuts)
sso.hmc <- as.shinystan.tmb(hmc)
launch_shinystan(sso.hmc)
sso.rwm <- as.shinystan.tmb(rwm)
launch_shinystan(sso.rwm)

## Some quick profiling
Rprof()
nuts <- run_mcmc(mvnd.obj, iter=2*Niter, init=init, chains=3, alg="NUTS",
                 covar=diag(2), thin=2, upper=c(4,4), lower=c(-4,-4))
Rprof(NULL)
head(summaryRprof()$by.self)


xx <- extract_samples(nuts)
plot(xx)


## Quick experiment for looking at typical set from iid Z as dimension
## grows.
n <- 5000
xx <- ldply(2^(1:10), function(D){
  x <- matrix(rnorm(D*n), ncol=D)
  dist <- apply(x, 1, function(y) (sum(y^2)))
  cbind(D=D, IQR=IQR(dist))
}
)
plot(xx)


m <- 'simple'
m <- 'catage'
setwd(main.dir)
setwd(file.path('models', m))
## compile and run it
system(paste('admb_original',m))
system(m)
## res <- read_pars(m)
## par.names <- names(res$coefficients[1:res$npar])

L.vec <- floor(seq(2, 50, len=10))
eps.vec <- exp(seq(log(.01), log(.2), len=10))
iter <- 2e3
seed <- sample(1:1e3, size=1)
results <- ldply(L.vec, function(L) {
    ldply(eps.vec, function(eps){
        print(paste(eps, L, sep=";"))
        if(file.exists(file.path(m,'.psv')))
            file.remove(file.path(m, '.psv'))
        time <- system.time(
            system(paste(m, '-noest -hybrid -mcmc', iter, '-hynstep', L, '-hyeps',
                         eps, '-mcseed', seed), ignore.stdout=TRUE))['elapsed']
        Sys.sleep(.1)
        psv <- read_psv(m, names=NULL)
        ## psv.array <- array(NA, dim=c(nrow(psv),1,ncol(psv)))
        ## rstan::monitor(sims=psv, warmup=0, probs=.5)
        ess <- effectiveSize(psv)
        accept.rate <- 1- sum(psv[-nrow(psv),1]== psv[-1,1])/iter
        x <- cbind(L, eps, accept.rate, min.ess=min(ess), time=as.numeric(time))
        x})})
results <- within(results, perf <- min.ess/time)

ggplot(results, aes(L, time, group=eps)) + geom_line()
ggplot(results, aes(eps, min.ess, group=L, color=factor(L))) + geom_line()
ggplot(results, aes(eps, perf, group=L, color=factor(L))) + geom_line()
ggplot(results, aes(L, eps, size=perf)) + geom_point()
ggplot(results, aes(eps, accept.rate, group=L, color=factor(L))) +
    geom_line() + ylim(0,1)


xx <- ldply(c(.8, .85, .9, .95), function(delta){
  ldply(1:3, function(seed){
  system(paste('catage -noest -mcmc 2000 -hmc -hyeps 1 -adapt_delta',
               delta, '-mcseed', seed),
         ignore.stdout=TRUE, ignore.stderr=TRUE )
  system('catage -mceval')
  cbind(seed=seed, delta=delta, read.csv('adaptation.csv'))})})

ggplot(xx, aes(iteration, epsvec, color=factor(seed))) + geom_line() +
  facet_wrap('delta') + scale_y_log10()

yy <- subset(xx, iteration==200/2)
ggplot(yy, aes(delta, epsvec)) + geom_point() + scale_y_log10()

zz <- ddply(xx, .(seed, delta), summarize, accept.rate=mean(accepted))

## make sure that when a divergence occurs that the parameters repeat

delta <- .8
seed <- 13
setwd('C:/Users/Cole/admb/examples/admb/catage/')
system('admb_hmc catage')
system(paste('catage -mcmc 2000 -hmc'),
       ignore.stdout=TRUE, ignore.stderr=TRUE )
adapt <- read.csv('adaptation.csv')
## fit <- read_pars('catage')
pars <- read_psv('catage')
pars2 <- array(0, dim=c(nrow(pars), 1, ncol(pars)))
pars2[,1,] <- as.matrix(pars)
str(adapt)
str(pars2)
dimnames(pars2) <- list(iterations=NULL, chains="chain:1", parameters=names(pars))
ss <- monitor(sims=pars2)
y <- vector("list", length=length(names(pars)))
names(y) <- names(pars)
z <- lapply(y, function(x) x=numeric(0))
sso2 <- shinystan:::shinystan(model_name='catage', param_names=names(pars),
                  param_dims=z, posterior_sample=pars2,
                  sampler_params=adapt,
                  summary=ss, n_chain=1, n_iter=nrow(pars),
                  n_warmup=nrow(pars)/2, model_code='NA',
                  misc=list(max_td=10, stan_method='sampling',
                            stan_algorithm='HMC',
                            sso_version=utils::packageVersion('shinystan')))

launch_shinystan(sso2)


library(shinystan)
library(rstan)
hmc.fit <- readRDS('hmc.fit.RDS')
nuts.fit <- readRDS('nuts.fit.RDS')

launch_shinystan(hmc.fit)
launch_shinystan(nuts.fit)
xx <- array(NA, dim=c(100,1,2))
xx[,1,] <- cbind(rnorm(100), rnorm(100))
sso <- as.shinystan(hmc.fit)
ss <- matrix(rnorm(20), ncol=10)
sp <- get_sampler_params(hmc.fit)
str(sp)
sp <- list(accept_stat__=runif(100), stepsize__=runif(100, .5,.9),
           int_time__=runif(100, 6,8), energy=rnorm(100, 50, 5))
pp <- cbind(runif(100), runif(100, .5,.9),
           runif(100, 6,8), rnorm(100, 50, 5))
attr(pp, 'dimnames') <-
  list(NULL, c("accept_stat__", "stepsize__", "int_time__", "energy__"))
sso <- shinystan:::shinystan(model_name='test', param_names=c('x','y'),
                             param_dims=list(x=1,y=1), posterior_sample=xx,
                             sampler_params=list(pp),
                             summary=ss, n_chain=1, n_iter=nrow(xx),
                             n_warmup=50, model_code='NA',
                             misc=list(max_td=10, stan_method='sampling',
                                       stan_algorithm='NUTS',
                                       sso_version=utils::packageVersion('shinystan')))
launch_shinystan(sso)
