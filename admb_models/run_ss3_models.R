devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

## Investigate performance differences between algorithms and settings

reps <- 4                        # chains/reps to run
tt <- 100 # thin rate for RWM
td <- 10
iter <- 500
warmup <- (iter/2)
hh <- 10                           # hours to run
d <- 'cod_fast'
m <- 'cod_fast'
ad <- .9                                # adapt_delta
source('template.R')

cov0 <- get.admb.cov(d)

sfStop()
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) mle$est[1:N])
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=3, sigma=covar)))
inits <- NULL
eps <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()



fit.rwm.mle <- sample_admb(m, iter, warmup=warmup, init=inits, chains=reps,
                          thin=tt,  parallel=TRUE, dir=d, algorithm='RWM')
rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
## Examine pairs in bounded space.
posterior <- extract_samples(fit.rwm.mle, inc_lp=TRUE)
ess <- fit.rwm.mle$ess
pars <- names(sort(ess))[1:20]
pairs_admb(posterior=posterior, mle=fit.rwm.mle$mle, pars=pars)
launch_shinyadmb(fit.rwm.mle)

## Now convert to unbounded and check
mle <- fit.rwm.mle$mle
unbounded$lp__ <- posterior$lp__
names(unbounded) <- names(posterior)
mle$cov <- fit.rwm.mle$covar.est
mle$se <- sqrt(diag(mle$cov))
mle$cor <- cov2cor(mle$cov)
mle$est <- apply(unbounded, 2, mean)
pairs_admb(posterior=unbounded, mle=mle, pars=pars, diag='trace')

## Try to get NUTS workign for snowcrab
sfStop()
iter <- 20
warmup <- iter/2
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.unbounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=3, sigma=covar)))
inits <- lapply(1:reps, function(i) mle$est[1:N])
inits <- NULL
eps <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.nuts.unit <- sample_admb(m, iter, warmup=warmup, init=inits, chains=1,
                            thin=1,  parallel=F, dir=d, algorithm='NUTS',
                            control=list(metric=NULL, stepsize=eps))

df <- scan('snowcrab/gradient.dat', what=character(), skip=1)
df2 <- as.data.frame(matrix(df, ncol=3, byrow=TRUE), stringsAsFactors=FALSE)
df2[,2] <- as.numeric(df2[,2])
df2[,3] <- as.numeric(df2[,3])
points(df2[,3], col=2)


launch_shinyadmb(fit.nuts.unit )

rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
## Examine pairs in bounded space.
posterior <- extract_samples(fit.rwm.mle, inc_lp=TRUE)
ess <- fit.rwm.mle$ess
pars <- names(sort(ess))[1:20]
pairs_admb(posterior=posterior, mle=fit.rwm.mle$mle, pars=pars)
launch_shinyadmb(fit.rwm.mle)

## Now convert to unbounded and check
mle <- fit.rwm.mle$mle
unbounded$lp__ <- posterior$lp__
names(unbounded) <- names(posterior)
mle$cov <- fit.rwm.mle$covar.est
mle$se <- sqrt(diag(mle$cov))
mle$cor <- cov2cor(mle$cov)
mle$est <- apply(unbounded, 2, mean)
pairs_admb(posterior=unbounded, mle=mle, pars=pars, diag='trace')


rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
mle <- r4ss::read.admbFit(file.path(d,m))
## Examine pairs in bounded space
pairs_admb(posterior=extract_samples(fit.rwm.mle), mle=mle,
           which.keep=NULL)
## Now convert to unbounded and check
mle$cov <- fit.rwm.mle$covar.est
mle$std <- sqrt(diag(mle$cov))
mle$cor <- cov2cor(mle$cov)
mle$est <- apply(unbounded, 2, mean)
pairs_admb(posterior=unbounded, mle=mle,
           which.keep=NULL)

