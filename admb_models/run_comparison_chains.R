## devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(rstan)
library(adnuts)
library(snowfall)
td <- 10
## Investigate performance differences between algorithms and settings

reps <- 3                        # chains/reps to run

iter <- 1000; warmup <- (iter/10)
m <- d <- 'halibut2'
ad <- .9                                # adapt_delta
source('template.R')

reps <- 6                        # chains/reps to run
iter <- 1000; warmup <- (iter/10)
m <- d <- 'cod'
ad <- .9                                # adapt_delta
source('template.R')


## Old code to look at unbounded space
rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
## Examine pairs in bounded space.
posterior <- extract_samples(fit.rwm.mle, inc_lp=TRUE)
ess <- fit.rwm.mle$ess
pars <- names(sort(ess))[1:8]
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


