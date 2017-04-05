setwd('admb_models/')
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

reps <- 6                        # chains/reps to run
td <- 12
iter <- 2000
warmup <- iter/2
tt <- 100 # thin rate
hh <- 10                           # hours to run
d <- 'catage'
m <- 'catage'
ad <- .9                                # adapt_delta
source('template.R')

fit.nuts.mle <- sample_admb(m, 2000, init=NULL, chains=1, dir=d)
rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
unbounded <- read.csv(file.path(d,'unbounded.csv'), head=FALSE)
mle <- r4ss::read.admbFit(file.path(d,m))
## Examine pairs in bounded space
pairs_admb(posterior=extract_samples(fit.nuts.dense), mle=mle,
           which.keep=NULL)
## Now convert to unbounded and check
mle$cov <- fit.nuts.mle$covar.est
mle$std <- sqrt(diag(mle$cov))
mle$cor <- cov2cor(mle$cov)
mle$est <- apply(unbounded, 2, mean)
pairs_admb(posterior=unbounded, mle=mle,
           which.keep=NULL)

fit.rwm.mle <- sample_admb(m, 2000*100, init=NULL, thin=100, chains=1, dir=d, algorithm='RWM')
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

