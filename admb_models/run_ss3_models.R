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

rotated <- read.csv(file.path(d,'rotated.csv'), head=FALSE)
apply(rotated, 2, var)
xx <- cor(rotated)
plot(sort((abs(xx[lower.tri(xx, FALSE)]))), type='h')
post <-
mle <- r4ss::read.admbFit(file.path(d,m))
pairs_admb(posterior=extract_samples(fit.nuts.dense), mle=mle,
           which.keep=1:15)
pairs_admb(posterior=extract_samples(fit.rwm.mle), mle=mle,
           which.keep=1:15)
