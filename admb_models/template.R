library(ggplot2)
library(adnuts)
library(snowfall)
library(shinystan)
library(plyr)

## Run model with hbf=1 to get right covariance matrix and MLEs
setwd(d)
system(paste('admb', m))
system(paste(m, '-hbf 1 -nox'))
setwd('..')

sfStop()
## Draw inits from MVT using MLE and covar
inits <- NULL
eps <- NULL
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()


## Run NUTS for different mass matrices
fit.nuts.mle <-
  sample_admb(m, iter=iter, init=inits, algorithm='NUTS',
               parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, metric=NULL,
                           adapt_delta=ad))
## covar.diag <- diag(x=diag(fit.nuts.mle$covar.est))
## fit.nuts.diag <-
##   sample_admb(m, iter=iter, init=inits, algorithm='NUTS',
##               parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
##               control=list(max_treedepth=td, stepsize=eps,
##                            metric=covar.diag, adapt_delta=ad))
covar.dense <- fit.nuts.mle$covar.est
fit.nuts.dense <-
  sample_admb(m, iter=iter, init=inits, algorithm='NUTS',
               parallel=TRUE, chains=reps, warmup=warmup, dir=d, cores=reps,
              control=list(max_treedepth=td, stepsize=eps, extra.args='-noest',
                           metric=covar.dense, adapt_delta=ad))

## Now run RWM but using a thinning rate similar to NUTS so the time is
## roughly equivalent.
tt <- 4*floor(mean(extract_sampler_params( fit.nuts.mle)$n_leapfrog__))
## Rerun model with hbf=0
setwd(d)
system(paste(m, '-nox'))
setwd('..')

fit.rwm.mle <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt,
              parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=NULL),
              algorithm='RWM')
covar.diag <- diag(x=diag(fit.rwm.mle$covar.est))
## fit.rwm.diag <-
##   sample_admb(m, iter=tt*iter, init=inits,  thin=tt,
##               parallel=TRUE, chains=reps, warmup=tt*warmup,
##               dir=d, cores=reps, control=list(metric=covar.diag),
##               algorithm='RWM')
covar.dense <- fit.rwm.mle$covar.est
fit.rwm.dense <-
  sample_admb(m, iter=tt*iter, init=inits, thin=tt,
              parallel=TRUE, chains=reps, warmup=tt*warmup,
              dir=d, cores=reps, control=list(metric=covar.dense),
              algorithm='RWM')

## Gather adaptation and performance metrics
ff <- function(labels, ...){
  fits <- list(...)
  if(length(fits)!=length(labels)) stop("aa")
  do.call(rbind,  lapply(1:length(fits), function(i){
    x <- fits[[i]]$sampler_params
    data.frame(m=labels[i], chain=1:length(x),
               eps=as.numeric(do.call(rbind, lapply(x, function(l) tail(l[,2],1)))),
               divergences=as.numeric(do.call(rbind, lapply(x, function(l) sum(l[-(1:fits[[i]]$warmup),5])))),
               accept_prob=as.numeric(do.call(rbind, lapply(x,
                          function(l) mean(l[-(1:fits[[i]]$warmup),1])))),
               nsteps=as.numeric(do.call(rbind, lapply(x,
                          function(l) mean(l[-(1:fits[[i]]$warmup),4])))))
  }))}
adaptation <- ff(c("mle", "dense"), fit.nuts.mle, fit.nuts.dense)
stats.nuts.mle <- with(fit.nuts.mle, data.frame(alg='nuts', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.mle <- data.frame(alg='nuts', m='mle', efficiency=min(stats.nuts.mle$n_eff)/sum(fit.nuts.mle$time.total))
## stats.nuts.diag <- with(fit.nuts.diag, data.frame(alg='nuts', m='diag', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
## perf.nuts.diag <- data.frame(alg='nuts', m='diag', efficiency=min(stats.nuts.diag$n_eff)/sum(fit.nuts.diag$time.total))
stats.nuts.dense <- with(fit.nuts.dense, data.frame(alg='nuts', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.nuts.dense <- data.frame(alg='nuts', m='dense', efficiency=min(stats.nuts.dense$n_eff)/sum(fit.nuts.dense$time.total))
stats.rwm.mle <- with(fit.rwm.mle, data.frame(alg='rwm', m='mle', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.mle <- data.frame(alg='rwm', m='mle', efficiency=min(stats.rwm.mle$n_eff)/sum(fit.rwm.mle$time.total))
## stats.rwm.diag <- with(fit.rwm.diag, data.frame(alg='rwm', m='diag', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
## perf.rwm.diag <- data.frame(alg='rwm', m='diag', efficiency=min(stats.rwm.diag$n_eff)/sum(fit.rwm.diag$time.total))
stats.rwm.dense <- with(fit.rwm.dense, data.frame(alg='rwm', m='dense', time.total=sum(time.total), rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE)))
perf.rwm.dense <- data.frame(alg='rwm', m='dense', efficiency=min(stats.rwm.dense$n_eff)/sum(fit.rwm.dense$time.total))
stats.all <- rbind(stats.nuts.mle, stats.nuts.dense, stats.rwm.mle, stats.rwm.dense)
stats.all[,c('mean', 'se_mean', 'sd', 'X50.')] <- NULL
stats.all <- ddply(stats.all, .(alg, m), mutate, perf=(n_eff)/time.total)
stats.long <- reshape2::melt(stats.all, c('alg', 'm'))
perf.all <- rbind(perf.nuts.mle, perf.nuts.dense, perf.rwm.mle, perf.rwm.dense)
adaptation.long <- reshape2::melt(adaptation, c('m', 'chain'))

## Quick plots
ggwidth <- 7
ggheight <- 5
g <- ggplot(adaptation.long, aes(y=value, x=m)) + geom_point(alpha=.5) +
  facet_wrap('variable', scales='free')
ggsave(paste0('plots/', d, '_adaptation.png'),g, width=ggwidth, height=ggheight, units='in')
g <- ggplot(stats.long, aes(y=value, x=m, color=alg)) + geom_jitter(alpha=.5) +
  facet_wrap('variable', scales='free')
ggsave(paste0('plots/', d, '_stats.png'),g, width=ggwidth, height=ggheight, units='in')
g <- ggplot(perf.all, aes(m, efficiency, color=alg)) + geom_point() + scale_y_log10()
ggsave(paste0('plots/', d, '_perf.png'),g, width=ggwidth, height=ggheight, units='in')

## Save fits
saveRDS(list(fit.nuts.mle=fit.nuts.mle, #fit.nuts.diag=fit.nuts.diag,
             fit.nuts.dense=fit.nuts.dense, fit.rwm.mle=fit.rwm.mle,
             #fit.rwm.diag=fit.rwm.diag,
             fit.rwm.dense=fit.rwm.dense),
        file=paste0('results/', d,'_fits.RDS'))
