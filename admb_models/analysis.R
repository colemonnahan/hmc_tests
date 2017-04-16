## This file is where I do the analysis and explore the MCMC fits from the
## already run models
library(adnuts)
library(shinystan)
library(rstan)
library(plyr)

plot.ess <- function(rwm, nuts){
  model <- rwm$model
  x <- rwm$ess/sum(rwm$time.total);
  y <- nuts$ess/sum(nuts$time.total)
  temp <- range(c(x, y,0))
  png(paste0('plots/ess_comparison_',model, '.png'), width=7, height=3,
  units='in', res=500)
  par(mfrow=c(1,3))
  plot(x=x, y=y, xlim=temp, ylim=temp, xlab='RWM', ylab='NUTS')
  abline(0,1)
  col1 <- gray(.7)
  barplot(sort(x), ylim=temp, main='RWM', col=col1, border=col1)
  barplot(sort(y), ylim=temp, main='NUTS', col=col1, border=col1)
  dev.off()
}

n.slow <- 8 # number of parameters to show in pairs plot

cod.rwm <- readRDS('results/long_rwm_cod_fast.RDS')
cod.post <- extract_samples(cod.rwm, inc_lp=TRUE)
chain <- rep(1:dim(cod.rwm$samples)[2], each=dim(cod.rwm$samples)[1]-cod.rwm$warmup)
slow <- names(sort(cod.rwm$ess))[1:n.slow]
png('plots/pairs.cod.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(cod.post, mle=cod.rwm$mle, chains=chain, pars=slow);dev.off()
cod.nuts <- readRDS('results/long_nuts_cod_fast.RDS')
cod.post <- extract_samples(cod.nuts, inc_lp=TRUE)
divs <- extract_sampler_params(cod.nuts)$divergent__
slow <- names(sort(cod.nuts$ess))[1:n.slow]
png('plots/pairs.cod.nuts.png', width=7, height=5, units='in', res=500)
pairs_admb(cod.post, mle=cod.nuts$mle, pars=slow, chains=chain, divergences=divs);dev.off()
plot.ess(cod.rwm, cod.nuts)
## launch_shinyadmb(cod.rwm)
## launch_shinyadmb(cod.nuts)

halibut.rwm <- readRDS('results/long_rwm_halibut.RDS')
halibut.post <- extract_samples(halibut.rwm, inc_lp=TRUE)
slow <- names(sort(halibut.rwm$ess))[1:n.slow]
png('plots/pairs.halibut.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut.post, mle=halibut.rwm$mle, pars=slow);dev.off()
halibut.nuts <- readRDS('results/long_nuts_halibut.RDS')
halibut.post <- extract_samples(halibut.nuts, inc_lp=TRUE)
divs <- extract_sampler_params(halibut.nuts)$divergent__
slow <- names(sort(halibut.nuts$ess))[1:n.slow]
png('plots/pairs.halibut.nuts.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut.post, mle=halibut.nuts$mle, pars=slow, divergences=divs);dev.off()
plot.ess(halibut.rwm, halibut.nuts)
## launch_shinyadmb(halibut.rwm)
## launch_shinyadmb(halibut.nuts)
##
halibut2.rwm <- readRDS('results/long_rwm_halibut2.RDS')
halibut2.post <- extract_samples(halibut2.rwm, inc_lp=TRUE)
slow <- names(sort(halibut2.rwm$ess))[1:n.slow]
png('plots/pairs.halibut2.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.rwm$mle, pars=slow);dev.off()
halibut2.nuts <- readRDS('results/long_nuts_halibut2.RDS')
halibut2.post <- extract_samples(halibut2.nuts, inc_lp=TRUE)
divs <- extract_sampler_params(halibut2.nuts)$divergent__
slow <- names(sort(halibut2.nuts$ess))[1:n.slow]
png('plots/pairs.halibut2.nuts.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.nuts$mle, pars=slow, divergences=divs);dev.off()
plot.ess(halibut2.rwm, halibut2.nuts)
## launch_shinyadmb(halibut2.rwm)
## launch_shinyadmb(halibut2.nuts)

hake.rwm <- readRDS('results/long_rwm_hake.RDS')
hake.post <- extract_samples(hake.rwm, inc_lp=TRUE)
slow <- names(sort(hake.rwm$ess))[1:n.slow]
png('plots/pairs.hake.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(hake.post, mle=hake.rwm$mle, pars=slow);dev.off()
hake.nuts <- readRDS('results/long_nuts_hake.RDS')
hake.post <- extract_samples(hake.nuts, inc_lp=TRUE)
divs <- extract_sampler_params(hake.nuts)$divergent__
slow <- names(sort(hake.nuts$ess))[1:n.slow]
png('plots/pairs.hake.nuts.png', width=7, height=5, units='in', res=500)
pairs_admb(hake.post, mle=hake.nuts$mle, pars=slow, divergences=divs);dev.off()
plot.ess(hake.rwm, hake.nuts)
## launch_shinyadmb(hake.rwm)
## launch_shinyadmb(hake.nuts)


snowcrab.rwm <- readRDS('results/long_rwm_2016sc.RDS')
snowcrab.post <- extract_samples(snowcrab.rwm, inc_lp=TRUE)
chain <- rep(1:dim(snowcrab.rwm$samples)[2], each=dim(snowcrab.rwm$samples)[1]-snowcrab.rwm$warmup)
slow <- names(sort(snowcrab.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.snowcrab.rwm.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=snowcrab.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
## I found this manually by looking at par file
hitbounds <- sort(c(293,309, 310:312, 323, 324,331,6,268, 269, 284, 286, 287))
png('plots/pairs.snowcrab.rwm.hitbounds.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=snowcrab.rwm$mle, chain=chain, diag='trace', pars=hitbounds);dev.off()
Rhat <- names(sort(snowcrab.rwm$Rhat, TRUE))[1:n.slow]
png('plots/pairs.snowcrab.rwm.Rhat.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=snowcrab.rwm$mle, chain=chain, diag='trace', pars=Rhat);dev.off()
## launch_shinyadmb(snowcrab.rwm)
## launch_shinyadmb(snowcrab.nuts)

tanner.rwm <- readRDS('results/long_rwm_tanner.RDS')
tanner.post <- extract_samples(tanner.rwm, inc_lp=TRUE)
chain <- rep(1:dim(tanner.rwm$samples)[2], each=dim(tanner.rwm$samples)[1]-tanner.rwm$warmup)
slow <- names(sort(tanner.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.tanner.rwm.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner.post, mle=tanner.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
## I found this manually by looking at par file
hitbounds <- c(1,7,33, 67,70)
png('plots/pairs.tanner.rwm.hitbounds.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner.post, mle=tanner.rwm$mle, chain=chain, diag='trace', pars=hitbounds);dev.off()
Rhat <- names(sort(tanner.rwm$Rhat, TRUE))[1:n.slow]
png('plots/pairs.tanner.rwm.Rhat.png', width=7, height=5, units='in', res=500)
pairs_admb(tanner.post, mle=tanner.rwm$mle, chain=chain, diag='trace', pars=Rhat);dev.off()
plot.ess(tanner.rwm, tanner.nuts)


## Look at which parameter MLE vs posterior variances are different
var.poster <- apply(extract_samples(snowcrab.rwm),2, var)
var.mle <- diag(snowcrab.rwm$mle$cov)[1:snowcrab.rwm$mle$nopar]
plot(log10(var.mle), log10(var.poster)); abline(0,1)

all.fits <- list(cod.rwm, cod.nuts, halibut.rwm, halibut.nuts, hake.rwm,
                 hake.nuts)
perf.wide <- ldply(all.fits, function(x){
  data.frame(model=x$model, alg=x$algorithm, runtime=sum(x$time.total),
        minESS=min(x$ess), maxRhat=max(x$Rhat))})
g <- ggplot(perf.wide, aes(model, maxRhat, color=alg)) + geom_point()
ggsave('plots/maxRhat_comparison.png', g, width=7, height=5)
perf.wide$perf <- with(perf.wide, minESS/runtime)
g <- ggplot(perf.wide, aes(model, y=log(perf), color=alg)) + geom_point()
ggsave('plots/efficienty_comparison.png', g, width=7, height=5)
temp <- reshape2::melt(subset(perf.wide, select=-c(maxRhat)), c('model', 'alg'))
perf.long <- reshape2::dcast(temp, model+variable~alg)
g <- ggplot(perf.long, aes(x=log10(RWM), y=log10(NUTS), color=model)) + geom_point() +
  geom_abline(slope=1)+ facet_wrap('variable', scales='free') + coord_equal()
ggsave('plots/perf_comparison.png', g, width=7, height=5)

