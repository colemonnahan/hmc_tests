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

n.slow <- 10 # number of parameters to show in pairs plot

cod.rwm <- readRDS('results/long_rwm_cod_fast.RDS')
cod.post <- extract_samples(cod.rwm)
slow <- names(sort(cod.rwm$ess))[1:n.slow]
png('plots/pairs.cod.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(cod.post, mle=cod.rwm$mle, pars=slow);dev.off()
cod.nuts <- readRDS('results/long_nuts_cod_fast.RDS')
cod.post <- extract_samples(cod.nuts, inc_lp=TRUE)
divs <- extract_sampler_params(cod.nuts)$divergent__
slow <- names(sort(cod.nuts$ess))[1:n.slow]
png('plots/pairs.cod.nuts.png', width=7, height=5, units='in', res=500)
pairs_admb(cod.post, mle=cod.nuts$mle, pars=slow, divergences=divs);dev.off()
plot.ess(cod.rwm, cod.nuts)
## launch_shinyadmb(cod.rwm)
## launch_shinyadmb(cod.nuts)

halibut.rwm <- readRDS('results/long_rwm_halibut.RDS')
halibut.post <- extract_samples(halibut.rwm)
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

hake.rwm <- readRDS('results/long_rwm_hake.RDS')
hake.post <- extract_samples(hake.rwm)
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

