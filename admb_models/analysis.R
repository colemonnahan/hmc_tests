## This file is where I do the analysis and explore the MCMC fits from the
## already run models
library(adnuts)
library(shinystan)
library(rstan)
library(plyr)

plot.uncertainties <- function(fit1, fit2=NULL, xlims, ylims){
  n <- NROW(fit1$dq)
  png(paste0('plots/uncertainties_', fit1$model, '.png'), units='in', width=7,
             height=3, res=300)
  par(mfrow=c(1,n), mar=c(5,1,1,1 ), oma=c(0,0, 0, 0))
  for(i in 1:n){
    ii <- fit1$dq$dq[i]
    xx <- fit1$dq.post[,ii]
    hist(xx, freq=FALSE, xlim=xlims[[i]], ylim=ylims[[i]], col=gray(.8),
  border=gray(.8), breaks=50, xlab=fit1$dq$dq[i], main=NA, ylab=NA)
    abline(v=mean(xx), col=2)
    lines(x <- seq(min(xlims[[i]]), max(xlims[[i]]), len=10000), y=dnorm(x, fit1$dq[i, 'mle'], fit1$dq[i,'se']))
    abline(v=fit1$dq[i, 'mle'], col=1)
  }
  dev.off()
}
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
plot.improvement <- function(fit1, fit2){
  ## xx <- rbind(data.frame(model=fit1$model, ess=fit1$ess,
  ##                  perf=fit1$ess/sum(fit1$time.total)),
  ##       data.frame(model=fit2$model, ess=fit2$ess,
  ##                  perf=fit2$ess/sum(fit2$time.total)) )
  ## xx$par <- row.names(xx)
  ## levels(xx$model) <- c("Original", "Fixed")
  ## ggplot(xx, aes(x=model, y=ess)) +geom_violin() + scale_y_log10()
  png(paste0('plots/ess_improvement_',fit1$model, '.png'), width=3, height=5,
      units='in', res=500)
  vioplot::vioplot(log10(fit1$ess), log10(fit2$ess), names=c("Original", "Fixed"))
  mtext(fit1$model, line=1, cex=1.5)
  mtext("log10(ESS)", side=2, line=2.5, cex=1.25)
  dev.off()
}

n.slow <- 5 # number of parameters to show in pairs plot

cod.rwm <- readRDS('results/pilot_rwm_cod.RDS')
cod.post <- extract_samples(cod.rwm, inc_lp=TRUE)
cod.post <- cbind(cod.post, cod.rwm$dq.post)
mle <- cod.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(cod.rwm$samples)[2], each=dim(cod.rwm$samples)[1]-cod.rwm$warmup)
slow <- c(names(sort(cod.rwm$ess))[1:n.slow], names(cod.rwm$dq.post), 'lp__')
png('plots/pairs.cod.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(cod.post, mle=mle, chains=chain, pars=slow, diag='trace')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(cod.rwm),2, var)
var.mle <- diag(cod.rwm$mle$cov)[1:cod.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.', m, '.png'), g, width=7, height=5)
## Compare estimates of DQs
xlims <- list(c(0, 2e9), c(0, 1.5e09))
ylims <- list(c(0, 2.5e-9), c(0, 3.5e-08))
plot.uncertainties(cod.rwm, xlims=xlims, ylims=ylims)


hake.rwm <- readRDS('results/pilot_rwm_hake.RDS')
hake.post <- extract_samples(hake.rwm, inc_lp=TRUE)
hake.post <- cbind(hake.post, hake.rwm$dq.post)
mle <- hake.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(hake.rwm$samples)[2], each=dim(hake.rwm$samples)[1]-hake.rwm$warmup)
slow <- c(names(sort(hake.rwm$ess))[1:n.slow], names(hake.rwm$dq.post), 'lp__')
png('plots/pairs.hake.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(hake.post, mle=mle, chains=chain, pars=slow, diag='trace')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(hake.rwm),2, var)
var.mle <- diag(hake.rwm$mle$cov)[1:hake.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.', m, '.png'), g, width=7, height=5)
## Compare estimates of DQs
xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
plot.uncertainties(hake.rwm, xlims=xlims, ylims=ylims)


n.slow <- 13
canary.rwm <- readRDS('results/pilot_rwm_canary.RDS')
canary.post <- extract_samples(canary.rwm, inc_lp=TRUE)
#canary.post <- cbind(canary.post, canary.rwm$dq.post)
mle <- canary.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(canary.rwm$samples)[2], each=dim(canary.rwm$samples)[1]-canary.rwm$warmup)
slow <- c(names(sort(canary.rwm$ess))[1:n.slow], names(canary.rwm$dq.post), 'lp__')
png('plots/pairs.canary.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(canary.post, mle=mle, chains=chain, pars=slow, diag='acf')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(canary.rwm),2, var)
var.mle <- diag(canary.rwm$mle$cov)[1:canary.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.canary.png'), g, width=7, height=5)
## ## Compare estimates of DQs
## xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
## ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
## plot.uncertainties(canary.rwm, xlims=xlims, ylims=ylims)
canary2.rwm <- readRDS('results/pilot_rwm_canary2.RDS')
canary2.post <- extract_samples(canary2.rwm, inc_lp=TRUE)
#canary2.post <- cbind(canary2.post, canary2.rwm$dq.post)
mle <- canary2.rwm$mle
## Run mcsave and get generated quantities
chain <- rep(1:dim(canary2.rwm$samples)[2], each=dim(canary2.rwm$samples)[1]-canary2.rwm$warmup)
slow <- c(names(sort(canary2.rwm$ess))[1:n.slow], names(canary2.rwm$dq.post), 'lp__')
png('plots/pairs.canary2.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(canary2.post, mle=mle, chains=chain, pars=slow, diag='acf')
dev.off()
## Look at which parameter MLE vs posterior variances are different
var.post <- apply(extract_samples(canary2.rwm),2, var)
var.mle <- diag(canary2.rwm$mle$cov)[1:canary2.rwm$mle$nopar]
vars <- data.frame(post=var.post, mle=var.mle)
g <- ggplot(vars, aes(x=log10(mle), log10(post))) + geom_point(alpha=.7) +
  geom_abline(slope=1) + xlab("MLE Variance") + ylab("Posterior Variance")
ggsave(paste0('plots/vars.canary2.png'), g, width=7, height=5)
plot.improvement(canary.rwm, canary2.rwm)
## ## Compare estimates of DQs
## xlims <- list(c(0, 7e6), c(0, 1.5), c(0, 3e6))
## ylims <- list(c(0, 6e-7), c(0, 3),c(0, 2e-6))
## plot.uncertainties(canary2.rwm, xlims=xlims, ylims=ylims)


halibut.rwm <- readRDS('results/long_rwm_halibut.RDS')
halibut.post <- extract_samples(halibut.rwm, inc_lp=TRUE)
slow <- names(sort(halibut.rwm$ess))[1:n.slow]
png('plots/pairs.halibut.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut.post, mle=halibut.rwm$mle, pars=slow);dev.off()
halibut2.rwm <- readRDS('results/long_rwm_halibut2.RDS')
chain <- rep(1:dim(halibut.rwm$samples)[2], each=dim(halibut.rwm$samples)[1]-halibut.rwm$warmup)
halibut2.post <- extract_samples(halibut2.rwm, inc_lp=TRUE)
png('plots/pairs.halibut2.rwm.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.rwm$mle, pars=slow);dev.off()
recdev2 <- names(halibut2.post)[4:34][1:15]
png('plots/pairs.halibut2.rwm.recdev2.png', width=7, height=5, units='in', res=500)
pairs_admb(halibut2.post, mle=halibut2.rwm$mle, diag='trace',chain=chain, pars=recdev2);dev.off()
plot.improvement(halibut.rwm, halibut2.rwm)
## launch_shinyadmb(halibut2.rwm)
## launch_shinyadmb(halibut2.nuts)


n.slow <- 16
snowcrab.rwm <- readRDS('results/pilot_rwm_snowcrab.RDS')
snowcrab.post <- extract_samples(snowcrab.rwm, inc_lp=TRUE)
chain <- rep(1:dim(snowcrab.rwm$samples)[2], each=dim(snowcrab.rwm$samples)[1]-snowcrab.rwm$warmup)
slow <- names(sort(snowcrab.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.snowcrab.rwm.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=snowcrab.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
## I found this manually by looking at par file
hitbounds <- sort(c(293,309, 310:312, 323, 324,331,6,268, 269, 284, 286, 287))
png('plots/pairs.snowcrab.rwm.hitbounds.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab.post, mle=snowcrab.rwm$mle, chain=chain, diag='trace', pars=hitbounds);dev.off()
snowcrab2.rwm <- readRDS('results/pilot_rwm_snowcrab2.RDS')
snowcrab2.post <- extract_samples(snowcrab2.rwm, inc_lp=TRUE)
chain <- rep(1:dim(snowcrab2.rwm$samples)[2], each=dim(snowcrab2.rwm$samples)[1]-snowcrab2.rwm$warmup)
#slow <- names(sort(snowcrab2.rwm$ess, FALSE))[1:n.slow]
png('plots/pairs.snowcrab2.rwm.slow.png', width=7, height=5, units='in', res=500)
pairs_admb(snowcrab2.post, mle=snowcrab2.rwm$mle, chain=chain, diag='trace', pars=slow);dev.off()
plot.improvement(snowcrab.rwm, snowcrab2.rwm)
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

