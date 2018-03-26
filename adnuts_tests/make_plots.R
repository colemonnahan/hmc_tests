## This file makes some exploratory plots given the empirical and simulated
## results loaded in the workspace

g <- ggplot(simulated, aes(Npar, y=median.efficiency, color=platform)) +
  geom_line() + facet_wrap("model", scales='free') + scale_y_log10() +
  geom_pointrange(aes(ymin=lwr.efficiency, ymax=upr.efficiency))
ggsave('plots/efficiency_simulated.png', g, width=ggwidth, height=ggheight)

g <- ggplot(simulated, aes(Npar, y=time.total, color=platform)) +
  geom_point() + facet_wrap("model", scales='free') + scale_y_log10()
ggsave('plots/runtime_simulated.png', g, width=ggwidth, height=ggheight)

g <- ggplot(simulated, aes(Npar, y=log10(minESS), color=platform)) +
  geom_point() + facet_wrap("model", scales='free') + scale_y_log10()
ggsave('plots/ESS_simulated.png', g, width=ggwidth, height=ggheight)

## ## show iidz vs zdiag to test adaptation of mass matrix
## x2 <- subset(simulated, model %in% c('iidz', 'zdiag'))
## x3 <- dcast(x2, formula=platform+Npar~model, fun=median, value.var='median.efficiency')
## x3$ratio <- x3$zdiag/x3$iidz
## ggplot(x3, aes(log2(Npar), y=ratio, color=platform)) + geom_line()
## ggplot(x2, aes(log2(Npar), y=median.efficiency, group=model, color=model)) + geom_line()+
##   facet_wrap("platform") + scale_y_log10()

g <- ggplot(empirical, aes(minESS, y=minESS.coda, group=platform, color=platform))  +
  scale_x_log10()+scale_y_log10()+
    geom_point() + geom_abline(intercept=0,slope=1)+
      facet_wrap('model', scales='fixed')
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)

g <- ggplot(empirical, aes(platform, efficiency, group=seed)) +
  geom_line( alpha=.5) +
  facet_wrap("model") + scale_y_log10() +
  geom_line(aes(platform, median.efficiency), col='red', lwd=2)
ggsave('plots/empirical_efficiency.png', g, width=ggwidth, height=ggheight)

g <- ggplot(empirical, aes(platform, y=100*(minESS/Nsims)))  +
  geom_point()+  facet_wrap('model') +
  #theme(axis.text.x = element_text(angle = 90)) +
  ylab("% ESS")
ggsave('plots/ESS_percentages.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(platform, y=log(time.total)))  +
  geom_point()+ facet_wrap('model') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/runtime.png', g, width=ggwidth, height=ggheight)
g <- ggplot(empirical, aes(platform, y=Rhat.max))  +
  geom_point()+ facet_wrap('model') +
    theme(axis.text.x = element_text(angle = 90))
ggsave('plots/Rhat.png', g, width=ggwidth, height=ggheight)


## Some adaptation plots
g0 <- ggplot(adapt_empirical, aes(x=platform, color=platform)) +
   facet_wrap("model", scales='free_y') #+ scale_y_log10()

g <- g0+ geom_point(aes(y=eps.final), alpha=.5) + ylab("Step size")
ggsave('plots/adapt_eps.png', g, width=ggwidth, height=ggheight)

g <- g0+ geom_point(aes(y=delta.mean), alpha=.5) + ylab("Acceptance Rate")
ggsave('plots/adapt_delta.png', g, width=ggwidth, height=ggheight)
g <- g0+ geom_jitter(aes(y=nsteps.mean), alpha=.5) + ylab("Mean Steps")
ggsave('plots/adapt_nsteps.png', g, width=ggwidth, height=ggheight)
g <- g0+ geom_point(aes(y=max_treedepths), alpha=.5) + ylab("# Max treedepths")
ggsave('plots/adapt_max_treedepths.png', g, width=ggwidth, height=ggheight)
g <- g0+ geom_point(aes(y=100*ndivergent/Nsims), alpha=.5) + ylab('% Divergences')
ggsave('plots/adapt_ndivergent.png', g, width=ggwidth, height=ggheight)

## g <- ggplot(adapt_empirical, aes(log(eps.final), log10(nsteps.mean),
##   color=model, shape=platform)) + geom_point()
## ggsave('plots/adapt_eps_vs_nsteps.png', g, width=ggwidth, height=ggheight)

