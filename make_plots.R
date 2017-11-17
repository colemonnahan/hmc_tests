m <- c('ss_logistic', 'wildflower_nc', 'redkite', 'growth_nc', 'swallows', 'mvnc','mvnd')
## g <- ggplot(subset(empirical, platform!='jags' & model %in% m)) +
##     geom_point(aes(delta.target, log(samples.per.time))) + facet_wrap('model', scales='free_y')
## ggsave('plots/optimal_delta.png', g, width=ggwidth, height=ggheight)

g <- ggplot(simulated, aes(Npar, y=median.efficiency, color=platform)) +
  geom_line() + facet_wrap("model") + scale_y_log10() +
  geom_pointrange(aes(ymin=lwr.efficiency, ymax=upr.efficiency))
ggsave('plots/efficiency_simulated.png', g, width=ggwidth, height=ggheight)

g <- ggplot(simulated, aes(Npar, y=time.total, color=platform)) +
  geom_point() + facet_wrap("model") + scale_y_log10()
ggsave('plots/runtime_simulated.png', g, width=ggwidth, height=ggheight)

g <- ggplot(simulated, aes(Npar, y=log10(minESS), color=platform)) +
  geom_point() + facet_wrap("model") + scale_y_log10()
ggsave('plots/ESS_simulated.png', g, width=ggwidth, height=ggheight)

g <- ggplot(empirical, aes(minESS, y=minESS.coda, group=platform, color=platform))  +
  scale_x_log10()+scale_y_log10()+
    geom_point() + geom_abline(intercept=0,slope=1)+
      facet_wrap('model', scales='fixed')
ggsave('plots/ESS_comparison.png', g, width=ggwidth, height=ggheight)

ggplot(empirical, aes(platform, efficiency, group=seed)) + geom_line() +
  facet_wrap("model") + scale_y_log10()

g <- ggplot(empirical, aes(platform, y=100*(minESS/Nsims)))  +
  geom_point()+ ylim(0,100) + facet_wrap('model') +
    theme(axis.text.x = element_text(angle = 90)) + ylab("% ESS")
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
g <- ggplot(adapt_empirical, aes(model, log(eps.final), color=platform)) + geom_point(alpha=.5)
ggsave('plots/adapt_eps.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, delta.mean, color=platform)) + geom_point(alpha=.5)
ggsave('plots/adapt_delta.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, log10(nsteps.median), color=platform)) + geom_jitter(alpha=.5)
ggsave('plots/adapt_nsteps.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, max_treedepths, color=platform)) + geom_point(alpha=.5)
ggsave('plots/adapt_max_treedepths.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(model, 100*ndivergent/Nsims, color=platform)) +
    geom_point(alpha=.5) + ylab('Percent diverged transitions')
ggsave('plots/adapt_ndivergent.png', g, width=ggwidth, height=ggheight)
g <- ggplot(adapt_empirical, aes(log(eps.final), log10(nsteps.mean),
  color=model, shape=platform)) + geom_point()
ggsave('plots/adapt_eps_vs_nsteps.png', g, width=ggwidth, height=ggheight)

