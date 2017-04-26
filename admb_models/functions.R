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
  library(vioplot)
  png(paste0('plots/ess_improvement_',fit1$model, '.png'), width=3, height=5,
      units='in', res=500)
  vioplot(log10(fit1$ess), log10(fit2$ess), names=c("Original", "Fixed"))
  mtext(fit1$model, line=1, cex=1.5)
  mtext("log10(ESS)", side=2, line=2.5, cex=1.25)
  dev.off()
}
