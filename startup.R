## Load libraries, functions, and global variables
## devtools::install('C:/Users/Cole/adnuts')
## devtools::document('C:/Users/Cole/adnuts')
library(adnuts)
library(coda)
library(ggplot2)
library(plyr)
library(rstan)
library(reshape2)
library(TMB)
library(R2admb)
library(shinystan)
ggwidth <- 8
ggheight <- 5


#' Run Stan, TMB, and ADMB versions of the same model with NUTS.
#'
#' @param model Character string for model name
#' @param seeds A vector of seeds to run across (sequentially) for each
#' algorithm.
#' @param inits A list of lists of randomly drawn independent starting
#' values, of the same length as seeds
#' @param Nout The number of resulting iterations after thinning and warmup
#' @param Nthin The number to thin, defaults to 1.
#' @param delta A vector of target acceptance rates. Defaults to 0.8.
#' @param sink.console Whether to sink console output to trash file to cleanup
#' console output. Defaults to TRUE. Makes it easier to see the progress on
#' the console.
#' @param max_treedepth The maximum times the NUTS algorithm can double
#' before exiting that iteration. Default is 10 in Stan and this function.
#' @param metric A vector of metrics to run across for HMC and NUTS. Must
#' be one of c("unit_e", "diag_e", "dense_e").
#' @return A list of two data frames. adapt is the adaptive results from
#' Stan, and perf is the performance metrics for each run.
run.chains <- function(obj.stan, obj.tmb, model, covar, seeds, Nout, Nthin=1, delta=.8,
                       metric='diag', data, inits, pars,
                       sink.console=FALSE, max_treedepth=12){
  if(Nthin!=1) stop('this probably breaks if Nthin!=1')
  Niter <- 2*Nout*Nthin
  Nwarmup <- Niter/2
  ind.warmup <- 1:Nwarmup              # index of samples, excluding warmup
  ind.samples <- (Nwarmup+1):Niter     # index of warmup samples
  perf.list <- list()
  adapt.list <- list()
  k <- 1
  ## if(sink.console){
  ##   sink(file='trash.txt', append=FALSE, type='output')
  ##   on.exit(sink())
  ## }

  for(seed in seeds){
    message(paste('==== Starting seed',seed, 'at', Sys.time()))
    inits.seed <- list(inits[[which(seed==seeds)]])
    set.seed(seed)
    for(idelta in delta){
      for(imetric in metric){
        ## Crazy way to convert imetric into arguments for the different
        ## platforms.
        N <- length(as.vector(unlist(inits.seed)))
        if(imetric=='unit')
          M <- diag(N)
        else if(imetric=='diag'){
          M <- matrix(0, N, N)
          diag(M) <- diag(covar)
          }
        else if(imetric=='dense')
          M <- covar
        else stop("Invalid metric option supplied")

        fit.stan <-
          stan(fit=obj.stan, data=data, iter=Niter,
               warmup=Nwarmup, chains=1, thin=Nthin, algorithm='NUTS',
               init=inits.seed, seed=seed, par=pars,
               control=list(adapt_engaged=TRUE, adapt_delta=idelta,
                            metric=paste0(imetric,'_e'), max_treedepth=max_treedepth))
        sims.stan <- extract(fit.stan, permuted=FALSE)
        perf.stan <- data.frame(monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5))
        Rhat.stan <- with(perf.stan, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
        adapt <- as.data.frame(get_sampler_params(fit.stan, inc_warmup=FALSE))
        adapt.list[[k]] <-
          data.frame(platform='stan', seed=seed,
                     Npar=dim(sims.stan)[3]-1,
                     Nsims=dim(sims.stan)[1],
                     delta.mean=mean(adapt$accept_stat__),
                     delta.target=idelta,
                     eps.final=tail(adapt$stepsize__,1),
                     max_treedepths=sum(adapt$treedepth__>max_treedepth),
                     ndivergent=sum(adapt$divergent__),
                     nsteps.mean=mean(adapt$n_leapfrog__),
                     nsteps.median=median(adapt$n_leapfrog__),
                     nsteps.sd=sd(adapt$n_leapfrog__),
                     metric=imetric)
        perf.list[[k]] <-
          data.frame(platform='stan',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=tail(adapt$stepsize__,1),
                     Npar=dim(sims.stan)[3]-1,
                     time.warmup= get_elapsed_time(fit.stan)[1],
                     time.total= sum(get_elapsed_time(fit.stan)),
                     minESS=min(perf.stan$n_eff),
                     medianESS=median(perf.stan$n_eff),
                     Nsims=dim(sims.stan)[1],
                     minESS.coda=min(effectiveSize(as.data.frame(sims.stan[,1,]))),
                     Rhat.stan)
        k <- k+1
        ## Start of TMB run
        fit.tmb <-
          sample_tmb(obj=obj.tmb, iter=Niter, warmup=Nwarmup, chains=1, thin=Nthin,
                     init=inits.seed, control=list(metric=M,
          adapt_delta=idelta, max_treedepth=max_treedepth))
        ## saveRDS(fit.tmb, file=paste('fits/tmb_', metric, idelta, seed,'.RDS', sep='_'))
        sims.tmb <- fit.tmb$samples[-(1:Nwarmup),,, drop=FALSE]
        perf.tmb <- data.frame(monitor(sims=sims.tmb, warmup=0, print=FALSE, probs=.5))
        Rhat.tmb <- with(perf.tmb, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
        adapt <- as.data.frame(fit.tmb$sampler_params[[1]])
        adapt.list[[k]] <-
          data.frame(platform='tmb', seed=seed,
                     Npar=dim(sims.tmb)[3]-1,
                     Nsims=dim(sims.tmb)[1],
                     delta.mean=mean(adapt$accept_stat__),
                     delta.target=idelta,
                     eps.final=tail(adapt$stepsize__,1),
                     max_treedepths=sum(adapt$treedepth__>max_treedepth),
                     ndivergent=sum(adapt$divergent__),
                     nsteps.mean=mean(adapt$n_leapfrog__),
                     nsteps.median=median(adapt$n_leapfrog__),
                     nsteps.sd=sd(adapt$n_leapfrog__),
                     metric=imetric)
        perf.list[[k]] <-
          data.frame(platform='tmb',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=tail(adapt$stepsize__,1),
                     Npar=dim(sims.tmb)[3]-1,
                     time.warmup=fit.tmb$time.warmup,
                     time.total=fit.tmb$time.total,
                     minESS=min(perf.tmb$n_eff),
                     medianESS=median(perf.tmb$n_eff),
                     Nsims=dim(sims.tmb)[1],
                     minESS.coda=min(effectiveSize(as.data.frame(sims.tmb[,1,]))),
                     Rhat.tmb)
        k <- k+1
        fit.admb <- sample_admb(dir='admb', model=model, iter=Niter, warmup=Nwarmup,
                                   init=inits.seed, thin=Nthin,
          control=list(metric=M, max_treedepth=max_treedepth, adapt_delta=idelta))
        sims.admb <- fit.admb$samples[-(1:Nwarmup),,, drop=FALSE]
        perf.admb <- data.frame(monitor(sims=sims.admb, warmup=0, print=FALSE, probs=.5))
        Rhat.admb <- with(perf.admb, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
        adapt <- as.data.frame(fit.admb$sampler_params[[1]])
        adapt.list[[k]] <-
          data.frame(platform='admb', seed=seed,
                     Npar=dim(sims.admb)[3]-1,
                     Nsims=dim(sims.admb)[1],
                     delta.mean=mean(adapt$accept_stat__),
                     delta.target=idelta,
                     eps.final=tail(adapt$stepsize__,1),
                     max_treedepths=sum(adapt$treedepth__>max_treedepth),
                     ndivergent=sum(adapt$divergent__),
                     nsteps.mean=mean(adapt$n_leapfrog__),
                     nsteps.median=median(adapt$n_leapfrog__),
                     nsteps.sd=sd(adapt$n_leapfrog__),
                     metric=imetric)
        perf.list[[k]] <-
          data.frame(platform='admb',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=tail(adapt$stepsize__,1),
                     Npar=dim(sims.admb)[3]-1,
                     time.warmup=fit.admb$time.warmup,
                     time.total=fit.admb$time.total,
                     minESS=min(perf.admb$n_eff),
                     medianESS=median(perf.admb$n_eff),
                     Nsims=dim(sims.admb)[1],
                     minESS.coda=min(effectiveSize(as.data.frame(sims.admb[,1,]))),
                     Rhat.admb)
        k <- k+1
        ## ADMB default RWM algorithm
        fit.admb.rwm <-
          sample_admb(dir='admb', model=model, iter=10*Niter, warmup=Nwarmup,
                      init=inits.seed, chains=1, thin=10*Nthin,
                      control=list(metric=M,algorithm='RWM'))
        sims.admb.rwm <- fit.admb.rwm$samples[-(1:Nwarmup),,, drop=FALSE]
        perf.admb.rwm <- data.frame(monitor(sims=sims.admb.rwm, warmup=0, print=FALSE, probs=.5))
        Rhat.admb.rwm <- with(perf.admb.rwm, data.frame(Rhat.min=min(Rhat), Rhat.max=max(Rhat), Rhat.median=median(Rhat)))
        adapt <- as.data.frame(fit.admb.rwm$sampler_params[[1]])
        perf.list[[k]] <-
          data.frame(platform='admb.rwm',
                     seed=seed, delta.target=idelta, metric=imetric,
                     eps.final=NA,
                     Npar=dim(sims.admb.rwm)[3]-1,
                     time.warmup=fit.admb.rwm$time.warmup,
                     time.total=fit.admb.rwm$time.total,
                     minESS=min(perf.admb.rwm$n_eff),
                     medianESS=median(perf.admb.rwm$n_eff),
                     Nsims=dim(sims.admb.rwm)[1],
                     minESS.coda=min(effectiveSize(as.data.frame(sims.admb.rwm[,1,]))),
                     Rhat.admb.rwm)
        k <- k+1
        ## End of ADMB model
      } # End loop over metric
    }   # end loop over adapt_delta
    rm(fit.stan, sims.stan, perf.stan, adapt)
    rm(fit.tmb, sims.tmb, perf.tmb)
    rm(fit.admb, sims.admb, perf.admb)
    rm(fit.admb.rwm, sims.admb.rwm, perf.admb.rwm)
  }
  perf <- do.call(rbind.fill, perf.list)
  perf$efficiency <- perf$minESS/perf$time.total
  adapt <- do.call(rbind.fill, adapt.list[!ldply(adapt.list, is.null)])
  perf$model <- adapt$model <- model
  return(invisible(list(adapt=adapt, perf=perf)))
}

#' Verify models and then run empirical tests across delta
fit.empirical <- function(obj.stan, obj.tmb, model, pars, covar, inits, data, seeds,
                          delta, model.stan, Nout,  metric,
                          Nthin=1, sink.console=FALSE, ...){
    ## Now rerun across gradient of acceptance rates and compare to JAGS
  message('Starting empirical runs')
  ## For now using diagonal to more closely match what Stan is doing.
  covar2 <- diag(x=diag(covar))
    results.empirical <-
      run.chains(obj.stan=obj.stan, obj.tmb=obj.tmb, model=model, seeds=seeds,
                 Nout=Nout, covar=covar2,
                 metric=metric, delta=delta, data=data,
                 Nthin=Nthin, inits=inits, pars=pars,
                 sink.console=sink.console, ...)
    with(results.empirical, plot.empirical.results(perf, adapt))
    write.csv(file=paste0(m, '_adapt_empirical.csv'), results.empirical$adapt)
    write.csv(file=paste0(m, '_perf_empirical.csv'), results.empirical$perf)
}

#' Make plots comparing the performance of simulated data for a model.
#'
#' @details Makes plots with Npar (model complexity) on the x-axis and a
#' variety of things on the y-axis.
#' @param perf A data frame containing model performance data, as returned
#' by run.chains
#' @param adapt A data frame containing adaptation information from Stan,
#' as returned by run.chains.
#'
#' @return Nothing. Makes plots in local folder of performance comparisons
#' and adaptation results
plot.simulated.results <- function(perf, adapt){
    model.name <- as.character(perf$model[1])
    perf.long <- melt(perf,
                      id.vars=c('model', 'platform', 'seed', 'Npar',
                                'Nsims', 'metric', 'cor'),
                      measure.vars=c('time.total', 'minESS', 'efficiency'))
    perf.long <- ddply(perf.long, .(platform, cor, Npar, variable), mutate,
                       mean.value=mean(value))
    g <- ggplot(perf.long, aes(Npar, log10(value), group=platform, color=platform)) +
        geom_point()+
            geom_line(data=perf.long, aes(Npar, log10(mean.value))) +
                facet_grid(variable~cor, scales='free_y') + ggtitle("Performance Comparison")
    ggsave(paste0('plots/', model.name, '_perf_simulated.png'), g, width=ggwidth, height=ggheight)
    adapt$pct.divergent <- with(adapt, ndivergent/Nsims)
    adapt$pct.max.treedepths <- with(adapt, max_treedepths/Nsims)
    adapt.long <- melt(adapt,
                      id.vars=c('model', 'platform', 'seed', 'Npar',
                                'Nsims', 'cor'),
                      measure.vars=c('delta.mean', 'eps.final',
                                     'pct.divergent', 'nsteps.mean'))
    adapt.long <- ddply(adapt.long, .(platform, cor, Npar, variable), mutate,
                       mean.value=mean(value))
    g <- ggplot(adapt.long, aes(Npar, value, group=platform, color=platform)) +
        geom_point() + geom_line(data=adapt.long, aes(Npar, mean.value)) +
                facet_grid(variable~cor, scales='free_y') + ggtitle("Adaptation Comparison")
    ggsave(paste0('plots/',model.name, '_adapt_simulated.png'), g, width=ggwidth, height=ggheight)
}

#' Make plots comparing the performance of empirical data for a model.
#'
#' @details Makes plots with delta (acceptance rate) on the x-axis and a
#' variety of things on the y-axis.
#' @param perf A data frame containing model performance data, as returned
#' by run.chains
#' @param adapt A data frame containing adaptation information from Stan,
#' as returned by run.chains.
#'
#' @return Nothing. Makes plots in local folder of performance comparisons
#' and adaptation results
plot.empirical.results <- function(perf, adapt){
  model.name <- as.character(perf$model[1])
  perf.long <-
    melt(perf, id.vars=c('platform', 'seed', 'delta.target', 'metric'),
         measure.vars=c('time.total', 'minESS', 'efficiency'))
  perf.long$seed2 <- as.factor(with(perf.long, paste(seed, metric, sep="_")))
  g <- ggplot(perf.long, aes(platform, log10(value), group=seed2, color=metric))+
    geom_line() + geom_point() + facet_grid(variable~., scales='free_y')
  ggsave(paste0('plots/',model.name, '_perf_empirical.png'), g, width=ggwidth, height=ggheight)
  adapt.long <- melt(adapt, id.vars=c('platform', 'seed', 'delta.target', 'metric'),
                     measure.vars=c('eps.final', 'delta.mean', 'nsteps.mean'))
  adapt.long$seed2 <- as.factor(with(adapt.long, paste(seed, metric, sep="_")))
  g <- ggplot(adapt.long, aes(platform, value, group=seed2, color=metric))+
    geom_line() + facet_grid(variable~., scales='free_y')
  ggsave(paste0('plots/',model.name, '_adapt_empirical.png'), g, width=ggwidth, height=ggheight)
}


#' Verify the models are the same (coding errors) by plotting QQ plots of
#' each parameter between Stan and JAGS.
#'
#' @param sims.stan A data frame of MCMC samples from a single chain of
#' a Stan run, using extract(fit.stan).
#' @param sims.tmb A data frame of MCMC samples from a single chain of a
#' TMB run.
#' @param perf.platforms A wide data frame with platform, variable and
#' values for Rhat and n_eff, one for each variable of a chain for each
#' platform. As created by verify.model.
#' @return ggplot object invisibly. Also makes plot in folder 'plots' in
#' current working directory.
plot.model.comparisons <- function(sims.stan, sims.tmb, sims.admb, perf.platforms=NULL){
  ## Clean up names so they match exactly
  names(sims.stan) <- gsub('\\.', '', x=names(sims.stan))
  names(sims.tmb) <- gsub('\\.', '', x=names(sims.tmb))
  names(sims.admb) <- gsub('\\.', '', x=names(sims.admb))
  sims.stan$lp__ <- sims.tmb$lp__ <- sims.admb$lp__ <- NULL
  par.names <- names(sims.tmb)
  sims.stan <- sims.stan[,par.names]
  ## Massage qqplot results into long format for ggplot
  qq <- ldply(1:length(par.names), function(i){
    tmb <- as.data.frame(qqplot(sims.tmb[,i], sims.stan[,i], plot.it=FALSE))
    admb <- as.data.frame(qqplot(sims.admb[,i], sims.stan[,i], plot.it=FALSE))
    return(rbind(cbind(par=par.names[i], platform='tmb',tmb),
            cbind(par=par.names[i], platform='admb',admb)))})
  ## Since can be too many parameters, break them up into pages. Stolen
  ## from
  ## http://stackoverflow.com/questions/22996911/segment-facet-wrap-into-multi-page-pdf
  noVars <- length(par.names)
  noPlots <- 25
  plotSequence <- c(seq(0, noVars-1, by = noPlots), noVars)
  ## pdf('plots/model_comparison_qqplots.pdf', onefile=TRUE,
  ## width=ggwidth,height=ggheight)
  png('plots/model_comparison_qqplots%02d.png', units='in', res=500,
      width=ggwidth,height=ggheight)
  for(ii in 2:length(plotSequence)){
    start <- plotSequence[ii-1] + 1;   end <- plotSequence[ii]
    tmp <- subset(qq, par %in% par.names[start:end])
    g <- ggplot(tmp, aes(x,y, color=platform))+ geom_point(alpha=.5) +
      geom_abline(slope=1, col='red') + facet_wrap('par',
                                                   scales='free', nrow=5) + xlab('tmb')+ ylab('stan')
    ## theme(axis.text.x=element_blank(), axis.text.y=element_blank())
    g <- g+ theme(text=element_text(size=7))
    print(g)
  }
  dev.off()
  if(!is.null(perf.platforms)){
    g <- ggplot(perf.platforms, aes(platform, value)) +
      facet_wrap('variable', scales='free')
    g <- g + if(nrow(perf.platforms)>50) geom_violin() else geom_point()
    ggsave('plots/model_comparison_convergence.png', g, width=9, height=5)
  }
  return(NULL)
}
#' Run models with thinning to get independent samples which are then used
#' to verify the posteriors are the same, effectively checking for bugs
#' between models before doing performance comparisons
#'
#' @param admb.columns Columns of the ADMB sample outputs to
#'   exponentiate. THis is needed b/c ADMB needs to add jacobian manually
#'   for bounded (0, Inf) parameters. Thus exponentiate these columns to
#'   match TMB and Stan.
verify.models <- function(obj.stan, obj.tmb, model, covar, pars, inits, data, Nout, Nthin,
                          sink.console=TRUE, dir=NULL, admb.columns=NULL, ...){
  message('Starting independent runs')
  ## if(sink.console){
  ##   sink(file='trash.txt', append=FALSE, type='output')
  ##   on.exit(sink())
  ## }
  Niter <- 2*Nout*Nthin
  Nwarmup <- Niter/2
  fit.stan <- stan(fit=obj.stan, data=data, iter=Niter, chains=1,
                   thin=Nthin, init=inits)
  sims.stan <- extract(fit.stan, permuted=FALSE)
  perf.stan <- data.frame(rstan::monitor(sims=sims.stan, warmup=0, print=FALSE, probs=.5))
  fit.tmb <- sample_tmb(obj=obj.tmb, iter=Niter,
               warmup=Nwarmup, chains=1, thin=Nthin,
               init=inits, control=list(metric=covar), ...)
  sims.tmb <- fit.tmb$samples[-(1:fit.tmb$warmup),,,drop=FALSE]
  perf.tmb <- data.frame(rstan::monitor(sims=sims.tmb, warmup=0, print=FALSE, probs=.5))
  fit.admb <- sample_admb(dir=dir, model=model, iter=Niter,
               warmup=Nwarmup, chains=1, thin=Nthin,
               init=inits, control=list(metric=covar))
  sims.admb <- fit.admb$samples[-(1:fit.admb$warmup),,,drop=FALSE]
  if(!is.null(admb.columns))
    sims.admb[,,admb.columns] <- exp(sims.admb[,,admb.columns])
  perf.admb <- data.frame(rstan::monitor(sims=sims.admb, warmup=0, print=FALSE, probs=.5))
  perf.platforms <- rbind(cbind(platform='tmb',perf.tmb),
                          cbind(platform='admb',perf.admb),
                          cbind(platform='stan',perf.stan))
  perf.platforms <- melt(perf.platforms, c('Rhat', 'n_eff'), id.vars='platform')
  plot.model.comparisons(
    sims.stan=as.data.frame(sims.stan[,1,]),
    sims.tmb=as.data.frame(sims.tmb[,1,]),
    sims.admb=as.data.frame(sims.admb[,1,]),
                         perf.platforms)
  ## Save independent samples for use in intial values later.
  sims.ind <- data.frame(sims.stan[,1,])
  sims.ind <- sims.ind[, names(sims.ind) != 'lp__']
  saveRDS(sims.ind, file='sims.ind.RDS')
}
make.trace <- function(df.thinned, model, string){
    nrows <- ceiling(sqrt(ncol(df.thinned)))
    png(paste0('plots/', model, '.trace.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df.thinned)){
        plot(df.thinned[,i], type='l', axes=FALSE,
             ylim=range(df.thinned[,i]), col=rgb(0,0,0,.5)); box()
        title(names(df.thinned)[i], line=-1)
    }
    dev.off()
}
make.acf <- function(df, model, string){
    nrows <- ceiling(sqrt(ncol(df)))
    png(paste0('plots/', model, '.acf.', string,'.png'), width=9, height=6,
        units='in', res=500)
    par(mfrow=c(nrows,nrows), mar=.1*c(1,1,1,1))
    for(i in 1:ncol(df)) {
        acf(df[,i], axes=FALSE);box()
        title(names(df)[i], line=-1)
    }
    dev.off()
}

## Growth model functions
sample.vbgf <- function(ages, Linf, k,  t0, sigma.obs){
    lengths <- Linf*(1-exp(-k*(ages-t0)))
    loglengths <- log(lengths)+ rnorm(n=length(lengths), mean=0, sd=sigma.obs)
    data.frame(ages=ages, loglengths=loglengths)
}
sample.ages <- function(n.ages, t0, Ntime) {sample((t0+1):Ntime, size=n.ages, replace=FALSE)}
sample.lengths <- function(Nfish, n.ages, logLinf.mean, logLinf.sigma,
                           logk.mean, logk.sigma, sigma.obs, t0){
    Linf.vec <- exp(logLinf.mean + rnorm(n=Nfish, 0, sd=logLinf.sigma))
    k.vec <- exp(logk.mean +rnorm(n=Nfish, mean=0, sd=logk.sigma))
    dat <- ldply(1:Nfish, function(i)
        cbind(fish=i, sample.vbgf(ages=sample.ages(n.ages, t0=t0, Ntime=Ntime),
              Linf=Linf.vec[i], k=k.vec[i], sigma.obs=sigma.obs, t0=t0)))
   return( dat)
}
ss_logistic.traj <- function(r, K, num.years, sd.catch, prop.caught,
                             years.obs, sd.obs, sd.process, F, plot=TRUE){
    catches <- trajectory <- rep(0,num.years)
    u <- rnorm(n=num.years, mean=0, sd=sd.process)
    trajectory[1] <- K + u[1]
    for( yr in 2: num.years){
        Nt <- trajectory[yr-1]
        Ct <- catches[yr-1] <- Nt*(1-exp(-ifelse(Nt<K/10, 0, F)))
        trajectory[yr] <- (Nt+r*Nt*(1- (Nt/K))-Ct) * exp(u[yr]-sd.process^2/2)
        if( trajectory[yr]<=0 ) { trajectory[yr] <- NA; break}
    }
    ## generate lognormal samples
    log.pop.obs <-
        rnorm(n=length(years.obs),
              mean=log(trajectory[years.obs]),sd=sd.obs) -sd.obs^2/2
    if(plot){
        plot(1:num.years, ylim=c(0,K*1.1),y= trajectory, type='l')
        points(1:num.years, catches, type='h')
        points(years.obs, exp(log.pop.obs), col='red')
    }
    return(list(N=num.years, catches=catches, u=u, logcpue=log.pop.obs))
    ## return(trajectory)
}
clean.TMB.files <- function(model.path=getwd()){
  o <- paste0(model.path,'.o')
  dll <- paste0(model.path, '.dll')
  tryCatch(expr=dyn.unload(dynlib(model.path)), error=function(cond){x=3})
  if(file.exists(dll)) trash <- file.remove(dll)
  if(file.exists(o)) trash <- file.remove(o)
}
