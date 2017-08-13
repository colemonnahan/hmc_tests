## Generate simulated data

adapt.list <- perf.list <- list()
k <- 1
for(i in seq_along(Npar.vec)){
  Npar <- Npar.vec[i]
  message(paste("======== Starting Npar=", Npar))

  ## Update data, inits and model objects
  data <- list(n=Npar, x=rep(0, Npar))
  inits <- lapply(1:length(seeds), function(i) list(mu=rnorm(Npar)))
  obj.stan <- stan(fit=obj.stan, data=data, iter=100,
                   chains=1, init=list(inits[[1]]), verbose=TRUE)
  obj.tmb <- MakeADFun(data=data, parameters=inits[[1]], DLL=m)
  setwd('admb')
  write.table(x=unlist(data), file=paste0(m,'.dat'), row.names=FALSE,
                                          col.names=FALSE )
  system(paste(m, '-hbf')); setwd('..')
  ## Run with updated model size
  temp <- run.chains(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m,
                     data=data, inits=inits, pars=pars, metric='diag', seeds=seeds,
                     Nout=Nout, Nthin=1, delta=delta, useRWM=TRUE)
  adapt.list[[k]] <- cbind(cor=j,temp$adapt)
  perf.list[[k]] <- cbind(cor=j,temp$perf)
  ## Save them as we go in case it fails
  perf <-  do.call(rbind, perf.list)
  adapt <- do.call(rbind, adapt.list)
  plot.simulated.results(perf, adapt)
  write.csv(x=perf, file=paste0(m,'_perf_simulated.csv'))
  write.csv(x=adapt, file=paste0(m,'_adapt_simulated.csv'))
  rm(temp)
  k <- k+1
}



