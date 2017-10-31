## Generate simulated data
m <- 'growth'
adapt.list <- perf.list <- list()
k <- 1
for(i in seq_along(Npar.vec)){
  Npar <- Npar.vec[i]
  message(paste("======== Starting Npar=", Npar))

  ## Update data, inits and model objects
  temp <- growth_setup(N=Npar, seed=Npar)
  data <- temp$data
  inits <- temp$inits
  obj.tmb <- MakeADFun(data=data, parameters=inits(), DLL=m)
  setwd('admb')
  write.table(x=unlist(data), file=paste0(m,'.dat'), row.names=FALSE,
                                          col.names=FALSE )
  system(paste(m))
  setwd('..')
  ## Run with updated model size
  temp <- run.chains(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m,
                     data=data, inits=inits, pars=pars, metric='diag', seeds=seeds,
                     Nout=1000, Nthin=1)
  adapt.list[[k]] <- cbind(temp$adapt)
  perf.list[[k]] <- cbind(temp$perf)
  ## Save them as we go in case it fails
  perf <-  do.call(rbind, perf.list)
  adapt <- do.call(rbind, adapt.list)
  plot.simulated.results(perf, adapt)
  write.csv(x=perf, file=paste0(m,'_perf_simulated.csv'))
  write.csv(x=adapt, file=paste0(m,'_adapt_simulated.csv'))
  rm(temp)
  k <- k+1
}



