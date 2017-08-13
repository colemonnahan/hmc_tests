## Generate simulated data
stop("No data generation setup")

adapt.list <- perf.list <- list()
k <- 1
j <- 0
for(i in seq_along(Npar.vec)){
  Npar <- Npar.vec[i]
  message(paste("======== Starting Npar=", Npar))
  ## Reproducible data since seed set inside the function
  set.seed(115)

  ## Update data, inits and model objects
  data <- NULL
  inits <- NULL # needs to be list of length seeds

  obj.stan <- stan(fit=obj.stan, data=data, iter=100,
                   chains=1, init=list(inits[[1]]), verbose=FALSE)
  obj.tmb <- MakeADFun(data=data, parameters=inits[[1]])
  setwd('admb')
  write.table(x=unlist(data), file=paste0(m,'.dat'), row.names=FALSE,
                                          col.names=FALSE )
  system(m); setwd('..')
  ## Run with updated model size
  temp <- run.chains(obj.stan=obj.stan, obj.tmb=obj.tmb, model=m,
                     inits=inits, pars=pars, data=data, metric='diag',
                     seeds=seeds, Nout=Nout, Nthin=1, delta=delta)
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



