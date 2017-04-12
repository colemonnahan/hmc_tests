## Generate simulated data and run across vector of Npars

adapt.list <- perf.list <- list()
k <- 1
for(i in seq_along(Npar.vec)){
  Npar <- Npar.vec[i]
  message(paste("======== Starting Npar=", Npar))
  for(j in cor.vec){
    ## Reproducible data since seed set inside the function
    set.seed(115)
    ## Update data, inits and model objects
    df <- max(Npar.vec)
    if(j==0){
      covar <- diag(Npar)
    } else if(j==1){
      covar <- rWishart(n=1, df=df, Sigma=diag(df))[1:Npar,1:Npar,1]
    } else {
      stop("Invalid cor value in mvnd")
    }
    data <- list(Npar=Npar, covar=covar, x=rep(0, len=Npar))
    inits <- lapply(1:length(seeds), function(i)
      list(mu=as.vector(mvtnorm::rmvnorm(n=1, mean=rep(0, len=Npar), sigma=covar))))

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
                       covar=covar,
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
}


