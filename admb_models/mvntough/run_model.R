## Generate random covariance matrix with bounds
setwd('mvntough')
Npar <- 25
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
covar <- cov2cor(covar)
lwr <- -3
upr <- 3
data <- list(Npar=Npar, covar=covar, x=rep(0,Npar), lwr=lwr, upr=upr)
pars <- 'mu'
write.table(x=unlist(data), file='mvntough.dat', row.names=FALSE, col.names=FALSE)
setwd('..')


