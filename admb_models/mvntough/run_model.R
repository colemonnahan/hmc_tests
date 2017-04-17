## Generate random covariance matrix with bounds

setwd('mvntough')

Npar <- 100
covar <- rWishart(n=1, df=Npar, Sigma=diag(Npar))[,,1]
covar <- cov2cor(covar)
#covar <- diag(Npar)
lwr <- -5
upr <- 5
data <- list(Npar=Npar, covar=covar, lwr=lwr, upr=upr)
pars <- 'mu'
write.table(x=unlist(data), file='mvntough.dat', row.names=FALSE, col.names=FALSE)

setwd('..')


