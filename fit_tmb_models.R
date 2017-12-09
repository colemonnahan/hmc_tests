## Script to do MLE fits for the example models

## Working directory needs to be set to the folder of this file before
## proceeding.

## Note: these models include priors, and use exponentiation in the
## template with a Jacobian adjustment instead of external bounds of
## (0,Inf) for hypervariance parameters. Thus, the SD parameters are really
## in log space. See template for more details. Should we turn the priors
## off for MLE?

source("startup.R")

## swallows model
data <- swallows_setup()$data
inits <- swallows_setup()$inits
compile('models/swallows/swallows.cpp')
dyn.load('models/swallows/swallows')
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'),
            DLL='swallows')
fit.swallows <- with(obj, nlminb(par, fn, gr))
mcmc.swallows <- tmbstan(obj, iter=2000, chains=3)

## wildf model
data <- wildf_setup()$data
inits <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'),
            DLL='wildf')
fit.wildf <- with(obj, nlminb(par, fn, gr))
mcmc.wildf <- tmbstan(obj, iter=2000, chains=3)
mcmc.wildf.la <- tmbstan(obj, iter=2000, chains=3, laplace=TRUE)

temp <- extract(mcmc.wildf, permuted=FALSE)
x1 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
temp <- extract(mcmc.wildf.la, permuted=FALSE)
x2 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
pars <- dimnames(x2)[[2]][1:9]
post <- data.frame(rbind(x1[,pars], x2[,pars]))
model <- as.factor(rep(c("normal", "LA"), each=3000))
png("pairs_wildf_LA.png", width=7, height=5, units='in', res=200)
pairs(post, col=model, pch='.')
dev.off()

png("qqplots_wildf_LA.png", width=7, height=5, units='in', res=200)
par(mfrow=c(3,3), mar=c(3,3,.5,.5), oma=c(2,2,2,0))
for(i in 1:9){
  qqplot(x=x1[,i], y=x2[,i], main=NA, xlab=NA, ylab=NA)
  mtext(pars[i])
  ## qqline(x=x1[,i], y=x2[,i])
  abline(a=0, b=1, lwd=2, col='red')
  }
mtext("Full Integration", side=1, outer=TRUE)
mtext("LA Integration", side=2, outer=TRUE)
dev.off()

## Quickly rerun to plot MCMC vs MLE for fixed effects
mcmc.wildf2 <- sample_tmb(obj, iter=5000, warmup=1000, init=inits, chains=3, cores=3,
                          parallel=TRUE, path='models/wildf', control=list(adapt_delta=.59)
fit.wildf <- with(obj, nlminb(par, fn, gr))
## Have to manually modify the list so that pairs_admb can use it. It's not
## designed for tmb output.
sd <- sdreport(obj)
se <- sqrt(diag(sd$cov.fixed))
cor <- sd$cov.fixed / se %o% se
mle <- list(par.names=dimnames(mcmc.wildf2$samples)[[3]][1:9],
            est=fit.wildf$par, se=se, cor=cor)
mcmc.wildf2$mle <- mle
png('pairs_wildf.png', width=6, height=4, units='in', res=500)
pairs_admb(mcmc.wildf2, pars=mcmc.wildf2$mle$par.names)
dev.off()
