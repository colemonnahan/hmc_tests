## Script to demonstrate testing the Laplace approximation of a TMB
## model. Verisons are run with and without the LA turned on for the random
## effects, and differences in posteriors for fixed effects are compared.

## Working directory needs to be set to the folder of this file before
## proceeding.

## Note: these models include priors, and use exponentiation in the
## template with a Jacobian adjustment instead of external bounds of
## (0,Inf) for hypervariance parameters. Thus, the SD parameters are really
## in log space. See model files for more details.

source("startup.R")
options(mc.cores = 4) ## parallel TMB runs with tmbstan

## swallows model
data <- swallows_setup()$data
inits.swallows <- swallows_setup()$inits
compile('models/swallows/swallows.cpp')
dyn.load('models/swallows/swallows')
obj.swallows <-
  MakeADFun(data=data, parameters=inits.swallows(),
            random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'),
            DLL='swallows')
## wildf model
set.seed(321)
data <- wildf_setup()$data
inits.wildf <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj.wildf <- MakeADFun(data=data, parameters=inits.wildf(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'), DLL='wildf')
## MLEs with timing
time.swallows.mle <- system.time(fit.swallows <- with(obj.swallows, nlminb(par, fn, gr)))[3]
time.wildf.mle <- system.time(fit.wildf <- with(obj.wildf, nlminb(par, fn, gr)))[3]


## Run default settings to test for performance
mcmc.swallows <- tmbstan(obj.swallows, init=inits.swallows,
                         chains=4, control=list(adapt_delta=.9))
mcmc.swallows.la <- tmbstan(obj.swallows, init=inits.swallows, laplace=TRUE,
                            chains=4, control=list(adapt_delta=.9))
mcmc.wildf <- tmbstan(obj.wildf, init=inits.wildf,
                         chains=4, control=list(adapt_delta=.9))
mcmc.wildf.la <- tmbstan(obj.wildf, init=inits.wildf, laplace=TRUE,
                            chains=4, control=list(adapt_delta=.9))
## Calculate performance of the LA tests
calc.perf <- function(fit, model, version){
  ## Calculate minESS/time for a stanfit obj
  minESS <- min(monitor(extract(fit, permuted=FALSE),
                        print=FALSE)[,'n_eff'])
  ## Time assumed to be the sum of sampling + warmup for all chains
  time <- sum(get_elapsed_time(fit))
  return(data.frame(model=model, version=version, time=time,
                    ESS=minESS,efficiency=minESS/time))
}
x1 <- data.frame(model='swallows', version='MLE', time=time.swallows.mle, ESS=NA, efficiency=NA)
x2 <- calc.perf(mcmc.swallows, 'swallows', 'Full Bayesian')
x3 <- calc.perf(mcmc.swallows.la, 'swallows', 'Laplace')
x4 <- data.frame(model='wildflower', version='MLE', time=time.wildf.mle, ESS=NA, efficiency=NA)
x5 <- calc.perf(mcmc.wildf, 'wildflower', 'Full Bayesian')
x6 <- calc.perf(mcmc.wildf.la, 'wildflower', 'Laplace')
table.la.perf <- rbind(x1,x2,x3,x4,x5,x6)
write.csv(table.la.perf, file='results/table.la.perf.csv')
rm(x1,x2,x3,x4,x5,x6)

## Re-run long, thinned chains to ensure mixing is relatively close for
## comparing marginal distributions. Based on table.la.perf I'm thinning
## the versions at different rates. THe LA ones don't need as much and take
## longer to run.
wm <- 1000 # warmup
mcmc.swallows <- tmbstan(obj.swallows, iter=1000*10+wm, warmup=wm, thin=10,
                         init=inits.swallows,
                         chains=4, control=list(adapt_delta=.9))
saveRDS(mcmc.swallows, file='results/mcmc.swallows.RDS')
mcmc.swallows.la <- tmbstan(obj.swallows, iter=1000*4+wm, warmup=wm, thin=4,
                            init=inits.swallows,
                            chains=4, laplace=TRUE, control=list(adapt_delta=.9))
saveRDS(mcmc.swallows.la, file='results/mcmc.swallows.la.RDS')
mcmc.wildf <- tmbstan(obj.wildf, iter=1000*10+wm, warmup=wm, thin=10,
                      init=inits.wildf,
                      chains=4, control=list(adapt_delta=.9))
saveRDS(mcmc.wildf, file='results/mcmc.wildf.RDS')
mcmc.wildf.la <- tmbstan(obj.wildf, iter=1000*2+wm, warmup=wm, thin=2,
                         init=inits.wildf, control=list(adapt_delta=.9),
                         chains=4, laplace=TRUE)
saveRDS(mcmc.wildf.la, file='results/mcmc.wildf.la.RDS')

## Read the results in and process them
mcmc.wildf <- readRDS('results/mcmc.wildf.RDS')
mcmc.wildf.la <- readRDS('results/mcmc.wildf.la.RDS')
mcmc.swallows <- readRDS('results/mcmc.swallows.RDS')
mcmc.swallows.la <- readRDS('results/mcmc.swallows.la.RDS')

rbind(calc.perf(mcmc.swallows, 'swallows', 'Full Bayesian'),
      calc.perf(mcmc.swallows.la, 'swallows', 'Laplace'),
      calc.perf(mcmc.wildf, 'wildflower', 'Full Bayesian'),
      calc.perf(mcmc.wildf.la, 'wildflower', 'Laplace'))

## Make quick pairs plots
x1 <- as.matrix(mcmc.swallows)
x2 <- as.matrix(mcmc.swallows.la)
stopifnot(nrow(x1)==nrow(x2))
post <- data.frame(rbind(x1[,dimnames(x2)[[2]]], x2[,dimnames(x2)[[2]]]))
post$model <- as.factor(rep(c("normal", "LA"), each=nrow(x1)))
saveRDS(post, file='results/swallows.la.post.RDS')
png("plots/pairs_swallows_LA.png", width=7, height=5, units='in', res=200)
## Randomize order to prevent overplotting
ind <- sample(1:nrow(post), size=nrow(post))
post <- post[ind,]
pairs(post[,1:10], col=post$model, pch=16, cex=.5)
dev.off()
## Same for wildf
x1 <- as.matrix(mcmc.wildf)
x2 <- as.matrix(mcmc.wildf.la)
stopifnot(nrow(x1)==nrow(x2))
post <- data.frame(rbind(x1[,dimnames(x2)[[2]]], x2[,dimnames(x2)[[2]]]))
post$model <- as.factor(rep(c("normal", "LA"), each=nrow(x1)))
saveRDS(post, file='results/wildf.la.post.RDS')
png("plots/pairs_wildf_LA.png", width=7, height=5, units='in', res=200)
## Randomize order to prevent overplotting
ind <- sample(1:nrow(post), size=nrow(post))
post <- post[ind,]
pairs(post[,1:9], col=post$model, pch=16, cex=.5)
dev.off()


