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
inits <- swallows_setup()$inits
compile('models/swallows/swallows.cpp')
dyn.load('models/swallows/swallows')
obj <-
  MakeADFun(data=data, parameters=inits(),
            random=c('fameffphi_raw', 'fameffp_raw', 'yeareffphi_raw'),
            DLL='swallows')
time.swallows.mle <- system.time(fit.swallows <- with(obj, nlminb(par, fn, gr)))[3]
## Run long, thinned chains to ensure mixing is relatively close
th <- 5 # thin rate
wm <- 1000 # warmup
mcmc.swallows <- tmbstan(obj, iter=1000*th+wm, warmup=wm, thin=th,
                         chains=4, control=list(adapt_delta=.9))
saveRDS(mcmc.swallows, file='results/mcmc.swallows.RDS')
mcmc.swallows.la <- tmbstan(obj, iter=1000*th+wm, warmup=wm, thin=th,
                            chains=4, laplace=TRUE, control=list(adapt_delta=.9))
saveRDS(mcmc.swallows.la, file='results/mcmc.swallows.la.RDS')

## wildf model
set.seed(321)
data <- wildf_setup()$data
inits <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj <- MakeADFun(data=data, parameters=inits(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'), DLL='wildf')
time.wildf.mle <- system.time(fit.wildf <- with(obj, nlminb(par, fn, gr)))[3]
## Run long, thinned chains to ensure mixing is relatively close
th <- 5 # thin rate
wm <- 1000 # warmup
mcmc.wildf <- tmbstan(obj, iter=1000*th+wm, warmup=wm, thin=th, chains=4)
saveRDS(mcmc.wildf, file='results/mcmc.wildf.RDS')
mcmc.wildf.la <- tmbstan(obj, iter=1000*th+wm, warmup=wm, thin=th, chains=4, laplace=TRUE)
saveRDS(mcmc.wildf.la, file='results/mcmc.wildf.la.RDS')


## Read the results in and process them
mcmc.wildf <- readRDS('results/mcmc.wildf.RDS')
mcmc.wildf.la <- readRDS('results/mcmc.wildf.la.RDS')
mcmc.swallows <- readRDS('results/mcmc.swallows.RDS')
mcmc.swallows.la <- readRDS('results/mcmc.swallows.la.RDS')

ess.swallows <- monitor(extract(mcmc.swallows, permuted=FALSE), print=FALSE)[,'n_eff']
ess.swallows.la <- monitor(extract(mcmc.swallows.la, permuted=FALSE), print=FALSE)[,'n_eff']
ess.wildf <- monitor(extract(mcmc.wildf, permuted=FALSE), print=FALSE)[,'n_eff']
ess.wildf.la <- monitor(extract(mcmc.wildf.la, permuted=FALSE), print=FALSE)[,'n_eff']
min(ess.swallows)
min(ess.swallows.la)
min(ess.wildf)
min(ess.wildf.la)

## Make quick pairs plots
temp <- extract(mcmc.swallows, permuted=FALSE)
x1 <- do.call(rbind, lapply(1:4,function(i)  temp[,i,]))
temp <- extract(mcmc.swallows.la, permuted=FALSE)
x2 <- do.call(rbind, lapply(1:4,function(i)  temp[,i,]))
stopifnot(nrow(x1)==nrow(x2))
pars <- dimnames(x2)[[2]]
post <- data.frame(rbind(x1[,pars], x2[,pars]))
model <- as.factor(rep(c("normal", "LA"), each=nrow(x1)))
saveRDS(cbind(model,post), file='results/post.swallows.RDS')
png("plots/pairs_swallows_LA.png", width=7, height=5, units='in', res=200)
## Randomize order to prevent overplotting
ind <- sample(1:nrow(post), size=nrow(post))
post <- post[ind,]
pairs(post[ind,1:10], col=model[ind], pch=16, cex=.5)
dev.off()
## Same for wildf
temp <- extract(mcmc.wildf, permuted=FALSE)
x1 <- do.call(rbind, lapply(1:4,function(i)  temp[,i,]))
temp <- extract(mcmc.wildf.la, permuted=FALSE)
x2 <- do.call(rbind, lapply(1:4,function(i)  temp[,i,]))
stopifnot(nrow(x1)==nrow(x2))
pars <- dimnames(x2)[[2]][1:9]
post <- data.frame(rbind(x1[,pars], x2[,pars]))
model <- as.factor(rep(c("normal", "LA"), each=nrow(x1)))
saveRDS(cbind(model,post), file='results/post.wildf.RDS')
png("plots/pairs_wildf_LA.png", width=7, height=5, units='in', res=200)
## Randomize order to prevent overplotting
ind <- sample(1:nrow(post), size=nrow(post))
post <- post[ind,]
pairs(post[ind,], col=model[ind], pch=16, cex=.5)
dev.off()

## Calculate performance of the LA tests
calc.perf <- function(fit){
  ## Calculate minESS/time for a stanfit obj
  minESS <- min(monitor(extract(fit, permuted=FALSE),
                        print=FALSE)[,'n_eff'])
  time <- max(rowSums(get_elapsed_time(fit)))
  return(minESS/time)
}

## Time is the longest running chain, both warmup + sampling
time.swallows <- max(rowSums(get_elapsed_time(mcmc.swallows)))
time.swallows.la <- max(rowSums(get_elapsed_time(mcmc.swallows.la)))
time.wildf <- max(rowSums(get_elapsed_time(mcmc.wildf)))
time.wildf.la <- max(rowSums(get_elapsed_time(mcmc.wildf.la)))
xx <- cbind(c(time.wildf.mle, time.wildf, time.wildf.la),
      c(NA, calc.perf(mcmc.wildf), calc.perf(mcmc.wildf.la)),
      c(time.swallows.mle, time.swallows, time.swallows.la),
      c(NA, calc.perf(mcmc.swallows), calc.perf(mcmc.swallows.la)))

xx
