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
options(mc.cores = 3) ## parallel TMB runs with tmbstan

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
th <- 10 # thin rate
wm <- 1000 # warmup
mcmc.swallows <- tmbstan(obj, iter=2000*th+wm, warmup=wm, thin=th, chains=3)
## Time is the longest running chain, both warmup + sampling
time.swallows <- max(rowSums(get_elapsed_time(mcmc.swallows)))
mcmc.swallows.la <- tmbstan(obj, iter=2000*th+wm, warmup=wm, thin=th, chains=3, laplace=TRUE)
time.swallows.la <- max(rowSums(get_elapsed_time(mcmc.swallows.la)))

temp <- extract(mcmc.swallows, permuted=FALSE)
x1 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
temp <- extract(mcmc.swallows.la, permuted=FALSE)
x2 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
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

## wildf model
data <- wildf_setup()$data
inits <- wildf_setup()$inits
compile('models/wildf/wildf.cpp')
dyn.load('models/wildf/wildf')
obj <- MakeADFun(data=data, parameters=inits(),
            random=c('yearInterceptEffect_raw',
                     'plantInterceptEffect_raw',
                     'plantSlopeEffect_raw'),
            DLL='wildf')
time.wildf.mle <- system.time(fit.wildf <- with(obj, nlminb(par, fn, gr)))[3]
## Run long, thinned chains to ensure mixing is relatively close
th <- 10 # thin rate
wm <- 1000 # warmup
mcmc.wildf <- tmbstan(obj, iter=2000*th+wm, warmup=wm, thin=th, chains=3)
time.wildf <- max(rowSums(get_elapsed_time(mcmc.wildf)))
mcmc.wildf.la <- tmbstan(obj, iter=2000*thin+wm, warmup=wm, thin=th, chains=3, laplace=TRUE)
time.wildf.la <- max(rowSums(get_elapsed_time(mcmc.wildf.la)))
stopifnot(nrow(x1)==nrow(x2))
temp <- extract(mcmc.wildf, permuted=FALSE)
x1 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
temp <- extract(mcmc.wildf.la, permuted=FALSE)
x2 <- do.call(rbind, lapply(1:3,function(i)  temp[,i,]))
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

