## Run long chains for individual models using RWM.
devtools::install('C:/Users/Cole/adnuts')
library(shinystan)
library(adnuts)
library(snowfall)

reps <- 6 # chains to run in parallel
hh <- 24 # hours to run

sfStop()
d <- m <- 'cod_fast'
thin <- 100
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM', extra.args=' -ainp cod_fast.par')

saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=3, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS',
              control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'hake'
thin <- 1000
iter <- 4000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))

sfStop()
d <- m <- 'tanner'
thin <- 1
iter <- 1000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
             parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
               parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'halibut'
thin <- 1000
iter <- 1000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.nuts <- sample_admb(m, iter=iter, init=inits,  thin=1,
              parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.95))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))
d <- m <- 'halibut2'
inits <- NULL
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
i <- 2
temp <- sample_admb(m, iter=100*(2^i), init=inits, par.names=par.names, thin=1,
               parallel=TRUE, chains=reps, warmup=25*(2^i),
              dir=d, cores=reps, algorithm='NUTS',
              control=list(metric='unit', max_treedepth=5))
fit.nuts <- sample_admb(m, iter=iter, init=inits, thin=1,
              parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS',
              control=list(adapt_delta=.8, metric='unit'))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))




sfStop()
d <- 'snowcrab'; m <- '2016sc'
thin <- 1000
iter <- 1000
warmup <- iter/4
mle <- get_m(paste0(d,'/',m))
N <- mle$nopar
par.names <- paste0(1:N, "_", mle$names[1:N])
## Draw inits from MVT using MLE and covar
covar <- get.admb.cov(d)$cov.bounded
inits <- lapply(1:reps, function(i) mle$est[1:N]+as.vector(mvtnorm::rmvt(n=1, df=50, sigma=covar)))
inits <- NULL
inits <- lapply(1:reps, function(i) mle$est[1:N])
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, par.names=par.names, thin=thin,
              duration=hh*60, parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
fit.rwm <- readRDS('results/long_rwm_2016sc.RDS')
## This model has all sorts of bounds issues which derails NUTS b/c the
## initial value is inf due to differences in bound functions. So pick a
## more reasonable place to start and bootstrap up to a good mass matrix.
mle.adjusted <- mle$est
## selsmo10ind(6) and (22) are right at upper bound of -0.001
mle.adjusted[c(293,309)] <- -.01
## selsmo09ind(1) (2) and (3) are  at lower bound of -4
mle.adjusted[310:312] <- -3.9
#selsmo09ind (14) and (15) and (22)  at upper bound of -.001
mle.adjusted[c(323, 324,331)] <- -.01
mle.adjusted[6] <- 1.1
mle.adjusted[268:269] <- .9
mle.adjusted[284] <- 57
mle.adjusted[286:287] <- .9
inits <- lapply(1:reps, function(i) mle.adjusted)
temp <- list(covar.est=diag(N))
for(i in 1:3){
temp <- sample_admb(m, iter=100*(2^i), init=inits, par.names=par.names, thin=1,
               parallel=TRUE, chains=reps, warmup=25*(2^i),
              dir=d, cores=reps, algorithm='NUTS',
              control=list(metric=temp$covar.est, max_treedepth=5))
}
fit.nuts <- sample_admb(m, iter=iter, init=inits, par.names=par.names, thin=1,
              parallel=TRUE, chains=reps, warmup=warmup,
              dir=d, cores=reps, algorithm='NUTS',
              control=list(metric=temp$covar.est, adapt_delta=.9))
saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))


sfStop()
d <- m <- 'canary'
thin <- 1000
iter <- 1000
warmup <- iter/4
mle <- r4ss::read.admbFit(paste0(d,'/',m))
N <- mle$nopar
inits <- lapply(1:reps, function(i) mle$est[1:N])
temp <- c(0.667932, 9.04907, 60.0444, 0.128736, 0.109042, -1.35737, -0.114137,
  0.224431, 0.0420364, -1.18266, 1.10623, 0.542608, -0.00392423,
  -0.00354755, -0.00311344, -0.00272063, -0.0024948, -0.0024558,
  -0.00249081, -0.0024877, -0.00199517, -0.00115078, 0.000488141,
  0.000310569, -0.000304327, -0.00207686, -0.00518243, -0.00984975,
  -0.0159893, -0.0241392, -0.0341713, -0.0459476, -0.0604696, -0.0681047,
  -0.0825857, -0.0970526, -0.0994629, -0.0940686, -0.0676346, -0.0254009,
  0.04378, 0.0791615, 0.157366, 0.177877, 0.205543, 0.245198, 0.204598,
  -0.130925, -0.236923, -0.258321, -0.317246, -0.320321, -0.050177,
  0.0520908, 0.0119062, -0.320327, -0.392566, -0.369134, -0.0138708,
  -0.22056, -0.501851, -0.347695, 0.0639402, -0.37856, 0.143644, 0.0431916,
  -0.209958, -0.0459135, -0.114176, 0.124875, -0.0869618, -0.116457,
  -0.158047, -0.0363377, 0.427613, 0.078282, 0.180554, -0.0681416,
  -0.246816, -0.355915, -0.656147, 0.278408, 0.550664, 1.00967, -0.0336861,
  0.880991, 0.260892, 1.48685, -0.507539, -0.0795356, 1.09422, -0.113079,
  -0.211892, 0, -0.049647, -0.0521108, -0.0546626, -0.0573545, -0.060259,
  -0.0634185, -0.0668194, -0.070451, -0.0744166, -0.0786126, -0.0832015,
  -0.087443, -0.0918219, -0.0959709, -0.100133, -0.104205, -0.107923,
  -0.111061, -0.113767, -0.116248, -0.118316, -0.125002, -0.131235,
  -0.140237, -0.156747, -0.184097, -0.224445, -0.23219, -0.314068,
  -0.34967, -0.364081, -0.30802, -0.282229, -0.282187, -0.307438,
  -0.570209, -0.384907, -0.321725, -0.320619, -0.463101, -0.478543,
  -0.449394, -0.478213, -0.603201, -0.362558, 0.0307809, 0.350292,
  -0.136224, -0.581776, -0.327881, -0.10776, -0.0211126, 0.351443,
  0.186848, 0.30852, 0.319298, 0.466284, 0.14346, 0.518141, 0.0542307,
  0.0319405, 1.3121, 0.248275, 0.00419739, 0.883068, -0.0967021,
  0.00891414, 0.460654, -0.270442, -0.302727, 1.3948, -0.0421638, 1.99587,
  1.19578, -0.286724, 2.44091, -0.441265, 0.320963, 0.174421, -0.209732,
  0.317509, 0, 8.03694, -0.126693, -0.132112, -0.137838, -0.14397,
  -0.150527, -0.157492, -0.164792, -0.172239, -0.179036, -0.18566,
  -0.192412, -0.201062, -0.208986, -0.216168, -0.223253, -0.230588,
  -0.238119, -0.245058, -0.251484, -0.256436, -0.26096, -0.24543,
  -0.230374, -0.200049, -0.14984, -0.0436263, 0.0795042, 0.0269414,
  0.33706, 0.426508, 0.450238, 0.146275, -0.0287136, -0.0192518, 0.149992,
  0.967852, 0.252482, -0.200305, -0.290326, 0.243858, 0.373395, 0.301688,
  0.231147, 0.174543, -0.307435, 0.38436, 0.876661, -1.08124, 0.331093,
  -0.172887, -0.616866, 0.479902, -0.171624, -0.142624, 0.160288,
  0.0245147, 0.314947, 0.0480755, 0.179072, 0.241305, -0.158219, 0.362694,
  0.26286, 0.00445095, 0.37694, -0.228295, -0.151966, -0.285628, 0.226225,
  -0.00956329, 0.14378, -0.382053, -0.961338, -0.461078, 0.203337,
  -0.941728, -0.0460424, 0.285251, -0.480272, -0.667028, -0.778699,
  -0.404558, 0, 0, 0.513864, 0.390702, 0.985528, 48.8625, 4.29864,
  0.0118053, 1.18783, 43.23, 5.00815, 7.32338, 4.59544, 32.3961, 3.75467,
  3.40737, -1.66833, 45.9363, 2.99353, 0.00642784, 0.605525, 13.017,
  4.51006, 8.93911, 2.28129, 53.9796, 6.55668, 0.424443, 0.556848, 51.0633,
  5.66468, 4.31027, -0.25247, 43.6794, 43.2385, 4.25186, 5.03605, 1.58579,
  1.58691, 0.188873, -0.778762, 40.6513, 4.31715, 6.61523, 4.21022)
inits <- lapply(1:reps, function(i) temp)
sfInit(parallel=TRUE, cpus=reps)
sfExportAll()
fit.rwm <- sample_admb(m, iter=iter*thin, init=inits, thin=thin,
              parallel=TRUE, chains=reps, warmup=warmup*thin,
              dir=d, cores=reps, algorithm='RWM')
launch_shinyadmb(fit.rwm)
saveRDS(fit.rwm, file=paste0("results/long_rwm_", m, ".RDS"))
## fit.nuts <- sample_admb(m, iter=iter, init=inits,  thin=1,
##               parallel=TRUE, chains=reps, warmup=warmup,
##               dir=d, cores=reps, algorithm='NUTS', control=list(adapt_delta=.95))
## saveRDS(fit.nuts, file=paste0("results/long_nuts_", m, ".RDS"))

