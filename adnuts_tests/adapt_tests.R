## Test that step size and mass matrix adaptation is consistent across TMB
## and ADMB. Start with simple MVN model.
## devtools::install("../adnuts")
source("startup.R")
set.seed(235)
Npar <- 2
corr <- matrix(.8, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
#covar <- diag(2)
covar2 <- diag(x=diag(covar))
covar.inv <- solve(covar)
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
compile(file='mvnd.cpp')
dyn.load(dynlib('mvnd'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd')
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd -hbf 1')
setwd('..')
stan.fit <- stan(file='mvnd.stan', data=data, chains=1, iter=1000)
extract_adapt <- function(fit){
  ldply(1:length(fit$sampler_params), function(i){
    ind <- -(1:fit$warmup)
    xx <- fit$sampler_params[[i]]
    cbind(chain=i, stepsize=tail(xx[ind,2],1), nsteps.median= median(xx[ind,4]),
          nsteps.mean=mean(xx[ind,4]),
          accept.prob=mean(xx[ind,1]),
          stepsize.final=head(xx[ind,2],1))})
}

## Run across fixed step size to see if matches, then with adaptation to
## see if step size is similar
set.seed(2115)
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')
devtools::load_all('C:/Users/Cole/adnuts')
iter <- 1100
chains <- 1
td <- 12
eps <- 1.21
rm(admb, tmb, stan2)
init <- sapply(1:chains, function(x) list(mu=c(0,0)))
stan2 <- stan(fit=stan.fit, data=data, chains=chains, iter=iter,
              init=list(init),
              control=list(metric='unit_e', stepsize=eps,
                           adapt_engaged=FALSE, max_treedepth=td))
admb <- sample_admb(dir=dir, model=model, iter=iter, init=init,
                        chains=chains, control=list(stepsize=eps, max_treedepth=td,
                                                    metric=diag(2)))
tmb <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, warmup=iter-100,
                 control=list(stepsize=eps, metric=diag(2), max_treedepth=td))
x1 <- data.frame(admb$sampler_params[[1]],platform='admb')
x2 <- data.frame(tmb$sampler_params[[1]], platform='tmb')
x3 <- data.frame(get_sampler_params(stan2)[[1]],platform='stan')
plot(x2$n_leapfrog__)
plot(x2$accept_stat__)
plot(x2$stepsize__)
plot(tmb$samples[,1,2])
nevals <- cbind(rbind(x1,x2,x3),i=1:iter)
## ggplot(nevals, aes(i, log2(n_leapfrog__), color=platform)) + geom_jitter(alpha=.5)
## ggplot(nevals, aes(log2(n_leapfrog__))) +
##   geom_histogram() + facet_wrap('platform')
os <- .3
setwd('../../')
png('plots/stan_tmb_admb_comparison.png', width=9, height=6.5, units='in', res=300)
par(mfrow=c(2,2), mar=c(4,2.5,.5,.5), mgp=c(1.5,.5,0), oma=c(0,0,1.5,0))
plot(sort(x1[,4])+2*os, ylim=c(0,1+max(nevals$n_leapfrog__)),
     ylab='n_leapfrog', xlab='Sorted iteration')
points(sort(x2[,4]), col=2)
points(sort(x3[,4])-2*os, col=3)
plot(sort(x1[,1]), ylab='acceptance probability', , xlab='Sorted iteration')
points(sort(x2[,1]), col=2)
points(sort(x3[,1]), col=3)
legend('bottomright', legend=c('admb', 'tmb', 'stan'), pch=16, col=1:3)
plot(sort(x1[,3])+os, ylim=c(0, 1+max(nevals$treedepth__)),
     ylab='treedepths', xlab='Sorted iteration')
points(sort(x2[,3]), col=2)
points(sort(x3[,3])-os, col=3)
admb.samples <- extract_samples(admb)
tmb.samples <- extract_samples(tmb)
stan.samples <- data.frame(extract(stan2, permuted=FALSE)[,1,])
qqplot(stan.samples[,1], admb.samples[,1], main=NA, ylab=NA, xlab='Stan')
mtext("QQ Plot", line=-2)
x <- qqplot(stan.samples[,1], tmb.samples[,1], plot=FALSE)
points(x[[1]], x[[2]], col=2)
mtext("Trajectory info for MVN w/ fixed eps", line=0, outer=TRUE)
dev.off()
##ggplot(nevals, aes(treedepth__, n_leapfrog__, color=platform)) + geom_jitter()

## --------------- Same test but with step size adaptation
set.seed(2115)
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')
## devtools::load_all('C:/Users/Cole/adnuts')
iter <- 500
warmup <- 400
chains <- 20
td <- 14
eps <- NULL
rm(admb, tmb, stan2)
init <- sapply(1:chains, function(x) list(mu=c(0,0)))
stan2 <- stan(fit=stan.fit, data=data, chains=chains, iter=iter, warmup=warmup,
              control=list(metric='unit_e',
                            max_treedepth=td))
admb <- sample_admb(dir=dir, model=model, iter=iter, init=init, warmup=warmup,
                        chains=chains, control=list(stepsize=eps, max_treedepth=td,
                                                    metric=diag(2)))
tmb <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, warmup=warmup,
                 control=list(stepsize=eps, metric=diag(2), max_treedepth=td))
x1 <- data.frame(do.call(rbind, admb$sampler_params),platform='admb',i=1:iter)
x2 <- data.frame(do.call(rbind, tmb$sampler_params),platform='tmb',i=1:iter)
x3 <- data.frame(do.call(rbind, get_sampler_params(stan2)),platform='stan',i=1:iter)
nevals <- cbind(rbind(x1,x2,x3))
setwd('../../')
g <- ggplot(subset(nevals, i>5), aes(i, stepsize__)) + ylab("Log10 Stepsize")+
  geom_point(alpha=.1, size=.25) + facet_wrap('platform')+ scale_y_log10() +
  ggtitle("Adaptation over 20 replicates for iid bivariate Normal") + xlab("Iteration")
ggsave('plots/stepsize_adaptation.png', g, width=7, height=4)

## Compare effect of mass matrix on a MVN model. I can't pass a matrix to
## Stan so cannot compare easily here.
setwd('C:/Users/Cole/hmc_tests/models/mvnd')
Npar <- 5
corr <- matrix(.98, nrow=Npar, ncol=Npar)
diag(corr) <- 1
se <- seq(.01, .5, len=Npar)
covar <- corr * (se %o% se)
covar2 <- diag(x=diag(covar))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=rep(0, Npar)), DLL='mvnd')
chains <- 15
init <- lapply(1:chains, function(x) rnorm(Npar, sd=se*5))
iter <- 1000
tmb1 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=diag(Npar)))
tmb2 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=covar2))
tmb3 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=covar))
adapt.tmb1 <- cbind(form=1, extract_adapt(tmb1))
adapt.tmb2 <- cbind(form=2, extract_adapt(tmb2))
adapt.tmb3 <- cbind(form=3, extract_adapt(tmb3))
adapt.tmb <- rbind(adapt.tmb1, adapt.tmb2, adapt.tmb3)
admb1 <- sample_admb(dir=dir, model, iter=iter, init=init, chains=chains, control=list(metric=diag(Npar)))
admb2 <- sample_admb(dir=dir, model, iter=iter, init=init, chains=chains, control=list(metric=covar2))
admb3 <- sample_admb(dir=dir, model, iter=iter, init=init, chains=chains, control=list(metric=covar))
adapt.admb1 <- cbind(form=1, extract_adapt(admb1))
adapt.admb2 <- cbind(form=2, extract_adapt(admb2))
adapt.admb3 <- cbind(form=3, extract_adapt(admb3))
adapt.admb <- rbind(adapt.admb1, adapt.admb2, adapt.admb3)
yy1 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar='unit', m='admb', as.data.frame(admb1$sampler_params[[i]])))
xx1 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar='unit', m='tmb', as.data.frame(tmb1$sampler_params[[i]])))
yy2 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar='diag', m='admb', as.data.frame(admb2$sampler_params[[i]])))
xx2 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar='diag', m='tmb', as.data.frame(tmb2$sampler_params[[i]])))
yy3 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter,  covar='dense',m='admb', as.data.frame(admb3$sampler_params[[i]])))
xx3 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar='dense', m='tmb', as.data.frame(tmb3$sampler_params[[i]])))
zz <- rbind(yy1,yy2,yy3, xx1,xx2,xx3)
zz <- zz[sample(1:nrow(zz), size=nrow(zz)),]
adapt <- rbind(cbind(name='tmb', adapt.tmb),
               cbind(name='admb', adapt.admb))
adapt$nsteps.median <- NULL
adapt.long <- reshape2::melt(adapt, c("name", "form", "chain"))
adapt.long$form <- factor(adapt.long$form)
levels(adapt.long$form) <- c("Unit", "Diag", "Dense")
## Effective sample sizes should be the same too
x <- data.frame(m='admb', form='unit', ess=admb1$ess)
y <- data.frame(m='admb', form='diag', ess=admb2$ess)
z <- data.frame(m='admb', form='dense', ess=admb3$ess)
library(rstan)
xx <- data.frame(m='tmb', form='unit', ess=monitor(tmb1$samples, warmup=tmb1$warmup)[,'n_eff'])
yy <- data.frame(m='tmb', form='diag', ess=monitor(tmb2$samples, warmup=tmb1$warmup)[,'n_eff'])
zz <- data.frame(m='tmb', form='dense', ess=monitor(tmb3$samples, warmup=tmb1$warmup)[,'n_eff'])
ess.all <- rbind(x,y,z,xx,yy,zz)
ess.all$ess.par <- row.names(ess.all)
ess.all <- ess.all[grep('mu', x=row.names(ess.all)),]
setwd('../..')

## Make final plots. These should match. YES!!!
g <- ggplot(subset(zz, iter>5), aes(iter, log10( stepsize__), group=seed, color=m)) +
  facet_grid(covar~., scales='free') + geom_point(alpha=.15, size=.5)
ggsave("plots/adapt_tests_stepsize.png", g, width=7, height=4)
g <- ggplot(adapt.long, aes(form, value, color=name)) + geom_point(alpha=.5, position = position_jitter(w = 0.3, h = 0)) +
  facet_wrap('variable', scales='free') + scale_y_log10()
ggsave("plots/adapt_tests.png", g, width=7, height=5)
g <- ggplot(ess.all, aes(form, ess, color=m)) +
  geom_point(alpha=.5, position=position_jitter(w = 0.3, h = 0)) #+ scale_y_log10
ggsave("plots/adapt_tests_ess.png", g, width=7, height=5)

### --------------------------------------------------
## Test mass matrix adaptation
## Run across fixed step size to see if matches, then with adaptation to
## see if step size is similar
devtools::install("c:/Users/Cole/adnuts")
set.seed(2115)
dir <- 'admb'; model <- 'mvnd'
chains <- 5
seeds <- 1:chains
cores <- 5
td <- 12
eps <- NULL
init <- sapply(1:chains, function(x) list(mu=c(0,0)))
setwd('admb')
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd -hbf 1')
setwd('..')


## Loop across warmup lengths and see how estmiated step size changes
output <- ldply(100*2^(0:5), function(warmup){
  ldply(c(TRUE, FALSE), function(am) {
    tmb <- sample_tmb(mvnd.obj, iter=warmup+200, seeds=seeds, init=init, chains=chains, warmup=warmup,
                      parallel=TRUE, cores=cores,
                      control=list(stepsize=eps, max_treedepth=td, adapt_mass=am))
    tmb.avg.steps <- sapply(1:chains, function(i)
      mean(tmb$sampler_params[[i]][-(1:warmup),4]))
    tmb.eps <- sapply(1:chains, function(i)
      tail(tmb$sampler_params[[i]][,2],n=1))
    x1 <- data.frame(platform='tmb', warmup=warmup, am=am, seeds=seeds,
                     runtime=tmb$time.total, steps=tmb.avg.steps, eps=tmb.eps)
    mm <- if(am) NULL else 'unit'
    admb <- sample_admb('mvnd', path='admb', iter=warmup+200, seeds=seed, init=init, chains=chains, warmup=warmup,
                      parallel=TRUE, cores=cores,
                      control=list(metric=mm, adapt_mass=am))
    admb.avg.steps <- sapply(1:chains, function(i)
      mean(admb$sampler_params[[i]][-(1:warmup),4]))
    admb.eps <- sapply(1:chains, function(i)
      tail(admb$sampler_params[[i]][,2],n=1))
    x2 <- data.frame(platform='admb', warmup=warmup, am=am, seeds=seeds,
                     runtime=admb$time.total, steps=admb.avg.steps, eps=admb.eps)
    return(rbind(x1,x2))
    })
  })

output.long <- melt(output, id.vars=c('warmup', 'am', 'seeds', 'platform'))
output.long <- subset(output.long, variable!='runtime')
g <- ggplot(output.long, aes(log2(warmup), y=value, color=am, shape=platform)) +
  geom_jitter(alpha=.75, height=0, width=.2) + facet_grid('variable~.', scales='free')
g <- g + scale_y_log10()
setwd('C:/Users/Cole/hmc_tests/')
ggsave('plots/mass_matrix_adaptation.png', g, width=7, height=3.5)

