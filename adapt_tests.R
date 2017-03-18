## Test that step size and mass matrix adaptation is consistent across TMB
## and ADMB. Start with simple MVN model.
## devtools::install("../adnuts")
source("startup.R")
set.seed(235)
Npar <- 2
corr <- matrix(.98, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
covar <- diag(2)
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
system('mvnd')
setwd('..')
stan.fit <- stan(file='mvnd.stan', data=data, chains=1, iter=10)
extract_adapt <- function(fit){
  ldply(1:length(fit$sampler_params), function(i){
    ind <- -(1:fit$warmup)
    xx <- fit$sampler_params[[i]]
    cbind(chain=i, stepsize=tail(xx[ind,2],1), nsteps.median= median(xx[ind,4]),
          nsteps.mean=mean(xx[ind,4]),
          accept.prob=mean(xx[ind,1]),
          stepsize0=head(xx[ind,2],1))})
}

## Run across fixed step size to see if matches
dir <- 'admb'; model <- 'mvnd'
setwd(dir)
write.table(x=c(Npar, covar), file='mvnd.dat', row.names=FALSE, col.names=FALSE)
system('admb mvnd')
system('mvnd')
setwd('..')
devtools::load_all('C:/Users/Cole/adnuts')
iter <- 1000
chains <- 1
td <- 4
eps <- 1.312/10
rm(admb, tmb, stan2)
init <- sapply(1:chains, function(x) list(mu=c(0,0)))
stan2 <- stan(fit=stan.fit, data=data, chains=chains, iter=iter,
              init=list(init),
              control=list(metric='unit_e', stepsize=eps,
                           adapt_engaged=FALSE, max_treedepth=td))
admb <- sample_admb(dir=dir, model=model, iter=iter, init=init,
                        chains=chains, control=list(stepsize=eps, max_treedepth=td,
                                                    metric=diag(2)))
tmb <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains,
                 control=list(stepsize=eps, metric=diag(2), max_treedepth=td))
x1 <- data.frame(admb$sampler_params[[1]],platform='admb')
x2 <- data.frame(tmb$sampler_params[[1]], platform='tmb')
x3 <- data.frame(get_sampler_params(stan2)[[1]],platform='stan')
nevals <- cbind(rbind(x1,x2,x3),i=1:iter)
## ggplot(nevals, aes(i, log2(n_leapfrog__), color=platform)) + geom_jitter(alpha=.5)
## ggplot(nevals, aes(log2(n_leapfrog__))) +
##   geom_histogram() + facet_wrap('platform')
par(mfrow=c(2,2))
plot(sort(x1[,4])+.1, ylim=c(0,1+max(nevals$n_leapfrog__)), ylab='n_leapfrog')
points(sort(x2[,4]), col=2)
points(sort(x3[,4])-.1, col=3)
plot(sort(x1[,1]), ylab='acceptance probability')
points(sort(x2[,1]), col=2)
points(sort(x3[,1]), col=3)
legend('bottomright', legend=c('admb', 'tmb', 'stan'), pch=16, col=1:3)
plot(sort(x1[,3])+.1, ylim=c(0, 1+max(nevals$treedepth__)), ylab='treedepths')
points(sort(x2[,3]), col=2)
points(sort(x3[,3])-.1, col=3)
admb.samples <- extract_samples(admb)
tmb.samples <- extract_samples(tmb)
stan.samples <- data.frame(extract(stan2, permuted=FALSE)[,1,])
qqplot(stan.samples[,1], admb.samples[,1])
x <- qqplot(stan.samples[,1], tmb.samples[,1], plot=FALSE)
points(x[[1]], x[[2]], col=2)


ggplot(nevals, aes(treedepth__, n_leapfrog__, color=platform)) + geom_jitter()
sum(!2^(x1$treedepth__-1)-1 < x1$n_leapfrog__ &
 x1$n_leapfrog__ <= 2^(x1$treedepth__)-1)
sum(!2^(x2$treedepth__-1)-1 < x2$n_leapfrog__ &
 x2$n_leapfrog__ <= 2^(x2$treedepth__)-1)
sum(!2^(x3$treedepth__-1)-1 < x3$n_leapfrog__ &
 x3$n_leapfrog__ <= 2^(x3$treedepth__)-1)

## launch_shinystan_admb(admb)
## launch_shinystan_admb(tmb)
## launch_shinystan(stan2)


chains <- 6
init <- lapply(1:chains, function(x) rnorm(2, sd=se*1))
iter <- 500
tmb1 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=diag(2)))
tmb2 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=covar2))
tmb3 <- sample_tmb(mvnd.obj, iter=iter, init=init, chains=chains, control=list(metric=covar))
adapt.tmb1 <- cbind(form=1, extract_adapt(tmb1))
adapt.tmb2 <- cbind(form=2, extract_adapt(tmb2))
adapt.tmb3 <- cbind(form=3, extract_adapt(tmb3))
adapt.tmb <- rbind(adapt.tmb1, adapt.tmb2, adapt.tmb3)
admb1 <- sample_admb(dir=dir, model, iter=iter, init=init, chains=chains, control=list(metric=diag(2)))
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
ggplot(subset(zz, iter>5), aes(iter, log10( stepsize__), group=seed, color=m)) +
  facet_grid(covar~., scales='free') + geom_point(alpha=.15, size=.5)

adapt <- rbind(cbind(name='tmb', adapt.tmb),
               cbind(name='admb', adapt.admb))
adapt.long <- reshape2::melt(adapt, c("name", "form", "chain"))
ggplot(adapt.long, aes(form, value, color=name)) + geom_jitter() +
  facet_wrap('variable', scales='free')
ggplot(adapt, aes(form, stepsize, color=name)) + geom_jitter()
ggplot(adapt, aes(form, nsteps.median, color=name)) + geom_jitter()
ggplot(adapt, aes(form, accept.prob, color=name)) + geom_jitter()
ggplot(adapt, aes(form, log(stepsize0), color=name)) + geom_jitter()
launch_shinystan_tmb(tmb1)
## launch_shinystan_tmb(tmb2)
## launch_shinystan_tmb(tmb3)
launch_shinystan_admb(admb1)
launch_shinystan_admb(admb2)
## launch_shinystan_admb(admb3)



dir="C:/Users/Cole/hmc_tests/models/catage"
model='catage'
x <- sample_admb(dir=dir, model=model, iter=20000,
                   chains=1, eps=.2, max_treedepth=10, thin=100)
x <- sample_admb(dir=dir, model=model, iter=20000,
                   chains=1, eps=.2, max_treedepth=10, thin=100)
launch_shinystan_admb(x)
setwd(dir)
system('admb catage')
system('catage -nohess -mcmc 10 -nuts -mcseed 5')
adapt <- read.csv("adaptation.csv")
pars <- read_psv(model)
fit <- read_admb(model)
