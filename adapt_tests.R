
## Test that step size and mass matrix adaptation is consistent across TMB
## and ADMB. Start with simple MVN model.
set.seed(235)
Npar <- 2
corr <- matrix(.98, nrow=2, ncol=2)
diag(corr) <- 1
se <- c(1, .10)
covar <- corr * (se %o% se)
covar <- diag(2)
covar2 <- diag(x=se^2)
covar.inv <- solve(covar)
samples <- mvtnorm::rmvnorm(n=1e5, sigma=covar)
apply(samples, 2, var)
setwd('C:/Users/Cole/hmc_tests/')
compile(file='models/mvnd/mvnd_tmb.cpp')
dyn.load(dynlib('models/mvnd/mvnd_tmb'))
data <- list(covar=covar, Npar=Npar, x=rep(0, len=Npar))
mvnd.obj <- MakeADFun(data=data, parameters=list(mu=c(0,0)), DLL='mvnd_tmb')
model.path <- 'models/mvnd'
model.name <- 'mvnd'
write.table(x=c(2, covar), file='models/mvnd/mvnd.dat', row.names=FALSE, col.names=FALSE)
setwd('models/mvnd')
system('admb mvnd')
system('mvnd')
setwd('../..')
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
iter <- 1000
chains <- 1
eps <- .5
init <- lapply(1:chains, function(x) rnorm(2, sd=se*1))
admb <- run_admb_mcmc(model.path, model.name, iter=iter, init=init,
                       chains=chains, covar=diag(2), eps=eps)
tmb <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains,
                 covar=diag(2), eps=eps)
par(mfrow=c(2,2))
plot(admb$sampler_params[[1]][,4])
points(tmb$sampler_params[[1]][,4], col=2)
plot(admb$sampler_params[[1]][,1])
points(tmb$sampler_params[[1]][,1], col=2)
plot(admb$sampler_params[[1]][,2])
points(tmb$sampler_params[[1]][,2], col=2)
admb.samples <- extract_samples(admb)
tmb.samples <- extract_samples(tmb)
qqplot(admb.samples[,1], tmb.samples[,1])
## launch_shinystan_admb(admb)
## launch_shinystan_admb(tmb)


chains <- 6
init <- lapply(1:chains, function(x) rnorm(2, sd=se*1))
iter <- 500
tmb1 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=diag(2))
tmb2 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar2)
tmb3 <- run_mcmc(mvnd.obj, iter=iter, init=init, chains=chains, covar=covar)
adapt.tmb1 <- cbind(form=1, extract_adapt(tmb1))
adapt.tmb2 <- cbind(form=2, extract_adapt(tmb2))
adapt.tmb3 <- cbind(form=3, extract_adapt(tmb3))
adapt.tmb <- rbind(adapt.tmb1, adapt.tmb2, adapt.tmb3)
admb1 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=diag(2))
admb2 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=covar2)
admb3 <- run_admb_mcmc(model.path, model.name, iter=iter, init=init, chains=chains, covar=covar)
adapt.admb1 <- cbind(form=1, extract_adapt(admb1))
adapt.admb2 <- cbind(form=2, extract_adapt(admb2))
adapt.admb3 <- cbind(form=3, extract_adapt(admb3))
adapt.admb <- rbind(adapt.admb1, adapt.admb2, adapt.admb3)

yy1 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar=1, m='admb', as.data.frame(admb1$sampler_params[[i]])))
xx1 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar=1, m='tmb', as.data.frame(tmb1$sampler_params[[i]])))
yy2 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar=2, m='admb', as.data.frame(admb2$sampler_params[[i]])))
xx2 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar=2, m='tmb', as.data.frame(tmb2$sampler_params[[i]])))
yy3 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter,  covar=3,m='admb', as.data.frame(admb3$sampler_params[[i]])))
xx3 <- ldply(1:chains, function(i)
  cbind(seed=i, iter=1:iter, covar=3, m='tmb', as.data.frame(tmb3$sampler_params[[i]])))
zz <- rbind(yy1,yy2,yy3, xx1,xx2,xx3)
ggplot(zz, aes(iter, log( stepsize__), group=seed)) +
  facet_grid(m~covar) + geom_point(alpha=.15, size=.5)

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



model.path="C:/Users/Cole/hmc_tests/models/catage"
model.name='catage'
x <- run_admb_mcmc(model.path=model.path, model.name=model.name, iter=20000,
                   chains=1, eps=.2, max_treedepth=10, thin=100)
x <- run_admb_mcmc(model.path=model.path, model.name=model.name, iter=20000,
                   chains=1, eps=.2, max_treedepth=10, thin=100)
launch_shinystan_admb(x)
setwd(model.path)
system('admb catage')
system('catage -nohess -mcmc 10 -nuts -mcseed 5')
adapt <- read.csv("adaptation.csv")
pars <- read_psv(model.name)
fit <- read_admb(model.name)
