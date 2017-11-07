## This model was doing weird things with adapt_delta=.8. The mass
## adaptation was maybe interacting with step size adapt and going
## haywire. Need to come back and diagnose this issue.
compile(paste0(m, '.cpp'))
dyn.load(paste0(m))
obj.tmb <- MakeADFun(data=data, parameters=inits(), DLL=m)
opt <- nlminb(obj.tmb$par, obj.tmb$fn, obj.tmb$gr)
set.seed(1)
ii <- inits()
ii$sigmayearphi <- 1
ii$a <- rep(1.5,17)
test2 <- tmbstan(obj.tmb, iter=500, warmup=250, init=list(ii), seed=1,
                 chains=1, control=list(adapt_engaged=TRUE))
vars <- apply(as.data.frame(extract(test2)), 2, var)[-178]

test <- sample_tmb(obj.tmb, iter=1000, warmup=900, init=list(ii),
                    chains=1, seeds=10,
                   control=list(metric=NULL, adapt_delta=.8, stepsize=NULL,
                                max_treedepth=7))

p <- extract_samples(test, TRUE,TRUE)
plot(p$sigmayearphi)
plot(p$sigmaphi)
abline(v=c(125,225,425,875))
plot(p$lp__)
with(p, plot(sigmayearphi, lp__))
with(p, plot(sigmaphi, lp__))
with(p, plot(sigmap, lp__))
with(p, plot(exp(sigmayearphi), p[,'a[1]']))

sp <- extract_sampler_params(test, TRUE)
ggplot(sp, aes(iteration, log(stepsize__), color=energy__)) + geom_point()
