## Some code to profile and try to improve TMB speed.

library(TMB)
library(adnuts)
TMB::runExample("simple")
## Remake obj without random components
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
                 random=NULL, DLL="simple", silent=TRUE)

## Note that this model does not use any box constraints. These can be
## passed to sample_tmb just like with nlminb. Also, there are no explicit
## priors in this model, which can be influential particularly for the
## logsd params.

## init can be a function, or a list of lists, or NULL to use starting
## values in obj.
init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)
Rprof()
tmb <- sample_tmb(obj=obj, iter=2000, init=init, chains=1,
                          algorithm='NUTS')
Rprof(NULL)
summaryRprof()$by.self[1:15,]

library(microbenchmark)

op <- ldply(c(2, 10, 50, 100, 200, 500, 1000, 2000, 5000), function(i) {
  x <- rnorm(i)
  fn <- function(x) sum(x^2)
  gr <- function(x) -2*x
  M <- diag(i)
  chd <- t(chol(M))               # lower triangular Cholesky decomp.
  ## Redefine these functions
  fn2 <- function(theta) fn(chd %*% theta)
  gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
  M <- rep(1,i)
  ## Redefine these functions
  fn3 <- function(theta) fn(M * theta)
  gr3 <- function(theta) gr(M * theta)
  cbind(i=i, microbenchmark(
    fn=fn(x), fn2=fn2(x), fn3=fn3(x),
    times=500L))
})

library(ggplot2) #nice log plot of the output
op2 <- ddply(op, .(expr, i), summarize, mean=mean(time))
qplot(x=i, y=mean, data=op2, colour=factor(expr), geom='line') + scale_y_log10()
