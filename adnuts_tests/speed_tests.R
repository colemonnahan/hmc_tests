## Some code to profile and try to improve TMB speed.
rm(list=ls())
devtools::install('C:/Users/Cole/adnuts')

library(TMB)
library(adnuts)
TMB::runExample("simple")
## Remake obj without random components
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
                 random=NULL, DLL="simple", silent=TRUE)
## init can be a function, or a list of lists, or NULL to use starting
## values in obj.
init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)

## Note that this model does not use any box constraints. These can be
## passed to sample_tmb just like with nlminb. Also, there are no explicit
## priors in this model, which can be influential particularly for the
## logsd params.
fn <- obj$fn; gr <- obj$gr
y.cur <- unlist(init())
## Make sure these match whether mass matrix is vector or matrix
M1 <- diag(abs(rnorm(length(unlist(init())))))
r1 <- adnuts:::rotate_space(fn=fn, gr=gr, M=M1, y.cur=y.cur)
M2 <- diag(M1)
r2 <- adnuts:::rotate_space(fn=fn, gr=gr, M=M2, y.cur=y.cur)
head(r1$gr2(r1$x.cur))
head(r2$gr2(r2$x.cur))
r1$fn2(r1$x.cur)
r2$fn2(r2$x.cur)
str(diag(r1$chd))
str(r2$chd)
x <- r1$x.cur
library(microbenchmark)
op <- microbenchmark(fn=fn(x), fn1=r1$fn(x), fn2=r2$fn2(x),
               gr=gr(x), gr1=r1$gr(x), gr2=r2$gr2(x),
               times=5000L)
plot(op, log='y')

## Now profile over NUTS
Rprof()
tmb <- sample_tmb(obj=obj, iter=500000, init=init, chains=3, thin=1000,
                          algorithm='RWM')
Rprof(NULL)
summaryRprof()$by.self[1:15,]



## library(microbenchmark)
## library(plyr)
## op <- ldply(c(2, 10, 50, 100, 200, 500, 1000, 2000, 5000),
##             function(i) {
##   x <- rnorm(i)
##   fn <- function(x) sum(x^2)
##   gr <- function(x) -2*x
##   M <- diag(i)
##   chd <- t(chol(M))               # lower triangular Cholesky decomp.
##   ## Redefine these functions
##   fn2 <- function(theta) fn(chd %*% theta)
##   gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
##   M <- rep(1,i)
##   ## Redefine these functions
##   fn3 <- function(theta) fn(M * theta)
##   gr3 <- function(theta) gr(M * theta)
##   cbind(i=i, microbenchmark(
##     fn=fn(x), fn2=fn2(x), fn3=fn3(x),
##     times=500L))})
## library(ggplot2) #nice log plot of the output
## op2 <- ddply(op, .(expr, i), summarize, median=median(time))
## qplot(x=log10(i), y=median, data=op2, colour=factor(expr), geom='line') + scale_y_log10()
