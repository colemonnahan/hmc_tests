library(R2admb)
library(shinystan)

N <- 2
covar <- diag(N)
covar <- matrix(rWishart(n=1, df=N, Sigma=diag(N)), nrow=N)
write.table(x=c(N, covar), file='mvn.dat', , row.names=FALSE,
            col.names=FALSE)

getADMBCovariance <- function(){
    ## This function reads in all of the information contained in the
    ## admodel.cov file. Some of this is needed for relaxing the
    ## covariance matrix, and others just need to be recorded and
    ## rewritten to file so ADMB "sees" what it's expecting.
    filename <- file("admodel.cov", "rb")
    on.exit(close(filename))
    num.pars <- readBin(filename, "integer", 1)
    cov.vec <- readBin(filename, "numeric", num.pars^2)
    cov <- matrix(cov.vec, ncol=num.pars, nrow=num.pars)
    hybrid_bounded_flag <- readBin(filename, "integer", 1)
    scale <- readBin(filename, "numeric", num.pars)
    result <- list(num.pars=num.pars, cov=cov,
                   hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
    return(result)
}

getADMBCovariance()

system('mvn -mcmc 1000000 -mcsave 100')
rwm <- read_psv('mvn')
system('mvn -mcmc 10000 -hybrid -hynstep 100 -hyeps .5')
hmc <- read_psv('mvn')

pairs(rwm, pch='.')
pairs(hmc, pch='.')
chains <- array(data=NA, dim=c(nrow(rwm),2, ncol(rwm)))
chains[,1,] <- as.matrix(rwm)
chains[,2,] <- as.matrix(hmc)


launch_shinystan(as.shinystan(chains))
launch_shinystan(as.shinystan(hmc2))
