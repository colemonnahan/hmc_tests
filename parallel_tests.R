library(snowfall)
## ## devtools::document('C:/Users/Cole/admbtools/')
## devtools::install('C:/Users/Cole/admbtools/')
## library(admbtools)
Ncores <- 2
sfStop()
sfInit( parallel=TRUE, cpus=Ncores )

sfExportAll()
xx <- sfLapply(1:Ncores, sample_admb_parallel, dir='admb', iter=2000,
                           model='mvnd', init=rep(list(c(0,0)),1))
out <- sample_admb_parallel(1, dir, model='mvnd', iter=2000, init=list(rnorm(2)))
out <- sample_admb(chains=2, dir, model='mvnd', iter=2000, init=rep(list(rnorm(2)),2))
