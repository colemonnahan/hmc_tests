## This file loads in the results from all models, empirical and simulated,
## and prepares them for plotting. DFs are returned in the global
## workspace, and two .csv files empirical.csv and simulated.csv are
## written to the results folder.
setwd(main.dir)



empirical <- ldply(list.files('results', pattern='perf_empirical'), function(i)
    read.csv(paste0('results/',i)))
empirical$platform <- factor(empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
empirical <- droplevels(subset(empirical, !model%in%c('iidz', 'zdiag')))

## normalize by maximum run time across delta.target values
empirical <-
    ddply(empirical, .(platform, model, delta.target), mutate,
          mean.efficiency=mean(efficiency),
          sd.efficiency=sd(efficiency),
          median.efficiency=quantile(efficiency, probs=.5),
          lwr.efficiency=min(efficiency),
          upr.efficiency=max(efficiency))
## A separate DF for the medians
empirical.median <-
    ddply(empirical, .(model, platform, Npar), summarize,
          median.efficiency=median(efficiency))
empirical.median.wide <- dcast(subset(empirical.median), model+Npar~platform,
                              value.var='median.efficiency')
## Make table of relative efficiencies (relative to stan)
rnd <- 3
empirical.perf.table <- within(empirical.median.wide, {
  stan.re <- round(stan/stan,rnd)
  tmbstan.re <- round(tmbstan/stan,rnd)
  tmb.re <- round(tmb/stan,rnd)
  admb.re <- round(admb/stan,rnd)})
empirical.perf.table <- empirical.perf.table[, c(1, 10,9,8,7)]
print(empirical.perf.table)

## Adaptation
adapt_empirical<- ldply(list.files('results', pattern='adapt_empirical'), function(i)
    read.csv(paste0('results/',i)))
adapt_empirical$platform <- factor(adapt_empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
adapt_empirical <- droplevels(subset(adapt_empirical, !model%in%c('iidz', 'zdiag')))

## Write file to use for making figures
write.table(empirical, file='results/empirical.csv', sep=',',
            row.names=FALSE, col.names=TRUE)


### ---------- Simulated results
simulated <- ldply(list.files('results', pattern='perf_simulated'), function(i)
    read.csv(paste0('results/',i)))
simulated$platform <- factor(simulated$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
simulated <-
    ddply(simulated, .(platform, model, Npar), mutate,
          mean.efficiency=mean(efficiency),
          sd.efficiency=sd(efficiency),
          median.efficiency=quantile(efficiency, probs=.5),
          lwr.efficiency=min(efficiency),
          upr.efficiency=max(efficiency))

## Write them to file
write.table(simulated, file='results/simulated.csv', sep=',',
            row.names=FALSE, col.names=TRUE)

## Check that everything worked properly among platforms and models
temp <- subset(empirical, platform=='tmb' & seed ==1,
             select=c(model, Nsims))
print(temp)
temp <- ddply(empirical, .(model, platform), summarize, reps=length(seed))
print(dcast(temp, model~platform, value.var='reps'))

