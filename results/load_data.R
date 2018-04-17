## This file loads in the results from all models, empirical and simulated,
## and prepares them for plotting. DFs are returned in the global
## workspace, and two .csv files empirical.csv and simulated.csv are
## written to the results folder.

empirical <- ldply(list.files('results', pattern='perf_empirical'), function(i)
    read.csv(paste0('results/',i)))
## reorder platforms
empirical$platform <- factor(empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
## Drop "tmb" platform which is sample_tmb from adnuts package.
#empirical <- droplevels(subset(empirical, platform != 'tmb'))

## normalize by maximum run time across delta.target values
empirical <-
    ddply(empirical, .(platform, model), mutate,
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

## Make table of relative efficiencies (relative to stan). This table is in
## the paper.
rnd <- 3
empirical.perf.table <- within(empirical.median.wide, {
  stan.re <- round(stan/stan,rnd)
  tmbstan.re <- round(tmbstan/stan,rnd)
  admb.re <- round(admb/stan,rnd)})
empirical.perf.table <- empirical.perf.table[, c(1,2, 8,7,6)]
print(empirical.perf.table)

## Adaptation
adapt_empirical<- ldply(list.files('results', pattern='adapt_empirical'), function(i)
    read.csv(paste0('results/',i)))
adapt_empirical$platform <- factor(adapt_empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
adapt_empirical <- droplevels(subset(adapt_empirical, !model%in%c('iidz', 'zdiag')))

## Write file to use for making figures and  make_plots.R
write.table(empirical, file='results/empirical.csv', sep=',',
            row.names=FALSE, col.names=TRUE)


### ---------- Simulated results
simulated <- ldply(list.files('results', pattern='perf_simulated'), function(i)
    read.csv(paste0('results/',i)))
simulated$platform <- factor(simulated$platform, levels= c("stan", "tmbstan", "admb", "tmb"))
## Drop "tmb" platform which is sample_tmb from adnuts package.
simulated <- droplevels(subset(simulated, platform != 'tmb'))
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
temp <- ddply(empirical, .(model, platform), summarize, reps=length(seed))
print(dcast(temp, model~platform, value.var='reps'))
#ddply(simulated, .(model, platform, Npar), summarize, reps=length(seed))
