setwd(main.dir)

## Loop through and get the maximum correlation of each model for merging
## into the main results.
## cor.table <- ldply(list.files('models'), function(i) {
##               xx <- readRDS(file.path('models', i, 'sims.ind.RDS'))
##               cortemp <- cor(xx)
##               max.cor <- max(abs(cortemp[lower.tri(cortemp)]))
##               median.cor <- median(abs(cortemp[lower.tri(cortemp)]))
##               data.frame(model=i, max.cor=max.cor, median.cor=median.cor)
##             })
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


empirical <- ldply(list.files('results', pattern='perf_empirical'), function(i)
    read.csv(paste0('results/',i)))
empirical$platform <- factor(empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))

## normalize by maximum run time across delta.target values
empirical <-
    ddply(empirical, .(platform, model, delta.target), mutate,
          mean.efficiency=mean(efficiency),
          sd.efficiency=sd(efficiency),
          median.efficiency=quantile(efficiency, probs=.5),
          lwr.efficiency=min(efficiency),
          upr.efficiency=max(efficiency))
## empirical <- merge(x=empirical, y=cor.table, by='model')
empirical.means <-
    ddply(empirical, .(platform, model, Npar), summarize,
          mean.efficiency=mean(efficiency))
empirical.means.wide <- dcast(subset(empirical.means), model+Npar~platform,
                              value.var='mean.efficiency')
## empirical.means.wide$stan_re <- with(empirical.means.wide, round(stan/tmb, 2))
simulated <- ldply(list.files('results', pattern='perf_simulated'), function(i)
    read.csv(paste0('results/',i)))
simulated <-
    ddply(simulated, .(platform, model, delta.target, Npar), mutate,
          mean.efficiency=mean(efficiency),
          sd.efficiency=sd(efficiency),
          median.efficiency=quantile(efficiency, probs=.5),
          lwr.efficiency=min(efficiency),
          upr.efficiency=max(efficiency))
## ## Select Stan modles with default delta.target level
## growth <-
##     subset(simulated, model %in%
##                c('growth','growth_t','growth_nc','growth_nct'))
## growth$model <- as.character(growth$model)
## growth$centered <- 'centered'
## growth$centered[growth$model %in% c('growth_nc', 'growth_nct')] <- 'noncentered'
## growth$normal <- 'normal'
## growth$normal[growth$model %in% c('growth_t', 'growth_nct')] <- 'student-t'
## growth.means <- ddply(growth, .(platform, model, Npar, normal, centered), summarize,
##                       mean.efficiency=round(mean(efficiency),2))
## growth.means.wide <- dcast(subset(growth.means), model+Npar+centered+normal~platform, value.var='mean.efficiency')
## growth.means.wide$stan_re <- with(growth.means.wide, round(stan.nuts/jags, 2))
mvn <- subset(simulated, model %in% c('mvnd', 'mvnc'))
mvn.means <- ddply(mvn, .(platform, model, Npar, cor), summarize,
                      mean.efficiency=round(mean(efficiency),2))

## summarize adaptation info
adapt_empirical<- ldply(list.files('results', pattern='adapt_empirical'), function(i)
    read.csv(paste0('results/',i)))
adapt_empirical$platform <- factor(adapt_empirical$platform, levels= c("stan", "tmbstan", "admb", "tmb"))

## Write them to file
write.table(empirical, file='results/empirical.csv', sep=',',
            row.names=FALSE, col.names=TRUE)
write.table(simulated, file='results/simulated.csv', sep=',',
            row.names=FALSE, col.names=TRUE)
## write.table(file='results/table_growth.csv', x=growth.means.wide, sep=',',
##             row.names=FALSE, col.names=TRUE)
## write.table(file='results/table_cor.csv', x=cor.table, sep=',',
##             row.names=FALSE, col.names=TRUE)
## write.table(file='results/table_perf.csv', x=empirical.means.wide, sep=',',
##             row.names=FALSE, col.names=TRUE)

## Check that everything worked properly among platforms and models
print(subset(empirical, platform=='tmb' & seed ==seeds[1],
             select=c(model, Nsims)))
print(ddply(empirical, .(model, platform), summarize, reps=length(seeds)))
