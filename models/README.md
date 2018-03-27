Models used in the study.
=================

See the paper for more information and references. Inside each folder is a
Stan model (.stan), TMB model (.cpp), and folder containing the ADMB
version (.tpl). If appropriate, there is also a data.RDS file containing
the real data. This must be read in with R's `loadRDS` function. Most of
them are lists. There is also a "plots" folder that contains exploratory
plots of each model, including performance metrics and verification that
each model is identical (via qqplots).


* _zdiag_ is a normal model with variable variances, but no
   correlations. Data are simulated. Provides testing of adaptation with
   widely different scales on a simple model.

* _growth_ is a Von Bertalannfy non-linear somatic growth model. Data are
   simulated. Provides a scalable non-linear mixed effects model.

* _sslog_ is a logistic state-space fisheries model. Real data.

* _swallows_ is a mark-recpature model that estimates state-space survival
   and detection with environmental covariates and three random effect
   components for a total of 5 fixed effects and 172 random effects. Real
   data.

* _wildf_ (wildflower) is binomial generalized linear mixed effects model
  for flowering success of a perennial plant. 
