Performance tests for no-U-turn sampling for TMB and ADMB vs. Stan 
=========

This repository contains code, data, and models to recreate the analysis in
the paper:

Monnahan CC, Kristensen K (2018) No-U-turn sampling for fast Bayesian
inference in ADMB and TMB: Introducing the adnuts and tmbstan R
packages. PLoS ONE 13(5): e0197954. https://doi.org/10.1371/journal.pone.0197954

The repository evolved from some previous testing code which is left for
now and should be ignored. The plots folder contains extra performance
plots, the models folder the models and data, and the results folder saved
output from running the R scripts.

The paper has three components with scripts to recreate them:

Demonstration of features
-----

This can recreated by executing the file run_demo.R and shows basic
functionality of the two packages using two real models.

Laplace approximation checks
------

The paper includes a component where a TMB model is run with and without
the Laplace approximation turned on, in order to test the accuracy of this
on the fixed effects. The file run_laplace.R recreates this analysis. Note
that this takes a long time to run.

Performance testing
-------
The supplementary material has a more thorough exploration of peformance
between the three software platforms: Stan, TMB, and ADMB, using their
respective R packages (rstan, tmbstan, adnuts). This is run on two
simulated models with increasing dimensionality (zdiag and growth) as well
as three "real" models. The file run_analysis.R recreates this part. It
takes a long time and is not setup to run in parallel. 

