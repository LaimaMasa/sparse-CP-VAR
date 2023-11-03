# Sparse change point vector autoregression, or sparse CP-VAR(p)

## Overview

The sparse CP-VAR(p) model is developed using several sophisticated statistical building blocks:
- VAR(p) model large multivariate autoregression model denoted by parameters n - number of time series, and p - number of lags in the time series. The total number of parameters in such a model is (n × p + 1) × n,
- Sparsity, or parameter penalization technique that shrinks irrelevant parameter values to zero,
- Change or break points, which identify the point in time at which the parameters of the model change drastically,
- MCMC with Gibbs Sampling is an approach to perform Monte Carlo Markov Chain sampling using a special case of the Metropolis-Hastings algorithm, where we sample parameter values from fully conditional distributions.

This model is implemented with: 
- R software,
- JAGS (Just Another Gibbs Sampling), an open-source statistical inference tool for hierarchical Bayesian methods, accessed via the R package runjags.

The original model was developed in [Dufays et al. 2021]. The coded implementation provided here is part of my master's thesis, and contains several deviations from the original model.
