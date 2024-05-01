# `R/paramedic`: Predicting Absolute and Relative Abundance by Modeling Efficiency to Derive Intervals and Concentrations <img src="man/figures/paramedic-logo.png" align="right" width="165px"/>

[![Build Status](https://travis-ci.org/statdivlab/paramedic.svg?branch=master)](https://travis-ci.org/statdivlab/paramedic)
[![codecov](https://codecov.io/gh/statdivlab/paramedic/branch/master/graph/badge.svg?token=GnLFG7QNsh)](https://codecov.io/gh/statdivlab/paramedic)
[![License: BSD-3](https://img.shields.io/badge/License-BSD--3--Clause-yellow)](https://opensource.org/licenses/BSD-3-Clause)

**Software author:** [Brian Williamson](https://bdwilliamson.github.io/)

**Methodology authors:** [Brian Williamson](https://bdwilliamson.github.io/), [Jim Hughes](http://faculty.washington.edu/jphughes/) and [Amy Willis](http://statisticaldiversitylab.com/team/amy-willis)

------------------------------

## Introduction

`paramedic` is a R package for estimating microbial concentration. `paramedic` uses information from 16S count data (compositional data on all taxa) and absolute data on a subset of taxa (e.g., qPCR or flow cytometry) to estimate the absolute abundance of all taxa. The method accounts for differing taxon detection efficiencies between the two methods, and produces prediction and confidence intervals as well as point estimates of the absolute abundance of all taxa. Check out [the paper](https://onlinelibrary.wiley.com/doi/10.1111/biom.13503) (also available [here](https://pubmed.ncbi.nlm.nih.gov/34048057/)) for more details.

------------------------------

## R installation

You may install a stable release of `paramedic` from GitHub via  [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) by running the following code (replace `v0.1.0` with the tag for the specific release you wish to install, and only install `devtools` if needed):

```r
install.packages("devtools")
devtools::install_github(repo = "statdivlab/paramedic@v0.1.0")
```

You may install a development release of `paramedic` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) by running the following code:

```r
install.packages("devtools")
devtools::install_github(repo = "statdivlab/paramedic")
```

------------------------------

## Example

This example shows how to use `paramedic` in a simple setting with data inspired by the University of Washington (UW) Sexually Transmitted Infections Cooperative Research Center (STICRC). For more examples and detailed explanation, please see the [vignette](https://statdivlab.github.io/paramedic/articles/introduction_to_paramedic.html).

```r
## load required functions and libraries
library("rstan")
library("paramedic")

## read in the data
data(example_16S_data)
data(example_qPCR_data)

# sample hyperparameter values
sigma_beta <- 1
sigma_Sigma <- 1
alpha_sigma <- 2
kappa_sigma <- 1

## -------------------------------------------------------------
## Estimate concentrations for the first 7 taxa
## -------------------------------------------------------------
## this is a small number of iterations, only for illustration
## also, shows how to use control parameters for rstan::stan
stan_mod <- run_paramedic(W = example_16S_data[, 1:8], V = example_qPCR_data,
                      sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                      alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                      n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                      control = list(adapt_delta = 0.85, max_treedepth = 15),
                      verbose = FALSE)
stan_mod_summ <- summary(stan_mod$stan_fit, probs = c(0.025, 0.975))$summary
stan_mod_samps <- extract(stan_mod$stan_fit)

## -------------------------------------------------------------
## Extract posterior estimates for the taxon missing abolute abundance data
## -------------------------------------------------------------
posterior_summaries <- extract_posterior_summaries(stan_mod_summ, stan_mod_samps,
                                                   taxa_of_interest = 7, mult_num = 1,
                                                   level = 0.95, interval_type = "wald")

## -------------------------------------------------------------
## Print out the posterior mean efficiencies and concentrations
## -------------------------------------------------------------
posterior_summaries$estimates
posterior_summaries$est_efficiency
```

------------------------------

## Citation

After using the `paramedic` package, please cite the following:

```
  @article{williamson2020,
    author={Williamson, BD and Hughes, JP and Willis, AD},
    title={A multi-view model for relative and absolute microbial abundances},
    journal={Biometrics},
    year={2021},
    note={doi.org/10.1111/biom.13503}
  }
```

------------------------------

## Issues

We use Stan to fit the hierarchical model to the data. Read any warning messages returned by the algorithm carefully, as these can help diagnose convergence issues. Carefully choosing initialisation values can speed up convergence.

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/statdivlab/paramedic/issues). If your feature request relates to the hierarchical models implemented in `paramedic`, please read [the hierarchical models vignette](https://statdivlab.github.io/paramedic/articles/hierarchical_model_specification.html) prior to filing an issue or creating a pull request.

------------------------------
