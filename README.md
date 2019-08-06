# paramedic: Predicting Absolute and Relative Abundance by Modeling Efficiency to Derive Intervals and Concentrations

[![License: BSD-3](https://img.shields.io/badge/License-BSD--3--Clause-yellow)](https://opensource.org/licenses/BSD-3-Clause)

**Author:** Brian Williamson

------------------------------

## Introduction


------------------------------

## Issues

If you encounter any bugs or have any specific feature requests, please [file an issue](https://github.com/statdivlab/paramedic/issues).

------------------------------

## R installation

You may install a stable release of `paramedic` from GitHub via  [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) by running the following code (you may replace `v0.0.1` with the tag for the specific release you wish to install):

```r
## install.packages("devtools") # only run this line if necessary
devtools::install_github(repo = "statdivlab/paramedic@v0.0.1")
```

You may install a development release of `paramedic` from GitHub via [`devtools`](https://www.rstudio.com/products/rpackages/devtools/) by running the following code:

```r
## install.packages("devtools") # only run this line if necessary
devtools::install_github(repo = "statdivlab/paramedic")
```

------------------------------

## Example

This example shows how to use `paramedic` in a simple setting with simulated data. For more examples and detailed explanation, please see the [vignette](vignettes/introduction_to_paramedic.Rmd).

```r
## load required functions and libraries
library("rstan")
library("paramedic")

## -------------------------------------------------------------
## problem setup
## -------------------------------------------------------------
## set up the data
q <- 3
q_obs <- 2
n <- 10
## ------------------------------------------
## set up hyperparameters for data generation
## ------------------------------------------
Sigma <- diag(1, nrow = 10, ncol = 10)
m_min <- 10000
m_max <- 100000
set.seed(4747)
beta_init <- rnorm(q, 0, sqrt(50))
## order by abundance
beta <- beta_init[order(beta_init, decreasing = TRUE)]
## create the model parameters
e <- rep(1, q)
m <- sample(m_min, m_max, n)
log_mu <- rnorm(n, beta, diag(Sigma))
mu <- exp(log_mu)
## ------------------------------------------
## create the observed data
## ------------------------------------------
V <- matrix(NA, nrow = n, ncol = q_obs)
for (i in 1:q_obs) {
  V[, i] <- rpois(n, mu[, i])    
}
W <- matrix(NA, nrow = n, ncol = q)
for (i in 1:n) {
  p_i <- e*mu[i, ]/sum(e*mu[i, ])
  W[i, ] <- rmultinom(n = 1, size = m[i], prob = p_i)    
}
full_data <- cbind(W, V)

## -------------------------------------------------------------
## preliminary step: process the data
## -------------------------------------------------------------
processed_data <- process_data(full_data = full_data, rel_inds = 1:q,
                                abs_inds = (q + 1):(q + 1 + q_obs),
                                abs_plus_rel_inds = 1:q_obs,
                                regex_thr = NA, regex_abs = "", llod = 0,
                                m_min = 1000, div_num = 1000)
qpcr <- processed_data$absolute
br16s <- processed_data$relative

## -------------------------------------------------------------
## Run the Stan algorithm that models varying efficiency
## -------------------------------------------------------------
## this is a small number of iterations, only for illustration
## also, shows how to use control parameters for rstan::stan
stan_mod <- run_paramedic(W = br16s, V = qpcr,
                      stan_model = "src/stan_files/variable_efficiency.stan",
                      n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                      params_to_save = c("mu", "Sigma", "beta", "e"),
                      control = list(adapt_delta = 0.85, max_treedepth = 15),
                      verbose = FALSE)
stan_mod_summ <- summary(stan_mod, probs = c(0.025, 0.975))$summary
stan_mod_samps <- extract(stan_mod)

## -------------------------------------------------------------
## Extract posterior estimates for the taxon missing abolute abundance data
## -------------------------------------------------------------
posterior_summaries <- extract_posterior_summaries(stan_mod_summ, stan_mod_samps,
                                                   taxa_of_interest = q, mult_num = 1,
                                                   level = 0.95, interval_type = "wald")

## -------------------------------------------------------------
## Print out the posterior mean efficiencies and concentrations
## -------------------------------------------------------------
posterior_summaries$estimates
posterior_summaries$est_efficiency
```
