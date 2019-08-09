## test that wald intervals work!
context("Test wald intervals")

## set up the data
library("testthat")
library("rstan")
library("paramedic")
data(example_16S_data)
data(example_qPCR_data)

## run paramedic (with a *very* small number of iterations, for illustration only)
mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
                                n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                control = list(max_treedepth = 15))
## get model summary
mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
## get samples
mod_samps <- rstan::extract(mod)
test_that("wald intervals work", {
  summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, taxa_of_interest = 1:3, mult_num = 1, level = 0.95, interval_type = "wald")
  ## check that I have proper dimension of intervals
  expect_equal(dim(summs$pred_intervals), c(dim(example_16S_data)[1], 2, length(1:3)))
  ## check that the intervals cover the observed qPCR for taxon 1
  expect_equal(mean(summs$pred_intervals[,,1][,1] <= example_qPCR_data[, 1] & summs$pred_intervals[,,1][,2] >= example_qPCR_data[, 1]), 0.95, tolerance = 0.1)
})
