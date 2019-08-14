## test that extracting posterior summaries works!
context("Test extracting posterior summaries")

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
test_that("extracting posterior summaries works", {
  summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, taxa_of_interest = 1:3, mult_num = 1, level = 0.95, interval_type = "wald")
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1; same for 2
  expect_equal(colMeans(summs$estimates)[1], mean(example_qPCR_data$Gardnerella.vaginalis), tolerance = 100)
  expect_equal(colMeans(summs$estimates)[2], mean(example_qPCR_data$Lactobacillus.crispatus), tolerance = 100)
})
