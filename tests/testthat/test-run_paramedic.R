## test that run_paramedic works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
data(example_16S_data)
data(example_qPCR_data)


## run paramedic (with a *very* small number of iterations, for illustration only)
## also only on the first 10 taxa
test_that("paramedic works", {
  mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
                                  n_iter = 200, n_burnin = 150, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15))
  mod_centered <- paramedic::run_paramedic_centered(W = example_16S_data[, 1:10], V = example_qPCR_data,
                                                    n_iter = 200, n_burnin = 150, n_chains = 1, stan_seed = 4747,
                                                    control = list(max_treedepth = 15))
  ## get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  mod_summ_c <- rstan::summary(mod_centered, probs = c(0.025, 0.975))$summary
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), mean(example_qPCR_data[, 1]), tolerance = 100)
  expect_equal(mean(mod_summ_c[grepl("mu", rownames(mod_summ_c)) & grepl(",1]", rownames(mod_summ_c)), 1]), mean(example_qPCR_data[, 1]), tolerance = 100)
})
