## test that run_paramedic works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
data(example_16S_data)
data(example_qPCR_data)


## run paramedic (with a *very* small number of iterations, for illustration only)
test_that("paramedic works", {
  mod <- paramedic::run_paramedic(W = example_16S_data, V = example_qPCR_data,
                                  stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  params_to_save = c("mu", "Sigma", "beta", "e"))
  ## get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), mean(example_qPCR_data[, 1]), tolerance = 0.3)
})
