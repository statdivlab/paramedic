## test that run_paramedic works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
data(simple_example_data)
q <- 3
q_obs <- 2
## process the data
processed_data <- paramedic::process_data(full_data = simple_example_data, rel_inds = 1:q,
                               abs_inds = (q + 1):(q + q_obs),
                               abs_plus_rel_inds = 1:q_obs,
                               regex_thr = "", regex_abs = "_cps", llod = 0,
                               m_min = 1000, div_num = 1)
## run paramedic (with a *very* small number of iterations, for illustration only)
test_that("paramedic works", {
  mod <- paramedic::run_paramedic(W = processed_data$br, V = processed_data$qpcr,
                              stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                              params_to_save = c("mu", "Sigma", "beta", "e"))
  ## get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  ## check that dimensions are correct
  expect_equal(dim(mod_summ), c(34, 7))
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), mean(processed_data$qpcr[, 1]), tolerance = 0.3)
})
