## test that quantile intervals work!
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
mod <- paramedic::run_paramedic(W = processed_data$br, V = processed_data$qpcr,
                            stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                            params_to_save = c("mu", "Sigma", "beta", "e"))
## get model summary
mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
## get samples
mod_samps <- rstan::extract(mod)
test_that("quantile intervals work", {
  summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, q = q, taxa_of_interest = 1:3, mult_num = 1, level = 0.95, interval_type = "quantile")
  ## check that I have proper dimension of intervals
  expect_equal(dim(summs$pred_intervals), c(dim(processed_data$br)[1], 2, length(1:3)))
  ## check that the intervals cover the observed qPCR for taxon 1
  expect_equal(mean(summs$pred_intervals[,,1][,1] <= processed_data$qpcr[, 1] & summs$pred_intervals[,,1][,2] >= processed_data$qpcr[, 1]), 0.95, tolerance = 0.1)
})
