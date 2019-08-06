## test that extracting posterior summaries works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
data(simple_example_data)
q <- 3
q_obs <- 2
## process the data
processed_data <- paramedic::process_data(full_data = simple_example_data, br_inds = 1:q,
                                          qpcr_inds = (q + 1):(q + q_obs),
                                          pcr_plus_br_inds = 1:q_obs,
                                          regex_thr = "", regex_cps = "_cps", llod = 0,
                                          m_min = 1000, div_num = 1)
## run paramedic (with a *very* small number of iterations, for illustration only)
mod <- paramedic::paramedic(W = processed_data$br, V = processed_data$qpcr, q = q, q_obs = q_obs,
                            stan_model = "src/stan_files/variable_efficiency.stan", n_iter = 30, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                            params_to_save = c("mu", "Sigma", "beta", "e"))
## get model summary
mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
## get samples
mod_samps <- rstan::extract(mod)
test_that("extracting posterior summaries works", {
  summs <- extract_posterior_summaries(stan_mod = mod_summ, stan_samps = mod_samps, q = q, taxa_of_interest = 1:3, mult_num = 1, level = 0.95, interval_type = "wald")
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1; same for 2
  expect_equal(colMeans(summs$estimates)[1], mean(processed_data$qpcr[, 1]), tolerance = 0.3)
  expect_equal(colMeans(summs$estimates)[2], mean(processed_data$qpcr[, 2]), tolerance = 0.3)
})
