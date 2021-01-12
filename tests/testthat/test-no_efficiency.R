# test that no_efficiency works!
# set up the data
library("testthat")
library("rstan")
library("dplyr")
data(example_16S_data)
data(example_qPCR_data)

# hyperparameter values
sigma_beta <- 5
sigma_Sigma <- 5
alpha_sigma <- 0
kappa_sigma <- 0
alpha_phi <- 0
beta_phi <- 0

# run efficiency-naive (with a *very* small number of iterations, for illustration only)
# also only on the first 10 taxa
test_that("no-efficiency works", {
  expect_warning(mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data, 
                                                 sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,  
                                                 alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                 alpha_phi = alpha_phi, beta_phi = beta_phi,
                                                 n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                                 control = list(max_treedepth = 15)))
  # get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  # check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})
# test that efficiency-naive works with covariates
X <- cbind(example_qPCR_data[, 1], rbinom(dim(example_qPCR_data)[1], 1, prob = 0.6))
test_that("no-efficiency works with covariates", {
  expect_warning(mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data, X = X,
                                                 sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,  
                                                 alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                 alpha_phi = alpha_phi, beta_phi = beta_phi,
                                                 n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                                 control = list(max_treedepth = 15)))
  # get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  # check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})