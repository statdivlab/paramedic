context("Test the no-efficiency models (only present to replicate results from the paramedic manuscript)")

## test that no_efficiency works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
library("dplyr")
data(example_16S_data)
data(example_qPCR_data)

# hyperparameter values
sigma_beta <- 5
sigma_Sigma <- 5
alpha_sigma <- 2
kappa_sigma <- 1

## run no_efficiency (with a *very* small number of iterations, for illustration only)
## also only on the first 10 taxa
test_that("no-efficiency works", {
  expect_warning(mod <- paramedic::no_efficiency(W = example_16S_data[, 1:10], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  ## get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 0.5, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})
# test that no_efficiency works with covariates
X <- cbind(example_qPCR_data[, 1], rbinom(dim(example_qPCR_data)[1], 1, prob = 0.6))
test_that("no-efficiency works with covariates", {
  expect_warning(mod <- paramedic::no_efficiency(W = example_16S_data[, 1:10], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma, X = X,
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  ## get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})

## check to see that errors/warnings work for no_efficiency
test_that("errors and warnings for no_efficiency work", {
  ## check to make sure that if W and V have different numbers of rows, we stop
  expect_error(paramedic::no_efficiency(W = example_16S_data[1, 1:10], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  ## expect error if q < q_obs
  expect_error(paramedic::no_efficiency(W = example_16S_data[1, 1:3], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                        n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                        control = list(max_treedepth = 15)))
  ## check to make sure that both W and V have as first column the sample IDs
  wrong_sample_id_name <- example_16S_data %>%
    mutate(SampleID = sample_id) %>%
    select(SampleID, names(example_16S_data)[2:dim(example_16S_data)[2]], -sample_id)
  expect_error(paramedic::no_efficiency(W = wrong_sample_id_name[, 1:10], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                        n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                        control = list(max_treedepth = 15)))
  ## check what happens if rows are scrambled
  scrambled_rows <- example_16S_data %>%
    arrange(desc(sample_id))
  expect_warning(paramedic::no_efficiency(W = scrambled_rows[, 1:10], V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                          n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                          control = list(max_treedepth = 15)))
  ## check what happens if columns are scrambled
  scrambled_cols <- example_16S_data %>%
    select(sample_id, Lactobacillus.iners, Gardnerella.vaginalis, Lactobacillus.crispatus,
           Lactobacillus.jensenii:Lactobacillus.gasseri)
  expect_warning(paramedic::no_efficiency(W = scrambled_cols, V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                          n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                          control = list(max_treedepth = 15)))
  ## make sure that code with rows and columns scrambled works
  scrambled_rows_and_cols <- example_16S_data %>%
    arrange(desc(sample_id)) %>%
    select(sample_id, Lactobacillus.iners, Gardnerella.vaginalis, Lactobacillus.crispatus,
           Lactobacillus.jensenii:Lactobacillus.gasseri)
  expect_warning(mod <- paramedic::no_efficiency(W = scrambled_rows_and_cols, V = example_qPCR_data, sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,                                                  
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 0.5, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})
