# test that run_paramedic works!
# set up the data
library("testthat")
library("rstan")
library("dplyr")
data(example_16S_data)
data(example_qPCR_data)

# hyperparameter values
sigma_beta <- 5
sigma_Sigma <- 5
alpha_sigma <- 2
kappa_sigma <- 1
alpha_phi <- 0
beta_phi <- 0

# run paramedic (with a *very* small number of iterations, for illustration only)
# also only on the first 10 taxa
test_that("paramedic works", {
  expect_warning(mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
                                                 sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                 alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                 alpha_phi = alpha_phi, beta_phi = beta_phi,
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(adapt_delta = 0.9, max_treedepth = 15)))
  # get model summary
  mod_summ <- summary(mod, probs = c(0.025, 0.975))$summary
  # check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), tolerance = 1, 
               scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})
# test that paramedic works with covariates
X <- cbind(example_qPCR_data[, 1], rbinom(dim(example_qPCR_data)[1], 1, prob = 0.6))
test_that("paramedic works with covariates", {
  expect_warning(mod <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data, X = X,
                                                 sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                 alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                 alpha_phi = alpha_phi, beta_phi = beta_phi,
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(adapt_delta = 0.9, max_treedepth = 15)))
  # get model summary
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  # check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})

test_that("centered-paramedic works", {
  expect_warning(mod_centered <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data,
                                                                   sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                                   alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                          alpha_phi = alpha_phi, beta_phi = beta_phi,
                                                          centered = TRUE,
                                                    n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                                    control = list(adapt_delta = 0.9, max_treedepth = 15)))
  mod_summ_c <- rstan::summary(mod_centered, probs = c(0.025, 0.975))$summary
  expect_equal(mean(mod_summ_c[grepl("mu", rownames(mod_summ_c)) & grepl(",1]", rownames(mod_summ_c)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})

test_that("centered-paramedic works with covariates", {
  expect_warning(mod_centered <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data, X = X,
                                                                   sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                                   alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                          alpha_phi = alpha_phi, beta_phi = beta_phi, centered = TRUE,
                                                    n_iter = 50, n_burnin = 30, n_chains = 1, stan_seed = 4747,
                                                    control = list(adapt_delta = 0.9, max_treedepth = 15)))
  mod_summ_c <- rstan::summary(mod_centered, probs = c(0.025, 0.975))$summary
  expect_equal(mean(mod_summ_c[grepl("mu", rownames(mod_summ_c)) & grepl(",1]", rownames(mod_summ_c)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})

test_that("negative-binomial paramedic works", {
  expect_warning(mod_negbin <- paramedic::run_paramedic(W = example_16S_data[, 1:10], V = example_qPCR_data, X = X,
                                                        sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                        alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                        alpha_phi = 3, beta_phi = 10,
                                                        n_iter = 50, n_burnin = 30, n_chains = 1, stan_seed = 4747,
                                                        control = list(adapt_delta = 0.9, max_treedepth = 15)))
  mod_summ_nb <- rstan::summary(mod_negbin, probs = c(0.025, 0.975))$summary
  expect_equal(mean(mod_summ_nb[grepl("mu", rownames(mod_summ_nb)) & grepl(",1]", rownames(mod_summ_nb)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})

# generate technical replicates
K <- 2
set.seed(12345)
W_diff <- replicate(ncol(example_16S_data) - 1, runif(n = nrow(example_16S_data), min = 0, max = 100))
V_diff <- replicate(ncol(example_qPCR_data) - 1, runif(n = nrow(example_qPCR_data), min = 0, max = 100))
W_replicate <- example_16S_data
W_replicate[, -1] <- W_replicate[, -1] + as.integer(W_diff)
W <- bind_rows(example_16S_data %>% mutate(batch = 1), W_replicate %>% mutate(batch = 2))
V_replicate <- example_qPCR_data
V_replicate[, -1] <- V_replicate[, -1] + as.integer(V_diff)
V <- bind_rows(example_qPCR_data %>% mutate(batch = 1), V_replicate %>% mutate(batch = 2))
test_that("paramedic with batches works", {
  expect_warning(mod_batches <- paramedic::run_paramedic(W = W[, c(1:10, 435)], V = V, k = K,
                                                         sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                         alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                         n_iter = 50, n_burnin = 30, n_chains = 1, stan_seed = 4747,
                                                         control = list(adapt_delta = 0.9, max_treedepth = 15)))
  mod_summ_batch <- rstan::summary(mod_batches, probs = c(0.025, 0.975))$summary
  e_mod <- mod_summ_batch[grepl("e", rownames(mod_summ_batch)) & !grepl("beta", rownames(mod_summ_batch)), ]
  expect_equal(nrow(e_mod), K * (ncol(W[, c(1:10, 435)]) - 2) + 1, 
               tolerance = 1)
})


# check to see that errors/warnings work for run_paramedic
test_that("errors and warnings for run_paramedic work", {
  # check to make sure that if W and V have different numbers of rows, we stop
  expect_error(paramedic::run_paramedic(W = example_16S_data[1, 1:10], V = example_qPCR_data,
                                        sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                        alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                        alpha_phi = alpha_phi, beta_phi = beta_phi,
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  # expect error if q < q_obs
  expect_error(paramedic::run_paramedic(W = example_16S_data[1, 1:3], V = example_qPCR_data,
                                        sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                        alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                        alpha_phi = alpha_phi, beta_phi = beta_phi,
                                        n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                        control = list(max_treedepth = 15)))
  # check to make sure that both W and V have as first column the sample IDs
  wrong_sample_id_name <- example_16S_data %>%
    mutate(SampleID = sample_id) %>%
    select(SampleID, names(example_16S_data)[2:dim(example_16S_data)[2]], -sample_id)
  expect_error(paramedic::run_paramedic(W = wrong_sample_id_name[, 1:10], V = example_qPCR_data,
                                        sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                        alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                        alpha_phi = alpha_phi, beta_phi = beta_phi,
                                        n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                        control = list(max_treedepth = 15)))
  # check what happens if rows are scrambled
  scrambled_rows <- example_16S_data %>%
    arrange(desc(sample_id))
  expect_warning(paramedic::run_paramedic(W = scrambled_rows[, 1:10], V = example_qPCR_data,
                                          sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                          alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                          alpha_phi = alpha_phi, beta_phi = beta_phi,
                                          n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                          control = list(max_treedepth = 15)))
  # check what happens if columns are scrambled
  scrambled_cols <- example_16S_data %>%
    select(sample_id, Lactobacillus.iners, Gardnerella.vaginalis, Lactobacillus.crispatus,
           Lactobacillus.jensenii:Lactobacillus.gasseri)
  expect_warning(paramedic::run_paramedic(W = scrambled_cols, V = example_qPCR_data,
                                          sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                          alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                          alpha_phi = alpha_phi, beta_phi = beta_phi,
                                          n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                          control = list(max_treedepth = 15)))
  # make sure that code with rows and columns scrambled works
  scrambled_rows_and_cols <- example_16S_data %>%
    arrange(desc(sample_id)) %>%
    select(sample_id, Lactobacillus.iners, Gardnerella.vaginalis, Lactobacillus.crispatus,
           Lactobacillus.jensenii:Lactobacillus.gasseri)
  expect_warning(mod <- paramedic::run_paramedic(W = scrambled_rows_and_cols, V = example_qPCR_data,
                                                 sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma,
                                                 alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma,
                                                 alpha_phi = alpha_phi, beta_phi = beta_phi,
                                  n_iter = 35, n_burnin = 25, n_chains = 1, stan_seed = 4747,
                                  control = list(max_treedepth = 15)))
  mod_summ <- rstan::summary(mod, probs = c(0.025, 0.975))$summary
  expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
               mean(example_qPCR_data$Gardnerella.vaginalis), 
               tolerance = 1, scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})