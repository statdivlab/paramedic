## test that predict_paramedic works!
## set up the data
library("testthat")
library("rstan")
library("paramedic")
library("dplyr")
data(example_16S_data)
data(example_qPCR_data)

# set up folds
set.seed(1234)
folds <- sample(1:2, size = nrow(example_16S_data), replace = TRUE)

# hyperparameter values
sigma_beta <- 5
sigma_Sigma <- 5
alpha_sigma <- 2
kappa_sigma <- 1
alpha_phi <- 0
beta_phi <- 0

## run paramedic (with a *very* small number of iterations, for illustration only)
## also only on the first 10 taxa
test_that("predictions from paramedic work", {
    expect_warning(
        mod <- paramedic::run_paramedic(
            W = example_16S_data[folds == 1, 1:10], 
            V = example_qPCR_data[folds == 1, ], 
            sigma_beta = sigma_beta, sigma_Sigma = sigma_Sigma, 
            alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma, 
            alpha_phi = alpha_phi, beta_phi = beta_phi, n_iter = 35, 
            n_burnin = 25, n_chains = 1, stan_seed = 4747, 
            control = list(adapt_delta = 0.9, max_treedepth = 15)
            )
        )
    ## get posterior distributions
    ext_mod <- rstan::extract(mod)
    beta0_post <- ext_mod$beta_0
    sigma_post <- ext_mod$Sigma
    sigmae_post <- ext_mod$sigma_e
    ## obtain predictions
    expect_warning(
        pred_mod <- paramedic::predict_paramedic(
            W = example_16S_data[folds == 2, 1:10], 
            V = example_qPCR_data[folds == 2, ], 
            beta_0 = beta0_post, Sigma = sigma_post, sigma_e = sigmae_post,
            n_iter = 100, 
            n_burnin = 25, n_chains = 1, stan_seed = 4747
        )
    )
    ## check that mean of mu for taxon 1 is "close" to mean qPCR for taxon 1
    expect_equal(mean(mod_summ[grepl("mu", rownames(mod_summ)) & grepl(",1]", rownames(mod_summ)), 1]), 
                 mean(example_qPCR_data$Gardnerella.vaginalis), tolerance = 1, 
                 scale = mean(example_qPCR_data$Gardnerella.vaginalis))
})