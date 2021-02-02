# test that predict_paramedic works!
# set up the data
library("testthat")
library("rstan")
library("dplyr")
data(example_16S_data)
data(example_qPCR_data)

# set up folds
set.seed(1234)
folds <- sample(1:2, size = nrow(example_16S_data), 
                prob = c(0.5, 0.5), replace = TRUE)

# hyperparameter values
sigma_beta <- 5
sigma_Sigma <- 5
alpha_sigma <- 2
kappa_sigma <- 1
alpha_phi <- 0
beta_phi <- 0

test_that("predictions from paramedic work (using posterior_predict.paramedic)", {
    set.seed(4747)
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
    set.seed(1234)
    pp <- posterior_predict(mod, 
                            W = example_16S_data[folds == 2, 1:10], 
                            V = example_qPCR_data[folds == 2, ],
                            alpha_sigma = alpha_sigma, kappa_sigma = kappa_sigma, 
                            alpha_phi = alpha_phi, beta_phi = beta_phi)
    # average over draws
    V_means <- apply(pp$V, c(1, 2), mean)
    expect_equal(colMeans(V_means)[1], 
                 mean(example_qPCR_data$Gardnerella.vaginalis[folds == 2]),
                 tolerance = 1, 
                 scale = mean(example_qPCR_data$Gardnerella.vaginalis[folds == 2]))
})

# run paramedic (with a *very* small number of iterations, for illustration only)
# also only on the first 10 taxa
test_that("predictions from paramedic work (using new stan model)", {
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
    # get posterior distributions
    ext_mod <- rstan::extract(mod)
    beta0_post <- ext_mod$beta_0
    sigma_post <- ext_mod$Sigma
    sigmae_post <- ext_mod$sigma_e
    # obtain predictions
    pred_mod <- paramedic::predict_paramedic(
            W = example_16S_data[folds == 2, 1:10], 
            V = example_qPCR_data[folds == 2, ], 
            beta_0 = beta0_post, Sigma = sigma_post, sigma_e = sigmae_post,
            n_iter = 100, 
            n_burnin = 25, n_chains = 1, stan_seed = 4747
    )
    # check accuracy on V; note that 
    # dimension 1: number of iterations (new)
    # dimension 2: number of samples from posterior on old dataset
    # dimension 3: sample size of new dataset
    # dimension 4: number of taxa
    ext_pred <- rstan::extract(pred_mod)
    V_pred <- apply(ext_pred$V, c(3, 4), mean)
    all_v_means <- colMeans(V_pred)
    expect_equal(all_v_means[1], 
                 mean(example_qPCR_data[folds == 2, ]$Gardnerella.vaginalis), 
                 tolerance = 1, 
                 scale = all_v_means[1])
})
