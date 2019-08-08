## test that processing data works!
## set up the data
library("testthat")
library("paramedic")

## create some data
## set up the data
q <- 3
q_obs <- 2
n <- 10
## ------------------------------------------
## set up hyperparameters for data generation
## ------------------------------------------
Sigma <- diag(1, nrow = 10, ncol = 10)
m_min <- 10000
m_max <- 100000
set.seed(4747)
beta_init <- rnorm(q, 0, sqrt(50))
## order by abundance
beta <- beta_init[order(beta_init, decreasing = TRUE)]
## create the model parameters
e <- rep(1, q)
m <- sample(m_min, m_max, n)
log_mu <- rnorm(n, beta, diag(Sigma))
mu <- exp(log_mu)
## ------------------------------------------
## create the observed data
## ------------------------------------------
V <- matrix(NA, nrow = n, ncol = q_obs)
for (i in 1:q_obs) {
    V[, i] <- rpois(n, mu[, i])    
}
W <- matrix(NA, nrow = n, ncol = q)
for (i in 1:n) {
    p_i <- e*mu[i, ]/sum(e*mu[i, ])
    W[i, ] <- rmultinom(n = 1, size = m[i], prob = p_i)    
}
full_data <- cbind(W, V)

## process the data
test_that("processing works", {
  processed_data <- paramedic::process_data(full_data = full_data, rel_inds = 1:q,
                                 abs_inds = (q + 1):(q + q_obs),
                                 abs_plus_rel_inds = 1:q_obs,
                                 regex_thr = "", regex_abs = "_cps", llod = 0,
                                 m_min = 1000, div_num = 1)
  expect_equal(dim(processed_data$relative)[2], q)
  expect_equal(dim(processed_data$absolute)[2], q_obs)
})
