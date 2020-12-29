    // sample size
    int<lower=1> N;
    // num taxa with observed absolute abundance
    int<lower=1> q_obs;
    // overall num taxa
    int<lower=1> q;
    // read counts
    int<lower=1>[N] M;
    // num covariates
    int<lower=0> d;
    // feature matrix (test data)
    matrix[N,d] X;
    // posterior distributions on sigma_e
    real[q] sigma_e;
    // posterior distributions on beta_0, beta_1, Sigma
    real[q] beta_0;
    matrix[d,q] beta_1;
    real[q] Sigma;
    // 0 for both = efficiency-naive model
    // otherwise, fit varying-efficiency model
    real<lower=0> alpha_sigma;
    real<lower=0> kappa_sigma;
    // 0 for both = Poisson model
    // otherwise, fit negative binomial
    real<lower=0> alpha_phi;
    real<lower=0> beta_phi;
