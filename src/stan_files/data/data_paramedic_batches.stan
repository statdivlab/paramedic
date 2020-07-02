    // sample size
    int<lower=1> N;
    // num taxa with observed absolute abundance
    int<lower=1> q_obs;
    // overall num taxa
    int<lower=1> q;
    // num covariates
    int<lower=0> d;
    // num batches
    int<lower=1> K;
    // observed absolute abundance
    int<lower=0> V[N,q_obs];
    // observed rel. abundance
    int<lower=0> W[N,q,K];
    // feature matrix
    matrix[N,d] X;
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    // 0 for both = efficiency-naive model
    // otherwise, fit varying-efficiency model
    real<lower=0> alpha_sigma;
    real<lower=0> kappa_sigma;
    // 0 for both = Poisson model
    // otherwise, fit negative binomial
    real<lower=0> alpha_phi;
    real<lower=0> beta_phi;
    // 0 = no additional studies
    real<lower=0> sigma_xi;
