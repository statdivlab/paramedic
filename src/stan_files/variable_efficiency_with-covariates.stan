data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> p;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,p] X;
}
parameters{
    vector[q] log_mu_tilde[N];
    vector[q] beta_0;
    row_vector[p] beta[q];
    vector<lower=0>[q] Sigma;
    vector<lower=0>[q] e;
    real<lower=0> sigma_e;
}
transformed parameters{
    vector<lower=0>[q] mu[N];
    linear_predictor = beta_0 + X * beta
    for (j in 1:N){
        mu[j] = exp(beta_0 + X[j] * beta + Sigma .* log_mu_tilde[j]);
    }
}
model {
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    real alpha_sigma;
    real kappa_sigma;

    sigma_beta = sqrt(50.0);
    sigma_Sigma = sqrt(50.0);
    alpha_sigma = 2;
    kappa_sigma = 1;

    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    beta_1 ~ normal()
    Sigma ~ lognormal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    e ~ lognormal(0, sqrt(sigma_e));

    for (j in 1:N){
        log_mu_tilde[j] ~ normal(0, 1);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial((e .* mu[j])/sum(e .* mu[j]));
    }
}
