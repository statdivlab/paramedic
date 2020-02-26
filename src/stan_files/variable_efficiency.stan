data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
}
parameters{
    vector[q] log_mu_tilde[N];
    vector[q] beta;
    vector<lower=0>[q] Sigma;
    vector<lower=0>[q] e;
    real<lower=0> sigma_e;
}
transformed parameters{
    matrix<lower=0>[N,q] mu;
    mu = exp(beta + log_mu_tilde * Sigma);
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
    beta ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    e ~ lognormal(0, sqrt(sigma_e));

    for (j in 1:N){
        log_mu_tilde[j] ~ normal(0, 1);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial((e .* mu[j])/sum(e .* mu[j]));
    }
}
