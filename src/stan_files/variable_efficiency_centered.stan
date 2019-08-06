data{
    int<lower=0> N;
    int<lower=0> q_obs;
    int<lower=0> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
}
parameters{
    vector[q] log_mu[N];
    vector[q] beta;
    row_vector<lower=0>[q] Sigma;
    vector<lower=0>[q] e;
    real<lower=0> sigma_e;
}
transformed parameters{
    vector<lower=0>[q] mu[N];
    vector<lower=0>[q] p[N];
    mu = exp(log_mu);
    for (j in 1:N){
        p[j] = (e .* mu[j])/sum(e .* mu[j]);
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

    // model
    beta ~ normal(0, sigma_beta);
    Sigma ~ lognormal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    e ~ lognormal(0, sqrt(sigma_e));

    for (j in 1:N){
        log_mu[j] ~ normal(beta, Sigma);
        V[j] ~ poisson(mu[j,1:q_obs]);
        W[j] ~ multinomial(p[j]);
    }
}
