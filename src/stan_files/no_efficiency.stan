data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N];
    // second-level hyperparameters
    vector[q] beta_0;
    vector[q] log_Sigma;
}
transformed parameters{
    vector[q] p[N];
    vector[q_obs] log_mu_v[N];
    vector[q] log_mu[N];

    for (i in 1:N) {
        log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
        p[i] = softmax(log_mu[i]);
        log_mu_v[i] = head(log_mu[i], q_obs);
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    for (i in 1:N){
        log_mu_tilde[i] ~ std_normal();
        V[i] ~ poisson_log(log_mu_v[i]);
        W[i] ~ multinomial(p[i]);
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] Sigma;

    mu = exp(log_mu);
    Sigma = exp(log_Sigma);
}
