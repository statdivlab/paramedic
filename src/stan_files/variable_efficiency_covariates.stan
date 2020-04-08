data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> p;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,p] X;
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
    real alpha_sigma;
    real kappa_sigma;
}
parameters{
    vector[q] log_mu_tilde[N];
    vector[q] beta_0;
    matrix[p,q] beta_1;
    vector[q] log_Sigma;
    vector[q] log_e;
    real<lower=0> sigma_e;
}
transformed parameters{
    vector[q] log_mu[N];
    for (i in 1:N){
        log_mu[i] = beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i];
    }
}
model {
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    sigma_e ~ inv_gamma(alpha_sigma, kappa_sigma);
    log_e ~ normal(0, sqrt(sigma_e));

    for (j in 1:q) {
        beta_1[:,j] ~ std_normal();
    }

    for (i in 1:N){
        log_mu_tilde[i] ~ std_normal();
        V[i] ~ poisson_log(head(log_mu[i],q_obs));
        W[i] ~ multinomial(softmax(log_e + log_mu[i]));
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] e;
    vector[q] Sigma;
    mu = exp(log_mu);
    e = exp(log_e);
    Sigma = exp(log_Sigma);
}
