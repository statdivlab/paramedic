data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> d;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,d] X;
    // hyperparameters
    real hyper_sigma_beta;
    real hyper_sigma_Sigma;
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N];
    // second-level hyperparameters
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
    // third-level hyperparameters
    vector[q] mu_beta;
    vector[q] sigma_beta;
    vector[q] mu_sigma;
    vector[q] sigma_Sigma;
}
transformed parameters{
    simplex[q] p[N];
    vector[q_obs] log_mu_v[N];
    for (i in 1:N){
        p[i] = softmax(beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i]);
        log_mu_v[i] = head(beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i], q_obs);
    }
}
model {
    // hierarchical model
    mu_beta ~ std_normal();
    mu_sigma ~ std_normal();
    sigma_beta ~ normal(hyper_sigma_beta, 1);
    sigma_Sigma ~ normal(hyper_sigma_Sigma, 1);
    beta_0 ~ normal(mu_beta, exp(sigma_beta));
    log_Sigma ~ normal(mu_sigma, exp(sigma_Sigma));

    for (j in 1:q) {
        beta_1[:,j] ~ std_normal();
    }

    for (i in 1:N){
        log_mu_tilde[i] ~ std_normal();
        V[i] ~ poisson_log(log_mu_v[i]);
        W[i] ~ multinomial(p[i]);
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] Sigma;
    for (i in 1:N) {
        mu[i] = exp(beta_0 + (X[i] * beta_1)' + exp(log_Sigma) .* log_mu_tilde[i]);
    }
    Sigma = exp(log_Sigma);
}
