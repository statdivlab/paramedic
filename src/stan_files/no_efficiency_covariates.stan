data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=0> d;
    int<lower=0> V[N,q_obs];
    int<lower=0> W[N,q];
    matrix[N,d] X;
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
}
parameters{
    // first-level parameters
    vector[q] log_mu_tilde[N];
    // second-level hyperparameters
    vector[q] beta_0;
    matrix[d,q] beta_1;
    vector[q] log_Sigma;
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
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

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
