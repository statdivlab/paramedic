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
    vector[q] log_mu_tilde[N];
    vector[q] beta_0;
    vector[q] log_Sigma;
}
model {
    vector[q] log_mu[N];
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    for (i in 1:N){
        log_mu_tilde[i] ~ std_normal();
        log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
        V[i] ~ poisson_log(head(log_mu[i],q_obs));
        W[i] ~ multinomial(softmax(log_mu[i]));
    }
}
generated quantities{
    vector[q] mu[N];
    vector[q] Sigma;

    for (i in 1:N) {
        mu[i] = exp(beta_0 + exp(log_Sigma) .* log_mu_tilde[i]);
    }
    Sigma = exp(log_Sigma);
}
