data{
    int<lower=1> N;
    int<lower=1> q_obs;
    int<lower=1> q;
    int<lower=1, upper=N*q> ii_obs[N*q_obs];
    int<lower=1, upper=N*q> ii_mis[N*q - N*q_obs];
    int<lower=1, upper=q> jj[N*q];
    int<lower=0> V_obs[N*q_obs];
    int<lower=0> W[N,q];
    // hyperparameters
    real sigma_beta;
    real sigma_Sigma;
}
transformed data {
    int<lower=1> n_obs = N*q_obs;
    int<lower=1> n_mis = N*q - N*q_obs;
    int<lower=0> n_tot = n_obs + n_mis;
}
parameters{
    int V_mis[n_mis];
    vector log_mu_tilde[n_tot];
    vector[q] beta_0;
    vector[q] log_Sigma;
}
transformed parameters{
    int V[n_tot];
    V[ii_obs] = V_obs;
    V[ii_mis] = V_mis;
}
model {
    vector log_mu[n_tot];
    // hierarchical model
    beta_0 ~ normal(0, sigma_beta);
    log_Sigma ~ normal(0, sigma_Sigma);

    log_mu_tilde ~ std_normal();

    for (i in 1:N){
        log_mu[i] = beta_0 + exp(log_Sigma) .* log_mu_tilde[i];
        W[i] ~ multinomial(softmax(log_mu[i]));
    }
    V ~ poisson_log(log_mu)
}
generated quantities{
    vector[q] mu[N];
    vector[q] Sigma;

    for (i in 1:N) {
        mu[i] = exp(beta_0 + exp(log_Sigma) .* log_mu_tilde[i]);
    }
    Sigma = exp(log_Sigma);
}
